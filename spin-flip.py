#!/usr/bin/python3

import os
os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"]      = "8"
os.environ["VECLIB_NUM_THREADS"]   = "8"
os.environ["QULACS_NUM_THREADS"]   = "8"
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from openfermion.chem import MolecularData
from openfermion.transforms import get_fermion_operator, jordan_wigner, bravyi_kitaev
from openfermion.linalg import get_sparse_operator
from openfermion.ops import FermionOperator
from openfermionpyscf import run_pyscf
from pyscf import fci, mcscf
from qulacs import QuantumState
from qulacs import QuantumCircuit
from qulacs.observable import create_observable_from_openfermion_text
from qulacs.gate import H, RX
from itertools import combinations as comb
from itertools import product as prod
from scipy.optimize import minimize
from math import pi
import re
import py3Dmol  # For visualization of molecules and orbitals (don't use in this program)
import pandas as pd  # For numerics
pd.options.display.float_format = "{:,.3f}".format  # {:, .3f} with space causes error!

def create_active_hf_state(mo_occ, nele, nqubits, print_state=None):
    hf_state_bin = ""
    ncore_orb = int(np.sum(mo_occ) - nele) //2
    count_ncore_orb = 0  # spatial orb
    count_nqubits = 0
    for occ in mo_occ:
        if count_nqubits < nqubits +ncore_orb *2:
            count_nqubits += 2
        else:
            break

        if count_ncore_orb < ncore_orb:
            count_ncore_orb += 1
            continue

        if int(occ) == 0:
            hf_state_bin += "00"
        elif int(occ) == 1:
            hf_state_bin += "10"  # considered as alpha spin
        elif int(occ) == 2:
            hf_state_bin += "11"
        else:
            raise ValueError("Invalid occ number detected!!!")

    hf_state_bin = hf_state_bin[::-1]
    hf_state_bin = "0b" + hf_state_bin
    state = QuantumState(nqubits)
    state.set_zero_state()
    state.set_computational_basis(int(hf_state_bin, 0))
    if print_state:
        print(f"HF state in binary = {hf_state_bin}, where 1 is occ, 0 is unocc, and consider this by rotating it 90º counterclockwise!")
        return state, hf_state_bin

    return state


def kron_N(hf_state_bin):
    hf_state_bin = hf_state_bin[2:]    # example: 0b0011 -> 0011
    hf_state_bin = hf_state_bin[::-1]  # example:   0011 -> 1100
    bra0 = np.array([[1,0]])
    bra1 = np.array([[0,1]])
    tmp = bra1 if hf_state_bin[0] == "1" else bra0
    for i in hf_state_bin[1:]:
        op = bra1 if i == "1" else bra0
        tmp = np.kron(tmp, op)
    return tmp


def decode_to_state(eigenvec, threshold=1e-6):
    bit_vec = np.where(eigenvec > threshold, True, False)
    digit = 0
    length = len(eigenvec)
    while length != 1:
        digit += 1
        length //= 2
    op = False
    for indx, bit in enumerate(bit_vec):
        if bit:
            op_str = "      +" if op else "|Psi> ="
            config = str(bin(indx))[2:]
            config = "0" *(digit -len(config)) +config
            print(f"{op_str} {eigenvec[indx]:14.10f} |{config}>")
            op = True

    return None


def decode_to_states(eigenenergies, eigenvecs, states=[0]):
    for istate in states:
        print(f"\n\n=== Electronic Configuration: {istate} ===\n")
        print(f"Energy: {eigenenergies[istate]} [Eh]")
        decode_to_state(eigenvecs[:,istate])

    return None



def create_excitations(hf_state_bin, delta_sz_from_hf_state):
    hf_state_bin = hf_state_bin[2:]    # example: 0b0011 -> 0011
    hf_state_bin = hf_state_bin[::-1]  # example:   0011 -> 1100

    sz_list = np.array([+1/2 if (i %2 == 0) else -1/2 for i in range(len(hf_state_bin))])
    spinocc_list  = np.array([i for i, val in enumerate(hf_state_bin) if val == "1"])
    print(spinocc_list)
    spinvirt_list  = np.array([i for i, val in enumerate(hf_state_bin) if val == "0"])
    print(spinvirt_list)
    singles = [sorted([r, p])  # reorder i > j
               for r, p in list(prod(spinocc_list, spinvirt_list))
               if (sz_list[p] -sz_list[r] == delta_sz_from_hf_state)]

    spinocc_comb  = list(comb(spinocc_list, 2))
    spinvirt_comb = list(comb(spinvirt_list, 2))

    doubles = [sorted([s, r, q, p])  # reorder i > j > k > l
               for (s, r), (q, p) in list(prod(spinocc_comb, spinvirt_comb))
               if (sz_list[p] +sz_list[q] -sz_list[r] -sz_list[s] == delta_sz_from_hf_state)]

    return singles, doubles



def judge_initial_sym(active_orb_sym, hf_state_bin, D2h_table):
    hf_state_bin = hf_state_bin[2:]    # example: 0b0011 -> 0011
    hf_state_bin = hf_state_bin[::-1]  # example:   0011 -> 1100
    initial_sym = "ag"  # just for the temporary symmetry for initailization

    for indx, val in enumerate(hf_state_bin):
        if val == "1":
            initial_sym = D2h_table.loc[initial_sym, active_orb_sym[indx //2]]
    return initial_sym



def judge_D2h_excitation(singles, doubles, active_orb_sym, initial_sym, D2h_table):
    new_singles = []
    new_doubles = []
    for excitation in singles:
        sym = initial_sym
        for indx in excitation:
            sym = D2h_table.loc[sym, active_orb_sym[indx //2]]
        if sym == "ag":  # because target ground state symmetry is "ag"
            new_singles.append(excitation)

    for excitation in doubles:
        sym = initial_sym
        for indx in excitation:
            sym = D2h_table.loc[sym, active_orb_sym[indx //2]]
        if sym == "ag":  # because target ground state symmetry is "ag"
            new_doubles.append(excitation)

    return new_singles, new_doubles



def t1_circuit(param, r, p, pauli_list, cnot_connects, circuit):  # Eq. (1)-(2)
    assert re.search("[^XY]", pauli_list) == None, "pauli matrices X or Y are acceptable."

    gate_r1 = H(r) if pauli_list[0] == "X" else RX(r, pi/2)  # else is "Y"
    gate_p1 = H(p) if pauli_list[1] == "X" else RX(p, pi/2)

    circuit.add_gate(gate_r1)
    circuit.add_gate(gate_p1)

    for cnot_connect in cnot_connects:
        circuit.add_CNOT_gate(cnot_connect[0], cnot_connect[1])  # control, target

    circuit.add_RZ_gate(p, param)

    for cnot_connect in reversed(cnot_connects):
        circuit.add_CNOT_gate(cnot_connect[0], cnot_connect[1])  # control, target

    gate_r2 = H(r) if pauli_list[0] == "X" else RX(r, -pi/2)  # else is "Y"
    gate_p2 = H(p) if pauli_list[1] == "X" else RX(p, -pi/2)

    circuit.add_gate(gate_p2)
    circuit.add_gate(gate_r2)

    return None



def add_fermionic_single(param, single, circuit):  # only 1 parameter
    # Excitation r -> p
    r = single[0]
    p = single[1]

    cnot_connects = [[i, i+1] for i in range(r, p)]

    # paramの定義がpennylaneと異なる
    t1_circuit(-param, r, p, "YX", cnot_connects, circuit)  # Eq. (1)
    t1_circuit( param, r, p, "XY", cnot_connects, circuit)  # Eq. (2)

    return None


from math import pi
from qulacs.gate import H, RX
import re

def t2_circuit(param, s, r, q, p, pauli_list, cnot_connects, circuit):  # Eq. (3)-(10)
    assert re.search("[^XY]", pauli_list) == None, "pauli matrices X or Y are acceptable."

    gate_s1 = H(s) if pauli_list[0] == "X" else RX(s, pi/2)  # else is "Y"
    gate_r1 = H(r) if pauli_list[1] == "X" else RX(r, pi/2)
    gate_q1 = H(q) if pauli_list[2] == "X" else RX(q, pi/2)
    gate_p1 = H(p) if pauli_list[3] == "X" else RX(p, pi/2)

    circuit.add_gate(gate_s1)
    circuit.add_gate(gate_r1)
    circuit.add_gate(gate_q1)
    circuit.add_gate(gate_p1)

    for cnot_connect in cnot_connects:
        circuit.add_CNOT_gate(cnot_connect[0], cnot_connect[1])  # control, target

    circuit.add_RZ_gate(p, param)

    for cnot_connect in reversed(cnot_connects):
        circuit.add_CNOT_gate(cnot_connect[0], cnot_connect[1])  # control, target

    gate_s2 = H(s) if pauli_list[0] == "X" else RX(s, -pi/2)  # else is "Y"
    gate_r2 = H(r) if pauli_list[1] == "X" else RX(r, -pi/2)
    gate_q2 = H(q) if pauli_list[2] == "X" else RX(q, -pi/2)
    gate_p2 = H(p) if pauli_list[3] == "X" else RX(p, -pi/2)

    circuit.add_gate(gate_p2)
    circuit.add_gate(gate_q2)
    circuit.add_gate(gate_r2)
    circuit.add_gate(gate_s2)

    return None



def add_fermionic_double(param, double, circuit):  # only 1 parameter
    # Excitation s,r -> q,p
    s = double[0]
    r = double[1]
    q = double[2]
    p = double[3]

    cnot_occ =  [[i, i+1] for i in range(s, r)]
    cnot_virt = [[i, i+1] for i in range(q, p)]
    cnot_connects = cnot_occ +[[r, q]] +cnot_virt

    # different definition from pennylane, but this is correct
    t2_circuit(-param/4.0, s, r, q, p, "XXYX", cnot_connects, circuit)  # Eq. (3)
    t2_circuit(-param/4.0, s, r, q, p, "YXYY", cnot_connects, circuit)  # Eq. (4)
    t2_circuit(-param/4.0, s, r, q, p, "XYYY", cnot_connects, circuit)  # Eq. (5)
    t2_circuit(-param/4.0, s, r, q, p, "XXXY", cnot_connects, circuit)  # Eq. (6)
    t2_circuit( param/4.0, s, r, q, p, "YXXX", cnot_connects, circuit)  # Eq. (7)
    t2_circuit( param/4.0, s, r, q, p, "XYXX", cnot_connects, circuit)  # Eq. (8)
    t2_circuit( param/4.0, s, r, q, p, "YYYX", cnot_connects, circuit)  # Eq. (9)
    t2_circuit( param/4.0, s, r, q, p, "YYXY", cnot_connects, circuit)  # Eq. (10)

    return None



def UCCSD_circuit(params, singles=None, doubles=None):
    circuit = QuantumCircuit(nqubits)
    if singles != None:
        for idx, single in enumerate(singles):
            add_fermionic_single(params[idx], single, circuit)
    if doubles != None:
        for idx, double in enumerate(doubles):
            add_fermionic_double(params[idx +len(singles)], double, circuit)  # setting t2 parameters after t1 parameters

    return circuit


D2h = [["ag",  "b1g", "b2g", "b3g", "au",  "b1u", "b2u", "b3u"],
       ["b1g", "ag",  "b3g", "b2g", "b1u", "au",  "b3u", "b2u"],
       ["b2g", "b3g", "ag",  "b1g", "b2u", "b3u", "au",  "b1u"],
       ["b3g", "b2g", "b1g", "ag",  "b3u", "b2u", "b1u", "au" ],
       ["au",  "b1u", "b2u", "b3u", "ag",  "b1g", "b2g", "b3g"],
       ["b1u", "au",  "b3u", "b2u", "b1g", "ag",  "b3g", "b2g"],
       ["b2u", "b3u", "au",  "b1u", "b2g", "b3g", "ag",  "b1g"],
       ["b3u", "b2u", "b1u", "au",  "b3g", "b2g", "b1g", "ag" ]]
D2h_table = pd.DataFrame(D2h,
                         index=["ag", "b1g", "b2g", "b3g", "au", "b1u", "b2u", "b3u"],
                         columns=["ag", "b1g", "b2g", "b3g", "au", "b1u", "b2u", "b3u"])
print("D2h symmetry:\n", D2h_table)


# equilibrium positional data from https://cccbdb.nist.gov/expgeom2x.asp
dist_list = [1.0, 1.3264, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

# {dist, [hf_ene, casci[0], vqe_ene, len(ene_hist)]}
pes_dict = {}

basis                  = "sto-3g"  #basis set
charge                 = 0         #total charge for the molecule
multiplicity           = 3         #spin multiplicity 2S +1 (Select 1 or 3)
delta_sz_from_hf_state = -1         # select from -2, -1, 0, 1, 2 at the moment

for dist in dist_list:
    print(f"\n########################################\n"
          f"### ditance is now {dist:<8f} angstrom ###\n"
          f"########################################")

    if dist <= 1.5:
        active_orb_sym = ["ag", "b1u", "b2u", "b3u", "ag", "b1u"]
    elif dist >= 2.0:
        active_orb_sym = ["ag", "b1u", "ag",  "b2u", "b3u", "b1u"]
    else:
        raise ValueError("Invalid dist detected!")

    # ------------------------------------------------------
    # xyz coordinates for atoms (angstrom)
    geometry = [("Be", (0.0, 0.0,  0.0)),
                ("H",  (0.0, 0.0,  dist)),
                ("H",  (0.0, 0.0, -dist))]
    # ------------------------------------------------------

    molecule = MolecularData(geometry, basis, multiplicity, charge)
    molecule = run_pyscf(molecule, run_scf=True, run_fci=False, run_ccsd=False, verbose=False)
    print("HF energy: {} (Hartree)".format(molecule.hf_energy))
    print("FCI energy: {} (Hartree)".format(molecule.fci_energy))

    nele = molecule.n_electrons
    norb = molecule.n_orbitals              # number of spatial orbitals
    nspinorb = nqubits = molecule.n_qubits  # number of spin orbitals
    print(f"nele = {nele}, norb = {norb}, nspinorb = nqubits = {nqubits}")


    # molecule._pyscf_data["scf"].analyze(verbose=5)
    # print(molecule._pyscf_data["scf"].__dict__)
    mol          = molecule._pyscf_data["mol"]
    mf           = molecule._pyscf_data["scf"]
    mo_occ_table = pd.DataFrame({"Energy": mf.mo_energy, "Occupancy": mf.mo_occ})
    print(mo_occ_table)

    # Setting for active electron & spinorbitals
    occ    = [0]            # spatial orbitals
    active = [1,2,3,4,5,6]  # spatial orbitals
    nele = 4
    nspinorb = nqubits = 12
    print(f"active nele = {nele}, active nspinorb = nqubits = {nqubits}")

    jw_hamiltonian     = jordan_wigner(get_fermion_operator(molecule.get_molecular_hamiltonian(occupied_indices=occ,
                                                                                               active_indices  =active)))
    qulacs_hamiltonian = create_observable_from_openfermion_text(str(jw_hamiltonian))
    jw_mat             = get_sparse_operator(jw_hamiltonian)



    print(f"MO occupancy in HF = {mf.mo_occ}")

    # Create HF state
    hf_state, hf_state_bin = create_active_hf_state(mf.mo_occ, nele, nspinorb, print_state=True)
    # Check HF energy
    hf_ene = qulacs_hamiltonian.get_expectation_value(hf_state)


    HFbra = kron_N(hf_state_bin)
    HFket = HFbra.T
    jw_matrix = get_sparse_operator(jw_hamiltonian)
    hf_ene_from_mat = np.real(HFbra.dot(jw_matrix.dot(HFket)))[0][0]
    print(f"Check HF energy:\n"
          f" pyscf           = {molecule.hf_energy:14.10f}\n"
          f" qulacs          = {hf_ene:14.10f}\n"
          f" hf_ene_from_mat = {hf_ene_from_mat:14.10f}")

    eigenenergies, eigenvecs = np.linalg.eigh(jw_matrix.toarray().real)


    decode_to_states(eigenenergies, eigenvecs, list(range(1)))

    print("\n\n=== CASCI Result ===\n")
    mycas = mcscf.CASCI(mf, len(active), (nele//2, nele//2))  # alpha=2, beta=2 in 6 spatial orbitals
    casci = mycas.kernel(verbose=False)
    print(f"casci energy = casci[0]")


    #hf_state_bin = "0b000000010111"
    singles, doubles = create_excitations(hf_state_bin, delta_sz_from_hf_state)
    print(f"singles = {singles}\n\ndoubles = {doubles}")

    initial_sym = judge_initial_sym(active_orb_sym, hf_state_bin, D2h_table)
    print(f"\ninitial state symmetry = {initial_sym}")
    singles, doubles = judge_D2h_excitation(singles, doubles, active_orb_sym, initial_sym, D2h_table)
    print(f"Symmetry reduced singles = {singles}\n\nSymmetry reduced doubles = {doubles}")



    def compute_energy(params):  # singles, doubles, state are working as global variables
        #hf_state = create_active_hf_state(np.array([2., 2., 1., 1., 0., 0., 0.]), nele, nqubits)
        hf_state = create_active_hf_state(np.array([2., 2., 2., 0., 0., 0., 0.]), nele, nqubits)
        #hf_state = create_active_hf_state(mf.mo_occ, nele, nqubits)        # HF state initialization each time
        circuit = UCCSD_circuit(params, singles=singles, doubles=doubles)  # create quantum cirquit
        circuit.update_quantum_state(hf_state)                             # apply quantum cirquit to state

        return qulacs_hamiltonian.get_expectation_value(hf_state)




    ene_hist = []
    if multiplicity == 1:
        params = np.zeros(len(singles +doubles))
    elif multiplicity == 3:
        params = np.full(len(singles +doubles), 0.6)
        #params = np.random.normal(0, np.pi, len(singles +doubles))
    ene_hist.append(compute_energy(params))

    method = "BFGS"
    options = {"disp": True, "maxiter": 1000, "gtol": 1e-10}

    #method = "COBYLA"
    #options = None

    opt = minimize(compute_energy,
                   params,
                   method=method,
                   callback=lambda params: ene_hist.append(compute_energy(params)),
                   options=options,
                   jac="2-point",
                   tol=1e-10)

    vqe_ene = compute_energy(opt.x)
    print(f"energy history [Eh]: {ene_hist}")
    print(f"VQE energy [Eh]: {vqe_ene}")
    print(f"VQE parameter: {opt.x}")
    pes_dict[dist] = [hf_ene, casci[0], vqe_ene, len(ene_hist)]


with open(f"result_multi{multiplicity}.dat", mode="w") as f:
    f.write("distance [ang]   HF energy [Eh]   CASCI energy [Eh]   VQE energy [Eh]   VQE total iteration\n")
    for dist in dist_list:
        f.write(f"{dist:<8f}         {pes_dict[dist][0]:<14.10f}   {pes_dict[dist][1]:<14.10f}      {pes_dict[dist][2]:<14.10f}    {pes_dict[dist][3]:<5d}\n")

print("\nNormal termination of the VQE calculation.\n")
