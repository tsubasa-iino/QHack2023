# spin-flip [QHack2023]
A project by team Qtest



To compute the ground state energy of BeH2 molecule, we have implemented sym-UCCSD<sup>[1]</sup> quantum cirquit (D2h molecular symmetry) using pyscf, openfermion, qulacs library.
[1] https://arxiv.org/abs/2109.02110


The use of symmetry reduces the number of excitation operators in UCCSD, resulting in reducing the depth of quantum cirquit greatly.


In singlet case:
the number of one-electron excitations is 16 -> 4
the number of two-electron excitations is 76 -> 14
See multi1.log for the details.


In triplet case (spin-flip):
the number of one-electron excitations is 15 -> 3
the number of two-electron excitations is 75 -> 13
See multi3.log for the details.



In addition, we implemented a new idea VQE algorithm: spin-flip<sup>[2]</sup>.
Spin-flip tequnique has been developed mainly by Anna I. Krylov group in quantum chemical calculation community.
[2] https://www.sciencedirect.com/science/article/abs/pii/S0009261401002871


The main idea is based on the fact that triplet state (Ms=1) shows essentially single reference character.
Thus, that state could be described by Hartree-Fock method with high-accuracy.
For computing ground state energy, the excitation changing up spin to down spin is applied to the triplet HF reference state: spin-flip.
We used spin-flip technique with sym-UCCSD targeting Ag ground state.


To verify above ideas (sym-UCCSD and spin-flip), we compute Ag ground state energy of BeH2 molecule with some Be-H bond distances: [1.0, 1.3264, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0] (angstrom) using sto-3g basis and active space (4e, 6o) -> 12 qubits


The main code of this project: spin-flip.py
You can use this program by: \$ python3 spin-flip.py > logfile
The logfile for singlet sym-UCCSD: multi1.log
The logfile for triplet spin-flip sym-UCCSD: multi1.log
The obtained energies are summarized in: result_multi1.dat and result_multi3.dat
PES.jpg is obtained by: \$ python3 plot.py result_multi1.dat result_multi3.dat


To make it easy to read our program code, we also made jupyter-notebook code: spin-flip.ipynb
In this notebook, we also checked MO symmetry and attributed each MO to D2h symmetry (Maybe this check can be automated by pyscf).


The obtained potential energy curve is shown in: PES.jpg
Unfortunately, the spin-flip technique converged some incomprehensible high-energy state or HF state.
Therefore, further studies are needed to investigate spin-flip properties in VQE.
However, in singlet case, sym-UCCSD gave almost the same energy with CASCI energy in all bond distances, implying the use of symmetry is very promising technique for further research in VQE.
