# spin-flip [QHack2023]
A project by team Qtest<br>
<br>
<br>
To compute the ground state energy of BeH2 molecule, we have implemented sym-UCCSD<sup>[1]</sup> quantum circuit (D2h molecular symmetry) using pyscf, openfermion, qulacs library.<br>
[1] https://arxiv.org/abs/2109.02110<br>
<br>
<br>
The use of symmetry reduces the number of excitation operators in UCCSD, resulting in reducing the depth of quantum circuit greatly.<br>
<br>
<br>
In singlet case:<br>
the number of one-electron excitations is 16 -> 4<br>
the number of two-electron excitations is 76 -> 14<br>
See multi1.log for the details.<br>
<br>
<br>
In triplet case (spin-flip):<br>
the number of one-electron excitations is 15 -> 3<br>
the number of two-electron excitations is 75 -> 13<br>
See multi3.log for the details.<br>
<br>
<br>
In addition, we implemented a new idea VQE algorithm: spin-flip<sup>[2]</sup>.<br>
Spin-flip tequnique has been developed mainly by Anna I. Krylov group in quantum chemical calculation community.<br>
[2] https://www.sciencedirect.com/science/article/abs/pii/S0009261401002871<br>
<br>
<br>
The main idea is based on the fact that triplet state (Ms=1) shows essentially single reference character.<br>
Thus, that state could be described by Hartree-Fock method with high-accuracy.<br>
For computing ground state energy, the excitation changing up-spin to down-spin is applied to the triplet HF reference state: spin-flip.<br>
We used spin-flip technique with sym-UCCSD targeting Ag ground state.<br>
<br>
<br>
To verify above ideas (sym-UCCSD and spin-flip), we compute Ag ground state energy of BeH2 molecule with some Be-H bond distances: [1.0, 1.3264, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0] (angstrom) using sto-3g basis and active space (4e, 6o) -> 12 qubits<br>
<br>
<br>
The main code of this project: spin-flip.py<br>
You can use this program by: \$ python3 spin-flip.py > logfile<br>
The logfile for singlet sym-UCCSD: multi1.log<br>
The logfile for triplet spin-flip sym-UCCSD: multi3.log<br>
The obtained energies are summarized in: result_multi1.dat and result_multi3.dat<br>
PEC.jpg is obtained by: \$ python3 plot.py result_multi1.dat result_multi3.dat<br>
<br>
<br>
To make it easy to read our program code, we also made jupyter-notebook code: spin-flip.ipynb (https://nbviewer.org/github/tsubasa-iino/spin-flip/blob/main/spin-flip.ipynb)<br>
In this notebook, we also checked each MO symmetry and attributed it to D2h symmetry (Maybe this check can be automated by pyscf). The generated cube files are stored in ./cmo/ directory.<br>
<br>
<br>
The obtained potential energy curve is shown in: PEC.jpg<br>
Unfortunately, the spin-flip technique converged some incomprehensible high-energy state or HF state.<br>
Therefore, further studies are needed to investigate spin-flip properties in VQE.<br>
However, in singlet case, sym-UCCSD gave energy values close to CASCI energies in most of bond distances, implying the use of symmetry is very promising technique for further research in VQE.<br>
