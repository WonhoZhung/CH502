import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from molecule import Atom, Molecule
from basis import STO_3G
from HF import run, nuc_repulsion
from utils import A_LIST, D_LIST
from tqdm import tqdm

x_list, y_list = [], []
for d in tqdm([0.01*x for x in range(50,450)]):

    atoms = [
           Atom(1, [0., 0., 0.        ]),
           Atom(1, [0., 0., d/0.529177])
    ]
    mol = Molecule(atoms)
    
    R_list = mol.coords
    basis_set = [STO_3G(A_LIST, D_LIST, R) for R in R_list]

    _, _, e, E = run(basis_set, mol)
    nuc_rep = nuc_repulsion(mol)

    x_list.append(d)
    y_list.append(E+nuc_rep)

plt.plot(x_list, y_list)
plt.xlabel("Distance (Å)")
plt.ylabel("Dissociation Energy (Hartree)")
plt.title("Dissociation_Curve_H2")
plt.savefig("Dissociation_Curve_H2.png")
