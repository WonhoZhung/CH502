import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from molecule import Atom, Molecule
from basis import STO_3G
from HF import run, nuc_repulsion

a_list = [3.42525 , 0.623914, 0.168855]
d_list = [0.154329, 0.535328, 0.444635]

x_list, y_list = [], []
for d in [0.1*x for x in range(7, 22)]:

    atoms = [
           Atom(1, [0., 0., 0.]),
           Atom(1, [0., 0., d ])
    ]
    mol = Molecule(atoms)
    
    R_list = mol.coords
    basis_set = [STO_3G(a_list, d_list, R) for R in R_list]

    _, _, e = run(basis_set, mol)
    nuc_rep = nuc_repulsion(mol)

    x_list.append(d)
    y_list.append(e[0]+nuc_rep)

plt.plot(x_list, y_list)
plt.savefig("Dissociation_Curve_H2.png")
