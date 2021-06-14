from basis import STO_3G
from HF import run, nuc_repulsion
from molecule import xyz_to_mol
from utils import A_LIST, D_LIST

import sys


# Read the xyz file to get molecule
H2 = xyz_to_mol(sys.argv[1])

# Define the basis for H2

R_list = H2.coords
basis_set = [STO_3G(A_LIST, D_LIST, R) for R in R_list]

# SCF
try: P, C, e, E = run(basis_set, H2)
except Exception as e: print(e); exit()

# Nuclear Repulsion
nuc_rep = nuc_repulsion(H2)

print()
print(f"### Total Energy:: {E+nuc_rep:.5f}")
print(f"### Electronic Energy:: {E:.5f}")
print(f"### Nuclear Repulsion:: {nuc_rep:.5f}")
print(f"### Orbital Energy::\t{e[0]:.5f}\t{e[1]:.5f}")
print(f"### Density Matrix::\n{P}")

