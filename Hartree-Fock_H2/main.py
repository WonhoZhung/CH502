from basis import STO_3G
from HF import run, nuc_repulsion
from molecule import xyz_to_mol

import sys


# Read the xyz file to get molecule
H2 = xyz_to_mol(sys.argv[1])

# Define the basis for H2

#----------------------------------------------------------------------
# Basis Set Exchange
# Version v0.8.13
# https://www.basissetexchange.org
#----------------------------------------------------------------------
#   Basis set: STO-3G
# Description: STO-3G Minimal Basis (3 functions/AO)
#        Role: orbital
#     Version: 1  (Data from Gaussian09)
#----------------------------------------------------------------------
# BASIS "ao basis" PRINT
# #BASIS SET: (3s) -> [1s]
# H    S
#     0.3425250914E+01       0.1543289673E+00
#     0.6239137298E+00       0.5353281423E+00
#     0.1688554040E+00       0.4446345422E+00
# END

a_list = [3.42525 , 0.623914, 0.168855]
d_list = [0.154329, 0.535328, 0.444635]
R_list = H2.coords
basis_set = [STO_3G(a_list, d_list, R) for R in R_list]

# SCF
try: P, C, e = run(basis_set, H2)
except Exception as e: print(e); exit()

# Nuclear Repulsion
nuc_rep = nuc_repulsion(H2)

print()
print(f"### Nuclear Repulsion:: {nuc_rep:.5f}")
print(f"### Orbital Energy::\t{e[0]:.5f}\t{e[1]:.5f}")
print(f"### Density Matrix::\n{P}")
#print(C)

