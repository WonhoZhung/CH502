import numpy as np
import sys

class Atom():
    
    def __init__(self, Z, coord):
        """
        params
        Z:: atom charge
        coord:: atom coordinate
        """
        self.Z = Z
        self.coord = coord

        
class Molecule():
  
    def __init__(self, atoms):
        """
        params
        atoms:: atoms consisting molecule
        """
        self.atoms = atoms
        
        self.N = len(atoms)
        self.Zs = [a.Z for a in atoms]
        self.coords = [a.coord for a in atoms]


def xyz_to_mol(xyz_fn):
    with open(xyz_fn, 'r') as f: lines = [l.strip() for l in f.readlines()]
    charge_dict = {'H': 1}
    atom_list = []
    for l in lines[2:]:
        split = l.split()
        a = Atom(split[0], [float(s) for s in split[1:]])
    return Molecule(atom_list)
        
# Define H2
# H2 = xyz_to_mol(sys.argv[1])
