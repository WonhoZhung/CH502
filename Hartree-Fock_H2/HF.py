import numpy as np
from integrals import S, T, V, multi_electron


def nuc_repulsion(molecule):
    retval = 0.
    for i in range(molecule.N):
        for j in range(i+1, molecule.N):
            if i == j: continue
            atom_i, atom_j = molecule.atoms[i], molecule.atoms[j]
            retval += \
                atom_i.Z*atom_j.Z/np.linalg.norm(atom_i.coord-atom_j.coord)
    return retval

def calc_matrices(basis_set, molecule):
    """
    params
    basis_set:: list of basis object
    molecule:: molecule object
    """
    # Number of basis function
    N = len(basis_set)

    # Initialize matrices
    kinetic = np.zeros((N, N))
    potential = np.zeros((N, N))
    overlap = np.zeros((N, N))
    multi = np.zeros((N, N, N, N))
    for i, bi in enumerate(basis_set):
        for gi, di in bi.gs:
            for j, bj in enumerate(basis_set):
                for gj, dj in bj.gs:
                    kinetic[i,j] += di*dj*T(gi, gj)
                    overlap[i,j] += di*dj*S(gi, gj)
                    for z, coord in zip(molecule.Zs, molecule.coords):
                        potential[i,j] += di*dj*V(gi, gj, z, coord)
                    for k, bk in enumerate(basis_set):
                        for gk, dk in bk.gs:
                            for l, bl in enumerate(basis_set):
                                for gl, dl in bl.gs:
                                    multi[i,j,k,l] += \
                                            multi_electron(gi, gj, gk, gl)
                                    multi[i,j,k,l] *= \
                                            di*dj*dk*dl
    return kinetic, potential, overlap, multi

def calc_error(p_old, p):
    """
    Calculate RMSE between p_old and p
    """
    return np.linalg.norm(p_old - p)

def orthogonalize(m):
    """
    Symmetric orthogonalization of matrix m
    """
    _, U = np.linalg.eig(m)
    m_diag = np.dot(U.T, np.dot(m, U))
    m_diag = np.diag(np.diagonal(m_diag)**-0.5)
    X = np.dot(U, np.dot(m_diag, U.T))
    return X

def run(basis_set, molecule, thr=1e-9, max_iter=1000):
    """
    Run the SCF Iteration
    """

    # Print Input Info
    print(f"### Molecule Info\n")
    print(molecule)

    # Number of basis function
    N = len(basis_set)

    # Initialize P
    P = np.zeros((N, N))

    # Initialize error
    error = 1e9

    # Initialize counter
    cnt = 0

    # Start iteration
    print(f"### Start SCF")
    while error > thr:
        # Back-up P
        P_old = P

        # Calculate matrices
        kinetic, potential, overlap, multi = \
                calc_matrices(basis_set, molecule)
        H_core = kinetic + potential
        G_ee = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        G_ee[i,j] += P[k,l]*(multi[i,j,k,l]-0.5*multi[i,l,k,j])
        F = H_core + G_ee
        
        # Solving the Roothan equation
        # F' = X^TFX
        X = orthogonalize(overlap) 
        F_prime = np.dot(X.T, np.dot(F, X))

        # F'C' = C'e
        # Eigenvectors:: C'
        # Eigenvalues:: e
        eigenvalues, eigenvectors = np.linalg.eig(F_prime)
        e = eigenvalues[eigenvalues.argsort()]
        C_prime = eigenvectors[eigenvalues.argsort()]

        # C = XC'
        C = np.dot(X, C_prime)

        # Calculate new P
        P = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                P[i,j] += 2*C[i,0]*C[j,0]

        error = calc_error(P_old, P)
        
        # Print log
        print(f"### {cnt}th iteration:: error - {error:.6f}")
        print(f"### Converged:: {error < thr}")

        cnt += 1
        
        assert cnt < max_iter, "### Iteration exceeds MAX counts"

    print(f"### SCF Converged!!!")
    return P, C, e


if __name__ == "__main__":

    from molecule import Atom, Molecule
    from basis import STO_3G
    atoms = [
           Atom(1, [0., 0., 0.]),
           Atom(1, [0., 0., 1.])
    ]
    mol = Molecule(atoms)
    
    a_list = [3.42525 , 0.623914, 0.168855]
    d_list = [0.154329, 0.535328, 0.444635]
    R_list = mol.coords
    basis_set = [STO_3G(a_list, d_list, R) for R in R_list]

    t, v, s, ee = calc_matrices(basis_set, mol)
    print(t)
    print(v)
    print(s)
    print(ee)

