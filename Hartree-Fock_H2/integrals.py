import numpy as np
from math import pi, exp
from scipy.special import erf


"""
Reference:: Modern Quantum Chemistry by Szabo & Ostlund
"""

def _EF(t):
    if t == 0: return 1.
    return (0.5*(pi/t)**0.5)*erf(t**0.5)

def S(g1, g2):
    """
    Calculate (g1|g2)
    """
    a, R_a = g1.a, g1.R
    b, R_b = g2.a, g2.R
    norm = np.linalg.norm(R_a - R_b)**2
    retval = (pi/(a+b))**1.5
    retval *= (4*a*b/pi**2)**0.75
    retval *= exp(-norm*a*b/(a+b))
    return retval

def T(g1, g2):
    """
    Calculate (g1|-0.5*\nabla^2|g2)
    """
    a, R_a = g1.a, g1.R
    b, R_b = g2.a, g2.R
    norm = np.linalg.norm(R_a - R_b)**2
    retval = (pi/(a+b))**1.5
    retval *= a*b/(a+b)
    retval *= 3-2*norm*a*b/(a+b)
    retval *= (4*a*b/pi**2)**0.75
    retval *= exp(-norm*a*b/(a+b))
    return retval

def V(g1, g2, Z_c, R_c):
    """
    Calculate (g1|-Z_c/R_c|g2)
    """
    a, R_a = g1.a, g1.R
    b, R_b = g2.a, g2.R
    _, p, R_p = g1.prod(g2)
    norm = np.linalg.norm(R_a - R_b)**2
    norm2 = np.linalg.norm(R_p - R_c)**2
    retval = -2*pi*Z_c/(a+b)
    retval *= (4*a*b/pi**2)**0.75
    retval *= exp(-norm*a*b/(a+b))
    retval *= _EF(norm2*(a+b))
    return retval

def two_electron(g1, g2, g3, g4):
    """
    Calculate (g1,g2|g3,g4)
    """
    a, R_a = g1.a, g1.R
    b, R_b = g2.a, g2.R
    c, R_c = g3.a, g3.R
    d, R_d = g4.a, g4.R
    _, p, R_p = g1.prod(g2)
    _, q, R_q = g3.prod(g4)
    norm = np.linalg.norm(R_a - R_b)**2
    norm2 = np.linalg.norm(R_c - R_d)**2
    norm3 = np.linalg.norm(R_p - R_q)**2
    retval = 2*pi**2.5/((a+b)*(c+d)*(a+b+c+d)**0.5)
    retval *= (4*a*b/pi**2)**0.75
    retval *= (4*c*d/pi**2)**0.75
    retval *= exp(-norm*a*b/(a+b))
    retval *= exp(-norm2*c*d/(c+d))
    retval *= _EF(norm3*(a+b)*(c+d)/(a+b+c+d))
    return retval

