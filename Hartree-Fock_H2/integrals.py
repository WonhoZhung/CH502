import numpy as np
from math import pi, exp

def S(g1, g2):
    """
    Calculate (g1|g2)
    """
    a, R_a = g1.a, g1.R
    b, R_b = g2.a, g2.R
    norm = sum((R_a - R_b)**2)
    retval = (pi/(a+b))**1.5
    retval *= exp(-norm*a*b/(a+b))
    return retval

def T(g1, g2):
    """
    Calculate (g1|-0.5*\nabla^2|g2)
    """
    a, R_a = g1.a, g1.R
    b, R_b = g2.a, g2.R
    norm = sum((R_a - R_b)**2)
    retval = (pi/(a+b))**1.5
    retval *= a*b/(a+b)
    retval *= 3-2*norm*a*b/(a+b)
    retval *= exp(-norm*a*b/(a+b))
    return retval

def V(g1, g2, Z_c, R_c):
    """
    Calculate (g1|-Z_c/R_c|g2)
    """
    a, R_a = g1.a, g1.R
    b, R_b = g2.a, g2.R
    p, R_p = g1.prod(g2)
    norm = sum((R_a - R_b)**2)

    return

def multi(g1, g2, g3, g4):
    """
    caculate (g1,g2|g3,g4)
    """
    return
