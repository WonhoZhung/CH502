import numpy as np

class Gaussian():
  
    def __init__(self, a, R):
        """
        params
        a:: Gaussian Exponents
        d:: Gaussian Centers
        """
        self.a = a
        self.R = R
        
    def prod(self, g1, g2):
        """
        inputs
        g1, g2:: Gaussian class
        
        output
        g:: Gaussian class, product of g1 & g2
        """
        p = g1.a + g2.a
        R_p = (g1.a*g1.R+g2.a*g2.R)/p
        g = Gaussian(p, R_p)
        return g
      
        
class STO_3G():
  
    def __init__(self, as, Rs, ds):
        """
        params
        as:: list of a
        Rs:: list of R
        ds:: list of d, contraction coefficients
        """
        assert len(as) == len(Rs) == len(ds)
        self.as = as
        self.Rs = Rs
        self.ds = ds
        
        self.gs = [(Gaussian(a, R), d) for a, R, d in zip(as, Rs, ds)]
        self.N = len(self.gs)
      
