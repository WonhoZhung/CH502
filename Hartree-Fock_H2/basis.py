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
        
    def prod(self, g2):
        """
        inputs
        self, g2:: Gaussian class
        
        output
        g:: Gaussian class, product of self & g2
        """
        p = self.a + g2.a
        R_p = (self.a*self.R+g2.a*g2.R)/p
        g = Gaussian(p, R_p)
        return g
      
        
class STO_3G():
  
    def __init__(self, a_list, R_list, d_list):
        """
        params
        a_list:: list of a
        R_list:: list of R
        d_list:: list of d, contraction coefficients
        """
        assert len(a_list) == len(R_list) == len(d_list)
        self.a_list = a_list
        self.R_list = R_list
        self.d_list = d_list
        
        self.gs = [(Gaussian(a, R), d) for a, R, d in \
                zip(a_list, R_list, d_list)]
        self.N = len(self.gs)
      
