# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:34:08 2023

@author: Css
"""

from ideal_gas import ideal_gas
import numpy as np
from scipy.integrate import quad

class component:
    name : str  = 'NoName'
    Mw : float
    Hf : float
    T0 : float = 298.15
    Cp_coeff = [0.0, 0.0, 0.0, 0.0, 0.0]
    mu : float
    
    def __init__(self, name, Mw, Cp, Hf, mu):
        self.name = name 
        self.Mw = Mw                                                        # kg/kmol
        self.Cp_coeff = Cp                                                  # j/kmol.K
        self.Hf = Hf                                                        # j/kmol
        self.EOS = ideal_gas()                                              # kmol/m^3
        self.mu = mu                                                        # pa.sec or kg/m.sec
     
    def mass_density(self, P, T):
        return self.EOS.molar_density(P, T) * self.Mw                       # kg/m^3
    
    def molar_density(self, P, T):
        return self.EOS.molar_density(P, T)                                 # kmol/m^3
    
    def viscosity(self, T):
        return self.mu                                                      # pa.sec
    
    def heat_capacity(self,T):
        return np.sum(self.Cp_coeff * np.array([T**i for i in range(5)]))   # j/kmol.K
    
    def Hs(self, T):
        return (quad(self.heat_capacity, self.T0, T)[0])                    # j/kmol
    
    def Ha(self, T):
        return self.Hf + self.Hs(T)                                         # j/kmol