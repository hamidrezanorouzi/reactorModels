# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:37:30 2023

@author: Css
"""

import numpy as np
from component import component

class component_set:
    
    def __init__(self, components):
        self.components = list(components)
        
    def append_component(self, comp : component):
        self.components.append(comp)
        
    def num_components(self):
        return len(self.components)
    
    def mixture_Mw(self, y):
        return np.sum([y[i] * self.components[i].Mw for i in range(self.num_components())])                 # kg/kmol
        
    def mixture_Cp(self, y, T):
        return np.sum([y[i] * self.components[i].heat_capacity(T) for i in range(self.num_components())])   # j/kmol.K
    
    def mixture_molar_conc(self, y, P, T):
        return [y[i] * self.components[i].molar_density(P, T) for i in range(self.num_components())]        # kmol/m^3
    
    def mixture_viscosity(self, y, T):
        return np.sum([self.mass_fraction(y)[i] * self.components[i].viscosity(T) for i in range(self.num_components())])          # kg/m.sec or pa.sec
    
    def mixture_mass_density(self, y, P, T):
        return np.sum ([self.mass_fraction(y)[i] * component.mass_density(P, T) for i, component in enumerate(self.components)])   # kg/m^3
    
    def mass_fraction(self, y):
        return [y[i] * self.components[i].Mw / self.mixture_Mw(y) for i in range(self.num_components())]