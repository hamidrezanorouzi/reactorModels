# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:35:24 2023

@author: Css
"""

R = 8314  

class ideal_gas:
    
    def molar_density(self, P, T):
        return P/(R*T)                              #kmol/m^3