import numpy as np
import math
from ideal_gas import R

class reaction:
    A: float
    E: float
    
    def __init__(self, A, E, components, reactionDict):
        self.A = A                                           # m^3/kmol.sec or 1/sec
        self.E = E                                           # j/kmol
        self.components = components
        self.nu = np.zeros(len(components))
        self.coeff = np.zeros(len(components))
        for i in range(len(components)):
            nuCoeff = reactionDict.get(components[i].name, [0,0])
            self.nu[i]= nuCoeff[0]
            self.coeff[i]=nuCoeff[1]
            
    def k(self, T):                                           # m^3/kmol.sec or 1/sec
        return (self.A * math.exp(-self.E / (R * T)))         # Arrhenius equation
            
    def rate_of_reaction(self, C, T):                         # kmol/m^3.sec
        return (self.k(T) * (np.prod( C ** self.coeff)))

    def net_production(self, C, T):
        return self.rate_of_reaction(C, T) * self.nu           # kmol/m^3.sec
    
    def Hrxn(self, T):
        hrxn = np.array([component.Ha(T) for component in self.components])
        return np.sum(hrxn * self.nu)                          # j/kmol
