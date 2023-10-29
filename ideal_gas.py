R = 8314                                           #j/kmol-k or pa-m^3/kmol-k

class ideal_gas:
    
    def molar_density(self, P, T):
        return P/(R*T)                              #kmol/m^3
