# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:48:43 2023

@author: Css
"""


import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt



class reactor:
    reactorLength : float
    tubeDiameter : float
    particleDiameter : float
    prosity : float    
    roughness : float
    componentsFlowRate : list = []
    Tempreture : list = []
    length : list = []
    pressure : list = []
    componentsMassFlowRate : list = []
   
    
    def __init__(self, reactorLength, tubeDiameter, reactions, components, particleDiameter, prosity, roughness):
        self.reactorLength = reactorLength                          # m
        self.tubeDiameter = tubeDiameter                            # m
        self.tubeArea = np.pi*tubeDiameter**2/4                     # Reactor cross-sectional area (m^2)
        self.reactions = reactions
        self.components = components
        self.particleDiameter = particleDiameter                    # m
        self.prosity = prosity
        self.roughness = roughness
        self.rQ = self.__Q

    def volume(self, z):
        return self.At*z                                            # Reactor volume (m^3)
    
    def __Q(self, z):                                                # w/m^2
        return 0.0;
    
    def __mass_balance(self, C, T):                                                                             # kmol/m
        dFdz = np.zeros_like(C)
        dFdz = self.reactions.all_net_productions(C, T) * self.tubeArea
        return dFdz
    
    def __energy_balance(self, F, y, C, T, z):                                                                   # K/m
        dTdz = (((self.rQ(z)*np.pi*self.tubeDiameter) - ((np.sum(self.reactions.Hrxn_dot_rate(T, C)) *self.tubeArea))) / (np.sum(self.components.mixture_Cp(y, T)* np.sum(F))))
        return dTdz
    
    def __pressure_drop(self, y, P, T, z):                                                                       # pa/m
        v = self.M0 / (self.tubeArea * self.components.mixture_mass_density(y, P, T))                             # m/s
        Re = (self.components.mixture_mass_density(y, P, T) * v * self.tubeDiameter) / (self.components.mixture_viscosity(y, T))
        f = (-1.8 * math.log((6.9/Re) + ((self.roughness/3.7)**1.11)))**(-2)
        term1 = (150*self.components.mixture_viscosity(y, T) *(1-self.prosity)**2*v)/(self.particleDiameter**2*self.prosity**3)
        term2 = (1.75*self.components.mixture_mass_density(y, P, T) *(1-self.prosity)*v**2)/(self.particleDiameter * self.prosity)
        term3 = (f * self.components.mixture_mass_density(y, P, T) * (v)**2) / (2 * self.tubeDiameter)
        dpdz = -(term1 + term2 + term3)
        return dpdz
    
    def __ode_system(self, z, x):
        F = x[:self.components.num_components()]
        T = x[self.components.num_components()]
        P = x[self.components.num_components()+1]
        y = F / np.sum(F)
        C = self.components.mixture_molar_conc(y, P, T)
        dFdz = self.__mass_balance(C, T)
        dTdz = self.__energy_balance(F, y, C, T, z)
        dpdz = self.__pressure_drop(y, P, T, z)
        return np.concatenate((dFdz, [dTdz, dpdz]))
    
    def solve(self, F0, Ti, P0, reactQ = None):
        if reactQ is None:
            self.rQ = self.__Q
        else:
            self.rQ = reactQ
            
        massF0 =[]
        for i in range(self.components.num_components()):
            massF0.append(F0[i]*self.components.components[i].Mw)
        self.M0 = np.sum(massF0)
        
        x0 = np.concatenate((F0, np.array([Ti, P0])))
        z = np.linspace(0, self.reactorLength)
        sol = solve_ivp(self.__ode_system, (z[0], z[-1]), x0, method='LSODA', rtol=1e-3)

        F = sol.y[:self.components.num_components()]  
        T = sol.y[self.components.num_components()]  
        P = sol.y[self.components.num_components()+1]
        z = sol.t
        M = []
        for i in range(self.components.num_components()):
            M.append(F[i]*self.components.components[i].Mw)
        self.componentsFlowRate = F
        self.Tempreture = T
        self.length = z
        self.pressure = P
        self.componentsMassFlowRate = M
        return np.concatenate((F, M, [T, P]))
    
    def plot(self):
        for i in range(self.components.num_components()):
            plt.plot(self.length, self.componentsFlowRate[i], 'm')
            plt.xlabel('Reactor length (m)')
            plt.ylabel('{} (kmol/sec)'.format(self.components.components[i].name))        
            plt.show()  
        plt.plot(self.length, self.Tempreture, 'm')
        plt.ylabel('Temperature (K)')
        plt.xlabel('Reactor length (m)')
        plt.show()
        plt.plot(self.length, self.pressure, 'm')
        plt.ylabel('Pressure (pa)')
        plt.xlabel('Reactor length (m)')
        plt.show()
        
    def profile(self):
        for j in range(self.components.num_components()):
            print("\t\t\t{} (kmol/sec)".format(self.components.components[j].name))
            for i in range(len(self.length)):
                print("{:.3f}\t\t{:.9f}".format(self.length[i], self.componentsFlowRate[j][i]))
            print("\n")
        print("\t\t\tTempreture (K)")
        for i in range(len(self.length)):
            print("{:.3f}\t\t{:.9f}".format(self.length[i], self.Tempreture[i]))
        print("\n")
        print("\t\t\tPressure (pa)")
        for i in range(len(self.length)):
            print("{:.3f}\t\t{:.9f}".format(self.length[i], self.pressure[i]))
        print("\n") 
        
    def outlet(self):
        out = (len(self.length)-1)
        F =[]
        for j in range(self.components.num_components()):
            F.append(self.componentsFlowRate[j][out])
        M =[]
        for j in range(self.components.num_components()):
            F.append(self.componentsMassFlowRate[j][out])
        T = self.Tempreture[out]
        P = self.pressure[out]
        return np.concatenate((F, M, [T, P]))
    
