
#تغییرات دادم دیتاهای ورودی رو

import numpy as np
from component import component
from component_set import component_set
from reaction import reaction
from reaction_set import reaction_set
from reactor import reactor

#define constants and parameters

F0 = [0.02087, 0, 0, 0, 0, 0, 0, 0, 0.0139335091]       # initial flow rate (kmol/sec)
# massF0 = 0.8783982648                              # kg/sec
Ti = 953                                            # initial tempreture(K)
P0 =303000                                          # initial pressure(pa)                   
                     
#list of component      
C2H6 = component ('C2H6', 30.06904  , [33643.99276, 127.6521301, -0.045480849, 7.29825E-06, -4.34961E-10], -83.852e+06, 2.60811056228900e-005)
C2H4 = component ('C2H4', 28.05316 , [33188.05475, 87.15891277, -0.030904916, 4.94374E-06, -2.94008E-10], 52.50e+06, 2.73856837201462e-005)
H2   = component ('H2'  , 2.01588  , [24383.84618, 6.872418638, -0.001217189, 1.28119E-07, -5.72672E-12], 0, 2.43657377632359e-005)
C3H8 = component ('C3H8', 44.09562  , [55436.86764, 171.3683459, -0.061237019, 9.84716E-06, -5.87761E-10], -104.68e+06, 2.48855589655182e-005)
CH4  = component ('CH4' , 16.04246 , [13597.76674, 83.84069976, -0.028011217, 4.44765E-06, -2.62322E-10], -74.60e+06, 2.72690307579320e-005)
C3H6 = component ('C3H6', 42.07974 , [50205.79329, 135.4881823, -0.04839834, 7.78138E-06, -4.64422E-10], 20.00e+06, 2.59244520821420e-005)
C2H2 = component ('C2H2', 26.03728  , [38733.13882, 40.60528892, -0.013371304, 2.05334E-06, -1.15237E-10], 228.2e+06, 2.89902738599355e-005)
C4H6 = component ('C4H6', 54.09044  , [133032.4296, 32.57634006, 0.009507536, -1.72869E-06, 6.29963E-11], 110.00e+06, 2.50134586143583e-005)
steam = component ('steam', 18.0153 , [30.092e+03, 6.832514, 6.793435e-03, 2.53448e-6, 0], -241.818e+06,3.57342710996558e-005)
allComponents = [C2H6, C2H4, H2, C3H8, CH4, C3H6, C2H2, C4H6, steam]

compSet = component_set(allComponents)
     
#list of reaction
reaction1 = reaction(4.652e13, 272.83864e+06, allComponents, {'C2H6':[-1,1], 'C2H4':[1,0], 'H2':[1,0]})
reaction2 = reaction(9.902744695e8, 137.9447384e+06, allComponents, {'C2H4':[-1,1], 'H2':[-1,1], 'C2H6':[1,0]})
reaction3 = reaction(3.85e11, 273.006e+06, allComponents, {'C2H6':[-2,1], 'C3H8':[1,0], 'CH4':[1,0]})
reaction4 = reaction(9.814e08 , 154.47328e+06, allComponents, {'C3H6':[-1,1], 'C2H2':[1,0], 'CH4':[1,0]})
reaction5 = reaction(6.03404188e4, 29.70901722e+06, allComponents, {'C2H2':[-1,1], 'CH4':[-1,1], 'C3H6':[1,0]})
reaction6 = reaction(1.026e12, 172.63184e+06, allComponents, {'C2H2':[-1,1], 'C2H4':[-1,1], 'C4H6':[1,0]})
reaction7 = reaction(7.083e13, 252.83912e+06, allComponents, {'C2H4':[-1,1], 'C2H6':[-1,1], 'C3H6':[1,0], 'CH4':[1,0]})
allReactions = [reaction1,reaction2, reaction3, reaction4, reaction5, reaction6, reaction7]

reactSet = reaction_set(allReactions)

def flux(z):
    return (96e+03 - (85.91e+03 * z/95) + (42.955e+03 * ((z/95)**2)))        # w/m^2
        
plug = reactor(95, 0.108, reactSet, compSet, 0.02, 1, 0.07)

solve = plug.solve(F0, Ti, P0, flux)

outlets = plug.outlet()
print(outlets)

# plots = plug.plot()

# profiles = plug.profile()

# validation with article
EC2H6 = (0.00723 -outlets[0])/0.00723*100
EC2H4 =(0.01114 - outlets[1])/0.01114*100
ECH4 = (0.00186 - outlets[4])/0.00186 *100
EC3H6 = (0.00007 - outlets[5])/0.00007 *100
ETEMP = (1077.086 - outlets[18])/1077.086*100
print(EC2H6 , EC2H4, ECH4, EC3H6, ETEMP)



