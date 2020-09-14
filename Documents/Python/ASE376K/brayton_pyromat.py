import pyromat as pyro
import numpy as np
import matplotlib.pyplot as plt

pyro.config['unit_pressure'] = 'kPa'
pyro.config['unit_temperature'] = 'K'
pyro.config['unit_energy'] = 'kJ'

air = pyro.get('ig.air')

T1 = 300    # kelvin
T3 = 1600

p1 = 50     # kPa
p2 = p1*50
p3 = p2
p4 = p1
R = 287; g = 1.4
cp = (R*g)/(g - 1)

s1 = 0# air.s(T1,p1)
T2 = air.T_s(s=s1,p=p2)
s2 = s1
s3 = air.s(T3,p3)
s4 = s3
T4 = air.T_s(s=s4,p=p4)#
'''
Now we will build isentropic compresion line. 
This is a vertical line. It consists of coordinates (s1,T1), (s1,T2) and (s3,T3), (s3,T4)
API Documentaion
https://steemit.com/utopian-io/@emirfirlar/how-to-build-thermodynamic-cycle-with-pyromat
'''

T = np.linspace(T2,T3)
plt.plot(air.s(T=T,p=p2),T,'r',linewidth=1.5)
T = np.linspace(T1,T4)
plt.plot(air.s(T=T,p=p1),T,'r',linewidth=1.5)
plt.plot([s1,s1],[T1,T2],'r',linewidth=1.5)
plt.plot([s3,s3],[T3,T4],'r',linewidth=1.5)
plt.xlabel('Entropy, s (kJ/kg/K)')
plt.ylabel('Temperature, T (K)')
plt.grid('on')
plt.title('Brayton Cycle T-s Graphic')
plt.show()

