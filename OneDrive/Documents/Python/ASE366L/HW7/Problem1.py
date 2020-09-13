import numpy as np
import math # 你好
import matplotlib.pyplot as plt
from ASE366L import beta, Functions as func, TimeConversions as TC, ValladoAlg
from ASE366L.ODE45 import RungeKutta as RK

RE = 6378.1363 # km
mu_Earth = 398600.4415 # km^3/s^2
J2 = 0.0010826267
J3 = -0.0000025327


#def Perturb(r0, v0, T0, dT, tF, mu, muCelest, RE, J2, J3, cD, A_m):
# OrbitPropvals_J2 = RK.Perturb(rIJK[0, :], vIJK[0, :], 0, 10/(60*24), 5, mu, muSun=0, muMoon=0, RE=0, J2=0, J3=-0.0000025327, cD=0, A_m=0.0)
# OrbitPropvals_J2 = np.delete(OrbitPropvals_J2, 720, axis=0)
# rD_J2 = np.sqrt((OrbitPropvals_J2[:, 1]) ** 2 + (OrbitPropvals_J2[:, 2]) ** 2 + (OrbitPropvals_J2[:, 3]) ** 2)
# tt = OrbitPropvals_J2[:, 0]*(24*60)

# Problem 1
#r = np.array([6092.032, 2487.062, 2487.062])    # test data; units: km
r = np.array([0.009, 4595.737, 4595.731]) # Problem 1 data; units:km
dY = RK.get_pertAcc_with_J2andJ3(r, mu_Earth, RE, J2, J3)
print('This is r [km]', dY)

# Problem 2
#d = RK.get_pertAcc_with_aSRP(2)

r = np.arange(64000, 900000, 10)
#r = np.array([1000000000, 1, 2])   # test
#atm = RK.get_atmDensity(r, RE)     # test
atm = np.zeros(len(r))

for ii in range(len(r)):
    atm[ii] = RK.get_atmDensity(r[ii], RE)

plt.figure()

plt.loglog(r, atm)
plt.savefig("out.png")
plt.title('Density as a Function of r')
plt.ylabel(' Atmospheric Density [kg/ km^3]', size=10)
plt.xlabel('r [km]', size=10)

plt.show()