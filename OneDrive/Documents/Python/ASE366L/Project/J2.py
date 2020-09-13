import numpy as np
import math # 你好
import matplotlib.pyplot as plt
from ASE366L import beta, Functions as func, TimeConversions as TC, ValladoAlg
from ASE366L.ODE45 import RungeKutta as RK

# Calculate Cart Values using Orb to Cart methods #
# Will be using Newton Raphson to calculate nu over time, then convert Orbital elements with changing nu to Cart Elm
a = 26560 # km
ed = 0.001
i = 50# deg
LitOmeg = 0 # deg
mu = func.Earth.mu

# Given March 1st, 2020, 12:00:00 UTC
tp = TC.CalToJD(2020, 3, 1, 12, 0, 0)   # units: days, converts given time in JD_UTC
t = np.asarray(np.arange(0, 5, 10/(24*60)))   # units: days
n = math.sqrt(mu / (a ** 3))            # units km3/s2
x = np.zeros((len(t), 6))
epoch = tp
# Propagate nu #
nu = ValladoAlg.ellip_prop(t, n, ed, tp)  # gets theta in rad

# Parameters for computing J2
RE = 6378.1363 # km
J2 = 0.0010826267 # no units
dT = t #np.array((t + tp) - tp)
rD = np.zeros(len(t))
BigOmeg = 0

x = np.zeros((len(t), 6))
rIJK = np.zeros((len(t), 3))
vIJK = np.zeros((len(t), 3))
rIJKDim = np.zeros((len(t)))
for ii in range(np.size(t)):
    x[ii] = np.array([a, ed, i, LitOmeg, BigOmeg, nu[ii]])  # gets rIJk, vIJk in km
    rIJK[ii], vIJK[ii] = beta.beta(mu, np.squeeze(np.asarray(x[ii])))
    rIJK[ii] = np.squeeze(np.asarray(rIJK[ii]))
    vIJK[ii] = np.squeeze(np.asarray(vIJK[ii])) # Orb to Cart Values
    rIJKDim[ii] = np.linalg.norm(rIJK[ii])  # Orb to Cart Values



# x = np.array([a, ed, i, LitOmeg, BigOmeg, nu[1]])  # gets rIJk, vIJk in km
# rIJK, vIJK= beta.beta(mu, np.squeeze(np.asarray(x)))
# rIJK = np.squeeze(np.asarray(rIJK))
# vIJK = np.squeeze(np.asarray(vIJK))
# rMag = np.linalg.norm(rIJK)  # Orb to Cart Values
# Propagate Orbit with J2 #
#def Perturb(r0, v0, T0, dT, tF, mu, muCelest, RE, J2, J3, cD, A_m):
OrbitPropvals_J2 = RK.Perturb(rIJK[0, :], vIJK[0, :], 0, 10/(60*24), 5, epoch, mu, muSun=0, muMoon=0, RE=0, J2=0, J3=-0.0000025327, cD=0, A_m=0.0)
OrbitPropvals_J2 = np.delete(OrbitPropvals_J2, 720, axis=0)
rD_J2 = np.sqrt((OrbitPropvals_J2[:, 1]) ** 2 + (OrbitPropvals_J2[:, 2]) ** 2 + (OrbitPropvals_J2[:, 3]) ** 2)
tt = OrbitPropvals_J2[:, 0]*(24*60)


plt.figure()
plt.suptitle('Position Error with respect to time [J2 minus Non Perturbed]')
plt.subplot(3, 1, 1)
plt.plot(tt*(24*60), rD_J2)#- OrbitPropvals[:, 3]), linewidth=0.8)
# plt.subplot(3, 1, 2)
# plt.plot(tt*(24*60), OrbitPropvals_J2[:, 2])#- OrbitPropvals[:, 3]), linewidth=0.8)
# plt.subplot(3, 1, 3)
# plt.plot(tt*(24*60), OrbitPropvals_J2[:, 3])#- OrbitPropvals[:, 3]), linewidth=0.8)
# plt.ylabel('J2 [DU]', size=10)
# plt.xlabel('TU [min]', size=16)

plt.show()