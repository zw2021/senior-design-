import numpy as np
import math # 你好
import matplotlib.pyplot as plt
from ASE366L import beta, Functions as func, TimeConversions as TC, ValladoAlg
from ASE366L.ODE45 import RungeKutta as RK

# Calculate Cart Values using Orb to Cart methods #
# Will be using Newton Raphson to calculate nu over time, then convert Orbital elements with changing nu to Cart Elm
a = 26560 # km
ed = 0.001
i = 50 # deg
BigOmeg = 0 # deg
LitOmeg = 0 # deg
mu = func.Earth.mu

# Given March 1st, 2020, 12:00:00 UTC
tp = TC.CalToJD(2020, 3, 1, 12, 0, 0)   # units: days, converts given time in JD_UTC
t = np.asarray(np.arange(0, 5, 10/(24*60)))   # units: days
n = math.sqrt(mu / (a ** 3))            # units km3/s2
x = np.zeros((len(t), 6))
rIJK = np.zeros((len(t), 3))
vIJK = np.zeros((len(t), 3))
rIJKDim = np.zeros((len(t)))

nu = ValladoAlg.ellip_prop(t, n, ed, tp)  # gets theta in rad

for ii in range(np.size(t)):
    x[ii] = np.array([a, ed, i, LitOmeg, BigOmeg, nu[ii]])  # gets rIJk, vIJk in km
    rIJK[ii], vIJK[ii] = beta.beta(mu, np.squeeze(np.asarray(x[ii])))
    rIJK[ii] = np.squeeze(np.asarray(rIJK[ii]))
    vIJK[ii] = np.squeeze(np.asarray(vIJK[ii])) # Orb to Cart Values
    rIJKDim[ii] = np.linalg.norm(rIJK[ii])  # Orb to Cart Values

# Calculate Cart Values using Orbit Propagator methods #
# Using RungeKutta
OrbitPropvals = RK.RungeKutta(rIJK[1, :], vIJK[1, :], 0, 10/(60*24), 5, mu)
OrbitPropvals = np.delete(OrbitPropvals, 720, axis=0)


# Grab magnitude of r values for propagator
tt = OrbitPropvals[:, 0]*(24*60)
rD = np.sqrt((OrbitPropvals[:, 1]) ** 2 + (OrbitPropvals[:, 2]) ** 2 + (OrbitPropvals[:, 3]) ** 2)  # Position

Error_r = rD - rIJKDim

plt.figure()
plt.suptitle('2 body Position Error with respect to time [numerically propagated minus non perturbed]')
plt.subplot(3, 1, 1)
plt.plot(tt, OrbitPropvals[:, 1], linewidth=.5)#(OrbitPropvals[:, 2] - rIJK[:, 1]))
plt.subplot(3, 1, 2)
plt.plot(tt, OrbitPropvals[:, 2])#(OrbitPropvals[:, 2] - rIJK[:, 1]))
plt.subplot(3, 1, 3)
plt.plot(tt, OrbitPropvals[:, 3])#(OrbitPropvals[:, 2] - rIJK[:, 1]))
plt.ylabel(' Magnitude of position [DU]', size=10)
plt.xlabel('TU [min]', size=16)

plt.show()