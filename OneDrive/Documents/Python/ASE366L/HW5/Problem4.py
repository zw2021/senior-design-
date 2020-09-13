from ASE366L import ValladoAlg as VA
from ASE366L import TimeConversions as TC
from ASE366L import Functions as func
import numpy as np
from matplotlib import pyplot as plt
from ASE366L import beta, alpha
from ASE366L.ODE45 import RungeKutta as RK
from mpl_toolkits.mplot3d import Axes3D

a = 70000   # km
e = 0.01
i = 20 # deg
BigOmeg = 180 # deg
LitOmeg = 25 # deg
nu = 0 # deg
mu = func.Earth.mu
muCelest = func.Sun.mu
x = np.array([a, e, i, LitOmeg, BigOmeg, nu])  # matrix of orbital elements
rIJK, vIJK = beta.beta(mu, x)
rIJK = np.squeeze(np.asarray(rIJK)); vIJK = np.squeeze(np.asarray(vIJK))
JD_TDB = 2451545
times = np.arange(JD_TDB, 2451545 + 365.25, .5)
r_Sat_GCRF = np.zeros((len(times), 6))

epoch = 2451545 # In time TDB
T0 = 0
tF = (365.25) *86400# Time Rep: UTC
dT = (.5)*86400
r0 = rIJK
v0 = vIJK
time_series = np.arange(T0, tF, dT)
print(np.shape(time_series))
Cart = RK.RungeKuttaPertRV(r0, v0, T0, dT, tF, mu, muCelest, epoch)  # j2000 frame, units km
print(np.shape(Cart))
times = np.array(Cart[:,0])
r_Cart = np.transpose(([Cart[:,1], Cart[:,2], Cart[:,3]]))
v_Cart = np.transpose(([Cart[:,4], Cart[:,5], Cart[:,6]]))
OrbElm = np.zeros((len(r_Cart), 6))

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot3D(r_Cart[:, 0], r_Cart[:, 1], r_Cart[:,2])
# plt.show()

print(mu)
for ii in range(len(r_Cart)):
   # print(type(r_Cart[ii, :]))
   # r_Cart[ii] = np.matmul(r_Cart[ii], Q)
   # r_Cart[ii] = np.squeeze(np.asarray(r_Cart[ii]))
   # r_Cart[ii] = func.KMtoAU(r_Cart[ii])
   OrbElm[ii] = alpha.cto(mu, np.squeeze(np.asarray(r_Cart[ii])), np.squeeze(np.asarray(v_Cart[ii])))

plt.figure()
plt.suptitle('r')
plt.subplot(3, 1, 1)
plt.plot(time_series, r_Cart[:, 0])
plt.subplot(3, 1, 2)
plt.plot(time_series, r_Cart[:, 1])
plt.subplot(3, 1, 3)
plt.plot(time_series, r_Cart[:, 2])

plt.suptitle('v')
plt.subplot(3, 1, 1)
plt.plot(time_series, v_Cart[:, 0])
plt.subplot(3, 1, 2)
plt.plot(time_series, v_Cart[:, 1])
plt.subplot(3, 1, 3)
plt.plot(time_series, v_Cart[:, 2])
plt.show()

plt.suptitle('Orbital Elements: a, e, i')
plt.subplot(3, 1, 1)
plt.plot(time_series, OrbElm[:, 0])
plt.ylabel('Semi Major Axis [km]', size=5)

plt.subplot(3, 1, 2)
plt.plot(time_series, OrbElm[:, 1])
plt.ylabel('Eccentricity', size=5)

plt.subplot(3, 1, 3)
plt.plot(time_series, OrbElm[:, 2])
plt.ylabel('Inclination [Deg]', size=5)

plt.figure()
plt.suptitle('Orbital Elements: w and RAAN')
plt.subplot(2, 1, 1)
plt.plot(time_series, OrbElm[:, 3])
plt.ylabel('RAAN', size=5)

plt.subplot(2, 1, 2)
plt.plot(time_series, OrbElm[:, 4])
plt.ylabel('Argument of Periapsis [Deg]', size=5)
plt.xlabel('Time [days]', size=10)

plt.show()
print('End of Computation')
