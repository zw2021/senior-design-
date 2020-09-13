from ASE366L import ValladoAlg as VA
from ASE366L import Functions as func
import numpy as np
from matplotlib import pyplot as plt

# test data #
JD_TDB = 2451545
mu = func.Sun.mu
r_XYZ, v_XYZ = VA.PlanetRV(mu, JD_TDB)
r_XYZ = np.matmul(VA.J2000toGCRF(r_XYZ), (r_XYZ))
print(r_XYZ)

# Problem 3 Part a#
print('This is Problem 3, Part a:')
JD_TDB = 2457793
mu = func.Earth.mu
r_XYZ, v_XYZ = VA.PlanetRV(mu, JD_TDB)
r_XYZ = np.matmul(VA.J2000toGCRF(r_XYZ), (r_XYZ))
r_XYZ = np.squeeze(np.asarray(r_XYZ))
print('GCRF position of Sun in at 2457794 TDB [AU]',r_XYZ)

# Problem 3 Part b #
JD_TDB = 2457700
times = np.arange(JD_TDB, JD_TDB + 365.25, 1)
mu = func.Earth.mu
r_XYZ = np.zeros((len(times), 3))
v_XYZ = np.zeros((len(times), 3))
r_XYZDim = np.zeros((len(times)))

Q = VA.J2000toGCRF(1)
for ii in range(len(times)):
    r_XYZ[ii], v_XYZ[ii] = VA.PlanetRV(mu, times[ii])
    r_XYZ[ii] = np.matmul(r_XYZ[ii], Q)
    r_XYZ[ii] = np.squeeze(np.asarray(r_XYZ[ii]))
    r_XYZ[ii] = func.KMtoAU(r_XYZ[ii])
    r_XYZDim[ii] = np.linalg.norm((r_XYZ[ii]))

plt.figure()
tt = np.arange(0, 365.25, 1)
plt.suptitle('Position of Sun as a function of time')
plt.subplot(3, 1, 1)
plt.plot(tt, r_XYZ[:, 0])
plt.ylabel('rx [AU]', size=10)

plt.subplot(3, 1, 2)
plt.plot(tt, r_XYZ[:, 1])
plt.ylabel('ry [AU]', size=10)

plt.subplot(3, 1, 3)
plt.plot(tt, r_XYZ[:, 2])
plt.ylabel('rz [AU]', size=10)
plt.xlabel('Time [Days]')

# Problem 1 Part c #
r_Sun_XYZ = np.zeros((len(times), 3))
v_Sun_XYZ = np.zeros((len(times), 3))
delta = np.zeros((len(times), 3))
deltaDim = np.zeros((len(times)))
r_Sun_XYZDim = np.zeros((len(times)))

times = np.arange(JD_TDB, JD_TDB+365.25, 1)

for ii in range(len(times)):
    r_Sun_XYZ[ii] = VA.sun(times[ii])   # units: km , MOD frame
   # print('this is rsun in mod!:',r_Sun_XYZ[ii])
    r_Sun_XYZ[ii] = VA.MODt0GCRF(times[ii], r_Sun_XYZ[ii])  # km, GCRF frame
   # print('this is rsun in gcrf!:',r_Sun_XYZ[ii])
    r_Sun_XYZ[ii] = np.squeeze(np.asarray(r_XYZ[ii]))
    r_Sun_XYZ[ii] = func.KMtoAU(r_Sun_XYZ[ii])  # units AU, GCRF frame
    r_Sun_XYZDim[ii] = np.linalg.norm(r_Sun_XYZ[ii])
    deltaDim[ii] = r_Sun_XYZDim[ii] - r_XYZDim[ii]   # units AU, GCRF frame
    #deltaDim[ii] = np.linalg.norm(delta[ii])

plt.figure()
tt = np.arange(0, 365.25, 1)
plt.suptitle('Difference in Position')
plt.plot(tt, (delta))
plt.ylabel('r [AU]', size=10)

plt.show()
print('End of Computation')




