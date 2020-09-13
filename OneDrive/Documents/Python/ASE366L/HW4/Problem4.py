import numpy as np
import math
from skyfield.api import load
from ASE366L.ValladoAlg import sun, MODt0GCRF
from ASE366L.Functions import deltatimes, DegtoRad
import matplotlib.pyplot as plt

# Test Data #
# r_vec = np.array([26493118, -132755488, -57556547]) # position of the Sun relative to the Earth at 2451545 UTC in GCRF
# tt = 2451545.0
# test = sun(tt)
# #test = MODt0GCRF(tt, test)
# print(test) # AU to km: Value *1.496*1e8

# Problem 4 #
print('This is Problem 4a:')
t =  2451545.0 # UTC
# convert UTC to UT1 #
r_MOD = sun(t)
r_GCRF = MODt0GCRF(t, r_MOD)
print('Position in GCRF Frame [km]:', np.multiply(r_GCRF, 149597870)) # AU to km: Value *1.496*1e8
print('Position in GCRF Frame [AU]:', r_GCRF) # AU to km: Value *1.496*1e8

print('This is Problem 4b')
t = 2457793.5 + 365.25 # UTC
times = np.arange(2457793.5, t, 1)
rr_MOD = np.zeros((len(times), 3))
rr_GCRF = np.zeros((len(times), 3))

for ii in range(len(times)):
    rr_MOD[ii] = sun(times[ii])
    rr_GCRF[ii] = MODt0GCRF(times[ii], rr_MOD[ii])

tt = np.arange(0, 366, 1)

# plot out 4b #
plt.figure()
plt.subplot(3, 1, 1)
plt.suptitle('Position Vector of the Sun')
plt.plot(tt, rr_GCRF[:,0])
plt.ylabel('r_i [AU]', size=10)

plt.subplot(3, 1, 2)
plt.plot(tt, rr_GCRF[:,1])
plt.ylabel('r_j  [AU]', size=10)

plt.subplot(3, 1, 3)
plt.plot(tt, rr_GCRF[:,2])
plt.ylabel('r_k  [AU]', size=10)
plt.xlabel('Time [days]', size=10)
plt.show()

# print('End of Computation')


