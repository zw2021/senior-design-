import numpy as np
import math
from ASE366L.ODE45 import CowellsMethod as CM
from ASE366L.Functions import Sun, Moon
import matplotlib.pyplot as plt

# Problem 4, part 1: considering perturbations due to moon
rd0 = 6600          # initial satellite position (given)
rdF = 50000         # Final satellite position (given)

r_SM = np.arange(rd0, rdF, 1)
rVec_Sat = np.zeros((len(r_SM), 3))
aTB = np.zeros((len(r_SM), 3))
aTBD = np.zeros(len(r_SM))

for ii in range(np.size(r_SM)):
    rVec_Earth = np.array([Moon.r, 0, 0])
    rVec_Sat[ii] = np.array([(Moon.r - r_SM[ii]), 0, 0])
    aTB[ii] = CM.Pert(Moon.mu, rVec_Sat[ii], rVec_Earth)
    aTBD[ii] = np.linalg.norm(aTB[ii], 2)

rr_SM = np.arange(rd0, rdF, 1)
rrVec_Sat = np.zeros((len(rr_SM), 3))
aTB_Sun = np.zeros((len(rr_SM), 3))
aTBD_Sun = np.zeros(len(rr_SM))

for jj in range(np.size(r_SM)):
    rVec_Earth_Sun = np.array([Sun.r, 0, 0])
    rrVec_Sat[jj] = np.array([(Sun.r - r_SM[jj]), 0, 0])
    aTB_Sun[jj] = CM.Pert(Sun.mu, rrVec_Sat[jj], rVec_Earth_Sun)
    aTBD_Sun[jj] = np.linalg.norm(aTB_Sun[jj], 2)

plt.figure()
plt.suptitle('Magnitude of Third Body perturbations')
plt.plot(r_SM, aTBD, label='Third Body Pert due to Moon')
plt.plot(r_SM, aTBD_Sun, label = 'Third Body Pert due to Sun')
plt.ylabel('Perturbation due to Sun  [km^2]', size=10)
plt.xlabel('Position r [km]', size=10)
plt.legend(loc="upper left")
plt.show()

print('End Computation')