import numpy as np
import math
from ASE366L.ODE45 import RungeKutta as RK, CowellsMethod as CM
from ASE366L import alpha, beta
from ASE366L.ValladoAlg import sun, MODt0GCRF
from ASE366L.Functions import Sun
from matplotlib import pyplot as plt

# Given Conditions of satellite's position in Orbital Elements #
mu = 398600 # km3/s2
a = 500000 # km
ed = 0.3
i = 50 # deg
BigOmeg = 350
LitOmeg = 90
nu = 0

x = np.array([a, ed, i, LitOmeg, BigOmeg, nu])  # matrix of orbital elements
a, b  = beta.beta(mu, x)    # initial position of the satellite = a; initial velocity of the satellite = b
r0 = np.array([ 39066.61396713, 221557.77754582, 268115.55509164])
v0 = np.array([-1.19827848e+00, 2.11288827e-01, 5.70743697e-17])

T0 = 2451545 # TAI add 37 to convert from UTC to TAI
tF = (365.25) * 86400 # Time Rep: UTC
dT = (.5)* 86400
T0 = 0
#d = rungekuttapert(r0, v0, T0, dT, tF, mu, muCelest, JD_TAI, )
Cart = RK.RungeKuttaPert(r0, v0, T0, dT, tF, mu, Sun.mu)

r_Cart = np.transpose(np.matrix([Cart[:,1], Cart[:,2], Cart[:,3]]))
v_Cart = np.transpose(np.matrix([Cart[:,4], Cart[:,5], Cart[:,6]]))
OrbElm = np.zeros((len(r_Cart),6))

h = np.cross(r_Cart, v_Cart)
for ii in range(len(r_Cart)):
    #print(type(r_Cart[ii, :]))
    OrbElm[ii] = alpha.cto(mu, np.squeeze(np.asarray(r_Cart[ii, :])), np.squeeze(np.asarray(v_Cart[ii, :])))

try:
    plt.figure()

    plt.suptitle('Changes in Angular Momentum')
    plt.subplot(2, 1, 1)
    plt.plot(np.arange(T0, tF, dT), h[:, 0])
    plt.ylabel('h_i [m2/s]', size=10)

    plt.subplot(2, 1, 2)
    plt.plot(np.arange(T0, tF, dT), h[:, 1])
    plt.ylabel('h_j [m2/s]', size=10)
    plt.xlabel('Time [days]', size=10)

    # plt.suptitle('Orbital Elements')
    # plt.subplot(5, 1, 1)
    # plt.plot(np.arange(T0, tF, dT), OrbElm[:, 0])
    # plt.ylabel('Semi Major Axis [km]', size=5)
    #
    # plt.subplot(5, 1, 2)
    # plt.plot(np.arange(T0, tF, dT), OrbElm[:, 1])
    # plt.ylabel('Eccentricity', size=5)
    #
    # plt.subplot(5, 1, 3)
    # plt.plot(np.arange(T0, tF, dT), OrbElm[:, 2])
    # plt.ylabel('Inclination [Deg]', size=5)
    #
    # plt.subplot(5, 1, 4)
    # plt.plot(np.arange(T0, tF, dT), OrbElm[:, 3])
    # plt.ylabel('RAAN', size=5)
    #
    # plt.subplot(5, 1, 5)
    # plt.plot(np.arange(T0, tF, dT), OrbElm[:, 4])
    # plt.ylabel('Argument of Periapsis [Deg]', size=5)
    # plt.xlabel('Time [days]', size=10)

    plt.show()
except Exception as e:
    print(e)

print('End of Computation')