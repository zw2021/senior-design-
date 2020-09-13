from skyfield.api import load
from ASE366L.HW3.scratch import ERA, JDtoUTT
from ASE366L.Functions import s_prime_Generator, WRPN_MatrixGenerator, R3Generator
import numpy as np

# see website for documentation: https://rhodesmill.org/skyfield/time.html
ts = load.timescale() # external function to convert JD to UT1
t1 = ts.utc(2020, 2, 14, 21)
t1_JD_UT1 = t1.ut1
# t1_JD_UT1 = t1._utc_float() - converts given term to utc

# Problem 6 Part a
print('Problem 6:')
print('This is JD_UT1 [seconds]:', t1_JD_UT1)

# Problem 6 Part b
theta_ERA = ERA(t1_JD_UT1)
#theta_ERA = 99.35274668*np.pi/180
print('This is theta_ERA [radians]:', theta_ERA)

# Problem 6 Part c
xp = 0.41075 # units: arc seconds
yp = 0.377800 # units: arc seconds
X = 396.8722 # units: arc seconds
Y = -1.481681 # units: arc seconds
theta_ERA = theta_ERA # Calculated Previously from Above
TUTT = JDtoUTT(t1.tt)
s_prime = s_prime_Generator(TUTT)
s = -0.000904
W,R,PN =  WRPN_MatrixGenerator(s, s_prime, xp, yp, X, Y, theta_ERA)
print('This is W:\n', W)
print('This is R:\n', R)
print('This is PN:\n', PN)

Q_ITRF_GCRF = np.matmul(PN, np.matmul(R,W))
print('This Q_ITRF_GCRF:\n', Q_ITRF_GCRF)


# Problem 7
print('Problem 7:')
r_ITRF = np.array([-742.845, -5463.244, 3196.066])

# Part a
GMST = 99.35274668
GMST = -1* GMST *(np.pi/180)
r_GCRF = np.matmul(R3Generator(GMST), np.transpose(r_ITRF))
print('This is r_GCRF using R3(GMST) [km]:', r_GCRF)

r_GCRF_better = np.matmul( Q_ITRF_GCRF, np.transpose(r_ITRF) )
print('This is r_GCRF using Q_ITRF_GCRF [km]:', r_GCRF_better)

Vec_error = np.array([r_GCRF_better - r_GCRF])
print('This is the error difference between the two previous vectors [km]:', Vec_error)
print('This is the magnitude of ther error:', np.linalg.norm(Vec_error))
print('End Computation')