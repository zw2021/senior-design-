import numpy as np
import math
from ASE366L.ODE45.RungeKutta import RungeKutta

# Problem 5 #
mu = 398600.4415 # units: km^3/s^2
AT = 37 # units: seconds
# Operator A information #
r = np.array([2030, 18638.0, 3707]) # units: km
v = np.array([-4.5386, 0.52049, .10353]) # units: km/s
TT = 2458165.4375 # units: seconds
# Operator B information #
rr = np.array([1715.7, 18671.5, -3713.7])    # units: km
vv = np.array([-4.546, 0.446796, 0.088875])  # units: km/s
UTC_OperatorB = 2458165.4375 # units: seconds
# Time Constraints #
UTC_Initial_OperatorA = TT - AT - 32.184
UTC_Final_OperatorA = 2458165.4375 # seconds
Delta_Time_UTC = UTC_Final_OperatorA - UTC_Initial_OperatorA # seconds
# Start Propagator #
times = np.arange( 0, Delta_Time_UTC, 1 )

dT = 0.01 # seconds
tF = Delta_Time_UTC
T0 = 0 # seconds
dd = RungeKutta( r, v, T0, dT, tF, mu)
print(dd)
print('End Computation')
# They're not the same!!! sadness...