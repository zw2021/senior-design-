#!/usr/local/bin/python
'''

  software for propagating an ODE for position and velocity state of a
  3-D orbit

  Be sure no arguments are loaded in alpha so that the compiler will not
  become confused
'''


################################################################################
#                     I M P O R T     L I B R A R I E S
################################################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import ode
import math
from ASE366L import alpha
from ASE366L.HW2.problem1 import myFunc

################################################################################
#                     S E T U P    L I B R A R I E S
################################################################################

np.set_printoptions(precision=15)
mpl.rcParams['figure.figsize'] = (8.0, 6.0)

################################################################################
#                   E X P O R T E D    F U N C T I O N S:
################################################################################
def OrbitProp( t, Y, mu ) :
  '''
  Here, we have the description of the dynamics to use for propagating
  the state
  '''

  dY = np.zeros([6,1])
  r = np.math.sqrt(Y[0] ** 2 + Y[1] ** 2 + Y[2] ** 2)
  dY[0] = Y[3]
  dY[1] = Y[4]
  dY[2] = Y[5]

  dY[3] = -(mu / r ** 3) * Y[0]  # velocity/ acceleration
  dY[4] = -(mu / r ** 3) * Y[1]
  dY[5] = -(mu / r ** 3) * Y[2]

  return dY

################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main() :

  a = 26000  # km
  ed = 0.72  # less than one, so we know that we can use elliptical function
  i = 75
  BigOmeg = 90
  LitOmeg = -90

  tf = 86400 * 3
  #t = np.arange(0, tf, 20)
  #tt = np.range(0, tf, 20)
  mu = 398600.4415  # km^3/s^2
  n = math.sqrt(mu / (a ** 3))
  OrbToCart, theta, Mf = myFunc(a, ed, i, LitOmeg, BigOmeg, mu, n, tf)

  #  Set the mu value used in the propagation
  mu = 398600.4415

  #  Setup the initial conditions for Orbit Propagator: Hard Coded from Results
  T0  = 0.0
  r0 = np.array([1.884202648346351e+03,  3.303972977775553e-13, - 7.031940015384417e+03])
  dr0 = np.array([4.404255738577948e-16, 9.704371673898565e+00, 5.739737826996621e-16])
  #  Group the initial conditions into a single vector
  Y0 = np.concatenate( [ r0, dr0 ], axis=0 )

  #  Set the times at which we want the solution.
  T = 3*(86400)
  times = np.arange( 0, 3*T, 20 )
  tF = T
  dT = 20

  # Define integrator, set tolerances and the initial state
  def derivFcn(times,Y):
    return OrbitProp(times, Y, mu )

  rv = ode( derivFcn )

  #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
  #  rtol and atol are the relative and absolute tolerances, respectively
  rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20 )
  rv.set_initial_value(Y0, T0)

  # Define output array
  output = []
  output.append(np.insert(Y0, 0, T0))

  # Run the integrator and populate output array with positions and velocities
  while rv.successful() and rv.t < tF:
    rv.integrate( rv.t + dT )
    output.append( np.insert( rv.y, 0, rv.t ) )

  output = np.array(output)

  # Convert the output a numpy array for later use
  output = np.array(output)
  rD = np.sqrt((output[:, 1])**2 + (output[:, 2])**2 + (output[:,3])**2)  # Position
  vD = np.sqrt((output[:, 4])**2 + (output[:, 5])**2 + (output[:,6])**2)  # Velocity

  Et = .5 * ((vD**2)) - mu / (rD)
  Et0 = .5*((np.linalg.norm(dr0))**2) - mu/(np.linalg.norm(r0))
  Z = Et - Et0

  #  Let's plot the results:
  # Plot Change in Energy

  plt.figure()
  plt.suptitle('Change in Energy With Respect to Time')
  plt.plot((output[:, 0]) / 86400, Z)  # output[:,0] is the time time
  plt.xlabel('TU')
  plt.ylabel('MU*DU^2/TU^2', size=16)

  plt.show()

# Problem 5
  ############################################################
  output = np.delete(output, 12960, 0)

  # grab the r and v vector values only (probably a more elegant solution - develop later)
  OrbPropVals = np.transpose(np.array([output[:, 1], output[:, 2], output[:, 3], output[:, 4], output[:, 5], output[:, 6]]))

  print('Problem 5:')

  # Error = Propagated Orbital Values - Values Converted from Orb to Cart
  Error_rv = OrbPropVals[12959,:] - OrbToCart[12959,:]
  #print(Error_rv)

  print('Error in r at Tf:', Error_rv[0:3:1], 'Magnitude of Error Vector in R', np.linalg.norm(Error_rv[0:3:1],2))
  print('Error in v at Tf:', Error_rv[3:6:1], 'Magnitude of Error Vector in V', np.linalg.norm(Error_rv[3:6:1],2))


# Problem 6
####################################################################33

  #error = np.zeros([len(OrbToCart),6])
  #for ii in range(len(OrbPropVals)):
    #error[:,ii] = abs(OrbPropVals[:,ii]) - abs(OrbToCart[:,ii]) # initial position, r0

  error_rx = np.array((OrbPropVals[:,0]) - (OrbToCart[:,0])) # initial position, r0
  error_ry = np.array((OrbPropVals[:,1]) - (OrbToCart[:,1])) # initial position, r0
  error_rz = np.array((OrbPropVals[:,2]) - (OrbToCart[:,2])) # initial position, r0
  error_vx = np.array((OrbPropVals[:,3]) - (OrbToCart[:,3])) # initial position, r0
  error_vy = np.array((OrbPropVals[:,4]) - (OrbToCart[:,4])) # initial position, r0
  error_vz = np.array((OrbPropVals[:,5]) - (OrbToCart[:,5])) # initial position, r0

  # Plot error components in position
  plt.figure()
  plt.suptitle('Position Error with respect to time [numerically propagated minus truth]')
  plt.subplot( 3, 1, 1 )
  plt.plot( (output[:,0])/3600, np.transpose(error_rx))
  plt.ylabel( 'X Component [DU]', size=10)
  plt.xlabel( 'TU', size=16 )

  plt.subplot( 3, 1, 2 )
  plt.plot( (output[:,0])/3600, np.transpose(error_ry))
  plt.ylabel( 'Y Component [DU]', size=10 )
  plt.xlabel( 'TU', size=16 )

  plt.subplot( 3, 1, 3 )
  plt.plot( (output[:,0])/3600, np.transpose(error_rz))
  plt.ylabel( 'Z Component [DU]', size=10 )
  plt.xlabel( 'TU', size=16 )

  plt.show()

  #  Plot error components in velocity
  plt.figure()
  plt.suptitle('Velocity Error with respect to time [numerically propagated minus truth]')
  plt.subplot(3, 1, 1)
  plt.plot((output[:, 0]) / 3600, np.transpose(error_vx))
  plt.ylabel('X Component [DU/TU]', size=10)
  plt.xlabel('TU', size=16)

  plt.subplot(3, 1, 2)
  plt.plot((output[:, 0]) / 3600, np.transpose(error_vy))
  plt.ylabel('Y Component [DU/TU]', size=10)
  plt.xlabel('TU', size=16)

  plt.subplot(3, 1, 3)
  plt.plot((output[:, 0]) / 3600, np.transpose(error_vz))
  plt.ylabel('Z Component [DU/TU]', size=10)
  plt.xlabel('TU', size=16)

  plt.show()
  return

if __name__ == "__main__":
  main()
