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
  '''
  '''

  #  Set the mu value used in the propagation
  mu = 398600.4415

  #  As part of the initial solution we need the amplitude of the oscillation
  #  and the phase offset phi
  T0  = 0.0
  RE = 6378.1363
  J2 = .00108248

  #  Setup the initial conditions for Orbit Propagator

  #r0 = np.array([-5282.628, -4217.988, 296.511])
  #dr0 = np.array([-4.541972, 5.885228, 2.043105])
  r0 = np.array([1.884202648346351e+03,  3.303972977775553e-13, - 7.031940015384417e+03])
  dr0 = np.array([4.404255738577948e-16, 9.704371673898565e+00, 5.739737826996621e-16])

  #  Group the initial conditions into a single vector
  Y0 = np.concatenate( [ r0, dr0 ], axis=0 )

  # Convert from Cartesian to Orbital Elements using Alpha Function
  xx = alpha.cto(mu, r0, dr0)
  sma = xx[0] # semi major axis

  #  Set the times at which we want the solution.
  #T = 2*np.pi*math.sqrt(a**3/mu)
  T = 3*(86400)
  times = np.arange( 0, 3*T, 20 )
  tF = 3*T
  dT = 20

  # Define integrator, set tolerances and the initial state
  def derivFcn(t,Y):
    return OrbitProp(t, Y, mu )

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

  #  Convert the output a numpy array for later use
  output = np.array(output)
  rD = np.sqrt((output[:, 1])**2 + (output[:, 2])**2 + (output[:,3])**2)  # Position
  vD = np.sqrt((output[:, 4])**2 + (output[:, 5])**2 + (output[:,6])**2)  # Velocity

  # Find Acceleration #
  a = np.zeros((len(rD),6))
  tt = np.column_stack((output[:,0]))
  YY = np.column_stack((output[:,1], output[:,2], output[:,3], output[:,4], output[:,5], output[:,6]))

  for i in range(0, len(rD)):#range(np.size(a)):
    #a[i] = OrbitProp(tt[i], YY[i],  mu)
    a[i] = -(mu / rD[i]** 3) * rD[i]

  R = np.column_stack((output[:,1], output[:,2], output[:,3]))
  V = np.column_stack((output[:,4], output[:,5], output[:,6]))

  # Calculate Change in Energy
  Et = .5*(vD**2) - mu/(rD)
  Et0 = .5*((np.linalg.norm(dr0))**2) - mu/(np.linalg.norm(r0))
  Z = Et - Et0

  # Calculate Orbital Elements
  dd = np.zeros((len(rD),6))
  for i in range(len(rD)):
    dd[i,:]= alpha.cto(mu, R[i,:], V[i,:])

  #  Let's plot the result

  # Plot Position
  plt.figure()
  plt.suptitle('Position, Velocity, Acceleration as Function of Time')
  plt.subplot( 3, 1, 1 )
  plt.plot((output[:,0])/3600, rD)   # output[:,0] is the time time
  plt.xlim([ 0, 5 ] )
  plt.ylabel( 'DU', size=16 )

  # Plot Velocity
  plt.subplot( 3, 1, 2 )
  plt.plot( (output[:,0])/3600, vD )
  plt.xlim([ 0, 5 ] )
  plt.ylabel( 'DU/TU', size=16 )
  plt.xlabel( 'TU', size=16 )

  # Plot Acceleration

  plt.subplot(3, 1, 3)
  plt.plot((output[:, 0]) / 3600, a)
  plt.xlim([0, 5])
  plt.ylabel('DU/TU^2', size=16)
  plt.xlabel('TU', size=16)

  plt.show()

  # Plot Change in Energy
  plt.figure()
  plt.suptitle('Change in Energy With Respect to Time')
  plt.plot((output[:, 0]) / 86400, Z)  # output[:,0] is the time time
  plt.xlim([0, 0.2])
  plt.ylabel('MU*DU^2/TU^2', size=16)

  plt.show()

  # Plot Orbital Elements
  plt.figure()
  plt.suptitle('Orbital Elements with respect to time')
  plt.subplot( 6, 1, 1 )
  plt.plot((output[:,0])/3600, dd[:,0])   # output[:,0] is the time time
  plt.xlim([ 0, 5 ] )
  plt.ylabel( 'Semi Major Axis [DU]', size=6 )

  plt.subplot(6, 1, 2)
  plt.plot((output[:, 0]) / 3600, dd[:, 1])  # output[:,0] is the time time
  plt.xlim([0, 5])
  plt.ylabel('Eccentricity', size=6)

  plt.subplot(6, 1, 3)
  plt.plot((output[:, 0]) / 3600, dd[:, 2])  # output[:,0] is the time time
  plt.xlim([0, 5])
  plt.ylabel('Inclination [rad]', size=6)

  plt.subplot(6, 1, 4)
  plt.plot((output[:, 0]) / 3600, dd[:, 3])  # output[:,0] is the time time
  plt.xlim([0, 5])
  plt.ylabel('Arg of Periapsis [rad]', size=6)

  plt.subplot(6, 1, 5)
  plt.plot((output[:, 0]) / 3600, dd[:, 4])  # output[:,0] is the time time
  plt.xlim([0, 5])
  plt.ylabel('RAAN [rad]', size=6)

  plt.subplot(6, 1, 6)
  plt.plot((output[:, 0]) / 3600, dd[:, 5])  # output[:,0] is the time time
  plt.xlim([0, 5])
  plt.ylabel('True Anomaly [rad]', size=6)

  plt.show()


  return

if __name__ == "__main__":
  main()
