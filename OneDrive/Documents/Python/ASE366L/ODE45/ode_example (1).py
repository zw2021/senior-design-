#!/usr/local/bin/python
'''

  Example software for propagating an ODE 


  Example software to propagate the position and velocity state for a
  1-D spring-mass system


'''

__author__ = 'Brandon A. Jones'
__version__ = '1'
__date__ = 'March 3, 2016'

################################################################################
#                     I M P O R T     L I B R A R I E S
################################################################################


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import ode

################################################################################
#                     S E T U P    L I B R A R I E S
################################################################################

np.set_printoptions(precision=15)
mpl.rcParams['figure.figsize'] = (8.0, 6.0)


################################################################################
#                   E X P O R T E D    F U N C T I O N S:
################################################################################
def harmoscillator( t, Y, kmratio ) :
  '''
  Here, we have the description of the dynamics to use for propagating
  the state
  '''

  dY = np.zeros( Y.shape )

  dY[0] = Y[1]

  dY[1] = -kmratio*Y[0]

  return dY

################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main() :
  '''
  '''

  def derivFcn( t, x ) :
    return harmoscillator( t, x, kmratio )

  #  Set the (k/m) value used in the propagation
  kmratio = 1.0

  #  As part of the initial solution we need the amplitude of the oscillation
  #  and the phase offset phi
  T0  = 0.0;
  A   = 1.34;
  phi = np.pi/3;

  #  Set the times at which we want the solution.
  times = np.linspace( 0, 20, 101 )
  tF = 20.0
  dT = 0.2

  #  Setup the initial conditions for the harmonic oscillator.  In general,
  #  we can't do this, but here we use the analytic solution.
  x0 = A*np.cos(times[0]*np.sqrt(kmratio)+phi)
  dx0 = -A*np.sqrt(kmratio)*np.sin(np.sqrt(kmratio)*T0+phi)

  #  Group the initial conditions into a single vector
  Y0 = np.array( [ x0, dx0 ] )

  # Define integrator, set tolerances and the initial state
  rv = ode( derivFcn )

  #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
  #  rtol and atol are the relative and absolute tolerances, respectively
  rv.set_integrator('dopri5', rtol=1e-6, atol=1e-6 )
  rv.set_initial_value( Y0, T0)

  # Define output array
  output = []
  output.append(np.insert(Y0, 0, T0))

  # Run the integrator and populate output array with positions and velocities
  while rv.successful() and rv.t < tF:
    rv.integrate( rv.t + dT )
    output.append( np.insert( rv.y, 0, rv.t ) )

  #  Convert the output a numpy array for later use
  output = np.array(output)

  #  Let's plot the result

  plt.figure()

  plt.subplot( 2, 1, 1 )
  plt.plot( output[:,0], output[:,1] )
  plt.xlim([ 0, 20 ] )
  plt.ylabel( 'DU', size=16 )

  plt.subplot( 2, 1, 2 )
  plt.plot( output[:,0], output[:,2] )
  plt.xlim([ 0, 20 ] )
  plt.ylabel( 'DU/TU', size=16 )
  plt.xlabel( 'TU', size=16 )

  plt.show()

  return


if __name__ == "__main__":
  main()
