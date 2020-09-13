import numpy as np
import math
from ASE366L import beta
from ASE366L.ODE45 import ode_noPert

def ellip_prop_initial(t, n, ed, tp):

    M = n * (t  - tp)

    # Computation for Newton Raphson Method #
    if M < np.pi:
        x_r_old = M + ed / 2
    elif M > np.pi:
        x_r_old = M - ed / 2

    f = lambda x: x - ed * np.sin(x) - M  # Function
    df = lambda x: 1 - ed * np.cos(x)  # Derivative of Function

    # Constraints fo Iteration
    e_a = 1  # Error
    TOL = 10 ** -10  # Max Error Tolerance
    MAX = 100  # Max Iterations (what's a good max?)
    itn = 0  # Initial Iteration

    while (e_a > TOL) & (itn <= MAX):
        itn = itn + 1
        x_r_new = x_r_old - (f(x_r_old)) / df(x_r_old)

        e_a = np.abs(x_r_new - x_r_old)
        x_r_old = x_r_new
        E = x_r_new

    theta = 2 * math.atan2(np.math.sqrt((1 - ed) / (1 + ed)) * np.math.tan(E / 2), np.math.sqrt(1 - ed))
    theta = theta*(180/np.pi)
    return theta

def ellip_prop(t, n, ed, tp):
    theta = np.zeros(len(t))
    M = np.zeros(len(t))
    E = np.zeros(len(t))
    x_r_old = np.zeros(len(t))
    x_r_new = np.zeros(len(t))

    for i in range(0,len(t),1):# range(len(t)):
        M[i] = n*(t[i]-tp)

        # Computation for Newton Raphson Method #
        if M[i] < np.pi:
            x_r_old[i] = M[i] + ed/2
        elif M[i] > np.pi:
            x_r_old[i] = M[i] - ed/2

        f = lambda x: x - ed*np.sin(x) - M[i]    # Function
        df = lambda x: 1 - ed*np.cos(x)          # Derivative of Function

        # Constraints fo Iteration
        e_a = 1             # Error
        TOL = 10**-10       # Max Error Tolerance
        MAX = 100           # Max Iterations (what's a good max?)
        itn = 0             # Initial Iteration

        while (e_a > TOL) & (itn <= MAX):
            itn = itn + 1
            x_r_new[i] = x_r_old[i] - (f(x_r_old[i]))/df(x_r_old[i])

            e_a = np.abs(x_r_new[i] - x_r_old[i])
            x_r_old[i] = x_r_new[i]
            E[i] = x_r_new[i]

        theta[i] = 2*math.atan2(math.sqrt((1+ed))*math.tan(E[i]/2), math.sqrt(1-ed))
    return theta

# Main Function
#########################################################################

def myFunc(a, ed, i, LitOmeg, BigOmeg, mu, n, tf):

    # Problem 1
    print('Problem 1:')
    tp = 0  # sec
    print('Time of Periapsis, [seconds]:', tp)

    t = np.array(0)
    thetaNot = ellip_prop_initial(t, n, ed, tp)
    print('At t = 0 [seconds], Initial True Anomaly [radians]:', thetaNot)

    OrbT0 = np.array([a, ed, i, LitOmeg, BigOmeg, thetaNot])
    CartT0 = beta.beta(mu, OrbT0)
    print('r_IJK:', CartT0[0:3:1])  # Orbital Elements at t = 0
    print('v_IJK:', CartT0[3:6:1])  # Orbital Elements at t = 0

    # Problem 2
    print('Problem 2:')
    t = np.arange(0, tf, 20)

    E0 = 2*math.atan2(np.math.sqrt(1-ed)*np.tan(thetaNot/2), np.math.sqrt(1+ed))
    print('At t = 0, Initial Eccentric Anomaly [radians]:', E0)

    M_not = E0 - ed*np.sin(E0)
    print('At t = 0, Initial Mean Anomaly [radians]:', M_not)

    theta = ellip_prop(t, n, ed, tp)    # gets theta in Rad
    thetaF = theta[12959] # 38879 is the index theta @ tf
    print('At t = tf, True Anomaly [deg]:', thetaF*(180/np.pi))

    EF = 2*math.atan2(np.math.sqrt(1-ed)*np.tan(thetaF/2), np.math.sqrt(1+ed))
    Ef = EF * (180/np.pi) # convert rad to deg
    print('At t = tf, Eccentric Anomaly [deg]:', Ef)

    MF = EF - ed*np.sin(EF)
    Mf = MF * (180/np.pi)   # convert rad to deg
    print('At t = tf, Mean Anomaly [deg]:', Mf)

    # Problem 3
    print('Problem 3:')

    # Compute "Truth Position": Orb to Cart for each iteration
    OrbToCart = np.zeros([len(theta),6])
    OrbElm = np.zeros([len(theta),6])

    for ii in range(np.size(theta)):
        # beta(a, e, i, BigOmeg, LitOmeg, theta)
        OrbElm[ii,:] = np.array([a, ed, i, LitOmeg, BigOmeg, theta[ii]])    # matrix of orbital elements
        OrbToCart[ii,:] = np.array(beta.beta(mu, OrbElm[ii,:]))                 # matrix for convtering orbital elements into Cartesian components

    OrbTf = np.array([a, ed, i, LitOmeg, BigOmeg, theta[12959]])    # Orbital Elements at Tf, indexed at 12959
    CartTf = np.array(beta.beta(mu,OrbTf))                  # Cartesian Elements at Tf

    print('At t = tf, r_IJK [km]:', CartTf[0:3:1])
    print('At t = tf, v_IJK [km/s]:', CartTf[3:6:1])

    # Problem 4
    print('Problem 4: See Output Graphs')

    return  OrbToCart, theta, Mf

