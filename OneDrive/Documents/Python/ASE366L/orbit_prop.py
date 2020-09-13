#this software uses Newton Raphson to propagate an orbit over time#
import numpy as np
import math
import matplotlib.pyplot as plt

# Propagator Function to Calculate Mean Anomaly using Newton Raphson Method for Elliptical Orbits
##################################################################################################

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

        theta[i] = 2*math.atan2(np.math.sqrt((1-ed)/(1+ed))*np.math.tan(E[i]/2), np.math.sqrt(1-ed))
    return theta

# Propagator Function to Calculate Mean Anomaly using Newton Raphson Method for Hyperbolic Orbits
##################################################################################################
def hyper_prop(t, n, ed, tp):
    theta = np.zeros(len(t))
    M = np.zeros(len(t))
    F = np.zeros(len(t))
    x_r_old = np.zeros(len(t))
    x_r_new = np.zeros(len(t))

    for i in range(0, len(t), 1):  # range(len(t)):
        M[i] = n * (t[i] - tp)

        # Computation for Newton Raphson Method #
        if M[i] < np.pi:
            x_r_old[i] = M[i] + ed / 2
        elif M[i] > np.pi:
            x_r_old[i] = M[i] - ed / 2

        f = lambda x: ed*math.sinh(x) - x - M[i]    # Anonymous Functions to Compute NR
        df = lambda x: ed*math.cosh(x) - 1

        # Constraints fo Iteration
        e_a = 1             # Error
        TOL = 10**-10       # Max Error Tolerance
        MAX = 100           # Max Iterations (what's a good max?)
        itn = 0             # Initial Iteration

        while (e_a > TOL) & (itn <= MAX):
            itn = itn + 1
            x_r_new[i] = x_r_old[i] - (f(x_r_old[i])) / df(x_r_old[i])

            e_a = np.abs(x_r_new[i] - x_r_old[i])
            x_r_old[i] = x_r_new[i]
            F[i] = x_r_new[i]

        F[i] = x_r_new[i]
        theta[i] = 2*math.atan2((np.math.sqrt(1+ed))*np.tanh(F[i]/2), np.math.sqrt(ed-1))

    return theta

# Main Function
##################################################################################################
#r = [0, 39991, -52979] # these r and v arrays show a hyperbolic orbit
#v = [0, 0, 3.5511]
#t = range(-20000, 20000, 500)
# = [-17130, 39, -784] # these r and v arrays show an eccentric orbit
#v = [-0.251, -2.827, 3.7540]
#t = range(0, 40000, 500)

mu = 398600.4415
RE = 6378.1363

# Calculate eccentricity vector and magnitude
h = np.cross(r,v)
hd = np.linalg.norm(h)
rd = np.linalg.norm(r)
e = (np.cross(v,h)/mu) - r/rd
ed = np.linalg.norm(e)

# Calculate True Anomaly
z = (1/ed)*((hd**2/(rd*mu))-1)
theta = np.arccos(z)
theta = np.real(theta)
if np.dot(r,v)<0:
    theta = 2*np.pi-theta

if ed < 1:
    print('Eccentricity, e shows an elliptical orbit:', ed)
    a = rd*(1+ed*np.cos(theta))/(1-ed**2)
    n = np.math.sqrt(mu/a**3)

    E0 = 2*math.atan2(np.math.sqrt(1-ed)*np.tan(theta/2), np.math.sqrt(1+ed))
    print('Eccentric Anomaly at t = 0, E0 [radians]:', E0)

    M_not = E0 - ed*np.sin(E0)
    print('At t = 0, M_0 [radians]:', M_not)

    tp = -1*M_not/n
    print('Time of Periapsis, [seconds]:', tp)
    xx = ellip_prop(t, n, ed, tp)
    #print(theta[20])
    print('At t = ___ [seconds], True Anomaly is [radians]', xx)

elif ed > 1:
    print('Eccentricity, e shows a hyperbolic orbit:', ed)

    F0 = 2*math.atanh((np.math.sqrt((ed-1)/(ed+1)))*np.tan(theta/2))
    print('Eccentric Anomaly at t = 0, F0 [radians]:', F0)

    M_not = -F0 + ed*np.sinh(F0)
    print('At t = 0, M_0 [radians]:', M_not)

    a = rd*(1+ed*np.cos(theta))/(1-ed**2)
    n = np.math.sqrt(mu/np.abs(a**3))
    tp = -M_not/n

    print('Time of Periapsis, [seconds]:', tp)
    yy = hyper_prop(t, n, ed, tp)

    print('At t = ___ [seconds], True Anomaly is [radians]', yy)

# print(' At t [seconds], True Anomaly [radians]:', theta) # change as question asks