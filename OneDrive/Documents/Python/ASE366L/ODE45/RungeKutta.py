import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import ode
import math
from ASE366L import ValladoAlg

def RungeKutta(r0, v0, T0, dT, tF, mu):
    Y0 = np.concatenate([r0, v0], axis=0)
    t = np.arange(T0,tF, dT)

    # Orbit Propagator that ignores Perturbations
    def OrbitProp( t, Y, mu ):

      dY = np.zeros([6,1])
      r = np.math.sqrt(Y[0] ** 2 + Y[1] ** 2 + Y[2] ** 2)
      dY[0] = Y[3]
      dY[1] = Y[4]
      dY[2] = Y[5]

      dY[3] = -(mu / r ** 3) * Y[0]  # velocity/ acceleration
      dY[4] = -(mu / r ** 3) * Y[1]
      dY[5] = -(mu / r ** 3) * Y[2]

      return dY

    def derivFcn(t, Y):
        return OrbitProp(t, Y, mu)

    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y0, T0)

    # Define output array
    output = []
    output.append(np.insert(Y0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    #  Convert the output a numpy array for later use
    output = np.array(output)
    rD = np.sqrt((output[:, 1]) ** 2 + (output[:, 2]) ** 2 + (output[:, 3]) ** 2)  # Position
    vD = np.sqrt((output[:, 4]) ** 2 + (output[:, 5]) ** 2 + (output[:, 6]) ** 2)  # Velocity

    return output

# Orbit Propagator that accounts for 2Body, J2, J3, Atmospheric Drag, SRP
def Perturb(r0, v0, T0, dT, tF, epoch, mu, muSun, muMoon, RE, J2, J3, cD, A_m):
    Y0 = np.concatenate([r0, v0], axis=0)
    t = np.arange(T0, tF, dT)
    dY = np.zeros([6, 1])
    vRel = np.zeros([3, 1])

    w = np.asarray([0, 0, 0.0000729211585530]) # See Vallado, very last page unit: rad/ solar sec
    # Orbit Propagator that ignores Perturbations
    def OrbitProp( t, epoch, Y, mu, muSun, muMoon, RE, J2, J3, cD, A_m):
        #   UNCOMMENT IF TOLD TO PROPAGATE IN TAI   #
        # Delta_UT1 = 0
        # Leap_Seconds = 32
        # JD_UT1 = JD_TAI - Leap_Seconds + Delta_UT1 + t/86400

        r = np.math.sqrt(Y[0] ** 2 + Y[1] ** 2 + Y[2] ** 2)

        # Account for 3rd Body Perturbations due to the Sun Using RV Code#
        SunPos, SunVeloc = (ValladoAlg.PlanetRV(muSun, epoch + t/86400))# J200 Frame, units: km

        # UNCOMMENT IF YOU USE SUN FUNCTION TO LOCATE POSITION OF SUN #
        #SunPos = np.multiply(ValladoAlg.sun(t/86400 + T0), 149597870)
        # Q = ValladoAlg.J2000toGCRF(1)
        # SunPos = np.matmul(Q, np.transpose(SunPos))
        # SunPos = np.squeeze(np.asarray(SunPos))

        #Sunpos = np.linalg.norm(Sunpos, 2)
        r_sat_sun = SunPos - Y[0:3]
        r_sat_sun = np.linalg.norm(r_sat_sun, 2)

        ax_Sun = muSun*(((SunPos[0]- Y[0]) / (r_sat_sun ** 3)) - SunPos[0]/(np.linalg.norm(SunPos, 2))**3)
        ay_Sun = muSun*(((SunPos[1]- Y[1]) / (r_sat_sun ** 3)) - SunPos[1]/(np.linalg.norm(SunPos, 2))**3)
        az_Sun = muSun*(((SunPos[2]- Y[2]) / (r_sat_sun** 3)) - SunPos[2]/(np.linalg.norm(SunPos, 2))**3)

        # Account for 3rd Body Perturbations due to the moon #
        MoonPos = ValladoAlg.getMoonPos_wrt_Earth(t/86400 + T0) # conver tht time to UT_TDB!!!

        r_sat_moon = MoonPos - Y[0:3]
        r_sat_moon = np.linalg.norm(r_sat_moon, 2)

        ax_Moon = muMoon* (((MoonPos[0] - Y[0]) / (r_sat_moon ** 3)) - MoonPos[0] / (np.linalg.norm(MoonPos, 2)) ** 3)
        ax_Moon = muMoon* (((MoonPos[1] - Y[1]) / (r_sat_moon ** 3)) - MoonPos[1] / (np.linalg.norm(MoonPos, 2)) ** 3)
        ax_Moon = muMoon* (((MoonPos[2] - Y[2]) / (r_sat_moon ** 3)) - MoonPos[2] / (np.linalg.norm(MoonPos, 2)) ** 3)

        # Calculate V relative for atmospheric drag computation #
        temp = np.asarray([Y[0], Y[1], Y[2]])
        var = np.cross(w, temp)
        vRel[0] = Y[0] - var[0]
        vRel[1] = Y[1] - var[1]
        vRel[2] = Y[2] - var[2]
        vRelMag = np.math.sqrt(vRel[0] ** 2 + vRel[1] ** 2 + vRel[2] ** 2)

        # Calculate atmospheric drag #
            # CHANGE ONCE YOU FIND OUT HOW TO FIND RHO BASED ON THE ALTITUDE - MEO IS GREATER THAN 1000...
        rho = get_atmDensity(r, RE)
        ax_Drag = -1 / 2 * (cD * A_m) * rho * vRelMag * vRel[0]
        ay_Drag = -1 / 2 * (cD * A_m) * rho * vRelMag * vRel[1]
        az_Drag = -1 / 2 * (cD * A_m) * rho * vRelMag * vRel[2]

        dY[0] = Y[3]
        dY[1] = Y[4]
        dY[2] = Y[5]

        # Cowell's Formulation; includes J2 and J3 Perturbed Acceleration
        dY[3] = -(mu / r ** 3) * Y[0] + ((-3 *J2 *mu *RE**2 *Y[0]) / (2 *r**5)) *(1 - ((5 *Y[2]**2) / (r **2)))\
                + ((-5 *J3 *mu *RE**3 *Y[0]) / (2 *r**7)) *(3 *Y[2] - ((7 *Y[2]**3) / (r **2))) \
                + ax_Drag + ax_Sun

        dY[4] = -(mu / r ** 3) * Y[1] + ((-3 *J2 *mu *RE**2 *Y[1]) / (2 *r**5)) *(1 - ((5 *Y[2]**2) / (r **2)))\
                + ((-5 *J3 *mu *RE**3 *Y[1]) / (2 *r**7)) *(3 *Y[2] - ((7 *Y[2]**3) / (r **2))) \
                + ay_Drag + ay_Sun

        dY[5] = -(mu / r ** 3) * Y[2] + ((-3 *J2 *mu *RE**2 *Y[2]) / (2 *r**5)) *(3 - ((5 *Y[2]**2) / (r **2)))\
                + ((-5 *J3 *mu *RE**3 ) / (2 *r**7)) *(6 *Y[2]**2 - ((7 *Y[2]**4) / (r **2)) - ((3 *r**2) / (5))) \
                + az_Drag + ay_Sun

        return dY

    def derivFcn(t, Y):
        return OrbitProp(t, epoch, Y, mu, muSun, muMoon, RE, J2, J3, cD, A_m)

    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y0, T0)

    # Define output array
    output = []
    output.append(np.insert(Y0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    #  Convert the output a numpy array for later use
    output = np.array(output)


    rD = np.sqrt((output[:, 1]) ** 2 + (output[:, 2]) ** 2 + (output[:, 3]) ** 2)  # Position
    vD = np.sqrt((output[:, 4]) ** 2 + (output[:, 5]) ** 2 + (output[:, 6]) ** 2)  # Velocity

    return output

# uses Sun algorithm to compute position of Sun relative to planet
def RungeKuttaPert(r0, v0, T0, dT, tF, mu, muCelest):
    Y0 = np.concatenate([r0, v0], axis=0)
    t = np.arange(T0, tF, dT)
    dY = np.zeros([6, 1])

    # Orbit Propagator that ignores Perturbations
    def OrbitProp( t, Y, mu, muCelest):
        # Delta_UT1 = 0
        # Leap_Seconds = 32
        # JD_UT1 = JD_TAI - Leap_Seconds + Delta_UT1 + t/86400

        dY = np.zeros([6, 1])
        r = np.math.sqrt(Y[0] ** 2 + Y[1] ** 2 + Y[2] ** 2)
        SunPos = np.multiply(ValladoAlg.sun(t/86400 + T0), 149597870)
        Q = ValladoAlg.J2000toGCRF(1)
        SunPos = np.matmul(Q, np.transpose(SunPos))
        SunPos = np.squeeze(np.asarray(SunPos))

        #Sunpos = np.linalg.norm(Sunpos, 2)
        r_sat_sun = SunPos - Y[0:3]
        r_sat_sun = np.linalg.norm(r_sat_sun, 2)
        dY[0] = Y[3]
        dY[1] = Y[4]
        dY[2] = Y[5]

        a1 = muCelest*(((SunPos[0]- Y[0]) / (r_sat_sun ** 3)) - SunPos[0]/(np.linalg.norm(SunPos, 2))**3)
        a2 = muCelest*(((SunPos[1]- Y[1]) / (r_sat_sun ** 3)) - SunPos[1]/(np.linalg.norm(SunPos, 2))**3)
        a3 = muCelest*(((SunPos[2]- Y[2]) / (r_sat_sun** 3)) - SunPos[2]/(np.linalg.norm(SunPos, 2))**3)

        dY[3] = -(mu / r ** 3) * Y[0] + a1
        dY[4] = -(mu / r ** 3) * Y[1] + a2
        dY[5] = -(mu / r ** 3) * Y[2] + a3
        return dY

    def derivFcn(t, Y):
        return OrbitProp(t, Y, mu, muCelest)

    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y0, T0)

    # Define output array
    output = []
    output.append(np.insert(Y0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    #  Convert the output a numpy array for later use
    output = np.delete(output, 731, axis=0)
    output = np.array(output)


    rD = np.sqrt((output[:, 1]) ** 2 + (output[:, 2]) ** 2 + (output[:, 3]) ** 2)  # Position
    vD = np.sqrt((output[:, 4]) ** 2 + (output[:, 5]) ** 2 + (output[:, 6]) ** 2)  # Velocity

    return output

# uses RVPlanet algorithm (Algorithm #33) to compute position of Sun relative to planet
def RungeKuttaPertRV(r0, v0, T0, dT, tF, mu, muCelest, epoch):
    Y0 = np.concatenate([r0, v0], axis=0)
    t = np.arange(T0, tF, dT)
    dY = np.zeros([6, 1])

    # Orbit Propagator that ignores Perturbations
    def OrbitProp( t, Y, mu, muCelest):
        # Delta_UT1 = 0
        # Leap_Seconds = 32
        # JD_UT1 = JD_TAI - Leap_Seconds + Delta_UT1 + t/86400

        dY = np.zeros([6, 1])
        r = np.math.sqrt(Y[0] ** 2 + Y[1] ** 2 + Y[2] ** 2)
        #SunPos, SunVeloc = ValladoAlg.PlanetRV(muCelest, t);
        #SunPos, SunVeloc = ValladoAlg.PlanetRV(t/86400+T0, t)
        SunPos, SunVeloc = (ValladoAlg.PlanetRV(muCelest, epoch + t/86400))#, 149597870)

        #Sunpos = np.linalg.norm(Sunpos, 2)
        r_sat_sun = SunPos - Y[0:3]
        r_sat_sun = np.linalg.norm(r_sat_sun, 2)

        dY[0] = Y[3]
        dY[1] = Y[4]
        dY[2] = Y[5]

        dY[3] = -(mu / r ** 3) * Y[0] + muCelest*(((SunPos[0]- Y[0]) / (r_sat_sun ** 3)) - SunPos[0]/(np.linalg.norm(SunPos, 2))**3)
        dY[4] = -(mu / r ** 3) * Y[1] + muCelest*(((SunPos[1]- Y[1]) / (r_sat_sun ** 3)) - SunPos[1]/(np.linalg.norm(SunPos, 2))**3)
        dY[5] = -(mu / r ** 3) * Y[2] + muCelest*(((SunPos[2]- Y[2]) / (r_sat_sun** 3)) - SunPos[2]/(np.linalg.norm(SunPos, 2))**3)
        # print("HULLLO")
        return dY

    def derivFcn(t, Y):
        return OrbitProp(t, Y, mu, muCelest)

    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y0, T0)

    # Define output array
    output = []
    output.append(np.insert(Y0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    #  Convert the output a numpy array for later use
    output = np.array(output)
    output = np.delete(output, 731, axis=0)

    rD = np.sqrt((output[:, 1]) ** 2 + (output[:, 2]) ** 2 + (output[:, 3]) ** 2)  # Position
    vD = np.sqrt((output[:, 4]) ** 2 + (output[:, 5]) ** 2 + (output[:, 6]) ** 2)  # Velocity

    return output   # position as a function of time in j2000 frame, units AU

# Helper function for HW 7 to initialize propagator
# gets perturbed Acceleration with J2 and J3
def get_pertAcc_with_J2andJ3(r, mu, RE, J2, J3):
    Y = r # set Y = r to make copy paste for RK easier.
    r = np.linalg.norm(r, 2)

    dY = np.asarray(np.zeros((3,1)))
    #-(mu / r ** 3) * Y[0] +
    # -(mu / r ** 3) * Y[1] +
    # -(mu / r ** 3) * Y[2] +
    dY[0] =  ((-3 * J2 * mu * RE ** 2 * Y[0]) / (2 * r ** 5)) * (
                1 - ((5 * Y[2] ** 2) / (r ** 2))) \
            + ((-5 * J3 * mu * RE ** 3 * Y[0]) / (2 * r ** 7)) * (3 * Y[2] - ((7 * Y[2] ** 3) / (r ** 2)))

    dY[1] =  ((-3 * J2 * mu * RE ** 2 * Y[1]) / (2 * r ** 5)) * (
                1 - ((5 * Y[2] ** 2) / (r ** 2))) \
            + ((-5 * J3 * mu * RE ** 3 * Y[1]) / (2 * r ** 7)) * (3 * Y[2] - ((7 * Y[2] ** 3) / (r ** 2)))


    dY[2] =  ((-3 * J2 * mu * RE ** 2 * Y[2]) / (2 * r ** 5)) * (
                3 - ((5 * Y[2] ** 2) / (r ** 2))) \
            + ((-5 * J3 * mu * RE ** 3) / (2 * r ** 7)) * (
                        6 * Y[2] ** 2 - ((7 * Y[2] ** 4) / (r ** 2)) - ((3 * r ** 2) / (5)))

    dY = np.asarray(dY)
    return dY

def get_pertAcc_with_aSRP(r_Earth_sun, r_sat_Earth):
    r_sat_sun = r_Earth_sun - r_Earth_sun
    pSRP = 4.57*1e-6
    gamma = 1.4 # can I assume this???

    aSRP = -1*pSRP*gamma*((Cr*A)/m)*(r_sat_sun/np.linalg.norm(r_sat_sun))

    return aSRP
def test(d):

    if d > 1 and d == 2:
        d = 2
    if d < 1:
        d = 1
    return d

def get_atmDensity(r, RE):
    r = np.linalg.norm(r)
    h_ellp = r - RE

    if h_ellp > 0 and h_ellp < 25:
        h0 = 0; rho0 = 1.225; H = 7.249
    if h_ellp > 25 and h_ellp < 30:
        h0 = 25; rho0 = 3.899*1e-2; H = 6.349
    if h_ellp > 30 and h_ellp < 40:
        h0 = 30; rho0 = 1.774*1e-2; H = 6.682
    if h_ellp > 40 and h_ellp < 50:
        h0 = 40; rho0 = 3.972*1e-3; H = 7.554
    if h_ellp > 50 and h_ellp < 60:
        h0 = 50; rho0 = 1.057*1e-3; H = 8.382
    if h_ellp > 60 and h_ellp < 70:
        h0 = 60; rho0 = 3.206*1e-4; H = 7.714
    if h_ellp > 70 and h_ellp < 80:
        h0 = 70; rho0 = 8.770*1e-5; H = 6.549
    if h_ellp > 80 and h_ellp < 90:
        h0 = 80; rho0 = 1.905*1e-5; H = 5.799
    if h_ellp > 90 and h_ellp < 100:
        h0 = 90; rho0 = 3.396*1e-6; H = 5.382
    if h_ellp > 100 and h_ellp < 110:
        h0 = 100; rho0 = 5.297*1e-7; H = 5.877
    if h_ellp > 110 and h_ellp < 120:
        h0 = 110; rho0 = 9.661*1e-8; H = 7.263
    if h_ellp > 120 and h_ellp < 130:
        h0 = 120; rho0 = 2.438*1e-8; H = 9.473
    if h_ellp > 130 and h_ellp < 140:
        h0 = 130; rho0 = 8.484*1e-9; H = 12.636
    if h_ellp > 140 and h_ellp < 150:
        h0 = 140; rho0 = 3.845*1e-9; H = 16.149
    if h_ellp > 150 and h_ellp < 180:
        h0 = 150; rho0 = 2.070*1e-9; H = 22.523
    if h_ellp > 180 and h_ellp < 200:
        h0 = 180; rho0 = 5.464*1e-10; H = 29.74;
    if h_ellp > 200 and h_ellp < 250:
        h0 = 200; rho0 = 2.789*1e-10; H = 37.105
    if h_ellp > 250 and h_ellp < 300:
        h0 = 250; rho0 = 7.248*1e-10; H = 45.546
    if h_ellp > 300 and h_ellp < 350:
        h0 = 300; rho0 = 2.418 * 1e-11; H = 53.628
    if h_ellp > 350 and h_ellp < 400:
        h0 = 350; rho0 =  9.518 * 1e-12; H = 53.298
    if h_ellp > 400 and h_ellp < 450:
        h0 = 400; rho0 = 3.725 * 1e-12; H = 58.515
    if h_ellp > 400 and h_ellp < 450:
        h0 = 450; rho0 = 1.585*1e-12; H = 60.828
    if h_ellp > 500 and h_ellp < 600:
        h0 = 500; rho0 = 6.967*1e-13; H = 63.822
    if h_ellp > 600 and h_ellp < 700:
        h0 = 600; rho0 = 1.454*1e-13; H = 71.835
    if h_ellp > 700 and h_ellp < 800:
        h0 = 700; rho0 = 3.614*1e-14; H = 88.667
    if h_ellp > 800 and h_ellp < 900:
        h0 = 800; rho0 = 1.170*1e-14; H = 124.64
    if h_ellp > 900 and h_ellp < 1000:
        h0 = 900; rho0 = 5.245*1e-15; H = 181.05
    if h_ellp > 1000:
        h0 = 1000; rho0 = 3.019*1e-15; H = 268

    atmDensity = rho0*math.exp(-1*((h_ellp) - h0)/(H))

    return atmDensity