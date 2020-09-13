import numpy as np
import math
from ASE366L import beta
from ASE366L.TimeConversions import JDonetoUTCone
from ASE366L.Functions import DegtoRad, RadtoDeg, R3Generator, R2Generator, R1Generator, AUtoKM

# Newton Raphson Method from SpaceD, outputs nu (or theta) as rad
# t is a time array; n, ed, tp as int values
def ellip_prop(t, n, ed, tp):
    theta = np.zeros(len(t))
    M = np.zeros(len(t))
    E = np.zeros(len(t))
    x_r_old = np.zeros(len(t))
    x_r_new = np.zeros(len(t))

    t = tp + t

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

        theta[i] = 2*math.atan2(math.sqrt((1+ed))*math.tan(E[i]/2), math.sqrt(1-ed)) # return theta as rad
    return theta

# Algorithm 2 pg 65: Kep Eqtn delta M, e to E
# Takes M as Rad
def KepEqntoE(M, e):

    if (-1*np.pi < M) or (M > np.pi):
        E = M - e
    else:
        E = M + e

    # E_n = E
    tol = 100*1e-12 # set tolerances

    # set while loop #
    E_new= E + (M - E + e * np.sin(E)) / (1 - e * np.cos(E))
    ii = 0

    while not abs(E_new - E) < tol:
        E_new = E + (M - E + e*np.sin(E))/(1-e*np.cos(E))
        E = E_new
        ii += 1
        # print("END OUR SUFFERING", ii)
    return  E_new    # returns E_new in Radians

# Algorithm 6: E to nu pg 77
def AnomtoNu(e, E):

    if e < 1:
        nu = 2*math.atan2(math.sqrt((1+e))*np.tan(E/2), math.sqrt(1-e))
        #nu = math.asin((np.sin(E)*math.sqrt(1-e**2))/(1-e*np.cos(E)))
        #EF = 2 * math.atan2(np.math.sqrt(1 - ed) * np.tan(thetaF / 2), np.math.sqrt(1 + ed))
    return nu

# Algorithm 10: COE2RV pg. 118
# outputs r_IJK and v_IJK in km
def COE2RV(mu, p, e, i, BigOmeg, LitOmeg, nu):#, u, LongTrue, OmegSquig):

   # r_IJK, v_IJK = COE2RV(mu, p, e, i, BigOmeg, LitOmeg, nu)#, u, LongTrue, OmegSquig)

    # if (e == 0) and (i == 0):
    #     LitOmeg = 0; BigOmeg = 0;
    #     nu = LongTrue
    # if (e == 0):
    #     LitOmeg = 0; nu = LongTrue
    # if (i == 0):
    #     BigOmeg = 0; LitOmeg = OmegSquig

    r_pqw = np.array([(p*np.cos(nu)/(1+e*np.cos(nu))), ((p*np.sin(nu))/(1+e*np.cos(nu))), 0])
    v_pqw = np.array([-1*(math.sqrt(mu/p)*np.sin(nu)), (math.sqrt(mu/p)*(e+np.cos(nu))), 0])
    Q = np.matmul(R3Generator(-1*BigOmeg), np.matmul(R1Generator(-1*i), R3Generator(-1*LitOmeg)))

    r_IJK = np.matmul(Q, np.transpose(r_pqw))
    v_IJK = np.matmul(Q, np.transpose(v_pqw))

    r_IJK = np.squeeze(np.asarray(r_IJK))
    v_IJK = np.squeeze(np.asarray(v_IJK))

    return r_IJK, v_IJK # returns r_IJK and v_IJK in km

# Algorithm 29: Sun(JD => position vec of Sun) See Vallado pg 279
def sun(JD_UT1):
    # returns position of Sun in MOD frame in km
    T_UT1 =  (JD_UT1 - 2451545)/36525 # Vallado notation TUT1 means time in UT1
    MeanLongSun = DegtoRad(280.46 + 36000.771*T_UT1)  # Calculate Mean Longitude of Sun in rad

    T_TDB = T_UT1
    MeanSun = DegtoRad(357.52772333 + 35999.05034*(T_TDB)) # Calculate Mean Anomaly of Sun in rad

    EccLong = MeanLongSun + (DegtoRad((1.914666471)*np.sin(MeanSun) + 0.019994643*np.sin(2*MeanSun))) # Calculate ecliptic Longitude
    r_sun = 1.000140612 - ( 0.016708617*np.cos(MeanSun) ) - (0.000139589*np.cos(2*MeanSun))
    Eps = DegtoRad(23.439291 - 0.0130042*(T_TDB) )   # Short for 'Epsilon': Obliquity of the Ecliptic


    # gets position of the sun in the Time of Day (TOD) frame, but we will assume that MOD (Mean of Date)
    # MOD was created with reference to the J2000 frame
    r_vec_Sun = np.array([r_sun*np.cos(EccLong), r_sun*np.cos(Eps)*np.sin(EccLong), r_sun*np.sin(Eps)*np.sin(EccLong)])

    return r_vec_Sun

def MODt0GCRF(t_UTC, r_MOD):
    # Function requires t_TT to be in arc seconds
    T_TT = (t_UTC - 2451545) / 36525
    zeta = (1/3600)*(2306.2181*T_TT + (0.30188*(T_TT**2)) + (0.017998*(T_TT**3)))
    theta = (1/3600)*(2004.3109*T_TT - 0.42665*(T_TT**2) - 0.041833*(T_TT**3))
    z = (1/3600)*(2306.2181*T_TT + 1.09468*(T_TT**2) + 0.018203*(T_TT**3))

    zeta = DegtoRad(zeta)
    theta = DegtoRad(theta)
    z = DegtoRad(z)

    R3_zeta = R3Generator(zeta)
    R2_theta = R2Generator(-1*theta)
    R3_z = R3Generator(z)

    Q_MOD_GCRF = np.matmul(R3_zeta, np.matmul(R2_theta, R3_z))
    r_GCRF = np.matmul(Q_MOD_GCRF, np.transpose(r_MOD))
    r_GCRF = np.squeeze(np.asarray(r_GCRF))
    return r_GCRF

def J2000toGCRF(r_MOD):
    alpha = 1/3600*(0.0146)
    alpha = DegtoRad(alpha)
    zeta = 1/3600*(-0.16617)
    zeta = DegtoRad(zeta)
    nu = 1/3600*(-0.0068192)
    nu = 1/3600*(nu)

    Q_J2000_GCRF = np.matmul(R3Generator(-alpha), np.matmul(R2Generator(-zeta), R1Generator(nu)))

    return Q_J2000_GCRF

# Algorith 33 from Vallado, pg 296 4th ed
# Outputs location of Earth relative to the Sun in J2000 frame, units: km
# Assumes JTDB is the input into algorithm
def PlanetRV(Planet, JD_TDB):
    T_TT = (JD_TDB - 2451545)/ 36525   # Included in Vallado but simplified for class
    T_TDB = T_TT

   # T_TDB = JD_TDB
    mu = Planet

    # coefficients on PG 1046 from D4: Planetary Ephemerides
    a = 1.000001018 # AU
    e = 0.01670862 - 0.000042037*T_TDB - 0.0000001236*T_TDB**2 + 0.00000000004*T_TDB**3 # AU
    i = 0.0000000 + 0.0130546*T_TDB - 0.00000931*T_TDB**2 - 0.000000034*T_TDB**3 # Deg
    # Check errata #
    BigOmeg = 174.873174 - 0.2410908*T_TDB + 0.00015026*T_TDB**2 + 0.000000478*T_TDB**3 # Deg
    OmegSquig = 102.937348 + 0.322557*T_TDB + 0.00015026*T_TDB**2 + 0.000000478*T_TDB**3    # Deg
    MeanLong = 100.466449 + 35999.3728519*T_TDB - 0.00000568*T_TDB**2 + 0.00000000000*T_TDB**3 # Deg

    i = DegtoRad(i); BigOmeg = DegtoRad(BigOmeg); # convert all angle values into radians
    OmegSquig = DegtoRad(OmegSquig); MeanLong = DegtoRad(MeanLong);


    a = AUtoKM(a)   # convert value to km
    M = MeanLong - OmegSquig    # units: rad
    LitOmeg = OmegSquig - BigOmeg   # units: rad

    # Algorithm 2, pg 65 4th ed Vallado #
    E = KepEqntoE(M,e)  # units: rad
    nu = AnomtoNu(e, E)
    #a = AUtoKM(a)
    x = np.array([a, e, i, LitOmeg, BigOmeg, nu])  # matrix of orbital elements
    p = a*(1 - e**2)  # p is in km
    # check for singularities: can we assume that we have none?#
    # x = np.array([a, ed, i, LitOmeg, BigOmeg, theta[ii]])  # matrix of orbital elements

    # r, v = beta.beta(mu, x)
    # r = np.squeeze(np.asarray(r)); v = np.squeeze(np.asarray(v))

    # if (i == 0) and (e == 0):
    #     u = np.arccos(np.dot((n / np.linalg.norm(n)), r) / np.linalg.norm(r))
    #     if r[2] < 0:
    #         u = 2 * np.pi - u
    # elif e > 0:
    #     u = LitOmeg + nu
    #
    # if e == 0:
    #     OmegSquig = np.arccos(np.dot(e, [1, 0, 0]) / np.linalg.norm(e))
    #     if e[1] < 0:
    #         OmegSquig = 2 * np.pi - OmegSquig
    #
    # elif e != 0:
    #     OmegSquig = 0
    #
    # if i == 0:
    #     LongTrue = np.arccos(np.dot(r, [1, 0, 0]) / np.linalg.norm(r))
    #     if r[1] < 0:
    #         LongTrue = (2 * np.pi) - LongTrue;
    #
    # elif i != 0:
    #     LongTrue = 0


    r_IJK, v_IJK = COE2RV(mu, p, e, i, BigOmeg, LitOmeg, nu)#, u, LongTrue, OmegSquig)

    Eps = 23.439291 # deg
    Eps = DegtoRad(Eps)
    Q = R1Generator(-1* Eps)

    r_XYZ = np.matmul(Q, np.transpose(r_IJK))
    v_XYZ = np.matmul(Q, np.transpose(v_IJK))

    r_XYZ = np.squeeze(np.asarray(r_XYZ))   # units: km
    v_XYZ = np.squeeze(np.asarray(v_XYZ))   # units: km

    return r_XYZ, v_XYZ # returns r and v in J200 frame, units: km


def getMoonPos_wrt_Earth(JD_TDB):
    '''
    Get the position of the Moon with respect to the Earth based on the
    2010 Astronomical Almanac algorithms.  Returns position vector in
    J2000, geocentric, equatorial frame.
    '''

    #  Convenience functions to make copy/paste from MATLAB easier
    def sind(arg):
        return np.sin(arg *np.pi/180.0)

    def cosd(arg):
        return np.cos(arg * np.pi/180.0)

    # Convert the input time into Julian Centuries
    T_TDB = (JD_TDB - 2451545.0) / 36525.0
    D2R_ = np.pi / 180.0

    long_ecliptic = 218.32 + 481267.8813 * T_TDB \
                    + 6.29 * sind(134.9 + 477198.85 * T_TDB) \
                    - 1.27 * sind(259.2 - 413335.38 * T_TDB) \
                    + 0.66 * sind(235.7 + 890534.23 * T_TDB) \
                    + 0.21 * sind(269.9 + 954397.70 * T_TDB) \
                    - 0.19 * sind(357.5 + 35999.05 * T_TDB) \
                    - 0.11 * sind(186.6 + 966404.05 * T_TDB)  # deg

    lat_ecliptic = 5.13 * sind(93.3 + 483202.03 * T_TDB) \
                   + 0.28 * sind(228.2 + 960400.87 * T_TDB) \
                   - 0.28 * sind(318.3 + 6003.18 * T_TDB) \
                   - 0.17 * sind(217.6 - 407332.20 * T_TDB)  # deg

    parralax = 0.9508 + 0.0518 * cosd(134.9 + 477198.85 * T_TDB) \
               + 0.0095 * cosd(259.2 - 413335.38 * T_TDB) \
               + 0.0078 * cosd(235.7 + 890534.23 * T_TDB) \
               + 0.0028 * cosd(269.9 + 954397.70 * T_TDB)  # deg

    #  Use the unwrap() function to get an angle in the range [-pi,pi]
    long_ecliptic = np.unwrap([long_ecliptic * D2R_])[0]
    lat_ecliptic = np.unwrap([lat_ecliptic * D2R_])[0]
    parralax = np.unwrap([parralax * D2R_])[0]

    obliquity = 23.439291 - 0.0130042 * T_TDB  # deg
    obliquity = obliquity * D2R_

    # Pre-compute cosine and sine values that will be needed at least twice
    cos_lag_ecl = np.cos(lat_ecliptic)
    sin_lon_ecl = np.sin(long_ecliptic)
    sin_lat_ecl = np.sin(lat_ecliptic)
    cos_obliq = np.cos(obliquity)
    sin_obliq = np.sin(obliquity)

    # ------------- calculate moon position vector ----------------
    pos_j2000 = np.zeros((3,))
    pos_j2000[0] = cos_lag_ecl * np.cos(long_ecliptic)
    pos_j2000[1] = cos_obliq * cos_lag_ecl * sin_lon_ecl - sin_obliq * sin_lat_ecl
    pos_j2000[2] = sin_obliq * cos_lag_ecl * sin_lon_ecl + cos_obliq * sin_lat_ecl

    #  Generate position in km where Re = 6378.1363 km
    pos_j2000 *= 6378.1363 / np.sin(parralax)

    #  Return position
    return pos_j2000