import numpy as np
import math

def R3Generator(theta):
    a = theta
    R3 = np.array([[(np.cos(a)), (np.sin(a)), (0)], \
                 [(-np.sin(a)), (np.cos(a)), (0)], \
                 [0, 0, 1]])
    R3 = np.matrix(R3)
    return R3

def R2Generator(theta):
    a = theta
    R2 = np.array([[(np.cos(a)), (0), (-np.sin(a))], \
                 [(0), (1), (0)], \
                 [(np.sin(a)), (0), (np.cos(a))]])
    R2 = np.matrix(R2)
    return R2

def R1Generator(theta):
    a = theta
    R1 = np.array([[(1), (0), (0)], \
                 [(0), (np.cos(a)), (np.sin(a))], \
                 [(0), (-np.sin(a)), (np.cos(a))]])
    R1 = np.matrix(R1)
    return R1

def DMStoRad(Deg,ArcMin, ArcSec):
    alpha = (Deg + ArcMin/60 + ArcSec/3600)*np.pi/(180)
    return alpha

def RadtoDMS(alpha):
    Temp = alpha*(180/np.pi)
    Deg = np.trunc(Temp)
    ArcMin = np.trunc((Temp - Deg)*60)
    ArcSec = 3600*(Temp - Deg - 1/60*(ArcMin))
    return Deg, ArcMin, ArcSec

# Functtion that calculates s_prime given Century in TT
def s_prime_Generator(TUTT):
    s_prime = -0.000047*(TUTT)
    return s_prime

def ArcmintoRad(a):
    a = a*(1 / 15) * (1 / 60) * (1 / 4) * (np.pi / 180)
    return a

def WRPN_MatrixGenerator(s, s_prime, xp, yp, X, Y, theta_ERA):

    s = ArcmintoRad(s)
    s_prime = ArcmintoRad(s_prime)
    xp = ArcmintoRad(xp)
    yp = ArcmintoRad(yp)
    s_prime = -1*s_prime

    R3_s = R3Generator(s_prime)
    R2_xp = R2Generator(xp)
    R1_yp = R1Generator(yp)
    W = np.matmul(R3_s, np.matmul(R2_xp, R1_yp))

    theta_ERA = -1*theta_ERA
    R = R3Generator(theta_ERA)

    a = 0.5 + .125*(X**2 + Y**2)
    a = ArcmintoRad(a)
    X = ArcmintoRad(X)
    Y = ArcmintoRad(Y)

    PN = np.array([[(1 - a*(X**2)), (-a*X*Y), (X)],\
                  [(-a*X*Y), (1- (a*Y**2)), (Y)],\
                  [(-X), (-Y), (1-(a*(X**2+Y**2)))]])
    R3_s = R3Generator(s)
    PN = np.matmul(PN, R3_s)

    return W,R,PN

# Function to convert AU to km
def AUtoKM(a):
    b = a*149597870.7
    return b
# Function to convert km to AU
def KMtoAU(b):
    a = b/149597870.7
    return a
# Function to convert Deg to Radians
def DegtoRad(a):
    b = a*(np.pi/180)
    return b

# Function to convert Deg to Radians
def RadtoDeg(a):
    b = a*(180/np.pi)
    return b

class Earth:
    mu = 398600.4418 # units km3/s2

class Moon:
    r = 384.4 * 1e3  # units km; Distance from sun to Earth's center of mass
    mu = 4.902*1e3 # units km3/s2

class Sun:
    r = 149.6*1e6 # units km; Distance from sun to Earth's center of mass
    mu = 1.327*1e11 # units km3/s2

class deltatimes:
    ut1 = -.200316 # from hw3
    AT = 37
