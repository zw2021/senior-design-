import numpy as np
import math

######
# Function that converts calendar date to Julian Date
def CalToJD(yr, mon, day, hr, min ,sec):
    JD = 367*yr - np.floor(7/4*(yr + np.floor((mon+9)/12))) + np.int(275/9*mon)\
        + day + 1721013.5 + 1/24*(hr + 1/60*(min + sec/60))
    return JD

# Function that converts Julian Century to UT1 Century
def JDtoUTC(JD):
    TUTC = (JD - 2451545)/36525 # Vallado notation TUT1 means time in UT1
    return TUTC

# Function that converts Julian Century to UTT Century
def JDtoUTT(JD):
    TUTT = (JD - 2451545)/36525 # Vallado notation TUT1 means time in UT1
    return TUTT

# Function that converts TUT1 to sidereal greenwich mean sidereal time (theta_GMST)
def TUT1toGMST(TUT1):

    GMST = 67310.54841 + ((876600*3600) + 8640184.812866)*TUT1\
                 + 0.093104*(TUT1**2) - 6.2*( 10e-6 )*(TUT1**3)
    return GMST

# Function that makes JD to modified JD
def JDtoMJD(JD):
    MJD = JD - 2400000.5
    return MJD

# Converts JD time to theta Era time
def ERA(t1_JD_UT1):
    theta_ERA = 2*np.pi*(0.7790572732640+1.00273781191135448*(t1_JD_UT1-2451545))
    theta_ERA = np.mod(theta_ERA, 2*np.pi)
    return theta_ERA

#####

# Problem 2 Check
yr = 2020
mon = 2; day = 14; hr = 21 # 3:00 PM + 6 or 15:00 PM + 6?;
min = 0; sec = 0

# yr = 1992
# mon = 8; day = 20; hr = 12
# min = 14; sec = 0

JD = CalToJD(yr, mon, day, hr, min, sec)    # Calculate JD in UTC
MJD = JDtoMJD(JD)   # Calculate resultant JD into MJD for higher precision
TUTC = JDtoUTC(JD) # Find TUTC from MJD in UTC
TUT1 = (TUTC) - 0.200316    # Take resultant TUTC and add corresponding deltaUTC to transform into TUT1
GMST = TUT1toGMST(TUTC)

t= np.mod(GMST, 86400)  # Get GMST in terms of 86400
tt = t/240              # Convert GMST from seconds to degrees
print('Problem 2:', tt)


######## Problem 3
# Gregorian Date to JD
yr = 2020
mon = 7; day = 4; hr = 14
min = 47; sec = 42

JD_2 = CalToJD(yr, mon, day, hr, min, sec)    # Calculate JD in UTC
TUTC = JDtoUTC(JD)
print('Problem')
###########################
