import numpy as np
import math

######
# Function that converts calendar date to Julian Date
def CalToJD(yr, mon, day, hr, min ,sec):

    JD = 367 * yr - np.floor(7 / 4 * (yr + np.floor((mon + 9) / 12))) + np.int(275 / 9 * mon) \
            + day + 1721013.5 + 1 / 24 * (hr + 1 / 60 * (min + sec / 60))

    return JD

# Function that converts Julian Date to UT1
def JDtoUTC(JD):
    TUTC = (JD - 2451545)/36525 # Vallado notation TUT1 means time in UT1
    return TUTC

# Function that converts Julian Date_UT1 to Time_UT1
def JDonetoUTCone(JD):
    TUTC = (JD - 2451545)/36525 # Vallado notation TUT1 means time in UT1
    return TUTC

# Function that converts TUT1 to sidereal greenwich mean sidereal time (theta_GMST)
def TUT1toGMST(TUT1):

    GMST = 67310.54841 + ((876600*3600) + 8640184.812866)*TUT1\
                 + 0.093104*(TUT1**2) - 6.2*( 10e-6 )*(TUT1**3)
    return GMST

# Function that makes JD to modified JD
def JDtoMJD(JD):
    MJD = JD - 2400000.5
    return MJD

# Algorithm 22: Julian Date to Gregorian Date
# def JDtoGregDate(JD):
#
#     T_1900 = (JD - 2415019.5)/(365.25)
#     Year = 1900 + np.floor(T_1900)
#     LeapYears = np.floor(1/4*(Year - 1900 - 1))
#     Days = (JD - 2415019.5) - (LeapYears + (Year - 1900)*365)
#     Days = (JD - 2415019.5) - (LeapYears + (Year - 1900)*365)
#
#     if Days < 1:
#         Year = Year - 1
#         LeapYears = np.floor((Year - 1900 -1)*1/4)
#
#     if (Year % 4) == 0:
#         Days = (JD - 2415019.5) - (LeapYears + 365*(Year - 1900))
#
#     LMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
#         LMonth[2] = 29
#
#     DOY = np.floor(Days)
#     while (np.sum(LMonth)+1) > DOY:
#         Mon = len(LMonth)#ask how to change this to number of months being summed
#         Day = DOY - LMonth*np.sum(LMonth)
#         T = (Days - DOY)*24
#         h = np.floor(T)
#         min = np.floor(60*(T - h))
#         s = (T - h - min/60)*3600
#         DOY += 1
#
#     return Year, Mon, Day, h, min, s

# def Lmonth(a):
#     if a == 1:
#         d = mon1 + mon2 + mon3 + mon3
#         lenMon = np.sum(d)
#     return lenMon