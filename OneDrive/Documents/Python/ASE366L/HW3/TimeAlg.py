import numpy as np
import math


# Algorithim for years 1900 - 2100
# ROUGH DRAFT CODE! SEE PG 202, Vallado 4th edition

def JD_Conv(JD,DOY):
    T_1900 = (JD-2415019.5)/365.25
    Year = 1900 + np.floor(T_1900)
    LeapYears = np.floor((Year - 1900 -1)*(.25))
    Days = (JD - 2415019.5) - ( (Year - 1900)*365 + LeapYears)

    LMonth = np.array(len(Days)) # array that contains the number of days in a year

    if Days < 1:
        Year = Year - 1
        LeapYears = np.floor((Year - 1900 - 1)*.25)
        Days = (JD - 2415019.5) - ((Year - 1900) * 365 + LeapYears)
    if (Year % 4) == 0:
        LMonth[2] = 29      # Accounts for Leap Years - AKA changing February

    DOY = np.floor(DOY)

    while LMonth +1 > DOY:
        Mon = len(LMonth)
        Day = DOY - LMonth
        T = (Days - DOY)*24
        h = np.floor(T)
        min = np.floor((T-h)*60)
        s = (T - h - min/60)*3600

    # Convert DOY to Month and Day (See Vallado 3.6.4)
    D = Days - DOY
    # Convert D (day fraction) to hour, minutes, seconds

    # Output: Year, Month, Day, Hour, Minute Seconds
    t = 3.555559
    xx = np.floor(t)
    return xx


# Algorithm for converting DOY to Month and Day
# Vallado 3.6.4
JD = np.array([0])
DOY =2
xx = JD_Conv(JD, DOY)
print(xx)