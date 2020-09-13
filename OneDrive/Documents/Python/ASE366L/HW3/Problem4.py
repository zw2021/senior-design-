import numpy as np
import math

# Helper Function
###############
def TTtoUTC(TT,DeltaAT):
    UTC = TT - DeltaAT - 32.184
    return UTC

# Main Function
################

TT = 2458165.4375 # units: seconds
DeltaAt = 37 # units: seconds

UTC = TTtoUTC(TT, DeltaAt)
