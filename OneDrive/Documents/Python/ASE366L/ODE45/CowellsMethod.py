import numpy as np
import math
from ASE366L.ODE45 import RungeKutta as RK


def Pert(mu, rVec_Sat, rVec_Earth):

    rD_Earth = np.linalg.norm(rVec_Earth, 2)
    rD_Sat = np.linalg.norm(rVec_Sat, 2)

    aTB = np.array([mu*((rVec_Sat/rD_Sat**3) - (rVec_Earth/(rD_Earth**3)))])

    # return acceleration due to third body perturbations
    return aTB