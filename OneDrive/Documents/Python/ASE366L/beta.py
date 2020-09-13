import numpy as np
import math # 你好
from ASE366L.Functions import R3Generator, R1Generator

def beta(mu, x):
#x = np.array([a, ed, i, LitOmeg, BigOmeg, theta[ii]])  # matrix of orbital elements

    a = x[0]
    e = x[1]
    i = x[2]
    LitOmeg = x[3]
    BigOmeg = x[4]
    theta = x[5]

    i = i*np.pi/180
    BigOmeg = BigOmeg*np.pi/180
    LitOmeg = LitOmeg*np.pi/180
    #theta = theta*np.pi/180 # UNCOMMENT THETA ALREADY FED AS RAD IN PREVIOUS HOMEWORK

    # Compute position and velocity vectors in pqw frame #
    p = a*(1-math.pow(e,2))
    r_pqw = (p/(1+e*np.cos(theta)))*np.array([np.cos(theta), np.sin(theta), 0])
    v_pqw = (math.sqrt(mu/p))*np.array([-np.sin(theta), (e+np.cos(theta)), 0])

    # Convert position and velocity vectors in ijk frame #
    # Find Q_ijk_pqw: converting from inertial ref to pqw ref #
    R3_LitOmeg = R3Generator(LitOmeg)
    R3_BigOmeg = R3Generator(BigOmeg)
    R1_i = R1Generator(i)

    Q_ijk_pqw = np.matmul(R3_LitOmeg, np.matmul(R1_i,R3_BigOmeg))
    Q_pqw_ijk = np.transpose(Q_ijk_pqw)

    # Matrix Q from ijk to pqw

    rIJK = np.matmul(Q_pqw_ijk, np.transpose(r_pqw))
    vIJK = np.matmul(Q_pqw_ijk, np.transpose(v_pqw))

    #value = np.concatenate((rIJK, vIJK))
    return rIJK, vIJK

# Main Function
#####################################################
# mu = 1  # DU3/TU2
# a = 1.4; e = 0.55069; i = 117.14;
# LitOmeg = 177.95; BigOmeg = 268.84; theta = 188.18
#
# x = np.array([a, e, i, LitOmeg, BigOmeg, theta])
# bb = beta(mu,x)
# print(bb)
#print('r_IJK:', bb[0])
#print('v_IJK:', bb[1])