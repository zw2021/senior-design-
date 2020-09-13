import math
import numpy as np
from numpy import matrix

# Prob 2
##################
#r1 = np.array([-2, 1, 1/2])# DU
#v1 = np.array([-.4, -.2, 1/2])# DU/TU
#h = np.cross(r1,v1)
#WW_hat = h/np.linalg.norm(h); #print(WW_hat)
#TT_hat = v1/np.linalg.norm(v1); #print(TT_hat)
#NN_hat = np.cross(TT_hat,WW_hat); #print(NN_hat)
#WW_hat_2 = np.cross(TT_hat,TT_hat); #print(WW_hat_2/np.linalg.norm(WW_hat_2))

# Prob 3
##################

r1 = np.array([-2, 1, 1/2]) # DU
v1 = np.array([-.4, -.2, 1/2])# DU/TU
R_hat = r1/np.linalg.norm(r1); #print(r1/np.linalg.norm(r1))
h1 = np.cross(r1,v1); W_hat = (h1/np.linalg.norm(h1))
S_hat = np.cross(W_hat,R_hat); print(S_hat)


r2 = np.array([-2, .9, .51]) # DU
v2 = np.array([-.39, -.21, .4]) # DU/TU

r21 = r2 - r1
print('r21:', r21)

Q_RSW_ijk = (R_hat, S_hat, W_hat)
Q_ijk_RSW = np.transpose(R_hat, S_hat, W_hat)

print('Q_ijk_RSW:\n', Q_ijk_RSW)

r21_RSW = np.matmul(Q_ijk_RSW, np.transpose(r21))
print('r21_RSW:', r21_RSW)

