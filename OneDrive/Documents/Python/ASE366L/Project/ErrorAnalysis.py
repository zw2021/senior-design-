import numpy as np
import math # 你好
import matplotlib.pyplot as plt

# 4 range measurement residual equation #
# only takes in xi, yi, zi array of size 1 by 4 or 4 by 1
# for higher order of residuals need to ask Dr. Jones
def CalcDOP(xi, x, yi, y, zi, z):

    R = np.zeros((len(xi), 1))
    A = np.asmatrix(np.zeros((len(xi), 4)))

    for ii  in range(len(xi)):
        R[ii] = math.sqrt((xi[ii] - x) ** 2 + (yi[ii] - y) ** 2 + (zi[ii] - z) ** 2)
        A[ii, :] = [(xi[ii] - x) / R[ii], (yi[ii] - y) / R[ii], (zi[ii] - z) / R[ii], -1]

    Q = np.transpose(np.transpose(A)*A)
    PDOP = math.sqrt((Q[0, 0] ** 2) + (Q[1, 1] ** 2) + (Q[2, 2]) ** 2)
    TDOP = math.sqrt((Q[3, 3]**2))
    GDOP = math.sqrt( PDOP ** 2 + TDOP ** 2)

    return  PDOP, TDOP, GDOP

# Test #
xi = np.array([1, 2, 3, 4])
yi = np.array([1, 2, 3, 4])
zi = np.array([1, 2, 3, 4])
x = 1; y = 2; z = 3
A, Q, B = CalcDOP(xi, x, yi, y, zi, z)

print(A)
print(Q)
print(B)