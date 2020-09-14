import numpy as np
import matplotlib.pyplot as plt


T1 = 300    # kelvin
T3 = 1600

p1 = 50     # kPa
p2 = p1*50
p3 = p2
p4 = p1
R = 287; g = 1.4
cp = (R*g)/(g - 1)
T4 = 513.139

s1 = 0
T2 = 917.363
s2 = cp*np.log(T2/T1)-R*np.log(p2/p1)
s3 = s2 + cp*np.log(T3/T2)-R*np.log(p3/p2)
s4 = s3
x = [s1,s2,s3,s4,s1]
y = [T1,T2,T3,T4,T1]
plt.plot(x, y, color='green', linestyle='dashed', linewidth = 1,
         marker='o', markerfacecolor='blue', markersize=6)
plt.annotate("T1,s1", (x[0], y[0]), weight='bold', fontsize=12, color='red')
plt.annotate("T2,s2", (x[1], y[1]), weight='bold', fontsize=12, color='red')
plt.annotate("T3,s3", (x[2], y[2]), weight='bold', fontsize=12, color='red')
plt.annotate("T4,s4", (x[3], y[3]), weight='bold', fontsize=12, color='red')

plt.xlabel('Entropy, s (kJ/kg/K)')
plt.ylabel('Temperature, T (K)')
plt.grid('on')
plt.title('Brayton Cycle T-s Graphic')
plt.show()

