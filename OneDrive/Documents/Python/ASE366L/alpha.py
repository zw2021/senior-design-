import math
import numpy as np

# cto = 'cartesian to orbital'
# cto function converts cartesian elements to orbital elements
######################################################################
def cto(mu, r, v):
    # Compute Intermediary Values #
    z = (np.linalg.norm(v)**2)/2 - mu/(np.linalg.norm(r))  # Compute for Specific Energy
    # try:
    h = np.cross(r, v)  # Compute for Angular Momentum
    # except Exception as e:
    #     h = np.cross(np.transpose([r[0], r[1], r[2]]), np.transpose([v[0], v[1], v[2]]))
    #     print("Went to the backup")
    e = 1/mu*(np.cross(v,h)+(-mu*(r/np.linalg.norm(r))))
    ed = np.linalg.norm(e)
    n = np.cross([0,0,1], h)

    # Computations #
    if z < 0:
        a = -1*mu/(2*z)
    elif z == 0:
        print('a = infinity')
    elif z > 0:
        a = mu/(2*z)

    theta = np.arccos((np.dot(r,e))/(np.linalg.norm(r)*ed))
    i = np.arccos((np.dot(h,[0,0,1])/np.linalg.norm(h)))
    BigOmeg = np.arccos(n[0]/np.linalg.norm(n))
    LitOmeg = np.arccos(np.dot(n,e)/(np.linalg.norm(n)*np.linalg.norm(e)))

    # Angle Ambiguity Checks #
    if np.dot(r,v)<0:
        theta = 2*np.pi-theta

    if n[1] < 0:
        BigOmeg = 2*np.pi-BigOmeg

    if e[2] < 0:
        LitOmeg = 2*np.pi - LitOmeg

    # Singularity Checks #

    if (i == 0 and ed == 0):
        u = np.arccos(np.dot((n/np.linalg.norm(n)),r)/np.linalg.norm(r))
        if r[2] < 0:
            u = 2*np.pi-u
    elif ed > 0:
        u = LitOmeg + theta

    if ed == 0:
        OmegBar = np.arccos(np.dot(e, [1,0,0])/np.linalg.norm(e))
        if e[1] < 0:
            OmegBar = 2*np.pi - OmegBar

    if i == 0:
        l = np.arccos(np.dot(r, [1,0,0])/np.linalg.norm(r))
        if r[1] < 0:
            l = 2*np.pi-l;


    if (ed < math.pow(10,-12)) and ( i < math.pow(10,-12)):
        BigOmeg = 0
        LitOmeg = 0
        theta = l
    elif ed<10**-12:
        LitOmeg = 0
        theta = u
    elif i < 10^-12:
        BigOmeg = 0
        LitOmeg = OmegBar

    # Convert Rad to Deg #
    i = i*180/np.pi
    BigOmeg = BigOmeg*180/np.pi
    LitOmeg = LitOmeg*180/np.pi
    theta = theta*180/np.pi

    ans = [a,ed,i,BigOmeg,LitOmeg,theta] #a, RAN, theta
    #print('ans', ans)
    return ans

# Main Function
############################################################################
#mu =1;
#r = [-(math.sqrt(2)), math.sqrt(2), 0]; v = [math.sqrt(2)/-4, math.sqrt(2)/-4, 1/2]
#r = np.array([1.73008, 1.25982, -0.07740]); v = np.array([-0.22351, -0.65058, 0.07880])

#xx = cto(mu, r, v)
#print(xx)

#r = [-(2**1/2), 2**1/2, 0]  # DU
#v = [-(2**1/2)/4, -(2**1/2)/4, 1/2]  # DU/TU

