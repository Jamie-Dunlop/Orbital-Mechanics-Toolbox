#Anna Middleton 26/10/2016
#Converting Keplerian to cartesian coordinates
#Solving Kepler's equation for eccentric anomaly
#in order to forecast position and velocity vectors


#def kepler(ecc, M):
import math
import numpy as np
import matplotlib.pyplot as plt
import Constants
import Main
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse, Circle

# fig = plt.figure()
# ax = fig.gca(projection='3d')

def rad(X):
    return ((X * math.pi) / 180)

def deg(X):
    return ((X / math.pi) * 180)

mu = Constants.mu
e = Main.e
i = rad(Main.i)
omega = rad(Main.omega)
RAAN = rad(Main.RAAN)

mean_motion = ((Main.Mean_motion) * 2 * math.pi) / 86400 #rad/s

a = (Constants.mu / mean_motion ** 2) ** (1/3)

#Error tolerance
etol = 1e-8

#Period
T = math.pi * 2 * math.sqrt(a ** 3 / mu)
print ("period", T, "secs")

t = 0
step = 0.1
r_resultsx = np.empty((0,1))
r_resultsy = np.empty((0,1))
r_resultsz = np.empty((0,1))
vsat_resultsx = np.empty((0,1))
vsat_resultsy = np.empty((0,1))
vsat_resultsz = np.empty((0,1))
try:
    while t <= T:

        #Mean anomaly
        M = math.sqrt(mu/a**3) * t

        #Initial guess for Eccentric anomaly
        if M < math.pi:
            E = M - e/2
        else:
            E = M + e/2

        #Iterates Kepler's equation until error is less than tolerance
        error = 1
        try:
            while abs(error) > etol:
                error = (M - E + e * math.sin(E))/ (1 - e*math.cos(E))
                E = E + error

        except KeyboardInterrupt:
            print('interrupted!')

        #True anomaly
        V = 2 * math.atan(math.sqrt((1+e) / (1-e) ) * math.tan(E/2))

        #Rotation matrix

        R11 = math.cos(omega) * math.cos(RAAN) - math.sin(omega) * math.cos(i) * math.sin (RAAN)
        R21 = math.cos(omega) * math.sin(RAAN) + math.sin(omega) * math.cos(i) * math.cos(RAAN)
        R31 = math.sin(omega) * math.sin(i)

        R12=-math.sin(omega)*math.cos(RAAN)-math.cos(omega)*math.cos(i)*math.sin(RAAN);
        R22=-math.sin(omega)*math.sin(RAAN)+math.cos(omega)*math.cos(i)*math.cos(RAAN);
        R32=math.cos(omega)*math.sin(i);

        R13=math.sin(i)*math.sin(RAAN);
        R23=-math.sin(i)*math.cos(RAAN);
        R33=math.cos(i);

        if e==1:
            p=2*a
        else:
            p = a*(1-e**2)

        #Radius
        r=(p)/(1 + (e * math.cos(V)))
        xp=r*math.cos(V)
        yp=r*math.sin(V)
        wom_dot = math.sqrt(mu*p)/r**2

        #velocity
        r_dot = math.sqrt(mu / p) * e * math.sin(V)
        vxp=r_dot*math.cos(V)-r*math.sin(V)*wom_dot
        vyp=r_dot*math.sin(V)+r*math.cos(V)*wom_dot

        t = t + step

        r_resultsx = R11*xp+R12*yp
        r_resultsy = R21*xp+R22*yp
        r_resultsz = R31*xp+R32*yp

        vsat_resultsx = R11*vxp+R12*vyp
        vsat_resultsy = R21*vxp+R22*vyp
        vsat_resultsz = R31*vxp+R32*vyp

except KeyboardInterrupt:
    print('interrupted!')

print('r', r_resultsx, r_resultsy, r_resultsz)
print('v', vsat_resultsx, vsat_resultsy, vsat_resultsz)

# #Plot Earth sphere
# u = np.linspace(0, 2 * np.pi, 100)
# vearth = np.linspace(0, np.pi, 100)
# x = 6371000* np.outer(np.cos(u), np.sin(vearth))
# y = 6371000* np.outer(np.sin(u), np.sin(vearth))
# z = 6371000 * np.outer(np.ones(np.size(u)), np.cos(vearth))
# ax.plot_surface(x, y, z, color='g')
#
# #Plot Orbit
# ax.plot(r_resultsx, r_resultsy, r_resultsz, color='g')
# ax.legend()
# #ax.axis('equal')
# max_range = np.array([r_resultsx.max()-r_resultsx.min(), r_resultsy.max()-r_resultsy.min(), r_resultsz.max()-r_resultsz.min()]).max() / 2.0
#
# mid_x = (r_resultsx.max()+r_resultsx.min()) * 0.5
# mid_y = (r_resultsy.max()+r_resultsy.min()) * 0.5
# mid_z = (r_resultsz.max()+r_resultsz.min()) * 0.5
# ax.set_xlim(mid_x - max_range, mid_x + max_range)
# ax.set_ylim(mid_y - max_range, mid_y + max_range)
# ax.set_zlim(-1e7, mid_z + max_range)
# plt.show()
