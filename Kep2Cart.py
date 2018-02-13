#Anna Middleton 26/10/2016
#Converting Keplerian to cartesian coordinates
#Solving Kepler's equation for eccentric anomaly
#in order to forecast position and velocity vectors


#def kepler(ecc, M):
import math
import time
import numpy as np
import matplotlib.pyplot as plt
import Constants
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse, Circle
t0 = time.time()
fig = plt.figure()
ax = fig.gca(projection='3d')

def rad(X):
    return ((X * math.pi) / 180)

def deg(X):
    return ((X / math.pi) * 180)

#Test values
# a = 26559000 #m
# e = 0.704482
# i = rad(63.1706)
# omega = rad(281.646)
# RAAN = rad(206.346)
mu = Constants.mu
#Geo-Orbit
# a = 42169440 #m
e = 0.0005304
i = rad(277.9340)
omega = rad(86.2686)
RAAN = rad(355.7087)

mean_motion = (1.00272536) #revs/day
mean_motion1 = mean_motion / (24*60*60) #revs/secs
mean_motion2 = mean_motion1 * 2 * math.pi

a = (Constants.mu / mean_motion2 ** 2) ** (1/3)

print('a', a)



#Error tolerance
etol = 1e-8
#period
T = math.pi * 2 * math.sqrt(a ** 3 / mu)
print ("period", T, "secs")
t = 0
step = 1
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

        ####print ("Error" ,error)
        ####print ("Eccentric anomaly", deg(E)) #2.71

        #true anomaly
        V = 2 * math.atan(math.sqrt((1+e) / (1-e) ) * math.tan(E/2))
        # print ("True anomaly" , deg(V))

        # rotation matrix
        u = omega + V

        R11 = math.cos(omega) * math.cos(RAAN) - math.sin(omega) * math.cos(i) * math.sin (RAAN)
        R21 = math.cos(omega) * math.sin(RAAN) + math.sin(omega) * math.cos(i) * math.cos(RAAN)
        R31 = math.sin(omega) * math.sin(i)
        

        R12=-math.sin(omega)*math.cos(RAAN)-math.cos(omega)*math.cos(i)*math.sin(RAAN);
        R22=-math.sin(omega)*math.sin(RAAN)+math.cos(omega)*math.cos(i)*math.cos(RAAN);
        R32=math.cos(omega)*math.sin(i);

        R13=math.sin(i)*math.sin(RAAN);
        R23=-math.sin(i)*math.cos(RAAN);
        R33=math.cos(i);


        ####print ("Rotation Matrix", R)
        if e==1:
            p=2*a
        else:
            p = a*(1-e**2)


        #radius
        r=(p)/(1 + (e * math.cos(V)))
        xp=r*math.cos(V)
        yp=r*math.sin(V)
        ####print ("r", r_comp/1000, "km")
        wom_dot = math.sqrt(mu*p)/r**2
        #velocity
        r_dot = math.sqrt(mu / p) * e * math.sin(V)
        vxp=r_dot*math.cos(V)-r*math.sin(V)*wom_dot
        vyp=r_dot*math.sin(V)+r*math.cos(V)*wom_dot
        #print ("v", v)

        #v_x = math.sqrt(mu/(a * (1-e**2))) * math.sin(V)
        #v_y = math.sqrt(mu/(a * (1-e**2))) * (e + math.cos(V))
        #v_mag = math.sqrt(v_x**2 + v_y**2)
        ####print ("v", (v_mag * R)/1000)
        t = t + step
        r_resultsx = R11*xp+R12*yp
        r_resultsy = R21*xp+R22*yp
        r_resultsz = R31*xp+R32*yp

        vsat_resultsx = R11*vxp+R12*vyp;
        vsat_resultsy = R21*vxp+R22*vyp;
        vsat_resultsz = R31*vxp+R32*vyp;

        ####print(t)
        #for k in range(1):
        # r_resultsx = np.append(r_resultsx, r_comp[0])
        # r_resultsy = np.append(r_resultsy, r_comp[1])
        # r_resultsz = np.append(r_resultsz, r_comp[2])
        #
        # vsat_resultsx = np.append(vsat_resultsx, vsat_comp[0])
        # vsat_resultsy = np.append(vsat_resultsy, vsat_comp[1])
        # vsat_resultsz = np.append(vsat_resultsz, vsat_comp[2])

except KeyboardInterrupt:
    print('interrupted!')

print ("True anomaly" , deg(V))

print('r', r_resultsx, r_resultsy, r_resultsz)
print('v', vsat_resultsx, vsat_resultsy, vsat_resultsz)
t1 = time.time()
total = t1-t0
print('Time',total)

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
