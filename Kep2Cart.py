#Anna Middleton 26/10/2016
#Converting Keplerian to cartesian coordinates
#Solving Kepler's equation for eccentric anomaly
#in order to forecast position and velocity vectors


#def kepler(ecc, M):
import math
import time
import numpy as np
import matplotlib.pyplot as plt
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
a = 26559000 #m
e = 0.704482
i = rad(63.1706)
omega = rad(281.646)
RAAN = rad(206.346)
mu = 3.985992e14

#Error tolerance
etol = 1e-8
#period
T = math.pi * 2 * math.sqrt(a**3/mu)
####print ("period", T, "secs")
t = 0
step = 1
resultsx = np.empty((0,1))
resultsy = np.empty((0,1))
resultsz = np.empty((0,1))
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
        ####print ("True anomaly" , deg(V))

        # rotation matrix
        u = omega + V

        R1 = math.cos(u) * math.cos(RAAN) - math.sin(u) * math.cos(i) * math.sin (RAAN)
        R2 = math.cos(u) * math.sin(RAAN) + math.sin(u) * math.cos(i) * math.cos(RAAN)
        R3 = math.sin(u) * math.sin(i)
        R = np.matrix([[R1], [R2], [R3]])

        ####print ("Rotation Matrix", R)

        #radius
        r=(a*(1 - e**2))/(1 + (e * math.cos(V)))
        r_comp = r * R
        ####print ("r", r_comp/1000, "km")

        #velocity
        #v = math.sqrt(mu * ((2/r) - (1/a)))
        #v_comp = v * R
        #print ("v", v_comp)

        v_x = math.sqrt(mu/(a * (1-e**2))) * math.sin(V)
        v_y = math.sqrt(mu/(a * (1-e**2))) * (e + math.cos(V))
        v_mag = math.sqrt(v_x**2 + v_y**2)
        ####print ("v", (v_mag * R)/1000)
        t = t + step
        ####print(t)
        #for k in range(1):
        resultsx = np.append(resultsx, r_comp[0])
        resultsy = np.append(resultsy, r_comp[1])
        resultsz = np.append(resultsz, r_comp[2])

except KeyboardInterrupt:
    print('interrupted!')

print(resultsx, resultsy, resultsz)
t1 = time.time()
total = t1-t0
print('Time',total)

#Plot Earth sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 6371000* np.outer(np.cos(u), np.sin(v))
y = 6371000* np.outer(np.sin(u), np.sin(v))
z = 6371000 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='b')

#Plot Orbit
ax.plot(resultsx, resultsy, resultsz, color='r')
ax.legend()
#ax.axis('equal')
max_range = np.array([resultsx.max()-resultsx.min(), resultsy.max()-resultsy.min(), resultsz.max()-resultsz.min()]).max() / 2.0

mid_x = (resultsx.max()+resultsx.min()) * 0.5
mid_y = (resultsy.max()+resultsy.min()) * 0.5
mid_z = (resultsz.max()+resultsz.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)
#ax.autoscale(enable=True, axis='both', tight=None)
plt.show()
