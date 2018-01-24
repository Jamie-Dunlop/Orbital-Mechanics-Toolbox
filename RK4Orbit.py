#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse, Circle
fig = plt.figure()
ax = fig.gca(projection='3d')
t0 = 0
tf = 864000 #s
mu = 3.986004418e14 #m something s something
r0 = np.array([6871008, 0, 0]) #m  42164000 - geo distance
rdot0 = np.array([0, math.sqrt(mu/np.linalg.norm(r0)), 0]) #m/s
h = 0.5
Area = 5 #m^2
BChigh = 1
BClow = 12
mass = 50 #kg
Rearth = 6371008 #m
n = (tf-t0)//h
n = int(n)

T = 2 * math.pi * math.sqrt (np.linalg.norm(r0) ** 3 / mu)

Cd = mass / ( BClow * Area)

x=[]
y=[]
z=[]
r = [[]]
r[0] = r0
r.insert(1,r0)
rmod = [[]]
rdot = [[]]
rdot[0] = rdot0
rdot.insert(1, rdot0)

t = np.linspace(t0,tf,n)

def Rho(r):
    alt = r - Rearth
    T = -131.21 + 0.00299 * alt
    p = 2.488 * ((T + 273.1) / (216.6)) ** -11
    return p / (0.2869 * ( T + 273.1 ))

def f(X):
    return ((-mu) * (X)) / np.linalg.norm(X) ** 3 - (0.5 * Rho(np.linalg.norm(r[j-1])) * np.linalg.norm(rdot[j-1]) ** 2 * Area * Cd) / mass

for j in range (1, n):

        k1rdot = f(r[j-1])
        k1r = rdot[j-1]

        # k1rdot = (h/2) * f(r[j-1])
        # k1r = (h/2) * rdot[j-1]

        # print('k1rdot', k1rdot)
        # print('k1r', k1r)

        k2rdot =f(r[j-1] + (k1r))
        k2r = rdot[j-1] + k1rdot

        # k2rdot =(h/2) * f(r[j-1] + (k1r))
        # k2r = (h/2) * (rdot[j-1] + k1rdot)

        # print('k2rdot', k2rdot)
        # print('k2r', k2r)

        k3rdot = f(r[j-1] + (k2r))
        k3r = rdot[j-1] + k2rdot

        # k3rdot = h * f(r[j-1] + (k2r))
        # k3r = h * (rdot[j-1] + k2rdot)

        # print('k3rdot', k3rdot)
        # print('k3r', k3r)

        k4rdot = f(r[j-1] + (k3r))
        k4r = rdot[j-1] + k3rdot

        # k4rdot = (h/2) * f(r[j-1] + (k3r))
        # k4r = (h/2) * (rdot[j-1] + k3rdot)

        # print('k4rdot', k4rdot)
        # print('k4r', k4r)

        rdot[j] = rdot[j-1] + ((h/6) * (k1rdot + 2 * k2rdot + 2 * k3rdot + k4rdot))
        r[j] = r[j-1] + ((h/6) * (k1r + 2*k2r + 2*k3r + k4r))
        rmod[j-1] = np.linalg.norm(r[j])
        # rdot[j] = rdot[j-1] + ((1/3) * (k1rdot + (2 * k2rdot) + k3rdot + k4rdot))
        # r[j] = r[j-1] + ((1/3) * (k1r + (2 * k2r) + k3r + k4r))

        r.insert(j, r[j])
        rmod.insert(j, rmod[j-1])
        rdot.insert(j, rdot[j])
        x.insert(j, r[j-1][0])
        y.insert(j, r[j-1][1])
        z.insert(j, r[j-1][2])

        if np.linalg.norm(r[j]) < Rearth:
            break

# # print ('r', r)
# # print('rdot', rdot)
# # print (x, y)
# plt.plot(x,y)
# plt.subplot(2,1,1)
# plt.plot(x,y)
# plt.xlabel("X-position (m)")
# plt.ylabel("Y-position (m)")
# plt.title("RK4 - Orbit")
# plt.grid()
#
# plt.subplot(2,1,2)
# plt.plot(t,rmod)
# plt.xlabel("Time (s)")
# plt.ylabel("Radius (m)")
# plt.title("Radisu change over time")
# plt.grid()
#
# plt.show()

# print ('r', r)
# print('rdot', rdot)
# print (x, y)

#Plot Orbit
ax.plot(x, y, z, color = 'r')
ax.set_xlabel("X-position (m)")
ax.set_ylabel("Y-position (m)")
ax.set_zlabel("Z-position (m)")
plt.title("RK4 - Orbit")
plt.grid()

#Plot Earth sphere
u = np.linspace(0, 2 * np.pi, 100)
vearth = np.linspace(0, np.pi, 100)
xe = 6371008 * np.outer(np.cos(u), np.sin(vearth))
ye = 6371008 * np.outer(np.sin(u), np.sin(vearth))
ze = 6371008 * np.outer(np.ones(np.size(u)), np.cos(vearth))
ax.plot_surface(xe, ye, ze, color='b')

ax.set_xlim(-8e6, 8e6)
ax.set_ylim(-8e6, 8e6)
ax.set_zlim(-8e6, 8e6)


# xmod = np.linalg.norm(x)
# ymod = np.linalg.norm(y)
# zmod = np.linalg.norm(z)

# max_range = np.array([xmod.max()-xmod.min(), ymod.max()-ymod.min(), zmod.max()-zmod.min()]).max() / 2.0
#
# mid_x = (xmod.max()+xmod.min()) * 0.5
# mid_y = (ymod.max()+ymod.min()) * 0.5
# mid_z = (zmod.max()+zmod.min()) * 0.5
# ax.set_xlim(mid_x - max_range, mid_x + max_range)
# ax.set_ylim(mid_y - max_range, mid_y + max_range)
# ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.show()
