#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import Main
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse, Circle
from Six_orbital_elements import dot, cross, mod, orb_elems
import Constants
import Density

fig = plt.figure()
ax = fig.gca(projection='3d')

T = 2 * math.pi * math.sqrt (np.linalg.norm(Main.r0) ** 3 / Constants.mu)

n = (Main.tf-Main.t0)//Main.h
n = int(n)

Cd = Main.mass / ( Main.BClow * Main.Area)

x=[]
y=[]
z=[]
r = [[]]
r[0] = Main.r0
r.insert(1,Main.r0)
rmod = [[]]
rdot = [[]]
rdot[0] = Main.rdot0
rdot.insert(1, Main.rdot0)

alist=[]
elist=[]
ilist=[]
omegalist=[]
RAANlist=[]
Vlist=[]

t = np.linspace(Main.t0,Main.tf,n)

def f(X):
    Gravity = ((-Constants.mu) * (X)) / np.linalg.norm(X) ** 3 #monopole gravity model?
    Drag = - (0.5 * Density.Density1(np.linalg.norm(r[j-1])) * np.linalg.norm(rdot[j-1]) ** 2 * Main.Area * Cd) / Main.mass
    return  Gravity + Drag

for j in range (1, n):

        k1rdot = f(r[j-1])
        k1r = rdot[j-1]

        k2rdot =f(r[j-1] + (k1r))
        k2r = rdot[j-1] + k1rdot

        k3rdot = f(r[j-1] + (k2r))
        k3r = rdot[j-1] + k2rdot

        k4rdot = f(r[j-1] + (k3r))
        k4r = rdot[j-1] + k3rdot

        rdot[j] = rdot[j-1] + ((Main.h/6) * (k1rdot + 2 * k2rdot + 2 * k3rdot + k4rdot))
        r[j] = r[j-1] + ((Main.h/6) * (k1r + 2*k2r + 2*k3r + k4r))
        rmod[j-1] = np.linalg.norm(r[j-1])

        rdot.insert(j, rdot[j])
        r.insert(j, r[j])
        rmod.insert(j, rmod[j-1])

        x.insert(j, r[j-1][0])
        y.insert(j, r[j-1][1])
        z.insert(j, r[j-1][2])

        a = orb_elems(r[j], rdot[j], Constants.mu)
        alist.insert(j, a[0])
        elist.insert(j, a[1])
        ilist.insert(j, a[2])
        omegalist.insert(j, a[3])
        RAANlist.insert(j, a[4])
        Vlist.insert(j, a[5])

        if np.linalg.norm(r[j]) < Constants.Rearth:
            t = t[:j]
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
# ax.plot(x, y, z, color = 'r')
# ax.set_xlabel("X-position (m)")
# ax.set_ylabel("Y-position (m)")
# ax.set_zlabel("Z-position (m)")
# plt.title("RK4 - Orbit")
# plt.grid()

#Plot Earth sphere
# u = np.linspace(0, 2 * np.pi, 100)
# vearth = np.linspace(0, np.pi, 100)
# xe = 6371008 * np.outer(np.cos(u), np.sin(vearth))
# ye = 6371008 * np.outer(np.sin(u), np.sin(vearth))
# ze = 6371008 * np.outer(np.ones(np.size(u)), np.cos(vearth))
# ax.plot_surface(xe, ye, ze, color='b')
#
# ax.set_xlim(-8e6, 8e6)
# ax.set_ylim(-8e6, 8e6)
# ax.set_zlim(-8e6, 8e6)

# plt.show()

plt.subplot(3,3,1)
plt.plot(t,alist)
plt.xlabel("Time (s)")
plt.ylabel("Semi-Major Axis (m)")
plt.title("a")
plt.grid()

plt.subplot(3,3,2)
plt.plot(t,elist)
plt.xlabel("Time (s)")
plt.ylabel("Eccentricity")
plt.title("e")
plt.grid()

plt.subplot(3,3,3)
plt.plot(t,ilist)
plt.xlabel("Time (s)")
plt.ylabel("Inclination")
plt.title("i")
plt.grid()

plt.subplot(3,3,4)
plt.plot(t,omegalist)
plt.xlabel("Time (s)")
plt.ylabel("omega")
plt.title("omega")
plt.grid()

plt.subplot(3,3,5)
plt.plot(t,RAANlist)
plt.xlabel("Time (s)")
plt.ylabel("RAAN")
plt.title("RAAN")
plt.grid()

plt.subplot(3,3,6)
plt.plot(t,Vlist)
plt.xlabel("Time (s)")
plt.ylabel("V")
plt.title("RAAN")
plt.grid()

plt.show()
