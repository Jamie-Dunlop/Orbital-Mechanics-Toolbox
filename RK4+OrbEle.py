#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse, Circle
from Six_orbital_elements import dot, cross, mod, orb_elems
import Constants

fig = plt.figure()
ax = fig.gca(projection='3d')

r0 = np.array([7371008, 0, 0]) #m  42164000 - geo distance
rdot0 = np.array([0, math.sqrt(Constants.mu/np.linalg.norm(r0)), 0]) #m/s
h = 0.5
Area = 5 #m^2
BChigh = 1
BClow = 12
mass = 50 #kg

T = 2 * math.pi * math.sqrt (np.linalg.norm(r0) ** 3 / Constants.mu)

NumOrbits = 5
t0 = 0
tf = NumOrbits * T #s

n = (tf-t0)//h
n = int(n)

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

alist=[]
elist=[]
ilist=[]
omegalist=[]
RAANlist=[]
Vlist=[]

t = np.linspace(t0,tf,n)

print(t)

#US atmospheric density model only applicable when alt > 25000m
def Rho(r):
    alt = r - Constants.Rearth
    Temp = -131.21 + 0.00299 * alt
    Pres = 2.488 * ((Temp + 273.1) / (216.6)) ** -11 #should that be 11.388?
    return Pres / (0.2869 * ( Temp + 273.1 ))


def f(X):
    Gravity = ((-Constants.mu) * (X)) / np.linalg.norm(X) ** 3 #monopole gravity model?
    Drag = - (0.5 * Rho(np.linalg.norm(r[j-1])) * np.linalg.norm(rdot[j-1]) ** 2 * Area * Cd) / mass
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

        rdot[j] = rdot[j-1] + ((h/6) * (k1rdot + 2 * k2rdot + 2 * k3rdot + k4rdot))
        r[j] = r[j-1] + ((h/6) * (k1r + 2*k2r + 2*k3r + k4r))
        rmod[j-1] = np.linalg.norm(r[j-1])

        rdot.insert(j, rdot[j])
        r.insert(j, r[j])
        rmod.insert(j, rmod[j-1])

        x.insert(j, r[j-1][0])
        y.insert(j, r[j-1][1])
        z.insert(j, r[j-1][2])

        a = orb_elems(r[j], rdot[j], Constants.mu)
        alist.insert(j-1, a[0])
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
