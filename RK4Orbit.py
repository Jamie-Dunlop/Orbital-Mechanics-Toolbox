#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import Main
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse, Circle
import Constants
import Density

fig = plt.figure()
ax = fig.gca(projection='3d')
n = (Main.tf-Main.t0)//Main.h
n = int(n)

#Period
T = 2 * math.pi * math.sqrt (np.linalg.norm(Main.r0) ** 3 / Constants.mu)

x=[]
y=[]
z=[]
r = [[]]
r[0] = Main.r0
r.insert(1,Main.r0)
rmod = [[]]
rdot = [[]]
rdot[0] = Main.rdot0
rmod[0] = np.linalg.norm(Main.r0)
rdot.insert(1, Main.rdot0)

t = np.linspace(Main.t0,Main.tf,n)

#Force Model
def f(X):
    Gravity = ((-Constants.mu) * (X)) / np.linalg.norm(X) ** 3 #monopole gravity model?
    Drag = - (0.5 * Main.DensityModel(np.linalg.norm(r[j-1])) * np.linalg.norm(rdot[j-1]) ** 2 * Main.Area * Main.Cd) / Main.mass
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
        rmod[j-1] = np.linalg.norm(r[j])

        r.insert(j, r[j])
        rmod.insert(j, rmod[j-1])
        rdot.insert(j, rdot[j])
        x.insert(j, r[j-1][0])
        y.insert(j, r[j-1][1])
        z.insert(j, r[j-1][2])

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
plt.figure(1)
ax.plot(x, y, z, color = 'r')
ax.set_xlabel("X-position (m)")
ax.set_ylabel("Y-position (m)")
ax.set_zlabel("Z-position (m)")
plt.title("RK4 - Orbit")
plt.grid()

#Plot Earth sphere
u = np.linspace(0, 2 * np.pi, 100)
vearth = np.linspace(0, np.pi, 100)
xe = Constants.Rearth * np.outer(np.cos(u), np.sin(vearth))
ye = Constants.Rearth * np.outer(np.sin(u), np.sin(vearth))
ze = Constants.Rearth * np.outer(np.ones(np.size(u)), np.cos(vearth))
ax.plot_surface(xe, ye, ze, color='b')

ax.set_xlim(-8e6, 8e6)
ax.set_ylim(-8e6, 8e6)
ax.set_zlim(-8e6, 8e6)

plt.figure(2)
plt.plot(t,rmod)
plt.xlabel("Time (s)")
plt.ylabel("Radius (m)")
plt.title("Radius change over time")
plt.grid()



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
