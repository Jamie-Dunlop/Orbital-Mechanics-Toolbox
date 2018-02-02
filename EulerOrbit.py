import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse, Circle
import Constants

fig = plt.figure()
ax = fig.gca(projection = '3d')

#Inital conditions
t0 = 0
tf = 86400
#r0 = 42164000 #m
r0 = np.array([42164000, 0, 0]) #m  42164000 - geo distance
#rdot0 = 3070 #m/s
rdot0 = np.array([0, math.sqrt(Constants.mu / np.linalg.norm(r0)), 0]) #m/s
#Acceleration
T = 2 * math.pi * math.sqrt(np.linalg.norm(r0) ** 3 / Constants.mu)
#Number of steps
n = 86400000

deltat = (tf-t0)/n
# r = np.zeros([n])
# rdot = np.zeros([n])
# xdot = np.zeros([n])
# ydot = np.zeros([n])
# zdot = np.zeros([n])
# # theta = np.zeros([n])
# x = np.zeros([n])
# y = np.zeros([n])
# z = np.zeros([n])
t = np.linspace(t0,tf,n)

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

# rdot[0] = rdot0
# r[0] = r0
# x[0] = r0
# y[0] = 0

# f = 0

for j in range(1,n):
    rdot[j] = deltat * (-Constants.mu * r[j-1] // np.linalg.norm(r[j-1]) ** 3) + rdot[j-1]
    r[j] = deltat * (rdot[j-1]) + r[j-1]

    rmod[j-1] = np.linalg.norm(r[j-1])

    # xdot[j] = (deltat*(-Constants.mu/r0**2)*math.cos(theta[j]) + xdot[j-1])
    # x[j] = deltat*(xdot[j-1]) + x[j-1]
    # ydot[j] = (deltat*(-Constants.mu/r0**2)*math.sin(theta[j]) + ydot[j-1])
    # y[j] = deltat*(ydot[j-1]) + y[j-1]
    # theta[j] = (t[j] / T) * 2 * math.pi - f * 2 * math.pi
    # if theta[j] > 2*math.pi:
    #     f = f + 1
    # x[j] = r[j]*math.cos(theta[j])
    # y[j] = r[j]*math.sin(theta[j])

    rdot.insert(j, rdot[j])
    r.insert(j, r[j])
    rmod.insert(j, rmod[j-1])

    x.insert(j, r[j-1][0])
    y.insert(j, r[j-1][1])
    z.insert(j, r[j-1][2])

    if np.linalg.norm(r[j]) < Constants.Rearth:
        t = t[:j]
        break



    # if theta[j] < f*2*math.pi:
    #     theta[j] = theta[j]
    #     print (theta[j])
    #     x[j] = r[j]*math.cos(theta[j])
    #     y[j] = r[j]*math.sin(theta[j])
    #     print('f',f)
    # else:
    #     theta[j] = theta[j] - f*2*math.pi
    #     print (theta[j])
    #     x[j] = r[j]*math.cos(theta[j])
    #     y[j] = r[j]*math.sin(theta[j])
    #     f = f + 1
    #     print ('f', f)

# for j in range(n):
#     print (r[j])

# u = np.linspace(0, 2 * np.pi, 100)
# vearth = np.linspace(0, np.pi, 100)
# xe[u] = 6371000*math.cos(u)
# ye[u] = 6371000*math.sin(u)
# plt.plot(xe,ye)

# print(r)
# print(rdot)
# print(z)

# plt.plot(x,y)
# plt.xlabel("X-Position")
# plt.ylabel("Y-position")
# plt.title("Orbit")
# plt.grid()
# plt.show()

#Plot Orbit
ax.plot(x, y, z, color = 'r')
ax.set_xlabel("X-position (m)")
ax.set_ylabel("Y-position (m)")
ax.set_zlabel("Z-position (m)")
plt.title("Euler - Orbit")
plt.grid()

# Plot Earth sphere
u = np.linspace(0, 2 * np.pi, 100)
vearth = np.linspace(0, np.pi, 100)
xe = Constants.Rearth * np.outer(np.cos(u), np.sin(vearth))
ye = Constants.Rearth * np.outer(np.sin(u), np.sin(vearth))
ze = Constants.Rearth * np.outer(np.ones(np.size(u)), np.cos(vearth))
ax.plot_surface(xe, ye, ze, color='b')

ax.set_xlim(-8e7, 8e7)
ax.set_ylim(-8e7, 8e7)
ax.set_zlim(-8e7, 8e7)

plt.show()
