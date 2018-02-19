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

No_of_orbits = int(Main.tf/T)

print('Period', T)
x= np.zeros([n])
y= np.zeros([n])
z= np.zeros([n])
xdot= np.zeros([n])
ydot= np.zeros([n])
zdot= np.zeros([n])
x[0] = Main.r0[0]
y[0] = Main.r0[1]
z[0] = Main.r0[2]
xdot[0] = Main.rdot0[0]
ydot[0] = Main.rdot0[1]
zdot[0] = Main.rdot0[2]
r = [[]]
rmod = np.zeros([n])
rdot = [[]]
diff = np.zeros([No_of_orbits])
Orbit_no = np.arange(1, No_of_orbits+1)
# rdot[0] = Main.rdot0
rmod[0] = np.linalg.norm(Main.r0)
# rdot.insert(1, Main.rdot0)

state = np.array([Main.r0,Main.rdot0])
t = np.linspace(Main.t0,Main.tf,n)

#Force Model including gravity and drag
def Accel(R,V):
    Gravity = ((-Constants.mu) * (R)) / np.linalg.norm(R) ** 3 #monopole gravity model?
    Drag = - (0.5 * Main.DensityModel(np.linalg.norm(R)) * np.linalg.norm(V) * Main.AreaH * Main.Cd*V) / Main.mass
    # print(Drag)
    return Gravity + Drag

#Gets position and velocity from state vector and calculates acceleration from Accel
def Orbit(t, state):
    pos, vel = state
    return np.array([vel, Accel(pos,vel)])

#Runge-Kutta 4 Integrator
def rk4(x, h, y, f):
    k1 = h * f(x, y)
    k2 = h * f(x + 0.5*h, y + 0.5*k1)
    k3 = h * f(x + 0.5*h, y + 0.5*k2)
    k4 = h * f(x + h, y + k3)
    return x + h, y + (k1 + 2*(k2 + k3) + k4)/6.0

#Starts RK4 and stores position and veloctiy values obtained from RK4
timed = 0
for j in range(1,n):
    timed, state = rk4(timed, Main.h, state, Orbit)
    x[j] = state[0][0]
    y[j] = state[0][1]
    z[j] = state[0][2]
    xdot[j] = state[1][0]
    ydot[j] = state[1][1]
    zdot[j] = state[1][2]
    r = [x,y,z]
    rdot = [xdot,ydot,zdot]
    vmod = math.sqrt(xdot[j]**2+ydot[j]**2+zdot[j]**2)
    rmod[j] = math.sqrt(x[j]**2+y[j]**2+z[j]**2)
    # print('Rmod',rmod[j])
    # print('Rmod',rmod[jj])
    if rmod[j] < Constants.Rearth:
        t = t[:j]
        # rmod = rmod[:j]
        break

## Prints Rmod at the start and end of an orbit for error check
for jj in range(0,int(Main.tf/T)):
    # print('R at start of orbit',jj+1,'>>',rmod[int((jj)*T/Main.h)],'m')
    # print('R at end of orbit',jj+1,'>>',rmod[int((jj+1)*T/Main.h)],'m')
    diff[jj] = rmod[int((jj+1)*T/Main.h)]-rmod[int((jj)*T/Main.h)]
    # print('Difference from R at start >>', diff)

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
print('name',Main.name)
# Plot Orbit
plt.figure(1)
ax.plot(x, y, z, color = 'r')
ax.set_xlabel("X-position (m)")
ax.set_ylabel("Y-position (m)")
ax.set_zlabel("Z-position (m)")
plt.title("RK4 - Orbit of {}".format(Main.name))
plt.grid()

# Plot Earth sphere
u = np.linspace(0, 2 * np.pi, 100)
vearth = np.linspace(0, np.pi, 100)
xe = Constants.Rearth * np.outer(np.cos(u), np.sin(vearth))
ye = Constants.Rearth * np.outer(np.sin(u), np.sin(vearth))
ze = Constants.Rearth * np.outer(np.ones(np.size(u)), np.cos(vearth))
ax.plot_surface(xe, ye, ze, color='b')

ax.set_xlim(-7e6, 7e6)
ax.set_ylim(-7e6, 7e6)
ax.set_zlim(-7e6, 7e6)

plt.figure(2)
plt.plot(t,rmod)
plt.xlabel("Time (s)")
plt.ylabel("Radius (m)")
plt.title("Radius change over time {}".format(Main.name))
plt.grid()

plt.figure(3)
plt.plot(Orbit_no, diff)
plt.xlabel("Orbit Number")
plt.ylabel("Difference in Radius (m)")
plt.title("Radius difference between orbits {}".format(Main.name))
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
