import math
import numpy as np
import matplotlib.pyplot as plt

#Inital conditions
t0 = 0
tf = 864000
r0 = 42164000 #m
mu = 3.986004418e14 #m something s something
rdot0 = 3070 #m/s
#Acceleration
T = 2*math.pi*math.sqrt(r0**3/mu)
#Number of steps
n = 10000

Area = 5 #m^2
BChigh = 1
BClow = 12
mass = 50 #kg
Rearth = 6371008 #m

h = r0 - Rearth

T = -131.21 + 0.00299 * h

p = 2.488 * ((T + 273.1) / (216.6)) ** -11.388

Rho = p / (0.2869 * ( T + 273.1 ))

V = np.sqrt ( mu / r0 )

Cd = mass / ( BChigh * Area)

Drag = 0.5 * Rho * V ** 2 * Area * Cd

deltat = (tf-t0)/n
r = np.zeros([n])
rdot = np.zeros([n])
xdot = np.zeros([n])
ydot = np.zeros([n])
theta = np.zeros([n])
x = np.zeros([n])
y = np.zeros([n])
t = np.linspace(t0,tf,n)

rdot[0] = rdot0
r[0] = r0
x[0] = r0
y[0] = 0
f = 0
for j in range(1,n):
    rdot[j] = deltat*(Drag / mass) + rdot[j-1]
    r[j] = deltat*(rdot[j-1]) + r[j-1]

    theta[j] = (t[j]/T)*2*math.pi - f*2*math.pi
    if theta[j] > 2*math.pi:
        f = f + 1
    x[j] = r[j]*math.cos(theta[j])
    y[j] = r[j]*math.sin(theta[j])



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
print(x)
print(y)
plt.plot(x,y)
plt.xlabel("X-Position")
plt.ylabel("Y-position")
plt.title("Orbit")
plt.grid()
plt.show()
