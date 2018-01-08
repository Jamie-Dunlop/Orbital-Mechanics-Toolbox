import math
import numpy as np
import matplotlib.pyplot as plt

#Inital conditions
t0 = 0
tf = 86400
r0 = 42164000 #m
mu = 3.986004418e14 #m something s something
rdot0 = 3070000 #m/s
#Acceleration
T = 2*math.pi*math.sqrt(r0**3/mu)
#Number of steps
n = 100000

deltat = (tf-t0)/n
r = np.zeros([n])
rdot = np.zeros([n])
theta = np.zeros([n])
x = np.zeros([n])
y = np.zeros([n])
t = np.linspace(t0,tf,n)

rdot[0] = rdot0
r[0] = r0
for j in range(1,n):
    rdot[j] = deltat*(-mu/r[j-1]**2) + rdot[j-1]
    r[j] = deltat*(rdot[j-1]) + r[j-1]
    theta[j] = (t[j]/T)*2*math.pi
    x[j] = r[j]*math.cos(theta[j])
    y[j] = r[j]*math.sin(theta[j])

for j in range(n):
    print (t[j],r[j])

plt.plot(x,y)
plt.xlabel("X-Position")
plt.ylabel("Y-position")
plt.title("Orbit")
plt.grid()
plt.show()
