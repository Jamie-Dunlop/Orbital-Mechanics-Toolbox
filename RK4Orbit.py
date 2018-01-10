#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt
t0 = 0
tf = 8640
r0 = 42164000 #m
mu = 3.986004418e14 #m something s something
rdot0 = 3070000 #m/s
h = 1
n = (tf-t0)//h
n = int(n)
T = 2*math.pi*math.sqrt(r0**3/mu)

r = np.zeros([n])
rdot = np.zeros([n])
theta = np.zeros([n])
x = np.zeros([n])
y = np.zeros([n])
t = np.linspace(t0,tf,n)

rdot[0] = rdot0
r[0] = r0
x[0] = r0
y[0] = 0

def f(Rad):
    return -mu/ Rad**2

for j in range (1, n):

        k1rdot = f(r[j-1])
        k1r = rdot[j-1]

        k2rdot =f(rdot[j-1] + (k1rdot * (h/2)))
        k2r = r[j-1] + k1r * (h/2)

        k3rdot = f(rdot[j-1] + (k2rdot * (h/2)))
        k3r = r[j-1] + k2r * (h/2)

        k4rdot = f(rdot[j-1] + (k3rdot * (h)))
        k4r = r[j-1] + k3r * (h)

        rdot[j] = (rdot[j-1] + ((h/6) * (k1rdot + 2*k2rdot + 2*k3rdot + k4rdot)))
        r[j] = (r[j-1] + ((h/6) * (k1r + 2*k2r + 2*k3r + k4r)))

        theta[j] = (t[j]/T)*2*math.pi

        print (r[j])

        x[j] = r[j]*math.cos(theta[j])
        y[j] = r[j]*math.sin(theta[j])

print (r)
print (x, y)
plt.plot(x, y)
plt.xlabel("X-position (m)")
plt.ylabel("Y-position (m)")
plt.title("RK4 - Orbit")
plt.grid()
plt.show()
