#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt

t0 = 0
tf = 86400 #s
mu = 3.986004418e14 #m something s something
r0 = np.array([42164000, 0, 0]) #m
rdot0 = np.array([0, 3070, 0]) #m/s
h = 1
n = (tf-t0)//h
n = int(n)

T = 2 * math.pi * math.sqrt (np.linalg.norm(r0) ** 3 / mu)

x=[]
y=[]
z=[]
r = [[]]
r[0] = r0
r.insert(1,r0)

rdot = [[]]
rdot[0] = rdot0
rdot.insert(1, rdot0)

t = np.linspace(t0,tf,n)

def f(X):
    return (mu) * (X) / np.linalg.norm(X) ** 3

for j in range (1, n):

        k1rdot = f(r[j-1])
        k1r = rdot[j-1]

        # print('k1rdot', k1rdot)
        # print('k1r', k1r)

        k2rdot =f(r[j-1] + (k1r * (h/2)))
        k2r = rdot[j-1] * k1rdot * (h/2)

        # print('k2rdot', k2rdot)
        # print('k2r', k2r)

        k3rdot = f(r[j-1] + (k2r * (h/2)))
        k3r = rdot[j-1] * k2rdot * (h/2)

        # print('k3rdot', k3rdot)
        # print('k3r', k3r)

        k4rdot = f(r[j-1] + (k3r * (h)))
        k4r = rdot[j-1] * k3rdot * (h)

        # print('k4rdot', k4rdot)
        # print('k4r', k4r)

        rdot[j] = rdot[j-1] + ((h/6) * (k1rdot + 2 * k2rdot + 2 * k3rdot + k4rdot))
        r[j] = r[j-1] + ((h/6) * (k1r + 2*k2r + 2*k3r + k4r))

        r.insert(j, r[j])
        rdot.insert(j, rdot[j])
        x.insert(j, r[j-1][0])
        y.insert(j, r[j-1][1])
        z.insert(j, r[j-1][2])



print (r)
print (x, y)
plt.plot(x, y)
plt.xlabel("X-position (m)")
plt.ylabel("Y-position (m)")
plt.title("RK4 - Orbit")
plt.grid()
plt.show()
