#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt

t0 = 0
tf = 864000 #s
mu = 3.986004418e14 #m something s something
r0 = np.array([32164000, 0, 0]) #m  42164000 - geo distance
rdot0 = np.array([0, 3070, 0]) #m/s
h = 1
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

        # rdot[j] = rdot[j-1] + ((1/3) * (k1rdot + (2 * k2rdot) + k3rdot + k4rdot))
        # r[j] = r[j-1] + ((1/3) * (k1r + (2 * k2r) + k3r + k4r))

        r.insert(j, r[j])
        rdot.insert(j, rdot[j])
        x.insert(j, r[j-1][0])
        y.insert(j, r[j-1][1])
        z.insert(j, r[j-1][2])

        if np.linalg.norm(r[j]) < Rearth:
            break



# print ('r', r)
# print('rdot', rdot)
# print (x, y)
plt.plot(x, y)
plt.xlabel("X-position (m)")
plt.ylabel("Y-position (m)")
plt.title("RK4 - Orbit")
plt.grid()
plt.show()
