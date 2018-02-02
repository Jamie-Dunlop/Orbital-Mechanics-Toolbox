#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt
import Constants

t0 = 0
tf = 1000
theta0 = 1
L = 9.798
thetadot0 = 0
h = 0.001
n = (tf-t0)//h
n = int(n)
theta = np.zeros([n])
thetadot = np.zeros([n])
thetadot[0] = thetadot0
theta[0] = theta0
t = np.linspace(t0,tf,n)

def f(x):
    return (-Constants.g / L) * x

for j in range (1, n):

        k1thetadot = f(theta[j-1])
        k1theta = thetadot[j-1]

        k2thetadot =f(theta[j-1] + (k1thetadot * (h/2)))
        k2theta = thetadot[j-1] + k1theta * (h/2)

        k3thetadot = f(theta[j-1] + (k2thetadot * (h/2)))
        k3theta = thetadot[j-1] + k2theta * (h/2)

        k4thetadot = f(theta[j-1] + (k3thetadot * (h)))
        k4theta = thetadot[j-1] + k3theta * (h)

        thetadot[j] = (thetadot[j-1] + ((h/6) * (k1thetadot + 2*k2thetadot + 2*k3thetadot + k4thetadot)))
        theta[j] = (theta[j-1] + ((h/6) * (k1theta + 2*k2theta + 2*k3theta + k4theta)))

print (theta)
plt.plot(t, theta)
plt.show()
