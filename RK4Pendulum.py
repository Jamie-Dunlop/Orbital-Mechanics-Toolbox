#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt
t0 = 0
tf = 100
theta0 = 0.5
g = 9.81
L = 9.81
thetadot0 = 0
n = 1000
h = (tf-t0)/n
theta = np.zeros([n])
thetadot = np.zeros([n])
thetadot[0] = thetadot0
theta[0] = theta0
t = np.linspace(t0,tf,n)

for j in range (1, n):

        k1thetadot = h * ((-g/L)*theta[j-1])
        k1theta = h * thetadot[j-1]

        k2thetadot =(thetadot[j-1] + (k1thetadot * (h/2)))
        k2theta = theta[j-1] + k1theta * (h/2)

        k3thetadot = (thetadot[j-1] + (k2thetadot * (h/2)))
        k3theta = theta[j-1] + k2theta * (h/2)

        k4thetadot = (thetadot[j-1] + (k2thetadot * (h)))
        k4theta = theta[j-1] + k2theta * (h)

        thetadot[j] = (thetadot[j] + ((1/6) * (k1thetadot + 2*k2thetadot + 2*k3thetadot + k4thetadot)))
        theta[j] = (theta[j] + ((1/6) * (k1theta + 2*k2theta + 2*k3theta + k4theta)))

print (theta)
plt.plot(t, theta)
plt.show()
