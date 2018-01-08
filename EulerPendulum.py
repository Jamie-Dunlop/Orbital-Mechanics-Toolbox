import math
import numpy as np
import matplotlib.pyplot as plt

#Inital conditions
t0 = 0
tf = 100
theta0 = 0.5
g = 9.81
L = 9.81
thetadot0 = 0
#Acceleration
#acc = (-g/L)*theta
#Number of steps
n = 1000

deltat = (tf-t0)/n
theta = np.zeros([n])
thetadot = np.zeros([n])
t = np.linspace(t0,tf,n)

thetadot[0] = thetadot0
theta[0] = theta0
for j in range(1,n):
    thetadot[j] = deltat*((-g/L)*theta[j-1]) + thetadot[j-1]
    theta[j] = deltat*(thetadot[j-1]) + theta[j-1]

for j in range(n):
    print (t[j],theta[j])

plt.plot(t,theta)
plt.show()
