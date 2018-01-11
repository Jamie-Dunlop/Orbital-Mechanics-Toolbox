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
h = 0.001

n = (tf-t0)/h
n = int(n)
theta = np.zeros([n])
thetadot = np.zeros([n])
t = np.linspace(t0,tf,n)

thetadot[0] = thetadot0
theta[0] = theta0
for j in range(1,n):
    thetadot[j] = h*((-g/L)*theta[j-1]) + thetadot[j-1]
    theta[j] = h*(thetadot[j-1]) + theta[j-1]

# for j in range(n):
#     print (t[j],theta[j])

plt.plot(t,theta)
plt.grid()
plt.show()
