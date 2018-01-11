import math
import numpy as np
import matplotlib.pyplot as plt

#Inital conditions
t0 = 0
tf = 50
x0 = -0.8
m = 1
k = 2
c = 0.4
xdot0 = 0
#Acceleration
#acc = (-g/L)*theta
#Number of steps
h = 0.001

n = (tf-t0)/h
n = int(n)
xe = np.zeros([n])
xedot = np.zeros([n])
t = np.linspace(t0,tf,n)

xedot[0] = xdot0
xe[0] = x0
for j in range(1,n):
    xedot[j] = deltat*((-c*xedot[j-1]-k*xe[j-1])/m) + xedot[j-1]
    xe[j] = deltat*(xedot[j-1]) + xe[j-1]

for j in range(n):
    print (t[j],xe[j])

plt.plot(t,xe)
axes = plt.gca()
axes.set_xlim([0, 50])
axes.set_ylim([-0.8,0.8])
plt.xlabel("Time (seconds)")
plt.ylabel("X-position (m)")
plt.title("Euler Spring Mass Damper System")
plt.grid()
plt.show()
