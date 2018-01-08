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
n = 1000

deltat = (tf-t0)/n
x = np.zeros([n])
xdot = np.zeros([n])
t = np.linspace(t0,tf,n)

xdot[0] = xdot0
x[0] = x0
for j in range(1,n):
    xdot[j] = deltat*((-c*xdot[j-1]-k*x[j-1])/m) + xdot[j-1]
    x[j] = deltat*(xdot[j-1]) + x[j-1]

for j in range(n):
    print (t[j],x[j])

plt.plot(t,x)
plt.xlim(1)
plt.xlabel("Time (seconds)")
plt.ylabel("X-position (m)")
plt.title("Spring Mass Damper System")
plt.grid()
plt.show()
