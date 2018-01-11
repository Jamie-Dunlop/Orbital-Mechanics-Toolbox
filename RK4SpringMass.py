#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt
t0 = 0
tf = 50
x0 = -0.8
m = 1
k = 2
c = 0.4
xdot0 = 0
h = 0.001
n = (tf-t0)//h
n = int(n)
x = np.zeros([n])
xdot = np.zeros([n])
t = np.linspace(t0,tf,n)

xdot[0] = xdot0
x[0] = x0

def f(a,b):
    return (-c*a-k*b)/m

for j in range (1, n):

        k1xdot = f(xdot[j-1], x[j-1])
        k1x = xdot[j-1]

        k2xdot =f(xdot[j-1] + (k1xdot * (h/2)), x[j-1] + (k1xdot * (h/2)))
        k2x = xdot[j-1] + k1x * (h/2)

        k3xdot = f(xdot[j-1] + (k2xdot * (h/2)), x[j-1] + (k2xdot * (h/2)))
        k3x = xdot[j-1] + k2x * (h/2)

        k4xdot = f(xdot[j-1] + (k3xdot * (h)), x[j-1] + (k3xdot * (h)))
        k4x = xdot[j-1] + k3x * (h)

        xdot[j] = (xdot[j-1] + ((h/6) * (k1xdot + 2*k2xdot + 2*k3xdot + k4xdot)))
        x[j] = (x[j-1] + ((h/6) * (k1x + 2*k2x + 2*k3x + k4x)))

print (x)
plt.plot(t, x)
axes = plt.gca()
axes.set_xlim([0, 50])
axes.set_ylim([-0.8,0.8])
plt.xlabel("Time (seconds)")
plt.ylabel("X-position (m)")
plt.title("RK4 Spring Mass Damper System")
plt.grid()
plt.show()
