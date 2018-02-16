import math
import numpy as np
import matplotlib.pyplot as plt
import Constants

#Inital conditions
t0 = 0
tf = 10
theta0 = 0.1
L = 9.81
thetadot0 = 0
h = 0.01
n = (tf-t0)//h
n = int(n)
t = np.linspace(t0,tf,n)
theta = np.zeros([n])
thetar = np.zeros([n])
thetadotr = np.zeros([n])
thetae = np.zeros([n])
thetadote = np.zeros([n])
print('t',t)
#Exact Solution
for j in range(0,n):
    theta[j] = theta0*math.cos(math.sqrt(Constants.g/L)*t[j])
print ('Real',theta)

#Euler
thetadote[0] = thetadot0
thetae[0] = theta0

for j in range(1,n):
    thetadote[j] = h*((-Constants.g/L)*thetae[j-1]) + thetadote[j-1]
    thetae[j] = h*(thetadote[j-1]) + thetae[j-1]
print('Euler',thetae)
#RK4
thetadotr[0] = thetadot0
thetar[0] = theta0

def f(x):
    return (-Constants.g / L) * x

for j in range (1, n):

        k1thetadot = f(thetar[j-1])
        k1theta = thetadotr[j-1]

        k2thetadot =f(thetar[j-1] + (k1thetadot * (h/2)))
        k2theta = thetadotr[j-1] + k1theta * (h/2)

        k3thetadot = f(thetar[j-1] + (k2thetadot * (h/2)))
        k3theta = thetadotr[j-1] + k2theta * (h/2)

        k4thetadot = f(thetar[j-1] + (k3thetadot * (h)))
        k4theta = thetadotr[j-1] + k3theta * (h)

        thetadotr[j] = (thetadotr[j-1] + ((h/6) * (k1thetadot + 2*k2thetadot + 2*k3thetadot + k4thetadot)))
        thetar[j] = (thetar[j-1] + ((h/6) * (k1theta + 2*k2theta + 2*k3theta + k4theta)))


ExactError = ((theta-theta))
print('Exact Error', ExactError)
EulerError = ((thetae-theta))
print('Euler Error', EulerError)
RungeError =((thetar-theta))
print('Runge-Kutta Error',RungeError)

plt.figure(1)
plt.plot(t,theta, 'r-')
plt.plot(t,thetae, 'y-')
plt.plot(t,thetar, 'b-')
axes = plt.gca()
axes.set_xlim([0, tf])
axes.set_ylim([-0.8,0.8])
plt.xlabel("Time (seconds)")
plt.ylabel("X-position (m)")
plt.title("Exact Pendulum Solution")
plt.grid()

plt.figure(2)
plt.plot(t,ExactError, 'r-')
plt.plot(t,EulerError, 'y-')
plt.plot(t,RungeError, 'b-')
axes = plt.gca()
axes.set_xlim([0, tf])
axes.set_ylim([-0.01,0.01])
plt.xlabel("Time (seconds)")
plt.ylabel("Percentage Error")
plt.title("Error Comparison Pendulum")
plt.grid()
plt.show()
