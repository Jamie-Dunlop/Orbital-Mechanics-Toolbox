import math
import numpy as np
import matplotlib.pyplot as plt

t0 = 0
tf = 500
x0 = -0.8
m = 1
k = 2
c = 0.4
xdot0 = 0
h = 0.001
n = (tf-t0)//h
n = int(n)
t = np.linspace(t0,tf,n)
x = np.zeros([n])
xr = np.zeros([n])
xdot = np.zeros([n])
xe = np.zeros([n])
xedot = np.zeros([n])

Natfreq = math.sqrt(k/m)
Critdamp = c/(2*Natfreq*m)
Dampfreq = Natfreq * math.sqrt(1-Critdamp**2)
A = x0
B = (x0 * Critdamp) / (math.sqrt(1-Critdamp**2))

#Exact Solution
for j in range(1,n):
    x[j] = math.exp(-Critdamp*Natfreq*t[j])*(A*math.cos(Dampfreq*t[j]) + B*math.sin(Dampfreq*t[j]))
print (x)

#Euler
xedot[0] = xdot0
xe[0] = x0

for j in range(1,n):
    xedot[j] = h*((-c*xedot[j-1]-k*xe[j-1])/m) + xedot[j-1]
    xe[j] = h*(xedot[j-1]) + xe[j-1]

#RK4
xdot[0] = xdot0
xr[0] = x0

def f(a,b):
    return (-c*a-k*b)/m

for j in range (1, n):

        k1xdot = f(xdot[j-1], xr[j-1])
        k1x = xdot[j-1]

        k2xdot =f(xdot[j-1] + (k1xdot), xr[j-1] + (k1x))
        k2x = xdot[j-1] + k1xdot

        k3xdot = f(xdot[j-1] + (k2xdot), xr[j-1] + (k2x))
        k3x = xdot[j-1] + k2xdot

        k4xdot = f(xdot[j-1] + (k3xdot), xr[j-1] + (k3x))
        k4x = xdot[j-1] + k3xdot

        xdot[j] = (xdot[j-1] + ((h/6) * (k1xdot + 2*k2xdot + 2*k3xdot + k4xdot)))
        xr[j] = (xr[j-1] + ((h/6) * (k1x + 2*k2x + 2*k3x + k4x)))


ExactError = ((x-x)/x)*100
print('Exact Error', ExactError)
EulerError = ((xe-x)/x)*100
print('Euler Error', EulerError)
RungeError =((xr-x)/x)*100
print('Runge-Kutta Error',RungeError)

plt.figure(1)
plt.plot(t,x, 'r-')
plt.plot(t,xe, 'y-')
plt.plot(t,xr, 'b-')
axes = plt.gca()
axes.set_xlim([0, 500])
axes.set_ylim([-0.8,0.8])
plt.xlabel("Time (seconds)")
plt.ylabel("X-position (m)")
plt.title("Exact Spring Mass Damper System")
plt.grid()

plt.figure(2)
plt.plot(t,ExactError, 'r-')
plt.plot(t,EulerError, 'y-')
plt.plot(t,RungeError, 'b-')
plt.xlabel("Time (seconds)")
plt.ylabel("Percentage Error")
plt.title("Error Comparison Spring Mass Damper")
plt.grid()
plt.show()
