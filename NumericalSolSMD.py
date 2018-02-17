import math
import numpy as np
import matplotlib.pyplot as plt

t0 = 0
tf = 50
x0 = -0.8
m = 1
k = 1
c = 0.05
xdot0 = 0
h = 0.01
n = (tf-t0)//h
n = int(n)
state = np.array([-0.8, 0])
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
x[0] = x0
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

def damped_spring(t, state):
    pos, vel = state
    stiffness = k
    damping = c
    return np.array([vel, -stiffness*pos - damping*vel])

def rk4(x, h, y, f):
    k1 = h * f(x, y)
    k2 = h * f(x + 0.5*h, y + 0.5*k1)
    k3 = h * f(x + 0.5*h, y + 0.5*k2)
    k4 = h * f(x + h, y + k3)
    return x + h, y + (k1 + 2*(k2 + k3) + k4)/6.0

timed = 0
xr[0] = x0
for jr in range(1,n):
    timed, state = rk4(timed, h, state, damped_spring)
    xr[jr] = state[0]


print('state',xr)
ExactError = x-x
print('Exact Error', ExactError)
EulerError = xe-x
print('Euler Error', EulerError)
RungeError =xr-x
print('Runge-Kutta Error',RungeError)

plt.figure(1)
plt.plot(t,x, 'r-')
plt.plot(t,xe, 'y-')
plt.plot(t,xr, 'b-')
axes = plt.gca()
axes.set_xlim([0, tf])
axes.set_ylim([-0.8,0.8])
plt.xlabel("Time (seconds)")
plt.ylabel("X-position (m)")
plt.title("Exact Spring Mass Damper System")
plt.grid()

plt.figure(2)
plt.plot(t,ExactError, 'r-')
plt.plot(t,EulerError, 'y-')
plt.plot(t,RungeError, 'b-')
axes = plt.gca()
axes.set_xlim([0, tf])
axes.set_ylim([-0.1,0.1])
plt.xlabel("Time (seconds)")
plt.ylabel("Absolute Error")
plt.title("Error Comparison Spring Mass Damper")
plt.grid()
plt.show()
