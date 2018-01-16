import math
import numpy as np
import matplotlib.pyplot as plt

h = np.linspace(25000, 600000, 1000000) #m

Area = 5 #m^2
BChigh = 1
BClow = 12
mass = 50 #kg
mu = 3.986004418e14 #m something s something
Rearth = 6371008 #m

R = h + Rearth
print(R)

T = -131.21 + 0.00299 * h

p = 2.488 * ((T + 273.1) / (216.6)) ** -11.388

Rho = p / (0.2869 * ( T + 273.1 ))

V = np.sqrt ( mu / R )

Cd = mass / ( BChigh * Area)

Drag = 0.5 * Rho * V ** 2 * Area * Cd

print(Drag)
plt.plot(Drag, h)
plt.xlabel("Drag Force (N)")
plt.ylabel("Altitude (m)")
plt.title("Drag Force Vs Altitude")
plt.grid()
plt.show()
