#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
#import Kep2Cart or correct terminology

fx = []
fy = []

def rk4(f, x, y, h):

    for i in range (0, 200):

        k1 = h * f(x, y)
        k2 = h * f(x + (h/2), y + (k1/2))
        k3 = h * f(x + (h/2), y + (k2/2))
        k4 = h * f(x + h, y + k3)

        fx.append(x + h)
        fy.append(y + ((1/6) * (k1 + 2*k2 + 2*k3 + k4)))

        x = fx[i]
        y = fy[i]

    return fx, fy

def f(x, y):
    return x * math.sqrt(y)

fx, fy = rk4(f, 0, 1, 0.1)

arr = np.array(fx)

print(arr, fy)
