#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse, Circle

fig = plt.figure()
ax = fig.gca(projection = '3d')

fv = []
fr = []
ft = []

r_vec = np.array([-10515.45, -5235.37, 49.17]) #m
v_vec = np.array([-2.10305, -4.18146, -5.563290]) #km/s

r = np.linalg.norm(r_vec)

mu = 3.985992e14 #units

ax = (-mu * r_vec[0]) / (r ** 3)
ay = (-mu * r_vec[1]) / (r ** 3)
az = (-mu * r_vec[2]) / (r ** 3)

a = np.array([ax, ay, az])

def rk4(a, r_vec, v_vec, t, h):

    #frx = np.empty((0,1))
    #fry = np.empty((0,1))
    #frz = np.empty((0,1))

    for i in range (0, t):

        k1v = a * r_vec
        k1r = v_vec
        k2v = a * (r_vec + (k1r * (h/2)))
        k2r = v_vec * k1v * (h/2)
        k3v = a * (r_vec + (k2r * (h/2)))
        k3r = v_vec * k2v * (h/2)
        k4v = a * (r_vec + (k3r * h))
        k4r = v_vec * k3v * h

        ft.append(t + h)
        fr.append(r_vec + ((1/6) * (k1r + 2*k2r + 2*k3r + k4r)))
        fv.append(v_vec + ((1/6) * (k1v + 2*k2v + 2*k3v + k4v)))

        r = fr[i]
        v = fv[i]

        #frx = np.append(frx, r_vec[0])
        #fry = np.append(fry, r_vec[1])
        #frz = np.append(frz, r_vec[2])

    return fr, fv

fr, fv = rk4(a, r_vec, v_vec, 1000, 1)

print(fr)
#ax.plot(fr)
#plot.show()
