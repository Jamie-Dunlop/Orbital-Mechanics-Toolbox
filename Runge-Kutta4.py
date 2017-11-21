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

r_0 = np.array([-10515.45, -5235.37, 49.17]) #m
v_0 = np.array([-2.10305, -4.18146, -5.563290]) #km/s

r = np.linalg.norm(r_0)

mu = 3.985992e14 #units

a = []

def rk4(a, r, r_vec, v_vec, t, h):

    #frx = np.empty((0,1))
    #fry = np.empty((0,1))
    #frz = np.empty((0,1))

    for i in range (0, t):

        if i == 0:

            r_scal = np.linalg.norm(r)

            a_x = (-mu * r_vec[0]) / (r_scal ** 3)
            a_y = (-mu * r_vec[1]) / (r_scal ** 3)
            a_z = (-mu * r_vec[2]) / (r_scal ** 3)

            a = np.array([a_x, a_y, a_z])

            k1v = a * r_vec
            k1r = v_vec
            k2v = a * (r_vec + (k1r * (h/2)))
            k2r = v_vec * k1v * (h/2)
            k3v = a * (r_vec + (k2r * (h/2)))
            k3r = v_vec * k2v * (h/2)
            k4v = a * (r_vec + (k3r * h))
            k4r = v_vec * k3v * h

        else:

            r_scal = np.linalg.norm(r)

            #index_x = 0 + (3)
            #index_y = 1 + (3)
            #index_z = 2 + (3)

            a_x = (-mu * r[i, 0]) / (r_scal ** 3)
            a_y = (-mu * r[i, 1]) / (r_scal ** 3)
            a_z = (-mu * r[i, 2]) / (r_scal ** 3)

            a = np.array([a_x, a_y, a_z])

            k1v = a * r
            k1r = v
            k2v = a * (r + (k1r * (h/2)))
            k2r = v * k1v * (h/2)
            k3v = a * (r + (k2r * (h/2)))
            k3r = v * k2v * (h/2)
            k4v = a * (r + (k3r * h))

            print(r)

        ft.append(t + h)
        fr.append(r_vec + ((1/6) * (k1r + 2*k2r + 2*k3r + k4r)))
        fv.append(v_vec + ((1/6) * (k1v + 2*k2v + 2*k3v + k4v)))

        r = fr[i]
        v = fv[i]

        #frx = np.append(frx, r_vec[0])
        #fry = np.append(fry, r_vec[1])
        #frz = np.append(frz, r_vec[2])

    return fr, fv

fr, fv = rk4(a, r, r_0, v_0, 3, 1)

print(fr)
#ax.plot(fr)
#plot.show()
