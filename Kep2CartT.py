#Anna Middleton 26/10/2016
#Converting Keplerian to cartesian coordinates
#Solving Kepler's equation for eccentric anomaly
#in order to forecast position and velocity vectors


def kepler(e,i,omega,RAAN,Mean_motion):
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    import Constants
    # from Main import e,i,omega,RAAN,Mean_motion
    import sys
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.patches import Ellipse, Circle

    # fig = plt.figure()
    # ax = fig.gca(projection='3d')

    def rad(X):
        return (list((i * math.pi) / 180 for i in X))

    def deg(X):
        return ((X / math.pi) * 180)

    mu = Constants.mu
    e = e
    i = rad(i)
    omega = rad(omega)
    RAAN = rad(RAAN)

    mean_motion = list((i * 2 * math.pi) / 86400 for i in Mean_motion) #rad/s

    a = list((mu / i ** 2) ** (1/3) for i in mean_motion)

    #Error tolerance
    etol = 1e-8

    #Period
    T = list(math.pi * 2 * math.sqrt(i ** 3 / mu) for i in a)
    index = 0
    r_resultsx = np.empty((len(T),1))
    r_resultsy = np.empty((len(T),1))
    r_resultsz = np.empty((len(T),1))
    vsat_resultsx = np.empty((len(T),1))
    vsat_resultsy = np.empty((len(T),1))
    vsat_resultsz = np.empty((len(T),1))
    V = np.empty((len(T),1))
    for j in range(0,len(T)):
        print('Simulation {} complete...'.format(index+1))
        t = 0
        step = 0.1
        try:
            while t <= T[index]:

                #Mean anomaly
                M = math.sqrt(mu/a[index]**3) * t

                #Initial guess for Eccentric anomaly
                if M < math.pi:
                    E = M - e[index]/2
                else:
                    E = M + e[index]/2

                #Iterates Kepler's equation until error is less than tolerance
                error = 1
                try:
                    while abs(error) > etol:
                        error = (M - E + e[index] * math.sin(E))/ (1 - e[index]*math.cos(E))
                        E = E + error

                except KeyboardInterrupt:
                    print('interrupted!')

                #True anomaly
                V[index] = 2 * math.atan(math.sqrt((1+e[index]) / (1-e[index]) ) * math.tan(E/2))

                #Rotation matrix

                R11 = math.cos(omega[index]) * math.cos(RAAN[index]) - math.sin(omega[index]) * math.cos(i[index]) * math.sin (RAAN[index])
                R21 = math.cos(omega[index]) * math.sin(RAAN[index]) + math.sin(omega[index]) * math.cos(i[index]) * math.cos(RAAN[index])
                R31 = math.sin(omega[index]) * math.sin(i[index])

                R12=-math.sin(omega[index])*math.cos(RAAN[index])-math.cos(omega[index])*math.cos(i[index])*math.sin(RAAN[index]);
                R22=-math.sin(omega[index])*math.sin(RAAN[index])+math.cos(omega[index])*math.cos(i[index])*math.cos(RAAN[index]);
                R32=math.cos(omega[index])*math.sin(i[index]);

                R13=math.sin(i[index])*math.sin(RAAN[index]);
                R23=-math.sin(i[index])*math.cos(RAAN[index]);
                R33=math.cos(i[index]);

                if e[index]==1:
                    p=2*a[index]
                else:
                    p = a[index]*(1-e[index]**2)

                #Radius
                r=(p)/(1 + (e[index] * math.cos(V[index])))
                xp=r*math.cos(V[index])
                yp=r*math.sin(V[index])
                wom_dot = math.sqrt(mu*p)/r**2

                #velocity
                r_dot = math.sqrt(mu / p) * e[index] * math.sin(V[index])
                vxp=r_dot*math.cos(V[index])-r*math.sin(V[index])*wom_dot
                vyp=r_dot*math.sin(V[index])+r*math.cos(V[index])*wom_dot

                t = t + step

                r_resultsx[index] = R11*xp+R12*yp
                r_resultsy[index] = R21*xp+R22*yp
                r_resultsz[index] = R31*xp+R32*yp

                vsat_resultsx[index] = R11*vxp+R12*vyp
                vsat_resultsy[index] = R21*vxp+R22*vyp
                vsat_resultsz[index] = R31*vxp+R32*vyp

        except KeyboardInterrupt:
            print('interrupted!')

        index += 1
        # sys.exit("Error message")
    # print('r', r_resultsx, r_resultsy, r_resultsz)
    # print('v', vsat_resultsx, vsat_resultsy, vsat_resultsz)
    return([r_resultsx, r_resultsy, r_resultsz],[ vsat_resultsx, vsat_resultsy, vsat_resultsz],V)

# r0 , rdot0, V = kepler()

# #Plot Earth sphere
# u = np.linspace(0, 2 * np.pi, 100)
# vearth = np.linspace(0, np.pi, 100)
# x = 6371000* np.outer(np.cos(u), np.sin(vearth))
# y = 6371000* np.outer(np.sin(u), np.sin(vearth))
# z = 6371000 * np.outer(np.ones(np.size(u)), np.cos(vearth))
# ax.plot_surface(x, y, z, color='g')
#
# #Plot Orbit
# ax.plot(r_resultsx, r_resultsy, r_resultsz, color='g')
# ax.legend()
# #ax.axis('equal')
# max_range = np.array([r_resultsx.max()-r_resultsx.min(), r_resultsy.max()-r_resultsy.min(), r_resultsz.max()-r_resultsz.min()]).max() / 2.0
#
# mid_x = (r_resultsx.max()+r_resultsx.min()) * 0.5
# mid_y = (r_resultsy.max()+r_resultsy.min()) * 0.5
# mid_z = (r_resultsz.max()+r_resultsz.min()) * 0.5
# ax.set_xlim(mid_x - max_range, mid_x + max_range)
# ax.set_ylim(mid_y - max_range, mid_y + max_range)
# ax.set_zlim(-1e7, mid_z + max_range)
# plt.show()
