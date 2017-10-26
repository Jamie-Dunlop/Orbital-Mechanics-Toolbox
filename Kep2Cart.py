#Anna Middleton 26/10/2016
#Converting Keplerian to cartesian coordinates
#Solving Kepler's equation for eccentric anomaly
#in order to forecast position and velocity vectors


#def kepler(ecc, M):
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def rad(X):
    return ((X * math.pi) / 180)

def deg(X):
    return ((X / math.pi) * 180)

#Test values
a = 26559000 #m
e = 0.704482
i = rad(63.1706)
omega = rad(281.646)
RAAN = rad(206.346)
mu=3.985992e14

#Error tolerance
etol = 1e-8
#period
T = math.pi * 2 * math.sqrt(a**3/mu)
#print ("period", T, "secs")
t = 0
step = 1
try:
    while t <= T:


        #Mean anomaly
        M = math.sqrt(mu/a**3) * t

        #Initial guess for Eccentric anomaly
        if M < math.pi:
            E = M - e/2
        else:
            E = M + e/2

        #Iterates Kepler's equation until error is less than tolerance
        error = 1
        try:
            while abs(error) > etol:
                error = (M - E + e * math.sin(E))/ (1 - e*math.cos(E))
                E = E + error
        except KeyboardInterrupt:
            print('interrupted!')

        print ("Error" ,error)
        print ("Eccentric anomaly", deg(E)) #2.71

        #true anomaly
        V = 2 * math.atan(math.sqrt((1+e) / (1-e) ) * math.tan(E/2))
        print ("True anomaly" , deg(V))

        # rotation matrix
        u = omega + V

        R1 = math.cos(u) * math.cos(RAAN) - math.sin(u) * math.cos(i) * math.sin (RAAN)
        R2 = math.cos(u) * math.sin(RAAN) + math.sin(u) * math.cos(i) * math.cos(RAAN)
        R3 = math.sin(u) * math.sin(i)
        R = np.matrix([[R1], [R2], [R3]])

        print ("Rotation Matrix", R)

        #radius
        r=(a*(1 - e**2))/(1 + (e * math.cos(V)))
        r_comp = r * R
        print ("r", r_comp/1000, "km")

        #velocity
        #v = math.sqrt(mu * ((2/r) - (1/a)))
        #v_comp = v * R
        #print ("v", v_comp)

        v_x = math.sqrt(mu/(a * (1-e**2))) * math.sin(V)
        v_y = math.sqrt(mu/(a * (1-e**2))) * (e + math.cos(V))
        v_mag = math.sqrt(v_x**2 + v_y**2)
        print ("v", (v_mag * R)/1000)
        t = t + step
        print(t)
except KeyboardInterrupt:
    print('interrupted!')


#fig = plt.figure()
#ax = Axes3D(fig)
#X_AXIS = np.arange(-1000000000, 1000000000, 10000000)
#Y_AXIS = np.arange(-1000000000, 1000000000, 10000000)
#Axes3D.autoscale_view(tight=None, scalex=True, scaley=True)
#ax.plot(r_comp[0], r_comp[1])
#plt.show()

ax = fig.gca(projection='3d')
#x = (-20, 20, 1)
#y = (-20, 20, 1)
#z = (-20, 20, 1)

#x_data = [-2959.5096, -1918.8406, 14909.782, 25153.0555, 16005.3093]
#y_data = [2405.9574, -18190.5537, -15084.9452, -5024.1555, 6722.9885]
#z_data = [-6859.6046, 30545.9435, 39809.5606, 30971.9153, 2132.0202]

x_data = [-2959.50962062, -11221.1414687, -11067.62129986, -8637.37019003, -5424.26938833, -1918.84057056, 1661.74817554, 5200.23696173, 8622.4404704, 11874.11811165, 14909.78220786, 17686.22646908, 20157.86009325, 22272.34059496, 23965.3470621, 25153.05546333, 25719.78882302, 25495.42571228, 24209.09833631, 21379.87158701, 16005.30928656, 5567.69841626]
y_data = [2405.95742113, -6972.35729419, -12431.58317202, -15539.17669177, -17308.86852879, -18190.5537088, -18422.31902398, -18147.05073813, -17458.33413722, -16421.57304422, -15084.94518873, -13485.68927331, -11654.10685544, -9616.48024251, -7397.67340042, -5024.15554193, -2528.60535715, 41.34610452, 2603.17452904, 4984.26727277, 6722.98851835, 6007.96977784]
z_data = [-6859.60463185, 2507.61352178, 12314.84908931, 19953.20525759, 25907.99950516, 30545.94345596, 34098.32242557, 36715.41080771, 38497.92789351, 39514.16195435, 39809.56059697, 39412.19168574, 38335.7397503, 36580.84392565, 34135.12004223, 30971.91530367, 27047.62975847, 22297.3841026, 16629.72348835, 9928.48858631, 2132.02015351, -5759.47623658]

#x_data = [r_comp[0]]
#y_data = [r_comp[1]]
#z_data = [r_comp[2]]

ax.plot3D(x_data, y_data, z_data, label='Orbit Path')
ax.legend()

plt.show()
