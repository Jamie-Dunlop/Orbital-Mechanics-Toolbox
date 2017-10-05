#Jamie Dunlop 02/10/2017
#Starting with position and velocity vectors, the six orbital elements
#that fully define an orbit are calculated
import math
import numpy as np
#Dot product function
def dot(A, B):
    return (A[0] * B [0] + A[1] * B[1] + A[2] * B[2])

#Cross product function
def cross(A, B):
    return [(A[1] * B[2] - A[2] * B[1]), (A[2] * B[0] - A[0] * B[2]), (A[0] * B[1] - A[1] * B[0])]

#Modulus of vector function
def mod(vector):
    sum = 0
    for i in vector:
        sum += i ** 2
    return math.sqrt(sum)


def orb_elems(r, v, mu):

    #Orbital energy
    epsilon = 0.5 * mod(v) ** 2 - mu/mod(r)
    #Semi-major axis
    a = -mu / (2 * epsilon)

    #Angular Momentum vector
    H = cross(r, v)
    #Eccenticity vector
    e_v = (1 / mu) * (cross(v, H) - (mu * np.array(r))/mod(r))
    e = mod(e_v)

    #Unit vector along direction of ascending nodes
    k = [0,0,1] #Unit vector in the inertial frame
    n = np.array((cross(k, H)))/mod(cross(k, H))

    #Right ascension of the ascending node
    RAAN = math.atan((H[0] / mod(H)) / -(-H[1] / mod(H))) * 180 / math.pi


    #Inclination
    i = (math.acos((dot(k, H))/mod(H))) * 180 / math.pi

    #Argument of perigee
    omega = math.acos((dot(n, e_v))/mod(e_v))
    if dot(e_v, k) < 0:
        omega = 2 * math.pi - omega

    #True anomaly
    V = math.acos((dot(e_v, r))/(e * mod(r)))
    if dot(r, v) < 0:
       V = 2 * math.pi - V

    print("Orbit energy: ", epsilon)
    print("Angular Momentum vector: ", H)
    print()
    print()
    print("Semi-major axis: ", a)
    print("Eccenticity: ", e)
    print("Inclination: ", i, "deg")
    print("Argument of the perigee: ", omega * 180 / math.pi, "deg")
    print("True anomaly: ", V * 180 / math.pi, "deg")
    print("RAAN: ", RAAN, "deg")

    return

def vectors(Vel, height, gamma):
    v_tan = Vel * math.sin(gamma * math.pi / 180)
    v_rad = Vel * math.cos(gamma * math.pi / 180)
    r_x = height * math.cos(115 * math.pi / 180)
    r_y = height * math.sin(115 * math.pi / 180)

    print ("v_tan", v_tan)
    print ("v_rad", v_rad)
    print ("r_x", r_x)
    print ("r_y", r_y)


    return v_tan, v_rad, r_x, r_y

Vel = 7605
gamma = 0.45
height = 495

v_tan, v_rad, r_x, r_y = vectors(Vel, height, gamma)

#orb_elems([r_x, r_y, -1.405667119e-4], [v_tan, v_rad, -4.05058737e-5], 3.986004e14)
orb_elems([0, -7000, -14000], [2, 1, -0.5], 398600) #Molniya orbit 
# r = [4.1852e7, 6.2778e7, 10.463e7]
# v = [2.5936e4, 5.1872e4, 0]
# mu = 1.40812
