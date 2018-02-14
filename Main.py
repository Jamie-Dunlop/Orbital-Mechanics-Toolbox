#Jamie Dunlop and Co 2/2/2018
#Main Simulation Script
import time
import Constants
import RK4Orbit
# import RK4OrbEle
import numpy as np
import math
import Six_orbital_elements
import Density

#Integration properties
h = 5 #Time step
t0 = 0  #Starting time seconds
tf = 1000000 #time or number of orbits
#Satellite properties
Area = 0.0612    #Wetted area m^2
AreaH = 0.7276   #High wetted area m^2
Cd = 2.147
mass = 50   #Mass of Satellite

#Position and velocity
def StateVec():
    # r0 = np.array([3169751.48119611, 5583111.16079282, -41650245.77267506]) #m
    r0 = np.array([Constants.Geo+77500, 0, 0])
    # rdot0 = np.array([-3058.91916149, 257.54617660, -200.43184936]) #m/s
    rdot0 = np.array([0, math.sqrt(Constants.mu / np.linalg.norm(r0)), 0])
    return (r0,rdot0)

#Orbital elements
def OrbElm():
    e = 0.0005304
    i = 277.9340 #degrees
    omega = 86.2686 #degrees
    RAAN = 355.7087 #degrees
    Mean_motion = 1.00272536 #revolutions per day

    return (e,i,omega,RAAN,Mean_motion)

# Selection of density model to be used
# Density1 = US Model       Density2 =
DensityModel = Density.Density1

# Obtain State Vectors for Satellite
(r0,rdot0) = StateVec()
# (e,i,omega,RAAN,Mean_motion) = OrbElm()
# import Kep2Cart
# (r0, rdot0) = Kep2Cart

#Call the relevant scripts
RK4Orbit
