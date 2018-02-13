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
start = time.time()

#Integration properties
h = 0.1 #Time step
t0 = 0  #Starting time seconds
tf = 100000 #time or number of orbits
#Satellite properties
Area = 0.0612    #Wetted area m^2
AreaH = 0.7276   #High wetted area m^2
Cd = 2.147
mass = 50   #Mass of Satellite

#Position and velocity
def StateVec():
    r0 = np.array([3169751.48119611, 5583111.16079282, -41650245.77267506]) #m
    rdot0 = np.array([-3058.91916149, 257.54617660, -200.43184936]) #m/s
    return (r0,rdot0)



#Orbital elements
def OrbElm():
    a = 0
    e = 0
    i = 0
    omega = 0
    RAAN = 0
    V = 0 #true anomaly

    return (a,e,i,omega,RAAN,V)

# Selection of density model to be used
# Density1 = US Model       Density2 =
DensityModel = Density.Density1

# Obtain State Vectors for Satellite
(r0,rdot0)=StateVec()

#Call the relevant scripts
RK4Orbit

end = time.time()
print(end - start)
