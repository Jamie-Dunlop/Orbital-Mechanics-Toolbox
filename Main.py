#Jamie Dunlop and Co 2/2/2018
#Main Simulation Script

import Constants
import RK4Orbit
# import RK4+OrbEle
import numpy as np
import math
import Six_orbital_elements


#Integration properties
h = 10  #Time step
t0 = 0  #Starting time seconds
tf = 86400 #time or number of orbits

#Satellite properties
Area = 5    #Wetted area
AreaH =10   #High wetted area
BChigh = 1  #High Baliistic coefficient
BClow = 12  #Low Baliistic coefficient
mass = 50   #Mass of Satellite

#Position and velocity
def StateVec():
    r0 = np.array([7371008, 0, 0]) #m
    rdot0 = np.array([0, math.sqrt(Constants.mu/np.linalg.norm(r0)), 0]) #m/s

    return (r0,rdot0)

#Orbital elements
def OrbElm():
    a =0
    e =0
    i =0
    omega =0
    RAAN =0
    V = 0#true anomaly

    return (a,e,i,omega,RAAN,V)




(r0,rdot0)=StateVec()
RK4Orbit
