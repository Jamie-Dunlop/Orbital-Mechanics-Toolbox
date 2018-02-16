#Jamie Dunlop and Co 2/2/2018
#Main Simulation Script
import time
import Constants
import RK4Orbittest
# import RK4OrbEle
import numpy as np
import math
# import Six_orbital_elements
import Density
import csv

for h in (1,10):
    #Integration properties
    h = h #Time step
    t0 = 0  #Starting time seconds
    tf = 200000 #time or number of orbits
    #Satellite properties
    Area = 0.0612    #Wetted area m^2
    AreaH = 0.7276   #High wetted area m^2
    Cd = 2.147
    mass = 50   #Mass of Satellite

    #Position and velocity
    def StateVec():
        # r0 = np.array([3169751.48119611, 5583111.16079282, -41650245.77267506]) #m
        r0 = np.array([Constants.Geo, 0, 0])
        # rdot0 = np.array([-3058.91916149, 257.54617660, -200.43184936]) #m/s
        rdot0 = np.array([0, math.sqrt(Constants.mu / np.linalg.norm(r0)), 0])
        return (r0,rdot0)

    #Orbital elements
    with open('TLE_Datatest.txt') as inf:
        reader = csv.reader(inf,delimiter=' ')
        second_col = list(zip(*reader))

    with open('NORAD_Satellite_Codes.txt') as inf:
        reader = csv.reader(inf,delimiter=' ')
        col = list(zip(*reader))

    def OrbElm():
        e = float(list(second_col[4])[0])/10000000
        i = float(list(second_col[2])[0]) #degrees
        omega = float(list(second_col[5])[0]) #degrees
        RAAN = float(list(second_col[3])[0]) #degrees
        Mean_motion = float(list(second_col[7])[0]) #revolutions per day
        Norad = list(second_col[1])[0] #NORAD ID number of Satellite
        return (e,i,omega,RAAN,Mean_motion)

    # Selection of density model to be used
    # Density1 = US Model       Density2 =
    DensityModel = Density.Density1

    # Obtain State Vectors for Satellite
    # (r0,rdot0) = StateVec()
    (e,i,omega,RAAN,Mean_motion) = OrbElm()
    import Kep2Cart
    r0 = Kep2Cart.r0
    rdot0 = Kep2Cart.rdot0


    #Call the relevant scripts
    RK4Orbittest
