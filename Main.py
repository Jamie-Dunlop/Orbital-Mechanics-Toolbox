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

#Integration properties
h = 10 #Time step
t0 = 0  #Starting time seconds
tf = 86400 #time or number of orbits
#Satellite properties
Area = 0.0612    #Wetted area m^2
AreaH = 0.7276   #High wetted area m^2
Cd = 2.147
mass = 80  #Mass of Satellite
info = {}
#Position and velocity
def StateVec():
    # r0 = np.array([3169751.48119611, 5583111.16079282, -41650245.77267506]) #m
    r0 = np.array([Constants.Rearth+500000, 0, 0])
    # rdot0 = np.array([-3058.91916149, 257.54617660, -200.43184936]) #m/s
    rdot0 = np.array([0, math.sqrt(Constants.mu / np.linalg.norm(r0)), 0])
    name = 'Satellite'
    return (r0,rdot0,name)

#Orbital elements
with open('TLE_Data.txt') as inf1:
    reader = csv.reader(inf1,delimiter=' ')
    second_col = list(zip(*reader))

with open('NORAD_Satellite_Codes.txt') as inf2:
    reader = csv.reader(inf2,delimiter=',')
    col = list(zip(*reader))

def OrbElm():
    linenum1 = 1
    e = np.zeros([1])
    i = np.zeros([1])
    omega = np.zeros([1])
    RAAN = np.zeros([1])
    Mean_motion = np.zeros([1])
    Norad = np.zeros([1])
    for line in open('TLE_Data.txt'):
        if line.startswith('1'):
            continue
        else:
            e[linenum1-1] = float(second_col[4][linenum1])/10000000
            i[linenum1-1] = float(second_col[2][linenum1]) #degrees
            omega[linenum1-1] = float(second_col[5][linenum1]) #degrees
            RAAN[linenum1-1] = float(second_col[3][linenum1]) #degrees
            Mean_motion[linenum1-1] = float(second_col[7][linenum1]) #revolutions per day
            Norad[linenum1-1] = second_col[1][linenum1] #NORAD ID number of Satellite

            linenum2 = 0
            for line in open('NORAD_Satellite_Codes.txt'):
                linenum2 = linenum2
                if str(int(Norad[linenum1-1])) in line:
                    name = list(col[1])[linenum2]
                linenum2 += 1
            linenum1 += 1
        return (e,i,omega,RAAN,Mean_motion,name)

# def OrbElm():
#     # for line in open('TLE_Data.txt'):
#     #     if line.startswith('1'):
#     #         return line
#     #     else:
#     #         print(list(second_col[4]))
#         e = float(second_col[4][0])/10000000
#         print(e)
#         i = float(second_col[2][0]) #degrees
#         omega = float(second_col[5][0]) #degrees
#         RAAN = float(second_col[3][0]) #degrees
#         Mean_motion = float(second_col[7][0]) #revolutions per day
#         Norad = second_col[1][0] #NORAD ID number of Satellite
#         linenum = 0
#         for line in open('NORAD_Satellite_Codes.txt'):
#             linenum = linenum
#             if Norad in line:
#                 name = list(col[1])[linenum]
#             linenum +=1
#         return (e,i,omega,RAAN,Mean_motion,name)

# Selection of density model to be used
# Density1 = US Model       Density2 =
DensityModel = Density.Density1

# Obtain State Vectors for Satellite
# (r0,rdot0,name) = StateVec()
(e,i,omega,RAAN,Mean_motion,name) = OrbElm()
print('name',name)
import Kep2Cart
r0 = Kep2Cart.r0
rdot0 = Kep2Cart.rdot0


#Call the relevant scripts
RK4Orbittest
