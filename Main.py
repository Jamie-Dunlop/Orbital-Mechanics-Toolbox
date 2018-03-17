#Jamie Dunlop and Co 2/2/2018
#Main Simulation Script
import time
import os
import subprocess
import sys
import Constants
import RK4Orbittest
# import RK4OrbEle
import numpy as np
import matplotlib.pyplot as plt
import math
# import Six_orbital_elements
import Density
import csv
os.system('CLS')

#Initialisation of defining satallite properties and simulation times
h = 10      #Time step (seconds)
t0 = 0      #Starting time (seconds)
tf = 86400 #Finish time (seconds)
#Satellite properties
Area = 0.037    #Nominal area m^2
AreaL = 0.02    #Low area m^2
AreaH = 0.195   #High area m^2
Cd = 2.147      #Drag coefficient
mass = 5        #Mass of Satellite

# Selection of density model to be used
DensityModel = Density.Density1 #US Standard Atmosphere

############## Single Sat obtained from manually entered state vectors ##################

#Position and velocity
def StateVec():
    # r0 = np.array([3169751.48119611, 5583111.16079282, -41650245.77267506]) #m
    r0 = np.array([Constants.Rearth+400000, 0, 0])
    # rdot0 = np.array([-3058.91916149, 257.54617660, -200.43184936]) #m/s
    rdot0 = np.array([0, math.sqrt(Constants.mu / np.linalg.norm(r0)), 0])
    name = 'Satellite'
    return (r0,rdot0,name)

############## Multiple Sats obtained from TLE Data ##################

#Opens Flock2e_TLE_Data.txt and transfers data to 'second_col' variable. File then closed.
with open('Flock2e_TLE_Data.txt') as inf1:
    reader = csv.reader(inf1,delimiter=' ')
    second_col = list(zip(*reader))

#Opens NORAD_Satellite_Codes.txt and transfers data to 'col' variable. File then closed.
with open('NORAD_Satellite_Codes.txt') as inf2:
    reader = csv.reader(inf2,delimiter=',')
    col = list(zip(*reader))

#Orbital elements
def OrbElm():
    #Initialises empty lists for required elements
    e = []
    i = []
    omega = []
    RAAN = []
    Mean_motion = []
    Norad = []
    names = []
    # Fetches orbital elements from Flock2e_TLE_Data.txt for the first (Max(range)-1) satellites
    for b in range (1,3):
        linenum1 = (2*b-1)
        e.append(float(second_col[4][linenum1])/10000000)
        i.append(float(second_col[2][linenum1])) #degrees
        omega.append(float(second_col[5][linenum1])) #degrees
        RAAN.append(float(second_col[3][linenum1])) #degrees
        Mean_motion.append(float(second_col[7][linenum1])) #revolutions per day
        Norad.append(second_col[1][linenum1]) #NORAD ID number of Satellite
        linenum2 = 0
        # Fetches the name of 2019 satellites from NORAD_Satellite_Codes.txt and appends to list. Later used to identify name
        # of chosen satellites in above loop from the NORAD number obtained from TLE number.
        for line in open('NORAD_Satellite_Codes.txt'):
            linenum2 = linenum2
            if str(int(Norad[b-1]))in line:
                names.append(list(col[1])[linenum2])
            linenum2 += 1
    return (e,i,omega,RAAN,Mean_motion,names)

############## Single Sat ##################
# def SingleOrbElm():
#     # for line in open('Flock2e_TLE_Data.txt'):
#     #     if line.startswith('1'):
#     #         return line
#     #     else:
#     #         print(list(second_col[4]))
#         e = float(second_col[4][1])/10000000
#         print(e)
#         i = float(second_col[2][1]) #degrees
#         omega = float(second_col[5][1]) #degrees
#         RAAN = float(second_col[3][1]) #degrees
#         Mean_motion = float(second_col[7][1]) #revolutions per day
#         Norad = second_col[1][1] #NORAD ID number of Satellite
#         linenum = 0
#         for line in open('NORAD_Satellite_Codes.txt'):
#             linenum = linenum
#             if Norad in line:
#                 name = list(col[1])[linenum]
#             linenum +=1
#         return (e,i,omega,RAAN,Mean_motion,name)
################################################

# Obtain State Vectors for Satellite
# (r0,rdot0,name) = StateVec()

######Single########
# (e,i,omega,RAAN,Mean_motion,name) = OrbElm()
#
# print('e', e,i,omega,RAAN,Mean_motion,name)
# import Kep2Cart
# r0 = Kep2Cart.r0
# rdot0 = Kep2Cart.rdot0
#########################

######Multi########
(e,i,omega,RAAN,Mean_motion,names) = OrbElm()
TrueAnomaly = []
print('Please standby, running simulations for the following: {}'.format(names))
import Kep2CartT
r0vec,rdot0vec,V = Kep2CartT.kepler(e,i,omega,RAAN,Mean_motion)
for b in range (0,len(r0vec[0])):
    r0temp = np.array([r0vec[0][b],r0vec[1][b],r0vec[2][b]]).tolist() #arrays to lists
    r0 = [y for x in r0temp for y in x] #List in a list to list
    print('r0',r0)
    print('r0 mod',np.linalg.norm(r0)-Constants.Rearth)

    rdot0temp = np.array([rdot0vec[0][b],rdot0vec[1][b],rdot0vec[2][b]]).tolist()
    rdot0 = [y for x in rdot0temp for y in x]
    print('rdot0',rdot0)
    print('rdot0 mod',np.linalg.norm(rdot0))

    name = names[b]
    
##################################

## Call the relevant scripts

    (V,t) = RK4Orbittest.wholefile(r0,rdot0,t0,h,tf,name,Area,AreaH,AreaL,Cd,mass,DensityModel,b)
    TrueAnomaly.append(V) #True anomaly at each step for every satellite run
    print('V',TrueAnomaly[b])

plt.plot(t,TrueAnomaly[0],'r-')
plt.plot(t,TrueAnomaly[1],'b-')
plt.xlabel("Time (s)")
plt.ylabel("True Anomaly (radians)")
plt.title("True Anomaly vs Time of {} satellies".format(len(r0vec[0])))
plt.grid()
plt.show()
