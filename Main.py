#Jamie Dunlop and Co 2/2/2018
#Main Simulation Script
import time
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

#Integration properties
h = 10 #Time step
t0 = 0  #Starting time seconds
tf = 86400 #time or number of orbits
#Satellite properties
Area = 0.037
AreaL = 0.02    #Wetted area m^2
AreaH = 0.195   #High wetted area m^2
Cd = 2.147
mass = 5  #Mass of Satellite
#Position and velocity
def StateVec():
    # r0 = np.array([3169751.48119611, 5583111.16079282, -41650245.77267506]) #m
    r0 = np.array([Constants.Rearth+400000, 0, 0])
    # rdot0 = np.array([-3058.91916149, 257.54617660, -200.43184936]) #m/s
    rdot0 = np.array([0, math.sqrt(Constants.mu / np.linalg.norm(r0)), 0])
    name = 'Satellite'
    return (r0,rdot0,name)

with open('Flock2e_TLE_Data.txt') as inf1:
    reader = csv.reader(inf1,delimiter=' ')
    second_col = list(zip(*reader))

with open('NORAD_Satellite_Codes.txt') as inf2:
    reader = csv.reader(inf2,delimiter=',')
    col = list(zip(*reader))

#### Multiple Sats######
#Orbital elements
def OrbElm():
    e = []
    i = []
    omega = []
    RAAN = []
    Mean_motion = []
    Norad = []
    names = []
    for b in range (1,3):
        linenum1 = (2*b-1)
        e.append(float(second_col[4][linenum1])/10000000)
        i.append(float(second_col[2][linenum1])) #degrees
        omega.append(float(second_col[5][linenum1])) #degrees
        RAAN.append(float(second_col[3][linenum1])) #degrees
        Mean_motion.append(float(second_col[7][linenum1])) #revolutions per day
        Norad.append(second_col[1][linenum1]) #NORAD ID number of Satellite
        linenum2 = 0
        for line in open('NORAD_Satellite_Codes.txt'):
            linenum2 = linenum2
            if str(int(Norad[b-1]))in line:
                names.append(list(col[1])[linenum2])
            linenum2 += 1
    return (e,i,omega,RAAN,Mean_motion,names)

##### Single Sat#########
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

# Selection of density model to be used
# Density1 = US Model       Density2 =
DensityModel = Density.Density1

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
print('e', e,i,omega,RAAN,Mean_motion,names)
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
    print(name)

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
