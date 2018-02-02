#Jamie Dunlop and Co 2/2/2018
#Main Simulation Script

import Constants
import RK4+OrbEle
import Six_orbital_elements
import Kep2Cart

#Position and velocity
def StateVec
    r0 = np.array([7371008, 0, 0]) #m
    rdot0 = np.array([0, math.sqrt(Constants.mu/np.linalg.norm(r0)), 0]) #m/s

    return (r0,rdot0)

#Orbital elements
def OrbElm
    a =
    e =
    i =
    omega =
    RAAN =
    V = #true anomaly

    return (a,e,i,omega,RAAN,V)

#Satellite properties
AreaH =
AreaL =
CdH =
CdL =
Mass =


#Integration properties
h = 0.01
t0 = 0
tf = #time or number of orbits
