#Jamie Dunlop 02/10/2017
#Solving Kepler's equation for eccentric anomaly
#in order to forecast position and velocity vectors


#def kepler(ecc, M):

import math
M = 3.0491
ecc = 0.4
#Error tolerance
etol = 1e-8

#Initial guess for Eccentric anomaly
if M < math.pi:
    E = M - ecc/2
else:
    E = M + ecc/2

#Iterates Kepler's equation until error is less than tolerance
error = 1
try:
    while abs(error) > etol:
        error = (M - E + ecc * math.sin(E))/ (1 - ecc*math.cos(E))
        E = E + error
except KeyboardInterrupt:
    print('interrupted!')

print (error)
print (E)
#    return
