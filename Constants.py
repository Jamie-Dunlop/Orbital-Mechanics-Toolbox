#https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html

MassEarth = 5.9723e24 #kg https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
GravConst = 6.67384e-11 #m^3/kg*s^2 CODATA
g = 9.798 #m/s^2
J2 = 1082.63e-6

Equatorial_radius = 6378137 #m
Polar_radius = 6356752 #m
Volumetric_mean_radius = 6371008 #m

Rearth = 6371008 # m

Geo = 42240435.255265124 #m Geostationary Radius calculated from T = 2*pi*sqrt(a^3/mu)

mu = MassEarth * GravConst
