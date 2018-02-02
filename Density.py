import Constants
#US atmospheric density model only applicable when alt > 25000m
def Density1(r):
    alt = r - Constants.Rearth
    Temp = -131.21 + 0.00299 * alt
    Pres = 2.488 * ((Temp + 273.1) / (216.6)) ** -11.388
    return Pres / (0.2869 * ( Temp + 273.1 ))

#NRLMSISE-00
def Density2(r):
    
