#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math

def rk4(f, x0, y0, x1, h):

    x[0] = x = x0
    y[0] = y = y0

    for i in range (0, n + 1):

        k1 = h * f(x, y)
        k2 = h * f(x + (h/2), y + (k1/2))
        k3 = h * f(x + (h/2), y + (k2/2))
        k4 = h * f(x + h, y + k3)

        x[i] = x = x0 + (i * h)
        y[i] = y = y0 + ((1/6) * (k1 + 2*k2 + 2*k3 + k4))

    return x, y

def f(x, y):
    #return (5 * x**2 - y) / math.exp(x + y)
    return x * math.sqrt(y)

x, y = rk4(f, 0, 1, 20, 0.1)
for x, y in list (zip (x, y)) [::10]:

    # within %number.decimalf
    # 'decimal' controls the number of decimal places
    # 'number' controls the spacing between printed numbers when printed
    # can write %+10.1e to give in scientific or %10.1f to give in decimal

    print (" %4.1f %10.3f %12.3f" % (x, y, x * math.sqrt(y)))
