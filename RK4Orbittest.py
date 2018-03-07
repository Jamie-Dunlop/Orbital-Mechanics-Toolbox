#Chris Lowther 06/11/2017
#Runge-Kutta 4th order

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse, Circle
import Constants
import Density
from Six_orbital_elements import dot, cross, mod, orb_elems

#Force Model including gravity and drag
def Accel(t,R,V,Area,AreaH,AreaL,Cd,mass,DensityModel,b,tf):
    if b == 1:
        if t > tf/4:
            Area = AreaH
        else:
            Area = Area
    Gravity = ((-Constants.mu) * (R)) / np.linalg.norm(R) ** 3 #monopole gravity model?
    Drag = - (0.5 * DensityModel(np.linalg.norm(R))*np.linalg.norm(V) * Area * Cd*V) / mass
    return Gravity + Drag

#Gets position and velocity from state vector and calculates acceleration from Accel
def Orbit(t, state,Area,AreaH,AreaL,Cd,mass,DensityModel,b,tf):
    pos, vel = state
    return np.array([vel, Accel(t,pos,vel,Area,AreaH,AreaL,Cd,mass,DensityModel,b,tf)])

#Runge-Kutta 4 Integrator
def rk4(x, h, y, f,Area,AreaH,AreaL,Cd,mass,DensityModel,b,tf):
    k1 = h * f(x, y,Area,AreaH,AreaL,Cd,mass,DensityModel,b,tf)
    k2 = h * f(x + 0.5*h, y + 0.5*k1,Area,AreaH,AreaL,Cd,mass,DensityModel,b,tf)
    k3 = h * f(x + 0.5*h, y + 0.5*k2,Area,AreaH,AreaL,Cd,mass,DensityModel,b,tf)
    k4 = h * f(x + h, y + k3,Area,AreaH,AreaL,Cd,mass,DensityModel,b,tf)
    return x + h, y + (k1 + 2*(k2 + k3) + k4)/6.0


def wholefile(r0,rdot0,t0,h,tf,name,Area,AreaH,AreaL,Cd,mass,DensityModel,b):

  fig = plt.figure()
  ax = fig.gca(projection='3d')
  n = (tf-t0)//h
  n = int(n)

  #Period
  T = 2 * math.pi * math.sqrt (np.linalg.norm(r0) ** 3 / Constants.mu)

  No_of_orbits = int(tf/T)
  limit = 7e6 #Graph limits for plotting
  print('Period', T)
  x= np.zeros([n])
  y= np.zeros([n])
  z= np.zeros([n])
  V = np.zeros([n])
  xdot= np.zeros([n])
  ydot= np.zeros([n])
  zdot= np.zeros([n])
  x[0] = r0[0]
  y[0] = r0[1]
  z[0] = r0[2]
  xdot[0] = rdot0[0]
  ydot[0] = rdot0[1]
  zdot[0] = rdot0[2]
  r = [[]]
  rmod = np.zeros([n])
  rdot = [[]]
  diff = np.zeros([No_of_orbits])
  Apogee = np.zeros([No_of_orbits])
  # diff[0] = 0
  Orbit_no = np.arange(1, No_of_orbits+1)
  # rdot[0] = Mainrdot0
  rmod[0] = np.linalg.norm(r0)
  # rdot.insert(1, Mainrdot0)

  state = np.array([r0,rdot0])
  t = np.linspace(t0,tf,n)

  #Starts RK4 and stores position and veloctiy values obtained from RK4
  timed = 0
  Count = 1
  for j in range(1,n):
      timed, state = rk4(timed, h, state, Orbit,Area,AreaH,AreaL,Cd,mass,DensityModel,b,tf)
      # if timed < Count*T:
      #     print(Count)
      x[j] = state[0][0]
      y[j] = state[0][1]
      z[j] = state[0][2]
      xdot[j] = state[1][0]
      ydot[j] = state[1][1]
      zdot[j] = state[1][2]
      r = [x,y,z]
      rdot = [xdot,ydot,zdot]
      vmod = math.sqrt(xdot[j]**2+ydot[j]**2+zdot[j]**2)
      rmod[j] = math.sqrt(x[j]**2+y[j]**2+z[j]**2)

      #True Anomaly at every step from Six_orbital_elements
      a = orb_elems([x[j],y[j],z[j]], [xdot[j],ydot[j],zdot[j]], Constants.mu)
      V[j-1] =  a[5]

      if rmod[j] < Constants.Rearth:
          t = t[:j]
          rmod = rmod[:j]
          No_of_orbits = int(t[j-1]/T)
          diff = np.zeros([No_of_orbits])
          n = (t[j-1]-t0)//h
          n = int(n)
          x= x[:j]
          y= y[:j]
          z= z[:j]
          Orbit_no = np.arange(1, No_of_orbits+1)
          print('***WARNING: Satellite has crashed!***')
          break
          # Apogee[(Count-1)] = max(rmod)
          # Count +=1
          # if Count == 15:
          #     break
          # print('Apogee',Apogee)

  ## Prints Rmod at the start and end of an orbit for error check
  for jj in range(0,No_of_orbits):
  # for jj in range(0,int(Maintf/T)):
      # print('R at start of orbit',jj+1,'>>',rmod[int((jj)*T/Mainh)],'m')
      # print('R at end of orbit',jj+1,'>>',rmod[int((jj+1)*T/Mainh)],'m')
      diff[jj] = rmod[int((jj)*T/h)]-rmod[0]
      # print('Difference from R at start >>', diff)

  print('First', rmod[0])
  # print('last',rmod[n-1])
  # print('Ratio',((rmod[int((Maintf/Mainh)-1)]-rmod[int(Maintf/(2*Mainh))])/(rmod[int(Maintf/(2*Mainh))]-rmod[0])))
  # # print ('r', r)
  # # print('rdot', rdot)
  # # print (x, y)
  # plt.plot(x,y)
  # plt.subplot(2,1,1)
  # plt.plot(x,y)
  # plt.xlabel("X-position (m)")
  # plt.ylabel("Y-position (m)")
  # plt.title("RK4 - Orbit")
  # plt.grid()
  #
  # plt.subplot(2,1,2)
  # plt.plot(t,rmod)
  # plt.xlabel("Time (s)")
  # plt.ylabel("Radius (m)")
  # plt.title("Radisu change over time")
  # plt.grid()
  #
  # plt.show()

  # print ('r', r)
  # print('rdot', rdot)
  # print (x, y)
  print('name',name)
  # Plot Orbit
  plt.figure(1)
  ax.plot(x, y, z, color = 'r')
  ax.set_xlabel("X-position (m)")
  ax.set_ylabel("Y-position (m)")
  ax.set_zlabel("Z-position (m)")
  plt.title("RK4 - Orbit of {}".format(name))
  plt.grid()

  # Plot Earth sphere
  u = np.linspace(0, 2 * np.pi, 100)
  vearth = np.linspace(0, np.pi, 100)
  xe = Constants.Rearth * np.outer(np.cos(u), np.sin(vearth))
  ye = Constants.Rearth * np.outer(np.sin(u), np.sin(vearth))
  ze = Constants.Rearth * np.outer(np.ones(np.size(u)), np.cos(vearth))
  ax.plot_surface(xe, ye, ze, color='b')

  ax.set_xlim(-limit, limit)
  ax.set_ylim(-limit, limit)
  ax.set_zlim(-limit, limit)
  # ax.set_aspect('equal')

  plt.figure(2)
  plt.plot(t,(rmod-Constants.Rearth)/1000)
  plt.xlabel("Time (s)")
  plt.ylabel("Altitude (km)")
  plt.title("Altitude vs Time of {}".format(name))
  plt.grid()

  plt.figure(3)
  plt.plot(Orbit_no, diff)
  plt.xlabel("Orbit Number")
  plt.ylabel("Difference in Radius (m)")
  plt.title("Radius Difference from Original Orbit of {}".format(name))
  plt.grid()

  plt.figure(4)
  plt.plot(x,y)
  plt.xlabel("X (m)")
  plt.ylabel("Y (m)")
  plt.title("2D Plot of {}".format(name))
  plt.grid()
  plt.axes().set_aspect('equal')



  # xmod = np.linalg.norm(x)
  # ymod = np.linalg.norm(y)
  # zmod = np.linalg.norm(z)

  # max_range = np.array([xmod.max()-xmod.min(), ymod.max()-ymod.min(), zmod.max()-zmod.min()]).max() / 2.0
  #
  # mid_x = (xmod.max()+xmod.min()) * 0.5
  # mid_y = (ymod.max()+ymod.min()) * 0.5
  # mid_z = (zmod.max()+zmod.min()) * 0.5
  # ax.set_xlim(mid_x - max_range, mid_x + max_range)
  # ax.set_ylim(mid_y - max_range, mid_y + max_range)
  # ax.set_zlim(mid_z - max_range, mid_z + max_range)

  plt.show()
  return V,t
