#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 12:03:40 2020

@author: stephen
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

a = 0.25
b = 4
F = 6  #28 #6
G = 0.8 #0.8 #1.1

timestep = 0.01
nstep = 15000
nTries = 7
transient_buffer = 200

def Xdot(X, Y, Z, t):
    return - Y**2 - Z**2 - a*X + a*F
 
def Ydot(X, Y, Z, t):
    return X*Y - b*X*Z  - Y + G

def Zdot(X, Y, Z, t):
    return b*X*Y + X*Z - Z


x1plot = []
x2plot = []
x3plot = []
# x4plot = []

x1BackPlot = []
x2BackPlot = []
x3BackPlot = []
# x4BackPlot = []

x1FrontPlot = []
x2FrontPlot = []
x3FrontPlot = []
# x4FrontPlot = []

### Perform 4th Order Runge-Kutta Integration:

def gle4_rk4(x1, x2, x3, t):

    for i in range(nstep - 1):
         t[i+1] = t[i] + timestep
    
         ka1 = Xdot( x1[i] , x2[i] , x3[i] , t[i] )
         ka2 = Xdot( ( (x1[i]) + (timestep/2)*(ka1) ) , ( (x2[i]) + (timestep/2)*(ka1) ) , ( (x3[i]) + (timestep/2)*(ka1) ) , ( (t[i]) + (timestep/2) ) )
         ka3 = Xdot( ( (x1[i]) + (timestep/2)*(ka2) ) , ( (x2[i]) + (timestep/2)*(ka2) ) , ( (x3[i]) + (timestep/2)*(ka2) ) , ( (t[i]) + (timestep/2) ) )
         ka4 = Xdot( ( (x1[i]) + (timestep)*(ka3) ) , ( (x2[i]) + (timestep)*(ka3) ), ( (x3[i]) + (timestep)*(ka3) ), ( (t[i]) + (timestep) ) )
    
         kb1 = Ydot( x1[i] , x2[i] , x3[i] , t[i] )
         kb2 = Ydot( ( (x1[i]) + (timestep/2)*(kb1) ) , ( (x2[i]) + (timestep/2)*(kb1) ) , ( (x3[i]) + (timestep/2)*(kb1) ) , ( (t[i]) + (timestep/2) ) )
         kb3 = Ydot( ( (x1[i]) + (timestep/2)*(kb2) ) , ( (x2[i]) + (timestep/2)*(kb2) ) , ( (x3[i]) + (timestep/2)*(kb2) ) , ( (t[i]) + (timestep/2) ) )
         kb4 = Ydot( ( (x1[i]) + (timestep)*(kb3) ) , ( (x2[i]) + (timestep)*(kb3) ), ( (x3[i]) + (timestep)*(kb3) ) , ( (t[i]) + (timestep) ) )
    
         kc1 = Zdot( x1[i] , x2[i] , x3[i] , t[i] )
         kc2 = Zdot( ( (x1[i]) + (timestep/2)*(kc1) ) , ( (x2[i]) + (timestep/2)*(kc1) ) , ( (x3[i]) + (timestep/2)*(kc1) ) , ( (t[i]) + (timestep/2) ) )
         kc3 = Zdot( ( (x1[i]) + (timestep/2)*(kc2) ) , ( (x2[i]) + (timestep/2)*(kc2) ) , ( (x3[i]) + (timestep/2)*(kc2) ) , ( (t[i]) + (timestep/2) ) )
         kc4 = Zdot( ( (x1[i]) + (timestep)*(kc3) ) , ( (x2[i]) + (timestep)*(kc3) ), ( (x3[i]) + (timestep)*(kc3) ) , ( (t[i]) + (timestep) ) )
    
         # kd1 = x4dot( x1[i] , x2[i] , x3[i] , x4[i] , t[i] )
         # kd2 = x4dot( ( (x1[i]) + (timestep/2)*(kd1) ) , ( (x2[i]) + (timestep/2)*(kd1) ) , ( (x3[i]) + (timestep/2)*(kd1) ) , ( (x4[i]) + (timestep/2)*(kd1) ) , ( (t[i]) + (timestep/2) ) )
         # kd3 = x4dot( ( (x1[i]) + (timestep/2)*(kd2) ) , ( (x2[i]) + (timestep/2)*(kd2) ) , ( (x3[i]) + (timestep/2)*(kd2) ) , ( (x4[i]) + (timestep/2)*(kd2) ) , ( (t[i]) + (timestep/2) ) )
         # kd4 = x4dot( ( (x1[i]) + (timestep)*(kd3) ) , ( (x2[i]) + (timestep)*(kd3) ), ( (x3[i]) + (timestep)*(kd3) ) , ( (x4[i]) + (timestep)*(kd3) ) , ( (t[i]) + (timestep) ) )
    
         x1[i+1] = x1[i] + (timestep/6)*(ka1 + 2*ka2 + 2*ka3 + ka4)
         x2[i+1] = x2[i] + (timestep/6)*(kb1 + 2*kb2 + 2*kb3 + kb4)
         x3[i+1] = x3[i] + (timestep/6)*(kc1 + 2*kc2 + 2*kc3 + kc4)
         # x4[i+1] = x4[i] + (timestep/6)*(kd1 + 2*kd2 + 2*kd3 + kd4)
    
         if ( i >= transient_buffer ):
              x1plot.append(x1[i+1])
              x2plot.append(x2[i+1])
              x3plot.append(x3[i+1])
              # x4plot.append(x4[i+1])
              if Xdot(x1[i+1], x2[i+1], x3[i+1], t) < 0:
                  x1BackPlot.append(x1[i+1])
                  x2BackPlot.append(x2[i+1])
                  x3BackPlot.append(x3[i+1])
                  # x4BackPlot.append(x4[i+1])
              elif Xdot(x1[i+1], x2[i+1], x3[i+1], t) > 0:
                  x1FrontPlot.append(x1[i+1])
                  x2FrontPlot.append(x2[i+1])
                  x3FrontPlot.append(x3[i+1])
                  # x4FrontPlot.append(x4[i+1])

    
    BackPlot = [[] for i in range(3)]
    FrontPlot = [[] for i in range(3)]
    MasterPlot = [[] for i in range(3)]
    
    BackPlot[0] = x1BackPlot
    BackPlot[1] = x2BackPlot
    BackPlot[2] = x3BackPlot
    # BackPlot[3] = x4BackPlot
    
    FrontPlot[0] = x1FrontPlot
    FrontPlot[1] = x2FrontPlot
    FrontPlot[2] = x3FrontPlot
    # FrontPlot[3] = x4FrontPlot
    
    MasterPlot[0] = x1plot
    MasterPlot[1] = x2plot
    MasterPlot[2] = x3plot
    # MasterPlot[3] = x4plot
    
    return BackPlot, FrontPlot, MasterPlot


### Execute nTries different trajectories with initial conditions defined below:

BackPlots = [[] for i in range(nTries)]
FrontPlots = [[] for i in range(nTries)]
MasterPlots = [[] for i in range(nTries)]

## Initial Conditions
# x1s = [1]
# x2s = [1]
# x3s = [1]
# # x4s = [1]
# for idx in range(12):
    
                
    
# x1s = [2.4, 2.4, 0.5, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# x2s = [1.0, 0.0, 0.5, 1.0, 2.4, 0.0, 1.0, 2.4, 0.9, 1.0, 0.0, 0.0]
# x3s = [0.0, 1.0, 0.5, 1.0, 0.0, 2.4, 2.4, 1.0, 0.0, 0.0, 0.5, 0.0]
# x4s = [0.3, 0.3, 0.3, 0.0, 0.3, 0.3, 0.1, 0.1, 0.1, 0.0, 0.5]

x1s = [1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0]
x2s = [1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
x3s = [1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0]

backCounter0 = 0
frontCounter0 = 0
masterCounter0 = 0


plt.figure(1)
plt.figure(2)
plt.figure(3)

for nTry in range(nTries):

    x1 = np.zeros(nstep)
    x2 = np.zeros(nstep)
    x3 = np.zeros(nstep)
    # x4 = np.zeros(nstep)
    t = np.zeros(nstep)

    x1[0] = x2s[nTry]#2.4 #0.0
    x2[0] = x2s[nTry]
    x3[0] = x3s[nTry]
    # x4[0] = x2s[nTry]
    t[0] = 0.0
    
    # [ BackPlots[nTry], FrontPlots[nTry], MasterPlots[nTry] ] = gle4_rk4(x1, x2, x3, x4, t)
    [ BackPlots[nTry], FrontPlots[nTry], MasterPlots[nTry] ] = gle4_rk4(x1, x2, x3, t)

    backCounter0 += len(BackPlots[nTry][0])
    frontCounter0 += len(FrontPlots[nTry][0])
    masterCounter0 += len(MasterPlots[nTry][0])
    
    # plt.figure(1)
    # plt.plot(BackPlots[nTry][1], BackPlots[nTry][2], 'b.', FrontPlots[nTry][1], FrontPlots[nTry][2], 'r.')
    # plt.xlabel('y')
    # plt.ylabel('z')
    
    plt.figure(2)
    plt.plot(BackPlots[nTry][0], BackPlots[nTry][1], 'b.', FrontPlots[nTry][0], FrontPlots[nTry][1], 'r.')
    plt.xlabel('x')
    plt.ylabel('y')
    
    fig3 = plt.figure(3)
    # ax = plt.axes(projection='3d')
    sub3 = fig3.add_subplot(1, 1, 1, projection='3d')
    # ax.plot3D(BackPlots[nTry][0], BackPlots[nTry][1], BackPlots[nTry][2], BackPlots[nTry][3], 'b.')
    # ax.plot3D(FrontPlots[nTry][0], FrontPlots[nTry][1], FrontPlots[nTry][2], FrontPlots[nTry][3], 'r.')
    sub3.plot3D(BackPlots[nTry][0], BackPlots[nTry][1], BackPlots[nTry][2], 'b.')
    sub3.plot3D(FrontPlots[nTry][0], FrontPlots[nTry][1], FrontPlots[nTry][2], 'r.')
    sub3.set_xlabel('x')
    sub3.set_ylabel('y')
    sub3.set_zlabel('z')
    plt.text(0.5, 0.5, 'a = ' + str(a) + ', b = ' + str(b) + ', F = ' + str(F) + ', G = ' + str(G), horizontalalignment='left', verticalalignment='bottom', transform=sub3.transAxes)
    plt.show()


