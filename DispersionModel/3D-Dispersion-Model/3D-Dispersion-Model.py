#!/usr/bin/env python
# Copyright (C) 2009-2022 Xuanpeng Zhao.

# @file    runner.py
# @author  Xuanpeng Zhao
# @date    2022-06-30

from __future__ import absolute_import
from __future__ import print_function
from audioop import mul
import enum
from sqlite3 import Timestamp
from tkinter import E, ROUND

from pkg_resources import ContextualVersionConflict
from MOVESTAR import movestar as ms
import os
import math
import sys
import optparse
import json
import random
import time
import queue
from optparse import OptionParser
import pandas as pd
# we need to import python modules from the $SUMO_HOME/tools directory
if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:
    sys.exit("please declare environment variable 'SUMO_HOME'")
import sumolib
from sumolib import checkBinary  # noqa
import traci # noqa
import matplotlib.pyplot as plt  # noqa
from sumolib.visualization import helpers  # noqa
import numpy as np
from copy import deepcopy


breathingRate = 17/24/3600
def predefinedPedestrainPara(pedType):
     
    if pedType == "adult":
        return [1.2, breathingRate]
    if pedType == "child":
        return [0.6, breathingRate]

def fading2D(conc, fading, threshold):
    result = conc * (1 - fading)
    result[result < threshold] = 0
    return result

def expanding2D(conc, gridSize, expanding, timeStep, windX, windY):
    L_e = (timeStep * expanding) * gridSize + gridSize
    L_e_y_plus = ((timeStep * expanding) / 2) * gridSize + windY * timeStep
    L_e_x_plus = ((timeStep * expanding) / 2) * gridSize + windX * timeStep
    L_e_y_minus = ((timeStep * expanding) / 2) * gridSize - windY * timeStep
    L_e_x_minus = ((timeStep * expanding) / 2) * gridSize - windX * timeStep
 
    if (gridSize - ((timeStep * expanding) / 2) * gridSize < abs(windY * timeStep)) or \
       (gridSize - ((timeStep * expanding) / 2) * gridSize < (abs(windX * timeStep))):
        print("Warning: Time Step OR Wind Speed Too Large")
     
    topLeft = max(min(L_e_y_plus, gridSize), 0) * max(min(L_e_x_minus, gridSize), 0)
    top = max(min(L_e_y_plus, gridSize), 0) * max(min(L_e_x_minus + gridSize, gridSize), 0)
    topRight =  max(min(L_e_y_plus, gridSize), 0) * max(min(L_e_x_plus, gridSize), 0)
    midLeft = max(min(L_e_y_minus + gridSize, gridSize), 0) * max(min(L_e_x_minus, gridSize), 0)
    midRight = max(min(L_e_y_minus + gridSize, gridSize), 0) * max(min(L_e_x_plus, gridSize), 0)
    botLeft =  max(min(L_e_y_minus, gridSize), 0) * max(min(L_e_x_minus, gridSize), 0)
    bot = max(min(L_e_y_minus, gridSize), 0) * max(min(L_e_x_minus + gridSize, gridSize), 0)
    botRight = max(min(L_e_y_minus, gridSize), 0) * max(min(L_e_x_plus, gridSize), 0)
    mid = 0
    result = np.array([topLeft, top, topRight, midLeft, mid, midRight, botLeft, bot, botRight])
    result = result / (L_e ** 2)
    result[4] = 1 - sum(result)
    
    result = result * conc
    return result

def dispersion(timeGap, concLast, X, Y, Z, gridSize, windEnvX, windEnvY, vehicleInfo, emissionRate, expanding, fading, decay, fadingThreshold, bound, globalFlag):

 
    # Initial
    concCurr = np.zeros((int (X / gridSize) + 1, int (Y / gridSize) + 1, len(Z)))
    # Expanding
    for iz in range(len(Z)):
        [row,col]  = np.nonzero(concLast[:, :, iz])

        if row.size > 0:
            for j in range(row.size):      
       
                neighbors = expanding2D(concLast[row[j], col[j], iz], gridSize, expanding, timeGap, windEnvX[row[j], col[j], iz], windEnvY[row[j], col[j], iz])
              
                concCurr[row[j], col[j], iz] += neighbors[4] 
                
                if (row[j] - 1 >= 0):
                    concCurr[row[j] - 1, col[j], iz] += neighbors[1]
                    if (col[j] + 1 < int (Y / gridSize)):
                        concCurr[row[j] - 1, col[j] + 1, iz] += neighbors[2]
            
                    if (col[j] - 1 >= 0):
                        concCurr[row[j] - 1, col[j] - 1, iz] += neighbors[0]
            
                
                if (row[j] + 1 < int (X / gridSize)):
                    concCurr[row[j] + 1, col[j], iz] += neighbors[7]
                    if (col[j] + 1 < int (Y / gridSize)):
                        concCurr[row[j] + 1, col[j] + 1, iz] += neighbors[8]
                    
                    if (col[j] - 1 >= 0):
                        concCurr[row[j] + 1, col[j] - 1, iz] += neighbors[6]

                if (col[j] + 1 < int (Y / gridSize)):
                    concCurr[row[j], col[j] + 1, iz] += neighbors[5]
                
                if (col[j] - 1 > 0):
                    concCurr[row[j], col[j] - 1, iz] += neighbors[3] 
         
    # New source
    for veh_id in vehicleInfo:
        pos = vehicleInfo[veh_id][traci.constants.VAR_POSITION]     
        er = emissionRate[veh_id]
        #print(int ((pos[0] - bound[0][0] + X/1.1*0.05)/ gridSize))
        #print(int ((pos[1] - bound[0][1] + Y/1.1*0.05)/ gridSize))
        emissionPosX = round ((pos[0] - bound[0][0] + X/1.1*0.05)/ gridSize)
        emissionPosY = round ((pos[1] - bound[0][1] + Y/1.1*0.05)/ gridSize)
        if emissionPosX >= 0 and emissionPosX < concCurr.shape[0] and emissionPosY >= 0 and emissionPosY < concCurr.shape[1]  :
            globalFlag += 1
            for iz in range(len(Z)):
                # CO(g),HC(g),NOx(g),PM2.5(g),PM10(g),Energy(KJ),CO2(g)
                #print(er[6])
                #print(er)
                concCurr[emissionPosX, emissionPosY, iz] += er[2] * math.exp((-(Z[iz] - 0.3)**2 / (2 * decay**2)) ) * timeGap / (gridSize**2)
        #else:
        #    print("pos: " + str(pos))
 
         
    # Concentration Fading
    concCurr = fading2D(concCurr, fading, fadingThreshold)

    concLast = concCurr
    
    return concLast, globalFlag
 
def dispersion2(timeGap, concLast, X, Y, Z, gridSize, vehicleInfo, emissionRate, bound, windEnvX, windEnvY, Kx, Ky, Kz):
   
    # Kx = gridSize * 2
    # Ky = gridSize * 2
    # sigw = 10
    # Kz = Z * sigw
    #concLast[40, 21, 0] += 0.0001
    for veh_id in vehicleInfo:
        pos = vehicleInfo[veh_id][traci.constants.VAR_POSITION]     
        er = emissionRate[veh_id]
        #print(int ((pos[0] - bound[0][0] + X/1.1*0.05)/ gridSize))
        #print(int ((pos[1] - bound[0][1] + Y/1.1*0.05)/ gridSize))
        posX = (pos[0] - bound[0][0] + X/1.1*0.05)
        posY = (pos[1] - bound[0][1] + Y/1.1*0.05)
        emissionPosX = round (posX / gridSize)
        emissionPosY = round (posY / gridSize)
        
        if emissionPosX >= 0 and emissionPosX < concLast.shape[0] and emissionPosY >= 0 and emissionPosY < concLast.shape[1]:
            #print(er[6] / 3600 * timeGap / (gridSize**2*Z[-1]))
            ix = int ()
            w1 = 1 / ((posX - emissionPosX * gridSize)**2+(posY-emissionPosY * gridSize)**2)
            w2 = 1 / ((posX - (emissionPosX + 1) * gridSize)**2+(posY-emissionPosY * gridSize)**2)
            w3 = 1 / ((posX - (emissionPosX + 1) * gridSize)**2+(posY - (emissionPosY + 1) * gridSize)**2)
            w4 = 1 / ((posX - emissionPosX * gridSize)**2+(posY-(emissionPosY + 1) * gridSize)**2)
            # w1~w4
            W=w1+w2+w3+w4
            rate = er[2]
            concLast[emissionPosX, emissionPosY, 0] += rate * w1/W * timeGap / (gridSize**2*(Z[0])) # change concCurr to concLast
            concLast[emissionPosX+1, emissionPosY, 0] += rate * w2/W * timeGap / (gridSize**2*(Z[0])) # change concCurr to concLast
            concLast[emissionPosX+1, emissionPosY+1, 0] += rate * w3/W * timeGap / (gridSize**2*(Z[0])) # change concCurr to concLast
            concLast[emissionPosX, emissionPosY+1   , 0] += rate * w4/W * timeGap / (gridSize**2*(Z[0])) # change concCurr to concLast
        else:
            traci.vehicle.remove(veh_id, 2)
            print("pos: " + str(pos))
    ##### TEST ONLY  ##############
    #emissionPosX = round ((1900 - bound[0][0] + X/1.1*0.05)/ gridSize)
    #emissionPosY = round ((1300 - bound[0][1] + Y/1.1*0.05)/ gridSize)
    #concLast[emissionPosX, emissionPosY, 0] += 0.00248 / 3600 * timeGap  / (gridSize**2*(Z[1]- Z[0]))
    ###########################
    for iz in range(len(Z)):
        concLast[:, :, iz] = Advect(windEnvX[:, :, iz], windEnvY[:, :, iz], timeGap, concLast[:, :, iz], X, Y, gridSize)
    for ix in range(int (X / gridSize) + 1):
        for iy in range(int (Y / gridSize) + 1):
            
            concLast[ix, iy, :] = Diffuse_Z(concLast[ix, iy, :], Z, Kz, timeGap)
        
    for iz in range(len(Z)):
        for iy in range(int (Y / gridSize) + 1):
            concLast[:, iy, iz] = Diffuse(concLast[:, iy, iz], X, gridSize, Kx, timeGap)
        for ix in range(int (X / gridSize) + 1):
            concLast[ix, :, iz] = Diffuse(concLast[ix, :, iz], Y, gridSize, Ky, timeGap)
    
    
    return concLast
def Diffuse_Z(V, Z, K, timeGap):
    e = np.zeros(len(Z))
    g = np.zeros(len(Z))
    h = np.zeros(len(Z))
    f = np.zeros(len(Z))
    for iz in range(1, len(Z) - 1):
        dzp = Z[iz + 1] - Z[iz]
        dzm = Z[iz] - Z[iz - 1]
        dz = (dzp + dzm) / 2
        Kzp = (K[iz + 1] + K[iz]) / 2
        Kzm = (K[iz] + K[iz - 1]) / 2

        alphap = Kzp * timeGap / (dzp * dzm)
        alpham = Kzm * timeGap / (dz * dzm)
        e[iz] = -alpham
        g[iz] = -alphap
        f[iz] = 1 + alpham + alphap
        h[iz] = V[iz]
    if len(Z) >= 3:
        alpha1 = 2 * K[1] * timeGap / ((Z[1] - Z[0])*(Z[2] - Z[0]))
        alpha2 = 2 * K[-2] * timeGap / ((Z[-1] - Z[-2])*(Z[-1] - Z[-3]))
        e[0] = 0.0
        f[0] = 1.0 + alpha1
        g[0] = - alpha1
        h[0] = V[0]
        e[-1] = - alpha2
        f[-1] = 1 + alpha2
        g[-1] = 0.0
        h[-1] = V[-1]
    else:
        e[0] = 0.0
        f[0] = 1.0
        g[0] = 0
        h[0] = V[0]
        e[-1] = 1.0
        f[-1] = -1.0
        g[-1] = 0.0
        h[-1] = 0.0
    Vn = Tridiag_Solver_ADM(e, f, g, h)  
     
     
    return Vn
def Diffuse(V, X, gridSize, K, timeGap):
    e = np.zeros((int (X / gridSize) + 1))
    g = np.zeros((int (X / gridSize) + 1))
    h = np.zeros((int (X / gridSize) + 1))
    f = np.zeros((int (X / gridSize) + 1))
    beta = K * timeGap
    for ix in range(1, int (X / gridSize)):
        dzp = gridSize
        dzm = gridSize
        dz = (dzp + dzm) / 2
        alphap = beta / (dzp * dzm)
        alpham = beta / (dz * dzm)
        e[ix] = -alpham
        g[ix] = -alphap
        f[ix] = 1 + alpham + alphap
        h[ix] = V[ix]
    alpha1 = 2 * K * timeGap / (gridSize*gridSize*2)
    alpha2 = 2 * K * timeGap / (gridSize*gridSize*2)
    e[0] = 0.0
    f[0] = 1.0 + alpha1
    g[0] = - alpha1
    h[0] = V[0]
    e[-1] = - alpha2
    f[-1] = 1 + alpha2
    g[-1] = 0.0
    h[-1] = V[-1]

    Vn = Tridiag_Solver_ADM(e, f, g, h)
     
    return Vn
    
def Businger_ADM(z, Met):
    kappa = 0.35
    zi=Met["zi"]
    turbz=Met["turbz"]
    L=Met["Lmon"]
    ustar=Met["ustar"]
    z = min(z, 0.1*zi)
    if z <= turbz:
        z = turbz
    psi = z / L
    if psi > 0:
        phih = 0.74 * (1 + 6.3 * psi)
    else:
        phih = 0.74 / (1 - 9 * psi)**0.5
    if z <= zi:
        Kz = kappa * ustar * z / phih
    else:
        Kz = 0.01
    return Kz
def Compute_Ustar_ADM(Met):
    kappa=0.4
    
    z0=Met["z0"]
    zi=Met["zi"]
    L=Met["Lmon"]
    dh=5*z0; 
    
    Wspd=Met["wspd"]
    
    z=Met["zref"]
    if L>0.0:
       psi1=-17*(1-math.exp(-0.29*(z-dh)/L))
       psi2=-17*(1-math.exp(-0.29*z0/L))
    else:
       x1=(1-16*(z-dh)/L)**0.25; x2=(1-16*z0/L)**0.25
       psi1=2*math.log((1+x1)/2)+math.log((1+x1*x1)/2)-2*math.atan(x1)+math.pi/2
       psi2=2*math.log((1+x2)/2)+math.log((1+x2*x2)/2)-2*math.atan(x2)+math.pi/2
    Met["ustar"]=kappa*Wspd/(math.log((z-dh)/z0)-psi1+psi2)
    return Met
def Similarity_Wind_Airport_ADM(z, Met):
    kappa=0.4
    z0=Met["z0"]
    zi=Met["zi"]
    Lmon=Met["Lmon"]
    dh=5*z0
    ust=Met["ustar"] 
    z=min(z,0.1*zi)
    z=max(z,1.1*(z0+dh))
    if Lmon>0.0: 
       psi1=-17*(1-math.exp(-0.29*(z-dh)/Lmon))
       psi2=-17*(1-math.exp(-0.29*z0/Lmon))
    else:
       x1=(1-16*(z-dh)/Lmon)**0.25; x2=(1-16*z0/Lmon)**0.25
       psi1=2*math.log((1+x1)/2)+math.log((1+x1*x1)/2)-2*math.atan(x1)+math.pi/2
       psi2=2*math.log((1+x2)/2)+math.log((1+x2*x2)/2)-2*math.atan(x2)+math.pi/2
    windSpeed=ust*(math.log((z-dh)/z0)-psi1+psi2)/kappa
    return windSpeed
def Initial(Z, Met):
    Met=Compute_Ustar_ADM(Met)
    windSpeed = np.zeros(len(Z))
    wd=(Met["wdir"])
    u =  np.zeros(len(Z))
    v =  np.zeros(len(Z))
    Kz = np.zeros(len(Z))
    for iz in range(len(Z)):
        Kz[iz] = Businger_ADM(Z[iz], Met)
        windSpeed[iz] = Similarity_Wind_Airport_ADM(Z[iz], Met)
         
        u[iz] = windSpeed[iz]*np.cos(wd)
        v[iz] = windSpeed[iz]*np.sin(wd)
    print(windSpeed)
    return Kz, u, v
def Tridiag_Solver_ADM(e, f, g, h):
    x = np.zeros(len(f))
    for k in range(1, len(f)):
        mult = e[k] / f[k-1]
        f[k] -= mult * g[k-1]
        h[k] -= mult * h[k-1]
    x[-1] = h[-1] / f[-1]
    for k in range(len(f) - 2, -1, -1):
        x[k] = (h[k] - g[k] * x[k+1]) / f[k]
    
     
    return x
def Advect(windX, windY, timeStep, concLast, X, Y, gridSize):
    concCurr = np.zeros((np.size(concLast, 0), np.size(concLast, 1)))
    for ix in range(int (X / gridSize) + 1):
        for iy in range(int (Y / gridSize) + 1):
             
            Xn = ix * gridSize - windX[ix, iy] * timeStep
            if round (Xn / gridSize) > int (X / gridSize) + 1 or round (Xn / gridSize) < 0:
                continue

            Yn = iy * gridSize - windY[ix, iy] * timeStep
            if round (Yn / gridSize) > int (Y / gridSize) + 1 or round (Yn / gridSize) < 0:
                continue
            i = round (Xn / gridSize)#min(max(0, int (Xn / gridSize) ), int (X / gridSize) + 1)
            j = round (Yn / gridSize)#min(max(0, int (Yn / gridSize) ), int (Y / gridSize) + 1)

           

            p = (Xn - i * gridSize) / (gridSize)
            q = (Yn - j * gridSize) / (gridSize)
            
            imax=int (X / gridSize)
            jmax=int (Y / gridSize)
            ip1=min(imax,i+1)
            jp1=min(jmax,j+1)
            if ip1 == imax: 
                concCurr[ix, iy] = 0
            elif jp1 == jmax: 
                concCurr[ix, iy] = 0
            else:
                Vb = p * concLast[ip1,j]+(1-p)*concLast[i,j]
                Vt = p * concLast[ip1,jp1]+(1-p)*concLast[i,jp1]
                concCurr[ix, iy] = q*Vt+(1-q)*Vb

            """
            p = max(0, Xn - ix * gridSize) / (gridSize)
            q = max(0, Yn - iy * gridSize) / (gridSize)
      
            V1 = concLast[i, j]
            V2 = concLast[i+1, j]
            V3 = concLast[i, j+1]
            V4 = concLast[i+1, j+1]

            concLast[i, j] = (V1 * (1-p) + V2 * p) * (1 - q) + (V4 * (1-p) + V3 * p) * q
            """
    return concCurr

 
    

def run(helpers, fig, bound, dispersionModel):

   

    np.set_printoptions(threshold=sys.maxsize)
    # helpers might get the distance between origin to the most left bottom corner 
    
    """execute the TraCI control loop"""
    # GUI settings
    traci.gui.setSchema('View #0', "real world")
   
    globalFlag = 0
    # initialize
    timeGap =  1 # s
    #print(round(2*timeGap*10/(2-timeGap*1.1)))
    gridSize = 5#round(2*timeGap*10/(2-timeGap*1.1))  # meter
    print("gridSize ", gridSize)
    #timeGap = 1
    net_X = bound[1][0] - bound[0][0]
    net_Y = bound[1][1] - bound[0][1]
    X = net_X + 0.1*net_X#round(net_X)# * gridSize) # meter
    Y = net_Y + 0.1*net_Y#round(net_Y)# * gridSize) # meter
    #Z = np.array([0.3, 0.6, 1.2, 2.4, 4.8, 9.6])
    Z = np.array([0.3, 1.2, 4.8, 19.2, 76.8, 200]) 
    #, 2.4, 4.8, 9.6, 19.2, 38.4, 76.8
    #Z = np.array([2, 3, 4.5, 6.75, 10.125])#[1, 1.5, 2.25, 3.375]


    offsetX = 0
    offsetY = 0
        
    # for i in range(gridSize):

    #     if len(range(round(bound[0][1] - 0.05*net_Y), round(bound[1][1] - i + 0.05*net_Y), gridSize)) == round (Y / gridSize):
    #         offsetY = i
    #         break
    # for j in range(gridSize):
    
    #     if len(range(round(bound[0][0] - 0.05*net_X), round(bound[1][0] - j + 0.05*net_X), gridSize)) == round (X / gridSize):
    #         offsetX = j
                
    #         break

    x, y = np.meshgrid(range(round(bound[0][0] - 0.05*net_X), round(bound[1][0] - offsetX + 0.05*net_X), gridSize), range(round(bound[0][1] - 0.05*net_Y), round(bound[1][1] - offsetY + 0.05*net_Y), gridSize))
    
    print("np.size(x, 0): ", np.size(x, 0))
    print("int (X / gridSize) + 1: ", int (X / gridSize) + 1)
    
    # Inhale
        

    
    # function
    expanding = 0.9
    fading = 0.055
    decay = 0.48
    fadingThreshold = 1.0e-15

    concCurr = np.zeros((int (X / gridSize) + 1, int (Y / gridSize) + 1, len(Z)))
    concLast = concCurr
    
    
    # Wind
    windSpeed = 10#2   # m/s
    windDirec = -3/4 * np.pi  #-(WINDDIREC(meas) / 360 * 2 * pi + pi/2) - pi/2 #
    # Wind X:1, Y:2


    windEnvX = np.ones((np.size(concCurr, 0), np.size(concCurr, 1), len(Z))) * windSpeed * np.cos(windDirec)
    windEnvY = np.ones((np.size(concCurr, 0), np.size(concCurr, 1), len(Z))) * windSpeed * np.sin(windDirec)
    
    Met = {}
    Met["wspd"] = windSpeed
    Met["wdir"] = -windDirec
    #Met["ustar"] = 0.3
    Met["Lmon"] = -200
    Met["zi"] = 200
    Met["z0"] = 0.1
    Met["zref"] = 5
    Met["turbz"] = 0.5
    if dispersionModel == "3D-grid":
        Kz, u, v = Initial(Z, Met)
        print("Kz: ", Kz)
        
        windEnvX = np.ones((np.size(concCurr, 0), np.size(concCurr, 1), len(Z))) * np.array(v)
        windEnvY = np.ones((np.size(concCurr, 0), np.size(concCurr, 1), len(Z))) * np.array(u)
    vehicleSpdAccQueue = {}
    failtoDispatched = []
    emissionRate = {}
    existingTime = {}
    mMOVESTAR = ms.MOVESTAR()
    inhaleMass = {} 

    Log_pedestrain = []
    Log_vehicle = []
    Log_concentration = []
    while traci.simulation.getMinExpectedNumber() > 0:
        
        try:
            traci.simulationStep()
        except Exception:
            #print(Log_concentration)
            #print(Log_pedestrain)
            #print(Log_vehicle)
            
            print("Error connection closed", Exception)
        
         
        for veh_id in traci.simulation.getDepartedIDList():
            traci.vehicle.subscribe(veh_id, [traci.constants.VAR_POSITION,
                                             traci.constants.VAR_SPEED, 
                                             traci.constants.VAR_ACCELERATION])
        
            
        vehicleInfo = traci.vehicle.getAllSubscriptionResults()
        #t1 = time.time()
        for veh_id in vehicleInfo:
        
            if veh_id not in existingTime: 
                existingTime[veh_id] = 0
            else:
                existingTime[veh_id] += timeGap
            if existingTime[veh_id] % (1/timeGap) == 0:
               
                if veh_id not in vehicleSpdAccQueue:
                    vehicleSpdAccQueue[veh_id] = {}
                    vehicleSpdAccQueue[veh_id]['speed'] = queue.Queue(3)
                    vehicleSpdAccQueue[veh_id]['accel'] = queue.Queue(3)
                currSpeed = vehicleInfo[veh_id][traci.constants.VAR_SPEED]
                if not vehicleSpdAccQueue[veh_id]['speed'].full():
                    vehicleSpdAccQueue[veh_id]['speed'].put(currSpeed)
                else:
                    vehicleSpdAccQueue[veh_id]["speed"].get()
                    vehicleSpdAccQueue[veh_id]['speed'].put(currSpeed)
                q = vehicleSpdAccQueue[veh_id]['speed'].qsize()
                if q > 1:
                    currAcc = vehicleSpdAccQueue[veh_id]['speed'].queue[q - 1] - vehicleSpdAccQueue[veh_id]['speed'].queue[q - 2]
                else:
                    currAcc = vehicleInfo[veh_id][traci.constants.VAR_ACCELERATION]
                if not vehicleSpdAccQueue[veh_id]['accel'].full():
                    vehicleSpdAccQueue[veh_id]['accel'].put(currAcc)
                else:
                    vehicleSpdAccQueue[veh_id]['accel'].get()
                    vehicleSpdAccQueue[veh_id]['accel'].put(currAcc)
                er = mMOVESTAR.movestar(np.array(vehicleSpdAccQueue[veh_id]['speed'].queue), Acc=np.array(vehicleSpdAccQueue[veh_id]['accel'].queue))#[100, 100, 100, 100, 100, 100, 100, 100,100]#[10, 10, 10, 10, 10, 10, 10, 10,10]#   
                 
                emissionRate[veh_id] = list(np.array(er) / 3600)
                # print(er)
        #t2 = time.time()-t1
        #print("time passed: " + str(t2))
        for veh_id in traci.simulation.getArrivedIDList():
            if veh_id in vehicleSpdAccQueue:
                vehicleSpdAccQueue.pop(veh_id)
                existingTime.pop(veh_id)
            if veh_id in emissionRate:
                emissionRate.pop(veh_id)
        timeStart = time.time()
        if dispersionModel == "2D-grid":
            concCurr, globalFlag = dispersion(timeGap, concLast, X, Y, Z, gridSize, windEnvX, windEnvY, vehicleInfo, emissionRate, expanding, fading, decay, fadingThreshold, bound, globalFlag)
        if dispersionModel == "3D-grid":
            concCurr = dispersion2(timeGap, concLast, X, Y, Z, gridSize, vehicleInfo, emissionRate, bound, windEnvX, windEnvY, 44, 44, Kz)
        concLast = concCurr
        deltaTime = timeStart - time.time()
        #print(deltaTime)
        # region Health intale calculation
     
         
        for ped_id in traci.simulation.getDepartedPersonIDList():
            traci.person.subscribe(ped_id, [traci.constants.VAR_POSITION,
                                            traci.constants.VAR_SPEED])
        personInfo = traci.person.getAllSubscriptionResults()
        for ped_id in personInfo:
            #print(traci.person.getTypeID(ped_id)) #getTypeID
            [height, breathingRate] = predefinedPedestrainPara(traci.person.getTypeID(ped_id))
            if traci.person.getStage(ped_id, 0).description != 'driving':
                # Computes inhale
                # print(traci.person.getStage(ped_id, 0))
                # np.where( (height / Z <= 1) == True)[0][0] is the index in Z axis
                pos = traci.person.getPosition(ped_id)
                inhaleCurr = concCurr[int ((pos[0] - bound[0][0] + X/1.1*0.05)/ gridSize), int ((pos[1] - bound[0][1] + Y/1.1*0.05)/ gridSize),  np.where( (height / Z <= 1) == True)[0][0]]
                if ped_id not in inhaleMass:
                    inhaleMass[ped_id] = 0
             
                inhaleMass[ped_id] += inhaleCurr * breathingRate * 1.0e6
                #print(inhaleMass[ped_id])
              
        for ped_id in traci.simulation.getArrivedPersonIDList():
            if ped_id in inhaleMass:
                inhaleMass.pop(ped_id)
 
    
        #endregion

        # region Curbside pick up
        reservations = traci.person.getTaxiReservations(1)
        reservation_ids = [r.id for r in reservations]
         
        #print(reservation_ids)
        try:
            fleet = traci.vehicle.getTaxiFleet(0)
            fleet = list(fleet)
        except (traci.exceptions.FatalTraCIError):
            print("No unoccupied taxi-fleet!")
        
        for i, val in enumerate(reservation_ids):
            succussDispatched = False
            print("fleet: ", fleet)
            print(reservation_ids)
            print(reservations)
            for j, val in enumerate(fleet):
                try:
                    traci.vehicle.dispatchTaxi(fleet[j], [reservation_ids[i]])
                    succussDispatched = True
                    print("Depatched taxi" + fleet[j] + " for " + reservation_ids[i])
                    fleet.pop(j)
                    break
                except:
                    print("Fail to depatched taxi")
                     
            if succussDispatched == False and reservation_ids[i] not in failtoDispatched:
                failtoDispatched.append(reservation_ids[i])
        failtoDispatchedCopy = failtoDispatched.copy()
        if failtoDispatchedCopy:
            for i, val in enumerate(failtoDispatchedCopy):
                for j, val in enumerate(fleet):
                    try:
                        traci.vehicle.dispatchTaxi(fleet[j], [failtoDispatchedCopy[i]])
                        print("Depatched taxi" + fleet[j] + " for " + failtoDispatchedCopy[i])
                        failtoDispatched.pop(i)
                        fleet.pop(j)
                        break
                    except:
                        print("Fail to depatched taxi")
        #endregion
        # dict append need "deepcopy()"
        
        Log_concentration.append(deepcopy(concCurr))
        Log_pedestrain.append(deepcopy([personInfo, inhaleMass]))
        Log_vehicle.append(deepcopy([vehicleInfo, emissionRate]))
        
         
        if traci.simulation.getTime() % 1200 == 0:
            with open("logs5\Log_concentration.npy","wb") as file:
                np.save(file, Log_concentration)
            with open("logs5\Log_pedestrain.json","w") as file:
                json.dump(Log_pedestrain, file)
            with open("logs5\Log_vehicle.json","w") as file:
                json.dump(Log_vehicle, file)

        #print(np.max(concCurr[:, :, 1]))
        #"""
        plt.ion()
        # print(np.max(concCurr[:, :, 0]))
        
        levels = np.arange(0, max(np.max(concCurr[:, :, 1]), 1e-20), max(np.max(concCurr[:, :, 1]), 1e-20) / 100)
        cs = plt.contourf(x, y, concCurr[:, :, 1].T, levels, cmap=plt.get_cmap('nipy_spectral_r'))
        cbar = fig.colorbar(cs)    
        helpers.plotNet(net, {}, {}, options)
        plt.draw()
        plt.pause(0.00001)
        fig.clear()
        #"""

        if traci.simulation.getTime() == traci.simulation.getEndTime():
            break

    traci.close()
    sys.stdout.flush()
   


def get_options(args=None):
    optParser = optparse.OptionParser()
    optParser.add_option("--nogui", action="store_true",
                         default=False, help="run the commandline version of sumo")
    optParser.add_option("-n", "--net", dest="net", metavar="FILE",
                         default="grid.net.xml", help="Defines the network to read")
    optParser.add_option("-c", "--color", dest="color",
                         default='r', help="Defines the dot color")            
    optParser.add_option("-w", "--width", dest="width",
                         type="float", default=20, help="Defines the width of the dots")
    optParser.add_option("--edge-width", dest="defaultWidth",
                         type="float", default=1, help="Defines the edge width")
    optParser.add_option("--edge-color", dest="defaultColor",
                         default='k', help="Defines the edge color")
    helpers.addInteractionOptions(optParser)
    helpers.addPlotOptions(optParser)
    options, remaining_args = optParser.parse_args(args=args)
    

    return options


# this is the main entry point of this script
if __name__ == "__main__":
    options = get_options()
    fig, ax = helpers.openFigure(options)
    ax.set_aspect("equal", None, 'C')
    net = sumolib.net.readNet(options.net)
    print(net.getBBoxXY())
    helpers.plotNet(net, {}, {}, options)
    
    #plt.draw()
    #plt.pause(0.001)
    
    
    # this script has been called from the command line. It will start sumo as a
    # server, then connect and run
    if options.nogui:
        sumoBinary = checkBinary('sumo')
    else:
        sumoBinary = checkBinary('sumo-gui')
    dispersionModel = "3D-grid"
    # this is the normal way of using traci. sumo is started as a
    # subprocess and then the python script connects and runs
    traci.start([sumoBinary, "-c", "dispersion.sumocfg", "-n", options.net,
                             "--tripinfo-output", "tripinfo.xml"])
    run(helpers, fig, net.getBBoxXY(), dispersionModel)
    #fig, ax = plt.subplots()#figure()
    #dispersion(helpers, fig)
    
