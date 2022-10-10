#!/usr/bin/env python
from copy import deepcopy
from difflib import diff_bytes
import json
import numpy as np
import matplotlib.pyplot as plt
from MOVESTARziran.movestar import movestar
from MOVESTAR import movestar as ms
from runner import Diffuse
import optparse
import sumolib
from sumolib.visualization import helpers  # noqa
busStops =[[2191, 1300]]#[[1892, 1405], [1829,1317], [1815,1293], [1870,1253], [2244,1314], [2191,1300], [1895,1258], [1870,1375]] #[[2191, 1300]]#[[2191, 1300]]#[[1892, 1405], [1829,1317], [1815,1293], [1870,1253], [2244,1314], [2191,1300], [1895,1258], [1870,1375]] #[[2191, 1300]]#
bound = [(1074.49, 880.15), (2677.65, 1718.8)]
net_X = bound[1][0] - bound[0][0]
net_Y = bound[1][1] - bound[0][1]
X = net_X + 0.1*net_X#round(net_X)# * gridSize) # meter
Y = net_Y + 0.1*net_Y#round(net_Y)# * gridSize) # meter
timeGap = 1
gridSize =  round(2*timeGap*10/(2-timeGap*1.1))  # meter
def analyzeConcentration(Log_concentration, Log_concentration2):
    plt.figure()
    concAtBusStops = {}
    for conc in Log_concentration:
        for i in range(len(busStops)):
            emissionPosX = round ((busStops[i][0] - bound[0][0] + X/1.1*0.05)/ gridSize)
            emissionPosY = round ((busStops[i][1] - bound[0][1] + Y/1.1*0.05)/ gridSize)
            if i not in concAtBusStops:
                concAtBusStops[i] = []
            concAtBusStops[i].append(conc[emissionPosX, emissionPosY, 0])
    
    concAtBusStops2 = {}
    for conc in Log_concentration2:
        for i in range(len(busStops)):
            emissionPosX = round ((busStops[i][0] - bound[0][0] + X/1.1*0.05)/ gridSize)
            emissionPosY = round ((busStops[i][1] - bound[0][1] + Y/1.1*0.05)/ gridSize)
            if i not in concAtBusStops2:
                concAtBusStops2[i] = []
            concAtBusStops2[i].append(conc[emissionPosX, emissionPosY, 0])
    pltBusStop2 = []
    pltBusStop = []
    for i in concAtBusStops:
        pltBusStop.append(plt.plot(range(len(concAtBusStops[i])), concAtBusStops[i]))
        print("Bus Stop " + str(i) + ": " + str(np.sum(concAtBusStops[i])/3600))
    for i in concAtBusStops2:
        pltBusStop2.append(plt.plot(range(len(concAtBusStops2[i])), concAtBusStops2[i]))
        print("Bus Stop " + str(i) + ": " + str(np.sum(concAtBusStops2[i])/3600))
    plt.title("Concentration at Bus Station")
    plt.axis([0, 3600, 0, 1E-6])
    plt.xlabel('Time steps')
    plt.ylabel('PM2.5 Concentration (g/m^3)')
    plt.legend(["Centralized Pickup", "Distributed Pickup"])
    #plt.legend(['0', '1', '2', '3', '4', '5', '6', '7'])

def analyzeInhale(Log_pedestrain, name, pltSize, timeInterval):
    plt.figure()
    inhaleTotal = {}
    for steps in Log_pedestrain:
        for peds in steps[1]:
            # print(peds[1][peds])
            if peds not in inhaleTotal:
                inhaleTotal[peds] = []
            inhaleTotal[peds].append(deepcopy(steps[1][peds]))
        # print(inhaleTotall)
    
    averageInhale = 0
    inhaleTotallEnd = []
    for peds in inhaleTotal:
        inhaleTotal[peds] = np.array(inhaleTotal[peds])
        inhaleTotal[peds] = inhaleTotal[peds] #/ 44 / 41
        #print(peds)
        #if int(peds) <= 139 and int(peds) > 137:
        plt.plot(range(int(peds)*timeInterval, int(peds)*timeInterval + len(inhaleTotal[peds])), inhaleTotal[peds])
        inhaleTotallEnd.append(inhaleTotal[peds][-1])
        averageInhale += inhaleTotal[peds][-1] / len(inhaleTotal)
        #print(len(inhaleTotall[peds]))
    plt.title(name)
    plt.axis(pltSize)
    plt.xlabel('Time steps')
    plt.ylabel('Inhale Mass (μg)')
    print("Total Peds: ", len(inhaleTotal))
    print("Maximum Inhale Mass: ", np.max(inhaleTotallEnd))
    print("Minimum Inhale Mass: ", np.min(inhaleTotallEnd))
    print("Median Inhale Mass: ", np.median(inhaleTotallEnd))
    print("Average Inhale Mass: ", averageInhale)
    

def analyzeVehicle(Log_vehicle):
    mMOVESTAR = ms.MOVESTAR()
    c2022_speed = []
    c2022_acc = [0]
    emissionRateTotal = {}
    travelDistanceTotal = {}
    for steps in Log_vehicle:
        for vehicles in steps[1]:
            if vehicles not in emissionRateTotal:
                emissionRateTotal[vehicles] = []
            emissionRateTotal[vehicles].append(deepcopy(steps[1][vehicles][3])) # 6: CO2, 3: PM2.5, 2: NOx,
            if vehicles not in travelDistanceTotal:
                travelDistanceTotal[vehicles] = 0
            travelDistanceTotal[vehicles] += steps[0][vehicles]["64"]
            if 'c1139' in vehicles:
                c2022_speed.append(steps[0][vehicles]["64"])
               
    print(sum(c2022_speed))
    emissionRateSumTaxi = 0
    emissionRateSumPassenger = 0
    totalDistanceTaxi = 0
    totalDistancePassenger = 0
    for vehicles in emissionRateTotal:
         
        if 'c' in vehicles:
            emissionRateSumPassenger += np.sum(emissionRateTotal[vehicles])    
            totalDistancePassenger += travelDistanceTotal[vehicles]
            # if np.sum(emissionRateTotal[vehicles]) > 100:
            #     print(vehicles, ": ", np.sum(emissionRateTotal[vehicles]))
            # if 'c2202' in vehicles:
            #     print("Emis: ", np.array((emissionRateTotal[vehicles])) *3600  )
            #     print("Total emis: ", sum(np.array((emissionRateTotal[vehicles])))   )
            #     print("Total time: ", len((emissionRateTotal[vehicles])))
            #     print("Total dist: ", travelDistanceTotal[vehicles])
            #     print("Avg Speed: ", travelDistanceTotal[vehicles] / 1609.34 / (len((emissionRateTotal[vehicles])) / 3600))
            #     print("Speed: ", np.array(c2022_speed) )#4.73090779e+00 9.45306715e+00 1.28569759e+01
            #     speed = np.array(c2022_speed)# * 2.23693629
            #     #print("Acc: ", np.array([0, c2022_speed]))# - np.array([c2022_speed, 0]))
            #     er = mMOVESTAR.movestar(speed[1:4])
            #     print(er)
        else:
            emissionRateSumTaxi += np.sum(emissionRateTotal[vehicles]) 
            totalDistanceTaxi += travelDistanceTotal[vehicles]
    
    totalDistanceTaxi /= 1609.34 # meter to mile
    totalDistancePassenger /= 1609.34 # meter to mile
    print("Total Emission of Passenger Vehicles: ", emissionRateSumPassenger)
    print("Total Emission of Taxis: ", emissionRateSumTaxi)
    print("Total Distance of Taxis: ", totalDistanceTaxi + totalDistancePassenger)
    print("Total Emission: ", emissionRateSumTaxi + emissionRateSumPassenger)
    print("VMT: ", (emissionRateSumTaxi + emissionRateSumPassenger) / (totalDistanceTaxi + totalDistancePassenger))
    #print("Avg Speed: ", )
def compareConcentration(conc1, conc2):
    diffMax = 0
    j = 5
    x = round((busStops[j][0] - bound[0][0] + X/1.1*0.05)/ gridSize)
    y = round ((busStops[j][1] - bound[0][0] + X/1.1*0.05)/ gridSize)
    for i in range(np.size(conc1, 0)):
        diff = abs(conc1[i] - conc2[i])
        # diffSum = diff[x, y, 2]
        diffSum = sum(sum(diff[:, :, 2]))
        if diffSum >= diffMax:# and i < 1800:
            bestDiff = i
            diffMax = diffSum
     
     
    print(x, y)
    print(diff[x, y, 2])
    print(bestDiff)
    print(diffMax)
     
    # diff = conc1[bestDiff] 
    # plt.figure()
    # levels = np.arange(0, max(np.max(diff[:, :, 2]), 0.01), max(np.max(diff[:, :, 2]), 0.01) / 100)
    # plt.contourf(range(np.size(diff, 0)), range(np.size(diff, 1)),  diff[:, :, 2].T, levels, cmap=plt.get_cmap('nipy_spectral_r'))
    # plt.draw()
    # diff = conc2[bestDiff]
    # plt.figure()
    # #levels = np.arange(0, max(np.max(diff[:, :, 2]), 0.01), max(np.max(diff[:, :, 2]), 0.01) / 100)
    # plt.contourf(range(np.size(diff, 0)), range(np.size(diff, 1)),  diff[:, :, 2].T, levels, cmap=plt.get_cmap('nipy_spectral_r'))
    # plt.draw()


    options = get_options()
    fig, ax = helpers.openFigure(options)
    #plt.ion()
    diff = conc1[bestDiff] - conc2[bestDiff]
    
    conc3 = diff[:, :, 2]
    # print(np.size(conc3))
    # aaa = np.tile(conc3, 3)
    # print(np.size(aaa, 1))

     
    # conc3 = conc3.repeat(22,axis=1)
    # conc3 = conc3.repeat(22,axis=0)
 



    #plt.figure()
    levels = np.arange(0, max(np.max(conc3), 0.000000001), max(np.max(conc3), 0.000000001) / 100)
    cs = plt.contourf(range(np.size(conc3, 0)), range(np.size(conc3, 1)),  conc3.T, levels, cmap=plt.get_cmap('nipy_spectral_r'))
    cbar = fig.colorbar(cs)
    # net = sumolib.net.readNet(options.net)
    # helpers.plotNet(net, {}, {}, options)
    plt.draw()
    plt.pause(1000)
def get_options(args=None):
    optParser = optparse.OptionParser()
    optParser.add_option("--nogui", action="store_true",
                         default=False, help="run the commandline version of sumo")
    optParser.add_option("-n", "--net", dest="net", metavar="FILE",
                         default="ChicagoAve.net.xml", help="Defines the network to read")
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
if __name__ == "__main__":

    # with open("logs2\/singleCO2\Log_concentration.npy","rb") as file:
    #     # size of Log_concentration is [time, gridsize, gridsize, len(Z)]
    #     Log_concentration_single = np.load(file)

    with open("logs2\/11_cent_PM\Log_concentration.npy","rb") as file:
        # size of Log_concentration is [time, gridsize, gridsize, len(Z)]
        Log_concentration_cent = np.load(file)
    with open("logs2\/11_dist_PM\Log_concentration.npy","rb") as file:
        # size of Log_concentration is [time, gridsize, gridsize, len(Z)]
        Log_concentration_dist = np.load(file)
    # with open("logs\/venkyAtEachEdge\Log_pedestrain.json","r") as file:
    #     Log_pedestrain = file.read()
    #     Log_pedestrain_venkyAtEachEdge = json.loads(Log_pedestrain)
    # with open("logs\/venkyAtBusStop\Log_pedestrain.json","r") as file:
    #     Log_pedestrain = file.read()
    #     Log_pedestrain_venkyAtBusStop = json.loads(Log_pedestrain)



    # with open("logs2\/11_cent_PM\Log_pedestrain.json","r") as file:
    #     Log_pedestrain = file.read()
    #     Log_pedestrain_atBusStop = json.loads(Log_pedestrain)
    # with open("logs2\/11_dist_PM\Log_pedestrain.json","r") as file:
    #     Log_pedestrain = file.read()
    #     Log_pedestrain_atEachEdge = json.loads(Log_pedestrain)

    # with open("logs2\/11_cent_NOx\Log_pedestrain.json","r") as file:
    #     Log_pedestrain = file.read()
    #     Log_pedestrain_atBusStop = json.loads(Log_pedestrain)
    # with open("logs2\/11_dist_NOx\Log_pedestrain.json","r") as file:
    #     Log_pedestrain = file.read()
    #     Log_pedestrain_atEachEdge = json.loads(Log_pedestrain)



    # with open("logs2\/1_dist_NOx/Log_vehicle.json","r") as file:
    #     Log_vehicle = file.read()
    #     Log_vehicle_atEachEdge = json.loads(Log_vehicle)
    # with open("logs2\/11_dist_NOx/Log_vehicle.json","r") as file:
    #     Log_vehicle = file.read()
    #     Log_vehicle_atBusStop = json.loads(Log_vehicle)

    # compareConcentration(Log_concentration_cent, Log_concentration_dist)
    #analyzeConcentration(Log_concentration_single) 
    analyzeConcentration(Log_concentration_cent, Log_concentration_dist)
    # analyzeConcentration(Log_concentration_dist)
    # analyzeInhale(Log_pedestrain_atBusStop, "Scenario 1: Centralized Pickup", [0, 3600, 0, 0.045], 15) #  PM: 0.045, NOx: 5
    # analyzeInhale(Log_pedestrain_atEachEdge, "Scenario 2: Distributed Pickup", [0, 3600, 0, 0.045], 15)
    # analyzeInhale(Log_pedestrain_atBusStop, "Scenario 1: Centralized Dispatching", [0, 3600, 0, 5], 15)
    # analyzeInhale(Log_pedestrain_atEachEdge, "Scenario 2: Distributed Dispatching", [0, 3600, 0, 5], 15)
    # analyzeVehicle(Log_vehicle_atBusStop)
    # analyzeVehicle(Log_vehicle_atEachEdge)
     
    plt.show()