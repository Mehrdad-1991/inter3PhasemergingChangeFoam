#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
# MIT License
# 
# Copyright (c) 2026 - Universität Rostock 
#                       Fakultät für Maschinenbau und Schiffstechnik
#                       Lehrstuhl für Modellierung und Simulation
#                       
# Author: Mehrdad Kazemi
# Date: September 2021 (validated)
# Address: Albert-Einstein-Str. 2, 18059 Rostock, Deutschland
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
################################################################################
from pathlib import Path
from typing import Counter
import numpy as np
from math import sqrt
import csv
import os
from sys import path
from pathlib import Path
import glob
import shutil
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import size
from numpy.core.numeric import count_nonzero
from numpy.lib.function_base import append
from numpy.lib.shape_base import apply_along_axis
##################################################
fig = plt.figure(figsize=(14, 10))
variablesp = ['time','p']
variablesu = ['time','ux' , 'uy' , 'uz']
#################################################
averagingtimestart = 15*0.006698565
averagingtimeend = 55*0.006698565
################################################
merging_items = ['U', 'p'] 

base_dir = Path("..")
for ii in range(1,12):
    pressure_file = base_dir / "postProcessing" / f"tap{ii}" / f"tap{ii}pressure.dat"
    if pressure_file.exists(): 
        os.remove(str(pressure_file))
    velocity_file = base_dir / "postProcessing" / f"tap{ii}" / f"tap{ii}velosity.dat"
    if velocity_file.exists():             
        os.remove(str(velocity_file))
    resultpressure = open(str(pressure_file),'w+')
    resultpressure.write("%s\n" %variablesp)
    resultvelocity = open(str(velocity_file),'w+') 
    resultvelocity.write("%s\n" %variablesu)
    #########clearing lists############################
    ftpcurrent = []
    utcurrent = []
    velosity = [[]for i in range(3)]
    pressure = []
    averagevelocity = [0]*4
    averagepressure = 0.0
    ftpcurrentlast = 0.0
    utcurrentlast = 0.0
    ######sorting times #############################
    path = str(base_dir / "postProcessing" / f"tap{ii}")
    dirname = [f for f in os.listdir(path) if not f.endswith('.dat')]
    merging_times = []*len(dirname)
    for lines in dirname: 
        if float(lines) % 1.0 == 0 : merging_times.append(int(lines))
        else: merging_times.append(float(lines))
    merging_times.sort()   
    ##################################################
    for times in merging_times:            
        for items in merging_items:
            #if 'Cavitaiton' in simulation or 'cavitaiton' in simulation: 
            #    if items =='p': items = 'taprgh'
            if os.path.isfile(str(base_dir / "postProcessing" / f"tap{ii}" / str(times) / items)):
                sourcedata = open(str(base_dir / "postProcessing" / f"tap{ii}" / str(times) / items), 'r')
                for line in sourcedata:
                    if line[0] != '#':
                        line = line.replace('\n', ' ')
                        line = line.replace('(', ' ')
                        line = line.replace(')', ' ')
                        line = line.replace('[', ' ')
                        line = line.replace(']', ' ')                              
                        line = line.replace(',', ' ')
                        line = [float(s) for s in line.split()]                     
                        if items == 'p' or items == 'taprgh':
                            if line[0] > averagingtimestart and line[0] < averagingtimeend:
                                if line[0] > ftpcurrentlast:                                                                            
                                    ftpcurrent.append(line[0])
                                    ftpcurrentlast = line[0]
                                    pressure.append(line[1])
                                    resultpressure.write("%s %s\n" %(line[0],line[1]))
                        if items == 'U':
                            if line[0] > averagingtimestart and line[0] < averagingtimeend:
                                if line[0] > utcurrentlast:                                  
                                    utcurrent.append(line[0])
                                    utcurrentlast = line[0]
                                    resultvelocity.write("%s " %line[0])
                                    for i in range(1,4):
                                        velosity[i-1].append(line[i])
                                        resultvelocity.write("%s " %line[i])
                                    resultvelocity.write('\n') 
                sourcedata.close() 
    ####################writing pressure result file##################
    resultpressure.write("%s,%s,%s,%s\n" % ('#Averaging time is from ', ftpcurrent[0] , ' to ' , ftpcurrent[-1] ))
    averagepressure = np.average(pressure)
    resultpressure.write("%s, %s\n" %("#Aceraged pressure ", averagepressure))
    maxpressure = np.max(pressure)
    minpressure = np.min(pressure)
    resultpressure.write("%s, %s, %s ,%s\n" %("#Maximum pressure ", maxpressure, "Minimum pressure ", minpressure))
    resultpressure.close()        
    ####################writing velocity result file##################
    resultvelocity.write("%s,%s,%s,%s\n" % ('#Averaging time is from ',utcurrent[0] ,' to ',utcurrent[-1] ))
    resultvelocity.write("%s\n" %variablesu[1:4]) 
    resultvelocity.write("#Aceraged Velocity ")             
    for i in range(3):
        averagevelocity[i] = np.average(velosity[i])
    resultvelocity.write("%s\n" %averagevelocity)
    resultvelocity.close()
                      
    #########clearing lists############################
    resultpressure.close()
    resultpressure = open(str(pressure_file),'r+')
    primpressure2 = 0.0
    for line in resultpressure:
        if line[0] != '#' and  line[0] != '[':
            line = line.replace('\n', ' ')
            line = line.replace('(', ' ')
            line = line.replace(')', ' ')
            line = line.replace(',', ' ')
            line = [float(s) for s in line.split()]  
            primpressure2 = ((line[1])- (averagepressure))**2 + (primpressure2)
    resultpressure.write("%s, %s\n" %("#Pressure Prime in power 2 average ", primpressure2/len(pressure)))
    resultpressure.close()
                      
#########clearing lists############################
pressureall = []
first_tap = True
for ii in range(1,13):
    pressure_file = base_dir / "postProcessing" / f"tap{ii}" / f"tap{ii}pressure.dat"
    if pressure_file.exists():
        sourcedata = open(str(pressure_file), 'r')
        counter = 0
        for line in sourcedata:
            counter += 1
            if line[0] != '#' and line[0] != '[':
                line = line.replace('\n', ' ')
                line = line.replace('(', ' ')
                line = line.replace(')', ' ')
                line = line.replace('[', ' ')
                line = line.replace(']', ' ')                              
                line = line.replace(',', ' ')
                line = [float(s) for s in line.split()]
                if first_tap:
                    pressureall.append([line[0]])
                else:
                    pressureall[counter-1].append(line[-1])
        sourcedata.close()
        if first_tap:
            pressureall.insert(0, ['time'])
            first_tap = False
        pressureall[0].append('P'+str(ii)+'')           
with open(str(base_dir / "postProcessing" / "ALL_POINTS_pressure.csv"), 'w', newline='') as dat_data:	
    writer = csv.writer(dat_data)
    writer.writerows(pressureall)
dat_data.close()