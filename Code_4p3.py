#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import analogdef as ana
import math
#
#
MaxNumberOfDevices=100
DevType=[0*i for i in range(MaxNumberOfDevices)]
DevLabel=[0*i for i in range(MaxNumberOfDevices)]
DevNode1=[0*i for i in range(MaxNumberOfDevices)]
DevNode2=[0*i for i in range(MaxNumberOfDevices)]
DevNode3=[0*i for i in range(MaxNumberOfDevices)]
DevModel=[0*i for i in range(MaxNumberOfDevices)]
DevValue=[0*i for i in range(MaxNumberOfDevices)]
Nodes=[]
k=1.3823e-23 
Temperature=300
NumberOfPoints=1000
#
#
modeldict=ana.readmodelfile('models.txt')
ICdict={}
Plotdict={}
Printdict={}
Optdict={}
#
#
DeviceCount=ana.readnetlist('netlist_4p3.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#    
DeviceCount=DeviceCount+1
MatrixSize=DeviceCount+len(Nodes) 
#
#
STA_matrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_rhs=[0 for i in range(MatrixSize)]
#
#
NumberOfNodes=len(Nodes)
for i in range(DeviceCount-1):
    if DevType[i]=='capacitor' or DevType[i]=='inductor':
        DevValue[i]*=(0+1j)
    if DevType[i] == 'resistor' or DevType[i] == 'inductor':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=-DevValue[i]        
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=1
        STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-1
        STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
    if DevType[i]=='capacitor':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
        STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=-DevValue[i]
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=+DevValue[i]
    if DevType[i]=='VoltSource':
        if DevNode1[i]!= '0':
            STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=1
            STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        if DevNode2[i] != '0':
            STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-1
            STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=0
        STA_rhs[NumberOfNodes+i]=0
    if DevType[i]=='CurrentSource':
        if DevNode1[i] != '0' :
            STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        if DevNode2[i] != '0' :
            STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
        STA_rhs[NumberOfNodes+i]=0        
        
val=[[0 for j in range(NumberOfPoints)] for i in range(DeviceCount)]
freqpnts=[i*1e8 for i in range(NumberOfPoints)]
for NoiseSource in range(DeviceCount):
    if DevType[NoiseSource]=='resistor':
        if DevNode1[NoiseSource] != '0' :
            STA_matrix[Nodes.index(DevNode1[NoiseSource])][NumberOfNodes+DeviceCount-1]=1
        if DevNode2[NoiseSource] != '0' :
            STA_matrix[Nodes.index(DevNode2[NoiseSource])][NumberOfNodes+DeviceCount-1]=-1
        STA_matrix[NumberOfNodes+DeviceCount-1][NumberOfNodes+DeviceCount-1]=1
#
#        
        for iter in range(NumberOfPoints):
            omega=iter*1e8*2*3.14159265
            for i in range(DeviceCount):
                if DevType[i]=='capacitor':
                    if DevNode1[i] != '0' : 
                        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=DevValue[i]*omega
                    if DevNode2[i] != '0' :
                        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-DevValue[i]*omega
                if DevType[i]=='inductor':
                    STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]*omega
                if DevType[i]=='resistor' and i==NoiseSource:
                    STA_rhs[NumberOfNodes+DeviceCount-1]=math.sqrt(4*k*Temperature/DevValue[i])
            sol=np.matmul(np.linalg.inv(STA_matrix),STA_rhs)
            val[NoiseSource][iter]=abs(sol[Nodes.index('out')])
        if DevNode1[NoiseSource] != '0' :
            STA_matrix[Nodes.index(DevNode1[NoiseSource])][NumberOfNodes+DeviceCount-1]=0
        if DevNode2[NoiseSource] != '0' :
            STA_matrix[Nodes.index(DevNode2[NoiseSource])][NumberOfNodes+DeviceCount-1]=0
TotalNoiseSpectrum=[0 for i in range(NumberOfPoints)]
for NoiseSource in range(DeviceCount):
    if DevType[NoiseSource]=='resistor':
        for i in range(NumberOfPoints):
            TotalNoiseSpectrum[i]+=abs(val[NoiseSource][i])*abs(val[NoiseSource][i])

plt.plot(freqpnts,TotalNoiseSpectrum)            
plt.title('Noise Power vs frequency')
plt.xlabel('frequency [Hz]')
plt.ylabel('Noise Power [V^2/Hz]')