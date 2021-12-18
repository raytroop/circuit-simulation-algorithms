#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import sys
import numpy
import matplotlib.pyplot as plt
import analogdef as ana
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
#
#
modeldict=ana.readmodelfile('models.txt')
ICdict={}
Plotdict={}
Printdict={}
Optionsdict={}
SetupDict={}
SimDict={}
#
#
DeviceCount=ana.readnetlist('netlist_4p2.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#    
MatrixSize=DeviceCount+len(Nodes) 
#
#
STA_matrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_rhs=[0 for i in range(MatrixSize)]
sol=[0 for i in range(MatrixSize)]
#
#
NumberOfNodes=len(Nodes)

for i in range(DeviceCount):
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
        STA_rhs[NumberOfNodes+i]=DevValue[i]
    if DevType[i]=='CurrentSource':
        if DevNode1[i] != '0' :
            STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        if DevNode2[i] != '0' :
            STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
        STA_rhs[NumberOfNodes+i]=0 
    if DevType[i]=='transistor':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=1/DevValue[i]
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode3[i])]=-1/DevValue[i]
        STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        STA_matrix[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1
#
#           
val=[[0 for i in range(100)] for j in range(MatrixSize)]
freqpnts=[0 for i in range(100)]
for iter in range(100):
    omega=iter*1e8*2*3.14159265
    for i in range(DeviceCount):
        if DevType[i]=='capacitor':
            if DevNode1[i] != '0' : 
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=DevValue[i]*omega
            if DevNode2[i] != '0' :
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-DevValue[i]*omega
        if DevType[i]=='inductor':
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]*omega
    STA_inv=numpy.linalg.inv(STA_matrix)
    sol=numpy.matmul(STA_inv,STA_rhs)
    freqpnts[iter]=iter*1e8
    for j in range(MatrixSize):
        val[j][iter]=abs(sol[j])
        
ana.plotdata(Plotdict,NumberOfNodes,freqpnts,val,Nodes)
if len(Printdict)> 0:
    ana.printdata(Printdict,NumberOfNodes,freqpnts,val,Nodes)        
plt.title('Voltage vs frequency')
plt.xlabel('frequency [Hz]')
plt.ylabel('|Voltage| [V]')  
plt.show