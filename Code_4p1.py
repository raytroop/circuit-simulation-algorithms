#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy as np
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
Optdict={}
#
#
DeviceCount=ana.readnetlist('netlist_4p1.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#    
MatrixSize=DeviceCount+len(Nodes) 
#
#
STA_matrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_rhs=[0 for i in range(MatrixSize)]
#
#
NumberOfNodes=len(Nodes)
for i in range(DeviceCount):
    if DevType[i]=='capacitor' or DevType[i]=='inductor':
        DevValue[i]*=(0+1j)
    STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=-DevValue[i]        
    if DevNode1[i] != '0' :
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=1
        STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
    if DevNode2[i] != '0' :
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-1
        STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
    if DevType[i]=='capacitor':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
        if DevNode1[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=-DevValue[i]
        if DevNode2[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=+DevValue[i]
    if DevType[i]=='VoltSource':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=0
        STA_rhs[NumberOfNodes+i]=DevValue[i]
#
#        
val=[0 for i in range(100)]        
for iter in range(100):
    omega=iter*1e9*2*3.14159265
    for i in range(DeviceCount):
        if DevType[i]=='capacitor':
            if DevNode1[i] != '0' : 
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=DevValue[i]*omega
            if DevNode2[i] != '0' :
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-DevValue[i]*omega
        if DevType[i]=='inductor':
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]*omega
    STA_inv=np.linalg.inv(STA_matrix)
    sol=np.matmul(STA_inv,STA_rhs)
    val[iter]=abs(sol[2])

plt.title('Voltage vs frequency')
plt.xlabel('frequency [GHz]')
plt.ylabel('|Voltage| [V]')  
plt.plot(val)