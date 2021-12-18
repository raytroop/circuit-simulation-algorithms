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
Ports=['vp1','vp2']
PortNodes=['in','out']
NPorts=2
#
#
modeldict=ana.readmodelfile('models.txt')
ICdict={}
Plotdict={}
Printdict={}
Optionsdict={}
#
#
DeviceCount=ana.readnetlist('netlist_4p6.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
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
#
#   
val=[[[0 for i in range(1000)] for j in range(NPorts)] for k in range(NPorts)]
for port in range(NPorts):      
    for iter in range(1000):
        omega=iter*1e8*2*3.14159265
        for i in range(DeviceCount):
            if DevType[i]=='capacitor':
                if DevNode1[i] != '0' : 
                    STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=DevValue[i]*omega
                if DevNode2[i] != '0' :
                    STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-DevValue[i]*omega
            if DevType[i]=='inductor':
                STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]*omega
            if DevLabel[i]==Ports[port]:
                print('Exciting port:',DevLabel[i])
                STA_rhs[NumberOfNodes+i]=2
            else:
                STA_rhs[NumberOfNodes+i]=0
        STA_inv=np.linalg.inv(STA_matrix)
        sol=np.matmul(STA_inv,STA_rhs)
        for j in range(NPorts):
            print('Sniffing port: ',PortNodes[j])
            if port != j :
                val[port][j][iter]=20*math.log10(abs(sol[Nodes.index(PortNodes[j])]))
            else :
                val[port][j][iter]=20*math.log10(abs(sol[Nodes.index(PortNodes[j])]-1))
plt.plot(val[0][0])
plt.plot(val[0][1])
plt.plot(val[1][0])
plt.plot(val[1][1])
plt.title('S-parameters vs frequency')
plt.xlabel('frequency [100MHz]')
plt.ylabel('S-par [dB]')  
plt.show
        
        

    

    

    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
