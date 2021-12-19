#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy as np
import matplotlib.pyplot as plt
import math
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
Optionsdict['deltaT']=1e-12
Optionsdict['NIterations']=200
#
#
DeviceCount=ana.readnetlist('netlist_4p7.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#
MatrixSize=DeviceCount+len(Nodes)
STA_matrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_rhs=[0 for i in range(MatrixSize)]
sol=[0 for i in range(MatrixSize)]
solm1=[0 for i in range(MatrixSize)]
#
#
deltaT=Optionsdict['deltaT']
NIterations=int(Optionsdict['NIterations'])
#
#
NumberOfNodes=len(Nodes)
if len(ICdict)>0:
    for i in range(len(ICdict)):
        for j in range(NumberOfNodes):
            if Nodes[j]==ICdict[i]['NodeName']:
                sol[j]=ICdict[i]['Value']
                print('Setting ',Nodes[j],' to ',sol[j])
#
#
for i in range(DeviceCount):
    if DevType[i] != 'VoltSource' and DevType[i] != 'CurrentSource':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=-DevValue[i]
    if DevNode1[i] != '0' :
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=1
        STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
    if DevNode2[i] != '0' :
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-1
        STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
    if DevType[i]=='capacitor':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
        if DevNode1[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=-3/2.0*DevValue[i]/deltaT
        if DevNode2[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=3/2.0*DevValue[i]/deltaT
        if DevNode1[i] != '0' : STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        if DevNode2[i] != '0' : STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
        if DevNode1[i] != '0' and DevNode2[i] != '0':
            STA_rhs[NumberOfNodes+i]=DevValue[i]/deltaT*(-2*(sol[Nodes.index(DevNode1[i])]-sol[Nodes.index(DevNode2[i])])+1/2*(solm1[Nodes.index(DevNode1[i])]-solm1[Nodes.index(DevNode2[i])]) )
        if DevNode1[i] == '0':
            STA_rhs[NumberOfNodes+i]=DevValue[i]/deltaT*(-2*(-sol[Nodes.index(DevNode2[i])])+1/2*(-solm1[Nodes.index(DevNode2[i])] ))
        if DevNode2[i] == '0':
            STA_rhs[NumberOfNodes+i]=DevValue[i]/deltaT*(-2*(sol[Nodes.index(DevNode1[i])])+1/2*(solm1[Nodes.index(DevNode1[i])] ))
    if DevType[i]=='inductor':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
        if DevNode1[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=-2/3*deltaT/DevValue[i]
        if DevNode2[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=2/3*deltaT/DevValue[i]
        if DevNode1[i] != '0' : STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        if DevNode2[i] != '0' : STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
        STA_rhs[NumberOfNodes+i]=4/3*sol[NumberOfNodes+i]-1/3*solm1[NumberOfNodes+i]
    if DevType[i]=='VoltSource':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=0
        STA_rhs[NumberOfNodes+i]=ana.getSourceValue(DevValue[i],0)
#
#
val=[[0 for i in range(NIterations)] for j in range(MatrixSize)]
timeVector=[0 for i in range(NIterations)]
for iter in range(NIterations):
    SimTime=iter*deltaT
    STA_inv=np.linalg.inv(STA_matrix)
    solm1=sol[:]
    sol=np.matmul(STA_inv,STA_rhs)
    timeVector[iter]=SimTime
    for j in range(MatrixSize):
        val[j][iter]=sol[j]
    for i in range(DeviceCount):
        if DevType[i]=='capacitor':
            if DevNode1[i] != '0' and DevNode2[i] != '0':
                STA_rhs[NumberOfNodes+i]=DevValue[i]/deltaT*(-2*(sol[Nodes.index(DevNode1[i])]-sol[Nodes.index(DevNode2[i])])+1/2*(solm1[Nodes.index(DevNode1[i])]-solm1[Nodes.index(DevNode2[i])]) )
            if DevNode1[i] == '0':
                STA_rhs[NumberOfNodes+i]=DevValue[i]/deltaT*(-2*(-sol[Nodes.index(DevNode2[i])])+1/2*(-solm1[Nodes.index(DevNode2[i])] ))
            if DevNode2[i] == '0':
                STA_rhs[NumberOfNodes+i]=DevValue[i]/deltaT*(-2*(sol[Nodes.index(DevNode1[i])])+1/2*(solm1[Nodes.index(DevNode1[i])] ))
        if DevType[i]=='inductor':
            STA_rhs[NumberOfNodes+i]=4/3*sol[NumberOfNodes+i]-1/3*solm1[NumberOfNodes+i]
        if DevType[i]=='VoltSource':
            STA_rhs[NumberOfNodes+i]=ana.getSourceValue(DevValue[i],SimTime)
ana.plotdata(Plotdict,NumberOfNodes,timeVector,val,Nodes)