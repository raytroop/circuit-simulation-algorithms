#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy
import matplotlib.pyplot as plt
import math
import sys
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
vkmax=0
#
#
modeldict=ana.readmodelfile('models.txt')
ICdict={}
Plotdict={}
Printdict={}
Optionsdict={}
Optionsdict['reltol']=1e-2
Optionsdict['iabstol']=1e-7
Optionsdict['vabstol']=1e-2
Optionsdict['lteratio']=2
Optionsdict['deltaT']=1e-12
Optionsdict['NIterations']=200
Optionsdict['GlobalTruncation']=True
#
#
DeviceCount=ana.readnetlist('netlist_4p7.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#    
MatrixSize=DeviceCount+len(Nodes)
STA_matrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_rhs=[0 for i in range(MatrixSize)]
sol=[0 for i in range(MatrixSize)]
solold=[0 for i in range(MatrixSize)]
solm1=[0 for i in range(MatrixSize)]
solm2=[0 for i in range(MatrixSize)]
#
#
deltaT=Optionsdict['deltaT']
NIterations=int(Optionsdict['NIterations'])
GlobalTruncation=Optionsdict['GlobalTruncation']
PointLocal=not GlobalTruncation
reltol=Optionsdict['reltol']
iabstol=Optionsdict['iabstol']
vabstol=Optionsdict['vabstol']
lteratio=Optionsdict['lteratio']
#
#
NumberOfNodes=len(Nodes)
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
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=0
        if DevNode1[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=+1
        if DevNode2[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-1
        if DevNode1[i] != '0' : STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=+1
        if DevNode2[i] != '0' : STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
        STA_rhs[NumberOfNodes+i]=sol[Nodes.index(DevNode1[i])]-sol[Nodes.index(DevNode2[i])]+deltaT*sol[NumberOfNodes+i]/DevValue[i]
    if DevType[i]=='inductor':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
        if DevNode1[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=0
        if DevNode2[i] != '0' : STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=0
        if DevNode1[i] != '0' : STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        if DevNode2[i] != '0' : STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
        STA_rhs[NumberOfNodes+i]=sol[NumberOfNodes+i]+(sol[Nodes.index(DevNode1[i])]-sol[Nodes.index(DevNode2[i])])*deltaT/DevValue[i]
    if DevType[i]=='VoltSource':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=0
        STA_rhs[NumberOfNodes+i]=ana.getSourceValue(DevValue[i],0)
#
#        
val=[[0 for i in range(NIterations)] for j in range(MatrixSize)]
timeVector=[0 for i in range(NIterations)]
PredMatrix=[[0 for i in range(2)] for j in range(2)]
Predrhs=[0 for i in range(2)]
for iter in range(NIterations):
    SimTime=iter*deltaT       
    STA_inv=numpy.linalg.inv(STA_matrix)
    solm2=[solm1[i] for i in range(MatrixSize)]
    solm1=[solold[i] for i in range(MatrixSize)]
    solold=[sol[i] for i in range(MatrixSize)]
    sol=numpy.matmul(STA_inv,STA_rhs)
    for node in range(NumberOfNodes):
        vkmax=max(vkmax,abs(sol[node]))
    timeVector[iter]=SimTime
    for j in range(MatrixSize):
        val[j][iter]=sol[j]
    for i in range(DeviceCount):
        if DevType[i]=='capacitor':
            STA_rhs[NumberOfNodes+i]=sol[Nodes.index(DevNode1[i])]-sol[Nodes.index(DevNode2[i])]+deltaT*sol[NumberOfNodes+i]/DevValue[i]
        if DevType[i]=='inductor':
            STA_rhs[NumberOfNodes+i]=sol[NumberOfNodes+i]+(sol[Nodes.index(DevNode1[i])]-sol[Nodes.index(DevNode2[i])])*deltaT/DevValue[i]
        if DevType[i]=='VoltSource':
            STA_rhs[NumberOfNodes+i]=ana.getSourceValue(DevValue[i],SimTime)
    if iter>2:
        LTEConverged=True
        for i in range(NumberOfNodes):
            tau1=(timeVector[iter-2]-timeVector[iter-3])
            tau2=(timeVector[iter-1]-timeVector[iter-3])
            PredMatrix[0][0]=tau2
            PredMatrix[0][1]=tau2*tau2
            PredMatrix[1][0]=tau1
            PredMatrix[1][1]=tau1*tau1
            Predrhs[0]=solold[i]-solm2[i]
            Predrhs[1]=solm1[i]-solm2[i]
            Predsol=numpy.matmul(numpy.linalg.inv(PredMatrix),Predrhs)
            vpred=solm2[i]+Predsol[0]*(SimTime-timeVector[iter-3])+Predsol[1]*(SimTime-timeVector[iter-3])*(SimTime-timeVector[iter-3])
            if PointLocal:
                for node in range(NumberOfNodes):
                    vkmax=max(abs(sol[node]),abs(solm1[node]))
                    print('Is vkmax correct here?')
                    if abs(vpred-sol[i])> lteratio*(vkmax*reltol+vabstol):
                        LTEConverged=False
            elif GlobalTruncation:
                for node in range(NumberOfNodes):
                    if abs(vpred-sol[i])> lteratio*(vkmax*reltol+vabstol):
                        LTEConverged=False
            else:
                print('Error: Unknown truncation error')
                sys.exit()
        if not LTEConverged:
            print('LTE NOT converging, change time step')
            sys.exit(0)

ana.plotdata(Plotdict,NumberOfNodes,timeVector,val,Nodes)    