#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy as np
import scipy.linalg as slin
import matplotlib.pyplot as plt
import math
import sys
import analogdef as ana

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
Optdict={}
Optdict['reltol']=1e-2
Optdict['iabstol']=1e-7
Optdict['vabstol']=1e-2
Optdict['lteratio']=2
Optdict['MaxTimeStep']=1e-11
Optdict['FixedTimeStep']='False'
Optdict['GlobalTruncation']='True'
Optdict['deltaT']=3e-13
Optdict['MaxSimulationIterations']=200000
Optdict['MaxSimTime']=1e-8
Optdict['ThreeLevelStep']='True'
Optdict['method']='trap'
Optdict['CheckLTE']='True'
#
#
'''
5.4 Transient Nonlinear Simulation
As a result, we have Code 5.6 in Sect. 5.6.6. This code will now be used for all our subsequent simulations using transient technique
'''
# DeviceCount=ana.readnetlist('netlist_5p7.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
# DeviceCount=ana.readnetlist('netlist_5p8.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
# DeviceCount=ana.readnetlist('netlist_5p9.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
# DeviceCount=ana.readnetlist('netlist_5p10.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
# DeviceCount=ana.readnetlist('netlist_5p11.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
# DeviceCount=ana.readnetlist('netlist_5p12.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
DeviceCount=ana.readnetlist('netlist_5p13.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#
reltol=Optdict['reltol']
iabstol=Optdict['iabstol']
vabstol=Optdict['vabstol']
lteratio=Optdict['lteratio']
MaxTimeStep=Optdict['MaxTimeStep']
FixedTimeStep=(Optdict['FixedTimeStep']=='True')
GlobalTruncation=(Optdict['GlobalTruncation']=='True')
deltaT=Optdict['deltaT']
MaxSimulationIterations=int(Optdict['MaxSimulationIterations'])
ThreeLevelStep=(Optdict['ThreeLevelStep']=='True')
MaxSimTime=Optdict['MaxSimTime']
PointLocal=not GlobalTruncation
method=Optdict['method']
CheckLTE=(Optdict['CheckLTE']=='True')
#
#
NumberOfNodes=len(Nodes)
NumberOfCurrents=DeviceCount
MatrixSize=DeviceCount+len(Nodes)
Jacobian=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
Jac_inv=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
Spare=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_matrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_rhs=[0 for i in range(MatrixSize)]
STA_nonlinear=[0 for i in range(MatrixSize)]
f=[0 for i in range(MatrixSize)]
SetupDict={}
SetupDict['NumberOfNodes']=NumberOfNodes
SetupDict['NumberOfCurrents']=NumberOfCurrents
SetupDict['DeviceCount']=DeviceCount
SetupDict['Nodes']=Nodes
SetupDict['DevNode1']=DevNode1
SetupDict['DevNode2']=DevNode2
SetupDict['DevNode3']=DevNode3
SetupDict['DevValue']=DevValue
SetupDict['DevType']=DevType
SetupDict['DevModel']=DevModel
SetupDict['MatrixSize']=MatrixSize
SetupDict['Jacobian']=Jacobian
SetupDict['STA_matrix']=STA_matrix
SetupDict['STA_rhs']=STA_rhs
SetupDict['STA_nonlinear']=STA_nonlinear
SetupDict['FixedTimeStep']=FixedTimeStep
SetupDict['method']=method
SetupDict['GlobalTruncation']=GlobalTruncation
SetupDict['PointLocal']=PointLocal
SetupDict['vkmax']=vkmax
SetupDict['Vthermal']=1.38e-23*300/1.602e-19
SetupDict['reltol']=reltol
SetupDict['iabstol']=iabstol
SetupDict['vabstol']=vabstol
SetupDict['lteratio']=lteratio
SetupDict['MaxTimeStep']=MaxTimeStep
#
#
sol=[0 for i in range(MatrixSize)]
solm1=[0 for i in range(MatrixSize)]
solm2=[0 for i in range(MatrixSize)]
soltemp=[0 for i in range(MatrixSize)]
SimDict={}
SimDict['deltaT']=deltaT
SimDict['ThreeLevelStep']=ThreeLevelStep
SimDict['sol']=sol
SimDict['solm1']=solm1
SimDict['solm2']=solm2
SimDict['soltemp']=soltemp
SimDict['f']=f
#
#
if len(ICdict)>0:
    for i in range(len(ICdict)):
        for j in range(NumberOfNodes):
            if Nodes[j]==ICdict[i]['NodeName']:
                sol[j]=ICdict[i]['Value']
                solm1[j]=ICdict[i]['Value']
                solm2[j]=ICdict[i]['Value']
                print('Setting ',Nodes[j],' to ',sol[j])
#
#
ana.build_SysEqns(SetupDict, SimDict, modeldict)
#
f=np.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
#
#
val=[[0 for i in range(MaxSimulationIterations)] for j in range(MatrixSize)]
vin=[0 for i in range(MaxSimulationIterations)]
timeVector=[0 for i in range(MaxSimulationIterations)]
TotalIterations=0
iteration=0
SimTime=0
NewtonIter=0
NewtonConverged=False
LTEIter=0
Converged=False
for i in range(MatrixSize):
    soltemp[i]=sol[i]
#
#
while SimTime<MaxSimTime and iteration<MaxSimulationIterations:
    if iteration%100==0:
        print("Iter=",iteration,NewtonIter,LTEIter,deltaT,vkmax,SimTime)
    SimTime=SimTime+deltaT
    NewtonIter=0
    LTEIter=0
    ResidueConverged=False
    UpdateConverged=False
    NewtonConverged=False
    Converged=False
    while (not NewtonConverged) and NewtonIter<50:
        NewtonIter=NewtonIter+1
        for i in range(MatrixSize):
            STA_nonlinear[i]=0
        ana.update_SysEqns(SimTime, SetupDict, SimDict, modeldict)
        SimDict['f']=f=np.matmul(STA_matrix,soltemp)-STA_rhs+STA_nonlinear
        ResidueConverged=ana.DidResidueConverge(SetupDict, SimDict)
#
        for i in range(MatrixSize):
            for j in range(MatrixSize):
                Jacobian[i][j]=STA_matrix[i][j]
        ana.build_Jacobian(SetupDict, SimDict, modeldict)
        UpdateConverged=ana.DidUpdateConverge(SetupDict, SimDict)
        SolutionCorrection=np.matmul(np.linalg.inv(Jacobian),f)
        for i in range(MatrixSize):
            soltemp[i]=soltemp[i]-SolutionCorrection[i]
        NewtonConverged=ResidueConverged and UpdateConverged
#
    LTEConverged, MaxLTERatio=ana.DidLTEConverge(SetupDict, SimDict, iteration, LTEIter, NewtonConverged, timeVector, SimTime, SolutionCorrection)
    if not CheckLTE:
        LTEConverged=True
    if not LTEConverged:
        NewtonIter=0
        for i in range(MatrixSize):
            soltemp[i]=sol[i]
#
    deltaT, iteration, SimTime, Converged=ana.UpdateTimeStep(SetupDict, SimDict, LTEConverged, NewtonConverged, val, iteration, NewtonIter, MaxLTERatio, timeVector, SimTime)
    TotalIterations=TotalIterations+NewtonIter

    SimDict['deltaT']=deltaT
    if Converged:
        for i in range(MatrixSize):
            sol[i]=soltemp[i]
        for node in range(NumberOfNodes):
            vkmax=max(vkmax,abs(sol[node]))
            SetupDict['vkmax']=vkmax
    if deltaT<1e-15:
        print('Warning: Timestep too short: ',deltaT)
        sys.exit(0)

reval=[[0 for i in range(iteration)] for j in range(MatrixSize)]
retime=[0 for i in range(iteration)]
logvalue=[0 for i in range(iteration)]
for i in range(iteration):
    for j in range(MatrixSize):
        reval[j][i]=val[j][i]
    retime[i]=timeVector[i]

ana.plotdata(Plotdict,NumberOfNodes,retime,reval,Nodes)
if len(Printdict)> 0:
    ana.printdata(Printdict,NumberOfNodes,retime,reval,Nodes)
print('TotalIterations ',TotalIterations)