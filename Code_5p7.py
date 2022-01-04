#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy as np
import sys
import analogdef as ana
import matplotlib.pyplot as plt

'''
 shooting method, supporting capacitors as the dynamic element (inductor NOT supported).
'''
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
Optionsdict['method']='be'
Optionsdict['Period']=1e-9
#
#
# DeviceCount=ana.readnetlist('netlist_5p7.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
DeviceCount=ana.readnetlist('netlist_5p14.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#
MatrixSize=DeviceCount+len(Nodes)
NumberOfNodes=len(Nodes)
NumberOfCurrents=DeviceCount
Jacobian=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
Jac_inv=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
Spare=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_matrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_rhs=[0 for i in range(MatrixSize)]
STA_nonlinear=[0 for i in range(MatrixSize)]
CapMatrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
FrachetMatrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
for i in range(MatrixSize):
    FrachetMatrix[i][i]=1
IdentityMatrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
for i in range(MatrixSize):
    IdentityMatrix[i][i]=1
f=[0 for i in range(MatrixSize)]
maxSol=[0 for i in range(10)]
#
deltaT=Optionsdict['deltaT']
NIterations=int(Optionsdict['NIterations'])
GlobalTruncation=Optionsdict['GlobalTruncation']
PointLocal=not GlobalTruncation
reltol=Optionsdict['reltol']
iabstol=Optionsdict['iabstol']
vabstol=Optionsdict['vabstol']
lteratio=Optionsdict['lteratio']
method=Optionsdict['method']
Period=Optionsdict['Period']
#
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
SetupDict['method']=method
SetupDict['reltol']=reltol
SetupDict['iabstol']=iabstol
SetupDict['vabstol']=vabstol
SetupDict['lteratio']=lteratio
SetupDict['GlobalTruncation']=GlobalTruncation
SetupDict['PointLocal']=PointLocal
SetupDict['vkmax']=vkmax
#
#
deltaT=Optionsdict['deltaT']
sol=np.zeros(MatrixSize)
solm1=np.zeros(MatrixSize)
soltemp=np.zeros(MatrixSize)
solInit=np.zeros(MatrixSize)
SimDict={}
SimDict['deltaT']=deltaT
SimDict['sol']=sol
SimDict['solm1']=solm1
SimDict['solInit']=solInit
SimDict['soltemp']=soltemp
SimDict['f']=f
#
#
ana.build_SysEqns(SetupDict, SimDict, modeldict)
f=np.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
#
#
Npnts=int(Period/deltaT)
val=[[0 for i in range(Npnts)] for j in range(MatrixSize)]
vin=[0 for i in range(20)]
for ShootingIter in range(10):
    for i in range(MatrixSize):
        sol[i]=solInit[i]
    timepnts=[i*deltaT for i in range(Npnts)]
    for iter in range(Npnts):
        if iter%100==0:
            print("Iter=",iter)
        SimTime=iter*deltaT
        for i in range(MatrixSize):
            soltemp[i]=sol[i]
        NewIter=500
        ResidueConverged=False
        UpdateConverged=False
        NewtonConverged=False
        NewtonIter=0
    #
    #
        while not NewtonConverged and NewtonIter<50:
            NewtonIter=NewtonIter+1

            for i in range(MatrixSize):
                STA_nonlinear[i]=0
            ana.update_SysEqns(SimTime, SetupDict, SimDict, modeldict)
            SimDict['f']=f=np.matmul(STA_matrix,soltemp)-STA_rhs+STA_nonlinear

            ResidueConverged=ana.DidResidueConverge(SetupDict, SimDict)
    #
    #
            for i in range(MatrixSize):
                for j in range(MatrixSize):
                    Jacobian[i][j]=STA_matrix[i][j]
            SetupDict['Jacobian']=Jacobian
            ana.build_Jacobian(SetupDict, SimDict, modeldict)
            UpdateConverged=ana.DidUpdateConverge(SetupDict, SimDict)
            Jacobian=SetupDict['Jacobian']
            SolutionCorrection=np.matmul(np.linalg.inv(Jacobian),f)
            for i in range(MatrixSize):
                soltemp[i]=soltemp[i]-SolutionCorrection[i]

            NewtonConverged=ResidueConverged and UpdateConverged
        if NewtonConverged:
            for i in range(MatrixSize):
                solm1[i]=sol[i]
            for node in range(NumberOfNodes):
                vkmax=max(vkmax,abs(sol[node]))
                SetupDict['vkmax']=vkmax
            f=np.matmul(STA_matrix,soltemp)-STA_rhs+STA_nonlinear
            for i in range(MatrixSize):
                sol[i]=soltemp[i]
                val[i][iter]=sol[i]
            #
            #
            for i in range(DeviceCount):
                if DevType[i]=='capacitor':
                    CapMatrix[Nodes.index(DevNode1[i])][Nodes.index(DevNode2[i])]=DevValue[i]
                    CapMatrix[Nodes.index(DevNode2[i])][Nodes.index(DevNode1[i])]=DevValue[i]
                elif DevType[i]=='inductor':
                    print('Error: this shooting code only implements capacitors as dynamic elements\n')
                    sys.exit(1)
            FrachetMatrix=np.matmul(np.linalg.inv(Jacobian),np.matmul(CapMatrix,FrachetMatrix))/deltaT
        else:
            print('Newtoniteration did not converge.\nexiting ...\n')
            sys.exit(0)
    #
    #
    maxSol[ShootingIter]=max(abs(sol-solInit))
    #
    #
    solInit=solInit+np.matmul(np.linalg.inv(IdentityMatrix-FrachetMatrix),(sol-solInit))
ana.plotdata(Plotdict,NumberOfNodes,timepnts,val,Nodes)
plt.figure(2)
plt.plot(maxSol)
plt.xlabel('Iteration')
plt.ylabel('Max Difference IC')
plt.show()
#
if len(Printdict)> 0:
    ana.printdata(Printdict,NumberOfNodes,timepnts,val,Nodes)