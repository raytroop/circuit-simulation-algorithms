#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy
import matplotlib.pyplot as plt
import math
import analogdef as ana

"""
Newton-Raphson algorithm with simplest MOS model I=K*(Vg-Vs)^2, i.e. CMOS model1
"""

#
#
DeviceCount=0
MaxNumberOfDevices=100
DevType=[0*i for i in range(MaxNumberOfDevices)]
DevLabel=[0*i for i in range(MaxNumberOfDevices)]
DevNode1=[0*i for i in range(MaxNumberOfDevices)]
DevNode2=[0*i for i in range(MaxNumberOfDevices)]
DevNode3=[0*i for i in range(MaxNumberOfDevices)]
DevValue=[0*i for i in range(MaxNumberOfDevices)]
DevModel=[0*i for i in range(MaxNumberOfDevices)]
Nodes=[]
#
modeldict=ana.readmodelfile('models.txt')
ICdict={}
Plotdict={}
Printdict={}
Optdict={}
Optdict['MaxNewtonIterations']=int(5)
#
#
DeviceCount=ana.readnetlist('netlist_5p1.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#
NumberOfNodes=len(Nodes)
MatrixSize=DeviceCount+len(Nodes)
Jacobian=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
Jac_inv=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
Spare=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_matrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_rhs=[0 for i in range(MatrixSize)]
STA_nonlinear=[0 for i in range(MatrixSize)]
f=[0 for i in range(MatrixSize)]
#
#
sol=[0 for i in range(MatrixSize)]
solm1=[0 for i in range(MatrixSize)]
#
#
if len(ICdict)>0:
    for i in range(len(ICdict)):
        for j in range(NumberOfNodes):
            if Nodes[j]==ICdict[i]['NodeName']:
                sol[j]=ICdict[i]['Value']
                print('Setting ',Nodes[j],' to ',sol[j])
#
#
for i in range(DeviceCount):
    if DevType[i] != 'transistor':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=-DevValue[i]
        if DevNode1[i] != '0' :
            STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=1
            STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        if DevNode2[i] != '0' :
            STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-1
            STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
        if DevType[i]=='capacitor':
            # Do nothing since DC sim
            STA_rhs[NumberOfNodes]=STA_rhs[NumberOfNodes]
        if DevType[i]=='inductor':
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=0
            STA_rhs[NumberOfNodes+i]=0
        if DevType[i]=='VoltSource':
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=0
            STA_rhs[NumberOfNodes+i]=ana.getSourceValue(DevValue[i],0)
        if DevType[i]=='CurrentSource':
            if DevNode1[i] != '0' :
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=0
                STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=0
            if DevNode2[i] != '0' :
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=0
                STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=0
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
            STA_rhs[NumberOfNodes+i]=ana.getSourceValue(DevValue[i],0)
            if DevNode1[i] != '0' and DevNode2[i]!='0':
                STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
                STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
            elif DevNode2[i] != '0' :
                STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
            elif DevNode1[i] != '0' :
                STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
    if DevType[i]=='transistor':
        lambdaT=ana.findParameter(modeldict,DevModel[i],'lambdaT')
        VT=ana.findParameter(modeldict,DevModel[i],'VT')
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=0
        STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode3[i])]=0
        STA_matrix[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1
        VG=sol[Nodes.index(DevNode2[i])]
        VS=sol[Nodes.index(DevNode3[i])]
        Vgs=VG-VS
        if DevModel[i][0]=='p':
            Vgs=-Vgs
        STA_nonlinear[NumberOfNodes+i]=Vgs**2

f=numpy.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
#
#
NewIter=int(Optdict['MaxNewtonIterations'])
val=[[0 for i in range(NewIter+1)] for j in range(MatrixSize)]
for j in range(MatrixSize):
    val[j][0]=sol[j]
Iteration=[i for i in range(NewIter+1)]
for Newtoniter in range(NewIter):
    for i in range(MatrixSize):
        STA_nonlinear[i]=0
    for i in range(DeviceCount):
        if DevType[i]!='transistor':
            if DevType[i]=='capacitor':
                STA_rhs[NumberOfNodes+i]=STA_rhs[NumberOfNodes+i]
            if DevType[i]=='inductor':
                STA_rhs[NumberOfNodes+i]=0
            if DevType[i]=='VoltSource':
                STA_rhs[NumberOfNodes+i]=ana.getSourceValue(DevValue[i],0)
            if DevType[i]=='CurrentSource':
                STA_rhs[NumberOfNodes+i]=ana.getSourceValue(DevValue[i],0)
        if DevType[i]=='transistor':
                VG=sol[Nodes.index(DevNode2[i])]
                VS=sol[Nodes.index(DevNode3[i])]
                Vgs=VG-VS
                if DevModel[i][0]=='p':
                    Vgs=-Vgs
                STA_nonlinear[NumberOfNodes+i]=Vgs**2
    f=numpy.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
#
#
    for i in range(MatrixSize):
        for j in range(MatrixSize):
            Jacobian[i][j]=STA_matrix[i][j]
    for i in range(DeviceCount):
        if DevType[i]=='transistor':
            Jacobian[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]
            VG=sol[Nodes.index(DevNode2[i])]
            VS=sol[Nodes.index(DevNode3[i])]
            Vgs=VG-VS
            if DevModel[i][0]=='p':
                PFET=-1
                Vgs=-Vgs
            else:
                PFET=1
            Jacobian[NumberOfNodes+i][Nodes.index(DevNode2[i])]=2*PFET*Vgs
            Jacobian[NumberOfNodes+i][Nodes.index(DevNode3[i])]=-2*PFET*Vgs
    sol=sol-numpy.matmul(numpy.linalg.inv(Jacobian),f)
    Jac_inv=numpy.linalg.inv(Jacobian)
    for j in range(MatrixSize):
        val[j][Newtoniter+1]=sol[j]

ana.plotdata(Plotdict,NumberOfNodes,Iteration,val,Nodes)
ana.printdata(Printdict,NumberOfNodes,Iteration,val,Nodes)