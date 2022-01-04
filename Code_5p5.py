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
Homotopy Methods - gmin Stepping
"""
#
#
def f_NL(STA_matrix, STA_rhs, STA_nonlinear, solution):
    return numpy.matmul(STA_matrix,solution)-STA_rhs+STA_nonlinear

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
#
#
# Fig.5.13 -> Fig.5.14
DeviceCount=ana.readnetlist('netlist_5p6.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
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
deltaT=1e-12
sol=[0 for i in range(MatrixSize)]
solm1=[0 for i in range(MatrixSize)]
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
            # Do nothing
            STA_rhs[NumberOfNodes]=STA_rhs[NumberOfNodes]
        if DevType[i]=='inductor':
            # For DC we treat this as a voltage source with V=0
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=0
            STA_rhs[NumberOfNodes+i]=0
        if DevType[i]=='VoltSource':
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=0
            STA_rhs[NumberOfNodes+i]=DevValue[i]
        if DevType[i]=='CurrentSource':
            if DevNode1[i] != '0' :
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=0
                STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=0
            if DevNode2[i] != '0' :
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=0
                STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=0
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
            STA_rhs[NumberOfNodes+i]=DevValue[i]
            if DevNode1[i] != '0' and DevNode2[i]!='0':
                STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
                STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
            elif DevNode2[i] != '0' :
                STA_matrix[Nodes.index(DevNode2[i])][NumberOfNodes+i]=-1
            elif DevNode1[i] != '0' :
                STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
    if DevType[i]=='transistor':
        if DevModel[i][0]=='p':
            PFET=1
        else:
            PFET=1
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]
        if DevNode1[i] != '0' :
            STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=0
            STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        if DevNode3[i] != '0' :
            STA_matrix[NumberOfNodes+i][Nodes.index(DevNode3[i])]=0
            STA_matrix[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1
        if DevNode1[i] != '0' :
            if DevNode2[i] != '0' and DevNode3[i] != '0':
                STA_nonlinear[NumberOfNodes+i]=(sol[Nodes.index(DevNode2[i])]-sol[Nodes.index(DevNode3[i])])**2
            if DevNode2[i] == '0':
                STA_nonlinear[Nodes.index(DevNode1[i])]=sol[Nodes.index(DevNode3[i])]**2
                STA_nonlinear[NumberOfNodes+1]=(sol[Nodes.index(DevNode3[i])])**2
            if DevNode3[i] == '0':
                STA_nonlinear[Nodes.index(DevNode1[i])]=sol[Nodes.index(DevNode2[i])]**2
        if DevNode3[i] != '0' :
            if DevNode2[i] != '0':
                STA_nonlinear[NumberOfNodes+i]=(sol[Nodes.index(DevNode2[i])]-sol[Nodes.index(DevNode3[i])])**2
            else:
                STA_nonlinear[Nodes.index(DevNode3[i])]=-sol[Nodes.index(DevNode3[i])]**2
#
f=numpy.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
#
#
SimTime=0
NewIter=15
val=[0 for i in range(NewIter)]
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
                if DevLabel[i]=='vinp':
                        STA_rhs[NumberOfNodes+i]=2*math.sin(2*math.pi*1e9*SimTime)+1
                if DevLabel[i]=='vinn':
                        STA_rhs[NumberOfNodes+i]=-2*math.sin(2*math.pi*1e9*SimTime)+1
        if DevType[i]=='transistor':
            if DevNode1[i] != '0' :
                if DevNode2[i] != '0' and DevNode3[i] != '0':
                    STA_nonlinear[NumberOfNodes+i]=(sol[Nodes.index(DevNode2[i])]-sol[Nodes.index(DevNode3[i])])**2
                if DevNode2[i] == '0':
                    STA_nonlinear[NumberOfNodes+i]=sol[Nodes.index(DevNode3[i])]^2
                if DevNode3[i] == '0':
                    STA_nonlinear[NumberOfNodes+i]=sol[Nodes.index(DevNode2[i])]^2
    f=numpy.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
#
#
    for i in range(MatrixSize):
        for j in range(MatrixSize):
            Jacobian[i][j]=STA_matrix[i][j]
    for i in range(DeviceCount):
        if DevType[i]=='transistor':
            Jacobian[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]
            if DevNode1[i] != '0' :
                if DevNode2[i] != '0' and DevNode3[i] != '0':
                    Jacobian[NumberOfNodes+i][Nodes.index(DevNode2[i])]=2*PFET*(sol[Nodes.index(DevNode2[i])]-sol[Nodes.index(DevNode3[i])])
                    Jacobian[NumberOfNodes+i][Nodes.index(DevNode3[i])]=-2*PFET*(sol[Nodes.index(DevNode2[i])]-sol[Nodes.index(DevNode3[i])])
                    Jacobian[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
                    Jacobian[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1
                elif DevNode2[i] == '0':
                    Jacobian[NumberOfNodes+i][Nodes.index(DevNode3[i])]=-1
                    Jacobian[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
                    Jacobian[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1
                elif DevNode3[i] == '0':
                    Jacobian[NumberOfNodes+i][Nodes.index(DevNode2[i])]=1
                    Jacobian[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
                    Jacobian[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1
    sol=sol-numpy.matmul(numpy.linalg.inv(Jacobian),f)
    Jac_inv=numpy.linalg.inv(Jacobian)
    val[Newtoniter]=sol[2]

plt.plot(val)
plt.show()











































