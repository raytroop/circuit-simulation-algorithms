#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy as np
import analogdef as ana

"""
Newton-Raphson algorithm with three regions MOS model, i.e. CMOS model2
1) cut off region, when Vgs - Vt < 0:
 Id=0
2) triode region, when Vds < Vgs - Vt:
 Id=2K((Vgs-Vt)Vds-0.5*Vds^2)
3) saturation region
 Id=K(Vgs-Vt)^2*(1+lamda*Vds)
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

modeldict=ana.readmodelfile('models.txt')
ICdict={}
Plotdict={}
Printdict={}
Optdict={}
Optdict['MaxNewtonIterations']=int(5)
#
#
DeviceCount=ana.readnetlist('netlist_5p1.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
# netlist_5p3 is simple PMOS follower
#DeviceCount=ana.readnetlist('netlist_5p3.txt',modeldict,ICdict,Plotdict,Printdict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
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
        VD=sol[Nodes.index(DevNode1[i])]
        VG=sol[Nodes.index(DevNode2[i])]
        VS=sol[Nodes.index(DevNode3[i])]
        Vgs=VG-VS
        Vds=VD-VS
        if DevModel[i][0]=='p':
            Vds=-Vds
            Vgs=-Vgs
        if Vds < Vgs-VT :
            STA_nonlinear[NumberOfNodes+i]=2*((Vgs-VT)*Vds-0.5*Vds**2)
        else :
            STA_nonlinear[NumberOfNodes+i]=(Vgs-VT)**2*(1+lambdaT*Vds)

#
f=np.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
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
                lambdaT=ana.findParameter(modeldict,DevModel[i],'lambdaT')
                VT=ana.findParameter(modeldict,DevModel[i],'VT')
                VD=sol[Nodes.index(DevNode1[i])]
                VG=sol[Nodes.index(DevNode2[i])]
                VS=sol[Nodes.index(DevNode3[i])]
                Vgs=VG-VS
                Vds=VD-VS
                if DevModel[i][0]=='p':
                    Vds=-Vds
                    Vgs=-Vgs
                if Vgs<VT:
                    STA_nonlinear[NumberOfNodes+i]=1e-5
                elif Vds < Vgs-VT:
                    STA_nonlinear[NumberOfNodes+i]=2*((Vgs-VT)*Vds-0.5*Vds**2)
                else :
                    STA_nonlinear[NumberOfNodes+i]=(Vgs-VT)**2*(1+lambdaT*Vds)
    f=np.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
#
#
    for i in range(MatrixSize):
        for j in range(MatrixSize):
            Jacobian[i][j]=STA_matrix[i][j]
    for i in range(DeviceCount):
        if DevType[i]=='transistor':
            lambdaT=ana.findParameter(modeldict,DevModel[i],'lambdaT')
            VT=ana.findParameter(modeldict,DevModel[i],'VT')
            Jacobian[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]
            VD=sol[Nodes.index(DevNode1[i])]
            VG=sol[Nodes.index(DevNode2[i])]
            VS=sol[Nodes.index(DevNode3[i])]
            Vgs=VG-VS
            Vds=VD-VS
            Vgd=VG-VD
            if DevModel[i][0]=='p':
                PFET=-1
                Vgs=-Vgs
                Vds=-Vds
                Vgd=-Vgd
            else:
                PFET=1
            if Vgs<VT :
                Jacobian[NumberOfNodes+i][Nodes.index(DevNode1[i])]=PFET*1e-1
                Jacobian[NumberOfNodes+i][Nodes.index(DevNode2[i])]=PFET*1e-1
                Jacobian[NumberOfNodes+i][Nodes.index(DevNode3[i])]=-PFET*1e-1
                Jacobian[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
                Jacobian[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1
            elif Vds <= Vgs-VT:
                Jacobian[NumberOfNodes+i][Nodes.index(DevNode1[i])]=PFET*2*(Vgd-VT)
                Jacobian[NumberOfNodes+i][Nodes.index(DevNode2[i])]=PFET*2*Vds
                Jacobian[NumberOfNodes+i][Nodes.index(DevNode3[i])]=-PFET*2*(Vgs-VT)
                Jacobian[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
                Jacobian[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1
            else :
                Jacobian[NumberOfNodes+i][Nodes.index(DevNode1[i])]=PFET*lambdaT*(Vgs-VT)**2
                Jacobian[NumberOfNodes+i][Nodes.index(DevNode2[i])]=PFET*2*(Vgs-VT)*(1+lambdaT*Vds)
                Jacobian[NumberOfNodes+i][Nodes.index(DevNode3[i])]=PFET*(-2*(Vgs-VT)*(1+lambdaT*Vds)-lambdaT*(Vgs-VT)**2)
                Jacobian[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
                Jacobian[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1
    sol=sol-np.matmul(np.linalg.inv(Jacobian),f)
    Jac_inv=np.linalg.inv(Jacobian)
    for j in range(MatrixSize):
        val[j][Newtoniter+1]=sol[j]

ana.plotdata(Plotdict,NumberOfNodes,Iteration,val,Nodes)
ana.printdata(Printdict,NumberOfNodes,Iteration,val,Nodes)