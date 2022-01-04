#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy as np
from numpy import linalg
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import colors as mcolors
import math
import analogdef as ana

'''
 harmonic balance simulation implementation with PAC
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
#
#
modeldict=ana.readmodelfile('models.txt')
ICdict={}
Plotdict={}
Printdict={}
Optionsdict={}
SetupDict={}
Optionsdict['NHarmonics']=32
Optionsdict['Period']=1e-9
Optionsdict['PAC']='False'
Optionsdict['MaxNewtonIterations']=15
Optionsdict['iabstol']=1e-7
#
# DeviceCount=ana.readnetlist('netlist_5p7.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
# DeviceCount=ana.readnetlist('netlist_5p15.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
DeviceCount=ana.readnetlist('netlist_5p16.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#
NHarmonics=Optionsdict['NHarmonics']
Period=Optionsdict['Period']
run_PAC=(Optionsdict['PAC']=='True')
NSamples=2*(NHarmonics-1)
TotalHarmonics=2*NHarmonics-1
NumberOfNodes=len(Nodes)
MatrixSize=(DeviceCount+len(Nodes))*TotalHarmonics
Jacobian=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
Jac_inv=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
Spare=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_matrix=[[0 for i in range(MatrixSize)] for j in range(MatrixSize)]
STA_rhs=[0 for i in range(MatrixSize)]
STA_nonlinear=[0 for i in range(MatrixSize)]
f=[0 for i in range(MatrixSize)]
Template=[[0 for i in range(TotalHarmonics)] for j in range(TotalHarmonics)]

Jacobian_Offset=int(TotalHarmonics/2)
omegak=[0 for i in range(TotalHarmonics)]
HarmonicsList=[0 for i in range(TotalHarmonics)]
Samples=[i*Period/NSamples for i in range(NSamples)]
#
#
SetupDict['NumberOfNodes']=NumberOfNodes
SetupDict['DeviceCount']=DeviceCount
SetupDict['DevValue']=DevValue
SetupDict['TotalHarmonics']=TotalHarmonics
SetupDict['Jacobian_Offset']=Jacobian_Offset
SetupDict['NSamples']=NSamples
SetupDict['DevNode1']=DevNode1
SetupDict['DevNode2']=DevNode2
SetupDict['DevNode3']=DevNode3
SetupDict['DevType']=DevType
SetupDict['DevLabel']=DevLabel
SetupDict['DevModel']=DevModel
SetupDict['Nodes']=Nodes
SetupDict['STA_matrix']=STA_matrix
SetupDict['Jacobian']=Jacobian
SetupDict['STA_rhs']=STA_rhs
SetupDict['STA_nonlinear']=STA_nonlinear
SetupDict['omegak']=omegak
iabstol=Optionsdict['iabstol']

sol=np.zeros(MatrixSize)+1j*np.zeros(MatrixSize)
SolutionCorrection=np.zeros(MatrixSize)+1j*np.zeros(MatrixSize)
TransistorOutputTime=[0 for i in range(NSamples)]
TransistorOutputTimeDerivative=[0 for i in range(TotalHarmonics)]
TransistorOutputFreq=[0 for i in range(TotalHarmonics)]
Jlkm=[0 for i in range(TotalHarmonics)]
Jlko=[0 for i in range(TotalHarmonics)]
Vg=[0 for i in range(TotalHarmonics)]
Vs=[0 for i in range(TotalHarmonics)]
Vd=[0 for i in range(TotalHarmonics)]
VgTime=[0 for i in range(NSamples)]
VsTime=[0 for i in range(NSamples)]
VdTime=[0 for i in range(NSamples)]
gm=[0 for i in range(TotalHarmonics)]
Vp=[0 for i in range(TotalHarmonics)]
Vn=[0 for i in range(TotalHarmonics)]
VpTime=[0 for i in range(NSamples)]
VnTime=[0 for i in range(NSamples)]
IOscFilterSpec=[0 for i in range(TotalHarmonics)]
IOscFilter=[0 for i in range(NSamples)]
#
SimDict={}
SimDict['sol']=sol
#
#
for i in range(2*NHarmonics-1):
    HarmonicsList[i]=i
for i in range(TotalHarmonics):
    omegak[i]=(1-NHarmonics+i)*2*math.pi/Period

ana.build_SysEqns_HB(SetupDict, SimDict, modeldict)

f=np.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
#
#
TimePnts=[i*Period/NSamples for i in range(NSamples)]
NewIter=int(Optionsdict['MaxNewtonIterations'])
Newtoniter=0
while Newtoniter < NewIter and abs(max(f)) > iabstol:
    print('NewtonIteration :',Newtoniter,abs(max(f)))
    #
    #
    for i in range(MatrixSize):
        STA_nonlinear[i]=0
    ana.update_SysEqns_HB(SetupDict, SimDict, modeldict)
    f=np.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear

#
#
    for i in range(MatrixSize):
        for j in range(MatrixSize):
            Jacobian[i][j]=STA_matrix[i][j]
    ana.build_Jacobian_HB(SetupDict, SimDict, modeldict)

    SolutionCorrection=np.matmul(np.linalg.inv(Jacobian),f)
    for i in range(MatrixSize):
        sol[i]=sol[i]-SolutionCorrection[i]
    Newtoniter=Newtoniter+1

try:
    for j in range(TotalHarmonics):
        # always save node 'outp' as Vntime
        Vn[j]=sol[Nodes.index('outp')*TotalHarmonics+j]
except:
    for j in range(TotalHarmonics):
        # guess output node is named as 'out' instead of 'outp'
        Vn[j]=sol[Nodes.index('out')*TotalHarmonics+j]
VnTime=ana.idft(Vn,TotalHarmonics)

if run_PAC:
    _, axs = plt.subplots(nrows=2)
else:
    _, axs = plt.subplots(nrows=1)

try:
    axs[0].plot(TimePnts,VnTime)
    axs[0].set_title('outp Voltage vs time')
    axs[0].set_xlabel('time [s]')
    axs[0].set_ylabel('outp Voltage [V]')
except:
    axs.plot(TimePnts,VnTime)
    axs.set_title('Voltage vs time')
    axs.set_xlabel('time [s]')
    axs.set_ylabel('Voltage [V]')

if run_PAC:
    STA_rhs=[0 for i in range(MatrixSize)]
    NUM_FREQ_STEP = 100
    NUM_HARMONIC_PAC = 4
    val=[[0 for i in range(NUM_FREQ_STEP)] for j in range(NUM_HARMONIC_PAC)]
    freqx=[[0 for i in range(NUM_FREQ_STEP)] for j in range(NUM_HARMONIC_PAC)]
    for iter in range(NUM_FREQ_STEP):
        omega=iter*1e6*2*math.pi
        print('Frequency sweep(Hz):',iter*1e6)
        for i in range(DeviceCount):
            for row in range(TotalHarmonics):
                if DevType[i]=='capacitor':
                    if DevNode1[i] != '0' :
                        Jacobian[(NumberOfNodes+i)*TotalHarmonics+row][Nodes.index(DevNode1[i])*TotalHarmonics+row]=1j*(omegak[row]+(np.sign(omegak[row])+(omegak[row]==0))*omega)*DevValue[i]
                    if DevNode2[i] != '0' :
                        Jacobian[(NumberOfNodes+i)*TotalHarmonics+row][Nodes.index(DevNode2[i])*TotalHarmonics+row]=-1j*(omegak[row]+(np.sign(omegak[row])+(omegak[row]==0))*omega)*DevValue[i]
                if DevType[i]=='inductor':
                    Jacobian[(NumberOfNodes+i)*TotalHarmonics+row][(NumberOfNodes+i)*TotalHarmonics+row]=-1j*(omegak[row]+(np.sign(omegak[row])+(omegak[row]==0))*omega)*DevValue[i]
                if DevType[i]=='CurrentSource':
                    if DevLabel[i]=='i1':
                        STA_rhs[(NumberOfNodes+i)*TotalHarmonics+row]=1*(row==Jacobian_Offset)
                    else:
                        STA_rhs[(NumberOfNodes+i)*TotalHarmonics+row]=-(row==Jacobian_Offset)
        sol=np.matmul(np.linalg.inv(Jacobian),STA_rhs)
        val[0][iter]=abs(sol[6*TotalHarmonics+Jacobian_Offset])
        val[1][iter]=abs(sol[6*TotalHarmonics+Jacobian_Offset+1])
        val[2][iter]=abs(sol[6*TotalHarmonics+Jacobian_Offset+2])
        val[3][iter]=abs(sol[6*TotalHarmonics+Jacobian_Offset+3])

        freqx[0][iter] = (0/Period + iter*1e6)
        freqx[1][iter] = (1/Period + iter*1e6)
        freqx[2][iter] = (2/Period + iter*1e6)
        freqx[3][iter] = (3/Period + iter*1e6)

    lines = [list(zip(freqx[i], val[i])) for i in range(4)]
    colors = [mcolors.to_rgba(c) for c in plt.rcParams['axes.prop_cycle'].by_key()['color']]
    lc = mc.LineCollection(lines, linewidths=2, colors=colors)
    axs[1].add_collection(lc)
    axs[1].autoscale()
    axs[1].set_title('PAC simulation')
    axs[1].set_xlabel('frequency [Hz]')
    axs[1].set_ylabel('TransConductance [S] vs Harmonic')

plt.show()
