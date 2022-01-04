#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy as np
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt
import math
import analogdef as ana

'''
 PNOISE implementation
'''

##
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
Optionsdict['NHarmonics']=16
Optionsdict['Period']=1/5032661878.243104
Optionsdict['PNOISE']='False'
Optionsdict['iabstol']=1e-11
#
#
DeviceCount=ana.readnetlist('netlist_5p17.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
#
#
NHarmonics=int(Optionsdict['NHarmonics'])
Period=Optionsdict['Period']
run_PNOISE=(Optionsdict['PNOISE']=='True')
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
f=np.zeros(MatrixSize)+1j*np.zeros(MatrixSize)
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
#
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

for i in range(DeviceCount):
    if DevType[i] == 'oscfilter':
        OscFilterIndex=i
    if(DevLabel[i] == 'vinp'):
        StimulusIndex=i
ana.build_SysEqns_HB(SetupDict, SimDict, modeldict)
#
#
drivenAmp = []
drivenIrms = []
for AmpIndex in range(1):
    STA_rhs[(NumberOfNodes+StimulusIndex)*TotalHarmonics+Jacobian_Offset+1]=.23577+float (AmpIndex)/10000
    STA_rhs[(NumberOfNodes+StimulusIndex)*TotalHarmonics+Jacobian_Offset-1]=.23577+float (AmpIndex)/10000
    NewIter=15
    Newtoniter=0
    f[0]=1j
    while Newtoniter < NewIter and max(abs(f)) > iabstol:
        print('NewtonIteration :',Newtoniter,max(abs(f)))
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

    f=np.matmul(STA_matrix,sol)-STA_rhs+STA_nonlinear
    for j in range(TotalHarmonics):
        Vp[j]=sol[3*TotalHarmonics+j]
        Vn[j]=sol[4*TotalHarmonics+j]
        IOscFilterSpec[j]=sol[(NumberOfNodes+OscFilterIndex)*TotalHarmonics+j]
    VnTime=ana.idft(Vn,TotalHarmonics)
    VpTime=ana.idft(Vp,TotalHarmonics)
    IOscFilter=ana.idft(IOscFilterSpec,TotalHarmonics)

    drivenAmp_i = STA_rhs[(NumberOfNodes+StimulusIndex)*TotalHarmonics+Jacobian_Offset+1]
    drivenIrms_i= np.real(np.sqrt(np.mean(IOscFilter**2)))
    drivenAmp.append(drivenAmp_i)
    drivenIrms.append(drivenIrms_i)

    print('rms ',drivenIrms_i, drivenAmp_i)

FIGNUM = 0
plt.figure(FIGNUM)
FIGNUM += 1
plt.plot(abs(VnTime))
plt.title('Steady state harmonics')
plt.xlabel('time [s]')
plt.ylabel('Output Amplitude [V]')

if len(drivenAmp) > 1:
    plt.figure(FIGNUM)
    FIGNUM += 1
    plt.plot(drivenAmp, drivenIrms)
    plt.title('residue current through the flter vs driving voltage amplitude')
    plt.xlabel('Driven Amplitude [V]')
    plt.ylabel('Driven Current, rms[A]')

if run_PNOISE:
    #
    #
    Jacobian=SetupDict['Jacobian']
    Jacobian[(NumberOfNodes+OscFilterIndex)*TotalHarmonics+Jacobian_Offset-1][(NumberOfNodes+OscFilterIndex)*TotalHarmonics+Jacobian_Offset-1]=-1e18
    Jacobian[(NumberOfNodes+OscFilterIndex)*TotalHarmonics+Jacobian_Offset+1][(NumberOfNodes+OscFilterIndex)*TotalHarmonics+Jacobian_Offset+1]=-1e18
    #
    #
    STA_rhs=[0 for i in range(MatrixSize)]
    val=[[0 for i in range(100)] for j in range(4)]
    for iter in range(100):
        omega=iter*1e6*2*math.pi
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
                    STA_rhs[(NumberOfNodes+i)*TotalHarmonics+row]=.5*(row==Jacobian_Offset+1)+.5*(row==Jacobian_Offset-1)
        sol=np.matmul(np.linalg.inv(Jacobian),STA_rhs)
        val[0][iter]=abs(sol[3*TotalHarmonics+Jacobian_Offset])
        val[1][iter]=abs(sol[3*TotalHarmonics+Jacobian_Offset+1])
        val[2][iter]=abs(sol[3*TotalHarmonics+Jacobian_Offset+2])
        val[3][iter]=20*math.log10(val[1][iter])
    plt.figure(FIGNUM)
    plt.plot(val[3])
    plt.title('PNOISE output of simple VCO analysis')
    plt.xlabel('Offset Freq from 1st Harmonic [MHz]')
    plt.ylabel('20*log10(pnoise) [arbitrary units]')

plt.show()