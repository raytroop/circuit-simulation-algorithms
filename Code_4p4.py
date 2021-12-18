#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 22:33:04 2019

@author: mikael
"""
import numpy as np
import matplotlib.pyplot as plt
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
FreqStep=3e7
#
#
modeldict=ana.readmodelfile('models.txt')
ICdict={}
Plotdict={}
Printdict={}
Optionsdict={}
#
#
DeviceCount=ana.readnetlist('netlist_4p5.txt',modeldict,ICdict,Plotdict,Printdict,Optionsdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNumberOfDevices)
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
        STA_rhs[NumberOfNodes+i]=0 # No AC source, so put to zero
    if DevType[i]=='transistor':
        STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=1
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=1/DevValue[i]
        STA_matrix[NumberOfNodes+i][Nodes.index(DevNode3[i])]=-1/DevValue[i]
        STA_matrix[Nodes.index(DevNode1[i])][NumberOfNodes+i]=1
        STA_matrix[Nodes.index(DevNode3[i])][NumberOfNodes+i]=-1    
#
#
print('Setting up stab run 1')
for i in range(DeviceCount):
    if DevType[i]=='VoltSource':
        if DevLabel[i]=='vstab':
            STA_rhs[NumberOfNodes+i]=1
            VElabel=DevNode1[i]
            StabProbeIndex=i
            print('Found stability probe')
        else:
            STA_rhs[NumberOfNodes+i]=0    
#
#           
val=[[0 for i in range(100)] for j in range(MatrixSize)]
freqpnts=[0 for i in range(100)]
D=[0+0j for i in range(100)]
B=[0+0j for i in range(100)]
for iter in range(100):
    omega=iter*FreqStep*2*3.14159265
    for i in range(DeviceCount):
        if DevType[i]=='capacitor':
            if DevNode1[i] != '0' : 
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=DevValue[i]*omega
            if DevNode2[i] != '0' :
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-DevValue[i]*omega
        if DevType[i]=='inductor':
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]*omega
    STA_inv=np.linalg.inv(STA_matrix)
    sol=np.matmul(STA_inv,STA_rhs)
    freqpnts[iter]=iter*FreqStep
    for j in range(MatrixSize):
        val[j][iter]=abs(sol[j])
        if j<NumberOfNodes:
            if Nodes[j]==VElabel:
                D[iter]=sol[j] 
        if j==NumberOfNodes+StabProbeIndex:
            B[iter]=sol[j]
            

print('Setting up stab run 2')                
for i in range(DeviceCount):
    if DevType[i]=='VoltSource':
        STA_rhs[NumberOfNodes+i]=0    
    if DevType[i]=='CurrentSource':
        if DevLabel[i]=='istab':
            STA_rhs[NumberOfNodes+i]=1
            print('Found stability current probe')
#
#           
val=[[0 for i in range(100)] for j in range(MatrixSize)]
freqpnts=[0 for i in range(100)]
C=[0+0j for i in range(100)]
A=[0+0j for i in range(100)]
for iter in range(100):
    omega=iter*FreqStep*2*3.14159265
    for i in range(DeviceCount):
        if DevType[i]=='capacitor':
            if DevNode1[i] != '0' : 
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode1[i])]=DevValue[i]*omega
            if DevNode2[i] != '0' :
                STA_matrix[NumberOfNodes+i][Nodes.index(DevNode2[i])]=-DevValue[i]*omega
        if DevType[i]=='inductor':
            STA_matrix[NumberOfNodes+i][NumberOfNodes+i]=DevValue[i]*omega
    STA_inv=np.linalg.inv(STA_matrix)
    sol=np.matmul(STA_inv,STA_rhs)
    freqpnts[iter]=iter*FreqStep
    for j in range(MatrixSize):
        val[j][iter]=abs(sol[j])
        if j<NumberOfNodes:
            if Nodes[j]==VElabel:
                C[iter]=sol[j]
        if j==NumberOfNodes+StabProbeIndex:
            A[iter]=-sol[j] 

T=[0+0j for i in range(100)]                            
magdB=[0 for i in range(100)]
phasedegree=[0 for i in range(100)]
for iter in range(100):
    T[iter]=(2*(A[iter]*D[iter]-B[iter]*C[iter])-A[iter]+D[iter])/(2*(B[iter]*C[iter]-A[iter]*D[iter])+A[iter]-D[iter]+1)
    magdB[iter]=20*np.log10(np.abs(T[iter]))
    if (180/np.pi*np.arctan(np.imag(T[iter])/np.real(T[iter])))>=0:
        phasedegree[iter]=180-180/np.pi*np.arctan(np.imag(T[iter])/np.real(T[iter]))
    else:
        phasedegree[iter]=180-(180+180/np.pi*np.arctan(np.imag(T[iter])/np.real(T[iter])))

for i in range(100-1):
    if magdB[i]*magdB[i+1]<0:
        print('Phasemargin: ',phasedegree[i])
        print('UGF ',freqpnts[i])

for i in range(100-1):
    if phasedegree[i]*phasedegree[i+1]<0:
        print('Gainmargin: ',-magdB[i])
        
#
plt.xscale('log')
plt.title('Gain/Phase vs Frequency')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [degrees], Gain [dB]')
#
plt.plot(freqpnts,magdB,label='Gain')
plt.plot(freqpnts,phasedegree,label='Phase')
plt.subplot(111).legend(loc='upper center', bbox_to_anchor=(0.8, 0.97), shadow=True)
plt.show()
#
if len(Printdict)> 0:
    ana.printdata(Printdict,NumberOfNodes,freqpnts,val,Nodes)        
