#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:04:22 2020

@author: mikael
"""
import sys
import re
import math
import matplotlib.pyplot as plt
import numpy
from scipy.fftpack import fft, ifft

Vthermal=1.38e-23*300/1.602e-19
TINY=1e-5





def readnetlist(netlist,modeldict,ICdict,Plotdict,Writedict,Optiondict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MaxNDevices):
    try:
        myfile=open(netlist,'r')
    except:
        print('netlist file',netlist,' not found')
        sys.exit()
    DeviceCount=0
    if len(modeldict)==0:
        print('Warning: model dictionary is empty!')
    line=myfile.readline()
    while line !='' :
        DevType[DeviceCount]='empty'
        if line[0]=='*':
            print('comment')
        if line[0]=='v':
            print('VoltSource')
            DevType[DeviceCount]='VoltSource'
        if line[0]=='i':
            print('CurrSource')
            DevType[DeviceCount]='CurrentSource'
        if line[0]=='r':
            print('resistor')
            DevType[DeviceCount]='resistor'
        if line[0]=='f':
            print('Special Oscillator filter found')
            DevType[DeviceCount]='oscfilter'
        if line[0]=='l':
            print('inductor')
            DevType[DeviceCount]='inductor'
        if line[0]=='c':
            print('capacitor')
            DevType[DeviceCount]='capacitor'
        if line[0]=='m':
            print('transistor')
            DevType[DeviceCount]='transistor'
        if line[0]=='q':
            print('bipolar')
            DevType[DeviceCount]='bipolar'
        if re.split(' ',line)[0]=='.ic':
            print('Initial Condition Statement')
            lineSplit=re.split(' ',line)
            for i in range(len(lineSplit)-1):
                ConditionSplit=re.split('\(|\)|=|\n',re.split(' ',line)[i+1])
                if len(ConditionSplit)>2:
                    try:
                        ICdict[i]={}
                        ICdict[i]['NodeName']=ConditionSplit[1]
                        ICdict[i]['Value']=float(ConditionSplit[3])
                    except:
                        print('Syntax Error in .ic statement')
                        sys.exit()
                else:
                    print('Warning: Odd characters in IC statement \'',ConditionSplit,'\'')
        if re.split(' ',line)[0]=='.plot':
            print('Plot Statement')
            lineSplit=re.split(' ',line)
            for i in range(len(lineSplit)-1):
                ConditionSplit=re.split('\(|\)|=|\n',re.split(' ',line)[i+1])
                if len(ConditionSplit)>2:
                    try:
                        Plotdict[i]={}
                        Plotdict[i]['NodeName']=ConditionSplit[1]
                    except:
                        print('Syntax Error in .plot statement')
                        sys.exit()
                else:
                    print('Warning: Odd characters in .plot statement \'',ConditionSplit,'\'')
        if re.split(' ',line)[0]=='.write':
            print('Write Statement')
            lineSplit=re.split(' ',line)
            Writedict[0]={}
            Writedict[0]['filename']=lineSplit[1]
            for i in range(len(lineSplit)-2):
                ConditionSplit=re.split('\(|\)|=|\n',re.split(' ',line)[i+2])
                if len(ConditionSplit)>2:
                    try:
                        Writedict[i+1]={}
                        Writedict[i+1]['NodeName']=ConditionSplit[1]
                    except:
                        print('Syntax Error in .write statement')
                        sys.exit()
                else:
                    print('Warning: Odd characters in .write statement \'',ConditionSplit,'\'')
        if re.split(' ',line)[0]=='.options':
            print('Option Statement')
            lineSplit=re.split(' ',line)
            for i in range(len(lineSplit)-1):
                ConditionSplit=re.split('=|\n',re.split(' ',line)[i+1])
                if len(ConditionSplit)>=2:
                    try:
                        Optiondict[ConditionSplit[0]]=float(ConditionSplit[1])
                    except:
                        try:
                            Optiondict[ConditionSplit[0]]=ConditionSplit[1]
                        except:
                            print('Syntax Error in .options statement')
                            sys.exit()
                else:
                    print('Warning: Odd characters in .options statement \'',ConditionSplit,'\'')
            # overrdie "float"
            if 'NHarmonics' in Optiondict:
                Optiondict['NHarmonics'] = int(Optiondict['NHarmonics'])
        if DevType[DeviceCount]!='empty':
            if DevType[DeviceCount] != 'transistor' and DevType[DeviceCount] != 'bipolar':
                try:
                    DevLabel[DeviceCount]=line.split(' ')[0]
                except:
                    print('Syntax Error in line:',line)
                    sys.exit();
                try:
                    DevNode1[DeviceCount]=line.split(' ')[1]
                except:
                    print('Syntax Error in line:',line)
                    sys.exit()
                if DevType[DeviceCount] != 'VoltSource' and DevType[DeviceCount] != 'CurrentSource' and DevNode1[DeviceCount]=='0':
                    print('Error: Node \'0\' only allowed for independent sources')
                    print('line',line)
                    sys.exit()
                try:
                    DevNode2[DeviceCount]=line.split(' ')[2]
                except:
                    print('Syntax Error in line:',line)
                    sys.exit()
                if DevType[DeviceCount] != 'VoltSource' and DevType[DeviceCount] != 'CurrentSource' and DevNode2[DeviceCount]=='0':
                    print('Error: Node \'0\' only allowed for independent sources')
                    print('line',line)
                    sys.exit()
                try:
                    DevValue[DeviceCount]=float(line.split(' ')[3])
                except:
                    print('Value is not a number')
                    if DevType[DeviceCount] != 'VoltSource' and DevType[DeviceCount] != 'CurrentSource':
                        sys.exit(0)
                    srcdict={}
                    try:
                        DevValue[DeviceCount]=re.split((' |\('),line)[3]
                    except:
                        print('Syntax Error in line:',line)
                        sys.exit();
                    srcdict[0]={}
                    srcdict[0]['type']=DevValue[DeviceCount]
                    if DevValue[DeviceCount]=='pwl':
                        DoneReadingPoints=False
                        pnt=1
                        while not DoneReadingPoints:
                            try:
                                TimePnt=float(re.split((' |\('),line)[2+pnt*2])
                            except:
                                DoneReadingPoints=True
                            if not DoneReadingPoints:
                                srcdict[pnt]={}
                                srcdict[pnt]['time']=TimePnt
                                try:
                                    SrcPnt=float(re.split((' |\(|\)'),line)[3+pnt*2])
                                except:
                                    print('Syntax Error in line:',line)
                                    sys.exit();
                                srcdict[pnt]['value']=SrcPnt
                                pnt=pnt+1
                    if DevValue[DeviceCount]=='sin':
#                        srcdict[1]={}
                        try:
                            srcdict['Offset']=float(re.split((' |\('),line)[4])
                        except:
                            print('Syntax Error in line:',line)
                            sys.exit();
                        try:
                            srcdict['Amplitude']=float(re.split((' |\('),line)[5])
                        except:
                            print('Syntax Error in line:',line)
                            sys.exit();
                        try:
                            srcdict['Freq']=float(re.split((' |\('),line)[6])
                        except:
                            print('Syntax Error in line:',line)
                            sys.exit();
                        try:
                            srcdict['TDelay']=float(re.split((' |\('),line)[7])
                        except:
                            print('Syntax Error in line:',line)
                            sys.exit();
                        try:
                            srcdict['Theta']=float(re.split((' |\(|\)'),line)[8])
                        except:
                            print('Syntax Error in line:',line)
                            sys.exit();
                    DevValue[DeviceCount]=srcdict
                if DevNode1[DeviceCount] not in Nodes and DevNode1[DeviceCount]!='0':
                    Nodes.append(DevNode1[DeviceCount])
                if DevNode2[DeviceCount] not in Nodes and DevNode2[DeviceCount]!='0':
                    Nodes.append(DevNode2[DeviceCount])
            else:
                try:
                    DevLabel[DeviceCount]=line.split(' ')[0]
                except:
                    print('Syntax Error in line:',line)
                    sys.exit();
                try:
                    DevNode1[DeviceCount]=line.split(' ')[1]
                except:
                    print('Syntax Error in line:',line)
                    sys.exit();
                if DevNode1[DeviceCount]=='0':
                    print('Error: Node \'0\' only allowed for independent sources')
                    print('line',line)
                    sys.exit()
                try:
                    DevNode2[DeviceCount]=line.split(' ')[2]
                except:
                    print('Syntax Error in line:',line)
                    sys.exit();
                if DevNode2[DeviceCount]=='0':
                    print('Error: Node \'0\' only allowed for independent sources')
                    print('line',line)
                    sys.exit()
                try:
                    DevNode3[DeviceCount]=line.split(' ')[3]
                except:
                    print('Syntax Error in line:',line)
                    sys.exit();
                if DevNode3[DeviceCount]=='0':
                    print('Error: Node \'0\' only allowed for independent sources')
                    print('line',line)
                    sys.exit()
                try:
                    DevModel[DeviceCount]=line.split(' ')[4]
                    DevModel[DeviceCount]=DevModel[DeviceCount].rstrip('\n')
                    modelIndex=findmodelIndex(modeldict,DevModel[DeviceCount])
                except:
                    print('Syntax Error in line4:',line)
                    sys.exit();
                DevValue[DeviceCount]=-1/modeldict[modelIndex]['K']#-1/K
                if DevNode1[DeviceCount] not in Nodes and DevNode1[DeviceCount]!='0':
                    Nodes.append(DevNode1[DeviceCount])
                if DevNode2[DeviceCount] not in Nodes and DevNode2[DeviceCount]!='0':
                    Nodes.append(DevNode2[DeviceCount])
                if DevNode3[DeviceCount] not in Nodes and DevNode3[DeviceCount]!='0':
                    Nodes.append(DevNode3[DeviceCount])
            DeviceCount+=1
            if DeviceCount>=MaxNDevices:
                print('Too many devices in the netlist: Max is set to ',MaxNDevices)
                sys.exit()
        line=myfile.readline()
    return DeviceCount

def readmodelfile(filename):
    modeldict={}
    index=0
    modelfile=open(filename,'r')
    line=modelfile.readline()
    while line != '' :
        modeldict[index]={}
        name=line.split(' ')[0]
        print('Reading model ',name)
        modeldict[index]['modelName']=name
        for i in range(3):
            dum=line.split(' ')[i+1]
            try:
                dum.index("=")
            except:
                print('Syntax error in model file, line ',line)
                sys.exit()
            Parname=dum.split('=')[0]
            try:
                ParValue=float(dum.split('=')[1])
            except:
                print('Syntax error: Parameter',Parname,' value is not a number',dum)
            modeldict[index][Parname]=ParValue
        index=index+1
        line=modelfile.readline()
    return modeldict

def findmodelIndex(modeldict,name):
    for i in range(len(modeldict)):
        if modeldict[i]['modelName'] == name:
            return i
    print('model name ',name,' is not found in modelfile' )
    sys.exit()


def getSourceValue(DevValue,SimTime):
    if type(DevValue)==float:
        return DevValue
    if type(DevValue)==dict:
        if DevValue[0]['type']=='sin':
            A=DevValue['Amplitude']
            freq=DevValue['Freq']
            Offset=DevValue['Offset']
            TDelay=DevValue['TDelay']
#            Theta=DevValue['Theta']
            return Offset+A*math.sin(freq*2*math.pi*(SimTime-TDelay))
        if DevValue[0]['type']=='pwl':
            TimeIndex=1
            while SimTime >= DevValue[TimeIndex]['time'] and TimeIndex<len(DevValue)-1:
                TimeIndex=TimeIndex+1
            if SimTime>=DevValue[len(DevValue)-1]['time']:
                return DevValue[len(DevValue)-1]['value']
            else:
                PrevTime=DevValue[TimeIndex-1]['time']
                NextTime=DevValue[TimeIndex]['time']
                PrevValue=DevValue[TimeIndex-1]['value']
                NextValue=DevValue[TimeIndex]['value']
                return (NextValue-PrevValue)*(SimTime-PrevTime)/(NextTime-PrevTime)+PrevValue


def findParameter(modeldict,modelname,parameterName):
    for i in range(len(modeldict)):
        if modeldict[i]['modelName']==modelname:
            try:
                return modeldict[i][parameterName]
            except:
                print('Error: Parameter ',parameterName,' not found in ',modeldict[i])


def setupDicts(SimDict,SetupDict,Optdict,DevType,DevValue,DevLabel,DevNode1,DevNode2,DevNode3,DevModel,Nodes,MatrixSize,Jacobian,STA_matrix,STA_rhs,STA_nonlinear,sol,solm1,solm2,f):
    SetupDict['NumberOfNodes']=10
    SetupDict['NumberOfCurrents']=10
    SetupDict['DeviceCount']=10
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
    SetupDict['Vthermal']=1.38e-23*300/1.602e-19
    Optdict['reltol']=1e-3
    Optdict['iabstol']=1e-7
    Optdict['vabstol']=1e-6
    Optdict['lteratio']=2
    Optdict['MaxTImeStep']=1e-11
    Optdict['FixedTimeStep']='False'
    Optdict['GlobalTruncation']='True'
    Optdict['deltaT']=3e-13
    Optdict['MaxSimulationIterations']=200000
    Optdict['MaxSimTime']=1e-8
    Optdict['MaxNewtonIter']=5
    SimDict['deltaT']=1e-12
    SimDict['sol']=sol
    SimDict['solm1']=solm1
    SimDict['solm2']=solm2
    SimDict['f']=f

def plotdata(Plotdict,NumberOfNodes,retime,reval,Nodes):
    if len(Plotdict)> 0:
        ax = plt.subplot(111)
        for j in range(NumberOfNodes):
            for i in range(len(Plotdict)):
                if Plotdict[i]['NodeName']==Nodes[j]:
                        ax.plot(retime, reval[j], label=Nodes[j])
        plt.title('Voltage vs time')
        ax.legend(loc='upper center', bbox_to_anchor=(1.2, 0.97), shadow=True, ncol=2)
        plt.xlabel('time [s]')
        plt.ylabel('Voltage [V]')
        plt.show()

def printdata(Printdict,NumberOfNodes,retime,reval,Nodes):
    if len(Printdict)> 0:
        fp=open(Printdict[0]['filename'],"w+")
        fp.write('time ')
        for i in range(len(Printdict)-1):
            fp.write('%s ' % Printdict[i+1]['NodeName'])
        fp.write('\n')
        for i in range(len(retime)):
            fp.write("%g " % retime[i])
            for j in range(NumberOfNodes):
                for k in range(len(Printdict)-1):
                    if Printdict[k+1]['NodeName']==Nodes[j]:
                        fp.write("%g " % reval[j][i])
            fp.write('\n')
        fp.close()


def build_SysEqns(SetupDict, SimDict, modeldict):
    DeviceCount=SetupDict['DeviceCount']
    DevType=SetupDict['DevType']
    deltaT=SimDict['deltaT']
    for i in range(DeviceCount):
        if DevType[i]=='resistor':
            build_SysEqn_resistor(i, SetupDict)
        if DevType[i]=='capacitor':
            build_SysEqn_capacitor(deltaT, i, SetupDict, SimDict)
        if DevType[i]=='inductor':
            build_SysEqn_inductor(deltaT, i, SetupDict, SimDict)
        if DevType[i]=='VoltSource':
            build_SysEqn_VSource(i, SetupDict)
        if DevType[i]=='CurrentSource':
            build_SysEqn_ISource(i, SetupDict)
        if DevType[i]=='transistor':
            build_SysEqn_MOS(i, SetupDict, SimDict, modeldict)
        if DevType[i]=='bipolar':
            build_SysEqn_bipolar(i, SetupDict, SimDict, modeldict)

def update_SysEqns(SimTime, SetupDict, SimDict, modeldict):
    DeviceCount=SetupDict['DeviceCount']
    DevType=SetupDict['DevType']
    deltaT=SimDict['deltaT']
    for i in range(DeviceCount):
        if DevType[i]=='capacitor':
            update_SysEqn_capacitor(deltaT, i, SetupDict, SimDict)
        if DevType[i]=='inductor':
            update_SysEqn_inductor(deltaT, i, SetupDict, SimDict)
        if DevType[i]=='VoltSource':
            update_SysEqn_VSource(i, SimTime, SetupDict)
        if DevType[i]=='CurrentSource':
            update_SysEqn_ISource(i, SimTime, SetupDict)
        if DevType[i]=='transistor':
            update_SysEqn_MOS(i, SetupDict, SimDict, modeldict)
        if DevType[i]=='bipolar':
            update_SysEqn_bipolar(i, SetupDict, SimDict, modeldict)

def build_Jacobian(SetupDict, SimDict, modeldict):
    DeviceCount=SetupDict['DeviceCount']
    DevType=SetupDict['DevType']
    for i in range(DeviceCount):
        if DevType[i]=='transistor':
            build_Jacobian_MOS(i, SetupDict, SimDict, modeldict)
        if DevType[i]=='bipolar':
            build_Jacobian_bipolar(i, SetupDict, SimDict, modeldict)

def build_Jacobian_HB(SetupDict, Simdict, modeldict):
    DeviceCount=SetupDict['DeviceCount']
    DevType=SetupDict['DevType']
    for i in range(DeviceCount):
        if DevType[i]=='transistor':
            build_Jacobian_MOS_HB(i, SetupDict, Simdict, modeldict)

def update_SysEqns_HB(SetupDict, SimDict, modeldict):
    DeviceCount=SetupDict['DeviceCount']
    DevType=SetupDict['DevType']
    TotalHarmonics=SetupDict['TotalHarmonics']
    for i in range(DeviceCount):
        for row in range(TotalHarmonics):
            if DevType[i]=='transistor':
                update_SysEqn_MOS_HB(i, row, SetupDict, SimDict, modeldict)
            if DevType[i]=='bipolar':
                print('Error: Harmonic Balance for Bipolar Transistors not implemented')
                sys.exit(0)


def build_SysEqns_HB(SetupDict, SimDict, modeldict):
    DeviceCount=SetupDict['DeviceCount']
    DevType=SetupDict['DevType']
    TotalHarmonics=SetupDict['TotalHarmonics']
    for i in range(DeviceCount):
        for row in range(TotalHarmonics):
            if DevType[i]=='resistor':
                build_SysEqn_resistor_HB(i, row, SetupDict)
            if DevType[i] == 'oscfilter':
                build_SysEqn_oscfilter_HB(i, row, SetupDict)
            if DevType[i]=='capacitor':
                build_SysEqn_capacitor_HB(i, row, SetupDict, SimDict)
            if DevType[i]=='inductor':
                build_SysEqn_inductor_HB(i, row, SetupDict, SimDict)
            if DevType[i]=='VoltSource':
                build_SysEqn_VSource_HB(i, row, SetupDict)
            if DevType[i]=='CurrentSource':
                build_SysEqn_ISource_HB(i, row, SetupDict)
            if DevType[i]=='transistor':
                build_SysEqn_MOS_HB(i, row, SetupDict, SimDict, modeldict)
            if DevType[i]=='bipolar':
                print('Error: Harmonic Balance for Bipolar Transistors not implemented')
                sys.exit(0)


def build_SysEqn_oscfilter_HB(DeviceNr, row, SetupDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    TotalHarmonics=SetupDict['TotalHarmonics']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    Jacobian_Offset=SetupDict['Jacobian_Offset']
    if row==Jacobian_Offset+1 or row==Jacobian_Offset-1:
        OscFilterValue=DevValue[DeviceNr]
    else:
        OscFilterValue=1e18
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-OscFilterValue
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row]=1
    STA_matrix[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row]=-1
    STA_matrix[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1


def build_SysEqn_resistor_HB(DeviceNr, row, SetupDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    TotalHarmonics=SetupDict['TotalHarmonics']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-DevValue[DeviceNr]
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row]=1
    STA_matrix[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row]=-1
    STA_matrix[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1

def build_SysEqn_capacitor_HB(DeviceNr, row, SetupDict, SimDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    TotalHarmonics=SetupDict['TotalHarmonics']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    omegak=SetupDict['omegak']
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row]=1j*omegak[row]*DevValue[DeviceNr]
    STA_matrix[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row]=-1j*omegak[row]*DevValue[DeviceNr]
    STA_matrix[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1

def build_SysEqn_inductor_HB(DeviceNr, row, SetupDict, SimDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    TotalHarmonics=SetupDict['TotalHarmonics']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    omegak=SetupDict['omegak']
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1j*omegak[row]*DevValue[DeviceNr]
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row]=1
    STA_matrix[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row]=-1
    STA_matrix[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1

def build_SysEqn_VSource_HB(DeviceNr, row, SetupDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    TotalHarmonics=SetupDict['TotalHarmonics']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevLabel=SetupDict['DevLabel']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_rhs=SetupDict['STA_rhs']
    Jacobian_Offset=SetupDict['Jacobian_Offset']
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=0
    if DevNode1[DeviceNr] != '0' :
        STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row]=1
        STA_matrix[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1
    if DevNode2[DeviceNr] != '0' :
        STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row]=-1
        STA_matrix[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1
    if DevLabel[DeviceNr] != 'vinp' and DevLabel[DeviceNr] != 'vinn':
        STA_rhs[(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=getSourceValue(DevValue[DeviceNr],0)*(row==Jacobian_Offset)
    if(DevLabel[DeviceNr] == 'vinp'):
        STA_rhs[(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=.5*((row==Jacobian_Offset+1)+(row==Jacobian_Offset-1))
    if(DevLabel[DeviceNr] == 'vinn'):
        STA_rhs[(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-.5*((row==Jacobian_Offset+1)+(row==Jacobian_Offset-1))
    if(DevLabel[DeviceNr] == 'vin'):
        STA_rhs[(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-.02*((row==Jacobian_Offset+1)+(row==Jacobian_Offset-1))+0.2*(row==Jacobian_Offset)

def build_SysEqn_ISource_HB(DeviceNr, row, SetupDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    TotalHarmonics=SetupDict['TotalHarmonics']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_rhs=SetupDict['STA_rhs']
    Jacobian_Offset=SetupDict['Jacobian_Offset']
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1
    STA_rhs[(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=getSourceValue(DevValue[DeviceNr],0)*(row==Jacobian_Offset)
    if DevNode1[DeviceNr] != '0' and DevNode2[DeviceNr]!='0':
        STA_matrix[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1
        STA_matrix[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1
    elif DevNode2[DeviceNr] != '0' :
        STA_matrix[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1
    elif DevNode1[DeviceNr] != '0' :
        STA_matrix[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1

def build_SysEqn_MOS_HB(DeviceNr, row, SetupDict, SimDict, modeldict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    TotalHarmonics=SetupDict['TotalHarmonics']
    NSamples=SetupDict['NSamples']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevNode3=SetupDict['DevNode3']
    DevModel=SetupDict['DevModel']
    Nodes=SetupDict['Nodes']
    sol=SimDict['sol']
    STA_matrix=SetupDict['STA_matrix']
    STA_nonlinear=SetupDict['STA_nonlinear']
    lambdaT=findParameter(modeldict,DevModel[DeviceNr],'lambdaT')
    VT=findParameter(modeldict,DevModel[DeviceNr],'VT')
    K=findParameter(modeldict,DevModel[DeviceNr],'K')
    Vg=[0 for i in range(TotalHarmonics)]
    Vs=[0 for i in range(TotalHarmonics)]
    Vd=[0 for i in range(TotalHarmonics)]
    TransistorOutputTime=[0 for i in range(NSamples)]
    TransistorOutputFreq=[0 for i in range(TotalHarmonics)]
    for j in range(TotalHarmonics):
        Vg[j]=sol[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+j]
        Vs[j]=sol[Nodes.index(DevNode3[DeviceNr])*TotalHarmonics+j]
        Vd[j]=sol[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+j]
    TransistorOutputTime=TransistorModel(idft(Vg,TotalHarmonics),idft(Vs,TotalHarmonics),idft(Vd,TotalHarmonics),NSamples, K, VT, lambdaT)
    TransistorOutputFreq=dft(TransistorOutputTime,NSamples)
    STA_matrix[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1
    STA_matrix[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1
    STA_matrix[Nodes.index(DevNode3[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1
    STA_nonlinear[(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=TransistorOutputFreq[row]

def update_SysEqn_MOS_HB(DeviceNr, row, SetupDict, SimDict, modeldict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    TotalHarmonics=SetupDict['TotalHarmonics']
    NSamples=SetupDict['NSamples']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevNode3=SetupDict['DevNode3']
    DevModel=SetupDict['DevModel']
    Nodes=SetupDict['Nodes']
    sol=SimDict['sol']
    STA_nonlinear=SetupDict['STA_nonlinear']
    lambdaT=findParameter(modeldict,DevModel[DeviceNr],'lambdaT')
    VT=findParameter(modeldict,DevModel[DeviceNr],'VT')
    K=findParameter(modeldict,DevModel[DeviceNr],'K')
    STA_nonlinear=SetupDict['STA_nonlinear']
    Vg=[0 for i in range(TotalHarmonics)]
    Vs=[0 for i in range(TotalHarmonics)]
    Vd=[0 for i in range(TotalHarmonics)]
    TransistorOutputTime=[0 for i in range(NSamples)]
    TransistorOutputFreq=[0 for i in range(TotalHarmonics)]
    for j in range(TotalHarmonics):
        Vg[j]=sol[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+j]
        Vs[j]=sol[Nodes.index(DevNode3[DeviceNr])*TotalHarmonics+j]
        Vd[j]=sol[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+j]
    TransistorOutputTime=TransistorModel(idft(Vg,TotalHarmonics),idft(Vs,TotalHarmonics),idft(Vd,TotalHarmonics),NSamples, K, VT, lambdaT)
    TransistorOutputFreq=dft(TransistorOutputTime,NSamples)
    STA_nonlinear[(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=TransistorOutputFreq[row]


def build_SysEqn_resistor(DeviceNr, SetupDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_matrix[NumberOfNodes+DeviceNr][NumberOfNodes+DeviceNr]=-DevValue[DeviceNr]
    STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=1
    STA_matrix[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
    STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=-1
    STA_matrix[Nodes.index(DevNode2[DeviceNr])][NumberOfNodes+DeviceNr]=-1

def build_SysEqn_capacitor(deltaT, DeviceNr, SetupDict, SimDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_rhs=SetupDict['STA_rhs']
    method=SetupDict['method']
    sol=SimDict['sol']
    solm1=SimDict['solm1']
    deltaT=SimDict['deltaT']
    STA_matrix[NumberOfNodes+DeviceNr][NumberOfNodes+DeviceNr]=1
    STA_matrix[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
    STA_matrix[Nodes.index(DevNode2[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    if method=='trap':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-2.0*DevValue[DeviceNr]/deltaT
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=2.0*DevValue[DeviceNr]/deltaT
        STA_rhs[NumberOfNodes+DeviceNr]=-2*DevValue[DeviceNr]/deltaT*(sol[Nodes.index(DevNode1[DeviceNr])]-sol[Nodes.index(DevNode2[DeviceNr])])-sol[NumberOfNodes+DeviceNr]
    elif method=='gear2':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-3.0/2.0*DevValue[DeviceNr]/deltaT
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=3.0/2.0*DevValue[DeviceNr]/deltaT
        STA_rhs[NumberOfNodes+DeviceNr]=DevValue[DeviceNr]/deltaT*(-2*(sol[Nodes.index(DevNode1[DeviceNr])]-sol[Nodes.index(DevNode2[DeviceNr])])+1/2*(solm1[Nodes.index(DevNode1[DeviceNr])]-solm1[Nodes.index(DevNode2[DeviceNr])]) )
    elif method=='be':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-1.0*DevValue[DeviceNr]/deltaT
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=1.0*DevValue[DeviceNr]/deltaT
        STA_rhs[NumberOfNodes+DeviceNr]=-DevValue[DeviceNr]/deltaT*(sol[Nodes.index(DevNode1[DeviceNr])]-sol[Nodes.index(DevNode2[DeviceNr])])
    else:
        print('Warning: unknown integration method',method)

def update_SysEqn_capacitor(deltaT, DeviceNr, SetupDict, SimDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_rhs=SetupDict['STA_rhs']
    method=SetupDict['method']
    sol=SimDict['sol']
    solm1=SimDict['solm1']
    if method=='trap':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-2.0*DevValue[DeviceNr]/deltaT
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=2.0*DevValue[DeviceNr]/deltaT
        STA_rhs[NumberOfNodes+DeviceNr]=-2*DevValue[DeviceNr]/deltaT*(sol[Nodes.index(DevNode1[DeviceNr])]-sol[Nodes.index(DevNode2[DeviceNr])])-sol[NumberOfNodes+DeviceNr]
    elif method=='gear2':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-3.0/2.0*DevValue[DeviceNr]/deltaT
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=3.0/2.0*DevValue[DeviceNr]/deltaT
        STA_rhs[NumberOfNodes+DeviceNr]=DevValue[DeviceNr]/deltaT*(-2*(sol[Nodes.index(DevNode1[DeviceNr])]-sol[Nodes.index(DevNode2[DeviceNr])])+1/2*(solm1[Nodes.index(DevNode1[DeviceNr])]-solm1[Nodes.index(DevNode2[DeviceNr])]) )
    elif method=='be':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-1.0*DevValue[DeviceNr]/deltaT
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=1.0*DevValue[DeviceNr]/deltaT
        STA_rhs[NumberOfNodes+DeviceNr]=-DevValue[DeviceNr]/deltaT*(sol[Nodes.index(DevNode1[DeviceNr])]-sol[Nodes.index(DevNode2[DeviceNr])])
    else:
        print('Warning: unknown integration method',method)

def build_SysEqn_inductor(deltaT, DeviceNr, SetupDict, SimDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_rhs=SetupDict['STA_rhs']
    method=SetupDict['method']
    sol=SimDict['sol']
    solm1=SimDict['solm1']
    STA_matrix[NumberOfNodes+DeviceNr][NumberOfNodes+DeviceNr]=1
    STA_matrix[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
    STA_matrix[Nodes.index(DevNode2[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    if method=='trap':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-deltaT/DevValue[DeviceNr]/2
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=deltaT/DevValue[DeviceNr]/2
        STA_rhs[NumberOfNodes+DeviceNr]=sol[NumberOfNodes+DeviceNr]+deltaT*(sol[Nodes.index(DevNode1[DeviceNr])]-sol[Nodes.index(DevNode2[DeviceNr])])/(2*DevValue[DeviceNr])
    elif method=='gear2':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-2/3*deltaT/DevValue[DeviceNr]
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=2/3*deltaT/DevValue[DeviceNr]
        STA_rhs[NumberOfNodes+DeviceNr]=4/3*sol[NumberOfNodes+DeviceNr]-1/3*solm1[NumberOfNodes+DeviceNr]
    elif method=='be':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-deltaT/DevValue[DeviceNr]
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=deltaT/DevValue[DeviceNr]
        STA_rhs[NumberOfNodes+DeviceNr]=sol[NumberOfNodes+DeviceNr]
    else:
        print('Warning: unknown integration method',method)


def update_SysEqn_inductor(deltaT, DeviceNr, SetupDict, SimDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_rhs=SetupDict['STA_rhs']
    method=SetupDict['method']
    sol=SimDict['sol']
    solm1=SimDict['solm1']
    deltaT=SimDict['deltaT']
    if method=='trap':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-deltaT/DevValue[DeviceNr]/2
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=deltaT/DevValue[DeviceNr]/2
        STA_rhs[NumberOfNodes+DeviceNr]=sol[NumberOfNodes+DeviceNr]+deltaT*(sol[Nodes.index(DevNode1[DeviceNr])]-sol[Nodes.index(DevNode2[DeviceNr])])/(2*DevValue[DeviceNr])
    elif method=='gear2':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-2/3*deltaT/DevValue[DeviceNr]
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=2/3*deltaT/DevValue[DeviceNr]
        STA_rhs[NumberOfNodes+DeviceNr]=4/3*sol[NumberOfNodes+DeviceNr]-1/3*solm1[NumberOfNodes+DeviceNr]
    elif method=='be':
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=-deltaT/DevValue[DeviceNr]
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=deltaT/DevValue[DeviceNr]
        STA_rhs[NumberOfNodes+DeviceNr]=sol[NumberOfNodes+DeviceNr]
    else:
        print('Warning: unknown integration method',method)



def build_SysEqn_VSource(DeviceNr, SetupDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_rhs=SetupDict['STA_rhs']
    if DevNode1[DeviceNr] != '0' :
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=1
        STA_matrix[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
    if DevNode2[DeviceNr] != '0' :
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=-1
        STA_matrix[Nodes.index(DevNode2[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    STA_matrix[NumberOfNodes+DeviceNr][NumberOfNodes+DeviceNr]=0
    STA_rhs[NumberOfNodes+DeviceNr]=getSourceValue(DevValue[DeviceNr],0)

def update_SysEqn_VSource(DeviceNr, SimTime, SetupDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    STA_rhs=SetupDict['STA_rhs']
    STA_rhs[NumberOfNodes+DeviceNr]=getSourceValue(DevValue[DeviceNr],SimTime)

def build_SysEqn_ISource(DeviceNr, SetupDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_rhs=SetupDict['STA_rhs']
    if DevNode1[DeviceNr] != '0' :
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=0
        STA_matrix[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=0
    if DevNode2[DeviceNr] != '0' :
        STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=0
        STA_matrix[Nodes.index(DevNode2[DeviceNr])][NumberOfNodes+DeviceNr]=0
    STA_matrix[NumberOfNodes+DeviceNr][NumberOfNodes+DeviceNr]=1
    STA_rhs[NumberOfNodes+DeviceNr]=getSourceValue(DevValue[DeviceNr],0)
    if DevNode1[DeviceNr] != '0' and DevNode2[DeviceNr]!='0':
        STA_matrix[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
        STA_matrix[Nodes.index(DevNode2[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    elif DevNode2[DeviceNr] != '0' :
        STA_matrix[Nodes.index(DevNode2[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    elif DevNode1[DeviceNr] != '0' :
        STA_matrix[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1

def update_SysEqn_ISource(DeviceNr, SimTime, SetupDict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    STA_rhs=SetupDict['STA_rhs']
    STA_rhs[NumberOfNodes+DeviceNr]=getSourceValue(DevValue[DeviceNr],SimTime)

def build_SysEqn_MOS(DeviceNr, SetupDict, SimDict, modeldict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevNode3=SetupDict['DevNode3']
    DevModel=SetupDict['DevModel']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_nonlinear=SetupDict['STA_nonlinear']
    sol=SimDict['sol']
    lambdaT=findParameter(modeldict,DevModel[DeviceNr],'lambdaT')
    VT=findParameter(modeldict,DevModel[DeviceNr],'VT')
    STA_matrix[NumberOfNodes+DeviceNr][NumberOfNodes+DeviceNr]=DevValue[DeviceNr]
    STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=0
    STA_matrix[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
    STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode3[DeviceNr])]=0
    STA_matrix[Nodes.index(DevNode3[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    VD=sol[Nodes.index(DevNode1[DeviceNr])]
    VG=sol[Nodes.index(DevNode2[DeviceNr])]
    VS=sol[Nodes.index(DevNode3[DeviceNr])]
    Vgs=VG-VS
    Vds=VD-VS
    if DevModel[DeviceNr][0]=='p':
        Vds=-Vds
        Vgs=-Vgs
    if Vds < Vgs-VT :
        STA_nonlinear[NumberOfNodes+DeviceNr]=2*((Vgs-VT)*Vds-0.5*Vds**2)
    else :
        STA_nonlinear[NumberOfNodes+DeviceNr]=(Vgs-VT)**2*(1+lambdaT*Vds)

def update_SysEqn_MOS(DeviceNr, SetupDict, SimDict, modeldict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevNode3=SetupDict['DevNode3']
    DevModel=SetupDict['DevModel']
    Nodes=SetupDict['Nodes']
    STA_nonlinear=SetupDict['STA_nonlinear']
    soltemp=SimDict['soltemp']
    lambdaT=findParameter(modeldict,DevModel[DeviceNr],'lambdaT')
    VT=findParameter(modeldict,DevModel[DeviceNr],'VT')
    VD=soltemp[Nodes.index(DevNode1[DeviceNr])]
    VG=soltemp[Nodes.index(DevNode2[DeviceNr])]
    VS=soltemp[Nodes.index(DevNode3[DeviceNr])]
    Vgs=VG-VS
    Vds=VD-VS
    if DevModel[DeviceNr][0]=='p':
        Vds=-Vds
        Vgs=-Vgs
    if Vgs<VT:
        STA_nonlinear[NumberOfNodes+DeviceNr]=1e-5
    elif Vds < Vgs-VT:
        STA_nonlinear[NumberOfNodes+DeviceNr]=2*((Vgs-VT)*Vds-0.5*Vds**2)
    else :
        STA_nonlinear[NumberOfNodes+DeviceNr]=(Vgs-VT)**2*(1+lambdaT*Vds)


def build_Jacobian_MOS(DeviceNr, SetupDict, SimDict, modeldict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevNode3=SetupDict['DevNode3']
    DevModel=SetupDict['DevModel']
    DevValue=SetupDict['DevValue']
    Nodes=SetupDict['Nodes']
    Jacobian=SetupDict['Jacobian']
    soltemp=SimDict['soltemp']
    lambdaT=findParameter(modeldict,DevModel[DeviceNr],'lambdaT')
    VT=findParameter(modeldict,DevModel[DeviceNr],'VT')
    Jacobian[NumberOfNodes+DeviceNr][NumberOfNodes+DeviceNr]=DevValue[DeviceNr]
    VD=soltemp[Nodes.index(DevNode1[DeviceNr])]
    VG=soltemp[Nodes.index(DevNode2[DeviceNr])]
    VS=soltemp[Nodes.index(DevNode3[DeviceNr])]
    Vgs=VG-VS
    Vds=VD-VS
    Vgd=VG-VD
    if DevModel[DeviceNr][0]=='p':
        PFET=-1
        Vgs=-Vgs
        Vds=-Vds
        Vgd=-Vgd
    else:
        PFET=1
    if Vgs<VT :
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=TINY
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=TINY
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode3[DeviceNr])]=-2*TINY
        Jacobian[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
        Jacobian[Nodes.index(DevNode3[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    elif Vds <= Vgs-VT:
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=PFET*2*(Vgd-VT)
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=PFET*2*Vds
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode3[DeviceNr])]=-PFET*2*(Vgs-VT)
        Jacobian[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
        Jacobian[Nodes.index(DevNode3[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    else :
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=PFET*lambdaT*(Vgs-VT)**2
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=PFET*2*(Vgs-VT)*(1+lambdaT*Vds)
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode3[DeviceNr])]=PFET*(-2*(Vgs-VT)*(1+lambdaT*Vds)-lambdaT*(Vgs-VT)**2)
        Jacobian[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
        Jacobian[Nodes.index(DevNode3[DeviceNr])][NumberOfNodes+DeviceNr]=-1


def build_SysEqn_bipolar(DeviceNr, SetupDict, SimDict, modeldict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevNode3=SetupDict['DevNode3']
    DevModel=SetupDict['DevModel']
    Nodes=SetupDict['Nodes']
    STA_matrix=SetupDict['STA_matrix']
    STA_nonlinear=SetupDict['STA_nonlinear']
    soltemp=SimDict['soltemp']
    VEarly=findParameter(modeldict,DevModel[DeviceNr],'Early')
    STA_matrix[NumberOfNodes+DeviceNr][NumberOfNodes+DeviceNr]=DevValue[DeviceNr]
    STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=0
    STA_matrix[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
    STA_matrix[NumberOfNodes+DeviceNr][Nodes.index(DevNode3[DeviceNr])]=0
    STA_matrix[Nodes.index(DevNode3[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    VC=soltemp[Nodes.index(DevNode1[DeviceNr])]
    VB=soltemp[Nodes.index(DevNode2[DeviceNr])]
    VE=soltemp[Nodes.index(DevNode3[DeviceNr])]
    Vbe=VB-VE
    Vce=VC-VE
    if Vbe < 0 :
        STA_nonlinear[NumberOfNodes+DeviceNr]=0
    else :
        STA_nonlinear[NumberOfNodes+DeviceNr]=math.exp(Vbe/Vthermal)*(1+Vce/VEarly)

def build_Jacobian_MOS_HB(DeviceNr, SetupDict, SimDict, modeldict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    TotalHarmonics=SetupDict['TotalHarmonics']
    Jacobian_Offset=SetupDict['Jacobian_Offset']
    NSamples=SetupDict['NSamples']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevNode3=SetupDict['DevNode3']
    DevModel=SetupDict['DevModel']
    Nodes=SetupDict['Nodes']
    Jacobian=SetupDict['Jacobian']
    sol=SimDict['sol']
    lambdaT=findParameter(modeldict,DevModel[DeviceNr],'lambdaT')
    VT=findParameter(modeldict,DevModel[DeviceNr],'VT')
    K=findParameter(modeldict,DevModel[DeviceNr],'K')
    Vg=[0 for i in range(TotalHarmonics)]
    Vs=[0 for i in range(TotalHarmonics)]
    Vd=[0 for i in range(TotalHarmonics)]
    gm=[0 for i in range(TotalHarmonics)]
    for j in range(TotalHarmonics):
        Vg[j]=sol[Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+j]
        Vs[j]=sol[Nodes.index(DevNode3[DeviceNr])*TotalHarmonics+j]
        Vd[j]=sol[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+j]
    gm=TransistorModel_dIdVg(idft(Vg,TotalHarmonics),idft(Vs,TotalHarmonics),idft(Vd,TotalHarmonics),NSamples,K, VT, lambdaT)
    go=TransistorModel_dIdVd(idft(Vg,TotalHarmonics),idft(Vs,TotalHarmonics),idft(Vd,TotalHarmonics),NSamples,K, VT, lambdaT)
    Jlkm=dft(gm,NSamples)
    Jlko=dft(go,NSamples)
    for j in range(TotalHarmonics):
        Jlkm[j]=Jlkm[j]+TINY*(j==Jacobian_Offset)
        Jlko[j]=Jlko[j]+TINY*(j==Jacobian_Offset)
    for row in range(TotalHarmonics):
        Jacobian[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1
        Jacobian[Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=1
        Jacobian[Nodes.index(DevNode3[DeviceNr])*TotalHarmonics+row][(NumberOfNodes+DeviceNr)*TotalHarmonics+row]=-1
        for col in range(TotalHarmonics):
            if(col-row+Jacobian_Offset>=0 and col-row+Jacobian_Offset<TotalHarmonics):
                Jacobian[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode1[DeviceNr])*TotalHarmonics+col]=Jlko[col-row+Jacobian_Offset]
                Jacobian[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode2[DeviceNr])*TotalHarmonics+col]=Jlkm[col-row+Jacobian_Offset]
                Jacobian[(NumberOfNodes+DeviceNr)*TotalHarmonics+row][Nodes.index(DevNode3[DeviceNr])*TotalHarmonics+col]=-Jlkm[col-row+Jacobian_Offset]-Jlko[col-row+Jacobian_Offset]


def update_SysEqn_bipolar(DeviceNr, SetupDict, SimDict, modeldict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevNode3=SetupDict['DevNode3']
    DevModel=SetupDict['DevModel']
    Nodes=SetupDict['Nodes']
    STA_nonlinear=SetupDict['STA_nonlinear']
    soltemp=SimDict['sol']
    VEarly=findParameter(modeldict,DevModel[DeviceNr],'Early')
    VC=soltemp[Nodes.index(DevNode1[DeviceNr])]
    VB=soltemp[Nodes.index(DevNode2[DeviceNr])]
    VE=soltemp[Nodes.index(DevNode3[DeviceNr])]
    Vbe=VB-VE
    Vce=VC-VE
    if Vbe<0:
        STA_nonlinear[NumberOfNodes+DeviceNr]=0
    else :
        STA_nonlinear[NumberOfNodes+DeviceNr]=math.exp(Vbe/Vthermal)*(1+Vce/VEarly)

def build_Jacobian_bipolar(DeviceNr, SetupDict, SimDict, modeldict):
    NumberOfNodes=SetupDict['NumberOfNodes']
    DevValue=SetupDict['DevValue']
    DevNode1=SetupDict['DevNode1']
    DevNode2=SetupDict['DevNode2']
    DevNode3=SetupDict['DevNode3']
    DevModel=SetupDict['DevModel']
    Nodes=SetupDict['Nodes']
    Jacobian=SetupDict['Jacobian']
    soltemp=SimDict['soltemp']
    VEarly=findParameter(modeldict,DevModel[DeviceNr],'Early')
    Jacobian[NumberOfNodes+DeviceNr][NumberOfNodes+DeviceNr]=DevValue[DeviceNr]
    VC=soltemp[Nodes.index(DevNode1[DeviceNr])]
    VB=soltemp[Nodes.index(DevNode2[DeviceNr])]
    VE=soltemp[Nodes.index(DevNode3[DeviceNr])]
    Vbe=VB-VE
    Vce=VC-VE
    Vbc=VB-VC
    if Vbe<=0 :
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=1e-5
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=1e-5
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode3[DeviceNr])]=-1e-5
        Jacobian[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
        Jacobian[Nodes.index(DevNode3[DeviceNr])][NumberOfNodes+DeviceNr]=-1
    else :
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode1[DeviceNr])]=math.exp(Vbe/Vthermal)/VEarly
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode2[DeviceNr])]=math.exp(Vbe/Vthermal)*(1+Vce/VEarly)/Vthermal
        Jacobian[NumberOfNodes+DeviceNr][Nodes.index(DevNode3[DeviceNr])]=(-math.exp(Vbe/Vthermal)/VEarly-math.exp(Vbe/Vthermal)*(1+Vce/VEarly)/Vthermal)
        Jacobian[Nodes.index(DevNode1[DeviceNr])][NumberOfNodes+DeviceNr]=1
        Jacobian[Nodes.index(DevNode3[DeviceNr])][NumberOfNodes+DeviceNr]=-1

def DidResidueConverge(SetupDict, SimDict ):
    NumberOfNodes=SetupDict['NumberOfNodes']
    NumberOfCurrents=SetupDict['NumberOfCurrents']
    STA_matrix=SetupDict['STA_matrix']
    soltemp=SimDict['soltemp']
    f=SimDict['f']
    reltol=SetupDict['reltol']
    iabstol=SetupDict['iabstol']
    ResidueConverged=True
    node=0
    while ResidueConverged and node<NumberOfNodes:
        MaxCurrent=0
        for current in range(NumberOfCurrents):
            MaxCurrent=max(MaxCurrent,abs(STA_matrix[node][NumberOfNodes+current]*(soltemp[NumberOfNodes+current])))
        if f[node] > reltol*MaxCurrent+iabstol:
            ResidueConverged=False
        node=node+1
    return ResidueConverged


def DidUpdateConverge(SetupDict, SimDict ):
    NumberOfNodes=SetupDict['NumberOfNodes']
    Jacobian=SetupDict['Jacobian']
    soltemp=SimDict['soltemp']
    f=SimDict['f']
    reltol=SetupDict['reltol']
    vabstol=SetupDict['vabstol']
    PointLocal=SetupDict['PointLocal']
    GlobalTruncation=SetupDict['GlobalTruncation']
    vkmax=SetupDict['vkmax']
    SolutionCorrection=numpy.matmul(numpy.linalg.inv(Jacobian),f)
    UpdateConverged=True
    if PointLocal:
        for node in range(NumberOfNodes):
            vkmax=max(abs(soltemp[node]),abs(soltemp[node]-SolutionCorrection[node]))
            if abs(SolutionCorrection[node])>vkmax*reltol+vabstol:
                UpdateConverged=False
    elif GlobalTruncation:
        for node in range(NumberOfNodes):
            if abs(SolutionCorrection[node])>vkmax*reltol+vabstol:
                UpdateConverged=False
    else:
        print('Error: Unknown truncation error')
        sys.exit()
    return UpdateConverged

def DidLTEConverge(SetupDict, SimDict, iteration, LTEIter, NewtonConverged, timeVector, SimTime, SolutionCorrection):
    NumberOfNodes=SetupDict['NumberOfNodes']
    sol=SimDict['sol']
    soltemp=SimDict['soltemp']
    solm1=SimDict['solm1']
    solm2=SimDict['solm2']
    PointLocal=SetupDict['PointLocal']
    GlobalTruncation=SetupDict['GlobalTruncation']
    lteratio=SetupDict['lteratio']
    vkmax=SetupDict['vkmax']
    reltol=SetupDict['reltol']
    vabstol=SetupDict['vabstol']
    LTEConverged=True
    MaxLTERatio=0
    PredMatrix=[[0 for i in range(2)] for j in range(2)]
    Predrhs=[0 for i in range(2)]
    if iteration>200 and NewtonConverged:
        LTEIter=LTEIter+1
        for i in range(NumberOfNodes):
            tau1=(timeVector[iteration-2]-timeVector[iteration-3])
            tau2=(timeVector[iteration-1]-timeVector[iteration-3])
            PredMatrix[0][0]=tau2
            PredMatrix[0][1]=tau2*tau2
            PredMatrix[1][0]=tau1
            PredMatrix[1][1]=tau1*tau1
            Predrhs[0]=sol[i]-solm2[i]
            Predrhs[1]=solm1[i]-solm2[i]
            Predsol=numpy.matmul(numpy.linalg.inv(PredMatrix),Predrhs)
            vpred=solm2[i]+Predsol[0]*(SimTime-timeVector[iteration-3])+Predsol[1]*(SimTime-timeVector[iteration-3])*(SimTime-timeVector[iteration-3])
            if PointLocal:
                for node in range(NumberOfNodes):
                    vkmax=max(abs(soltemp[node]),abs(soltemp[node]-SolutionCorrection[node]))

                    if abs(vpred-soltemp[i])> lteratio*(vkmax*reltol+vabstol):
                        LTEConverged=False
                    else:
                        MaxLTERatio=max(abs(vpred-soltemp[i])/(vkmax*reltol+vabstol),MaxLTERatio)
            elif GlobalTruncation:
                for node in range(NumberOfNodes):
                    if abs(vpred-soltemp[i])> lteratio*(vkmax*reltol+vabstol):
                        LTEConverged=False

            else:
                print('Error: Unknown truncation error')
                sys.exit()
    return LTEConverged, MaxLTERatio


def UpdateTimeStep(SetupDict, SimDict, LTEConverged, NewtonConverged, val, iteration, NewtonIter, MaxLTERatio, timeVector, SimTime):
    MatrixSize=SetupDict['MatrixSize']
    FixedTimeStep=SetupDict['FixedTimeStep']
    MaxTimeStep=SetupDict['MaxTimeStep']
    soltemp=SimDict['soltemp']
    sol=SimDict['sol']
    solm1=SimDict['solm1']
    solm2=SimDict['solm2']
    deltaT=SimDict['deltaT']
    ThreeLevelStep=SimDict['ThreeLevelStep']
    Converged=NewtonConverged and LTEConverged
    if Converged:
        if not NewtonConverged:
            print('Some trouble converging, skipping')
            NewtonConverged=True
        for i in range(MatrixSize):
            solm2[i]=solm1[i]
            solm1[i]=sol[i]
        if iteration > -1:
            for j in range(MatrixSize):
                val[j][iteration]=soltemp[j]
        timeVector[iteration]=SimTime
        iteration=iteration+1
        if not FixedTimeStep:
            if ThreeLevelStep:
                if 0.9<MaxLTERatio<1.0:
                    deltaT=deltaT/1.1
                else:
                    if MaxLTERatio<0.1:
                        deltaT=1.01*deltaT
            else:
                deltaT=1.001*deltaT
            deltaT=min(deltaT,MaxTimeStep)
    else:
        if FixedTimeStep:
            if iteration>100:
                if not NewtonConverged:
                    print('Newton failed to converge',NewtonIter)
                if not LTEConverged:
                    print('LTE failed to converge',NewtonIter)
                sys.exit()
        else:
            SimTime=max(SimTime-deltaT,0)
            deltaT=deltaT/1.1
    return deltaT, iteration, SimTime, Converged

def dft(Samples,N):
    sol=[0+0j for i in range(N+1) ]
    y=fft(Samples)
    for i in range(int(N/2)):
        sol[i]=y[int(N/2)+i]/N
    for i in range(int(N/2)+1):
        sol[int(N/2)+i]=y[i]/N
        # DC??
    return sol

def idft(Vk, N):
    y=[0 for i in range(N-1)]
    for i in range(int(N/2)+1):
        y[i]=Vk[int(N/2)+i]*(N-1)
    for i in range(int(N/2)-1):
        y[int(N/2)+i+1]=Vk[i+1]*(N-1)
    return ifft(y,N-1)

def TransistorModel(Vg, Vs, Vd, N, K, VT, lambdaT):
    Id=[0 for i in range(N)]
    for i in range(N):
        Vgs=Vg[i]-Vs[i]
        Vds=Vd[i]-Vs[i]
        if Vgs<VT:
            Id[i]=K*0
        elif Vds < Vgs-VT :
            Id[i]=2*K*((Vgs-VT)*Vds-0.5*Vds**2)
        else:
            Id[i]=K*(Vgs-VT)**2*(1+lambdaT*Vds)
    return Id
def TransistorModel_dIdVg(Vg, Vs, Vd, N, K, VT, lambdaT):
    gmm=[0 for i in range(N)]
    for i in range(N):
        Vgs=Vg[i]-Vs[i]
        Vds=Vd[i]-Vs[i]
        if Vgs<VT:
            gmm[i]=K*0
        elif Vds < Vgs-VT :
            gmm[i]=K*2*Vds
        else:
            gmm[i]=2*K*(Vgs-VT)*(1+lambdaT*Vds)
    return gmm
def TransistorModel_dIdVd(Vg, Vs, Vd, N, K, VT, lambdaT):
    goo=[0 for i in range(N)]
    for i in range(N):
        Vgs=Vg[i]-Vs[i]
        Vds=Vd[i]-Vs[i]
        Vgd=Vg[i]-Vd[i]
        if Vgs<VT:
            goo[i]=K*0
        elif Vds < Vgs-VT :
            goo[i]=2*K*(Vgd-VT)
        else:
            goo[i]=K*lambdaT*(Vgs-VT)**2
    return goo
def TimeDerivative(inp,wk,N):
    deriv=[0+0j for i in range(N)]
    for i in range(N):
        deriv[i]=inp[i]*wk[i]*1j
    return deriv