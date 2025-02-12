
#***********************************************************************
#DD package, data collection and analysis of 2D electronic spectra
#Copyright (C) 2016, 2017  Jan Alster (Charles Univesity, Prague)
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>
#***********************************************************************

# -*- coding: utf-8 -*-

try:
    from . import libPyPopot as Popot   
except SystemError:
    import libPyPopot as Popot   
    
import time
Popot.seed((int)((time.time()%100)*10000))
import numpy as np
from PyQt5 import  QtCore

try:
    from .helpers import icmCoef # w[1/cm] = w[rad/fs] * icmCoef
except SystemError:
    from helpers import icmCoef

import time

def fit_gaussian(x, y, maxfev=100, weights=1):
    dimension = 3
    
    #def lbound(index):
        #return -30
    #def ubound(index):
        #return 30

    lbound=lambda index: [x[0], 0.9*min(y), 1e-10][index]
    ubound=lambda index: [x[-1], 1.1*max(y), abs(x[-1]-x[0])][index]
    
    amplitude=y.max()-y.min()
    
    def stop(fitness, epoch):
        return bool(fitness <= 0.0001*amplitude**2) | (epoch >= maxfev)
        #return (fitness <= 100) | (epoch >= maxfev)
        
    def evaluate(param):
        gauss=param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2)
        return np.sum(((gauss-y)*weights)**2)/(len(gauss)-1.)
        
    #Initialize libPyPopot
    algo = Popot.Stochastic_PSO_2006(dimension, lbound, ubound, stop, evaluate)

    #Run until the stop criterion is met
    algo.run(0)
    #return algo.bestParticule(), algo.getEpoch(), algo.bestFitness()
    param=algo.bestParticule()
    return param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2), tuple(param), algo.getEpoch(), algo.bestFitness()

def fit_gaussian0(x, y, maxfev=100, weights=1):
    dimension = 4
    
    lbound=lambda index: [x[0], min(y), 1e-10, min(y)][index]
    ubound=lambda index: [x[-1], max(y), abs(x[-1]-x[0]), max(y)][index]
    
    def stop(fitness, epoch):
        return (fitness <= 100) | (epoch >= maxfev)
        
    def evaluate(param):
        gauss=param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2)+param[3]
        return np.sum(((gauss-y)*weights)**2)/(len(gauss)-1.)
        
    #Initialize libPyPopot
    algo = Popot.Stochastic_PSO_2006(dimension, lbound, ubound, stop, evaluate)

    #Run until the stop criterion is met
    algo.run(0)
    #return algo.bestParticule(), algo.getEpoch(), algo.bestFitness()
    param=algo.bestParticule()
    return param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2)+param[3], tuple(param), algo.getEpoch(), algo.bestFitness()


def fit_lorentzian(x, y, maxfev=100, weights=1):
    dimension = 3
    
    #def lbound(index):
        #return -30
    #def ubound(index):
        #return 30

    lbound=lambda index: [x[0], min(y), 1e-10][index]
    ubound=lambda index: [x[-1], max(y), abs(x[-1]-x[0])][index]
    
    def stop(fitness, epoch):
        return (fitness <= 100) | (epoch >= maxfev)
        
    def evaluate(param):
        lorentz=param[1]/(((x-param[0])/(0.5*param[2]))**2+1)
        return np.sum(((lorentz-y)*weights)**2)/(len(lorentz)-1.)
        
    #Initialize libPyPopot
    algo = Popot.Stochastic_PSO_2006(dimension, lbound, ubound, stop, evaluate)

    #Run until the stop criterion is met
    algo.run(0)
    #return algo.bestParticule(), algo.getEpoch(), algo.bestFitness()
    param=algo.bestParticule()
    return param[1]/(((x-param[0])/(0.5*param[2]))**2+1), tuple(param), algo.getEpoch(), algo.bestFitness()    

def fit_lorentzian0(x, y, maxfev=100, weights=1):
    dimension = 4
    
    lbound=lambda index: [x[0], min(y), 1e-10, min(y)][index]
    ubound=lambda index: [x[-1], max(y), abs(x[-1]-x[0]), max(y)][index]
    
    def stop(fitness, epoch):
        return (fitness <= 100) | (epoch >= maxfev)
        
    def evaluate(param):
        lorentz=param[1]/(((x-param[0])/(0.5*param[2]))**2+1)+param[3]
        return np.sum(((lorentz-y)*weights)**2)/(len(lorentz)-1.)
        
    #Initialize libPyPopot
    algo = Popot.Stochastic_PSO_2006(dimension, lbound, ubound, stop, evaluate)

    #Run until the stop criterion is met
    algo.run(0)
    #return algo.bestParticule(), algo.getEpoch(), algo.bestFitness()
    param=algo.bestParticule()
    return param[1]/(((x-param[0])/(0.5*param[2]))**2+1)+param[3], tuple(param), algo.getEpoch(), algo.bestFitness()    


#models for GA - we can allow
# only populations
# also frequencices with +,- with same amplitude and lifetime
# also frequencices with +,- with different amplitude and same lifetime
# also frequencices with +,- with same amplitude and different lifetime
# also frequencices with +,- with different amplitude and lifetime
# fixing of some parameters

#we will get [(lifetime, fixed, frequency type, frequency, fixed)]
#[(4.0, False, u'\xd8', 0.0, False), (6.0, False, u'\xd8', 0.0, False), (3.0, True, u'\xb1', 0.0, False)]

#this will be slow...
def C(rates, frequencies, freqTypes, t2):
    temp=[]
    for r, f, ft in zip(rates, frequencies, freqTypes):
        if ft=='Ø':   
            temp.append(-r)
        #TODO: check signs here
        if ft=='+' or ft=='±':
            temp.append(-r+f*1j)
        if ft=='-' or ft=='±':
            temp.append(-r-f*1j)
        pass
    return np.exp(np.outer(t2, temp))

#take values from params and update rates and w which are not fixed
def update(params, rates, w, freqTypes, rFix, wFix):
    for i, freqType in enumerate(freqTypes):
        if not rFix[i]:
            rates[i]=params.pop(0)
        if freqType!='Ø' and not wFix[i]:
            w[i]=params.pop(0)
        pass
    pass

#scan input list and create list of params that are not fixed
def divide(inputList):
    params=[]
    rates=[]
    w=[]
    rFix=[]
    wFix=[]
    types=[]
    
    for lifetime, lFix, freqType, freq, fFix in inputList:
        rFix.append(lFix)
        rates.append(1./lifetime)
        w.append(freq/icmCoef)
        wFix.append(fFix)
        types.append(freqType)
        if not lFix:
            params.append(rates[-1])
        if freqType!='Ø' and not fFix:
            params.append(w[-1])
        pass

    #TODO: get bounds for these params
    return params, rates, w, rFix, wFix, types

def combine(params, inputList):
    for i in range(len(inputList)):
        lifetime, lFix, freqType, freq, fFix=inputList[i]
        if not lFix:
            lifetime=1./params.pop(0)
        if freqType!='Ø' and not fFix:
            freq=params.pop(0)*icmCoef
        inputList[i]=(lifetime, lFix, freqType, freq, fFix)
        pass
    pass
    
def model(data, t2, inputList):
    params, rates, w, rFix, wFix, types=divide(inputList)
    c=C(rates, w, types, t2)
    cp=np.linalg.pinv(c)
    return np.tensordot(np.dot(c, cp), data, 1)

def DASandRes(data, t2, inputList):
    params, rates, w, rFix, wFix, types=divide(inputList)
    c=C(rates, w, types, t2)
    cp=np.linalg.pinv(c)
    return np.tensordot(cp, data, 1), np.tensordot(np.eye(c.shape[0])-np.dot(c, cp), data, 1)
    


def setupGA(data, t2, inputList, maxfev=100):
    params, rates, w, rFix, wFix, types=divide(inputList)

    if len(params)==0:
        import warnings
        warnings.warn('Nothing to fit')
        return None
    
    dimension=len(params)

    #popot has a quirk that it reevaluates the best position of each particle, which is not needed here and cannot be switched off from python
    #storing the values will speed things up twice, hopefully it will not cause other problems
    store={}
    
    def evaluate(p):
        key=tuple(p)
        if key in store:
            return store[key]
        
        update(p, rates, w, types, rFix, wFix)
        c=C(rates, w, types, t2)
        cp=np.linalg.pinv(c)
        R=np.tensordot(np.eye(c.shape[0])-np.dot(c, cp), data, 1) #(I-C(C+))Y == Y - C(C+)Y == Y - M
        res=np.sum(np.abs(R), dtype=np.float64)
        store[key]=res
        return res

    def stop(fitness, epoch):
        return (fitness <= 100) | (epoch >= maxfev)

    #bounds for lifetimes and frequencies
    lbound=lambda index: 1./(10*t2[-1])
    ubound=lambda index: 10./(t2[1]-t2[0])

    print("setupGA: starting params", params)
    print("bounds", lbound(0), ubound(0))
    algo = Popot.Stochastic_PSO_2006(dimension, lbound, ubound, stop, evaluate) 
    print("setting starting particule")
    algo.changeParticule(params)

    print(t2)
    return algo #so that user can run it several times, or start/stop at will



class SwarmThread(QtCore.QThread):
    report=QtCore.pyqtSignal(list, int, float)
    finalReport=QtCore.pyqtSignal(list, int, float)
    def __init__(self, parent, data, t2, inputList, reportEachNth=10, maxSteps=1000):
        super(SwarmThread, self).__init__(parent)
        self._data=data
        self._t2=t2
        self._inputList=inputList
        self.reportEachNth=reportEachNth
        self._maxSteps=maxSteps
        pass

    def run(self):
        algo=setupGA(self._data, self._t2, self._inputList)
        counter=0
        while not self.isInterruptionRequested() and counter<self._maxSteps:
            #next step
            algo.step()
            
            #emit results (or emit results every N steps)
            counter+=1
            if (counter % self.reportEachNth) == 0:
                self.report.emit(algo.bestParticule(), counter, algo.bestFitness())
            pass
        self.finalReport.emit(algo.bestParticule(), counter, algo.bestFitness())
        #TODO: report on fitness and counter too
    pass
