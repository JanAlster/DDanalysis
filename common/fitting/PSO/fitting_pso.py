
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

#globalni
import time
import numpy as np

from PyQt5 import  QtCore, QtWidgets

if __name__=="__main__":
    import sys
    sys.path.append("../../..")

#lokalni
try:
    from .pso import PSO
except (ModuleNotFoundError, SystemError):
    from pso import PSO

try:
    from . import icmCoef # w[1/cm] = w[rad/fs] * icmCoef
except (ModuleNotFoundError, SystemError):
    from common import icmCoef
except ImportError:
    print("fitting_pso.py: blabol")
    from common import icmCoef


#obecne fitovani
def fit(func, target, lbound, ubound, start=None, maxfev=100, weights=1, args=tuple(), kw={}, **PSOkw):
    #TODO: try to pop omega, phip, etc from kw
    #print("fitting_pso::fit", *args, **kw)
    target = np.asarray(target)
    N = (max(target.size-1, 1))  #prevent divide-by-zero in case of just one point targer
    if np.iscomplex(target).any():
        def evaluate(params):
            return np.sum(np.abs((func(params, *args, **kw)-target)*weights)**2)/N 
    else:
        def evaluate(params):
            return np.sum(((func(params, *args, **kw)-target)*weights)**2)/N
        
    swarm = PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8, **PSOkw)
    
    if start is not None:
        swarm.set_particle(0, start)
        print("starting particle value", swarm.obj(start))
        #Todo: possibly a bug, if particle 0 happens to be the best in the swarm, it could lead to problems (at minimum we are replacing a better solution which might not be scanned around now)
        # but likely not an issue
    
    return swarm.run()

#jednoduche funkce
#TODO: better bounds - allow more for amplitude
def fit_gaussian(x, y, maxfev=300, weights=1):
    """
    Fit zero based gaussian
    """
    lbound=[x[0], min(y), 1e-10]
    ubound=[x[-1], max(y), abs(x[-1]-x[0])]
    
    def evaluate(param):
        gauss=param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2)
        return np.sum(((gauss-y)*weights)**2)/(len(gauss)-1.)

    swarm=PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8)
        
    param, ff = swarm.run()
    return param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2), tuple(param), swarm.it, ff

def fit_2_gaussian(x, y, maxfev=300, weights=1):
    """
    Fit zero based gaussian
    """
    lbound = [x[0], min(y), 1e-10, x[0], min(y), 1e-10]
    ubound = [x[-1], max(y), abs(x[-1] - x[0]), x[-1], max(y), abs(x[-1] - x[0])]

    def evaluate(param):
        gauss = param[1] * np.exp(-2.7725887222397811 * ((x - param[0]) / param[2]) ** 2)
        gauss2 = param[4] * np.exp(-2.7725887222397811 * ((x - param[3]) / param[5]) ** 2)
        return np.sum(((gauss+gauss2 - y) * weights) ** 2) / (len(gauss) - 1.)

    swarm = PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8)

    param, ff = swarm.run()
    return (param[1] * np.exp(-2.7725887222397811 * ((x - param[0]) / param[2]) ** 2)
           + param[4] * np.exp(-2.7725887222397811 * ((x - param[3]) / param[5]) ** 2),
            tuple(param), swarm.it, ff)


def fit_gaussian0(x, y, maxfev=300, weights=1):
    """
    Fit gaussian with constant background
    """
    lbound=[x[0], min(y), 1e-10, min(y)]
    ubound=[x[-1], max(y), abs(x[-1]-x[0]), max(y)]
    
    def evaluate(param):
        gauss=param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2)+param[3]
        return np.sum(((gauss-y)*weights)**2)/(len(gauss)-1.)

    swarm=PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8)
        
    param, ff = swarm.run()
    return param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2)+param[3], tuple(param), swarm.it, ff


def fit_lorentzian(x, y, maxfev=100, weights=1):
    lbound=[x[0], min(y), 1e-10]
    ubound=[x[-1], max(y), abs(x[-1]-x[0])]
    
    def evaluate(param):
        lorentz=param[1]/(((x-param[0])/(0.5*param[2]))**2+1)
        return np.sum(((lorentz-y)*weights)**2)/(len(lorentz)-1.)

    swarm=PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8)
        
    param, ff = swarm.run()
    return param[1]/(((x-param[0])/(0.5*param[2]))**2+1), tuple(param), swarm.it, ff

def fit_lorentzian0(x, y, maxfev=100, weights=1):
    lbound=[x[0], min(y), 1e-10, min(y)]
    ubound=[x[-1], max(y), abs(x[-1]-x[0]), max(y)]
    
    def evaluate(param):
        lorentz=param[1]/(((x-param[0])/(0.5*param[2]))**2+1)+param[3]
        return np.sum(((lorentz-y)*weights)**2)/(len(lorentz)-1.)

    swarm=PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8)
        
    param, ff = swarm.run()
    return param[1]/(((x-param[0])/(0.5*param[2]))**2+1)+param[3], tuple(param), swarm.it, ff


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
#we might do this faster considering that freqTypes do not change
def _innerC(rates, frequencies, freqTypes):
    temp=[]
    for r, f, ft in zip(rates, frequencies, freqTypes):
        if ft=='Ø':   
            temp.append(-r)
        #Todo: check signs here
        """
        according to my Lund lecture from 11.6.2014
        for model with expression OM_N(w1, w3) * exp(-1j*wN*t2-kN*t2)
        we have positive wN for Feyman digram with larger left side in waiting time part
        (and I hope that that lecture has it correct)
        """
        if ft=='+' or ft=='±':
            temp.append(-r-f*1j)
        if ft=='-' or ft=='±':
            temp.append(-r+f*1j)
        pass
    return temp

#for testing
innerC = _innerC    

def C(rates, frequencies, freqTypes, t2):
    temp=_innerC(rates, frequencies, freqTypes)
    # ~ print(temp)
    #this should include Heaviside(t2)
    # ~ print(np.ma.masked_less(t2, 0))
    # ~ print(np.outer(np.ma.masked_less(t2, 0), temp))
    return np.exp(np.ma.outer(np.ma.masked_less(t2, 0), temp)).filled(0)

#-----------------------------------------------------------------------
#fukce pro IRF (see GA for details)

#exp_conv, Nexp_conv
#see GA/fitovanic/fitovani.c

#pocita M=sum_i(conv(p[2+N+i]*exp(-p[2+i]*t)*theta(t), norm*exp(-(t+p[0])**2/(2*d**2))))
#kde
#norm=1/(d*(2*pi)**0.5)
#d=p[1]/(2*(2*ln(2))**0.5)

#i.e. sum of convolutions of exp decay with gaussian IRF

#seems like numpy does not have erfc and we need it (and do not want to depend on scipy)
# and we need the complex version
#this is actually bad, python itself does not have a complex version of erf :(
# there is one in scipy, but I do not want to depend on scipy it is too large 

#we could use mpmath (not happy about new dependence)
#  and is arbitrary precision, not native float, could be slow?
#  scipy implements erf in fortran
"""
If you use mpmath in your research, please cite it! In BibTeX format, the following entry can be used:

@manual{mpmath,
  key     = {mpmath},
  author  = {Fredrik Johansson and others},
  title   = {mpmath: a {P}ython library for arbitrary-precision floating-point arithmetic (version 0.18)},
  note    = {{\tt http://mpmath.org/}},
  month   = {December},
  year    = {2013},
}

This might render as:

    Fredrik Johansson and others. mpmath: a Python library for arbitrary-precision floating-point arithmetic (version 0.18), December 2013. http://mpmath.org/.

Preferably, change the version number to the version you used, and use the date of its release.
"""

"""
alternatively, I can do python binding for this
https://jugit.fz-juelich.de/mlz/libcerf
https://packages.debian.org/jessie/libs/libcerf1

it is another dependence but it could be fast and should be fairly straightforfard?

but it does not have binary package for Windows and I want to avoid compilation
"""

#seems like the best approach is actually to depend on scipy :(
import math
import cmath

def _RADAs(x):
    val=1./x
    inc=val
    n=0
    dx=-0.5*x**-2
    while n<=100:
        n+=1
        inc=(2*n-1)*dx*inc
        val+=inc
        if (np.abs(inc)<1e-15).all():
            break
    return val

#these are included just for reference from GA, do not use    
def exp_conv(x, p, res=None):
    flip=False
    if x[0]>x[-1]:
        flip=True
        x=x[::-1]
        
    d=p[1]*0.42466090014400952136
    p2d=p[2]*d
    a=(p2d-(x+p[0])/d)*0.7071067811865475244
    n=len(x)
    if res==None:
        res=np.zeros_like(x)
        
    A=(p2d-np.array([6, -10])/0.7071067811865475244)*d-p[0]
    i1,i2 = x.searchsorted(A)

    if i2>i1 or i2<len(x):
        res[i1:]=p[3]*np.exp(p[2]*((d*p2d*0.5-p[0])-x[i1:])) #nebo jen pro x[i1:]

    if i1>0:
        res[:i1]=(p[3]*0.5*0.56418958354775628695)*np.exp(-((x[:i1]+p[0])/(2**0.5*d))**2)*_RADAs(a[:i1])
    if i2>i1:
        erfc = np.array([math.erfc(it) for it in a[i1:i2]])
        res[i1:i2]*=0.5*erfc
        
    return res[::-1] if flip else res
    
#these are included just for reference from GA, do not use    
def Nexp_conv(x, p, res=None):
    x=np.array(x, dtype=np.float64)
    #TODO: this will work only for ndarray
    flip=False
    if x[0]>x[-1]:
        flip=True
        x=x[::-1]

    n=len(x)
    if res==None:
        res=np.zeros_like(x)
    else:
        res[:]=0

    d=p[1]*0.42466090014400952136

    N=(len(p)-2)/2
        
    x0=x+p[0]
    x0d=x0/d
    for j in range(N):
        p2d=p[2+j]*d
        a=(p2d-x0d)*0.7071067811865475244

        A=(p2d-np.array([6, -10])/0.7071067811865475244)
        i1,i2 = x0d.searchsorted(A)

        if i1>0:
            res[:i1]+=(p[2+N+j]*0.5*0.56418958354775628695)*np.exp(-0.5*x0d[:i1]**2)*_RADAs(a[:i1])

        if i2>i1:
            erfc=np.array([math.erfc(it) for it in a[i1:i2]])
            res[i1:i2]+=0.5*erfc*p[2+N+j]*np.exp(p2d*((p2d*0.5)-x0d[i1:i2]))

        if i2<len(x):
            res[i2:]+=p[2+N+j]*np.exp(p2d*((p2d*0.5)-x0d[i2:]))

        pass

    return res[::-1] if flip else res

#however here we need equivalent of C, i.e. not sum but array of such convolutions
#TODO: not sure if convolution with complex exp is the same as for real ones in GA
#  looking on wolfram|alpha evaluation of "integrate from 0 to infty exp(-ax^2+b x+c)dx" it does not matter if coeficient "b" is complex
# but condition handling possibly needs to be adjusted

#NOTE: originaly I have wanted to have a possibility to not use IRF if scipy is not available (it is pain to have a dependency just for one function)
# but the effect on data is so big that it makes IRF mandatory
from scipy.special import erfc

def CIRF(rates, frequencies, freqTypes, t2, IRFOffset, IRFFWHM, res=None):
    # ~ print("\nCIRF start", rates, frequencies, freqTypes, t2, IRFOffset, IRFFWHM)
    temp = _innerC(rates, frequencies, freqTypes)
    temp = np.asarray(temp)
    # ~ print("\t", temp)
    # ~ print()
    
    x = np.asarray(t2, dtype=np.float64)

    #this is likely needed because we use searchsorted
    flip=False
    if x[0]>x[-1]:
        flip=True
        x=x[::-1]

    n = len(x)
    N = len(temp)
    if res is None:
        res = np.zeros((N, n), dtype=temp.dtype)
    else:
        res[:,:]=0

    d = IRFFWHM*0.42466090014400952136

    x0d = (x+IRFOffset)/d
    
    for j in range(N):
        # ~ p2d = -temp[j].conjugate()*d #done: this calculation assumes positive real(rate) but _innerC already includes the minus; using minus here will change sign of frequency component too (thats why we have conjugate) -> redo to sign in calculations from GA to work with innerC better
        p2d = -temp[j]*d #like this it should work (see simulation at the end of the file)
        a = (p2d-x0d)*0.7071067811865475244

        A=(p2d-np.array([6, -10])/0.7071067811865475244)
        i1,i2 = x0d.searchsorted(A)

        

        if i1>0:
            res[j, :i1]+=(0.5*0.56418958354775628695)*np.exp(-0.5*x0d[:i1]**2)*_RADAs(a[:i1])

        if i2>i1:
            # ~ erfc=np.array([erfc(it) for it in a[i1:i2]])
            # ~ print("erfc", erfc(a[i1:i2]))
            # ~ print("a", a[i1:i2])
            # ~ print("p2d", p2d)
            # ~ print("x0d", x0d[i1:i2])
            res[j, i1:i2]+=0.5*erfc(a[i1:i2])*np.exp(p2d*((p2d*0.5)-x0d[i1:i2]))

        if i2<len(x):
            res[j, i2:]+=np.exp(p2d*((p2d*0.5)-x0d[i2:]))

        # ~ print("CIRF", j, i1, i2, res[j])
        pass

    return res[:,::-1].T if flip else res.T #we need .T to be compatible with C TODO: probably could be handled better
#-----------------------------------------------------------------------

#take values from params and update rates and w which are not fixed
def update(params, rates, w, freqTypes, rFix, wFix):
    pos=0
    for i, freqType in enumerate(freqTypes):
        if not rFix[i]:
            rates[i]=params[pos]
            pos+=1
        if freqType!='Ø' and not wFix[i]:
            w[i]=params[pos]
            pos+=1
        pass
    pass

def updateIRF(params, rates, w, freqTypes, rFix, wFix, IRF):
    pos=0
    for i, freqType in enumerate(freqTypes):
        if not rFix[i]:
            rates[i]=params[pos]
            pos+=1
        if freqType!='Ø' and not wFix[i]:
            w[i]=params[pos]
            pos+=1
        pass
    if not IRF[1]:
        IRF[0] = params[pos]
        pos += 1
    if not IRF[3]:
        IRF[2] = params[pos]
        pos += 1
    pass

#~ def update(params, rates, w, freqTypes, rFix, wFix):
    #~ for i, freqType in enumerate(freqTypes):
        #~ if not rFix[i]:
            #~ rates[i]=params.pop(0)
        #~ if freqType!='Ø' and not wFix[i]:
            #~ w[i]=params.pop(0)
        #~ pass
    #~ pass


#scan input list and create list of params that are not fixed
def divide(inputList):
    params=[]
    rates=[]
    w=[]
    rFix=[]
    wFix=[]
    types=[]
    paramsTypes = ""
    # ~ print("inputList", inputList)
    for lifetime, lFix, freqType, freq, fFix in inputList:
        rFix.append(lFix)
        rates.append(1./lifetime)
        w.append(freq/icmCoef)
        wFix.append(fFix)
        types.append(freqType)
        if not lFix:
            params.append(rates[-1])
            paramsTypes += "r"
        if freqType!='Ø' and not fFix:
            params.append(w[-1])
            paramsTypes += "w"
        pass

    #bounds for these params cannot be determined here, because those depend on t2
    return params, rates, w, rFix, wFix, types, paramsTypes

def divideIRF(inputList, inputIRF):
    params=[]
    rates=[]
    w=[]
    rFix=[]
    wFix=[]
    types=[]
    paramsTypes = ""
    IRF = [it for it in inputIRF] #copy 
    
    for lifetime, lFix, freqType, freq, fFix in inputList:
        rFix.append(lFix)
        rates.append(1./lifetime)
        w.append(freq/icmCoef)
        wFix.append(fFix)
        types.append(freqType)
        if not lFix:
            params.append(rates[-1])
            paramsTypes += "r"
        if freqType!='Ø' and not fFix:
            params.append(w[-1])
            paramsTypes += "w"
        pass
    
    if not IRF[1]:
        params.append(IRF[0])
        paramsTypes += "W" #FWHM
    
    if not IRF[3]:
        params.append(IRF[2])
        paramsTypes += "O" #Offset
    

    #bounds for these params cannot be determined here, because those depend on t2
    return params, rates, w, rFix, wFix, types, paramsTypes, IRF


def combine(params, inputList):
    pos=0
    for i in range(len(inputList)):
        lifetime, lFix, freqType, freq, fFix=inputList[i]
        if not lFix:
            lifetime=1./params[pos]
            pos+=1
        if freqType!='Ø' and not fFix:
            freq=params[pos]*icmCoef
            pos+=1
        inputList[i]=(lifetime, lFix, freqType, freq, fFix)
        pass
    pass

def combineIRF(params, inputList, IRF):
    pos=0
    for i in range(len(inputList)):
        lifetime, lFix, freqType, freq, fFix=inputList[i]
        if not lFix:
            lifetime=1./params[pos]
            pos+=1
        if freqType!='Ø' and not fFix:
            freq=params[pos]*icmCoef
            pos+=1
        inputList[i]=(lifetime, lFix, freqType, freq, fFix)
        pass

    if not IRF[1]:
        IRF[0] = params[pos] #FWHM
        pos += 1 
    
    if not IRF[3]:
        IRF[2] = params[pos]
        pos += 1
 
    pass


def combineCopy(params, inputList):
    pos=0
    R = []
    for i in range(len(inputList)):
        lifetime, lFix, freqType, freq, fFix=inputList[i]
        if not lFix:
            lifetime=1./params[pos]
            pos+=1
        if freqType!='Ø' and not fFix:
            freq=params[pos]*icmCoef
            pos+=1
        R.append((lifetime, lFix, freqType, freq, fFix))
        pass
    return R

#~ def combine(params, inputList):
    #~ for i in range(len(inputList)):
        #~ lifetime, lFix, freqType, freq, fFix=inputList[i]
        #~ if not lFix:
            #~ lifetime=1./params.pop(0)
        #~ if freqType!='Ø' and not fFix:
            #~ freq=params.pop(0)*icmCoef
        #~ inputList[i]=(lifetime, lFix, freqType, freq, fFix)
        #~ pass
    #~ pass
        
def model(data, t2, inputList):
    params, rates, w, rFix, wFix, types, paramsTypes = divide(inputList)
    c=C(rates, w, types, t2)
    cp=np.linalg.pinv(c)
    return np.tensordot(np.dot(c, cp), data, 1)

def DASandRes(data, t2, inputList):
    params, rates, w, rFix, wFix, types, paramsTypes = divide(inputList)
    c=C(rates, w, types, t2)
    cp=np.linalg.pinv(c)
    return np.tensordot(cp, data, 1), np.tensordot(np.eye(c.shape[0])-np.dot(c, cp), data, 1)
    
def DAS(data, t2, inputList):
    params, rates, w, rFix, wFix, types, paramsTypes = divide(inputList)
    c=C(rates, w, types, t2)
    cp=np.linalg.pinv(c)
    return np.tensordot(cp, data, 1)

def DASIRF(data, t2, inputList, IRFOffset, IRFFWHM):
    params, rates, w, rFix, wFix, types, paramsTypes = divide(inputList)
    c=CIRF(rates, w, types, t2, IRFOffset, IRFFWHM)
    cp=np.linalg.pinv(c)
    return np.tensordot(cp, data, 1)

def DASandResIRF(data, t2, inputList, IRFOffset, IRFFWHM):
    params, rates, w, rFix, wFix, types, paramsTypes = divide(inputList)
    c=CIRF(rates, w, types, t2, IRFOffset, IRFFWHM)
    cp=np.linalg.pinv(c)
    return np.tensordot(cp, data, 1), np.tensordot(np.eye(c.shape[0])-np.dot(c, cp), data, 1)

def setupGA(data, t2, inputList, maxfev=100):
    params, rates, w, rFix, wFix, types, paramsTypes = divide(inputList)

    if len(params)==0:
        import warnings
        warnings.warn('Nothing to fit')
        return None
    
    #TODO: storing is transferred from popot version, check if it is needed here (try to count how many times store is accessed for reading)
    #popot has a quirk that it reevaluates the best position of each particle, which is not needed here and cannot be switched off from python
    #storing the values will speed things up twice, hopefully it will not cause other problems
    store={}
    
    dataT = data.transpose(2,1,0)
    c = C(rates, w, types, t2)
    E = np.eye(c.shape[0])
    
    def evaluate(p):
        key=tuple(p)
        if key in store:
            return store[key]
        
        update(p, rates, w, types, rFix, wFix)
        c=C(rates, w, types, t2)
        cp=np.linalg.pinv(c)
        #for optimization attempts see sandbox/test_GA_performance.py

        #R=np.tensordot(np.eye(c.shape[0])-np.dot(c, cp), data, 1) #(I-C(C+))Y == Y - C(C+)Y == Y - M
        R = np.dot(dataT, (E-np.dot(c, cp)).T) #should be about 4% faster than tensordot (not much optimization)
        #other optimization would require limiting number of datapoints (do not consider areas with no signal)
        res = np.sum(np.abs(R)).real
        store[key]=res
        return res

    #bounds for lifetimes and frequencies
    dt2 = np.diff(t2).min() #but not zero
    if dt2==0:
        ttt = np.unique(t2)
        dt2 = np.diff(ttt).min()
    lb = {}
    ub = {}
    #r bounds
    #lower bound is theoretically 0, because we can have infinite lifetime (constant signal)
    lb["r"] = 0.
    #upper bound should be such that rate influences at least two frames, i.e. exp(-dt2*r)>=alpha -> r = -ln(alpha)/dt2
    alpha = 1e-4
    ub["r"] = -np.log(alpha)/dt2
    #w bounds
    #w lower is 0 rad/fs, but then it will be basically simple DAS
    lb["w"] = 0
    #w upper is given by Nyquist (we can do twice)
    ub["w"] = 2*np.pi/dt2 #if you change this, please adjust lMaxFreq calculation in tabGA
    
    #~ lbound=[1./(10*t2[-1])]*len(params)
    #~ ubound=[10./(t2[1]-t2[0])]*len(params)
    lbound = [lb[it] for it in paramsTypes]
    ubound = [ub[it] for it in paramsTypes]

    print("setupGA: starting params", params)
    print("bounds", lbound, ubound)

    swarm = PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8)
    print("setting starting particle")
    swarm.set_particle(0, params)
    swarm.store = store #might be usefull for analysis, i.e. to take a look where the swarm looked
    # ~ print(t2)
    return swarm #so that user can run it several times, or start/stop at will


def setupGAIRF(data, t2, inputList, IRF, maxAmplitude=None, maxfev=100, limit_to_start=None):
    params, rates, w, rFix, wFix, types, paramsTypes, IRF = divideIRF(inputList, IRF)

    params = np.asarray(params)

    if len(params)==0:
        import warnings
        warnings.warn('Nothing to fit')
        return None
    
    #TODO: storing is transferred from popot version, check if it is needed here (try to count how many times store is accessed for reading)
    #popot has a quirk that it reevaluates the best position of each particle, which is not needed here and cannot be switched off from python
    #storing the values will speed things up twice, hopefully it will not cause other problems
    store={}
    
    dataT = data.transpose(2,1,0)
    c = CIRF(rates, w, types, t2, IRF[2], IRF[0])
    E = np.eye(c.shape[0])
    
    def evaluate(p):
        key=tuple(p)
        if key in store:
            return store[key]
        
        #print("evaulate pre", IRF)
        updateIRF(p, rates, w, types, rFix, wFix, IRF)
        #print("evaluate after", IRF)
        c=CIRF(rates, w, types, t2, IRF[2], IRF[0])
        try:
            cp=np.linalg.pinv(c)
        except:
            print("WARNING: setupGAIRF evaluate pinv failed")
            print("c", c)
            print("t2", t2)
            print("rates", rates)
            print("w", w)
            print("types", types)
            print("IRF", IRF)
            raise
        #for optimization attempts see sandbox/test_GA_performance.py

        #R=np.tensordot(np.eye(c.shape[0])-np.dot(c, cp), data, 1) #(I-C(C+))Y == Y - C(C+)Y == Y - M
        R = np.dot(dataT, (E-np.dot(c, cp)).T) #should be about 4% faster than tensordot (not much optimization)
        #other optimization would require limiting number of datapoints (do not consider areas with no signal)
        res = np.sum(np.abs(R)).real
        store[key]=res
        return res

    #bounds for lifetimes and frequencies
    dt2 = np.diff(t2).min() #but not zero
    if dt2==0:
        ttt = np.unique(t2)
        dt2 = np.diff(ttt).min()
    lb = {}
    ub = {}
    #r bounds
    #lower bound is theoretically 0, because we can have infinite lifetime (constant signal)
    # however that might cause round off errors for UI, where inf cannot be displayed
    # so limit to 10*max(t2)
    lb["r"] = 0.1/t2.max()
    #upper bound should be such that rate influences at least two frames, i.e. exp(-dt2*r)>=alpha -> r = -ln(alpha)/dt2
    # with IRF two frames are not a good heuristics, because kinetic will be spread over whole pulse
    # but maybe we can do two frames in plain exp?
    if maxAmplitude is not None:
        alpha = 2 #lets allow components with their max amplitude at t2=0fs of alpha*maxAmplitude
        #this does not work properly if (selected) t2[0] is <=0fs, so limit it to dt2
        ub["r"] = -np.log(np.abs(data[0]).max()/(alpha*maxAmplitude))/max(t2[0], dt2)
    else:
        alpha = 1e-4
        ub["r"] = -np.log(alpha)/dt2
    ub["r"] = 100 #for testing
    #w bounds
    #w lower is 0 rad/fs, but then it will be basically simple DAS
    lb["w"] = 0
    #w upper is given by Nyquist (we can do twice - do not do twice, it will look like 0rad/fs)
    ub["w"] = 1.3*np.pi/dt2 #if you change this, please adjust lMaxFreq calculation in tabGA
    

    lb["O"] = -50. #lower offset
    ub["O"] = 50.  #upper offset
    lb["W"] = 0.01 #lower width
    ub["W"] = 100. #upper width
    
    
    #~ lbound=[1./(10*t2[-1])]*len(params)
    #~ ubound=[10./(t2[1]-t2[0])]*len(params)
    lbound = [lb[it] for it in paramsTypes]
    ubound = [ub[it] for it in paramsTypes]

    if limit_to_start is not None:
        lbound = params * (1-limit_to_start)
        ubound = params * (1+limit_to_start)

    print("setupGAIRF: starting params", params)
    print("bounds", lbound, ubound)

    swarm = PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8)
    print("setting starting particle")
    swarm.set_particle(0, params)
    #for testing
    #swarm.set_particle(1, [0.656743829471214, 0.0005576455505932945])
    #swarm.set_particle(2, [0.6578947368421053,0.0005576467307960407])
    swarm.store = store #might be usefull for analysis, i.e. to take a look where the swarm looked
    # ~ print(t2)
    return swarm #so that user can run it several times, or start/stop at will


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
        import cProfile
        pr = cProfile.Profile()
        pr.enable()

        swarm = setupGA(self._data, self._t2, self._inputList)
        if swarm is None: #if there is nothing to fit
            self.finalReport.emit([], 0, 0)
            return
            
        counter = 0
        while not self.isInterruptionRequested() and counter<self._maxSteps:
            #next step
            swarm.step()
            
            #emit results (or emit results every N steps)
            counter+=1
            if (counter % self.reportEachNth) == 0:
                self.report.emit(list(swarm.g), counter, swarm.fg)
            pass
        self.finalReport.emit(list(swarm.g), counter, swarm.fg)
        pr.disable()
        print()
        print()
        print()
        print("SWARM PROFILE REPORT")
        pr.print_stats()        
        print()
        print()
    pass

class SwarmThreadIRF(QtCore.QThread):
    report=QtCore.pyqtSignal(list, int, float)
    finalReport=QtCore.pyqtSignal(list, int, float)
    def __init__(self, parent, data, t2, inputList, IRF, reportEachNth=10, maxSteps=1000, maxAmplitude=None):
        super(SwarmThreadIRF, self).__init__(parent)
        self._data=data
        self._t2=t2
        self._inputList=inputList
        self.reportEachNth=reportEachNth
        self._maxSteps=maxSteps
        self._IRF = IRF
        self._maxAmplitude = maxAmplitude
        pass

    
    def run(self):
        try:
            import cProfile
            pr = cProfile.Profile()
            pr.enable()

            swarm = setupGAIRF(self._data, self._t2, self._inputList, self._IRF, self._maxAmplitude)
            if swarm is None: #if there is nothing to fit
                self.finalReport.emit([], 0, 0)
                return
                
            counter = 0
            while not self.isInterruptionRequested() and counter<self._maxSteps:
                #next step
                swarm.step()
                
                #emit results (or emit results every N steps)
                counter+=1
                if (counter % self.reportEachNth) == 0:
                    self.report.emit(list(swarm.g), counter, swarm.fg)
                pass
            self.finalReport.emit(list(swarm.g), counter, swarm.fg)
            pr.disable()
            print()
            print()
            print("SWARM FINAL CHECK")
            fg = swarm.obj(swarm.g)
            print(swarm.g, swarm.fg, fg)
            if fg != swarm.fg:
                QtWidgets.QMessageBox.warning(None, "WARNING", "PSO warning, reported fitness does not match final params\n (see stdout for details)")
            print()
            print()
            print()
            print("SWARM PROFILE REPORT")
            pr.print_stats()        
            print()
            print()
        except:
            print("WARNING: SwarmThreadIRF.run failed")
            #QtWidgets.QMessageBox.warning(None, "WARNING", "PSO warning, SwarmThreadIRF.run failed")
            import traceback
            traceback.print_exc()
            self.finalReport.emit([], 0, 0)
    pass

class SwarmThreadIRFBootstrap(QtCore.QThread):
    # ~ report=QtCore.pyqtSignal(list, int, float)
    # ~ finalReport=QtCore.pyqtSignal(list, int, float)
    finalReport=QtCore.pyqtSignal()

    def __init__(self, parent, filepath, data, t2, w1, w3, LO, inputList, IRF,
                 bootstrap_runs=100, bootstrap_reduction_factor=0.7, metadata=None,
                 reportEachNth=10, maxSteps=1000, maxAmplitude=None):
        super(SwarmThreadIRFBootstrap, self).__init__(parent)
        self._data=data
        self._filepath = filepath
        self._t2=t2
        self._inputList=inputList
        self.reportEachNth=reportEachNth
        self._maxSteps=maxSteps
        self._IRF = IRF
        self._maxAmplitude = maxAmplitude
        self._bootstrap = bootstrap_runs
        self._bootstrap_reduction_factor = bootstrap_reduction_factor
        self._metadata = metadata
        self._metadata["w1"] = w1
        self._metadata["w3"] = w3
        self._metadata["LO"] = LO
        pass

    
    def run(self):
        try:
            import cProfile
            pr = cProfile.Profile()
            pr.enable()

            import shelve
            with shelve.open(self._filepath) as storage:

                # it is possible that self._filepath already exists and contains previous run (bootstrap can take long
                #  time to compute, and it makes sense to split it to multiple sessions)
                # we should check that input is the same so that it makes sense to add new bootstrap runs to it

                if "data" in storage:
                    stored_data = storage["data"]
                    assert (stored_data[0] == self._data).all()
                    assert (stored_data[1] == self._t2).all()
                else:
                    storage["data"] = (self._data, self._t2)

                # FIXME: this will fail because numpy.ndarray cannot test with assert without .all()
                #  but we do not know which items in metadata need this :(
                #  numpy.testing.assert_allclose will work on single values or arrays, but not on strings
                if "start" in storage:
                    assert storage["start"] == (self._inputList, self._IRF)
                else:
                    storage["start"] = (self._inputList, self._IRF)

                if "metadata" in storage:
                    # allow extending metadata
                    stored_metadata = storage["metadata"]
                    for item in self._metadata:
                        if item in stored_metadata:
                            np.testing.assert_allclose()
                            assert stored_metadata[item] == self._metadata[item]
                        else:
                            stored_metadata[item] = self._metadata[item]
                    storage["metadata"] = stored_metadata
                else:
                    storage["metadata"] = self._metadata
                
                N = len(self._t2)
                Nshort = int(self._bootstrap_reduction_factor*N)
                
                for i in range(self._bootstrap):
                    #select a random subset of data
                    I = np.random.permutation(N)[:Nshort]
                    I.sort()
                    print("bootstrap selection", I)
                    print("\t\t", self._data.shape, self._t2.shape)

                    swarm = setupGAIRF(self._data[I], self._t2[I], self._inputList, self._IRF, self._maxAmplitude,
                                       limit_to_start=0.1)  # TODO: test limit_to_start
                    if swarm is None: #if there is nothing to fit
                        continue
                        
                    counter = 0
                    while not self.isInterruptionRequested() and counter<self._maxSteps:
                        #next step
                        swarm.step()
                        
                        #emit results (or emit results every N steps)
                        counter+=1
                        pass

                    # ~ self.finalReport.emit(list(swarm.g), counter, swarm.fg)
                    key = ",".join(str(it) for it in I)
                    while key in storage:
                        key += "_repeat"
                    storage[key] = (list(swarm.g), counter, swarm.fg)
            pr.disable()
                    
            print()
            print()
            print("BOOTSTRAP SWARM PROFILE REPORT")
            pr.print_stats()         # todo: save to file
            print()
            print()
        except:
            print("WARNING: SwarmThreadIRFBootstrap.run failed")
            #QtWidgets.QMessageBox.warning(None, "WARNING", "PSO warning, SwarmThreadIRF.run failed")
            import traceback
            traceback.print_exc()
        self.finalReport.emit()
    pass


def setupLA(data, t2, inputList, maxfev=100):
    params, rates, w, rFix, wFix, types, paramsTypes = divide(inputList)

    if len(params)==0:
        import warnings
        warnings.warn('Nothing to fit')
        return None
    
    #TODO: storing is transferred from popot version, check if it is needed here (try to count how many times store is accessed for reading)
    #popot has a quirk that it reevaluates the best position of each particle, which is not needed here and cannot be switched off from python
    #storing the values will speed things up twice, hopefully it will not cause other problems
    store={}
    storeCount = [0]
    
    #in LA version, data is single dimension kinetic
    dataT = data
    c = C(rates, w, types, t2)
    E = np.eye(c.shape[0])
    
    def evaluate(p):
        key=tuple(p)
        if key in store:
            storeCount[0] += 1
            return store[key]
        
        update(p, rates, w, types, rFix, wFix)
        c=C(rates, w, types, t2)
        cp=np.linalg.pinv(c)
        #for optimization attempts see sandbox/test_GA_performance.py

        #R=np.tensordot(np.eye(c.shape[0])-np.dot(c, cp), data, 1) #(I-C(C+))Y == Y - C(C+)Y == Y - M
        R = np.dot(dataT, (E-np.dot(c, cp)).T) #should be about 4% faster than tensordot (not much optimization)
        # TODO: we can add here some additional cost for amplitudes, to discourage overinflated and pair-compensated amplitudes
        #other optimization would require limiting number of datapoints (do not consider areas with no signal)
        res = np.sum(np.abs(R)).real
        store[key]=res
        return res

    def amplitudes(p):
        update(p, rates, w, types, rFix, wFix)
        c=C(rates, w, types, t2)
        cp=np.linalg.pinv(c)

        R = np.dot(dataT, cp.T)
        return R


    #bounds for lifetimes and frequencies
    dt2 = np.diff(t2).min() #but not zero
    if dt2==0:
        ttt = np.unique(t2)
        dt2 = np.diff(ttt).min()
    lb = {}
    ub = {}
    #r bounds
    #lower bound is theoretically 0, because we can have infinite lifetime (constant signal)
    lb["r"] = 0.  # do not do that, it makes more sense to say the lifetime is >>max(t2)
    lb["r"] = 1/3/t2.max()
    #upper bound should be such that rate influences at least two frames, i.e. exp(-dt2*r)>=alpha -> r = -ln(alpha)/dt2
    alpha = 1e-4
    ub["r"] = -np.log(alpha)/dt2
    #w bounds
    #w lower is 0 rad/fs, but then it will be basically simple DAS
    lb["w"] = 0
    #w upper is given by Nyquist (we can do twice)
    ub["w"] = 2*np.pi/dt2
    
    #~ lbound=[1./(10*t2[-1])]*len(params)
    #~ ubound=[10./(t2[1]-t2[0])]*len(params)
    lbound = [lb[it] for it in paramsTypes]
    ubound = [ub[it] for it in paramsTypes]

    print("setupLA: starting params", params)
    print("bounds", lbound, ubound)

    swarm = PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8)
    print("setting starting particle")
    swarm.set_particle(0, params)
    swarm.store = store #might be usefull for analysis, i.e. to take a look where the swarm looked
    swarm.storeCount = storeCount[0]
    swarm.amplitudes = amplitudes
    # ~ print(t2)
    return swarm #so that user can run it several times, or start/stop at will


def setupCA(data, t2, inputList, maxfev=100):
    params, rates, w, rFix, wFix, types, paramsTypes = divide(inputList)

    if len(params) == 0:
        import warnings
        warnings.warn('Nothing to fit')
        return None

    # TODO: storing is transferred from popot version, check if it is needed here (try to count how many times store is accessed for reading)
    # popot has a quirk that it reevaluates the best position of each particle, which is not needed here and cannot be switched off from python
    # storing the values will speed things up twice, hopefully it will not cause other problems
    store = {}
    storeCount = [0]

    # in CA version, data is single column (i.e. all data excited at single frequency)
    dataT = data.T
    c = C(rates, w, types, t2)
    E = np.eye(c.shape[0])

    def evaluate(p):
        key = tuple(p)
        if key in store:
            storeCount[0] += 1
            return store[key]

        update(p, rates, w, types, rFix, wFix)
        c = C(rates, w, types, t2)
        cp = np.linalg.pinv(c)
        # for optimization attempts see sandbox/test_GA_performance.py

        # R=np.tensordot(np.eye(c.shape[0])-np.dot(c, cp), data, 1)
        # (I-C(C+))Y == Y - C(C+)Y == Y - M
        R = np.dot(dataT, (E - np.dot(c, cp)).T)  # should be about 4% faster than tensordot (not much optimization)

        # TODO: we can add here some additional cost for amplitudes, to discourage overinflated and pair-compensated amplitudes
        # other optimization would require limiting number of datapoints (do not consider areas with no signal)
        A = np.dot(dataT, cp.T)
        res = np.mean(np.abs(R)).real + np.mean(np.abs(A))*1e-5 # add amplitude cost to prevent too large amplitudes

        # M = np.dot(c, A.T).T
        # R2 = dataT - M
        print("CA fitting amplitude ccost test", np.mean(np.abs(R)).real, np.mean(np.abs(A)), A)
        # print(R[:5, :5])
        # print(R2[:5, :5])
        # import pylab as pl
        # pl.figure()
        # pl.subplot(121)
        # pl.pcolormesh(dataT.real)
        # pl.colorbar()
        # pl.subplot(122)
        # print(A)
        # pl.pcolormesh(M.real.T)
        # pl.colorbar()
        # pl.show()


        store[key] = res
        return res

    def amplitudes(p):
        update(p, rates, w, types, rFix, wFix)
        c = C(rates, w, types, t2)
        cp = np.linalg.pinv(c)

        R = np.dot(dataT, cp.T)
        print("CA fitting amplitude test", R)
        return R

    # bounds for lifetimes and frequencies
    dt2 = np.diff(t2).min()  # but not zero
    if dt2 == 0:
        ttt = np.unique(t2)
        dt2 = np.diff(ttt).min()
    lb = {}
    ub = {}
    # r bounds
    # lower bound is theoretically 0, because we can have infinite lifetime (constant signal)
    lb["r"] = 0.  # do not do that, it makes more sense to say the lifetime is >>max(t2)
    lb["r"] = 1 / 3 / t2.max()
    # upper bound should be such that rate influences at least two frames, i.e. exp(-dt2*r)>=alpha -> r = -ln(alpha)/dt2
    alpha = 1e-4
    ub["r"] = -np.log(alpha) / dt2
    # w bounds
    # w lower is 0 rad/fs, but then it will be basically simple DAS
    lb["w"] = 0
    # w upper is given by Nyquist (we can do twice)
    ub["w"] = 2 * np.pi / dt2

    # ~ lbound=[1./(10*t2[-1])]*len(params)
    # ~ ubound=[10./(t2[1]-t2[0])]*len(params)
    lbound = [lb[it] for it in paramsTypes]
    ubound = [ub[it] for it in paramsTypes]

    print("setupCA: starting params", params)
    print("bounds", lbound, ubound)

    swarm = PSO(evaluate, lbound, ubound, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev, minstep=1e-8, minfunc=1e-8)
    print("setting starting particle")
    swarm.set_particle(0, params)
    swarm.store = store  # might be usefull for analysis, i.e. to take a look where the swarm looked
    swarm.storeCount = storeCount[0]
    swarm.amplitudes = amplitudes
    # ~ print(t2)
    return swarm  # so that user can run it several times, or start/stop at will


class LASwarmThread(QtCore.QThread):
    report = QtCore.pyqtSignal(float) #just percent of job done
    finalReport=QtCore.pyqtSignal(object, object)
    def __init__(self, parent, data, t2, inputList, maxSteps=1000, singleInputList=None, skipFilled=None):
        #if you give singleInputList (point-by-point staring parameters) it will be prefered, but swarm will use inputList as default value in case of problems
        super().__init__(parent)
        self._data=data
        self._t2=t2
        self._inputList=inputList
        self._maxSteps=maxSteps
        self._singleInputList = singleInputList
        self._skipFilled = skipFilled
        pass

    def run(self):
        import cProfile
        pr = cProfile.Profile()
        pr.enable()

        sh = self._data.shape

        if len(sh)!=3 or (np.asarray(sh)==0).any():
            self.finalReport.emit(None)
            return
            
        LA = [[None for j in range(sh[2])] for i in range(sh[1])]
        A = [[None for j in range(sh[2])] for i in range(sh[1])]
        
        storeCount = 0
        
        for i in range(sh[1]):
            for j in range(sh[2]):
                if self._skipFilled is not None and self._skipFilled[i][j] is not None:
                    LA[i][j] = self._skipFilled[i][j]
                    continue
                
                ttt = self._inputList if self._singleInputList is None else self._singleInputList[i][j]
                if ttt is None: ttt = self._inputList
                print("-"*60)
                swarm = setupLA(self._data[:,i,j], self._t2, ttt, self._maxSteps)
                if swarm is None: continue
                    
                # ~ counter = 0
                # ~ while not self.isInterruptionRequested() and counter<self._maxSteps:
                    # ~ #next step
                    # ~ swarm.step()
                    
                    # ~ #emit results (or emit results every N steps)
                    # ~ counter+=1
                    # ~ pass
                swarm.run() #this can auto stop, if some fit goodness criteria are met
                #Todo: we need to make sure that those criteria are well selected
                counter = swarm.it

                LA[i][j] = combineCopy(swarm.g, ttt)
                A[i][j] = swarm.amplitudes(swarm.g)
                print("\t\t", sh, i, j, LA[i][j], swarm.g, swarm.fg, counter)
                print("\t\tamplitudes", A[i][j])
                self.report.emit((i*sh[2]+j)/(sh[1]*sh[2])*100)
                storeCount += swarm.storeCount
                if self.isInterruptionRequested(): break
                pass
            if self.isInterruptionRequested(): break
            pass
        
        print("PSO store used", storeCount,"times")
        self.finalReport.emit(LA, A)
        pr.disable()
        print()
        print()
        print("-"*60)
        print("LA SWARM PROFILE REPORT")
        pr.print_stats()        
        print("-"*60)
        print()
        print()
    pass


class CASwarmThread(QtCore.QThread):
    report = QtCore.pyqtSignal(float)  # just percent of job done
    finalReport = QtCore.pyqtSignal(object, object)

    def __init__(self, parent, data, t2, inputList, maxSteps=1000, singleInputList=None, skipFilled=None):
        # if you give singleInputList (point-by-point staring parameters) it will be prefered, but swarm will use inputList as default value in case of problems
        super().__init__(parent)
        self._data = data
        self._t2 = t2
        self._inputList = inputList
        self._maxSteps = maxSteps
        self._singleInputList = singleInputList
        self._skipFilled = skipFilled
        pass

    def run(self):
        import cProfile
        pr = cProfile.Profile()
        pr.enable()

        sh = self._data.shape # should be [t2, w1, w3]

        if len(sh) != 3 or (np.asarray(sh) == 0).any():
            self.finalReport.emit(None)
            return

        CA = [None for i in range(sh[1])]
        A = [[None for i in range(sh[2])] for i in range(sh[1])]

        storeCount = 0

        for i in range(sh[1]):
            if self._skipFilled is not None and self._skipFilled[i] is not None:
                CA[i] = self._skipFilled[i]
                continue

            ttt = self._inputList if self._singleInputList is None else self._singleInputList[i][j]
            if ttt is None: ttt = self._inputList
            print("-" * 60)
            swarm = setupCA(self._data[:, i], self._t2, ttt, self._maxSteps)
            if swarm is None: continue

            # ~ counter = 0
            # ~ while not self.isInterruptionRequested() and counter<self._maxSteps:
            # ~ #next step
            # ~ swarm.step()

            # ~ #emit results (or emit results every N steps)
            # ~ counter+=1
            # ~ pass
            swarm.run()  # this can auto stop, if some fit goodness criteria are met
            # Todo: we need to make sure that those criteria are well selected
            counter = swarm.it

            CA[i] = combineCopy(swarm.g, ttt)
            # reshape this: A[i] = swarm.amplitudes(swarm.g)
            for j, a in enumerate(swarm.amplitudes(swarm.g)):
                A[i][j] = a # we need to make it list of arrays
            print("\t\t", sh, i, CA[i], swarm.g, swarm.fg, counter)
            print("\t\tamplitudes", A[i])
            self.report.emit((i / sh[1]) * 100)
            storeCount += swarm.storeCount
            if self.isInterruptionRequested(): break
            pass

        print("PSO store used", storeCount, "times")
        self.finalReport.emit(CA, A)
        pr.disable()
        print()
        print()
        print("-" * 60)
        print("LA SWARM PROFILE REPORT")
        pr.print_stats()
        print("-" * 60)
        print()
        print()

    pass


if __name__=="__main__":
    #test CIRF
    from pqr import pqr2 as pqr
    pqr.init()
    fw = pqr.FigureWidget()
    p = fw.plot()
    
    rates = [0.01, 0.01]
    frequencies = [3.14, -0.25] #1200./5308] #750icm to rad/fs conversion
    freqTypes = ['Ø', '+']
    
    t2 = np.arange(-50, 300, 0.1)

    eC = C(rates, frequencies, freqTypes, t2).T
    print(eC)
    eCIRF = CIRF(rates, frequencies, freqTypes, t2, 0, 25).T
    
    
    p.addSeriesPlotter(color="red").setData(t2, eC[0].real)
    p.addSeriesPlotter(p.x1, p.y1, color="darkred").setData(t2, eCIRF[0].real)
    # ~ p.addSeriesPlotter(p.x1, p.y1, color="darkred", style="..").setData(t2, eCIRF[0].real-eC[0].real)

    p.addSeriesPlotter(p.x1, p.y1, color="green").setData(t2, eC[1].real)
    p.addSeriesPlotter(p.x1, p.y1, color="darkgreen").setData(t2, eCIRF[1].real)
    
    p.zoomToData()
    
    """
    there is a problem 
    positive frequency in terms of density matrix and feynman diagrams means
    
    for model with expression OM_N(w1, w3) * exp(-1j*wN*t2-kN*t2)
    we have positive wN for Feyman digram with larger left side in waiting time part
    
    however numpy.fft.fft will place this at -wN
    we need to use numpy.fft.ifft to place it at wN
    
    and CIRF calculates wrong phase???
    """
    eCIRF = CIRF(rates, frequencies, freqTypes, t2, 0, 10).T
    # ~ eCIRF *= 10
    fw = pqr.FigureWidget()
    p = fw.plot()
    w2 = np.arange(len(t2))*(2*np.pi/(len(t2)*(t2[-1]-t2[0])/(len(t2)-1.)))
    p.addSeriesPlotter(color="red").setData(w2, np.abs(np.fft.ifft(eC[1])))
    p.addSeriesPlotter(p.x1, p.y1, color="blue").setData(w2, np.abs(np.fft.ifft(eCIRF[1])))
    p.addSeriesPlotter(p.x1, p.y1, color="green").setData(w2, np.abs(np.fft.ifft(np.exp(-1j*frequencies[1]*t2))))
    p.zoomToData()
    
    p = fw.addPlot(2,1)
    p.addSeriesPlotter(color="red").setData(t2, eC[1].real)
    p.addSeriesPlotter(p.x1, p.y1, color="red", style="..").setData(t2, eC[1].imag)
    p.addSeriesPlotter(p.x1, p.y1, color="blue").setData(t2, eCIRF[1].real)
    p.addSeriesPlotter(p.x1, p.y1, color="blue", style="..").setData(t2, eCIRF[1].imag)
    p.addSeriesPlotter(p.x1, p.y1, color="green").setData(t2, np.exp(-1j*frequencies[1]*t2).real)
    p.addSeriesPlotter(p.x1, p.y1, color="green", style="..").setData(t2, np.exp(-1j*frequencies[1]*t2).imag)
    p.zoomToData()
    
    # ~ OM_N(w1, w3) * exp(-1j*wN*t2-kN*t2)

    #compare different lifetimes shift with respect to 0fs
    fw = pqr.FigureWidget()
    p = fw.plot()        
    p.setup(x1={}, y1={})
    for it in [1, 10, 25, 100]:
        eCIRF = CIRF(rates, frequencies, freqTypes, t2, 0, it).T
        p.addSeriesPlotter(p.x1, p.y1, name=str(it)+"  "+str(t2[eCIRF[0].real.argmax()]), color=pqr.color()).setData(t2, eCIRF[0].real)
    p.zoomToData()
    pqr.show()

    """
    IMPORTANT TODO: this is seriously bad, if we include convolution with IRF, any fast oscillation is *strongly* suppressed
    no wonder out oscillations are very weak
    for 750/cm (typical Chl a vibration) the suppression is about 1/3 with FWHM 25fs pulse (which is roughly what we have)
    oscillations go to almost zero for 1200/cm vibration (thankfully we do not observe that that often)
    however for 0.4rad/fs spectral bandwith we can see potentially up to 0.4rad/fs oscillations (electronic coherences)
     - but that will not be visible
     - probably any oscillation that has period comparable or larger than FWHM will not be visible
     - for 25fs FWHM that is 0.25rad/fs, i.e. 1300/cm
     - we basically need FWHM less than optical period corresponding to laser bandwidth, but that is 15fs and that is not possible with out setup :(
    
    IMPORTANT TODO: if convolution with IRF is included, and we start at t2=0fs, there will be fast rising component, simply because of the convolution (but it might have been interpreted as system component by mistake - fortunately nothing is published yet)
     - it helps if we cut first 40fs (for FWHM 25fs) as after that the effect is negligible
    """
