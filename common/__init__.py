import logging

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
import numpy as np

"""
FFT axes
t3 = np.arange(len(w3))*(2*np.pi/(len(w3)*(w3[-1]-w3[0])/(len(w3)-1.)))
w1 = np.arange(len(t1))*(2*np.pi/(len(t1)*(t1[-1]-t1[0])/(len(t1)-1.)))
"""
def FFTaxis(x):
    N = len(x)
    return np.arange(N)*(2*np.pi/(N*(x[-1]-x[0])/(N-1.)))

def normalize(x):
    return x/x.max()

#might come in handy
def subdivide(a, N=1):
    res = np.empty(len(a)*2-1, dtype=a.dtype)
    res[::2] = a
    res[1::2] = 0.5*(a[1:]+a[:-1]) 
    return res if N==1 else subdivide(res, N-1)


def moving_average(a, n=3) :
    """https://stackoverflow.com/a/14314054"""
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

c = 299.792485 #nm/fs

unit_to_nm={'nm': lambda x: x, 'um': lambda x:1000.*x, '1/cm': lambda x:1e7/x, 'rad/fs': lambda x:np.pi*2*299.792485/x}
materials={'water':([0.5684027565, 0.1726177391, 0.02086189578, 0.1130748688],[5101.829712, 18211.53963, 26207.22293, 10697927.21]),
'air':([5.792105e-8, 1.67917e-9],[0.0002380185, 5.7362e-5]),
'fused silica':([0.6961663, 0.4079426, 0.8974794],[4679.14826, 13512.0631, 97934002.5]),
'bk7':([1.03961212, 0.231792344, 1.01046945],[6000.69867, 20017.9144, 103560653.]),
'sapphire ordinary':([1.4313493, 0.65054713, 5.3414021],[5279.9261, 14238.2647, 325017834.]),
'sapphire extraordinary':([1.5039759, 0.55069141, 6.5927379],[5480.41129, 14799.4281, 402895140.]),
'soda lime glass':[1.5130, -0.003169*1e-6, 0.003962*1e6]}

#Done: check n definitions by refractionindex www - checks out OK
def n(x, unit="nm", material="air"):
    ret=np.ones_like(x, dtype=float)
    if material=='vacuum':
        pass
    elif material=='air':
        x=unit_to_nm[unit](np.asarray(x))**-2.
        for b,c in zip(*materials['air']): ret+=b/(c-x)
    elif material=="soda lime glass":
        x=unit_to_nm[unit](np.asarray(x))**2
        a, b, c = materials[material]
        ret = a + b*x + c/x
    else:
        x=unit_to_nm[unit](np.asarray(x))**-2.
        for b,c in zip(*materials[material]): ret+=b/(1.-c*x)
        ret=ret**0.5
    return ret

def phase_shift(w, thickness, material='fused silica', material2='air'):
    return np.exp(1j*thickness/0.299792458*w*(n(w,'raf/fs',material)-n(w,'rad/fs',material2)))

icmCoef=1./(2*np.pi*2.99792e-5) # w[1/cm] = w[rad/fs] * icmCoef (icmCoef = 5308.8455693245762318168517963293 fs/rad/cm)

#frequency axis for FT
# ~ wfine = np.arange(len(tfine))*(2*np.pi/(len(tfine)*(tfine[-1]-tfine[0])/(len(tfine)-1.)))

def cosWindow(x, start, end, edge):
    ttt = np.empty_like(x, dtype=float)
    s, se, ee, e = x.searchsorted([start, start+edge, end-edge, end])
    ttt[:s] = -np.pi
    ttt[s:se] = np.pi*(x[s:se]-start)/edge-np.pi
    ttt[se:ee] = 0
    ttt[ee:e] = np.pi*(x[ee:e]-(end-edge))/edge
    ttt[e:] = np.pi
    return np.cos(ttt)/2+0.5

def sigmoidWindow(x, center, width, steepness):
    x=np.asarray(x)

    width*=0.5

    t1=(np.exp(-width/steepness)+1.)**-2 if steepness!=0 else 1.
    try:
        t2=1./(np.exp(-(x-center+width)/steepness)+1)
    except RuntimeWarning:
        print("RuntimeWarning")
        print(-(x-center+width)/steepness)
        raise
    
    try:
        t3=1./(np.exp(-(width+center-x)/steepness)+1)
    except RuntimeWarning:
        print("RuntimeWarning")
        print(-(width+center-x)/steepness)
        raise

    window=t2*t3/t1

    #~ window[np.where(np.logical_or(window<1e-100, np.isnan(window)))]=0.
    window[np.isnan(window)] = 0
    window[window<1e-100] = 0

    return window

def cosWindow(x, start, end, edge):
    #this assumes x ascending
    x = np.asarray(x)
    y = np.empty_like(x)
    
    #sanitize
    if end<start: start, end = end, start
    edge = min(edge, (end-start)/2)
    
    #parts of window
    s, se, ee, e = x.searchsorted((start, start+edge, end-edge, end))
    
    y[:s] = -np.pi
    y[s:se] = (x[s:se]-start)/edge*np.pi-np.pi
    y[se:ee] = 0
    y[ee:e] = (x[ee:e]-(end-edge))/edge*np.pi
    y[e:] = np.pi
    return 0.5*(np.cos(y)+1)

# phase unwrap
def midunwrap(a, axis=-1):
    N = a.shape[axis]//2
    t = a.take(N, axis=axis)
    a = np.unwrap(a, axis=axis)
    s = np.array(a.shape)
    s[axis] = 1
    t -= a.take(N, axis=axis)
#    print(a.shape, axis, t.shape, s)
    a += t.reshape(s) #but we likely need some magic here
    return a

def midunwrap2d(a):
    return midunwrap(midunwrap(a, axis=0))

#pretty print of angle
def rad2dms(angle):
    angle=angle/np.pi*180
    d=int(angle)
    angle=(angle-d)*60
    m=int(angle)
    s=(angle-m)*60
    return d, m, s
    
def rad2dms_s(angle):
    return "{}° {}' {:.2f}''".format(*rad2dms(angle))

#pretty print of numbers
# TODO: does this also use superscript "-"?
reg = "-0123456789"
sup = "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"
tr = "".maketrans(reg, sup)


def formatValue(x):
    """
    pretty print mean value with std range (single sigma)
    """
    mx = x.mean()
    sx = x.std()
    
    # ~ print(sx)
    if sx == 0:
        o = 0
    else:
        o = int(np.floor(np.log10(sx)))

    o = -o +1
    if o >= -3:    
        t = f"{np.round(mx, o):.{max(0, o)}f}±{np.round(sx, o):.{max(0, o)}f}" #single sigma
        # ~ t = f"${np.round(mx, o):.{max(0, o)}f}\\pm{np.round(sx, o):.{max(0, o)}f}$" #single sigma
        # ~ t = fr"\num{{{np.round(mx, o):.{max(0, o)}f}\pm{np.round(sx, o):.{max(0, o)}f}}}"
        # ~ t = fr"{np.round(mx, o):.{max(0, o)}f}+-{np.round(sx, o):.{max(0, o)}f}"
        # ~ t = f"{mx}({sx})"
    else:
        # we need exp format
        ttt = 10**(-o)
        t = fr"{np.round(mx, o)/ttt:.0f}±{np.round(sx, o)/ttt:.0f} x10{str(-o).translate(tr)}"
        print("format to order of ten", mx, sx, t)
        # ~ t = fr"(&{np.round(mx, o)/ttt:.0f}& \pm ${np.round(sx, o)/ttt:.0f}$) 10^{-o}"
        # ~ t = fr"\num{{{np.round(mx, o)/ttt:.0f} \pm {np.round(sx, o)/ttt:.0f}e{-o}}}"
        # ~ t = fr"{np.round(mx, o)/ttt:.0f} +- {np.round(sx, o)/ttt:.0f}e{-o}"
        

    
    # ~ t = f"{np.round(mx, -(o-1)):.{max(0, -o+1)}f}±{np.round(sx, -(o-1)):.{max(0, -o+1)}f}" #single sigma
    # ~ t = f"{mx}({sx})"
    return t, mx, sx


#catch exception decorator
def catchException(function):
    """
    Calling code has to be able to handle None as return value (best if function does not return value).
    """
    def wrapper(*args, **kw):
        try:
            return function(*args, **kw)
        except:
            print("-"*60)
            print("catchException decorator has caught exception for function", function)
            import traceback
            traceback.print_exc()
            print("-"*60)
    wrapper.__name__='catchException_'+function.__name__
    return wrapper


#TODO: we could do this usable with "with ...: code block, which would guarantee the cursor is restored in case of crash inside
#override cursor decorator
from PyQt5 import QtWidgets, QtGui, QtCore
def busyCursor(function):
    logging.debug("Checking if Qt is running")
    if QtCore.QCoreApplication.instance() is None:
        return function

    def wrapped(self, *args, **kw):
        try:
            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
            res = function(self, *args, **kw)
            QtWidgets.QApplication.restoreOverrideCursor()
            return res
        except:
            QtWidgets.QApplication.restoreOverrideCursor()
            raise
    wrapped.__name__='busyCursor_'+function.__name__
    return wrapped

#delayed update decorator
from collections import OrderedDict
class LastUpdatedOrderedDict(OrderedDict):
    'Store items in the order the keys were last added'

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        self.move_to_end(key)
        
_delayedUpdateStore = LastUpdatedOrderedDict() #will store function which were delayed, but only the last call
_delayedUpdateLock = 0

def delayUpdate():
    global _delayedUpdateLock
    _delayedUpdateLock += 1

#TODO: this is not safe in the way that it stores inner functions of some objects to be called later
#   if the object itself is destroyed (or should be destroyed) while the update is locked
#   it will lead to undefined behaviour
# do not decorate anything that is not persistend during update lock
#TODO: this will be a global lock, perhaps we might need to use several lock in parallel?
def releaseUpdate():
    global _delayedUpdateLock
    if _delayedUpdateLock==0:
        print("WARNING: releaseUpdate (common delayedUpdate decorator) - update is not delayed, but release was called")
    else:
        _delayedUpdateLock -= 1
        global _delayedUpdateStore
        while _delayedUpdateLock == 0 and len(_delayedUpdateStore)>0:
            #proceed with update
            #we should be carefull with lock, because any call here might lock update again
            function, (args, kw) = _delayedUpdateStore.popitem(False)
            function(*args, **kw)
    pass

#to be used as decorator for function calls that should be delayed under some conditions (controlled by delayUpdate, releaseUpdate)
#NOTE that decorated functions must return None (or rather that it is impossible to get the return value)
def delayedUpdate(function):
    def wrapped(*args, **kw):
        global _delayedUpdateLock
        global _delayedUpdateStore
        if _delayedUpdateLock>0:
            _delayedUpdateStore[function] = (args, kw)
            return None
        else:
            #we are not locked, but function could be scheduled to be called by us
            #e.g. we have several functions scheduled to be called during delayUpdate lock
            # and after release, first one triggers one of the others - it is no longer needed to call it second time
            # ~ if function in _delayedUpdateStore:
                # ~ del _delayedUpdateStore[function]
            # ~ return function(*args, **kw)
            
            #alternatively we can keep order of calling
            if function in _delayedUpdateStore:
                _delayedUpdateStore[function] = (args, kw)
            else:
                return function(*args, **kw)
    wrapped.__name__='delayedUpdate_'+function.__name__
    return wrapped
    

   
if __name__=="__main__":
#    print((n(700, "nm", "fused silica")-n(700, "nm", "air"))*4e6 / c)
    for i in [500, 600, 700, 800]:
        print(2*n(i)*1e6/118568/299.792485*17787900)
