
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

icmCoef=1./(2*np.pi*2.99792e-5) # w[1/cm] = w[rad/fs] * icmCoef

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

    window[np.where(np.logical_or(window<1e-100, np.isnan(window)))]=0.

    return window

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

#override cursor decorator
from PyQt5 import QtWidgets, QtGui, Qt
def busyCursor(function):
    def wrapped(self, *args, **kw):
        try:
            QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.Qt.WaitCursor))
            res = function(self, *args, **kw)
            QtWidgets.QApplication.restoreOverrideCursor()
        except:
            QtWidgets.QApplication.restoreOverrideCursor()
            raise
    wrapped.__name__='busyCursor_'+function.__name__
    return wrapped
   
if __name__=="__main__":
#    print((n(700, "nm", "fused silica")-n(700, "nm", "air"))*4e6 / c)
    for i in [500, 600, 700, 800]:
        print(2*n(i)*1e6/118568/299.792485*17787900)
