
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

scheme_lund={-1: '#190093', -0.75: '#7C07A9', -0.5: '#da70d6', -0.25:'#1e90ff', 0.: '#FFFFFF', 0.25: '#00FB00', 0.50: '#FFFA00', 0.75: '#FF0000', 1: '#A52A2A'}
scheme_hotAndCold={-1:'#0000FF', 0.:'#FFFFFF', 1: '#FF0000'}

import numpy as np

class Colormap(object):
    def __init__(self, scheme=scheme_lund):
        self.scheme=scheme
        self.points=sorted(scheme.keys())
        colors=np.array([int(scheme[it][1:],16) for it in self.points])
        #print [np.base_repr(it, base=16) for it in colors]
        self.B=colors.astype(np.uint8)
        #print [np.base_repr(it, base=16) for it in self.B]
        colors=np.right_shift(colors, 8)
        self.G=colors.astype(np.uint8)
        #print [np.base_repr(it, base=16) for it in self.G]
        colors=np.right_shift(colors, 8)
        self.R=colors.astype(np.uint8)
        #print [np.base_repr(it, base=16) for it in self.R]
        
        self.dR=[(float(self.R[i+1])-self.R[i])/(self.points[i+1]-self.points[i]) for i in range(len(self.points)-1)]
        self.dG=[(float(self.G[i+1])-self.G[i])/(self.points[i+1]-self.points[i]) for i in range(len(self.points)-1)]
        self.dB=[(float(self.B[i+1])-self.B[i])/(self.points[i+1]-self.points[i]) for i in range(len(self.points)-1)]
        pass

    def __call__(self, data, vmin=None, vmax=None, symmetrical=False):
        #data to -1, 0, 1
        try:
            iter(data)
        except:
            is_scalar = True
            data=np.array([data])
        else:
            is_scalar=False
            data=np.array(data)


        if vmin==None: vmin = data.min()
        if vmax==None: vmax = data.max()
        if symmetrical:
            vmax=max(abs(vmax), abs(vmin))
            vmin=-vmax

        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result=zeros_like(data, dtype=np.float32)
            pass
        else:
            vmin = float(vmin) #just in case it is int
            vmax = float(vmax)

            if vmin>=0 or vmax<=0:
                if vmin>=0:
                    result=(data/vmax).astype(np.float32)
                else:
                    result=(data/abs(vmin)).astype(np.float32)
                pass
            else:
                result=(data/vmax).astype(np.float32)
                result[np.where(result<0)]/=abs(vmin)/vmax
                pass

        #normalized data to RGB
        R=np.zeros_like(result, dtype=np.uint8)
        G=np.zeros_like(result, dtype=np.uint8)
        B=np.zeros_like(result, dtype=np.uint8)
        for i in range(len(self.points)-1):
            indices=np.where(np.bitwise_and(result>=self.points[i], result<self.points[i+1]))
            if len(indices)>0:
                ttt1=result-self.points[i]
                R[indices]=(self.R[i]+self.dR[i]*ttt1).astype(np.uint8)[indices]
                G[indices]=(self.G[i]+self.dG[i]*ttt1).astype(np.uint8)[indices]
                B[indices]=(self.B[i]+self.dB[i]*ttt1).astype(np.uint8)[indices]
        
        indices=np.where(result==self.points[-1])
        R[indices]=self.R[-1]
        G[indices]=self.G[-1]
        B[indices]=self.B[-1]
        
        argb=np.ones(result.shape, dtype=np.uint32)*255
        argb=np.left_shift(argb,8)
        argb+=R
        argb=np.left_shift(argb,8)
        argb+=G
        argb=np.left_shift(argb,8)
        argb+=B
        
        if is_scalar: argb = argb[0]
        return argb
    pass




def dataToRGB(data, colormap='lund'):
    rr=np.zeros(data.shape, dtype=np.uint8)
    argb=np.ones(data.shape, dtype=np.uint32)*255
    rr=((1-np.abs(data)/np.abs(data).max())*255).astype(np.uint8)
    #print
    #print rr[0,0], np.base_repr(rr[0,0], base=16)
    #print argb[0,0], np.base_repr(argb[0,0], base=16)
    argb=np.left_shift(argb,8)
    #print argb[0,0], np.base_repr(argb[0,0], base=16)
    argb+=rr
    #print np.base_repr(argb[0,0], base=16)
    argb=np.left_shift(argb,8)
    #print np.base_repr(argb[0,0], base=16)
    argb+=rr #as gg
    #print np.base_repr(argb[0,0], base=16)
    argb=np.left_shift(argb,8)
    #print np.base_repr(argb[0,0], base=16)
    argb+=rr #as bb
    #print np.base_repr(argb[0,0], base=16)
    return argb
