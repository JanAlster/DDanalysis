
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

def interpolate(x,n):
    """interpolate x to n points, keep derivation at edges"""
    if np.issubdtype(x.dtype, complex):
        print("complex")
        return np.fft.ifft(np.fft.fft(x), n)
    else:
        print("real")
        return np.fft.irfft(np.fft.rfft(x),n).astype(x.dtype)*(n/x.size)
        
def interpolate(x,n):
    """interpolate x to n points, keep derivation at edges"""
    return np.fft.irfft(np.fft.rfft(x), n).astype(x.dtype)*(n/x.size)

# def interpolate(y,M,x=None, axis=-1):
#     """advanced interpolator uses FT for interpolation and avoids ripples in case of jumps on edges"""
#     N=y.shape[axis]
#     if N==M or M<=0:
#         return y if x==None else (y,x)
#     
#     t1=arange(N) #or maybe arange(N+1)?
#     te=exp(-t1/(N/10))
#     
#     eshape=y.shape
#     eshape[axis]*=3
#     
#     Y=zeros(eshape, dtype=y.dtype)
#     
#     pos=tuple([slice(None, None) for i in eshape])
#     pos[axis]=slice(N, 2*N)
#     
#     Y[pos]=y
#
#
#     
#     (y[0]/((y[0]-y[1])*(t1+1))*te)[::-1]
#     ((y[-1]-y[-2])*(t1+1))/y[-1]*te
#     pass

def interpolate(y, M, x=None, axis=0):
    #TODO: efficiency
    #TODO: this cannot handle work/2019_solvents/sandbox/test_dataTTW2.py
    """advanced interpolator uses FT for interpolation and avoids ripples in case of jumps on edges
    1D y
    (note that this will uniformly add more points to y; optionally calculate new point positions if given positions of original points (x))
    """
    #swap axes
    if axis!=0:
        y = np.swapaxes(y, axis, 0)
    
    #interpolate first axis
    shape = y.shape
    N = shape[0]
    y = y.reshape((N, -1))
        
    if N==M or M<=0:
        return y if x==None else (y,x)

    MM = int((M*N)/(N-1)) #MM is the number we have to interpolate to to have M points in the cropped range (so that we do not overextend the original range)
    
    t1 = np.arange(N, dtype=np.float32)*N/(N-1) #or maybe arange(N+1)?
    te = np.exp(-t1/(N/30.))
    
    ttt = list(y.shape)
    ttt[0] = 3*N
    Y = np.zeros(ttt, dtype=y.dtype)
    
    print(Y.shape, Y[0:N].shape)
    print(y[0].shape)
    print(np.outer(t1+1, y[0]-y[1]).shape)
    
    #TODO: see frameListModel.updateUnphased for maybe a better version of extension
    
    Y[0:N]=((y[0]+np.outer(t1+1, y[0]-y[1])).T*te).T[::-1]
    Y[N:N*2]=y
    Y[2*N:]=((np.outer(t1+1, y[-1]-y[-2])+y[-1]).T*te).T
    
    FY = np.fft.fft(Y, axis=0)
    #get FY to MM size
    ttt = list(FY.shape)
    ttt[0] = 3*MM
    FFY = np.zeros(ttt, dtype=FY.dtype)
    
    hMM = int((3*MM+1.1)/2)
    hN = int((3*N+1.1)/2)
    
    if N>MM:
        #cut middle
        FFY[:hMM]=FY[:hMM]
        FFY[-hMM+1:]=FY[-hMM+1:]
    else:
        #add zeroes
        FFY[:hN]=FY[:hN]
        FFY[-hN+1:]=FY[-hN+1:]
        pass
        
    Y=np.fft.ifft(FFY, axis=0)
    #cut appropriate piece
    koef=(float(MM)/N)
    res=Y[MM:MM+M].astype(y.dtype)*koef
    print(MM,MM+M)

    shape = list(shape)
    shape[0]=M
    res = res.reshape(shape)
    if axis!=0:
        res = np.swapaxes(res, axis, 0)

    
    if x is None:
        return res
    else:
        dx=(x[-1]-x[0])/(N-1)
        X=np.r_[x[0]-N*dx:x[0]-N*dx+(3*MM-1)*dx/koef:dx/koef]
        return res, X[MM:MM+M]


#what I also want is (non-linear - use numpy.interp for that) interpolarion to arbitrary positions
def interpolateA(x, y, xnew):
    """
    take points (x, y) and estimate what would be the point of the same curve on xnew coordinates
    
    
    """
    raise NotImplementedError("This is WIP, does not work yet.")
    
    
    if x is None:
        x = np.arange(y.shape[-1])
        
    #we will use FT, but we must prevent waves from sharp edges
    #basically the simples way is 
    #smallest step which is not zero
    if False:
        ttt = np.diff(x)
        dx = ttt.min()
        if dx == 0:
            #TODO: we need to select larger (in case of repeating values in x)
            raise NotImplementedError
        
        ttt2 = np.diff(xnew)
        dx2 = ttt2.min()
        if dx2 == 0:
            #TODO: we need to select larger (in case of repeating values in x)
            raise NotImplementedError
        
        dx = min(dx, dx2) #this should give as range of FT axis

    #TODO: extend for N-D y
    print(x)
    print(y)
    if y[0]!=y[-1]: #there is a jump on edge (FT does not like that)
        A = 4 #best will work with 1, but that might be too large
        E = len(y)//A
            
        xe = np.arange(E, dtype=np.float32) #TODO: not sure why this redefines t1, t1s is used for w1 calculation

        #I do not want to depend on other modules here
        def cosWindow(x, start, end, edge):
            x = np.asarray(x)
            ttt = np.empty_like(x, dtype=float)
            s, se, ee, e = x.searchsorted([start, start+edge, end-edge, end])
            ttt[:s] = -np.pi
            ttt[s:se] = np.pi*(x[s:se]-start)/edge-np.pi
            ttt[se:ee] = 0
            ttt[ee:e] = np.pi*(x[ee:e]-(end-edge))/edge
            ttt[e:] = np.pi
            return np.cos(ttt)/2+0.5
            
        # ~ s = sigmoidWindow(range(2*E), E, E/2, E/20)
        s = cosWindow(range(2*E), 0, 2*E-1, E)

        shape = [it for it in y.shape]
        shape[-1] += 2*E
        ye = np.empty(shape, dtype=y.dtype)

        shape1 = list(y.shape)
        shape1[-1] = 1
        
        shapeE = list(y.shape)
        shapeE[-1] = E
        
        ettt = y[1]-y[0]
        ye[...,:E] = ((y[...,0]-ettt*E).reshape(shape1) + np.outer(ettt, xe).reshape(shapeE)) * s[:E]
        ye[...,E:-E] = y
        ye[...,-E:] = (y[...,-1].reshape(shape1) + np.outer(y[...,-1]-y[...,-2], xe+1).reshape(shapeE))*s[E:]
        
        y = ye
        
        #extend x too (should not matter how)
        dx = np.diff(x).mean()
        x = np.concatenate((xe*dx+(x[0]-E*dx), x, x[-1]+dx+xe*dx))
        # ~ return x, y
    
    print("extended")
    print(x)
    print(y)

    xnew_orig = xnew
    xnew = np.array(x)
    xnew[:len(xnew_orig)] = xnew_orig
        
    #TODO: this might not work as intended if x is nonuniform, or has repeating values
    dw = (2*np.pi/(len(x)*(x[-1]-x[0])/(len(x)-1.)))
    dw2 = (2*np.pi/(len(xnew)*(xnew.max()-xnew.min())/(len(xnew)-1.)))
    # ~ w = np.arange(N)*dw*dw2
    w = np.sort(np.concatenate([np.arange(len(x))*dw, np.arange(len(xnew))*dw2])) #combining x and xnew FT conjugate axes; this should make this pseudo-FT work as good as possible

    T = np.dot(np.exp(-1j*np.outer(xnew, w)), np.exp(1j*np.outer(w, x)))/(len(x)+len(xnew))
    print("T")
    print(T)
    print("ones")
    print(np.dot(T, np.ones_like(y)))
    return np.dot(T, y).real[:len(xnew_orig)]

def interpolate2D(y,M,x=None):
    #TODO: efficiency
    #TODO: what about axes???
    #TODO: looks like this is not finished
    """advanced interpolator uses FT for interpolation and avoids ripples in case of jumps on edges
    2D y
    """
    N=y.shape[0]
    
    if N==M or M<=0:
        return y if x==None else (y,x)

    MM=int((M*N)/(N-1)) #MM is the number we have to interpolate to to have M points in the cropped range (so that we do not overextend the original range)
    
    t1=np.arange(N, dtype=np.float32)*N/(N-1) #or maybe arange(N+1)?
    te=np.exp(-t1/(N/30.))
    
    Y=np.zeros((3*N,y.shape[1]), dtype=y.dtype)
    Y[0:N]=(np.outer(te,y[0])+np.outer((t1+1)*te,y[0]-y[1]))[::-1]
    Y[N:N*2]=y
    Y[2*N:]=np.outer((t1+1)*te,y[-1]-y[-2])+np.outer(te,y[-1])

    FY=np.fft.fft(Y, axis=0)
    #get FY to MM size
    FFY=np.zeros((3*MM,y.shape[1]), dtype=FY.dtype)
    
    hMM=int((3*MM+1.1)/2)
    hN=int((3*N+1.1)/2)
    
    if N>MM:
        #cut middle
        FFY[:hMM]=FY[:hMM]
        FFY[-hMM+1:]=FY[-hMM+1:]
    else:
        #add zeroes
        FFY[:hN]=FY[:hN]
        FFY[-hN+1:]=FY[-hN+1:]
        pass
        
    Y=np.fft.ifft(FFY, axis=0)

    #cut appropriate piece
    koef=(float(MM)/N)
    res=Y[MM:MM+M].astype(y.dtype)*koef
    
    if x==None:
        return res
    else:
        dx=(x[-1]-x[0])/(N-1)
        X=np.r_[x[0]-N*dx:x[0]-N*dx+(3*MM-1)*dx/koef:dx/koef]
        return res, X[MM:MM+M]
    
if __name__=='__main__':
    x=np.arange(10, dtype=np.float32)
    x=np.array([0,0,0,1.,2.,3.,4.,5.,0,0,0])
    #x=np.exp(-((np.arange(-2, 2, 0.1)/0.4)**2))
    x=np.exp(-(((np.arange(-2, 2, .1)-1.5)/0.4)**2))
    
    import sys
    sys.path.append("/home/araigne/Code/DD")
    sys.path.append("/home/alster/code/DD")
    # ~ from PQCP import PQCP
    from pqr import pqr2 as pqr
    pqr.init()
    
    f = pqr.FigureWidget()
    p = f.plot()
    X = list(range(len(x)))
    p.addSeriesPlotter().setData(X, x)

    #xx=interpolate(x, 100)
    print(X)
    #pl.plot(np.r_[0:len(x)-1:len(xx)*1j], xx)
    xx=np.fft.ifft(x)
    xx=np.r_[xx[:xx.size//2],np.zeros(70-xx.size),xx[xx.size//2:]]
    xx=np.fft.fft(xx)
    print(np.r_[0:len(x)-1:len(xx)*1j])
    # ~ PQCP.plot(np.r_[0:len(x)-1:len(xx)*1j], xx, color="red")
    p.addSeriesPlotter(p.x1, p.y1,color="red").setData(np.r_[0:len(x)-1:len(xx)*1j], xx)
    
    Y,XX = interpolate(np.array([x, x]), 70, X, axis=1)
    print(XX)
    for line in Y:
        # ~ PQCP.plot(X,line, color="green")
        p.addSeriesPlotter(p.x1, p.y1, color="green").setData(XX, line)
    
    
    Xnew = X + np.random.randn(len(X))*(0.1*np.diff(X).mean())
    Xnew = X
    Xnew = X[::2]
    # ~ Xnew = np.asarray(X) +0.1
    xnew = interpolateA(X, x, Xnew)
    # ~ xxx, yyy = interpolateA(X, x, Xnew)
    p.addSeriesPlotter(p.x1, p.y1,color="purple").setData(Xnew, xnew)
    # ~ p.addSeriesPlotter(p.x1, p.y1,color="purple").setData(xxx, yyy)
    
    p.zoomToData()
    
    #test
    #ones
    #centered gaussian
    #centered sine
    #line (sharp jump on edge)
    #jagged line (sharp jumps in the middle)
    
    #N-D data
    
    
    pqr.show()
    pass

    
