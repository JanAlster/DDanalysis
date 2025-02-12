
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

#NOTE: there is a possible bug - once we have _u and _a and clear intermediates, we can still setup new _p or _P
#   which will not be able to get new _u and _a (because other things are missing)
#   but will be used in model calculation - i.e. wrong results
#   possible fix is to keep copy of the last _p, _P, _u used for _a calculation
#   but I will probably leave it as it is and whoever decides to abuse things will have to live with it

import numpy as np
from numbers import Number

class Spline(object):
    #~ #degree
    #~ p = 3
    #~ #control points - will be fitted (for 1d curve these are heights of basis functions)
    #~ n = int(len(y)/5) #number of control points, this is also the number of basis functions of degree p
    #~ #n = len(y)
    #~ #100 is not enough for sin(x**2), len(y) will describe it perfectly (it hits every knot by definition), probably looks awful on finer grid

    #~ #knot vector
    #~ m = n + p + 1 #number of knots, also number (+1) of basis functions of degree 0
    #~ #knot vector u should span x
    #~ #but only u[p:-p] will have p basis functions on each knot spacing (i.e. u[:p] and u[-p:] will not be fitted by arbitrary/best p degree polynom, but by something forced by inner points)

    #~ N = len(y)

    #~ P=4

    #~ #choose knot grid
    #~ indices = np.r_[0:len(x)-1:1j*(m-2*P)].astype(np.int)
    #~ #note that start and end needs k repetitions to fit curve to the edges
    #~ indices = np.concatenate(([indices[0]]*P, indices, [indices[-1]]*P)).astype(np.int)
    #~ #print(indices)
    #~ u = x[indices] #knot vector
    #~ uy = y[indices] #values at know points, useful only for plotting
    
    #we can setup knot vector once we know x range and number of internal knots (and degree)
    #we can setup basis functions once we know the whole x axis (likely to be huge and full of zeros)
    #we can calculate control point values once we know y data and smoothing parameter
    
    #Todo: this can be optimized, as most of B is zeros
    
    #Spline can handle multidimensional y with some limitations (basically only setup and model are tested, math might work although is not tested and inverse will probably not work)
    # - also not tested for more than 2D y
    
    def __init__(self, *args, **kw):
        super().__init__()
        #it init we need to know just the x range, number on (inner) knots (will be equispaced) and degree
        self._x = None #data
        self._y = None #data
        self._n = None #number of inner knots
        self._p = None #degree of spline
        self._P = None #multiplicity of edge knots, reduces edge effects (should be at least p+1)
        self._l = None #smoothing parameter
        
        #structures
        self._B = None #basis functions
        self._BB = None
        self._DD = None
        self._I = None
        self._u = None
        self._a = None #control points of bspline that fits the data

        if len(args)>0 and isinstance(args[0], Spline):
            other = args[0]
            self._a = np.array(other._a) if other._a is not None else None #copy
            self._u = np.array(other._u) if other._u is not None else None #copy
            self._P = other._P
            self._p = other._p
            self._I = np.array(other._I) if other._I is not None else None
        else:
            self.setup(*args, **kw)
        pass
        
    def setup(self, x=None, y=None, n=None, p=None, P=None, l=None, extendMin=None, extendMax=None):
        if x is None and y is None and n is None and p is None and P is None and l is None: return #from default __init__
        
        #changing x, n, P, or p needs new B
        #changing x, y needs new sort
        #changing just y needs to sort y, it has to be the same order as old one (i.e. one-on-one correspondence with old x, not the one sorted here)
        #changing n, P, or p needs new D
        #changing l, B or D needs new inverse
        #changing y or inverse needs new model
        
        #extendMin, extendMax are cheat variables which can extend range of Spline beyond range of x (it will be poorly defined there so do it at your own risk)
        # - they are evaluated every time the Basis is recalculated
        # - they are not stored (preferably should not be used so the use is not convenient, this might change in the future)
        
        calculateBasis = False
        calculateD = False
        calculateI = False
        calculateA = False
        
        if x is None and y is None:
            x = self._x
            y = self._y
        elif x is None:
            #new y, has to be the same order as old y
            y = y[self._sind]
            self._y = y
            calculateA = True
        elif y is None:
            #just new x - unset y, cannot perform spline fit
            #but possibly we could (in principle) somehow duplicate the object and use different y values for same x
            sind = x.argsort()
            x = x[sind]
            self._x = x
            self._sind = sind
            self._y = None #assuming old y is not compatible with new x
            
            calculateBasis = True
        else:
            #we need sorted data
            #this could be slow and is not in place, but is easy to code
            sind = x.argsort()
            x = x[sind]
            y = y[sind]
            #we will have local copy of data
            self._x = x
            self._y = y
            self._sind = sind 
            
            calculateBasis = True
            
        if p is None:
            p = self._p
        else:
            self._p = p
            calculateBasis = True
            calculateD = True
            
        if n is None:
            n = self._n
        else:
            self._n = n
            calculateBasis = True
            calculateD = True
        
        if P is None:
            P = self._P
        else:
            self._P = P
            calculateBasis = True
            calculateD = True
        
        if l is None:
            l = self._l
        else:
            self._l = l
            calculateI = True

        #print("Spline.setup", x is None, y is None, p, P, l, n)

        #calculations
        if calculateBasis and not (p is None or P is None or n is None or x is None):
            calculateI = True
            #print("\tcalculateBasis")
            xmin = x.min()
            xmax = x.max()

            #we can avoid problems with inclusion of last value if we move knots just a little bit out from xmin, xmax
            dx = 2*(xmax-xmin)/n #setting dx to 2 actually helps a bit to reduce edge effects (together with P = 40 ~ 60)
            inner = np.r_[min(extendMin, xmin-dx) if extendMin is not None else xmin-dx:max(extendMax, xmax+dx) if extendMax is not None else xmax+dx:n*1j]
            u = np.concatenate(([inner[0]]*P, inner, [inner[-1]]*P)) #knot vector
            self._u = u
            
            #this might be usefull
            uind = x.searchsorted(u)
            
            m = len(u) #number of knots - also number (+1) of basis functions of degree 0
            #nn = m - p - 1 #number of control points (will be fitted from data) - also number of basis functions of degree self._p
            #print("u", len(u), n+2*P, "nn", m - p - 1, 2*P + n -p -1)
            N = len(x)
            
            #base function matrix
            B = np.zeros((N, m-1))

            for j in range(P, m-1-P): #the other will be 0
                B[uind[j]:uind[j+1],j] = 1

            def Bnext(Bprev, k):
                du = [u[j+k]-u[j] for j in range(m-k)]
                Bk = np.zeros((N, m-1-k))
                for j in range(0, m-1-k):
                    if du[j] != 0:
                        Bk[uind[j]:uind[j+k],j] += (x[uind[j]:uind[j+k]]-u[j])/du[j] * Bprev[uind[j]:uind[j+k],j]
                    if du[j+1] != 0:
                        Bk[uind[j+1]:uind[j+k+1],j] += (u[j+1+k]-x[uind[j+1]:uind[j+1+k]])/du[j+1] * Bprev[uind[j+1]:uind[j+k+1], j+1]
                    pass
                pass
                
                return Bk
            
            #NOTE that Bks and Ds are really sparse matrices, but I do not want to go there (might as well install whole scipy)
            for ik in range(1, p+1):
                #print("B degree", ik)
                B = Bnext(B, ik)
                
            self._B = B
            self._BB = np.dot(B.T, B)
            pass

        if calculateD and not (n is None or P is None or p is None):
            #print("\tcalculateD", n, P, p)
            calculateI = True

            nn = 2*P + n -p -1            
            D = np.diff(np.diff(np.diff(np.eye(nn+3)))) #3 for number of diff
            self._DD = np.dot(D.T, D) 
            
        if calculateI and not (self._BB is None or self._DD is None or self._B is None or l is None):
            #print("\tcalculateI")
            calculateA = True
            self._I = np.dot(np.linalg.inv(self._BB+l*self._DD), self._B.T)
            
        if calculateA and not (self._I is None or self._y is None):
            #print("\tcalculateA")
            self._a = np.dot(self._I, self._y) #this is what we finally need to evaluate the spline
        pass

    def ready(self):
        return not (self._u is None or self._p is None or self._a is None or self._P is None)

    def __iadd__(self, other):
        #only for number other - different type will lead to undefined behaviour
        #we need _I for this, so cannot be done after clearIntermediates
        #todo: this could be extended to vector other
        self._a+= np.sum(self._I, axis=1)*other
        return self

    def __imul__(self, other):
        #only for number other - different type will lead to undefined behaviour
        self._a *= other
        return self

    #Note: this can get called with all sorts of garbage from UI
    #Note that this does not respect x type, it will return numpy.float64
    def model(self, x=None):
        """
        with x==None returns (x, spline(x)) on original x used for spline setup
        otherwise returns just spline(x)
        """
        if x is None:
            #this is model on original x used for setup
            #TODO: it is quite likely that the original dataset self._x, self._y has multiple y for one x
            #   i.e. there are repeating columns in B
            #   this is needed to allow easy fit via linalg
            #   however it also means we have duplicates in values returned here, where is no need for them
            return self._x, np.dot(self._B, self._a)
        else:
            #check data type
            #print("Spline.model", x, Number)
            isNumber = isinstance(x, Number)
            #print("\t", isinstance(x, Number))
            if isinstance(x, list): x = np.array(x)
            isNumberArray = isinstance(x, np.ndarray) and x.dtype.kind in "iuf"
            #note that we do not accept list
            if not (isNumber or isNumberArray):
                raise TypeError("Spline.model does not accept data of type", x.__class__)

            #check data value
            if isNumber and not (x>=self._u[0] and x<=self._u[-1]):
                raise ValueError("Spline.model: value", x, "out of spline range", self._u[0], self._u[-1])
            if isNumberArray and not np.logical_and(x>=self._u[0], x<=self._u[-1]).all():
                raise ValueError("Spline.model: value range", x.min(), x.max(), "out of spline range", self._u[0], self._u[-1])

            #model on arbitrary x
            #TODO? remove invalid x
            u = self._u
            p = self._p
            i = u.searchsorted(x)
            if isNumberArray:
                i[np.where(i==0)] = self._P+1
                #i[np.where(i==len(u))] = len(u)-self._P-1  #TODO: but this likely means attempt at extrapolation
                k = i-1
                aa = np.zeros((p+1,p+1, len(x))+self._a.shape[1:])
                for i in range(len(k)):
                    aa[0, :, i] = self._a[k[i]-p:k[i]+1]
                pass
            else:
                if i==0: i = self._P+1
                #if i==len(u): i = len(u)-self._P-1 #TODO: but this likely means attempt at extrapolation
                k = i-1
                aa = np.zeros((p+1,p+1))
                aa[0] = self._a[k-p:k+1]
                pass
            
            for r in range(1, p+1):
                for i in range(-p+r, 1):
                    alfa = (x-u[k+i])/(u[k+i+p+1-r]-u[k+i])
                    aa[r, i+p] = ((1-alfa) * aa[r-1, i-1+p].T).T + (alfa * aa[r-1, i+p].T).T
            
            return aa[r, -1]
        pass

    #Note: this can get called with all sorts of garbage from UI
    def inverse(self, y, inversePrecision=1e-8):
        #print("Spline.inverse", y)
        """
        This works only for monotonous splines and we also assume that they are mostly linear (our use case).
        """
        #TODO: select inversePrecision properly
        #it should be possible to analytically invert B-spline, but I will not go there
        #find position x such that self.model(x)==y
        
        #check data type
        isNumber = isinstance(y, Number)
        if isinstance(y, list): y = np.array(y)
        isNumberArray = isinstance(y, np.ndarray) and y.dtype.kind in "iuf"
        #note that we do not accept list
        if not (isNumber or isNumberArray):
            raise TypeError("Spline.inverse does not accept data of type", y.__class__)

        #check data value
        xmin, xmax = self._u[0], self._u[-1]
        ymin, ymax = self.model([xmin, xmax])
        if isNumber and not (y>=ymin and y<=ymax):
            raise ValueError("Spline.inverse: value", y, "out of spline range", ymin, ymax)
        if isNumberArray and not np.logical_and(y>=ymin, y<=ymax).all():
            raise ValueError("Spline.inverse: value range", y.min(), y.max(), "out of spline range", ymin, ymax)

        #find inverse (assuming mostly linear spline)
        if isNumberArray:
            xmin = np.ones_like(y)*xmin
            xmax = np.ones_like(y)*xmax
            ymin = np.ones_like(y)*ymin
            ymax = np.ones_like(y)*ymax
        
            while True: #TODO: I am slightly nervous about the infinite loop
                x = xmin + (y-ymin)*(xmax-xmin)/(ymax-ymin) 
                yy = self.model(x)#this should not discard any x
                #yy == y: we are done
                #yy < y: y0 = yy
                #yy > y: y1 = yy
                
                if (np.abs(yy-y)<inversePrecision).all(): 
                    #print("Spline.inverse end", x)
                    return x
                i = np.where(yy<y)
                ymin[i] = yy[i]
                xmin[i] = x[i]
                
                i = np.where(yy>y)
                ymax[i] = yy[i]
                xmax[i] = x[i]
            #TODO: it is not efficient calculating model on those positions where yy ~ y
            pass
        else:
            while True: #TODO: I am slightly nervous about the infinite loop
                #print((xmin, xmax), (ymin, ymax))
                x = xmin + (y-ymin)*(xmax-xmin)/(ymax-ymin)
                yy = self.model(x)#this should not discard any x
                
                #for scalars
                if np.abs(yy-y)<inversePrecision: 
                    #print("Spline.inverse end", x)
                    return x
                if yy<y:     
                    ymin = yy
                    xmin = x
                if yy>y:
                    ymax = yy
                    xmax = x
                pass
            pass
        pass

    def knots(self):
        return self._u, self.model(self._u)
        
    def limits(self):
        #Todo: this will crash if spline is not ready yet
        return [self._u[0], self._u[-1]]

    def clearIntermediates(self):
        #remove (memory heavy) intermediate structures
        #keep only _a, _P, _p, _u - this allows to calculate spline at arbitrary x
        del self._B
        del self._BB
        del self._DD
        del self._I
        del self._x
        del self._y
        #_n and _l are not needed, but also are not memory heavy

        self._x = None #data
        self._y = None #data
        
        #structures
        self._B = None #basis functions
        self._BB = None
        self._DD = None
        self._I = None
        pass

    #only save the minimum needed information: _a, _u, _p, _P
    #save/loaded spline cannot be changed without setting it up from scratch (apart from scalar scaling __imul__)
    def save(self, filepath):
        #use a binary pickle
        import pickle
        with open(filepath, "wb") as f:
            data = (self._p, self._P, self._u, self._a)
            pickle.dump(data, f) #note that in Python3 this will use protocol version 3 by default which is not compatible with Python 2.x
        pass
        
    def load(self, filepath):
        import pickle
        with open(filepath, "rb") as f:
            self.clearIntermediates()
            self._p, self._P, self._u, self._a = pickle.load(f)
        pass

    def saveData(self):
        return (self._p, self._P, self._u, self._a)
        
    def loadData(self, data):
        self.clearIntermediates()
        self._p, self._P, self._u, self._a = data
    pass
