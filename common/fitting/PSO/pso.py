
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
#copyright 2015 Jan Alster
#based on pyswarm.pso
#which is Copyright: 2013 Abraham Lee, released under BSD license (or any other the author approves, just ask!)
#nicmene to je nejspis primy prepis algorithmu z wikipedie

#TODO: initialize the random seed   
#TODO: fitovani vahovych parametru pro danou modelovou funkci   
   
import numpy as np

import threading

class PSO(object):
    def __init__(self, func, lb, ub, ieqcons=[], f_ieqcons=None, args=(), kwargs={}, 
        swarmsize=None, omega=0.5, phip=0.5, phig=0.5, maxiter=1000, 
        minstep=1e-8, minfunc=1e-8, debug=False):
        """
        Perform a particle swarm optimization (PSO)
       
        Parameters
        ==========
        func : function
            The function to be minimized
        lb : array
            The lower bounds of the design variable(s)
        ub : array
            The upper bounds of the design variable(s)
       
        Optional
        ========
        ieqcons : list
            A list of functions of length n such that ieqcons[j](x,*args) >= 0.0 in 
            a successfully optimized problem (Default: [])
        f_ieqcons : function
            Returns a 1-D array in which each element must be greater or equal 
            to 0.0 in a successfully optimized problem. If f_ieqcons is specified, 
            ieqcons is ignored (Default: None)
        args : tuple
            Additional arguments passed to objective and constraint functions
            (Default: empty tuple)
        kwargs : dict
            Additional keyword arguments passed to objective and constraint 
            functions (Default: empty dict)
        swarmsize : int
            The number of particles in the swarm (Default: 10)
        omega : scalar
            Particle's inertia weight (Default: 0.5)
        phip : scalar
            Particle's cognitive weight (Default: 0.5)
        phig : scalar
            Swarm's social weight (Default: 0.5)
        maxiter : int
            The maximum number of iterations for the swarm to search (Default: 1000)
        minstep : scalar
            The minimum stepsize of swarm's best position before the search
            terminates (Default: 1e-8)
        minfunc : scalar
            The minimum change of swarm's best objective value before the search
            terminates (Default: 1e-8)
        debug : boolean
            If True, progress statements will be displayed every iteration
            (Default: False)
       
        """
       
        assert len(lb)==len(ub), 'Lower- and upper-bounds must be the same length'
        assert hasattr(func, '__call__'), 'Invalid function handle'
        lb = np.array(lb)
        ub = np.array(ub)
        assert np.all(ub>lb), 'All upper-bound values must be greater than lower-bound values'
       
        vhigh = np.abs(ub - lb)
        vlow = -vhigh
        
        # Check for constraint function(s) #########################################
        obj = lambda x: func(x, *args, **kwargs)
        if f_ieqcons is None:
            if not len(ieqcons):
                if debug:
                    print('No constraints given.')
                cons = lambda x: np.array([0])
            else:
                if debug:
                    print('Converting ieqcons to a single constraint function')
                cons = lambda x: np.array([y(x, *args, **kwargs) for y in ieqcons])
        else:
            if debug:
                print('Single constraint function given in f_ieqcons')
            cons = lambda x: np.array(f_ieqcons(x, *args, **kwargs))
        self.cons=cons            
            
        # Initialize the particle swarm ############################################
        D = len(lb)  # the number of dimensions each particle has
        S = swarmsize if swarmsize is not None else 5*D #automate swarm size if not specified
        x = np.random.rand(S, D)  # particle positions
        v = np.zeros_like(x)  # particle velocities
        p = np.zeros_like(x)  # best particle positions
        fp = np.zeros(S)  # best particle function values
        g = []  # best swarm position
        fg = 1e100  # artificial best swarm position starting value
        
        for i in range(S):
            # Initialize the particle's position
            x[i, :] = lb + x[i, :]*(ub - lb)
       
            # Initialize the particle's best known position
            p[i, :] = x[i, :]
           
            # Calculate the objective's value at the current particle's
            fp[i] = obj(p[i, :])
           
            # If the current particle's position is better than the swarm's,
            # update the best swarm position
            if fp[i]<fg and self.is_feasible(p[i, :]):
                fg = fp[i]
                g = p[i, :].copy()
           
            # Initialize the particle's velocity
            v[i, :] = vlow + np.random.rand(D)*(vhigh - vlow)
            pass
            
        print("PSO.__init__ best from random positions: g", g, "fg", fg)
            
        self.it = 1
        
        self.D=D
        self.S=S
        
        self.v=v
        self.omega=omega
        self.phip=phip
        self.p=p
        self.x=x
        self.phig=phig
        self.fg=fg
        self.g=g
        self.fp=fp
        
        self.lb=lb
        self.ub=ub
        self.obj=obj
        
        def objThread(x, storage, index):
            storage[index] = obj(x)
        self.objThread = objThread
        
        self.maxiter=maxiter
        self.minfunc=minfunc
        self.minstep=minstep
        
        self.debug=debug

        pass

    def set_particle(self, i, x):
        assert(self.is_feasible(x))
        
        #TODO: ensure that we do not replace the best particle - or if so that we update best position
        self.x[i, :] = x
        fx = self.obj(x)
        
        print("PSO: setting particle", x, "with goodness", fx)
        
        #check if it is a better position
        if fx < self.fp[i]:
            self.p[i,:]=x
            self.fp[i]=fx
            
            #check if it is a global best
            if fx < self.fg:
                self.fg=fx
                self.g=x.copy()
                pass
            pass
        pass
       

    def is_feasible(self, x):
        check = np.all(self.cons(x)>=0)
        return check

        
    def run(self, steps=None, maxiter=None):
        #TODO: run specific number of steps, override maxiter
        """    
        Returns
        =======
        g : array
            The swarm's best known position (optimal design)
        f : scalar
            The objective value at ``g``
            
        """
        #TODO: maybe make it local here
        while self.it<=self.maxiter:
            dx, dfx = self.step() #could be None if no particle improved

            #TODO: here is slight conceptual problem
            # - if we start from fairly good position (e.g. repeating already performed fit)
            #   then we might not be able to update to better position
            #   and we will get only None, None from self.step
            #   so we cannot stop the iteration short of running out of steps
            #   but the whole fit is pointless and a waste of effort
            # - if we cannot move, we cannot say if we are close enough to minimum to satisfy stopping criteria
            #   even though they are in fact satisfied the whole time

            # Compare swarm's best position to current particle's position
            # (Can only get here if constraints are satisfied)
            if dfx!=None:
                if self.debug:
                    print('New best for swarm at iteration %d:'%self.it, self.g, self.fg)

                if self.minfunc!=None and dfx<=self.minfunc:
                    if self.debug:
                        print('Stopping search: Swarm best objective change less than:', self.minfunc)
                    return self.g, self.fg
                elif self.minstep!=None and dx<=self.minstep:
                    if self.debug:
                        print('Stopping search: Swarm best position change less than:', self.minstep)
                    return self.g, self.fg
                pass
                
            if self.debug:
                print('Best after iteration %d'%self.it, self.g, self.fg)
                
            # ~ print("Iteration", self.it, self.fg)    
            self.it += 1
            pass
            
        if self.debug:
            print('Stopping search: maximum iterations reached -->', self.maxiter)
        
        if self.g is []:
            print('Warning: No feasible point found')
        return self.g, self.fg
    
    
    def step(self):
        """
        Makes a step with each particle.
        
         Returns
        =======
        sx : array
            The swarm's best improved position in this step (or None if no particle improved)
        sfx : scalar
            The objective value at ``g`` (or None if no particle improved)
        """
        
        rp = np.random.uniform(size=(self.S, self.D))
        rg = np.random.uniform(size=(self.S, self.D))
        #Todo: np.random.rand(S, D)
        
        sfx = None
        sx = None
        
        #Todo: this might be done not in loop perhaps...
        #~ for i in range(self.S):
            #~ # Update the particle's velocity
            #~ #TODO: this does not work if we have only one reference point
            #~ #   or possibly if func gives unexpected return value (shape)
            #~ self.v[i, :] = self.omega*self.v[i, :] + self.phip*rp[i, :]*(self.p[i, :] - self.x[i, :]) + \
                      #~ self.phig*rg[i, :]*(self.g - self.x[i, :])
            #~ # Update the particle's position, correcting lower and upper bound 
            #~ # violations, then update the objective function value
            #~ x = self.x[i, :] + self.v[i, :]
            #~ mark1 = x<self.lb
            #~ mark2 = x>self.ub
            #~ x[mark1] = self.lb[mark1]
            #~ x[mark2] = self.ub[mark2]
            #~ self.x[i, :] = x
            #~ fx = self.obj(x) #this is typically a bottleneck, note that my objective function saves already calculated values, so asking again should be fast if some particle stops moving
            #~ 
            #~ # Compare particle's best position (if constraints are satisfied)
            #~ if fx<self.fp[i] and self.is_feasible(self.x[i, :]):
                #~ self.p[i, :] = self.x[i, :].copy()
                #~ self.fp[i] = fx
                #~ 
                #~ if sfx==None or fx<sfx: 
                    #~ sfx=fx
                    #~ sx=self.x[i,:].copy()
                #~ pass
            #~ pass
        
        # Update the particle's velocity
        self.v[:] = self.omega*self.v + self.phip*rp*(self.p - self.x) + self.phig*rg*(self.g - self.x)
        # Update the particle's position, correcting lower and upper bound 
        # violations, then update the objective function value
        fx = np.empty((self.S,))
        T = []
        for i in range(self.S):    
            x = self.x[i, :] + self.v[i, :]
            mark1 = x<self.lb
            mark2 = x>self.ub
            x[mark1] = self.lb[mark1]
            x[mark2] = self.ub[mark2]
            self.x[i, :] = x
            #fx = self.obj(x) #this is typically a bottleneck, note that my objective function saves already calculated values, so asking again should be fast if some particle stops moving
            T.append(threading.Thread(target=self.objThread, args=(x, fx, i)))
            T[-1].start()
        
        for i in range(self.S):
            T[i].join()
        
        for i in range(self.S):    
            # Compare particle's best position (if constraints are satisfied)
            if fx[i]<self.fp[i] and self.is_feasible(self.x[i, :]):
                self.p[i, :] = self.x[i, :].copy()
                self.fp[i] = fx[i]
                
                if sfx is None or fx[i]<sfx: 
                    sfx=fx[i]
                    sx=self.x[i,:].copy()
                pass
            pass
        
        
        #update global best
        if sfx is not None and sfx<self.fg: #there was an improvement
            dx = np.sqrt(np.sum((self.g-sx)**2))
            dfx = np.abs(self.fg - sfx)

            self.g=sx
            self.fg=sfx

            return dx, dfx #return step size

        return None, None #no improvement
    pass

def metafit(func, lbound, ubounds, params=None, average=1, maxiter=100):
    #try to find good omega, phip, phig to fit func(params)
    #i.e. we know the solution and try to find weights so that swarm arives to the solution
    #alternatively use (several) random params from between lbounds, ubounds
    if params==None:
        #choose params randomly (average over several choices)
        raise NotImplementedError
    
    #is it too cheeky to use PSO for metafit?
    #is it too cheeky to use metafit to find weights for PSO used in metafit?
    
    def metafunc(weights, lbounds=lbound, ubounds=ubounds, params=params):
        omega, phip, phig = weights
        ret=np.zeros(average)
        for i in range(average):
#            if params==None:
#                solution=lbounds + np.random.rand(D)*(ubounds - lbounds)
#            else:
#                solution=params
            swarm=PSO(func, lbounds, ubounds, omega=omega, phip=phip, phig=phig)
            swarm.run()
            #check the results
            #i.e. closeness of swarm.g to params, and swarm.it
            #create fitness from this
            #TODO: add number of iterations to fitness (i.e. speed of fit)
            #Todo: add result spread to fitness (i.e. stability of fit) - uz je castecne zahrnuta diky prumerovani
            ret[i] = np.sum((swarm.g-params)**2)**0.5/len(params)
            pass
        ret=ret.mean()
        print("metafunc", ret, *weights)
        return ret
        
    #TODO: maybe the bounds can be different?
    metaubounds=[1,1,1.]
    metalbounds=[0,0,0.]
    metaswarm=PSO(metafunc, metalbounds, metaubounds, omega=0.9, phip=0.8, phig=0.2, maxiter=maxiter) 
    #get the weights for fitting func with PSO
    return metaswarm.run()



def test_metafit():
    def f(param, x):
        return param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2)       
        
    #select params
    solution=[5.3, 23.1, 3.4]
    x=np.arange(-5, 10, .1)
    y=f(solution, x)
    
    lbound=[x[0], min(y), 1e-10]
    ubound=[x[-1], max(y), abs(x[-1]-x[0])]
    
    #Todo: this could be done by metafit
    def evaluate(param):
        return np.sum((f(param, x)-y)**2)/(len(y)-1.)
    
    mf = metafit(evaluate, lbound, ubound, solution, average=10, maxiter=500)
    print("metafit", mf)
    #metafit (array([ 0.69230046,  0.31876719,  1.        ]), 1.6810983333638489e-07)  - maxiter==1000, no averaging
    #metafit (array([ 0.53522652,  1.        ,  1.        ]), 3.4220394985592628e-06) - maxiter 100, average 3
    #metafit (array([ 0.89963513,  0.40104498,  0.65064116]), 2.8699523877615324e-06) - maxiter 100, average 3
    #metafit (array([ 0.58855188,  1.        ,  1.        ]), 2.3289163598623693e-06) - maxiter 100, average 3
    #jde videt, ze je to hodne nestabilni, zkusim average=5, maxiter=500
    #metafit (array([ 0.66813029,  0.57323875,  0.95387921]), 3.6852337795412415e-06)
    #metafit (array([ 0.7682925 ,  0.45187318,  0.96342436]), 3.6972703535702845e-06)
    #metafit (array([ 0.74110786,  0.        ,  0.97681726]), 7.8853406209880581e-06)
    #metafit (array([ 0.79866493,  0.32857793,  1.        ]), 3.3411647463476432e-06)
    #zkusim average=10, maxiter=500
    #metafit (array([ 0.65643923,  0.24458033,  0.96691738]), 5.6107063106889458e-06)
    #metafit (array([ 0.83317526,  0.14090742,  0.7079887 ]), 5.9109478640175721e-06)
    #metafit (array([ 0.6878308 ,  0.84994993,  1.        ]), 5.8684500193699068e-06)
    #metafit (array([ 0.75545307,  0.29480582,  0.79363519]), 6.5135745598214495e-06)
    #tohle neni moc stabilni...
    pass

def test(maxfev=100):
    dimension = 3
    p0=[10., 5., 3.]
    print("real params", p0)
    x=np.arange(-5, 25, .1)
    y=p0[1]*np.exp(-2.7725887222397811*((x-p0[0])/p0[2])**2)
    
    
    lbound=[x[0], min(y), 1e-10]
    ubound=[x[-1], max(y), abs(x[-1]-x[0])]
    print("lb", lbound)
    print("ub", ubound)
    
    def evaluate(param):
        gauss=param[1]*np.exp(-2.7725887222397811*((x-param[0])/param[2])**2)
        return np.sum((gauss-y)**2)/(len(gauss)-1.)
        
        
    amplitude=y.max()-y.min()
    #bool(fitness <= 0.0001*amplitude**2)
    
    #Initialize libPyPopot
    #algo = Popot.Stochastic_PSO_2006(dimension, lbound, ubound, stop, evaluate)
    swarm=PSO(evaluate, lbound, ubound, swarmsize=12, omega=0.9, phip=0.8, phig=0.2, maxiter=maxfev,         minstep=1e-8, minfunc=1e-8, debug=False) #hodne zalezi na vahovych parametrech


    import time
    from fitting import fit_gaussian

    t0=time.clock()
    f, ff = swarm.run()
    t1=time.clock()
    pf=fit_gaussian(x,y, maxfev)
    t2=time.clock()
    print('pso', t1-t0, 'popot', t2-t1)

    print()
    print("fit", f, ff, "\n", f-p0, "\n", swarm.it, np.sum((f-p0)**2)**0.5/len(p0))
    print("popot", pf[1], pf[3], "\n", np.array(pf[1])-p0, "\n", np.sum((np.array(pf[1])-p0)**2)**0.5/len(p0))
    import pylab as pl
    pl.plot(x,y, marker='x')
    pl.plot(x, f[1]*np.exp(-2.7725887222397811*((x-f[0])/f[2])**2))
    pl.plot(x, pf[0])
    pl.show()

if __name__=='__main__':
    #test(500)
    test_metafit()

