"""
PQR package, Qt-based plotting widgets 
Copyright (C) 2019  Jan Alster 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""

#~ #Todo?
#~ we can have defaultParams part of class definition
  #~ -specifying which attributes are aprams and should be part of automatic save/load 
  #~ -and their default values
  #~ - base class init needs ability to collect defaultParams from whole descendant chain
  #~ - base class init will create simple attribute params (just setattr(self, name, defaultvalue))
#~ derived class init might redefine the param to be property or external property
#~ as it is now, params and their default values are instance dependend

#we can separate params (attributes with automatic save/load) and properties
# we can treat params as class defined (with default values)
#  and maybe provide some convinience for defining properties, but is should be straightforward to that without help
#  possibly some interface for Qt signals and params groups

"""
#small problem with getattr
#   setattr sets the value of property
#   getattr returns the getter, not the value

that is actually worse
 - that was a bug in my code, trying to set property on instance, but it has to be done on class
     i.e.   addParam(name, value, Class.getter, Class.setter)   and inside setattr(Class, name, property(Class.getter, Class.setter)
 - but I cannot make an external params as property
 because I want to tie *instance* of UI elment (like QComboBox) which is not even created or known on Param-derived class declaration to one of the class attributes
 
 I guess easiest is to make (automatically or manually if needed) getter and setter for everything and make params accessible only through them - not properties
 
 i.e.
 addParam(name, defValue, value)
 will
  - create self.__param_name
  - create getter self.name
  - create setter self.setName
  - will be basically just a way how to not write trivial getters and setter by hand
    - register params for auto save/load under name with created getter/setter

 addParam(name, defValue, value, getter, setter)
  - will not create internal _param_name or getter/setter
  - will register params under name with user supplied getter/setter
  - params value will be stored somewhere arbitrary, where getter and setter can get them
  - getter/setter can be self methods, or methods of external objects
"""

"""
now I have a design issue
params were originaly intended for easy access to things like Plotter line width or color
  - which are things that maybe should not be saved, because those might be hardcoded
  
saveload elements were designed for UI elements, which are inherently not hardcoded and their state cannot be retrieved after program termination

but there are things like contours colors which belong to both groups

I might add a option specifieng if params should be included for save/load (which might change as I expose more things to UI)
 to have the same interface for setting which are and which are not exposed to UI
 
there is also possibility that params do not describe reality completely - e.g. I might want to (or allow to) pass arbitrary QPen for plot line visuals, which will noy be completely described by width, color, pattern (I do not want to recreate Qt interface, but I also want to make it easy to specify "normal" pens for lines without need of creating QPen)
 - however such a pen will not be described in save/load, although its width might be changed later by UI

"""


"""
another thing is preferences to set params 
 individual setter - which are needed
 setup
 
 it is likely crazy to want from individual setter (which might belong to completely unrelated classes) to call setup
  - I might do a wrapper to do that
  - but it will be more readable to have relavant pieces of code in per param methods than in one big blob
  
 however setup might (if possible) block generation of signal, or raise some flags, to prevent multiple redraws when setting more params at once - .e.g if we want to set line width and color of pen, it is only needed to be generated once
 
 we can do param group setters? setup could aggregate all params of one group and send it to group setter
 individual setter would just forward it to group setter
 
 but that can be inplemented later

"""
"""       
    we have to create the register before creating the attributes / params
      - will it work if it is a class attribute?
      - likely have to call base class init directly, before creating / registering params
      - we need to call Qt base class (multiple inheritance) init directly to create UI elements getters/setters
      - Params base class cannot initialize params values from init **kw or call setup (even descendant setup), because at that point param are likely not yet included - that will be handled inside derived class init after call to base class init (unles both register and params are class attributes defined from class definition - which needs at minimum decorator params which will not work for dynamically loaded UI)
      - params can be initialized to init kw value during cretion (just pass kw as kw and not **kw)
      - problem is that we cannot separate kw params for this class and for ancestor class, because during call of ancestor init we do not yet know this class params - ancestor init needs to be able to handle arbitrary number of unknown kw 
      and we have no way how to check that all kw were used/handled
       - we can require that all params will have kw explicitly in init
    
    setup is covenient interface to set multiple values at once (prefereably without emiting signals twice)  
    
    we might even have group of params with their own notification signals 
        (original SaveLoadElement included uninitialized params and test of groups if all params in group were initialized)
    but that was not good, it is probably better to have a default value

    if Qt derived we can also have value change notification signals (but that could be handled by custom setter)
    """
        #there is problem with order - I assume that params will be created in derived class
        # however the initialization (i.e. from __init__ **kw) will be done in base class
        # then if we call super().__init__ before creating the params, they will not be initialized from **kw because they do not exists yet
        # but it might be inconvenient (esp. with multiple inheritance) to call super().__init__ after params definition in derived class __init__
        # it would be best to create params during derived class definition (not __init__)
        # or somehow declare them for base class init to be able to defin them 
        #    which is the way with class attribute paramsDefaults, but that cannot specify setters easily :(
        #    and those should be active before base class init initialization from **kw
        #   (I do not want to initialize params twice - once defaults and second time from base class init **kw)
        # the way would be to use decorator properties in derived class, but that is not convenient syntax
        # and it cannot easily create "simple" attributes


    #~ paramsDefaults = {}
    #TODO?: not sure if this can handle inheriting the params properly
    # we can use super().paramsDefaults to iterate through inheritance
    # we probably should also call super().setup in setup (esp. if reimplemented in derived class)
    #~ question is if super().__init__ can see original (descendant) class paramsDefaults
    #~ and we cannot use this interface to change param to property
    
    #TODO?: params should be set only through setup
    # but we might do a property that would call self.setup(param=value) upon self.param=value
    # or setup might be a convinient way how to change multiple parameters without several signal emits

Unset = object() #unique value for testing if param value is set
from collections import OrderedDict
class Param(object):
    __slots__=("get", "set", "preSet", "onSet", "onGet", "group")
    def __init__(self, get=None, set=None, preSet=None, onSet=None, onGet=None, group=None):
        self.get = get
        self.set = set
        self.preSet = preSet
        self.onSet = onSet
        self.onGet = onGet
        self.group = group
        #~ self.value = value
    pass

class ParamGroup(object):
    __slots__=("onSet", "params", "counter")
    def __init__(self, onSet=None):
        self.onSet = onSet
        self.params = []
        self.counter = 0
    pass

class Params(object):
    Unset = Unset
    """
    unified params handling
    - param is a property-like instance attribute (managed via custom Params.__getattr__ and Params.__setattr__)
    - can have either internal data storage, or external data storage accessible via getter and setter
    - allows defining response callback onGet and onSet
    - optional automatic save and load 
    - (todo) groups of params with group onGet and onSet
        - when setting values via global setter Params.setParams, group response will be triggered only once even 
        - (??) might even do some timer triggers to only trigger group response only after some time after last single params setter of the group is invoked (so you do not have to use global setter and still only get one group response when setting individual group params in close succession)
    
    - getter should work with no arguments, setter with just one - the new value - for global set and retrieve to work, individual setter and getters can take whatever parameters
    
    Interface
    register new param
    - addParam(name, defaultValue=Unset, getter=Unset, setter=Unset, onGet=Unset, onSet=Unset, save=False)
        - name - has to be valid python indentifier and must not start with "_".
          Also it might overwrite _name so do not use that simultaneously with param name. Access to the params is property-like, i.e. self.name for both reading and writing.
        - value will be set to param upon creation (Params unique object Unset is used to test if value is specified as it can be None). You should give some value, unless it is stored on external object (and Params does not have to allocate any data storage).
        - getter / setter  - interface for data storage
             - if both are given, they take complete control of data storage
             - if only one is given, assume internal data storage self._name, the one given should use that, but might modify its value
             - if none is given, use data storage on self.__params
        - onGet / onSet    - response for getter or setter, cannot modify stored value directly, but can invoke getter / setter (beware of infinite recursion though)
        - save  - will be automaticaly included in Params.saveSettings and Params.loadSettings
    Typically, only things that cannot be reconstructed from code (typically UI interface) should be saved and loaded (saveload is optional for that reason).
    
    global getter / setter
     - params(names) - to retrieve {name:value, ...} for all params
     - setParams({name:value, ...} - to set specified params (TODO: "as one step" in case these lead to signals)
     
    automatic saving and loading
      - saveSettings
      - loadSettings
      Can be overloaded in derived classes to include other things. Call Params version to include params.
   """

    """
    there are actually two separate concepts that can be solved using getter setter and we need to cover both
     - data storage: can be external and we need getter and setter to access this storage (i.e. storing parameters on UI elements - it is good if instance attribute is connected to UI element and the most certain way how to accomplish this is if the attribute takes its value from UI element, no need for synchronization of two data storages; on the other hand acess time is probably slower, so it might not be good for often queried parameters)
     - response to data change (could be also e.g. access counter on getter) - best ensured if this is a part of setter/getter, but that is not compatible with external data storage, because external object does not care about response of this instance
     
     we can do more attributes of params: getter, setter (these two will handle data storage), onGet, onSet (these two will handle response to data change and cannot directly change the data storage)
    
    
    we will still have problems if e.g. instance uses setter to check param values (as it might need to change the value) and it wants to use external data storage (because that also need a setter) and I might want to install the external storage from outside of the instance
     - but that would need to have 
         - callback before setting to check the value and possibly modify
         - setter
         - callback after setting the value (which cannot modify it anymore) to react to the change
    """
    
    """
    another thing
     - when change UI, i.e. load items to QComboBox, the attached param attribute has changed value, but it will not trigger onSet or preSet, because the new value was not input through the attribute, but externaly 
      - i.e. if I change value of the param in the external data storage in another way, the param attribute will not know it :(
      - I still need to trigger change response manually
    since that is the case, wouldn't it be better to have internal data storage for each param and just keep them synchronized?
     - esp UI elements will be changed in another way (ie. by user) - that is their purpose
    
    """

    """
    for now, params should not be instances, but only values - something that can survive pickle dump and load cycle
    """
    
    """
    we might need a way to set parameters without triggering onSet, or at least delaying it
    use case: we have several parameters which are all needed to calculate something (possibly all have the same onSet, but not necesarily)
    if we set the first one, the others are not set yet, they do not even exist (so callback cannot test their values) and we have problems
    
    it is the same for setup
     
    basically we need a way to delay calling callbacks (and aggregate callbacks if called more than once during the waiting time)
    
    #TODO: another problem is that linesH will create lines, but not position them
        # initially right after creation we cannot position the lines because separationH is not defined
        # but after it is and we change number of lines, they should be position as well and not wait for other impulses
        
        #possibly we can create groups in another way: it will be defined after all parameters are set and will trigger if any is changed (obviously cannot be triggered when some parameter is not set)
        # or we can define conditional onSet, e.g. we can set the color of lines on change of color, but only if some other parameters are set
        #but it is crazy to do all that just to facilitate __init__ (once all params are set, there are not problems)
        
        # only how to set color on ancestor class and have it processed in descendant only after it is ready?
        
        # we can leave params values unset in __init__ and call (overloaded) setup at __init__ end?
        # i.e. in Desc.__init__(**kw_params_setup) call super().__init__() without params kw, that will just create params but not assign value
        # (i.e. Anc.__init__ will call setup, but without any values and it will not do anything)
        # Decs.__init__ can call setup after it is ready using all params kw (even those for Anc) and its setup will then pass them to Anc.setup 
            
    """

    def __init__(self, *args, **kw):
        self.__params = OrderedDict()
        self.__paramsToSave = []
        self.__groups = OrderedDict()
        self.__flagSetup = False
        #~ print(**kw)
        if kw.pop("startSetup", False):
            self.startSetup()
            
        if len(kw)>0:
            print("Params.__init__: WARNING likely unused kw", kw, "for", self)
        super().__init__(*args, **kw)  #send the rest to object, it should complain if there are unused kw
        pass
    
    #between calls to startSetup and stopSetup invocations of onSet will be withheld and only processed (each only once, in original order but with the last one value) with stopSetup
    def startSetup(self):
        self.__flagSetup = True
        self.__callbackQueue = OrderedDict()
    
    def stopSetup(self):
        self.__flagSetup = False
        for it in self.__groups.values():
            if it.counter>0:
                it.counter = 0
                if it.onSet is not None: it.onSet()
        #only callback the last one
        for it in self.__callbackQueue:
            it(self.__callbackQueue[it])
        del self.__callbackQueue
    
    def __getattr__(self, name):
        #~ print("Params.__getattr__", name)
        if not name.startswith("_") and name in self.__params:
            p = self.__params[name]
            #there is no preGet
            if p.get is not None:
                value = p.get()
            else:
                #assuming internal data storage
                value = self.__getattribute__("_"+name+"Param")
            
            if p.onGet is not None: p.onGet(value)
            return value
        raise AttributeError("Params.__getattr__:", name, "is not a param of", self)
    
    def __setattr__(self, name, value):
        if not name.startswith("_") and hasattr(self, "_Params__params") and name in self.__params:
            #~ print("Params.__setattr__", name, value)
            p = self.__params[name]
            
            if p.get is None and p.set is None and p.onGet is None and p.onSet is None and p.preSet is None:
                #params is not managed at all (apart from possible automatic saving) and can be changed to regular attribute
                return object.__setattr__(self, name, value)
            
            if p.preSet is not None: value = p.preSet(value) #give chance to change the value before set
            if p.set is not None:
                p.set(value) #this should not change the value, just store it
            else:
                #we operate on internal data storage
                self.__setattr__("_"+name+"Param", value)
            #Todo: we might store the old value and pass it to onSet    
            self._afterSet(name)
            #~ if p.onSet is not None: p.onSet(value) #allow response to value change    
            #~ if p.group is not None: 
                #~ if self.__flagSetup:
                    #~ self.__groups[p.group].counter += 1
                #~ else:
                    #~ self.__groups[p.group].onSet()
            return 
        return object.__setattr__(self, name, value)
    
    def _afterSet(self, name, value=Unset):
        #this is here as a response for notification that param stored in external data storage (e.g. UI element) was changed bypassing the property-like attribute settings (which is the function of UI after all)
        #but we want to trigger after setter responses as if the value was set on attribute
        p = self.__params[name]
        if value is Unset: value = getattr(self, name) #this is for response for UI, we need to read the new value
        if p.onSet is not None: 
            if self.__flagSetup:
                self.__callbackQueue[p.onSet] = value
            else:
                p.onSet(value) #allow response to value change    
        if p.group is not None: 
            if self.__flagSetup:
                self.__groups[p.group].counter += 1
            else:
                self.__groups[p.group].onSet()
        pass
    
    def __delattr__(self, name):
        if not name.startswith("_") and name in self.__params:
            p = self.__params[name]
            print("Params.__delattr__", name, p)
            if p.get is None and p.set is None and p.onGet is None and p.onSet is None and p.preSet is None:
                #params was released as normal attribute
                object.__delattr__(self, name)
            elif not (p.get is not None and p.set is not None):
                #we have used internal data storage
                object.__delattr__(self, "_"+name+"Param")
            
            del self.__params[name]
            return
        return object.__delattr__(self, name)
    
    def addParamGroup(self, name, onSet=None):
        self.__groups[name] = ParamGroup(onSet)
    
    def addParam(self, name, value=Unset, getter=Unset, setter=Unset, onGet=Unset, onSet=Unset, preSet=Unset, group=Unset, save=False, changeSignal=None):
        #Note that using value=Unset and trying to read value of param without external data storage (or otherwise initialized) will result in exception as the attribute is only created when some value is given.
        
        #Todo: we might distinguish between getter or setter passed as None and Unset
        #  one can do automatic handling, other can do read-only (or set-only??) param
        #  which does not make sense, because params should change
        #  but possibly they might be changed by something external
        
        #actually if we have no getter or setter name can live directly in self
        #but we have to remove it if it is redefined
        #and we might want to keep original callbacks if we redefine params (set it to None to overwrite the old ones with nothing)
        if name in self.__params:
            if value is Unset: #if param is redefined, try to keep its value
                value = getattr(self, name)
            p = self.__params[name]
            if getter is Unset: getter = p.getter
            if setter is Unset: setter = p.setter
            if onGet is Unset: onGet = p.onGet
            if onSet is Unset: onSet = p.onSet
            if preSet is Unset: preSet = p.preSet
            if group is Unset: group = p.group
            delattr(self, name)
        
        if getter is Unset: getter = None
        if setter is Unset: setter = None
        if onGet is Unset: onGet = None
        if onSet is Unset: onSet = None
        if preSet is Unset: preSet = None
        if group is Unset: group = None
        
        self.__params[name] = Param(getter, setter, preSet, onSet, onGet, group)#getter, setter, onGet, onSet
        if group is not None: self.__groups[group].params.append(name)
        if value is not Unset: setattr(self, name, value)
        if save: self.__paramsToSave.append(name)
        if changeSignal is not None: changeSignal.connect(lambda *args: print("changeSignal", name, args) or self._afterSet(name, *args)) #this is Qt "ready"
        pass
        

    #global setter
    def setParams(self, **kw):
        #NOTE that order of setting params is random
        #only change present parameters, leave other intact
        #TODO: (setup should raise global flag so that only one of each change signal is triggered) - see setup
        for name in kw: 
            if name in self.__params: setattr(self, name, kw[name])
        pass

    #test
    def setup(self, **kw):
        self.startSetup()
        self.setParams(**kw)
        self.stopSetup()

    #getter
    def params(self, keys=None):
        if keys is None: keys = self.__params.keys()
        res = OrderedDict()
        for name in keys: res[name] = getattr(self, name) 
        return res
    
    def saveSettings(self):
        """
        overload this to record params or other settings of members in derived classes
        compatibility with SaveLoadElement
        """
        #TODO: we have a problem here - we cannot save instances as those cannot be recovered when loaded
        #   typically QColor()
        #   Param might need to have separate getter for saving instances (i.e. if getter returns an instance, we will need saveGetter(instance) that will give something serializable which setter (or loadSetter) can digest
        #    Param is getting quite complex ...
        return self.params(self.__paramsToSave)
        
    def loadSettings(self, settings):
        """
        overload this to load params or other settings of members in derived classes
        compatibility with SaveLoadElement
        """
        print("Params.loadSettings", self, settings)
        #self.setParams(**settings)
        #keep order of params in settings (if it is ordered sequence)
        for name in settings: 
            if name in self.__params: setattr(self, name, settings[name])
        
        pass
    pass

if __name__=="__main__":
    
    #Param use showcase
    class Helper:
        def __init__(self):
            self._value = None
            
        def set(self, value):
            print("Helper.set", self, self._value, value)
            self._value = value
        
        def get(self):
            print("Helper.get", self, self._value)
            return self._value

    class ShowcaseParams(Params):
        def __init__(self, p1=1, p2=2, p3=5, notparam=None, **kw):
            super().__init__(**kw)
            self.addParam("p1", p1) #simple param
            self.addParam("p2", p2, self._getP2, self._setP2) # property param
            self._p3object = Helper()
            self.addParam("p3", p3, self._p3object.get, self._p3object.set, save=True) #external object param
            
            self._notparam = notparam
            pass
            
        def _setP2(self, value):
            print("ShowcaseParams._setP2", self, value)
            self._p2 = value
            
        def _getP2(self):
            print("ShowcaseParams._getP2", self, self._p2)
            return self._p2
            
    class Level2(ShowcaseParams):
        def __init__(self, o1=1, p1=5, **kw):
            super().__init__(**kw)
            #~ self.addParam("p1", p1, self._getP1, self._setP1)
            self.addParam("o1", o1)
            
        def _setP1(self, value):
            print("Level2._setP1", self, self._p1, value)
            self._p1 = value
            
        def _getP1(self):
            print("Level2._getP1", self, self._p1)
            return self._p1
            
    l2 = Level2()
    print("\t save", l2.saveSettings())
    print("\t load", l2.loadSettings({"o1":5, "p1":8, "p2":4, "p3":45}))
    print("\t params", l2.params())
    print()
    print(dir(l2))
    print("o1", l2.o1)
    
    #test if setattr uses property setter - it does
    #~ class OurClass:

        #~ def __init__(self, a):
            #~ print("init")
            #~ self.OurAtt = a

        #~ @property
        #~ def OurAtt(self):
            #~ print("getter")
            #~ return self.__OurAtt

        #~ @OurAtt.setter
        #~ def OurAtt(self, val):
            #~ print("setter")
            #~ if val < 0:
                #~ self.__OurAtt = 0
            #~ elif val > 1000:
                #~ self.__OurAtt = 1000
            #~ else:
                #~ self.__OurAtt = val

        #~ def test(self, value):
            #~ print("test")
            #~ setattr(self, "OurAtt", value)

    #~ x = OurClass(10)
    #~ print(x.OurAtt)
    #~ x.test(15)
    #~ x.OurAtt = 10
    #~ print("getattr", getattr(x, "OurAtt"))

    #test if super() init uses descendant method - it does
    #~ class A:
        #~ def __init__(self):
            #~ print("A.init")
            #~ self.m()
        #~ def m(self):
            #~ print("A.m")
            
    #~ class B(A):
        #~ def m(self):
            #~ print("B.m")
            #~ super().m()
            
    #~ B()

    if False:
        #can function remove elements from kw?
        def test2(d):
            if "a" in d: 
                print("removing a", d["a"])
                del d["a"]
            
        def test(**kw):
            print("pred", kw)
            test2(kw)
            print("po", kw)
            
        test(a=5, b=4)
    
    if False:
        #Todo: include this if params are to be specified as class {}
        #test if base class init sees class attributes of derived class
        
        class A(object):
            a = {"A":1}
            def __init__(self):
                print("A.__init__", self, self.a)
                print("\t\t aggregate", self.get())
                
            def get(self):
                d = {}
                print(self.__class__.__mro__)
                for it in self.__class__.__mro__:
                    try:
                        d.update(it.a)
                    except AttributeError:
                        print("class", it ,"is not part of Params chain")
                return d
        
        class B(A):
            #~ a = {"B":2}
            pass
        
        
        class C(B):
            a = {"C":3}    

        class D(C):
            a = {"D":4}    
        
        b = B()
        c = C()
        d = D()
                    
    
