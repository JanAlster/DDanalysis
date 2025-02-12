
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

#seems like this needs another overhaul
# - we need constant covertors (that cannot be modified)
# - we need memory convertors, that do not depend on file
# - we need way how to store convertors
#
# for example
# - in place simple range axis for plot needs a convertor (that does nothing basically)
# - convertors are to be stored inside data, then we do not have a file, only byte stream to load (and it cannot be modified later)
# - if convertor contains a custom class as a part of params, i.e. Spline in DS convertors and that is pickled
#      if we then change the class definition or move the class somewhere else in the file hiearchy
#      unpickle will/can fail
#
# using generic data container, we do not know which axitype/convertor will be used apriory, so we need to pickle the object and go with that in 
#    saving the generic data container
#
# it would be better to save the convertors into the data (wasting some space, but allowing transition to other places)
#    (as long as the classes are available though)
#    there is likely no way how to make a generic convertor describing an arbitrary transformation
#
# -> if convertor is not tied to a file, we need it to be immutable itself (or at least you can only add units)
#       so we cannot load/setup multiple times
#       if you want to change, you need to have a new convertor object, the data files will of course hold the original one


#TODO: in case of mising hardware and not available calibration files, even default ones are *not* created and DD crashes when trying to save
# with Convertor._filepath not set up
# - possibly not an issue with constant Convertor modification (not set _filepath ignores saving attempts and adding units)

#we can do a Convertor class that would basically be a convertor (redesign - position and convertor were originally just for DS)
#it would work on top of a convertor data file (each instance connected to a real file) so that there is a guarantee that the file exists
#   - the file could even be in another session, but I would probably prefere to copy it over in that case, so that session are self contained
#it would provide uniform interface for using calibrations in datafiles
# - filename
# - base unit
# - units and conversions (but it does *not* have a selected unit) - or rather parameters to create them


#we can save some parameters for each unit, but the derived classes have to specify how the conversion should look like given those parameters
# - seems like it will have to be hardcoded what units the Convertor can handle
# - you need to be backward compatible when changing code here

# - demand that all convertors are valid on predefined range (basically any value we can get from the experiment)


#there is a tricky part - wedge calibrations come in pairs and it is nonsense to use one for one wedge and different one for the other wedge
#   but each axis have its own convertor (and the idea was to save the parametes (in this case the splines) into convertors data files) 
# - we will have to leave it to the user not to be stupid about it ...
# - on the other hand if we want to say that some axis is e.g. a time delay - that would depend on *both* wedges positions (as the delay is the difference between them)
# - even more rDS should align all four beams (or just three as one distance is fixed)
# - anyway we should save position of all DSs even if only one is scanned
# => this is for rDS to handle


#to use this
#create Convertor and set it up using either
# - load(filepath) - will connect Convertor to existing filepath and load from the file (error if filepath does not contain at least base units)
# - setup(baseUnit, limit, filepath) - setup required elements, if filepath exists, it will be loaded
# afterward you can add new units (will be immediately saved to filepath), but you cannot modify already loaded ones
# you can reset the Convertor by calling load or setup again

#adding new units is split to two steps (adding a group with some parameters and adding individual units with conversion functions) because
#   conversion functions cannot be easily saved to calibration files - paramets can

import os
import os.path as op
import pickle
#TODO: use ordered dict for params and values (to keep order of units)
#TODO: maybe prevent overriding calibration files at all (always create a new calibration file, even when adding new units)
class Convertor(object):
    """
    Usage:
     - first create the convertor "convertor = Convertor()" (or use derived class)
     - initialize the convertor with either "convertor.load(filepath)" or "convertor.setup(baseunit, limits, filepath)"
     - (optionally, but doesn't make much sense not to do this) add units based on some parameters "convertor.addUnits(unitgroup, parameters)" - only derived classes know how to handle their own units
     - then it can be used to convert values to/from baseunits and added units
     
    Uninitialized convertor will behave as with only no-name base unit without limits (i.e. directly return the value)
    Note that convertor is bound to the file specified during initiallization and each unit group can only be added once. You have to reinitiallize the convertor with another file if you want to change things.
    """
    #reimplement in derived classes        
    def _addUnits(self, group, params):
        #reimplement in derived classes and call super()._addUnits() to catch base unit
        #you can set up more units from one "group" of parameters
        #example
        #unit2base = lambda x : x * params[0]
        #base2unit = lambda x : x / params[0]
        #self._addUnit(group+"0", unit2base, base2unit)
        #unit2base = lambda x : x * params[1]
        #base2unit = lambda x : x / params[1]
        #self._addUnit(group+"1", unit2base, base2unit)
        
        if group == "base":
            self._params["base"] = params
            base2base = lambda x : x
            self._addUnit(params["unit"], base2base, base2base, params["label"]) #we need this for QConvertor
        pass

    #outer interface
    #setup
    def setup(self, baseUnit, limits, filepath, label=""):
        #assert filepath valid and does not exists
        print("Convertor.setup", baseUnit, limits, filepath)
        print(os.path.exists(filepath), os.access(os.path.dirname(filepath), os.W_OK), os.path.dirname(filepath))
        assert((not os.path.exists(filepath)) and os.access(os.path.dirname(filepath), os.W_OK))
            
        self._filepath = filepath
        #reset
        self._reset()
        #setup basics
        self.addUnits("base", {"unit":baseUnit, "limits":limits, "label":label})
        pass
        
    def addUnits(self, group, params):
        print("Convertor.addUnits", group, params, self._filepath is None)
        #NOTE: do not use "base" as group
        if self._filepath is None: 
            print("Convertor.addUnits on constant covertor is ignored")
            return
        assert(group not in self._params) #do not allow changing of existing units
        self._params[group] = params
        self._save() #if filename is not setup by this point there will be problems
        self._addUnits(group, params)
        pass
    
    def baseUnit(self):
        return None if "base" not in self._params else self._params["base"]["unit"]
        
    def limits(self):
        return None if "base" not in self._params else self._params["base"]["limits"]
    
    def _isBaseUnit(self, unit):
        #check for aliases of baseUnits
        return unit is None or unit == "" or ("base" in self._params and unit == self._params["base"]["unit"])
        
    #TODO: this might come in handy if some new units has tighter limits than base, but it has to change calfile, so the correct way is probably to reset the convertor
    #def resetLimits(self, limits=None):
    #    self._params["base"]["limits"]=limits
        
    
    #~ def units(self, all=False):
        #~ if all:
            #~ return list(self._units.keys())+[self.baseUnit()]
        #~ else:
            #~ return list(self._units.keys())
        #~ pass
        
    def units(self):
        return list(self._units.keys())
        
    def label(self, unit):
        label = self._units[unit][3]
        if label is None:
            return "" if "base" not in self._params else self._params["base"]["label"]
        return label
        
    #conversion
    def V2U(self, valueBase, unit, check=False):
        if check:
            coerced, valueBase = self._check(valueBase, self.limits())
        
        if self._isBaseUnit(unit): 
            res = valueBase
        else:
            res = self._units[unit][1](valueBase)

        return res if not check else (coerced, res)

    def U2V(self, valueUnit, unit, check=False):
        if self._isBaseUnit(unit):
            if check:
                return self._check(valueUnit, self.limits())
            else:
                return valueUnit
        else:
            if check:
                coerced, valueUnit = self._check(valueUnit, self._units[unit][2])
            res = self._units[unit][0](valueUnit)
            return res if not check else (coerced, res)
        pass

    #SaveLoadElement        
    #this can be used as getter for SaveLoadElement
    def filepath(self):
        return self._filepath

    #this can be used as setter for SaveLoadElement
    def load(self, filepath, constant=False):
        try:    
            #in case filepath is opened file-like object
            params = pickle.load(filepath)
            self._filepath = None #in this case constant has to be True so we can ignore call "constant" option
        except ImportError:
            #TODO: old data saved as a part of some pickled files custom classes, we have to be able to find definitions
            # - this is a backward compatibility hack (which will not work on another PC)
            import sys
            sys.path.append(r"C:\Users\DD\Desktop\code\DD_v02")
            params = pickle.load(f)
        except:
            #filepath is normal string
            #assert filepath valid and exists - exception if not (not our problem)
            with open(filepath, "rb") as f:
                self._filepath = filepath if not constant else None
                print("Convertor.load", filepath)
                try:
                    params = pickle.load(f)
                except ImportError:
                    #TODO: old data saved as a part of some pickled files custom classes, we have to be able to find definitions
                    # - this is a backward compatibility hack (which will not work on another PC)
                    import sys
                    sys.path.append(r"C:\Users\DD\Desktop\code\DD_v02")
                    params = pickle.load(f)
                pass
            pass
        
        baseParams = params.pop("base")

        #reset - only if we actually have the baseParams
        self._reset()

        #"base" first
        #backward compatibility hack
        if "label" not in baseParams: baseParams["label"]=""

        #skip save() in addUnits
        self._params["base"] = baseParams
        self._addUnits("base", baseParams)

        for group in params.keys():
            p = params[group]
            self._params[group] = p
            self._addUnits(group, p)
        pass        

    #inner interface
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self._params = {}
        self._units = {}
        pass
    
    def _reset(self):
        self._params = {}
        self._units = {}
        pass    
    
    def _save(self):
        with open(self._filepath, "wb") as f:
            pickle.dump(self._params, f)
        pass
        
    def _saveHack(self):
        #NOTE: do not use this, it is intended for internal testing
        return pickle.dumps(self._params)

    def _addUnit(self, unit, unit2base, base2unit, label=None):
        #print("Convertor._addUnit", unit, unit2base, base2unit)
        #check that unit is not already present - that is a coding error - two parameter groups try to setup the same unit
        assert(unit not in self._units)
        limits = self.limits()
        unitLimits = None if limits is None else tuple(sorted((base2unit(limits[0]), base2unit(limits[1])))) #note that unit conversion might change axis direction (i.e. lower limit of base value can be upper limit of unit value)
        self._units[unit]=(unit2base, base2unit, unitLimits, label)
        pass
        
    #we can do without numpy here, but it will not be fast
    @staticmethod
    def _check(value, limits):
        #print("Convertor._check", value, limits)
        if limits is None: return False, value
        #import warnings
        try: 
            if (min(value)<limits[0] or max(value)>limits[1]):
                #warnings.warn("convertor: value ({}, {}) coerced to convertor limits ({}, {})".format(min(value), max(value), *limits))
                for i in range(len(value)):
                    value[i] = limits[0] if value[i]<limits[0] else limits[1] if value[i]>limits[1] else value[i]
                return True, value
            else:
                return False, value
        except TypeError:
            #not iterable
            if (value<limits[0] or value>limits[1]):
                #warnings.warn("convertor: value {} coerced to convertor limits ({}, {})".format(value, *limits))
                value = limits[0] if value<limits[0] else limits[1]
                return True, value
            else:
                return False, value
        pass
    pass

#TODO: we could generalize this to convertor with specified precision in base units (that should work well with precision in DoubleWidget)
#TODO: test this    
def ConvertorWithIntBase(Convertor):
    #conversions to base units are always rounded off to int
    #conversion
    import math
    def V2U(self, valueBase, unit, check=False):
        valueBase = round(valueBase,0)
        if check:
            coerced, ttt = self._check(valueBase, self.limits())
            if coerced:
                valueBase = math.ceil(ttt) if ttt>valueBase else math.floor(ttt)
            else:
                valueBase = ttt
        
        if self._isBaseUnit(unit): 
            res = valueBase
        else:
            res = self._units[unit][1](valueBase)

        return res if not check else (coerced, res)

    #Todo: tohle by slo asi zjednodusit
    def U2V(self, valueUnit, unit, check=False):
        if self._isBaseUnit(unit):
            valueBase = round(valueUnit,0)
            if check:
                coerced, ttt = self._check(valueBase, self.limits())
                if coerced:
                    valueBase = math.ceil(ttt) if ttt>valueBase else math.floor(ttt)
                else:
                    valueBase = ttt
                return (coerced, valueBase)
            else:
                return valueBase
        else:
            if check:
                coerced, valueUnit = self._check(valueUnit, self._units[unit][2])
            res = self._units[unit][0](valueUnit)
            valueBase = round(res,0)
            if check:
                coerced, ttt = self._check(valueBase, self.limits())
                if coerced:
                    valueBase = math.ceil(ttt) if ttt>valueBase else math.floor(ttt)
                else:
                    valueBase = ttt
                return (coerced, valueBase)
            else:
                return valueBase
        pass
    pass
    
class Quantity(object):
    """Simple class to store a value that can be converted to different units without changing.
    Otherwise just converting the value there and back might change it due to rounding off errors.
    """
    #Note that if your convertors can handle it, value can be either scalar or vector
    def __init__(self, convertor, value=None):
        #print("Quantity.__init__", convertor, value)
        super().__init__()
        self._convertor = convertor
        if value is not None:
            self.setValue(value) #use this for range checking
        else:
            self.coerced = False #flag indicating if last set value has been coerced to convertor limits
            self._value = None #no value set, no conversion possible
        pass
        
    def setConvertor(self, convertor): #do not abuse this - like changing the base units, etc.
        #print("Quantity.setConvertor", convertor)
        self._convertor = convertor
        self.setValue(self._value) #to check limits on new convertor
        pass
        
    def value(self, unit=None):
        #print("Quantity.value", unit)
        if unit is None or self._value is None:
            return self._value
        else:
            return self._convertor.V2U(self._value, unit)
        pass
        
    def setValue(self, value, unit=None):
        #print("Quantity.setValue", value, unit)
        self.coerced, self._value = self._convertor.U2V(value, unit, True) #go through convertor for range checking (even for base units)
        pass
    pass    
