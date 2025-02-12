
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

#TODO: consider nesting of SaveLoadElements

class SaveLoadElement(object):
    def __init__(self, *args):
        #print("SaveLoadElement.__init__", *args)
        super().__init__(*args)
        self._editables = []
        self._groups = {}
        self._callbacks = {}
        self._disabledGroups = set()
        self._set = {}
        #TODO: this could add each instance to an application wide list of cofigurable elements (so far I have to keep track in main application)
        pass
        
    def addEditable(self, name, getter, setter, connectToChanged=None, group=None, initialize=False):
        self._editables.append((name, getter, setter, group))
        if connectToChanged is not None:
            connectToChanged.connect(lambda *args: self.changed(*args, name=name, group=group))
        if group is not None:
            if group not in self._groups:
                self.addGroup(group)
            self._groups[group].append(name)
        self._set[name] = initialize #i.e. if initialize == True, current value of editable is deemed valid and set; if not, current value is not defined (this is a change after long time, I hope it works as intended :) - i.e. considering default values of ui elements as valid)
        pass

    #if only basic settings are applicable, addEditable is sufficient
    #if not, overload saveSettings and loadSettings, but do not forget to call these ones if you want to keep addEditable functionality
    def saveSettings(self, **kw):
        res={}
        for name, getter, setter, group in self._editables:
            res[name]=getter()
        return res
        
    def loadSettings(self, settings, **kw):
        #print("formTemplate.loadSettings", settings)
        #possibly change state here in derived classes
        for name, getter, setter, group in self._editables:
            #print("\t", name, setter)
            try:
                if name in settings: 
                    #TODO: this should call changed only once (setter might cause connectToChanged signal)
                    setter(settings[name])
                    self.changed(name=name, group=group)
            except Exception as e:
                print("SaveLoadElement cannot set", name," to value", settings[name])
                import traceback
                print(traceback.format_exc())
                pass
            pass
        pass

    def addGroup(self, name, changeCallback=None):
        #just for adding empty groups
        self._groups[name] = []
        if changeCallback is not None:
            self._callbacks[name] = changeCallback
        pass

    def disableGroup(self, group):
        """temporarily disable callbacks for given group"""
        self._disabledGroups.add(group)
    
    def enableGroup(self, group):
        """reenable previously disabled group"""
        if group in self._disabledGroups:
            self._disabledGroups.remove(group)

    def isSet(self, name):
        print(name, "is element", name in self._set, "is group", name in self._groups)
        if name in self._set: return self._set[name]
        print([(it, self._set[it]) for it in self._groups[name]])
        if name in self._groups: return len(self._groups[name])==0 or all(self._set[it] for it in self._groups[name])
        raise IndexError("Unknown SaveLoadElement", name)

    def unset(self, name=None):
        if name is None: 
            for it in self._set: self._set[it] = False
        else:
            self._set[name]=False
        pass

    #note that this does not work perfectly, some UI elements could be changed/set without triggering connectToChanged signal
    def changed(self, *args, name=None, group=None):
        #in case some derived class needs to react to changing editables
        #here we could e.g. save 'current' settings, but we might leave that to self._close
        #or mark self as edited/changed
        if name is not None:
            self._set[name] = True
            if group is not None and group in self._callbacks and group not in self._disabledGroups:
                if self.isSet(group):
                    self._callbacks[group]()
                pass
            pass
        pass        
    
    pass
