
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


from PyQt5 import QtCore

class UIGroup(QtCore.QObject):
    """similar to QActionGroup or QButtonGroup
    but aggregating UI elements 
    and triggering on any element change
    but only if all elements are set in UI (i.e. were added marked as set or were changed at least once)
    """
    #todo: 
    # keep track of signals and disconnect them on close
    # implement removing elements
    
    #problem
    # Qt UI maji vetsinou signaly na zmenu hodnoty
    # pokud ale je tady oaznaicm jako unset
    # a pak nastavim hodnotu na tu co uz maji, tak by se tady mely oznacit jako set, ale signal zmeny hodnoty nebude vygenerovany
    
    changed = QtCore.pyqtSignal()
    
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self._members = {}
        pass
    
    def add(self, element, signal, isSet=False):
        #print(dir(signal))
        #element = signal.__self__ #signal is assumed to be a bound method of UI element
        #well, pyqtSignal does not have __self__
        self._members[element] = isSet
        signal.connect(self._onChanged)
        pass
        
    def unsetAll(self):
        for it in self._members: self._members[it] = False
        pass
        
    def setAll(self):
        for it in self._members: self._members[it] = True
        pass
        
    def isAllSet(self):
        return len(self._members)==0 or all(self._members.values())
        
    def _onChanged(self, *args):
        #*args will be data of change signal, ignored here
        element = self.sender()
        self._members[element] = True
        if all(self._members.values()): self.changed.emit()
        pass
    pass    
    
