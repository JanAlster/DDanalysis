
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

from PyQt5 import QtWidgets, QtGui, QtCore

try:
    from common.conversion import Convertor, Quantity

    #TODO: position does not check value type, this can cause problems. be carefull what you send there
    #TODO: Qt uses locale decimal separator, float() does not, we should handle both ',' and '.'
    from common.widgets.DoubleWidget import DoubleWidget
except:
    from ..conversion import Convertor, Quantity
    from .DoubleWidget import DoubleWidget
    
#TODO: we could catch exceptions in case some convertor goes crazy
#TODO: changing convertor is not very polished (I only use it for setting convertors to UI form elements)

#different approach to DoubleQuantityEdit
# - assume that baseUnits are integer
# - we cannot set arbitrary value
# - during setting value from UI convert to base units and round off to nearest integer
# - change value to that corresponding to the rounded off base value
# - more general - simply change value to correspond to the base value (whatever the convertor decides to adjust base value to)
# - set sensible precision

class DoubleQuantityEdit(DoubleWidget):
    def __init__(self, convertorSelector, value=None, *args, **kw):
        editValue = kw.pop('editValue', True)
        super().__init__(*args, **kw)
        assert(convertorSelector is not None)
        
        self._quantity = Quantity(convertorSelector, value)
        convertorSelector.currentUnitChanged.connect(self._setCurrentUnit)
        self.valueChanged.connect(self._updateMasterValue) #DoubleWidget signal, we want to listen to it when the change comes from DoubleWidget._editor, not from DoubleQuantityEdit.setValue
        
        if not editValue:
            self.setReadOnly(True)
            
        self._setCurrentUnit(convertorSelector.currentUnit())
        pass
        
    def setConvertorSelector(self, convertorSelector):
        assert(convertorSelector is not None)
        self._quantity._convertor.currentUnitChanged.disconnect(self._setCurrentUnit)
        self._quantity.setConvertor(convertorSelector)
        convertorSelector.currentUnitChanged.connect(self._setCurrentUnit)
        self._setCurrentUnit(convertorSelector.currentUnit())
        pass
    
    def _updateMasterValue(self, value):
        #~ #print("DoubleQuantityEdit._updateMasterValue", value, self._currentUnit)
        #~ self._quantity.setValue(value, self._currentUnit)
        #~ #print("\tcoerced", self._quantity.coerced)
        #~ self._updateCurrentValue() #always update current value
        #~ if self._quantity.coerced:
            #~ #also change current value - master value was coerced
            #~ self.setStyleSheet("QLineEdit{background : rgba(255, 170, 0, 33%)}")
        #~ else:
            #~ self.setStyleSheet("")
        #~ pass
        self.setValue(value, self._currentUnit)
        pass
        
    def _updateCurrentValue(self):
        #print("DoubleQuantityEdit._updateCurrentValue start", self._value, self._currentUnit)
        self.valueChanged.disconnect(self._updateMasterValue)
        #problem is how to select precision
        # - under the assumption that convertor operates on integer base units, it should be easy to select
        # - however edit does not know that
        # - so precision should be given by convertor
        #TODO: using precision==4 and contract==True works nicely for DD user case, but is not general
        DoubleWidget.setValue(self, self.value(self._currentUnit), 4, contract=True)
        self.valueChanged.connect(self._updateMasterValue)
        #print("DoubleQuantityEdit._updateCurrentValue end")
        pass

    def _setCurrentUnit(self, unit):
        #print("DoubleQuantityEdit._setCurrentUnit", unit)
        self._currentUnit = unit
        self._updateCurrentValue()
        pass

    #we need to redefine value and setValue to behave as position
    def value(self, unit=None):
        return self._quantity.value(unit)
        
    def currentValue(self):
        return self.value(self._currentUnit) #should be the same as DoubleWidget.value(self) up to precision of display
        
    def setValue(self, value, unit=None):
        #print("DoubleQuantityEdit.setValue", value, unit)
        self._quantity.setValue(value, unit)
        if self._quantity.coerced:
            self.setStyleSheet("QLineEdit{background : rgba(255, 170, 0, 33%)}")
        else:
            self.setStyleSheet("")
        self._updateCurrentValue()
        #print("DoubleQuantityEdit.setValue end")
    pass


#TODO: this concept has a problem
# - I want to keep convertors with data (also wedge calibration dialog), but there is no need for UI element there
#   - data need access to convertors to be properly able to handle covnersion of data, not just axes
# - I also need to select units of these convertors once we show the data
#   - UI element should be part of the plot widget
# - but mostly it is nuisance to keep convertor and selector separate - need to propagate two elements from resources instead of one
#   - the would require to design selector that has all the outer interface used by other parts (and preferably none of convertor itself duplicated) - might work, if it only concerns current units
#   - convertor would have to notify of new units, which it cannot do since it is not Qt
#   - so such selector has no way how to tell new units were loaded to convertor
# - or I can make the data container a Qt object (e.g. Dialog) with selection buttons (and possibly pop it up/include in plot widget)
#   - the question is then how to select units in method calls - we probably need to interpret None as current units not as base units
#   - problem is that such Qt object cannot be easily unplugged from the rest of UI
#   - for this we need separate convertor and selector that can be easily connected/disconnected
#~ class ConvertorSelector(Convertor, QtWidgets.QComboBox):
    #~ currentUnitChanged = QtCore.pyqtSignal(str)
    
    #~ def __init__(self, *args, **kw):
        #~ QtWidgets.QComboBox.__init__(self, *args, **kw)
        #~ Convertor.__init__(self)
        #~ #TODO: not sure if super().__init__(*args, **kw) would work
        #~ self.setSizeAdjustPolicy(self.AdjustToContents)

        #~ self.currentIndexChanged[str].connect(self.currentUnitChanged)

        #~ self.setEnabled(False)
        #~ pass   

    #~ def _addUnit(self, unit, unit2base, base2unit):
        #~ super()._addUnit(unit, unit2base, base2unit)
        #~ self.setEnabled(True)
        #~ self.addItem(unit)
        #~ pass
    
    #~ def _reset(self):
        #~ self.currentIndexChanged[str].disconnect(self.currentUnitChanged)
        #~ self.clear()
        #~ self.currentIndexChanged[str].connect(self.currentUnitChanged)
        #~ super()._reset()
        #~ pass

    #~ def setCurrentUnit(self, unit):
        #~ if unit in self.units():
            #~ self.setCurrentText(unit) #TODO: will this trigger the signal?
        #~ else:
            #~ print("ConvertorSelector.selectUnit: requested unit", unit, "is not available")
        #~ pass

    #~ def currentUnit(self):
        #~ return self.currentText() #None or "" will be interpreted as base unit by Convertor
        #~ #unit = self.currentText()
        #~ #if unit == "": return None
        #~ #return unit

    #~ def C2V(self, data, check=False):
        #~ return self.U2V(data, self.currentUnit(), check)
        
    #~ def V2C(self, data, check=False):
        #~ return self.V2U(data, self.currentUnit(), check)
    #~ pass

#~ #problem with this approach is duplicity of SelectorFactory and StaticSelector
#~ # - every change has to be make twice
#~ # - resulting object are different (although they have the same "selector" interface)

class SelectorFactory(QtWidgets.QComboBox):
    #use for making Convertor to ConvertorSelector
    #class ConvertorSelector(SelectorFactory, Convertor):
    #   pass
    
    currentUnitChanged = QtCore.pyqtSignal(str)
    unitsChanged = QtCore.pyqtSignal()
    
    def __init__(self, *args, **kw):
        QtWidgets.QComboBox.__init__(self, *args, **kw)
        self.setSizeAdjustPolicy(self.AdjustToContents)

        self.currentIndexChanged[str].connect(self.currentUnitChanged)

        self.setEnabled(False)
        pass   

    def _addUnit(self, unit, unit2base, base2unit, label=None):
        #print("SelectorFactory._addUnit", unit, unit2base, base2unit)
        super()._addUnit(unit, unit2base, base2unit, label)
        self.setEnabled(True)
        self.addItem(unit)
        self.unitsChanged.emit()
        pass
    
    def _reset(self):
        self.currentIndexChanged[str].disconnect(self.currentUnitChanged)
        self.clear()
        self.currentIndexChanged[str].connect(self.currentUnitChanged)
        super()._reset()
        self.unitsChanged.emit()
        pass

    def setCurrentUnit(self, unit):
        if unit in self.units():
            self.setCurrentText(unit) #TODO: will this trigger the signal?
        else:
            print("ConvertorSelector.selectUnit: requested unit", unit, "is not available")
        pass

    def currentUnit(self):
        return self.currentText() #None or "" will be interpreted as base unit by Convertor
        #unit = self.currentText()
        #if unit == "": return None
        #return unit

    def C2V(self, data, check=False):
        return self.U2V(data, self.currentUnit(), check)
        
    def V2C(self, data, check=False):
        return self.V2U(data, self.currentUnit(), check)
        
    def copy(self):
        #Todo: this will work only if __bases__ = (SelectorFactory, ArbitraryConvertor)
        try:
            assert(self.__class__.__bases__[0] == SelectorFactory)
            n = self.__class__.__bases__[1]()
            n._params.update(self._params)
            n._units.update(self._units)
            return n
        except:
            print("SelectorFactory.copy: cannot remove Selector from copied instance")
            import traceback
            print(traceback.print_exc())
            #as a fall back return just a copy
            return super().copy()
        pass
        
    pass


class StaticSelector(QtWidgets.QComboBox):
    currentUnitChanged = QtCore.pyqtSignal(str)
    #only shows the units the convertor had at the time it was set and it cannot react to unit changes
    #so only use it for convertors that are finished and will not change
    def __init__(self, *args, **kw):
        QtWidgets.QComboBox.__init__(self, *args, **kw)
        self.setSizeAdjustPolicy(self.AdjustToContents)
        self.currentIndexChanged[str].connect(self.currentUnitChanged)
        self.setEnabled(False)
        self.convertor = None
        pass   

    def setConvertor(self, convertor, currentUnit=None):
        #it should be possible to setConvertor multiple times
        self.convertor = convertor
        self.currentIndexChanged[str].disconnect(self.currentUnitChanged)
        self.clear()
        self.addItems(convertor.units())
        self.currentIndexChanged[str].connect(self.currentUnitChanged)
        if currentUnit is None: currentUnit = convertor.baseUnit()
        self.setCurrentUnit(currentUnit)
        self.setEnabled(True)
        pass

    def setCurrentUnit(self, unit):
        if self.convertor is None: return
        if unit in self.convertor.units():
            self.setCurrentText(unit) #TODO: will this trigger the signal?
        else:
            print("StaticSelector.selectUnit: requested unit", unit, "is not available")
        pass

    def currentUnit(self):
        if self.convertor is None: return
        return self.currentText() #None or "" will be interpreted as base unit by Convertor
        #unit = self.currentText()
        #if unit == "": return None
        #return unit

    def currentLabel(self):
        if self.convertor is None: return
        return self.convertor.label(self.currentUnit())

    def C2V(self, data, check=False):
        if self.convertor is None: return
        return self.convertor.U2V(data, self.currentUnit(), check)
        
    def V2C(self, data, check=False):
        if self.convertor is None: return #TODO: not sure if this is a good idea, possibly allow to fail, but then we need some empty convertor for NDContainer.RangeAxis
        return self.convertor.V2U(data, self.currentUnit(), check)
    pass
    
    
"""
class QConvertor(Convertor, QtCore.QObject):
    unitAdded = QtCore.pyqtSignal(str)
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        pass
        
    def _addUnit(self, unit, unit2base, base2unit):
        super()._addUnit(unit, unit2base, base2unit)
        self.unitAdded.emit(unit)
        pass
    pass

class UnitSelector(QtWidgets.QComboBox):
    #unit selector for convertor
    unitChanged = QtCore.pyqtSignal(str)
    
    def __init__(self, *args, convertor=None, **kw):
        #it can happen that the first argument is parent of the QComboBox (automatic ui loader does this)
        super().__init__(*args, **kw)
        self.currentIndexChanged[str].connect(self.unitChanged)
        self.setConvertor(convertor)
        self.unitChanged.connect(self.slot)
        pass
        
    def slot(self, *args):
        print("UnitSelector.slot:", *args)
        pass
        
    def setConvertor(self, convertor):
        self._convertor = convertor
        self.currentIndexChanged[str].disconnect(self.unitChanged)     
        current = self.currentText()   
        print("UnitSelector.setConvertor", current, convertor)
        self.clear()
        self.currentIndexChanged[str].connect(self.unitChanged)
        
        if self._convertor is None:
            self.setEnabled(False)
        else:
            baseUnit = self._convertor.baseUnit()
            if baseUnit == None: 
                self.setEnabled(False)
            else:
                self.setEnabled(True)
                self.addItem(baseUnit)
                self.addItems(self._convertor.units())
                print("\tprevious", current, current in self._convertor.units())

            self._convertor.unitAdded.connect(self._updateConvertor)
            pass
        pass
    
    def convertor(self):
        return self._convertor
    
    def _updateConvertor(self, newUnit):
        self.setEnabled(True)
        self.addItem(newUnit)
        pass
        
    def selectUnit(self, unit):
        #TODO: or just use setCurrentText directly (it should change current index 
        #   I am just not sure how it reacts when text is not in list (probably does nothing)
        if self._convertor is not None and unit in self._convertor.units(all=True):
            self.setCurrentText(unit) #TODO: will this trigger the signal?
        else:
            print("UnitSelector.selectUnit: requested unit", unit, "is not available")
        pass
    pass
"""    
