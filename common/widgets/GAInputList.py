
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
from PyQt5 import QtWidgets, QtCore
#from PyQt5.uic import loadUiType

from .GrowingList import GrowingList

import os.path as op
#widgetUI, widgetClass = loadUiType(op.join(op.dirname(__file__),'GAItemWidget.ui'))

#class GAItemWidget(widgetClass, widgetUI):
class GAItemWidget(QtWidgets.QWidget):
    SomeValue=QtCore.pyqtSignal()
    NoValue=QtCore.pyqtSignal()
    
    def __init__(self, parent=None):
        super().__init__(parent)
        #self.setupUi(self)
        layout = QtWidgets.QHBoxLayout()
        self.dLifetime = QtWidgets.QDoubleSpinBox()
        self.dLifetime.setMaximum(999999999.000000)
        layout.addWidget(self.dLifetime)
        self.bLifetimeFix = QtWidgets.QCheckBox()
        layout.addWidget(self.bLifetimeFix)
        self.cbFrequencyType = QtWidgets.QComboBox()
        self.cbFrequencyType.addItems(["Ø", "+", "-", "±"])
        layout.addWidget(self.cbFrequencyType)
        self.dFrequency = QtWidgets.QDoubleSpinBox()
        self.dFrequency.setMaximum(999999999.000000)
        layout.addWidget(self.dFrequency)
        self.bFrequencyFix = QtWidgets.QCheckBox()
        layout.addWidget(self.bFrequencyFix)
        layout.setContentsMargins(0,0,0,0)
        self.setLayout(layout)
        

        self.dLifetime.editingFinished.connect(self.testValue)
        self.dFrequency.editingFinished.connect(self.testValue)
        self.bFrequencyFix.stateChanged.connect(self.testValue)
        self.bLifetimeFix.stateChanged.connect(self.testValue)
        self.cbFrequencyType.currentIndexChanged.connect(self.testValue)
        pass

    def testValue(self, *args):
        value=self.dLifetime.value()
        if value==0. : self.NoValue.emit()
        else: self.SomeValue.emit()

    def value(self):
        return self.dLifetime.value(), self.bLifetimeFix.isChecked(), self.cbFrequencyType.currentText(), self.dFrequency.value(), self.bFrequencyFix.isChecked()

    def setValue(self, lifetime, lFix=None, freqType=None, freq=None, fFix=None):
        try:
            lifetime, lFix, freqType, freq, fFix=lifetime #in case lifetime is tuple
        except: pass
        self.dLifetime.setValue(lifetime)
        self.bLifetimeFix.setChecked(lFix)
        self.cbFrequencyType.setCurrentText(freqType)#TODO: hopefully this will not allow text not present in combobox
        self.dFrequency.setValue(freq)
        self.bFrequencyFix.setChecked(fFix)
        pass
    pass

class GAInputList(QtWidgets.QScrollArea):
    valueChanged=QtCore.pyqtSignal()
    def __init__(self, *args):
        super(GAInputList, self).__init__(*args)
        widget=GrowingList(GAItemWidget, GAItemWidget.value, GAItemWidget.setValue, self)
        widget.widgetAdded.connect(self.ensureWidgetVisible)
        widget.valueChanged.connect(self.valueChanged)
        self.setWidget(widget)
        pass

    def value(self):
        return self.widget().value("value") #value is a getter of GAItemWidget

    def setValue(self, data, emit=True):
        widget=self.widget()
        widget.valueChanged.disconnect(self.valueChanged)
        ret=self.widget().setValue(data, "setValue") #setValue is a setter of GAItemWidget
        widget.valueChanged.connect(self.valueChanged)
        if emit: self.valueChanged.emit()
        return ret
    pass
    
