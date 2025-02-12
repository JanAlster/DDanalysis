#***********************************************************************
#DD package, data collection and analysis of 2D electronic spectra
#Copyright (C) 2019  Jan Alster (Charles Univesity, Prague)
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

#!/usr/bin/env python
# -*- coding: utf-8 -*-


#TODO: redesign tracer to be decorator, add tracer to mainPlot (or AddTracer to context menu of Data area)
#TODO: highlight item under mouse - give FigureScene top level item which will not be selectable by mouse and will be ignored by highlight (also can be send to top or bottom of z stack), and send it QPainterPath item.shape of item under mouse in mouse move envent, probably with timer so that it is not changed too often
#TODO: grid, diagonal grid, labels, ...


from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.uic import loadUi

import os
import os.path as op

import pickle

if __name__=="__main__":
    #create simple app containing just the dialog - preparation
    import sys
    sys.path.append(op.join(sys.path[0], "..")) #need to import from directory one above this file
    os.chdir("..")


from common.SaveLoadElement import SaveLoadElement
from common.interpolator import interpolate
#from common.UIGroup import UIGroup
from pqr import params
from pqr import pqr

import numpy as np

from collections import OrderedDict

def StyleSheet2QColor(stylesheet, key='color'):
    stylesheet=stylesheet.replace(' ','')
    i=stylesheet.find(key+':')+1+len(key)
    return QtGui.QColor(stylesheet[i:stylesheet.find(';', i)])

TRN2trn={'Total':'total', 'Rephasing':'rephasing', 'Non-Rephasing':'nonrephasing'}
trn2TRN={v: k for k, v in TRN2trn.items()}

PAIR={'Real': lambda x: x.real, 'Imag': lambda x: x.imag, 'Ampl': lambda x: np.abs(x), 'Phase': lambda x: np.angle(x), None:lambda x:x}

import pickle

class ControlBlob(params.Params, QtWidgets.QFrame): #Todo: possibly scroll area
    def __init__(self, plotMain, loadFromUI, *args, master=True):
        super().__init__(*args)
        self.plotMain = plotMain
        self._master = master

        #TODO: most of the things here are candidates for Params
        self.cbState = QtWidgets.QComboBox()
        self.cbFrame = QtWidgets.QComboBox()
        self.cbTRN = QtWidgets.QComboBox()
        self.cbTRN.addItems(["Total", "Rephasing", "Non-Rephasing"])
        self.cbPAIR = QtWidgets.QComboBox()
        self.cbPAIR.addItems(["Real", "Imag", "Ampl", "Phase"])
        self.iInterpolation1 = QtWidgets.QSpinBox()
        self.iInterpolation1.setMinimum(0)
        self.iInterpolation1.setMaximum(800)
        self.iInterpolation1.setSingleStep(100)
        self.iInterpolation3 = QtWidgets.QSpinBox()
        self.iInterpolation3.setMinimum(0)
        self.iInterpolation3.setMaximum(800)
        self.iInterpolation3.setSingleStep(100)
        self.bGlobalColors = QtWidgets.QCheckBox("Global")
        if master:
            self.pushAddContours = QtWidgets.QPushButton("+")
            self.pushAddContours.clicked.connect(self._addContours)
        else:
            self.pushAddContours = QtWidgets.QPushButton("-")
            self.pushAddContours.clicked.connect(self.parent()._removeContours)
        self.pushAddContours.setMaximumSize(20,20)

        #contour control
        self.bContoursVisible = QtWidgets.QCheckBox("Cont")
        self.pushContoursColor = QtWidgets.QPushButton() #color icon
        self.pushContoursColor.clicked.connect(self._setContoursColor)
        #~ px = QtGui.QPixmap(20,20)
        #~ px.fill(QtCore.Qt.black)
        #~ self.pushContoursColor.setIcon(QtGui.QIcon(px))
        self.dContoursMin = QtWidgets.QDoubleSpinBox()
        self.dContoursMax = QtWidgets.QDoubleSpinBox()
        self.iContoursLevels = QtWidgets.QSpinBox()
        self.iContoursLevels.setRange(2, 15)
        self.iContoursLevels.setValue(7)
        


        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(0,0,0,0)
        #dataset selection
        l0 = QtWidgets.QHBoxLayout()
        l0.addWidget(self.cbState, 1)
        l0.addWidget(self.cbFrame)
        l0.setContentsMargins(0,0,0,0)
        l0.addWidget(self.cbTRN)
        layout.addLayout(l0)
        
        #data selection and interpolation
        l1 = QtWidgets.QHBoxLayout()
        l1.setContentsMargins(0,0,0,0)
        l1.addWidget(self.cbPAIR)
        l1.addWidget(self.bGlobalColors)
        l1.addStretch(1)
        l1.addWidget(QtWidgets.QLabel("Interp"))
        l1.addWidget(self.iInterpolation1)
        l1.addWidget(self.iInterpolation3)
        l1.addWidget(self.bContoursVisible)
        l1.addWidget(self.pushContoursColor)
        l1.addWidget(self.dContoursMin)
        l1.addWidget(self.dContoursMax)
        l1.addWidget(self.iContoursLevels)
        l1.addWidget(self.pushAddContours)
        layout.addLayout(l1)
        
        self.setLayout(layout)
        
        p = self.plotMain.plot()
        if master:
            xaxis = p.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.bottom, label='ω₁ / cm⁻¹')
            yaxis = p.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.left, label='ω₃ / cm⁻¹')
            caxis = p.addAxis(None, pqr.AxisPosition.right, pqr.ColorAxisArea, globalZoom=False)
        else:
            xaxis = p.b1
            yaxis = p.l1
            caxis = p.r1
        self._ddPlotter = p.addDDPlotter(xaxis, yaxis, caxis, contours=True, colorPlot=master) #TODO: bilinear color scale, pqr.MZMTransform()) #TODO: I can do "typical" DD plot which uses these transforms and keeps track of correct stops and mnmx

        #Params right now cannot handle signals
        # esp group signal that would be perfect here
        #but they can do automatic saving
        #data selection
        self.addParam("state", getter=self.cbState.currentIndex, setter=self.cbState.setCurrentIndex, save=True, onSet=self.updateShow, changeSignal=self.cbState.currentIndexChanged)#we actually operate here with currentData, but that is not good for saving
        self.addParamGroup("dataBuffer", self._updateDataBuffer)
        self.addParam("TRN", getter=self.cbTRN.currentText, setter=self.cbTRN.setCurrentText, save=True, group="dataBuffer", changeSignal=self.cbTRN.currentTextChanged)
        self.addParam("PAIR", getter=self.cbPAIR.currentText, setter=self.cbPAIR.setCurrentText, save=True, group="dataBuffer", changeSignal=self.cbPAIR.currentTextChanged)
        self.addParam("interpolation1", getter=self.iInterpolation1.value, setter=self.iInterpolation1.setValue, save=True, group="dataBuffer", changeSignal=self.iInterpolation1.valueChanged)
        self.addParam("interpolation3", getter=self.iInterpolation3.value, setter=self.iInterpolation3.setValue, save=True, group="dataBuffer", changeSignal=self.iInterpolation3.valueChanged)
        #still data selection
        self.addParamGroup("plotMain", self.updatePlotMain)
        self.addParam("frame", getter=self.cbFrame.currentIndex, setter=self.cbFrame.setCurrentIndex, save=True, group="plotMain", changeSignal=self.cbFrame.currentIndexChanged)
        self.addParam("contoursMax", getter=self.dContoursMax.value, setter=self.dContoursMax.setValue, save=True, group="plotMain", changeSignal=self.dContoursMax.valueChanged)
        self.addParam("contoursMin", getter=self.dContoursMin.value, setter=self.dContoursMin.setValue, save=True, group="plotMain", changeSignal=self.dContoursMin.valueChanged)
        self.addParam("contoursLevels", getter=self.iContoursLevels.value, setter=self.iContoursLevels.setValue, save=True, group="plotMain", changeSignal=self.iContoursLevels.valueChanged)

        #TODO: I need to preserve order of params in loadSettings, otherwise I cannot guarantee that e-g- dContoursMin will not be reset in _updateDataBuffer if e.g. interpolation3 is loaded after it

        #visuals, these belong to _ddPlotter
        #~ self._ddPlotter.addParam("contoursColor", #Todo: it works in the old way, lets keep it for now
        #self.pushContoursColor - does not store QColor and cannot be used as data storage, unless you count color of installed icon, which is probably valid approach
        
        import re
        self._colorRE = re.compile("background-color:\s*([#\w]+)")

        self._ddPlotter.addParam("contoursColor", getter=lambda : self._colorRE.findall(self.pushContoursColor.styleSheet())[0], setter=self._setContoursColor, save=True)
        self._ddPlotter.addParam("contoursVisible", getter=self.bContoursVisible.isChecked, setter=self.bContoursVisible.setChecked, save=True, changeSignal=self.bContoursVisible.toggled)

        #~ self.bContoursVisible.setChecked(self._ddPlotter.contoursVisible)
        #~ self.bContoursVisible.toggled.connect(lambda x: self._ddPlotter.setup(contoursVisible=x))
        
        self._currentData = None #interpolated data to be shown on plotMain

        #~ self.cbTRN.currentTextChanged.connect(self._updateDataBuffer)
        #~ self.cbPAIR.currentTextChanged.connect(self._updateDataBuffer)
        #~ self.iInterpolation1.valueChanged.connect(self._updateDataBuffer)
        #~ self.iInterpolation3.valueChanged.connect(self._updateDataBuffer)

        #~ self.cbState.currentIndexChanged.connect(self.updateShow)
        #~ self.cbFrame.currentIndexChanged.connect(self.updatePlotMain)
        #~ self.dContoursMin.valueChanged.connect(self.updatePlotMain)
        #~ self.dContoursMax.valueChanged.connect(self.updatePlotMain)
        #~ self.iContoursLevels.valueChanged.connect(self.updatePlotMain)
        
        if not master:
            #we have to update self.controls.cbState with parent.cbState, also take data from it, but probably not store them separately
            self.cbState.setModel(self.parent().cbState.model())

        if loadFromUI:
            self.cbState.addItem("UI : Data", None)
            self.cbState.addItem("UI : DAS", None)
            self.cbState.addItem("UI : Res", None)
            self.cbState.addItem("UI : OM", None)
        pass
        
    def saveSettings(self):
        res = super().saveSettings() #this is self
        #add contour widgets
        for i in range(2, self.layout().count()):
            w = self.layout().itemAt(i).widget()
            key = "Contours"+str(i-2)
            res[key] = w.saveSettings()
        return res
    
    def loadSettings(self, settings):
        keys = list(settings.keys())
        keys = sorted([it for it in keys if it.startswith("Contours")])
        #remove all contours
        while self.layout().count()>2:
            w = self.layout().itemAt(2).widget()
            self.layout().removeWidget(w)
            w.setParent(None)
        
        for it in keys:
            w = self._addContours()
            w.loadSettings(settings.pop(it))
            
        super().loadSettings(settings)
    
    def _setContoursColor(self, color):
        if color is False or color is True:
            color = QtWidgets.QColorDialog.getColor(QtCore.Qt.black, self, "Select color for contours...")
        else:
            try:
                ccolor=QtGui.QColor(color)
                print("_setContoursColor: color", color, ccolor, ccolor.isValid(), ccolor.red(), ccolor.green(), ccolor.blue())
                color=ccolor
            except:
                import traceback
                print("_setContoursColor, failed, getting color from dialog")
                print(traceback.print_exc())
                color = QtWidgets.QColorDialog.getColor(QtCore.Qt.black, self, "Select color for contours...")
        if not color.isValid: return
        
        self.pushContoursColor.setStyleSheet("QPushButton { background-color: %s}" % color.name())
        #~ px = QtGui.QPixmap(20,20)
        #~ px.fill(color)
        #~ self.pushContoursColor.setIcon(QtGui.QIcon(px))
        self._ddPlotter._afterSet("contoursColor", color)
        #~ self._ddPlotter.setup(contoursColor=color)        
        pass
    
    def updateShow(self, *args):
        #TODO: this whole new approach is weird, it starts to get too ugly
        #   having property-like instance attributes with automatic saving is good
        #   but this whole connection to UI elements is heavy handed and cannot be precisely controled, e.g. here
        #~ self.cbFrame.currentIndexChanged.disconnect(self.updatePlotMain)
        self.cbFrame.clear() #this should not trigger updatePlotMain
        data = self.cbState.currentData()
        if data is not None:
            self.cbFrame.addItems([str(it) for it in data["a2"]])
        #~ self.cbFrame.currentIndexChanged.connect(self.updatePlotMain)
        
        self._updateDataBuffer()
        pass
        
        
    def _updateDataBuffer(self):
        _dataShow = self.cbState.currentData() #TODO: test what is the invalid value (is it None?)
        #~ print("_updateDataBuffer", self.cbState.currentText())
        if _dataShow is None: return

        trn = TRN2trn[self.cbTRN.currentText()]
        pair = self.cbPAIR.currentText()
        i1 = self.iInterpolation1.value()
        i3 = self.iInterpolation3.value()

        #getting the correct dataRange is not easy - we need to reflect PAIR selection
        # - i.e. max is not always from absolute value, but could be just the biggest real part
        # - but that could mean to apply PAIR transfrom on the whole dataset all the time :(
        #TODO: we could hold a copy of data (not interpolated) for every PAIR (more memory, but uninterpolated data are usually not that big)

        data = PAIR[pair](_dataShow[trn])
        a1 = _dataShow["a1"]
        a2 = _dataShow["a2"]
        a3 = _dataShow["a3"]
        
        dataRange = (data.min(), data.max()) #global range

        if i1>0:
            data, a1 = interpolate(data, i1, a1, axis=1)
        if i3>0:
            data, a3 = interpolate(data, i3, a3, axis=2)

        self._currentData = data, a1, a3, a2, dataRange

        #~ self.dContoursMin.valueChanged.disconnect(self.updatePlotMain)
        #~ self.dContoursMax.valueChanged.disconnect(self.updatePlotMain)
        self.dContoursMin.setRange(*dataRange)
        self.dContoursMax.setRange(*dataRange)
        s = (dataRange[1]-dataRange[0])/100.
        self.dContoursMin.setSingleStep(s)
        self.dContoursMax.setSingleStep(s)
        self.dContoursMin.setValue(dataRange[0])
        self.dContoursMax.setValue(dataRange[1])
        #~ self.dContoursMin.valueChanged.connect(self.updatePlotMain)
        #~ self.dContoursMax.valueChanged.connect(self.updatePlotMain)

        self.updatePlotMain()
        pass
        
    def updatePlotMain(self, *args):
        #TODO: main plot should keep zoom/range state of color scale when just changing current (or all the time?)
        #  rescale color scale only on updateBuffer
        print("updatePlotMain", self._currentData is None)
        if self._currentData is None: 
            #TODO: clear plot
            return

        #~ index = self.cbFrame.currentIndex()
        index = self.frame #via params
        print("\t", index)
        if index is -1:
            print('tabCompare.updatePlotMain: invalid current index')
            return

        data, a1, a3, a2, dataRange = self._currentData
        frame = data[index]

        if self.bGlobalColors.isChecked():
            a = max(abs(dataRange[0]), abs(dataRange[1]))
        else:
            a = max(abs(frame.min()), abs(frame.max()))
        dataRange = (-a, a)
        
        print("\tsetData")

        #for NumEnumTransform to work properly it has to be initialized with correct stops
        if self._master:
            self._ddPlotter._xaxis._transform.stops = a1
            self._ddPlotter._yaxis._transform.stops = a3
            #~ self._ddPlotter._xaxis._transform.stops = a3
            #~ self._ddPlotter._yaxis._transform.stops = a1
            #TODO: this should update all contours

        #TODO: bilinear color scale
        #~ self._ddPlotter.colorAxis()._transform.min = frame.min()
        #~ self._ddPlotter.colorAxis()._transform.max = frame.max()
        #~ print("\tmaster", self._master, "color", self._ddPlotter._createColorPlot, "contours", self._ddPlotter.contours)
        #~ print("\tdata", frame.T)
        self._ddPlotter.setup(contours=np.r_[self.dContoursMin.value():self.dContoursMax.value():1j*self.iContoursLevels.value()])
        self._ddPlotter.setData(frame.T, a1, a3, dataRange=(-a, a)) #this could reinitialize NumEnumTransforms, but it has no guarantee that these transforms are actually used
        #~ self._ddPlotter.setData(frame, a3, a1) #this could reinitialize NumEnumTransforms, but it has no guarantee that these transforms are actually used
        #~ print("\tcolor", self._ddPlotter._colorPlot, "contours", self._ddPlotter._contours)
        if self._master:
            self.plotMain.plot().zoomToData()
            #~ self._ddPlotter._xaxis.zoomToData(self._ddPlotter)
            #~ self._ddPlotter._yaxis.zoomToData(self._ddPlotter)
            #~ self._ddPlotter.colorAxis().setRange(-a, a)
        pass
        
    def _addContours(self, *args):
        #we need to create another contour widget and include it into layout
        w = ControlBlob(self.plotMain, False, self, master=False)
        self.layout().addWidget(w)
        print("addContours", w)
        return w
        
    def _removeContours(self, *args):
        print("removeContours", args)
        l = self.layout()
        w = self.sender().parent()
        l.removeWidget(w)
        w.setParent(None)
        
        ...
    pass
