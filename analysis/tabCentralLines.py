
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

# -*- coding: utf-8 -*-
from .tabTemplate import *

"""
calculate central lines of selected peak (peaks) and plot they params (slope, offset) against population time

work of selected data frames


no side plots for cuts
add ROI + fit central lines
side plots for params

show central lines directly in main plot

it will fit the interpolated data (+trn and pair)

"""

from PyQt5 import QtWidgets, QtCore

from common.widgets.GrowingList import GrowingList
from common import delayedUpdate, delayUpdate, releaseUpdate

import os.path as op
class ROIItemWidget(QtWidgets.QWidget):
    SomeValue=QtCore.pyqtSignal()
    NoValue=QtCore.pyqtSignal()
    
    def __init__(self, parent=None):
        super().__init__(parent)
        #self.setupUi(self)
        layout = QtWidgets.QHBoxLayout()
        def _create():
            ttt = QtWidgets.QDoubleSpinBox()
            ttt.setMaximum(99999.00)
            ttt.setVisible(False)
            ttt.setDecimals(0)
            layout.addWidget(ttt)
            return ttt
        self.dW1_1 = _create()
        self.dW1_2 = _create()
        self.dW3_1 = _create()
        self.dW3_2 = _create()
        
        self._button = QtWidgets.QPushButton("+")
        self._button.clicked.connect(self._trigger)
        layout.addWidget(self._button)
        
        layout.setContentsMargins(0,0,0,0)
        self.setLayout(layout)
        
        pass

    def _trigger(self, *args):
        if self.dW1_1.isVisible():
            self.NoValue.emit()
        else:
            self.dW1_1.setVisible(True)
            self.dW1_2.setVisible(True)
            self.dW3_1.setVisible(True)
            self.dW3_2.setVisible(True)
            self._button.setText("x")
            self.SomeValue.emit()
        pass

    def value(self):
        w11 = self.dW1_1.value()
        w12 = self.dW1_2.value()
        w31 = self.dW3_1.value()
        w32 = self.dW3_2.value()
        return min(w11, w12), max(w11, w12), min(w31, w32), max(w31, w32)

    def setValue(self, w11, w12=None, w31=None, w32=None):
        # ~ print("ROIItemWidget.setValue", w11, w12, w31, w32)
        try:
            w11, w12, w31, w32 = w11 #in case lifetime is tuple
        except: pass
        self.dW1_1.setValue(min(w11, w12))
        self.dW1_2.setValue(max(w11, w12))
        self.dW3_1.setValue(min(w31, w32))
        self.dW3_2.setValue(max(w31, w32))
        
        # ~ print("\t", self.dW1_1.isVisible())
        self._trigger()
        # ~ print("\t", self.dW1_1.isVisible())
        pass
    pass

class ROIInputList(QtWidgets.QScrollArea):
    valueChanged=QtCore.pyqtSignal()
    def __init__(self, *args):
        super(ROIInputList, self).__init__(*args)
        widget=GrowingList(ROIItemWidget, ROIItemWidget.value, ROIItemWidget.setValue, self)
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
    

class tabCentralLines(tabTemplateSmall):
    #TODO:can we supress dissection plots from tabTemplate init?
    # ~ outputChanged = QtCore.pyqtSignal(dict)
    def __init__(self, main, *args):
        super().__init__(main, __file__.replace("tabCentralLines.py", 'tabCentralLines.ui'), *args)
        self.pushCalculate.clicked.connect(self._calculate)
        
        self._roilines = []
        
        plot = self.plotSpectral.plot()
        plot.addDDPlotter()
        # ~ plot.setup(x1label="t2 / fs", l1label='ω₁ / cm⁻¹', r1label='ω₃ / cm⁻¹') #TODO: this has to wait for switch to pqr2
        
        plot = self.plotKinetic.plot()
        plot.addAxis(None, pqr.AxisPosition.B, label="t₂ / fs")
        plot.addAxis(None, pqr.AxisPosition.L, label='slope')
        plot.addAxis(None, pqr.AxisPosition.R, label='offset')
        
        
        self.cbKinetic.currentIndexChanged.connect(self.updateROIKinetic)
        self.cbSpectral.currentIndexChanged.connect(self.updateROISpectral)
        
        self.addEditable("ROIInput", self.listROIInput.value, self.listROIInput.setValue)
        pass

    def resetData(self):
        self._centralLines = None 
        super().resetData()

    def updateInput(self, data=None):
        #get data from outside
        #proccess up to the point where input from UI is needed
        #store to self._dataInput
        if data is None: return
        self._dataInput = data
        self.updateShow()

    def updatePlotMain(self, *args, plotMain=True):
        if plotMain: super().updatePlotMain(*args)
        #add ROIS and central lines
        if self._centralLines is None: return
            #note that we do not have centralLines for disabled frames
        
        plot = self.plotMain.plot()
        
        for it in self._roilines: #TODO: possibly clearGraphs will be enough (or put all these to "layer/group" and remove the layer
            plot.removePlotter(it)
            
        self._roilines = []
        self._kineticlines = []
        
        index = self.listFrames.current()
        
        for ROI in self._centralLines:
            w1n, w1x, w3n, w3x = ROI
            ttt = plot.addSeriesPlotter(plot.b1, plot.l1, color="black")
            ttt.setData([w1n, w1x, w1x, w1n, w1n],[w3n, w3n, w3x, w3x, w3n])
            self._roilines.append(ttt)
            
            ttt = self._centralLines[ROI]
            x = ttt["x"]
            y = ttt["y"][index]
            C = ttt["coefs"][index]
            
            ttt = plot.addSeriesPlotter(plot.b1, plot.l1, color="black")
            ttt.setData(x,y)
            self._roilines.append(ttt)

            
            ttt = plot.addSeriesPlotter(plot.b1, plot.l1, color="orange")
            from numpy.polynomial.polynomial import polyval
            ttt.setData(x, polyval(x, C))
            self._roilines.append(ttt)
            
            
        pass
            
    #spectral and kinetic plots are ROI specific
    
    
    @delayedUpdate
    def updateROIKinetic(self, index=-1):
        if index==-1: return
        ROI = self.cbKinetic.itemData(index)
        plot = self.plotKinetic.plot(0)
        plot.clearPlotters() 
        color = pqr.ColorCycle()
        if ROI is None: #all
            for ROI in self._centralLines:
                c = color()
                ttt = self._centralLines[ROI]
                x = ttt["key"]
                C = ttt["coefs"].T
                #TODO: offset should be relative to diagonal (or something...) - offset to zero x is meaningless
                plot.addSeriesPlotter(plot.x1, plot.r1, color=c, name=None, style="--").setData(x, C[0])
                plot.addSeriesPlotter(plot.x1, plot.l1, color=c, name=str(ROI)).setData(x, C[1])
        elif ROI in self._centralLines:
            ttt = self._centralLines[ROI]
            x = ttt["key"]
            C = ttt["coefs"].T
            plot.addSeriesPlotter(plot.x1, plot.r1, color="blue", name="offset", style="--").setData(x, C[0])
            plot.addSeriesPlotter(plot.x1, plot.l1, color="green", name="slope").setData(x, C[1])
        else:
            #this should not happen
            ...
        plot.zoomToData()
        pass
        
    @delayedUpdate
    def updateROISpectral(self, index=-1):
        if index==-1: return
        ROI = self.cbSpectral.itemData(index)
        plot = self.plotSpectral.plot(0)
        p = plot.plotter(0)
        if ROI in self._centralLines:
            data = self._centralLines[ROI]
            p.setData(data["y"].T, data["key"], data["x"])
            plot.zoomToData()
        else:
            p.setData([])#Todo: test this
        pass

    def _calculate(self, *args):
        if self._currentData is None: return
        
        #these are interpolated - i.e. what is shown in main plot
        data, a1, a3, a2, dataRange = self._currentData
        
        self._centralLines = {}
        
        ROIs = self.listROIInput.value()
        for ROI in ROIs:
            w1n, w1x, w3n, w3x = ROI
            
            i1n, i1x = a1.searchsorted((w1n, w1x))
            i3n, i3x = a3.searchsorted((w3n, w3x))
            
            _data = data[:,i1n:i1x, i3n:i3x]
            _a1 = a1[i1n:i1x]
            _a3 = a3[i3n:i3x]
            
            #calculate central lines
            self._centralLines[ROI] = {"x":_a1, "y":_a3[_data.argmax(axis=2)], "key": a2}
            
            #fit slope
            from numpy.polynomial.polynomial import polyfit
            self._centralLines[ROI]["coefs"] = polyfit(_a1, self._centralLines[ROI]["y"].T, 1).T
        
        oldKinetic = self.cbKinetic.currentData()
        oldSpectral = self.cbSpectral.currentData()
        
        delayUpdate() #this should prevent several calls of updateROI... when cb are cleared
        self.cbKinetic.clear()
        self.cbSpectral.clear()
        self.cbKinetic.addItem("all", None)
        for ROI in ROIs:
            self.cbKinetic.addItem(str(ROI), ROI)
            self.cbSpectral.addItem(str(ROI), ROI)
        
        index = self.cbKinetic.findData(oldKinetic)
        self.cbKinetic.setCurrentIndex(index if index>=0 else 0)

        index = self.cbSpectral.findData(oldSpectral)
        self.cbSpectral.setCurrentIndex(index if index>=0 else 0)
        
        releaseUpdate()
        
        self.updatePlotMain(plotMain=False)
        pass
        
    pass

