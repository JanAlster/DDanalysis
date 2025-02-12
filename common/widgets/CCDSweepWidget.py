
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

from PyQt5 import QtWidgets, QtCore, QtGui
import time

if __name__=="__main__":
    #create simple app containing just the widget - preparation
    import sys
    import os
    import os.path as op
    sys.path.append(op.join(sys.path[0], "..")) #need to import from directory one above this file
    os.chdir("..")
    print(sys.path)


from common.widgets import PluginPlot, DoubleWidget
from common.widgets.conversionWidgets import StaticSelector
from common.widgets.PluginPlot import QCP

"""
widget to whow CCDSweepSetData

combo box to select sweep from the set (populated from data sweep Ids)

individual sweep is 2D array with one axis CCD width (selectable unit px/nm/(rad/fs)) and the other scan axis
scan axis has one DS moving, the other two stationary
 combobox to select width axis unit
 combobox to select scan axis unit
 labels to show static postions (we do not need multiple axes as only one DS should move)
"""

#plot to show DD data cut (or integral)
# since unit conversion will make axes nonuniform (different step size) and QCP cannot handle that
# keep internal plot (axes) data in base units and only change labels (and tracer)
# but update data to new units (we are showing spectral density which depends on units)
# scan axis can have up to three value ranges (we can use both sides of plot and add one more), show the constant ones too (but not as axis)

WIP

class SweepPlot(PluginPlot.Saving, PluginPlot.Zoom, PluginPlot.Tracer, PluginPlot.DD, PluginPlot.Base):
    def __init__(self, *args, **kw):
        #turn off automatic labels  OR  overload setupTickVectors (auto and then change labels)
        self._axesTransform = {}
        self._ticks = {}
        self._labels = {}

        super().__init__(*args, **kw)
        self._inflateRangeCoef = 0.0 #PluginPlot.Zoom setting, tight fit of data
        
        ar = self.axisRect()
        
        self.xAxis.ticksRequest.connect(lambda : self._transformLabels(self.xAxis))
        self.yAxis.ticksRequest.connect(lambda : self._transformLabels(self.yAxis))
        
        self.xAxis.setAutoTickLabels(False)
        self.yAxis.setAutoTickLabels(False)
        
        self._pointLabel = QCP.QCPItemText(self)
        self._pointLabel.setText("")
        self._pointLabel.setPositionAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignBottom)
        self._pointLabel.setClipToAxisRect(False)
        self._pointLabel.setLayer("legend")
        self._pointLabel.setPadding(QtCore.QMargins(1,1,1,3))
        pos = self._pointLabel.position
        pos.setTypeY(pos.ptViewportRatio)
        pos.setTypeX(pos.ptAxisRectRatio)
        pos.setCoords(1,1)
        pass
        
    #overload DD.setData
    #def setData(self, A, x=None, y=None, xlabel=None, ylabel=None, rescale=True):
    def setData(self, A, scanAxes, w, cut=None, xlabel=None, ylabel=None, zlabel=None, rescale=True, sweepId=None):
        #TODO: labels
        pointLabels = []

        a1, a2, a3 = scanAxes
        try: 
            l1 = len(a1)
            scanAxis = a1
        except:
            l1 = 0
            if xlabel is not None: pointLabels.append((a1, xlabel[0]))
            
        try: 
            l2 = len(a2)
            scanAxis = a2
        except:
            l2 = 0
            if xlabel is not None: pointLabels.append((a2, xlabel[1]))

        try: 
            l3 = len(a3)
            scanAxis = a3
        except:
            l3 = 0
            if xlabel is not None: pointLabels.append((a3, xlabel[2]))

        self._axesTransform[self.yAxis] = scanAxis
        if xlabel is not None: self.yAxis.setLabel(xlabel[0])
        self.yAxis.setVisible(True)

        lw = len(w)
        self._axesTransform[self.xAxis] = w
        if ylabel is not None: self.xAxis.setLabel(ylabel)
        #Todo: naming convention is messed up here
        
        #add sweepId to label
        if sweepId is not None:
            pointLabels.append((sweepId, "sweep"))
        
        if len(pointLabels) > 0:
            #show these points in self._pointLabel
            ttt = []
            for value, label in pointLabels:
                value = str(value)
                if label is None:
                    text = value
                elif " / " in label:
                    text = label.replace(" / ", " "+value+" ") #Todo: only replace first occurence
                else:
                    text = label+" "+value
                ttt.append(text)
            self._pointLabel.setText(" | ".join(ttt))
            pass
        
        #plot data in base units, but intercept tick labels generation and change them
        PluginPlot.DD.setData(self, A, list(range(len(w))), list(range(len(scanAxis))), zlabel=zlabel, rescale=False)
        if rescale: self.rescaleAxes() #might do this directly in DD.setData if we do not do anything else after that
        
        #show point labels
        pass
        
    #overload PluginPlot.Tracer tracer position labels update
    def updateXYPosition(self, key, value):
        #transform key, value
        if self.xAxis in self._axesTransform:
            ikey = int(round(key))
            transform = self._axesTransform[self.xAxis]
            if ikey < 0: ikey = 0
            if ikey >= len(transform): ikey = len(transform)-1
            key = transform[ikey]

        if self.yAxis in self._axesTransform:
            ivalue = int(round(value))
            transform = self._axesTransform[self.yAxis]
            if ivalue < 0: ivalue = 0
            if ivalue >= len(transform): ivalue = len(transform)-1
            value = transform[ivalue]
            
        super().updateXYPosition(key, value)
        pass
    
    def _transformLabels(self, axis):
        #NOTE: apparently I cannot change number of ticks or we face crushes
        #TODO: that leads to multiple labels with the same value at different places (within round of to int) or at the same place (if we set ticks to int values, but keep their number)
        # - we need to tell QCP to use at least 1 step tick increment regardless of range zoom
        
        #assumptions 
        # - we have self._xbase and self._xunits
        # - we have autogenerated tick vector on _xbase space
        # we want to transform it to xunits space (1:1 mapping) 
        # - the transformation is reasonably smooth so we can interpolate if needed (but it should not be necessary if we operate on pixel space)
        # (we can even force ticks to be on pixel space if QCP tries to subdivide pixels)
        #axis = self.sender()
        #print("CutPlot._transformLabelsX: sender", axis, axis in self._axesTransform)
        ticks = axis.tickVector()
        #print("\tticks", ticks) #assuming list of floats
        if axis in self._axesTransform: #saveguard against showing axes before setting actual data
            #we operate on pixel space
            N = len(self._axesTransform[axis])
            nticks = [int(round(it)) for it in ticks]#prevent ticks outside of the data range (not defined in xunits)
            nticks = [it if it>=0 else 0 for it in nticks]
            nticks = [it if it<N else N-1 for it in nticks]
            
            #transform
            
            #TODO: I have a feeling that self._xbase.index(it) == it
            #print("\ttest")
            #print(ticks)
            #print([self._xbase.index(it) for it in ticks])
            ticks = [self._axesTransform[axis][it] for it in nticks] 
            
            #print("\tticks", ticks)
        
        #this is how QCP formats numbers, we want to do the same
        #mTickVectorLabels[i] = mParentPlot->locale().toString(mTickVector.at(i), mNumberFormatChar.toLatin1(), mNumberPrecision);
        toString = self.locale().toString
        nf = axis.numberFormat()[0]
        p = axis.numberPrecision()
        labels = [toString(float(it), nf, p) for it in ticks]
        #print("\tticks", axis.tickVector())
        #print("\tlabels", labels)
        axis.setTickVectorLabels(labels)
        self._labels[axis] = labels #just store to prevent crashes due to garbage collection
        pass
    pass

WIP 

class CCDScanWidget(QtWidgets.QWidget):
    #CutPlot + selection widgets + CCDScanData
    def __init__(self, *args):
        super().__init__(*args)
        
        self._data = None
        
        layout = QtWidgets.QVBoxLayout()
        
        self.plot = CutPlot(self)
        self._rescale = True
        layout.addWidget(self.plot, 1)
        
        cl = QtWidgets.QHBoxLayout()
        self.kind = QtWidgets.QComboBox()
        self.kind.addItems(["spe, spa", "spe, scan", "spa, scan"])
        self.kind.currentTextChanged.connect(self._changeKind)
        
        self.step = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.step.setMinimum(0)
        self.step.setMaximum(0)
        self._stepMemory = dict(zip(["spe, spa", "spe, scan", "spa, scan"], [0,0,0]))
        self.step.valueChanged.connect(self._replot)
            
       
        cl.addWidget(QtWidgets.QLabel("Type"))
        cl.addWidget(self.kind)
        cl.addWidget(QtWidgets.QLabel("Step"))
        cl.addWidget(self.step, 1)

        selectors = ["Spectral", "Height", "Amplitude", "DS1", "DS2", "DS3"] #need to correspond to CCDScanData.convertors
        self._selectors = {}
        for it in selectors:
            sel = StaticSelector()
            sel.currentUnitChanged.connect(self._replot)
            self._selectors[it] = sel
            label = QtWidgets.QLabel(it)
            label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
            cl.addWidget(label)
            cl.addWidget(sel)
        
        layout.addLayout(cl)
        self.setLayout(layout)
        pass

    def setData(self, data):
        #TODO: order of units in convertors is random (that follows from saving calibration files as dict)
        #Todo: unset previous data if needed
        self._lastScanStep = None
        self._data = data
        
        #connect convertors
        for it in self._selectors:
            sel = self._selectors[it]
            unit = sel.currentUnit()
            #disconnect is not needed
            sel.currentUnitChanged.disconnect(self._replot)
            sel.setConvertor(data.convertors[it], unit)
            sel.currentUnitChanged.connect(self._replot)
            pass
            
        self._stepMemory = dict(zip(["spe, spa", "spe, scan", "spa, scan"], [0,0,0]))    
        self._changeKind(self.kind.currentText()) #to set self.step
        pass
    
    def update(self, step):
        print("CCDScanWidget.update", step)
        #if showing cut with scan coordinate, only update internal data buffer of CutPlot by lines from last update step to this update step to minimize time issues
        # - probably will have to keep track of data range too (to keep color scale in check and prevent going over whole dataset each time
        
        #elif showing scanCut with step step, replot
        kind = self.kind.currentText()
        if kind == "spe, spa":
            if self.step.value() == step:
                #if the currently shown is updated, replot
                self._replot()
            elif self.step.value() == self._lastScanStep: #TODO: this branch does not work
                #if it is on the last one available before now, keep to the last one available now
                self.step.setValue(step)
                self._replot()
            return
        #TODO: the rest does not work
        elif kind == "spe, scan":
            cut = self._data.cutSpatial
        else:
            cut = self._data.cutSpectral
        
        scanUnits = [self._selectors[it].currentUnit() for it in ["DS1", "DS2", "DS3"]]
        spatialUnit = self._selectors["Height"].currentUnit()
        spectralUnit = self._selectors["Spectral"].currentUnit()
        
        scanRange = slice(self._lastScanStep, step+1)
        point, data, valueAxis, keyAxis = cut(step, scanUnits, spatialUnit, spectralUnit, scanRange=scanRange)
        
        #now we need to get data into cutPlot
        self.plot.updateData(data, scanRange)
        pass
        
    def _changeKind(self, kind):
        self._rescale = True
        
        #change range of step
        shape = self._data.shape()
        
        if kind == "spe, spa":
            m = shape[0]-1
        elif kind == "spe, scan":
            m = shape[1]-1
        else:
            m = shape[2]-1
        
        self.step.setMaximum(m)
        o = self.step.value()
        n = self._stepMemory[kind]
        if n == o:
            self._replot()
        else:
            self.step.setValue(n) #this will trigger replot
        pass
    
    def _replot(self, *args):
        print("CCDScanWidget._replot start", time.perf_counter())
        #readout units
        scanUnits = [self._selectors[it].currentUnit() for it in ["DS1", "DS2", "DS3"]]
        spatialUnit = self._selectors["Height"].currentUnit()
        spectralUnit = self._selectors["Spectral"].currentUnit()
        
        #select cut 
        kind = self.kind.currentText()
        step = self.step.value()
        self._stepMemory[kind] = step
        
        #labels
        spatialLabel = self._selectors["Height"].currentLabel()+" / "+spatialUnit
        spectralLabel = self._selectors["Spectral"].currentLabel()+" / "+spectralUnit
        scanLabel = [ self._selectors[name].currentLabel()+name[-1]+" / "+unit for name, unit in zip(["DS1", "DS2", "DS3"], scanUnits)]
        
        if kind == "spe, spa":
            cut = self._data.cutScan
            valueLabel = spatialLabel
            keyLabel = spectralLabel
            cutLabel = scanLabel
        elif kind == "spe, scan":
            cut = self._data.cutSpatial
            valueLabel  = scanLabel
            keyLabel = spectralLabel
            cutLabel = spatialLabel
        else:
            cut = self._data.cutSpectral
            valueLabel = scanLabel
            keyLabel = spatialLabel
            cutLabel = spectralLabel
        
        #get data
        point, data, valueAxis, keyAxis = cut(step, scanUnits, spatialUnit, spectralUnit)
        
        #send to self.plot
        #TODO: axis labels + units
        #print("CCDScanWidget._replot", step, data.shape, keyAxis, valueAxis)
        #print("\t", data)
        #TODO: this is terible, we are calling _replot with every step update to scan and it will reset all data to DD point by point, and find minmax and all that :(
        #   - time this and compare to data collection - also terrible - there are at least 3 experiment data updates for one replot (with full wedge scan)
        #   - at least I have to disconnect _replot from updates while it is plotting
        # - it laso depends on what is displayed - problems occur when big cut with scan dimension is shown as we are constantly setting and scanning a large array
        # - we need ability to update data from cut by only one line - and insert directly into QCPColorMapData - i.e. do not use _replot but _update kind of function
        self.plot.setData(data, keyAxis, valueAxis, cut=point, xlabel=keyLabel, ylabel=valueLabel, zlabel=self._selectors["Amplitude"].currentLabel(), cutlabel=cutLabel, rescale=self._rescale)
        self.plot.replot()
        self._rescale = False
        print("CCDScanWidget._replot end", time.perf_counter())
        pass
    pass

if __name__=="__main__":
    #create simple app containing just the dialog
    app = QtWidgets.QApplication(sys.argv)
    window = QtWidgets.QMainWindow()
    w = CutPlot(window)
    window.setCentralWidget(w)
    window.show()
    
    #data
    import numpy as np
    data = np.outer(np.arange(10), np.arange(15))
    x1 = np.arange(15)*2
    x2 = x1[::-1]
    x3 = x2 + 4.1
    y = np.arange(10)+23.1
    y2 = y[::-1]
    y3 = y2+4
    
    w.setData(data, (x1, x2, x3), (y, y2, y3))
    
    app.exec_()
    
