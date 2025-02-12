
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
we need a widget to show this data set 
 - possibly multiple axes if showing scan coordinate
 
we also need a lower dimensionality containers for e.g. data integrated over CCD height - this is equivalent to a cut, only the widget does not need so much selection tools
problem with transformation of data to different units is that QCP cannot handle non-uniform axes, esp. with 2D color maps (like stretching some cells of the colormap)
 - we can either interpolate data to uniform step - which I do not want to do
 - or keep the plot in uniform axes - i.e. CCD pixels (base units of convertors) and only convert what is shown on axes labels (and tracer)
 - of course handling of the scan axis has to be special too
"""

#TODO: IMPORTANT this is seriously not working, can lead to crashes easily

#plot to show DD data cut (or integral)
# since unit conversion will make axes nonuniform (different step size) and QCP cannot handle that
# keep internal plot (axes) data in base units and only change labels (and tracer)
# but update data to new units (we are showing spectral density which depends on units)
# scan axis can have up to three value ranges (we can use both sides of plot and add one more), show the constant ones too (but not as axis)

class CutPlot(PluginPlot.Saving, PluginPlot.Zoom, PluginPlot.Tracer, PluginPlot.DD, PluginPlot.Base):
    def __init__(self, *args, **kw):
        #turn off automatic labels  OR  overload setupTickVectors (auto and then change labels)
        self._axesTransform = {}
        self._ticks = {}
        self._labels = {}

        super().__init__(*args, **kw)
        self._inflateRangeCoef = 0.0 #PluginPlot.Zoom setting, tight fit of data
        
        ar = self.axisRect()
        self.yAxis3 = ar.addAxis(self.xAxis.atRight)
        self.xAxis3 = ar.addAxis(self.xAxis.atTop)
        
        self.xAxis3.setVisible(False)
        self.yAxis3.setVisible(False)
        
        self.xAxis2.setTickLabels(True)
        self.yAxis2.setTickLabels(True)
        
        self.xAxis.ticksRequest.connect(lambda : self._transformLabels(self.xAxis))
        self.xAxis2.ticksRequest.connect(lambda : self._transformLabels(self.xAxis2))
        self.xAxis3.ticksRequest.connect(lambda : self._transformLabels(self.xAxis3))
        self.yAxis.ticksRequest.connect(lambda : self._transformLabels(self.yAxis))
        self.yAxis2.ticksRequest.connect(lambda : self._transformLabels(self.yAxis2))
        self.yAxis3.ticksRequest.connect(lambda : self._transformLabels(self.yAxis3))
        
        self.xAxis.setAutoTickLabels(False)
        self.xAxis2.setAutoTickLabels(False)
        self.xAxis3.setAutoTickLabels(False)
        self.yAxis.setAutoTickLabels(False)
        self.yAxis2.setAutoTickLabels(False)
        self.yAxis3.setAutoTickLabels(False)
        
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
        
    #we need DD cut data, original axes, transformed axes, labels
    #overload DD.setData
    #def setData(self, A, x=None, y=None, xlabel=None, ylabel=None, rescale=True):
    def setData(self, A, x, y, cut=None, xlabel=None, ylabel=None, cutlabel=None, zlabel=None, rescale=True):
        #TODO: labels
        pointLabels = []
        #any of (but only one) x and y can be triplet of (scanDS1, scanDS2, scanDS3)
        #this can be difficult to recognize automatically
        #length of x is 3 and at least one of elements is iterable
        try:
            assert(len(x)==3)
            l1 = 0
            try: 
                l1 = len(x[0])
            except:
                pass
                
            l2 = 0
            try: 
                l2 = len(x[1])
            except:
                pass

            l3 = 0
            try: 
                l3 = len(x[2])
            except:
                pass
            
            #x is scan axis
            
            lx = max([l1,l2,l3])
            assert(lx>0)
            axes = [self.xAxis, self.xAxis2, self.xAxis3]
            #print("CutPlot.setData: x is scan axis, length", lx)
            
            if l1>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = x[0]
                if xlabel is not None: axis.setLabel(xlabel[0])
                axis.setVisible(True)
            else:
                #add x[0] scalar value to displaying element
                if xlabel is not None: pointLabels.append((x[0], xlabel[0]))
                pass
                
            if l2>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = x[1]
                if xlabel is not None: axis.setLabel(xlabel[1])
                axis.setVisible(True)
            else:
                #add x[1] scalar value to displaying element
                if xlabel is not None: pointLabels.append((x[1], xlabel[1]))
                pass
            
            if l3>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = x[2]
                if xlabel is not None: axis.setLabel(xlabel[2])
                axis.setVisible(True)
            else:
                #add x[2] scalar value to displaying element
                if xlabel is not None: pointLabels.append((x[2], xlabel[2]))
                pass
                
            #hide remaining axes
            for axis in axes:
                axis.setVisible(False)
                
            pass
        except:
            #x is normal axis
            lx = len(x)
            #print("CutPlot.setData: x is normal axis, length", lx)
            self._axesTransform[self.xAxis] = x
            if xlabel is not None: self.xAxis.setLabel(xlabel)
            #hide xAxis2, xAxis3
            self.xAxis2.setVisible(False)
            self.xAxis3.setVisible(False)
            pass

        try:
            assert(len(y)==3)
            l1 = 0
            try: 
                l1 = len(y[0])
            except:
                pass
                
            l2 = 0
            try: 
                l2 = len(y[1])
            except:
                pass

            l3 = 0
            try: 
                l3 = len(y[2])
            except:
                pass
            
            #y is scan axis
            
            ly = max([l1,l2,l3])
            assert(ly>0)
            axes = [self.yAxis, self.yAxis2, self.yAxis3]
            
            #print("CutPlot.setData: y is scan axis, length", ly)
            
            if l1>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = y[0]
                if ylabel is not None: axis.setLabel(ylabel[0])
                axis.setVisible(True)
            else:
                #add y[0] scalar value to displaying element
                if ylabel is not None: pointLabels.append((y[0], ylabel[0]))
                pass
                
            if l2>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = y[1]
                if ylabel is not None: axis.setLabel(ylabel[1])
                axis.setVisible(True)
            else:
                #add y[1] scalar value to displaying element
                if ylabel is not None: pointLabels.append((y[1], ylabel[1]))
                pass
            
            if l3>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = y[2]
                if ylabel is not None: axis.setLabel(ylabel[2])
                axis.setVisible(True)
            else:
                #add y[2] scalar value to displaying element
                if ylabel is not None: pointLabels.append((y[2], ylabel[2]))
                pass
                
            #hide remaining axes
            for axis in axes:
                axis.setVisible(False)
                
            pass
        except:
            #y is normal axis
            ly = len(y)
            #print("CutPlot.setData: y is normal axis, length", ly)
            self._axesTransform[self.yAxis] = y
            if ylabel is not None: self.yAxis.setLabel(ylabel)
            #hide yAxis2, yAxis3
            self.yAxis2.setVisible(False)
            self.yAxis3.setVisible(False)
            pass
        
        #cut point
        #could be scan point or single point
        if cut is not None:
            try:
                assert(len(cut)==3)
                #scan point
                #print("CutPlot.setData: cut is scan point")
                if cutlabel is None: cutlabel = ("","","")
                for it in zip(cut, cutlabel):
                    pointLabels.append(it)
            except:
                #single point
                #print("CutPlot.setData: cut is normal point")
                if cutlabel is None: cutlabel = ""
                pointLabels.append((cut, cutlabel))
                pass
        
        if len(pointLabels) > 0:
            #TODO: show these points in self._pointLabel
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
        
        #print("CutPlot.setData: data size", A.shape)
        #TODO: possibly flip data (and axes) if axes are descending
        
        #plot data in base units, but intercept tick labels generation and change them
        PluginPlot.DD.setData(self, A, list(range(lx)), list(range(ly)), zlabel=zlabel, rescale=False)
        if rescale: self.rescaleAxes() #might do this directly in DD.setData if we do not do anything else after that
        
        #show point labels
        
        pass
        
    def updateData(self, data, scanRange):
        print("CutPlot.updateData", data.shape, scanRange, scanRange.start)
        #insert data directly into DD colormap data
        #TODO: this will not update DD._cmorder (will not update plot in contrast colors)
        #also it is not elegant to access directly _cmdata
        #also DD show actually a *copy* of _cdmdata or _cmorder in DD.cm.data
        #to be fast we can put data directly into DD.cm.data, but it will not work if _cmorder is there (wrong scaling)
        # and it will be only temporary as any change (e.g. switching to/from contrast colors) will rest it 
        # to store it in all DD.cm.data and DD._cmdata and DD._cmorder is too much trouble from here
        # on the other hand, I dont see a need for DD to allow rewriting part of the data
        # the best way would be to make DD.setData fast, but that is probably impossible
        cmdata = self.cm.data()
        for j in range(data.shape[1]):
            for i in range(data.shape[0]):
                cmdata.setCell(j,i+scanRange.start, data[i,j])
                pass
            pass            
        self.cm.rescaleDataRange() #this should update min max from cmdata internal values (not that I believe that much)
        self.replot()
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
        self._labels[axis] = labels #just store to prevent crashs due to garbage collection
        pass
    pass


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
        if self._data is None: return
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
    
