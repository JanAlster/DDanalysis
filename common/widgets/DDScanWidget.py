
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

# TODO: move the app to utils
if __name__=="__main__":
    #create simple app containing just the widget - preparation
    import sys
    import os
    import os.path as op
    sys.path.append(op.join(sys.path[0], "..")) #need to import from directory one above this file
    os.chdir("..")
    print(sys.path)


from common.widgets import DoubleWidget
from common.widgets.conversionWidgets import StaticSelector
from pqr import pqr2 as pqr

"""
we need a widget to show this data set 
 
problem with transformation of data to different units is that QCP cannot handle non-uniform axes, esp. with 2D color maps (like stretching some cells of the colormap)
 - we can either interpolate data to uniform step - which I do not want to do
 - or keep the plot in uniform axes - i.e. CCD pixels (base units of convertors) and only convert what is shown on axes labels (and tracer)
 - of course handling of the scan axis has to be special too
"""

#plot to show DD data cut (or integral)
# since unit conversion will make axes nonuniform (different step size) and QCP cannot handle that
# keep internal plot (axes) data in base units and only change labels (and tracer)
# but update data to new units (we are showing spectral density which depends on units)

#TODO: for DDScanData only spectral axis can change, other two are fixed (and coherenceAxis is equispaced)

class CutPlot2(pqr.FigureWidget):
    # PQR  - this is a mess :)
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self._inflateRangeCoef = 0.0 #PluginPlot.Zoom setting, tight fit of data
        #turn off automatic labels  OR  overload setupTickVectors (auto and then change labels)
        self._axesTransform = {}
        self._ticks = {}
        self._labels = {}

        plot = self.plot()
        plot.addDDPlotter()

        # PQR
        # self.xAxis.ticksRequest.connect(lambda : self._transformLabels(self.xAxis))
        # self.yAxis.ticksRequest.connect(lambda : self._transformLabels(self.yAxis))
        
        # self.xAxis.setAutoTickLabels(False)
        # self.yAxis.setAutoTickLabels(False)
        
        # self._pointLabel = QCP.QCPItemText(self)
        # self._pointLabel.setText("Ahoj")
        # self._pointLabel.setPositionAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignBottom)
        # self._pointLabel.setClipToAxisRect(False)
        # self._pointLabel.setLayer("legend")
        # self._pointLabel.setPadding(QtCore.QMargins(1,1,1,3))
        # pos = self._pointLabel.position
        # pos.setTypeY(pos.ptViewportRatio)
        # pos.setTypeX(pos.ptAxisRectRatio)
        # pos.setCoords(1,1)
        
        
        pass
        
    #we need DD cut data, original axes, transformed axes, labels
    #overload DD.setData
    #def setData(self, A, x=None, y=None, xlabel=None, ylabel=None, rescale=True):
    def setData(self, A, x, y, cut=None, xlabel=None, ylabel=None, cutlabel=None, zlabel=None, rescale=True):
        plot = self.plot()
        #TODO: labels
        pointLabels = []

        #x is normal axis
        lx = len(x)
        self._axesTransform[plot.x1] = x
        if xlabel is not None: plot.x1.label = xlabel

        #y is normal axis
        ly = len(y)
        self._axesTransform[plot.y1] = y
        if ylabel is not None: plot.y1.label = ylabel
        
        #cut point
        #could be scan point or single point
        if cut is not None:
            #single point
            #print("CutPlot2.setData: cut is normal point")
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
            # PQR self._pointLabel.setText(" | ".join(ttt))  - possibly send to hover
            pass
        
        #print("CutPlot2.setData: data size", A.shape)
        #TODO: possibly flip data (and axes) if axes are descending
        
        #plot data in base units, but intercept tick labels generation and change them
        plot.plotter(0).setData(A) #PQR , zlabel=zlabel, rescale=False)
        # PQR if rescale: self.rescaleAxes() #might do this directly in DD.setData if we do not do anything else after that
        plot.zoomToData()

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
        self._labels[axis] = labels #just store to prevent crashs due to garbage collection
        pass
    pass


class DDScanWidget(QtWidgets.QWidget):
    #CutPlot2 + selection widgets + DDScanData
    def __init__(self, *args):
        super().__init__(*args)
        
        self._data = None
        
        layout = QtWidgets.QVBoxLayout()
        
        self.plot = CutPlot2(self)
        self._rescale = True
        layout.addWidget(self.plot, 1)
        
        cl = QtWidgets.QHBoxLayout()
        self.kind = QtWidgets.QComboBox()
        self.kind.addItems(["spe, coh", "spe, pop", "coh, pop"])
        self.kind.currentTextChanged.connect(self._changeKind)
        
        self.step = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.step.setMinimum(0)
        self.step.setMaximum(0)
        self._stepMemory = dict(zip(["spe, coh", "spe, pop", "coh, pop"], [0,0,0]))
        self.step.valueChanged.connect(self._replot)
            
       
        cl.addWidget(QtWidgets.QLabel("Type"))
        cl.addWidget(self.kind)
        cl.addWidget(QtWidgets.QLabel("Step"))
        cl.addWidget(self.step, 1)

        selectors = ["Spectral", "Amplitude"] #need to correspond to DDScanData.convertors
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
            
        self._stepMemory = dict(zip(["spe, coh", "spe, pop", "coh, pop"], [0,0,0]))    
        self._changeKind(self.kind.currentText()) #to set self.step
        pass
        
    def _changeKind(self, kind):
        self._rescale = True
        
        #change range of step
        shape = self._data.shape()
        
        if kind == "spe, coh":
            m = shape[0]-1
        elif kind == "spe, pop":
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
        #readout units
        spectralUnit = self._selectors["Spectral"].currentUnit()
        
        #select cut 
        kind = self.kind.currentText()
        step = self.step.value()
        self._stepMemory[kind] = step
        #TODO: memory step does not work properly
        
        #labels
        spectralLabel = self._selectors["Spectral"].currentLabel()+" / "+spectralUnit
        coherenceLabel = "Coherence time / fs"
        populationLabel = "Population time / fs"
        
        if kind == "spe, coh":
            cut = self._data.cutPopulation
            valueLabel = coherenceLabel
            keyLabel = spectralLabel
            cutLabel = populationLabel
        elif kind == "spe, pop":
            cut = self._data.cutCoherence
            valueLabel  = populationLabel
            keyLabel = spectralLabel
            cutLabel = coherenceLabel
        else:
            cut = self._data.cutSpectral
            valueLabel = populationLabel
            keyLabel = coherenceLabel
            cutLabel = spectralLabel
        
        #get data
        point, data, valueAxis, keyAxis = cut(step, spectralUnit)
        
        #send to self.plot
        #TODO: axis labels + units
        #print("DDScanWidget._replot", step, data.shape, keyAxis, valueAxis)
        #print("\t", data)
        self.plot.setData(data, keyAxis, valueAxis, cut=point, xlabel=keyLabel, ylabel=valueLabel, zlabel=self._selectors["Amplitude"].currentLabel(), cutlabel=cutLabel, rescale=self._rescale)
        self.plot.replot()
        self._rescale = False
        pass
    pass

if __name__=="__main__":
    #create simple app containing just the dialog
    app = QtWidgets.QApplication(sys.argv)
    window = QtWidgets.QMainWindow()
    w = CutPlot2(window)
    window.setCentralWidget(w)
    window.show()
    
    #data
    import numpy as np
    data = np.outer(np.arange(10), np.arange(15))
    x1 = np.arange(15)*2
    y = np.arange(10)+23.1
    
    w.setData(data, x1, y)
    
    app.exec_()
    
