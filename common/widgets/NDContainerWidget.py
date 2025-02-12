
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
    sys.path.append(op.join(op.dirname(__file__), "../..")) #need to import from directory one above this file
    os.chdir("..")
    print(sys.path)


#from common.widgets import PluginPlot, DoubleWidget
from common.widgets.conversionWidgets import StaticSelector
#from common.widgets.PluginPlot import QCP
from pqr import pqr

import numpy as np

"""
NDContainer is a N dimensional container

we need the ability to display such data

for now we will assume that N>1

we can only show 2D so we will select 
 horizontal axis
 vertical axis
 (and units for both)
and the rest of the axes needs to select one point
and plot will display such a cut through the data

ideally the other axes would select specific points by sliders for easy adjustment
(also units)

axis is defined by NDContainer.NDContainerAxis interface
(axisID, convertor (with units), label (for plots))

so we need control widget for each axis
(axisID label, unit selector, position slider, position value label)
   posibly check box indicating data should be integrated over this axis

control widgets for axes selected for plot should be hidden

plot also needs control widget
(V: axisID selector, units selector, <->, H:, axisID selector, units selector)
 plot will also use axis label
   (note that e.g. ScanAxis might need to put static subaxes positions into the axis label or somewhere else)

 units selector in plot control widget should be tied with selector in axis control widget for selected axis
  (or at least it should be updated after axis is deselected from the plot)

"""

"""
whole plotting is ugly
QCP will not be fast for DD plots so I can probalby give up on separating setData and updateData
(anyway has to redraw everything and setting axes will not be that slow)

at this point it would be better to keep it simple
 - it makes sense to separate axes creation and data setting for NDContainerWidget (creation of axesControl)
 - but not for CutPlot
"""

#plot to show DD data cut
# since unit conversion will make axes nonuniform (different step size) and QCP cannot handle that
# keep internal plot (axes) data in base units and only change labels (and tracer)
# but update data to new units (we are showing spectral density which depends on units)

#TODO: point labels in Axescontrolwidgets change size too much which resizes the slider, this is not nice when you try to drag the slider, we should fix the size

#tracer does not keep position in continuous experiment (then it does after zoom in zoom out) ???

class CutPlot(pqr.FigureWidget):
    #this can accommodate three sub axes for one axis
    def __init__(self, *args, **kw):
        #turn off automatic labels  OR  overload setupTickVectors (auto and then change labels)
        self._axesTransform = {}
        self._ticks = {}
        self._labels = {}

        super().__init__(*args, **kw)
        
        p = self.plot()

        #potentially one of the displayed data cut axes can be a triple axis (for that we have three axes) - position of delay stages
        self.yAxis = p.addAxis(None, pqr.AxisPosition.L)
        self.xAxis = p.addAxis(None, pqr.AxisPosition.B)
        self.yAxis2 = p.addAxis(None, pqr.AxisPosition.R)
        self.xAxis2 = p.addAxis(None, pqr.AxisPosition.T)
        self.yAxis3 = p.addAxis(None, pqr.AxisPosition.R)
        self.xAxis3 = p.addAxis(None, pqr.AxisPosition.T)
        
        self.xAxis3.setVisible(False)
        self.yAxis3.setVisible(False)
        
        #TODO: right now axes are ignored, we need to connect old QCP approach to pqr axes tranaforms (difficult thing is that we might need to swap axes)
#        self.xAxis.ticksRequest.connect(lambda : self._transformLabels(self.xAxis))
#        self.xAxis2.ticksRequest.connect(lambda : self._transformLabels(self.xAxis2))
#        self.xAxis3.ticksRequest.connect(lambda : self._transformLabels(self.xAxis3))
#        self.yAxis.ticksRequest.connect(lambda : self._transformLabels(self.yAxis))
#        self.yAxis2.ticksRequest.connect(lambda : self._transformLabels(self.yAxis2))
#        self.yAxis3.ticksRequest.connect(lambda : self._transformLabels(self.yAxis3))
        
        p.addDDPlotter(p.x1, p.y1)
        
        self._pointLabel = pqr.LabelGraphicsLayoutItem("")
#        self._pointLabel.setPositionAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignBottom)
#        self._pointLabel.setClipToAxisRect(False)
#        self._pointLabel.setLayer("legend")
#        self._pointLabel.setPadding(QtCore.QMargins(1,1,1,3))
#        pos = self._pointLabel.position
#        pos.setTypeY(pos.ptViewportRatio)
#        pos.setTypeX(pos.ptAxisRectRatio)
#        pos.setCoords(1,1)
        pass
        
    #we need DD cut data, original axes, transformed axes, labels
    #overload DD.setData
    #def setData(self, A, x=None, y=None, xlabel=None, ylabel=None, rescale=True):
    def setData(self, A, x, y, xlabel=None, ylabel=None, zlabel=None, rescale=True):
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
                if xlabel is not None: axis.label = xlabel[0]
                axis.setVisible(True)
            else:
                #add x[0] scalar value to displaying element
                if xlabel is not None: pointLabels.append((x[0], xlabel[0]))
                pass
                
            if l2>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = x[1]
                if xlabel is not None: axis.label = xlabel[1]
                axis.setVisible(True)
            else:
                #add x[1] scalar value to displaying element
                if xlabel is not None: pointLabels.append((x[1], xlabel[1]))
                pass
            
            if l3>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = x[2]
                if xlabel is not None: axis.label = xlabel[2]
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
            if xlabel is not None: self.xAxis.label = xlabel
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
                if ylabel is not None: axis.label = ylabel[0]
                axis.setVisible(True)
            else:
                #add y[0] scalar value to displaying element
                if ylabel is not None: pointLabels.append((y[0], ylabel[0]))
                pass
                
            if l2>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = y[1]
                if ylabel is not None: axis.label = ylabel[1]
                axis.setVisible(True)
            else:
                #add y[1] scalar value to displaying element
                if ylabel is not None: pointLabels.append((y[1], ylabel[1]))
                pass
            
            if l3>0:
                axis = axes.pop(0)
                self._axesTransform[axis] = y[2]
                if ylabel is not None: axis.label = ylabel[2]
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
            if ylabel is not None: self.yAxis.label = ylabel
            #hide yAxis2, yAxis3
            self.yAxis2.setVisible(False)
            self.yAxis3.setVisible(False)
            pass
        
        if len(pointLabels) > 0: #TODO: no need to show axis points in point label (NDContainerWidgets shows that), but we might show value under tracer
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
        
        #TODO: possibly flip data (and axes) if axes are descending
        
        #plot data in base units, but intercept tick labels generation and change them
        self.plot().plotter(0).setData(A) 
#        PluginPlot.DD.setData(self, A, list(range(lx)), list(range(ly)), zlabel=zlabel, rescale=rescale)
        pass

#ignore tracer for now    
#    #overload PluginPlot.Tracer tracer position labels update
#    def updateXYPosition(self, key, value):
#        #transform key, value
#        if self.xAxis in self._axesTransform:
#            ikey = int(round(key))
#            transform = self._axesTransform[self.xAxis]
#            if ikey < 0: ikey = 0
#            if ikey >= len(transform): ikey = len(transform)-1
#            key = transform[ikey]
#
#        if self.yAxis in self._axesTransform:
#            ivalue = int(round(value))
#            transform = self._axesTransform[self.yAxis]
#            if ivalue < 0: ivalue = 0
#            if ivalue >= len(transform): ivalue = len(transform)-1
#            value = transform[ivalue]
#            
#        super().updateXYPosition(key, value)
#        pass
    
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
        if axis in self._axesTransform: #saveguard against showing axes before setting actual data
            #we operate on pixel space
            N = len(self._axesTransform[axis])
            nticks = [int(round(it)) for it in ticks]#prevent ticks outside of the data range (not defined in xunits)
            nticks = [it if it>=0 else 0 for it in nticks]
            nticks = [it if it<N else N-1 for it in nticks]
            
            #transform
            
            #TODO: I have a feeling that self._xbase.index(it) == it
            ticks = [self._axesTransform[axis][it] for it in nticks] 
            
        
        #this is how QCP formats numbers, we want to do the same
        #mTickVectorLabels[i] = mParentPlot->locale().toString(mTickVector.at(i), mNumberFormatChar.toLatin1(), mNumberPrecision);
        toString = self.locale().toString
        nf = axis.numberFormat()[0]
        p = axis.numberPrecision()
        labels = [toString(float(it), nf, p) for it in ticks]
        axis.setTickVectorLabels(labels)
        self._labels[axis] = labels #just store to prevent crashs due to garbage collection
        pass
    pass


class NDContainerAxisControlWidget(QtWidgets.QFrame):
    currentUnitChanged = QtCore.pyqtSignal(str)
    currentPositionChanged = QtCore.pyqtSignal(int)
    def __init__(self, axis, *args):
        super().__init__(*args)
        self._axis = axis
        
        layout = QtWidgets.QHBoxLayout()
        
        #axisID label
        layout.addWidget(QtWidgets.QLabel(str(axis.axisID)))
        self.axisID = axis.axisID
        
        #unit selector that can work as selector for axis convertor
        self.selector = StaticSelector()
        if axis.convertor is not None:
            self.selector.setConvertor(axis.convertor)
        self.selector.currentUnitChanged.connect(self.currentUnitChanged)
        layout.addWidget(self.selector)
        
        #position slider
        self.step = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.step.setMinimum(0)
        self.step.setMaximum(axis.length()-1)
        self.step.valueChanged.connect(self.currentPositionChanged)
        self.step.valueChanged.connect(self._updateCurrentPositionLabel)
        layout.addWidget(self.step, 1)
        
        #position label (in units)
        self.lPosition = QtWidgets.QLabel("")
        layout.addWidget(self.lPosition)
        
        layout.setContentsMargins(0,0,0,0)
        self.setLayout(layout)
        
        #self.setFrameStyle(self.Panel)
        self._updateCurrentPositionLabel(self.position())
        pass
        
    def position(self):
        return self.step.value()
        
    def unit(self):
        return self.selector.currentUnit()
        
    def units(self):
        return [self.selector.itemText(i) for i in range(self.selector.count())]
        
    def _updateCurrentPositionLabel(self, step):
        self.lPosition.setText(self._axis.formatPosition(step, self.unit()))
        pass
    pass


class ComboBox(QtWidgets.QComboBox):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self._prevText = None
        self._thisText = None
        self.currentTextChanged.connect(self._slot)
        pass
    
    def _slot(self, text):
        self._prevText = self._thisText
        self._thisText = text
        pass
        
    def clear(self, *args):
        super().clear(*args)
        self._prevText = None
        self._thisText = None
        pass
    pass

#TODO:
# recent change to ScanAxis inner works make container widget display even non-moving axes (instead of just printing the static point)

#TODO:
#interplay of unit selection between HV and axesControl
#interplay of unit selection between self and master unit selectors (on resources)
#   not sure how to approach that as master unit selector are not guaranteed to have the same units as convertors in data
class NDContainerWidget(QtWidgets.QWidget):
    cutChanged = QtCore.pyqtSignal(tuple) #axes or data changed
    tracesChanged = QtCore.pyqtSignal(tuple) #all is the same, only tracer moved -> this should give V and H cut through cut on position of tracer (maybe even trace through the other dimensions at tracer position? i.e. population evolution for DD)
    """
    setData - set NDContainer to show in the widget, will adjust self controls to accommodate all axes
    updateData - replot data (i.e. new data available from experiment), must not change axes, data are stored in the original container set by setData
    """
    def __init__(self, *args):
        super().__init__(*args)
        
        self._data = None
        self._cutData = None #actual exported data shown in self.plot (extra copy for easy reporoting tracesChanged)
        
        layout = QtWidgets.QVBoxLayout()
        
        self.plot = CutPlot(self)
#        self.plot.tracerMoveIndex.connect(self._tracerMoveIndex)
        
        layout.addWidget(self.plot, 1)
        
        cl = QtWidgets.QHBoxLayout()
        cl.addWidget(QtWidgets.QLabel("V:"))
        self.cbVertical = ComboBox()
        self.cbVertical.currentTextChanged.connect(self._changeVertical) #this will trigger change of units, which will trigger axesChanged signal
        cl.addWidget(self.cbVertical)
        self.cbVerticalUnits = QtWidgets.QComboBox()
        self.cbVerticalUnits.currentTextChanged.connect(self.updateData) #_changeVerticalUnits)
        cl.addWidget(self.cbVerticalUnits)
        #flip is not needed, it is enough to select the other axis in the first selector (if it is the same as the axis in the other selector, axes will be flipped)
        #self.pushFlip = QtWidgets.QPushButton("<->")
        #self.pushFlip.clicked.connect(self._flipVH)
        #cl.addWidget(self.pushFlip)
        cl.addWidget(QtWidgets.QLabel("H:"))
        self.cbHorizontal = ComboBox()
        self.cbHorizontal.currentTextChanged.connect(self._changeHorizontal)
        cl.addWidget(self.cbHorizontal)
        self.cbHorizontalUnits = QtWidgets.QComboBox()
        self.cbHorizontalUnits.currentTextChanged.connect(self.updateData) #_changeHorizontalUnits)
        cl.addWidget(self.cbHorizontalUnits)
        cl.addStretch(1)

        #axis control widgets can only be added after we have data
        #   the same for filling cbVertical and Horizontal
        self._axesControls = {}
        
        self._replotLock = None


        layout.addLayout(cl)
        #layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(1)
        self.setLayout(layout)
        pass

    #setData will store the data buffer, update/create axesControls, etc. (optionally plot too)
    #    (data buffer might have all axes ready, but not data themselves)
    #updateData will plot relevant cut into self.plot
    def setData(self, data, plot=False, rescale=None):
        #unset previous data if needed
        self._replotLock = "data"
        prevH = None
        prevV = None
        prevVU = None
        prevHU = None
        prevData = self._data
        if self._data is not None:
            #(maybe?) we need to reset cbVertical, cbVerticalUnits, cbHorizontal, cbHorizontalUnits
            #we need to remove axis control widgets connected to previous data
            layout = self.layout()
            for axis in self._axesControls:
                A = self._axesControls[axis]
                layout.removeWidget(A)
                A.currentPositionChanged.disconnect(self.updateData)
                A.currentUnitChanged.disconnect(self.updateData)
                self._axesControls[axis].setParent(None)
                pass
            prevV = self.cbVertical.currentText()
            prevH = self.cbHorizontal.currentText()
            prevVU = self.cbVerticalUnits.currentText()
            prevHU = self.cbHorizontalUnits.currentText()
            self.cbVertical.clear() #this should not trigger signal, or signal handler must handle empty selection
            self.cbHorizontal.clear() #this should not trigger signal, or signal handler must handle empty selection
            self._axesControls = {}
            pass
        
        #set new data
        self._data = data #data is NDContainer (or compatible)
        
        axes = data.axes() #axes asixIDs in correct order, suitable for setting cbVertical, cbHorizontal and axis control widgets
        assert(len(axes)>=2)
        
        for axis in axes:
            control = NDContainerAxisControlWidget(data.axis(axis))
            control.currentPositionChanged.connect(self.updateData)
            control.currentUnitChanged.connect(self.updateData)
            self.layout().addWidget(control)
            self._axesControls[axis] = control

        #fill in cbVertical na cbHorizontal
        #that should also trigger filling unit selectors
        
        self.cbVertical.addItems(axes)
        self.cbHorizontal.addItems(axes)

        #we should set self._cutAxesChanged (mark if shown axes changed with new setData, used to auto-rescale only with new data size) 
        #to False if (cut axes do not change)
        # we have prevV and prevH
        # new data have them too
        # old and new axes are the same
        #to True otherwise
        self._cutAxesChanged = True

        if prevV is None or prevV not in axes or prevH not in axes:
            #select first two axes
            self.cbVertical.setCurrentText(axes[0]) #probably already is from addItems
            self.cbHorizontal.setCurrentText(axes[1])
        else:
            self.cbVertical.setCurrentText(prevV)
            self.cbHorizontal.setCurrentText(prevH)
            self.cbVerticalUnits.setCurrentText(prevVU)
            self.cbHorizontalUnits.setCurrentText(prevHU)
            #note that prevData cannot be None here
            try:
                if np.allclose(prevData.axis(prevV).values(prevVU), data.axis(prevV).values(prevVU)) and np.allclose(prevData.axis(prevH).values(prevHU), data.axis(prevH).values(prevHU)):
                    self._cutAxesChanged = False
            except:
                #in case of an error, like axes cannot be broadxcasted to the same shape
                pass
            pass
        self._replotLock = None
        
        if plot: self.updateData(rescale)
        pass
        
    def _changeVertical(self, text):
        print("_changeVertical", text, "H", self.cbHorizontal.currentText())
        if text=="" or text is None: return
        if self._replotLock is None: self._replotLock = "vertical" #this should prevent _changeHorizontal handler replotting, because we will replot once finished here
        
        #check that the two axes are not the same
        if self.cbHorizontal.currentText() == text:
            self.cbHorizontal.setCurrentText(self.cbVertical._prevText)
        
        if self.cbVertical._prevText is not None and self.cbVertical._prevText!=self.cbHorizontal.currentText():
            C = self._axesControls[self.cbVertical._prevText]
            print("_changeVertical enable", self.cbVertical._prevText)
            C.setEnabled(True)
            C.selector.setCurrentUnit(self.cbVerticalUnits.currentText())
            
        self._axesControls[text].setEnabled(False) #TODO: and hide
        
        #note that _replotLock should prevent _changeVerticalUnits from triggering replot
        self.cbVerticalUnits.clear()
        self.cbVerticalUnits.addItems(self._axesControls[text].units())
        self.cbVerticalUnits.setCurrentText(self._axesControls[text].unit())
        
        if self._replotLock is "vertical": self._replotLock = None
        self.updateData()
        pass
    
    def _changeHorizontal(self, text):
        print("_changeHorizontal", text, "V", self.cbVertical.currentText())
        if text=="" or text is None: return
        if self._replotLock is None: self._replotLock = "horizontal" #this should prevent _changeVertical handler replotting, because we will replot once finished here
        #check that the two axes are not the same
        if self.cbVertical.currentText() == text:
            self.cbVertical.setCurrentText(self.cbHorizontal._prevText)
        
        if self.cbHorizontal._prevText is not None and self.cbHorizontal._prevText!=self.cbVertical.currentText():
            C = self._axesControls[self.cbHorizontal._prevText]
            C.setEnabled(True)
            C.selector.setCurrentUnit(self.cbHorizontalUnits.currentText())
            
        self._axesControls[text].setEnabled(False) #TODO: and hide

        #note that _replotLock should prevent _changeHorizontalUnits from triggering replot
        self.cbHorizontalUnits.clear()
        self.cbHorizontalUnits.addItems(self._axesControls[text].units())
        self.cbHorizontalUnits.setCurrentText(self._axesControls[text].unit())
        
        if self._replotLock is "horizontal": self._replotLock = None
        self.updateData()
        pass
  
    def _flipVH(self, *args):
        self.cbVertical.setCurrentText(self.cbHorizontal.currentText()) #should be enough, rest is handled by _changeVertical
        pass

    def updateData(self, rescale=True):
        """
        call this on axes change
        replot axes and labels
        
        then replot data (probably call update)
        """
        print("NDContainerWidget.updateData")
        print("NDContainerWidget.updateData: rescale", rescale)
        if self._replotLock is not None: return
        if self._data is None: return
        print("NDContainerWidget.updateData: plotting")


        multiIndex = []
        for axisControl in self._axesControls:
            C = self._axesControls[axisControl]
            if C.isEnabled():
                multiIndex.append((C.axisID, C.position(), C.unit()))
        
        cut, points = self._data.multiCut(multiIndex)
        
        #plot axes
        vAxisID = self.cbVertical.currentText()
        vAxisUnits = self.cbVerticalUnits.currentText()
        
        hAxisID = self.cbHorizontal.currentText()
        hAxisUnits = self.cbHorizontalUnits.currentText()
        
        #we do not know if vAxisID is first or second in cut now
        
        if cut.axes()[0] != vAxisID:
            print("swapping axes")
            cut.swapAxes(0, 1)
            #after this first axis of cut should be vAxis, second hAxis and there should be no other
        data, axes = cut.export((vAxisUnits, hAxisUnits))
        #valueAxis je prvni rozmer data
        
        
        valueAxis, keyAxis = axes
        valueLabel = cut.axis(0).label
        keyLabel = cut.axis(1).label
        zlabel = cut._valueLabel
        
        self._cutData = data, valueAxis, keyAxis, valueLabel, keyLabel, zlabel
        self.cutChanged.emit(self._cutData)
        print("NDContainerWidget.updateData: rescale", rescale)

        if rescale=="onAxesChange":
            rescale = self._cutAxesChanged
        elif rescale is False:
            pass
        else:
            rescale = True
            
        if rescale: 
            #move tracer to the middle
            #TODO
#            self.plot.setTracerPosition(0.5, 0.5, "rect")
            pass
        
        print("NDContainerWidget.updateData: rescale", rescale)
        self.plot.setData(data, keyAxis, valueAxis, xlabel=keyLabel, ylabel=valueLabel, zlabel=zlabel, rescale=rescale)
#        self.plot.replot()
        
#        self._tracerMoveIndex()
        
        print("NDContainerWidget._replot end", time.perf_counter())
        pass
        
    def _tracerMoveIndex(self, ix=None, iy=None):
        if self._cutData is None: return
        if ix is None or iy is None:
            #triggered after unit/data change
            ix, iy = self.plot.tracerPosition(True)

        data, valueAxis, keyAxis, valueLabel, keyLabel, zlabel = self._cutData

        #this should never be outside of the data range
        if ix < 0: ix=0
        if ix >= data.shape[1]: ix = data.shape[1]-1
        if iy < 0: iy=0
        if iy >= data.shape[0]: iy = data.shape[0]-1
        
        #emit tracesChanged
        V = (valueAxis, data[:, ix], valueLabel, zlabel)
        H = (keyAxis, data[iy], keyLabel, zlabel)
        self.tracesChanged.emit((V, H))
        pass
    
    pass

if __name__=="__main__":
    #~ #create simple app containing just the dialog
    #~ app = QtWidgets.QApplication(sys.argv)
    #~ window = QtWidgets.QMainWindow()
    #~ w = CutPlot(window)
    #~ window.setCentralWidget(w)
    #~ window.show()
    
    #~ #data
    #~ import numpy as np
    #~ data = np.outer(np.arange(10), np.arange(15))
    #~ x1 = np.arange(15)*2
    #~ x2 = x1[::-1]
    #~ x3 = x2 + 4.1
    #~ y = np.arange(10)+23.1
    #~ y2 = y[::-1]
    #~ y3 = y2+4
    
    #~ w.setData(data, (x1, x2, x3), (y, y2, y3))
    
    #create simple app containing just the dialog
    app = QtWidgets.QApplication(sys.argv)
    window = QtWidgets.QWidget()
    w = NDContainerWidget(window)
    
    p = QtWidgets.QPushButton("Load data")
    
    layout = QtWidgets.QVBoxLayout(window)
    layout.addWidget(w)
    layout.addWidget(p)
    
    def load():
        name = QtWidgets.QFileDialog.getOpenFileName(window)
        if name[0] != "":
            w.setData(NDContainer.NDContainer.fromFile(name[0]))
    
    p.clicked.connect(load)
    
    window.show()
    
    #data
    import numpy as np
    import common.NDContainer as NDContainer
    
    data = np.outer(np.arange(10), np.arange(15))[:, :, None]* np.arange(10)
    print(data.shape)
    d = NDContainer.NDContainer("ahoj")
    d.addAxis(NDContainer.RangeAxis(10, "X", "xlabel"))
    d.addAxis(NDContainer.RangeAxis(15, "Y", "ylabel"))
    #d.addAxis(NDContainer.RangeAxis(5, "Z", "zlabel"))
    from common.conversion import Convertor
    C = Convertor()
    C.setup("px", (0, 100), None)
    d.addAxis(NDContainer.ScanAxis(10, np.arange(10), np.arange(3, 13), C, C, C, "Z", ["zlabel", "z2", "z3"]))
    d.setData(data)
    
    #w.setData(data, (x1, x2, x3), (y, y2, y3))
    w.setData(d)
    
    
    app.exec_()
    
