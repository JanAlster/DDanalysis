
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

#we need a view of 3D data
#do not store the data, keep them elsewhere and just not the reference
#use PluginPlots to show 2D cuts

from PyQt5 import QtWidgets, QtCore, QtGui
from widgets import PluginPlot, DoubleWidget

class MainPlot(PluginPlot.Saving, PluginPlot.FreeZoom, PluginPlot.DD, PluginPlot.Base):
    pass

#you can only set the labels at init

#TODO: connect axes to convertors and take units from convertors


class DDDPlot(QtWidgets.QFrame):
    def __init__(self, *args):
        super().__init__(*args)

        self._data = None
        self._axes = [None, None, None]
        self._labels = ["X", "Y", "Z"]
        self._lastIndices = [0,0,0]
        
        self._plot = MainPlot(self)
        #self._plot.addGraph()
        self._rescalePlot = False

        index = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        index.sliderMoved.connect(self._moveIndex)
        index.valueChanged.connect(self._changeIndex)
        self._index = index

        value = DoubleWidget.DoubleWidget() #TODO: this could be position, but needs to be reconnected to different convertor based on selected axis
        value.valueChanged.connect(self._changeValue)
        self._value = value

        cut = QtWidgets.QComboBox()
        cut.currentIndexChanged[int].connect(self._changeCut)
        cut.addItems(self._labels)
        self._cut = cut

        l = QtWidgets.QVBoxLayout()
        l.addWidget(self._plot, 1)
        c = QtWidgets.QHBoxLayout()
        c.addWidget(cut)
        c.addWidget(value)
        c.addWidget(index, 1)
        l.addLayout(c)
        
        self.setLayout(l)
        
        #for testing
        import numpy as np
        ax = np.arange(5)
        ay = np.arange(10)
        az = np.arange(20)
        data = np.ones((5, 10, 20)) * np.outer(ay, az)
        data = (data.T * ax).T
        self.setup(data=data, x=ax, y=ay, z=az, xlabel="ble", ylabel="haf", zlabel="huu")
        pass
    
    #we need to setup for every axis
    #   - label
    #   - all values (could go with setData)
    
    #we should also allow setting/updating just one axis
    #changing code does not need replot
    #typically we
    #   - set all labels
    #   - set data and all axes values
    #   - set one axis label and value (change axis unit)
    
    def setup(self, **settings):
        """valid parameters are
            - data - 3D data adarray
            - x, y, z - axes, must be of correct length ot match data.shape
            - xlabel, ylabel, zlabel - label (with units) of axes
            TODO: datalabel - label for data amplitude on colorscale
        """
        print("DDDPlot.setup", settings.keys())
        replot = False
        rescale = False

        #labels
        lkey = ["xlabel", "ylabel", "zlabel"]
        for i, lkey in enumerate(lkey):
            if lkey in settings: 
                self._labels[i] = settings.pop(lkey)
                self._cut.setItemText(i, self._labels[i])
                #replot if i is not the current cut (i.e. label on plot changed)
                if i == self._currentKey: 
                    self._plot.setup(xlabel = self._labels[i])
                    replot = True
                elif i == self._currentValue: 
                    self._plot.setup(ylabel = self._labels[i])
                    replot = True
                pass
            pass
            
        #data before axes
        newData = False
        if "data" in settings:
            #this is the only place where we can change dimension
            #TODO: sanity check
            self._data = settings.pop("data")
            newData = True
            #newData implicates replot and rescale
            
            #reset index memory
            self._lastIndices = [0,0,0] 
            pass
            
        #axes
        akey = ["x", "y", "z"]
        for i, akey in enumerate(akey):
            newAxis = False
            if akey in settings:
                axis = settings[akey]
                try:
                    assert(len(axis) == self._data.shape[i])
                    self._axes[i] = axis
                    newAxis = True
                except Exception as e:
                    print("DDDPlot.setup: cannot set axis", akey, "- action failed with:\n", e)
                    pass
                pass
            
            if newData and not newAxis:
                #check old axis
                try:
                    assert(len(self._axes[i]) == self._data.shape[i])
                    #leave old axis
                    pass
                except Exception as e:
                    print("DDDPlot.setup: cannot use old axis", akey, " for new data - action failed with:\n", e)
                    #create new index axis
                    self._axes[i] = None if self._data is None else list(range(self._data.shape[i]))
                    newAxis = True
                    pass
                    
            if newAxis:
                if i == self._currentCut:
                    if self._axes[self._currentCut] is not None:
                        self._value.setValue(self._axes[self._currentCut][self._index.value()])
                    pass
                else:
                    replot = True
                pass
            pass
            
        #after axes are sorted out
        if newData:
            #update self._index range
            self._index.setMaximum(0 if self._axes[self._currentCut] is None else len(self._axes[self._currentCut])-1)
            #update self._index position (reset)
            self.setIndex(0) #this will replot and rescale
            pass
            
        if replot:
            self._update(rescale)
        pass
    
    def clearData(self):
        self.setup(data=None, x=None, y=None, z=None) #TODO: this will likely fail
        pass
    
    def _update(self, rescale=False):
        if rescale:
            self._plot.rescaleAxes()
        self._plot.replot() #we might want to update axes labels even if there are not data yet
        pass
        
    def _moveIndex(self, index):
        #response to index slider movement
        print("DDDPlot._moveIndex", index, self._lastIndices, self._currentCut)
        self._lastIndices[self._currentCut] = index
        if self._data is not None:
            #self._plot.setData(self._data[self._indexFactory(index)], self._keyAxis, self._valueAxis)
            self._plot.setData(self._data[self._indexFactory(index)], self._axes[self._currentKey], self._axes[self._currentValue])
        self._update()
        pass
        
    def _changeIndex(self, index):
        #react to any change of index slider (movement or from outside)
        print("DDDPlot._changeIndex", index)
        if self._axes[self._currentCut] is not None:
            self._value.setValue(self._axes[self._currentCut][index])
        pass    
        
    def setIndex(self, index):
        #set index from outside
        print("DDDPlot.setIndex", index, self._lastIndices, self._currentCut)
        self._index.setValue(index) #this will limit index to valid range and update self._value (via valueChanged and _changeIndex)
        self._lastIndices[self._currentCut] = self._index.value()
        if self._data is not None:
            #self._plot.setData(self._data[self._indexFactory(index)], self._keyAxis, self._valueAxis)
            self._plot.setData(self._data[self._indexFactory(index)], self._axes[self._currentKey], self._axes[self._currentValue])
        self._update(True)
        pass
        
    def _changeValue(self, value):
        #change index to corespond to value (as well as possible)
        print("DDDPlot._changeValue", value)
        pass
        
    def _changeCut(self, index):
        print("DDDPlot._changeCut", index)
        #order is x, y, z
        self._currentCut = index
        if index == 0:
            self._indexFactory = lambda t: t
            self._currentKey = 2
            self._currentValue = 1
        elif index == 1:
            self._indexFactory = lambda t: (slice(None), t)
            self._currentKey = 2
            self._currentValue = 0
        elif index == 2:
            self._indexFactory = lambda t: (slice(None), slice(None), t)
            self._currentKey = 1
            self._currentValue = 0
        else:
            #complain
            raise NotImplementedError()
        
        self._plot.setup(xlabel = self._labels[self._currentKey], ylabel = self._labels[self._currentValue])
  
        self._index.setMaximum(0 if self._axes[self._currentCut] is None else len(self._axes[self._currentCut])-1)
        self.setIndex(self._lastIndices[index])
        pass
    pass
