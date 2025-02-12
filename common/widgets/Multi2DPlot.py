
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

#import PyQCustomPlot as QCP
from PyQt5 import QtCore

#TODO: needed for DDphasing

#TODO: does not work will PluginPlot.DD contrast colors

#from .Tracer2DPlot import Tracer2DPlot, QCP
# ~ from common.widgets.PluginPlot import Tracer2DPlot, QCP

#transition from QCP to pqr, possibly not the best approach :)
from pqr import pqr

class Multi2DPlot(pqr.FigureWidget):
    tracerMoveIndex = QtCore.pyqtSignal(int, int)
    def __init__(self, *args):
        super().__init__(*args, addFirstPlot=False)
        plot = self.addPlot(0,0)
        #we can keep multiple DDPlotters, but we need to keep track of multiple axes too
        self._data = {}
        self._currentData = None
        
        self.tracer = plot.addDecorator(pqr.TracerDecorator(plot))
        self.tracer.originChangedL.connect(self._emitIndex)
        # ~ self.tracerMove.connect(self._emitIndex)
        
    def dataNames(self):
        return list(self._data.keys())

    def data(self, name):
        return self._data[name]

    def addData(self, data, x, y, name, xlabel=None, ylabel=None):
        
        plot = self.plot()
        xaxis = plot.addAxis(pqr.NumEnumTransform(x), pqr.AxisPosition.B, label=str(xlabel)+name)
        yaxis = plot.addAxis(pqr.NumEnumTransform(y), pqr.AxisPosition.L, label=str(ylabel)+name)
        zaxis = plot.addAxis(None, pqr.AxisPosition.R, axisClass=pqr.ColorAxisArea, globalZoom=False)
        dd = plot.addDDPlotter(xaxis, yaxis, zaxis)
        dd.setData(data.T, x, y)

        self._data[name]={'cm':dd, 'x':x, 'y':y, 'xlabel':xlabel, 'ylabel':ylabel, 'data':data, "xaxis":xaxis, "yaxis":yaxis, "zaxis":zaxis} #TODO: it is not good to keep data twice (once in cmdata, second time in data)
        for it in ["cm", "xaxis", "yaxis", "zaxis"]:
            self._data[name][it].setVisible(False)
        
        if self._currentData==name: self.switchData(name)
        pass

    def removeData(self, name):
#         if name=='default':
#             raise IndexError('Default data cannot be removed')
        if name==self._currentData:
            raise IndexError('Current data cannot be removed')
        del self._data[name]
        #TODO: this should remove DDPlotter and axes too
        #use this sparingly
        pass
 
    def switchData(self, name):
        #print("Multi2DPlot.switchData", name)
        #hide present
        if self._currentData is not None:
            old = self._data[self._currentData]
            for it in ["cm", "xaxis", "yaxis", "zaxis"]:
                old[it].setVisible(False)
        
        new = self._data[name]
        for it in ["cm", "xaxis", "yaxis", "zaxis"]:
            new[it].setVisible(True)
        
        self.tracer.setAxes(new["xaxis"], new["yaxis"])
        
        self.plot().zoomToData()
        #TODO: possibly update tracer position

        self._currentData=name
        pass


    """this should not be needed for pqr, right?
    def plot(self, name=None, centerTracer=True):
        #print( "Multi2DPlot.plot", name, centerTracer, self._currentData)
        if name!=None and self._currentData!=name: self.switchData(name)

        #self.rescaleAxes(True)
        #print( self.xAxis.range().lower, self.xAxis.range().upper)
        #self.xAxis.rescale()
        #print( self.xAxis.range().lower, self.xAxis.range().upper)
        #print(self._data[name]["x"])
        #self.yAxis.rescale()
        #but not the other two
        #TODO: actually after change to PluginPlot it does not work without rescaleAxes (all of them) - probably zoom guard in Zoom (have to mark axis as pending rescale to temporarily disable the guard)
        #   but this will rescale badly window marker in plotRaw
        # ~ self.rescaleAxis(self.xAxis, True)
        # ~ self.rescaleAxis(self.yAxis, True)
        self.plot().zoomToData()
        
        #move cursor to the middle
        self.tracer.pos = 0.5, 0.5
        # ~ ttt=self.tracer.position
        # ~ if centerTracer:
            # ~ ttt.setType(self.tracer.position.ptAxisRectRatio)
            # ~ ttt.setCoords(0.5,0.5)
            # ~ ttt.setType(self.tracer.position.ptPlotCoords)
            # ~ pass
        # ~ self.replot()
        #this should also emit cursor move signals
        # ~ key=ttt.key()
        # ~ value=ttt.value()
        # ~ self.tracerMove.emit(key, value)
        pass
    """
    
    def _emitIndex(self, key, value):
        # ~ cmdata = self.cm.data()
        # ~ nki, nvi = cmdata.coordToCell(key,value)
        # ~ #ensure to be inside valid index range
        # ~ if nki >=0 and nki < cmdata.keySize() and nvi >= 0 and nvi < cmdata.valueSize():
            # ~ self.tracerMoveIndex.emit(nki, nvi)
        # ~ pass
        try:
            self.tracerMoveIndex.emit(int(key), int(value))
        except:
            import traceback
            print("Multi2DPlot._emitIndex: encountered exception", key, value)
            print(traceback.print_exc())
            pass

    pass

"""
old QCP version
class Multi2DPlot(Tracer2DPlot):
    tracerMoveIndex=QtCore.pyqtSignal(int, int)
    def __init__(self, *args):
        Tracer2DPlot.__init__(self, *args)
        self._data={'default':{'cm':self.cm.data(), 'x':[], 'y':[], 'xlabel':None, 'ylabel':None, 'data':[]}}
        self._currentData='default'
        
        self.tracerMove.connect(self._emitIndex)
        pass

    def dataNames(self):
        return list(self._data.keys())

    def data(self, name):
        return self._data[name]

    def addData(self, data, x, y, name, xlabel=None, ylabel=None):
        #create QCPColorMapData from data
        #print("Multi2DPlot.addData", data, x, y, name)
        cmdata=QCP.QCPColorMapData(data.shape[0], data.shape[1], QCP.QCPRange(x[0], x[-1]),QCP.QCPRange(y[0], y[-1]) )
        #TODO: do something about setting data one by one
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                cmdata.setCell(i,j,data[i,j])
                pass
            pass
        self._data[name]={'cm':cmdata, 'x':x, 'y':y, 'xlabel':xlabel, 'ylabel':ylabel, 'data':data} #TODO: it is not good to keep data twice (once in cmdata, second time in data)
        if self._currentData==name: self.switchData(name)
        pass

    def removeData(self, name):
#         if name=='default':
#             raise IndexError('Default data cannot be removed')
        if name==self._currentData:
            raise IndexError('Current data cannot be removed')
        del self._data[name]
        pass

    def switchData(self, name):
        #print("Multi2DPlot.switchData", name)
        self.cm.setData(self._data[name]['cm'], True) #I have to copy, otherwise it will delete the old one :(
        if self._data[name]['xlabel']!=None: self.xAxis.setLabel(self._data[name]['xlabel'])
        if self._data[name]['ylabel']!=None: self.yAxis.setLabel(self._data[name]['ylabel'])
        self.cm.rescaleDataRange(True)
        #print("\t", self._data[name]['cm'].keyRange().lower, self._data[name]['cm'].keyRange().upper)

        ttt=self.tracer.position
        key=ttt.key()
        if not self._data[name]['cm'].keyRange().contains(key):
            key=self._data[name]['cm'].keyRange().center()

        value=ttt.value()
        if not self._data[name]['cm'].valueRange().contains(value):
            value=self._data[name]['cm'].valueRange().center()

        ttt.setCoords(key, value)

        self._currentData=name
        pass

    def plot(self, name=None, centerTracer=True):
        #print( "Multi2DPlot.plot", name, centerTracer, self._currentData)
        if name!=None and self._currentData!=name: self.switchData(name)

        #self.rescaleAxes(True)
        #print( self.xAxis.range().lower, self.xAxis.range().upper)
        #self.xAxis.rescale()
        #print( self.xAxis.range().lower, self.xAxis.range().upper)
        #print(self._data[name]["x"])
        #self.yAxis.rescale()
        #but not the other two
        #TODO: actually after change to PluginPlot it does not work without rescaleAxes (all of them) - probably zoom guard in Zoom (have to mark axis as pending rescale to temporarily disable the guard)
        #   but this will rescale badly window marker in plotRaw
        self.rescaleAxis(self.xAxis, True)
        self.rescaleAxis(self.yAxis, True)

        #move cursor to the middle
        ttt=self.tracer.position
        if centerTracer:
            ttt.setType(self.tracer.position.ptAxisRectRatio)
            ttt.setCoords(0.5,0.5)
            ttt.setType(self.tracer.position.ptPlotCoords)
            pass
        self.replot()
        #this should also emit cursor move signals
        key=ttt.key()
        value=ttt.value()
        self.tracerMove.emit(key, value)
        pass
    
    def _emitIndex(self, key, value):
        cmdata = self.cm.data()
        nki, nvi = cmdata.coordToCell(key,value)
        #ensure to be inside valid index range
        if nki >=0 and nki < cmdata.keySize() and nvi >= 0 and nvi < cmdata.valueSize():
            self.tracerMoveIndex.emit(nki, nvi)
        pass

"""
