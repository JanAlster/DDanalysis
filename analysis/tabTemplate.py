
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


from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.uic import loadUi

from analysis import data_pipe
from common.SaveLoadElement import SaveLoadElement
from common import sigmoidWindow, fitting
from common.interpolator import interpolate
#from common.UIGroup import UIGroup
import os
import os.path as op
import numpy as np

import warnings

from pqr import pqr2 as pqr
#~ import PyQCustomPlot as QCP

#it will keep local copy of main.data
# with applied current PAIR, TRN, interpolation 
# and show that in timeline and dissection plots
# keeping tract of disabled frames
# of course, other might want to access this data ... 

#e.g. tabGA might want to use original main.data, original main.data with some frames disabled, or interpolated and manipulated data from tabData

#another possibility is to interpolate only data in plotMain (diccestion) and let timeline uninterpolated
# - I dont see a point of doing GA on interpolated data, we can just interpolate results of GA
# - then PAIR, TRN and interpolation can be pushed to plotMain, but not global colors ...

#sublimation
# - inputData (all TRN)
# - showData (interpolated, PAIR)
# - outputData
# tabData     (input: original data, show: original data interpolated, output: original data (selection))
# tabGA       (input: original data, show: DAS interpolated, output: model)
# tabResidues (input: original data, model, show: residues interpolated, output: residues (selection))
# tabFreq     (input: residues, show: OM interpolated, output: )


def StyleSheet2QColor(stylesheet, key='color'):
    stylesheet=stylesheet.replace(' ','')
    i=stylesheet.find(key+':')+1+len(key)
    return QtGui.QColor(stylesheet[i:stylesheet.find(';', i)])

def StyleSheet2Color(stylesheet, key='color'):
    stylesheet=stylesheet.replace(' ','')
    i=stylesheet.find(key+':')+1+len(key)
    return stylesheet[i:stylesheet.find(';', i)]

def midunwrap(a, axis=-1):
    N = a.shape[axis]//2
    t = a.take(N, axis=axis)
    a = np.unwrap(a, axis=axis)
    s = np.array(a.shape)
    s[axis] = 1
    t -= a.take(N, axis=axis)
#    print(a.shape, axis, t.shape, s)
    a += t.reshape(s) #but we likely need some magic here
    return a

def midunwrap2d(a):
    return midunwrap(midunwrap(a, axis=0))


TRN2trn={'Total':'total', 'Rephasing':'rephasing', 'Non-Rephasing':'nonrephasing'}
trn2TRN={v: k for k, v in TRN2trn.items()}


# testing
figure = None


# def _AmplToNoiseOld(x):
#     global figure, plotter
#     if figure is None:
#         figure = pqr.FigureWidget()
#         plot = figure.plot()
#         plotter = plot.addDDPlotter()
#         figure.show()
#
#     ax = np.abs(x)
#     # percentile does not work, because different size of spectrum might be covered by noise, we need to figure out
#     # the average noise level somehow
#     plotter.setData(ax)
#     print(ax)
#     return ax / np.percentile(ax, 2)

# def _AmplToNoise(x):
#     print("_AmplToNoise", x.shape)
#     N = 2
#     sr = np.empty((x.shape[1] - N, x.shape[2] - N), float)
#     si = np.empty((x.shape[1] - N, x.shape[2] - N), float)
#     S = []
#     res = np.abs(x)
#     for k in range(x.shape[0]):
#         data = x[k]
#         for i in range(sr.shape[0]):
#             for j in range(sr.shape[1]):
#                 sr[i, j] = data.real[i:i + N, j:j + N].std()
#                 si[i, j] = data.imag[i:i + N, j:j + N].std()
#
#         def filter(y, M=1):
#             for i in range(M):
#                 c, edges = np.histogram(y, bins=20, range=(np.nanmin(y), np.nanmax(y)))
#                 ci = c.argmax()
#                 my = 0.5 * (edges[ci] + edges[ci + 1])
#                 sy = np.nanmean(y)
#                 y[y > my + 1.5 * sy] = np.NaN
#
#         filter(sr, 5)
#         filter(si, 5)
#
#         s = (np.nanmean(sr)**2+np.nanmean(si)**2)**0.5
#         print("_AmplToNoise", k, s)
#         res[k] /= s
#     return res


#TODO: unwrap phase from maximum amplitude
PAIR={'Real': lambda x: x.real,
      'Imag': lambda x: x.imag,
      'Ampl': lambda x: np.abs(x),
      'Phase': lambda x: np.angle(x),
      None:lambda x:x,
      "PhaseA": lambda x: (np.angle(x), np.abs(x)),
      "PhaseU": lambda x: midunwrap2d(np.angle(x)),
      "PhaseUA": lambda x: (midunwrap2d(np.angle(x)), np.abs(x)),
      "Ampl/Noise": data_pipe.amplitude_to_noise,
      }

class tabTemplateSmall(QtWidgets.QWidget, SaveLoadElement):
    outputChanged = QtCore.pyqtSignal(dict)
    def __init__(self, main, form, *args): 
        super().__init__(*args)
        loadUi(form, self)

        self.boss = main #link to other things
        self.resetData() #reset/prepare data buffers
        
        #plotting
        #windowing main plot (first plot is added automatically)
        plot = self.plotMain.plot()
        xaxis = plot.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.B, label='ω₁ / cm⁻¹', tickFormat="{:.0f}")
        yaxis = plot.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.L, label='ω₃ / cm⁻¹', tickFormat="{:.0f}")
        xaxis2 = plot.addAxis(pqr.AOverTransform(1e7, xaxis._transform), pqr.AxisPosition.top, label="λ₁ / nm")
        yaxis2 = plot.addAxis(pqr.AOverTransform(1e7, yaxis._transform), pqr.AxisPosition.right, label="λ₃ / nm")
        
        xaxis.zoomRangeLChanged.connect(xaxis2._setZoomRange)
        xaxis2.zoomRangeLChanged.connect(xaxis._setZoomRange)
        yaxis.zoomRangeLChanged.connect(yaxis2._setZoomRange)
        yaxis2.zoomRangeLChanged.connect(yaxis._setZoomRange)

        plotter = plot.addDDPlotter(xaxis, yaxis)
        tracer = plot.addDecorator(pqr.TracerDecorator(plot, xaxis, yaxis))
        diagonal = plot.addSeriesPlotter(xaxis, yaxis, width=0.5, name="diagonal")
        
        #~ plotter.contours = 12
        self.cbPAIR.clear()
        self.cbPAIR.addItems([it for it in PAIR.keys() if it is not None])
        
        #setup SaveLoadElement(s)
        #NOTE that order here is important, because setting UI elements will trigger calculation of intermediate data automatically
        #no UI elements affect internal raw data (in template)
        self.addGroup("show", self.updateShow)
        #these will recalculate (TRN, PAIR, interpolation) data buffer for UI elements
        self.addGroup("buffer", self._updateDataBuffer) 
        self.addEditable("TRN", self.cbTRN.currentText, self.cbTRN.setCurrentText, self.cbTRN.currentTextChanged, group="buffer", initialize=True)
        self.addEditable("PAIR", self.cbPAIR.currentText, self.cbPAIR.setCurrentText, self.cbPAIR.currentTextChanged, group="buffer", initialize=True)
        self.addEditable('I1', self.iInterpolation1.value, self.iInterpolation1.setValue, self.iInterpolation1.valueChanged, group="buffer", initialize=True)
        self.addEditable('I3', self.iInterpolation3.value, self.iInterpolation3.setValue, self.iInterpolation3.valueChanged, group="buffer", initialize=True)
        #UI elements controlling individual plots
        self.addEditable("globalColors", self.bGlobalColors.isChecked, self.bGlobalColors.setChecked, initialize=True)
        self.addEditable("fixColors", self.bFixColors.isChecked, self.bFixColors.setChecked, initialize=True)

        #inform downstream that data selection changed (UI triggered)
        self.addGroup("output", self.updateOutput)
        self.addEditable("disabledFrames", self.listFrames.disabled, self.listFrames.setDisabled, self.listFrames.disabledChanged, group="output", initialize=True)

        self.addEditable("currentFrame", self.listFrames.current, self.listFrames.setCurrent, initialize=True)


        #connect elements
        #these will update UI elements (with data from precalculated data buffer)
        self.bGlobalColors.stateChanged.connect(self.updatePlotMain)
        self.listFrames.currentChanged.connect(self.updatePlotMain)
        
        #Todo: possibly settings of plots                
        pass

    def disabled(self):
        return self.listFrames.disabled()
        
    def enabled(self):
        return self.listFrames.enabled()

    def saveSettings(self, saveData=False):
        settings = super().saveSettings() #those setup SaveLoadElement(s) above
        if saveData and self._dataShow is not None:
            #I want to save uninterpolated dataShow
            settings["data"] = self._dataShow
        return settings

    #saved data will not be loaded, they are used for something else
    #~ def loadSettings(self, settings):
        #~ print("tabTemplate.loadSettings")
        #~ super().loadSettings(settings) #those setup SaveLoadElement(s) above
        #~ #plots
        #~ pass
        
    def resetData(self):
        self._dataInput = None #input data, processed up to the point that does not depend on UI
        self._dataShow = None #internal data that are shown in main plot (raw data)
        self._dataOutput = None #output data 
        #above data are {"total": ndarray, "rephasing":..., "nonrephasing":..., "a1":..., "a2":..., "a3":...}
        self._currentData = None #internal data that are shown in main plot (interpolated, TRN, PAIR)
        pass  
        
        
    #for export dialog
    def dataExport(self):
        return self._dataShow
        
    def updateInput(self, data=None):
        #get data from outside
        #proccess up to the point where input from UI is needed
        #store to self._dataInput
        raise NotImplementedError

    def updateShow(self):
        #~ print("tabTemplate.updateShow", not self.isSet("show"), self._dataInput is None)
        #get input from UI needed to calculate data to show in plots, etc..
        #store ti self._dataShow
        #(assumption is that correct UI state is loaded before this is called)
        if not self.isSet("show"): return 
        #if no input is needed use the default implementation
        if self._dataInput is None: return #in case some relevant UI element is set before we have data
        self._dataShow = self._dataInput
        #print("\tsetting listFrames", self._dataShow["a2"])
        self.listFrames.setItems(self._dataShow["a2"], force=False) #to not reset the listFrame if we are changing UI params, but data remain the same
        self._updateDataBuffer()
        self.updateOutput()
        pass
        
    def updateOutput(self):
        #~ print("tabTemplate.updateOutput", not self.isSet("output"), self._dataShow is None)
        #should be triggered every time input data changes or relevant UI settings changes
        #get input from UI needed to calculate output data 
        #(mostly this is a selection from shown data)
        # if so, you can use the default implementation
        if not self.isSet("output"): return #you should use this in reimplementation
        if self._dataShow is None: return #in case some relevant UI element is set before we have data
        #limit output to enabled frames
        enabled = self.enabled()
        self._dataOutput = {"a1":self._dataShow["a1"], "a3":self._dataShow["a3"]}
        for it in ["total", "rephasing", "nonrephasing", "a2"]:
            try:
                self._dataOutput[it] = self._dataShow[it][enabled]
            except KeyError:
                # this is for backward compatibility with badly saved data from previous buggy version
                print("warning at tabData.updateShow: something went wrong")
                continue

        if "LO" in self._dataShow and self._dataShow["LO"] is not None:
            self._dataOutput["LO"] = self._dataShow["LO"][enabled]

        self.outputChanged.emit(self._dataOutput)
        pass
        
    def _updateDataBuffer(self):
        #most of the data proccessing before it can be shown is time consuming (interpolation, ...)
        # so it is best to precalculate the data buffer
        # and only select the correct part in updatePlot methods
        print("tabTemplateSmall._updateDataBuffer", not self.isSet("buffer"), self._dataShow is None)
        if not self.isSet("buffer"): return
        if self._dataShow is None: return         
        trn = TRN2trn[self.cbTRN.currentText()]
        pair = self.cbPAIR.currentText()
        i1 = self.iInterpolation1.value()
        i3 = self.iInterpolation3.value()


        #getting the correct dataRange is not easy - we need to reflect PAIR selection
        # - i.e. max is not always from absolute value, but could be just the biggest real part
        # - but that could mean to apply PAIR transfrom on the whole dataset all the time :(
        #TODO: we could hold a copy of data (not interpolated) for every PAIR (more memory, but uninterpolated data are usually not that big)

        data = PAIR[pair](self._dataShow[trn])
        a1 = self._dataShow["a1"]
        a2 = self._dataShow["a2"]
        a3 = self._dataShow["a3"]
        
        #handle phaseA and phaseU
        if pair=="PhaseA" or pair=="PhaseUA":
            data, temp = data
            #toto transformace by mohla byt trochu nelinearni, utlumit jen male maplitudy, kde se fazi neda verit
            #alpha = temp * (255/temp.reshape((temp.shape[0], -1)).max(axis=1).reshape((-1, 1, 1)))
            temp *= 255/temp.reshape((temp.shape[0], -1)).max(axis=1).reshape((-1, 1, 1))
            # ~ alpha = (temp > 50) * 255
            # ~ del temp
            alpha = temp
        else:
            alpha = None
        
        if pair.startswith("Phase"):
            #TODO: use colorcet linear perception colormaps here, allow user selection
            self.plotMain.plot().y3.colorMap = pqr.cmRGBCycle #this should be colorAxis
        else:
            self.plotMain.plot().y3.colorMap = pqr.cmDeepBlueRed
        
        #for dataRange we need to keep infs and nans out
        temp = data[np.isfinite(data)]
        if temp.size==0:
            dataRange = (-1, 1)
        else:
            dataRange = (temp.min(), temp.max()) #global range

        if i1>0:
            data, a1 = interpolate(data, i1, a1, axis=1)
        if i3>0:
            data, a3 = interpolate(data, i3, a3, axis=2)

        self._currentData = data, a1, a3, a2, dataRange, alpha

        self.updatePlotMain()
        pass
        
    def updatePlotMain(self, *args):
        #TODO: main plot should keep zoom/range state of color scale when just changing current (or all the time?)
        #  rescale color scale only on updateBuffer
        print("tabTemplateSmall.updatePlotMain", "no current data", self._currentData is None)
        if self._currentData is None: 
            print("tabTemplateSmall.updatePlotMain", "no current data")
            #TODO: clear plot
            return
        index = self.listFrames.current()
        print("\t\tindex", index)
        if index is None or index==-1:
            print("tabTemplateSmall.updatePlotMain", "bad index", index)
            warnings.warn(' is not yet calculated')
            return

        data, a1, a3, a2, dataRange, alpha = self._currentData
        if index>len(data): 
            print("WARNING: tabTemplateSmall.updatePlotMain - mismatched index and data")
            return
        
        frame = data[index]
        if alpha is not None: alpha = alpha[index].T
        # ~ print("\t\tframe", frame)
        # ~ print("\t\tframe all nan", np.all(np.isnan(frame)))


        # ~ print("\t", index, frame)
        if self.bGlobalColors.isChecked():
            a = max(abs(dataRange[0]), abs(dataRange[1]))
        else:
            #for dataRange we need to keep infs and nans out
            temp = frame[np.isfinite(frame)]
            if temp.size==0:
                a = 1
            else:
                a = max(abs(temp.min()), abs(temp.max()))
        dataRange = (-a, a)
        # ~ print("\t\tdata range", -a, a)
        #TODO: contours do not work properly from some time (not sure what happend), they display at 0,0 not at axis coordinates
        plot = self.plotMain.plot()
        plotter = plot.plotter(0)

        plotter._xaxis._transform.stops = a1
        plotter._yaxis._transform.stops = a3

        zoomDataRange = plotter._zaxis.linear2data(plotter._zaxis._zoomRange())
        
        plotter.setData(frame.T, a1, a3, dataRange=dataRange, alpha=alpha) #Todo: keep contours settings from previous (Traced2DPlot is originaly inteded for one time use, not interactive)
        #~ self.plotMain.replot()
        plot.zoomToData()

        if self.bFixColors.isChecked():
            plotter._zaxis._setZoomRange(*plotter._zaxis.data2linear(zoomDataRange))

        diagonal = plot.plotter(1)
        limits = (max(a1[0], a3[0]), min(a1[-1], a3[-1]))
        diagonal.setData(limits, limits)

        pass

    pass

class tabTemplate(tabTemplateSmall):
    def __init__(self, main, form, *args): 
        super().__init__(main, form, *args)
        #Dissection
        plot = self.plotCuts.plot()
        xaxis = plot.addAxis(None, pqr.AxisPosition.B, label='ω / cm⁻¹')
        xaxis2 = plot.addAxis(None, pqr.AxisPosition.T, label='ω / cm⁻¹')
        yaxis = plot.addAxis(None, pqr.AxisPosition.L)
        
        #Todo: diagonal and antidiagonal should be part of tracer
        #~ self.plotMain.addGraph().setPen(self.plotMain.tracer.pen())#diagonal tracer
        #~ self.plotMain.addGraph().setPen(self.plotMain.tracer.pen())#antidiagonal tracer
        
        plot.addSeriesPlotter(xaxis, yaxis, color=StyleSheet2Color(self.cbDisVertical.styleSheet())) #vertical
        plot.addSeriesPlotter(xaxis, yaxis, color=StyleSheet2Color(self.cbDisHorizontal.styleSheet())) #horizontal
        colorD = StyleSheet2QColor(self.cbDisDiagonal.styleSheet())
        colorA = StyleSheet2QColor(self.cbDisAntiDiagonal.styleSheet())

        plot.addSeriesPlotter(xaxis2, yaxis, color=colorD) #diagonal
        plot.addSeriesPlotter(xaxis2, yaxis, color=colorA) #antidiagonal
        
        plot.addSeriesPlotter(xaxis2, yaxis, color=colorD.darker(), style=QtCore.Qt.DashLine) #diagonal fit
        plot.addSeriesPlotter(xaxis2, yaxis, color=colorA.darker(), style=QtCore.Qt.DashLine) #antidiagonal fit
        

        self.addEditable("disHor", self.cbDisHorizontal.isChecked, self.cbDisHorizontal.setChecked, initialize=True)
        self.addEditable("disVer", self.cbDisVertical.isChecked, self.cbDisVertical.setChecked, initialize=True)
        self.addEditable("disDia", self.cbDisDiagonal.isChecked, self.cbDisDiagonal.setChecked, initialize=True)
        self.addEditable("disADia", self.cbDisAntiDiagonal.isChecked, self.cbDisAntiDiagonal.setChecked, initialize=True)
        self.addEditable("disDiaFit", self.cbDisDiagonalFit.currentText, self.cbDisDiagonalFit.setCurrentText, initialize=True)
        self.addEditable("disADiaFit", self.cbDisAntiDiagonalFit.currentText, self.cbDisAntiDiagonalFit.setCurrentText, initialize=True)

        self.listFrames.currentChanged.connect(self.updateDissectionSidePlots) 

        #connect elements
        self.cbDisVertical.clicked.connect(self.updateDissectionSidePlots)
        self.cbDisHorizontal.clicked.connect(self.updateDissectionSidePlots)
        self.cbDisDiagonal.clicked.connect(self.updateDissectionSidePlots)
        self.cbDisAntiDiagonal.clicked.connect(self.updateDissectionSidePlots)
        self.cbDisDiagonalFit.currentIndexChanged.connect(lambda x: self.updateDissectionSidePlots())
        self.cbDisAntiDiagonalFit.currentIndexChanged.connect(lambda x: self.updateDissectionSidePlots())

        self.cbDisDiagonal.toggled.connect(self.cbDisDiagonalFit.setEnabled)
        self.cbDisAntiDiagonal.toggled.connect(self.cbDisAntiDiagonalFit.setEnabled)
        
        tracer = self.plotMain.plot(0).decorator(0)
        tracer.originChangedL.connect(self.updateDissectionSidePlots) #if NumEnumTransform is used for axis, linear space should be index of axis NumEnumTransform.stops, which mightbe tied to data array index

        #Todo: possibly settings of plots                
        pass
    
    def _updateDataBuffer(self):
        super()._updateDataBuffer()
        self.updateDissectionSidePlots()
        pass

    def updateDissectionSidePlots(self, ix=None, iy=None):
        #~ print("tabTemplate.updateDissectionSidePlots", ix, iy, self._currentData is None)
        if self._currentData is None: 
            #TODO: clear plot
            return

        index = self.listFrames.current()
        if index is None: return

        if ix is None or iy is None:
            #~ ix,iy = self.plotMain.tracerPosition(True)
            #what we need here is originL value of tracer
            tracer = self.plotMain.plot().decorator(0)
            ix = tracer._xL
            iy = tracer._yL
            print("ix,iy", ix, iy)
            pass

        # FIXME: there is a possibility that tracer goes crazy and sends NaN here - we should restart it in that case
        ix = int(round(ix))
        iy = int(round(iy))
        
        """
        TODO:
        alternatively use:
        w1, w3 = tracer.origin
        
        trn = TRN2trn[self.cbTRN.currentText()]
        pair = self.cbPAIR.currentText()
        
        data = PAIR[pair](self._dataInput[trn])
        a1 = self._dataInput["a1"]
        a2 = self._dataInput["a2"]
        a3 = self._dataInput["a3"]
        
        ix = a1.searchsorted(w1)
        iy = a3.searchsorted(w3)        
        
        not sure which is better - I am not sure if linear space for DDPlotter is certain to be the index space for data array
         if we plot colors under different transform, it might be different
         - which is something that might happen once we are able to plot two DD on the same axes (that is likely to involve some interpolation)
        """
        
        data, a1, a3, a2, dataRange, alpha = self._currentData
        frame = data[index]
        if ix<0: ix = 0
        elif ix>=data.shape[1]: ix = data.shape[1]-1
        if iy<0: iy = 0
        elif iy>=data.shape[2]: iy=data.shape[2]-1

        try:
            #cuts
            #vertical
            graph=self.plotCuts.plot().plotter(0)
            if self.cbDisVertical.isChecked():
                graph.setVisible(True)
                graph.setData(a3, frame[ix])
            else:
                graph.setVisible(False)

            #horizontal
            graph=self.plotCuts.plot().plotter(1)
            if self.cbDisHorizontal.isChecked():
                graph.setVisible(True)
                graph.setData(a1, frame[:,iy])
            else:
                graph.setVisible(False)

            #diagonal
            graph=self.plotCuts.plot().plotter(2)
            #~ graph2=self.plotMain.plot().plotter(0)
            graph3=self.plotCuts.plot().plotter(4)
            if self.cbDisDiagonal.isChecked():
                graph.setVisible(True)
                X=a1[ix]
                Y=a3[iy]
                Xmin=max(X+a3[0]-Y, a1[0])
                Xmax=min(X+a3[-1]-Y, a1[-1])
                idia1=np.where(np.logical_and(a1>=Xmin, a1<=Xmax))[0]
                dia1=a1[idia1]
                dia3=Y-X+dia1
                #tttidia3=a3.searchsorted(dia3)
                idia3=np.interp(dia3, a3, list(range(len(a3)))).round(0).astype(np.int32)
                #print "diagonal index same", (tttidia3==idia3).all()
                #print tttidia3-idia3
                #they differ by one at cca half points
                dia=(dia1-X)*2.0**0.5 #TODO: this multiplication is suspicious
                diaY=[frame[idx,idy] for idx,idy in zip(idia1, idia3)]
                graph.setData(dia, diaY)

                #~ graph2.setVisible(True)
                #~ graph2.setData(a1[idia1], a3[idia3])

                weights=sigmoidWindow(dia, 0, 500., 15)

                fitType=self.cbDisDiagonalFit.currentText()
                if fitType!="Do not fit":
                    graph3.setVisible(True)
                    #fit peak nearest dia==0
                    if fitType=="Gaussian":
                        fit=fitting.fit_gaussian(dia, diaY, weights=weights)
                        self.labDisDiagonalFit.setText("μ: %.0f A: %.0f FWHM: %.0f" % fit[1])
                    elif fitType=="2 Gaussian":
                        fit = fitting.fit_2_gaussian(dia, diaY, weights=weights)
                        self.labDisDiagonalFit.setText("μ: %.0f A: %.0f FWHM: %.0f; μ: %.0f A: %.0f FWHM: %.0f" % fit[1])
                    elif fitType=="Lorentzian":
                        fit=fitting.fit_lorentzian(dia, diaY, weights=weights)
                        self.labDisDiagonalFit.setText("μ: %.0f A: %.0f FWHM: %.0f" % fit[1])
                    graph3.setData(dia, fit[0])
                    pass
                else:
                    graph3.setVisible(False)
                    self.labDisDiagonalFit.setText("")
                pass
            else:
                graph.setVisible(False)
                #~ graph2.setVisible(False)
                graph3.setVisible(False)

            #antidiagonal
            graph=self.plotCuts.plot().plotter(3)
            #~ graph2=self.plotMain.plot().plotter(1)
            graph3=self.plotCuts.plot().plotter(5)
            if self.cbDisAntiDiagonal.isChecked():
                graph.setVisible(True)
                X=a1[ix]
                Y=a3[iy]
                Xmin=max(X-a3[-1]+Y, a1[0])
                Xmax=min(X-a3[0]+Y, a1[-1])
                idia1=np.where(np.logical_and(a1>=Xmin, a1<=Xmax))[0]
                dia1=a1[idia1]
                dia3=Y+X-dia1
                #idia3=a3.searchsorted(dia3)
                idia3=np.interp(dia3, a3, list(range(len(a3)))).round(0).astype(np.int32)
                dia=(dia1-X)*2.0**0.5
                diaY=[frame[idx,idy] for idx,idy in zip(idia1, idia3)]
                graph.setData(dia, diaY)

                weights=sigmoidWindow(dia, 0, 500., 15)

                #~ graph2.setVisible(True)
                #~ graph2.setData(a1[idia1], a3[idia3])

                fitType=self.cbDisAntiDiagonalFit.currentText()
                if fitType!="Do not fit":
                    graph3.setVisible(True)
                    #fit peak nearest dia==0
                    if fitType=="Gaussian":
                        fit=fitting.fit_gaussian(dia, diaY, weights=weights) 
                    elif fitType=="Lorentzian":
                        fit=fitting.fit_lorentzian(dia, diaY, weights=weights) 
                    graph3.setData(dia, fit[0])
                    self.labDisAntiDiagonalFit.setText("μ: %.0f A: %.0f FWHM: %.0f" % fit[1])
                    pass
                else:
                    graph3.setVisible(False)
                    self.labDisAntiDiagonalFit.setText("")
                pass
            else:
                graph.setVisible(False)
                #~ graph2.setVisible(False)
                graph3.setVisible(False)
        except:
            import traceback
            print("tabTemplate.updateDissectionSidePlots: exception")
            print(traceback.print_exc())
            pass
            
        #~ self.plotCuts.rescaleAxes(True)
        self.plotCuts.plot().zoomToData()

        #~ self.plotCuts.replot()
        #~ self.plotMain.replot() #Todo? only call this if necesary
        return ix, iy
    pass
