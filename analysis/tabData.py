
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

#sublimation
# - inputData (all TRN)
# - showData (interpolated, PAIR)
# - outputData
# tabData     (input: original data, show: original data interpolated, output: original data (selection))
# tabGA       (input: original data, show: DAS interpolated, output: model)
# tabResidues (input: original data, model, show: residues interpolated, output: residues (selection))
# tabFreq     (input: residues, show: OM interpolated, output: )

#TODO: do not rescale Timeline axes on index move (e.g. when we want to look in detail at starting T)
# - only rescale when there are new data with different data range
# - improve Zoom of PluginPlot (see notes in PQCP.PluginPlot)

from .tabTemplate import *
from . import data_pipe

class tabData(tabTemplate):
    def __init__(self, main, *args):
        super().__init__(main, __file__.replace("tabData.py", 'tabData.ui'), *args)
        self._model = None

        #Dissection
        plot = self.plotTimeline.plot()
        xaxis = plot.addAxis(None, pqr.AxisPosition.B, label="t₂ / fs") #btw this is also candidate for NumEnumTransform, to avoid congested  start (but even better it should allow changing Transforms, which is not possible at the moment)
        plot.addSeriesPlotter(xaxis, name="data") #vertical
        plot.addSeriesPlotter(xaxis, plot.y1, name="DAS model", color="red")

        plot = self.plotTimeline2.plot()
        xaxis3 = plot.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.B, label="t₂ / fs", name="blackSheep") #actually it could share the transform with tabResidues
        plot.addSeriesPlotter(xaxis3, name="data") #vertical
        plot.addSeriesPlotter(xaxis3, plot.y1, name="DAS model", color="red")

        plot = self.plotLO.plot()
        plot.addAxis(None, pqr.AxisPosition.B, label='ω₃ / cm⁻¹')
        plot.addAxis(None, pqr.AxisPosition.L)
        
        plot = self.plotLOIntegral.plot()
        xaxis2 = plot.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.B, label="t₂ / fs")
        plot.addSeriesPlotter(xaxis2)
        # plot.y1.format = "b {:.2f}"
        plot.y1.tickFormat = "{:.2f}"

        #
        #~ self.plotTimeline.xAxis.rangeChanged.connect(self.plotLOIntegral.xAxis.setRange)
        xaxis3.zoomRangeLChanged.connect(xaxis2._setZoomRange) #TODO: such connection is not very good as it can lead to (re)cycling, multiple _updateTicks, _setZoomRange, etc. or maybe not
        xaxis2.zoomRangeLChanged.connect(xaxis3._setZoomRange)

        self.cbDisDiagonal.toggled.connect(self.pushTimelineD.setEnabled)
        self.cbDisAntiDiagonal.toggled.connect(self.pushTimelineA.setEnabled)
        
        
        #diagonal plot
        plot = self.plotDiagonal.plot()
        plot.setup(x1={"label":"t₂ / fs", "tickFormat":"{:.0f}", "tickNum":7}, y1={"label":'ω / cm⁻¹', "tickFormat":"{:.0f}", "tickNum":5}) #for now it will be the main diagonal, cannot choose any parallel
        plot.addDDPlotter(plot.x1, plot.y1)

        #add to context menu
        acShowDisabled = QtWidgets.QAction("Show disabled", self)
        acShowDisabled.setCheckable(True)
        acShowDisabled.triggered.connect(self.updateDissectionSidePlots)
        self.plotTimeline.addAction(acShowDisabled)
        self._acShowDisabled = acShowDisabled
        
        #correct LO and crop
        self.addEditable('CropW3Min', self.iCropW3Min.value, self.iCropW3Min.setValue, self.iCropW3Min.valueChanged, group="show", initialize=True)
        self.addEditable('CropW3Max', self.iCropW3Max.value, self.iCropW3Max.setValue, self.iCropW3Max.valueChanged, group="show", initialize=True)
        self.addEditable("CorrectLO", self.bCorrectLO.isChecked, self.bCorrectLO.setChecked, self.bCorrectLO.stateChanged, group="show", initialize=True)
        
        pass

    def _slot(self, *args):
        print("tabData._slot", *args)
        self.plotLOIntegral.xAxis.setRange(*args)
        self.plotLOIntegral.replot()

    def updateInput(self, data=None):
        print("tabData.updateInput")
        #~ if data is None:
            #~ data = self.boss.data._data["data"] #use the original buffer to save memory
        if data is None: return
        self._dataInput = data
        #print("tabData.updateInput: data", data)
        

        #~ will be done by updateShow
        #~ self._dataShow = data
        #~  self.listFrames.setItems(self._dataShow["a2"]) #TODO: this should not trigger changeCurrent
        
        self.updateShow()
        pass

    def updateShow(self):
        #(assumption is that correct UI state is loaded before this is called)
        if not self.isSet("show"): return 

        #if no input is needed use the default implementation
        if self._dataInput is None: return #in case some relevant UI element is set before we have data
        
        self._dataShow = data_pipe.apply_data(self._dataInput, self.bCorrectLO.isChecked(),
                                              self.iCropW3Min.value(), self.iCropW3Max.value())

        #update LO plots
        plot = self.plotLO.plot()
        plot.clearPlotters()
        print("tabData.unpdateInput, create LOs")
        #possibly outside can send data without LO, so load it directly from boss.data
        LO = self._dataShow["LO"]
        if LO is not None:
            for it in LO:
                plot.addSeriesPlotter(plot.x1, plot.y1).setData(self._dataShow["a3"], it)
            
            plot.zoomToData()
            
            LO = LO.sum(1)
            self._LOIntegral = LO/LO[0]*100
        else:
            self._LOIntegral = None


        self.listFrames.setItems(self._dataShow["a2"], force=False) #to not reset the listFrame if we are changing UI params, but 
        
        #data remain the same
        self._updateDataBuffer()
        self.updateOutput()
        pass
     
        
    def updateModel(self, data):
        self._model = data
        self.updateDissectionSidePlots()
        pass
    

    def _updateDataBuffer(self):
        super()._updateDataBuffer()
        
        #plotDiagonal with current interpolation
        data, a1, a3, a2, dataRange, alpha = self._currentData    
        
        wn = max(a1[0], a3[0])
        wx = min(a1[-1], a3[-1])

        idia1 = np.where(np.logical_and(a1>=wn, a1<=wx))[0]
        
        dia = a1[idia1]

        idia3 = np.interp(dia, a3, np.arange(len(a3))).round(0).astype(np.int32)

        #dia=(dia1-X)*2.0**0.5 #TODO: this multiplication is suspicious
        
        diagonal = data[:, idia1, idia3]
        p = self.plotDiagonal.plot()
        p.plotter(0).setData(diagonal.T, a2, dia)
        p.zoomToData()
        
    
    def updateDissectionSidePlots(self, ix=None, iy=None):
        super().updateDissectionSidePlots(ix, iy)
        print("tabData.updateDissectionSidePlots", ix, iy)
        if self._currentData is None:
            return

        newdata = ix is None and iy is None

        if ix is None or iy is None:
            #~ ix,iy = self.plotMain.tracerPosition(True)
            #what we need here is originL value of tracer
            tracer = self.plotMain.plot().decorator(0)
            ix = tracer._xL
            iy = tracer._yL
            print("ix,iy", ix, iy)
            pass
        ix = int(round(ix))
        iy = int(round(iy))
        
        data, a1, a3, a2, dataRange, alpha = self._currentData
        if ix<0: ix = 0
        elif ix>=data.shape[1]: ix = data.shape[1]-1
        if iy<0: iy = 0
        elif iy>=data.shape[2]: iy=data.shape[2]-1

        if newdata:
            self.plotTimeline2.plot().x1._transform.stops = a2
            self.plotLOIntegral.plot().x1._transform.stops = a2

        #timeline
        #TODO: decide if to show disabled frames in timeline
        y = data[:, ix, iy]
        if self._acShowDisabled.isChecked():
            self.plotTimeline.plot().plotter(0).setData(a2, y)
            self.plotTimeline2.plot().plotter(0).setData(a2, y)
            if self._LOIntegral is not None:
                self.plotLOIntegral.plot().plotter(0).setData(a2, self._LOIntegral)
        else:
            I = self.listFrames.enabled()
            self.plotTimeline.plot().plotter(0).setData(a2[I], y[I])
            self.plotTimeline2.plot().plotter(0).setData(a2[I], y[I])
            if self._LOIntegral is not None:
                self.plotLOIntegral.plot().plotter(0).setData(a2[I], self._LOIntegral[I])

        #model should not rescale axes (or at least not much, particularly for negative T)    
        self.plotTimeline.plot().plotter(1).setVisible(False)
        self.plotTimeline2.plot().plotter(1).setVisible(False)
        self.plotTimeline.plot().zoomToData(self.plotTimeline.plot().plotter(0)) #rescaleAxes(True)
        self.plotTimeline2.plot().zoomToData(self.plotTimeline2.plot().plotter(0), verbose=True) #rescaleAxes(True)
        self.plotLOIntegral.plot().zoomToData() #rescaleAxes()

        #model
        if self._model is not None:
            #problem je, ze model neni interpolovany
            trn = TRN2trn[self.cbTRN.currentText()]
            pair = self.cbPAIR.currentText()
            i1 = self.iInterpolation1.value()
            i3 = self.iInterpolation3.value()

            data = self._model[trn]

            if i1>0: ix = int(ix * data.shape[1]/float(i1))
            if i3>0: iy = int(iy * data.shape[2]/float(i3))
            
            if pair in ("PhaseA", "PhaseUA"):
                pair = "Phase"
            self.plotTimeline.plot().plotter(1).setData(self._model["a2"], PAIR[pair](data[:, ix, iy]))
            self.plotTimeline2.plot().plotter(1).setData(self._model["a2"], PAIR[pair](data[:, ix, iy]))
            self.plotTimeline.plot().plotter(1).setVisible(True)
            self.plotTimeline2.plot().plotter(1).setVisible(True)
            pass
        pass
    pass
