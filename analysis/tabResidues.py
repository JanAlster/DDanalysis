
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

from .tabTemplate import *
from common import busyCursor


class tabResidues(tabTemplate):
    def __init__(self, main, *args):
        super().__init__(main, __file__.replace("tabResidues.py",'tabResidues.ui'), *args)

        #Dissection
        plot = self.plotTimeline.plot()
        xaxis = plot.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.B, label="t₂ / fs") #btw this is also candidate for NumEnumTransform, to avoid congested  start (but even better it should allow changing Transforms, which is not possible at the moment)
        plot.addSeriesPlotter(xaxis, name="data") #vertical

        self.cbDisDiagonal.toggled.connect(self.pushTimelineD.setEnabled)
        self.cbDisAntiDiagonal.toggled.connect(self.pushTimelineA.setEnabled)
        self._timelineWidget = None
        self.pushTimelineA.clicked.connect(self._timelineA)
        
        #add to context menu
        acShowDisabled = QtWidgets.QAction("Show disabled", self)
        acShowDisabled.setCheckable(True)
        acShowDisabled.triggered.connect(self.updateDissectionSidePlots)
        self.plotTimeline.addAction(acShowDisabled)
        self._acShowDisabled = acShowDisabled
        pass

    def updateInput(self, data=None):
        #~ print("tabResidues.updateInput")
        #this does not depend on UI at all so we might as well calculate it here
        #TODO: or we might feed it residues directly from tabGA
        if data is None:
            data = self.boss.tabGA._dataOutput
            
        model = data
        data = self.boss.tabData._dataShow

        if data["total"].shape != model["total"].shape:
            print("WARNING: tabResidues.updateInput shape mismatch for data and model")
            print(data["total"].shape, model["total"].shape)
            #this can happen if data change (i.e. tabData crop) and model is not updated yet
            #TODO: however that case should be fixed by control of update events
            #  tabRes should get new data only after tabGA has updated model ready
            return 

        res = {"a1":data["a1"], "a2":data["a2"], "a3":data["a3"]}
        for it in ['total', 'rephasing', 'nonrephasing']:
            res[it] = data[it]-model[it]

        self._dataInput = res #optionaly we can store here the model (and modify updateOutput), but it does not matter
        #~ self._dataShow = res
        #~ self.listFrames.setItems(self._dataShow["a2"]) #TODO: this should preserve selection state when just updating GA analysis (but not when loading new data)
        self.updateShow()
        pass
       
    @busyCursor
    def _timelineA(self, *args):
        #display timeline of FWHM fits of antidiagonal at current tracer position for enabled times
        if self._currentData is None:
            return
        
        fitType = self.cbDisAntiDiagonalFit.currentText()
        if fitType == "Do not fit":
            return
             
        tracer = self.plotMain.plot().decorator(0)
        ix = tracer._xL
        iy = tracer._yL
        ix = int(round(ix))
        iy = int(round(iy))
        
        data, a1, a3, a2, dataRange = self._currentData
        if ix<0: ix = 0
        elif ix>=data.shape[1]: ix = data.shape[1]-1
        if iy<0: iy = 0
        elif iy>=data.shape[2]: iy=data.shape[2]-1

        I = self.listFrames.enabled()
        
        t2 = a2[I]
        frames = data[I]
        
        #TODO: check this
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
        
        FWHM = []
        A = []
        S = []
        for frame in frames:
            diaY=[frame[idx,idy] for idx,idy in zip(idia1, idia3)]
        
            weights=sigmoidWindow(dia, 0, 500., 15)
            if fitType=="Gaussian":
                fit=fitting.fit_gaussian(dia, diaY, weights=weights) 
            elif fitType=="Lorentzian":
                fit=fitting.fit_lorentzian(dia, diaY, weights=weights) 
                
            FWHM.append(fit[1][-1])
            A.append(fit[1][1])
            S.append(fit[1][0])
            pass
        
        FWHM = np.array(FWHM)
        
        #now create window to show plot
        if self._timelineWidget is not None:
            #Todo? destroy self._timelineWidget
            pass
            
        w = pqr.FigureWidget(figureTitle="Antidiagonal FWHM timeline")
        p = w.plot()
        xaxis = p.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.B, label="t₂ / fs")
        xaxis._transform.stops = t2
        p.addSeriesPlotter(xaxis, color="red", name="FWHM").setData(t2, FWHM)
        p.addSeriesPlotter(xaxis, color="green", name="Amplitude").setData(t2, A)
        p.addSeriesPlotter(xaxis, color="blue", name="Shift").setData(t2, S)
        p.zoomToData()
        w.resize(400,400)
        w.show()
        self._timelineWidget = w
        pass
    
        
    def updateDissectionSidePlots(self, ix=None, iy=None):
        super().updateDissectionSidePlots(ix, iy)
        #~ print("tabResidues.updateDissectionSidePlots", ix, iy)
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
            self.plotTimeline.plot().x1._transform.stops = a2

        #timeline
        #TODO: decide if to show disabled frames in timeline
        y = data[:, ix, iy]
        if self._acShowDisabled.isChecked():
            self.plotTimeline.plot().plotter(0).setData(a2, y)
        else:
            I = self.listFrames.enabled()
            self.plotTimeline.plot().plotter(0).setData(a2[I], y[I])

        #~ self.plotTimeline.rescaleAxes()
        self.plotTimeline.plot().zoomToData()

        #~ self.plotTimeline.replot()
        pass
    
    pass
