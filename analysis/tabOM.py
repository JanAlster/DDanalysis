
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

#TODO: add maximum line to global OM plot (i.e. not sum of icm frame, but its maximum)

#sublimation
# - inputData (all TRN)
# - showData (interpolated, PAIR)
# - outputData
# tabData     (input: original data, show: original data interpolated, output: original data (selection))
# tabGA       (input: original data, show: DAS interpolated, output: model)
# tabResidues (input: original data, model, show: residues interpolated, output: residues (selection))
# tabFreq     (input: residues, show: OM interpolated, output: )

from .tabTemplate import *
from common import fitting, icmCoef

class tabOM(tabTemplate):
    def __init__(self, main, *args):
        super().__init__(main, __file__.replace("tabOM.py",'tabResidues.ui'), *args)

        self._iOversample = QtWidgets.QSpinBox()
        self._iOversample.setMinimum(1)
        self._iOversample.valueChanged.connect(lambda x: self.updateInput())
        self._iExtend = QtWidgets.QSpinBox()
        self._iExtend.setMinimum(1)
        self._iExtend.valueChanged.connect(lambda x: self.updateInput())
        l = QtWidgets.QHBoxLayout()
        l.addWidget(QtWidgets.QLabel("oversample"))
        l.addWidget(self._iOversample)
        l.addWidget(QtWidgets.QLabel("extend"))
        l.addWidget(self._iExtend)
        self.verticalLayout_4.insertLayout(0, l)

        #Dissection
        plot = self.plotTimeline.plot()
        xaxis = plot.addAxis(None, pqr.AxisPosition.B, label='ω₂ / cm⁻¹') 
        plot.addSeriesPlotter(xaxis, name="data", style=(2,2))
        plot.addSeriesPlotter(xaxis, name="global", color="red", width=1)
        plot.addSeriesPlotter(xaxis, plot.y2, name="global2", color="green", width=1)
        
        pass

    def updateInput(self, data=None):
        print("tabOM.updateInput")
        #this does not depend on UI at all so we might as well calculate it here
        if data is None:
            data = self.boss.tabResidues._dataOutput

        res = data #enabled residues
        a1 = res["a1"]
        a3 = res["a3"]
        t2 = res["a2"]
        # ~ print("t2", t2)
        
        #now calculate frequency analysis of residues under the premise that some steps
        # are not allowed and therefore (actually even if not) time step is not equispaced
        # so FFT is not the correct choice
        # however if timesteps happen to be equispaced, we want results consistent with FFT

        #note that we expect *decaying* oscillatory traces so FT is not the best model anyway
        # (however these can be described by GA and here we need just a quick look)
        
        #in case of FFT
        # dw is given by max(t)
        # max(w) is given by dt
        
        #we can choose mean dt, or minimum dt
        #we can choose max(t) or len(t)*dt
        # probably good to oversample a bit
        dt = np.diff(t2).min() #but not zero
        N = len(t2) * self._iOversample.value()
        if dt==0:
            ttt = np.unique(t2)
            dt = np.diff(ttt).min()
            N = len(ttt) * self._iOversample.value()
        mt = N*dt
        dw = np.pi*2/mt * self._iExtend.value()
        w = np.arange(-(N//2)*dw, (N//2+.5)*dw, dw)
        core = np.exp(1j*np.outer(w, t2)) #Todo: check sign
        # ~ print()
        # ~ print(dt, N, dw, t2)
        # ~ print("tabOM w test", dt, N, dw, w[0], w[-1])
        # ~ print()

        
        OM = {"a1":a1, "a3":a3, "a2":w*icmCoef}
        for it in ['total', 'rephasing', 'nonrephasing']:
            OM[it] = np.tensordot(core, res[it], 1)

        self._dataInput = OM
        #~ self._dataShow = OM
        #~ self.listFrames.setItems(self._dataShow["a2"])
        self.updateShow()
        pass
        
            
    def updateDissectionSidePlots(self, ix=None, iy=None):
        super().updateDissectionSidePlots(ix, iy)
        print("tabOM.updateDissectionSidePlots", ix, iy)
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
            #~ ix,iy = self.plotMain.tracerPosition(True)
            #~ print("ix,iy", ix, iy)
            pass
        ix = int(round(ix))
        iy = int(round(iy))

        data, a1, a3, a2, dataRange, alpha = self._currentData
        if ix<0: ix = 0
        elif ix>=data.shape[1]: ix = data.shape[1]-1
        if iy<0: iy = 0
        elif iy>=data.shape[2]: iy=data.shape[2]-1

        #timeline
        #TODO: decide if to show disabled frames in timeline
        y = data[:, ix, iy]
        self.plotTimeline.plot().plotter(0).setData(a2, y)
        if newdata: 
            ttt = data.sum((1,2))
            self.plotTimeline.plot().plotter(1).setData(a2, ttt/ttt.max()) #plotting this is only needed on _currentData change (probably can be identified by ix and iy is None)
            ttt = (data/data.mean(0)).sum((1,2))
            self.plotTimeline.plot().plotter(2).setData(a2, ttt/ttt.max())
        
        self.plotTimeline.plot().zoomToData()
        #self.plotTimeline.yAxis.setRangeUpper(data.max()) #TODO: we cannot do this because Tracer plugin links axes together (which makes sense for DD, but not here) (this is outdated, this version does not use TRacer plugin, so it might be implemented again)
        #~ self.plotTimeline.replot()
        pass
    
    pass
