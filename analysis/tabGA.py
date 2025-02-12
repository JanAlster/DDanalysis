
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

#TODO: changing IRF ui during fit will trigger DAS update and whole recalculation; this slows down the fit significaly

#sublimation
# - inputData (all TRN)
# - showData (interpolated, PAIR)
# - outputData
# tabData     (input: original data, show: original data interpolated, output: original data (selection))
# tabGA       (input: original data, show: DAS interpolated, output: model)
# tabResidues (input: original data, model, show: residues interpolated, output: residues (selection))
# tabFreq     (input: residues, show: OM interpolated, output: )

#beside tabGA, dataShow need no input from UI
from .tabTemplate import *

class LimitsDialog(QtWidgets.QDialog):
    def __init__(self, *args):
        super().__init__(*args)
        
        # ~ self.setModal(False)
        
        l = QtWidgets.QFormLayout()

        self.lifetimeLow = QtWidgets.QSpinBox()
        self.lifetimeLow.setRange(-1, 10000000)
        self.lifetimeLow.setValue(-1) #-1 means unset
        l.addRow("lifetime Low", self.lifetimeLow)
        
        self.lifetimeHigh = QtWidgets.QSpinBox()
        self.lifetimeHigh.setRange(-1, 10000000)
        self.lifetimeHigh.setValue(-1)
        l.addRow("lifetime High", self.lifetimeHigh)

        #TODO: frequency limits
        #TODO: IRF limits

        #TODO: setup from data

        self.setLayout(l)
        pass
    pass

class tabGA(tabTemplate):
    def __init__(self, main, *args):
        super().__init__(main, __file__.replace("tabGA.py", 'tabGA.ui'), *args)

        plot = self.plotFitness.plot()
        xaxis = plot.addAxis(None, pqr.AxisPosition.B, label='swarm step')
        yaxis = plot.addAxis(None, pqr.AxisPosition.L, label='fitness', tickFormat="{:.3g}")
        plot.addSeriesPlotter(xaxis, yaxis)
        
        self.pushUpdateGA.clicked.connect(self.updateShow)
        self.pushSwarm.clicked.connect(self.swarm)
        self.pushBootstrap.clicked.connect(self.swarmBootstrap)
        
        self.addEditable("GAInput", self.listGAInput.value, self.listGAInput.setValue, self.listGAInput.valueChanged, group="show")
        self.addEditable("IRFOffset", self.dIRFOffset.value, self.dIRFOffset.setValue, self.dIRFOffset.valueChanged, group="show", initialize=True)
        self.addEditable("IRFFWHM", self.dIRFFWHM.value, self.dIRFFWHM.setValue, self.dIRFFWHM.valueChanged, group="show", initialize=True)
        self.addEditable("IRFOffsetFix", self.bIRFOffsetFix.isChecked, self.bIRFOffsetFix.setChecked, self.bIRFOffsetFix.stateChanged, group="show", initialize=True)
        self.addEditable("IRFFWHMFix", self.bIRFFWHMFix.isChecked, self.bIRFFWHMFix.setChecked, self.bIRFFWHMFix.stateChanged, group="show", initialize=True)
        
        
        self._figSSR = pqr.FigureWidget()
        self._figSSR.hide()
        self._figSSR.plot().setup(x1label = "t₂ / fs", y1label="SSR")
        
        self.pushSSR.clicked.connect(self._figSSR.show)
        
        plot = self.plotTimeline.plot()
        plot.setup(x1={"label":"t₂ / fs", "tickFormat":"{:.0f}", "tickNum":7, "transform":pqr.NumEnumTransform()},
                   y1={"label":'2DES / a.u.', "tickFormat":"{:.1f}", "tickNum":3})
        plot.addSeriesPlotter(plot.x1, plot.y1, color="red", name="real")
        plot.addSeriesPlotter(plot.x1, plot.y1, color="blue", name="imag")
        plot.addSeriesPlotter(plot.x1, plot.y1, color="green", name="abs")
        plot.addSeriesPlotter(plot.x1, plot.y1, color="black", name="real orig")
        
        
        self._limitsDia = LimitsDialog(self)
        self.pushLimits.clicked.connect(self._limitsDia.show)
        
        pass

    def updateInput(self, data=None):
        print("tabGA.updateInput")
        if data is None:
            data = self.boss.data._data["data"] #use the original buffer to save memory
        self._dataInput = data
        #update lMaxFreq (notifier label for max freq used in swarm fitting, see setupGA in fitting_pso.py)
        t2 = data["a2"]
        if len(t2)>1:
            dt2 = np.diff(t2).min() #but not zero
            if dt2==0:
                ttt = np.unique(t2)
                if len(ttt)>1:
                    dt2 = np.diff(ttt).min()
                else: #i.e. bleaching test run with just repeating one population time 
                    dt2 = np.NaN
            ub = 2*np.pi/dt2 #if you change this, please adjust lMaxFreq calculation in fitting
            # note that fitting now uses 1.3 * nyquist frequency to avoid aliasing to 0rad/fs (it is annoying to have 2*nyquist freq components masking as nonoscillatory ones)
            from common import icmCoef
            self.lMaxFreq.setText("Max frequency estimate {:.1f}/cm".format(ub*icmCoef))
        else:
            self.lMaxFreq.setText("Max frequency cannot be determined")
        self.updateShow()
        pass
        
    def updateOutput(self):
        print("tabGA.updateOutput", not self.isSet("output"), self._dataShow is None)
        if not self.isSet("output"): return
        if self._dataShow is None: return

        #use DAS (calculated using only enabled frames) to get model on all t2
        enabled = self.enabled() #which DAS use to build the model
        if len(enabled)==0: return
        #inputList = [self._dataShow["input"][it] for it in enabled]
        inputList = self._dataShow["input"]
        params, rates, w, rFix, wFix, types, paramsTypes = fitting.divide(inputList)
        #~ print(inputList) #[(2.0, False, 'Ø', 0.0, False), (150.57, False, 'Ø', 0.0, False), (8434.35, False, 'Ø', 0.0, False), (710.16, False, 'Ø', 0.0, False), (1000000.0, False, 'Ø', 0.0, False), (500.0, False, '±', 750.0, False)]
        #~ print(rates) # [0.5, 0.006641429235571495, 0.00011856278195711583, 0.0014081333783936016, 1e-06, 0.002]
        #~ print(types) # ['Ø', 'Ø', 'Ø', 'Ø', 'Ø', '±']
        #~ print(w) # [0.0, 0.0, 0.0, 0.0, 0.0, 0.14127365172074868]
        #~ print(enabled) # [0, 1, 2, 3, 4, 5, 6]
        #~ print(self._dataShow["a2"]) #[2.0, 150.57, 8434.35, 710.16, 1000000.0, (-500+750j), (-500-750j)]
        
        #~ #TODO: selection is nontrivial here, number of elements in each list is different (due to +- type)
        #~ # we probably need to build rates, types and w for fitting.C from scratch
        #~ # because +- type will not be selected as one
        #~ # and it does not matter for the model if we separate them to + and -, their lifetimes will be fixed
        #~ # w.append(freq/icmCoef)
        from common import icmCoef
        a2 = self._dataShow["a2"]
        rates = [1./a2[i].real for i in enabled]
        types = ['Ø' if a2[i].imag==0 else '+' if a2[i].imag>0 else '-' for i in enabled]
        w = [abs(a2[i].imag)/icmCoef for i in enabled]
        
        # ~ print()
        # ~ print(inputList)
        # ~ print(a2)
        # ~ print(rates)
        # ~ print(types)
        # ~ print(w)
        # ~ print()
        
        #~ print()
        #~ print(rates)
        #~ print(types)
        #~ print(w)
        #~ print()
        
        #~ rates = [rates[it] for it in enabled]
        #~ types = [types[it] for it in enabled]

        #TODO: improve model calculation at negative T (it should not be exp, but exp * H, or exp * H convolved with IRF)
        t2 = self.boss.tabData._dataInput["a2"]
        # ~ c = fitting.C(rates, w, types, t2) #kinetic profiles on whole t2
        
        c = fitting.CIRF(rates, w, types, t2, self.dIRFOffset.value(), self.dIRFFWHM.value())

        #for testing
        plot = self._figSSR.plot()
        plot.clearPlotters()
        plot.setup(legend=True)
        for C, r, t, w in zip(c.T, rates, types, w):
            plot.addSeriesPlotter(plot.x1, plot.y1, color=pqr.color(), name=str(r)+" "+t+" "+str(w)).setData(t2, C.real)
        plot.zoomToData()

        
        M = {"a1":self._dataShow["a1"], "a3":self._dataShow["a3"], "a2":t2}
        for it in ['total', 'rephasing', 'nonrephasing']:
            das = self._dataShow[it][enabled]
            M[it] = np.tensordot(c, das,1)
        
        self._dataOutput = M
        #notify model changed
        print("tabGA.updateOutput: emit outputChanged")
        self.outputChanged.emit(self._dataOutput)
        pass

    def updateShow(self):
        print("tabGA.updateShow", not self.isSet("show"), self._dataInput is None)
        if not self.isSet("show"):
            return #unset this flag when loading new UI state
        if self._dataInput is None: return 
        #calculate DAS after change of GAInput
        #or when we have new data to fit
        inputList = self.listGAInput.value()
        print("tabGA.updateGA", inputList)
        if len(inputList)==0:
            self.listFrames.setItems([])
            return

        DAS_t2=[]
        for lifetime, lFix, ft, freq, fFix in inputList:
            #this has to copy fitting.c
            if ft=='Ø':   
                DAS_t2.append(lifetime) #technically it should be -lifetime
            #TODO: check signs here
            assert(freq>=0.)
            if ft=='+' or ft=='±':
                DAS_t2.append(lifetime+freq*1j) #technically it should be -lifetime
            if ft=='-' or ft=='±':
                DAS_t2.append(lifetime-freq*1j) #technically it should be -lifetime       
            pass
        #here t2 is less obvious, it should be taken from input list, but might show frequencies as well

        #for DAS calculation use uninterpolated data, only enabled frames (but we calculate model on whole a2)

        DAS={'a2':DAS_t2, "a1":self._dataInput["a1"], "a3":self._dataInput["a3"], "input":inputList[:], "disabled":self.boss.tabData.disabled()} 

        for it in ['rephasing', 'nonrephasing']:
            DAS[it] = fitting.DASIRF(self._dataInput[it], self._dataInput["a2"], inputList, self.dIRFOffset.value(), self.dIRFFWHM.value())

        DAS['total'], R = fitting.DASandResIRF(self._dataInput["total"], self._dataInput["a2"], inputList, self.dIRFOffset.value(), self.dIRFFWHM.value())

        """
        #for testing
        with open("testSSR.pickle", "wb") as f:
            import pickle
            pickle.dump((DAS["total"], R, self._dataInput["total"], self._dataInput["a2"], inputList, self.dIRFOffset.value(), self.dIRFFWHM.value() ), f)
        """
        
        #here we can output sum of squared residues
        #see common/fitting for the metric used
        R = np.abs(R)
        res = np.sum(R).real
        self.lSSR.setText(f"Sum of abs(residues): {res:.4e}")
        
        R2 = np.sum(R.reshape((len(R), -1)), 1)
        plot = self._figSSR.plot()
        plot.clearPlotters()
        plot.addSeriesPlotter(plot.x1, plot.y1).setData(self._dataInput["a2"], R2)
        plot.zoomToData()


        self._dataShow = DAS
        self.listFrames.setItems(self._dataShow["a2"])

        #update listFrames
        self._updateDataBuffer() 
        self.updateOutput()
        pass        

    def updateDissectionSidePlots(self, ix=None, iy=None):
        super().updateDissectionSidePlots(ix, iy)
        """
        add plotting the timeline for the selected DAS component
        """
        
        #TODO: this is not very efficient as we repeat a lot what is done in super().updateDissctionSidePlots() :(
        
        print("tabGA.updateDissectionSidePlots", ix, iy)
        if self._currentData is None:
            return

        indices = self.listFrames.selectedItems()
        if len(indices) <= 1: #TODO: there is a problem, we react to currentChanged, but selection is not updated at that time; it is especially bad if we select just the current
            index = self.listFrames.current()
            if index is None: return
            indices = [index]

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
        frame = data[indices]
        if ix<0: ix = 0
        elif ix>=data.shape[1]: ix = data.shape[1]-1
        if iy<0: iy = 0
        elif iy>=data.shape[2]: iy=data.shape[2]-1

        #we want to show frame[ix, iy]* np.exp(a2[index]*t2)
        
        #get t2 from tabData - we likely want ALL data points even for disabled data frames
        t2 = self.boss.tabData._dataInput["a2"]
        t2 = self._dataInput["a2"] #selected frames
        #but make fine step otherwise we will get glitches in time
        dt2 = np.diff(t2).mean()
        t2fine = np.arange(t2[0], t2[-1], dt2)
        
        #convert DAS component to correct units
        from common import icmCoef
        a2sel = np.asarray(a2)[indices]
        rates = 1./a2sel.real
        types = ['Ø' if it==0 else '+' if it>0 else '-' for it in a2sel.imag]
        w = abs(a2sel.imag)/icmCoef

        #calculate model
        cfine = fitting.CIRF(rates, w, types, t2fine, self.dIRFOffset.value(), self.dIRFFWHM.value())
        c = fitting.CIRF(rates, w, types, t2, self.dIRFOffset.value(), self.dIRFFWHM.value())

        trace_fine = np.dot(cfine, frame[:, ix, iy])
        trace = np.dot(c, frame[:, ix, iy])

        pl = self.plotTimeline.plot()
        pl.x1._transform.stops = t2fine
        pl.plotter(0).setData(t2fine, trace_fine.real)
        pl.plotter(1).setData(t2fine, trace_fine.imag)
        pl.plotter(2).setData(t2fine, np.abs(trace_fine))
        pl.plotter(3).setData(t2, trace.real)
        pl.zoomToData()

        pass

    def swarmBootstrap(self):
        #TODO: possibly allow fitting all TRN at the same time
        trn = TRN2trn[self.cbTRN_GA.currentText()]        
        data = self._dataInput[trn]
        t2 = self._dataInput["a2"]
        w1 = self._dataInput["a1"]
        w3 = self._dataInput["a3"]
        LO = self._dataInput["LO"]
        
        #TODO: this will not save w1 and w3 axes with bootstrap
        
        self.pushBootstrap.setEnabled(False) #prevent running two swarms at the same time
        params = self.listGAInput.value()
        
        IRF = self.dIRFFWHM.value(), self.bIRFFWHMFix.isChecked(), self.dIRFOffset.value(), self.bIRFOffsetFix.isChecked()
        
        #for guessing rate limit (to prevent single frame components) we need max amplitude of the signal overall
        # regardless of what is selected for fitting
        # actually it should be max at t2=0fs, but lets assume that it is the max overall
        M = np.abs(self.boss.tabData._dataShow["total"]).max()
        
        #filepath to store data
        startFilename = self.boss.data._filename
        filename = QtWidgets.QFileDialog.getSaveFileName(filter="Bootstrap files (*.bootstrap)", directory=startFilename[:startFilename.rfind('.')])[0]
        if filename == "":
            return

        if not filename.endswith(".bootstrap"):
            filename += ".bootstrap"

        metadata = {"datafile": startFilename, "state": self.boss.data.lastState}
        
        #separate thread approach
        swarm = fitting.SwarmThreadIRFBootstrap(self, filename, data, t2, w1, w3, LO, params, IRF,
                                                metadata=metadata, maxAmplitude=M)
        print(params)
        print(swarm)

        # ~ self.pushStopSwarm.clicked.connect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        #Todo: interrup swarm if new data loaded or exit
        #self.pushLoad.clicked.connect(swarm.requestInterruption)#if new data are loaded, stop this
        #self.pushExit.clicked.connect(swarm.requestInterruption)#if program terminates, stop this #TODO: program can wait for this to terminate...
        # ~ swarm.report.connect(self.updateGAfromSwarm)#reports from swarm
        # ~ swarm.finalReport.connect(self.updateGAfromSwarmFinal)#reports from swarm
        swarm.finalReport.connect(lambda: self.pushBootstrap.setEnabled(True))#reports from swarm
        swarm.start()

        #response slot should do this and update UI
        #fitting.combine(result, inputList)
        #print inputList
        pass


    def swarm(self):
        #TODO: possibly allow fitting all TRN at the same time
        trn = TRN2trn[self.cbTRN_GA.currentText()]        
        data = self._dataInput[trn]
        t2 = self._dataInput["a2"]
        
        self.pushSwarm.setEnabled(False) #prevent running two swarms at the same time
        params = self.listGAInput.value()
        
        IRF = self.dIRFFWHM.value(), self.bIRFFWHMFix.isChecked(), self.dIRFOffset.value(), self.bIRFOffsetFix.isChecked()
        
        #for guessing rate limit (to prevent single frame components) we need max amplitude of the signal overall
        # regardless of what is selected for fitting
        # actually it should be max at t2=0fs, but lets assume that it is the max overall
        M = np.abs(self.boss.tabData._dataShow["total"]).max()
        
        #separate thread approach
        swarm=fitting.SwarmThreadIRF(self, data, t2, params, IRF, 1, maxAmplitude=M)
        print(params)
        print(swarm)

        self.plotFitness.plot().plotter(0).setData([],[])
        
        #for testing
        self._prevParams = None
        self._bestParams = None
        self._bestFitness = None
        
        self.pushStopSwarm.clicked.connect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        #Todo: interrup swarm if new data loaded or exit
        #self.pushLoad.clicked.connect(swarm.requestInterruption)#if new data are loaded, stop this
        #self.pushExit.clicked.connect(swarm.requestInterruption)#if program terminates, stop this #TODO: program can wait for this to terminate...
        swarm.report.connect(self.updateGAfromSwarm)#reports from swarm
        swarm.finalReport.connect(self.updateGAfromSwarmFinal)#reports from swarm
        swarm.start()

        #response slot should do this and update UI
        #fitting.combine(result, inputList)
        #print inputList
        pass

    def updateGAfromSwarm(self, params, step, fitness):
        if self._prevParams is None or self._prevParams!=params:
            self._prevParams = params
            inputList = self.listGAInput.value()
            inputIRF = [self.dIRFFWHM.value(), self.bIRFFWHMFix.isChecked(), self.dIRFOffset.value(), self.bIRFOffsetFix.isChecked()]
            print("pre", inputIRF)
            fitting.combineIRF(params, inputList, inputIRF)
            print("after", inputIRF)
            self.listGAInput.setValue(inputList, emit=False)
            print("change GAInput", inputList)
            
            
            _,  R = fitting.DASandResIRF(self._dataInput["total"], self._dataInput["a2"], inputList, self.dIRFOffset.value(), self.dIRFFWHM.value())
            R = np.abs(R)
            res = np.sum(R).real
            print(inputList, res)

            _,  R = fitting.DASandResIRF(self._dataInput["total"], self._dataInput["a2"], self.listGAInput.value(), self.dIRFOffset.value(), self.dIRFFWHM.value())
            R = np.abs(R)
            res = np.sum(R).real
            print(self.listGAInput.value(), res)


            #TODO: this must not emit (not possible with QDoubleSpinBox)
            # or do not listen (is editable in group show - so we need to disconenct it from the group or block the group)
            self.disableGroup("show")
            if not inputIRF[1]:
                self.dIRFFWHM.setValue(inputIRF[0])
            if not inputIRF[3]:
                self.dIRFOffset.setValue(inputIRF[2])
            self.enableGroup("show")
            
        if self._bestFitness is None or fitness<self._bestFitness:
            self._bestFitness = fitness
            self._bestParams = params
            
        #plot fitness
        self.plotFitness.plot().plotter(0).removeDataBefore(step-40)#TODO: update this 40 depending on reportEachNth parameter of swarm
        self.plotFitness.plot().plotter(0).addData([step], [fitness])
        self.plotFitness.plot().zoomToData() #rescaleAxes()
        #plot residues
        #only with final report
        pass

    def updateGAfromSwarmFinal(self, params, step, fitness):
        print("updateGAfromSwarmFinal: fitness", fitness, "params", params)
        print("\best", self._bestFitness, self._bestParams)
        if params == []:
            #this indicated an error
            self.pushSwarm.setEnabled(True)
            return
        
        self.updateGAfromSwarm(params, step, fitness)
        self.pushSwarm.setEnabled(True)
        self.updateShow()
        
        #TEST
        #there is a discrepancy between fitness reported and calculated sum of residues
        params = self.listGAInput.value()

        _,  R = fitting.DASandResIRF(self._dataInput["total"], self._dataInput["a2"], params, self.dIRFOffset.value(), self.dIRFFWHM.value())

        #here we can output sum of squared residues
        #see common/fitting for the metric used
        R = np.abs(R)
        res = np.sum(R).real
        print("\tsummed", res)
        pass

    pass
