
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

"""
This is column-wise global analysis based on assumpiton that peaks with same lifetimes should appear
at the same excitation frequency due to energy transfer (not all of them of course, but energy transfer is what
we are looking for).

WE do not need to change much from tabLA, only the model
Perhaps we can even inherit tabLA instead of copy/paste?

"""
import numpy as np

# copied from tabLA
#TODO: global colors should be pseudo frame specific - because we have different units for different pseudo frames
#TODO: persistent saving during fitting - could be directly to self._LA - although that is not persistent, we would need to save "LA" state with every iteration (or once a minute) to be persistent - we would need to be able to save state without changing current state
#TODO: lifetime selection has to allow for infinity, both in pseudo frame selection and in data range for colors
# perhaps just put one number and consider it division/border between pseudo frames (with the last one being larger than the boundary)
# but that will not work for frequency selections
# possibly also include "goodness of fit" pseudoframe
#TODO: change color of nan in QImage from black to some grey
#TODO: we can change pseudo frames to use 2d color map -> (saturated) color for lifetime and alpha (lightness) for amplitude
#    - but that would require redesign of colormaps
#    - on the other hand we could use that for phase maps for OM
# at very least we need amplitude (or raw data amplitude) contours over lifetime pseudoframes

# ~ TODO iam not so ssure that removing current seudo frame will not lead to roblems - buttongrou
 # ~ might not be uadate and current might oiny to bad index in calculated seudoframes

#TODO: maybe we can do LA as a correction to GA - that way we have give number of components, but still LA can suppress some by amplitude 
 
#TODO: muzeme zkusit semi-GA - drzet stejne parametry pro cely radek nebo sloupec ale menit je mezi sloupci - to by mohlo potlacit sum ale dovolit sledovat zmeny s excitaci 

# TODO: show amplitudes for single point kinetics
# TODO: somehow limit crazy values like amplitudes of overfitted single points, etc. (select point with lifetime in range and amplitude in range)


# from CA
# TODO: we do not save to state, so we cannot have multiple different CAs
# TODO: selected regions/selection rules are not saved to file
# TODO: add button to clear the CA buffer in case we want to start from nothing

#beside tabGA, dataShow need no input from UI
from .tabTemplate import *

#we need something to select what to show in listFrames
from PyQt5 import QtWidgets, QtCore

from common.widgets.GrowingList import GrowingList
from .tabLA import LAInputList, LAItemWidget


    
class tabCA(tabTemplate):
    #this needs to replace listFrames with something different, because LA does not give "frames" - we got bunch of lifetimes and we have to select an interval (preferably with only one value at each position
    #only calculate LA at specified region - easiest implementation is between two tracers
    # but do not forget what is already fitted if region is changed
    # allow fitting single point to fine tune
    #this will not play nicely with interpolation
    #   - changing interpolation will discard all already fitted LA
    # - we could interpolated only the results, but I am not sure how will that work with LA...
    # - if we interpolate we should interpolate source data, but that will blow up LA work
    
    #can we highjack editable I1 and I3 (group buffer) to warn user before changing interpolation?
    # or we can store LA for different interpolations
    def __init__(self, main, *args):
        super().__init__(main, __file__.replace("tabCA.py", 'tabLA.ui'), *args)

        plot = self.plotTimeline.plot()
        xaxis = plot.addAxis(None, pqr.AxisPosition.B, label="t₂ / fs")
        yaxis = plot.addAxis(None, pqr.AxisPosition.L, label='2DES / a.u.')
        plot.addSeriesPlotter(xaxis, yaxis, name="data") 
        plot.addSeriesPlotter(xaxis, yaxis, name="DAS model", color="red")
        
        # ~ self.addEditable("LAInput", self.listGAInput.value, self.listGAInput.setValue, self.listGAInput.valueChanged, group="show")

        #tabTemplate also connects currentChanged to updatePlotMain
        self.listFrames.currentChangedRaw.connect(self._updateShowRaw)
        self.listFrames.valueChanged.connect(self.updateShow)
        
        self.listLAInputSingle.valueChanged.connect(self.updateDissectionSidePlots2) #this should update the model
        """
        actually no
        - when tracer moves we should update listLAInputSingle to what is in self._LA
        - in updatedissectionPlots we should do the model by updateDissectionPlots
        - but change in listLAInput (by user in UI) shoudl trigger updateDissectionPlots
        -> we cannot use updateDissectionPlots to set listLAInputSingle to _LA
        - trace r move can trigger updateDissectionPlots two times - directly from tabTemplate and after listLAInputSingle changes value (but if the fitted params on new point are the same as on old point, change will not be triggered)
        
        we can bypass tabTemplate call of updateDissectionPlots to trigger move and alway call updateDissectionPlots via setting listLAInputSingle.setValue(..., emit=True) - emit of valueChanged is forced by emit=True
        (of course updateDissectionPlots will call super().updateDissectionPlots)
        
        or we can highjack updateDissectionPlots
        """
        
        
        
        """
        we have two inputs 
        - "global" listLAInput which should be used as starting point when whole region is fitted for the first time
        - single point listLAInputSingle which will show single point best fit params with cursor move and can be used as a starting point for single kinetic fitting (or simply to set the params for this point)
        (note that there is no need to set whole region to something, because that is not LA)
        
        we have four "fit" buttons
        - pushSwarm - this should start point-by-point fitting of selected region using listLAInput as starting point
        - pushSwarmEmpty - point-by-point fitting in selected region using listLAInput as starting point and skipping already fitted points
        - pushSwarmContinue - this should start point-by-point fitting of selected region using per point best fit params as starting point (continue fitting where it left off)
        - pushSwarmSingle - this should start single point kinetic fitting using listLAInputSingle as starting point (which should be the best previous fit params if jsut moved )
        
        """
        self.pushSwarm.clicked.connect(self.swarmFromStart)
        self.pushSwarmContinue.clicked.connect(self.swarmContinue)
        self.pushSwarmSingle.clicked.connect(self.swarmSingle)
        self.pushSwarmEmpty.clicked.connect(self.swarmEmpty)
        
        #add two more tracers to main plot to select region for fitting LA
        plot = self.plotMain.plot()
        self._tracer1 = plot.addDecorator(pqr.TracerDecorator(plot, plot.b1, plot.l1, color="orange"))
        self._tracer2 = plot.addDecorator(pqr.TracerDecorator(plot, plot.b1, plot.l1, color="orange"))
        # ~ self._tracer1.setPos(0,0)
        self._tracer1.setPos(0.4,0.4)
        # ~ self._tracer2.setPos(1,1)
        self._tracer2.setPos(0.6,0.6)
        timeaxis = plot.addAxis(None, pqr.AxisPosition.R, label='lifetime / fs', tickFormat="{:.0f}")
        self._column_lifetime_series = plot.addSeriesPlotter(plot.b1, timeaxis, color="lime")
        pass

    def saveSettings(self, saveData=False):
        #TODO: we also need to saveload LAinput lists, they are not automated in init so far
        settings = super().saveSettings() #those setup SaveLoadElement(s) above
        settings["CA"] = self._CA
        settings["A"] = self._A
        return settings

    #saved data will not be loaded, they are used for something else
    def loadSettings(self, settings):
        if "CA" in settings: self._CA = settings["CA"]
        if "A" in settings: self._A = settings["A"]
        super().loadSettings(settings) #those setup SaveLoadElement(s) above
         
    def resetData(self):
        self._CA = None
        self._A = None 
        self._lastPoint = None
        super().resetData()

    def swarmFromStart(self, *args):
        self.pushSwarm.clicked.disconnect(self.swarmFromStart)
        self.pushSwarmEmpty.setEnabled(False)
        self.pushSwarmSingle.setEnabled(False)
        self.pushSwarmContinue.setEnabled(False)

        trn = TRN2trn[self.cbTRN_GA.currentText()]        
        data = self._dataInput[trn]
        t2 = self._dataInput["a2"]
        
        params = self.listLAInput.value()
        
        #limit to region
        x1, y1 = self._tracer1.origin
        x2, y2 = self._tracer2.origin
        i1n, i1x = self._dataInput["a1"].searchsorted((min(x1, x2), max(x1, x2)))
        i3n, i3x = self._dataInput["a3"].searchsorted((min(y1, y2), max(y1, y2)))
        #and store this in case user changes the region during the fit
        self._region = i1n, i1x, i3n, i3x
        
        #separate thread approach
        swarm = fitting.CASwarmThread(self, data[:,i1n:i1x,i3n:i3x], t2, params, 1000)

        # todo: should not we disconnect these signals later on?
        self.pushSwarm.setText("Stop")
        self.pushSwarm.clicked.connect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        #Todo: interrup swarm if new data loaded or exit
        #TODO: program can wait for this to terminate...
        swarm.report.connect(self.updateLAFromStartProgress)
        swarm.finalReport.connect(self.updateLAFromStart)#reports from swarm
        swarm.start()
    
    def updateLAFromStartProgress(self, p):
        self.pushSwarm.setText("Stop with {:.0f}% done".format(p))
        
    def updateLAFromStart(self, paramsCA, A):
        print("updateLAFromStart")

        #TODO: this is probably not disconnected? self.pushSwarm.clicked.disconnect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        self.pushSwarm.setText("Fit Region from Start")
        self.pushSwarm.clicked.connect(self.swarmFromStart)
        
        i1n, i1x, i3n, i3x = self._region 
        self._CA[i1n:i1x] = paramsCA #only update fitted region
        self._A[i1n:i1x, i3n:i3x] = A #only update fitted region

        self.updateShow() #recalculate data to show

        self.pushSwarmSingle.setEnabled(True)
        self.pushSwarmContinue.setEnabled(True)
        self.pushSwarmEmpty.setEnabled(True)
        pass

    def swarmContinue(self, *args):
        self.pushSwarmContinue.clicked.disconnect(self.swarmContinue)
        self.pushSwarmEmpty.setEnabled(False)
        self.pushSwarmSingle.setEnabled(False)
        self.pushSwarm.setEnabled(False)

        trn = TRN2trn[self.cbTRN_GA.currentText()]        
        data = self._dataInput[trn]
        t2 = self._dataInput["a2"]
        
        params = self.listLAInput.value()
        
        #limit to region
        x1, y1 = self._tracer1.origin
        x2, y2 = self._tracer2.origin
        i1n, i1x = self._dataInput["a1"].searchsorted((min(x1, x2), max(x1, x2)))
        i3n, i3x = self._dataInput["a3"].searchsorted((min(y1, y2), max(y1, y2)))
        #and store this in case user changes the region during the fit
        self._region = i1n, i1x, i3n, i3x
        
        #separate thread approach
        swarm = fitting.CASwarmThread(self, data[:,i1n:i1x,i3n:i3x], t2, params, 1000, singleInputList=self._CA[i1n:i1x])

        self.pushSwarmContinue.setText("Stop")
        self.pushSwarmContinue.clicked.connect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        #Todo: interrup swarm if new data loaded or exit
        #TODO: program can wait for this to terminate...
        swarm.report.connect(self.updateLAContinueProgress)
        swarm.finalReport.connect(self.updateLAContinue)#reports from swarm
        swarm.start()
    
    def updateLAContinueProgress(self, p):
        self.pushSwarmContinue.setText("Stop with {:.0f}% done".format(p))
        
    def updateLAContinue(self, paramsCA, A):
        print("updateLAContinue")

        #TODO: this is probably not disconnected? self.pushSwarm.clicked.disconnect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        self.pushSwarmContinue.setText("Fit Region More")
        self.pushSwarmContinue.clicked.connect(self.swarmContinue)
        
        i1n, i1x, i3n, i3x = self._region 
        self._CA[i1n:i1x] = paramsCA #only update fitted region
        self._A[i1n:i1x, i3n:i3x] = A #only update fitted region

        self.updateShow() #recalculate data to show

        self.pushSwarmSingle.setEnabled(True)
        self.pushSwarm.setEnabled(True)
        self.pushSwarmEmpty.setEnabled(True)
        pass
        
    def swarmEmpty(self, *args):
        self.pushSwarmEmpty.clicked.disconnect(self.swarmEmpty)
        self.pushSwarm.setEnabled(False)
        self.pushSwarmSingle.setEnabled(False)
        self.pushSwarmContinue.setEnabled(False)

        trn = TRN2trn[self.cbTRN_GA.currentText()]        
        data = self._dataInput[trn]
        t2 = self._dataInput["a2"]
        
        params = self.listLAInput.value()
        
        #limit to region
        x1, y1 = self._tracer1.origin
        x2, y2 = self._tracer2.origin
        i1n, i1x = self._dataInput["a1"].searchsorted((min(x1, x2), max(x1, x2)))
        i3n, i3x = self._dataInput["a3"].searchsorted((min(y1, y2), max(y1, y2)))
        #and store this in case user changes the region during the fit
        self._region = i1n, i1x, i3n, i3x
        
        #separate thread approach
        swarm = fitting.CASwarmThread(self, data[:,i1n:i1x,i3n:i3x], t2, params, 1000, skipFilled=self._CA[i1n:i1x])

        self.pushSwarmEmpty.setText("Stop")
        self.pushSwarmEmpty.clicked.connect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        #Todo: interrup swarm if new data loaded or exit
        #TODO: program can wait for this to terminate...
        swarm.report.connect(self.updateLAEmptyProgress)
        swarm.finalReport.connect(self.updateLAEmpty)#reports from swarm
        swarm.start()
    
    def updateLAEmptyProgress(self, p):
        self.pushSwarmEmpty.setText("Stop with {:.0f}% done".format(p))
        
    def updateLAEmpty(self, paramsCA, A):
        print("updateLAEmpty")

        #TODO: this is probably not disconnected? self.pushSwarm.clicked.disconnect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        self.pushSwarmEmpty.setText("Fit Empty from Start")
        self.pushSwarmEmpty.clicked.connect(self.swarmEmpty)
        
        i1n, i1x, i3n, i3x = self._region 
        self._CA[i1n:i1x] = paramsCA #only update fitted region
        self._A[i1n:i1x, i3n:i3x] = A #only update fitted region

        self.updateShow() #recalculate data to show

        self.pushSwarmSingle.setEnabled(True)
        self.pushSwarmContinue.setEnabled(True)
        self.pushSwarm.setEnabled(True)
        pass
        
    def swarmSingle(self, *args):
        if self._lastPoint is None:
            return 

        self.pushSwarmSingle.clicked.disconnect(self.swarmSingle)
        self.pushSwarmEmpty.setEnabled(False)
        self.pushSwarm.setEnabled(False)
        self.pushSwarmContinue.setEnabled(False)

        trn = TRN2trn[self.cbTRN_GA.currentText()]        
        data = self._dataInput[trn]
        t2 = self._dataInput["a2"]
        
        params = self.listLAInputSingle.value()
        
        #limit to point
        ix, iy = self._lastPoint
        self._region = ix, iy
        
        #separate thread approach
        swarm = fitting.CASwarmThread(self, data[:,ix:ix+1,iy:iy+1], t2, params, 1000)

        self.pushSwarmSingle.setText("Stop")
        self.pushSwarmSingle.clicked.connect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        #Todo: interrup swarm if new data loaded or exit
        #TODO: program can wait for this to terminate...
        swarm.report.connect(self.updateLASingleProgress)
        swarm.finalReport.connect(self.updateLASingle)#reports from swarm
        swarm.start()
    
    def updateLASingleProgress(self, p):
        self.pushSwarmSingle.setText("Stop with {:.0f}% done".format(p))
        
    def updateLASingle(self, paramsCA, A):
        print("updateLASingle")

        #TODO: this is probably not disconnected? self.pushSwarm.clicked.disconnect(swarm.requestInterruption)#stop when self.pushStopSwarm is pressed
        self.pushSwarmSingle.setText("Fit Single Point")
        self.pushSwarmSingle.clicked.connect(self.swarmSingle)
        
        ix, iy = self._region
        self._CA[ix:ix+1] = paramsCA #only update fitted region
        self._A[ix:ix+1, iy:iy+1] = A #only update fitted region

        self.updateShow() #recalculate data to show

        self.pushSwarm.setEnabled(True)
        self.pushSwarmContinue.setEnabled(True)
        self.pushSwarmEmpty.setEnabled(True)
        pass

    def updateInput(self, data=None):
        print("tabCA.updateInput")
        #TODO: is this a good place to resetData?
        self.resetData() #if we do not reset data somewhere, _updataShowRaw fill throw an exception for incompatible data
        
        if data is None:
            data = self.boss.data._data["data"] #use the original buffer to save memory
        self._dataInput = data
        #this will destroy all LA, and that could be annoying - perhaps we should save special "LA backup" state after every LA fit
        #LA storage is same grid as data, with arbitrary length list at each position (listLAInput params basically)
        #   - should it have amplitudes too? or can it be calculated similarly to DAS?
        if data is None:
            self._CA = None
            self._A = None
            self.listFrames.setItems([], True)
        else:
            self._CA = np.asarray([None for it in data["a1"]])
            self._A = np.asarray([[None for it in data["a3"]] for it in data["a1"]])
            self.listFrames.setItems(data["a2"], True)
        #Todo: keep tracer positions here (relative to frame)? or absolute to data coordinates?
        #  it does make sense to keep data coordinates if new data have similar range (i.e. tracers are in new data range)
        #  otherwise it would be better to reset to some starting position
        self.updateShow()
        self.listFrames.setCurrent(0)
        pass
        
    def updateOutput(self):
        #no output
        pass

    def _updateShowRaw(self, index):
        #this means that the raw data pseudo frame (which should be the first frame of self._currentData)
        # is changed (new data or user selected different raw frame to show)
        #_currentData should be updated (but no need to calculate other pseudo frames
        print("tabCA._updateShowRaw", index)
        if self._dataInput is None: return 
        if self._dataShow is None: return
        raw = self.listFrames.rawIndex()
        if raw<0: return
        
        for it in ["total", "rephasing", "nonrephasing"]:
            self._dataShow[it][0] = self._dataInput[it][raw]
        self._updateDataBuffer() 
        self.updateOutput()
        

    def updateShow(self):
        print("tabCA.updateShow", not self.isSet("show"), self._dataInput is None)
        if not self.isSet("show"): return #unset this flag when loading new UI state        
        if self._dataInput is None: return 
        
        #here we do not use _dataInput (raw data)
        # but fitted single point kinetics (when and where available)
        
        #go through listFrames (LA version) and calculate "frames"
        # - first frame is slected rawData frame for easy region selection
        # - each of other componets have two frames - amplitudes of selected components (ie. local DAS) and their lifetimes (which are constant for DAS but spectraly dependend here)
        
        #output should be {"total":[frame0, frame1, ...], "rephasing":[f0, f1, ...], "nonrephasing":[f0,f1], "a1":a1, "a3":a3}
        raw = self.listFrames.rawIndex()
        if raw<0: return
        
        C = self.listFrames.count()
        
        TT = [self._dataInput["total"][raw]]
        RR = [self._dataInput["rephasing"][raw]]
        NN = [self._dataInput["nonrephasing"][raw]]
        
        sh = TT[0].shape
        
        P = self.listFrames.value()
        
        if len(P) != (C-1)//3:
            print("WARNING: tabCA.updateShow - unexpected number of pseudo-frames", len(P), C, (C-1)//3)
        
        for K in range((C-1)//3):
            #now we have to keep order of filling T,R,N to what listFrame(LAInputList) expects, otherwise there will be problems
            #in other tabs, the tab is the one to set items in listFrames, but here it is the user
            typ, n, m = P[K]
            A = np.full(sh, np.nan) #TODO:  we need separate total, rephasing and nonrephasing amplitudes
            T = np.full(sh, np.nan)
            W = np.full(sh, np.nan)
            for i in range(sh[0]):
                params = self._CA[i]
                if params is None:
                    continue
                #it would be probably best to internaly change all "+-" to "+" and "-" to avoid confusion with amplitude calculations
                ttt = []
                for lifetime, lFix, freqType, freq, freqFix in params:
                    if freqType == "±":
                        ttt.append((lifetime, lFix, "+", freq, freqFix))
                        ttt.append((lifetime, lFix, "-", freq, freqFix))#this is the order fitting_pso uses to split the component
                    else:
                        ttt.append((lifetime, lFix, freqType, freq, freqFix))

                for j in range(sh[1]):
                    print("\t\tpoint", i, j, params, n, m)

                    ampl = self._A[i,j]
                    if ampl is None:
                        continue

                    if len(ttt)==len(ampl):
                        params = [(a,*b) for a,b in zip(ampl, ttt)]
                    else:
                        print("tabCA.updateShow: _CA and _A do not match", self._CA[i], ampl)
                        params = [(np.nan,*b) for b in ttt]
                    
                    print("\t\tmodified params", params)
                    #now check if some (and only one) matches the criterion
                    a = np.nan 
                    t = np.nan
                    w = np.nan
                    for amplitude, lifetime, lFix, freqType, freq, freqFix in params:
                        if typ=="t" and freqType=="Ø" and n<=lifetime and lifetime<=m:
                            if t is np.nan:
                                t = lifetime
                                a = amplitude
                                # do not set w, it does not make sense for lifetime components (note that you cannot select oscillatory components based on their lifetime)
                            else:
                                t = -np.inf
                                a = np.inf
                                break #we have more that one component - possibly mark as inf? to distinguish from nan
                        elif typ=="w":
                            if freqType=="+" and n<=freq and freq<=m: #"+", "-", "±"
                                if t is np.nan:
                                    t = lifetime
                                    # ~ a = now we have a problem since we do not store amplitudes
                                    a = amplitude
                                    w = freq
                                else:
                                    t = -np.inf
                                    a = np.inf
                                    w = np.inf
                                    break #we have more that one component - possibly mark as inf? to distinguish from nan
                            elif freqType=="-" and n<=-freq and -freq<=m: #"+", "-", "±"
                                if t is np.nan:
                                    t = lifetime
                                    # ~ a = now we have a problem since we do not store amplitudes
                                    a = amplitude
                                    w = -freq
                                else:
                                    t = -np.inf
                                    a = np.inf
                                    w = np.inf
                                    break #we have more that one component - possibly mark as inf? to distinguish from nan
                    print("\t\ta", a, "t", t, "w", w)
                    A[i,j] = a
                    T[i,j] = t
                    W[i,j] = w
            print("\tadding a pseudoframe")
            TT.append(A)#should be AT
            TT.append(T)
            TT.append(W)
            RR.append(A)#shoudl be AR
            RR.append(T)
            RR.append(W)
            NN.append(A)#shoudl be AN
            NN.append(T)
            NN.append(W)
            #TODO: it is wastefull to repeat T and W but whatevevr
        
        self._dataShow = {"total":np.asarray(TT), "rephasing":np.asarray(RR), "nonrephasing":np.asarray(NN), "a1":self._dataInput["a1"], "a2":None, "a3":self._dataInput["a3"]}
        
        print("\t\ttotal pseudo frames all nan", [np.all(np.isnan(it)) for it in self._dataShow["total"]])
        self._updateDataBuffer() 
        self.updateOutput()
        pass        

    def updatePlotMain(self, *args):
        #todo: test projection

        super().updatePlotMain(*args)

        # we need to repeat the checks and get the data, but focus on time frames [raw, a, t, w, a, t, w, a, t, w, ...]
        #  so if index > 0, the nearest 2+3*N
        if self._currentData is None:
            return
        index = self.listFrames.current()
        if index is None or index == -1:
            return

        data, a1, a3, a2, dataRange, alpha = self._currentData
        if index > len(data):
            print("WARNING: tabTemplateSmall.updatePlotMain - mismatched index and data")
            return

        if index==0:
            self._column_lifetime_series.setVisible(False)
        else:
            index = ((index-1)//3)*3 + 2
            frame = np.nanmean(data[index], axis=1) # should be the same for all w3 for CA
            indices = np.isfinite(frame)
            self._column_lifetime_series.setData(a1[indices], frame[indices])
            self._column_lifetime_series.setVisible(True)
            self._column_lifetime_series._yaxis.zoomToData()


    def updateDissectionSidePlots(self, ix=None, iy=None):
        print("tabCA.updateDissectionSidePlots")
        res = super().updateDissectionSidePlots(ix, iy)
        if res is None: return
        ix, iy = res
        self._lastPoint = ix, iy
        #use this to set listLAInputSingle value
        # and let it call updateDissectionSidePlots2
        params = self._CA[ix]
        if params is None:
            self.listLAInputSingle.setValue([], True)
        else:
            self.listLAInputSingle.setValue(params, True)
        pass

    def _storeSingle(self, *args):
        #we have to find where to put the params
        if self._lastPoint is None: return
        ix, iy = self._lastPoint
        self._CA[ix] = self.listLAInputSingle.value()
        pass

    def updateDissectionSidePlots2(self):
        # ~ print("tabLA.updateDissectionSidePlots", ix, iy)
        #we need different approach here, since we do not show _currentData in plotTimeline
        #  we show _dataInput, but that also should be transformed TRN and PAIR ?
        if self._dataInput is None:
            return

        if self._lastPoint is None: return
        ix, iy = self._lastPoint

        tracer = self.plotMain.plot().decorator(0)
        w1, w3 = tracer.origin
        
        trn = TRN2trn[self.cbTRN.currentText()]
        pair = self.cbPAIR.currentText()
        
        data = PAIR[pair](self._dataInput[trn])
        a1 = self._dataInput["a1"]
        a2 = self._dataInput["a2"]
        a3 = self._dataInput["a3"]
        
        # ~ ix = a1.searchsorted(w1)
        # ~ iy = a3.searchsorted(w3)
        
        # ~ if ix<0: ix = 0
        # ~ elif ix>=data.shape[1]: ix = data.shape[1]-1
        # ~ if iy<0: iy = 0
        # ~ elif iy>=data.shape[2]: iy=data.shape[2]-1

        #timeline
        y = data[:, ix, iy]
        self.plotTimeline.plot().plotter(0).setData(a2, y)

        #model should not rescale axes (or at least not much, particularly for negative T)    
        self.plotTimeline.plot().plotter(1).setVisible(False)
        self.plotTimeline.plot().zoomToData(self.plotTimeline.plot().plotter(0)) #rescaleAxes(True)

        #model
        # ~ if self._LA[ix, iy] is not None:
            # ~ params = self._LA[ix, iy]

            # ~ #A = self._A[ix, iy]
            # ~ #lets try
            # ~ print("tabLA: model", params)
            # ~ model = fitting.model(y, a2, params)
            
            # ~ self.plotTimeline.plot().plotter(1).setData(a2, model) #not so sure this is exactly what LA calculates during fitting
            # ~ self.plotTimeline.plot().plotter(1).setVisible(True)
            # ~ pass
        
        params = self.listLAInputSingle.value()
        print("tabCA: model", params)
        
        if len(params)>0: #TODO: not sure if we should calculate model from PAIR(data) or rather do PAIR(model from data)
            model = fitting.model(y, a2, params)
            
            self.plotTimeline.plot().plotter(1).setData(a2, model) #not so sure this is exactly what LA calculates during fitting
            self.plotTimeline.plot().plotter(1).setVisible(True)
            pass
        pass
    pass

    pass
