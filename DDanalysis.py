
#***********************************************************************
#DD package, data collection and analysis of 2D electronic spectra
#Copyright (C) 2019 Jan Alster (Charles Univesity, Prague)
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

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This version uses pqr for images so it will require no external things beside python3, numpy and PyQt
"""


#NEXT STEP - check data calculation
#TODO: axis labels, remember axis ranges when changing rawShowType
#TODO: absolute units of 2D plots
#TODO: units - conversion to 1/cm
#TODO: add unit coversion to interpolator?
#TODO: add undo/redo (saving more states) to GA fitting (see QUndoCommand)
#TODO: DAS should store which 2D were disabled during calculation (and maybe all other parameters)

#TODO: (WIP) add state changed ability to SaveLoadElement
 
#TODO: add ability to hold (and compare) multiple datasets
 
#TODO: if we only load data (without analysis state, e.g. fresh from DDphasing)
# - we will not trigger UI set and plots will not update
# - also loading state that does not change saveloadelements values seems to not trigger plot update 
 
#TODO: add Notes (expandable/hideable side panel)
# - common for whole analysis in case more datasets are loaded
# - this would imply another format/file for saving analysis
# - probably best to continues saving GA into the original data files (.DD), but possibly also have analysis session file (which datasets are open, global notes)
# - we can have global notes for DD file, which will include notes from phasing

#TODO: open filename to window title (if more datasets open, add dataset name to State menu and show active dataset in window title / or analysis session filename)

#TODO: persistent marks for DAS analysis - energy levels shown on all main plots (switchable)
#Todo: diagonal mark for all main plots

#TODO: countour comparison, comparison tab (two main plots side by side), with loading additional datasets (one can be active for GA, others can be selected for comparison contours)
# if we should have more datasets open at the same time (possible DD are not that big), it would make sense to have calculaations part of DData, not UI

#TODO: add lockable and hideable overlay items like displaced diagonals, grids, etc. to help identify peaks, measure peak distance s etc

#TODO: saving state with data followed by saving state without data might lead to uncorrelated params and data

#TODO: correlation of GA params might be doable from swarm trace (not the same as bootstrap, but readily available)

#TODO: possible segfault, removing component from GA input by setting zero and leaving focus via press of TAB (not enter)

#TODO: GA has problems with fitting random noise rather than weak signals, perhaps we can adjust residues to add a measure of "continuous area"
# i.e. reject components with very noisy spectra and prefer components with peaks (continuous areas with similar amplitude and phase)

#TODO: there is discrepany in data export, context menu save DD will use pqr to directly save data from plot, ie. _currentData.
# DDanalysis main Export dialog will use _dataShow, which might be different. And it uses transposed data.
# This needs to be unified (not so sure which one of them is the bad guy here).
# Export dialog does PAIR and intepolation on its own. We need a unified approach so that the output is guaranteed to be the same.
# Right now, we get different scale from save DD and export.

#TODO: allow direct export of whole matrix (data, residues, OMs - does not make much sense for GA, but should work too) as .npy, not only sliced to 2D pieces in T

# TODO: add control to swarm to specify how wide the range of param space should be, so that the swarm does not start so far off
#  the starting location (or at least is more focused on the starting location if needed)

"""
TODO: on the problem with GA focusing on very fast kinetics to fit first few frames completely instead
of taking a bit of all frames

and related problme of weirdly scaled DAS amplitudes

actually it is not unexpected that DAS amplitudes scale with rate and fast rates leading to very large amplitudes
this all shoudl be connected to badly chosen t0

A*exp(-k*(t-t0)) = A*exp(-k*t)*exp(k*t0) = B(A, k, t0) * exp(-k*t)
We see that DAS amplitudes will be scaled differently depending on their rate and the global t0.
Of course we want t0 to be zero and that should lead to correct amplitudes, but
we do not know the correct experimental t2. 

I need to do some simulations and select a reasonable bounds for t0 and lowest allowed t0.

on LHCII_3 DD8sub dataset GA
if we change offset t0 but keep all others the same
the amplitude of total DAS changes depending on rate
           t0 /fs   0   10   -10    -50   
lifetime/fs
1.87              557776   6227657   84155   67477
313               258      266       250     220
8481              376      376       376     383

This is based on data with initial t2 part removed. But we cannot fit the pulse overlap region anyway,
because it is fraud with non-resonant signals which GA cannot model (yet).

Current version uses convolution with IRF so the situation will be slightly different, but if we ignore pulse overlap
region it will not be that different.
"""

"""
Current idea is to separate data processing and analysis.
Original Data is output from DDphasing. DDanalysis then adds (many versions of) analysis params (like global analysis lifetimes) into the same file.
Optionaly it can also store resulting spectra for later comparison.  (this is not perfect since it can enlarge data file significantly; unfortunately data calculation is tightly bound to UI and cannot be easily separated - me and my experimental software)


we need three stages for all tabs
- load data (possibly only Data/Timeline tab, data is stored by main window)
- set parameters (UI elements)
- make calculations, which will possibly trigger other changes

probably we need to keep order in which tabs update

I will have to solve the specifics as it comes :(
"""


from PyQt5.QtWidgets import QApplication
from PyQt5.uic import loadUiType

import sys
import os.path as op
import pickle

from PyQt5 import QtGui, QtWidgets, QtCore
from common import sigmoidWindow 
from analysis.DDdata2 import DDdata
from analysis.tabData import tabData
from analysis.tabGA import tabGA
from analysis.tabResidues import tabResidues
from analysis.tabOM import tabOM
from analysis.diaExport import ExportDialog
from analysis.tabCompare import tabCompare
from analysis.tabCentralLines import tabCentralLines
from analysis.tabLA import tabLA
from analysis.tabCA import tabCA

import numpy as np

from common import fitting
from common.interpolator import interpolate2D as interpolate  # TODO: interpolator needs more work
import time
import warnings

#for testing
def slot(*args, **kw):
    print("slot",time.clock(),":", args, kw)
    pass

class DDAMainWindow(QtWidgets.QMainWindow):
    closing = QtCore.pyqtSignal()
    def __init__(self, parent=None):
        super(DDAMainWindow, self).__init__(parent)

        #data container
        self.data = DDdata()

        #create menu
        # - File
        #    - Load data (possibly will also prompt selection of analysis state if some is present in the datafile)
        #    - Exit
        # - State  (only available if some datafile is loaded)
        #    - Save state (save analysis state into current state in the datafile - or triggers Save state as
        #    - Save state as (save analysis state under a new name into the datafile)
        #    ---------------
        #    - (list all available analysis states saved in the datefile, load upon click)
        menuFile = self.menuBar().addMenu("&File")
        menuFile.addAction("&Load data").triggered.connect(self.loadData)
        acCorrectLO = menuFile.addAction("Correct LO")
        menuFile.addAction("&Exit").triggered.connect(self.close)
        
        menuView = self.menuBar().addMenu("&View")
        acViewSigSignal = menuView.addAction("Significant Signal")
        acViewSigSignal.setCheckable(True)
        
        
        menuState = self.menuBar().addMenu("&State")
        menuState.addAction("&Save state").triggered.connect(self._saveState)
        menuState.addAction("&Save state as").triggered.connect(self._saveStateAs)
        menuState.addAction("&Save state with data").triggered.connect(self._saveStateWithData)
        menuState.addAction("&Save state with data as").triggered.connect(self._saveStateWithDataAs)
        #Remove state sub menu
        menuRemove = menuState.addMenu("&Remove")
        menuRemove.triggered.connect(self._removeState)
        menuState.addSeparator()
        self._menuStateStates = QtWidgets.QActionGroup(self)
        self._menuStateStates.triggered.connect(self._changeState)
        menuState.setEnabled(False)
        self._menuState = menuState
        self._menuRemove = menuRemove
        
        menuExport = self.menuBar().addAction("&Export")
        menuExport.setEnabled(False)
        self._menuExport = menuExport

        #create widgets
        central=QtWidgets.QWidget()
        layout=QtWidgets.QHBoxLayout()
        #tabs (user input/output)
        tabWidget=QtWidgets.QTabWidget()

        
        #observe (possibly manipulate, for now just disable some frames from fitting) data
        self.tabData = tabData(self)
        tabWidget.addTab(self.tabData, "Data")

        #central lines
        self.tabCentralLines = tabCentralLines(self)
        tabWidget.addTab(self.tabCentralLines, "Central Lines")
        self.tabData.outputChanged.connect(self.tabCentralLines.updateInput)

        self.tabLA = tabLA(self)
        tabWidget.addTab(self.tabLA, "LA")
        self.tabData.outputChanged.connect(self.tabLA.updateInput)

        self.tabCA = tabCA(self)
        tabWidget.addTab(self.tabCA, "CA")
        self.tabData.outputChanged.connect(self.tabCA.updateInput)

        #make GA and calculate model
        self.tabGA = tabGA(self)
        tabWidget.addTab(self.tabGA, "Global Analysis")
        self.tabData.outputChanged.connect(self.tabGA.updateInput)
        self.tabGA.outputChanged.connect(self.tabData.updateModel) #backchannel

        #show residues after GA (preparation for frequency analysis, probably will be mostly ignored)
        self.tabResidues = tabResidues(self)
        tabWidget.addTab(self.tabResidues, "Residues")
        self.tabGA.outputChanged.connect(self.tabResidues.updateInput)

        #calculate oscillation maps
        self.tabOM = tabOM(self)
        tabWidget.addTab(self.tabOM, "Oscillation Maps")
        self.tabResidues.outputChanged.connect(self.tabOM.updateInput)
        
        #compare datasets
        self.tabCompare = tabCompare(self)
        tabWidget.addTab(self.tabCompare, "Compare")
        
        #notes, connect from all tabs that have them
        self.tabData.pteNotes.textChanged.connect(self._notesChanged)
        self.tabLA.pteNotes.textChanged.connect(self._notesChanged)
        self.tabCA.pteNotes.textChanged.connect(self._notesChanged)
        self.tabGA.pteNotes.textChanged.connect(self._notesChanged)
        self.tabResidues.pteNotes.textChanged.connect(self._notesChanged)
        self.tabOM.pteNotes.textChanged.connect(self._notesChanged)

        self._editTimer = QtCore.QTimer()
        self._editTimer.setSingleShot(True)
        self._editTimer.timeout.connect(self._saveNotes)        
        
        #we might design tabs as dockable widgets, however they are likely to be large
        # - we can do that easily after if need be
        
        self.tabWidget = tabWidget
        self._tabs=[("t"+''.join(tabWidget.tabText(i).split()), tabWidget.widget(i)) for i in range(tabWidget.count())]

        layout.addWidget(tabWidget,1)
        layout.setContentsMargins(0,0,0,0)
        central.setLayout(layout)
        self.setCentralWidget(central)
        #order of disconnect is probably important
        #experiments first
        #for it in _experiments:
        #    self.closing.connect(it.requestInterruption)

        self.exportDialog=ExportDialog(self)
        self.exportDialog.leFolder.setText(QtCore.QDir.currentPath())
        menuExport.triggered.connect(self.export)
           
        self._acCorrectLO = acCorrectLO
        self._defaultState = self.getState()
        #~ print("default state")
        #~ print(self._defaultState)

        
        acCorrectLO.triggered.connect(self._correctLO)
        acCorrectLO.setCheckable(True)
        
        pass

    def _notesChanged(self, *args):
        #distribute the same to all tabs
        s = self.sender()
        t = s.toPlainText()

        for it in [self.tabData, self.tabLA, self.tabCA, self.tabGA, self.tabResidues, self.tabOM]:
            if it.pteNotes is not s:
                it.pteNotes.textChanged.disconnect(self._notesChanged)
                it.pteNotes.setPlainText(t)
                it.pteNotes.textChanged.connect(self._notesChanged)
        
        self._editTimer.start(1000) #start the timer to save the notes
        
    def _saveNotes(self, *args):
        #TODO: saving notes with every change can be incredibly slow on network files, with data stored in states, DD files can be several 10MB it would be best to just update the file, but how to do that?
        self.data.setNotes(self.tabData.pteNotes.toPlainText())
        self.data.save()

    def export(self):
        print("DDanalysis.export2")
        return self.exportDialog.exec_()

    def _correctLO(self, *agrs):
        if self.data.hasData():
            self._loadData()
        pass
           
    def loadData(self):
        print("load")
        # TODO: allow loading .bootstrap file and open the original file stored in its metadata so that bootstrap can be run under the same conditions
        filename=QtWidgets.QFileDialog.getOpenFileName(filter="DD files (*.DD)")[0]
        print(filename)
        if filename=="": return

        #here we should turn UI to default state
        if self.data.hasData():
            self.setState(self._defaultState)

        self.data.load(filename)
        if not self.data.hasData():
            print("could not load data", filename)
            #~ QtWidgets.QApplication.instance().setApplicationDisplayName('DD Analysis')
            self.setWindowTitle("")
            return
        
        print()
        print("setting up application display name", QtWidgets.QApplication.instance())
        #~ QtWidgets.QApplication.instance().setApplicationDisplayName('DD Analysis '+filename)
        self.setWindowTitle(filename)
 
        #update states
        print("lastState", self.data.lastState)
        self._menuState.setEnabled(True)
        lastStateAc = self._updateMenuState()
        print("loadData: states", self.data.listStates())
        self._menuExport.setEnabled(True)

        self._loadData()

        #load last state
        if self.data.lastState is not None:
            #select last state
            self._changeState(lastStateAc, True)
        else:
            #TODO? load some defaults, defaults could be hardcoded into tabs
            print("no state, possibly blank UI")
            ...
        
        #notes (are not part of state)
        self.tabData.pteNotes.setPlainText(self.data.notes()) #should be distributed to all tabs, also will save notes back to data, but whatever
        #we need to trigger change signal, to load notes into all tabs, but stop saving notes back to datafile
        self._editTimer.stop() #this might not be fast enough, but it should work most of the time
        
    def _loadData(self):
        #data first
        T,R,N,a1,a3,a2 = self.data.data("all")
        LO = self.data.LO()
        if LO is None:
            data = {"total":T, "rephasing":R, "nonrephasing":N, "a1":a1, "a2":a2, "a3":a3, "LO":None}
        else:
            for it in LO:
                print(it[0:5], it[-5:])
            if self._acCorrectLO.isChecked() and False: #TODO: temp. disabled before correctLO in tabData is tested, remove after
                LOi = LO.sum(1)
                #note that DDphasing already divides data by LO**0.5 in one direction (likely w3)
                # so v1 makes more sense than v2, although LO ~ E**2 and 2DES ~ E**3
                # also relative peak amplitudes are not corrected for laser spectral profile in w1 direction
                # and detected LO goes through the sample
                
                # v1 
                cor = LOi.reshape((-1,1,1))/LOi.max()
                data = {"total":T/cor, "rephasing":R/cor, "nonrephasing":N/cor, "a1":a1, "a2":a2, "a3":a3, "LO":LO/LOi.reshape((-1,1))*LOi.max()}
                # v2 
                #~ cor = ((LOi/LOi[0])**(3/2)).reshape((-1,1,1))
                #~ data = {"total":T/cor, "rephasing":R/cor, "nonrephasing":N/cor, "a1":a1, "a2":a2, "a3":a3, "LO":LO/LOi.reshape((-1,1))}
                # v3
                
                #~ print("_loadData", T.shape, LO.shape)
                #~ #T.shape = t2, w1, w3
                #~ #  w3 is the same as for LO
                #~ #  but w1 is not, we need to interpolate LO on w1 axis to calculate v2
                #~ LO1 = np.array([np.interp(a1, a3, it) for it in LO])
                #~ cor = np.array([np.outer(it1, it2) for it1,it2 in zip(LO1, LO)])
                #~ data = {"total":T/cor, "rephasing":R/cor, "nonrephasing":N/cor, "a1":a1, "a2":a2, "a3":a3, "LO":LO/LO}
            else:
                data = {"total":T, "rephasing":R, "nonrephasing":N, "a1":a1, "a2":a2, "a3":a3, "LO":LO}

        print("DDanalysis2.loadData: calling tabData.updateInput")
        self.tabData.updateInput(data)
        pass
     
    def _updateMenuState(self):
        #remove everything after the separator (self._menuStateBookmark)
        for ac in self._menuStateStates.actions():
            self._menuState.removeAction(ac)
            self._menuStateStates.removeAction(ac)
            ac.setParent(None)
            del ac
        
        for rac in self._menuRemove.actions():
            self._menuRemove.removeAction(rac)
            rac.setParent(None)
            del rac
            
        #add item per every state
        #items should be checkable and unique in the group
        # - either none or only one can be selected at the same time
        states = self.data.listStates()
        lastStateAc = None
        lastState = self.data.lastState
        for name in states:
            ac = self._menuState.addAction(name)
            ac.setCheckable(True)
            self._menuStateStates.addAction(ac)
            #test if there is data (this is weird, just getting the state from data will mark it as last state)
            # so basicaly I cannot check the state ...
            # problem is that both "data in state" and "last state" are concept of analysis UI (saved in data for convenience) and not data
            S = self.data.state(name)
            if ('tData' in S and "data" in S['tData']) or ('tGlobalAnalysis' in S and "data" in S['tGlobalAnalysis']) or ('tResidues' in S and "data" in S['tResidues']) or ('tOscillationMaps' in S and "data" in S['tOscillationMaps']):
                f = ac.font()
                f.setWeight(f.DemiBold)
                ac.setFont(f)
            
            if name == lastState:
                lastStateAc = ac
                ac.setChecked(True)
            #add to remove menu
            rac = self._menuRemove.addAction(name)
            rac.setData(ac)
        self.data.lastState = lastState
        return lastStateAc
        
    def _removeState(self, removeAc):
        stateAc = removeAc.data()
        name = stateAc.text()
        self._menuRemove.removeAction(removeAc)
        removeAc.setParent(None)
        del removeAc
        self._menuState.removeAction(stateAc)
        self._menuStateStates.removeAction(stateAc)
        stateAc.setParent(None)
        del stateAc
        self.data.removeState(name)
        self.data.save()
        pass
        
    def _saveState(self, *args):
        self._saveStateCommon(self.data.lastState)
        pass
        
    def _saveStateAs(self, *args):
        self._saveStateCommon()
        pass

    def _saveStateWithData(self, *args):
        self._saveStateCommon(self.data.lastState, saveData=True)
        pass
        
    def _saveStateWithDataAs(self, *args):
        self._saveStateCommon(saveData=True)
        pass
        
    def _saveStateCommon(self, name=None, **kw):
        if name is None:
            name, ok = QtWidgets.QInputDialog.getText(self, "Save state as...", "Save analysis state as: ")
            if not ok: return
       
        state = self.getState(**kw)
        self.data.setState(name, state)
        self.data.save() #should save under original filename
        self._updateMenuState() #Todo: it is not necessary to rebuild the whole menu every time 
        pass        
        
    def _changeState(self, stateAc, force=False):
        print("DDanalysis2._changeState:", stateAc.text())
        #triggered by menuState state items
        #(optionaly) ask to save current state (if there is some change?)
        #load new state
        if force or (stateAc.text() != self.data.lastState):
            stateAc.setChecked(True)
            self.setState(self.data.state(stateAc.text()))
        pass

    def getState(self, **kw):
        #iterate over SaveLoadElements present (Todo: this might be automated by SaveLoadElement)
        state = self.saveSettings()
        for name, element in self._tabs:
            print(name)
            state[name] = element.saveSettings(**kw)
            pass
        return state    
    
    def setState(self, state):
        print("DDanalysis3.setState")
        #mark UI state of all tabs as undefined (this will prevent update of UI until we are ready for it)
        for name, element in self._tabs:
            print(name, element)
            if name in state: element.unset()

        # here we have a problem
        #  since values are stored in UI (a bad idea, really), we do not keep the original default value
        #  if we then unset everything and load data that do not contain correct params
        #  (this will happen if you add a new tab and load data from previous version)
        #  the params will remain unset
        # admittedly, this does not happen often, but it can make unrecoverable error (if some UI elements are disabled
        #  as is the case e.g. in tabCA)
        # as a workaround, do not unset elements for which we do not have values in stata
        #  this will leave the previous (not necasarily the default) values, but it should be better than leaving
        #  the params unset with no way of changing them

        for name, element in self._tabs:
            print("\n\tname", name)
            if name in state: element.loadSettings(state[name]) 
            #loading UI state should trigger appropriate update(s) in correct order
            pass

        self.loadSettings(state) #TODO: loading main window size does not go well with tiling window manager

        #for name, element in self._tabs:
        #    element.update()
        #    pass
        pass
        

    #--reimplemented from SaveLoadElement-------------------------------------
    #settings specific to main UI window
    def saveSettings(self):
        #what will be saved to session config
        settings = {}
        t = self.size()
        #~ settings['window size'] = t.width(), t.height() #this leads to undesirabel behavior when transferring data to another computer
        settings["tab"] = self.tabWidget.currentIndex()
        #settings["notes"] = self.eNote.toPlainText()
        settings["corLO"] = self._acCorrectLO.isChecked()
        return settings

    def loadSettings(self, settings):
        #~ if 'window size' in settings: self.resize(*tuple(settings['window size']))
        if "tab" in settings: self.tabWidget.setCurrentIndex(settings["tab"])
        self._acCorrectLO.setChecked(settings.pop("corLO", False))
        #if "notes" in settings: self.eNote.setPlainText(settings["notes"])
        pass
    #-----------------------------------------------------------------------    

    #notify others that we are closing
    def closeEvent(self, ev):
        print("DDAMainWindow.closeEvent")
        #save most settings to session file
        #self.closeSession()

        #saving last session info into global configs
        #config=Config(CONFIG_PATH)
        #config['last session']=self._path
        #QtWidgets.QApplication.instance().config.save() #last session is set up during __init__ but do not save until closing
        
        self.closing.emit() #however do not emit this if you prevent closing window here
        return super(DDAMainWindow, self).closeEvent(ev)
    pass
    
    
class DDAApplication(QtWidgets.QApplication):
    def __init__(self, *args):
        super().__init__(*args)
        self.setApplicationDisplayName('DD Analysis')

        #load default (global) config  (before window)
        #self.config = Config(CONFIG_PATH) #this is hardcoded

        self._window = DDAMainWindow()
        self._window.show()
        pass
    
    def sessionPath(self):
        #return self.config.get("last session", None)
        return
        
    def run(self):
        return self.exec_()

        

if __name__ == "__main__":
    import faulthandler
    faulthandler.enable()
    app = DDAApplication(sys.argv)
    sys.exit(app.run())
    
