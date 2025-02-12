
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

#TODO: add saved flag

# from PyQt5 import QtCore, QtWidgets
import sys

from PyQt5 import QtWidgets

if sys.version_info<(3,8): #3.8 has new picke protocol, pickle5 is backport to 3.7 and lower
    import pickle5 as pickle
else:
    import pickle
import numpy as np
from common.interpolator import interpolate2D as interpolate  # TODO: interpolator needs more work
import warnings
from common import busyCursor

PAIR={'Real': lambda x: x.real, 'Imag': lambda x: x.imag, 'Ampl': lambda x: np.abs(x), 'Phase': lambda x: np.angle(x), None:lambda x:x}

DEBUG=False

#somewhat changed design from original
# - save multiple UI/analysis states within the datafile to allow multiple versions of analysis
# - do not keep calculated datasets, only original data and analysis parameters
# - do we even need it? thgis could be handled directly by app
# - we can just to this as an access wrapper for data from PP
# - on the other hand we might store here data from individual parts of analysis (just not store them to datafile)
#   - we need some central place to keep track, or do we?

# class DDdata(QtCore.QObject):
#     """Encapsulation of DD data, with Qt signals on data change."""
    #data from DDphasing (these will not be touched by DDanalysis, but are kept to see how the data were phased)
    # originalConfigFile
    # data   : total, rephasing, nonrephasing, t2, w1, w3, LO
    # params : filename, windowing, cropping, phasing, PP
    
    #states of DDanalysis as output of SaveLoadElement (basically DDdata does not care about this)
    
    #signals
    #dataChanged  = QtCore.pyqtSignal() #when selected data change (the one plotted on main plot)
        

class DDdata:
    # version without Qt, it seems that Qt signals were not used anyway
    def __init__(self, *args):
        super(DDdata, self).__init__(*args)
        self._data = {}
        self._filename = None
        pass

    def __contains__(self, x):
        return x in self._data

    @busyCursor
    def load(self, filename):
        #test if not pickled file
        try:
            with open(filename, 'rb') as ff: #do not forget binary flag when opening
                data=pickle.load(ff, encoding='latin1') #encoding needed to load numpy arrays (this is python2 to python3 pickle incompatibility)
                print (data.keys())
                pass
            pass
        except:
            print("unknown data format")
            import traceback
            traceback.print_exc()
            data = {}
            raise
            pass

        if "states" not in data:
            data["states"] = {}
            self.lastState = None
        else:
            self.lastState = data["lastState"]

        if "notes" not in data:
            data["notes"] = ""

        if False: #reverse old Lund data
            data["data"]["total"] *= -1
            data["data"]["rephasing"] *= -1
            data["data"]["nonrephasing"] *= -1

        self._filename = filename
        self._data = data
        pass

    def save(self, filename=None):
        if filename is None:
            filename = self._filename
        if filename is None:
            return
        self._data["lastState"] = self.lastState
        try:
            with open(filename, 'wb') as ff: #do not forget binary flag when opening
                pickle.dump(self._data, ff, -1)
                pass
        except:
            import traceback
            # FIXME: unresolved import
            QtWidgets.QMessageBox.warning(None, "Data save failed", "Data save failed with following exception:\n"+traceback.format_exc())
        pass

    def hasData(self):
        return "data" in self._data and "states" in self._data


    #data access
    def value(self, *args):
        temp = self._data
        for it in args:
            temp = temp[it]
        return temp

    #bit more explicit than value
    def data(self, trn, pair=None):
        data = self._data["data"]
        a1 = data["w1"]
        a3 = data["w3"]
        a2 = data["t2"]
        if trn=="all":
            return PAIR[pair](data["total"]), PAIR[pair](data["rephasing"]), PAIR[pair](data["nonrephasing"]), a1, a3, a2
        else:
            return PAIR[pair](data[trn]), a1, a3, a2
        pass

    def LO(self):
        return self._data["data"]["LO"] if "LO" in self._data["data"] else None

    def total(self, pair=None):
        return PAIR[pair](self._data["data"]["total"])

    def rephasing(self, pair=None):
        return PAIR[pair](self._data["data"]["rephasing"])

    def nonrephasing(self, pair=None):
        return PAIR[pair](self._data["data"]["nonrephasing"])

    
    #state access
    def listStates(self):
        return list(self._data["states"].keys())

    def state(self, state):
        self.lastState = state
        return self._data["states"][state]
        
    def setState(self, state, data):
        self.lastState = state
        self._data["states"][state] = data
        
    def setNotes(self, notes):
        self._data["notes"] = notes
        
    def notes(self):
        return self._data["notes"]
        
    def removeState(self, state):
        if self.lastState == state: self.lastState=None
        del self._data["states"][state]
        pass

    pass
