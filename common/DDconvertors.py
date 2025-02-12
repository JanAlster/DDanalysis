
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
try:
    from common import n, c
except:
    from . import n, c

try:
    from common.conversion import Convertor
    from common.widgets.conversionWidgets import SelectorFactory #ConvertorSelector
except:
    from .conversion import Convertor
    from .widgets.conversionWidgets import SelectorFactory #ConvertorSelector
    
import numpy as np

## CCD convertors ######################################################
#TODO: CCD convertors could use a default hardcoded (not connected to any file) convertor
# - as it is we need a (empty) calibration file for both CCD height and CCD amplitude
class CCDHeightConvertor(Convertor):
    pass
    
class CCDAmplitudeConvertor(Convertor):
    pass

class SpectralConvertor(Convertor):
    """CCD spectral calibration"""
    def _addUnits(self, group):
        super()._addUnits(group)
        #add mm
        params = self._params[group]
        if group == "calpol":
            #params is calibration polynom x[nm] = np.polyval(params, x[px])
            self._params[group] = params
            #unit, unit2base, base2unit
            self._addUnit("nm", NotImplemented, lambda x: np.polyval(params, x), "Wavelength")
            self._addUnit("rad/fs", NotImplemented, lambda x: 2*np.pi*c/np.polyval(params, x), "Frequency")
        pass

class CCDHeightConvertorSelector(SelectorFactory, CCDHeightConvertor):
    pass
    
class CCDAmplitudeConvertorSelector(SelectorFactory, CCDAmplitudeConvertor):
    pass
    
class SpectralConvertorSelector(SelectorFactory, SpectralConvertor):
    pass


## DS convertors #######################################################
class M406Convertor(Convertor):
    def _addUnits(self, group):
        super()._addUnits(group)
        #add mm
        params = self._params[group]
        if group == "counts":
            #params is [number counts/mm, number zero position in base units]
            #backward compatibility
            try:
                len(params)
            except TypeError:
                params = [params, 0]
            print("M406Convertor._addUnits", group, params)
            self._params[group] = params
            #unit, unit2base, base2unit
            self._addUnit("mm", lambda x: x*params[0]+params[1], lambda x: (x-params[1])/params[0])
            koef = -params[0] / (2*n(600)/c*1e6) #2x air space: dt = 2*x*n/c, n_air is assumed flat enough to ignore wavelength dependence, minus is to flip axis direction
            self._addUnit("fs", lambda x: x*koef+params[1], lambda x: (x-params[1])/koef, "Delay")
        pass

class M406ConvertorSelector(SelectorFactory, M406Convertor):
    pass

class M112Convertor(Convertor):
    def _addUnits(self, group):
        super()._addUnits(group)
        params = self._params[group]
        #add mm
        if group == "mm":
            #params is number nm/imp
            self._params[group] = params
            self._addUnit("mm", lambda x: x/(params*1e-6), lambda x: x*(params*1e-6))#6.88 nm/imp

        #add fs (simple)
        if group == "fs":
            #params is number fs/imp
            self._params[group] = params
            self._addUnit("fs (simple)", lambda x: x/params, lambda x: x*params, "Delay")
            
        #add fs (wedge cal)
        if group == "wedgecal":
            #params is spline
            #Todo: how to save spline? - we could pickle the whole thing, but it might lead to problems later on (cannot be unpickled without class definition and if something changes...)
            # - yeh I just got this problem :(
            self._params[group] = params
            self._addUnit("fs", lambda x: params.inverse(x), lambda x: params.model(x), "Delay")
            
        pass

class M112ConvertorSelector(SelectorFactory, M112Convertor):
    pass
