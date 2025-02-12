"""

bit more flexible container for CCD data

what I want is an unified way how to export data with proper conversion of axes (Convertors)
if we do not have an universal container, we will need to have that code at multiple places which is bad


CCD frame - single frame
CCD scan - multiple frames with arbitrary scan axis
    - scan axis will be hard to convert if it is completely arbitrary
    - scan axis can be also restricted to one DS axis
CCD scan cuts / projections / averages
    - i.e. average over CCD height


in genaral

arbitrary number of dimensions (in reality probably 2-4)
arbitrary type of dimensions

each dimension has its axis (values and Convertor and label (to be shown in plots) and maybe ID(for easy access))

we need ability to reserve space beforehand and fill in data in steps

we need ability to create cuts and export data

convertors should be saved in data not by reference
"""

import numpy as np
import pickle
import os.path as op

#from .conversion import Convertor
#TODO: can this replace DDScanData too?

#note that axis objects need to be pickleable
#TODO: not sure how convertors go with that
# - there will be problems if class definitions are updated, we will not be able to load older data anymore
"""
If that's the case… to make this work, you'd need to do some trickery. 
First you get the source for the old class definition, and when you grab the raw pickle 
you need to change the reference from the existing class to match the path to the code 
for the old version of the class. This should be in clear text (even in HIGHEST_PROTOCOL), 
so it shouldn't be an issue to grab and edit that part of the pickle string. You then would 
be able to unpickle the old objects, but they'd be pointing to the old class definition. 
A "format converter" would be needed to convert the old objects you'd be unpicking to the new 
class instance objects -- basically create new class instances that take the relevant state 
from the old class instances.
"""
# one way would be to not change the class definition and make a new version of the class
# (and change code everywhere else to use the new version)
# (use name alias to not have to change the code elsewhere)

#to keep this general, we need to allow multiple convertors per axis
# - because ScanAxis has multiple convertors
# - otherwise we would need to create UI elements here which I do not want
# - do not require units to be always list, if derived axis only uses one (which should be the majority of cases) allow it to keep scalar arguments (python does not care)
# - only NDContrainerAxis.convertors need to be vector
# - or rather make a single combination convertor for ScanAxis? that makes more sense

#OR we can save using derived classes which will have fixed axes types, so that we do not have to pickle
# instances directly

# FIXME: interface for the containter is flawed, we might need to remove background from CCD signal before
#  exporting, but  there is not way how to do that (see 2022_phycob/baldet_test2.py)

class NDContainerAxis_v0(object):
    """
    base class and minimal interface for NDContainer axes
    """
    def __init__(self, convertor, axisID, label="index"):
        self.axisID = axisID
        self.label = label
        self.convertor = convertor
        pass
        
    def length(self):
        raise NotImplementedError
        return N
        
    def convert(self, unit, iAxis, data): #TODO: this might return label updated with units too
        """
        convert axis to "unit" units
        if needed manipulate data in "iAxis" dimension 
        
        return axis expressed in units and manipulated data (or original if manipulation is not needed)
        """
        raise NotImplementedError
        return uAxis, modData
            
    def cut(self, index, unit):
        """
        information needed for elimination of axis from data
        data are to be cut at this axis at position index
        
        return position of axis expressed in units, data scaling factor if applicable (1 otherwise)
        """
        raise NotImplementedError
        return uAxisValue, dataScale

    def values(self, unit):
        """
        convert axis to "unit" units
        (will not manipulate data even if needed for correct conversion - use self.convert for that)
        return axis expressed in units
        """
        raise NotImplementedError
        return uAxis
        
    def formatPosition(self, index, unit):
        raise NotImplementedError
    pass
NDContainerAxis = NDContainerAxis_v0


class RangeAxis_v0(NDContainerAxis):
    def __init__(self, N, *arg, **kw):
        self._N = int(N)
        super().__init__(None, *arg, **kw)
        pass
        
    def length(self):
        return self._N
        
    def convert(self, unit, i, data):
        return np.arange(self._N), data
            
    def cut(self, index, unit):
        return index, 1
        
    def values(self, unit):
        return np.arange(self._N)
        
    def formatPosition(self, i, unit):
        return str(i)
    pass
RangeAxis = RangeAxis_v0
    
class CCDAxis_v0(NDContainerAxis):
    #this is basically spectral density type of axis, dividing each value by size of corresponding interval (in whatever units)
    #however this can lead to descending axis and negative spectral density
    #TODO: using abs for calculation of interval size should solve that, but is not entirely consistent :(
    
    #TODO: there is a problem here - with consistent definition of spectral density returned in scanData, we would get negative spectral density for spectral axis
    #   expressed in rad/fs - because the axis is descending (for standart setup of hardware with axcending wavelengths with ascending pixels on CCD)
    #   however, the FT transformations in calculation of wedge calibration (wedge thickness from interferrogram) implicitely assume ascending(incresing) rad/fs axis
    #   therefore we need to flip the axis somewhere and also make the spectral density positive
    
    #where to do that?
    # CCDAxis can be specificaly asked for single value in whatever units at given index, should we flip the axis there too? only for some units?
    # should we make spectral density positive for all units, knowing it is somewhat wrong?
    
    #compromise here is to have spectral density positive not matter the axis direction
    # and keep axis direction consistent with base units (i.e. descending rad/fs in our case)
    # flipping will be done in UI if needed and in calculation code as needed
    # (we cannot assume here when it will be needed)
    # just keep in mind that spectral density returned does not depend on axis direction
    # (note that this is not consistent with CCDScanData which also flipped the axis)
    
    def __init__(self, ROI, convertor, *arg, **kw):
        print("args", convertor, *arg, **kw)
        super().__init__(convertor, *arg, **kw)
        self._ROI = ROI #so it can be saved exactly
        N = ROI[1]/ROI[2]
        
        self._edges = ROI[0] - 0.5 + np.arange(N+1) * ROI[2]
        self._centers = ROI[0] -.5 + ROI[2]/2 + np.arange(N) * ROI[2]
        # self._edges = ROI[0] - 0.5 + np.arange(N+1)
        # self._centers = ROI[0] + np.arange(N)
        pass
    
    def convert(self, unit, i, data):
        edges = self.convertor.V2U(self._edges, unit)
        dedges = np.diff(edges)
        centers = self.convertor.V2U(self._centers, unit)
        
        #this need to divide in the right dimension
        dim_array = np.ones((data.ndim,), int)
        dim_array[i] = -1
        data = data / np.abs(dedges).reshape(dim_array) 
        return centers, data
    
    def cut(self, index, unit):
        edges = self.convertor.V2U(self._edges[index:index+2], unit)
        dedges = np.diff(edges)
        center = self.convertor.V2U(self._centers[index], unit)
        return center, 1./np.abs(dedges[0])
        
    def values(self, unit):
        return self.convertor.V2U(self._centers, unit)
        
    def length(self):
        return len(self._centers)
        
    def formatPosition(self, index, unit):
        center = self.convertor.V2U(self._centers[index], unit)
        return "{}{}".format(center, unit)
    pass
CCDAxis = CCDAxis_v0
    
class CombinationConvertor(object):
    def __init__(self, convertors):
        self._convertors = convertors
        #combination of all
        from itertools import product
        self._units = [",".join(it) for it in product(*(it.units() for it in self._convertors))]
        self._baseUnit = ",".join((it.baseUnit() for it in self._convertors))
        pass
    
    def convertors(self):
        return self._convertors
    
    def __getitem__(self, key):
        return self._convertors[key]

    def baseUnit(self):
        return self._baseUnit
        
    def units(self):
        return self._units
    pass
    
class ScanAxis_v0(NDContainerAxis):
    """
    at least one and up to three DS can move at the same time
    DS that are static are represented by scalar value (position)
    """
    def __init__(self, scanAxis1, scanAxis2, scanAxis3, convertor1, convertor2, convertor3, *args, **kw):
        super().__init__(CombinationConvertor([convertor1, convertor2, convertor3]), *args, **kw)
        #note that all scanAxes should be either same length array or a scalar
        # arrays of different lengths will lead to problems
        # we can sanitize the array lengths here to prevent problems
        #find scanAxes lengths
        def nlen(a):
            try:
                return len(a)
            except:
                return 0
        
        Ns = [nlen(it) for it in (scanAxis1, scanAxis2, scanAxis3)]
        N = np.max(Ns)
        if N==0:
            N = 1
        
        def sanitize(axis, aN):
            if aN==0:
                ttt = np.array([axis]*N)
            elif aN<N:
                axis = np.array(axis)
                ttt = np.ones((N,), dtype=axis.dtype) * axis[-1]
                ttt[:aN] = axis
            else:
                ttt = np.array(axis)
            return ttt
        
        scanAxis1 = sanitize(scanAxis1, Ns[0])
        scanAxis2 = sanitize(scanAxis2, Ns[1])
        scanAxis3 = sanitize(scanAxis3, Ns[2])
        
        self._originalAxesLengths = Ns #in case we need this info later
        self._scanAxes = (scanAxis1, scanAxis2, scanAxis3)
        
        self._length = len(self._scanAxes[0])
        pass
        
    def length(self):
        return self._length
        
    #no cut factor, this is not density axis
    
    def convert(self, units, i, data):
        units = units.split(",")
        axes = [convertor.V2U(axis, unit) for axis, unit, convertor in zip(self._scanAxes, units, self.convertor)]
        return axes, data

    def cut(self, index, units):
        units = units.split(",")
        axes = [convertor.V2U(axis[index], unit) for axis, unit, convertor in zip(self._scanAxes, units, self.convertor)]
        return axes, 1

    def values(self, unit):
        units = units.split(",")
        axes = [convertor.V2U(axis, unit) for axis, unit, convertor in zip(self._scanAxes, units, self.convertor)]
        return axes
        
    def formatPosition(self, index, units):
        units = units.split(",")
        axes = [convertor.V2U(axis[index], unit) for axis, unit, convertor in zip(self._scanAxes, units, self.convertor)]
        return ", ".join(["{}{}".format(*it) for it in zip(axes, units)])
    pass
ScanAxis = ScanAxis_v0

class NDContainer_v0(object): #TODO: possibly add convertor for data?
    _extensions = {} #registered subclass extensions to be used for fromFile
    
    def __init__(self, valueLabel, valueConvertor=None, *args, **kw):
        self._axes = []
        self._axesIDs = {}
        self._data = None
        
        self._valueLabel = valueLabel
        self._valueConvertor = valueConvertor
        
        #TODO: process args, kw to fill axes and data
        pass
        
    #first add all axes (or load from file)
    # loading into already filled CCDContainer is undefined
    def addAxis(self, containerAxis): 
        assert(self._data is None) #no adding dimensions after data added/space reserved
        axisID = containerAxis.axisID
        assert(axisID not in self._axesIDs) #each axis in container has to have unique ID
        self._axesIDs[axisID] = len(self._axes)
        self._axes.append(containerAxis)
        pass
        
    #then fill in data 
    #either directly
    def setData(self, data):
        #axes has to be specified before this
        #data needs to have the right dimension
        data = np.array(data) #TODO: force copy?
        print(data.shape, tuple(axis.length() for axis in self._axes))
        assert(data.shape == tuple(axis.length() for axis in self._axes))
        self._data = data
        pass

    #or step-by-step
    def reserveDataSpace(self, fill=None, *args, **kw): #allow passing dtype and order to np.empty
        self._data = np.empty([axis.length() for axis in self._axes], *args, **kw)
        if fill is not None:
            self._data.fill(fill)
        #TODO: mark (un)filled data
        pass

    def _normAxis(self, axisID):
        return self._axesIDs[axisID] if axisID in self._axesIDs else axisID

    def setDataStep(self, axisID, axisIndex, data):
        """
        input:
         - axisID - either axisID or axis index specifying which axis to fill at
         - axisIndex - where to fill on that axis
         - data - with proper dimension
        
        data space has to be reserved first
        no checks of overwrite
        """
        iAxis = self._normAxis(axisID)
        s = [slice(None)]*len(self._axes)
        s[iAxis] = axisIndex
        self._data[s] = data
        pass

    #after you can manipulate data
    def cut(self, axisID, axisIndex, units, label=None):
        """
        input
         - axisID - either axisID or axis index specifying which axis to cut at
         - axisIndex - where to cut on that axis
         
         This will eliminate this axis from resulting object (one dimension less).
         Note that this might require manipulation of the data (e.g. for density type axis).
         So units has to be specified.
         
         Also should return axis point at axisIndex in units.
         
        """
        iAxis = self._normAxis(axisID)
        cut = NDContainer(label if label is not None else self._valueLabel, self._valueConvertor)
        for j in range(len(self._axes)):
            if j==iAxis: continue
            cut.addAxis(self._axes[j])
            
        s = [slice(None)]*len(self._axes)
        s[iAxis] = axisIndex
        data = self._data[s] #TODO: is this view or copy? I need it to be a copy
        
        cutPointInUnits, cutScaleFactor = self._axes[iAxis].cut(axisIndex, units)
        if cutScaleFactor!=1:
            data *= cutScaleFactor
            
        cut.setData(data)
        #problem is that we might need to preserve something from the axis even if we cut it
        # e.g. pixel width of CCD which need to be converted to units and data divided by it
        # by using cut this transformation should go from axis to data value convertor
        # or we can fix the units for the cut here
        # :( apart from this it seemed nicely general concept 
            
        return cut, cutPointInUnits
        
    def multiCut(self, multiIndex, label=None):
        print("multiCut", multiIndex, label)
        """
        multiIndex = [(axisID, axisIndex, units), ...]
        
        specify axis, position, unit for arbitrary number of axes
        and eliminate those axes from data
        i.e. make cut through data at those points
        keep just unspecified axes
        """
        s = [slice(None)]*len(self._axes)
        S = 1.
        skipAxes = []
        cutPoints = []
        for axisID, axisIndex, units in multiIndex:
            iAxis = self._normAxis(axisID)
            s[iAxis] = axisIndex
            cutPointInUnits, cutScaleFactor = self._axes[iAxis].cut(axisIndex, units)
            S *= cutScaleFactor
            skipAxes.append(iAxis)
            cutPoints.append(cutPointInUnits)
            pass
            
        cut = NDContainer(label if label is not None else self._valueLabel, self._valueConvertor)
        for j in range(len(self._axes)):
            if j in skipAxes: continue
            cut.addAxis(self._axes[j])
        data = self._data[s] #TODO: is this view or copy (this should be basic slicing and therefore a view)? I need it to be a copy - maybe we can use both versions
        if S!=1:
            data *= S
        cut.setData(data)
        return cut, cutPoints #same order as multiIndex
        
    def average(self, axisID, units, label=None):
        #TODO: possibly units = None by defaults and average over base units?
        iAxis = self._normAxis(axisID)
        #export along the axisID and average
        uaxis, data = self._axes[iAxis].convert(units, iAxis, self._data)
        data = data.mean(iAxis)
        
        average = NDContainer(label if label is not None else self._valueLabel, self._valueConvertor)
        for j in range(len(self._axes)):
            if j==iAxis: continue
            average.addAxis(self._axes[j])
        average.setData(data)
        print("NDContainer.average", self.axes(), average.axes())
        return average
        
    def truncate(self, axisID, axisIndexFrom, axisIndexTo):
        ...
        
    def export(self, units):
        """
        input:
         - units for all axes
         
        needs to be able to manipulate data with axis unit conversion
        
        returns something that does not depend on convertors
        """
        data = self._data
        axes = []
        
        print(len(units), len(self._axes), self.axes())
        assert(len(units)==len(self._axes))
        for i in range(len(self._axes)):
            uaxis, data = self._axes[i].convert(units[i], i, data)
            axes.append(uaxis)
        
        return data, axes

    def swapAxes(self, axisID1, axisID2):
        iAxis1 = self._normAxis(axisID1)
        iAxis2 = self._normAxis(axisID2)
        
        print(self._axes, self._data.shape)
        
        self._axes[iAxis1], self._axes[iAxis2] = self._axes[iAxis2], self._axes[iAxis1]
        print("swapAxes")
        print(self._data.shape)
        self._data = self._data.swapaxes(iAxis1, iAxis2)
        print(self._data.shape)
        
        pass

    #or store for later
    def save(self, filename):
        #cannot easily save arbitrary axes, so this is reserved for overloaded classes only
        raise NotImplementedError
        
    def _sanitizeFilename(self, filename):
        #this expexts that the derived class has attribute "extension" and we will
        #check if the filename conforms to that extension
        if not filename.endswith(self.extension):
            return filename+self.extension
        else:
            return filename

    @classmethod
    def fromFile(cls, filename):    
        """
        NDContainer cannot load data directly, because it does not know what data structure to expect. However it can call fromFile for registered subclasses. Note that returned instance is a subclass of NDContainer.
        """
        #cannot easily save arbitrary axes, so this is reserved for overloaded classes only
        ext = op.splitext(filename)[1]
        if ext in cls._extensions:
            return cls._extensions[ext].fromFile(filename)
        raise ValueError("NDContainer_v0.fromFile cannot open", filename, "; its extension is not a registered subclass extension.")
    
    @classmethod
    def registerExtension(cls, subclass):
        print("NDContainer_v0.registerExtension", cls, subclass)
        cls._extensions[subclass.extension] = subclass
        return subclass
    
    def axes(self):
        return [axis.axisID for axis in self._axes]
        
    def axis(self, axisID):
        iAxis = self._normAxis(axisID)
        return self._axes[iAxis]
    pass
NDContainer = NDContainer_v0

#TODO: maybe we can put these to separate file so that NDContainer is self-contained and DD stuff is separate
from .DDconvertors import CCDHeightConvertor, SpectralConvertor, M112Convertor, M406Convertor
import io

#TODO: add CCD ampltide convertors (pass to NDContainer.__init__)
@NDContainer.registerExtension
class CCDFrameContainer(NDContainer):
    extension = ".NDframe"
    def __init__(self, roi, spatialConvertorSetup, spectralConvertorSetup, data=None, dtype=float):
        """
        [spatial|spectra]ConvertorSetup should be one of
          - filepath to stored convertor
          - file-like object filled with stored convertor
        """
        super().__init__("intensity / counts")
        self._setup(roi, spatialConvertorSetup, spectralConvertorSetup, data=data, dtype=dtype)
        pass
    
    def _setup(self, roi, spatialConvertorSetup, spectralConvertorSetup, data=None, dtype=float):
        self._roi = roi
        spatialConvertor = CCDHeightConvertor()
        spatialConvertor.load(spatialConvertorSetup, True)
        self.addAxis(CCDAxis(roi[3:], spatialConvertor, "spa", "row"))
        spectralConvertor = SpectralConvertor()
        spectralConvertor.load(spectralConvertorSetup, True)
        self.addAxis(CCDAxis(roi[:3], spectralConvertor, "spe", "wavelength"))
        if data is None:
            self.reserveDataSpace(dtype=dtype)
        else:
            self.setData(data)
        pass

    @classmethod
    def fromFile(cls, filename):
        print("CCDFrameContainer.fromFile")
        with open(filename, "rb") as f:
            s = pickle.load(f)
        
        spectralIO = io.BytesIO(s["spectral"])
        spatialIO = io.BytesIO(s["spatial"])
        
        return cls(s["roi"], spatialIO, spectralIO, s["data"])
    
    def save(self, filename):
        filename = self._sanitizeFilename(filename)
        s = {}
        s["data"] = self._data
        s["roi"] = self._roi
        s["spectral"] = self.axis("spe").convertor.dumpConfig()
        s["spatial"] = self.axis("spa").convertor.dumpConfig()
        
        with open(filename, "wb") as f:
            pickle.dump(s, f)
        pass
    pass

@NDContainer.registerExtension
class CCDMultiframeContainer(NDContainer):
    extension = ".NDmframe"
    def __init__(self, count, roi, spatialConvertorSetup, spectralConvertorSetup, data=None, dtype=float):
        """
        [spatial|spectra]ConvertorSetup should be one of
          - filepath to stored convertor
          - file-like object filled with stored convertor
        """
        super().__init__("intensity / counts")
        self._setup(count, roi, spatialConvertorSetup, spectralConvertorSetup, data=data, dtype=dtype)
        pass
    
    def _setup(self, count, roi, spatialConvertorSetup, spectralConvertorSetup, data=None, dtype=float):
        self._roi = roi
        self._count = count
        self.addAxis(RangeAxis(count, "N", "frames"))
        spatialConvertor = CCDHeightConvertor()
        spatialConvertor.load(spatialConvertorSetup, True)
        self.addAxis(CCDAxis(roi[3:], spatialConvertor, "spa", "row"))
        spectralConvertor = SpectralConvertor()
        spectralConvertor.load(spectralConvertorSetup, True)
        self.addAxis(CCDAxis(roi[:3], spectralConvertor, "spe", "wavelength"))
        if data is None:
            self.reserveDataSpace(dtype=dtype)
        else:
            self.setData(data)
        pass

    @classmethod
    def fromFile(cls, filename):
        with open(filename, "rb") as f:
            s = pickle.load(f)
        
        spectralIO = io.BytesIO(s["spectral"])
        spatialIO = io.BytesIO(s["spatial"])
        
        return cls(s["count"], s["roi"], spatialIO, spectralIO, s["data"])
    
    def save(self, filename):
        filename = self._sanitizeFilename(filename)
        s = {}
        s["data"] = self._data
        s["count"] = self._count
        s["roi"] = self._roi
        s["spectral"] = self.axis("spe").convertor.dumpConfig()
        s["spatial"] = self.axis("spa").convertor.dumpConfig()
        
        with open(filename, "wb") as f:
            pickle.dump(s, f)
        pass
    pass

#replacement of CCDScanData
@NDContainer.registerExtension
class CCDScanContainer(NDContainer):
    extension = ".NDscan"
    def __init__(self, roi, scanAxis1, scanAxis2, scanAxis3, spatialConvertorSetup, spectralConvertorSetup, scanConvertor1Setup, scanConvertor2Setup, scanConvertor3Setup, data=None, dtype=float):
        """
        [spatial|spectra]ConvertorSetup should be one of
          - filepath to stored convertor
          - file-like object filled with stored convertor
        """
        super().__init__("intensity / counts")
        print("CCDScanContainer.__init__", scanAxis1, scanAxis2, scanAxis3)
        self._setup(roi, scanAxis1, scanAxis2, scanAxis3, spatialConvertorSetup, spectralConvertorSetup, scanConvertor1Setup, scanConvertor2Setup, scanConvertor3Setup, data=data, dtype=dtype)
        print("\t", self.axis("pop").length())
        pass
    
    def _setup(self, roi, scanAxis1, scanAxis2, scanAxis3, spatialConvertorSetup, spectralConvertorSetup, scanConvertor1Setup, scanConvertor2Setup, scanConvertor3Setup, data=None, dtype=float):
        self._roi = roi
        self._scanAxes = scanAxis1, scanAxis2, scanAxis3
        scanConvertor1 = M112Convertor()
        scanConvertor1.load(scanConvertor1Setup, True)
        scanConvertor2 = M112Convertor()
        scanConvertor2.load(scanConvertor2Setup, True)
        scanConvertor3 = M406Convertor()
        scanConvertor3.load(scanConvertor3Setup, True)
        self.addAxis(ScanAxis(scanAxis1, scanAxis2, scanAxis3, scanConvertor1, scanConvertor2, scanConvertor3, "pop", ("Population 1", "Population 2", "Population 3")))
        
        spatialConvertor = CCDHeightConvertor()
        spatialConvertor.load(spatialConvertorSetup, True)
        self.addAxis(CCDAxis(roi[3:], spatialConvertor, "spa", "row"))
        
        spectralConvertor = SpectralConvertor()
        spectralConvertor.load(spectralConvertorSetup, True)
        self.addAxis(CCDAxis(roi[:3], spectralConvertor, "spe", "wavelength"))
        
        if data is None:
            self.reserveDataSpace(dtype=dtype)
        else:
            self.setData(data)
        pass
        
    @classmethod
    def fromFile(cls, filename):
        with open(filename, "rb") as f:
            s = pickle.load(f)
        
        spectralIO = io.BytesIO(s["spectral"])
        spatialIO = io.BytesIO(s["spatial"])
        
        scans = s["scan"]
        scan1IO = io.BytesIO(scans[0])
        scan2IO = io.BytesIO(scans[1])
        scan3IO = io.BytesIO(scans[2])
        
        sA = s["scanAxes"]
        return cls(s["roi"], sA[0], sA[1], sA[2], spatialIO, spectralIO, scan1IO, scan2IO, scan3IO, s["data"])
    
    def save(self, filename):
        filename = self._sanitizeFilename(filename)
        s = {} #TODO: this could be persistent dict, always saved to a file
        s["data"] = self._data
        s["roi"] = self._roi
        s["scanAxes"] = self._scanAxes
        s["spectral"] = self.axis("spe").convertor.dumpConfig()
        s["spatial"] = self.axis("spa").convertor.dumpConfig()
        s["scan"] = [it.dumpConfig() for it in self.axis("pop").convertor.convertors()] #ScanAxis uses CombinationConvertor
        
        with open(filename, "wb") as f:
            pickle.dump(s, f)
        pass
    pass

#TODO: possibly mapping from file extension to correct class?
