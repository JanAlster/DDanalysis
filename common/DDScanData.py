
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

#just for testing
import os
import os.path as op

if __name__=="__main__":
    #create simple app containing just the dialog - preparation
    import sys
    sys.path.append(op.join(sys.path[0], "..")) #need to import from directory one above this file
    os.chdir("..")
    print(sys.path)


#the rest

import numpy as np
from common.DDconvertors import CCDAmplitudeConvertor, SpectralConvertor
import pickle


"""
DD scan data container

we have basically 4D data - coherence axis, spectral axis, population axis and 1D CCD amplitude

spectral axis is in px and can be converted to nm, 1/cm, rad/fs, etc. if we have spectral calibration of CCD
amplitude of CCD is in counts and can be converted to photons, W, etc if we have intensity calibration of CCD (we will not have it at first)
coherence axis is in fs and cannot be converted (rather DDScan cannot be made without proper calibration available, DS1 and DS2)
population axis is in fs and cannot be converted (rather DDScan cannot be made without proper calibration available, DS3)
    
this will be hardcoded for this use, not flexible!
"""
import zipfile
import io 

class DDScanData(object):
    """
    data = np.array([population, coherence, spectral])
    
    how is spectral calibration aligned with pixels - centers or edges?
     - 0 is center of first pixel, etc.
     - center of last pixel is N-1 (N = data.shape[2])
     - since pixels have constant width, pixel edges are at [-0.5, 0.5, 1.5, ..., N-1.5, N-0.5]
     
    save calibration files and make own convertors - even if app convertor changes, data will keep the original calibration used to collect these data
     - in case you mess up the calibration, there will be no way how to change it later!
     - Todo: might do an external tool for that 
     - convertors are set and will not change!
     - Todo? might even store convertor parameters inside own datafile instead of links to calibration files ... it would be self-contained, but possibly wasting space and not in line with Convertor API
     
    allow online saving of data with each setDataStep if filepath is provided in setup 
     
    we also have to handle unit conversion of data
     - CCD measures value proportional to number of photons striking a pixel - integral over spectral and height range
     - what we want to show is more likely spectral density, but the units can change
     - CCD frame axes therefore have one more value than data (they represent edges of pixels, data the integral value inside)
     - we have to recalculate the value every time we want to export
     
    """
    
    def setup(self, data, ROI, coherenceAxisSetup, populationAxis, spectralCalibration, amplitudeCalibration, filepath=None, note=""):
        #CCD frame axes
        self._ROI = ROI #so it can be saved exactly
        wN = ROI[1]/ROI[2]
        self._wEdges = ROI[0] - 0.5 + np.arange(wN+1)
        self._wCenters = ROI[0] + np.arange(wN)

        self._coherenceAxisSetup = coherenceAxisSetup #From, To, Step
        self._coherenceAxis = np.arange(coherenceAxisSetup[0], coherenceAxisSetup[1]+0.5*coherenceAxisSetup[2], coherenceAxisSetup[2]) #overshoot To to prevent rounding off errors
        print("DDScanData.setup: coherence axis:", self._coherenceAxis)
        
        self._populationAxis = populationAxis
        
        #data
        # - assumption is that data are either set directly in setup (and not changed afterwards) or set frame-by-frame later using setDataStep
        if data is None:
            #initialize empty data arrray of correct size
            self._data = np.zeros((len(populationAxis), len(self._coherenceAxis), wN))
            self._dataSet = np.zeros(len(populationAxis), dtype=np.bool_)
        else:
            self._data = data
            self._dataSet = np.ones(len(populationAxis), dtype=np.bool_)
                        
        #CCD convertors
        self._spectralConvertor = SpectralConvertor()
        self._spectralConvertor.load(spectralCalibration)
        
        self._amplitudeConvertor = CCDAmplitudeConvertor()
        self._amplitudeConvertor.load(amplitudeCalibration)
        
        self.convertors = {"Spectral":self._spectralConvertor, "Amplitude":self._amplitudeConvertor}
        
        self._filepath = filepath
        #TODO: check the filepath for correct extension
        self._note = note
        
        if filepath is not None:
            self.save()
        pass
        
    def setDataStep(self, step, frame):
        #in case you do not setup data before (using DDScanData as a container during scan acquisition)
        self._data[step] = frame
        self._dataSet[step] = True
        if self._filepath is not None:
            #looks like there is no "easy" way to incrementally save ndarray to file
            # - we can write the whole thing each time - but it can be probably large and that would slow things down
            # - or we can do another file format (not pickle)
            # - possibly npz, which is supposed to be just a zip archive - easy (?) to insert new frames as separate npy files into the archive
            # - this does not keep data as one block when loading (probably could live with that)
            # - but we also need to include "header" pickle with the rest of the things to save
            # - all this is mostly to prevent data loss in case of crash 
            # - as a bonus we will avoid saving all the zeros
            # try to use zipfile
            # - assumption is that filepath is already saved (header, convertor calibration files, etc.) and we only fill in the data frames
            #with zipfile.ZipFile(self._filepath, "a") as z:
            #    z.writestr("f"+str(int(step))+".npy", frame.tobytes())
            #    pass
            self._saveFrame(step)
            pass
        pass
    
    #save
    #TODO: save/load filepaths could be relative to session (depends on use case though)
    # rather relative to filename
    
    def _saveFrame(self, i):
        with zipfile.ZipFile(self._filepath, "a") as z:
            #z.writestr("f"+str(int(i))+".npy", self._data[i].astype(np.float).tobytes())
            b = io.BytesIO()
            np.save(b, self._data[i])
            z.writestr("f"+str(int(i))+".npy", b.getvalue())
            pass
        
    def save(self, filepath=None):
        if filepath is None:
            filepath = self._filepath
        
        if filepath is None:
            #TODO: complain
            print("DDScanData.save: filepath is not specified")
            return
        
        #use zipfile
        # - save a (human readable) pickle header 
        # - save each scan step frame as npy
        # - also copy convertor calibration files (should be small enough not to matter and the zip will be self-contained)
        # - do not use compression
        # - possibly save note as ZipFile.comment
        with zipfile.ZipFile(filepath, "w") as z:
            #header
            
            h = {}
            h["spectral"] = self._ROI
            h["populations"] = self._populationAxis
            h["coherence"] = self._coherenceAxisSetup
            h["calSpe"] = self._spectralConvertor.filepath()
            h["calAmp"] = self._amplitudeConvertor.filepath()
            
            z.writestr("header.pickle", pickle.dumps(h, 0))
            
            #convertor calibration files
            z.write(h["calSpe"], "calSpe.pickle")
            z.write(h["calAmp"], "calAmp.pickle")
            
            #comment
            comment = self._note.encode("utf-8", "replace")
            z.comment = comment
            #in case comment is not displayed by file manager
            #TODO: I feel that there should be a difference between those two, but there is no apparent one on my test case
            z.writestr("note.txt", self._note)
            z.writestr("comment.txt", comment)
            
            #data
            # - is empty?
            # - here we are overwriting the file anyway
            # - but it is not necessary to write out non-initialized data
            for i in range(len(self._data)):
                #if self._dataSet[i]:
                #    z.writestr("f"+str(int(i))+".npy", self._data[i].tobytes())
                #pass
                if self._dataSet[i]: self._saveFrame(i)
            pass
        pass
        
    def load(self, filename):
        with zipfile.ZipFile(filename, "r") as z:
            with z.open("header.pickle", "r") as hf:
                h = pickle.load(hf) #TODO: possibly pickle.loads(z.read("header.pickle"))
                pass    

            print("DDScanData.load: header", h)
            
            #data - agregate f#.npy files
            names = z.namelist()
            print("DDScanData.load: names", names)
            data = []
            for i in range(len(h["populations"])):
                name = "f"+str(int(i))+".npy"
                if name in names:
                    b = io.BytesIO(z.read(name))
                    data.append(np.load(b))
                else:
                    #TODO: complain
                    print("DDScanData.load: cannot find frame", i, "in scan archive", filename)
                    break
                pass   
            data = np.array(data)                 
                        
            
            #use original calibration files here - convertors likely cannot live on top of files inside zip archive (we will need constant convertor that loads parameters file and does not need it anymore (also cannot change afterwards)
                        
            #TODO: solve save/load todo above and delete this
            #need to hack this to be able to open moved session
            #assuming calibration files at the same place as filename (later could be different)
            import ntpath
            import os.path as op
            dirname = op.dirname(filename)
            h["calSpe"] = op.join(dirname, ntpath.basename(h["calSpe"]))
            h["calAmp"] = op.join(dirname, ntpath.basename(h["calAmp"]))
            
            #note
            note = str(z.comment, "utf-8", "replace")
            
            self.setup(data, h["spectral"], h["coherence"], h["populations"], h["calSpe"], h["calAmp"], note=note)
                
            pass
        pass
        
    #cuts
    #TODO: handle convertor flipping axis (from ascending to descending) - divide by diff will even change sign of data!
    # - return cut position, cut data (as density), vertical axis, horizontal axis
    #TODO: I want to use the same display widget as for CCDScanData, but the interface uses different names (and arguments) than here
    def cutPopulation(self, i, spectralUnit=None):
        #spectral
        cut = self._data[i]
        
        wEdges = self._spectralConvertor.V2U(self._wEdges, spectralUnit)
        dwEdges = np.diff(wEdges)
        wCenters = self._spectralConvertor.V2U(self._wCenters, spectralUnit)
        
        cut = cut/dwEdges
        
        #amplitude unit is amplitudeUnit/spectralUnit (originally is also /spatialUnit, but that should be integrated over)
        return self._populationAxis[i], cut, self._coherenceAxis, wCenters
        
    def cutSpectral(self, i, spectralUnit=None):
        #spatial/scan cut
        cut = self._data[:,:,i]

        #single point spectrum
        wEdges = self._spectralConvertor.V2U(self._wEdges[i:i+2], spectralUnit)
        dwEdges = np.diff(wEdges)
        wCenters = self._spectralConvertor.V2U(self._wCenters[i], spectralUnit)

        cut = cut / dwEdges[0]
        
        #amplitude unit is amplitudeUnit/spectralUnit
        return wCenters, cut, self._populationAxis, self._coherenceAxis
        
    def cutCoherence(self, i, spectralUnit=None):
        #spectral/scan cut
        cut = self._data[:, i]
        
        wEdges = self._spectralConvertor.V2U(self._wEdges, spectralUnit)
        dwEdges = np.diff(wEdges)
        wCenters = self._spectralConvertor.V2U(self._wCenters, spectralUnit)

        cut = cut / dwEdges

        if dwEdges.mean() < 0:
            wCenters = wCenters[::-1]
            cut = -cut[:,::-1]
        
        #amplitude unit is amplitudeUnit/spectralUnit
        return self._coherenceAxis[i], cut, self._populationAxis, wCenters

    def export(self, spectralUnit=None, background=None):
        wEdges = self._spectralConvertor.V2U(self._wEdges, spectralUnit)
        dwEdges = np.diff(wEdges)
        wCenters = self._spectralConvertor.V2U(self._wCenters, spectralUnit)
        
        if background is not None:
            if background == "auto":
                h, b = np.histogram(self._data, range(590, 610))
                background = b[h.argmax()]
                print("DDScanData.export applying correction for automatically detected background", background)
            else:
                print("DDScanData.export applying correction for background", background)
            data = self._data - background
        else:
            data = np.array(self._data) #copy
        
        data /= dwEdges
        
        if dwEdges.mean() < 0:
            wCenters = wCenters[::-1]
            data = -data[:,:,::-1]

        return data, self._populationAxis, self._coherenceAxis, wCenters
    
    def shape(self):
        return self._data.shape
        
    pass

#test    
if __name__=="__main__":
    d = DDScanData()
    pop = np.arange(10)
    coh = [-50, 50. ,1]
    d.setup(None, [0, 1600, 1], coh, pop, "/home/araigne/Code/DD_v02/data/testSessionStyryl9/calspe_680nm_Xe_B4_2016-09-02_12-09-38.pickle", "/home/araigne/Code/DD_v02/data/testSessionStyryl9/defaultAmplitudeCalibration.pickle", "test.DDscan", "Ahoj Fíbí")
    d.setDataStep(0, np.ones((101,1600)))
    
    #, ROI, coherenceAxisSetup, populationAxis, spectralCalibration, amplitudeCalibration, filepath=None, note=""):
    
    b = DDScanData()
    b.load("test.DDscan")
    print(b, b._note, b._coherenceAxis)
    print(b._data)
    
