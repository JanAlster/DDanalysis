# ***********************************************************************
# DD package, data collection and analysis of 2D electronic spectra
# Copyright (C) 2016, 2017  Jan Alster (Charles Univesity, Prague)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
# ***********************************************************************

# just for testing
import os
import os.path as op
from typing import Type

if __name__ == "__main__":
    # create simple app containing just the dialog - preparation
    import sys

    sys.path.append(op.join(sys.path[0], ".."))  # need to import from directory one above this file
    os.chdir("..")
    print(sys.path)

# TODO: IMPORTANT - how is frame naming by index compatible with random order population time measurement?

# TODO: save chopper filter reference and DS convertors reference too (in case of "Save full dataset" option)

# TODO: IMPORTANT - save here open zipfile as "w" which might destroy already stored data
# - it is not supposed to access already created archive, but there ware some confusion that might have allowed that
# - this should be checked (possibly change the mode to "a")

# TODO: overwriting file leads to duplicate entries in zip file

# TODO: seems like we need to store convertors inside the datafiles, I too often move the data out of session folder

# the rest

import numpy as np

try:  # never sure which version will be needed (this is invoked from all kinds of scripts in all kinds of places)
    from common.DDconvertors import CCDAmplitudeConvertor, SpectralConvertor, Convertor
except:
    from .DDconvertors import CCDAmplitudeConvertor, SpectralConvertor, Convertor
import pickle

"""
DD scan data container

we have basically 4D data - coherence axis, spectral axis, population axis and 1D CCD amplitude

spectral axis is in px and can be converted to nm, 1/cm, rad/fs, etc. if we have spectral calibration of CCD
amplitude of CCD is in counts and can be converted to photons, W, etc if we have intensity calibration of CCD (we will not have it at first)
coherence axis is in fs and cannot be converted (rather DDScan cannot be made without proper calibration available, DS1 and DS2)
population axis is in fs and cannot be converted (rather DDScan cannot be made without proper calibration available, DS3)
    
this will be hardcoded for this use, not flexible!

v2 is able to store multiple datasets (signal, LO, scatterings, ...)
"""
import zipfile
import io


# if we just want to look
def readHeader(filename):
    with zipfile.ZipFile(filename, "r") as z:
        z.debug = 3
        # ~ print(z.namelist())
        # ~ print(filename)
        try:
            with z.open("header.pickle", "r") as hf:
                h = pickle.load(hf)  # TODO: possibly pickle.loads(z.read("header.pickle"))
                pass
        except:
            h = pickle.loads(z.read("header.pickle"))

        # ~ print("DDScanDataRich.load: header", h)
    return h


class DDScanDataRich(object):
    """
    data = {"dataset": np.array([population, coherence, spectral]), ... }
    
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

    def setup(self, data, ROI, coherenceAxisSetup, populationAxis, spectralCalibration, amplitudeCalibration,
              filepath=None, note="", dataKeys=None, defaultSet="S123"):
        # CCD frame axes
        self._ROI = ROI  # so it can be saved exactly
        wN = int(round(ROI[1] / ROI[2], 0))
        # todo: TEST under assumption that line calibration is done on full ROI
        self._wEdges = ROI[0] - 0.5 + np.arange(wN + 1) * ROI[2]
        self._wCenters = ROI[0] -0.5 + ROI[2]/2 + np.arange(wN) * ROI[2]
        # self._wEdges = ROI[0] - 0.5 + np.arange(wN + 1)
        # self._wCenters = ROI[0] + np.arange(wN)

        self._coherenceAxisSetup = coherenceAxisSetup  # From, To, Step
        self._coherenceAxis = np.arange(coherenceAxisSetup[0], coherenceAxisSetup[1] + 0.5 * coherenceAxisSetup[2],
                                        coherenceAxisSetup[2])  # overshoot To to prevent rounding off errors
        print("DDScanDataRich.setup: coherence axis:", self._coherenceAxis)

        self._populationAxis = populationAxis

        # data
        # - assumption is that data are either set directly in setup (and not changed afterwards) or set frame-by-frame later using setDataStep
        if data is None:
            # initialize empty data arrray of correct size
            print(wN, dataKeys)
            self._data = {it: np.zeros((len(populationAxis), len(self._coherenceAxis), wN)) for it in dataKeys}
            self._dataSet = np.zeros(len(populationAxis),
                                     dtype=np.bool_)  # mark which (population time) frames are set (contain data)
        else:
            self._data = data
            self._dataSet = np.ones(len(populationAxis), dtype=np.bool_)
        self._defaultSet = defaultSet

        # CCD convertors
        if isinstance(spectralCalibration, Convertor):
            self._spectralConvertor = spectralCalibration
        else:
            self._spectralConvertor = SpectralConvertor()
            self._spectralConvertor.load(spectralCalibration, memory=True)

        if isinstance(amplitudeCalibration, Convertor):
            self._amplitudeConvertor = amplitudeCalibration
        else:
            self._amplitudeConvertor = CCDAmplitudeConvertor()
            self._amplitudeConvertor.load(amplitudeCalibration, memory=True)

        # test:
        # print("original", self._spectralConvertor.V2U(ROI[0] + np.arange(wN), "rad/fs")*5308)
        # print("binned", self._spectralConvertor.V2U(self._wCenters, "rad/fs")*5308)

        self.convertors = {"Spectral": self._spectralConvertor, "Amplitude": self._amplitudeConvertor}

        self._filepath = filepath
        # TODO: check the filepath for correct extension
        self._note = note

        # subs tep temporal storage
        self._subStep = None
        self._subSet = None
        self._subData = None
        self._wN = wN

        if filepath is not None:
            self.save()

        return self._coherenceAxis

    # try not to abuse these
    # - they expect orderly flow of data, keep the coherence axis, etc..
    def setDataSubStep(self, step, subStep, data):
        if self._subStep is not None and step != self._subStep:
            # complain
            print(
                "DDScanDataRich.setDataSubStep: starting new sub step storage with unfinished old one (will be thrown away)")
            self._subStep = None  # throw away old substep

        if self._subStep is None:
            # starting new substep storage
            self._subStep = step
            self._subSet = np.zeros(len(self._coherenceAxis), dtype=np.bool_)
            self._subData = {it: np.empty((len(self._coherenceAxis), self._wN)) for it in self._data}

        # save sub step
        self._subSet[subStep] = True
        for it in data:
            if it in self._subData: self._subData[it][subStep] = data[it]

        # note that this will consider the sub step finished as soon as all sub steps are filled disregarding which datasets are filled (or even if any dataset is filled completely)
        # so please fill all datasets at once for every substep
        # this is not trying to be foolproof
        if self._subSet.all():
            self.setDataStep(self._subStep, self._subData)
            self._subStep = None
            self._subSet = None
            self._subData = None
        pass

    def setDataStep(self, step, frame):
        # in case you do not setup data before (using DDScanData as a container during scan acquisition)
        # self._data[step] = frame
        for it in frame:
            if it in self._data: self._data[it][step] = frame[it]
        self._dataSet[step] = True
        if self._filepath is not None:
            # looks like there is no "easy" way to incrementally save ndarray to file
            # - we can write the whole thing each time - but it can be probably large and that would slow things down
            # - or we can do another file format (not pickle)
            # - possibly npz, which is supposed to be just a zip archive - easy (?) to insert new frames as separate npy files into the archive
            # - this does not keep data as one block when loading (probably could live with that)
            # - but we also need to include "header" pickle with the rest of the things to save
            # - all this is mostly to prevent data loss in case of crash 
            # - as a bonus we will avoid saving all the zeros
            # try to use zipfile
            # - assumption is that filepath is already saved (header, convertor calibration files, etc.) and we only fill in the data frames
            # with zipfile.ZipFile(self._filepath, "a") as z:
            #    z.writestr("f"+str(int(step))+".npy", frame.tobytes())
            #    pass
            with zipfile.ZipFile(self._filepath, "a") as z:
                self._saveFrame(step, z)
            pass
        pass

    def datasets(self):
        return list(self._data.keys())

    # save
    # TODO: save/load filepaths could be relative to session (depends on use case though)
    # rather relative to filename

    def _saveFrame(self, i, z):
        print("DDScanDataRich._saveFrame", i, self._filepath)
        for it in self._data:
            b = io.BytesIO()
            np.save(b, self._data[it][i])
            z.writestr(str(it) + "_" + str(int(i)) + ".npy", b.getvalue())
        pass

    def save(self, filepath=None):
        if filepath is None:
            filepath = self._filepath

        if filepath is None:
            # TODO: complain
            print("DDScanDataRich.save: filepath is not specified")
            return

        # use zipfile
        # - save a (human readable) pickle header 
        # - save each scan step frame as npy
        # - also copy convertor calibration files (should be small enough not to matter and the zip will be self-contained)
        # - do not use compression
        # - possibly save note as ZipFile.comment
        with zipfile.ZipFile(filepath, "a") as z:  # TODO: check the code and select proper mode here "w"/"a"
            # header

            h = {}
            h["spectral"] = self._ROI
            h["populations"] = self._populationAxis
            h["coherence"] = self._coherenceAxisSetup
            h["calSpe"] = self._spectralConvertor.filepath()
            h["calAmp"] = self._amplitudeConvertor.filepath()
            h["dataSets"] = list(self._data.keys())
            h["defaultSet"] = self._defaultSet

            z.writestr("header.pickle", pickle.dumps(h, 0))

            # convertor calibration files
            # TODO: this is ugly, we need this for constant convertors (remastering) but it is not elegant at all (the convertor system got out of hand now)
            z.writestr("calSpe.pickle", self._spectralConvertor._saveHack())
            z.writestr("calAmp.pickle", self._amplitudeConvertor._saveHack())
            # z.write(h["calSpe"], "calSpe.pickle")
            # z.write(h["calAmp"], "calAmp.pickle")

            # comment
            comment = self._note.encode("utf-8", "replace")
            z.comment = comment
            # in case comment is not displayed by file manager
            # TODO: I feel that there should be a difference between those two, but there is no apparent one on my test case
            z.writestr("note.txt", self._note)
            z.writestr("comment.txt", comment)

            # data
            # - is empty?
            # - here we are overwriting the file anyway
            # - but it is not necessary to write out non-initialized data
            for i in range(len(self._dataSet)):
                # if self._dataSet[i]:
                #    z.writestr("f"+str(int(i))+".npy", self._data[i].tobytes())
                # pass
                if self._dataSet[i]: self._saveFrame(i, z)
            pass
        pass

    def addNote(self, note):
        self._note += note
        if self._filepath is not None:
            with zipfile.ZipFile(self._filepath, "a") as z:
                comment = self._note.encode("utf-8", "replace")
                z.comment = comment
        pass

    def setNote(self, note):
        self._note = note
        if self._filepath is not None:
            with zipfile.ZipFile(self._filepath, "a") as z:
                comment = self._note.encode("utf-8", "replace")
                z.comment = comment
        pass


    def load(self, filename, frames=None, datasets=None):
        with zipfile.ZipFile(filename, "r") as z:
            z.debug = 3
            print(z.namelist())
            print(filename)
            try:
                with z.open("header.pickle", "r") as hf:
                    h = pickle.load(hf)  # TODO: possibly pickle.loads(z.read("header.pickle"))
                    pass
            except:
                h = pickle.loads(z.read("header.pickle"))

            print("DDScanDataRich.load: header", h)

            # data - agregate f#.npy files
            names = z.namelist()
            print("DDScanDataRich.load: names", names)
            data = {}
            loaded = []  # which frames are loaded, TODO: this will fail horribly if not all datasets are available for the frame
            framesToLoad = range(len(h["populations"])) if frames is None else frames
            if datasets is None: datasets = h["dataSets"]
            for key in datasets:
                dataset = []
                for i in framesToLoad:
                    name = str(key) + "_" + str(int(i)) + ".npy"
                    if name in names:
                        b = io.BytesIO(z.read(name))
                        dataset.append(np.load(b))
                        if i not in loaded: loaded.append(i)
                    else:
                        # TODO: complain
                        print("DDScanDataRich.load: cannot find frame", key, i, "in scan archive", filename)
                        # ~ break
                    pass
                data[key] = np.array(dataset)

            populations = h["populations"][loaded]

            # use original calibration files here - convertors likely cannot live on top of files inside zip archive (we will need constant convertor that loads parameters file and does not need it anymore (also cannot change afterwards)

            try:
                _spectralConvertor = SpectralConvertor()
                b = io.BytesIO(z.read("calSpe.pickle"))
                _spectralConvertor.load(b)
            except:
                _spectralConvertor = h["calSpe"]

            try:
                _amplitudeConvertor = CCDAmplitudeConvertor()
                b = io.BytesIO(z.read("calAmp.pickle"))
                _amplitudeConvertor.load(b)
            except:
                _amplitudeConvertor = h["calAmp"]

            # TODO: solve save/load todo above and delete this
            # need to hack this to be able to open moved session
            # assuming calibration files at the same place as filename (later could be different)
            # import ntpath
            # import os.path as op
            # dirname = op.dirname(filename)
            # h["calSpe"] = op.join(dirname, ntpath.basename(h["calSpe"]))
            # h["calAmp"] = op.join(dirname, ntpath.basename(h["calAmp"]))

            # note
            note = str(z.comment, "utf-8", "replace")

            # compatibility with old data
            if "defaultSet" not in h: h["defaultSet"] = "S123"

            self.setup(data, h["spectral"], h["coherence"], populations, _spectralConvertor, _amplitudeConvertor, None,
                       note=note, defaultSet=h["defaultSet"])  # filepath=None to prevent saving in setup

            pass
        pass

    # cuts
    # TODO: handle convertor flipping axis (from ascending to descending) - divide by diff will even change sign of data!
    # - return cut position, cut data (as density), vertical axis, horizontal axis
    # TODO: I want to use the same display widget as for CCDScanData, but the interface uses different names (and arguments) than here
    def cutPopulation(self, i, spectralUnit=None, dataset=None):
        # spectral
        if dataset is None: dataset = self._defaultSet
        cut = self._data[dataset][i]

        wEdges = self._spectralConvertor.V2U(self._wEdges, spectralUnit)
        dwEdges = np.diff(wEdges)
        wCenters = self._spectralConvertor.V2U(self._wCenters, spectralUnit)

        cut = cut / dwEdges

        # amplitude unit is amplitudeUnit/spectralUnit (originally is also /spatialUnit, but that should be integrated over)
        return self._populationAxis[i], cut, self._coherenceAxis, wCenters

    def cutSpectral(self, i, spectralUnit=None, dataset=None):
        # spatial/scan cut
        if dataset is None: dataset = self._defaultSet
        cut = self._data[dataset][:, :, i]

        # single point spectrum
        wEdges = self._spectralConvertor.V2U(self._wEdges[i:i + 2], spectralUnit)
        dwEdges = np.diff(wEdges)
        wCenters = self._spectralConvertor.V2U(self._wCenters[i], spectralUnit)

        cut = cut / dwEdges[0]

        # amplitude unit is amplitudeUnit/spectralUnit
        return wCenters, cut, self._populationAxis, self._coherenceAxis

    def cutCoherence(self, i, spectralUnit=None, dataset=None):
        # spectral/scan cut
        if dataset is None: dataset = self._defaultSet
        cut = self._data[dataset][:, i]

        wEdges = self._spectralConvertor.V2U(self._wEdges, spectralUnit)
        dwEdges = np.diff(wEdges)
        wCenters = self._spectralConvertor.V2U(self._wCenters, spectralUnit)

        cut = cut / dwEdges

        if dwEdges.mean() < 0:
            wCenters = wCenters[::-1]
            cut = -cut[:, ::-1]

        # amplitude unit is amplitudeUnit/spectralUnit
        return self._coherenceAxis[i], cut, self._populationAxis, wCenters

    def export(self, spectralUnit=None, background=None, dataset=None):
        wEdges = self._spectralConvertor.V2U(self._wEdges, spectralUnit)
        dwEdges = np.diff(wEdges)
        wCenters = self._spectralConvertor.V2U(self._wCenters, spectralUnit)

        if dataset is None: dataset = self._defaultSet

        if background is not None:
            if background == "auto":
                h, b = np.histogram(self._data[dataset], range(-610,
                                                               610))  # Todo: we need to go from 0 here since 1 is used for PP data - it shows that this autoprocedure will fail if background is not in expected range
                # range starting from -610 is an attempt to sidestep a bug causing background subtracted twice for some data (I do not want to remeasure the data, new data shoudl be fine)
                # IMPORTANT TODO: this is not good
                #   - for multiframe measurement the apparent background is not 600
                #   - for broad LO we might not be able to get a background estimate
                if np.all(h == 0):
                    print("WARNING! DDScanData.export cannot autodetect background")
                    background = 0
                else:
                    background = b[h.argmax()]
                    print("DDScanData.export applying correction for automatically detected background", background)
            else:
                print("DDScanData.export applying correction for background", background)
            data = self._data[dataset] - background
        else:
            data = np.array(self._data[dataset])  # copy

        data /= dwEdges

        if dwEdges.mean() < 0:
            wCenters = wCenters[::-1]
            data = -data[:, :, ::-1]

        return data, self._populationAxis, self._coherenceAxis, wCenters

    def shape(self):
        return (len(self._populationAxis), len(self._coherenceAxis), self._wN)

    pass


class Dataset:
    def __init__(self, shape, indices, data=None, dtype: Type = float):
        self.indices = indices  # read-only
        if data is None:
            # initialize empty data array of correct size
            self.data = np.zeros(shape, dtype=dtype)
            # mark which (population time) frames are set (contain data)
            self._dataSet = np.zeros(shape[0], dtype=bool)
        else:
            self.data = data  # reuse memory, but be careful
            self._dataSet = np.ones(shape[0], dtype=bool)

    def frame_is_set(self, i) -> bool:
        return self._dataSet[i]

    def set_frame(self, i, data):
        assert not self._dataSet[i]
        self.data[i] = data
        self._dataSet[i] = True


class DDScanDatasets:
    """
    Similar to DDScanDataRich, but with a variable size of individual datasets.

    Object itself will store the metadata (like convertors) for reuse by the datasets.

    File storage as a zipfile.

    Individual frames are stored as they come (which need not be in order of population axis).

    The easiest approach is probably to reuse CCDFrameContainer; there will be some duplication, but
     hopefully not too much.
    """

    def __init__(self, ROI, coherenceAxis, populationAxis, spectralCalibration, amplitudeCalibration,
                 filepath=None, note="", dataKeys=None, defaultSet="B24"):
        """

        :param data:
        :param ROI: x_start, x_width, x_step, y_start, y_height, y_step
        :param coherenceAxisSetup:
        :param populationAxis:
        :param spectralCalibration:
        :param amplitudeCalibration:
        :param filepath:
        :param note:
        :param dataKeys:
        :param defaultSet:
        """
        # CCD frame axes
        self._ROI = ROI  # so it can be saved exactly
        wN = int(round(ROI[1] / ROI[2], 0))
        # todo: TEST under assumption that line calibration is done on full ROI
        self._wEdges = ROI[0] - 0.5 + np.arange(wN + 1) * ROI[2]
        self._wCenters = ROI[0] - 0.5 + ROI[2] / 2 + np.arange(wN) * ROI[2]
        # self._wEdges = ROI[0] - 0.5 + np.arange(wN + 1)
        # self._wCenters = ROI[0] + np.arange(wN)

        self._coherenceAxis = coherenceAxis
        self._populationAxis = populationAxis

        # data
        self._datasets = {}

        # CCD convertors
        if isinstance(spectralCalibration, Convertor):
            self._spectralConvertor = spectralCalibration
        else:
            self._spectralConvertor = SpectralConvertor()
            self._spectralConvertor.load(spectralCalibration, memory=True)

        if isinstance(amplitudeCalibration, Convertor):
            self._amplitudeConvertor = amplitudeCalibration
        else:
            self._amplitudeConvertor = CCDAmplitudeConvertor()
            self._amplitudeConvertor.load(amplitudeCalibration, memory=True)

        self.convertors = {"Spectral": self._spectralConvertor, "Amplitude": self._amplitudeConvertor}

        self._filepath = filepath
        # TODO: check the filepath for correct extension
        self._note = note

        self._defaultSet = defaultSet

        # subs tep temporal storage
        self._recorded_population = None  # which frame from populationAxis is being recorder
        self._recorded_coherences = None  # which indices from coherenceAxis are already recorded
        self._sub_data = None  # recorded data

        if filepath is not None:
            self.save()

        pass

    def addDataset(self, key: str, indices, data=None, dtype: Type = float):
        """
        Will prepare and/or store data to a single dataset.
        A dataset is accessible under key.
        Each dataset has a unique number of spectra for each population time and coherence time pair
        of DDScanDatasets instance. This axis corresponds to the original indices of rows on the CCDreadout frame
        for this particular signal channel.

        :param key:
        :param indices:
        :param data:
        :param dtype: Datatype to reserve if data is None
        :return:
        """

        assert key not in self._datasets, "DDScanDatasets.addDataset: cannot overwrite an existing dataset"

        # a dataset is a (could be dict, it should not matter; could be a dataclass) indices, data
        self._datasets[key] = Dataset((len(self._populationAxis), len(self._coherenceAxis),
                                       len(indices), len(self._wCenters)),
                                      indices, data, dtype)

    def datasets(self):
        return list(self._datasets.keys())

    def raw(self, dataset):
        """
        Gives direct access to stored data. Modify at your own risk.

        :param dataset:
        :return:
        """
        dataset = self._datasets[dataset]
        return dataset.data, self._populationAxis, self._coherenceAxis, dataset.indices, self._wCenters

    def export(self, spectralUnit=None, background=None, dataset=None):
        wEdges = self._spectralConvertor.V2U(self._wEdges, spectralUnit)
        dwEdges = np.diff(wEdges)
        wCenters = self._spectralConvertor.V2U(self._wCenters, spectralUnit)

        if dataset is None: dataset = self._defaultSet

        if background is not None:
            if background == "auto":
                h, b = np.histogram(self._datasets[dataset].data, range(-610,
                                                               610))  # Todo: we need to go from 0 here since 1 is used for PP data - it shows that this autoprocedure will fail if background is not in expected range
                # range starting from -610 is an attempt to sidestep a bug causing background subtracted twice for some data (I do not want to remeasure the data, new data shoudl be fine)
                # IMPORTANT TODO: this is not good
                #   - for multiframe measurement the apparent background is not 600
                #   - for broad LO we might not be able to get a background estimate
                if np.all(h == 0):
                    print("WARNING! DDScanData.export cannot autodetect background")
                    background = 0
                else:
                    background = b[h.argmax()]
                    print("DDScanData.export applying correction for automatically detected background", background)
            else:
                print("DDScanData.export applying correction for background", background)
            data = self._datasets[dataset].data - background
        else:
            data = np.array(self._datasets[dataset])  # copy

        data /= dwEdges

        if dwEdges.mean() < 0:
            wCenters = wCenters[::-1]
            data = -data[:, :, :, ::-1]

        return data, self._populationAxis, self._coherenceAxis, self._datasets[dataset].indices, wCenters

    def save(self, filepath=None):
        if filepath is None:
            filepath = self._filepath

        if filepath is None:
            # TODO: complain
            print("DDScanDatasets.save: filepath is not specified")
            return

        # use zipfile
        # - save a (human readable) pickle header
        # - save each scan step frame as npy
        # - also copy convertor calibration files (should be small enough not to matter and the zip will be self-contained)
        # - do not use compression
        # - possibly save note as ZipFile.comment
        with zipfile.ZipFile(filepath, "a") as z:  # TODO: check the code and select proper mode here "w"/"a"
            # header

            h = {}
            h["spectral"] = self._ROI
            h["populations"] = self._populationAxis
            h["coherence"] = self._coherenceAxis  # TODO: this is not backward compatible
            h["calSpe"] = self._spectralConvertor.filepath()
            h["calAmp"] = self._amplitudeConvertor.filepath()
            h["dataSets"] = list(self._datasets.keys())
            h["defaultSet"] = self._defaultSet
            h["indices"] = {key: self._datasets[key].indices for key in self._datasets}

            z.writestr("header.pickle", pickle.dumps(h, 0))

            # convertor calibration files
            # TODO: this is ugly, we need this for constant convertors (remastering) but it is not elegant at all (the convertor system got out of hand now)
            z.writestr("calSpe.pickle", self._spectralConvertor._saveHack())
            z.writestr("calAmp.pickle", self._amplitudeConvertor._saveHack())
            # z.write(h["calSpe"], "calSpe.pickle")
            # z.write(h["calAmp"], "calAmp.pickle")

            # comment
            comment = self._note.encode("utf-8", "replace")
            z.comment = comment
            # in case comment is not displayed by file manager
            # TODO: I feel that there should be a difference between those two, but there is no apparent one on my test case
            z.writestr("note.txt", self._note)
            z.writestr("comment.txt", comment)

            # data
            # - is empty?
            # - here we are overwriting the file anyway
            # - but it is not necessary to write out non-initialized data

            for i in range(len(self._populationAxis)):
                self._saveFrame(i, z)
        pass

    def _saveFrame(self, i, z):
        print("DDScanDataset._saveFrame", i, self._filepath)
        for key in self._datasets:
            if not self._datasets[key].frame_is_set(i):
                return

        for key in self._datasets:
            b = io.BytesIO()
            np.save(b, self._datasets[key].data[i])
            z.writestr(str(key) + "_" + str(int(i)) + ".npy", b.getvalue())
        pass

    def shape(self, key):
        return len(self._populationAxis), len(self._coherenceAxis), len(self._datasets[key].indices), len(
            self._wCenters)

    @classmethod
    def from_file(cls, filename, frames=None, datasets=None):
        with zipfile.ZipFile(filename, "r") as z:
            z.debug = 3
            print(z.namelist())
            print(filename)
            try:
                with z.open("header.pickle", "r") as hf:
                    h = pickle.load(hf)  # TODO: possibly pickle.loads(z.read("header.pickle"))
                    pass
            except:
                h = pickle.loads(z.read("header.pickle"))

            print("DDScanDatasets.from_file: header", h)

            # data - aggregate f#.npy files
            names = z.namelist()

            data = {}

            if frames is None:
                frames = list(range(len(h["populations"])))

            if datasets is None:
                datasets = h["dataSets"]

            loaded_frames = {key: [] for key in datasets}

            for key in datasets:
                dataset_data = []
                for i in frames:
                    name = str(key) + "_" + str(int(i)) + ".npy"
                    if name in names:
                        b = io.BytesIO(z.read(name))
                        dataset_data.append(np.load(b))
                        loaded_frames[key].append(i)
                    else:
                        # TODO: complain
                        print("DDScanDataset.load: cannot find frame", key, i, "in scan archive", filename)
                        # ~ break
                    pass
                data[key] = dataset_data

            # only use frames which are loaded for all datasets
            loaded = []
            for it in frames:
                if all(it in loaded_frames[key] for key in datasets):
                    loaded.append(it)

            for key in datasets:
                data[key] = np.array([data[key][it] for it in loaded])

            populations = h["populations"][loaded]

            # use original calibration files here - convertors likely cannot live on top of files inside zip archive
            # (we will need constant convertor that loads parameters file and does not need it anymore
            # (also cannot change afterwards)

            try:
                _spectralConvertor = SpectralConvertor()
                b = io.BytesIO(z.read("calSpe.pickle"))
                _spectralConvertor.load(b)
            except:
                _spectralConvertor = h["calSpe"]

            try:
                _amplitudeConvertor = CCDAmplitudeConvertor()
                b = io.BytesIO(z.read("calAmp.pickle"))
                _amplitudeConvertor.load(b)
            except:
                _amplitudeConvertor = h["calAmp"]

            # TODO: solve save/load todo above and delete this
            # need to hack this to be able to open moved session
            # assuming calibration files at the same place as filename (later could be different)
            # import ntpath
            # import os.path as op
            # dirname = op.dirname(filename)
            # h["calSpe"] = op.join(dirname, ntpath.basename(h["calSpe"]))
            # h["calAmp"] = op.join(dirname, ntpath.basename(h["calAmp"]))

            # note
            note = str(z.comment, "utf-8", "replace")

        # compatibility with old data
        if "defaultSet" not in h: h["defaultSet"] = "B24"

        res = cls(h["spectral"], h["coherence"], populations, _spectralConvertor, _amplitudeConvertor, None,
                  note=note, defaultSet=h["defaultSet"])

        for key in data:
            res.addDataset(key, h["indices"][key], data[key])
        return res

    def addNote(self, note):
        self._note += note
        if self._filepath is not None:
            with zipfile.ZipFile(self._filepath, "a") as z:
                comment = self._note.encode("utf-8", "replace")
                z.comment = comment
        pass

    def setNote(self, note):
        self._note = note
        if self._filepath is not None:
            with zipfile.ZipFile(self._filepath, "a") as z:
                comment = self._note.encode("utf-8", "replace")
                z.comment = comment
        pass

    def setDataStep(self, step, frame):
        """
        Store aggregated data for a single frame to datasets (and filepath).

        :param step:
        :param frame:
        :return:
        """
        # in case you do not setup data before (using DDScanData as a container during scan acquisition)
        # self._data[step] = frame
        assert all(key in frame for key in self._datasets)

        for key in self._datasets:
            self._datasets[key].set_frame(step, frame[key])

        if self._filepath is not None:
            # looks like there is no "easy" way to incrementally save ndarray to file
            # - we can write the whole thing each time - but it can be probably large and that would slow things down
            # - or we can do another file format (not pickle)
            # - possibly npz, which is supposed to be just a zip archive - easy (?) to insert new frames as separate npy files into the archive
            # - this does not keep data as one block when loading (probably could live with that)
            # - but we also need to include "header" pickle with the rest of the things to save
            # - all this is mostly to prevent data loss in case of crash
            # - as a bonus we will avoid saving all the zeros
            # try to use zipfile
            # - assumption is that filepath is already saved (header, convertor calibration files, etc.) and we only fill in the data frames
            # with zipfile.ZipFile(self._filepath, "a") as z:
            #    z.writestr("f"+str(int(step))+".npy", frame.tobytes())
            #    pass
            with zipfile.ZipFile(self._filepath, "a") as z:
                self._saveFrame(step, z)
            pass
        pass

    # try not to abuse these
    # - they expect orderly flow of data, keep the coherence axis, etc..
    def setDataSubStep(self, population_step, coherence_step, data):
        """
        Start recording a new frame, allocate space for data for individual datasets.
        We could use the space allocated in Datasets, but the intention is to move the data only after the whole
        frame is finished to preserve data coherence.

        :param step: - population axis index
        :param subStep: - coherence axis index
        :param data: {key:(dataset.indices) x (spectral) for key in datasets}
        :return: None
        """

        assert all(key in data for key in self._datasets)

        if self._recorded_population is not None and population_step != self._recorded_population:
            # complain
            print(
                "DDScanDatasets.setDataSubStep: starting new sub step storage with unfinished old one (will be thrown away)")
            self._recorded_population = None  # throw away old substep

        if self._recorded_population is None:
            # starting new substep storage
            self._recorded_population = population_step
            self._recorded_coherences = np.zeros(len(self._coherenceAxis), dtype=np.bool_)
            self._sub_data = {key: np.empty((len(self._coherenceAxis), len(self._datasets[key]), len(self._wCenters)),
                                            dtype=self._datasets[key].data.dtype)
                              for key in self._datasets}

        # save sub step
        self._recorded_coherences[coherence_step] = True
        for key in data:
            self._sub_data[key][coherence_step] = data[key]

        # note that this will consider the sub step finished as soon as all sub steps are filled disregarding which
        # datasets are filled (or even if any dataset is filled completely)
        # so please fill all datasets at once for every substep
        # this is not trying to be foolproof
        if self._recorded_coherences.all():
            self.setDataStep(self._recorded_population, self._sub_data)
            self._recorded_coherences = None
            self._recorded_population = None
            self._sub_data = None
        pass
