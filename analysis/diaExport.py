
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

from PyQt5.QtWidgets import QApplication
from PyQt5.uic import loadUiType

import numpy as np
import os.path as op

from analysis import data_pipe
from common.interpolator import interpolate

from PyQt5 import QtWidgets

from pqr import pqr

#~ exportDialogUI, exportDialogClass = loadUiType('analysis/exportDialog.ui')
exportDialogUI, exportDialogClass = loadUiType(__file__.replace("diaExport.py", 'exportDialog.ui'))
class ExportDialog(exportDialogClass, exportDialogUI):
    def __init__(self, parent=None):
        exportDialogClass.__init__(self, parent)
        self.setupUi(self)
        self.bScopeCurrent.toggled.connect(self._updateScope)
        self.bScopeSelected.toggled.connect(self._updateScope)
        self.bScopeAll.toggled.connect(self._updateScope)
        
        plot = self.plot.plot()
        xaxis = plot.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.B, label='ω₁ / cm⁻¹', tickFormat="{:.0f}", tickNum=4)
        yaxis = plot.addAxis(pqr.NumEnumTransform(), pqr.AxisPosition.L, label='ω₃ / cm⁻¹', tickFormat="{:.0f}", tickNum=4)
        xaxis2 = plot.addAxis(pqr.AOverTransform(1e7, xaxis._transform), pqr.AxisPosition.top, label="λ₁ / nm", tickFormat="{:.0f}")
        yaxis2 = plot.addAxis(pqr.AOverTransform(1e7, yaxis._transform), pqr.AxisPosition.right, label="λ₃ / nm", tickFormat="{:.0f}")
        
        xaxis.zoomRangeLChanged.connect(xaxis2._setZoomRange)
        xaxis2.zoomRangeLChanged.connect(xaxis._setZoomRange)
        yaxis.zoomRangeLChanged.connect(yaxis2._setZoomRange)
        yaxis2.zoomRangeLChanged.connect(yaxis._setZoomRange)

        plotter = plot.addDDPlotter(xaxis, yaxis)
        
        #self._plot = QCP.QCustomPlot(None) #or PluginPlot DD+Base
        #self._cm = QCP.QCPColorMap(self._plot.xAxis, self._plot.yAxis)
        #self._plot.addPlottable(self._cm)
        #self._cmdata = self._cm.data()
        self._shownDatasets = None
        self.iShow.valueChanged.connect(self._updateShownFrame)
        self.cbShow.currentTextChanged.connect(self._updateShownDataset)
        self.pushShow.clicked.connect(self._updateShown)
        self.pushAdjustToFrame.clicked.connect(self._updateAspectFrame)
        self.pushAdjustToData.clicked.connect(self._updateAspectData)
        
        self.iWidth.valueChanged.connect(self._updatePlotWidth)
        self.iHeight.valueChanged.connect(self._updatePlotHeight)
        self._updatePlotSize()
        
        self.siAxis1.valueChanged.connect(self._updateInterpolationWidth)
        self.siAxis3.valueChanged.connect(self._updateInterpolationHeight)
        
        self.bExportASCII.toggled.connect(self._updateExport)
        self.bExportImage.toggled.connect(self._updateExport)
        self.gbImageSettings.setVisible(self.bExportImage.isChecked())
        self.plot.setVisible(self.bExportImage.isChecked())
        self.adjustSize()
        #self._plot.setVisible(self.bExportImage.isChecked())
        
        self.pushExport.clicked.connect(self._export)
        
        self._bDatasets = [self.bData, self.bDAS, self.bRes, self.bFreq]
        main = self.parent()
        self._tabs = [main.tabData, main.tabGA, main.tabResidues, main.tabOM, None]
        self._names = ["DD", "DAS", "R", "OM"]
        pass

    def folder(self):
        fn=QtWidgets.QFileDialog.getExistingDirectory(self, directory=self.leFolder.text())
        if fn!="":
            self.leFolder.setText(fn)
        pass

    def scope(self):
        if self.bScopeCurrent.isChecked(): return "current"
        if self.bScopeSelected.isChecked(): return "selected"
        if self.bScopeAll.isChecked(): return "all"
        pass
    
    def datasets(self):
        return [(self._names[i],self._tabs[i]) for i in range(len(self._bDatasets)) if self._bDatasets[i].isChecked() and self._bDatasets[i].isEnabled()]

    def _updateInterpolationWidth(self, width):
        if self.bKeepAspectRatioInterpolation.isChecked():
            #we need data aspect ratio 
            a = self._dataAspectRatio[0]
        
            #adjust height too
            self.siAxis3.valueChanged.disconnect(self._updateInterpolationHeight)
            self.siAxis3.setValue(int(round(a*width, 0)))
            self.siAxis3.valueChanged.connect(self._updateInterpolationHeight)
        pass

    def _updateInterpolationHeight(self, height):
        if self.bKeepAspectRatioInterpolation.isChecked():
            #we need data aspect ratio 
            a = self._dataAspectRatio[0]

            #adjust height too
            self.siAxis1.valueChanged.disconnect(self._updateInterpolationWidth)
            self.siAxis1.setValue(int(round(height/a, 0)))
            self.siAxis1.valueChanged.connect(self._updateInterpolationWidth)
        pass


    def _updatePlotWidth(self, width):
        if self.bKeepAspectRatio.isChecked():
            #adjust height too
            self.iHeight.valueChanged.disconnect(self._updatePlotHeight)
            #we have to get the old size from plot
            # ~ v = self.plot.viewport()
            v = self.plot.size() #TODO: this is just a quess, there is not proper transition to pqr yet
            self.iHeight.setValue(int(round(float(v.height())/v.width()*width, 0)))
            self.iHeight.valueChanged.connect(self._updatePlotHeight)
        self._updatePlotSize()
        pass

    def _updatePlotHeight(self, height):
        if self.bKeepAspectRatio.isChecked():
            #adjust height too
            self.iWidth.valueChanged.disconnect(self._updatePlotWidth)
            #we have to get the old size from plot
            # ~ v = self.plot.viewport()
            v = self.plot.size() #TODO: this is just a quess, there is not proper transition to pqr yet
            self.iWidth.setValue(int(round(float(v.width())/v.height()*height, 0)))
            self.iWidth.valueChanged.connect(self._updatePlotWidth)
        self._updatePlotSize()
        pass

    def _updatePlotSize(self, *args):
        #print("ExportDialog._updatePlotSize", args, self.plot.size(), self.plot.minimumSize(), self.plot.maximumSize())
        self.plot.setMinimumSize(self.iWidth.value(), self.iHeight.value())
        self.plot.setMaximumSize(self.iWidth.value(), self.iHeight.value())
        print("ExportDialog._updatePlotSize", self.pos(), self.size())
        self.adjustSize()
        print("ExportDialog._updatePlotSize", self.pos(), self.size())
        #TODO: after this the dialog moves a bit, not sure where in code it occurs?
        pass

    def _updateExport(self, checked):
        if checked:
            self.gbImageSettings.setVisible(self.bExportImage.isChecked())
            self.plot.setVisible(self.bExportImage.isChecked())
            print("ExportDialog._updateExport", self.pos(), self.size())
            self.adjustSize()
            print("ExportDialog._updateExport", self.pos(), self.size())
            #self._plot.setVisible(self.bExportImage.isChecked())
            #print("ExportDialog._updateExport", self._plot.size())
        pass

    def _updateScope(self, checked):
        #setup check for selections
        #scope
        # - current forces dataset (and trn/pair?)
        # - selection forces dataset
        # - all allows to choose dataset
        #if checked: #we will also get the unchecked signals, which we can ignore
        #    if self.bScopeCurrent.isChecked() or self.bScopeSelected.isChecked():
        #       self.gbDataset.setEnabled(False) #Todo: possibly choose current dataset from UI (and remember what was selected here before?)
        #    else:
        #        self.gbDataset.setEnabled(True)
        #    pass
        pass

    def _getDatasets(self):
        print("_getDatasets")
        scope = self.scope()
        #TODO: limit to ZoomedIn region?
        # - not sure if interpolate only the zoomed in region or the whole frame (probably the former)

        #we need data, a2, a1, a3, ai1, ai3
        res = {}
        for name, tab in self.datasets():
            dataset = tab.dataExport()
            if scope == "current":
                indices = [tab.listFrames.current()]
            elif scope == "selected":
                indices = tab.listFrames.enabled()
            else:
                indices = slice(None)
            
            a1 = dataset["a1"]
            a3 = dataset["a3"]
            try:
                a2 = dataset["a2"][indices]
            except:
                a2 = [dataset["a2"][it] for it in indices]

            data = dataset["total"][indices]
            #Todo? use PAIR from DDdata.py
            if self.bTR.isChecked(): res[name+"TR"] = data.real, a1,a2, a3
            if self.bTI.isChecked(): res[name+"TI"] = data.imag, a1, a2, a3
            if self.bTA.isChecked(): res[name+"TA"] = np.abs(data), a1, a2, a3
            if self.bTP.isChecked(): res[name+"TP"] = np.angle(data), a1, a2, a3
            if self.bTN.isChecked(): res[name+"TN"] = data_pipe.amplitude_to_noise(data), a1, a2, a3

            data = dataset["rephasing"][indices]
            if self.bRR.isChecked(): res[name+"RR"] = data.real, a1, a2, a3
            if self.bRI.isChecked(): res[name+"RI"] = data.imag, a1, a2, a3
            if self.bRA.isChecked(): res[name+"RA"] = np.abs(data), a1, a2, a3
            if self.bRP.isChecked(): res[name+"RP"] = np.angle(data), a1, a2, a3
            if self.bRN.isChecked(): res[name+"RN"] = data_pipe.amplitude_to_noise(data), a1, a2, a3

            data = dataset["nonrephasing"][indices]
            if self.bNR.isChecked(): res[name+"NR"] = data.real, a1, a2, a3
            if self.bNI.isChecked(): res[name+"NI"] = data.imag, a1, a2, a3
            if self.bNA.isChecked(): res[name+"NA"] = np.abs(data), a1, a2, a3
            if self.bNP.isChecked(): res[name+"NP"] = np.angle(data), a1, a2, a3
            if self.bNN.isChecked(): res[name + "NN"] = data_pipe.amplitude_to_noise(data), a1, a2, a3
            pass

        #interpolation
        ai1 = self.siAxis1.value() #interpolation
        ai3 = self.siAxis3.value() #interpolation

        if ai1>0 or ai3>0:
            for name in res:
                data, a1, a2, a3 = res[name]
                if ai1>0:
                    data, a1 = interpolate(data, ai1, a1, axis=1)
                if ai3>0:
                    data, a3 = interpolate(data, ai3, a3, axis=2)
                
                res[name] = data, a1, a2, a3
                pass
            pass
              
        #TODO: normalization and global colors does not work really global - it only works on selected subset (indices here)
        #normalization
        for name in res:
            data, a1, a2, a3 = res[name]
            
            if name[-1] == "P": #different behaviour for phase data  (no normalization, fixed data range)
                dataRange = (-np.pi, np.pi)
            else:
                if self.bNormalizeZScale.isChecked():
                    data = data/np.abs(data).max()
                    
                dataRange = (data.min(), data.max()) if self.bGlobalZScale.isChecked() else None
            res[name] = data, a1, a2, a3, dataRange
            pass

        print("\t", list(res.keys()))
        return res

    def _showImage(self, frame, a1, a3, **kw):
        plot = self.plot.plot()
        print("_showImage", a1, a3, frame, kw)
        plotter = plot.plotter(0)
        plotter._xaxis._transform.stops = a1
        plotter._yaxis._transform.stops = a3
        plotter.setData(frame.T, a1, a3, **kw)
        plot.zoomToData()

        # ~ plot = self.plot
        # ~ print(plot.viewport(), plot.axisRect().rect())
        # ~ w=float(a1.size)/plot.axisRect().rect().width()
        # ~ h=float(a3.size)/plot.axisRect().rect().height()
        # ~ #width and height params times scale param give you new viewport (size of image)
        # ~ #plot.saveJpg(basename % (t2[0],), 320,240, 2*max(1., min(w,h)))
        # ~ print("_showImage", w, h, 2*max(1., min(w,h)))
        pass


    def _updateShownFrame(self, i):
        print("_updateShownFrame", i)
        dataset = self.cbShow.currentText()
        if dataset not in self._shownDatasets: return
        
        data, a1, a2, a3, dataRange = self._shownDatasets[dataset]
        frame = data[i]
        a2 = a2[i]
        
        #plot
        self._showImage(frame, a1, a3, dataRange=dataRange)
        
        pass
        
    def _updateShownDataset(self, dataset):
        print("_updateShownDataset")
        if dataset not in self._shownDatasets: return
        self.iShow.setMaximum(len(self._shownDatasets[dataset][0])-1)
        #trigger frame update Todo: if range reset did not
        self._updateShownFrame(self.iShow.value())
        pass
    
    def _updateShown(self):
        print("_updateShown")
        self._shownDatasets = self._getDatasets()
        
        prev = self.cbShow.currentText(), self.iShow.value()
        
        self.cbShow.clear()
        self.cbShow.addItems(list(self._shownDatasets.keys()))
        if prev[0] in self._shownDatasets:
            self.cbShow.setCurrentText(prev[0])
            self.iShow.setValue(prev[1])
            pass
        pass
        
    def _updateAspectFrame(self):
        #we assume that all frames are the same size
        #we adjust plot size to reflect well frame size 
        # - either viewport of plot has the pixel dimension of frame (or of interpolation) - results in square pixels
        # - or viewport has aspect ratio of axes dimensions - results in same size x and y step
        #lets keep this simple
        # ~ ar = self.plot.axisRect()
        ar = self.plot.plot().size() #TODO: this is just a quess, there is not proper transition to pqr yet
        
        #self._updateShown()
        #data, a1, a2, a3, dataRange = self._shownDatasets[list(self._shownDatasets.keys())[0]]
        #a = float(len(a3)) / len(a1) #this is the aspect ratio we want the axisRect to have
        a = self._dataAspectRatio[0]
        #but we are setting viewport (whole plot size)
        #so the question to what set the viewport (plot) size to have axisRect size with given aspect ratio
        
        arh = ar.height()
        narh = int(round(ar.width() * a,0))
        h = self.iHeight.value()
        nh = h + narh-arh
        self.iHeight.valueChanged.disconnect(self._updatePlotHeight)
        self.iHeight.setValue(nh)
        self.iHeight.valueChanged.connect(self._updatePlotHeight)
        self._updatePlotSize()
        pass
        
    def _updateAspectData(self):
        #we assume that all frames are the same size
        #we adjust plot size to reflect well frame size 
        # - either viewport of plot has the pixel dimension of frame (or of interpolation) - results in square pixels
        # - or viewport has aspect ratio of axes dimensions - results in same size x and y step
        #lets keep this simple
        # ~ ar = self.plot.axisRect()
        ar = self.plot.plot().size() #TODO: this is just a quess, there is not proper transition to pqr yet

        #self._updateShown()
        #data, a1, a2, a3, dataRange = self._shownDatasets[list(self._shownDatasets.keys())[0]]
        #a = abs(a3[0]-a3[-1]) / abs(a1[0]-a1[-1])
        a = self._dataAspectRatio[1]

        arh = ar.height()
        narh = int(round(ar.width() * a,0))
        h = self.iHeight.value()
        nh = h + narh-arh
        self.iHeight.valueChanged.disconnect(self._updatePlotHeight)
        self.iHeight.setValue(nh)
        self.iHeight.valueChanged.connect(self._updatePlotHeight)
        self._updatePlotSize()
        pass

    def exec_(self):
        print("ExportDialog.exec_")
        #TODO: set plot size according to data aspect ratio (at least the first time around)
        # - try to keep tho original size
        
        
        #dis/enable UI to reflect which datasets are available
        for i in range(len(self._bDatasets)):
            self._bDatasets[i].setEnabled(self._tabs[i].dataExport() is not None)
        
        #by default select presently selected tab
        for it in self._bDatasets:
            it.setChecked(False)
        tab = self.parent().tabWidget.currentIndex()
        # FIXME: tab listed in self._bDatasets do NOT correspond to tabs in main window anymore, so either add all new tabs to self._bDatasets (name and tab links in __init__) or remove the dependency on currentIndex()
        if tab>=0: self._bDatasets[tab].setChecked(True)
        
        self._currentInUI = self._tabs[tab]
        #TODO: select also current TRN, PAIR

        #pre calculate data aspect ratio
        #any date really, all should have the same dimensions
        data = self._tabs[0].dataExport()
        if data is not None:
            a1 = data["a1"]
            a3 = data["a3"]
            self._dataAspectRatio = float(len(a3)) / len(a1), float(abs(a3[0]-a3[-1])) / abs(a1[0]-a1[-1])
        else:
            self._dataAspectRation = None, None

        res = super().exec_()
        #self._plot.setVisible(False)
        return res
        
    def _export(self, *args):
        #export as images or ascii
        
        print("Exporting")
        try:
            basedir = self.leFolder.text()
            
            def normname(base, T, i, ext):
                if T is None:
                    return base+"."+ext
                if T.imag != 0:
                    return base+f"_{i}_{int(T.real):d}fs_{int(T.imag):d}icm."+ext
                return base+f"_{i}_{int(T.real):d}."+ext
            
            datasets = self._getDatasets()
            if self.bExportASCII.isChecked():
                for name in datasets:
                    data, a1, a2, a3, dataRange = datasets[name]
                    for i in range(len(data)):
                        frame = data[i]
                        T = a2[i]
                        #filename = op.join(basedir, name+"_%d.txt" % (T,))
                        filename = op.join(basedir, normname(name, T, i, "txt"))
                        e = np.empty((frame.shape[0]+1,frame.shape[1]+1), dtype=frame.dtype)
                        e[1:,1:] = frame
                        e[0, 0] = T.real
                        e[0, 1:] = a3
                        e[1:,0] = a1
                        np.savetxt(filename, e, '%.10e')
                        pass
                    pass
                pass
            elif self.bExportNpyFrames.isChecked():
                for name in datasets:
                    data, a1, a2, a3, dataRange = datasets[name]
                    for i in range(len(data)):
                        frame = data[i]
                        T = a2[i]
                        #filename = op.join(basedir, name+"_%d.npy" % (T,))
                        filename = op.join(basedir, normname(name, T, i, "npy"))
                        e = np.empty((frame.shape[0] + 1, frame.shape[1] + 1), dtype=frame.dtype)
                        e[1:, 1:] = frame
                        e[0, 0] = T.real
                        e[0, 1:] = a3
                        e[1:, 0] = a1
                        np.save(filename, e)
                        pass
                    #whole matrix
                    #we need data to be 3D array first in case it is not (one frame export)
                    if len(data.shape)==2:
                        data = data.reshape(1, data.shape[0], data.shape[1])
    #TODO: tady

                    pass
                pass
            elif self.bExportNpySingle.isChecked():
                for name in datasets:
                    data, a1, a2, a3, dataRange = datasets[name]
                    filename = op.join(basedir, normname(name, None, None, "npz"))
                    np.savez(filename, data=data, a1=a1, a2=a2, a3=a3)
                pass
            elif self.bExportImage.isChecked():
                format = self.cbType.currentText()
                #w = self.iWidth.value()
                #h = self.iHeight.value()
                s = self.iScale.value()
                #dpi = self.iResolution.value() #this is not supported by my version of QCP (I hate it when they do new version and I do not have the old API documentation)
                
                # ~ antialiased = self.bAntialiased.isChecked()
                # ~ if antialiased:
                    # ~ old = self.plot.antialiasedElements() 
                    # ~ self.plot.setAntialiasedElements(QCP.QCP.aeAll)
                
                
                for name in datasets:
                    data, a1, a2, a3, dataRange = datasets[name]
                    for i in range(len(data)):
                        frame = data[i]
                        T = a2[i]
                        #filename = op.join(basedir, name+"_%d." % (T,))+format
                        filename = op.join(basedir, normname(name, T, i, format))
            
                        #plot
                        self._showImage(frame, a1, a3, dataRange=dataRange)
                        #TODO: we might do self plot as SavingDDPlot and use saving to save the image (but that would have to be reworked too - it asks user for filename and is just a quick draft anyway)
                        #TODO: does not have to replot probably
                        #TODO: use string .format, adjust unit to dataset
                        #TODO: image saving is not propertly trnasitioned to pqr yet, this needs testing
                        if format == "jpg":
                            self.plot.saveImage(filename, subpixel=s) #0,0 to keep widget size (what you see is what you get)
                        elif format == "png":
                            self.plot.saveImage(filename, subpixel=s) #0,0 to keep widget size (what you see is what you get)
                        elif format == "pdf":
                            self.plot.savePdf(filename) #True to avoid hairlines
                        pass
                    pass
                
                # ~ if antialiased:
                    # ~ self.plot.setAntialiasedElements(old)
                pass
        except:
            print("WARNING: Export failed")
            import traceback
            traceback.print_exc()
        pass
    pass

