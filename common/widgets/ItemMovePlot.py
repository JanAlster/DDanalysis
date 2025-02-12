
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

# -*- coding: utf-8 -*-

from PyQt5 import QtCore, QtWidgets

from pqr import pqr2 as pqr
from pqr.pqr2 import FigureWidget

#TODO: this is transition from QCP to pqr, probably not the best way to handle things
class ModPlot(pqr.Plot):
    itemMoved=QtCore.pyqtSignal(pqr.OriginDecorator, QtCore.QPointF)
    def addDecorator(self, decorator):
        decorator.originChanged.connect(lambda x, y: self.itemMoved.emit(decorator, QtCore.QPointF(x, y)))
        return super().addDecorator(decorator)

class ItemMovePlot(FigureWidget):
    itemMoved=QtCore.pyqtSignal(pqr.OriginDecorator, QtCore.QPointF)

    def __init__(self, *args):
        super().__init__(addFirstPlot=False)
        # ~ plot = self.addPlot(0, 0, plotClass=ModPlot)
        #problem is that QCP version output delta movement and this gives absolute position
        #also units mightbe different
        
    def addPlot(self, *args, plotClass=ModPlot, **kw):
        plot = super().addPlot(*args, plotClass=ModPlot, **kw)
        plot.itemMoved.connect(self.itemMoved)
        return plot
    
"""
odl QCP version

#from .SavingPlot import SavingPlot, QCP
from common.widgets.PluginPlot import SavingPlot, QCP

class ItemMovePlot(SavingPlot):
    itemMoved=QtCore.pyqtSignal(QCP.QCPAbstractItem, QtCore.QPointF)
    #cursorMoveIndex=QtCore.pyqtSignal(int, int)

    def __init__(self, *args):
        super(ItemMovePlot, self).__init__(*args)
        
        self.mouseMoveTimer=QtCore.QElapsedTimer()
        self.mouseMoveTimer.start()#Todo?: only start this in mouse press event? and stop in mouse release event? but possibly it is not needed
        self._movingItem=None
        self._lastMousePos=None
        pass
    
    def mousePressEvent(self, mpe):
        #print "ItemMovePlot.mousePressEvent"
        if mpe.button()==QtCore.Qt.LeftButton:
            item=self.itemAt(mpe.pos())
            if item and item.selectable():
                print("ItemMovePlot moving item:",item)
                item.setSelected(True)
                self._movingItem=item
                self._lastMousePos=mpe.pos()
                mpe.accept()
                self.replot()
                return
            pass
        return super(ItemMovePlot, self).mousePressEvent(mpe)

    
    def mouseMoveEvent(self, mme):
        if self.mouseMoveTimer.elapsed()<20: #otherwise this has the potential to suffocate the GUI
            return QCP.QCustomPlot.mouseMoveEvent(self, mme)
            
        if self._movingItem!=None:
            #move the item by delta mouse
            mousePos=mme.pos()
            
            #do not move if pos is outside of plot
            if self.axisRect().rect().contains(mousePos, True):
                delta=mousePos-self._lastMousePos
            
                for it in self._movingItem.positions():
                    t=it.type()
                    it.setType(it.ptAbsolute)
                    it.setPixelPoint(it.pixelPoint()+delta)
                    it.setType(t)

                self._lastMousePos=mousePos
                self.itemMoved.emit(self._movingItem, delta)
                pass
                        
            #accept and report
            mme.accept()
            self.replot()
            self.mouseMoveTimer.restart()
            return
        
        return super(ItemMovePlot, self).mouseMoveEvent(mme)

        
    def mouseReleaseEvent(self, mre):
        if mre.button()==QtCore.Qt.LeftButton:
            if self._movingItem!=None:
                self._movingItem.setSelected(False)
                self._movingItem=None
                self._lastMousePos=None
                mre.accept()
                self.replot()
                return
            pass
        return super(ItemMovePlot, self).mouseReleaseEvent(mre)

    pass
    
"""    
    
