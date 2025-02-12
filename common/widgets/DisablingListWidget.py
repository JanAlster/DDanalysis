
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

from PyQt5 import QtWidgets, QtGui, QtCore

class DisablingListWidget(QtWidgets.QFrame):
    currentChanged=QtCore.pyqtSignal(int)
    disabledChanged=QtCore.pyqtSignal(int)  #index of changed item or -1 if more are changed
    def __init__(self, *args):
        super().__init__(*args)
        
        #we need button and listwidget
        layout = QtWidgets.QVBoxLayout()
        self._pushToggle = QtWidgets.QPushButton("Toggle")
        self._list = QtWidgets.QListWidget(self)
        self._list.setSelectionMode(self._list.ExtendedSelection)
        self._list.itemChanged.connect(self._itemChanged)
        self._list.currentRowChanged.connect(self.currentChanged)
        self._list.itemDoubleClicked.connect(self._toggleSelected)
        
        self._pushToggle.clicked.connect(self._toggleSelected)
        
        layout.addWidget(self._pushToggle)
        layout.addWidget(self._list, 1)
        
        self.setLayout(layout)
        pass
    
    #TODO: ensure self._list has the size of model, do not show scrollbars, remove stretch from tabData.ui
    #~ QAbstractScrollArea::setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff )
    #~ QAbstractScrollArea::setVerticalScrollBarPolicy( Qt::ScrollBarAlwaysOff ) 
    
    def _slot(self, *args):
        print("slot", args)
        pass
    
    def _itemChanged(self, item=None):
        self.disabledChanged.emit(self._list.row(item) if item is not None else -1) #assuming the only change we can do is change checkstate
        
    def count(self):
        return self._list.count()
        
    def clear(self, emit=True):
        self._list.itemChanged.disconnect(self._itemChanged)
        self._list.clear()
        self._list.itemChanged.connect(self._itemChanged)
        if emit: self._itemChanged()
        pass
        
    def setItems(self, items, force=True):
        #if you change this consider changing tabLA.LAInputList too
        if not force:
            #only do this if items are different from current items
            if len(items)!=self._list.count():
                force = True
            else:
                for i in range(len(items)):
                    if str(items[i]) != self._list.item(i).text():
                        force = True
                        break
                    pass
        
        if force:
            self._list.itemChanged.disconnect(self._itemChanged)
            self._list.clear()
            self._list.addItems([str(it) for it in items])
            self._list.itemChanged.connect(self._itemChanged)
            self.setDisabled([]) #here it can emit disabledChanged (once)
            
            c = self._list.currentRow()
            if c>=0: self.currentChanged.emit(c)
        pass
        
    def current(self): 
        return self._list.currentRow()

    def enabled(self):
        return [i for i in range(self._list.count()) if self._list.item(i).checkState()==QtCore.Qt.Checked]

    def disabled(self):
        return [i for i in range(self._list.count()) if self._list.item(i).checkState()!=QtCore.Qt.Checked]

    def selectedItems(self):
        sel =  self._list.selectedItems()
        return [i for i in range(self._list.count()) if self._list.item(i) in sel]

    def _toggleSelected(self, *args):
        l = self._list.selectedItems()
        if len(l)>0:
            self._list.itemChanged.disconnect(self._itemChanged)
            for it in l:
                it.setCheckState(QtCore.Qt.Checked if it.checkState()==QtCore.Qt.Unchecked else QtCore.Qt.Unchecked)
            self._list.itemChanged.connect(self._itemChanged)
            self._itemChanged(l[0])
        pass

    def setDisabled(self, indices):
        #this will enable all others
        self._list.itemChanged.disconnect(self._itemChanged)
        for i in range(self._list.count()):
            self._list.item(i).setCheckState(QtCore.Qt.Checked if i not in indices else QtCore.Qt.Unchecked)
        self._list.itemChanged.connect(self._itemChanged)
        self._itemChanged() #to keep things simple we will not test if there was an actual change or if it was one or more elements, just signal disabledChanged once
        pass

    def setCurrent(self, index):
        if index==self.current(): return 
        self._list.setCurrentRow(index)
        pass

    #~ def setCurrent(self, index):
        #~ #print( "DisablingListWidget.setCurrent", index)
        #~ if index==self._current:
            #~ return
            #~ 
        #~ if self._current>=0:
            #~ #print("\t", self._current, self._layout.itemAt(self._current))
            #~ w=self._layout.itemAt(self._current).widget()
            #~ w.setStyleSheet("")
            #~ #clear old current
            #~ pass
        #~ self._current=index
        #~ #mark new current
        #~ if self._current>=0:
            #~ w=self._layout.itemAt(self._current).widget()
            #~ w.setStyleSheet("background-color: #97ffff;")
        #~ 
        #~ if self._current>=0: self.currentChanged.emit(self._current)
        #~ pass
    pass
    
    
class DisablingListWidgetOld(QtWidgets.QScrollArea):
    """ListWidget that can disable items by double click
    and return list of indices of enabled items

    it does not change current when double clicking
    current is highlighted
    """
    currentChanged=QtCore.pyqtSignal(int)
    disabledChanged=QtCore.pyqtSignal(int)
    
    def __init__(self, *args):
        super(DisablingListWidget, self).__init__(*args)
        self.timer=QtCore.QTimer()
        self.timer.setSingleShot(True)
        self.timer.setInterval(QtWidgets.QApplication.doubleClickInterval()+50)
        self.timer.timeout.connect(self.timeout)

        widget=QtWidgets.QWidget(self)
        self.setBackgroundRole(QtGui.QPalette.Base)
        layout=QtWidgets.QVBoxLayout()
        layout.setContentsMargins(3, 3, 3, 3)
        layout.setSizeConstraint(layout.SetFixedSize)
        widget.setLayout(layout)

        self.setWidget(widget)
        self._layout=layout
        self._current=-1
        self._future=-1
        pass

    def count(self):
        return self._layout.count()

    def clear(self):
        cl=self._layout
        #clear if not empty
        while cl.count()>0:
            it=cl.takeAt(0)
            it.widget().setParent(None)
            pass
        self._current=-1
        pass
    
    def setItems(self, items):
        self.clear()
        cl=self._layout
        for it in items:
            w=QtWidgets.QLabel(str(it))
            cl.addWidget(w)
            pass
        self.setCurrent(0 if len(items)>0 else None)
        pass

    def timeout(self):
        if self._future>=0: self.setCurrent(self._future)
        pass
        

    def mousePressEvent(self, mpe):
        #print "DisablingListWidget.mousePressEvent", mpe.pos(), self.verticalScrollBar().value()
        w=self.widget().childAt(mpe.x(), mpe.y()+self.verticalScrollBar().value())
        if w is not None:
            self._future=self._layout.indexOf(w)
            self.timer.start()
        return super(DisablingListWidget, self).mousePressEvent(mpe)
        
    def mouseDoubleClickEvent(self, dce):
        w=self.widget().childAt(dce.x(), dce.y()+self.verticalScrollBar().value())
        if w is not None:
            w.setEnabled(not w.isEnabled())
            self.disabledChanged.emit(self._layout.indexOf(w))
            pass
        self.timer.stop()
        dce.accept()
        pass

    def keyPressEvent(self, kpe):
        print("DisablingListWidget.keyPressEvent", kpe.key())
        key = kpe.key()
        #moving in list so smaller T are more up
        if key == QtCore.Qt.Key_Up and self._current > 0:
            self.setCurrent(self._current-1)
        elif key == QtCore.Qt.Key_Down and self._current < self.count()-1:
            self.setCurrent(self._current+1)
        else:
            return super().keyPressEvent(kpe)
        pass

    def current(self): 
        return self._current

    def enabled(self):
        return [i for i in range(self._layout.count()) if self._layout.itemAt(i).widget().isEnabled()]

    def disabled(self):
        return [i for i in range(self._layout.count()) if not self._layout.itemAt(i).widget().isEnabled()]

    def setDisabled(self, indices):
        #done: shouldn't this enable all others?
        #for i in indices:
        #    assert i>=0 and i<self._layout.count()
        #    self._layout.itemAt(i).widget().setEnabled(False)
        for i in range(self._layout.count()):
            self._layout.itemAt(i).widget().setEnabled(i not in indices)
        pass

    def setCurrent(self, index):
        #print( "DisablingListWidget.setCurrent", index)
        if index==self._current:
            return
            
        if self._current>=0:
            #print("\t", self._current, self._layout.itemAt(self._current))
            w=self._layout.itemAt(self._current).widget()
            w.setStyleSheet("")
            #clear old current
            pass
        self._current=index
        #mark new current
        if self._current>=0:
            w=self._layout.itemAt(self._current).widget()
            w.setStyleSheet("background-color: #97ffff;")
        
        if self._current>=0: self.currentChanged.emit(self._current)
        pass
    pass
    
if __name__=="__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    app = QApplication(sys.argv)
    class TempMainWindow(QtWidgets.QMainWindow):
        def __init__(self, parent=None):
            super().__init__(parent)
            self._central=DisablingListWidget2()
            self.setCentralWidget(self._central)
            self._central.setItems(["ahoj", "nazdar", "1", "2", "3", "4", "5"])
            pass
        pass
        
    window=TempMainWindow()
    window.show()
    res=app.exec_()
    sys.exit(res)    
