
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

#Done?: work in progress - it is causing Illegal instruction sometimes :( - disconnect in removeWidget possibly helped

from PyQt5 import QtWidgets, QtCore

class GrowingList(QtWidgets.QFrame):
    """
    Groups widgetClass widgets and keep one empty as the last one. 
    Widgets have to emit SomeValue signal when it changes from empty and NoValue when it changes to empty (it will be removed).
    """
    widgetAdded=QtCore.pyqtSignal(QtWidgets.QWidget)
    #widgetRemoved=QtCore.pyqtSignal(QtWidgets.QWidget)
    valueChanged=QtCore.pyqtSignal()
    
    def __init__(self, widgetClass, getter=None, setter=None, *args):
        super(GrowingList, self).__init__(*args)
        self._widgetClass=widgetClass
        #test getter/setter (todo? however shouldn't we test it on instance? they could be created dynamically)
        #getter/setter should be name of widgetClass method, or widgetClass.method
        if hasattr(getter, '__call__'):
            if hasattr(widgetClass, getter.__name__):
                self._getter=getter
            else:
                self._getter=None
        else: 
            try:  
                self._getter=getattr(widgetClass, getter)
            except:
                self._getter=None
            
        if hasattr(setter, '__call__'):
            if hasattr(widgetClass, setter.__name__):
                self._setter=setter
            else:
                self._setter=None
        else: 
            try:  
                self._setter=getattr(widgetClass, setter)
            except:
                self._setter=None
            
        layout=QtWidgets.QVBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSizeConstraint(layout.SetFixedSize)
        self.setLayout(layout)
        self._signalsToResend=[]
        self.newWidget()
        pass
    
    #def registerSignalToResend(self, signal):
        #self.signalsToResend.append(signal)
        ##retroactive
        #l=self.layout()
        #for i in range(l.count()):
            #wid=l.itemAt(i).widget()
            #self.connect(wid, signal, self, signal)
            #pass
        #pass
    
    def count(self):
        return self.layout().count()-1 #one is not filled (hopefully this update will not break anything)
    
    def newWidget(self):
        """append a new Limit widget to the end"""
        #print "GrowingList.newWidget",
        print("GrowingList.newWidget")
        layout=self.layout()
        wid=self._widgetClass()
        #print wid
        
        wid.SomeValue.connect(self.newWidget)
        wid.SomeValue.connect(self.valueChanged)
        #for it in self.signalsToResend:
            #self.connect(wid, it, self, it)
            ##import sys
            ##self.connect(wid, it, lambda : sys.__stdout__.write(str(wid)+' signal '+str(it)+'\n'))
            #pass
        
        if layout.count()>=1:
            #print "disconnect",layout.itemAt(layout.count()-1).widget()
            prevwid=layout.itemAt(layout.count()-1).widget()
            prevwid.SomeValue.disconnect(self.newWidget) #only last widget will cause addition on a new one after it changes from NoLimit
            prevwid.NoValue.connect(self.removeWidgetCall) #only the last one should not be removed, but that will never emit NoLimit (since it already is NoLimit)
        layout.addWidget(wid)
        wid.show() #otherwise wid position will not be update here and now (ad we want that)
        self.widgetAdded.emit(wid)
        pass
    
    def removeWidgetCall(self):
        wid=self.sender()
        self.removeWidget(wid)
        pass
    
    def removeWidget(self, wid):
        #print "remove widget"
        #remove wid from layout
        wid.NoValue.disconnect()#seems like it is needed to avoid segfaults
        wid.SomeValue.disconnect()
        self.layout().removeWidget(wid)
        wid.setParent(None) #this should also delete it 
        del wid 
        #self.widgetRemoved.emit()
        #print "remove widget end"
        self.valueChanged.emit()
        pass
    
    def removeAllWidgets(self):
        for i in range(self.count()-1):
            self.removeWidget(self.widget(0))
            pass
        pass
        
    def widget(self, index):
        if index>=self.layout().count():
            return None
        return self.layout().itemAt(index).widget()
    
    def indexOf(self, widget):
        return self.layout().indexOf(widget)
    
    def value(self, getter):
        if self._getter==None: raise UserError
        data=[]
        l=self.layout()
        for i in range(l.count()-1):#last one is empty
            wid=l.itemAt(i).widget()
            data.append(self._getter(wid))
            pass
        return data

    def setValue(self, data, setter):
        if self._setter==None: raise UserError
        #match the list length to number of data
        l=self.layout()
        for i in range(l.count()-1-len(data)):
            self.removeWidget(self.widget(0))#let the last one intact
            pass
        
        for i in range(len(data)-l.count()+1):
            self.newWidget()
            pass
        
        assert(len(data)==l.count()-1)
        for i in range(len(data)):
            wid=l.itemAt(i).widget()
            self._setter(wid, data[i])
            pass
        pass
    pass

class SortableGrowingList(GrowingList):
    """Growing list with ability to sort its widgets upon request. You must provide method sortValue(widget) providing some sortable (sorted list) value of widgets."""
    def __init__(self, widgetClass, sortValue, *args):
        super(SortableGrowingList, self).__init__(widgetClass, *args)
        self._sortValue=sortValue
        pass
    
    def sort(self, ascending=True):
        """sort widgets in list according to sortValue (the last empty one would remain last)"""
        l=self.layout()
        m={}
        for i in range(l.count()-1):
            w=self.widget(i)
            m[self._sortValue(w)]=w
            pass
        a=list(m.keys())
        b=sorted(a)
        if not ascending:
            b=b[::-1]
        if a!=b:
            last=self.widget(l.count()-1)
            for i in range(l.count()):
                l.removeItem(l.itemAt(0))
                pass
            for it in b:
                l.addWidget(m[it])
                pass
            l.addWidget(last)
            pass
        pass
