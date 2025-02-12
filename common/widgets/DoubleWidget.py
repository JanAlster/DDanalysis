
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

#Todo: maybe reusable widgets
#qtpython StretchTextEdit - text edit which adjusts its height to fit given text for given (from layout) width
#GA SDevEdit - ~ spin box showing number with standard deviation

#TODO: these should not change position of cursor when SlowLineEdit changes text (not sure which one of them does this, could be as far as checking the value with hardware)

from PyQt5 import QtWidgets, QtCore, QtGui

#I want line edit that would emit change signal not all the time or only when lossing focus (pressing enter), but when it is not edited for some (setable) time
#i.e. once user stops typing

class SlowLineEdit(QtWidgets.QLineEdit):
    textEditedSlow=QtCore.pyqtSignal(str)
    def __init__(self, *args):
        super(SlowLineEdit, self).__init__(*args)
        self._editTimer=QtCore.QTimer()
        self._editTimer.setSingleShot(True)
        self._editTimer.timeout.connect(self._timeout)
        self.textEdited.connect(self._textEdited)
        self.textChanged.connect(lambda x: self.updateGeometry())
        pass
        
    def _timeout(self, *args):
        print("SlowLineEdit._timeout")
        print("\t", "old", self._oldStyleSheet, "current", self.styleSheet())
        self.setStyleSheet(self._oldStyleSheet)
        self.textEditedSlow.emit(self._text)
        pass
        
    def _textEdited(self, text):
        print("SlowLineEdit._textEdited", text)
        self._text=text
        if not self._editTimer.isActive(): self._oldStyleSheet = self.styleSheet()
        print("\t", self._oldStyleSheet)
        self.setStyleSheet("QLineEdit {background : rgba(255, 165, 0, 10%)}")
        self._editTimer.start(1000)
        pass
    pass

    def minimumSizeHint(self):
        #TODO: this is just slightly adjusted copy paste from QLineEdit (I don't really understand what is going on)
        self.ensurePolished()

        fm=self.fontMetrics()
        m = self.textMargins() #not sure if it is the same thing as QLineEditPrivate.topmargin, etc.
        #h = fm.height() + max(2*self.verticalMargin, fm.leading()) + m.top() + m.bottom() #cannot access verticalMargin
        h = fm.height() +  fm.leading() + m.top() + m.bottom()
        w = max(fm.width('000'), fm.width(self.text())) + m.left() + m.right()
        w += 2 # cursor blinking space
        w += 5 #seems needed
        #print("SlowLineEdit.minimumSizeHint", self.text(), fm.width('000'), fm.width(self.text()), m.left(), m.right())

        opt=QtWidgets.QStyleOptionFrame() #V2()
        self.initStyleOption(opt)
        hint=QtCore.QSize(w, h)

        hint=self.style().sizeFromContents(QtWidgets.QStyle.CT_LineEdit, opt, hint, self)
        return hint
    
    def sizeHint(self):
        return self.minimumSizeHint()



#I want something alike QDoubleSpinBox, but without limits
#and with better sizing (i.e. prefered size changes according to number)
#and with step changing according to cursor position

class FixedPointValidator(QtGui.QValidator):
    def __init__(self, *args):
        super().__init__(*args)
        pass
        
    def validate(self, text, pos):
        #print("FixedPointValidator.validate", text, pos)
        try:
            text = text.replace(",", ".")
            value = float(text)
            if "e" in text.casefold():
                #now it has to change but the one to change it is FixedPointLineEdit, because it knows how to format the number from value
                pass
            return (self.Acceptable, text, pos)
        except Exception as e:
            if len(text) == 0:
                return (self.Intermediate, text, pos)
            if text == "-":
                return (self.Intermediate, text, pos)
            print("FixedPointValidator.validate failed:", e)
            return (self.Invalid, text, pos)
        pass
        
    def fixup(self, *args):
        print("FixedPointValidator.fixup", args)
        return super().fixup(*args)
    pass

class DoubleWidget(SlowLineEdit):
    #for editing floating numbers in fixed point notation
    #enforce python formating of numbers
    #TODO: could place small spaces every three digits (and skip them when moving cursor and eliminate tham when copying the number)
    #TODO: allow pasting scientific notation - you can, but it is not translated to fixed point before first stepBy and that has bad increment
    #do not try complex numbers, part of the code might handle it without failing, but overall the result is not defined
    
    #problem is that it will emit also all kinds of signals that are not the final result/value
    
    valueChanged=QtCore.pyqtSignal(float)
    #we want to emit valueChanged signal when 
    #   - text is edited and timeout occured and new text is a valid value 
    #   - when setValue changes the value
    #   - when arrow key or wheel calls stepBy and value is changed - this calls setValue so the path is the same as second case
    
    
    def __init__(self, *args):
        super().__init__(*args)
        validator=FixedPointValidator()
        self.setValidator(validator)
        self._value = 0 #last valid value that was entered to the editor
        #text of the editor should alway correspond to _value unless it is being edited, once editing is done it should be validated and either _value or text changed to match
        self._precision = 0
        self._minimum = None
        self._maximum = None
        self._updateTextFromValue()
        pass
    
    def _timeout(self):
        #print("DoubleWidget._timeout")
        #SlowLineEdit timeout has to be first to properly handle background changing during timeout
        super()._timeout() #note that this will still emit textEditedSlow with the wrong text, but we can leave that be
        #we have to check that entered text is a valid value
        self._validateTextAndSetValue()
        pass
    
    def stepBy(self, steps):
        if steps==0: return

        text=self.text()
        pos=self.cursorPosition()
        
        value = float(self._value) #copy
        
        #change value according to steps and cursor position
        dpos=text.find(".")
        
        if dpos==-1 : dpos=len(text)
        
        spos=1 if text[0]=='+' or text[0]=='-' else 0

        if pos<spos:
            #do not change if positioned at '+' or '-' at the start
            pass
        elif pos<=dpos:
            value+=steps*10**(dpos-pos)
        elif pos>dpos and pos<=len(text):
            value+=steps*10**(dpos+1-pos)
        else:
            value+=steps #this should be reasonable (e.g. wheel when pos is at the end)

        #print("DoubleWidget.stepBy", value, self._value)
        #we need to validate the range
        self._validateValueAndSetText(value, max(0,len(text)-dpos-1))
        newtext = self.text()
        self.setCursorPosition(len(newtext)-len(text)+pos)#keep the same position from the end
        #no need to emit this - self.textEditedSlow.emit(newtext)
        pass
        
    def setMinimum(self, minimum):
        self._minimum = minimum
        if minimum is not None and self._value<minimum:
            self._validateValueAndSetText(minimum, self._precision) #TODO: self._precision might not be adequate after changing value (e.g. setting minimum of "0.01" to "0" valued DoubleWidget)
        pass
        
    def setMaximum(self, maximum):
        self._maximum = maximum
        if maximum is not None and self._value>maximum:
            self._validateValueAndSetText(maximum, self._precision)
        pass

    def setValue(self, value, precision = None, contract=False):
        """Be carefull: if precision is not given, DoubleWidget will guess it and round off the value with that precision.
        """
        assert(value is not None)
        self._validateValueAndSetText(value, precision, contract)
        pass
    
    def value(self):
        return self._value
        
    def precision(self):
        return self._precision

    #inner works
    def _updateTextFromValue(self):
        #print("DoubleWidget._updateTextFromValue", self._value, self._precision, "{:.{p}f}".format(self._value, p=self._precision))
        self.setText("{:.{p}f}".format(self._value, p=self._precision))
        pass
        
    def _adjustToRange(self, value, precision):
        if self._minimum is not None and value<self._minimum : return self._minimum
        if self._maximum is not None and value>self._maximum : return self._maximum
        return value
        
    def _setValue(self, value, precision, forceUpdateText=False, contract=False):
        #print("DoubleWidget._setValue", value, precision, forceUpdateText, self._value, self._minimum, self._maximum)
        if precision is None: 
            #precision = self._precision
            #we have to guess the correct precision
            #the value we hold is limited by the number of digits we show in editor, if someone tries to set value from outside
            #and does not set the precision we should reproduce the number as well as possible
            #which is basically presicion 15 - trailing zeroes (round off errors are at 1e-16)
            #lets try starting from 13
            #or just look for a large bunch of zeros
            #print(value)
            text = "{:.{p}f}".format(value, p=13).strip("0")
            dpos = text.find(".")
            precision = 0 if dpos==-1 else len(text)-dpos-1
            pass
        #adjust range
        adjusted = self._adjustToRange(value, precision)
        #update text if value is adjusted, or if update is forced
        updateText = forceUpdateText or adjusted is not value
        #emit signal only if value changed within precision
        emit = abs(self._value - adjusted) >= 0.5*(10**-precision) #TODO: this is a bit controvertial - adjusting text of edit might not change value - TODO: revert to old text in that case!
        #print("\t", emit, self._value, adjusted, abs(self._value - adjusted),  10**-precision)
        #print("DoubleWidget._setValue", value, precision, emit, self._value, adjusted, abs(self._value - adjusted), 10**-precision)
        self._value = adjusted
        self._precision = precision
        
        if contract:
            text = "{:.{p}f}".format(value, p=self._precision).strip("0")
            dpos = text.find(".")
            self._precision = 0 if dpos==-1 else len(text)-dpos-1
            
        
        if emit: self.valueChanged.emit(adjusted)
            
        if updateText: self._updateTextFromValue()
        pass
        
    def _validateValueAndSetText(self, value, precision=None, contract=False):
        try:
            self._setValue(value, precision, True, contract=contract) #force update text 
        except Exception as e:
            #do not change self._value
            #do not change text
            #TODO: complain a lot
            print("DoubleWidget._validateValueAndSetText", value, precision)
            import traceback
            print(traceback.format_exc())
            pass
        pass
        
    def _validateTextAndSetValue(self):
        #print("DoubleWidget._validateTextAndSetValue")
        try:
            #TODO: we could possibly here change decimal separator (instead of using validator) and replace scientific notation for fixed point
            text = self.text()
            value = float(text)
            #precision=sum(map(str.isdigit, text.lstrip('+-0').lstrip(",."))) #this calculates number of significant digits
            dpos=text.find(".")
            precision = 0 if dpos==-1 else len(text)-dpos-1
            #print("DoubleWidget._validateTextAndSetValue: setValue", value, precision)
            self._setValue(value, precision, False) #do not update text if not adjusted
        except Exception as e:
            #reset text to the last valid value
            print("DoubleWidget._validateTextAndSetValue: could not validate", e)
            self._updateTextFromValue()
        pass
        
        
    def wheelEvent(self, we):
        self.stepBy(we.angleDelta().y()/120.)
        pass

    def keyPressEvent(self, kpe):
        #filter out up/down
        if kpe.key()==QtCore.Qt.Key_Up:
            self.stepBy(1)
            kpe.accept()

        if kpe.key()==QtCore.Qt.Key_Down:
            self.stepBy(-1)
            kpe.accept()
        
        super().keyPressEvent(kpe)
        pass
    pass

class IntWidget(SlowLineEdit):
    #for editing floating numbers in fixed point notation
    #enforce python formating of numbers
    #TODO: could place small spaces every three digits (and skip them when moving cursor and eliminate tham when copying the number)
    #TODO: allow pasting scientific notation - you can, but it is not translated to fixed point before first stepBy and that has bad increment
    #do not try complex numbers, part of the code might handle it without failing, but overall the result is not defined
    
    #problem is that it will emit also all kinds of signals that are not the final result/value
    
    valueChanged=QtCore.pyqtSignal(int)
    #we want to emit valueChanged signal when 
    #   - text is edited and timeout occured and new text is a valid value 
    #   - when setValue changes the value
    #   - when arrow key or wheel calls stepBy and value is changed - this calls setValue so the path is the same as second case
    
    
    def __init__(self, *args):
        super().__init__(*args)
        validator=QtGui.QIntValidator()
        self.setValidator(validator)
        self._value = 0 #last valid value that was entered to the editor
        #text of the editor should alway correspond to _value unless it is being edited, once editing is done it should be validated and either _value or text changed to match
        self._minimum = None
        self._maximum = None
        self._updateTextFromValue()
        pass
    
    def _timeout(self):
        #we have to check that entered text is a valid value
        self._validateTextAndSetValue()
        return super()._timeout() #note that this will still emit textEditedSlow with the wrong text, but we can leave that be
    
    def stepBy(self, steps):
        if steps==0: return
        steps = int(steps)
        text=self.text()
        pos=self.cursorPosition()
        
        value = int(self._value) #copy
        
        #change value according to steps and cursor position
        dpos=len(text)
        
        spos=1 if text[0]=='+' or text[0]=='-' else 0

        if pos<spos:
            #do not change if positioned at '+' or '-' at the start
            pass
        elif pos<=dpos:
            value+=steps*10**(dpos-pos)
        else:
            value+=steps #this should be reasonable (e.g. wheel when pos is at the end)

        #we need to validate the range
        self._validateValueAndSetText(value)
        newtext = self.text()
        self.setCursorPosition(len(newtext)-len(text)+pos)#keep the same position from the end
        #no need to emit this - self.textEditedSlow.emit(newtext)
        pass
        
    def setMinimum(self, minimum):
        self._minimum = minimum
        if minimum is not None and self._value<minimum:
            self._validateValueAndSetText(minimum)
        pass
        
    def setMaximum(self, maximum):
        self._maximum = maximum
        if maximum is not None and self._value>maximum:
            self._validateValueAndSetText(maximum)
        pass
        
    def setValue(self, value):
        self._validateValueAndSetText(value)
    
    def value(self):
        return self._value

    #inner works
    def _updateTextFromValue(self):
        self.setText("{:d}".format(self._value))
        pass
        
    def _adjustToRange(self, value):
        if self._minimum is not None and value<self._minimum : return self._minimum
        if self._maximum is not None and value>self._maximum : return self._maximum
        return value

    def _setValue(self, value, forceUpdateText=False):
        #adjust range
        adjusted = self._adjustToRange(value)
        #update text if value is adjusted, or if update is forced
        updateText = forceUpdateText or adjusted is not value
        #emit signal only if value changed
        emit = self._value != adjusted
        self._value = adjusted
        if emit: self.valueChanged.emit(adjusted)
        if updateText: self._updateTextFromValue()
        pass
        
    def _validateValueAndSetText(self, value):
        try:
            self._setValue(value, True)
        except Exception as e:
            #do not change self._value
            #do not change text
            #TODO: complain a lot
            print("IntWidget._validateValueAndSetText", value)
            print("\terror:", e)
            pass
        pass
        
    def _validateTextAndSetValue(self):
        try:
            self._setValue(int(self.text()), False)
            pass
        except:
            #reset text to the last valid value
            self._updateTextFromValue()
        pass
        
        
    def wheelEvent(self, we):
        self.stepBy(we.angleDelta().y()/120.)
        pass

    def keyPressEvent(self, kpe):
        #filter out up/down
        if kpe.key()==QtCore.Qt.Key_Up:
            self.stepBy(1)
            kpe.accept()

        if kpe.key()==QtCore.Qt.Key_Down:
            self.stepBy(-1)
            kpe.accept()
        
        super().keyPressEvent(kpe)
        pass
    pass

   
if __name__=="__main__":
    def slot(*args):
        print("slot", *args)

    import sys
    app = QtWidgets.QApplication(sys.argv)
    window=QtWidgets.QMainWindow()
    #DoubleWidget
    #w=DoubleWidget()
    #w = IntWidget()
    w = IntWidget()
    w.setMinimum(5)
    w.setMaximum(10)
    #w.setMaximum(None)
    
    #SlowLineEdit
    #w=SlowLineEdit()
    #w.textEditedSlow.connect(slot)
    w.valueChanged.connect(slot)
    
    window.setCentralWidget(w)
    window.show()
    res=app.exec_()
    w.setParent(None)
    del w
    print("konec")
    sys.exit(res)
    
