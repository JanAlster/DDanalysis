#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
PQR package, Qt-based plotting widgets 
Copyright (C) 2019  Jan Alster 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""
import faulthandler
faulthandler.enable()

from PyQt5 import QtCore, QtWidgets, QtGui, QtSvg
from PyQt5.uic import loadUi

import sys
import os
import os.path as op

from collections import OrderedDict


PQR_VERBOSE = False
def vprint(*args, **kw):
    global PQR_VERBOSE
    if PQR_VERBOSE: print(*args, **kw)
"""
plotting widget

handling of non-equispaced axes

use cases

spectral plot with one horizontal axis in nm and another one (zooming the same) in some frequency units (obviously non-equispaced)

time plot of decay with very different time steps, but not suitable for logarithmic scale

spectral plot with common horizontal axis and independent (many) vertical axes (separate zooming, but possibly with common zero level)

the same axis (including zoom) on two different plots

polar axis?


inner works

axis can be mapped (non-linearly, but uniquely and reversibly) to another
simple case is pan/zoom 
another example is nm to icm conversion
yet another enumerated axis to pixels

AxisX is default horizontal axis of a plot, always has range (0,1) which will be stetched to whatever pixels are available for the widget
Axis1 is some axis with defined units and some range (xmin, xmax) which will be mapped
Axis1.parent = AxisX


in effect
 - axis stack should transform coordinates of data to coordinates of AxisX and AxisY (which in fact need not be real, it can be just a concept)
 
Axis plotter can be used/displayed on different plots (meaning it can have more than one "widget")
 
"""


"""
version  2

Axis: units to linear scale (not uniform, not linear transform, whatever inversible)
Data: e.g. XY scatter series

DataPlotter(Data, Axis, Axis) -> gritems on linear scale, adjustable graphics, possibly adjustable data (not needed to replot all)

AxisArea: part of plot, keeps zoom/pan range, works on Axis to get ticks values, converts Axis lineaer scale to (0, 1) DataArea range (scaling gritems)

cannot use gritems on multiple plots, because each gritem can be on scene only once and can be only on one scene
DataPlotter might keep tracks (or offer signals) to modify (visuals of) multiple generated gritems sets 
axis cannot convert directly to DAtaArea units, bacuase we would need to recreate gritems every time we pan or zoom, whereas it should be faster to just move and scale gritems

TODO:
 how in this tie two axes to zoom/pan together?
   - basically two AxisAreas need to have the same transformation from linear to dataarea, but not from units to linear (each can have different transform)
   - AxisArea can keep two transforms, one from units to linear and second from linear to dataarea
      - first will usually be identity
      - second will usually be Zoom and pan
      - zoom and pan can be tied together by sharing second transform, which will have internal settings for transformation shared too (can be even shared across multiple figures)
      AxisArea will need to keep track what is the data range shown on axis as second transform itself does not have access to this information; 
         it also has to translate user requests for axis range from units to params of Zoom transform
      if plotter first transform returns linear in (0,1) range and it should not be scaled/moved, the second transform can be identity
      the second transform has to be linear (otherwise moving and scaling gritems will not work) so it is limited to identity or pan and zoom
      
""" 

"""
TODO:
 - alpha channel for DD plots - use case: partial overlap of 2DES two color experiments with multiple color combinations
 
 - we have a problem, if we want to show two DDPlotters using NumEnumTransform to make things equispaced
    - if second is inside the area of first and they do not have matching axes, then it can be deformed (it will not be equispaced on the NumEnumTransform of the first)
    - if they are side by side, we can add stops from both, but what about the space between them?
    - what if there is partial overlap?
    - we can squeeze/stretch contours as needed, but not pixmap - very easily the pixel might not be on the correct position (which is what pqr is about, but with current implementation it only works for single DDPlotter)
    
    for now, I cannot plot two DDPlotter on one Plot correctly, and contours need to be limited to the range of original NumEnumTransform
    I am afraid that general solution requires pixmap with variable length pixels... (e.g. lots of RectItems or custom Painter)    
"""


#TODO: if we have two connected axes (ie. nm and freq) and we plot some data against one and some against the other, zoomToData will not work properly, because the two axes does not share their dataRanges, but they will influece each other zoom ranges
#   what is even worse is that order of axis zooming is arbitrary (dict) so it is not constant which axis will go last and will force its data range
#   intended behaviour in this case is to cover all data on all connected axes 
#   probably connecting axes should not be done just by linking zoomrange signals, but needs to be deeper

#TODO: we might specify positions as (x+y*1j) where x would be absolute distance in pixels and y*1j distance in plot area coordinates (i.e. 0.5j would be middle of plot area, 10+0.5j would be 10 pixels (or other units) from plot arrea middle - this would be good to position ticks and other elements which should not scale with plot size 
#    the absolute might be in fact relative to figure size (e.g. tick size should change with figure size)
# I should probably also make figure size independent from widget size (or underlying pixmap size - or the bigger from its width and height) so that text size scales properly regardless of widget resize (this should be possible without scaling dataarea, but the will make relative positioning more difficult and it will not work with cosmetic pens)
# anyway current approach is not optimal

#TODO: auto adjust format of ticks depending on smallest tick step

#TODO: allow axis settings via UI (tick number, tick format) and on exit (optionally) output setup string which will be able to recover this
#  that way it can be used even for simple scripts where you can copy paste this into the script easily

#TODO: sizing of multiple plots is weird, they do not get the same amount of space

#TODO: allow zoom to scene, i.e. enlarging selected plot (perhaps Ctrl+wheel?)

#TODO: DDPlotter with alpha should not use colormap with white color, otherwise it can be confusing

#TODO: linked tracers which move to the same data coordinate across multiple plots

import numpy as np

from enum import Enum

try:
    from params import Params, Unset
except:
    try:
        from pqr.params import Params, Unset
    except:
        from .params import Params, Unset
#Params should be used 
# for attributes that user should be able to change once the Figure is complete to tune the visuals
# only things directly tied to UI should be saved

#can we do global pen for every line?
globalPen = QtGui.QPen()
#problem is that what is already done will not update with change to this pen
# at least I think that Qt lines etc keep a copy of the pen, no? yes :(
# that means that we cannot dynamically change thickness of pens with widget resize and are stuck with cosmetic pens (which cannot show true line border so mouse clicking is not working)
# ~ globalPen.setCosmetic(True)
globalPen.setWidthF(2.)
#line thickness is not preserved when saving image with subpixel
#basically we need to determine correct line width before plotting anything
#and it will not scale nicely with the image
#maybe we could do something with the PlotArea instead? not use scaling so the we do no need cosmetic pens?
# but how to preserve coordinates mapping?
#we could register every line to have its pen changed when needed and redrawn, but that seems like overkill
# it is not unreasonable to set pen before plotting
# but it has to scale properly with change of widget size
# i.e. we need to set line width relative to widget size
#   graphics scene could handle scaling, but it must keep aspect ratio, or it has to replot
# we could keep line widths in relative size of min(width, height)
#   but we would have to redraw everything with new (pixel) line width with every widget size change
#   or we would need to scale heavily - always draw scene such that min(width, height) is 1; that would change the other scene size with widget resize, but it will keep line widths and there would be not need to change pens, but we would have to redraw (at least with every size change that changes aspect ratio)
# we can set scene (and its graphicswidget layout container) to specified size (in user units in which all other dimensions would be set) and use graphicsview.fitinview(scene.boundingrect, keepaspectratio) to handle resizing (will keep acpect ration, which is not that bad - we can even limit widget size to certain aspect ratio - that will not work very well with tiling wm - we could only set minimum scene size and calculate the other one acording to view aspect ratio - or we could update scene size with every view size/aspect change - that will be slower)
# with that  at least plot lines would not need to be cosmetic and would scale properly (width will be set in scene coordinates) - but they would have to be outside of PlotArea._area and would have to redrawn if AxisArea changes size  or could it ignore width change from just parent transform? but not scene transform?
#   something like tlabel.setFlag(tlabel.ItemIgnoresTransformations, True) - but that could ignore all transformations and I want the lines to scale size with AxisArea
#problem is that we need two types of units
#  - tick position is relative to AxisArea size and should scale with it 
#  - but tick length is relative to scene size and should not scale with AxisArea size
# - we could hardcode some things, but user might want to use this too (e.g. to place legend either to specific fixed position from plot edge (moving when the plot edge moves), or relative to plot size (e.g. in the middle))
#  - but how to handle that, and how to force Qt to use it...
#  - I could do (absolute + 1j*relative) to input the coordinates, but Qt objects cannot handle those
#   I would need to do my own objects that would just keep coordinates for Qt objects and do the transform to scene units
#  or do e.e. anchors (which would be in relative PlotArea units, will not be visible) and visible objects would be placed absolute distance from anchors (will have to remember the anchor or the anchor would have to remember all objects that use it, and move them if itself moves - then the object has to keep the absolute distance from anchor but also has to be drawn in scene coordinates - and I guess graphicslineitem cannot have end points each parented to different anchor - if that worked it would solve my problems somewhat)
# data area would have one anchor as origin of pan; or it could be left as is, there is not point of data lines to have different anchor.. and pan/zoom works fine as is (maybe line width will be a problem)

# line width of data lines should scale with scene scale (view scale) but not with zoom
#we should be moreorless able to preserve dataline width if it is cosmetic pen set to specific width, but the scene size (initial before anything is drawn) should be close to actuall used size, i.e. the viewport should not scale - because dataline width will not scale with viewport

#we could do setSceneSize(newsize, scaleParams=True) which wil rescale all default params like line width and tick size to new size (same proportion to old size)

#TOOD: adjust also font size to match scene size fraction
#   and force different labels to keep the size



#global settings
#all measures floating point precision and the same units
# if you put GS.width = 100 and assume it is 100mm, then tickMinorSpacing of 2 will be in mm too
# when saving, you can just scale it to whatever pixels you want (although default is GS.width pixels)
# (changes to settings will affect all subsequent drawing, it will not change what is already drawn)
class GS(object):
    def __init__(self):
        super().__init__()
        #this set is for screen (in px)
        # scene size 
        self.width = None
        self.height = None
        #this will tell FigureWidget to keep its scene on the same size as itself
        # but it will not play nice with saving images without showing the widget first - you have at least set the widget size before trying to save image
        # and it will not scale lines - line width is then in pixels

        # plot
        self.plotLinePen = QtGui.QPen() #axis baseline, major ticks
        self.plotLinePen.setWidthF(1)
        # self.plotLinePen.setColor(QtGui.QColor("red")) #you can change color
        #plot lines are assumed to be the same, so they have single pen; data lines tend to need individual pens so only width is uniform
        self.dataLineWidth = 1 #WARNING: this will not scale with the rest, you should put the width that will look good in the final resolution
        # also note that setting width of data line when plotting will use multiple of this value
        # also note that values > 1 might make plot harder to render

        self.axisLabelSize = 7 #font size for axis label
        self.plotTitleSize = 7 #font size for plot title

        # ticks
        self.tickMinorSpacing = 2 # from tick to tick labels
        self.tickMajorSpacing = 5 # from tick labels to axis label
        self.tickLength = 5 
        self.tickLabelSize = 7 #font size for tick labels
        #NOTE: tick label format is not part of GS and can be changed even after tick labels are drawn
        
        self.colorBarWidth = 10
        pass
    
    def loadPreset(self, preset):
        #TODO: this could be handled better, but should suffice for now
        if preset == "scaling":
            self.width = 600
            self.height = 500
        elif preset == "AIP1":
            #set for AIP single column
            #when saving image use saveImage(..., width=4016) to get 600dpi or saveImage(..., width=2008) to get 300dpi (should be properly scaled)
            """
            output
            black and white 600dpi

            color online 300dpi

            1col width 8.5cm
            2col width 17cm

            max height 21.1cm

            line width min 0.5pt (0.18mm)
            text height min 2.8mm

            EPS, TIFF, PDF, JPEG 

            17cm at 600dpi -> 4016px
            8.5cm at 600dpi -> 2008px

            17cm at 300dpi -> 2008px
            8.5cm at 300dpi -> 1004px

            """    
            self.width = 85 #mm
            self.height = 211 #mm, this is maximum height, should be updated to correct one
            
            self.plotLinePen.setWidthF(0.2) #mm
            
            self.dataLineWidth = 0.2*4016/85 #for 600dpi
            # ~ self.dataLineWidth = 0.2*2008/85 #for 300dpi
            
            self.tickLabelSize = 2.8 #mm
            self.tickLength = -2 #mm 
            self.tickMinorSpacing = 1 # from tick to tick labels
            self.tickMajorSpacing = 1 # from tick labels to axis label    
            
            self.axisLabelSize = 2.8 #mm
            self.plotTitleSize = 2.8 #mm
            
            self.colorBarWidth = 2
        elif preset=="AIP2":
            self.loadPreset("AIP1")
            self.width = 170 #mm
        else:
            print("WARNING: unknow preset", preset)
        pass
    pass

GS = GS()

#helper
class ColorCycle(object):
    """
    Simple way how to manage different colors for individual graphs.
    """
    _colorPool = ["blue", "red", "green", "purple", "olive", "magenta", "orange", "violet", "darksalmon", "crimson", "mediumblue", "limegreen", "navy", "darkyellow", "chocolate", "slateblue", "indianred", "teal", "lawngreen", "tomato", "sienna", "royalblue", "fuchsia"]
    def __init__(self, index=0, color_pool=None):
        super().__init__()
        if color_pool is not None:
            self._colorPool = color_pool
        self._index=index
        pass
        
    def __call__(self):
        color = self._colorPool[self._index]
        self._index = (self._index+1) % len(self._colorPool)
        return color
        
    def prev(self):
        return self._colorPool[self._index-1]
        
    previous = prev
    
    def next(self):
        return self._colorPool[(self._index+1) % len(self._colorPool)]
        
    def reset(self):
        self._index=0
        pass
    pass

color = ColorCycle() #for convenience

def cmap_color_cycle(color_map:"ColorMap", size):
    color_pool = color_map(np.r_[0:1.:size*1j])
    return ColorCycle(color_pool=color_pool)

#we might need to derive this from QObject to enable signals
class Transform(object):
    """
    this should handle transform from data space (non-uniform) to linear space
    can be linked to another Transform (linear space of this Transform is data space of parent Axis)
    
    Transform need not be tied to specific Plot or Figure, since it only describes the transformation of coordinates
    
    it can possibly keep default params for AxisArea (which will handle zoom/pan, ticks, labels, data range, etc)
    """

    def __init__(self, parent=None):
        self.parent = parent

    def data2linear(self, data, single=False):
        #go through the axis chain to calculate the ultimate values for plot view coordinates
        if single or self.parent is None:
            return self._data2linear(data)
        else:
            return self.parent.data2linear(self._data2linear(data))
        
    def linear2data(self, linear, single=False):
        if single or self.parent is None:
            return self._linear2data(linear)
        else:
            return self._linear2data(self.parent.linear2data(linear))
        
    def _data2linear(self, data):
        return data
        
    def _linear2data(self, linear):
        return linear
    pass

#Todo: this is not tested, I cannot think of a good use case
class ZoomTransform(Transform):
    """
    transforms its own data coordinates (in units) to units of parent (None parent indicates directly the plot view with range 0,1)
     - mapping is "min" to parentMin and "max" to parentMax (note that if parent is not plot view, then displayed min and max values of the axis will not be min and max)
     
    this should not be used for zooming data on DataArea
    only use this if you want to scale/offset two transforms 
    """
    #start and parentStart are aligned
    #scale is kept in range/parentRange ratio (basically we have too many parameters here, but it should be more convenient to setup)
    def __init__(self, start=0., range=1., parentStart=0., parentRange=1., *args, **kw):
        super().__init__(*args, **kw)
        self.start = start
        self.range = range
        self.parentStart = parentStart
        self.parentRange = parentRange
        pass

    #~ def setup(self, **kw):
        #~ vprint("ZoomAxis.setup", kw)
        #~ return super().setup(**kw)

    def zoom(self, relchange, point=None):
        #we need ability to relegate zoom to parent Axis
        #  i.e. EnumAxis should not zoom, it should have (hidden) ZoomAxis underneath
        
        #two kinds of zooming
        # zoom around point - keep specified point fixed and extend/shrink range on both sides
        # - now the point can be either mouse, or fixed to specific level (so that e.g. two vertical axes can have different scales  but common zero level)
        # it is probably resonable to keep the start point fixed
        
        #point is in d, if None use self.start
        #relchange is multiplier for self.range
        if point is None: point = self.start
        
        r = relchange * self.range
        start = point - relchange*(point-self.start)
        
        self.setup(start=start, range=r)
        pass

    def _data2linear(self, d, single=False):
        k = self.parentRange/self.range
        l = self.parentStart - k*self.start
        return k*d+l
        
    def _linear2data(self, p, single=False):
        k = self.range/self.parentRange
        l = self.start - k*self.parentStart 
        return k*p+l
    pass


#might be usefull 
class Reverse(Transform):
    def _data2linear(self, d):
        return -d
    def _linear2data(self, p):
        return -p
    pass

class MZMTransform(Transform):
    #WIP
    """
    Max-Zero-Min transform, bilinear transform that keeps 0data to 0linear and scales positive and negative separately
    """
    def __init__(self, min=-1, max=1., *args, **kw):
        super().__init__(*args, **kw)
        self.min = min
        self.max = max
        
    def _data2linear(self, d):
        D = np.full_like(d, self.max)
        D[np.where(d<0)] = -self.min
        return d/D
        
    def _linear2data(self, p):
        D = np.full_like(p, self.max)
        D[np.where(p<0)] = -self.min
        return p*D

class EnumTransform(Transform):
    #set of arbitrary elements assigned indexes which serve as units
    #note that elements are part of axis, not data (i.e. data shoudl contain element indices - although this might change)
    def __init__(self, elements = [], *args, **kw):
        super().__init__(*args, **kw)
        self.elements = elements
        
    def _data2linear(self, d):
        #here d should be integer
        return np.floor(d)+0.5
        
    def _linear2data(self, p):
        return np.floor(p)
    
    #TODO: this will need some special tick handling
    pass

class NumEnumTransform(Transform):
    #set of ordered numerical elements (but with no requirement on their spacing) which will behave as Enum, but can be linearly interpolated in between elements
    #stops should be ordered (ascending) np.array of numbers with arbitrary distances - at least 2 long
    
    def __init__(self, stops=np.array([0., 1.]), *args, **kw):
        super().__init__(*args, **kw)
        self.stops = stops
    
    def _data2linear(self, d):
        #for best performance d should be inside the range specified by stops
        N = len(self.stops)
        #anything at 0 shoudl be set to min or interpolated below
        #anything at N should be set to max or interpolated above
        #~ vprint("NumEnumTransform", self.stops)
        #~ vprint("\t", d)
        #~ vprint("\t", np.interp(d, self.stops, np.arange(N)))
        #~ return np.interp(d, self.stops, np.arange(N)/(N-1))
        return np.interp(d, self.stops, np.arange(N))

    def _linear2data(self, p):
        N = len(self.stops)
        #~ return np.interp(p, np.arange(N)/(N-1), self.stops)
        return np.interp(p, np.arange(N), self.stops)
        
    pass

#TODO: possibly do something similar to NumEnumTransform but with smooth interpolation between points, likely spline (if it could be kept inversible)

"""    
class ArbitraryTransform(Transform):
    defaultParams = {whatever params you would like}
    paramsDefaults.update(Axis.paramsDefaults)#for each parent class (this is slightly annoying)
    
    def _d2p(self, d, single=False):
        retrun transform(params, d)
        
    def _p2d(self, p, single=False):
        return inverseTransform(params, p)
        
    pass
"""

class AOverTransform(Transform):
    """
    d = A/p
    p = A/d
    """
    #TODO: do some trick near zero to prevent exceptions
    def __init__(self, A = 1., *args, **kw):
        super().__init__(*args, **kw)
        self.A = A
        
    def _data2linear(self, d):
        return self.A/np.array(d)
    def _linear2data(self, p):
        return self.A/np.array(p)
    pass

class PolarTransform(Transform):
    ...


"""
for plotting we need a widget (Figure widget)
probably containing QGraphicsView looking at figure scene containing plots, title, etc..

we want plots and stuff be organized into grid, so we need it to be a 

QGraphicsView     FigureWidget
QGraphicsScene    figure scene 
QGraphicsWidget::setLayout    containing one organizing widget 
QGraphicsGridLayout     everything will be put into this layout
 

Plot shoudl be an QGraphicsWidget, to be includable into the grid layout

it will have its own QGraphicsGridLayout
centerred around a separate QGraphicsView looking on data scene
    or we can make this dataview a simple QGraphicsWidget and just clip everything?
with AxesWidget all around 
"""

class Anchor(QtWidgets.QGraphicsObject): #cannot use QGraphicsItem since it is not QObject and cannot sent signals
    """
    position is defined as absolute+1j*relative
    where absolute is in scene (or item) coordinates and relative is scaled to map 0,1 to parent size
    
    Qt position is masked and set privately
    
    it has to keep track of parent resize to update its real position acording to self._x, self._y
    
    its parent probably has to be QGraphicsWidget
    """
    
    #TODO: we could allow for Anchor parent be another Anchor, from where it would start its position (basically Anchor would be a vector, from starting position (defaults to 0,0) wit hrealtive coord taken from common parent; but we can probably just add two anchors (or add offset to Anchor to get the same result))
    
    posChanged = QtCore.pyqtSignal(QtCore.QPointF)
    
    def __init__(self, x=0, y=0, parent=None, **kw):
        color = kw.pop("color", None)
        super().__init__(None, **kw)
        self._x = x
        self._y = y
        
        #we have to set parent item by hand, super().__init__ will not call setParentItem
        self.setParentItem(parent)
        self._updatePos()
        
        #cannot do this: 
        self.setFlag(self.ItemHasNoContents, True)# , because it will have children. Or is it related to this item only?
        #for debug
        # ~ size = 7
        # ~ self._rect = QtCore.QRectF(-size/2, -size/2, size, size)
        # ~ self.marker = QtWidgets.QGraphicsEllipseItem(self._rect,self)
        # ~ self.marker.setParentItem(self)
        # ~ if color is not None:
            # ~ self.marker.setBrush(QtGui.QBrush(QtGui.QColor(color)))
        pass
    
    def paint(self, *args):
        pass
    
    def boundingRect(self, *args):
        return self.childrenBoundingRect()
    
    def setParentItem(self, parentItem):
        # ~ print("Anchor.setParentItem", parentItem)
        old = self.parentItem()
        if old is not None:
            #unregister
            old.geometryChanged.disconnect(self._updatePos)
        if parentItem is not None:
            #register
            parentItem.geometryChanged.connect(self._updatePos)
        return super().setParentItem(parentItem)
    
    #TODO: this should be called setAnchorPos and setPos (to directly set item pos) should be disabled
    def setPos(self, x, y):
        self._x = x
        self._y = y
        self._updatePos()
        pass
    
    #this is questionable, it might be better to return calculated/actuall position
    # ~ def pos(self):
        # ~ return self._x, self._y 
    
    def anchorPos(self):
        return np.array([self._x, self._y])
        
    def setX(self, x):
        self._x = x
        self._updatePos()
        
    def setY(self, y):
        self._y = y
        self._updatePos()
    
    def _getParentSize(self):
        parent = self.parentItem()
        if parent is None:
            scene = self.scene()
            if scene is None:
                size = QtCore.QSizeF(1,1)
            else:
                size = scene.size()
        else:
            size = self.parentItem().size()
        return size
    
    def _updatePos(self, *args):
        # ~ vprint("Anchor._updatePos", args)
        size = self._getParentSize()
        super().setPos(self._x.real+self._x.imag*size.width(), self._y.real+self._y.imag*size.height())
        # ~ vprint("\t\tparent", self.parentItem(), self.parentItem().size(), super().pos())
        self.posChanged.emit(super().pos())
        pass
        
    def __getitem__(self, key):
        if key==0:
            return self._x
        elif key == 1:
            return self._y
        else:
            raise IndexError("Anchor.__getitem__: index", key, "is out of range")

    #this will create problems as intermediate anchors are children too and will not be destroyed, but are not kept in python
    # ~ def __add__(self, other):
        # ~ ox, oy = other
        # ~ return Anchor(self._x+ox, self._y+oy, self.parentItem())

    # ~ def __sub__(self, other):
        # ~ ox, oy = other
        # ~ return Anchor(self._x-ox, self._y-oy, self.parentItem())

    # ~ def __mul__(self, other):
        # ~ return Anchor(self._x*other, self._y*other, self.parentItem())
        
    # ~ def __truediv__(self, other):
        # ~ return Anchor(self._x/other, self._y/other, self.parentItem())
    
    # ~ __radd__ = __add__
        
    # ~ def __rsub__(self, other):
        # ~ ox, oy = other
        # ~ return Anchor(ox-self._x, oy-self._y, self.parentItem())
    
    # ~ __rmul__ = __mul__

#this is getting out of hand
        # maybe we should redesign axes to transform not to (0,1) but to externaly set (posmin, posmax) which will be updated during DataArea resize
class MovableRelativeAnchorRY(Anchor):
    """
    this can simulate DataArea._area to certain extend
    
    basically anything that tries to setPos (or setX, setY) will in fact set relative
    position within parent item, which will be updated when parent changes size
    
    pos should be within (0,1)
    
    reversed y
     i.e. qtpos should be as if y was 1-y
    
    """        
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw) #again this will call _updatePos but we have not yet prepared it
        self.setFlag(self.ItemIsMovable, True)
        #we need to keep correct _x and _y when we move pos by mouse
        #hopefully we will be the first to get change signal
        self.xChanged.connect(self._updateFromPos)
        self.yChanged.connect(self._updateFromPos)
    
    
    def _updatePos(self):
        #to prevent circular emiting
        self._flagSetting = True
        super()._updatePos() #this will change Qt pos
        self._flagSetting = False
    
    def setPos(self, x, y):
#        vprint("RelativeAnchorRY.setPos", x, y, self.pos())
        self._x = x*1j
        self._y = (1-y)*1j
        self._updatePos()
        pass
    
    def setX(self, x):
#        vprint("RelativeAnchorRY.setX", x, self.pos())
        self._x = x*1j
#        self._y = (1-self.y())*1j #theses are not updated automatically
        self._updatePos()

    def setY(self, y):
#        vprint("RelativeAnchorRY.setY", y, self.pos())
#        self._x = self.x()*1j #theses are not updated automatically
        self._y = (1-y)*1j
        self._updatePos()
        
#    def moveBy(self, dx, dy): #this is left as backdoor for Qtpos setting
#        self._x += dx*1j
#        self._y -= dy*1j
#        self._updatePos()
        
    #OriginDecorator needs to live on 0,1 space, but AnchorRectItem needs to live on Qt pos space
    def pos(self):
        return QtCore.QPointF(self._x.imag, 1-self._y.imag)
        
    def QtPos(self):
        return super().pos()
    
    def x(self):
#        w = self._getParentSize().width()
#        return super().x()/w if w>0 else 0
        return self._x.imag

    def y(self):
#        h = self._getParentSize().height()
#        return 1-super().y()/h if h>0 else 0
        return 1-self._y.imag

    def _updateFromPos(self, *args):
        #this depends on _x and _y matching internal pos
        #how ever it is possible that internal pos is changed in other way than setPos, seX, setY or moveBy (e.g.  mouse)
        #then call this
        if self._flagSetting: return
        
        pos = super().pos()
        size = self._getParentSize()
        w = size.width()
        h = size.height()
        x = pos.x()/w if w>0 else 0
        y = pos.y()/h if h>0 else 0
        self._x = x*1j
        self._y = y*1j
        self.posChanged.emit(super().pos())
#        vprint("MovableRelativeAnchorRY._updateFromPos", self, self._x, self._y, self.pos())
#


        
        
        
        
class AnchorLineItem(QtWidgets.QGraphicsLineItem):
    """I strongly doubt that graphicsLineItem can parent each point to different graphicsitem
    This one can.
    WIP: - this is likely to be more complex if Anchor do not share parentItem
    """
    def __init__(self, start, end):
        """
        start, end - Anchor from which position the line is drawn from/to.
        They should notify of position change to keep the line updated. 
        Does not have to be Anchor, but it might not work.
        """
        super().__init__(start.parentItem())
        self._start = start
        self._end = end
        try:
            self._start.posChanged.connect(self._update)
        except: #TOOD: which one is missing signal
            pass
        try:
            self._end.posChanged.connect(self._update)
        except: #TOOD: which one is missing signal
            pass
        if start.parentItem() != end.parentItem():
            print("AnchorLineItem.__init__ warning: start and end Anchor do not share parentItem. Things are likely to break.")
        self._update()

    def setLine(self, *args):
        #just to prevent QGraphicsLineItem from modiying the line directly
        print("AnchorLineItem.setLine warning: cannot set line directly")
        #TODO: possibly allow to change edge Anchors
        pass
    
    def _update(self, *args):
        s = self._start.pos()
        e = self._end.pos()
        super().setLine(s.x(),s.y(),e.x(),e.y())
    pass
        
class AnchorRectItem(QtWidgets.QGraphicsRectItem):
    """I strongly doubt that graphicsLineItem can parent each point to different graphicsitem
    This one can.
    WIP: - this is likely to be more complex if Anchor do not share parentItem
    """
    def __init__(self, topLeft, bottomRight):
        """
        start, end - Anchor from which position the line is drawn from/to.
        They should notify of position change to keep the line updated. 
        Does not have to be Anchor, but it might not work.
        """
        super().__init__(topLeft.parentItem())
        self._tl = topLeft
        self._br = bottomRight
        try:
            self._tl.posChanged.connect(self._update)
        except: #TOOD: which one is missing signal
            pass
        try:
            self._br.posChanged.connect(self._update)
        except: #TOOD: which one is missing signal
            pass
        if topLeft.parentItem() != bottomRight.parentItem():
            print("AnchorRectItem.__init__ warning: corner Anchors do not share parentItem. Things are likely to break.")
        self._update()

    def setRect(self, *args):
        #just to prevent QGraphicsLineItem from modiying the line directly
        print("AnchorRectItem.setRect warning: cannot set rect directly")
        #TODO: possibly allow to change edge Anchors
        pass
    
    def _update(self, *args):
        s = self._tl.QtPos()
        e = self._br.QtPos()
#        vprint("AnchorRectItem._update", s, e, e.x()-s.x(), e.y()-s.y())
        super().setRect(min(s.x(), e.x()),min(s.y(), e.y()),abs(e.x()-s.x()),abs(e.y()-s.y()))
    pass
            


class PlotArea(Params, QtWidgets.QGraphicsWidget):
    """
    convenience base class for pieces shown of PlotWidget
    - possible color background
    - has _area QGraphicsItem with coordinates scaled to (0,1)x(0,1)
        - so if you want to position something relative to area size (and not pixels), just make it child of _area
        - warning: all pens need to be cosmetic, or results will be probably ugly
        #TODO: we might just use _area to calculate coordinates of items in PlotArea coordinates (from input on (0,0,1,1) rect), but pens likely still need be cosmetic as we stretch plotters to achieve axis scaling
    """
    def __init__(self, *args, backgroundColor=None, name=None, **kw):
        vprint("PlotArea.__init__")
        super().__init__(*args, **kw)
        
        
        self._name = name
        
        if backgroundColor is not None:
            p = self.palette()
            p.setColor(p.Window, QtGui.QColor(backgroundColor))
            self.setPalette(p)
            self.setAutoFillBackground(True)

        # ~ self.setPreferredSize(1,1) #as far as I can tell, size of QGraphicsWidget does not correspond at all with children QGraphicsItems, so we will need to set it by hand
        #~ self.resize(1,1) #this has to be initialized, but cannot call resize direcly, because derived classes might not be ready for it before *their* __init__ finishes
        
        self.topLeft = Anchor(0, 0, self, color="red")
        self.topRight = Anchor(1j, 0, self, color="orange")
        self.bottomRight = Anchor(1j, 1j, self, color="green")
        self.bottomLeft = Anchor(0j, 1j, self, color="lime")
        vprint("\tPlotArea.bottomLeft.parentItem()", self.bottomLeft.parentItem())
        

class DataArea(PlotArea):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.setFlag(self.ItemIsMovable)
        # ~ self.setPreferredSize(1,1) #as far as I can tell, size of QGraphicsWidget does not correspond at all with children QGraphicsItems, so we will need to set it by hand
        #~ self.resize(1,1) #this has to be initialized, but cannot call resize direcly, because derived classes might not be ready for it before *their* __init__ finishes
        # ~ self.setMaximumSize(1,1)
        # ~ self.setTransform(QtGui.QTransform.fromScale(1, -1).translate(0, -1), False) 
        # ~ this is broken, we need to resize DataArea
        # ~ we need _area back
        
        self._area = QtWidgets.QGraphicsRectItem(0,0,1,1,self)
        p = self._area.pen()
        p.setCosmetic(True)
        p.setStyle(QtCore.Qt.NoPen)
        self._area.setPen(p)
        self.setPreferredSize(1,1)
        self._area.setTransform(QtGui.QTransform.fromScale(1, -1).translate(0, -1), False) 
        
        #Todo: hover could be optional
        self.setAcceptHoverEvents(True)
        pass
        
    def hoverMoveEvent(self, event):
        p = self.parentItem()
        p._hover(event.pos())
    
    def resizeEvent(self, ev):
        vprint("\t\tDataArea.resizeEvent")
        size = ev.newSize()
        self._area.setTransform(QtGui.QTransform.fromScale(size.width(), -size.height()).translate(0, -1), False) #this should keep local coordinates of self._area children in [0,0,1,1] spread across the whole data area, flip y axis
        #this does not work, it will scale the DataArea, we want to scale thing inside of it
        #or we would need to keep DataArea size to 1, 1
        """
        tady je problem
        pokud chci udrzet data na oblasti (0,0)-(1,1)
        tak musim tady skalovat
        ale to zrusi sirku pen, takze musi byt cosmetic, ale stejne nejde vyhledavat, protoze vyhledavaci alrogitmus nepocita s cosmetic pen
        
        data chci na 0-1 proto, aby se nemuselo vsechno prekreslovat, kdyz se zmeni velikost dataarea
        i kdyz ono by mozna nemuselo, jen se zmeni clip? a musel by se zmenit i range na osach, ale to zas nechci
        
        chci aby viditelne rozsah nebyl zavisly na velikosti data area, takze pri zmene velikosti dataarea musim preskalovat zobrazena data
        nebo prekreslit
        
        vybirani krivek mysi (a tedy vyhledavani car) je docela dulezite, pro right click - contex menu, takze nemuzu pouzit cosmetic pen
        
        takze si v podstate muzu vybrat
         - pouzit cosmetic pen a nemoct vybirat cary
         - nezachovat data range pri zmene velikosti dataarea - tohle nechci
         - vsechno prekreslit pri zmene data area
         
        pomoc Anchor jsem tenhle problem dostal z AxisArea
        ale v DataArea jeste porad pretrvava
         
        jedna moznost je nedavat area na 0-1, ale na neco bliz skutecne hodnote, treba 0-100 nebo 0-1000, pak se skalovani sirky pen tak moc neprojevi...
        - jenze kdyz budu hodne zoomovat, tak to nepomuze, cary budou tloustnout
        - z hlediska zachovani sirky car pri zoomovani muzu but pouzit cosmetic pen, nebo vsechno porad prekreslovat ...
        - ovsem pozor, mozna pomuze deviceTransform parametr z QGraphicsScene.itemAt
        
        v podstate musim pro DataArea zachovat _area
        a zvolit sirku pen pro cary v dataarea tak, aby byla hezka pri konecnem rozliseni, protoze se skaluje ...
        (problem s itemAt je castecne vyresen pouzitim velmi tenkeho pen pro pocitani tvaru SmoothPathItem, ktery ovsem taky skaluje, takze nebude fungovat stejne pro vsechna zvetseni)
        """
        
    def wheelEvent(self, ev):
        br = self.boundingRect()
        vprint("DataArea.wheelEvent", ev.pos(), ev.scenePos(), br, ev.pos().x()/br.width(), ev.pos().y()/br.height())
        p = self.parentItem()
        if p is not None:
            #TODO: if shift or alt is pressed modify ev.delta() for more or less zoom
            p._zoomDataArea(self.mapToItem(self._area, ev.pos()), ev.delta())
            # ~ x = ev.pos().x()
            # ~ y = ev.pos().y()
            # ~ br = self.boundingRect()
            # ~ p._zoomDataArea(QtCore.QPointF(x/br.width(), 1.-y/br.height()), ev.delta())
            # ~ p._zoomDataArea(ev.pos(), ev.delta())
        pass   
    
    def mousePressEvent(self, ev):
        #middle click (or doubleclick) to zoomToData
        if ev.buttons() == QtCore.Qt.MiddleButton:
            self.parentItem().zoomToData()
        else:
            super().mousePressEvent(ev)
        pass
   
    def mouseMoveEvent(self, ev):
        if ev.buttons() == QtCore.Qt.LeftButton:
            lp = self.mapToItem(self._area, ev.lastPos())
            p = self.mapToItem(self._area, ev.pos())
            self.parentItem()._panDataArea(lp, p)
        pass

    def contextMenuActions(self, ev):
        #specific for FigureScene, returns iterable of (name, callback) pairs which will be included in context menu
        return [("Add Note", lambda : self.addNote(ev)), ("Add Tracer", lambda : self.addTracer(ev)), ("Add Crosshair", lambda : self.addCrosshair(ev)), ("Add Diagonal", lambda : self.addDiagonal(ev)), ("Add Rect", lambda : self.addRect(ev))]
    
    def addNote(self, ev, *args):
        #~ vprint("DataArea.addNote", ev, args, ev.pos(), ev.scenePos())
        #~ self.parentItem().addNote(self.mapToItem(self._area, self.mapFromScene(ev.scenePos())))
        self.parentItem().addDecoratorFromDialog(TextDecorator, self.mapToItem(self._area, self.mapFromScene(ev.scenePos())), ev.screenPos())

    def addTracer(self, ev, *args):
        #~ self.parentItem().addTracer(self.mapToItem(self._area, self.mapFromScene(ev.scenePos())))
        self.parentItem().addDecoratorFromDialog(TracerDecorator, self.mapToItem(self._area, self.mapFromScene(ev.scenePos())), ev.screenPos())

    def addCrosshair(self, ev, *args):
        #~ self.parentItem().addTracer(self.mapToItem(self._area, self.mapFromScene(ev.scenePos())))
        self.parentItem().addDecoratorFromDialog(CrosshairDecorator, self.mapToItem(self._area, self.mapFromScene(ev.scenePos())), ev.screenPos())

    def addDiagonal(self, ev, *args):
        #~ self.parentItem().addTracer(self.mapToItem(self._area, self.mapFromScene(ev.scenePos())))
        self.parentItem().addDecoratorFromDialog(DiagonalDecorator, self.mapToItem(self._area, self.mapFromScene(ev.scenePos())), ev.screenPos())

    def addRect(self, ev, *args):
        #~ self.parentItem().addTracer(self.mapToItem(self._area, self.mapFromScene(ev.scenePos())))
        try:
            self.parentItem().addDecoratorFromDialog(RectDecorator, self.mapToItem(self._area, self.mapFromScene(ev.scenePos())), ev.screenPos())
        except:
            import traceback
            traceback.print_exc()


class AxisPosition(Enum): #TODO: Flag would be nice, but it is not available in my version of Python
    #values are chosen to allow easy comparison of positions: orthogonal == ((position1 * position2)<0), horizontal == (position>0), vertical == (position<0)
    #with aliases
    horizontal = H = h = X = x = 1
    bottom = B = b = 2
    top = T = t =3
    vertical = V = v = Y = y = -1
    left = L = l = -2
    right = R = r = -3
    
    def isHorizontal(self):
        return self.value>0
        
    def isVertical(self):
        return self.value<0
        
    def isOrtogonal(self, other):
        if self.__class__ is other.__class__:
            return (self.value * other.value) < 0
        return NotImplemented
    pass


"""
there will be multiple objects that are drawn
and have the same atributes (color, line width, ...)

so we might want to take that out and reuse

these need to be instance atributes
 but question is if we should have 
 line1.width = 5
 line1.color = "blue"
 
 line2.width = line1.width
 line2.color = line1.color
 
 (i.e. class Line(LineParams, ...))
 
 or
 
 line1.property.width = 5
 line1.property.color = "blue"
 
 line2.property = line1.property
 (i.e. Line().property = LineParams())
 
 or keep the properties as an object, but forward all individual 
 to using class methods? that sound horrible
 
 maybe we can have 
 all encompasing property
 properties and use that
 class Line(LineParams):
    ...
 
 line1 = Line()
 line1.width = 5
 line1.color="blue"
 
 line2=Line()
 line2.params = line1.params #will copy all params that are shared 
 
 
 keeping params as atribute will allow sharing params across different lines/markers/...
   and change params of already created objects... (but then we need to keep track of updating)
   
 
 
 another question is if we have something (e.g. SeriesPlotter)
  that is both Line and Marker, how will the conflicting params be handled?
  (shall the less important group of params be separated into atribute?)
  (or hardcode line and marker border the same params?)
  
  
  
  when some params changes, we need to have possibility of response (e.g. redraw)

    we have conflict of names with Params from old pqr ...
    (basically this is similar)
"""

"""
lets try this
class Attr

attrs are set in class, not on instance via __init__ addParam
can be grouped and group signal will be emited if any attr in group changes

basically attrs are python properties
with some automated creation  and groups 

so instead of 

class C:
   t = property(getT, setT, delT)
   
you would do

class C:
   t = Attr(group="...")
   
to get basic attr
or you could create full property in any way available
and do 

class C:
   t = property(...)
   
   hmm to nejde
   
   nechci duplikovat property uplne
   ale pokud nemam pristup k get a set tak nejsem schopny
   kontrolovat zmenu hodnoty - to totiz musi delat set
   
   pokud ale uzivatel bude chtit od get a set neco extra, tak narazim
   a to uz se dostavame k Params
   jenze Params je krok stranou, ty __getattr__ a podobne jsou hrozne
   bylo by dobre, kdyby Attr byly normalni property
   
   ale hodilo by se zachovat groups, onSet, onGet, mozna i preSet
   
   zasadni rozdil ale je, ze Params muzou vyuzivat get a set odjinud
   (tj z instance)
   zatimco property jen ze tridy (takze nemuzou vyuzivat treba UI, aspon ne primo)
   
   
   taky narazim na problem
   x = property(get, set)
   uvnitr property netusim jak se jmenuje x
   
   @property
   def x(self):
      ...
    uvnitr property sice vim, jak se jmenuje x, ale musim def x
    
    property("x", get, set)
     - musim __setattr__ coz je necitelne, navic je to na instanci

class Attr:
    Unset = object()
    
    @classmethod
    def basic(cls, group = None):
        value = cls.Unset
        
        def get(self):
            return self._value 
        
        def set(self, value):
            self._value = value
            #value change handling
            ale tady narazim, protoze vsechny takto vytoverene property budu sahat na stejne misto
            
        return property(get, set)
"""         
    

wrappertype = type(QtCore.QObject)
class MetaAttrs(wrappertype):
    def __new__(cls, future_class_name, future_class_parents, future_class_attrs):
        """
          Return a class object, changing some attributes to properties with signals
        """
        # ~ print("meta", cls, future_class_name, future_class_parents, future_class_attrs)
        
        _attrs = {} #list of attributes, use default class as starting point
        for parent in future_class_parents:
            if future_class_name!="Attrs" and issubclass(parent, Attrs):
                _attrs.update(**parent._attrs)
        
        new_attrs = {}
        for name, attr in future_class_attrs.items():
            if isinstance(attr, Attr):
                _attrs[name] = attr.default
                
                #attrChange signal
                signal = QtCore.pyqtSignal(object)
                new_attrs[name+"Changed"] = signal
        
                def get(self, name=name):
                    return self._attrs[name]
                
                def set(self, value, name=name, attr=attr): #we can save current values to defaults of kw
                    # ~ print("set", name, attr)
                    sv = self._attrs[name]

                    if sv != value:
                        self._attrs[name] = value
                        signal = getattr(self, name+"Changed")
                        # ~ print(signal, dir(signal))
                        # ~ print("test2", attr)
                        signal.emit(value) #this only works if order of classes to inherit is correct (Attrs first)
                        # ~ print("test2.5")
                        if "groupSignal" in attr.kw: 
                            #the signal object stored in Attr is unbound
                            # that cannot be emmited
                            # we can either pass name of the signal to Attr
                            # or find the object in future_class_attrs
                            for n, a in future_class_attrs.items():
                                if a is attr.kw["groupSignal"]:
                                    break
                            getattr(self, n).emit()
                        # ~ print("test3")
                
                new_attrs[name] = property(get, set)
            else:
                new_attrs[name] = attr
                
        new_attrs.update(_attrs = _attrs) #default values
        
        
        #derived class might want to add some attr to existing group
        #  in that case we should not create a new signal
        # ~ print(new_attrs)

        # let `type` do the class creation
        return super().__new__(cls, future_class_name, future_class_parents, new_attrs)
    

class Attr(object):
    """This will describe attributes to be processed by MetaAttrs"""
    def __init__(self, default, **kw):
        self.default = default
        self.kw = kw

# ~ class Attrs(QtCore.QObject, metaclass=MetaAttrs):
class Attrs(QtWidgets.QGraphicsObject, metaclass=MetaAttrs):
    # ~ _attrs = {}
    
    # ~ @classmethod
    # ~ def split(cls, kw):
        # ~ own = {}
        # ~ other = {}
        # ~ for it in kw:
            # ~ if it in cls._attrs:
                # ~ own[it] = kw[it]
            # ~ else:
                # ~ other[it] = kw[it]
        # ~ return own, other
    """
    doc?
    """
    
    
    """
    we cannot use classmethod for this because
    it needs cls as argument and that means that cls is defined and done at this point
    but Qt side will also be done and we cannot add more signals
    (unless we can update QMetaObject)
    so we need to create this without access to cls :(
    
    this means that we do not know the name of the attr
    and we cannot assign it to class directly (unless we use metaclasses, but that would not be readable at all)
    
    width, widthChanged = Attr.addBasic("width", lineParamsChanged)
      we can skip width, but then we need to generate some unique name for the instance to store the data
        this storage must not be class atr! and must be initialized before first set
        just name them in sequence and copy the starage of defaults (class) to instance
        where we can override it safely - but default values must not be mutable
        or we need to do deepcopy
    
    bad thing is that split will not know the names of the arguments
    
    staticmethod is also called after class object is done?
    this is bad...
    how do we store default values if so?
    
    propery and signal can be created by module scope function, but we need to store
    default value somewhere
    
    even if we do it module scope function and use
    width, widthChanged = addBasic(default, globalSignal, _attrs)
    _attrs is defined in parent class and is not accessible in base class definition
    anyway it should be specific to the derived class, not shared across all derived classes!
    so we need to create it in the derived class :(
    
    I was able to make the attr to property change via metaclass
    but it fails with segfault if it is mixed with QObject and OGraphicsObject
    """
    
    def __init__(self, **kw):
        super().__init__()
        #when we assign default values, we need to have the private data storage ready
        # because set tests the actual value
        # and we do not want to trigger signals here
        self._attrs = dict(self._attrs) #copy class atribute (defaults) to instance (data storage)
        self.setup(**kw)
        #test
        self.width = 15
    
    #TODO: as in Params, this should emit only one signal of every kind
    # maybe use common signal lock
    def setup(self, **kw):
        other = {}
        for it in kw:
            if hasattr(self, it):
                setattr(self, it, kw[it])
            else:
                other[it] = kw[it]
        return other

'''
class LineParams(QtCore.QObject):
    """
    all line properties
     - width, color, style
    
    takes care of setting pen
    """
    lineParamsChanged = QtCore.pyqtSignal()
    
    widthChanged = QtCore.pyqtSignal(float)
    
    @property
    def width(self):
        return self._width
    
    @width.setter
    def width(self, value):
        if self._width != value:
            self._width = value
            self.widthChanged.emit(value)
            self.lineParamsChanged.emit()
        pass

    colorChanged = QtCore.pyqtSignal(object)
    styleChanged = QtCore.pyqtSignal(object)

    
    def __init__(self, width = 1, color="black", style="-"):
        self.width = width
        self.color = color
        self.style = style
''' 
   
class LineParams(Attrs):
    """
    all line properties
     - width, color, style
    
    takes care of creating pen
    """
    
    lineParamsChanged = QtCore.pyqtSignal()
    
    width = Attr(1, groupSignal=lineParamsChanged)
    color = Attr("black", groupSignal=lineParamsChanged)
    style = Attr("solid", groupSignal=lineParamsChanged)
   
  
    # this could directly handle changes to generate appropriate pen
    
    def _getLinePen(self, widthScale = 1):
        pen = QtGui.QPen()

        pen.setColor(QtGui.QColor(self.color))

        pen.setWidthF(self.width * widthScale)
        #this is not unversal
        # due to scaling of data area
        # lines for SeriesPlotter use this: self._pen.setWidthF(width*GS.dataLineWidth)
        # use appropriate widthScale

        style = self.style
        if style in penStyles.values():
            pen.setStyle(style)
        elif style in penStyles:
            pen.setStyle(penStyles[style]) #this will set dash offset to zero
        else:
            pen.setDashPattern(style) #this will set style to Qt.Qt.CustomDashLine
        return pen
    
class MarkerParams(LineParams):
    """
    all marker properties
    - line properties are used as marker border
    - fill, fillColor (for shapes that have fillable space, like "o")
    - shape
    """
    ...
    
class TextParams:
    """
    color, font, font size, anchor? alignment? rotation?
    """
    ...


"""
ticks rozvaha (muzeme udelat Ticker, ktery to ude resit)
 - uniformly spread over axis (_tickArea, tickNum)
    - bad: not nice round numbers unless axis range is set accordingly
    - bad: does not match nicely with DDllotter edges
    - good: simple
 
 - pseudo uniform in area (tickNum) - adjust to nice numbers   
   - bad: has to adjust position with every zoom/pan
       (but probably can keep already used positions in data coordinates)
    - good: looks nice
       
 - uniform in data coordiantes
    - will be ugly for 
    
 - uniform in data frames (uniform in linear coordinates)
     - might not hit nice numbers
     - should be better for nonuniform data
 
 - specified positions
     - some might not be in visible range
     
     
(muzeme udelat Ticker, ktery to ude resit)
 - graphicsitem, ktery zaroven kresli tick - takze lze setShown
 - pole tick (to bude graphics item cary+label, ktery ma pozici v data coord, linear a area - pricemz jedna z nich bude fixni)
 - num (pocet tick), musi prepocitat pole tick pokud se zmeni
 - musi prepocitat pozice a label tick pokud se osa pohne (pan/zoom)
 
 Ticker.num - number of ticks, some could be out of range, or 
 Ticker.ticks - pole ticks
 Ticker.labels - pole label
 Ticker._anchors ? 
 Ticker.update - recalcualte positions and/or labels
"""

#TODO: Ticker needs to be also TextParams, which will be used for labels
class Ticker(LineParams, QtWidgets.QGraphicsObject):
    #num and length are basically on the same level as LineParams
    # or not, LineParams need redraw, but num needs update
    
    """
    @property
    def num(self):
        return self._num
    
    @num.setter
    def num(self, num):
        if num == self._num : return 
        num = max(0, num)
        self._num = num
        
    """
    
    def _updateNum(self, num) :   
        #recreate internal list of ticks/anchors/labels, ...
        N = len(self._anchors)
        
        #we can reuse old stuff
        if N>num:
            delAnchors = self._anchors[num:]
            self._anchors = self._anchors[:num]
            
            delLines = self._lines[num:]
            self._lines = self._lines[:num]

            delLabels = self._labels[num:]
            self._labels = self._labels[:num]
            
            #TODO: del all cleanly
        else:
            #add some more
            pen = self._getLinePen()
            
            for _ in range(N, num):
                tanchor = Anchor(0, 0, self.parentItem(), color="grey") #TODO: this will only work if Anchor can be repositioned
                self._anchors.append(tanchor)
                
                tend = self.length*self.normal
                tline = QtWidgets.QGraphicsLineItem(QtCore.QLineF(QtCore.QPointF(0,0), tend), tanchor)
                tline.setPen(pen)
                self._lines.append(tline)
            
                tlabel = QtWidgets.QGraphicsSimpleTextItem("", tanchor)
                f = tlabel.font()
                f.setPointSize(GS.tickLabelSize)
                tlabel.setFont(f)
                tlabel.setPos((tend if self.length>0 else QtCore.QPointF(0,0))+GS.tickMinorSpacing*self.normal)
                self._labels.append(tlabel)        
        
        self.update()
    
    num = Attr(5)
    length = Attr(5)
    format = Attr("{:.3g}")
    labelShift = Attr((0,0)) 
    normal = Attr(QtCore.QPointF(0,0))

    #format for labels
    #  but we might also want to use tickFormat from transform
    #  we also need to handle changing transform :(
    #  we should probably ask parent  for transform  every time and 
    #    let parent call updatePositions/updateLabels when transform changes

    
    def __init__(self, axis, length = GS.tickLength, width = GS.plotLinePen.width(), **kw):
        #note that length will override class default with actual GS value
       
        # ~ own, other = Attr.split(kw)
        
        LineParams.__init__(self)
        kw.update(length=length, width=width)
        other = self.setup(**kw)
        QtWidgets.QGraphicsItem.__init__(self, **other)
       
        self._anchors = []
        self._lines = []
        self._labels = []
        
        self.lineParamsChanged.connect(self.redraw)
        self.lengthChanged.connect(self.redrawLength)
        self.formatChanged.connect(self.update)
        
        # ~ self._num = None
        
        self.numChanged.connect(self._updateNum)
        
        self.setParentItem(axis) 
        axis.setTicker(self)
        
        
        self.labelShiftChanged.connect(self.update)
        
        
    
    def redraw(self):
        #visuals changed, redraw
        pen = self._getLinePen()
        for it in self._lines:
            it.setPen(pen)
        pass
        
    def redrawLength(self, *args):
        tend = self.length*self.normal
        
        line = QtCore.QLineF(QtCore.QPointF(0,0), tend)
        for it in self._lines:
            it.setLine(line)
        
        end = line.p2()
            
        for it in self._labels:
            it.setPos(end+GS.tickMinorSpacing*self.normal)
        pass
        
    def _getPositions(self):
        """
        calcualte position of tick anchors in area units (0,1)
        this might be target for overloading in derived classes
        """
        raise NotImplementedError
    
    def update(self, *args): #TODO: not sure if this does not conflict with some Qt name
        if self.parentItem() is None:
             return
             
        #TODO: possibly do prepareGeometryChange() for Qt
        
        positions = self._getPositions()
        
        startPos = self.parentItem().start.anchorPos()
        endPos = self.parentItem().end.anchorPos()
        
        tformat = self.format
        
        data = self.parentItem().area2data(positions)
        # ~ print(positions)
        # ~ print(self._anchors)
        # ~ print(self._labels)
        # ~ print(data)
        #move anchors
        for p, a, l, d in zip(positions, self._anchors, self._labels, data):
            #move anchor, hopefully this will move everything

            position = (1.-p)*startPos+p*endPos
            # ~ print(position)
            a.setPos(*position)
        
            a.setVisible(p>=0 and p<=1)
        
            #TODO: Transform might take precedence in generating tlabels string (e.g. EnumTransform can return strings, we might have time/date, etc.)
            #   we should at least take the default tickFormat from Transform (or override AxisArea default with that)
            l.setText(tformat.format(d))
            #move to anchor
            br = l.boundingRect()
            #top:  -0.5, -1 
            #bottom: -0.5, 0 
            #left:  -1, -0.5
            #right: 0, -0.5
            dx, dy = self.labelShift
            l.setTransform(QtGui.QTransform.fromTranslate(br.width()*dx, br.height()*dy), False)
            
        #TODO: this might be different ?
        self.parentItem()._updatePreferredSize()            

        # ~ print("Ticker.update", [it.pos() for it in self._anchors])
        pass

    def boundingRect(self):
        return self.childrenBoundingRect()
        
    def paint(self, painter, option, widget):
        pass    
    
    pass

class AreaUniformTicker(Ticker):
    def _getPositions(self):
        return np.r_[0.:1.:1j*self.num] #ticks in area units
    pass

class DataConstantTicker(Ticker):
    positions = Attr([])
    
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.positionsChanged.connect(self.update)
    
    def _getPositions(self):
        return self.parentItem().data2area(self.positions)

class AxisArea(PlotArea):
    """
    base class for axis areas
    - is tied to one Axis instance
    - is tied to one PlotGWidget 
    
    keeps
    - data ranges (it needs to keep track of all DataPlotters that use it
    - zoom and pan (mapping of DataPlotters space to (0,1) range of DataArea
    - units, ticks, labels, etc.
    
    #Todo? it would be nice to be able to switch transform, e.g. for timeline plot in DDanalysis switching between no transform and NumEnum
       but in current paradigm it would require resetting all data as DataPlotters do not store data apart from gritems and they use original transform to create gritems
    """
    
    
    
    axisChanged = QtCore.pyqtSignal()
    zoomRangeLChanged = QtCore.pyqtSignal(float, float)
    
    #TODO: this should be the other way round, derive class should override ancestor defaults
    #TODO: we might store default values simply just by keywords=default in __init__
    # but setup will not be able to iterate through that
    #tickFormat - "tickFormat".format(value) will be used to convert data coordinates to string (planned: unless transform provides its own label generation)
     
    #Todo: param specifying if the axis should be zoomed to data on data change, or keep zoom/pan
    def __init__(self, transform=None, position=AxisPosition.bottom, units=None, label=None, tickFormat="{:.1f}", tickNum=7, symmetric=False, *args, **kw):
        self.setTransform(transform)
        #~ kw["backgroundColor"] = "grey"
        super().__init__(*args, **kw)

        self._symmetric = symmetric	       

        if position == AxisPosition.top:
            self.start = self.bottomLeft
            self.end = self.bottomRight
            self.normal = QtCore.QPointF(0,-1)
            self.labelAnchor = Anchor(0.5j,0j, self, color="magenta")
            self._tickLabelShift = -0.5, -1
            self._labelShift = -0.5, 0
        elif position == AxisPosition.bottom:
            self.start = self.topLeft
            self.end = self.topRight
            self.normal = QtCore.QPointF(0,1)
            self.labelAnchor = Anchor(0.5j,1j, self, color="magenta")
            self._tickLabelShift = -0.5, 0
            self._labelShift = -0.5, -1
        elif position == AxisPosition.left:
            self.start = self.bottomRight
            self.end = self.topRight
            self.normal = QtCore.QPointF(-1,0)
            self.labelAnchor = Anchor(0j,0.5j, self, color="magenta")
            self._tickLabelShift = -1, -0.5
            self._labelShift = 0, 0.5
        elif position == AxisPosition.right:
            self.start = self.bottomLeft
            self.end = self.topLeft
            self.normal = QtCore.QPointF(1,0)
            self.labelAnchor = Anchor(1j, 0.5j, self, color="magenta")
            self._tickLabelShift = 0, -0.5
            self._labelShift = -1, 0.5
        else:
            #fall apart
            raise ValueError("AxisArea.__init__: Unknown position", position)

        self.position = position

        if self.position.isHorizontal():
            self._updatePreferredSize = self._updatePreferredSizeH
        else:
            self._updatePreferredSize = self._updatePreferredSizeV

        #zooming
        # - user request to set range shown on Axis to umn, umx (data units space, possibly None to keep previous), but not fix
        # - user request to *fix* shown range to umn, umx (data units space, possibly None to unfix, setter should allow to change just one edge)
        # - dataRange of any data can be registered (assuming that data uses this axis, but not check), lnm, lmx (linear space, basically size of gritems)
        #      - all of registered are taken into account for zoomToData
        # - zoom or pan - adjust existing range by scaling (with one point fixed) and/or moving

        #TODO: _zoomLm[n|x] could be params
        #   or not, because it is more an internal thing, we should keep user specified ranges in data units
        #   but actual zoom state should be stored, because there might be *no* user limits
        # problem is that user should not be concerned about linear space, but if there is no transform, it is actually the user space
        self._zoomLmn = 0.
        self._zoomLmx = 1. #possibly change to zoomLstart, zoomLrange as that is more used (although less intuitive)
        self._dataRange = None, None #min/max range of all present data on this axis (to allow zoom out), None means no data
        self._dataRanges = {}


        self._label = QtWidgets.QGraphicsSimpleTextItem("", self.labelAnchor)
        f = self._label.font()
        f.setPointSize(GS.axisLabelSize)
        self._label.setFont(f)

        self._baseline = AnchorLineItem(self.start, self.end)
        self._baseline.setPen(GS.plotLinePen)

        if self.position.isVertical():
            self._label.setRotation(-90)
            #~ self._label.setPos(0, 0.5)
        if self.position == AxisPosition.T:
            self._label.setY(1)
        # ~ self._label.setFlag(self._label.ItemIgnoresTransformations, True)

        self.ticks = None
        self._defaultTickerSettings = dict(normal=self.normal, labelShift=self._tickLabelShift, num=tickNum, format=tickFormat)
        AreaUniformTicker(self) #Ticker.__init__ will call setTicker (since we need to send self as axis to bind to)


        #params has to be so late, because their onSet will need stuff above
        self.addParam("units", units)
        self.addParam("label", label, onSet=self._updateLabel) #label could be stored directly in self._label (QGraphicsSimpleTextItem) and react to change signal of that, but whatever
        # ~ self.addParam("tickFormat", tickFormat, onSet=self._updateTicks)
        #TODO: add ticker to params (save/load)
       
       
        
        #from UI
        self.addParam("zoomRange", getter=self._zoomRange, setter=lambda r:self._setZoomRange(*r), save=True)
        
        
        
        
        self.setFlag(self.ItemIsMovable)
        #self._updatePreferredSize() #todo: not sure if it is not called from somewhere else
        pass

    def setTicker(self, ticker):
        if self.ticks is not None:
            #we assume here single ownership
            self.zoomRangeLChanged.disconnect(self.ticks.update)
            self.ticks.setParentItem(None)
            self.ticks.scene().removeItem(self.ticks)
            #we need to remove all children as well
            # problem is that tick lines and labels are children of tick anchors
            #  and anchors are children of axis, not ticker itself
            #TODO: move this to Ticker
            for it in self.ticks._anchors:
                it.setParentItem(None)
                it.scene().removeItem(it)
            
            
        ticker.setParentItem(self)
        self.ticks = ticker
        ticker.setup(**self._defaultTickerSettings)
        self.zoomRangeLChanged.connect(self.ticks.update)
        
        pass

    def setTransform(self, transform):
        """
        Sets new transform to the axis. But it will not affect what is already plotted.
        Could be usefull when redoing plot from scratch.
        WIP - be carefull
        """
        self._transform = transform #should be Transform-derived instance, but if it can quack...
        #TODO: something like
        self._tickLabelFactory = self._transform.tickLabelFactory if hasattr(self._transform, "tickLabelFactory") else lambda t: self.tickFormat.format(t)
        # but will it update tickFormat if it is changed?
        # similarly it can be done self.data2linear = self._transform.data2linear if self._transform is not None else lambda x:x
        # so that it need not be tested all the time        

    def setup(self, **kw):
        if "label" in kw: self.label = kw.pop("label")
        if "tickFormat" in kw: self.tickFormat = kw.pop("tickFormat")
        if len(kw)>0:
            print("WARNING: AxisArea.setup cannot (yet) process ", list(kw.keys()), " keywords")
        pass
        
    def _updateLabel(self, text):
        if text is None:
            self._label.setVisible(False)
        else:
            self._label.setText(text)
            self._label.setVisible(True)
            vprint("AxisArea._updateLabelPosition", self.geometry(), self._label.text())
            dx, dy = self._labelShift
            br = self._label.boundingRect()

            if self.position.isHorizontal():
                self._label.setPos(dx*br.width(), dy*br.height()) #relative to self.labelAnchor
                pass
            else:
                self._label.setPos(dx*br.height(), dy*br.width())
                pass

        self._updatePreferredSize() #and position of label normal to axis
        pass
        
    def _preferredHeight(self):
        #lets calculate size - we only care about one dimension
        #horizontal axis - height
        #  - some basic base
        h = 0
        
        #if we have ticks, tickLength
        #btw. if you set tickLength negative, it should drow tick into the data area
        # but it will also mess up tick label position
        #TODO: nevertheless it could be usefull and could be handled
        h += max(self.ticks.length, 0)
        
        #if we have tick labels, highest tick label
        h += max([0]+[it.boundingRect().height() for it in self.ticks._labels])+GS.tickMinorSpacing
            
        if self._label.isVisible():
            #if we have axis label, add its height
            h += self._label.boundingRect().height() + GS.tickMajorSpacing
        
        return h
    
    #_updatePreferredSize and _preferredHeight are separated to allow derived clases to adjust preferred size
    def _updatePreferredSizeH(self):
        h = self._preferredHeight()
        self.setPreferredHeight(h)
        self.setMinimumHeight(h)
        self.setMaximumHeight(h)
    
    def _preferredWidth(self):
        #vertical axis - width
        #  - some basic base
        w = 0

        #if we have ticks, tickLength
        #btw. if you set tickLength negative, it should draw tick into the data area
        # but it will also mess up tick label position
        #TODO: nevertheless it could be usefull and could be handled
        w += max(self.ticks.length, 0)
        
        
        #if we have tick labels, highest tick label
        w += max([0]+[it.boundingRect().width() for it in self.ticks._labels])+GS.tickMinorSpacing
            
        if self._label.isVisible():
            #if we have axis label, add its height
            w += self._label.boundingRect().height() + GS.tickMajorSpacing
        
        return w
    
    def _updatePreferredSizeV(self):
        w = self._preferredWidth()
        self.setPreferredWidth(w)
        self.setMinimumWidth(w)
        self.setMaximumWidth(w)
        #Todo: this could be done by sice policy
        

       
    def wheelEvent(self, ev):
        #~ vprint("AxisArea.wheelEvent", ev.pos(), ev.scenePos(), ev.delta(), self._name, self.mapToItem(self._area, ev.pos()))
        if self.position.isHorizontal():
            l = self.area2linear(ev.pos().x()/self.boundingRect().width())
        else:
            l = self.area2linear(1.-ev.pos().y()/self.boundingRect().height()) #(virtual) area unit space assumes 0,0 corner at bottom left (origin of axes)
        self.zoom(1.+0.03*np.sign(ev.delta()), l)
        pass   
    
    def mouseMoveEvent(self, ev):
        #~ vprint("AxisArea.mouseMoveEvent", "button Left", ev.buttons()==Qt.Qt.LeftButton)
        #super().mouseMoveEvent(ev) #only if you want to move the item
        vprint("\t mouse move", ev.pos().x(), ev.pos().y(), ev.lastPos().x()/self.boundingRect().width())        
        if ev.buttons() == QtCore.Qt.LeftButton:
            if self.position.isHorizontal():
                lpx = ev.lastPos().x()
                w = self.boundingRect().width()
                px = ev.pos().x()
                #TOOD: tohle asi neni dostatecny mapping, bude potreba zohlednit i pocatecni pozici, ne jen velikost
                #muzu zacit tim, ze si najdu v jakych jednotkach je ev.pos()
                #je potreba to namapovat na pozici v dataarea nebo axisarea v anchor space (tj na 0,1)
                d = self.area2linear(lpx/w) - self.area2linear(px/w)
            else:
                lpy = ev.lastPos().y()
                h = self.boundingRect().height()
                py = ev.pos().y()
                d = self.area2linear(1.-lpy/h) - self.area2linear(1.-py/h)
            self.pan(d)
            pass
        pass
    
    def mousePressEvent(self, ev):
        #middle click (or doubleclick) to zoomToData
        if ev.buttons() == QtCore.Qt.MiddleButton:
            self.zoomToData()
        else:
            super().mousePressEvent(ev)
        pass
    
    def data2linear(self, x):
        return self._transform.data2linear(x) if self._transform is not None else x
        
    def linear2data(self, x):
        return self._transform.linear2data(x) if self._transform is not None else x
        
    def area2linear(self, x):
        #~ (x-0)/(X-zoomMn) == (1-0)/(zoomMx-zoomMn)
        #~ return self._zoomRange[0] + x*(self._zoomRange[1]-self._zoomRange[0])
        return self._zoomLmn + x*(self._zoomLmx-self._zoomLmn)
        
    def linear2area(self, x):
        #~ return (x-self._zoomRange[0])/(self._zoomRange[1]-self._zoomRange[0])
        return (x-self._zoomLmn)/(self._zoomLmx-self._zoomLmn)
        
    def data2area(self, x):
        return self.linear2area(self.data2linear(x))
        
    def area2data(self, x):
        return self.linear2data(self.area2linear(x))
    
    def setDataRange(self, data, mn, mx):
        #~ vprint("AxisArea.setDataRange", data, mn, mx)
        #update data range to include mn mx (linear space)
        #problem is that one changed data range cannot easily shrink overall data range (it can easily expand it)
        #TODO: we can do updates to global data range here more effeciently than in _resetDataRange, but it will be more coding and there is no time now
        #TODO: can we afford to exchange limits? what about descending axes? (we might require that linear space is ascending)
        if (mn is not None) and (mx is not None) and mn>mx : mn, mx = mx, mn
        if data in self._dataRanges:
            if mn is None and mx is None and data in self._dataRanges:
                del self._dataRanges[data]
                return self._resetDataRange()
            else:
                omn, omx = self._dataRanges[data]
                gmn, gmx = self._dataRange
                self._dataRanges[data] = mn, mx
                #TODO: we only need to update global data range if mn or mx is outside of current data range or if their old values were the edges of data range (there is no guarantee that another data will not have the same data range, but that cannot be easily determined)
                #~ if (gmn is None and mn is not None) or (
                self._resetDataRange()
                return
        else:
            if mn is not None or mx is not None:
                self._dataRanges[data] = mn, mx
                #possibly we need to update
                self._resetDataRange()
                return
        #if we are set to adjust zoom/pan on change of data range
        ...
        
    def dataRange(self, data=None):
        if data is None:
            return self._dataRange
        else:
            return self._dataRanges[data]
        pass
        
        
    def _resetDataRange(self):
        #go through all data and recreate data range
        #this might be long if lot of data is present
        v = iter(self._dataRanges.values())
        try:
            gmn, gmx = next(v)
            for mn, mx in v:
                if mn<gmn: gmn = mn
                if mx>gmx: gmx = mx
        except StopIteration:
            gmn, gmx = None, None
        self._dataRange = gmn, gmx
        #~ vprint("AxisArea._resetDataRange", self._dataRange, self.linear2data(self._dataRange))
        #if we are set to adjust zoom/pan on change of data range
        ...
        
    
    def pz(self):
        #give pan and zoom coeffs needed to properly move and scale gritems to match selected zoom range
        #TODO: test this
        #~ vprint("AxisArea.pz", "_zoomRange", self._zoomRange, "range", self._zoomRange[1]-self._zoomRange[0])
        return self._zoomLmn, self._zoomLmx-self._zoomLmn
        
    
    def _zoomRange(self):
        return self._zoomLmn, self._zoomLmx
    
    def _setZoomRange(self, lmn, lmx):
        #~ vprint("AxisArea._setZoomRange", self._name, lmn, lmx)
        #lmn, lmx can be None, but will be ignored
        emit = False
        if lmn is not None and lmn!=self._zoomLmn: 
            self._zoomLmn = lmn
            emit = True
        if lmx is not None and lmx!=self._zoomLmx: 
            self._zoomLmx = lmx
            emit = True
        #sanitize this
        if self._zoomLmn>self._zoomLmx: 
            self._zoomLmn, self._zoomLmx = self._zoomLmx, self._zoomLmn
        if self._zoomLmn==self._zoomLmx:
            if self._zoomLmx == 0:
                self._zoomLmn = -0.05
                self._zoomLmx = 0.05
            else:
                self._zoomLmn = 0.95*self._zoomLmx
                self._zoomLmx = 1.05*self._zoomLmx
        if self._symmetric:
            self._zoomLmx = max(abs(self._zoomLmx), abs(self._zoomLmn))
            self._zoomLmn = -self._zoomLmx
        
        if emit: 
            #~ vprint("\t\temit zoomRangeLChanged")
            self.zoomRangeLChanged.emit(self._zoomLmn, self._zoomLmx)
        pass
        
    def setRange(self, umn, umx, fix=False):
        #user request
        self._setZoomRange(*self.data2linear((umn,umx)))
        #TODO: fixed (i.e. if you do not want to show all data)
        ...
    
    def pan(self, pan):
        vprint("AxisArea.pan", pan)
        #~ self._zoomRange = self._zoomRange[0]+pan, self._zoomRange[1]+pan
        #~ self.axisChanged.emit()
        self._setZoomRange(self._zoomLmn+pan, self._zoomLmx+pan)
        pass
        
    
    def zoom(self, relchange, point=None):
        #two kinds of zooming
        # zoom around point - keep specified point fixed and extend/shrink range on both sides
        # - now the point can be either mouse, or fixed to specific level (so that e.g. two vertical axes can have different scales  but common zero level)
        
        if point is None: point = self._zoomLmn
        
        r = relchange * (self._zoomLmx-self._zoomLmn)
        start = point - relchange*(point-self._zoomLmn)
        
        self._setZoomRange(start, start+r)
        pass
    
    def zoomToData(self, data=None, verbose=False):
        #TODO: (option) zoom to visible data only
        #~ vprint("AxisArea.zoomToData", self._name, data, verbose)
        if verbose: print("AxisArea.zoomToData", self._name, self._dataRange, None if data is None else self._dataRanges[data])
        if data is not None and data not in self._dataRanges:
            data = None
            # todo: we can send a plotter to Plot.zoomToData, which will forward it to all axes, but the plotter
            #  might not be tied to all axes - so if self does not know it, simply ignore it (at least for now)
            # it would be nice if it zoomed even unrelated axes only on the visible range, if the plotter

        vprint("AxisArea.zoomToData", self._name, self._dataRange, None if data is None else self._dataRanges[data])
        if data is None:
            self._setZoomRange(*self._dataRange)
        else:
            self._setZoomRange(*self._dataRanges[data])
        pass
    pass
    
class ColorMap(object):
    def __init__(self, *args):
        self._stops = np.empty((4,0), dtype=float)
        # [[s0, s1, s2],
        #  [r0, r1, r2],
        #  [g0, g1, g2],
        #  [b0, b1, b2]]
        if len(args)>0 : self.addStops(*args)
        pass
        
    def addStops(self, stops):
        """
        stops is either {s0:(r0, g0, b0), s1:(r1, g1, b1), ...}
        or [[s0, s1, s2], [r0, r1, r2], [g0, g1, g2], [b0, b1, b2]]
        (you can use [[s0, r0, g0, b0], [s1, r1, b1, g1], ...].T)
        
        where sN is position of stop on (0, 1) range
        and rN, gN, bN are (0, 255) values for RGB color of this position
        
        position in between will be linearly interpolated
        """
        
        if isinstance(stops, dict):
            sstops = np.empty((len(stops), 4), dtype=float)
            for i, s in enumerate(stops):
                sstops[i][0] = s
                try:
                    ttt =  stops[s].lstrip("#")
                    sstops[i][1:] = [int(ttt[i:i+2], 16) for i in (0, 2, 4)]
                except:
                    sstops[i][1:] = stops[s]
            stops = sstops.T
            
        stops = np.concatenate((self._stops, stops), axis=1)
        self._stops = stops[:, stops[0].argsort()]
        pass

    def color(self, data, alpha=None):
        stops = self._stops[0]
        # todo: is it possible for data to be a single value?
        #TODO: this might be a bottleneck
        #we need to add 255 as A for ARGB, even if the image format is RGB32, otherwise we will get rendering errors
        dt = np.dtype((np.int32, {'B':(np.uint8,0), 'G':(np.uint8,1),  'R':(np.uint8,2), 'A':(np.uint8,3)}))
        #NOTE: cannot print this type of array due to bug in numpy
        res = np.empty(data.shape, dtype=dt) #NOTE: cannot print this type of array due to bug in numpy
        #~ res = np.empty_like(data, dtype=dt) #NOTE: cannot use empty_like as that might create (probably) view of array which is not acceptable for QImage
        #TODO: this should handle values over top/bottom of range - easily done with interp left and right optional arguments
        #  - now they are the edge colors
        #TODO: this should handle nans - as it is now, they will be black
        res["A"] = 255 if alpha is None else alpha
        res["R"] = np.interp(data, stops, self._stops[1])
        res["G"] = np.interp(data, stops, self._stops[2])
        res["B"] = np.interp(data, stops, self._stops[3])

        # undefined data points in special color
        # TODO: make the bad value color a param of colormap
        # todo: similar for out of the range values
        bad_indices = np.isnan(data)
        if np.sum(bad_indices)>0:
            res["A"][bad_indices] = 128
            res["R"][bad_indices] = 255
            res["B"][bad_indices] = 0
            res["G"][bad_indices] = 255
        return res

    __call__ = color #just an alias
    pass

class ColorMapString(ColorMap):  #todo: this is a terrible name
    """Returns a #AARRGGBB string instead of int32
    Note that this will make ndarray of object which is possibly a bad idea."""
    def color(self, data, alpha=None):
        # todo: is it possible for data to be a single value?
        stops = self._stops[0]

        #TODO: this should handle values over top/bottom of range - easily done with interp left and right optional arguments
        #  - now they are the edge colors
        #TODO: this should handle nans - as it is now, they will be black

        def stringify(*args):
            return ("#" + ("{:0>2x}" * len(args))).format(*np.asarray(args).astype(np.uint8))


        vstringify = np.vectorize(stringify)


        return vstringify(255 if alpha is None else alpha,
                          np.interp(data, stops, self._stops[1]),
                          np.interp(data, stops, self._stops[2]),
                          np.interp(data, stops, self._stops[3])
                          )

    __call__ = color  # todo: this is also terrible repetition


#TODO: add more colormaps
cmBlueRed = ColorMap({0.:(0, 0, 255), 0.5:(255,255,255), 1.:(255,0,0)})
cmRed = ColorMap({0.:(0, 0, 0), 1.:(255,0,0)})
cmBlue = ColorMap({0.:(0, 0, 0), 1.:(0,0,255)})
cmGreen = ColorMap({0.:(0, 0, 0), 1.:(0,255,0)})
cmWhiteBlack = ColorMap({0.:(255, 255, 255), 1.:(0,0,0)})
cmDeepBlueRed = ColorMap({0.:(0, 0, 50), 0.25:(0, 0, 255), 0.5:(255,255,255), 0.75:(255,0,0), 1.:(50,0,0)})
cmRGBCycle = ColorMap({0.:(255, 0, 0), 1/3:(0,0,255), 2/3:(0,255,0),  1.:(255, 0, 0)})

import colorcet as cc
def _createColorMap(cetCM):
    x = np.r_[0:1.:1j*len(cetCM)]
    return ColorMap(dict(zip(x, np.asarray(cetCM)*255)))

# cmCetBlueRed = _createColorMap(cc.diverging_bwr_20_95_c54)
cmCetBlueGreen = _createColorMap(cc.diverging_bwg_20_95_c41)
cmCetBlueRed = _createColorMap(cc.diverging_bwr_20_95_c54)

colormaps = {"blue-red":cmBlueRed, "red":cmRed, "white-black": cmWhiteBlack, "dark blue-red": cmDeepBlueRed,
             "red-green-blue cycle": cmRGBCycle, "cet green-blue":cmCetBlueGreen, "cet blue-red":cmCetBlueRed}

class ColorAxisArea(AxisArea):
    """
    can also match area space to colors according to colorMap
    """    
    #TODO: params ensuring that coloraxis is zero symmetric (some colormap expect that)
    #TODO: add params ensuring that coloraxis has 0 linear fixed to 0.5 (but that would be basically the same as before since we cannot have two linear scales on one axis)
    #we can have bi-linear transform, why not??
    #   alternatively we can shift colormap "zero" point to 0 linear level to have + and - range different
    #~ paramsDefaults = {"colorMap":cmBlueRed, "colorLow":QtGui.QColor("black").rgb(), "colorHigh":QtGui.QColor("white").rgb()}
    #~ paramsDefaults.update(AxisArea.paramsDefaults)
    #TODO: this is not integrated well with AxisArea preferredSize calculation - but maybe it needs not to be?
    
    def __init__(self, *args, **kw):
        self._colorBarWidth = GS.colorBarWidth
        self._colorBar = QtWidgets.QGraphicsPixmapItem()
        self._levels = 100
        super().__init__(*args, **kw)
        self._colorBar.setParentItem(self)
        self._colorBar.setZValue(-1)
        #TODO: this needs adjusting
        if False:
            if self.position == AxisPosition.right:
                self._area.setX(self._colorBarWidth)
            elif self.position == AxisPosition.top:
                self._area.setY(-self._colorBarWidth)
            elif self.position == AxisPosition.bottom:
                self._area.setY(self._colorBarWidth)
            else: #left
                self._area.setX(-self._colorBarWidth)
        #bude stacit posunout Anchor? TODO: this needs testing
        # ~ vprint("self._colorBarWidth", self._colorBarWidth)
        # ~ vprint("labelAnchor", self.labelAnchor.anchorPos(), self.labelAnchor.pos())
        # ~ for a in [self.start, self.end, self.labelAnchor]+self._tickAnchors:
        for a in [self.start, self.end]:
            if self.position == AxisPosition.right:
                a.setX(a.anchorPos()[0]+self._colorBarWidth)
            elif self.position == AxisPosition.top:
                a.setY(a.anchorPos()[1]+self._colorBarWidth)
            elif self.position == AxisPosition.bottom:
                a.setY(a.anchorPos()[1]-self._colorBarWidth)
            else: #left
                a.setX(a.anchorPos()[0]-self._colorBarWidth)
        self.ticks.update() #this will move ticks anchors with start and end
        # ~ vprint("labelAnchor", self.labelAnchor.anchorPos(), self.labelAnchor.pos())    
      
      
        #~ self._colorBar.setTransformationMode(Qt.Qt.SmoothTransformation)
        
        #~ if self.position.isHorizontal():
            #~ self._colorBar.setTransform(QtGui.QTransform.fromScale(1./self._levels, 0.5), False)
        #~ else:
            #~ self._colorBar.setTransform(QtGui.QTransform.fromScale(0.5, 1./self._levels), False)

        #setting the color map will _updateColorBar
        # self.addParam("colorMap", cmDeepBlueRed, onSet=self._updateColorBar) #this is not yet exposed to UI so do not save
        self.addParam("colorMap", cmCetBlueRed, onSet=self._updateColorBar) #this is not yet exposed to UI so do not save
        #TODO: add colorLow and colorHigh - probably should be part of ColorMap itself, e.g. first and last point?
        self._updatePreferredSize()
        
        pass

    def _preferredHeight(self):
        return super()._preferredHeight()+self._colorBarWidth

    def _preferredWidth(self):
        return super()._preferredWidth()+self._colorBarWidth

    def contextMenuActions(self, ev):
        # specific for FigureScene, returns iterable of (name, callback) pairs which will be included in context menu
        return [("Change ColorMap ...", self._changeColorMap)]

    def _changeColorMap(self):
        """
        Display a dialog for user to select a new colormap
        :return:
        """
        selection, ok = QtWidgets.QInputDialog.getItem(None,  "Select color map", "Color map", list(colormaps.keys()),
                                                   editable=False)
        # print("TEST _changeColorMap selection", selection)
        if ok:
            self.colorMap = colormaps[selection]
            # TODO: this will not update existing plots

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        size = ev.newSize()
        if self.position == AxisPosition.top:
            self._colorBar.setTransform(QtGui.QTransform.fromTranslate(0, self._colorBarWidth).scale(ev.newSize().width()/self._levels, self._colorBarWidth), False)
        elif self.position == AxisPosition.right:
            self._colorBar.setTransform(QtGui.QTransform.fromTranslate(0, size.height()).scale(self._colorBarWidth, -ev.newSize().height()/self._levels), False)
        elif self.position == AxisPosition.bottom:
            self._colorBar.setTransform(QtGui.QTransform.fromScale(ev.newSize().width()/self._levels, self._colorBarWidth), False)
        else: #left
            self._colorBar.setTransform(QtGui.QTransform.fromTranslate(size.width()-self._colorBarWidth, size.height()).scale(self._colorBarWidth, -ev.newSize().height()/self._levels), False)
            pass

    def _updateColorBar(self, *args):
        data = np.array([np.r_[0:1.:1j*self._levels]])
        if self.position.isVertical():
            data = data.T
        self._temp = self.colors(data, False)
        im = QtGui.QImage(self._temp, data.shape[1], data.shape[0], QtGui.QImage.Format_RGB32)
        #~ im.bits() #TODO: for some reason unless I do this, it will output corrupted image for some image sizes ??? - likely caused by destroying of temporal data array created by self.colors after method ends. QImage requires that it survives for whole image life
        self._colorBar.setPixmap(QtGui.QPixmap.fromImage(im))
        pass
    
    def colors(self, data, transform=True, alpha=None):
        """
        data are in linear space of axis 
        zoom/pan to (0, 1)
        assign colors according to colorMap
        TODO: (values outside of range should have specific colors)
        """
        if transform: data = self.linear2area(data)
        return self.colorMap.color(data, alpha=alpha)


#WARNING: when deriving Plotters, esp. implementing Plotter._plot, all "top level" GraphicsItems has to be parented to Plotter itself, otherwise QGraphicsItemGroup.addToGroup will screw up item transform
#from this perspective it would be better to derive Plotter from QGraphicsItem directly, but then we need to implement paint
#also it will be better to have it derived from QGraphicsObject to have signals, which is usefull for e.g. tracer
# or we might have Decorator completely separate from Plotter and have that an object
class DataPlotter(QtWidgets.QGraphicsItemGroup, Params):
    """
    base class for data plotting elements
     - single data elements (e.g. line) on single plot (if you want the same data on another plot, you have to make another instance)
     - data is not stored (apart from generated gritems)
     - visuals can be changed
     - data can be reset (full replot)
     - can be setup without data
     - subclasses might allow altering data (without recreating all gritems)
     - (todo:) can be hidden/show
     - hold appropriate params to define drawing (line color, width, symbols, etc)

    - we need to inform relevant AxisArea of data range if it changes
     
    - this might be derived from QGraphicsItem directly, or possibly from graphics item group, then I would not need to keep track of replotting
            any change in included items should be directly applied to scene and if it is one item, it should be done at once
            
    TODO: we might save reference to original data (not copy), it will not waste too much memory and we could change axis transform
    typical use case would be DD plot with some data and contours from different sources, after than change picture data to some with different axes. If we use NumEnumTransforms, it should redefine Transform according to new axes, but that *has to* change present lines (contours, etc), preferably without replotting everything by hand
    """
    def __init__(self, plot, *args, **kw):
        self._plot = plot
        self._gritems = [] #it is better to store python wrapper object otherwise garbage collector can destroy C++ objects thinking they are no longer referenced (when they are referenced just in C++ insides and not python anymore)
        visible = kw.pop("visible", True)
        super().__init__(*args, **kw) #this will call setup and intialize all plotting parameters, adjusting gritems to look accordingly
        self.setVisible(visible)
        pass

    @property
    def plot(self): #should be read only, setup only through __init__
        return self._plot

    def _removeItem(self, item):
        self.removeFromGroup(item)
        item.setParentItem(None)
        item.scene().removeItem(item)
        self._gritems.remove(item)
        del item
        pass

    def setData(self, *args, **kw):
        if len(self._gritems) > 0:
            #TODO: all this should cause only one redraw (ideally including the _plot afterwards)
            #remove existing gritems from plot scene - destroy all children of self
            for it in self.childItems():
                self._removeItem(it)
            try:
                assert(self._gritems == [])
            except:
                print(self._gritems)
                raise
        
        #create gritems according to data
        self._gritems = self._setData(*args, **kw)
        for it in self._gritems:
            self.addToGroup(it)
        #TODO: we need to test what is the best for redrawing
         #   - directly parent individual gritems to self in _plot
         #   - or create a group (or single item) and parent it to self directly here (single update)
        #TODO: redraw
        pass


# if derived from QGraphicsItem, we need more reimplementations (we will see if things work out with ItemGroup first)
#~ NotImplementedError: QGraphicsItem.boundingRect() is abstract and must be overridden
#~ NotImplementedError: QGraphicsItem.paint() is abstract and must be overridden

    """
    #to subclass, reimplements these:
    def _plot(self, data): 
        create gritems using relevent visual parameters, applying all AxisArea.transforms, i.e. gritems are on linear space
        return [gritem1, gritem2, ...]
        
    def setup(self, **kw): 
        #this is reimplementation of Params.setup
        #it should directly apply relevant changes to gritems (note that there might be no gritems available)
        #we need to handle ours Params separately and not send them to super().setup
        #   we could let setting of the parameters to self attributes to Params.setup, but it would be cleaner to do it here, as it minimizes time when parameters of gritems and attributes are mismatched
            
        emit appropriate signals
        return super().setup(**kw)

    #if is is possible and desirable we can do 
    def adjustData(self, *args, **kw):
        #adjust data in some way (acting on existing gritems)
        #you can of course just create everything from scratch by setData
        pass
        
    def cleanup(self):
        #disconnect any signals, unset data ranges on Axes, etc.
        pass
    """
    
    pass
    
penStyles = {
    "":QtCore.Qt.NoPen, None:QtCore.Qt.NoPen,
    "solid":QtCore.Qt.SolidLine, "-":QtCore.Qt.SolidLine,
    "dash":QtCore.Qt.DashLine,"--":QtCore.Qt.DashLine,
    "dot":QtCore.Qt.DotLine, "..":QtCore.Qt.DotLine,
    "dashDot":QtCore.Qt.DashDotLine, "-.":QtCore.Qt.DashDotLine,
    "dashDotDot":QtCore.Qt.DashDotDotLine, "-..":QtCore.Qt.DashDotDotLine
    }
  
class XYPlotter(DataPlotter): 
    #TODO: clean up order of arguments of __init__ - esp. keywords and *args
    #   do this also in derived classes 
    #   (As it is it is perhaps not the best solution)
    def __init__(self, plot, xaxis, yaxis, *args, name=None, **kw):
        #xaxis and yaxis are AxisArea instances and should not be changed later
        self._xaxis = xaxis
        self._yaxis = yaxis
        super().__init__(plot, *args, **kw)
        self.addParam("name", name)
        if xaxis.position.isHorizontal():
            self._updatePZ = self._updatePZH
        else:
            self._updatePZ = self._updatePZV
        xaxis.zoomRangeLChanged.connect(self._updatePZ)
        yaxis.zoomRangeLChanged.connect(self._updatePZ)
        self._updatePZ() #initialize
        pass        
        
    def _updatePZH(self, *args):
        xp, xz = self._xaxis.pz()
        yp, yz = self._yaxis.pz()
        if xz==0 or yz==0:return 
        self.setTransform(QtGui.QTransform.fromScale(1./xz, 1./yz).translate(-xp, -yp), False)
        pass

    def _updatePZV(self, *args):
        xp, xz = self._xaxis.pz()
        yp, yz = self._yaxis.pz()
        if xz==0 or yz==0:return 
        self.setTransform(QtGui.QTransform.fromScale(1./yz, 1./xz).translate(-yp, -xp), False)
        pass
    
    def cleanup(self):
        self._xaxis.zoomRangeLChanged.disconnect(self._updatePZ)
        self._yaxis.zoomRangeLChanged.disconnect(self._updatePZ)
        self._xaxis.setDataRange(self, None, None)
        self._yaxis.setDataRange(self, None, None)
        pass
    pass

class SmoothPathItem(QtWidgets.QGraphicsPathItem):
    # ~ def __init__(self, *args, **kw):
        # ~ super().__init__(*args, **kw)
        # ~ self._pathItem = QtWidgets.QGraphicsPathItem()
        # ~ self._pathItem.setParentItem(self)
        # ~ p = self._pathItem.pen()
        # ~ p.setCosmetic(True)
        # ~ self._pathItem.setPen(p)
        # ~ pass
    
    def paint(self, painter, *args, **kw):
        #~ vprint("SmoothPathItem.paint", painter, args, kw)
        painter.setRenderHint(painter.Antialiasing) #TODO: this should preserve other flags
        return super().paint(painter, *args, **kw)
    
    def shape(self):
        #QGraphicsLineItem.shape uses its pen to create the shape
        #but it might not preserve cosmetic feature of the pen
        #ie. the shape can be scaled and the width will be changed
        #we might try to set the pen width temporarily here to offset this
        #self.setPen(p)  #this will trigger shape ... :( and infinity recursion
        
        #lets do it by hand and bypass super().shape() completely
        path = self.path()
        # ~ return path #probably too tight
        
        # ~ p = path.translated(0, 0.1) #this works but is not different from what is below
        # ~ p.addPath(path)
        # ~ return p
        
        if (path == QtGui.QPainterPath()):
            return path
            
        ps = QtGui.QPainterPathStroker()
        pen = self.pen()

        ps.setCapStyle(pen.capStyle())

        #we need to manipulate this to reverse scale
        #first lets try something small to see if it can limit itemAt detection
        # ~ ps.setWidth(pen.widthF() <= 0.0 ? penWidthZero : pen.widthF())
        vprint("\t\tSmoothPathItem.shape", self.parentItem().parentItem().scale())
        ps.setWidth(0.01) #OK this will limit with of the shape, however it is not connected to zooming
        #and unfortunately, shape is not called when DataArea is zoomed
        # so we cannot update the shape when zooming and we need to do that to be exact :(
        #we can try to set the width to match initial zoom level and hope that user will not zoom in/out too much
        #however we might even create these path before DataArea is scaled, yes
        #and anyway, scale is bad place to look at, we need transform()
        
        #recalculating all the shapes all the time will be a pain
        # but maybe it is enough to recalculate them when mouseclick is being resolved? i.e. right before scene.itemAt?

        ps.setJoinStyle(pen.joinStyle())
        ps.setMiterLimit(pen.miterLimit())
        p = ps.createStroke(path)
        # ~ self._pathItem.setPath(p)
        # ~ p.addPath(path) #TODO: not sure what this does? it is copy paste from Qt source
        return p

 
        return res

class SeriesPlotter(XYPlotter): 
    """
    - create visuals for sequence of (X, Y) points
    - is tied to specified xaxis and yaxis Axis instaces
    
    visuals can be
     - (todo) scatter points: marker params
     - linear segment line (todo: possibly spline or some other interpolations)
     - (todo) columns
    
    """
    #~ paramsDefaults = {"name":None, "color":"black", "width":1, "style":Qt.Qt.SolidLine, "fill":None}
    #Todo: implement fill param
    #TODO: markers
    #TODO: visible param
    #TODO: when I set pen of line and then change the pen, will the line change?
    #Todo: how to keep legend icon updated when series look change?
    def __init__(self, plot, xaxis, yaxis, color="black", width=1., style=QtCore.Qt.SolidLine, fill=None,*args,  **kw):
        super().__init__(plot, xaxis, yaxis, *args, **kw)
        self._pen = QtGui.QPen() 
        self._pen.setCosmetic(True)
        self.addParamGroup("pen", self._updateLines)
        self.addParam("color", color, getter=lambda : self._pen.color().name(), setter=lambda x:self._pen.setColor(QtGui.QColor(x)), group="pen") #unfortunately we cannot return QColor instance directly as that cannot be recovered when saved to settings
        self.addParam("width", width, getter=self._getPenWidth, setter=self._setPenWidth, group="pen")
        self.addParam("style", style, onSet=self._onSetStyle, group="pen")
        self.addParam("fill", fill)
        #setCapStyle
        #setDashOffset
        #setJoinStyle
        #setMiterLimit
        
        #these are candidates for param group and storing the thing directly in QPen with group getter setter
        pass

    def _setPenWidth(self, width):
        vprint("SeriesPlotter._setPenWidth", width, width*GS.dataLineWidth)
        self._pen.setWidthF(width*GS.dataLineWidth)
    
    def _getPenWidth(self):
        return self._pen.widthF()/GS.dataLineWidth

    def _saveData(self, *args):
        vprint("SeriesPlotter._saveData")
        filepath = QtWidgets.QFileDialog.getSaveFileName(None, "Save XY data", None, "numpy 2x1D array (*.npy)")[0]
        if filepath=="": return
        #we have to reconstruct data from path
        item = self._gritems[0]
        path = item.path()
        
        N = path.elementCount() #we basically assume that all that is presents is lineTo
        px = np.empty(N)
        py = np.empty(N)
        
        for i in range(N):
            e = path.elementAt(i)
            px[i] = e.x
            py[i] = e.y
        
        if self._xaxis.position.isVertical():
            px, py = py, px
        
        x = self._xaxis.linear2data(px)
        y = self._yaxis.linear2data(py)
        
        np.save(filepath, np.array([x,y]))
        
    def contextMenuActions(self, ev):
        #specific for FigureScene, returns iterable of (name, callback) pairs which will be included in context menu
        return [("Save series "+str(self.name), self._saveData)]
    
    
    def _onSetStyle(self, style):
        if style in penStyles.values():
            self._pen.setStyle(style)
        elif style in penStyles:
            self._pen.setStyle(penStyles[style]) #this will set dash offset to zero
        else:
            self._pen.setDashPattern(style) #this will set style to Qt.Qt.CustomDashLine
        self._updateLines()
        pass

    def _updateLines(self):
        if len(self._gritems) > 0:
            self._gritems[0].setPen(self._pen)

    def _setData(self, x, y=None, xIsLinear=False, yIsLinear=False):
        #~ item = QtWidgets.QGraphicsPathItem()
        item = SmoothPathItem()
        item.setParentItem(self)
        vprint("setting pen with width", self._pen.widthF())
        item.setPen(self._pen)
        
        if self.fill is not None:
            brush = QtGui.QBrush(QtGui.QColor(self.fill))
            item.setBrush(brush)
        
        if y is None:
            x,y = y,x
        
        if x is None:
            x = np.arange(len(y))
        
        if len(x)*len(y)>0:
            x = np.asarray(x)
            y = np.asarray(y)
            
            
            x = self._xaxis.data2linear(x) if not xIsLinear else x
            y = self._yaxis.data2linear(y) if not yIsLinear else y
            
            self._xaxis.setDataRange(self, min(x), max(x)) #if you want to delete this Plotter, send setDataRange(None, None, self) to axis
            self._yaxis.setDataRange(self, min(y), max(y))
            
            #TODO: xaxis is not guaranteed to be horizontal, it should be possible to plot lines to (yaxis, xaxis) instead of (xaxis, yaxis) (but not (xaxis, xaxis))
            #   _updatePZ should also reflect this
            # it might be enough to exchange x,y = y,x
            #TODO: this needs testing
            if self._xaxis.position.isVertical():
                x, y = y, x
            
            path = item.path()
            
            path.moveTo(x[0], y[0])
            for i in range(1, min(len(x), len(y))):
                path.lineTo(x[i], y[i])
            
            item.setPath(path)
        else:
            self._xaxis.setDataRange(self, None, None)
            self._yaxis.setDataRange(self, None, None)
        
        return [item]
    
    #this is for compatibility with QCP feature that I use
    #these both assume that you already have the items created and just edit the internal path
    def removeDataBefore(self, x):
        #remove all points that have x<x
        xlim = self._xaxis.data2linear(x)
        item = self._gritems[0]
        path = item.path()
        
        N = path.elementCount() #we basically assume that all that is presents is lineTo
        px = np.empty(N)
        py = np.empty(N)
        
        for i in range(N):
            e = path.elementAt(i)
            px[i] = e.x
            py[i] = e.y
        
        if self._xaxis.position.isVertical():
            px, py = py, px
            
        I = px>xlim
        self.setData(px[I], py[I])
        #this is likely very inneficient, but I do not want to waste time on this
        
    def addData(self, x, y):
        #append data to existing curve
        x = self._xaxis.data2linear(x)
        y = self._yaxis.data2linear(y)        
        item = self._gritems[0]
        path = item.path()

        try:
            xmn, xmx = self._xaxis.dataRange(self)
            self._xaxis.setDataRange(self, min(xmn, min(x)), max(xmx, max(x)))
            ymn, ymx = self._yaxis.dataRange(self)
            self._yaxis.setDataRange(self, min(ymn, min(y)), max(ymx, max(y)))
        except:
            #if we set up and empty path setData([],[]), we will have no dataRange on axes
            self._xaxis.setDataRange(self, min(x), max(x))
            self._yaxis.setDataRange(self, min(y), max(y))

        if self._xaxis.position.isVertical():
            x, y = y, x

        for i in range(0, min(len(x), len(y))):
            path.lineTo(x[i], y[i])
        
        item.setPath(path)
        
    def setup(self, **kw):
        #TODO: possibly allow setting QPen directly? but it will have to be kw separate from Params (which is possible, but Params will not reflect correct values ... (we could copy params from pen itself, but these params will not be to recreate complete QPen (it is out of question to recreate ever option of QPen as Param) - we could unset the params, but that might cause more problems than good
        #TODO: it would be good to trigger redraw only once per setup
        self.startSetup()
        if "color" in kw: self.color=kw["color"]
        if "width" in kw: self.width=kw["width"]
        if "style" in kw: self.style=kw["style"]
        self.stopSetup()
    
        if len(self._gritems) > 0:
            if "fill" in kw:
                self.fill = kw.pop("fill")
                brush = QtGui.QBrush(QtGui.QColor(self.fill))
                self._gritems[0].setBrush(brush)
                #Todo: more parameters for brush
                
            #TODO: possibly unset pen/brush with color=None (is done with style==None), fill=None
            pass
        return super().setup(**kw)
    pass

def calculateContours(A, levels):
    """
    return [[c1, c2, c3, ... at levels[0]], [d1, d2, d3, ... at levels[1], ...]
    """
    assert(np.isfinite(levels).all()) #we cannot calculate imposible contours
    
    CC = [] #levels
    
    for level in levels:
        #vprint("\nlevel", level)
        B = A>level
        imax = A.shape[0]-1
        jmax = A.shape[1]-1

        C = [] #all contours at this level
        
        #we iterate through squares, defined by upper left corner
        # for each edge we need to record if
        # - the edge was not yet evaluated
        # - there is no contour going through it
        # - contour is going through the edge at partial index x

        #note that nan in A indicates unknown value
        # contour cannot be on edge with unknown value (similar to array side)
        # therefore nan should result in False in BL and BU
        # ~ BL = B[:-1]!=B[1:] #left edge has contour
        # ~ BU = B[:,:-1]!=B[:, 1:] #upper edge has contour
        nB = np.logical_not(np.isnan(A))
        BL = np.logical_and(B[:-1]!=B[1:], np.logical_and(nB[:-1], nB[1:])) #vertices are on different side of level and neither is nan
        BU = np.logical_and(B[:,:-1]!=B[:,1:], np.logical_and(nB[:,:-1], nB[:,1:]))
        
        L = np.ones_like(BL, dtype=bool) #left edge not visited yet
        U = np.ones_like(BU, dtype=bool) #upper edge not visited yet
        
        def partial(level, bottom, top):
            # ~ x = (level-bottom) / (top-bottom)
            # ~ however things can get complicated in we have inf or nan
            # with inf as and end point, we cannot use linear interpolation to get partial position
            #both should not be inf, since we should have contour on the edge (and hopefully noone is going to plot inf contour)
            #   however we can have inf and -inf
            # with nan it is again hard because we do not know the value at all - we even do not know if there is contour on the edge - so we should not get nan here
            if np.isinf(bottom):
                if np.isinf(top):
                    #should be +inf, -inf
                    return 0.5
                else:
                    return 1.
            else:
                if np.isinf(top):
                    return 0.
                else:
                    return (level-bottom) / (top-bottom)
        
        def trace(i, j, edge):
            #trace contour starting at cell (i,j) entering from edge "edge"
            # in each cell test if the edge was visited before 
            #  if so we have completed the cycle, mark contour as cyclic and end
            #  if not, calculate partial index, append to contour, mark edge as visited and find exit edge
            contour = QtGui.QPolygonF()
            
            go = True
            
            while go:
                #record entry point 
                #find exit point (possibly record exit point at end)
                if edge == "left":
                    if L[i,j]: #edge not visited yet
                        # ~ x = (level-A[i, j])/(A[i+1,j]-A[i,j]) 
                        x = partial(level, A[i,j], A[i+1,j])
                        contour.append(QtCore.QPointF(i+x,j))
                        L[i,j] = False
                    
                        if j>=jmax: 
                            go = False #we are at edge and cannot go on
                        else:
                            up = BU[i,j]
                            down = BU[i+1,j]
                            if up and down:
                                #we have ambiguous case, decide by middle point (average)
                                up = (A[i:i+2,j:j+2].mean() > level) != B[i,j]
                                down = not up
                            
                            if up:
                                i = i-1
                                edge = "down"
                            elif down:
                                i = i+1
                                edge = "up"
                            elif BL[i,j+1]: #right
                                 j = j+1
                            else:
                                #there is no exit point??? probably can happen if we have nan (empty value) inside
                                go = False
                            pass
                        pass
                    else:
                        #edge already visited, mark contour as cyclic
                        contour.append(contour[0])
                        go = False
                    pass
                elif edge == "right": #after left we cannot get right, so we might as well use elif
                    if L[i,j+1]:
                        # ~ x = (level-A[i, j+1])/(A[i+1,j+1]-A[i,j+1]) 
                        x = partial(level, A[i, j+1], A[i+1,j+1])
                        contour.append(QtCore.QPointF(i+x,j+1))
                        L[i,j+1] = False
                        
                        if j<0: 
                            go = False #end contour if out of bounds
                        else:
                            up = BU[i,j]
                            down = BU[i+1,j]
                            if up and down:
                                #we have ambiguous case, decide by middle point (average)
                                down = (A[i:i+2,j:j+2].mean() > level) != B[i,j]
                                up = not down
                                
                            if up:
                                i = i-1
                                edge = "down"
                            elif down:
                                i = i+1
                                edge = "up"
                            elif BL[i,j]: #left
                                j = j-1
                            else:
                                go = False
                            pass
                        pass
                    else:
                        #edge already visited, mark contour as cyclic
                        contour.append(contour[0])
                        go = False
                        
                
                if edge == "up":
                    if U[i,j]:
                        # ~ x = (level-A[i, j])/(A[i,j+1]-A[i,j]) 
                        x = partial(level, A[i, j], A[i,j+1])
                        contour.append(QtCore.QPointF(i,j+x))
                        U[i,j] = False
                    
                        if i>=imax: 
                            go = False
                        else:
                            right = BL[i,j+1]
                            left = BL[i,j]
                            if right and left:
                                #we have ambiguous case, decide by middle point (average)
                                left = (A[i:i+2,j:j+2].mean() > level) != B[i,j]
                                right = not left
                                
                            if right:
                                j = j+1
                                edge = "left"
                            elif left:
                                j = j-1
                                edge = "right"
                            elif BU[i+1,j]: #down
                                i = i+1
                            else:
                                go = False
                            pass
                        pass
                    else:
                        #edge already visited, mark contour as cyclic
                        contour.append(contour[0])
                        go = False
                    pass
                elif edge == "down":
                    if U[i+1,j]:
                        # ~ x = (level-A[i+1, j])/(A[i+1,j+1]-A[i+1,j]) 
                        x = partial(level, A[i+1,j], A[i+1,j+1])
                        contour.append(QtCore.QPointF(i+1,j+x))
                        U[i+1,j] = False
                        
                        if i<0: 
                            go = False
                        else:
                            right = BL[i,j+1]
                            left = BL[i,j]
                            if right and left:
                                #we have ambiguous case, decide by middle point (average)
                                right = (A[i:i+2,j:j+2].mean() > level) != B[i,j]
                                left = not right
                                
                            if right:
                                j = j+1
                                edge = "left"
                            elif left:
                                j = j-1
                                edge = "right"
                            elif BU[i,j]: #up
                                i = i-1
                            else:
                                go = False
                            pass
                        pass
                    else:
                        contour.append(contour[0])
                        go = False
                #end contour if looped
                #~ if len(contour)>1 and contour[0]==contour[-1]: go=False #this should not be needed
                pass
            return contour
        
        j = 0
        for i in range(imax):
            if BL[i,j] and L[i,j]:
                #left edge has a contour point
                C.append(trace(i, j, "left"))
        j = jmax
        for i in range(imax):
            if BL[i,j] and L[i,j]:
                C.append(trace(i, j-1, "right"))
            pass
        i = 0
        for j in range(jmax):
            if BU[i,j] and U[i,j]:
                C.append(trace(i,j,"up"))
                pass
        i = imax
        for j in range(jmax):
            if BU[i,j] and U[i,j]:
                C.append(trace(i-1, j, "down"))
                pass
            
        for j in range(1,jmax-1):
            for i in range(1,imax-1):
                #test this block left and top edges
                if BL[i,j] and L[i,j]:
                    C.append(trace(i,j,"left"))
                    pass
                elif BU[i,j] and U[i,j]:
                    C.append(trace(i,j,"up"))
                    pass
                pass
            pass
        CC.append(C)        
        pass
        
    #print statistics
    # ~ N = 0
    # ~ NC = 0
    # ~ for C in CC:
        # ~ NC += 1
        # ~ for it in C:
            # ~ N+=it.size()    
    # ~ vprint("\tcalculateContours stats: levels", len(CC), "contours", NC, "points", N)
    
    # ~ vprint("\tcalculateContours")
    # ~ for C in CC:
        # ~ for it in C:
            # ~ vprint("\tcontour", it)
            # ~ for i in range(it.size()):
                # ~ vprint("\t\t", it.at(i))
    return CC
                
        
        
        

class DDPlotter(XYPlotter):
    """
    plot 2D arrays as image
    with possible contours, tracers, colorbar axis, etc
    probably needs to use NumEnum for axes otherwise cannot use same sized pixels (unless it can be guaranteed that axes are equispaced)
    
        Frankly I would let user do as he/she/it wishes
        use data range of supplied 
        then convert data to center of pixel kind
        (using NumEnum with stops at pixel edges should interpolate if that is desired)
        (using Enum will keep uninterpolated pixel position)
        if they mix things up ... they will have to live with it (I might adjust as it goes)
    """
    #~ paramsDefaults = {"colorPlot":True, "contours":False, "contoursWidth":1, "contoursColor":"black", "contoursStyle":"auto", "contoursVisible":True}
    #if colorPlot: try to plot color plot
    #if contours: try to plot contours
    #   controus = scalar - number of levels in dataRange
    #   contours = True, default to 10 levels (5% to 95%)
    #   contours = vector - plot contours on these levels
    # for now, contours will be plotted when setting data and cannot be changed (apart from visuals)
    
    def __init__(self, plot, xaxis, yaxis, colorAxis, colorPlot=True, contours=False, contoursWidth=1., contoursColor="black", contoursStyle="auto", contoursVisible=True, *args, **kw):
        #xaxis and yaxis are AxisArea instances and should not be changed later
        self._zaxis = colorAxis
        self._colorPlot = None
        self._contours = []
        self._alpha = None
        self._ddL = None

        super().__init__(plot, xaxis, yaxis, *args, **kw)
        self._zaxis.zoomRangeLChanged.connect(self._updateImage) 
        
        self._createColorPlot = colorPlot #only valid for the subsequent setData

        self.addParam("contours", contours)
        self.addParamGroup("pen", onSet=self._updateLines) #we might not have one pen for every contour line
        self.addParam("contoursWidth", contoursWidth, onSet=self._updateWidth, group="pen")
        self.addParam("contoursColor", contoursColor, onSet=self._updateColor, group="pen" )
        self.addParam("contoursStyle", contoursStyle, onSet=self._updateStyle, group="pen")
        self.addParam("contoursVisible", contoursVisible, onSet=self._updateVisible, group="pen")
        pass        
    

    def _updateColor(self, color):
        for it in self._contours:
            pen = it.pen()
            pen.setColor(QtGui.QColor(self.contoursColor))
            it.setPen(pen)

    def _updateStyle(self, style):
        if style=="auto": return #cannot change existing contours
        for it in self._contours:
            pen = it.pen()
            style = self.contoursStyle
            if style in penStyles.values():
                pen.setStyle(style)
            elif style in penStyles:
                pen.setStyle(penStyles[style]) #this will set dash offset to zero
            else:
                pen.setDashPattern(style) #this will set style to Qt.Qt.CustomDashLine
            it.setPen(pen)

    def _updateWidth(self, width):
        for it in self._contours:
            pen = it.pen()
            pen.setWidthF(width)
            it.setPen(pen)

    def _updateVisible(self, visible):
        for it in self._contours:
            it.setVisible(self.contoursVisible)
       
    def _updateLines(self):
        self.update()
        pass
       
    def colorAxis(self):
        return self._zaxis
        
    def cleanup(self):
        super().cleanup()
        self._zaxis.zoomRangeLChanged.disconnect(self._updateImage)
        self._zaxis.setDataRange(self, None, None)
        pass

    def setup(self, **kw):
        #TODO: it would be good to trigger redraw only once per setup
        #contours has to be set on per one basis as they might have different pens
        self.startSetup()
        if "contoursColor" in kw: self.contoursColor = kw.pop("contoursColor")
        if "contoursWidth" in kw: self.contoursWidth = kw.pop("contoursWidth")
        if "contoursStyle" in kw: self.contoursStyle = kw.pop("contoursStyle")
        if "contoursVisible" in kw: self.contoursVisible = kw.pop("contoursVisible")
        self.stopSetup()
        return super().setParams(**kw)

    #we need to save data to be able to change color maps, so we might as well allow data saving
    # TODO: ideally just the visible (zoomed-in) part
    # but for now everything, in linear space, without axes
    def _saveData(self, *args):
        vprint("DDPlotter._saveData")
        filepath, fil = QtWidgets.QFileDialog.getSaveFileName(None, "Save DD data", None, "DD + axes, numpy 2D array (*.npy);;DD, numpy 2D array (*.npy)")
        if filepath=="": return
        if fil== "DD + axes, numpy 2D array (*.npy)":
            #try to reconstruct pixel centers and compile larger matrix
            ttt = np.empty((self._ddL.shape[0]+1, self._ddL.shape[1]+1), dtype=self._ddL.dtype)
            ttt[1:,1:] = self._zaxis.linear2data(self._ddL)
            xe, xf = self._xaxis.dataRange(self)
            ye, yf = self._yaxis.dataRange(self)
            #dd.shape == len(y), len(x)
            xd = (xf-xe)/self._ddL.shape[1]
            x = (np.arange(self._ddL.shape[1])+0.5)*xd+xe
            yd = (yf-ye)/self._ddL.shape[0]
            y = (np.arange(self._ddL.shape[0])+0.5)*yd+ye
            ttt[1:,0] = self._yaxis.linear2data(y)
            ttt[0, 1:] = self._xaxis.linear2data(x)
            np.save(filepath, ttt)
        elif fil=="DD, numpy 2D array (*.npy)":
            np.save(filepath, self._ddL)
        pass
        
    def contextMenuActions(self, ev):
        #specific for FigureScene, returns iterable of (name, callback) pairs which will be included in context menu
        return [("Save DD data", self._saveData)]
    
    def dataL(self): #direct access to data (in z-axis linear scale)
        return self._ddL
        
    def _setData(self, dd, x=None, y=None, xAlignment=None, yAlignment=None, dataRange=None, alpha=None):
        """
        dd is 2D array of data such that
            dd.shape == len(y), len(x)
        or 
            dd.shape == len(y)-1, len(x)-1
            
        x and y are data unit coordinates of individual data points. Unlike for lines, both x and y should be monotonous (see below) and you should ensure that they are matched with supplied axes transform to have datapoints equispaced of linear space (e.g. EnumTransform or NumEnumTransform if not sure).

        Basically we do not need whole x and y. We only need edges and information if those are edges of pixels or centers of pixels and we assume that after transform, individual x would be equispaced.
        
        xAlignment can be
          None - pass whole x and edges will be determined automatically
          "center" - values are given for pixel; data range will be expanded 0.5px on both sides
          "leftEdge"  - values are given for left edges; data range will be expanded 1px to right
          "rightEdge" - values are given for right edges; data range will be expanded 1px to left
          "outerEdge" - values are directly data range; this corresponds to situation when len(x) == dd.shape[1]+1
          It is enough to just pass min and max, or rather edge values (which should be min and max after transform), for "center" or "leftEdge" or "rightEdge". You should (not need to) only pass two values for "outerEdge".
        
        self.xaxis is used to transform x (dtto for y)
        self.colorAxis is used to map data amplitude to color levels
        
        TODO: not sure what happens if axis is descending after transform, possibly the color map will not be reversed
           and the axis might be, but I am not sure
           also leftEdges and rightEdges might not work as intended
        """
        
        #~ vprint("DDPlotter._setData", self._zaxis)
        #~ vprint("data alignment", dd.shape, len(x), len(y), len(self._xaxis._transform.stops), len(self._yaxis._transform.stops))
        dd = np.asarray(dd)
        
        self._ddL = self._zaxis.data2linear(dd)
        self._alpha = alpha
        
        #we have to keep nans and infs outside of range determination
        temp = self._ddL.flat
        temp = temp[np.isfinite(temp)] #using numpy.nanmin is not enough since we have to disregard infs too
        #possibly use masked arrays, but that is likely the same underhood (can try timing it to see if it is perhaps faster):  np.ma.masked_invalid([1, 2, np.nan, np.NINF]).min()
        
        try:
            ddLmn = temp.min()
        except:
            ddLmn = 0
        try:
            ddLmx = temp.max()
        except:
            ddLmx = ddLmn + 1
        
        vprint("\tDDLmn, DDLmx", ddLmn, ddLmx)
        
        if dataRange is None:
            self._zaxis.setDataRange(self, ddLmn, ddLmx)
        else:
            self._zaxis.setDataRange(self, *self._zaxis.data2linear(dataRange))
        
        def edges(a, align, axis, N):
            if a is None:
                return axis.data2linear([0, N])

            if align is None:
                M = len(a)
                if M==N: align = "center"
                if M==N+1: align = "outerEdge"

            if align=="outerEdge":
                return axis.data2linear([a[0], a[-1]])
            elif align=="leftEdge":
                e, f = axis.data2linear([a[0], a[-1]])
                return [e, f+(f-e)/(N-1)]
            elif align=="rightEdge":
                e, f = axis.data2linear([a[0], a[-1]])
                return [e-(f-e)/(N-1), f]
            elif align=="center": #alternative to center alignment, delete the one above if you want this
                e, f = axis.data2linear([a[0], a[-1]])
                return [e-(f-e)/(N-1)/2, f+(f-e)/(N-1)/2]
            
            vprint("DDPlotter._setData: cannot determine edge positions.")
            vprint(a)
            vprint(align)
            vprint(axis)
            vprint("a", len(a), "N",  N)
            raise ValueError("DDPlotter._setData: cannot determine edge positions.")

        x = edges(x, xAlignment, self._xaxis, dd.shape[1])
        y = edges(y, yAlignment, self._yaxis, dd.shape[0])
        #~ vprint("DDPlotter._setData", x, y, dd.shape)
        #x and y are now outer edges in linear space where the pixmap should be
        #this is not ncesarily the same size as pixmap itself (unless transform spaces x steps to one unit)
        #  because we only require for transformed x to be equispaced, not going in steps of one
        #  so we might need to scale the pixmap a bit
        sx = abs(x[0]-x[1])/dd.shape[1]
        sy = abs(y[0]-y[1])/dd.shape[0]

        self._xaxis.setDataRange(self, min(x), max(x)) 
        self._yaxis.setDataRange(self, min(y), max(y))
        #if you want to delete this Plotter, send setDataRange(None, None, self) to axis

        self._colorPlot = None #DataPlotter.setData will remove _gritems, but not plotter specific members
        if self._createColorPlot: #TODO: DataPlotter.setData will remove all gritems, but that is not needed here, we just need to call self._updateImage, not create new pixmapItem
            item = QtWidgets.QGraphicsPixmapItem()
            item.setShapeMode(item.BoundingRectShape)
            #~ item.setTransformationMode(Qt.Qt.SmoothTransformation) #TODO: I do not want to do this, but FastTransformation sometimes leads to rendering errors when zooming (ie. scaling the image)
            self._colorPlot = item
            self._colorPlot.setParentItem(self)
            self._updateImage(item)
            if self._xaxis.position.isHorizontal():
                T = QtGui.QTransform.fromTranslate(x[0], y[0]).scale(sx, sy)#TODO: this might not work of x vertical
            else:
                T = QtGui.QTransform.fromTranslate(x[0], y[0]).scale(sy, sx)#TODO: this might not work of x vertical
            #~ vprint("colorPlot scaling", sx, sy, abs(x[0]-x[1]), abs(y[0]-y[1]), dd.shape)
            item.setTransform(T, False)
            pass
            
        self._contours = [] #DataPlotter.setData will remove _gritems, but not plotter specific members
        if (isinstance(self.contours, np.ndarray) and len(self.contours)>0) or self.contours:
            #~ T = QtGui.QTransform.fromTranslate(0.5,0.5).translate(x[0], y[0]).scale(sx, sy)#TODO: this might not work of x vertical
            T = QtGui.QTransform.fromTranslate(x[0], y[0]).scale(sx, sy)#TODO: this might not work of x vertical
            #get contours levels
            dd = (ddLmx-ddLmn)*0.05
            if self.contours is True:
                contours = np.r_[ddLmn+dd:ddLmx-dd:10j] 
            else:
                try:
                    contours = np.r_[ddLmn+dd:ddLmx-dd:1j*self.contours] 
                except (TypeError, ValueError):
                    # ~ contours = self.contours
                    contours = []
                    for it in self.contours:
                        if it.imag!=0:
                            contours.append(ddLmn+(ddLmx-ddLmn)*it.imag)
                        else:
                            contours.append(it)
                pass
            #~ vprint("generate contours", contours, x, y)
            if self._xaxis.position.isVertical():
                C = calculateContours(self._ddL, contours)
            else:
                C = calculateContours(self._ddL.T, contours)
            pen = QtGui.QPen()
            pen.setColor(QtGui.QColor(self.contoursColor))
            pen.setCosmetic(True) #cosmetic pen will have the same width regardless transform
            #setCapStyle
            #setDashOffset
            #setJoinStyle
            #setMiterLimit
            pen.setWidthF(self.contoursWidth)
        
        
            #I will not respect contourStyle here, or at least add another option - to distinguish countour level by dashes
            #here we can possibly also generate different colors for individual levels
            if self.contoursStyle=="auto":
                #we will set contour style on per contour level basis
                #going from solid line for 0 level to dashes at maximum
                # similar to minimum, but with dots inserted
                # 0: solid
                # positive, Np: [mx-i*dn/Np, 2] going from mx to mx-(1-1/Np)*dn (we could go with -i*dn/(Np-1) but that will not work with Np==1)
                # negative, Nn: [mx-i*dn/Nn, 2, 2, 2]
                ttt = np.array(contours)
                ttt.sort()
                positive = ttt[np.where(ttt>0)]
                negative = ttt[np.where(ttt<0)][::-1]
                dashes = {0.:None}
                Np = len(positive)
                Nn = len(negative)
                mx = 50.
                dn = 45.
                for i in range(Np):
                    dashes[positive[i]] = [int(mx-i*dn/Np), 3]
                for i in range(Nn):
                    dashes[negative[i]] = [int(mx-i*dn/Nn), 3, 3, 3]
                pens = []
                for it in contours:
                    pens.append(QtGui.QPen(pen))
                    if dashes[it] is not None:
                        pens[-1].setDashPattern(dashes[it])
                    else:
                        pens[-1].setStyle(QtCore.Qt.SolidLine)
            else:
                if self.contoursStyle in penStyles.values():
                    pen.setStyle(self.contoursStyle)
                elif self.contoursStyle in penStyles:
                    pen.setStyle(penStyles[self.contoursStyle]) #this will set dash offset to zero
                else:
                    pen.setDashPattern(self.contoursStyle) #this will set style to Qt.Qt.CustomDashLine
                pens = [pen]
            
            
            import itertools
            #TODO: contours are not always aligned with pixels - basically contours should go to pixel centers not edges
            #   right now the colorPlot should be 
            for cs, pen in zip(C, itertools.cycle(pens)):
                for c in cs:
                    #~ item = QtWidgets.QGraphicsPathItem()
                    item = SmoothPathItem()
                    item.setParentItem(self)
                    path = item.path()
                    #~ path.moveTo(*c[0])
                    #~ for i in range(1, len(c)):
                        #~ path.lineTo(*c[i])
                    #TODO? close path?
                    p = QtGui.QPolygonF(c)
                    p.translate(0.5,0.5)
                    path.addPolygon(p)
                    item.setPath(path)
                    item.setPen(pen)
                    item.setVisible(self.contoursVisible)
                    item.setTransform(T, False)
                    self._contours.append(item)
        return ([self._colorPlot] if self._colorPlot is not None else []) +self._contours
        
    def _updateImage(self, *args):
        #recreates self._gritems[0].pixmap when
        # - data changes
        # - color map changes
        # - color axis zoom changes
        try:
            if self._colorPlot is None: return 
            #Todo: test - this needs to handle situation where xaxis is not horizontal
            if self._xaxis.position.isVertical():
                self._temp = self._zaxis.colors(self._ddL.T, alpha=self._alpha)
                im = QtGui.QImage(self._temp, self._ddL.shape[0], self._ddL.shape[1], QtGui.QImage.Format_RGB32 if self._alpha is None else QtGui.QImage.Format_ARGB32)
            else:
                #~ vprint("DDPlotter._updateImage", self._zaxis.colors(self._ddL))
                self._temp = self._zaxis.colors(self._ddL, alpha=self._alpha)
                im = QtGui.QImage(self._temp, self._ddL.shape[1], self._ddL.shape[0], QtGui.QImage.Format_RGB32 if self._alpha is None else QtGui.QImage.Format_ARGB32)
            #~ im.bits() #TODO: for some reason unless I do this, it will output corrupted image for some image sizes ???
            # ccan it be caused by destroying of temporal data array returned by _zaxis.colors? QImage might expect its survival
            #~ del self._temp #yes very likely
            px = QtGui.QPixmap.fromImage(im)
            #~ vprint("DDPlotter._updateImage",self._colorPlot, self._colorPlot.isVisible(), self._colorPlot.isEnabled(), self._colorPlot.isObscured(self._colorPlot.boundingRect()))
            #~ px = QtGui.QPixmap(10000,1000)
            #~ px.fill(QtGui.QColor("red"))
            self._colorPlot.setPixmap(px)
            #~ del self._temp # it seems that even QPixmap uses the original data buffer, good
        except:
            import traceback
            print("DDPlotter._updateImage")
            print(traceback.print_exc())
        pass
    ...

"""
Marker
"""
class Marker(QtWidgets.QGraphicsItem):
    def __init__(self, color=None, movable=False, parent=None, **kw):
        self._moveCallback = None
        
        super().__init__(parent)
        #TODO: based on setting select color and shape
        size = 10
        self._rect = QtCore.QRectF(-size/2, -size/2, size, size) 
        self._drawing = QtWidgets.QGraphicsEllipseItem(self._rect,self)
        p = QtGui.QPen(GS.plotLinePen)
#        p.setCosmetic(True)
#        p.setWidthF(1)
        if color is not None: p.setColor(QtGui.QColor(color))
        self._drawing.setPen(p)
        
        
        
        if movable:
            self.setFlag(self.ItemIsMovable, True)
            self.setFlag(self.ItemSendsGeometryChanges, True)
            #I do not want to do this QGraphicsObject but we need a posibility of callback
            if movable is not True:
                self._moveCallback = movable
        pass

    
    def itemChange(self, change, value):
#        vprint("Marker.itemChange", change, value)
        if self._moveCallback is not None and change==self.ItemPositionHasChanged:
            self._moveCallback(value)
        return super().itemChange(change, value)
    
    def pen(self): #so that other can duplicate it
        return self._drawing.pen()
    
    def setColor(self, color):
        p = self._drawing.pen()
        p.setColor(color)
        self._drawing.setPen(p)

    #TODO: we should do proper bounding rect as Qt has bug that does not consider cosmetic pens
    def boundingRect(self):
        #~ res = self.childrenBoundingRect()
        #this will give (-1.6, -1.6, 3.2, 3.2) for (-.1, -.1, .2, .2) ellipse with 3 wide pen 
        #  (-0.6, -.6, 1.2, 1.2) for the same ellipse with 1 wide pen
        #   (-0.15, -0.15, .3, .3) for 0.1 wide pen (which is not visible ...)
        #... it likely does not consider cosmetic pens
        # which is apparent old issue with Qt and not likely to be resolved any time soon :(
        #~ vprint("Decorator.boundingRect", res)        
        return self._rect
    
    def paint(self, *args):
        #at this point uses drawing, but that could change
        pass

#class OriginDecorator(Params, QtWidgets.QGraphicsObject): #note that Params has to be first, because QGraphicsObject apparently also defines its own __getattr__
class OriginDecorator(Params, MovableRelativeAnchorRY): #note that Params has to be first, because QGraphicsObject apparently also defines its own __getattr__
    """
    It inherits "pos" position (setPos, setX, setY) from QGraphicsObject which is in parents coordinates and for us serve as DataArea coordinates on rect [0,0]-[1,1].
    It also has "origin" position in data coordinates (for this to work, x and y axis need to be set in __init__).
    
    graphicdobject position is local coordinates (0,0) mapped to parent coordinates
    decorator should be moveable by mouse, which in default implementation moves graphicsobject position?
    
    decorator should be parented to _area of dataArea (with space limited to 0,0,1,1 rect)
    
    it needs to know to which axes it is connected to to map this to data coordinates and zoom/pan with axes
    (or it might be axis-less and just keep position relative to data area)
    
    
    decorator has position in dataArea (0,0,1,1)
    and (if init with axes) also origin in data coordinates
    """
    
    """
    there is a discrepancy between Plotters and Decorators
     - decorators (in current implementation) place their lines primarily in parent coordinates (i.e. area coordinates minus area coordinates of origin (i.e. pos)) and are not stretched with axes zoom/pan
     - it is about GraphicsItem.pos and GraphicsItem.transform
     - plotters plot in linear space and stretch it to dataarea, whereas decorators take dataarea coordinates (edges) and if need be recalculate data coordinates directly to area (for lines, but not for origin)
     - it is a mess :(    
     
    """
    #TODO: origin is not preserved if DDplotter is changed (i.e. changing interpolation in DDanalysis)
    #   in that case, underlying NumEnumTransform changes, so origin likely keeps its NumEnum linear space position (index) but transform stops change externaly so it corresponds to different data value
    #   what we want to preserve in this situation is data coordinates
    
    #(partially done) add context menu "decorator options" which will trigger setup dialog (could be the same as creational one)
    #   dialog should likely be part of instance ((todo?)each decorator its own dialog? or class atribute?)
    #   (todo?) dialog could be used to store params (and decorator be modified on the run even when dialog is shown)
    #it has inherited xChanged(), yChanged() to notify of position change relative to dataArea if that is what you are interested in (if you zoom/pan it will not move in data space, only on data area)
    #we need also signal for change in data space
    #I will still call it origin to differentiate from pos (which is position in data area coordinates)
    originChanged = QtCore.pyqtSignal(float, float) #notify of origin movement in data coordinates
    originChangedL = QtCore.pyqtSignal(float, float) #notify of origin movement in linear coordinates (after transform)
    decoratorType = "origin"  #for identifying decorator class in saved settings
    dialogUI = None #this will be used for dialog
    def __init__(self, plot, xaxis=None, yaxis=None, color="black", *args, **kw):
        #~ vprint("Decorator.__init__", plot, xaxis, yaxis)
        self._plot = plot
        self._xaxis = None
        self._yaxis = None
        
        #problem now is that we need to divide **kw to those of params and those for other ancestors (QtWidgets.QGraphicsObject)
        # but we do not know params yet, they are defined later in __init__
        # so are we going back to cls.defaultparams dict? but that cannot specify onSet etc...
        # maintaining two lists of params seems tedious
        super().__init__(*args, **kw)
        
        self.setFlag(self.ItemHasNoContents, False) #Anchor sets this, we do not want it here
        self.setFlag(self.ItemClipsToShape, True) #otherwise it will clip to boundingRect which we do not want
        
        self._timer = QtCore.QTimer()
        self._timer.setSingleShot(True)
        self._timer.setInterval(100) #ms
        self._timer.timeout.connect(self._emitOriginChangedTimeout)
        
        #self.pos is in area units
        # but if we are connected to axes, we need to keep constant position in respect to data space (or linear)
        self._xL = None
        self._yL = None
        
        self.setAxes(xaxis, yaxis)
        if plot is not None: plot.addDecorator(self)
        
        

        self.marker = Marker(color=color, parent=self)

        self.addParam("origin", getter=self._getOrigin, setter=self._setOrigin, save=True) #in data coordinates
        self.addParam("color", color, onSet=self._updateColor, save=True)
        
        #the whole Param think is a step in a wrong direction, at least it is not compatible with Qt naming convention (param(), setParam())

        
        #Todo? make position (i.e. widget coordinates) a param too? it will have to be linked to origin
        pass

#    def contains(self, *args):
#        ttt = super().contains(*args)
#        vprint("OriginDecorator.contains", args, ttt, self.isClipped())
#        return ttt
#
    def shape(self, *args):
#        ttt = self.marker.shape()
#        vprint("OriginDecorator.shape", ttt.boundingRect())
#        return ttt
        return self.marker.shape()

#    def collidesWithPath(self, *args):
#        vprint("OriginDecorator.collidesWithPath", self, args)
#        res  = super().collidesWithPath(*args)
#        vprint("\t\t", res)
#        return res

    #testing, likely to cause problems
    #if this works, we can allow dialog to changes axes too
    def setAxes(self, xaxis=Unset, yaxis=Unset):
        if (xaxis is not Unset) and (self._xaxis is not xaxis):
            if self._xaxis is not None:
                self._xaxis.zoomRangeLChanged.disconnect(self._updateXA)
                self.xChanged.disconnect(self._updateXL) 
            
            self._xaxis = xaxis
            if xaxis is not None:
                xaxis.zoomRangeLChanged.connect(self._updateXA)
                self.xChanged.connect(self._updateXL) 
                self._updateXL()

        if (yaxis is not Unset) and (self._yaxis is not yaxis):
            if self._yaxis is not None:
                self._yaxis.zoomRangeLChanged.disconnect(self._updateYA)
                self.yChanged.disconnect(self._updateYL) 
            
            self._yaxis = yaxis
            if yaxis is not None:
                yaxis.zoomRangeLChanged.connect(self._updateYA)
                self.yChanged.connect(self._updateYL) 
                self._updateYL()

    def cleanup(self, *args):
        if self._xaxis is not None:
            self._xaxis.zoomRangeLChanged.disconnect(self._updateXA)
        if self._yaxis is not None:
            self._yaxis.zoomRangeLChanged.disconnect(self._updateYA)
        pass

    def centerOrigin(self, *args):
        #frankly you can just change pos to 0.5, 0.5
        #move origin to 0.5, 0.5 in area units (i.e. center of data area)
        #TODO: this will not work anymore, as self.pos() is not in area units now
        ...

    @property
    def plot(self): #should be read only, setup only through __init__
        return self._plot

    def saveSettings(self):
        #We need to save which axes is self.origin tied to. But we cannot save axes instances. We can only do axes positions in Plot. Also need to know what type of decorator it is (it might be distinguished by key in Plot settings, but I do not want to do a separate category for every type of decorator.
        res = super().saveSettings() #this should save all Params
        #~ res.update({"class":self.__class__, "xaxis":self._xaxis.plotKey, "yaxis":self._yaxis.plotKey})
        res.update({"class":self.decoratorType, "xaxis":self._xaxis.plotKey if self._xaxis is not None else None, "yaxis":self._yaxis.plotKey if self._yaxis is not None else None})
        #saving class as self.__class__ into pickle settings has a disadvantage of easily losing the class definition by e.g. renaming the class, changing module structure, etc. resulting in corrupted save files
        return res

    #for Qt
    def boundingRect(self):
#        res = self.childrenBoundingRect()
        #this will give (-1.6, -1.6, 3.2, 3.2) for (-.1, -.1, .2, .2) ellipse with 3 wide pen 
        #  (-0.6, -.6, 1.2, 1.2) for the same ellipse with 1 wide pen
        #   (-0.15, -0.15, .3, .3) for 0.1 wide pen (which is not visible ...)
        #... it likely does not consider cosmetic pens
        # which is apparent old issue with Qt and not likely to be resolved any time soon :(
#        vprint("Decorator.boundingRect", res, self)
        return self.childrenBoundingRect()

        
    def paint(self, painter, option, widget):
        #only draw children items
        #Todo: origin marker might be composed of children items or drawn here
        #considering the issue of cosmetic pens and children bounding rect, probably draw here
        #(or redesign _area idea, but I will not go there now)
        pass

    #internal
    def _updateColor(self, color):
        color = QtGui.QColor(color)
        self.marker.setColor(color)
        return color
    
    #origin handling
    def _getOrigin(self):
        xD = self._xaxis.linear2data(self._xL) if self._xaxis is not None else None
        yD = self._yaxis.linear2data(self._yL) if self._yaxis is not None else None
        return xD, yD
        
    def _setOrigin(self, origin):
        xD, yD = origin
        #we should change pos directly as that will set _xL, _yL
        #it is an error to try to set origin on non-existing axis
        # but we might allow passing None to existing axis, assuming we only want to change one dimension
        assert(not (xD is not None and self._xaxis is None))
        assert(not (yD is not None and self._yaxis is None))
        xA = self._xaxis.data2area(xD) if xD is not None else None
        yA = self._yaxis.data2area(yD) if yD is not None else None
        if xA is not None and yA is not None:
            self.setPos(xA, yA)
        elif xA is not None:
            self.setX(xA)
        elif yA is not None:
            self.setY(yA)
        #do nothing if both are None
        vprint("OriginDecorator._setOrigin", origin, self.pos())
        vprint("\t\t", self._yaxis.area2data(np.array([0,1])))
        pass
        
    #time this not to emit too often
    def _emitOriginChanged(self):
        #~ vprint("OriginDecorator._emitOriginChanged")
        #~ self._timer.start() #like so it will only emit 100ms after the last request (so after movement ends)
        #it is possible to emit e.g. every 100ms during the movement
        if not self._timer.isActive():
            self._timer.start()
    
    def _emitOriginChangedTimeout(self):
        #~ vprint("OriginDecorator._emitOriginChangedTimeout")
        self.originChanged.emit(*self.origin)
        self.originChangedL.emit(self._xL, self._yL)

    #connecting origin and pos
    def _updateXL(self):
        #position is change in area units, we need to update data coordinates
        #this should not be called unless we have the axis, so output _xL cannot be None
        old = self._xL
#        vprint("OriginDecorator._updateXL", self.x(), self._xaxis.area2linear(self.x()), old)
        self._xL = self._xaxis.area2linear(self.x())
        if old is None or abs(old-self._xL)>1e-9: self._emitOriginChanged()

    def _updateXA(self, *args):
        #axis moved and therefore data coordinates now correspond to different area coordinates
        #this should only be triggered by axis itself, so it exits
#        vprint("OriginDecorator._updateXA", self._xL, self._xaxis.linear2area(self._xL))
        self.setX(self._xaxis.linear2area(self._xL)) #preserve position in data coordinates

    def _updateYL(self):
        #position is change in area units, we need to update data coordinates
        old = self._yL
#        vprint("OriginDecorator._updateYL", self.y(), self._yaxis.area2linear(self.y()), old)        
#        try:
#            vprint("OriginDecorator._updateYL", self.pos(), self.origin)
#        except:
#            pass
        self._yL = self._yaxis.area2linear(self.y())
        if old is None or abs(old-self._yL)>1e-9: self._emitOriginChanged()

    def _updateYA(self, *args):
        #axis moved and therefore data coordinates now correspond to different area coordinates
        #this should only be triggered by axis itself, so it exits
#        vprint("OriginDecorator._updateYA", self._yL, self._yaxis.linear2area(self._yL))
        self.setY(self._yaxis.linear2area(self._yL)) #preserve position in data coordinates

    @classmethod
    def _aggregateDia(cls):
        res = []
        for it in cls.__mro__:
            try:
                res += [(it.decoratorType, it.dialogUI)]
            except AttributeError:
                print("class", it ,"is not part of Decorator chain")
            if it is OriginDecorator:
                break
        return res


    @classmethod
    def dialog(cls, plot, pos, defXAxis=None, defYAxis=None):
        """
        if you need to adjust the dialog you have three options:
        
        1) Set cls.dialogUI = "form.ui", it will be added to default dialog.
        
        2) Overload cls.dialog() completely and return the dialog you would like (it must include cbX, cbY, leX, leY, color, etc. that the default implementation of fromDialog expects).
        
        3) Adjust default dialog by hand
        @classmethod
        def dialog(cls, *args):
            dia = super().dialog(*args)
            dia.layout.addWidget(...)
            ...
            return dia
            
        In all cases adjust setupFromDialog() accordingly.
        
        Decorators created via code will not have self.dia and cannot be changed in UI context menu. If you want that, you have to create self.dia by calling cls.dialog() and filling up the dialog with correct values.
        """
        vprint("cls OriginDecorator.dialog")
        #~ dia = QtWidgets.QDialog()
        #~ loadUi(__file__.replace("pqr.py", "dia.ui"), dia)
        dia = loadUi(__file__.replace("pqr2.py", "dia.ui"))
        #just for testing
        for typ, it in cls._aggregateDia()[::-1]:
            if it is not None:
                gb = loadUi(__file__.replace("pqr2.py", it))
                dia.layout.insertWidget(dia.layout.count()-1, gb)
                setattr(dia, gb.objectName(), gb)
        
        #TODO: tie pos with dia.leX, dia.leY and current selection of cbX, cbY
        
        #add None to axes selections (decorator is not tied to data and can only be moved about in DataArea coordinates)
        dia.cbX.addItem("None", None)
        dia.cbY.addItem("None", None)

        #load axes to cb
        for i, ap in enumerate(plot.axes(AxisPosition.bottom)):
            label = "B"+str(i)+" "+str(ap.label)
            dia.cbX.addItem(label, ap)
        for i, ap in enumerate(plot.axes(AxisPosition.top)):
            label = "T"+str(i)+" "+str(ap.label)
            dia.cbX.addItem(label, ap)
        for i, ap in enumerate(plot.axes(AxisPosition.left)):
            label = "L"+str(i)+" "+str(ap.label)
            dia.cbY.addItem(label, ap)
        for i, ap in enumerate(plot.axes(AxisPosition.right)):
            label = "R"+str(i)+" "+str(ap.label)
            dia.cbY.addItem(label, ap)
        
        #select defXAxis and defYAxis in dia.cbX and dia.cbY if possible
        dia.cbX.setCurrentIndex(dia.cbX.findData(defXAxis))
        dia.cbY.setCurrentIndex(dia.cbY.findData(defYAxis))
        
        dia.color = QtCore.Qt.black
        dia.pushColor.setStyleSheet("QPushButton { background-color: %s}" % "black")
        def setColor(*args):
            dia.color = QtWidgets.QColorDialog.getColor(dia.color, dia)
            dia.pushColor.setStyleSheet("QPushButton { background-color: %s}" % dia.color.name())
        
        dia.pushColor.clicked.connect(setColor)
        dia.adjustSize()
        return dia

    def setupFromDialog(self):
        vprint("OriginDecorator.setupFromDialog")
        dia = self.dia
        self.color = dia.color
        self.marker.setVisible(dia.bShowMarker.isChecked())
        #TODO: dialog also allows changing axes...
        #TODO: dialog also allows changing origin

    @classmethod
    def fromDialog(cls, plot, pos, xaxis=None, yaxis=None, screenPos=None, *args, **kw): # plot, pos, xaxis=None, yaxis=None
        vprint("OriginDecorator.fromDialog", cls, args)

        dia = cls.dialog(plot, pos, xaxis, yaxis)
        #Todo: move dia to mouse click position (or same position as context menu)
        if screenPos is not None:
            dia.move(screenPos)
        
        res = dia.exec_()
        #Todo: add selected color to customColors
        #   check that selected color is not already included
        #   but that is not that easy as we need to check if there is space
        #   shift colors is there is not space
        if res:
            instance = cls(plot, dia.cbX.currentData(), dia.cbY.currentData())
            instance.dia = dia
            dia.cbX.setEnabled(False) #for now Decorators cannot change axes
            dia.cbY.setEnabled(False) #for now Decorators cannot change axes
            # ~ instance.setPos(pos) #TODO: this needs to be replaced by value readout from the dialog (see Todo in cls.dialog() )
            instance.setupFromDialog()
            #now we need to convert dia.leX and dia.leY to origin if possible (they might be in dataarea coordiantes if xaxis and/or yaxis is None)
            # note that dia.leX is QLineEdit because "data coordinates" of axis might be enum or basically whatever
            #  we might need AxisArea to have its own data point entry widget and input that into the dialog ...
            #  or limit ourselves to numbers for now (or just catch EnumTransform and modify default dialog)
            
            return instance
        else:
            return
    
    def contextMenuActions(self, ev):
        vprint("OriginDecorator.contextMenuAction")
        #~ res = [(self.decoratorType.capitalize()+" remove", self._menuRemove)]
        #~ if hasattr(self, "dia"):
            #~ res += [(self.decoratorType.capitalize()+" setup", self._menuSetup)]
        res = [(self.decoratorType.capitalize(), None), (" remove", self._menuRemove)]
        if hasattr(self, "dia"):
            res += [(" setup", self._menuSetup)]
        return res

    def _menuRemove(self, *args):
        self._plot.removeDecorator(self)
        
    def _menuSetup(self, *args):
        dia = self.dia
        res = dia.exec_()
        if res:
            self.setupFromDialog()
        pass

    """
    fromDialog
     - create dia
     - create instance (plot, xaxis, yaxis)
     - setup instance from dia
    
    context setup
     - create dia
     - setup self from dia
     
    but can we make it interactive?
    we cannot for fromDialog as instance is not created before dia is closed
    
    can we make it interactive for after setup adjustment?
    
    
    """
    pass
    
class TextDecorator(OriginDecorator):
    #TODO: we need to be able to set alignment of text around origin (likely override set origin and set position and set text)
    decoratorType = "text"
    #moveable text note with (visible) anchor in data units
    dialogUI = "text.ui"
    def __init__(self, *args, **kw):
        label = QtWidgets.QGraphicsSimpleTextItem("")
        f = label.font()
        f.setPointSize(GS.axisLabelSize)
        label.setFont(f)
        self.label = label
        
        super().__init__(*args, **kw)
        
        label.setParentItem(self)
        # ~ label.setFlag(label.ItemIgnoresTransformations, True)
        #~ self.addToGroup(label)
        
        self.addParam("text", "", onSet=self._updateText, save=True)
        self.marker.setVisible(False) #TODO: something makes the marker appear again
        pass
    
    #~ def _updateOrigin(self, origin):
        #~ #vprint("TextPlotter._updatePos", self, origin)
        #~ x, y = super()._updateOrigin(origin)        
        #~ self.label.setPos(x, y)     
    
    def _updateText(self, text):
        self.label.setText(text)
    
    def _updateColor(self, color):
        color = super()._updateColor(color)
        b = self.label.brush()
        b.setColor(color)
        self.label.setBrush(b)
        pass
    
    #~ @classmethod
    #~ def dialog(cls, *args):
        #~ vprint("TextDecorator.dialog")
        #~ dia = super().dialog(*args)
        #~ t = QtWidgets.QLineEdit()
        #~ dia.formLayout.addRow("Text", t)
        #~ dia.text = t
        #~ dia.adjustSize()
        #~ return dia

    def setupFromDialog(self):
        super().setupFromDialog()
        self.text = self.dia.text.leText.text()
    
#TODO: possibly merge Tracer, Crosshair and Diagonal into one
class TracerDecorator(OriginDecorator):
    """crosshair style tracer showing data unit position of Plot"""
    decoratorType = "tracer"
    dialogUI = "tracer.ui"
    #basically is connected to dataarea, but it will move with original xaxis and yaxis when those are zoom/panned
    def __init__(self, *args, **kw):
        #TODO: it we keep Origin as relative anchor, we need to recat this to anchorline, or directly draw
        #TODO: for this we really should have reasonable shape
        self._hLine = QtWidgets.QGraphicsLineItem(-1, 0., 1, 0.)
        self._vLine = QtWidgets.QGraphicsLineItem(0., -1, 0., 1)
        
        super().__init__(*args, **kw) #will setup some params, therefore all items referenced in _updateParam methods need to be created before super().__init__
        #we have a single point (origin) which lives in data linear space and can be transformed by axes
        # and two lines which operate on _area space, each taking one of origin x,y and the other (0,1)
        # possibly also diagonal (diagonal in data units space)
        
        #we need two lines which will not go out of Plot._da
#        self.setFlag(self.ItemSendsGeometryChanges)
        
        self._hLine.setParentItem(self)
        self._vLine.setParentItem(self)

        pen = self.marker.pen()
        self._hLine.setPen(pen)
        self._vLine.setPen(pen)
        
        self._hText = QtWidgets.QGraphicsSimpleTextItem()
        self._hText.setParentItem(self)
        # ~ self._hText.setFlag(self.ItemIgnoresTransformations)
        
        self._vText = QtWidgets.QGraphicsSimpleTextItem()
        self._vText.setParentItem(self)
        # ~ self._vText.setFlag(self.ItemIgnoresTransformations)
        
        self.yChanged.connect(self._updateVLine)#Todo: this is slightly weird as these are switched (X, y) - (v, h) but otherwise it will not work properly...
        self.xChanged.connect(self._updateHLine)
        #we also need to report the value
        #Todo: add context menu to show/hide hLine, vLine, origin, diagonal, value labels
        #Todo: allow to connect to other axes, with their own value labels
        
        self.setPos(0.5,0.5)
        
    def _updateVLine(self):
        # ~ print("TracerDecorator._updateVLine", self.pos(), self.origin, "line", (0, self.mapFromParent(0,0).y(), 0, self.mapFromParent(1000,1000).y()))
        # ~ print("\t", self.parentItem(), self.parentItem().size())
        h = self.parentItem().size().height()
        o = self.mapFromParent(0,0).y()
        self._vLine.setLine(0, o, 0, o+h)
        #self._vLine.setLine(0, self.mapFromParent(0,0).y(), 0, self.mapFromParent(1000,1000).y())
        ttt = self.mapFromParent(0.,0.1)
        self._vText.setPos(0, ttt.y())
        if self._yL is not None: self._hText.setText("{:.2f}".format(self._yaxis.linear2data(self._yL)))
    
    def _updateHLine(self):
        h = self.parentItem().size().width()
        o = self.mapFromParent(0,0).x()
        self._hLine.setLine(o, 0, o+h, 0)

        # ~ self._hLine.setLine(self.mapFromParent(0,0).x(), 0, self.mapFromParent(1,1).x(), 0)
        ttt = self.mapFromParent(0.01,0)
        self._hText.setPos(ttt.x(), 0)
        if self._xL is not None: self._vText.setText("{:.2f}".format(self._xaxis.linear2data(self._xL)))

    def _updateColor(self, color):
        color = super()._updateColor(color)
        pen = self._hLine.pen()
        pen.setColor(color)
        self._hLine.setPen(pen)
        self._vLine.setPen(pen)
        pass

    def setupFromDialog(self):
        super().setupFromDialog()
        tr = self.dia.tracer
        h = tr.bHorizontal.isChecked()
        v = tr.bVertical.isChecked()
        self._vLine.setVisible(v)
        self._vText.setVisible(v)
        self._hLine.setVisible(h)
        self._hText.setVisible(h)
        pass

class RectDecorator(OriginDecorator):
    """
    actually has two points which can be moved to form a rectangle
    """
    
    rectChanged = QtCore.pyqtSignal(QtCore.QRectF) #notify of (originx, originy, width, height) change in data coordinates
    rectChangedL = QtCore.pyqtSignal(QtCore.QRectF) #notify of origin movement in linear 


    decoratorType = "rect"
    dialogUI = None
    
    def __init__(self, *args, **kw):
        vprint("RectDecorator.__init__", args, kw)
        
#        self._theOther = Marker(movable=self._updateRectSize, color="red")
        self._theOther = OriginDecorator(None, color="red")
#        problem je origin by zachovaval pozici pri zoo,ovani 
#        coz marker nedela
#        a druhou stranu origin nemuze byt child tady ale plot ... 
#        ale v konstruktoru ani self nema nastavena parent, to se stane az pozdeji
#        i kdyz origin decorator potrebuje znat svu plot
        
        
        super().__init__(*args, **kw, startSetup=True) #will setup some params, therefore all items referenced in _updateParam methods need to be created before super().__init__
        #startSetup should delay that
        
        self._rectangle = AnchorRectItem(self, self._theOther)
        
#        self._theOther.originChangedL.connect(self._updateRectSize)
#        self.xChanged.connect(self._updateRectSize)
#        self.yChanged.connect(self._updateRectSize)
#        self._theOther.xChanged.connect(self._updateRectSize)
#        self._theOther.yChanged.connect(self._updateRectSize)
        #TODO: we should connect theOther to move with self (otpinaly with key modifier)
        self.setFlag(self.ItemSendsGeometryChanges, True)
        
#        self._rectangle.setParentItem(self)
        pen = self.marker.pen()
        self._rectangle.setPen(pen)        
        
#        self._theOther.originChangedL.connect(self._updateRect)
        
        
        
        self._timer2 = QtCore.QTimer()
        self._timer2.setSingleShot(True)
        self._timer2.setInterval(100) #ms
        self._timer2.timeout.connect(self._emitRectChangedTimeout)        

        self._theOther.setPos(self.x()+0.1, self.y()+0.1) #bad thing is that we do not have connected origin movements
        
        self.originChangedL.connect(self._emitRectChanged)
        self._theOther.originChangedL.connect(self._emitRectChanged)
        
        self.stopSetup() #the name is terrible , but it should let Params setup all withheld params
        pass

# I cannot even uncomment this, because if I do, it will crash when creating instance of RectDecorator, not output, no questions asked
        # cannot import whole Qt module, but we can live witouth it
    def itemChange(self, change, value):
        vprint("RectDecorator.itemChange", change, value)
        if change==self.ItemPositionChange:
            #value is new pos
            old = self.QtPos() #but we need Qt pos
            diff = value-old
            self._theOther.moveBy(diff.x(), diff.y())
            return value
        return super().itemChange(change, value)
        
    def _updateRectSize(self, *args):
#        vprint("RectDecorator._updateRect", args)
        vprint("RectDecorator._updateRectSize", self.pos(), self._theOther.pos())
        posOther = self._theOther.pos()
        self._rectangle.setRect(0,0, posOther.x()-self.x(), posOther.y()-self.y())
        self._emitRectChanged()

    def _emitRectChanged(self, *args):
        #~ vprint("OriginDecorator._emitOriginChanged")
        #~ self._timer.start() #like so it will only emit 100ms after the last request (so after movement ends)
        #it is possible to emit e.g. every 100ms during the movement
        if not self._timer2.isActive():
            self._timer2.start()
    
    def _emitRectChangedTimeout(self):
        #TODO: this should handle situation when we do not have axes
#        vprint(self.pos, self._theOther.pos)
        
        if self._xL is None or self._yL is None or self._theOther._xL is None or self._theOther._yL is None:
            return
        
        #now what should be the order of values?
        rectL = QtCore.QRectF(self._xL, self._yL, self._theOther._xL-self._xL, self._theOther._yL-self._yL).normalized()
        x, y  = self.origin
        ox, oy = self._theOther.origin
        rectD = QtCore.QRectF(x,y, ox-x, oy-y).normalized()
        
        self.rectChangedL.emit(rectL)
        self.rectChanged.emit(rectD)


    def boundingRect(self):
        #this is bad, beacuse it will not include the area where it was drawn before rectChange was scheduled, it should include _rectangle too
        return self.marker.boundingRect().united(self._theOther.marker.boundingRect()).normalized()

    def setParentItem(self, parent):
        res = super().setParentItem(parent)
        self._theOther.setParentItem(parent)
        return res

    def setAxes(self, xaxis=Unset, yaxis=Unset):
        res = super().setAxes(xaxis, yaxis)
        self._theOther.setAxes(xaxis, yaxis)
        return res

    def _updateColor(self, color):
        color = super()._updateColor(color)
        pen = self._rectangle.pen()
        pen.setColor(color)
        self._rectangle.setPen(pen)
        pass

    

class HorizontalDecorator(OriginDecorator):
    """crosshair style tracer showing data unit position of Plot"""
    decoratorType = "horizontal"
    dialogUI = None #TODO: adjust dialog
    #basically is connected to dataarea, but it will move with original xaxis and yaxis when those are zoom/panned
    def __init__(self, *args, **kw):
        self._hLine = QtWidgets.QGraphicsLineItem(-1, 0., 1, 0.)
        
        super().__init__(*args, **kw) #will setup some params, therefore all items referenced in _updateParam methods need to be created before super().__init__
        #we have a single point (origin) which lives in data linear space and can be transformed by axes
        # and two lines which operate on _area space, each taking one of origin x,y and the other (0,1)
        # possibly also diagonal (diagonal in data units space)
        
        #we need two lines which will not go out of Plot._da
        self.setFlag(self.ItemSendsGeometryChanges)
        
        self._hLine.setParentItem(self)

        pen = self.marker.pen()
        self._hLine.setPen(pen)
       
        self.xChanged.connect(self._updateHLine)
        #we also need to report the value
        self.setPos(0.5,0.5)
        
    def _updateHLine(self):
        self._hLine.setLine(self.mapFromParent(0,0).x(), 0, self.mapFromParent(1,1).x(), 0)

    def _updateColor(self, color):
        color = super()._updateColor(color)
        pen = self._hLine.pen()
        pen.setColor(color)
        self._hLine.setPen(pen)
        pass

    def setupFromDialog(self):
        super().setupFromDialog()
        tr = self.dia.tracer
        h = tr.bHorizontal.isChecked()
        self._hLine.setVisible(h)
        pass


class CrosshairDecorator(OriginDecorator):
    """
    similar to tracer, but has horizontal and vertical lines, separated by "spacing"
    """
    #TODO: will not update lines when zooming
    decoratorType = "crosshair"
    dialogUI = "crosshair.ui"
    def __init__(self, *args, **kw):
        #~ self._hLines = [QtWidgets.QGraphicsLineItem(-1, 0., 1, 0.) for it in range(-1, 2)]
        #~ self._vLines = [QtWidgets.QGraphicsLineItem(0., -1, 0., 1) for it in range(-1, 2)]

        self._hLines = []
        self._vLines = []
        color = kw.pop("color", "black")
        super().__init__(*args, color=Unset, **kw) #will setup some params, therefore all items referenced in _updateParam methods need to be created before super().__init__
        
        self.setFlag(self.ItemSendsGeometryChanges)
        
        self._hText = QtWidgets.QGraphicsSimpleTextItem()
        self._hText.setParentItem(self)
        # ~ self._hText.setFlag(self.ItemIgnoresTransformations)
        
        self._vText = QtWidgets.QGraphicsSimpleTextItem()
        self._vText.setParentItem(self)
        # ~ self._vText.setFlag(self.ItemIgnoresTransformations)
        
        self.yChanged.connect(self._updateLinesH)
        self.yChanged.connect(self._updateLinesV)
        self.xChanged.connect(self._updateLinesV)
        self.xChanged.connect(self._updateLinesH)
        #we also need to report the value
        #Todo: add context menu to show/hide hLine, vLine, origin, diagonal, value labels
        #Todo: allow to connect to other axes, with their own value labels

        #~ this will create _hLines and _vLines so it has to be before __init__ which sets the color for them
        #~ but it also has to be after __init__ otherwise params cannot be set
        self.addParam("linesH", 1, onSet=self._createLinesH, save=True)
        self.addParam("linesV", 1, onSet=self._createLinesV, save=True)
        self.addParam("separationH", 0.1, onSet=self._updateLinesH, save=True)
        self.addParam("separationV", 0.1, onSet=self._updateLinesV, save=True)
        self.color = color #using pqr.Unset for super().__init__ allows postponing setting color only after lines are ready
        #TODO: another problem is that linesH will create lines, but not position them
        # initially right after creation we cannot position the lines because separationH is not defined
        # but after it is and we change number of lines, they should be position as well and not wait for other impulses
        #possibly we can create groups in another way: it will be defined after all parameters are set and will trigger if any is changed (obviously cannot be triggered when some parameter is not set)
        # or we can define conditional onSet, e.g. we can set the color of lines on change of color, but only if some other parameters are set
        #but it is crazy to do all that just to facilitate __init__ (once all params are set, there are not problems)
        # only how to set color on ancestor class and have it processed in descendant only after it is ready?
        # we can leave params values unset in __init__ and call (overloaded) setup at __init__ end?
        # i.e. in Desc.__init__(**kw_params_setup) call super().__init__() without params kw, that will just create params but not assign value
        # (i.e. Anc.__init__ will call setup, but without any values and it will not do anything)
        # Decs.__init__ can call setup after it is ready using all params kw (even those for Anc) and its setup will then pass them to Anc.setup 
        #I can solve this with present Params by setting value of color in Anc.__init__ from kw and sending pqr.Unset to super().__init__ here, only sending the intended value after lines are setup
        #TODO: optimize order of calling in response to params change

        self.setPos(0.5,0.5)
        pass
        
    def _createLinesH(self, *args):
        for it in self._hLines:
            it.setParentItem(None)
            if it.scene() is not None: it.scene().removeItem(it)
        self._hLines = [QtWidgets.QGraphicsLineItem(-1, 0., 1, 0.) for it in range(1+2*self.linesH)]
        pen = self.marker.pen()
        for it in self._hLines:
            it.setParentItem(self)
            it.setPen(pen)
        pass
    
    def _createLinesV(self, *args):
        for it in self._vLines:
            it.setParentItem(None)
            if it.scene() is not None: it.scene().removeItem(it)
        self._vLines = [QtWidgets.QGraphicsLineItem(-1, 0., 1, 0.) for it in range(1+2*self.linesV)]
        pen = self.marker.pen()
        for it in self._vLines:
            it.setParentItem(self)
            it.setPen(pen)
        pass
        
    def _updateLinesH(self, *args):
        tl = self.mapFromParent(0,0)
        top = tl.y()
        left = tl.x()
        br = self.mapFromParent(1,1)
        bottom = br.y()
        right = br.x()
        sep = self.separationH
        
        if self._yL is not None:
            yD = self._yaxis.linear2data(self._yL)
            for i in range(-self.linesH, self.linesH+1):
                ttop = self._yaxis.data2area(yD+i*sep)-self.pos().y()
                self._hLines[i+self.linesH].setLine(left, ttop, right, ttop)
        else:
            for i in range(-self.linesH, self.linesH+1):
                self._hLines[i+self.linesH].setLine(left, i*sep, right, i*sep)
        
        ttt = self.mapFromParent(0.01,0)
        self._hText.setPos(ttt.x(), 0)
        if self._xL is not None: self._vText.setText("{:.2f}".format(self._xaxis.linear2data(self._xL)))

    def _updateLinesV(self, *args):
        tl = self.mapFromParent(0,0)
        top = tl.y()
        left = tl.x()
        br = self.mapFromParent(1,1)
        bottom = br.y()
        right = br.x()
        sep = self.separationV
        
        if self._xL is not None:
            xD = self._xaxis.linear2data(self._xL)
            for i in range(-self.linesV, self.linesV+1):
                tt = self._xaxis.data2area(xD+i*sep)-self.pos().x()
                self._vLines[i+self.linesV].setLine(tt, top, tt, bottom)
        else:
            for i in range(-self.linesV, self.linesV+1):
                self._vLines[i+self.linesV].setLine(i*sep, top, i*sep, bottom)
        
        ttt = self.mapFromParent(0.,0.1)
        self._vText.setPos(0, ttt.y())
        if self._yL is not None: self._hText.setText("{:.2f}".format(self._yaxis.linear2data(self._yL)))
        
    def _updateColor(self, color):
        color = super()._updateColor(color)
        pen = self.marker.pen()
        pen.setColor(color)
        for it in self._hLines: it.setPen(pen)
        for it in self._vLines: it.setPen(pen)
        pass
    
    @classmethod
    def dialog(cls, *args):
        #~ vprint("CrosshairDecorator.dialog")
        dia = super().dialog(*args)
        #we need to implement the logic
        ch = dia.crosshair
        def match(value):
            if ch.bSame.isChecked():
                ch.dH.valueChanged.disconnect(match)
                ch.dV.valueChanged.disconnect(match)
                ch.dV.setValue(value)
                ch.dH.setValue(value)
                ch.dH.valueChanged.connect(match)
                ch.dV.valueChanged.connect(match)
                
        ch.dH.valueChanged.connect(match)
        ch.dV.valueChanged.connect(match)
        return dia

    def setupFromDialog(self):
        ch = self.dia.crosshair
        self.linesH = ch.iH.value()
        self.linesV = ch.iV.value()
        self.separationH = ch.dH.value()
        self.separationV = ch.dV.value()
        h = ch.bH.isChecked()
        v = ch.bV.isChecked()
        for it in self._hLines: it.setVisible(h)
        for it in self._vLines: it.setVisible(v)
        super().setupFromDialog() #we can only update color after the lines are created during setting self.linesH and self.linesV
        pass

class DiagonalDecorator(OriginDecorator):
    """
    similar to tracer, but has diagonal and antidiagonal lines, separated by "spacing"
    """    
    #TODO: will not update lines when zooming
    decoratorType = "diagonal"
    dialogUI = "diagonal.ui"
    def __init__(self, *args, **kw):
        self._aLines = []
        self._dLines = []
        color = kw.pop("color", "black")
        super().__init__(*args, color=Unset, **kw) #will setup some params, therefore all items referenced in _updateParam methods need to be created before super().__init__
        
        self.setFlag(self.ItemSendsGeometryChanges)
       
        self.yChanged.connect(self._updateLines)
        self.xChanged.connect(self._updateLines)
        #we also need to report the value
        #Todo: add context menu to show/hide hLine, vLine, origin, diagonal, value labels
        #Todo: allow to connect to other axes, with their own value labels

        #~ this will create _hLines and _vLines so it has to be before __init__ which sets the color for them
        #~ but it also has to be after __init__ otherwise params cannot be set
        self.startSetup()
        self.addParam("linesA", 1, onSet=self._createLinesA, save=True)
        self.addParam("linesD", 1, onSet=self._createLinesD, save=True)
        self.addParam("separation", 0.1, onSet=self._updateLines, save=True)
        self.addParam("separationType", "XY", onSet=self._updateLines, save=True)
        self.addParam("points", 10, onSet=self._updateLines, save=True)

        self.color = color #using pqr.Unset for super().__init__ allows postponing setting color only after lines are ready
        self.stopSetup()
        
        #TODO: another problem is that linesH will create lines, but not position them
        # initially right after creation we cannot position the lines because separationH is not defined
        # but after it is and we change number of lines, they should be position as well and not wait for other impulses
        #possibly we can create groups in another way: it will be defined after all parameters are set and will trigger if any is changed (obviously cannot be triggered when some parameter is not set)
        # or we can define conditional onSet, e.g. we can set the color of lines on change of color, but only if some other parameters are set
        #but it is crazy to do all that just to facilitate __init__ (once all params are set, there are not problems)
        # only how to set color on ancestor class and have it processed in descendant only after it is ready?
        # we can leave params values unset in __init__ and call (overloaded) setup at __init__ end?
        # i.e. in Desc.__init__(**kw_params_setup) call super().__init__() without params kw, that will just create params but not assign value
        # (i.e. Anc.__init__ will call setup, but without any values and it will not do anything)
        # Decs.__init__ can call setup after it is ready using all params kw (even those for Anc) and its setup will then pass them to Anc.setup 
        #I can solve this with present Params by setting value of color in Anc.__init__ from kw and sending pqr.Unset to super().__init__ here, only sending the intended value after lines are setup
        #TODO: optimize order of calling in response to params change

        self.setPos(0.5,0.5)
        pass
        
    def _createLinesA(self, *args):
        vprint("Diagonal._createLinesA")
        for it in self._aLines:
            it.setParentItem(None)
            if it.scene() is not None: it.scene().removeItem(it)
        self._aLines = [QtWidgets.QGraphicsPathItem() for it in range(1+2*self.linesA)]
        pen = self.marker.pen()
        for it in self._aLines:
            it.setParentItem(self)
            it.setPen(pen)
        pass
    
    def _createLinesD(self, *args):
        for it in self._dLines:
            it.setParentItem(None)
            if it.scene() is not None: it.scene().removeItem(it)
        self._dLines = [QtWidgets.QGraphicsPathItem() for it in range(1+2*self.linesD)]
        pen = self.marker.pen()
        for it in self._dLines:
            it.setParentItem(self)
            it.setPen(pen)
        pass

    def _updateLines(self, *args):
        #~ vprint("DiagonalDecorator._updateCross")
        sep = self.separation
        sepT = self.separationType

        """
        we want to plot 
        x[i] = self._xL[data units] + i*vectX
        y[i] = self._yL[data units] + i*vectY
        
        with (most sense) vectX = vectY (to have diagonal)
        however the vector might be in principle whatever
        
        we want to distribute self.points number of points along the path (e.g. squeezing/stretching the vector if need be)
        
        if we do not have axis defined, then we should operate on DataArea units (what is the vector then is a mystery)
        
        ok so for now vectX,vectY = 1,1
        """
        x0 = self.pos().x()
        y0 = self.pos().y()
        
        if self._xaxis is not None:
            xmn, x, xmx = self._xaxis.area2data(np.array([0, x0, 1.]))
        else:
            xmn, x, xmx = 0., x0, 1.

        if self._yaxis is not None:
            ymn, y, ymx = self._yaxis.area2data(np.array([0, y0, 1.]))
        else:
            ymn, y, ymx = 0., y0, 1.
        
        for m in range(-self.linesD, self.linesD+1):
            #diagonal
            imn = max(xmn-x-m*sep, ymn-y)
            imx = min(xmx-x-m*sep, ymx-y)
            i = np.r_[imn:imx:1j*self.points]
            X = x+m*sep + i
            Y = y + i

            vprint("diagonal", m)
            vprint("X", X)
            vprint("dX", np.diff(X))
            vprint("Y", Y)
            vprint("dY", np.diff(Y))

            if self._xaxis is not None:
                X = self._xaxis.data2area(X)
            if self._yaxis is not None:
                Y = self._yaxis.data2area(Y)
            
            vprint()
            vprint("X", X)
            vprint("dX", np.diff(X))
            vprint("Y", Y)
            vprint("dY", np.diff(Y))
            vprint()
            vprint()
            
            
            path = QtGui.QPainterPath()
            path.moveTo(X[0]-x0, Y[0]-y0)
            for k in range(1, self.points):
                path.lineTo(X[k]-x0, Y[k]-y0)
            
            self._dLines[m+1].setPath(path)

        for m in range(-self.linesA, self.linesA+1):
            #antidiagonal
            imn = max(xmn-x-m*sep, y-ymn)
            imx = min(xmx-x-m*sep, y-ymx)
            i = np.r_[imn:imx:1j*self.points]
            X = x+m*sep + i
            Y = y - i
            
            if self._xaxis is not None:
                X = self._xaxis.data2area(X)
            if self._yaxis is not None:
                Y = self._yaxis.data2area(Y)
            
            path = QtGui.QPainterPath()
            path.moveTo(X[0]-x0, Y[0]-y0)
            for k in range(1, self.points):
                path.lineTo(X[k]-x0, Y[k]-y0)
            
            self._aLines[m+1].setPath(path)
        pass
        
    def _updateColor(self, color):
        color = super()._updateColor(color)
        pen = self.marker.pen()
        pen.setColor(color)
        for it in self._aLines: it.setPen(pen)
        for it in self._dLines: it.setPen(pen)
        pass
    
    def setupFromDialog(self):
        dd = self.dia.diagonal
        self.startSetup()
        self.linesA = dd.iA.value()
        self.linesD = dd.iD.value()
        self.separation = dd.dSep.value()
        self.separationType = dd.cbSep.currentText()
        super().setupFromDialog() #we can only update color after the lines are created during setting self.linesH and self.linesV
        self.stopSetup()
        d = dd.bD.isChecked()
        a = dd.bA.isChecked()
        for it in self._dLines: it.setVisible(d)
        for it in self._aLines: it.setVisible(a)
        pass

DecoratorClasses = {it.decoratorType:it for it in [OriginDecorator, TextDecorator, TracerDecorator, CrosshairDecorator, DiagonalDecorator, RectDecorator]}

#TODO: we need to allow to anchor legend to corners of plot
#  note that at creation we do not know size of DataArea
class Legend(QtWidgets.QGraphicsRectItem):
    _margin = 0.01
    _spacing = 0.005
    def __init__(self, plot, *args, **kw):
        self._plot = plot
        super().__init__(*args, **kw)
        self._entries = OrderedDict()
        p = self.pen()
        p.setCosmetic(True)
        p.setStyle(QtCore.Qt.NoPen)
        self.setPen(p) 
        self.setRect(0, 0, 0, 2*self._margin)
        self.setFlag(self.ItemIsMovable)
        pass

    def addSeries(self, seriesPlotter):
        #create gritems for it
        #also resize self
        #for now just name with seriesPlotter.color
        label = QtWidgets.QGraphicsSimpleTextItem(str(seriesPlotter.name))
        label.setParentItem(self)
        
        f = label.font()
        f.setPointSize(GS.axisLabelSize)
        label.setFont(f)
                
        # ~ label.setFlag(label.ItemIgnoresTransformations, True)
        b = label.brush()
        try:
            b.setColor(QtGui.QColor(seriesPlotter.color))
        except AttributeError: #TODO: clean up the interface after moving name form SeriesPlotter to XYPlotter
            b.setColor(QtGui.QColor(seriesPlotter.contoursColor))
        label.setBrush(b)
        self._entries[seriesPlotter] = label
        
        rl = label.mapRectToParent(label.boundingRect())
        r = self.rect()
        
        label.setPos(self._margin, r.bottom()-self._margin+self._spacing)
        
        r.setBottom(r.bottom()+rl.height()+self._spacing)
        if r.width()-2*self._margin < rl.width():
            r.setWidth(rl.width()+2*self._margin)
        
        
        self.setRect(r)
        vprint("Legend.addSeries: new rect", r)
        pass
    
    def removeSeries(self, seriesPlotter):
        if seriesPlotter not in self._entries:
            print("WARNING: Legend.removeSeries - series", seriesPlotter, "is not in Legend", self)
            return
        else:
            keys = self._entries.keys()
            width = 0
            for it in keys:
                if it is seriesPlotter:
                    label = self._entries.pop(seriesPlotter)
                    rl = label.mapRectToParent(label.boundingRect())
                    label.setParentItem(None)
                    break
                else:
                    label = self._entries[it]
                    rl = label.mapRectToParent(label.boundingRect())
                    width = max(width, rl.width())
            
            #now we have to move all labels below - entries is ordered dict so it should be possible
            for it in keys:
                self._entries[it].moveBy(0, -rl.height())
                label = self._entries[it]
                rl = label.mapRectToParent(label.boundingRect())
                width = max(width, rl.width())
                
            
            #and update master rect 
            r = self.rect()
            r.setBottom(r.bottom()-rl.height()-self._spacing)
            r.setWidth(width+2*self._margin)
            self.setRect(r)
        pass

#TODO: similar to Legend, we need the ability to anchor this to specific DataArea point
# 
class Label(QtWidgets.QGraphicsSimpleTextItem):
    def __init__(self, plot, *args, **kw):
        self._plot = plot
        super().__init__(*args, **kw)
        self.setFlag(self.ItemIsMovable)

        f = self.font()
        f.setPointSize(GS.axisLabelSize)
        self.setFont(f)
        pass

      
class LabelGraphicsLayoutItem(QtWidgets.QGraphicsLayoutItem):
    def __init__(self, text, *args, **kw):
        self._label = QtWidgets.QGraphicsSimpleTextItem()
        if "color" in kw:
            self._label.setBrush(QtGui.QBrush(QtGui.QColor(kw.pop("color"))))
        if "brush" in kw:
            self._label.setBrush(kw.pop("brush"))
        if "font" in kw:
            self._label.setFont(kw.pop("font"))
        if "fontSize" in kw:
            f = self._label.font()
            f.setPointSizeF(kw.pop("fontSize"))
            self._label.setFont(f)
        
        self._center = kw.pop("center", True)

        super().__init__(*args, **kw)

        #TODO: _background is more or less not needed (unless we want to do color background under the title, which I am not a fan of)
        #TODO: we might want to do alignment possible
        self._background = QtWidgets.QGraphicsRectItem(0,0,100,self._label.boundingRect().height())
        pen = self._background.pen()
        pen.setStyle(QtCore.Qt.NoPen)
        self._background.setPen(pen)
        self._label.setParentItem(self._background)
        self.setGraphicsItem(self._background)
        #~ self.setGraphicsItem(self._label)
        self.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Label)
        self.setText(text)
        pass 
    
    def setText(self, text=None):
        if text is None:
            self.setMaximumSize(0,0)
            self._label.setText("")
        else:
            self.setMaximumSize(-1, -1)
            self._label.setText(text)
        self.updateGeometry()
    
    def setGeometry(self, rect):
        #~ vprint("LabelGraphicsLayoutItem.setGeometry", rect, "text", self._label.text())
        #~ self.prepareGeometryChange()
        super().setGeometry(rect)
        self._background.setRect(0, 0, rect.width(), rect.height())
        self._background.setPos(rect.topLeft())
        #~ self._label.prepareGeometryChange()
        if self._center: self._label.setX(rect.width()/2-self._label.boundingRect().width()/2)
        #~ vprint("\t", self._label.pos(), self._label.boundingRect(), rect.width()/2-self._label.boundingRect().width()/2)
        pass
        
    #virtual QSizeF sizeHint(Qt::SizeHint which, const QSizeF &constraint = QSizeF()) const = 0
    def sizeHint(self, which, constraint):
        #~ vprint("LabelGraphicsLayoutItem.sizeHint", "text", self._label.text(), which, constraint, self._label.boundingRect().size(), self._label.boundingRect())
        #~ vprint("    maximum", self._label.parentItem(), QtCore.QSizeF(self._label.parentItem().boundingRect().width(), self._label.boundingRect().height()))
        #~ vprint("   parent", self._label.parentItem(), self._label.parentItem().boundingRect())
        #~ vprint("   parent", self._background.parentItem(), self._background.parentItem().boundingRect())
        if which==QtCore.Qt.MaximumSize:
            #~ return QtCore.QSizeF(self._label.parentItem().boundingRect().width(), self._label.boundingRect().height())
            return QtCore.QSizeF(self._background.parentItem().boundingRect().width(), self._label.boundingRect().height())
        else: # which==Qt.Qt.MinimumSize or which==Qt.Qt.PreferredSize or whatever
            return self._label.boundingRect().size()
        #Todo: what are constraints for?

    
class Plot(QtWidgets.QGraphicsWidget):
    #has placed axes
    #all data shown on the plot need to select one of its axes to be plotted against
    
    #has data view area, with [0,0,1,1] dimensions (stretched to fill up the available space)
    
    #possibly other elements, like title, legend etc..
    """
    has save/load interface similar to SaveLoadElement or Params
    keep in mind that save settings should only store UI defined values which are not known at compile time and can only be used on exactly the same Plot layout, with loaded data and added Plotters
    """
    def __init__(self, *args, **kw):
        super().__init__(*args)
        self._axes = {} #keep track of all AxisArea added to the Plot, store its globalZoom selection
        #Todo: this type of storage is not well readable, but if there is not other parameter we need
        #  we might do self._axes = [] and store globalZoom on axis itself

        #~ self._items = []
        self._plotters = []
        self._decorators = []
        self._legend = Legend(self)
        
        self._layout = QtWidgets.QGraphicsGridLayout()
        da = DataArea(name="data") #, backgroundColor="grey") #QtWidgets.QGraphicsWidget()
        da.setFlag(da.ItemClipsChildrenToShape, True)
        # ~ self.setAutoFillBackground(True)
        
        #~ self._legend.setParentItem(da._area)
        self._legend.setParentItem(da)
        
        #~ self._ttt = QtWidgets.QGraphicsRectItem(0.25,0.25,0.5,0.5)
        #~ self._ttt.setParentItem(da._area)
        
        self._da = da #and it should be scaled to whatever dimensions Plot can provide
        
        self._llabel = LabelGraphicsLayoutItem(kw.pop("title", None), fontSize=GS.plotTitleSize, center=False)  #TODO: will not be centered when not shown but directly saveImage 
        self._layout.addItem(self._llabel, 0,0, 1,3)
        
        self._lhover = Label(self)
        self._lhover.setParentItem(da)
        self._lhover.setPos(0, 0)
        self._lhover.setText("hover")
        
        #layouts for axes
        self._lt = QtWidgets.QGraphicsLinearLayout()
        self._lb = QtWidgets.QGraphicsLinearLayout()
        self._lr = QtWidgets.QGraphicsLinearLayout()
        self._ll = QtWidgets.QGraphicsLinearLayout()
        
        self._lt.setOrientation(QtCore.Qt.Vertical)
        self._lb.setOrientation(QtCore.Qt.Vertical)
        
        #TODO: add title self._layout.addItem(
        N = 1
        self._layout.addItem(self._da, N+1,1)
        self._layout.addItem(self._lt, N, 1)
        self._layout.addItem(self._lb, N+2,1)
        self._layout.addItem(self._lr, N+1, 2)
        self._layout.addItem(self._ll, N+1,0)
        
        # ~ self._layout.setSpacing(0)
        self._layout.setColumnSpacing(0, 0)
        self._layout.setColumnSpacing(1, 0)
        # ~ self._layout.setColumnSpacing(2, 0)
        self._layout.setRowSpacing(N, 0)
        self._layout.setRowSpacing(N+1, 0)
        # ~ self._layout.setRowSpacing(N+2, 0)
        self._layout.setContentsMargins(0,0,0,0)

        #layout stretch works in a strange way (?); documentation is sketchy, to make it work ~ as intended, stretch factor needs to be large (it seems it considers default stretch as 1, so it will try to stretch all columns
        #possibly we can override that using sizepolicy on individual widgets
        self._layout.setColumnStretchFactor(1, 1000)
        self._layout.setRowStretchFactor(N+1, 1000)
        #this stopped working so lets try to scale by changing size policy
        self._da.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        
        
        self.setLayout(self._layout)
        self.setup(**kw)
        pass

    def setup(self, **kw):
        """
        setup plot
        
        setup axes
         - name axes to create/modify by position around plot (t)op, (b)ottom, (l)efl, (r)ight or (x)horizontal, (y)vertical + number (1-9)
         - the behavior is slightly different based on whether the axis is already created or not
         - for example: 
              - "t1", "t2" on empty plot will create two axes on top of plot (the number is order of processing) - these will be accessible as Plot.t1 (or Plot.x1) and Plot.t2 (or Plot.x2)
              - "x1", "x2" alternates between "t" and "b" to keep number of axes the same
              - "t2" on empty plot will create single axis on top of plot (will be accessible as Plot.t1 or Plot.x1
              - "b1", "x1" will create only single axis from b1 params as Plot.b1 (or Plot.x1) that will be then modified by x1 params
              
              - if b1 and t1 are already created, "x1" refers to the one created first
        """
        if "title" in kw:
            self._llabel.setText(kw.pop("title"))
        
        if "legend" in kw:
            self._legend.setVisible(kw.pop("legend"))
        #have to be carefull, because axis added as b1 will be acesssible as x1 too, or possibly x2 if there was t1 added first
        # and e.g. color axis might cause some surprises
        # e.g. using b1label="b1" together with x1label="bbb" can have unexpected results, because it will not set up two axes (b1 will be created first and it will became x1, therefore x1 options will be then applied to just created b1)

        #process setup
        #axes
        #transform, where, axisClass=AxisArea, globalZoom=True, *args, **kw
        #two options: 
        #  "aN = {"transform":..., "label":...}"
        #  aNtransform = ..., aNlabel = ..., ...
        #where "a" gives position (where): "x" = horizontal, "y" = vertical, "b" = bottom, ...
        # and "N" gives order of calling addAxis, only 1-9
        # transform defaults to None
        #lets keep it simple :) although possibly slow
        for a in "btrlxy":
            for i in range(1,10):
                aN = a+str(i)
                has = False
                for it in kw:
                    if it.startswith(aN): 
                        has = True
                        break
                if has:
                    d = kw.pop(aN, {})
                    for it in list(kw.keys()):
                        if it.startswith(aN):
                            d[it[2:]] = kw.pop(it)
                    
                    #TODO: if user sends just e.g. b2 and b1 is not created, b1 will be created, not b2, because addAxis enumerates them automatically in the order of creation
                    #  shall we allow the use to name the axis as wanted?
                    # basically the number should represent order of processing of params
                    #   but if the axes are already created it represents order of axis creation 
                    if hasattr(self, aN):
                        getattr(self, aN).setup(**d)
                    else:
                        self.addAxis(d.pop("transform", None), AxisPosition[a], **d)
        pass
        

    """
    save / load of Plot
    Plot saves settings - i.e. things setup by UI, not data
    you need to set the plot up, plot all data and only after you can apply settings to get the UI to the same state
    plotters - data keeping elements, has to be create before loadSettings, setting consist of visuals (line colors, etc.)
    axes - position and connection to data has to be created before loadSettings, settings consist of pan/zoom (possilby label and visuals like fonts and colors)
    decorators - created through UI, but can be tied to axes (has to be loaded after axes are created), these can be created from settings alone
    """
    def resizeEvent(self, *args):
        # ~ vprint("\tPlot.resizeEvent", args, self.geometry(), self.parentItem().geometry())
        # ~ pol = self._da.sizePolicy()
        # ~ vprint("\t\tda", self._da.geometry(), self._da.preferredSize(), self._da.minimumSize(), self._da.maximumSize(), pol.horizontalPolicy(), pol.verticalPolicy(), pol.horizontalStretch(), pol.verticalStretch())
        # ~ try:
            # ~ vprint("\t\tx1", self.x1.geometry())
            # ~ pol = self.x1.sizePolicy()
            # ~ vprint("\t\tda", self.x1.geometry(), self.x1.preferredSize(), self.x1.minimumSize(), self.x1.maximumSize(), pol.horizontalPolicy(), pol.verticalPolicy(), pol.horizontalStretch(), pol.verticalStretch())
        # ~ except:
            # ~ pass
        # ~ try:
            # ~ vprint("\t\ty1", self.y1.geometry())
        # ~ except:
            # ~ pass
        # ~ try:
            # ~ vprint("\t\tx2", self.x2.geometry())
        # ~ except:
            # ~ pass
        # ~ try:
            # ~ vprint("\t\ty2", self.y2.geometry())
        # ~ except:
            # ~ pass
        self._llabel.updateGeometry() #because _llabel derives its maximum size from Plot geometry, we need to let it know it has changed
        return super().resizeEvent(*args)

    #~ def updateGeometry(self, *args):
        #~ vprint("\tPlot.updateGeometry", args, self.geometry())
        
        #~ return super().updateGeometry(*args)

    def setTitle(self, text):
        self._llabel.setText(text)

    def saveSettings(self):
        #~ vprint("Plot.saveSettings")
        res = {}
        plotters = []
        #plotters - has to be kept in order, must match already loaded data
        #~ for i in range(len(self._plotters)):
            #~ ttt = self._plotters[i].saveSettings()
            #~ vprint(i, self._plotters[i], ttt)
            #~ if len(ttt)>0 : res["plotter"+str(i)] = ttt
        for it in self._plotters:
            plotters.append(it.saveSettings())
        res["plotters"] = plotters
            
        #now we have problems because some plotters are completely UI generated (like notes)
        
        #DataArea should not have anything to save
        
        #axes
        axes = {ap.plotKey:ap.saveSettings() for ap in self._axes}
        #~ for i in range(self._lb.count()):
            #~ res["B"+str(i)] = self._lb.itemAt(i).saveSettings()
        #~ for i in range(self._lt.count()):
            #~ res["T"+str(i)] = self._lt.itemAt(i).saveSettings()
        #~ for i in range(self._ll.count()):
            #~ res["L"+str(i)] = self._ll.itemAt(i).saveSettings()
        #~ for i in range(self._lr.count()):
            #~ res["R"+str(i)] = self._lr.itemAt(i).saveSettings()
        res["axes"] = axes
            
        #decorators
        dec = []
        for it in self._decorators:
            dec.append(it.saveSettings())
        res["decorators"] = dec
        
        return res

    def loadSettings(self, settings):
        #plotters
        #~ for i in range(len(self._plotters)):
            #~ key = "plotter"+str(i)
            #~ if key in settings: self._plotters[i].loadSettings(settings[key])
        for pls, pl in zip(settings["plotters"], self._plotters):
            if len(pls)>0: pl.loadSettings(pls)
            
        #DataArea should not have anything to save
        
        #axes
        axes = settings["axes"]
        #~ vprint(axes)
        for i in range(self._lb.count()):
            key = "b"+str(i+1)
            #~ if key in axes: self._lb.itemAt(i).loadSettings(axes[key])
            self._lb.itemAt(i).loadSettings(axes[key])
        for i in range(self._lt.count()):
            key = "t"+str(i+1)
            #~ if key in axes: self._lt.itemAt(i).loadSettings(axes[key])
            self._lt.itemAt(i).loadSettings(axes[key])
        for i in range(self._lr.count()):
            key = "r"+str(i+1)
            #~ if key in axes: self._lr.itemAt(i).loadSettings(axes[key])
            self._lr.itemAt(i).loadSettings(axes[key])
        for i in range(self._ll.count()):
            key = "l"+str(i+1)
            #~ if key in axes: self._ll.itemAt(i).loadSettings(axes[key])
            self._ll.itemAt(i).loadSettings(axes[key])
            
        #decorators
        self.clearDecorators()
        for it in settings["decorators"]:
            xaxis = it.pop("xaxis")
            yaxis = it.pop("yaxis")
            if xaxis is not None: xaxis = getattr(self, xaxis)
            if yaxis is not None: yaxis = getattr(self, yaxis)
            
            try:
                dec = self.addDecorator(it["class"](self, xaxis, yaxis)) #TODO: this will only work on decorators derived from OriginDecorator
                it.pop("class")
            except:
                #~ dec = self.addDecorator(DecoratorClasses[it.pop("class")], xaxis, yaxis) #TODO: this will only work on decorators derived from OriginDecorator
                dec = self.addDecorator(DecoratorClasses[it.pop("class")](self, xaxis, yaxis)) #TODO: this will only work on decorators derived from OriginDecorator
            dec.loadSettings(it)
        pass

    def zoomToData(self, *args, exclude=[], **kw):
        #Todo: raise flag to postpone replots
        for axis in self._axes:
            #~ if self._axes[axis]: #globalzoom
            if axis not in exclude: axis.zoomToData(*args, **kw)
            
        #Todo: replot once
        pass

    def _zoomDataArea(self, daPos, delta):
        #forwarded wheel event from dataArea
        #zoom all axes
        #TODO: this should trigger only one update
        for ap in self._axes:
            if self._axes[ap]:
                #TODO: this duplicates a bit AxisArea.wheelEvent, same part of code shoudl be moved to a method?
                if ap.position.isHorizontal():
                    l = ap.area2linear(daPos.x())
                else:
                    l = ap.area2linear(daPos.y())
                ap.zoom(1.+0.03*np.sign(delta), l)
        pass

    def _panDataArea(self, lastPos, pos):
        for ap in self._axes:
            if self._axes[ap]:
                #TODO: this duplicates a bit AxisArea.mouseMoveEvent, same part of code shoudl be moved to a method?
                if ap.position.isHorizontal():
                    d = ap.area2linear(lastPos.x()) - ap.area2linear(pos.x())
                    #d *= 0.5 #TODO: this is hack I have no idea why DataArea pans twice in horizontal direction
                    #tady bude problem v okamziku, kdy mam dve spojene osy (jako treba vlnove delky a frekvence), a posunu tady obe dve
                    #to je konceptualni problem, musim spojene osy tady vyradit, ale nemam je jak poznat
                else:
                    d = ap.area2linear(lastPos.y()) - ap.area2linear(pos.y())
                # ~ vprint("Plot._panDataArea", ap.position.isHorizontal(), d)
                ap.pan(d)
        pass
     
    def _hover(self, pos):
        try:
            """make current position of mouse in data coordinates available"""
            x = pos.x()
            y = pos.y()
            
            #if _lhover is not anchored, we need to calculate relative position in DataArea
            s = self._da.size()
            x /= s.width()
            y = 1- y/s.height()
            
            t = ""
            
            for ax in self._axes:
                if isinstance(ax, ColorAxisArea): continue #this would be good too, but this needs access to data at give coordinates for DDplots (and this will be more complex)
                #TODO: an alternative approach would be to have an anchor/tracer follow the hover
                #TODO: switch hover on/off for relevant hoverEnter, hoveLeave
                #TODO: the whole hover coudl be part of DataArea instead of plot
                t += f"{ax.label if ax.label else ax.plotKey}:"+ \
                        ax.ticks.format.format(ax.area2data(x if ax.position.isHorizontal() else y))+", "
                       
                       
            self._lhover.setText(t[:-2])
        except:
            print("WARNING: _hover failed. This is not critical")
            import traceback
            traceback.print_exc()
            pass
    
    def axes(self, position):
        if position == AxisPosition.bottom:
            return [self._lb.itemAt(i) for i in range(self._lb.count())]

        if position == AxisPosition.top:
            return [self._lt.itemAt(i) for i in range(self._lt.count())]

        if position == AxisPosition.right:
            return [self._lr.itemAt(i) for i in range(self._lr.count())]

        if position == AxisPosition.left:
            return [self._ll.itemAt(i) for i in range(self._ll.count())]
        return
        

    def addAxis(self, transform, where, axisClass=AxisArea, globalZoom=True, *args, **kw):
        #this will create appropriate widget for it and add it to specified layout
        # data need to know about axis not the plotter widget
        # we need to create some shortcuts for easy access
        #on the other hand, axis plotter widget need to be accessible too, for setup ...
        #however we can have some hidden Axis instances included (e.g. as parents) to keep things tied together
        #probably return axis plotter object - that can be used for setup
        #note that it is not possible to remove axis
        
        if where == AxisPosition.horizontal:
            #look for which of top, bottom has more space and add there
            where = AxisPosition.bottom if (self._lb.count()<=self._lt.count()) else AxisPosition.top
        elif where == AxisPosition.vertical:
            #look for which of left, right has more space and add there
            where = AxisPosition.left if (self._ll.count()<=self._lr.count()) else AxisPosition.right
            
        ap = axisClass(transform, where, *args, **kw)
        #TODO: this adds left axes to right side (i.e. to plot not to edge), I guess the same goes for top
        if where == AxisPosition.top:
            self._lt.addItem(ap)
            setattr(self, "t"+str(self._lt.count()), ap)
            setattr(self, "x"+str(self._lt.count()+self._lb.count()), ap)
            ap.plotKey = "t"+str(self._lt.count())
        elif where == AxisPosition.bottom:
            self._lb.addItem(ap)
            setattr(self, "b"+str(self._lb.count()), ap)
            setattr(self, "x"+str(self._lt.count()+self._lb.count()), ap)
            ap.plotKey = "b"+str(self._lb.count())
        elif where == AxisPosition.left:
            self._ll.addItem(ap)
            setattr(self, "l"+str(self._ll.count()), ap)
            setattr(self, "y"+str(self._lr.count()+self._ll.count()), ap)
            ap.plotKey = "l"+str(self._ll.count())
        elif where == AxisPosition.right:
            self._lr.addItem(ap)
            setattr(self, "r"+str(self._lr.count()), ap)
            setattr(self, "y"+str(self._lr.count()+self._ll.count()), ap)
            ap.plotKey = "r"+str(self._lr.count())
        else:
            print("WARNING: Plot.addAxis - uknown position", where)
            return
            
        if "name" not in kw: ap._name = ap.plotKey
        self._axes[ap] = globalZoom
        #~ self._axes[axis] = ap
        #~ self._axesPositions[axis] = where
        return ap
        
    def plot(self, *args, **kw):
        #plot something
        # probably return relevant plotter object
        #this should be equivalent to 
        #  pl = self.addPlotter( ... )
        #  pl.setData(...)
        #or
        #  self.addPlotter(...).setData(...)  
        #possibly it might recognize what the plotter type should be
        pass
        
    def addDecoratorFromDialog(self, decoratorClass, pos, screenPos=None):
        #default axes
        dec = decoratorClass.fromDialog(self, pos, self.x1 if hasattr(self, "x1") else None, self.y1 if hasattr(self, "y1") else None, screenPos=screenPos)
#        if dec is not None: self.addDecorator(dec)

    def addDecorator(self, decorator):
        #this should be called by Decorator.__init__
        if decorator not in self._decorators: self._decorators.append(decorator)
#        decorator.setParentItem(self._da._area)
        decorator.setParentItem(self._da) #if decorators are derived from RelativeAnchor, they can be parented directly to _da
        return decorator
    #~ def addDecorator(self, decoratorClass, *args, **kw):
        #~ decorator = decoratorClass(self, *args, **kw)
        #~ self._decorators.append(decorator)
        #~ decorator.setParentItem(self._da._area)
        #~ return decorator
        
    def addDDPlotter(self, xaxis=None, yaxis=None, zaxis=None, **kw):
        if zaxis not in self._axes:
            zaxis = self.addAxis(zaxis, AxisPosition.right, axisClass=ColorAxisArea, globalZoom=False, **kw.pop("colorAxisSetup", {})) #TODO: this will add axis under self.yN series which is probably not good
        plotter = self.addXYPlotter(DDPlotter, xaxis, yaxis, zaxis, **kw)
        if plotter.name is not None: self._legend.addSeries(plotter)
        return plotter
        
    def addSeriesPlotter(self, xaxis=None, yaxis=None, **kw): #this is mostly equivalent to QCP.addGraph
        plotter = self.addXYPlotter(SeriesPlotter, xaxis, yaxis, **kw)
        if plotter.name is not None: self._legend.addSeries(plotter)
        return plotter
        
    def addXYPlotter(self, plotterClass, xaxis=None, yaxis=None, *args, **kw): #this is mostly equivalent to QCP.addGraph
        """
        xaxis (and yaxis) can be
         - None - a non-transformed AxisArea will be created, added to Plot and used
         - AxisArea instance which is already a part of the plot (for cenvenience theses can be accessed as Plot().x1 etc) - it will be used
            - using AxisArea which is not part of the Plot is undefined (likely to crash)
         - Transform instance - a transformed AxisArea will be created, added to the plot and used
            - Transform instances can be shared by multiple AxisArea in multiple Plots
        Axes will be added to top or bottom (xaxis) or 
        Using vertical axes for x data is not currently supported, but will be in the future. 
         
        """
        xok = xaxis in self._axes
        yok = yaxis in self._axes
        
        if xok!= yok:
            #one is added, one not - the not added one should go orthogonal to one added
            if xok:
                if xaxis.position.isHorizontal():
                    #y should be added to vertical position
                    yaxis = self.addAxis(yaxis, AxisPosition.vertical)
                else:
                    #y should be added to horizontal position
                    yaxis = self.addAxis(yaxis, AxisPosition.horizontal)
            else:
                if yaxis.position.isHorizontal():
                    #x should be added to vertical position
                    xaxis = self.addAxis(xaxis, AxisPosition.vertical)
                else:
                    #x should be added to horziontal position
                    xaxis = self.addAxis(xaxis, AxisPosition.horizontal)
                pass
            pass
        elif not xok:
            #neither is added - add both, x to horizontal, y to vertical - where is more space
            xaxis = self.addAxis(xaxis, AxisPosition.horizontal)
            yaxis = self.addAxis(yaxis, AxisPosition.vertical)
        else:
            #both are added already, check that they are ortogonal
            if not xaxis.position.isOrtogonal(yaxis.position):
                print("WARNING: Plot.addXYPlotter - specified axes are not orthogonal on plot", xaxis.position, yaxis.position)
                return
            pass

        #~ vprint("new plotter args, kw", args, kw)
        plotter = plotterClass(self, xaxis, yaxis, *args, **kw)
        self._plotters.append(plotter)
        # ~ plotter.setParentItem(self._da._area)
        plotter.setParentItem(self._da._area)
        
        return plotter
    
    def clearPlotters(self):
        for plotter in self._plotters:
            plotter.cleanup()
            plotter.setParentItem(None)
            plotter.scene().removeItem(plotter)
            self._legend.removeSeries(plotter)
        self._plotters = []
    
    def removePlotter(self, plotter):
        if plotter in self._plotters:
            plotter.cleanup()
            plotter.setParentItem(None)
            plotter.scene().removeItem(plotter)
            self._plotters.remove(plotter)
            self._legend.removeSeries(plotter)
            return True
        return False
       
    def plotter(self, index):
        return self._plotters[index] #TODO: allow named plotters
    
    def clearDecorators(self):
        for dec in self._decorators:
            dec.cleanup()
            dec.setParentItem(None)
            dec.scene().removeItem(dec)
        self._decorators = []
        
    def decorator(self, index):
        return self._decorators[index]

    def removeDecorator(self, decorator):
        if decorator in self._decorators:
            decorator.cleanup()
            decorator.setParentItem(None)
            decorator.scene().removeItem(decorator)
            self._decorators.remove(decorator)
            return True
        return False
        
    #~ def test(self, gritem):
        #~ #add gritem to data area
        #~ self._items.append(gritem)
        #~ gritem.setParentItem(self._da._area)
    pass

class FigureScene(QtWidgets.QGraphicsScene):
    #just for custom contextmenuevent
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self._lastImageFilter = None #TODO: can/should we do this class variable? (so all SavePlots will have it the same)
        
    def contextMenuEvent(self, ev):
        vprint("FigureScene.contextMenuEvent")
        #~ self.setContextMenuPolicy(Qt.Qt.ActionsContextMenu)
        #Todo: this should be part of Plot (each plot can have its own legend)
        #~ ac = QtWidgets.QAction("Legend", self)
        #~ ac.triggered.connect(self._setShowLegend)
        #~ ac.setCheckable(True)
        #~ self._acShowLegend = ac
        #~ vprint(ev.widget())
        #~ vprint(self.views())
        menu = QtWidgets.QMenu()
        
        #~ item = self.itemAt(ev.scenePos().toPoint(), self.views()[0].transform()) #it needs some deviceTransformation which I have no idea what it is - as it is done here, it will not work with multiple views
        item = self.itemAt(ev.scenePos().toPoint(), QtGui.QTransform()) #it needs some deviceTransformation which I have no idea what it is - as it is done here, it will not work with multiple views
        #just a note print(self.views()[0].transform()==QtGui.QTransform())  is True
#        print()
#        print()
#        print()
        print("items under mouse", self.items(ev.scenePos().toPoint(), QtCore.Qt.IntersectsItemShape))
#        print()
#        print()
#        print()

        #~ print("items under mouse, view", self.items(ev.scenePos().toPoint(), deviceTransform=self.views()[0].transform()))
        print("starting item chain")
        #TODO: it has very strange boundaries of items ... 
        #  as if the detection function ignores cosmetic pens, AxisArea graphicsLineItems extend way to much
        # possibly, I can use _area just to recalculate coordinates, but parent Items to PlotArea itself, there will be no need for cosmetic pens then; somewhat annoying, but if it works...
        #TODO: I can try to show highlight against boundingRect of top most item under mouse
        # another way is to change Z order of PlotAreas
        handled = []
        for item in self.items(ev.scenePos().toPoint(), QtCore.Qt.IntersectsItemShape):
            print(item, item in handled)
            while item:
                if item not in handled: 
                    handled.append(item)
                    # Try run items context if it has one
                    #~ go through item parent chain and add to menu
                    try:
                        #todo_ possibly add name of item as title of menu part (QMenu.addSection)
                        print(item)
                        acs = item.contextMenuActions(ev)
                        for name, callback in acs:
                            if callback is None:
                                menu.addSection(name)
                            else:
                                menu.addAction(name, callback)
                        menu.addSeparator()
                    except AttributeError:
                        pass
                    except:
                        import traceback
                        print("FigureScene.contextMenuEvent: item", item)
                        traceback.print_exc()
                        pass
                item = item.parentItem()
                
        acSaveImage = menu.addAction("Save Image", self.saveImage)
        #~ QAction *markAction = menu.addAction("Mark");
        selectedAction = menu.exec(ev.screenPos())
        #~ if selectedAction is acSaveImage:
            #~ self.saveImage()
    
        #Todo: also allow saving of individual plots - we can just draw the Plot widget
        #~ acSaveImage=QtWidgets.QAction("Save Image", self)
        #~ acSaveImage.triggered.connect(self.saveImage)
        #~ self.addAction(acSaveImage)

    def saveImage(self, filepath=None, width=None, height=None, *args, **kw):
        #if filepath is not writeable, ask for new one
        vprint("FigureScene.saveImage:", filepath, width, height)
        try:
            with open(filepath,'ab') as temp:
                pass
            #we can write to this file
        except:
            import traceback
            print("FigureScene.saveImage: cannot write to path", filepath)
            traceback.print_exc()
            print()
            
            #ask for filepath
            #~ filters = [ "PNG nice and big raster (*.png)", \
                    #~ "PNG nice and small raster (*.png)", \
                    #~ "JPEG small and dirty (*.jpg)",\
                    #~ "PDF nice and big vector (*.pdf)",\
                    #~ "PDF nice and big vector (no hairlines) (*.pdf)"] 
            filters = ["Images (*.png *.jpg *.bmp *.ppm *.xbm *.xpm)", "WIP PDF (*.pdf)", "PNG (*.png)", "JPEG (*.jpg)"]
            res = QtWidgets.QFileDialog.getSaveFileName(filter=";;".join(filters), directory=self.views()[0].window().windowFilePath(), initialFilter=self._lastImageFilter)
            #~ vprint(res)
            filepath = res[0]
            if filepath=="": return
        
            self._lastImageFilter = res[1]
            #file dialog under linux does not add extensions automatically :(
            name, ext = op.splitext(filepath)
            if len(ext)==0:
                i = self._lastImageFilter.find("*")+1
                filepath += self._lastImageFilter[i:i+4] #assuming only length 3 extensions used in filters above

        #TODO: allow some additional options like finer resolution than widget size, etc
        #   we can also render QGraphicsScene directly (probably will render whole scene and not just part visible in view, but that is the same for us)
        #   it will follow the size of the bitmap, *but* since I use cosmetic pens and text, these will not be scaled :(

        #decorate these with busy cursor, it can take a while :D
        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        NN = kw.pop("subpixel", 1) #not used anymore, but digest it for backward compatibility
        if filepath.endswith(".pdf"):
            self.savePdf(filepath, *args, **kw)
        else:
            sw = self.width()
            sh = self.height()
            
            if width is not None or height is not None:
                #we have at least one
                if width is None:
                    #calculate width 
                    width = height/sh*sw
                elif height is None:
                    height = width/sw*sh
            else:
                width = sw*NN
                height = sh*NN
            canvas = QtGui.QPixmap(int(width), int(height))
            canvas.fill()
            p = QtGui.QPainter(canvas)
            # ~ self.render(p, self.sceneRect(), QtCore.QRectF(0, 0, self.width(), self.height()))
            self.render(p)
            p.end()
            vprint("\t\tsave canvas", filepath)
            canvas.save(filepath, *args, **kw)
        QtWidgets.QApplication.restoreOverrideCursor()
        pass
    
    #TODO: this only works (almost) properly if the widget is shown before saving (i.e. run pqr.show() before widget.saveSvg)
    def saveSvg(self, filepath=None):
        #if filepath is not writeable, ask for new one
        try:
            with open(filepath,'ab') as temp:
                pass
            #we can write to this file
        except:
            #ask for filepath
            #~ filters = [ "PNG nice and big raster (*.png)", \
                    #~ "PNG nice and small raster (*.png)", \
                    #~ "JPEG small and dirty (*.jpg)",\
                    #~ "PDF nice and big vector (*.pdf)",\
                    #~ "PDF nice and big vector (no hairlines) (*.pdf)"] 
            res = QtWidgets.QFileDialog.getSaveFileName(filter="SVG files (*.svg)", directory=self.views()[0].window().windowFilePath())
            #~ vprint(res)
            filepath = res[0]
            if filepath=="": return
        
            self._lastImageFilter = res[1]
            #file dialog under linux does not add extensions automatically :(
            name, ext = op.splitext(filepath)
            if len(ext)==0:
                filepath += ".svg"
                i = self._lastImageFilter.find("*")+1
                filepath += self._lastImageFilter[i:i+4] #assuming only length 3 extensions used in filters above

        generator = QtSvg.QSvgGenerator()
        generator.setFileName(filepath)
        generator.setSize(QtCore.QSize(self.width(), self.height())) #looks like these are svg pt (unit), but the scene is not scaled
        #it is better when layout is activated prior save, but still the scale is not correct
        #~ generator.setViewBox(QtCore.QRect(0, 0, 400, 400));
        #~ generator.setTitle(tr("SVG Generator Example Drawing"));
        #~ generator.setDescription(tr("An SVG drawing created by the SVG Generator "
                                    #~ "Example provided with Qt."));
        painter = QtGui.QPainter()
        NN = 1
        painter.begin(generator)
        #~ self.render(painter, QtCore.QRectF(), QtCore.QRectF(0, 0, self.width()*NN, self.height()*NN))
        self.render(painter)
        painter.end()
        pass    
    
    def savePdf(self, filename, creator=None, pageLayout=None, pageMargins=None, pageOrientation=None, pageSize=None, resolution=None, title=None):
        #this is just conveninece, nothing stops you from rendering the widget whereever
        #TODO: this does not work perfectly (see e.g. test at the end of this file)
        canvas = QtGui.QPdfWriter(filename)
        if creator is not None: canvas.setCreator(creator)
        if pageLayout is not None: canvas.setPageLayout(pageLayout)
        if pageMargins is not None: canvas.setPageMargins(pageMargins)
        if pageOrientation is not None: canvas.setPageOrientation(pageOrientation)
        if pageSize is not None: canvas.setPageSize(pageSize)
        if resolution is not None:canvas.setResolution(resolution)
        if title is not None: canvas.setTitle(title)
        self.render(QtGui.QPainter(canvas))
        pass
        
        

class FigureWidget(QtWidgets.QGraphicsView):
    """
    has save/load interface similar to SaveLoadElement or Params
    keep in mind that save settings should only store UI defined values which are not known at compile time and can only be used on exactly the same Plot layout, with loaded data and added Plotters
    """
    def __init__(self, *args, addFirstPlot=True, figureTitle=None, floating=False, **kw):
        super().__init__(*args) #Qt stuff
        
        self._scene = FigureScene() #QtWidgets.QGraphicsScene()
        #~ scene->setItemIndexMethod(QGraphicsScene::NoIndex);
        width = kw.pop("width", GS.width)
        height = kw.pop("height", GS.height)
        if width is None or height is None:
            self._resizeScene = True
            width, height = 100, 100
        else:
            self._resizeScene = False
        
        self._scene.setSceneRect(0, 0, width, height)
        self._scene.changed.connect(self.updateScene)
        self.setScene(self._scene)
        #~ yourGraphicsView->setFixedSize(width+2*yourGraphicsView->frameWidth(), height+2*yourGraphicsView->frameWidth());
        #~ setCacheMode(CacheBackground);
        #~ setViewportUpdateMode(BoundingRectViewportUpdate);
        #~ setRenderHint(QPainter::Antialiasing);
        #~ setTransformationAnchor(AnchorUnderMouse);
        #~ scale(qreal(0.8), qreal(0.8));
        #~ setMinimumSize(400, 400);
        #~ setWindowTitle(tr("Elastic Nodes"));        
        
        self._widget = QtWidgets.QGraphicsWidget()
        self.layout = QtWidgets.QGraphicsGridLayout() #public layout for placing plots
        self.layout.setContentsMargins(0,0,0,0)
        self.layout.setSpacing(0) #TODO: expose to GS
        # ~ self._widget.setAutoFillBackground(True)
        
        self._llabel = LabelGraphicsLayoutItem(figureTitle, fontSize=14, center=False) 
        #~ if figureTitle is None: self._llabel.setMaximumSize(0,0) #Todo: set to invalid size to reenable sizeHint, alternatively hide elements inside of LabelGraphicsLayoutItem, or remove _llabel from _layout
        self._layout = QtWidgets.QGraphicsLinearLayout(QtCore.Qt.Vertical)
        self._layout.addItem(self._llabel)
        self._layout.addItem(self.layout)
        
        
        self._widget.setLayout(self._layout)
        self._widget.resize(width,height)
        self._scene.addItem(self._widget)
        self._scene._widget = self #TODO: just for testing
        
        #TODO: proper resizing
        #~ self.show() #for some reason has to be called before addPlot, otherwise plot will not show plotters (possilby resize thing?)
        #~ self.resize(300,300) #this is not enough
        #~ self._widget.show() #this is not enough
        
        self.setWindowTitle("Figure Widget" if figureTitle is None else str(figureTitle))
        if floating: self.setWindowFlags(Qt.Qt.Dialog) #I will try to make the window floating in tiling window managers (at least as an option) - use case: fixing window size for picture output
        #TODO: add Figure Title to figure itself
        self._plots = []
        #~ self._rects = [] #for updateScene testing 
        if addFirstPlot:
            vprint("adding first plot", **kw)
            self.addPlot(1, 1, **kw) #TODO: I will add first to 1,1 so that there is some space around, but it would be best if grid layout supported column/row insertions
        
        
        vprint("\tFigureWidget.show")
        self.show() #lets show widget by default
        self.updateGeometry() #TODO: try to force recalculation of title label item size
        self._widget.updateGeometry()
        pass

    def setSceneSize(self, size, size2=None):
        #we could ignore scene and do everything in _widget, since Plots are supposed to be organized in _widget.layout, but lets to it in scene, so that it is possible to ignore the layout if needed
        if size2 is None:
            #size is the smaller of width and height, keep aspect ratio of self (graphicsview widget)
            # this shoudl allow nice scaling 
            ...
        else:
            #set scene size to size, size2; when sel resizes, it should scale scale scene and keet its aspect ratio
            ...
        pass

    #~ def updateGeometry(self, *args):
        #~ vprint("\tFigureWidget.updateGeometry", args, self.geometry())
        #~ self._llabel.updateGeometry() #because _llabel derives its maximum size from Plot geometry, we need to let it know it has changed
        #~ return super().updateGeometry(*args)

    #~ def updateScene(self, rects):
        #~ #this will mark updated regions - for testing only
        #~ vprint("FigureWidget.updateScene", len(rects))
        #~ if len(rects)>0:
            #~ self._scene.changed.disconnect(self.updateScene)
            #~ for it in self._rects:
                #~ it.setParentItem(None)
                #~ self._scene.removeItem(it)
            #~ self._rects = []
            #~ for it in rects:
                #~ r = self._scene.addRect(it)
                #~ pen = r.pen()
                #~ pen.setCosmetic(True)
                #~ pen.setColor(Qt.Qt.green)
                #~ r.setPen(pen)
                #~ self._rects.append(r)
            #~ res =  super().updateScene(rects)
            #~ self._scene.changed.connect(self.updateScene)
            #~ self.repaint()
            #~ vprint("\tend1")
            #~ return res
        #~ vprint("\tend2")
        #~ return super().updateScene(rects)
    
    def saveSettings(self, **kw):
        #we want to save params of all plots, but these have to be somehow identified
        #loadSetting should not try to create Plots, only change 
        res = {}
        for i in range(len(self._plots)):
            ttt = self._plots[i].saveSettings()
            if len(ttt)>0: res["plot"+str(i)] = ttt
        return res
    
    def loadSettings(self, settings, *kw):
        for i in range(len(self._plots)):
            key = "plot"+str(i)
            if key in settings: self._plots[i].loadSettings(settings[key])
        pass
    
    def resizeEvent(self, ev):
        # ~ vprint("\tFigureWidget.resizeEvent", ev.size(), self.geometry())
        if self._resizeScene:
            s = ev.size()
            self._scene.setSceneRect(0, 0, s.width()-2*self.frameWidth(), s.height()-2*self.frameWidth())
            #~ self._rect.setRect(0+1, 0+1, s.width()-2*self.frameWidth()-2, s.height()-2*self.frameWidth()-2)
            #~ self._widget.resize(s.width(), s.height())
            #~ vprint("\tFigureWidget.resizeEvent call _widget.resize")
            self._widget.resize(s.width()-2*self.frameWidth(), s.height()-2*self.frameWidth())
        else:
            self.fitInView(self._scene.sceneRect(), QtCore.Qt.KeepAspectRatio)
        
        # ~ self._llabel.updateGeometry() #because _llabel derives its maximum size from Plot geometry, we need to let it know it has changed (Todo: this is weird, if you place it after _widget.resize, it will not update label size
        res = super().resizeEvent(ev)
        return res
    
    def addPlot(self, *args, plotClass=Plot, **kw):
        #plot can be basically any instance derived from QGraphicsWidget
        #todo: name will be used as key for settings storage
        plot = plotClass(**kw)
        self._plots.append(plot)
        self.layout.addItem(plot, *args)
        return plot
        
    def plot(self, index=-1):
        return self._plots[index]
    
    def saveImage(self, *args, **kw):
        if kw.pop("activateLayout", True): self.layout.activate() #in case image is saved before widget is shown
        #TODO: almost there, but this will not properly calculate space for figure title
        return self._scene.saveImage(*args, **kw)
    
    def saveSvg(self, *args, **kw):
        if kw.pop("activateLayout", True): self.layout.activate() #in case image is saved before widget is shown
        return self._scene.saveSvg(*args, **kw)
    #Todo: plotAt(self, mousePos)

    #can we delegate anything that we do not know about to the first plot? if there is only one plot then widget is basically just this plot
    
    pass
    

#interactive mode
# - has to create app if not created already
# - has to create Qt objects
# - has to display them
# there are two user cases
# - there is Qt app already and we want to just create additional window with plot
# - there is not Qt app

#These two should be minimal needed,
# but you have to keep all objects alive - reasigning variables might lead to garbage collector to destroy the windows (as it would destroy pqrApp if it was not global)
#For example
#  pqr.init() #interactive mode
#  w1 = pqr.FigureWidget(figureTitle=1) #will add 1,1 plot automatically
#  w1 = pqr.FigureWidget(figureTitle=2) #will add 1,1 plot automatically
#  pqr.show()
# will only display one figure window with title "2".
# Since this might (not)destroy some Qt internal C datastructures which python knows nothing about, this might lead to memory leaks, segfaults, or other problems.

pqrApp = None
def init(): #TODO: we should store figures internally in this mode so that they are not destroyed by garbage collector
    global pqrApp
    vprint("pqr.init")
    pqrApp = QtWidgets.QApplication(sys.argv)
    vprint("\tapp", pqrApp)
    
def show():
    global pqrApp
    vprint("pqr.show")
    vprint("\tapp", pqrApp)
    #TODO: test if any window is created and do not run if none is
    res = pqrApp.exec_()
    vprint("pqr.show done")
    return res

def plotDD(data, x=None, y=None, *args, **kw):
    #args are for figureWidget, kw are for Plot
    #data.shape == len(y), len(x) (or one less depending on whether axes show values on data values or between them (pixel center of pixel edges))
    f = FigureWidget(*args, **kw)
    p = f.plot()
    p.addDDPlotter().setData(data, x, y)
    p.zoomToData()
    return f

#TODO: (session) global settings for e.g. plot title fonts - see pqr2
#TODO: keep track of created widgets internally so that they will not be destroyed by garbage collector

if __name__=="__main__":
    app = QtWidgets.QApplication(sys.argv)
    GS.loadPreset("scaling")
    GS.tickLength = -3
    GS.colorBarWidth = 30

    window = FigureWidget(addFirstPlot=False, figureTitle="TEST FIGURE TITLE")
    #~ window.show()

    #tests
    if False:
        #two datasets with common x axis and individually scaled y axes
        #~ plot1 = Plot()
        #~ window.addPlot(plot1, 1, 1)
        plot1 = window.plot() #from addFirstPlot
        #~ plot1 = window.addPlot(0,0)
        print(plot1)
        x = [0.2, 0.5, .6, .3]
        y = [.1, .2, .3, .4]
        y2 = [10, 8, 2, 1]
        plot1.addSeriesPlotter(None, None, color="blue").setData(x, y) #this will add x1 and y1
        plot1.addSeriesPlotter(plot1.x1, None, color="red").setData(x, y2) 
        plot1.zoomToData()
        print("TEST", plot1._plotters)
        print("TEST", window.plot()._plotters)
        
    #two datasets with common x axis and individually scaled y axes with common level
    
    #dataset plotted against y,x
    if False:
        #two datasets with common x axis and individually scaled y axes
        plot1 = window.addPlot(3, 1)
        x = np.arange(10)
        y = np.exp(-(x-5)**2)
        
        plot1.addSeriesPlotter(None, None, color="blue").setData(x, y) #this will add x1 and y1
        plot1.addSeriesPlotter(plot1.y1, plot1.x1, color="red").setData(x, y) 
        plot1.zoomToData()
    
    #enum transformation - TODO: label handling
    
    #numenum transformation
    if False:
        plot1 = window.addPlot(2, 2)
        x = np.array([-0.1, 0, .1, 0.2, .3, .5, .8, 2, 3, 5, 10, 15, 30, 50, 200, 500, 1000, 2000, 5000, 10000, 20000])
        y = np.exp(-(x-5)/500)
        
        plot1.addSeriesPlotter(NumEnumTransform(stops=x), None, color="blue").setData(x, y) #this will add x1 and y1
        #actually stops need not be x
        x2 = [0.05, 0.7, 1, 7, 14, 100, 600, 3000, 30000] #however it *should* be inside stops, otherwise you will get this error (30000 is shown as 20000)
        y2 = x2
        plot1.addSeriesPlotter(plot1.x1, None, color="red").setData(x2, y2)
        plot1.zoomToData()
        
    #AOverTransformation, plot nm and icm with shared axis zoom
    if False:
        plot1 = window.addPlot(3,2)
        x = np.arange(300, 900)
        y = np.exp(-(x-400)**2/100) + np.exp(-(x-750)**2/200)
        plot1.addSeriesPlotter(None, None, color="blue").setData(x, y) #this will add x1 and y1
        at = plot1.addAxis(AOverTransform(A=1e7), AxisPosition.top) #TODO: add some kw to link zoom
        plot1.x1.zoomRangeLChanged.connect(at._setZoomRange)
        plot1.x2.zoomRangeLChanged.connect(plot1.x1._setZoomRange)
        plot1.x1.label = "Test x1"
        plot1.y1.label="TEST y1"
        plot1.x2.label="TEST x2 really long and long axis label"
        plot1.addAxis(None, AxisPosition.R, label="TEST RIGHT really long and long axis label")
        plot1.zoomToData()
        # ~ test = plot1.addDecorator(OriginDecorator(plot1, plot1.x1, plot1.y1))
        # ~ test.origin = 350, 0.5
        # ~ plot1._panDataArea(QtCore.QPointF(0,0), QtCore.QPointF(1,0))
    
    if False:
        p = window.addPlot(0,0)
        ttt = p.addSeriesPlotter()
        ttt.setData([0,1], [1,1])
        print("\t\t ttt", ttt)
        print("\t\t ttt", ttt._gritems[0])
        print("\t\t ttt", ttt._pen.widthF())
        print("\t\t ttt", ttt._gritems[0].shape())
    
    # ~ globalPen.setColor(QtGui.QColor("red"))
    
    if True:
        plot1 = window.addPlot(2,1, title="TEST Plot title")
        dd = plot1.addDDPlotter(contours=5)
        #~ data = np.array([[0,0,0], [1,1,1], [2,2,2], [3,3,3], [4,4,4], [5,5,10]])
        x = np.arange(0, 100+1e-12, 5)
        y = np.arange(0, 100+1e-12, 2)
        data = np.outer(np.exp(-(x-50)**2/100), np.exp(-(y-50)**2/1502))*15
        #~ dd.setData(data, y, x, xAlignment="leftEdge", yAlignment="rightEdge")
        
        plot1.addDecorator(RectDecorator(plot1, plot1.b1, plot1.l1, color="purple"))
        
        if True: #what if we add some nan or inf into data?
            #ok, we can get nans and infs from datarange calculation (slight problem is that both nan and inf appear the same (black) in image
            #however contours do not play nice if nan or inf gets involved (it can be fine if such pixel is far from contours)
            #  ui is very slow
            data[10,0] = np.nan
            # ~ data[0,0] = np.nan
            # ~ dd._zaxis.colorMap = cmBlueRed
            data[10,10] = -np.inf
            # ~ data[10,0] = np.inf
            # ~ data[0,0] = np.inf
            
            #normal
            # ~ generate contours [ 0.75   4.125  7.5   10.875 14.25 ] [-1.0, 101.0] [-2.5, 102.5]
            # ~ calculateContours stats: levels 5 contours 5 points 350
            
            #np.inf at 0,0 (outside of way on edge)
            # ~ generate contours [ 0.75   4.125  7.5   10.875 14.25 ] [-1.0, 101.0] [-2.5, 102.5]
            # ~ pqr.py:1677: RuntimeWarning: invalid value encountered in double_scalars
                # ~ x = (level-A[i, j])/(A[i+1,j]-A[i,j])
            # ~ pqr.py:1778: RuntimeWarning: invalid value encountered in double_scalars
                # ~ x = (level-A[i+1, j])/(A[i+1,j+1]-A[i+1,j])
            # ~ calculateContours stats: levels 5 contours 5 points 360
            
            #np.nan at 0,0 (outside of way on edge)
            # ~ generate contours [ 0.75   4.125  7.5   10.875 14.25 ] [-1.0, 101.0] [-2.5, 102.5]
            # ~ pqr.py:1645: RuntimeWarning: invalid value encountered in greater
                # ~ B = A>level
            # ~ calculateContours stats: levels 5 contours 5 points 350

            #np.inf at 10,10 - this makes problems
            # ~ generate contours [ 0.75   4.125  7.5   10.875 14.25 ] [-1.0, 101.0] [-2.5, 102.5]
            # ~ pqr.py:1778: RuntimeWarning: invalid value encountered in double_scalars
              # ~ x = (level-A[i+1, j])/(A[i+1,j+1]-A[i+1,j])
            # ~ pqr.py:1677: RuntimeWarning: invalid value encountered in double_scalars
              # ~ x = (level-A[i, j])/(A[i+1,j]-A[i,j])
            # ~ calculateContours stats: levels 5 contours 5 points 360

            #nan at 10,0 - this also make problems
            # ~ generate contours [ 0.75   4.125  7.5   10.875 14.25 ] [-1.0, 101.0] [-2.5, 102.5]
            # ~ pqr.py:1645: RuntimeWarning: invalid value encountered in greater
              # ~ B = A>level
                    # ~ calculateContours stats: levels 5 contours 5 points 353
            
            #but apparently the problem is not with number of contours generated
            #is it drawing?   inf at 10,10 draws some very weird contour
            # does it try to go to inf and draws all the way across the scene? (will be clipped by Plot so no way of knowing)
            # it seems UI has problems when there are nans in contour QPolygonF
            pass
        
        if True: # test alpha DDplotter
            alpha = np.outer(np.ones_like(x), np.arange(len(y))/len(y)*255 )
        else:
            alpha = None

        if True:
            DataConstantTicker( plot1.b1, positions=np.arange(-10, 80, 10))

        dd.setData(data, y, x, alpha=alpha)
        #~ dd.setData(data)
        plot1.zoomToData()
        print("colorMap", plot1.r1.colorMap)
        plot1.r1.label="Ahoj"
        
        tracer = plot1.addDecorator(TracerDecorator(plot1, plot1.b1, plot1.l1))
        #~ tracer.setData(50,50)
        tracer.origin = 50,50
        tracer.color = "orange"
        
    if False:
        #test color bar on different sides
        plot1 = window.addPlot(2,2, title="TEST Plot title")
        #~ plot1.addAxis(None, AxisPosition.top, ColorAxisArea, globalZoom=False)
        #~ za = plot1.addAxis(None, AxisPosition.right, ColorAxisArea, globalZoom=False)
        za = plot1.addAxis(None, AxisPosition.bottom, ColorAxisArea, globalZoom=False)
        za = plot1.addAxis(None, AxisPosition.bottom, ColorAxisArea, globalZoom=False)
        #~ za = plot1.addAxis(None, AxisPosition.left, ColorAxisArea, globalZoom=False)
        #~ za = plot1.addAxis(None, AxisPosition.left, ColorAxisArea, globalZoom=False)
        dd = plot1.addDDPlotter(contours=5, zaxis=za)
        #~ data = np.array([[0,0,0], [1,1,1], [2,2,2], [3,3,3], [4,4,4], [5,5,10]])
        x = np.arange(0, 100+1e-12, 5)
        y = np.arange(0, 100+1e-12, 2)
        data = np.outer(np.exp(-(x-50)**2/100), np.exp(-(y-50)**2/1502))*15
        #~ dd.setData(data, y, x, xAlignment="leftEdge", yAlignment="rightEdge")
        dd.setData(data, y, x)
        #~ dd.setData(data)
        plot1.zoomToData()
        

    #basics are done
    #~ next:
        #~ then proceed to 2D color plots
        #~ and contours
    #~ window.show()
    
    if False:
        #test scaling from not parent item
        parent = QtWidgets.QGraphicsRectItem(300, 300, 300, 200)
        window._scene.addItem(parent)
        
        area = QtWidgets.QGraphicsRectItem(0, 0, 1, 1, parent)
        p = area.pen()
        p.setCosmetic(True)
        p.setColor(QtGui.QColor("red"))
        area.setPen(p)
        area.setTransform(QtGui.QTransform.fromScale(300, 200), False)
        print(area.rect())
        print(area.mapRectToParent(area.rect()))
        print("scene", area.mapRectToScene(area.rect()))

        test = QtWidgets.QGraphicsRectItem(0.25, 0.25, .5, .5, parent)
        p = test.pen()
        p.setCosmetic(True)
        p.setColor(QtGui.QColor("green"))
        test.setPen(p)

        parent.setPos(500, 300)
        print("test", area.mapRectToParent(test.rect()))
        #in principle we can use area to calculate what test would be like in parent coordinates
        # and create test directly in parent, avoiding area massive scaling (and therefore cosmetic pens)
        # however this will not work well with changing window size, we will have to redraw/reposition everything


    #test save SVG
    if False:
        window.resize(500,500)
        window.show()
        window.saveSvg("test2.svg")

    if False:
        #test image saving
        # ~ window.resize(800,800)
        # ~ window.saveImage("test5.png")
        window.saveImage("testAIP.png", width=4016)
        # ~ window.saveImage("test4.png", width=900)
        # ~ window.saveImage("test2.pdf", subpixel=2) #will not work well
        

    res = app.exec_()
    print("settings")
    print(window.saveSettings())
    sys.exit(res)
