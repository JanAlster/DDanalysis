class DDPlotter(XYPlotter):
    """
    plot 2D arrays as image
    with possible contours, tracers, colorbar axis, etc
    probably needs to use NumEnum for axes otherwise cannot use same sized pixels (unless it can be guaranteed that axes are equispaced)
    
    __init__ needs to know against which axes of Plot to plot 
        but it does not know the data yet
        so it cannot gaurantee that data will be equispaced
        
        but to keep things simple, we want equispaced pixels, i.e. datapoints
        we can do that with NumEnum transform
        
        but there is not guarantee that supplied axes will use this transform
        
        we might just want from user axes positions (and plot) and create them here, but we might want to share some existing axes
        
        we might restrict user to use such data and transforms so that data are equispaced on linear space
           - otherwise results is undefined (which means that wrong values will be reported on axes positions)
           - there is not need to use strictly NumEnum  - you can generate data on 1/x and use AOverTransform
       
        we might subclass AxisArea to one where transform is set by first DataPlotter, but that will not help
        what if two DDPlotters want to share axes, but each is of different kind (described below)
        
        typically a CountourPLotter should share axes with DDPLotter
       
        otherwise I would need to implement creation of different sized pixels (might work with vector graphics, but creating big 2D array out of rectangles is probably stupid)
           - we might use smallest data spacing (or rather smallest common divisor), but that could lead to huge bitmaps
           - I also do not want to adjust bitmaps for every zoom level (although it is possible, esp. if interpolation is desired)
           - such things are for later
        
    data could be specified with axes positions that are 
       - in the middle of each pixel - edge pixels should be only half size (which might be hard to do, I will not likely make bitmaps twice the size just to do that)
                (e.g. data value is value at this point, but we want to plot everywhere and our resolution is limited)
                (this is equivalent to scatter plot with no connecting lines between points, where points are so big that they touch)
                (? if we zoom in so that pixels are big, should we interpolate?, should be report interpolated axis position?)
                 - it should be one pixel, one axis position, not interpolation of those
                 - this is actually more Enum transform use case 
       - on the edges of each pixel - x and y coordinates has to be one bigger than data array
                (e.g. data value is average (or sum) value in between interval)
                 - here we should probably interpolate axes positions
                 - so NumEnum seems more appropriate
        shall we allow switching between these two types within one instance? should they be separate classes?
       
        Frankly I would let user do as he/she/it wishes
        use data range of supplied 
        then convert data to center of pixel kind
        (using NumEnum with stops at pixel edges should interpolate if that is desired)
        (using Enum will keep uninterpolated pixel position)
        if they mix things up ... they will have to live with it (I might adjust as it goes)
    """

    def _plot(self, dd, x=None, y=None):
        """
        dd is 2D array of data such that
            dd.shape == len(y), len(x)
        or 
            dd.shape == len(y)-1, len(x)-1
            
        x and y are data unit coordinates of individual data points
        unlike for lines, both x and y should be monotnonous (ascending?)
        and you should ensure that they are matched with supplied axes transform to have datapoints equispaced of linear space
        (use EnumTransform or NumEnumTransform if not sure)
        
        z is used to map data amplitude to color levels, each level is assigned color by colormap
           color space is (0,1) divided into N discreet levels, each assigned a color (which might change if colormap changes)
           Z AxisArea will use its transform to map data amplitudes to linear space and zoom/pan to map resulting levels to colors
           
           this will not work so simply - normal axis uses zoom/pan to select a (Arbitrary) range from linear space 
                and map that to 0,1
            color axis should map not change image data, only colorTable
            which means that data level in linear space has to be mapped to (0,255) (8bit colorTable lookup index) 
                - actually this should respect global zrange if specified 
                - but than there might not be enough levels for data amplitudes and the plot will look ugly
            and color zoom maps part of colorMap range to those (0, 255) levels
            now the question is should the z transform data amplitudes, or colormap?
            
            frankly it should zoom data levels - i.e. mapping of data amplitudes in linearspace to (0, 255) indexes in color lookup table (here we might again reach limits)
            of finite resolution), but that would mean *changing* values of image pixels
               and it should be faster to change color than pixels
               
               we might combine things a bit - do image to 256 levels (scaled to its real data amplitudes range)
                 and scale colorMap to respect the global range
                 maybe have some offset from the global color map
                 otherwise it might not workout if colorAxis is shared between images with different data amplitude ranges
                 basically colorAxis cannot colorTable of image, it can only map data amplitdue linear space to colors
                    image has to updat its colorTable according to what colors correspond to its línear levels (there will be 256 of them, but need not be (0, 255)
               
           of course we can store original data (or their linear space transform), use RGB format, and with every zoom recalculate all colors
                with infinite resolution
                but likely much slower
                
        """
        #we also need a color map
        # but changing colormap might require storing original data :( unless we can invert color map and that would be slow and too complex
        # this might be solved by using color look up table and indexed image format, but we will not be able to do 
        # we will be limited to 8bit color map (should be enough)
        # we will not have alpha channel ...
        # 
        
        #dataAmplitudes -> z.transform -> dataAmpLinear -> select 256 levels in its range (equispaced)
           # - problem with this approach is that if we have data with 256 levels, we will not be able to see fine variations on those levels even if we zoom in
           
        """
        another approch would be to realocate those 256 levels when we zoom
            requiring recalculation of image with every zoom (of color map)
            but allowing fast switch of colormaps
            
        we have actually two scales here
         dataAmp - (0,255) - colorRange
         we might want to change two different things
           - map only part of dataAmp range to (0,255) 
           - map (0, 255) to only part of colorRange
           correspondence of dataAmpLinearSpace to colorRange is governed by (global) dataAmp range
           we mostly want to do the second, but first might be needed if we want to see some details
           *or* we might recalculate dataAmp levels by visible data area, but that would be slow
           *or* we might recalculate dataAmp levels if colorRange is below dataAmp levels
           
           we might have global 256 levels corresponding to actual dataAmp range
           if the colorMap is mapped to smaller range, we might limit ourselves to that range and recalculate levels
           
           so in case that global dataAmp range which colorMap is able to assign color is larger than image dataAmp range, 
           we will limit 256 levels to image dataAmp range and achieve best amp resolution for image
           (of course map a subset of colorMap range to those 256 levels)
           
           if the dataAmp range described by colorMap is smaller that image dataRange, we might recalculate levels to match only
           that range 
            - there is no point to have levels for which there is no color
            - it will improve amp resolution
        
        dataAmplitudes -> z.transform -> dataAmpLinear (stored) -> image dataAmpLinear range
        colorMap zoom and pan -> described dataAmpLinear range -> colormap dataAmpLinear range
        now select intersection of imageDataAmpLinear range and colormapDataAmpLinear range
           divide it into 256 (or 254 + less + more) levels
           create image by mapping levels to stored dataAmpLinear
           create colorTable by mapping levels to colors
           
           change colorTable if levels change or colorMap change

    from the point of view of speed of rendering, then creating QPixmap out of QImage is fastest for RGB32, or ARGB32_premultiplied, others are way slower
    #~ Format_ARGB32_Premultiplied      0.001488, 0.001535, 0.001386
    #~ Format_ARGB32                    0.2962
    #~ Format_RGB32                             0.001313, 0.001223
    #~ Format_Indexed8                  0.2962
    
    but we also need to consider speed of data creation - it is faster to calculate 256 RGB values than ~800x800 (for each point of pixmap)
    for RGB32 we will not have problems with discreet levels
    but my user case is not fast refresh of images (although that would be nice)
    we are likely to redo image only when changing data or changing color map
    
    I will just try to keep things simple for now
     - colorMap defines colors on range (0,1)
     - zaxis zoom and pan maps dataAmp linear space to this range
    
    on another point, what about displaying datasets with more points than pixels on target screen? we need to shrink the pixmap
      somehow - should we try to be smart? e.g. not to miss important values (like peaks?)
        """
        """
        ok, this also does not work 
         - because we will need to replot every timeaxis changes zoom/pan, but data is not available anymore (and it does not make sense to do that)
         - we need to be able to plot the data to linear scale (using Axis.d2p) and handle the zooming and panning differently, without need to 
         - we cannot zoom/pan the whole scene (at least not only), because individual lines might be scaled separately
        """

#TODO: we need to force scene to redraw, because it does not do that automatically 
#   i.e. it redraws contours, but not colorPlots, not sure why?
# if there are more Plots on scene, they can (Strangely) update parts of scene belonging to different plots
# it will not happen if PixmapItem is SmoothTransformed, but that is not what I want (at least not always)
#some forum suggested that using ARGB32 might help, but it is platform specific :( 
# possibly bug in Qt FastTransformation
# problem is obviously that during the transform big endian and little endian mix up somehow
# my original code creates ARGB with A = 0 starting from 24bits, with just R and B values
#    and it is displayed ok
#    but during some refreshes (even when zooming something on different part of scene) these values will be interpreted
#     as BGRA or something, leading to white
#     if I pass A = 255, then I get half pixels wrong - alternating white and black - but the other correct and stable
# I had to change input array for QImage, see ColorMap.color


#how to access QImage pixels (for reading, not sure if also for writing)
def _updateColorBar(self, *args):
    data = np.array([np.r_[0:1.:1j*self._levels]])
    if self.position.isVertical():
        data = data.T
    colors = self.colors(data, False)
    im = QtGui.QImage(colors, data.shape[1], data.shape[0], QtGui.QImage.Format_RGB32)
    #~ im = QtGui.QImage(np.zeros_like(colors), data.shape[1], data.shape[0], QtGui.QImage.Format_RGB32)
    #~ im.bits() #TODO: for some reason unless I do this, it will output corrupted image for some image sizes ???
    #~ im.fill(Qt.Qt.red)
    #~ record = [('b','u1'), ('g','u1'), ('r','u1'), ('unused','u1')] # reverse if big-endian
    #~ array = np.ndarray(data.T.shape, dtype=record, buffer=im.bits())
    w, h = data.shape
    s = im.bits().asstring(w * h * 4)
    arr = np.fromstring(s, dtype=np.uint8).reshape((w, h, 4)) 
    print(np.frombuffer(colors, dtype=np.uint8).reshape((w, h, 4)))
    print(arr)
    
    self._colorBar.setPixmap(QtGui.QPixmap.fromImage(im))
    pass
