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
# ~ import faulthandler
# ~ faulthandler.enable()

from PyQt5 import QtCore, QtWidgets, QtGui, QtSvg
from PyQt5.uic import loadUi


import sys
import os
import os.path as op

"""
I would like ot use QGraphicsAnchorLayout, but there is a problem
with how it assigns free spacing between the widgets. 

It is not documented and behaves in strange ways.

I need a better control over it. But I like ability to keep size and 
positions of selected widgets the same (in one dimension).

So I want to try this approach:
 - have a master widget with QGraphicsAnchorLayout
 - define AnchorEdges (could be simple QGraphicsLayoutItem) which keep
    track of distance between two edges and can serve as anchors for the 
    rest of the widgets
    - best regarded as a widget with fixed size with respect to qt layout
    - the size (in the dimension of interest - determined by attached edges)
       can be: - whatever, decided by qt layout
               - fixed pixel size
               - proportion of master widget size (or direct parent size) - regarding qt layout this is also fixed, but it changes when master widget size changes
               - actually it might be needed that this is larger to avoid problems with constraints
               - so we need alignment, to know where to add the extra space
               - but the extra added size should be as small as possible
               
               - proportion of extra distributed space (this is where qt layout goes boing)
               
               - we have to check if we do not add conflicting size to already connected edges
       
"""

class Port:
    def __init__(self, node=None):
        self.connections = []
        self.node = node

class Node:
    """
    serves as anchor point
    needs to have ability to calculate/cache position
    
    connections - list of all connections using this node as start
    """
    #TODO: possibly merge Nodes if we would want to remove widgets from layout
    
    def __init__(self, portIn, portOut):
        self._portIn = portIn
        self._portOut = portOut
        
        portIn.node = self
        portOut.node = self
        
        self._distance = None #distance from starting edge, to be filled by the layout
        pass
    
    @property
    def connections(self):
        return self._portOut.connections
    
    @property
    def backConnections(self):
        return self._portIn.connections
    
    @property
    def portIn(self):
        return self._portIn
    
    @property
    def portOut(self):
        return self._portOut
    
    def addConnection(self, connection):
        self._portOut.connections.append(connection)
        
    def addBackConnection(self, connection):
        self._portIn.connections.append(connection)
        
    def removeConnection(self, connection):
        self._portOut.connections.remove(connection)

    def removeBackConnection(self, connection):
        self._portIn.connections.remove(connection)
        
    def findPath(self, node):
        #try to find node in any of the connections
        for it in self.connections:
            if it.end is node or it.end.findPath(node):
                return True
        return False
    
    def paths(self):
        print("Node.paths", self, self.connections)
        if len(self.connections)==0:
            return [[]]
        else:
            #TODO: connections to the same end node with the same distance can be ignored (only one should be passed on?)
            return [[it]+path for it in self.connections for path in it.end.paths()]
    
    def split(self):
        """
        split iSSo 
        to  iSOo (self)  and  iOSo (new node)
        return the new node
        """
        nPortIn = Port()
        nPortOut = Port(self)
        other = Node(nPortIn, self._portOut) #this should also redirect node of _portOut to other
        self._portOut = nPortOut
        return other
        
    ...
from enum import Enum
class Alignment(Enum):
    start = 1
    center = 2
    end = 3

class Connection:
    """
    describes connection of Nodes, i.e. placement of a widget between two nodes
    
    start - starting node (left or top) (internally the port is saved in case the Node is split later)
    end - ending node  (right or bottom)
    widget - widget whose size is controlled by this connection
    
    size - desire size (fixed, fixed portion of whole, portion of extra space)
             the effective (resulting) size should be at least this (if there is not conflict)
             so if you do not want to give a size, use 0
        - whole int real should be fixed pixels
        - partial real <1 should be fixed portion of whole
        - imag should be portion of extra space
         (real a imag by melo jit kombinovat, real pak slouzi jako min size)
         size cannot be negative (this is not enforced, but it is likely to break things)
    effective size - size determined from all connections (i.e. a parallel route might enforce larger size)
    
    todo:
    alignment - where to place size in effective size interval if it is larger than size
    depth - how much below the master widget is the connected widget - serves as priority for determining
        conflicting effective sizes (note that widget need to have the final parent before adding it to the layout)
    """
    def __init__(self, start, end, widget, size):
        self.widget = widget
        self.size = size
        
        assert not end.findPath(start), "Cannot connect nodes in reverse"
        
        self._start = None
        self._end = None
        self.start = start
        self.end = end
        
        
    @property
    def start(self):
        return self._start.node

    @property
    def end(self):
        return self._end.node
        
    @start.setter
    def start(self, start):
        if self._start is not None:
            self._start.node.removeConnection(self)
            
        if isinstance(start, Node):
            self._start = start.portOut #we store acess through port in case the node splits
        else:
            self._start = start        
        self._start.node.addConnection(self)

    @end.setter
    def end(self, end):
        if self._end is not None:
            self._end.node.removeBackConnection(self)
        
        if isinstance(end, Node):
            self._end = end.portIn
        else:
            self._end = end
        #register with Nodes?
        self._end.node.addBackConnection(self)
        



class Part:
    """
    horizontal connection (left - right)
    vertical connection (top - bottom) 
    
    mozna zkratky
    left = horizontal.start
    atd.
    
    this might be hidden from widget, could be part of layout
    
    takze widget by nemela pristup ke connections, to by mohla delat jen master
    na druhou stranu potrebuju dat info kam se ma nova widget pridat
     ale to muzu delat pomoci enum, nemusi to byt objekt na widget
    
    
    v podstate chci asi misto master widget delat vlastni layout
    """
    ...

from enum import Enum
class Placement(Enum): #TODO: Flag would be nice, but it is not available in my version of Python
    #with aliases
    beforeStart = bS = 1 #will shift all parallel Parts
    afterStart = aS = 2  #will not shift parallel Parts
    beforeEnd = bE =3    #will not shift parallel Parts
    afterEnd = aE = 4    #will shift all parallel Parts
    parallel = same = 5  #will use the same nodes (create a parallel widget)

class PartsLayout(QtWidgets.QGraphicsLayout):
    def __init__(self):
        super().__init__()
        self.left = Node(Port(), Port())
        self.top = Node(Port(), Port())
        self.right = Node(Port(), Port())
        self.bottom = Node(Port(), Port())
        
        self._connectionsH = {} #widget -> connection map
        self._connectionsV = {} #widget -> connection map
        
       
    #we do not have itemAt, count, etc...
    
    """
    pridavat parts muzeme
    bud rozdelenim widget na dve a novou umistit vpravo nebo vlevo 
       (zbytek layout zustane stejny)
    nebo mezi dve sousedni widget tim se odsune i ostatni paralelni widgety
    nebo udanim dvou Node v kazdem smeru
    
    problem je najit sousedni widget pro druhy zpusob
      Node k tomu neni uzpusobena
      - bude to jedna z backconnections
      - museli bychom si pamatovat i hrany, ne jen Node (pozice hrany ve smeru)
      - mohlo by pomoc rozdeli node, to by ale node musela by z dvou casti (prichozi a odchozi)
       a connections by si muselo pamatovat jen prislusnou cast
       
       iAAo o11i iBBo
       (A, B jsou node, 11 je connection reprezentujici widget where)
       iAAo o22i iCCo o11i iBBo (add 22 after start)
       iACo o22i iCAo o11i iBBo (add 22 before start)
       iAAo o11i iCCo o22i iBBo (add 22 before end)
       iAAo o11i iBCo o22i iCBo (add 22 after end)
    """    
    def addPart(self, widget, size, whereH=None, howH=Placement.afterStart, whereV=None, howV=Placement.afterStart):
        """
        widget - what to add
        where - neighbouring widget where the new one should be added, None for master widget of the layout
        how - how it should be added (beforeStart, afterStart, beforeEnd, afterEnd)
        """
        #TODO: size and alignment shoudl be also H and V 
        #note that adding beforeStart on None should change all starting _connections!
        if whereH is None:
            if howH in (Placement.afterStart, Placement.beforeEnd, Placement.same):
                #we do not change border 
                conH = Connection(self.left, self.right, widget, size)
            elif howH == Placement.beforeStart:
                #if there is already a connection between self.left and self.right, 
                #  we need to create a new left
                #otherwise it is not needed
                if self.left.findPath(self.right):
                    oldLeft = self.left
                    self.left = Node(Port(), Port())
                    conH = Connection(self.left, oldLeft, widget, size)
                else:
                    conH = Connection(self.left, self.right, widget, size)
            elif howH == Placement.afterEnd:
                #see beforeStart
                if self.left.findPath(self.right):
                    oldRight = self.right
                    self.right = Node(Port(), Port())
                    conH = Connection(oldRight, self.right, widget, size)
                else:
                    conH = Connection(self.left, self.right, widget, size)
            else:
                raise ValueError("Unknow value for howH", howH)
        else:
            #now we need to keep the existing widget
            """
            iAAo o11i iBBo
            (A, B jsou node, 11 je connection reprezentujici widget where)
            iAAo o22i iCCo o11i iBBo (add 22 after start)
            iACo o22i iCAo o11i iBBo (add 22 before start)
            iAAo o11i iCCo o22i iBBo (add 22 before end)
            iAAo o11i iBCo o22i iCBo (add 22 after end)            
            """
            whereCon = self._connectionsH[whereH]
            A = whereCon.start
            B = whereCon.end

            if howH == Placement.beforeStart:
                C = A.split()
                conH = Connection(A, C, widget, size)
            elif howH == Placement.afterStart:
                #we need to modify whereCon, remove it from A
                C = Node(Port(), Port())
                whereCon.start = C
                conH = Connection(A, C, widget, size)
            elif howH == Placement.beforeEnd:
                C = Node(Port(), Port())
                whereCon.end = C
                conH = Connection(C, B, widget, size)
            elif howH == Placement.afterEnd:
                C = B.split()
                conH = Connection(B, C, widget, size)
            elif howH == Placement.same:
                conH = Connection(A, B, widget, size)
            else:
                raise ValueError("Unknow value for howH", howH)

        self._connectionsH[widget] = conH
        
        # the same for V
        if whereV is None:
            if howV in (Placement.afterStart, Placement.beforeEnd, Placement.same):
                #we do not change border 
                conV = Connection(self.top, self.bottom, widget, size)
            elif howV == Placement.beforeStart:
                #if there is already a connection between self.left and self.right, 
                #  we need to create a new left
                #otherwise it is not needed
                if self.top.findPath(self.bottom):
                    oldTop = self.top
                    self.top = Node(Port(), Port())
                    conV = Connection(self.top, oldTop, widget, size)
                else:
                    conV = Connection(self.top, self.bottom, widget, size)
            elif howV == Placement.afterEnd:
                #see beforeStart
                if self.top.findPath(self.bottom):
                    oldBottom = self.bottom
                    self.bottom = Node(Port(), Port())
                    conV = Connection(oldBottom, self.bottom, widget, size)
                else:
                    conV = Connection(self.top, self.bottom, widget, size)
            else:
                raise ValueError("Unknow value for howV", howV)
        else:
            #now we need to keep the existing widget
            """
            iAAo o11i iBBo
            (A, B jsou node, 11 je connection reprezentujici widget where)
            iAAo o22i iCCo o11i iBBo (add 22 after start)
            iACo o22i iCAo o11i iBBo (add 22 before start)
            iAAo o11i iCCo o22i iBBo (add 22 before end)
            iAAo o11i iBCo o22i iCBo (add 22 after end)            
            """
            whereCon = self._connectionsV[whereV]
            A = whereCon.start
            B = whereCon.end

            if howV == Placement.beforeStart:
                C = A.split()
                conV = Connection(A, C, widget, size)
            elif howV == Placement.afterStart:
                #we need to modify whereCon, remove it from A
                C = Node(Port(), Port())
                whereCon.start = C
                conV = Connection(A, C, widget, size)
            elif howV == Placement.beforeEnd:
                C = Node(Port(), Port())
                whereCon.end = C
                conV = Connection(C, B, widget, size)
            elif howV == Placement.afterEnd:
                C = B.split()
                conV = Connection(B, C, widget, size)
            elif howV == Placement.same:
                conV = Connection(A, B, widget, size)
            else:
                raise ValueError("Unknow value for howV", howV)

        self._connectionsV[widget] = conV
        
        #recalculate layout
        self.invalidate()
        pass
    
    def changeSize(self, widget, sizeH=None, sizeV=None):
        """
        In case that some widget wants to change its size
        
        cannot use Qt sizeHint, since we need complex sizes
        unless we can use custom sizehint types
        if that is possible, we can rely on updateGeometry
        """
        #TODO: we need something in place to prevent needless recalculation of layout
        # in case of 
        
        update = False
        
        if sizeH is not None:
            self._connectionsH[widget].size = sizeH
            update = True
        
        if sizeV is not None:
            self._connectionsV[widget].size = sizeV
            update = True
        
        if update:
            # ~ self.activate() #nebo invalidate, uvidime, jak se to bude chovat
            # ~ self.updateGeometry()
            self.invalidate()
            #this is annoying - setGeometry is called once and after a widget is added, it is not called again
            # only activate is called
            # as far as I can understand, invalidate/activate should call setGeometry in all cases
    
    #what should be reimplemented
    def setGeometry(self, rect):
        print("PartsLayout.setGeometry", rect)
        super().setGeometry(rect)
        
        """
        Pokud dobre chapu dokumentaci a zdrojaky, tak by setGeometry mela
        spocitat layout a umistit spravovane kousky.
        
        Ovsem predtim je mozne/nutne? volat super().setGeometry()
        ktera nastavni vnitrni uloziste (vyvolatelne pomoci geometry())
        ovsem hodnota muze byt upravena z rect na neco jineho podle min a max size
        
        v podstate bych tedy nemel pocitat s rect, ale s se self.geometry()
        """
        effectiveRect = self.geometry()
        width = effectiveRect.width()
        height = effectiveRect.height()
        left = effectiveRect.left()
        top = effectiveRect.top()
        
        #a ted muzu effectiveRect rozpocitat na jednotlive kousky a zavolat jejich setGeometry
                
        #calculate positions of all Nodes
        """
        try to find all paths from self.left (or self.top) as deep as possible (ideally to self.right / bottom)
        (basically it should be impossible to miss right or bottom, unless user defines his own Nodes)
        sort path based on depth
        if two paths of different depths result in different positions, higher paths get preference (not implemented)
        if two paths of the same depth results in different positions of Nodes we have a problem
        
        example (number in -()-> gives size of connection between nodes
        A -(1j)-> B -(1j)->C
        A -(1j)-> B -(2j)->C
        (left single widget size 1j, right two widgets above each other, one set to get 1j the other 2j
          all widgets the same depth)
        This is an impossible diagram (we do not know if width should be split to 2 or 3 parts). We might want to check in addPart if the newly added Part does not result in a conflict. But the would be hard if the left widget was added last.
        
        Maybe we can make a rule that desired size is the minimal, in that case the B->C size should be at least 1j or at least 2j, i.e. 2j.
        
        So if we have more than one path between two Nodes we should take the supreme. But it might be hard to calculate if imag proportions are involved. Is the proportion higher than a fixed size or not?
        
        We basically need to find all the possible paths. For each path calculate imag proportions. Then go back and compare them with respect to Node placement.
        Postupne prochaze Node od sousedu left a hledat, jestli je v cestach vice nez jedna cesta. Pokud ano, tak vybrat tu
        moznou, a pokud je nejaka, kde je vzdalenost mensi, tak je zahodit.
        Zafixovat Node a jit dal.
        """
        
        # helper functions
        def firstUnfixed(ppath):
            for it in ppath:
                if not ppath[it]["fixed"]: return it
            return None

        def calculateDistances(ppath, distance):
            s = 0
            for it in ppath:
                s += ppath[it]["size"]
                ppath[it]["distance"] = s
            
            #now calculate proportions to master distance
            #s is now the overall size for this path
            if s.imag >0:
                free = distance-s.real
                unit = free / s.imag
                for it in ppath:
                    ppath[it]["distance"] = ppath[it]["distance"].real + unit*ppath[it]["distance"].imag
                    ppath[it]["effectiveSize"] = ppath[it]["size"].real + unit*ppath[it]["size"].imag
            pass

        def calculatePaths(paths, width, endNode):
            print(paths)
            #reorganize to Node:distance maps
            ppaths = []
            for path in paths:
                if len(path)==0: continue
                
                m = {}# Node:{"fixed":False, "size":original_connection_size, "distance":calculated distance from starting edge}
                print("path", path)
                for connection in path:
                    print("connection", connection)
                    size = connection.size
                    #we can resolve proportion of master size, that will not change
                    if size.real < 1:
                        size.real *= width
                    m[connection.end] = {"fixed":False, "size":size}
                
                print(connection)
                assert connection.end is endNode, "Dangling path not connected to both edges of master widget are not allowed"
                
                calculateDistances(m, width)
                ppaths.append(m)
            
                #also reset distance of all nodes
                for it in m: it._distance = None
            
            #merge paths (enlarge effective distances as needed)
            
            finishedPaths = []
            
            while len(ppaths):
                #find first unfixed Node in each ppath
                # if it is not present in another ppath, fix it
                # if it is present in another path and it is not the first unfixed there, ignore it for now
                # if it is present in another path and it is the first unfixed in all such paths, resolve the conflict
                
                # continue till all Nodes in all paths are fixed (or paths are discarded)
                
                for m in ppaths:
                    f = firstUnfixed(m)
                    
                    if f is None:
                        #move m to finished
                        finishedPaths.append(m)
                        ppaths.remove(m)
                        break #we have modified iterator so lets end the cycle here
                    
                    #is the node present elsewhere?
                    conflicts = []
                    for o in ppaths:
                        if m is o: continue
                        if f in o: conflicts.append(o)
                    
                    if len(conflicts):
                        #are all conflicting Nodes the first unfixed on their paths?
                        if all(firstUnfixed(it) is f for it in conflicts):
                            #we can resolve the conflict (discard losing paths)
                            # select the path(s) with the largest Node distance
                            candidates = conflicts+[m]
                            sx = max(it[f]["distance"] for it in candidates)
                            for it in candidates:
                                if it[f]["distance"] < sx:
                                    #we need to make the distance to f equal to sx
                                    #questions is how to distribute the extra space between the already fixed
                                    # nodes of it? 
                                    # the simple is the since the other nodes are already fixed, let them be
                                    # and put the extra space to f
                                    
                                    # BUT it will look ugly in plots if the nodes in question are left side AxisAreas
                                    #   the one closest to DataArea will be larger that needed 
                                    
                                    # on the other hand, if the previous nodes are connected to something else
                                    #  we risk breaking an already resolved nodes :(
                                    
                                    # options are put the extra space to f
                                    #   distribute extra space equaly
                                    #   mark somehow, where the extra space should go - stretch like connection
                                    #  if we could not forget to place stretch conenction everywhere it should go
                                    #    it would be the best and fixed size could remain fixed
                                    #   but without it the conflict would remain unresolved
                                    #  we could make some stretch connection that would get the extra space preferentially
                                    #    if it is present during this situation
                                    #  stretch connection likely needs to have its end node for itself (only one in and one out
                                    #    connection) or maybe just one incomming is enough, there cannot be conflicts 
                                    #    and it can be freely moved
                                    #   if more stretches is present in conflicted region, distribute equally (or by value of the stretch)
                                    # it should be actually easy to do - if stretch is not accessible from outside, it cannot
                                    #   be used as anchor for others, and the node will remain one in one out
                                    
                                    #in first approximation, lets put the space to f
                                    it[f]["size"] = sx-it[f]["distance"]+it[f]["effectiveSize"]
                                    
                                    #recalculate distances
                                    calculateDistances(it, width)
                                else:
                                    #fix the Node
                                    it[f]["fixed"] = True
                            # here we will change ppaths, so we need to break the loop
                            break
                        else:
                            # we have to get to resolving the conflict later
                            continue 
                    else:
                        #fix this Node
                        m[f]["fixed"] = True
                
                    
            #now distances for individual nodes should match in all paths
            # all nodes should be visited
            for path in finishedPaths:
                for node in path:
                    d = path[node]["distance"]
                    if node._distance is None:
                        node._distance = d
                    else:
                        assert node._distance==d, "Conflicting paths in PartsLayout"
                assert d==width, "Path length is not correct"
            pass
            
        #resolve horizontal nodes
        calculatePaths(self.left.paths(), width, self.right)
        #resolve vertical nodes
        calculatePaths(self.top.paths(), height, self.bottom)
        
        #now all the nodes are positioned
        for widget in self._connectionsH: #assuming the all relevant pieces are present in both connectionsH and connectionsV
            conH = self._connectionsH[widget]
            conV = self._connectionsV[widget]
            widget.setGeometry(QtCore.QRectF(left+conH.start._distance, top+conV.start._distance,\
                   conH.end._distance-conH.start._distance, conV.end._distance-conV.start._distance))
        pass

    
    def sizeHint(self, which, constraint = QtCore.QSizeF()):
        print("PartsLayout.sizeHint", which, constraint)
        #we do not care about maximum size
        # we do have a minimum size probably
        # we do not have a preferred size, but we might want to assume that master widget is 600x600 px to calculate something
        #  resonable
        if which==QtCore.Qt.MaximumSize:
            return QtCore.QSizeF(-1,-1)
        return QtCore.QSizeF(600,600)
    
    def count(self):
        print("PartsLayout.count")
        return len(self._connectionsV)
    
    def itemAt(self, i):
        print("PartsLayout.itemAt", i)
        return None
    
    def removeAt(self, i):
        print("PartsLayout.itemAt", i)
        #do nothing :)
    
    def updateGeometry(self):
        print("PartsLayout.updateGeometry")
        super().updateGeometry()
        pass
    
    def activate(self):
        print("PartsLayout.activate", self.isActivated(), self.parentLayoutItem(), self.parentLayoutItem().isLayout())
        """
        super().activate will not trigger setGeometry
        if layout is already activated
        
        or layout does not have valid parentItem
            QGraphicsLayoutItem *parentItem = this;
            while (parentItem && parentItem->isLayout())
                parentItem = parentItem->parentLayoutItem();
            if (!parentItem)
                return;
            Q_ASSERT(!parentItem->isLayout());        
        
        
        how the ... am I suposed to deactivate the layout?
        """
        super().activate()
        print("after super activate")
        pass
    pass
    
"""
nova connections davame vzdy za nebo pred existujici Node
  resp pred nebo za existujici connection, protoze potrebujeme vedet, co rozdelit
  v podstate jsou 4 moznosti, pokud existuje spojeni A*B, muzeme vlozit novou (+)
     - pred A: C+A*B (puvodni spojeni zustane)
     - za A:   A+C*B (puvodni spojeni se musi zmenit)
     - pred B: A*C-B
     - za B:   A*B+C
  take je mozne vytvorit (paralelni - melo by byt na jinem "radku") dvou uz existujicich Node
     to pravdepodobne vyvola nutnost prepocitat efektivni velikosti spojeni
pokud je Node na kraji , tak musime vytovrit novou krajni

pri vypoctu je nutne zohlednit hierarchii widgetu
napr. pokud mame na radku dva ploty, ktere maji nastavenou portion of extra space 2:1
 (takze by se jejich velikost mela rozdeli 2:1)
 ale uvnitr maji AxisArea (fixed size) a DataArea (portion of extra space 1 a 1)
  (takze se snazi rozdelit extra space v pomeru 1:1 pro DataArea
  tak mame konflik zajmu
  
  pri vypoctu rozdeleni extra space musime postupovat od nejvyssich widget
  a jak jsou jednou pozice Node dany (v tomto pripade celkova velikost Plotu nemusela byt stanovena)
    tak s tim nizsi widget (vnorene) nesmi hybat
    
"""

class Figure(QtWidgets.QGraphicsWidget):
    """
    Main piece, has GraphicsAnchorLayout. Everything else is inside this.
    
    But we do not anchor directly to this layout.
    We anchor to AnchorEdge.
    
    Q: Can children be present in parents parent layout?
    """
    def __init__(self):
        super().__init__()
        self.layout = PartsLayout()
        self.layout.setContentsMargins(0,0,0,0)
        # ~ self.layout.setSpacing(0)
        self.setLayout(self.layout)
        
        if True: #testing
            p = self.palette()
            p.setColor(p.Window, QtGui.QColor("grey"))
            self.setPalette(p)        
            self.setAutoFillBackground(True)

    pass

class Part(QtWidgets.QGraphicsWidget):
    def __init__(self, *args, color="grey",  **kw):
        super().__init__(*args)
        self.layout = QtWidgets.QGraphicsAnchorLayout()
        self.setLayout(self.layout)

        # ~ self.setFlag(da.ItemClipsChildrenToShape, True)
        p = self.palette()
        self._color = color
        p.setColor(p.Window, QtGui.QColor(color))
        self.setPalette(p)        
        self.setAutoFillBackground(True)
        pass

if __name__=="__main__":
    pqrApp = QtWidgets.QApplication(sys.argv)
    
    view = QtWidgets.QGraphicsView()
    scene = QtWidgets.QGraphicsScene()
    
    view.setScene(scene)
    
    part = Figure()
    scene.addItem(part)
    part.resize(600,400)
    
    
    test = Part(color="red")
    test.setParentItem(part)
    print("pre addPart")
    part.layout.addPart(test, 200)
    print("post addPart")
    
    view.show()
    
    
    res = pqrApp.exec_()
