# -*- coding:utf-8 -*-

# ##### BEGIN LGPL LICENSE BLOCK #####
# GEOS - Geometry Engine Open Source
# http://geos.osgeo.org
#
# Copyright (C) 2011 Sandro Santilli <strk@kbt.io>
# Copyright (C) 2005 2006 Refractions Research Inc.
# Copyright (C) 2001-2002 Vivid Solutions Inc.
# Copyright (C) 1995 Olivier Devillers <Olivier.Devillers@sophia.inria.fr>
#
# This is free software you can redistribute and/or modify it under
# the terms of the GNU Lesser General Public Licence as published
# by the Free Software Foundation.
# See the COPYING file for more information.
#
# ##### END LGPL LICENSE BLOCK #####

# <pep8 compliant>

# ----------------------------------------------------------
# Partial port (version 3.7.0) by: Stephen Leger (s-leger)
#
# ----------------------------------------------------------

from .algorithms import (
    LineIntersector,
    PointLocator
    )
from .geomgraph import (
    GeometryGraphOperation,
    GeometryGraph,
    Node,
    NodeMap,
    Edge,
    Label,
    EdgeEndStar,
    EdgeEnd
    )
from .constants import (
    Position,
    Location,
    IntersectionMatrix
    )


class EdgeEndBuilder():
    """
     * Computes the geomgraph::EdgeEnd objects which arise
     * from a noded geomgraph::Edge.
    """
    def computeEdgeEnds(self, edgeEnds):
        # EdgeEnd
        l = []
        for edge in edgeEnds:
            self._computeEdgeEnds(edge, l)
        return l

    def _computeEdgeEnds(self, edge, edgeEnds):
        """
         * Creates stub edges for all the intersections in this
         * Edge (if any) and inserts them into the graph.
        """
        # EdgeIntersectionList
        eiList = edge.eiList
        # ensure that the list has entries for the first and last point of the edge
        eiList.addEndpoints()

        # no intersections, so there is nothing to do
        if eiList[0] == eiList[-1]:
            return

        i = 0
        # EdgeIntersectionList
        eiNext = eiList[i]
        i += 1
        eiPrev = None
        eiCurr = None
        while (True):
            eiPrev = eiCurr
            eiCurr = eiNext
            eiNext = None
            if eiList[i] is not eiList[-1]:
                eiNext = eiList[i]
                i += 1
            if eiCurr is not None:
                self._createEdgeEndForPrev(edge, edgeEnds, eiCurr, eiPrev)
                self._createEdgeEndForNext(edge, edgeEnds, eiCurr, eiNext)
            if eiCurr is None:
                break

    def _createEdgeEndForPrev(self, edge, l, eiCurr, eiPrev):
        """
         * Create a EdgeStub for the edge before the intersection eiCurr.
         * The previous intersection is provided
         * in case it is the endpoint for the stub edge.
         * Otherwise, the previous point from the parent edge will be the endpoint.
         *
         * eiCurr will always be an EdgeIntersection, but eiPrev may be null.
        """
        iPrev = eiCurr.segmentIndex
        if eiCurr.dist == 0.0:
            if iPrev == 0:
                return
            iPrev -= 1
        # Coordinate
        pPrev = edge._coords[iPrev]
        if eiPrev is not None and eiPrev.segmentIndex >= iPrev:
            pPrev = eiPrev._coord
        label = Label(edge.label)
        label.flip()
        ee = EdgeEnd(edge, eiCurr._coord, pPrev, label)
        l.append(ee)

    def _createEdgeEndForNext(self, edge, l, eiCurr, eiNext):
        """
         * Create a StubEdge for the edge after the intersection eiCurr.
         * The next intersection is provided
         * in case it is the endpoint for the stub edge.
         * Otherwise, the next point from the parent edge will be the endpoint.
         *
         * eiCurr will always be an EdgeIntersection, but eiNext may be null.
        """
        iNext = eiCurr.segmentIndex + 1
        # if there is no next edge there is nothing to do
        if iNext >= len(edge._coords) and eiNext is None:
            return
        # Coordinate
        pNext = edge._coords[iNext]
        # if the next intersection is in the same segment as the current, use it as the endpoint
        if eiNext is not None and eiNext.segmentIndex == eiCurr.segmentIndex:
            pNext = eiNext._coord
        ee = EdgeEnd(edge, eiCurr._coord, pNext, edge.label)
        l.append(ee)


class EdgeEndBundleStar(EdgeEndStar):
    """
     * An ordered list of EdgeEndBundle objects around a RelateNode.
     *
     * They are maintained in CCW order (starting with the positive x-axis)
     * around the node
     * for efficient lookup and topology building.
    """
    def __init__(self):
        EdgeEndStar.__init__(self)

    def insert(self, edgeEnd):
        """
         * Insert a EdgeEnd in order in the list.
         * If there is an existing EdgeStubBundle which is parallel, the EdgeEnd is
         * added to the bundle.  Otherwise, a new EdgeEndBundle is created
         * to contain the EdgeEnd.
        """
        edgeBundle = self.find(edgeEnd)
        if edgeBundle is None:
            edgeBundle = EdgeEndBundle(edgeEnd)
            self.insertEdgeEnd(edgeBundle)
        else:
            edgeBundle.insert(edgeEnd)

    def updateIM(self, im):
        """
         * Update the IM with the contribution for the EdgeStubs around the node.
        """
        for esb in self.edges:
            esb.updateIM(im)


class EdgeEndBundle(EdgeEnd):
    """
     * A collection of geomgraph::EdgeEnd objects which
     * originate at the same point and have the same direction.
    """
    def __init__(self, edgeEnd):
        EdgeEnd.__init__(self, edgeEnd.edge, edgeEnd.coord, edgeEnd.direction, edgeEnd.label)
        self._edgeEnds = []
        self.insert(edgeEnd)

    def insert(self, edgeEnd):
        self._edgeEnds.append(edgeEnd)

    def computeLabel(self, bnr):
        """
         * This computes the overall edge label for the set of
         * edges in this EdgeStubBundle.  It essentially merges
         * the ON and side labels for each edge.  These labels must be compatible
        """
        # create the label.  If any of the edges belong to areas,
        # the label must be an area label
        isArea = False

        for e in self._edgeEnds:
            if e.label.isArea():
                isArea = True
        if isArea:
            self.label = Label(Location.UNDEF, Location.UNDEF, Location.UNDEF)
        else:
            self.label = Label(Location.UNDEF)
        # compute the On label, and the side labels if present
        for i in range(2):
            self._computeLabelOn(i, bnr)
            if isArea:
                self._computeLabelSides(i)

    def updateIM(self, im):
        """
         * Update the IM with the contribution for the computed label for
         * the EdgeStubs.
        """
        Edge.updateIM(self, im)

    def _computeLabelOn(self, geomIndex, bnr):
        """
         * Compute the overall ON location for the list of EdgeStubs.
         *
         * (This is essentially equivalent to computing the self-overlay of
         * a single Geometry)
         *
         * edgeStubs can be either on the boundary (eg Polygon edge)
         * OR in the interiors (e.g. segment of a LineString)
         * of their parent Geometry.
         *
         * In addition, GeometryCollections use a algorithm::BoundaryNodeRule
         * to determine whether a segment is on the boundary or not.
         *
         * Finally, in GeometryCollections it can occur that an edge
         * is both
         * on the boundary and in the interiors (e.g. a LineString segment
         * lying on
         * top of a Polygon edge.) In this case the Boundary is
         * given precendence.
         *
         * These observations result in the following rules for computing
         * the ON location:
         *  - if there are an odd number of Bdy edges, the attribute is Bdy
         *  - if there are an even number >= 2 of Bdy edges, the attribute
         *    is Int
         *  - if there are any Int edges, the attribute is Int
         *  - otherwise, the attribute is NULL.
         *
        """
        # compute the ON location value
        boundaryCount = 0
        foundInterior = False

        for e in self._edgeEnds:
            loc = e.label.getLocation(geomIndex)
            if loc == Location.BOUNDARY:
                boundaryCount += 1
            if loc == Location.INTERIOR:
                foundInterior = True

        loc = Location.UNDEF
        if foundInterior:
            loc = Location.INTERIOR
        if boundaryCount > 0:
            loc = GeometryGraph.determineBoundary(bnr, boundaryCount)
        self.label.setLocation(geomIndex, loc)

    def _computeLabelSides(self, geomIndex):
        self._computeLabelSide(geomIndex, Position.LEFT)
        self._computeLabelSide(geomIndex, Position.RIGHT)

    def _computeLabelSide(self, geomIndex, side):
        """
         * To compute the summary label for a side, the algorithm is:
         *   FOR all edges
         *     IF any edge's location is INTERIOR for the side, side location = INTERIOR
         *     ELSE IF there is at least one EXTERIOR attribute, side location = EXTERIOR
         *     ELSE  side location = NULL
         *  <br>
         *  Note that it is possible for two sides to have apparently contradictory information
         *  i.e. one edge side may indicate that it is in the interiors of a geometry, while
         *  another edge side may indicate the exterior of the same geometry.  This is
         *  not an incompatibility - GeometryCollections may contain two Polygons that touch
         *  along an edge.  This is the reason for Interior-primacy rule above - it
         *  results in the summary label having the Geometry interiors on <b>both</b> sides.
        """
        for e in self._edgeEnds:
            if e.label.isArea():
                loc = e.label.getLocation(geomIndex, side)
                if loc == Location.INTERIOR:
                    self.label.setLocation(geomIndex, side, Location.INTERIOR)
                    return
                elif loc == Location.EXTERIOR:
                    self.label.setLocation(geomIndex, side, Location.EXTERIOR)
                    return


class RelateNode(Node):
    """
     * Represents a node in the topological graph used to compute spatial
     * relationships.
    """
    def __init__(self, coord, edges):
        Node.__init__(self, coord, edges)

    def updateIMFromEdges(self, im):
        """
         * Update the IM with the contribution for the EdgeEnds incident on this node.
        """
        # EdgeEndBundleStar
        self.edges.updateIM(im)

    def _computeIM(self, im):
        """
         * Update the IM with the contribution for this component.
         * A component only contributes if it has a labelling for both parent geometries
        """
        im.setAtLeastIfValid(self.label.getLocation(0), self.label.getLocation(1), 0)


class RelateNodeFactory():
    """
     * Used by the geomgraph::NodeMap in a RelateNodeGraph to create RelateNode objects.
    """
    def createNode(self, coords):
        return RelateNode(coords, EdgeEndBundleStar())


class RelateNodeGraph():

    """
     * Implements the simple graph of Nodes and geomgraph::EdgeEnd which is all that is
     * required to determine topological relationships between Geometries.
     *
     * Also supports building a topological graph of a single Geometry, to
     * allow verification of valid topology.
     *
     * It is <b>not</b> necessary to create a fully linked
     * PlanarGraph to determine relationships, since it is sufficient
     * to know how the Geometries interact locally around the nodes.
     * In fact, this is not even feasible, since it is not possible to compute
     * exact intersection points, and hence the topology around those nodes
     * cannot be computed robustly.
     * The only Nodes that are created are for improper intersections;
     * that is, nodes which occur at existing vertices of the Geometries.
     * Proper intersections (e.g. ones which occur between the interiors of
     * line segments)
     * have their topology determined implicitly, without creating a geomgraph::Node object
     * to represent them.
    """
    def __init__(self):
        self._nodes = NodeMap(RelateNodeFactory())

    @property
    def nodes(self):
        return self._nodes.values()

    def build(self, geomGraph):
        # compute nodes for intersections between previously noded edges
        self.computeIntersectionNodes(geomGraph, 0)

        """
         * Copy the labelling for the nodes in the parent Geometry.  These override
         * any labels determined by intersections.
        """
        self.copyNodesAndLabels(geomGraph, 0)

        """
         * Build EdgeEnds for all intersections.
        """
        eeBuilder = EdgeEndBuilder()
        # EdgeEnd
        eeList = eeBuilder.computeEdgeEnds(geomGraph.edges)

        self.insertEdgeEnds(eeList)

    def computeIntersectionNodes(self, geomGraph, argIndex):
        """
         * Insert nodes for all intersections on the edges of a Geometry.
         * Label the created nodes the same as the edge label if they do not
         * already have a label.
         * This allows nodes created by either self-intersections or
         * mutual intersections to be labelled.
         * Endpoint nodes will already be labelled from when they were inserted.
         *
         * Precondition: edge intersections have been computed.
        """
        edges = geomGraph.edges
        for edge in edges:
            eLoc = edge.label.getLocation(argIndex)
            eiL = edge.eiList
            for ei in eiL:
                # RelateNode
                n = self._nodes.addNode(ei.coord)

                if eLoc == Location.BOUNDARY:
                    n.setLabelBoundary(argIndex)
                else:
                    if n.label.isNull(argIndex):
                        n.setLabel(argIndex, Location.INTERIOR)

    def copyNodesAndLabels(self, geomGraph, argIndex):
        """
         * Copy all nodes from an arg geometry into this graph.
         * The node label in the arg geometry overrides any previously computed
         * label for that argIndex.
         * (E.g. a node may be an intersection node with
         * a computed label of BOUNDARY,
         * but in the original arg Geometry it is actually
         * in the interiors due to the Boundary Determination Rule)
        """
        nodes = geomGraph.nodes
        for node in nodes:
            newNode = self._nodes.addNode(node.coords)
            newNode.setLabel(argIndex, node.label.getLocation(argIndex))

    def insertEdgeEnds(self, edgeEnds):
        for edgeEnd in edgeEnds:
            self._nodes.addEdge(edgeEnd)


class RelateComputer():
    """
     * Computes the topological relationship between two Geometries.
     *
     * RelateComputer does not need to build a complete graph structure to compute
     * the IntersectionMatrix.  The relationship between the geometries can
     * be computed by simply examining the labelling of edges incident on each node.
     *
     * RelateComputer does not currently support arbitrary GeometryCollections.
     * This is because GeometryCollections can contain overlapping Polygons.
     * In order to correct compute relate on overlapping Polygons, they
     * would first need to be noded and merged (if not explicitly, at least
     * implicitly).
     *
    """
    def __init__(self, graphs):
        # GeometryGraph
        self.arg = graphs

        self.li = LineIntersector()
        self.ptLocator = PointLocator()

        # geomgraph::NodeMap
        self.nodes = NodeMap(RelateNodeFactory())

        # this intersection matrix will hold the results compute for the relate
        # IntersectionMatrix
        self.im = IntersectionMatrix()

        # geomgraph::Edge
        self.isolatedEdges = []
        # The intersection point found (if any)
        self.invalidPoint = None

    def computeIM(self):
        # since Geometries are finite and embedded in a 2-D space, the EE element must always be 2
        self.im.set(Location.EXTERIOR, Location.EXTERIOR, 2)
        # if the Geometries don't overlap there is nothing to do
        env0 = self.arg[0]._parentGeom.envelope
        env1 = self.arg[1]._parentGeom.envelope

        if not env0.intersects(env1):
            self.computeDisjointIM(self.im)
            return self.im

        # SegmentIntersector
        self.arg[0].computeSelfNodes(self.li, False)
        self.arg[1].computeSelfNodes(self.li, False)

        intersector = self.arg[0].computeEdgeIntersections(self.arg[1], self.li, False)

        self.computeIntersectionNodes(0)
        self.computeIntersectionNodes(1)

        """
         * Copy the labelling for the nodes in the parent Geometries.
         * These override any labels determined by intersections
         * between the geometries.
        """
        self.copyNodesAndLabels(0)
        self.copyNodesAndLabels(1)

        """
         * complete the labelling for any nodes which only have a
         * label for a single geometry
        """
        self.labelIsolatedNodes()

        """
         * If a proper intersection was found, we can set a lower bound
         * on the IM.
        """
        self.computeProperIntersectionIM(intersector, self.im)

        """
         * Now process improper intersections
         * (eg where one or other of the geometrys has a vertex at the
         * intersection point)
         * We need to compute the edge graph at all nodes to determine
         * the IM.
        """
        eeBuilder = EdgeEndBuilder()
        ee0 = eeBuilder.computeEdgeEnds(self.arg[0].edges)
        self.insertEdgeEnds(ee0)
        ee1 = eeBuilder.computeEdgeEnds(self.arg[1].edges)
        self.insertEdgeEnds(ee1)

        self.labelNodeEdges()

        """
         * Compute the labeling for isolated components.
         * Isolated components are components that do not touch any
         * other components in the graph.
         * They can be identified by the fact that they will
         * contain labels containing ONLY a single element, the one for
         * their parent geometry.
         * We only need to check components contained in the input graphs,
         * since isolated components will not have been replaced by new
         * components formed by intersections.
        """
        self.labelIsolatedEdges(0, 1)
        self.labelIsolatedEdges(1, 0)
        # update the IM from all components
        self.updateIM(self.im)
        return self.im

    def insertEdgeEnds(self, ee) -> None:
        for de in ee:
            self.nodes.add(de)

    def computeProperIntersectionIM(self, intersector, im) -> None:
        # If a proper intersection is found, we can set a lower bound on the IM.
        dimA = self.arg[0]._parentGeom.dimension
        dimB = self.arg[1]._parentGeom.dimension
        hasProper = intersector.hasProperIntersection
        hasProperInterior = intersector.hasProperInteriorIntersection
        # For Geometry's of dim 0 there can never be proper intersections.
        """
         * If edge segments of Areas properly intersect, the areas must properly overlap.
        """
        if dimA == 2 and dimB == 2:
            if hasProper:
                im.setAtLeast("212101212")

            """
             * If an Line segment properly intersects an edge segment of an Area,
             * it follows that the Interior of the Line intersects the Boundary of the Area.
             * If the intersection is a proper <i>interiors</i> intersection, then
             * there is an Interior-Interior intersection too.
             * Note that it does not follow that the Interior of the Line intersects the Exterior
             * of the Area, since there may be another Area component which contains the rest of the Line.
            """
        elif dimA == 2 and dimB == 1:
            if hasProper:
                im.setAtLeast("FFF0FFFF2")
            if hasProperInterior:
                im.setAtLeast("1FFFFF1FF")
        elif dimA == 1 and dimB == 2:
            if hasProper:
                im.setAtLeast("F0FFFFFF2")
            if hasProperInterior:
                im.setAtLeast("1F1FFFFFF")
            """
             * If edges of LineStrings properly intersect *in an interiors point*, all
             * we can deduce is that
             * the interiors intersect.  (We can NOT deduce that the exteriors intersect,
             * since some other segments in the geometries might cover the points in the
             * neighbourhood of the intersection.)
             * It is important that the point be known to be an interiors point of
             * both Geometries, since it is possible in a self-intersecting geometry to
             * have a proper intersection on one segment that is also a boundary point of another segment.
            """
        elif dimA == 1 and dimB == 1:
            if hasProperInterior:
                im.setAtLeast("0FFFFFFFF")

    def copyNodesAndLabels(self, argIndex: int) -> None:
        """
         * Copy all nodes from an arg geometry into this graph.
         * The node label in the arg geometry overrides any previously computed
         * label for that argIndex.
         * (E.g. a node may be an intersection node with
         * a computed label of BOUNDARY,
         * but in the original arg Geometry it is actually
         * in the interiors due to the Boundary Determination Rule)
        """
        nodes = self.arg[argIndex].nodes
        for node in nodes:
            newNode = self.nodes.addNode(node.coord)
            newNode.setLabel(argIndex, node.label.getLocation(argIndex))

    def computeIntersectionNodes(self, argIndex: int) -> None:
        """
         * Insert nodes for all intersections on the edges of a Geometry.
         * Label the created nodes the same as the edge label if they do not
         * already have a label.
         * This allows nodes created by either self-intersections or
         * mutual intersections to be labelled.
         * Endpoint nodes will already be labelled from when they were inserted.
        """
        edges = self.arg[argIndex].edges
        for edge in edges:
            loc = edge.label.getLocation(argIndex)
            eil = edge.eiList
            for ei in eil:
                node = self.nodes.addNode(ei.coord)
                if loc == Location.BOUNDARY:
                    node.setLabelBoundary(argIndex)
                else:
                    if node.label.isNull(argIndex):
                        node.setLabel(argIndex, Location.INTERIOR)

    def labelIntersectionNodes(self, argIndex: int) -> None:
        """
         * For all intersections on the edges of a Geometry,
         * label the corresponding node IF it doesn't already have a label.
         * This allows nodes created by either self-intersections or
         * mutual intersections to be labelled.
         * Endpoint nodes will already be labelled from when they were inserted.
        """
        edges = self.arg[argIndex].edges
        for edge in edges:
            loc = edge.label.getLocation(argIndex)
            eil = edge.eiList
            for ei in eil:
                node = self.nodes.find(ei.coord)
                if node.label.isNull(argIndex):
                    if loc == Location.BOUNDARY:
                        node.setLabelBoundary(argIndex)
                    else:
                        node.setLabel(argIndex, Location.INTERIOR)

    def computeDisjointIM(self, im) -> None:
        """
         * If the Geometries are disjoint, we need to enter their dimension and
         * boundary dimension in the Ext rows in the IM
        """
        ga = self.arg[0]._parentGeom
        if not ga.is_empty:
            im.set(Location.INTERIOR, Location.EXTERIOR, ga.dimension)
            im.set(Location.BOUNDARY, Location.EXTERIOR, ga.boundaryDimension)

        gb = self.arg[1]._parentGeom
        if not gb.is_empty:
            im.set(Location.EXTERIOR, Location.INTERIOR, gb.dimension)
            im.set(Location.EXTERIOR, Location.BOUNDARY, gb.boundaryDimension)

    def labelNodeEdges(self) -> None:
        nodes = sorted(list(self.nodes.values()), key=lambda n: (n.coord.x, n.coord.y))
        for node in nodes:
            node.edges.computeLabelling(self.arg)

    def updateIM(self, im) -> None:
        """
         * update the IM with the sum of the IMs for each component
        """
        for edge in self.isolatedEdges:
            edge.updateIM(im)
        nodes = sorted(list(self.nodes.values()), key=lambda n: (n.coord.x, n.coord.y))
        for node in nodes:
            node.updateIM(im)
            node.updateIMFromEdges(im)

    def labelIsolatedEdges(self, thisIndex: int, targetIndex: int) -> None:
        """
         * Processes isolated edges by computing their labelling and adding them
         * to the isolated edges list.
         * Isolated edges are guaranteed not to touch the boundary of the target
         * (since if they
         * did, they would have caused an intersection to be computed and hence would
         * not be isolated)
        """
        edges = self.arg[thisIndex].edges
        for edge in edges:
            if edge.isIsolated:
                self.labelIsolatedEdge(edge, targetIndex, self.arg[targetIndex]._parentGeom)
                self.isolatedEdges.append(edge)

    def labelIsolatedEdge(self, edge, targetIndex: int, target) -> None:
        """
         * Label an isolated edge of a graph with its relationship to the target
         * geometry.
         * If the target has dim 2 or 1, the edge can either be in the interiors
         * or the exterior.
         * If the target has dim 0, the edge must be in the exterior
        """
        # this won't work for GeometryCollections with both dim 2 and 1 geoms
        if target.dimension > 0:
            # since edge is not in boundary, may not need the full generality of PointLocator?
            # Possibly should use ptInArea locator instead?  We probably know here
            # that the edge does not touch the bdy of the target Geometry
            loc = self.ptLocator.locate(edge.coord, target)
            edge.label.setAllLocations(targetIndex, loc)
        else:
            edge.label.setAllLocations(targetIndex, Location.EXTERIOR)

    def labelIsolatedNodes(self) -> None:
        """
         * Isolated nodes are nodes whose labels are incomplete
         * (e.g. the location for one Geometry is null).
         * This is the case because nodes in one graph which don't intersect
         * nodes in the other are not completely labelled by the initial process
         * of adding nodes to the nodeList.
         * To complete the labelling we need to check for nodes that lie in the
         * interiors of edges, and in the interiors of areas.
        """
        nodes = self.nodes.values()
        for node in nodes:
            label = node.label
            if node.isIsolated:
                if label.isNull(0):
                    self.labelIsolatedNode(node, 0)
                else:
                    self.labelIsolatedNode(node, 1)

    def labelIsolatedNode(self, node, targetIndex: int) -> None:
        """
         * Label an isolated node with its relationship to the target geometry.
        """
        loc = self.ptLocator.locate(node.coord, self.arg[targetIndex]._parentGeom)
        node.label.setAllLocations(targetIndex, loc)


class RelateOp(GeometryGraphOperation):

    def __init__(self, g0, g1, boundaryNodeRule=None):
        GeometryGraphOperation.__init__(self, g0, g1, boundaryNodeRule)
        self.relateComp = RelateComputer([g0, g1])

    def getIntersectionMatrix(self):
        return self.relateComp.computeIM()

    @staticmethod
    def relate(a, b, boundaryNodeRule=None):
        """
         * Computes the geom::IntersectionMatrix for the spatial relationship
         * between two geom::Geometry objects, using a specified
         * Boundary Node Rule
         *
         * @param a a Geometry to test. Ownership left to caller.
         * @param b a Geometry to test. Ownership left to caller.
         * @param boundaryNodeRule the Boundary Node Rule to use.
         *
         * @return the IntersectionMatrix for the spatial relationship
         *         between the geometries. Ownership transferred.
        """
        relOp = RelateOp(a, b, boundaryNodeRule)
        return relOp.getIntersectionMatrix()
