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


from .constants import (
    logger,
    Envelope,
    GeometryCombiner,
    PolygonExtracter
    )
from .index_strtree import (
    STRtree,
    ItemsListItem
    )


class GeometryListHolder(list):
    """
     * Helper class holding Geometries, part of which are held by reference
     * others are held exclusively.
    """
    def __init__(self):
        list.__init__(self)


class CascadedUnion():
    """
     * Provides an efficient method of unioning a collection of Geometries
     *
     * This algorithm is more robust than the simple iterated approach
     * of repeatedly unioning each geometry to a result geometry.
    """

    """
     * The effectiveness of the index is somewhat sensitive
     * to the node capacity.
     * Testing indicates that a smaller capacity is better.
     * For an STRtree, 4 is probably a good number (since
     * this produces 2x2 "squares").
    """
    STRTREE_NODE_CAPACITY = 4

    def __init__(self, geoms=None):
        """
         * Creates a new instance to union
         * the given collection of {@link Geometry}s.
         *
         * @param geoms a collection of {@link Geometryal} {@link Geometry}s
         *        ownership of elements _and_ vector are left to caller.
        """
        self.geoms = geoms
        self._factory = None

    @staticmethod
    def union(polys):
        logger.debug("cascadedUnion")
        return CascadedUnion(polys)._union()

    def _union(self, geoms=None, start=None, end=None):
        """
         * Computes the union of a collection of {@link Geometry}s.
         *
         * @param geoms a collection of {@link Geometry}s.
         *        ownership of elements _and_ vector are left to caller.
         * @return the union of the input geometries
         * @return null if no input geometries were provided
        """
        if geoms is None:
            geoms = self.geoms

        if geoms is None:
            return None

        if start is not None and end is not None:
            geoms = [geoms[i] for i in range(start, end)]

        self._factory = geoms[0]._factory

        index = STRtree(CascadedUnion.STRTREE_NODE_CAPACITY)

        for geom in geoms:
            index.insert(geom.envelope, geom)

        itemTree = index.itemsTree()
        return self.unionTree(itemTree)

    def unionTree(self, geomTree):
        """
         * Recursively unions all subtrees in the list into single geometries.
         * The result is a list of Geometry's only
        """
        geoms = self.reduceToGeometries(geomTree)
        return self.binaryUnion(geoms)

    def binaryUnion(self, geoms):
        """
         * Unions a list of geometries
         * by treating the list as a flattened binary tree,
         * and performing a cascaded union on the tree.
         * @param geom GeometryListHolder
        """
        return self._binaryUnion(geoms, 0, len(geoms))

    def _binaryUnion(self, geoms, start, end):
        """
         * Unions a section of a list using a recursive binary union on each half
         * of the section.
         *
         * @param geoms GeometryListHolder
         * @param start
         * @param end
         * @return the union of the list section
        """
        if end - start < 2:
            return self.unionSafe(geoms[start], None)
        elif end - start == 2:
            return self.unionSafe(geoms[start], geoms[start + 1])
        else:
            mid = int((end + start) / 2)
            geom0 = self._binaryUnion(geoms, start, mid)
            geom1 = self._binaryUnion(geoms, mid, end)
            return self.unionSafe(geom0, geom1)

    def reduceToGeometries(self, geomTree):
        """
         * Reduces a tree of geometries to a list of geometries
         * by recursively unioning the subtrees in the list.
         *
         * @param geomTree a tree-structured list of geometries
         * @return a list of Geometrys
        """
        logger.debug("CascadedUnion.reduceToGeometries:\n %s", geomTree)
        geoms = []
        for item in geomTree:
            if item.t == ItemsListItem.item_is_list:
                geom = self.unionTree(item.l)
                geoms.append(geom)

            elif item.t == ItemsListItem.item_is_geometry:
                geoms.append(item.g)

            else:
                assert(0), "should never be reached"
        return geoms

    def unionSafe(self, geom0, geom1):
        """
         * Computes the union of two geometries,
         * either of both of which may be null.
         *
         * @param g0 a Geometry
         * @param g1 a Geometry
         * @return the union of the input(s)
         * @return null if both inputs are null
        """
        if geom0 is None and geom1 is None:
            return None
        if geom0 is None:
            return geom1.clone()
        if geom1 is None:
            return geom0.clone()

        return self.unionOptimized(geom0, geom1)

    def unionOptimized(self, geom0, geom1):
        env0 = geom0.envelope
        env1 = geom1.envelope

        if not env0.intersects(env1):
            return GeometryCombiner.combine(geom0, geom1)

        if geom0.numgeoms < 2 and geom1.numgeoms < 2:
            return self.unionActual(geom0, geom1)

        commonEnv = Envelope()
        env0.intersection(env1, commonEnv)

        return self.unionUsingEnvelopeIntersection(geom0, geom1, commonEnv)

    def unionUsingEnvelopeIntersection(self, geom0, geom1, env):
        """
         * Unions two geometries.
         * The case of multi geometries is optimized to union only
         * the components which lie in the intersection of the two geometry's
         * envelopes.
         * Geometrys outside this region can simply be combined with the union
         * result, which is potentially much faster.
         * This case is likely to occur often during cascaded union, and may also
         * occur in real world data (such as unioning data for parcels on
         * different street blocks).
         *
         * @param g0 a geometry
         * @param g1 a geometry
         * @param common the intersection of the envelopes of the inputs
         * @return the union of the inputs
        """
        disjointGeoms = []
        g0Int = self.extractByEnvelope(env, geom0, disjointGeoms)
        g1Int = self.extractByEnvelope(env, geom1, disjointGeoms)
        u = self.unionActual(g0Int, g1Int)
        disjointGeoms.append(u)

        return GeometryCombiner.combine(disjointGeoms)
    
    def _extractByEnvelope(self, env, geom, intersectingGeoms: list, disjointGeoms: list) -> None:
        for i in range(geom.ngeoms):
            g = geom.getGeometryN(i)
            if g.envelope.intersects(env):
                intersectingGeoms.append(g)
            else:
                disjointGeoms.append(g)
    
    def extractByEnvelope(self, env, geom, disjointGeoms: list):
        intersectingGeoms = []
        self._extractByEnvelope(env, geom, intersectingGeoms, disjointGeoms)
        return self._factory.buildGeometry(intersectingGeoms)

    def unionActual(self, geom0, geom1):
        """
         * Encapsulates the actual unioning of two polygonal geometries.
         *
         * @param g0
         * @param g1
         * @return
        """
        return geom0.union(geom1)

        
class CascadedPolygonUnion(CascadedUnion):
    """
     * Provides an efficient method of unioning a collection of
     * {@link Polygonal} geometries.
     * This algorithm is faster and likely more robust than
     * the simple iterated approach of
     * repeatedly unioning each polygon to a result geometry.
     *
     * The <tt>buffer(0)</tt> trick is sometimes faster, but can be less robust and
     * can sometimes take an exceptionally long time to complete.
     * This is particularly the case where there is a high degree of overlap
     * between the polygons.  In this case, <tt>buffer(0)</tt> is forced to compute
     * with <i>all</i> line segments from the outset,
     * whereas cascading can eliminate many segments
     * at each stage of processing.
     * The best case for buffer(0) is the trivial case
     * where there is <i>no</i> overlap between the input geometries.
     * However, this case is likely rare in practice.
    """
    
    """
     * The effectiveness of the index is somewhat sensitive
     * to the node capacity.
     * Testing indicates that a smaller capacity is better.
     * For an STRtree, 4 is probably a good number (since
     * this produces 2x2 "squares").
    """
    STRTREE_NODE_CAPACITY = 4
    
    def __init__(self, geoms):
        """
         * Creates a new instance to union
         * the given collection of {@link Geometry}s.
         *
         * @param geoms a collection of {@link Polygonal} {@link Geometry}s
         *        ownership of elements _and_ vector are left to caller.
        """
        CascadedUnion.__init__(self)
        
    def restrictToPolygons(self, geom):
        """
         * Computes a {@link Geometry} containing only {@link Polygonal} components.
         *
         * Extracts the {@link Polygon}s from the input
         * and returns them as an appropriate {@link Polygonal} geometry.
         *
         * If the input is already <tt>Polygonal</tt>, it is returned unchanged.
         *
         * A particular use case is to filter out non-polygonal components
         * returned from an overlay operation.
         *
         * @param g the geometry to filter
         * @return a Polygonal geometry
        """
        if geom.geometryTypeId == GeometryTypeId.GEOS_POLYGON:
            return geom
            
        polys = []
        PolygonExtracter.getPolygons(geom, polys)
        if len(polys) == 1:
            return polys[0].clone()
        
        newPolys = [poly.clone() for poly in polys]
        return self._factory.createMultiPolygon(newPolys)
        
    @staticmethod
    def union(geoms):
        """
         * Computes the union of
         * a collection of {@link Polygonal} {@link Geometry}s.
         *
         * @param polys a collection of {@link Polygonal} {@link Geometry}s.
        """
        polys = None
        try:
            iter(geoms)
            polys = []
            for poly in geoms:
                PolygonExtracter.getPolygons(geoms, polys)
        except TypeError:
            pass
            
        try:
            polys = []
            PolygonExtracter.getPolygons(geoms, polys)
        except:
            pass
        
        return CascadedPolygonUnion(polys)._union()
        
    def unionUsingEnvelopeIntersection(self, g0, g1, common):
        """
         * Unions two polygonal geometries, restricting computation
         * to the envelope intersection where possible.
         *
         * The case of MultiPolygons is optimized to union only
         * the polygons which lie in the intersection of the two geometry's
         * envelopes.
         * Polygons outside this region can simply be combined with the union
         * result, which is potentially much faster.
         * This case is likely to occur often during cascaded union, and may also
         * occur in real world data (such as unioning data for parcels on
         * different street blocks).
         *
         * @param g0 a polygonal geometry
         * @param g1 a polygonal geometry
         * @param common the intersection of the envelopes of the inputs
         * @return the union of the inputs
        """
        disjointPolys = []
        g0Int = self.extractByEnvelope(env, g0, disjointPolys)
        g1Int = self.extractByEnvelope(env, g1, disjointPolys)
        
        u = self.unionActual(g0Int, g1Int)
        
        if len(disjointPolys) == 0:
            return u
            
        polysOn = []
        polysOff = []
        self._extractGeomListByEnvelope(u.envelope, disjointPolys, polysOn, polysOff)
        if len(polysOn) == 0:
            disjointPolys.append(u)
            ret = GeometryCombiner.combine(disjointPolys)
        else:
            ret = GeometryCombiner.combine(disjointPolys)
            ret = self.unionActual(ret, u)
            
        return ret
    
    def _extractGeomListByEnvelope(self, env, geomList: list, intersectingPolys: list, disjointPolys: list) -> None:
        for g in geomList:
            if g.envelope.intersects(env):
                intersectingPolys.append(g)
            else:
                disjointPolys.append(g)
    
    def unionActual(self, g0, g1):
        """
         * Encapsulates the actual unioning of two polygonal geometries.
         *
         * @param g0
         * @param g1
         * @return
        """
        return self.restrictToPolygons(g0.union(g1))