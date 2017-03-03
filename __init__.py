# -*- coding:utf-8 -*-

#  ***** GPL LICENSE BLOCK *****
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  All rights reserved.
#  ***** GPL LICENSE BLOCK *****

bl_info = {
    'name': 'PolyLib',
    'description': 'Polygons detection from unordered splines',
    'author': 's-leger',
    'license': 'GPL',
    'deps': '',
    'version': (1, 0),
    'blender': (2, 7, 8),
    'location': 'View3D > Tools > Polygons',
    'warning': '',
    'wiki_url': '',
    'tracker_url': '',
    'link': '',
    'support': 'COMMUNITY',
    'category': '3D View'
    }
    
import os
import sys
import time
import math
import array


import shapely.geometry
import shapely.ops
import shapely.prepared
import shapely.speedups
import rtree.index
import numpy as np
from mathutils import Vector, Matrix
from mathutils.geometry import intersect_line_plane
import bpy
import bgl
import blf
from bpy_extras import view3d_utils
import bpy.utils.previews as iconsLib
from bpy.types import Operator, Panel, PropertyGroup
from bpy.props import StringProperty, FloatProperty, PointerProperty, EnumProperty, IntProperty, BoolProperty

icons_dict = {}

EPSILON = 1.0e-4

if shapely.speedups.available:
    shapely.speedups.enable()

"""
    basic bitarray for selections
"""
class BitArray():
    def __init__(self, bitSize, fill = 0):
        self.size = bitSize
        intSize = self.intSize(bitSize)
        if fill == 1:
            fill = 4294967295                                 # all bits set
        else:
            fill = 0                                      # all bits cleared
        self.bitArray = array.array('I')          # 'I' = unsigned 32-bit integer
        self.set_size(intSize, fill)
    def __str__(self):
        return str(self.list())
    def set_size(self, intSize, fill=0):
        self.bitArray.extend((fill,) * intSize)
    def intSize(self, bitSize):
        intSize = bitSize >> 5        
        if (bitSize & 31):   
            intSize += 1
        return intSize
    def test(self, bit_num):
        record = bit_num >> 5
        offset = bit_num & 31
        mask = 1 << offset
        return(self.bitArray[record] & mask)
    def set(self, bit_num):
        record = bit_num >> 5
        offset = bit_num & 31
        mask = 1 << offset
        self.bitArray[record] |= mask
    def clear(self, bit_num):
        record = bit_num >> 5
        offset = bit_num & 31
        mask = ~(1 << offset)
        self.bitArray[record] &= mask
    def toggle(self, bit_num):
        record = bit_num >> 5
        offset = bit_num & 31
        mask = 1 << offset
        self.bitArray[record] ^= mask
    def list(self):
        return [x for x in range(0, self.size) if self.test(x) > 0]
    def none(self):
        for i in range(0, len(self.bitArray)):
            self.bitArray[i] = 0
    def reverse(self):
        for i in range(0, len(self.bitArray)):
            self.bitArray[i] = 4294967295 ^ self.bitArray[i]
    def all(self):
        for i in range(0, len(self.bitArray)):
            self.bitArray[i] = 4294967295
    
"""
    test BitArray

    aa = BitArray(12)
    bb = BitArray(8)

    aa.set(0)
    aa.test(0)
    bb.set(7)
    bb.logical_or(aa)
    bb.test(0)
    bb.test(7)
    print("%s" % (bb))
"""  


# http://blenderscripting.blogspot.ch/2011/07/bgl-drawing-with-opengl-onto-blender-25.html
class GlDrawOnScreen():
    explanation_colour = (1.0, 1.0, 1.0, 0.7)
    line_colour = (1.0, 1.0, 1.0, 0.5)
    selected_colour = (1.0, 1.0, 0.0, 0.5)
    def String(self, text, x, y, size, colour):
        ''' my_string : the text we want to print
            pos_x, pos_y : coordinates in integer values
            size : font height.
            colour : used for definining the colour'''
        dpi, font_id = 72, 0 # dirty fast assignment
        bgl.glColor4f(*colour)
        blf.position(font_id, x, y, 0)
        blf.size(font_id, size, dpi)
        blf.draw(font_id, text)
    def _end(self):
        bgl.glEnd()
        bgl.glPopAttrib()
        bgl.glLineWidth(1)
        bgl.glDisable(bgl.GL_BLEND)
        bgl.glColor4f(0.0, 0.0, 0.0, 1.0)
    def _start_line(self, colour, width=2, style=bgl.GL_LINE_STIPPLE):
        bgl.glPushAttrib(bgl.GL_ENABLE_BIT)
        bgl.glLineStipple(1, 0x9999)
        bgl.glEnable(style)
        bgl.glEnable(bgl.GL_BLEND)
        bgl.glColor4f(*colour)
        bgl.glLineWidth(width)
        bgl.glBegin(bgl.GL_LINE_STRIP)
    def Line(self, x0, y0, x1, y1, colour, width=2):
        self._start_line(colour, width) 
        bgl.glVertex2i(x0, y0)
        bgl.glVertex2i(x1, y1)
        self._end()
    def Rectangle(self, x0, y0, x1, y1, colour, width=2):
        self._start_line(colour, width) 
        bgl.glVertex2i(x0, y0)
        bgl.glVertex2i(x1, y0)
        bgl.glVertex2i(x1, y1)
        bgl.glVertex2i(x0, y1)
        bgl.glVertex2i(x0, y0)
        self._end()
    def PolyLine(self, pts, colour, width=2):
        self._start_line(colour, width) 
        for pt in pts:
            x, y = pt
            bgl.glVertex2i(x, y)
        self._end()
    def Polygon(self, pts, colour, width=2, style=bgl.GL_LINE):
        self._start_line(colour, width, style) 
        #bgl.glPushAttrib(bgl.GL_ENABLE_BIT)
        #bgl.glEnable(bgl.GL_BLEND)
        #bgl.glColor4f(*colour)    
        #bgl.glBegin(bgl.GL_POLYGON)
        for pt in pts:
            x, y = pt
            bgl.glVertex2f(x, y)  
        self._end()
""" intersection of two segments outside segment (projection)
    c0 = extremite sur le segment courant
    c1 = intersection point on oposite segment
    id = oposite segment id
    t = param t on oposite segment 
    d = distance from ends to segment
    insert = do we need to insert the point on other segment
    use id, c1 and t to insert segment slices
""" 
class Prolongement():
    def __init__(self, c0, c1, id, t, d):
        self.length = c0.distance(c1)
        self.c0 = c0
        self.c1 = c1
        self.id = id
        self.t = t
        self.d = d
        
class SimplePoint():
    def __init__(self, co, precision=EPSILON):
        self.users = 0
        self.co = tuple(co)
        x,y,z = co
        self.shapeIds = []
        self.bounds = (x-precision, y-precision, x+precision, y+precision)
    """ vector from this point to another """    
    def vect(self, point):
        return np.subtract(point.co, self.co)
    """ euclidian distance between points """    
    def distance(self, point):
        return np.linalg.norm(self.vect(point))
    def add_user(self):
        self.users += 1
        
class SimpleSegment():
    def __init__(self, c0, c1, extend=EPSILON):
        self.available = True
        self.c0 = c0
        self.c1 = c1
        x0, y0, z0  = c0.co
        x1, y1, z1  = c1.co
        self.splits = []
        # this seg has an opposite
        self.opposite = False
        # source of opposite
        self.original = False
        self.bounds = (min(x0, x1)-extend, min(y0, y1)-extend, max(x0, x1)+extend, max(y0, y1)+extend)
    def vect_2d(self):
        v = self.c0.vect(self.c1)
        v[2] = 0
        return v
    def lerp(self, t):
        vect = self.c0.vect(self.c1)
        return np.add(self.c0.co, np.multiply(t, vect))
    """ _point_sur_segment 
        point: Point
        t: param t de l'intersection sur le segment courant
        d: distance laterale perpendiculaire
    """
    def _point_sur_segment(self, point):
        vect = self.c0.vect(self.c1)
        dp = point.vect(self.c0)
        dl = np.linalg.norm(vect)
        d = np.linalg.norm(np.cross(vect, dp))/dl
        t = -np.divide(np.dot(dp,vect),np.multiply(dl,dl))
        return d, t
    def point_sur_segment(self, segment):
        d, t = segment._point_sur_segment(self.c0)
        if d < EPSILON:
            if t > 0 and t < 1:
                segment.append_splits((t, self.c0))
        d, t = self._point_sur_segment(segment.c0)
        if d < EPSILON:
            if t > 0 and t < 1:
                self.append_splits((t, segment.c0))
        d, t = segment._point_sur_segment(self.c1)
        if d < EPSILON:
            if t > 0 and t < 1:
                segment.append_splits((t, self.c1))
        d, t = self._point_sur_segment(segment.c1)
        if d < EPSILON:
            if t > 0 and t < 1:
                self.append_splits((t, segment.c1))
    """ distance intersection extremite la plus proche
        t: param t de l'intersection sur le segment courant
        point: Point d'intersection
        return d: distance
    """
    def min_intersect_dist(self, t, point):
        if t > 0.5:
            return self.c1.distance(point)
        else:
            return self.c0.distance(point)        
    """ point_sur_segment return 
        p: point d'intersection
        u: param t de l'intersection sur le segment courant
        v: param t de l'intersection sur le segment segment
    """
    def intersect(self, segment):
        v2d = self.vect_2d()
        c2 = np.cross(segment.vect_2d(), (0,0,1)) 
        d  = np.dot(v2d, c2)
        if d == 0:
            # segments paralleles
            self.point_sur_segment(segment)
            return False, 0, 0, 0
        c1 = np.cross(v2d, (0,0,1))
        v3 = self.c0.vect(segment.c0)
        v3[2] = 0.0
        u = np.dot(c2,v3)/d
        v = np.dot(c1,v3)/d
        co = self.lerp(u)
        #print("pt:%s u:%f v:%f" % (co, u, v))
        return True, co, u, v
    """
    """
    def append_splits(self, split):
        if split not in self.splits:
            self.splits.append(split)
    def slice(self, d, t, point):
        if d > EPSILON:
            if t > 0.5:
                if point != self.c1:
                    self.append_splits((t, point))
            else:
                if point != self.c0:
                    self.append_splits((t,point))
    def add_user(self):
        self.c0.add_user()
        self.c1.add_user()
    def consume(self):
        self.available = False
        
class Tree():
    """
        Rtree wrapper
    """
    def __init__(self, extend=EPSILON):
        self.extend = extend
        self.idx = None
        self.geoms = []
    def init(self, extend=EPSILON):
        self.extend = extend
        self.geoms = []
        self.idx = rtree.index.Index()
    def build(self, geoms):
        t = time.time()
        self.init()
        self.geoms = geoms
        for i, geom in enumerate(geoms):
            self.idx.insert(i, geom.bounds)
        print("Tree.build() :%.2f seconds" % (time.time()-t))
    def insert(self, id, geom):
        self.idx.insert(id, geom.bounds)
    def newPoint(self, co):
        point = SimplePoint(co, self.extend)
        found = sorted(list(self.idx.intersection(point.bounds)))
        for id in found:
            return self.geoms[id]
        id = len(self.geoms)
        self.geoms.append(point)
        self.idx.insert(id, point.bounds)
        return point
    """
        allow "opposite" segments,
        those segments are not found by intersects
        and not stored in self.geoms
    """
    def newSegment(self, c0, c1):
        seg = SimpleSegment(c0, c1, self.extend)
        found  = sorted(list(self.idx.intersection(seg.bounds)))
        for id in found:
            if (self.geoms[id].c0 == c0 and self.geoms[id].c1 == c1):
                return self.geoms[id]
            if (self.geoms[id].c0 == c1 and self.geoms[id].c1 == c0):
                if not self.geoms[id].opposite:
                    self.geoms[id].opposite = seg
                    seg.original = self.geoms[id]
                return self.geoms[id].opposite
        id = len(self.geoms)
        self.geoms.append(seg)
        self.idx.insert(id, seg.bounds)
        return seg
    def intersects(self, geom):
        selection = list(self.idx.intersection(geom.bounds))
        count = len(selection)
        return count, sorted(selection)
    def delete(self, id, geom):
        self.idx.delete(id, geom.bounds)
        
"""
    Ensure uniqueness and fix precision issues by design
    
    Slice:
    1 insert new points on segments
    2 add segments and points for projections
    - dosent require tree
    3 traverse shapes and increase number of users for points and segments
    4 split shapes from end to start for points with more than 2 users
    4b consume segment so it is only used once 
    - may require a tree
    5 merge shapes ends and consume shapes so it is only used once 
    
    
    def __eq__(self, other): 
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False
    def __ne__(self, other):
        return not self.__eq__(other)
    
"""
class SimpleShape():
    """ implicit closed with last point 
        require p_tree and s_tree
    """
    def __init__(self, verts=[], closed=False):
        self.closed = closed
        self.segs = []
        self.create_segments(verts)
        self.available = True
        self.shapeId = []
    def create_segments(self, verts):
        global s_tree
        nbverts = len(verts)
        self.segs = [s_tree.newSegment(verts[v-1], verts[v]) for v in range(1, nbverts)]
    def merge(self, shape):
        if self.segs[-1].c1 == shape.segs[-1].c1 or self.segs[0].c0 == shape.segs[0].c0:
            shape.reverse()
        if self.segs[-1].c1 == shape.segs[0].c0:
            self.segs += shape.segs
        elif shape.segs[-1].c1 == self.segs[0].c0:
            self.segs = shape.segs + self.segs
        else:
            print("Shape merge failed %s %s %s %s" % (id(self), id(shape), self.shapeId, shape.shapeId))
        self.closed = bool(self.segs[0].c0 == self.segs[-1].c1)
    def nbsegs(self):
        return len(self.segs)
    def valid(self):
        return self.nbsegs() > 0
    def reverse(self):
        verts = [seg.c0 for seg in self.segs]
        verts.append(self.segs[-1].c1)
        verts = verts[::-1]
        self.create_segments(verts)
    def as_shapely(self):
        coords = [seg.c0.co for seg in self.segs]
        coords.append(self.segs[-1].c1.co)
        return shapely.geometry.LineString(coords)
    def as_spline(self, curve):
        coords = [seg.c0.co for seg in self.segs]
        if not self.closed:
            coords.append(self.segs[-1].c1.co)
        spline = curve.splines.new('POLY')
        spline.use_endpoint_u = False
        spline.use_cyclic_u   = self.closed
        spline.points.add(len(coords)-1)
        for i, coord in enumerate(coords):
            x,y,z = coord
            spline.points[i].co = (x,y,z,1)
    def from_spline(self, world_matrix, spline):
        global p_tree
        verts = []
        self.closed = spline.use_cyclic_u
        if spline.type == 'POLY':
            verts =  [p_tree.newPoint(world_matrix * p.co.to_3d()) for p in spline.points]
            if spline.use_cyclic_u:
                verts.append(verts[0])     
        if spline.type == 'BEZIER':
            verts = [p_tree.newPoint(world_matrix * p.co.to_3d()) for p in spline.bezier_points]
            if spline.use_cyclic_u:
                verts.append(verts[0])
        self.create_segments(verts)
    def from_shapely(self, geom):
        global p_tree
        verts =  [p_tree.newPoint(p) for p in list(geom.coords)]
        self.closed = geom.is_ring
        self.create_segments(verts)
    """
        slice shape into smaller parts at intersections
        if closed the remaining part at end becomes the first part
    """
    def slice(self, shapes):
        verts = []
        for seg in self.segs:
            if seg.available and not seg.original:
                seg.consume()
                verts.append(seg.c0)
                if seg.c1.users > 2:
                    verts.append(seg.c1)
                    shape = SimpleShape(verts, verts[0] == verts[-1])
                    shapes.append(shape)
                    verts = []
        if len(verts) > 0:
            verts.append(self.segs[-1].c1)
            shape = SimpleShape(verts, verts[0] == verts[-1])
            shapes.append(shape)    
    """
        add points from intersection data
    """
    def add_points(self):
        verts = []
        if len(self.segs) > 0:
            for seg in self.segs:
                verts.append(seg.c0)
                splits = sorted(seg.splits)
                for split in splits:
                    verts.append(split[1])
            verts.append(self.segs[-1].c1)
        self.create_segments(verts)
    """
        add users on segments and points
    """
    def set_users(self):
        for seg in self.segs:
            seg.add_user()
    def consume(self):
        self.available = False
        
class Io():
    """
        precision :
        scene.scale_length (factor blender unit > meter)
        1 = 1m  0.01 = 1cm  0.001 = 1mm   0.0001 = 1/10mm
    """ 
    def __init__(self, precision=4):
        self.precision = precision
    def ensure_array(self, obj):
        try:
            iterator = iter(obj)
        except TypeError:
            obj = [obj]
        return obj
    """ curve_as_line
        input curve as SimpleShape
        require p_tree and s_tree
    """
    def _curve_as_line(self, curve, lines):
        wM = curve.matrix_world
        for spline in curve.data.splines:
            line = SimpleShape()
            line.from_spline(wM, spline)
            lines.append(line)       
    def curves_as_line(self, objs, lines):
        objs = self.ensure_array(objs)
        t = time.time()
        for obj in objs:
            self._curve_as_line(obj, lines)
        print("Io.curves_as_line() :%.2f seconds" % (time.time()-t))
    """ _get_spline_co
        get vertex coords of spline as array
    """
    def _get_spline_co(self, wM, spline):
        pts = []
        if spline.type == 'POLY':
            pts =  [wM * p.co.to_3d() for p in spline.points]
            if spline.use_cyclic_u:
                pts.append(pts[0])     
        if spline.type == 'BEZIER':
            pts = [wM * p.co.to_3d() for p in spline.bezier_points]
            if spline.use_cyclic_u:
                pts.append(pts[0])
        return pts    
    """ curve_as_shapely
        input curve as LineString
    """
    def _curve_as_shapely(self, curve, lines):
        wM = curve.matrix_world
        for spline in curve.data.splines:
            pts = self._get_spline_co(wM, spline)
            line = shapely.geometry.LineString(pts)
            lines.append(line)       
    def curves_as_shapely(self, objs, lines):
        objs = self.ensure_array(objs)
        t = time.time()
        for obj in objs:
            self._curve_as_shapely(obj, lines)
        print("Io.curves_as_line() :%.2f seconds" % (time.time()-t))
    """ _geom_to_spline
        add spline to curve from geom
    """
    def _geom_to_spline(self, curve, geom, is_ring):
        co = list(geom.coords)
        len_co = len(co)
        is_ring = co[0] == co[-1]
        if is_ring:
            len_co-=1
        if len_co < 2:
            return
        if len(co[0]) < 3:
            for i in range(0, len_co):
                x,y = co[i]
                co[i] = (x, y, 0)
        spline = curve.splines.new('POLY')
        spline.use_endpoint_u = True
        spline.use_cyclic_u = is_ring
        spline.points.add(len_co-1)
        for i in range(0, len_co):
            x,y,z = co[i]
            spline.points[i].co = (x,y,z,1)
    """ _to_spline
        convert geom to spline
    """
    def _to_spline(self, curve, geom):
        if geom.geom_type == 'Point':
            return
        if geom.geom_type == 'Polygon':
            self._geom_to_spline(curve, geom.exterior, True)
            if geom.interiors is not None:
                for interior in geom.interiors:
                    self._geom_to_spline(curve, interior, True)
        elif 'Multi' in geom.geom_type:
            for g in geom.geoms:
                self._to_spline(curve,g)
        else:
            self._geom_to_spline(curve, geom, geom.is_ring)
    """ to_curve
        convert shapely geom to curve
        we may select some shape
    """ 
    def to_curve(self, geoms, name, dimensions = '3D'):
        t = time.time()
        geoms = self.ensure_array(geoms)
        scene = bpy.context.scene
        curve = bpy.data.curves.new(name, type = 'CURVE')
        curve.dimensions = dimensions
        for geom in geoms:
            self._to_spline(curve, geom)
        curve_obj = bpy.data.objects.new(name, curve) 
        scene.objects.link(curve_obj)
        print("Io.to_curve() :%.2f seconds" % (time.time()-t))
        return [curve_obj]
    """ shape_to_curve
        output shapes as curves
    """
    def shapes_to_curve(self, shapes, name, dimensions = '3D'):
        t = time.time()
        shapes = self.ensure_array(shapes)
        scene = bpy.context.scene
        curve = bpy.data.curves.new(name, type = 'CURVE')
        curve.dimensions = dimensions
        for shape in shapes:
            if shape is not None:
                shape.as_spline(curve)
        curve_obj = bpy.data.objects.new(name, curve) 
        scene.objects.link(curve_obj)
        print("Io.to_curve() :%.2f seconds" % (time.time()-t))
        return [curve_obj]
    """ to_curves
        convert shapely geoms to curves
    """ 
    def to_curves(self, geoms, name, dimensions = '3D'):
        t = time.time()
        geoms = self.ensure_array(geoms)
        scene = bpy.context.scene
        curve_objs = []
        for geom in geoms:
            curve = bpy.data.curves.new(name, type = 'CURVE')
            curve.dimensions = dimensions
            self._to_spline(curve, geom)
            curve_obj = bpy.data.objects.new(name, curve) 
            scene.objects.link(curve_obj)
            curve_objs.append(curve_obj)
        print("Io.to_curves() :%.2f seconds" % (time.time()-t))
        return curve_objs
    """
    """
    def as_shapely(self, shapes):
        return [shape.as_shapely() for shape in shapes]
    def _from_shapely(self, geom):
        shape = SimpleShape()
        shape.from_shapely(geom)
        return shape
    def from_shapely(self, geoms):
        return [self._from_shapely(geom) for geom in geoms]
    """ _geom_to_gl
        add gl from geom
    """
    def _position_2d_from_coord(self, context, coord):
        scene = context.scene
        region = context.region
        rv3d = context.region_data
        loc = view3d_utils.location_3d_to_region_2d(region, rv3d, coord)
        x, y = loc
        return [x, y]
    def _geom_to_coords(self, context, gl, geom, is_ring):
        poly = [self._position_2d_from_coord(context, co) for co in list(geom.coords)]
        if len(poly) > 1:
            gl.append(poly)
    """ _to_gl
        convert geom to gl
    """
    def geom_to_gl(self, context, gl, geom):
        if geom.geom_type == 'Point':
            return
        if geom.geom_type == 'Polygon':
            self._geom_to_coords(context, gl, geom.exterior, True)
            if geom.interiors is not None:
                for interior in geom.interiors:
                    self._geom_to_coords(context, gl, interior, True)
        elif 'Multi' in geom.geom_type:
            for g in geom.geoms:
                self.geom_to_gl(context, gl, g)
        else:
            self._geom_to_coords(context, gl, geom, geom.is_ring)
        
            
class Ops():
    def ensure_array(self, obj):
        try:
            iterator = iter(obj)
        except TypeError:
            obj = [obj]
        return obj
    """ union (shapely based)
        cascaded union
    """         
    def union(self, geoms):
        t = time.time()
        geoms = self.ensure_array(geoms)
        collection = shapely.geometry.GeometryCollection(geoms)
        union = shapely.ops.cascaded_union(collection)
        print("Ops.union() :%.2f seconds" % (time.time()-t))
        return union
    """ union2 (SimpleShape based) 
        cascaded union
        require p_tree and s_tree
    """     
    def union2(self, shapes, extend=0.001):
        split = self._split(shapes, extend=extend)
        union = self._merge(split)
        return union
    """ detect_polygons
    """ 
    def detect_polygons(self, geoms):
        print("Ops.detect_polygons()")
        t = time.time()
        result, dangles, cuts, invalids = shapely.ops.polygonize_full(geoms)
        print("Ops.detect_polygons() :%.2f seconds" % (time.time()-t))
        return result, dangles, cuts, invalids
    """ optimize 
    """ 
    def optimize(self, geoms, tolerance=0.001, preserve_topology=True):
        t = time.time()
        geoms = self.ensure_array(geoms)
        optimized = [geom.simplify(tolerance, preserve_topology) for geom in geoms]
        print("Ops.optimize() :%.2f seconds" % (time.time()-t))
        return optimized
    """ _split 
        detect intersections between segments and slice shapes according
        is able to project segment ends on closest segment
        require p_tree and s_tree
        TODO: build her own seek_tree to setup segs bounds according extend
    """
    def _intersection_point(self, d, t, point, seg):
        if d > EPSILON:
            return point
        elif t > 0.5:
            return seg.c1
        else:
            return seg.c0
    def _split(self, shapes, extend=0.01):
        global p_tree
        t = time.time()
        new_shapes = []
        segs = s_tree.geoms
        nbsegs = len(segs)
        it_start = [None for x in range(nbsegs)] 
        it_end   = [None for x in range(nbsegs)]
        for s, seg in enumerate(segs):
            count, idx = s_tree.intersects(seg)
            for id in idx:
                if id > s:
                    intersect, co, u, v = seg.intersect(segs[id])
                    if intersect:                    
                        point = p_tree.newPoint(co)
                        du = seg.min_intersect_dist(u, point)
                        dv = segs[id].min_intersect_dist(v, point)
                        # point intersection sur segment id
                        pt = self._intersection_point(dv, v, point, segs[id])
                        # print("s:%s id:%s u:%7f v:%7f du:%7f dv:%7f" % (s, id, u, v, du, dv))
                        different = du+dv > EPSILON
                        if u <= 0:
                            # prolonge segment s c0
                            if du < extend and different:
                                it = Prolongement(seg.c0, pt, id, v, du)
                                last = it_start[s]
                                if last is None or last.length > it.length:
                                    it_start[s] = it 
                        elif u < 1:
                            # intersection sur segment s
                            seg.slice(du, u, pt)
                        else:
                            # prolonge segment s c1
                            if du < extend and different:
                                it = Prolongement(seg.c1, pt, id, v, du)
                                last = it_end[s]
                                if last is None or last.length > it.length:
                                    it_end[s] = it
                        pt = self._intersection_point(du, u, point, seg)
                        if v <= 0:
                            # prolonge segment id c0
                            if dv < extend and different:
                                it = Prolongement(segs[id].c0, pt, s, u, dv)
                                last = it_start[id]
                                if last is None or last.length > it.length:
                                    it_start[id] = it 
                        elif v < 1:
                            # intersection sur segment s
                            segs[id].slice(dv, v, pt)
                        else:
                            # prolonge segment s c1
                            if dv < extend and different:
                                it = Prolongement(segs[id].c1, pt, s, u, dv)
                                last = it_end[id]
                                if last is None or last.length > it.length:
                                    it_end[id] = it
        for it in it_start:
            if it is not None:
                #print("it_start[%s] id:%s t:%4f d:%4f" % (s, it.id, it.t, it.d) )
                if it.t > 0 and it.t < 1:
                    segs[it.id].append_splits((it.t, it.c1))
                if it.d > EPSILON:
                    shape = SimpleShape([it.c0, it.c1])
                    shapes.append(shape)
        for it in it_end:            
            if it is not None:
                #print("it_end[%s] id:%s t:%4f d:%4f" % (s, it.id, it.t, it.d) )
                if it.t > 0 and it.t < 1:
                    segs[it.id].append_splits((it.t, it.c1))
                if it.d > EPSILON:
                    shape = SimpleShape([it.c0, it.c1])
                    shapes.append(shape)
        print("Ops.split() intersect :%.2f seconds" % (time.time()-t))
        t = time.time()
        for shape in shapes:
            shape.add_points()
        for shape in shapes:
            shape.set_users()
        for shape in shapes:
            shape.slice(new_shapes)
        print("Ops.split() slice :%.2f seconds" % (time.time()-t))
        return new_shapes
    """ _merge
        merge shapes ends 
        reverse use s_tree
        does not need tree as all:
        - set shape ids to end vertices
        - traverse shapes looking for verts with 2 shape ids
        - merge different shapes according
    """
    def _merge(self, shapes):
        t = time.time()
        merged = []
        nbshapes = len(shapes)
        for i, shape in enumerate(shapes): 
            shape.available = True
            shape.shapeId = [i]
            shape.segs[0].c0.shapeIds = []
            shape.segs[-1].c1.shapeIds = []
        for i, shape in enumerate(shapes):  
            shape.segs[0].c0.shapeIds.append(i)
            shape.segs[-1].c1.shapeIds.append(i)
        for i, shape in enumerate(shapes): 
            shapeIds = shape.segs[-1].c1.shapeIds
            if len(shapeIds) == 2:
                if shapeIds[0] in shape.shapeId:
                    s = shapeIds[1]
                else:
                    s = shapeIds[0]
                if shape != shapes[s]:
                    shape.merge(shapes[s])
                    shape.shapeId += shapes[s].shapeId
                    for j in shape.shapeId:
                        shapes[j] = shape
            shapeIds = shape.segs[0].c0.shapeIds
            if len(shapeIds) == 2:
                if shapeIds[0] in shape.shapeId:
                    s = shapeIds[1]
                else:
                    s = shapeIds[0]
                if shape != shapes[s]:
                    shape.merge(shapes[s])
                    shape.shapeId += shapes[s].shapeId
                    for j in shape.shapeId:
                        shapes[j] = shape
        for shape in shapes:
            if shape.available:
                shape.consume()
                merged.append(shape)
        print("Ops.merge() :%.2f seconds" % (time.time()-t))
        return merged        
    """ min_bounding_rect
        minimum area oriented bounding rect 
    """
    def min_bounding_rect(self, geom):
        # Compute edges (x2-x1,y2-y1)
        if geom.convex_hull.geom_type == 'Polygon':
            hull_points_2d = [list(coord[0:2]) for coord in list(geom.convex_hull.exterior.coords)]
        else:
            hull_points_2d = [list(coord[0:2]) for coord in list(geom.convex_hull.coords)]
        edges = np.zeros( (len(hull_points_2d)-1,2) ) # empty 2 column array
        for i in range( len(edges) ):
            edge_x = hull_points_2d[i+1][0] - hull_points_2d[i][0]
            edge_y = hull_points_2d[i+1][1] - hull_points_2d[i][1]
            edges[i] = [edge_x,edge_y]
        # Calculate edge angles   atan2(y/x)
        edge_angles = np.zeros( (len(edges)) ) # empty 1 column array
        for i in range( len(edge_angles) ):
            edge_angles[i] = math.atan2( edges[i,1], edges[i,0] )
        # Check for angles in 1st quadrant
        for i in range( len(edge_angles) ):
            edge_angles[i] = abs( edge_angles[i] % (math.pi/2) ) # want strictly positive answers
        # Remove duplicate angles
        edge_angles = np.unique(edge_angles)
        # Test each angle to find bounding box with smallest area
        min_bbox = (0, sys.maxsize, 0, 0, 0, 0, 0, 0) # rot_angle, area, width, height, min_x, max_x, min_y, max_y
        #print "Testing", len(edge_angles), "possible rotations for bounding box... \n"
        for i in range( len(edge_angles) ):
            # Create rotation matrix to shift points to baseline
            # R = [ cos(theta)      , cos(theta-PI/2)
            #       cos(theta+PI/2) , cos(theta)     ]
            R = np.array([ [ math.cos(edge_angles[i]), math.cos(edge_angles[i]-(math.pi/2)) ], [ math.cos(edge_angles[i]+(math.pi/2)), math.cos(edge_angles[i]) ] ])
            # Apply this rotation to convex hull points
            rot_points = np.dot(R, np.transpose(hull_points_2d) ) # 2x2 * 2xn
            # Find min/max x,y points
            min_x = np.nanmin(rot_points[0], axis=0)
            max_x = np.nanmax(rot_points[0], axis=0)
            min_y = np.nanmin(rot_points[1], axis=0)
            max_y = np.nanmax(rot_points[1], axis=0)
            # Calculate height/width/area of this bounding rectangle
            width = max_x - min_x
            height = max_y - min_y
            area = width*height
            # Store the smallest rect found first (a simple convex hull might have 2 answers with same area)
            if (area < min_bbox[1]):
                min_bbox = ( edge_angles[i], area, width, height, min_x, max_x, min_y, max_y )
        # Re-create rotation matrix for smallest rect
        angle = min_bbox[0]   
        R = np.array([ [ math.cos(angle), math.cos(angle-(math.pi/2)) ], [ math.cos(angle+(math.pi/2)), math.cos(angle) ] ])
        # min/max x,y points are against baseline
        min_x = min_bbox[4]
        max_x = min_bbox[5]
        min_y = min_bbox[6]
        max_y = min_bbox[7]
        # Calculate center point and project onto rotated frame
        center_x = (min_x + max_x)/2
        center_y = (min_y + max_y)/2
        center_point = np.dot( [ center_x, center_y ], R )
        if min_bbox[2] > min_bbox[3]:
            a = -math.cos(angle)
            b = math.sin(angle)
            w = min_bbox[2]
            h = min_bbox[3]
        else:
            a = -math.cos(angle+(math.pi/2))
            b = math.sin(angle+(math.pi/2))
            w = min_bbox[3]
            h = min_bbox[2]
        tM = Matrix([[a, b, 0, center_point[0]], [-b, a, 0, center_point[1]], [0,0,1,0], [0,0,0,1]])
        return tM, w, h
  
class Selectable():
    """ selectable"""
    def __init__(self):
        # shapely geometry to select from
        self.geoms = []
        # selection sets (bitArray)
        self.selections = []
        # selected objects on screen representation
        self.curves = []
        # Rtree to speedup region selections
        self.tree = None
        # BitArray ids of selected geoms 
        self.ba = None
        # Material to represent selection on screen
        self.mat = None
        self.io = Io()
        self.ops = Ops()
        self.gl = GlDrawOnScreen()
        self.action = None
        self.ready = False
    def build_display_mat(self, name, color=(1,1,0)):
        midx = bpy.data.materials.find(name)    
        if midx < 0:
            mat = bpy.data.materials.new(name)
            mat.use_object_color = True
            mat.diffuse_color = color
        else:
            mat = bpy.data.materials[midx]
        return mat
    def select_from(self, geoms):
        self.ready = True
        self.mat = self.build_display_mat("Selected", (1,1,0))
        self.geoms = geoms
        self.tree = Tree()
        self.tree.build(geoms)
        self.ba = BitArray(len(geoms))
    def _unselect(self, selection):
        t = time.time()
        for i in selection:
            self.ba.clear(i)
        print("Selectable._unselect() :%.2f seconds" % (time.time()-t))
    def _select(self, selection):
        t = time.time()
        for i in selection:
            self.ba.set(i)
        print("Selectable._select() :%.2f seconds" % (time.time()-t))
    def _position_3d_from_coord(self, context, coord):
        """Run this function on left mouse, execute the ray cast"""
        scene = context.scene
        region = context.region
        rv3d = context.region_data
        view_vector_mouse = view3d_utils.region_2d_to_vector_3d(region, rv3d, coord)
        ray_origin_mouse = view3d_utils.region_2d_to_origin_3d(region, rv3d, coord)
        loc = intersect_line_plane(ray_origin_mouse, ray_origin_mouse+view_vector_mouse, Vector((0,0,0)), Vector((0,0,1)), False)
        x, y, z = loc
        return (x, y, z)
    def _position_2d_from_coord(self, context, coord):
        scene = context.scene
        region = context.region
        rv3d = context.region_data
        loc = view3d_utils.location_3d_to_region_2d(region, rv3d, coord)
        x, y = loc
        return (int(x), int(y))
    def _contains(self, context, coord, event):
        t = time.time()
        point = self._position_3d_from_coord(context, coord)
        selection = []
        pt = shapely.geometry.Point(point)
        prepared_pt = shapely.prepared.prep(pt)
        count, gids = self.tree.intersects(pt)
        selection = [i for i in gids if prepared_pt.intersects(self.geoms[i])]
        print("Selectable._contains() :%.2f seconds" % (time.time()-t))
        if event.shift:
            self._unselect(selection)
        else:
            self._select(selection)
        self._draw(context)
    def _intersects(self, context, coord, event):
        t = time.time()
        c0 = self._position_3d_from_coord(context, coord)
        c1 = self._position_3d_from_coord(context, (coord[0], event.mouse_region_y))
        c2 = self._position_3d_from_coord(context, (event.mouse_region_x, event.mouse_region_y))
        c3 = self._position_3d_from_coord(context, (event.mouse_region_x, coord[1]))
        poly = shapely.geometry.Polygon([c0,c1,c2,c3])
        prepared_poly = shapely.prepared.prep(poly)
        count, gids = self.tree.intersects(poly)
        selection = [i for i in gids if prepared_poly.intersects(self.geoms[i])]
        print("Selectable._intersects() :%.2f seconds" % (time.time()-t))  
        if event.shift:
            self._unselect(selection)
        else:
            self._select(selection)
        self._draw(context)
    def _hide(self, context):
        t = time.time()
        if len(self.curves) > 0:
            try:
                for curve in self.curves:
                    data = curve.data
                    context.scene.objects.unlink(curve)
                    bpy.data.objects.remove(curve, do_unlink=True)                
                    if data is None:
                        return
                    name = data.name
                    if bpy.data.curves.find(name) > -1:
                        bpy.data.curves.remove(data, do_unlink=True)
            except:
                pass
            self.curves = [] 
        print("Selectable._hide() :%.2f seconds" % (time.time()-t))      
    def _draw(self, context):
        if self.ba is None:
            return
        print("Selectable._draw()")
        t = time.time()
        self._hide(context) 
        selection = [self.geoms[i] for i in self.ba.list()]
        if len(selection) > 1000:
            self.curves = self.io.to_curve(selection, 'selection', '3D')
        else:
            self.curves = self.io.to_curves(selection, 'selection', '2D')
        for curve in self.curves:
            curve.color = (1,1,0,1)
            if len(curve.data.materials) < 1:
                curve.data.materials.append(self.mat)
                curve.active_material = self.mat
            curve.select = True
        print("Selectable._draw() :%.2f seconds" % (time.time()-t))
    def store(self):
        if self.check():
            self.selections.append(self.ba)
            self.ba = BitArray(len(self.geoms))
    def recall(self):
        if self.check() and len(self.selections) > 0:
            self.ba = self.selections.pop()
    def check(self):
        return self.ready

class SelectLines(Selectable):
    """
        pick_tools actions
    """
    def _draw(self, context):
        """ override draw method """
        if self.ba is None:
            return
        print("SelectLines._draw()")
        t = time.time()
        self._hide(context) 
        selection = [self.geoms[i] for i in self.ba.list()]
        self.curves = self.io.to_curve(selection, 'selection', '3D')
        for curve in self.curves:
            curve.color = (1,1,0,1)
            curve.select = True
        print("SelectLines._draw() :%.2f seconds" % (time.time()-t))
    def init(self, pick_tool, context, action):
        # Post selection actions
        self.startPoint = (0,0)
        self.endPoint = (0,0)
        self.drag = False
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_callback, args, 'WINDOW', 'POST_PIXEL')
        self.action = action 
        self._draw(context)
        print("SelectLines.init()")
    def complete(self, context):
        print("SelectLines.complete()")
        t = time.time()
        self._hide(context)
        selection = [self.geoms[i] for i in self.ba.list()]
        if self.action == 'select' and len(selection) > 0:
            result = self.io.to_curve(selection, 'selection')
            result[0].select = True
            context.scene.objects.active = result[0]
        if self.action == 'union' and len(selection) > 0:
            shapes = self.io.from_shapely(selection)
            merged = self.ops._merge(shapes)
            union  = self.io.as_shapely(merged)
            #union  = self.ops.union(selection)
            resopt = self.ops.optimize(union)
            result = self.io.to_curve(resopt, 'union')
            result[0].select = True
            context.scene.objects.active = result[0]
        print("SelectLines.complete() :%.2f seconds" % (time.time()-t))
    def select(self, context, coord, event):
        #print("Selection.select() mouse %s %s" % (abs(event.mouse_region_x - coord[0]), abs(event.mouse_region_y - coord[1])))
        t = time.time()
        if abs(event.mouse_region_x - coord[0]) > 2 and abs(event.mouse_region_y - coord[1]) > 2:
            self._intersects(context, coord, event)
        else:
            self._contains(context, (event.mouse_region_x, event.mouse_region_y), event)
        print("SelectLines.select() :%.2f seconds" % (time.time()-t))  
    def keyboard(self, context, event):
        if event.type in {'A'}:
            if len(self.ba.list()) > 0:
                self.ba.none()
            else:    
                self.ba.all()
        elif event.type in {'I'}:
            self.ba.reverse()
        elif event.type in {'S'}:
            self.store()
        elif event.type in {'R'}:
            self.recall()
        self._draw(context)
    def modal(self, context, event):
        if event.type in {'I', 'A', 'S', 'R'} and event.value == 'PRESS':
            self.keyboard(context, event)
        elif event.type in {'RIGHTMOUSE', 'ESC'}:
            self.complete(context)
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            return {'FINISHED'}
        elif event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            self.drag = True
            self.startPoint = (event.mouse_region_x, event.mouse_region_y)
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            self.drag = False
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
            self.select(context, self.startPoint, event)
        elif event.type == 'MOUSEMOVE':
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        return {'RUNNING_MODAL'}
    def draw_callback(self, _self, context):
        self.gl.String("Select with left mouse, drag for areas, SHIFT to unselect", 10, 139, 10, self.gl.explanation_colour)
        self.gl.String("ESC or right click when done", 10, 126, 10, self.gl.explanation_colour)
        self.gl.String("A select all / none", 10, 113, 10, self.gl.explanation_colour)
        self.gl.String("I invert selection", 10, 100, 10, self.gl.explanation_colour)
        self.gl.String("S store selection", 10, 87, 10, self.gl.explanation_colour)
        self.gl.String("R retrieve selection", 10, 74, 10, self.gl.explanation_colour)
        if self.drag:
            x0, y0 = self.startPoint
            x1, y1 = self.endPoint
            self.gl.Rectangle(x0, y0, x1, y1, self.gl.line_colour)
         
class SelectPolygons(Selectable):
    def store_point(self, context, coord, event):
        self.last_point = self._position_3d_from_coord(context, (event.mouse_region_x, event.mouse_region_y))
    """
        pick_tools actions
    """
    def init(self, pick_tool, context, action):
        # Post selection actions
        self.object_location = None
        self.last_point = None
        self.selectMode = True
        self.startPoint = (0,0)
        self.endPoint = (0,0)
        self.drag = False
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_callback, args, 'WINDOW', 'POST_PIXEL')
        self.action = action 
        self._draw(context)
        print("SelectPolygons.init()")
    def complete(self, context):
        print("SelectPolygons.complete()")
        t = time.time()
        self._hide(context)
        selection = [self.geoms[i] for i in self.ba.list()]
        if self.action == 'select' and len(selection) > 0:
            result = self.io.to_curve(selection, 'selection')
            result[0].select = True
            context.scene.objects.active = result[0]
        if self.action == 'union' and len(selection) > 0:
            union  = self.ops.union(selection)
            resopt = self.ops.optimize(union)
            result = self.io.to_curve(resopt, 'union')
            result[0].select = True
            context.scene.objects.active = result[0]
        if self.action == 'object' and len(selection) > 0:
            tM, w, h  = self.object_location
            #pt = (self.last_point * tM.inverted())
            #if pt.y < 0: #reverse rotation
            #if pt.y > 0: #right
            poly = shapely.geometry.box(-w/2.0, -h/2.0, w/2.0, h/2.0, ccw=True)
            result = self.io.to_curve(poly, 'object')
            result[0].matrix_world = tM
            result[0].select = True
            self.ba.none()
            context.scene.objects.active = result[0]
        print("SelectPolygons.complete() :%.2f seconds" % (time.time()-t))
    def select(self, context, coord, event):
        #print("Selection.select() mouse %s %s" % (abs(event.mouse_region_x - coord[0]), abs(event.mouse_region_y - coord[1])))
        t = time.time()
        if abs(event.mouse_region_x - coord[0]) > 2 and abs(event.mouse_region_y - coord[1]) > 2:
            self._intersects(context, coord, event)
        else:
            self._contains(context, (event.mouse_region_x, event.mouse_region_y), event)
        print("SelectPolygons.select() :%.2f seconds" % (time.time()-t))  
    def keyboard(self, context, event):
        if event.type in {'A'}:
            if len(self.ba.list()) > 0:
                self.ba.none()
            else:    
                self.ba.all()
        elif event.type in {'I'}:
            self.ba.reverse()
        elif event.type in {'S'}:
            self.store()
        elif event.type in {'R'}:
            self.recall()
        elif event.type in {'O'}:
            self.selectMode = not self.selectMode
            sel  = [self.geoms[i] for i in self.ba.list()]
            geom = self.ops.union(sel)
            tM, w, h = self.ops.min_bounding_rect(geom)
            self.object_location = (tM, w, h)
            self.startPoint = self._position_2d_from_coord(context, tM.translation)
        self._draw(context)
    def modal(self, context, event):
        if event.type in {'I', 'A', 'S', 'R', 'O'} and event.value == 'PRESS':
            self.keyboard(context, event)
        elif event.type in {'RIGHTMOUSE', 'ESC'}:
            if self.action != 'object':
                self.complete(context)
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            return {'FINISHED'}
        elif event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            self.drag = True
            if self.selectMode:
                self.startPoint = (event.mouse_region_x, event.mouse_region_y)
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            self.drag = False
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
            if self.selectMode:
                self.select(context, self.startPoint, event)
            else:
                self.store_point(context, self.startPoint, event)
                self.complete(context)
            self.selectMode = True
        if event.type == 'MOUSEMOVE':
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        return {'RUNNING_MODAL'}
    def draw_callback(self, _self, context):
        """
        # draw selection using gl.
        # need cache
        selection = [self.geoms[i] for i in self.ba.list()]
        for geom in selection:
            gl =[]
            self.io.geom_to_gl(context, gl, geom)
            for poly in gl:
                self.gl.Polygon(poly, self.selected_colour)
        """
        if self.selectMode:
            self.gl.String("Select with left mouse, drag for areas, SHIFT to unselect", 10, 139, 10, self.gl.explanation_colour)
            self.gl.String("ESC or right click when done", 10, 126, 10, self.gl.explanation_colour)
            self.gl.String("A select all / none", 10, 113, 10, self.gl.explanation_colour)
            self.gl.String("I invert selection", 10, 100, 10, self.gl.explanation_colour)
            self.gl.String("S store selection", 10, 87, 10, self.gl.explanation_colour)
            self.gl.String("R retrieve selection", 10, 74, 10, self.gl.explanation_colour)
        else:
            self.gl.String("Pick a point on inside", 10, 139, 10, self.gl.explanation_colour)
            self.gl.String("ESC or right click to exit", 10, 126, 10, self.gl.explanation_colour)
        if self.drag:
            x0, y0 = self.startPoint
            x1, y1 = self.endPoint
            if self.selectMode:
                self.gl.Rectangle(x0, y0, x1, y1, self.gl.line_colour)
            else:
                self.gl.Line(x0, y0, x1, y1, self.gl.line_colour)
                 
""" Matrix centered on rectangular selection
    
    obj = C.object
    io = Io()
    lines = []
    io.curves_as_shapely(obj, lines)

    ops = Ops()
    tM, w, h = ops.min_bounding_rect(lines[0])

    C.object.matrix_world = tM

"""

class TOOLS_OP_PolyLib_Pick2DPolygons(Operator):
    bl_idname = "tools.poly_lib_pick_2d_polygons"
    bl_label = "Pick 2d"
    bl_description = "Pick polygons"
    bl_options = {'REGISTER', 'UNDO'}
    pass_keys = ['NUMPAD_0', 'NUMPAD_1', 'NUMPAD_3', 'NUMPAD_4',
                 'NUMPAD_5', 'NUMPAD_6', 'NUMPAD_7', 'NUMPAD_8',
                 'NUMPAD_9', 'MIDDLEMOUSE', 'WHEELUPMOUSE', 'WHEELDOWNMOUSE']
    action = StringProperty(name="action", default="select")
    @classmethod
    def poll(self, context):
        global select_polygons
        return 'select_polygons' in globals() and select_polygons.check()
    def modal(self, context, event):
        global select_polygons
        context.area.tag_redraw()
        if event.type in self.pass_keys:
            return {'PASS_THROUGH'}
        return select_polygons.modal(context, event)
    def invoke(self, context, event):
        global select_polygons
        if 'select_polygons' not in globals():
            self.report({'WARNING'}, "Use detect before")
            return {'CANCELLED'}
        elif context.space_data.type == 'VIEW_3D':
            select_polygons.init(self, context, self.action)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}
  
class TOOLS_OP_PolyLib_Pick2DLines(Operator):
    bl_idname = "tools.poly_lib_pick_2d_lines"
    bl_label = "Pick lines"
    bl_description = "Pick lines"
    bl_options = {'REGISTER', 'UNDO'}
    pass_keys = ['NUMPAD_0', 'NUMPAD_1', 'NUMPAD_3', 'NUMPAD_4',
                 'NUMPAD_5', 'NUMPAD_6', 'NUMPAD_7', 'NUMPAD_8',
                 'NUMPAD_9', 'MIDDLEMOUSE', 'WHEELUPMOUSE', 'WHEELDOWNMOUSE']
    action = StringProperty(name="action", default="select")
    @classmethod
    def poll(self, context):
        global select_lines
        return 'select_lines' in globals() and select_lines.check()
    def modal(self, context, event):
        global select_lines
        context.area.tag_redraw()
        if event.type in self.pass_keys:
            return {'PASS_THROUGH'}
        return select_lines.modal(context, event)
    def invoke(self, context, event):
        global select_lines
        if 'select_lines' not in globals():
            self.report({'WARNING'}, "Use detect before")
            return {'CANCELLED'}
        elif context.space_data.type == 'VIEW_3D':
            select_lines.init(self, context, self.action)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}
   
class TOOLS_OP_PolyLib_Detect(Operator):
    bl_idname = "tools.poly_lib_detect"
    bl_label = "Detect Polygons"
    bl_description = "Detect polygons from unordered splines"
    bl_options = {'REGISTER', 'UNDO'}
    extend = FloatProperty(name="extend", default=0.01, subtype='DISTANCE', unit='LENGTH', min=0)
    @classmethod
    def poll(self, context):
        return len(context.selected_objects) > 0 and context.object is not None and context.object.type == 'CURVE'
    def execute(self, context):
        print("Detect")
        global p_tree
        global s_tree
        global select_polygons
        global select_lines
        p_tree = Tree()
        s_tree = Tree()
        select_polygons = SelectPolygons()
        select_lines = SelectLines()
        t = time.time()
        lines = []
        objs = [obj for obj in context.selected_objects if obj.type == 'CURVE']
        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}
        for obj in objs:
            obj.select = False
        p_tree.init(0.5*EPSILON)
        s_tree.init(self.extend)
        # SimpleShape based union
        select_polygons.io.curves_as_line(objs, lines)
        union = select_polygons.ops.union2(lines, self.extend)
        geoms = select_polygons.io.as_shapely(union)
        # output select_lines here
        select_lines.select_from(geoms)
        # Shapely based union
        # select_polygons.io.curves_as_shapely(objs, lines)
        # geoms = select_polygons.ops.union(lines, self.extend)
        result, dangles, cuts, invalids = select_polygons.ops.detect_polygons(geoms)
        if len(invalids) > 0:
            errs = select_polygons.io.to_curve(invalids, "invalid_polygons")
            err_mat = select_polygons.build_display_mat("Invalid_polygon", (1,0,0))
            for curve in errs:
                # curve.data.bevel_depth = 0.02
                curve.color = (1,0,0,1)
                if len(curve.data.materials) < 1:
                    curve.data.materials.append(err_mat)
                    curve.active_material = err_mat
                curve.select = True
            self.report({'WARNING'}, str(len(invalids)) + " invalid polygons detected")
        select_polygons.select_from(result)
        print("Detect :%.2f seconds polygons:%s invalids:%s" % (time.time()-t, len(result), len(invalids)))
        return {'FINISHED'}
  
class TOOLS_OP_PolyLib_Offset(Operator):
    bl_idname = "tools.poly_lib_offset"
    bl_label = "Offset"
    bl_description = "Offset lines"
    bl_options = {'REGISTER', 'UNDO'}
    @classmethod
    def poll(self, context):
        return len(context.selected_objects) > 0 and context.object is not None and context.object.type == 'CURVE'
    def execute(self, context):
        wm = context.window_manager.poly_lib
        objs = [obj for obj in context.selected_objects if obj.type == 'CURVE']
        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}
        for obj in objs:
            obj.select = False
        io = Io()
        lines = []
        io.curves_as_shapely(objs, lines)
        offset = []
        for line in lines:
            res = line.parallel_offset(wm.offset_distance, wm.offset_side, resolution=wm.offset_resolution,  join_style=int(wm.offset_join_style), mitre_limit=wm.offset_mitre_limit) 
            offset.append(res)
        result = io.to_curve(offset, 'offset')
        result[0].select = True
        context.scene.objects.active = result[0]    
        return {'FINISHED'}
 
class TOOLS_OP_PolyLib_Simplify(Operator):
    bl_idname = "tools.poly_lib_simplify"
    bl_label = "Simplify"
    bl_description = "Simplify lines"
    bl_options = {'REGISTER', 'UNDO'}
    @classmethod
    def poll(self, context):
        return len(context.selected_objects) > 0 and context.object is not None and context.object.type == 'CURVE'
    def execute(self, context):
        wm = context.window_manager.poly_lib
        objs = [obj for obj in context.selected_objects if obj.type == 'CURVE']
        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}
        for obj in objs:
            obj.select = False
        io = Io()
        lines = []
        simple = []
        io.curves_as_shapely(objs, lines)
        for line in lines:
            res = line.simplify(wm.simplify_tolerance, preserve_topology=wm.simplify_preserve_topology) 
            simple.append(res)
        result = io.to_curve(simple, 'simplify')
        result[0].select = True
        context.scene.objects.active = result[0]    
        return {'FINISHED'}    
       
class TOOLS_OP_PolyLib_OutputPolygons(Operator):
    bl_idname = "tools.poly_lib_output_polygons"
    bl_label = "Output Polygons"
    bl_description = "Output all polygons"
    bl_options = {'REGISTER', 'UNDO'}
    @classmethod
    def poll(self, context):
        global select_polygons
        return 'select_polygons' in globals() and select_polygons.check()
    def execute(self, context):
        global select_polygons
        result = select_polygons.io.to_curve(select_polygons.geoms, 'polygons')
        result[0].select = True
        context.scene.objects.active = result[0]
        return {'FINISHED'}
 
class TOOLS_OP_PolyLib_OutputLines(Operator):
    bl_idname = "tools.poly_lib_output_lines"
    bl_label = "Output lines"
    bl_description = "Output all lines"
    bl_options = {'REGISTER', 'UNDO'}
    @classmethod
    def poll(self, context):
        global select_lines
        return 'select_lines' in globals() and select_lines.check()
    def execute(self, context):
        global select_lines
        result = select_lines.io.to_curve(select_lines.geoms, 'lines')
        result[0].select = True
        context.scene.objects.active = result[0]
        return {'FINISHED'}
 
class TOOLS_OP_PolyLib_Solidify(Operator):
    bl_idname = "tools.poly_lib_solidify"
    bl_label = "Extrude"
    bl_description = "Extrude all polygons"
    bl_options = {'REGISTER', 'UNDO'}
    @classmethod
    def poll(self, context):
        return len(context.selected_objects) > 0 and context.object is not None and context.object.type == 'CURVE'
    def execute(self, context):
        wm = context.window_manager.poly_lib
        objs = [obj for obj in context.selected_objects if obj.type == 'CURVE']
        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}
        for obj in objs:
            obj.data.dimensions = '2D'
            mod = obj.modifiers.new("Solidify", 'SOLIDIFY')    
            mod.thickness = wm.solidify_thickness
            mod.offset = 1.00
            mod.use_even_offset = True
            mod.use_quality_normals = True
        return {'FINISHED'}
 
class TOOLS_PT_PolyLib(Panel):
    bl_label = "Polygons"
    bl_idname = "TOOLS_PT_PolyLib"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Tools"
    def draw(self, context):
        layout = self.layout
        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        row.operator("tools.poly_lib_detect", icon_value=icons_dict["detect"].icon_id, text ='Detect').extend = context.window_manager.poly_lib.extend
        row.prop(context.window_manager.poly_lib, "extend")
        row = box.row(align=True)
        row.label(text="Polygons")
        row = box.row(align=True)
        row.operator("tools.poly_lib_pick_2d_polygons", icon_value=icons_dict["selection"].icon_id, text ='Select').action = 'select'
        row.operator("tools.poly_lib_pick_2d_polygons", icon_value=icons_dict["union"].icon_id, text ='Union').action = 'union'
        row.operator("tools.poly_lib_output_polygons", icon_value=icons_dict["polygons"].icon_id, text ='All')
        row.operator("tools.poly_lib_pick_2d_polygons", text ='Obj').action = 'object'
        row = box.row(align=True)
        row.label(text="Lines")
        row = box.row(align=True)
        row.operator("tools.poly_lib_pick_2d_lines", icon_value=icons_dict["selection"].icon_id, text ='Lines').action = 'select'
        row.operator("tools.poly_lib_pick_2d_lines", icon_value=icons_dict["union"].icon_id, text ='Union').action = 'union'
        row.operator("tools.poly_lib_output_lines", icon_value=icons_dict["polygons"].icon_id, text ='All')
        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        row.operator("tools.poly_lib_solidify")
        row.prop(context.window_manager.poly_lib, "solidify_thickness")
        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        row.operator("tools.poly_lib_simplify")
        row.prop(context.window_manager.poly_lib, "simplify_tolerance")
        row = box.row(align=True)
        row.prop(context.window_manager.poly_lib, "simplify_preserve_topology")
        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        row.operator("tools.poly_lib_offset")
        row = box.row(align=True)
        row.prop(context.window_manager.poly_lib, "offset_distance")
        row = box.row(align=True)
        row.prop(context.window_manager.poly_lib, "offset_side")
        row = box.row(align=True)
        row.prop(context.window_manager.poly_lib, "offset_resolution")
        row = box.row(align=True)
        row.prop(context.window_manager.poly_lib, "offset_join_style")
        row = box.row(align=True)
        row.prop(context.window_manager.poly_lib, "offset_mitre_limit")
        
class PolyLibParameters(PropertyGroup):
    bl_idname = 'tools.poly_lib_parameters'
    extend = FloatProperty(name="Extend", description="Extend to closest intersecting segment", default=0.01, subtype='DISTANCE', unit='LENGTH', min=0)
    offset_distance = FloatProperty(name="Distance", default=0.05, subtype='DISTANCE', unit='LENGTH', min=0)
    offset_side = EnumProperty(name="Side", default='left', items = [('left', 'Left', 'Left'),('right', 'Right', 'Right')])
    offset_resolution = IntProperty(name="Resolution", default=16)
    offset_join_style = EnumProperty(name="Style", default='2', items = [('1', 'Round', 'Round'),('2', 'Mitre', 'Mitre'),('3', 'Bevel', 'Bevel')])
    offset_mitre_limit = FloatProperty(name="Mitre limit", default=10.0, subtype='DISTANCE', unit='LENGTH', min=0)
    simplify_tolerance = FloatProperty(name="Tolerance", default=0.01, subtype='DISTANCE', unit='LENGTH', min=0)
    simplify_preserve_topology = BoolProperty(name = "Preserve topology", description="Preserve topology (fast without, but may introduce self crossing)", default = True)
    solidify_thickness = FloatProperty(name="Thickness", default=2.7, subtype='DISTANCE', unit='LENGTH', min=0)

def register():
    global icons_dict
    icons_dict = iconsLib.new()
    icons_dir = os.path.join(os.path.dirname(__file__), "icons")
    for icon in os.listdir(icons_dir):
        name, ext = os.path.splitext(icon)
        icons_dict.load(name, os.path.join(icons_dir, icon), 'IMAGE')
    bpy.utils.register_class(TOOLS_OP_PolyLib_Pick2DPolygons)
    bpy.utils.register_class(TOOLS_OP_PolyLib_Pick2DLines)
    bpy.utils.register_class(TOOLS_OP_PolyLib_OutputPolygons)
    bpy.utils.register_class(TOOLS_OP_PolyLib_OutputLines)
    bpy.utils.register_class(TOOLS_OP_PolyLib_Offset)
    bpy.utils.register_class(TOOLS_OP_PolyLib_Simplify)
    bpy.utils.register_class(TOOLS_OP_PolyLib_Detect)
    bpy.utils.register_class(TOOLS_OP_PolyLib_Solidify)
    bpy.utils.register_class(PolyLibParameters)
    bpy.types.WindowManager.poly_lib = PointerProperty(type = PolyLibParameters)
    bpy.utils.register_class(TOOLS_PT_PolyLib)

def unregister():
    global icons_dict
    iconsLib.remove(icons_dict)
    bpy.utils.unregister_class(TOOLS_OP_PolyLib_Pick2DPolygons)
    bpy.utils.unregister_class(TOOLS_OP_PolyLib_Pick2DLines)
    bpy.utils.unregister_class(TOOLS_OP_PolyLib_Detect)
    bpy.utils.unregister_class(TOOLS_OP_PolyLib_OutputPolygons)
    bpy.utils.unregister_class(TOOLS_OP_PolyLib_OutputLines)
    bpy.utils.unregister_class(TOOLS_OP_PolyLib_Offset)
    bpy.utils.unregister_class(TOOLS_OP_PolyLib_Simplify)
    bpy.utils.unregister_class(TOOLS_OP_PolyLib_Solidify)
    bpy.utils.unregister_class(TOOLS_PT_PolyLib)
    bpy.utils.unregister_class(PolyLibParameters)
    del bpy.types.WindowManager.poly_lib
 
if __name__ == "__main__":
    register()
