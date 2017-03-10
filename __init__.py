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
    'version': (1, 1),
    'blender': (2, 7, 8),
    'location': 'View3D > Tools > Polygons',
    'warning': '',
    'wiki_url': 'https://github.com/s-leger/blenderPolygons/wiki',
    'tracker_url': 'https://github.com/s-leger/blenderPolygons/issues',
    'link': 'https://github.com/s-leger/blenderPolygons',
    'support': 'COMMUNITY',
    'category': '3D View'
    }
    
import os
import sys
import time
import bpy
import bgl
import blf
import numpy as np
from math import cos, sin, pi, atan2

# let shapely import throw ImportError when missing
from shapely.geometry import (
    LineString, Polygon, GeometryCollection, box
)
import shapely.ops
import shapely.prepared
import shapely.speedups

from .bitarray import BitArray
from .pyqtree import _QuadTree
from mathutils import Vector, Matrix
from mathutils.geometry import intersect_line_plane
from bpy_extras import view3d_utils
from bpy.types import Operator, Panel, PropertyGroup
from bpy.props import StringProperty, FloatProperty, PointerProperty, EnumProperty, IntProperty, BoolProperty
from bpy.utils import previews as iconsLib

if shapely.speedups.available:
    shapely.speedups.enable()
    

icons_dict = {}

# precision 1e-4 = 0.1mm
EPSILON = 1.0e-4
# Qtree params
MAX_ITEMS = 10
MAX_DEPTH = 20

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

class Prolongement():
    """ intersection of two segments outside segment (projection)
        c0 = extremite sur le segment courant
        c1 = intersection point on oposite segment
        id = oposite segment id
        t = param t on oposite segment 
        d = distance from ends to segment
        insert = do we need to insert the point on other segment
        use id, c1 and t to insert segment slices
    """ 
    def __init__(self, c0, c1, id, t, d):
        self.length = c0.distance(c1)
        self.c0 = c0
        self.c1 = c1
        self.id = id
        self.t = t
        self.d = d
        
class Point():
    
    def __init__(self, co, precision=EPSILON):
        self.users = 0
        self.co = tuple(co)
        x,y,z = co
        self.shapeIds = []
        self.bounds = (x-precision, y-precision, x+precision, y+precision)
        
    def vect(self, point):
        """ vector from this point to another """    
        return np.subtract(point.co, self.co)
    
    def distance(self, point):
        """ euclidian distance between points """    
        return np.linalg.norm(self.vect(point))
    
    def add_user(self):
        self.users += 1
       
class Segment():
    
    def __init__(self, c0, c1, extend=EPSILON):
        
        self.c0 = c0
        self.c1 = c1
        self._splits = []
        
        self.available = True
        # ensure uniqueness when merge 
        
        self.opposite = False
        # this seg has an opposite
        
        self.original = False
        # source of opposite
        
        x0, y0, z0  = c0.co
        x1, y1, z1  = c1.co
        self.bounds = (min(x0, x1)-extend, min(y0, y1)-extend, max(x0, x1)+extend, max(y0, y1)+extend)
    
    @property
    def splits(self):
        return sorted(self._splits)
    
    @property
    def vect(self):
        """ vector c0-c1"""
        return np.subtract(self.c1.co, self.c0.co)
        
    @property
    def vect_2d(self):
        v = self.vect
        v[2] = 0
        return v
    
    def lerp(self, t):
        return np.add(self.c0.co, np.multiply(t, self.vect))
    
    def _point_sur_segment(self, point):
        """ _point_sur_segment 
            point: Point
            t: param t de l'intersection sur le segment courant
            d: distance laterale perpendiculaire
        """
        vect = self.vect
        dp = point.vect(self.c0)
        dl = np.linalg.norm(vect)
        d = np.linalg.norm(np.cross(vect, dp))/dl
        t = -np.divide(np.dot(dp,vect),np.multiply(dl,dl))
        if d < EPSILON:
            if t > 0 and t < 1:
                self._append_splits((t, point))
        
    def is_end(self, point):
        return point == self.c0 or point == self.c1
          
    def min_intersect_dist(self, t, point):
        """ distance intersection extremite la plus proche
            t: param t de l'intersection sur le segment courant
            point: Point d'intersection
            return d: distance
        """
        if t > 0.5:
            return self.c1.distance(point)
        else:
            return self.c0.distance(point)  
            
    def intersect(self, segment):
        """ point_sur_segment return 
            p: point d'intersection
            u: param t de l'intersection sur le segment courant
            v: param t de l'intersection sur le segment segment
        """
        v2d = self.vect_2d
        c2 = np.cross(segment.vect_2d, (0,0,1)) 
        d  = np.dot(v2d, c2)
        if d == 0:
            # segments paralleles
            segment._point_sur_segment(self.c0)
            segment._point_sur_segment(self.c1)
            self._point_sur_segment(segment.c0)
            self._point_sur_segment(segment.c1)
            return False, 0, 0, 0
        c1 = np.cross(v2d, (0,0,1))
        v3 = self.c0.vect(segment.c0)
        v3[2] = 0.0
        u = np.dot(c2,v3)/d
        v = np.dot(c1,v3)/d
        co = self.lerp(u)
        #print("pt:%s u:%f v:%f" % (co, u, v))
        return True, co, u, v
        
    def _append_splits(self, split):
        """
            append a unique split point 
        """
        if split not in self._splits:
            self._splits.append(split)
            
    def slice(self, d, t, point):
        if d > EPSILON:
            if t > 0.5:
                if point != self.c1:
                    self._append_splits((t, point))
            else:
                if point != self.c0:
                    self._append_splits((t,point))
                    
    def add_user(self):
        self.c0.add_user()
        self.c1.add_user()
        
    def consume(self):
        self.available = False

class Shape():
    """    
        Ensure uniqueness and fix precision issues by design
        implicit closed with last point 
        require p_tree and s_tree
    """
    
    def __init__(self, verts=[]):
        """
            @vertex: list of coords
        """
        self.available = True
        # Ensure uniqueness of shape when merging
    
        self._segs = []
        # Shape segments
    
        self.shapeId = []
        # Id of shape in shapes to keep a track of shape parts when merging
    
        self._create_segments(verts)
        
    def _create_segments(self, verts):
        global s_tree
        if 's_tree' not in globals():
            raise RuntimeError('Shape _create_segments require a global s_tree spacial index ')
        self._segs = list(s_tree.newSegment(verts[v], verts[v+1]) for v in range(len(verts)-1))
    
    @property
    def coords(self):
        coords = list(seg.c0.co for seg in self._segs)
        coords.append(self.c1.co)
        return coords
        
    @property
    def points(self):
        points = list(seg.c0 for seg in self._segs)
        points.append(self.c1)
        return points
    
    @property
    def c0(self):
        if not self.valid:
            raise RuntimeError('Shape does not contains any segments')
        return self._segs[0].c0
    
    @property
    def c1(self):
        if not self.valid:
            raise RuntimeError('Shape does not contains any segments')
        return self._segs[-1].c1
    
    @property
    def nbsegs(self):
        return len(self._segs)
    
    @property
    def valid(self):
        return self.nbsegs > 0
        
    @property
    def closed(self):
        return self.valid and bool(self.c0 == self.c1)
    
    def merge(self, shape):
        """ merge this shape with specified shape
            shapes must share at least one vertex
        """
        if not self.valid or not shape.valid:
            raise RuntimeError('Trying to merge invalid shape')
        if self.c1 == shape.c1 or self.c0 == shape.c0:
            shape._reverse()
        if self.c1 == shape.c0:
            self._segs += shape._segs
        elif shape.c1 == self.c0:
            self._segs = shape._segs + self._segs
        else:
            # should never happen
            raise RuntimeError("Shape merge failed {} {} {} {}".format(id(self), id(shape), self.shapeId, shape.shapeId))
    
    def _reverse(self):
        """
            reverse vertex order
        """
        verts = self.points[::-1]
        self._create_segments(verts)
                       
    def slice(self, shapes):
        """
            slice shape into smaller parts at intersections
        """
        if not self.valid:
            raise RuntimeError('Cant slice invalid shape')
        verts = []
        for seg in self._segs:
            if seg.available and not seg.original:
                seg.consume()
                verts.append(seg.c0)
                if seg.c1.users > 2:
                    verts.append(seg.c1)
                    shape = Shape(verts)
                    shapes.append(shape)
                    verts = []
        if len(verts) > 0:
            verts.append(self.c1)
            shape = Shape(verts)
            shapes.append(shape)    
           
    def add_points(self):
        """
            add points from intersection data
        """
        verts = []
        if self.nbsegs > 0:
            for seg in self._segs:
                verts.append(seg.c0)
                for split in seg.splits:
                    verts.append(split[1])
            verts.append(self.c1)
        self._create_segments(verts)
        
    def set_users(self):
        """
            add users on segments and points
        """
        for seg in self._segs:
            seg.add_user()
                   
    def consume(self):
        self.available = False
         
class Qtree(_QuadTree):
    """
    The top spatial index to be created by the user. Once created it can be
    populated with geographically placed members that can later be tested for
    intersection with a user inputted geographic bounding box. Note that the
    index can be iterated through in a for-statement, which loops through all
    all the quad instances and lets you access their properties.
    """
    def __init__(self, objs, extend=EPSILON, max_items=MAX_ITEMS, max_depth=MAX_DEPTH):
        """
            objs may be blender objects or shapely geoms
            extend: how mutch seek arround
        """
        self._extend = extend
        self._geoms = []
        x = []
        y = []
        if len(objs) > 0:
            if hasattr(objs[0], 'bound_box'):
                for obj in objs:
                    x.append(obj.bound_box[0][0])
                    x.append(obj.bound_box[6][0])
                    y.append(obj.bound_box[0][1])
                    y.append(obj.bound_box[6][1]) 
            elif hasattr(objs[0], 'bounds'):
                for geom in objs:
                    x0, y0, x1, y1 = geom.bounds
                    x.append(x0)
                    x.append(x1)
                    y.append(y0)
                    y.append(y1)
        else:
            raise Exception("Qtree require at least one object to initialize bounds")
        x0 = min(x)
        x1 = max(x)
        y0 = min(y)
        y1 = max(y)
        width, height = x1-x0, y1-y0
        midx, midy = x0+width/2.0, y0+height/2.0
        super(Qtree, self).__init__(midx, midy, width, height, max_items, max_depth)
        
    @property
    def nbgeoms(self):
        return len(self._geoms)
     
    def build(self, geoms):
        """
            Build a spacial index from shapely geoms
        """
        t = time.time()
        self._geoms = geoms
        for i, geom in enumerate(geoms):
            self._insert(i, geom.bounds)
        print("Qtree.build() :%.2f seconds" % (time.time()-t))
        
    def insert(self, id, geom):
        self._geoms.append(geom)
        self._insert(id, geom.bounds)
    
    def newPoint(self, co):
        point = Point(co, self._extend)
        count, found  = self.intersects(point)
        for id in found:
            return self._geoms[id]
        self.insert(self.nbgeoms, point)
        return point
    
    def newSegment(self, c0, c1):
        """
            allow "opposite" segments,
            those segments are not found by intersects
            and not stored in self.geoms
        """
        new_seg = Segment(c0, c1, self._extend)
        count, found  = self.intersects(new_seg)
        for id in found:
            old_seg = self._geoms[id]
            if (old_seg.c0 == c0 and old_seg.c1 == c1):
                return old_seg
            if (old_seg.c0 == c1 and old_seg.c1 == c0):
                if not old_seg.opposite:
                    old_seg.opposite = new_seg
                    new_seg.original = old_seg
                return old_seg.opposite
        self.insert(self.nbgeoms, new_seg)
        return new_seg
        
    def intersects(self, geom):
        selection = list(self._intersect(geom.bounds))
        count = len(selection)
        return count, sorted(selection)
       
class Io():
    
    @staticmethod
    def ensure_iterable(obj):
        try:
            iterator = iter(obj)
        except TypeError:
            obj = [obj]
        return obj
    
    # Conversion methods
    @staticmethod
    def _to_geom(shape):
        if not shape.valid:
            raise RuntimeError('Cant convert invalid shape to Shapely LineString')
        return shapely.geometry.LineString(shape.coords)
        
    @staticmethod
    def shapes_to_geoms(shapes):
        return [Io._to_geom(shape) for shape in shapes]
    
    @staticmethod
    def _to_shape(geom):
        global p_tree
        if "p_tree" not in globals():
            raise RuntimeError("geoms to shapes require a global p_tree spacial index")
        verts =  list(p_tree.newPoint(p) for p in list(geom.coords))
        shape = Shape(verts)
        return shape
    
    @staticmethod    
    def geoms_to_shapes(geoms):
        return [Io._to_shape(geom) for geom in geoms]

    # Input methods
    @staticmethod
    def _coords_from_spline(wM, spline):
        pts = []
        if spline.type == 'POLY':
            pts = [wM * p.co.to_3d() for p in spline.points]
            if spline.use_cyclic_u:
                pts.append(pts[0])     
        if spline.type == 'BEZIER':
            pts = [wM * p.co.to_3d() for p in spline.bezier_points]
            if spline.use_cyclic_u:
                pts.append(pts[0])
        return pts
    
    @staticmethod
    def _add_geom_from_curve(curve, geoms):
        wM = curve.matrix_world
        for spline in curve.data.splines:
            pts = Io._coords_from_spline(wM, spline)
            geom = shapely.geometry.LineString(pts)
            geoms.append(geom)
    
    @staticmethod
    def curves_to_geoms(curves, geoms=[]):
        """ 
            @curves : blender curves collection 
            Return shapely geometry 
        """
        curves = Io.ensure_iterable(curves)
        t = time.time()
        for curve in curves:
            Io._add_geom_from_curve(curve, geoms)
        print("Io.curves_as_line() :%.2f seconds" % (time.time()-t))
        return geoms            
    
    @staticmethod
    def _unique_coords_from_spline(wM, spline):
        global p_tree
        if "p_tree" not in globals():
            raise RuntimeError("shapes from curves require a global p_tree spacial index")
        pts = []
        if spline.type == 'POLY':
            pts =  [p_tree.newPoint(wM * p.co.to_3d()) for p in spline.points]
        if spline.type == 'BEZIER':
            pts = [p_tree.newPoint(wM * p.co.to_3d()) for p in spline.bezier_points]
        if spline.use_cyclic_u:
            pts.append(pts[0])
        return pts
        
    @staticmethod
    def _add_shape_from_curve(curve, shapes):
        wM = curve.matrix_world
        for spline in curve.data.splines:
            pts = Io._unique_coords_from_spline(wM, spline)
            shape = Shape(verts=pts)
            shapes.append(shape)
    
    @staticmethod
    def curves_to_shapes(curves, shapes=[]):
        """ 
            @curves : blender curves collection 
            Return simple shapes 
        """
        curves = Io.ensure_iterable(curves)
        t = time.time()
        for curve in curves:
            Io._add_shape_from_curve(curve, shapes)
        print("Io.curves_as_line() :%.2f seconds" % (time.time()-t))
        return shapes
        
    # Output methods
    @staticmethod
    def _add_spline(curve, geometry):
        coords = list(geometry.coords)
        spline = curve.splines.new('POLY')
        spline.use_endpoint_u = False
        spline.use_cyclic_u  = coords[0] == coords[-1]
        spline.points.add(len(coords)-1)
        for i, coord in enumerate(coords):
            x,y,z = coord
            spline.points[i].co = (x,y,z,1)
    
    @staticmethod
    def _as_spline(curve, geometry):
        """
            add a spline into a blender curve
            @curve : blender curve
        """
        if hasattr(geometry, 'exterior'):
            # Polygon
            Io._add_spline(curve, geometry.exterior)
            for geom in geometry.interiors:
                Io._add_spline(curve, geom)
        elif hasattr(geometry, 'geoms'):
            # Multi and Collections
            for geom in geometry.geoms:
                Io._as_spline(curve, geom)
        else:
            # LinearRing, LineString and Shape
            Io._add_spline(curve, geometry)
  
    @staticmethod
    def to_curve(scene, geoms, name, dimensions = '3D'):
        t = time.time()
        geoms = Io.ensure_iterable(geoms)
        curve = bpy.data.curves.new(name, type = 'CURVE')
        curve.dimensions = dimensions
        for geom in geoms:
            Io._as_spline(curve, geom)
        curve_obj = bpy.data.objects.new(name, curve) 
        scene.objects.link(curve_obj)
        curve_obj.select = True
        print("Io.to_curves() :%.2f seconds" % (time.time()-t))
        return curve_obj
    
    @staticmethod
    def to_curves(scene, geoms, name, dimensions = '3D'):
        geoms = Io.ensure_iterable(geoms)
        return [Io.to_curve(scene, geom, name, dimensions) for geom in geoms]
         
class ShapelyOps():
     
    @staticmethod
    def min_bounding_rect(geom):
        """ min_bounding_rect
            minimum area oriented bounding rect 
        """
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
        
    @staticmethod
    def detect_polygons(geoms):
        """ detect_polygons
        """ 
        print("Ops.detect_polygons()")
        t = time.time()
        result, dangles, cuts, invalids = shapely.ops.polygonize_full(geoms)
        print("Ops.detect_polygons() :%.2f seconds" % (time.time()-t))
        return result, dangles, cuts, invalids
        
    @staticmethod
    def optimize(geoms, tolerance=0.001, preserve_topology=True):
        """ optimize 
        """ 
        t = time.time()
        geoms = Io.ensure_iterable(geoms)
        optimized = [geom.simplify(tolerance, preserve_topology) for geom in geoms]
        print("Ops.optimize() :%.2f seconds" % (time.time()-t))
        return optimized
        
    @staticmethod
    def union(geoms):
        """ union (shapely based)
            cascaded union - may require snap before use to fix precision issues
            use union2 for best performances
        """         
        t = time.time()
        geoms = Io.ensure_iterable(geoms)
        collection = shapely.geometry.GeometryCollection(geoms)
        union = shapely.ops.cascaded_union(collection)
        print("Ops.union() :%.2f seconds" % (time.time()-t))
        return union
  
class ShapeOps():
    
    @staticmethod
    def union(shapes, extend=0.001):
        """ union2 (Shape based) 
            cascaded union
            require p_tree and s_tree
        """        
        split = ShapeOps.split(shapes, extend=extend)
        union = ShapeOps.merge(split)
        return union
        
    @staticmethod    
    def _intersection_point(d, t, point, seg):
        if d > EPSILON:
            return point
        elif t > 0.5:
            return seg.c1
        else:
            return seg.c0
    
    @staticmethod        
    def split(shapes, extend=0.01):
        """ _split 
            detect intersections between segments and slice shapes according
            is able to project segment ends on closest segment
            require p_tree and s_tree
        """
        global p_tree
        global s_tree
        t = time.time()
        new_shapes = []
        segs = s_tree._geoms
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
                        pt = ShapeOps._intersection_point(dv, v, point, segs[id])
                        # print("s:%s id:%s u:%7f v:%7f du:%7f dv:%7f" % (s, id, u, v, du, dv))
                        if u <= 0:
                            # prolonge segment s c0
                            if du < extend and not seg.is_end(pt):
                                it = Prolongement(seg.c0, pt, id, v, du)
                                last = it_start[s]
                                if last is None or last.length > it.length:
                                    it_start[s] = it 
                        elif u < 1:
                            # intersection sur segment s
                            seg.slice(du, u, pt)
                        else:
                            # prolonge segment s c1
                            if du < extend and not seg.is_end(pt):
                                it = Prolongement(seg.c1, pt, id, v, du)
                                last = it_end[s]
                                if last is None or last.length > it.length:
                                    it_end[s] = it
                        pt = ShapeOps._intersection_point(du, u, point, seg)
                        if v <= 0:
                            # prolonge segment id c0
                            if dv < extend and not segs[id].is_end(pt):
                                it = Prolongement(segs[id].c0, pt, s, u, dv)
                                last = it_start[id]
                                if last is None or last.length > it.length:
                                    it_start[id] = it 
                        elif v < 1:
                            # intersection sur segment s
                            segs[id].slice(dv, v, pt)
                        else:
                            # prolonge segment s c1
                            if dv < extend and not segs[id].is_end(pt):
                                it = Prolongement(segs[id].c1, pt, s, u, dv)
                                last = it_end[id]
                                if last is None or last.length > it.length:
                                    it_end[id] = it
        for it in it_start:
            if it is not None:
                #print("it_start[%s] id:%s t:%4f d:%4f" % (s, it.id, it.t, it.d) )
                if it.t > 0 and it.t < 1:
                    segs[it.id]._append_splits((it.t, it.c1))
                if it.d > EPSILON:
                    shape = Shape([it.c0, it.c1])
                    shapes.append(shape)
        for it in it_end:            
            if it is not None:
                #print("it_end[%s] id:%s t:%4f d:%4f" % (s, it.id, it.t, it.d) )
                if it.t > 0 and it.t < 1:
                    segs[it.id]._append_splits((it.t, it.c1))
                if it.d > EPSILON:
                    shape = Shape([it.c0, it.c1])
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
    
    @staticmethod
    def merge(shapes):
        """ merge
            merge shapes ends 
            reverse use s_tree
            does not need tree as all:
            - set shape ids to end vertices
            - traverse shapes looking for verts with 2 shape ids
            - merge different shapes according
        """
        t = time.time()
        merged = []
        nbshapes = len(shapes)
        for i, shape in enumerate(shapes): 
            shape.available = True
            shape.shapeId = [i]
            shape.c0.shapeIds = []
            shape.c1.shapeIds = []
        for i, shape in enumerate(shapes):  
            shape.c0.shapeIds.append(i)
            shape.c1.shapeIds.append(i)
        for i, shape in enumerate(shapes): 
            shapeIds = shape.c1.shapeIds
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
            shapeIds = shape.c0.shapeIds
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
 
class Selectable():
    """ selectable"""
    def __init__(self, geoms):
        # shapely geometry to select from
        self.geoms = geoms
        # selection sets (bitArray)
        self.selections = []
        # selected objects on screen representation
        self.curves = []
        # Rtree to speedup region selections
        self.tree = Qtree(geoms)
        self.tree.build(geoms)
        # BitArray ids of selected geoms 
        self.ba = BitArray(len(geoms))
        # Material to represent selection on screen
        self.mat = self.build_display_mat("Selected", (1,1,0))
        self.gl = GlDrawOnScreen()
        self.action = None
        self.ready = True
        
    def build_display_mat(self, name, color=(1,1,0)):
        midx = bpy.data.materials.find(name)    
        if midx < 0:
            mat = bpy.data.materials.new(name)
            mat.use_object_color = True
            mat.diffuse_color = color
        else:
            mat = bpy.data.materials[midx]
        return mat
    
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
        selection = [self.geoms[i] for i in self.ba.list]
        if len(selection) > 1000:
            self.curves = [Io.to_curve(context.scene, selection, 'selection', '3D')]
        else:
            self.curves = Io.to_curves(context.scene, selection, 'selection', '2D')
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
        selection = list(self.geoms[i] for i in self.ba.list)
        self.curves = [Io.to_curve(context.scene, selection, 'selection', '3D')]
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
        scene = context.scene
        selection = list(self.geoms[i] for i in self.ba.list)
        if self.action == 'select' and len(selection) > 0:
            result = Io.to_curve(scene, selection, 'selection')
            scene.objects.active = result
        if self.action == 'union' and len(selection) > 0:
            shapes = Io.geoms_to_shapes(selection)
            merged = ShapeOps.merge(shapes)
            union  = Io.shapes_to_geoms(merged)
            #union  = self.ops.union(selection)
            resopt = ShapelyOps.optimize(union)
            result = Io.to_curve(scene, resopt, 'union')
            scene.objects.active = result
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
            if len(self.ba.list) > 0:
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
        scene = context.scene
        self._hide(context)
        selection = list(self.geoms[i] for i in self.ba.list)
        if self.action == 'select' and len(selection) > 0:
            result = Io.to_curve(scene, selection, 'selection')
            scene.objects.active = result
        if self.action == 'union' and len(selection) > 0:
            union  = ShapelyOps.union(selection)
            resopt = ShapelyOps.optimize(union)
            result = Io.to_curve(scene, resopt, 'union')
            scene.objects.active = result
        if self.action == 'object' and len(selection) > 0:
            tM, w, h  = self.object_location
            #pt = (self.last_point * tM.inverted())
            #if pt.y < 0: #reverse rotation
            #if pt.y > 0: #right
            poly = shapely.geometry.box(-w/2.0, -h/2.0, w/2.0, h/2.0, ccw=True)
            result = Io.to_curve(scene, poly, 'object')
            result.matrix_world = tM
            self.ba.none()
            scene.objects.active = result
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
            if len(self.ba.list) > 0:
                self.ba.none()
            else:    
                self.ba.all()
        elif event.type in {'I'}:
            self.ba.reverse()
        elif event.type in {'S'}:
            self.store()
        elif event.type in {'R'}:
            self.recall()
        elif event.type in {'B'}:
            areas = [self.geoms[i].area for i in self.ba.list]
            area = max(areas)
            self.ba.none()
            for i, geom in enumerate(self.geoms):
                if geom.area > area:
                    self.ba.set(i)
        elif event.type in {'O'}:
            self.selectMode = not self.selectMode
            sel  = [self.geoms[i] for i in self.ba.list]
            geom = ShapelyOps.union(sel)
            tM, w, h = ShapelyOps.min_bounding_rect(geom)
            self.object_location = (tM, w, h)
            self.startPoint = self._position_2d_from_coord(context, tM.translation)
        self._draw(context)
        
    def modal(self, context, event):
        if event.type in {'I', 'A', 'S', 'R', 'O', 'B'} and event.value == 'PRESS':
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
        selection = [self.geoms[i] for i in self.ba.list]
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
            self.gl.String("B select bigger area than current", 10, 87, 10, self.gl.explanation_colour)
            self.gl.String("S store selection", 10, 74, 10, self.gl.explanation_colour)
            self.gl.String("R retrieve selection", 10, 61, 10, self.gl.explanation_colour)
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
        t = time.time()
        shapes = []
        objs = [obj for obj in context.selected_objects if obj.type == 'CURVE']
        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}
        for obj in objs:
            obj.select = False
        p_tree = Qtree(objs, extend=0.5*EPSILON)
        s_tree = Qtree(objs, extend=self.extend)
        # Shape based union
        Io.curves_to_shapes(objs, shapes)
        union = ShapeOps.union(shapes, self.extend)
        geoms = Io.shapes_to_geoms(union)
        # output select_lines here
        select_lines = SelectLines(geoms)
        # Shapely based union
        # select_polygons.io.curves_as_shapely(objs, lines)
        # geoms = select_polygons.ops.union(lines, self.extend)
        result, dangles, cuts, invalids = ShapelyOps.detect_polygons(geoms)
        select_polygons = SelectPolygons(result)
        if len(invalids) > 0:
            errs = Io.to_curve(context.scene, invalids, "invalid_polygons")
            err_mat = select_polygons.build_display_mat("Invalid_polygon", (1,0,0))
            for curve in errs:
                # curve.data.bevel_depth = 0.02
                curve.color = (1,0,0,1)
                if len(curve.data.materials) < 1:
                    curve.data.materials.append(err_mat)
                    curve.active_material = err_mat
                curve.select = True
            self.report({'WARNING'}, str(len(invalids)) + " invalid polygons detected")
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
        objs = list(obj for obj in context.selected_objects if obj.type == 'CURVE')
        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}
        for obj in objs:
            obj.select = False
        lines = []
        Io.curves_to_geoms(objs, lines)
        offset = []
        for line in lines:
            res = line.parallel_offset(wm.offset_distance, wm.offset_side, resolution=wm.offset_resolution,  join_style=int(wm.offset_join_style), mitre_limit=wm.offset_mitre_limit) 
            offset.append(res)
        result = Io.to_curve(context.scene, offset, 'offset')
        context.scene.objects.active = result    
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
        lines = []
        simple = []
        Io.curves_to_geoms(objs, lines)
        for line in lines:
            res = line.simplify(wm.simplify_tolerance, preserve_topology=wm.simplify_preserve_topology) 
            simple.append(res)
        result = Io.to_curve(context.scene, simple, 'simplify')
        context.scene.objects.active = result    
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
        result = Io.to_curve(context.scene, select_polygons, 'polygons')
        context.scene.objects.active = result
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
        result = Io.to_curve(context.scene, select_lines, 'lines')
        context.scene.objects.active = result
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
