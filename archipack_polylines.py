# -*- coding:utf-8 -*-

# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

# ----------------------------------------------------------
# Author: Stephen Leger (s-leger)
#
# ----------------------------------------------------------

import time
import bpy
import bgl
import numpy as np
from math import cos, sin, pi, atan2
import bmesh
from mathutils import Vector, Matrix
from mathutils.geometry import intersect_line_plane, interpolate_bezier
from bpy_extras import view3d_utils
from bpy.types import Operator, PropertyGroup
from bpy.props import (
    StringProperty,
    FloatProperty,
    PointerProperty,
    EnumProperty,
    IntProperty,
    BoolProperty
    )
from bpy.app.handlers import persistent
from .bitarray import BitArray
from .pyqtree import _QuadTree
from .materialutils import MaterialUtils
from .archipack_gl import (
    FeedbackPanel,
    GlCursorFence,
    GlCursorArea,
    GlLine,
    GlPolyline
    )
from .pygeos.op_polygonize import PolygonizeOp
from .pygeos.geom import GeometryFactory
from .pygeos.geom import Point as GeosPoint
from .pygeos.shared import (
    Coordinate,
    CoordinateSequence,
    Envelope,
    TopologyException
    )
from .pygeos.prepared import PreparedGeometryFactory
from .pygeos.op_polygonsunion import PolygonsUnionOp

import logging
logger = logging.getLogger("archipack")


# module shared
# precision 1e-4 = 0.1mm
EPSILON = 1.0e-4
# Qtree params
MAX_ITEMS = 10
MAX_DEPTH = 20

# module globals vars dict
vars_dict = {
    # keep track of shapely geometry selection sets
    'select_polygons': None,
    'select_lines': None,
    'select_points': None
    }


class Selectable(object):

    """ selectable shapely geoms """
    def __init__(self, geoms, coordsys):
        # selection sets (bitArray)
        self.selections = []
        # selected objects on screen representation
        self.curves = []
        # Rtree to speedup region selections
        self.tree = Qtree(coordsys)
        self.tree.build(geoms)
        # BitArray ids of selected geoms
        self.ba = BitArray(self.ngeoms)
        # Material to represent selection on screen
        self.mat = self.build_display_mat("Selected",
                color=bpy.context.user_preferences.themes[0].view_3d.object_selected)
        self.cursor_fence = GlCursorFence()
        self.cursor_fence.enable()
        self.cursor_area = GlCursorArea()
        self.feedback = FeedbackPanel()
        self.action = None
        self.store_index = 0
        self.gf = GeometryFactory()

    @property
    def coordsys(self):
        return self.tree.coordsys

    @property
    def geoms(self):
        return self.tree._geoms

    @property
    def ngeoms(self):
        return self.tree.ngeoms

    @property
    def nsets(self):
        return len(self.selections)

    def build_display_mat(self, name, color=(0.2, 0.2, 0)):
        mat = MaterialUtils.build_default_mat(name, color)
        mat.use_object_color = True
        mat.emit = 0.2
        mat.alpha = 0.2
        mat.game_settings.alpha_blend = 'ADD'
        return mat

    def _unselect(self, selection):
        t = time.time()
        for i in selection:
            self.ba.clear(i)
        print("Selectable._unselect() :%.2f seconds" % (time.time() - t))

    def _select(self, selection):
        t = time.time()
        for i in selection:
            self.ba.set(i)
        print("Selectable._select() :%.2f seconds" % (time.time() - t))

    def _position_3d_from_coord(self, context, coord):
        """return point in local input coordsys
        """
        region = context.region
        rv3d = context.region_data
        view_vector_mouse = view3d_utils.region_2d_to_vector_3d(region, rv3d, coord)
        ray_origin_mouse = view3d_utils.region_2d_to_origin_3d(region, rv3d, coord)
        loc = intersect_line_plane(ray_origin_mouse, ray_origin_mouse + view_vector_mouse,
                                   Vector((0, 0, 0)), Vector((0, 0, 1)), False)
        x, y, z = self.coordsys.invert * loc
        return Vector((x, y, z))

    def _position_2d_from_coord(self, context, coord):
        """ coord given in local input coordsys
        """
        region = context.region
        rv3d = context.region_data
        loc = view3d_utils.location_3d_to_region_2d(region, rv3d, self.coordsys.world * coord)
        x, y = loc
        return Vector((x, y))

    def _contains(self, context, coord, event):
        t = time.time()
        point = self._position_3d_from_coord(context, coord)
        selection = []
        pt = self.gf.createPoint(point)
        prepared_pt = PreparedGeometryFactory.prepare(pt)
        count, gids = self.tree.intersects(pt)
        selection = [i for i in gids if prepared_pt.intersects(self.geoms[i])]
        print("Selectable._contains() :%.2f seconds selected:%s" % (time.time() - t, len(selection)))
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
        cs = [self.gf.createCoordinate(pt) for pt in [c0, c1, c2, c3, c0]]
        ring = self.gf.createLinearRing(cs)
        if not ring.is_ccw:
            ring = self.gf.createLinearRing(list(reversed(cs)))
        poly = self.gf.createPolygon(ring)
        prepared_poly = PreparedGeometryFactory.prepare(poly)
        count, gids = self.tree.intersects(poly)

        if event.ctrl:
            selection = [i for i in gids if prepared_poly.contains(self.geoms[i])]
        else:
            selection = [i for i in gids if prepared_poly.intersects(self.geoms[i])]
        print("Selectable._intersects() :%.2f seconds selected:%s" % (time.time() - t, len(selection)))
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
                    if bpy.data.curves.find(name) > - 1:
                        bpy.data.curves.remove(data, do_unlink=True)
            except:
                pass
            self.curves = []
        print("Selectable._hide() :%.2f seconds" % (time.time() - t))

    def _draw(self, context):
        print("Selectable._draw() %s" % (self.coordsys.world))
        t = time.time()
        self._hide(context)
        selection = [self.geoms[i] for i in self.ba.list]
        if len(selection) > 1000:
            self.curves = [Io.to_curve(context.scene, self.coordsys, selection, 'selection', '3D')]
        else:
            self.curves = Io.to_curves(context.scene, self.coordsys, selection, 'selection', '2D')
        for curve in self.curves:
            curve.color = (1, 1, 0, 1)
            if len(curve.data.materials) < 1:
                curve.data.materials.append(self.mat)
                curve.active_material = self.mat
            curve.select = True
        print("Selectable._draw() :%.2f seconds" % (time.time() - t))

    def store(self):
        self.selections.append(self.ba.copy)
        self.store_index = self.nsets

    def recall(self):
        if self.nsets > 0:
            if self.store_index < 1:
                self.store_index = self.nsets
            self.store_index -= 1
            self.ba = self.selections[self.store_index].copy

    def select(self, context, coord, event):
        if abs(event.mouse_region_x - coord[0]) > 2 and abs(event.mouse_region_y - coord[1]) > 2:
            self._intersects(context, coord, event)
        else:
            self._contains(context, (event.mouse_region_x, event.mouse_region_y), event)

    def init(self, pick_tool, context, action):
        raise NotImplementedError("Selectable must implement init(self, pick_tool, context, action)")

    def keyboard(self, context, event):
        """ keyboard events modal handler """
        raise NotImplementedError("Selectable must implement keyboard(self, context, event)")

    def complete(self, context):
        raise NotImplementedError("Selectable must implement complete(self, context)")

    def modal(self, context, event):
        """ modal handler """
        raise NotImplementedError("Selectable must implement modal(self, context, event)")

    def draw_callback(self, _self, context):
        """ a gl draw callback """
        raise NotImplementedError("Selectable must implement draw_callback(self, _self, context)")


class SelectPoints(Selectable):

    def __init__(self, geoms, coordsys):
        super(SelectPoints, self).__init__(geoms, coordsys)

    def _draw(self, context):
        """ override draw method """
        print("SelectPoints._draw()")
        t = time.time()
        self._hide(context)
        selection = list(self.geoms[i] for i in self.ba.list)
        if len(selection) < 2:
            return
        gf = GeometryFactory()
        geom = gf.buildGeometry(selection)
        # geom = ShapelyOps.union(selection)
        # if geom is None:
        #    return
        self.curves = [Io.to_curve(context.scene, self.coordsys, geom.convex_hull, 'selection', '3D')]
        for curve in self.curves:
            curve.color = (1, 1, 0, 1)
            curve.select = True
        print("SelectPoints._draw() :%.2f seconds" % (time.time() - t))

    def init(self, pick_tool, context, action):
        # Post selection actions
        self.selectMode = True
        self.object_location = None
        self.startPoint = (0, 0)
        self.endPoint = (0, 0)
        self.drag = False
        self.feedback.instructions(context, "Select Points", "Click & Drag to select points in area", [
            ('SHIFT', 'deselect'),
            ('CTRL', 'contains'),
            ('A', 'All'),
            ('I', 'Inverse'),
            ('F', 'Create line around selection'),
            # ('W', 'Create window using selection'),
            # ('D', 'Create door using selection'),
            ('ALT+F', 'Create best fit rectangle'),
            ('R', 'Retrieve selection'),
            ('S', 'Store selection'),
            ('ESC or RIGHTMOUSE', 'exit when done')
        ])
        self.feedback.enable()
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_callback, args, 'WINDOW', 'POST_PIXEL')
        self.action = action
        self._draw(context)
        print("SelectPoints.init()")

    def complete(self, context):
        self.feedback.disable()
        self._hide(context)

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
        elif event.type in {'F'}:
            sel = [self.geoms[i] for i in self.ba.list]
            if len(sel) > 0:
                scene = context.scene
                gf = GeometryFactory()
                geom = gf.buildGeometry(sel)
                # geom = ShapelyOps.union(sel)
                if event.alt:
                    tM, w, h, poly, w_pts = ShapelyOps.min_bounding_rect(geom)
                    result = Io.to_curve(scene, self.coordsys, poly, 'points')
                    result.matrix_world = self.coordsys.world * tM
                    scene.objects.active = result
                else:
                    result = Io.to_curve(scene, self.coordsys, geom.convex_hull, 'points')
                    scene.objects.active = result
            self.ba.none()
            self.complete(context)
        elif event.type in {'W'}:
            sel = [self.geoms[i] for i in self.ba.list]
            if len(sel) > 0:
                scene = context.scene
                gf = GeometryFactory()
                geom = gf.buildGeometry(sel)
                if event.alt:
                    tM, w, h, poly, w_pts = ShapelyOps.min_bounding_rect(geom)
                    result = Io.to_curve(scene, self.coordsys, poly, 'points')
                    result.matrix_world = self.coordsys.world * tM
                    scene.objects.active = result
                else:
                    result = Io.to_curve(scene, self.coordsys, geom.convex_hull, 'points')
                    scene.objects.active = result
            self.ba.none()
            self.complete(context)
        elif event.type in {'D'}:
            sel = [self.geoms[i] for i in self.ba.list]
            if len(sel) > 0:
                scene = context.scene
                gf = GeometryFactory()
                geom = gf.buildGeometry(sel)
                # geom = ShapelyOps.union(sel)
                if event.alt:
                    tM, w, h, poly, w_pts = ShapelyOps.min_bounding_rect(geom)
                    result = Io.to_curve(scene, self.coordsys, poly, 'points')
                    result.matrix_world = self.coordsys.world * tM
                    scene.objects.active = result
                else:
                    result = Io.to_curve(scene, self.coordsys, geom.convex_hull, 'points')
                    scene.objects.active = result
            self.ba.none()
            self.complete(context)
        self._draw(context)

    def modal(self, context, event):
        if event.type in {'I', 'A', 'S', 'R', 'F'} and event.value == 'PRESS':
            self.keyboard(context, event)
        elif event.type in {'RIGHTMOUSE', 'ESC'}:
            self.complete(context)
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            return {'FINISHED'}
        elif event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            self.drag = True
            self.cursor_area.enable()
            self.cursor_fence.disable()
            self.startPoint = (event.mouse_region_x, event.mouse_region_y)
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            self.drag = False
            self.cursor_area.disable()
            self.cursor_fence.enable()
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
            self.select(context, self.startPoint, event)
        elif event.type == 'MOUSEMOVE':
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        return {'RUNNING_MODAL'}

    def draw_callback(self, _self, context):
        self.feedback.draw(context)
        self.cursor_area.set_location(context, self.startPoint, self.endPoint)
        self.cursor_fence.set_location(context, self.endPoint)
        self.cursor_area.draw(context)
        self.cursor_fence.draw(context)


class SelectLines(Selectable):

    def __init__(self, geoms, coordsys):
        super(SelectLines, self).__init__(geoms, coordsys)

    def _draw(self, context):
        """ override draw method """
        print("SelectLines._draw()")
        t = time.time()
        self._hide(context)
        selection = list(self.geoms[i] for i in self.ba.list)
        self.curves = [Io.to_curve(context.scene, self.coordsys, selection, 'selection', '3D')]
        for curve in self.curves:
            curve.color = (1, 1, 0, 1)
            curve.select = True
        print("SelectLines._draw() :%.2f seconds" % (time.time() - t))

    def init(self, pick_tool, context, action):
        # Post selection actions
        self.selectMode = True
        self.object_location = None
        self.startPoint = (0, 0)
        self.endPoint = (0, 0)
        self.drag = False
        self.feedback.instructions(context, "Select Lines", "Click & Drag to select lines in area", [
            ('SHIFT', 'deselect'),
            ('CTRL', 'contains'),
            ('A', 'All'),
            ('I', 'Inverse'),
            # ('F', 'Create lines from selection'),
            ('R', 'Retrieve selection'),
            ('S', 'Store selection'),
            ('ESC or RIGHTMOUSE', 'exit when done')
        ])
        self.feedback.enable()
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(
                self.draw_callback, args, 'WINDOW', 'POST_PIXEL')
        self.action = action
        self._draw(context)
        print("SelectLines.init()")

    def complete(self, context):
        print("SelectLines.complete()")
        t = time.time()
        self._hide(context)
        scene = context.scene
        selection = list(self.geoms[i] for i in self.ba.list)
        if len(selection) > 0:
            if self.action == 'select':
                result = Io.to_curve(scene, self.coordsys, selection, 'selection')
                scene.objects.active = result
            elif self.action == 'union':
                gf = GeometryFactory()
                coll = gf.buildGeometry(selection)
                # merged = gf.buildGeometry(coll.line_merge())
                # resopt = merged.simplify(tolerance=0.001, preserve_topology=False)
                resopt = gf.buildGeometry(coll.line_merge())
                """
                shapes = Io.geoms_to_shapes(selection)
                merged = ShapeOps.merge(shapes)
                union = Io.shapes_to_geoms(merged)
                # union = self.ops.union(selection)
                resopt = ShapelyOps.optimize(union)
                """
                result = Io.to_curve(scene, self.coordsys, resopt, 'union')
                scene.objects.active = result
        self.feedback.disable()
        print("SelectLines.complete() :%.2f seconds" % (time.time() - t))

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
            self.cursor_area.enable()
            self.cursor_fence.disable()
            self.startPoint = (event.mouse_region_x, event.mouse_region_y)
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            self.drag = False
            self.cursor_area.disable()
            self.cursor_fence.enable()
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
            self.select(context, self.startPoint, event)
        elif event.type == 'MOUSEMOVE':
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        return {'RUNNING_MODAL'}

    def draw_callback(self, _self, context):
        self.feedback.draw(context)
        self.cursor_area.set_location(context, self.startPoint, self.endPoint)
        self.cursor_fence.set_location(context, self.endPoint)
        self.cursor_area.draw(context)
        self.cursor_fence.draw(context)


class SelectPolygons(Selectable):

    def __init__(self, geoms, coordsys):
        super(SelectPolygons, self).__init__(geoms, coordsys)

    """
        pick_tools actions
    """
    def init(self, pick_tool, context, action):
        # Post selection actions
        self.need_rotation = False
        self.direction = 0
        self.object_location = None
        self.selectMode = True
        self.startPoint = (0, 0)
        self.endPoint = (0, 0)
        if action in ['select', 'union', 'rectangle']:
            self.feedback.instructions(context, "Select Polygons", "Click & Drag to select polygons in area", [
                ('SHIFT', 'deselect'),
                ('CTRL', 'contains'),
                ('A', 'All'),
                ('I', 'Inverse'),
                ('B', 'Bigger than current'),
                # ('F', 'Create  from selection'),
                ('R', 'Retrieve selection'),
                ('S', 'Store selection'),
                ('ESC or RIGHTMOUSE', 'exit when done')
            ])
        elif action == 'wall':
            self.feedback.instructions(context, "Select Polygons", "Click & Drag to select polygons in area", [
                ('SHIFT', 'deselect'),
                ('CTRL', 'contains'),
                ('A', 'All'),
                ('I', 'Inverse'),
                ('B', 'Bigger than current'),
                ('R', 'Retrieve selection'),
                ('S', 'Store selection'),
                ('ESC or RIGHTMOUSE', 'exit and build wall when done')
            ])
        elif action == 'window':
            self.feedback.instructions(context, "Select Polygons", "Click & Drag to select polygons in area", [
                ('SHIFT', 'deselect'),
                ('CTRL', 'contains'),
                ('A', 'All'),
                ('I', 'Inverse'),
                ('B', 'Bigger than current'),
                ('F', 'Create a window from selection'),
                ('ESC or RIGHTMOUSE', 'exit tool when done')
            ])
        elif action == 'door':
            self.feedback.instructions(context, "Select Polygons", "Click & Drag to select polygons in area", [
                ('SHIFT', 'deselect'),
                ('CTRL', 'contains'),
                ('A', 'All'),
                ('I', 'Inverse'),
                ('B', 'Bigger than current'),
                ('F', 'Create a door from selection'),
                ('ESC or RIGHTMOUSE', 'exit tool when done')
            ])
        self.gl_arc = GlPolyline((1.0, 1.0, 1.0, 0.5), d=3)
        self.gl_arc.width = 1
        self.gl_arc.style = bgl.GL_LINE_STIPPLE
        self.gl_line = GlLine(d=3)
        self.gl_line.colour_inactive = (1.0, 1.0, 1.0, 0.5)
        self.gl_line.width = 2
        self.gl_line.style = bgl.GL_LINE_STIPPLE
        self.gl_side = GlLine(d=2)
        self.gl_side.colour_inactive = (1.0, 1.0, 1.0, 0.5)
        self.gl_side.width = 2
        self.gl_side.style = bgl.GL_LINE_STIPPLE
        self.feedback.enable()
        self.drag = False
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(
                        self.draw_callback, args, 'WINDOW', 'POST_PIXEL')
        self.action = action
        self._draw(context)
        print("SelectPolygons.init()")

    def complete(self, context):
        print("SelectPolygons.complete()")
        t = time.time()
        scene = context.scene
        self._hide(context)
        selection = list(self.geoms[i] for i in self.ba.list)
        if len(selection) > 0:
            if self.action == 'select':
                result = Io.to_curve(scene, self.coordsys, selection, 'selection')
                scene.objects.active = result
            elif self.action == 'union':
                union = ShapelyOps.union(selection)
                # union = ShapelyOps.optimize(union)
                result = Io.to_curve(scene, self.coordsys, union, 'union')
                scene.objects.active = result
            elif self.action == 'wall':
                union = ShapelyOps.union(selection)
                # union = ShapelyOps.optimize(union)
                res = []
                z = context.window_manager.archipack_polylib.solidify_thickness
                Io.to_wall(scene, self.coordsys, union, z, 'wall', res)
                if len(res) > 0:
                    scene.objects.active = res[0]
                    if len(res) > 1:
                        bpy.ops.object.join()
                    bpy.ops.archipack.wall(z=z)

            elif self.action == 'rectangle':
                # currently only output a best fitted rectangle
                # over selection
                if self.object_location is not None:
                    tM, w, h, poly, w_rect = self.object_location
                    result = Io.to_curve(scene, self.coordsys, poly, 'rectangle')
                    result.matrix_world = self.coordsys.world * tM
                    scene.objects.active = result
                self.ba.none()
            elif self.action == 'window':
                if self.object_location is not None:

                    tM, w, h, l_pts, w_pts = self.object_location

                    if self.need_rotation:
                        rM = Matrix([
                            [-1, 0, 0, 0],
                            [0, -1, 0, 0],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1],
                        ])
                    else:
                        rM = Matrix()

                    if w > 1.8:
                        z = 2.2
                        altitude = 0.0
                    else:
                        z = 1.2
                        altitude = 1.0

                    bpy.ops.archipack.window(x=w, y=h, z=z, altitude=altitude, auto_manipulate=False)
                    result = context.object
                    result.matrix_world = self.coordsys.world * tM * rM
                    result.data.archipack_window[0].hole_margin = 0.02
                self.ba.none()
            elif self.action == 'door':
                if self.object_location is not None:

                    tM, w, h, l_pts, w_pts = self.object_location

                    if self.need_rotation:
                        rM = Matrix([
                            [-1, 0, 0, 0],
                            [0, -1, 0, 0],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1],
                        ])
                    else:
                        rM = Matrix()

                    if w < 1.5:
                        n_panels = 1
                    else:
                        n_panels = 2

                    bpy.ops.archipack.door(x=w, y=h, z=2.0, n_panels=n_panels,
                                direction=self.direction, auto_manipulate=False)
                    result = context.object
                    result.matrix_world = self.coordsys.world * tM * rM
                    result.data.archipack_door[0].hole_margin = 0.02
                self.ba.none()

        if self.action not in ['window', 'door']:
            self.feedback.disable()

        print("SelectPolygons.complete() :%.2f seconds" % (time.time() - t))

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
        elif event.type in {'F'}:

            sel = [self.geoms[i] for i in self.ba.list]
            if len(sel) > 0:
                if self.action == 'window':
                    self.feedback.instructions(context,
                        "Select Polygons", "Click & Drag to select polygons in area", [
                        ('CLICK & DRAG', 'Set window orientation'),
                        ('RELEASE', 'Create window'),
                        ('F', 'Return to select mode'),
                        ('ESC or RIGHTMOUSE', 'exit tool when done')
                    ])
                elif self.action == 'door':
                    self.feedback.instructions(context,
                        "Select Polygons", "Click & Drag to select polygons in area", [
                        ('CLICK & DRAG', 'Set door orientation'),
                        ('RELEASE', 'Create door'),
                        ('F', 'Return to select mode'),
                        ('ESC or RIGHTMOUSE', 'exit tool when done')
                    ])
                self.selectMode = not self.selectMode
                gf = GeometryFactory()
                geom = gf.buildGeometry(sel)
                # geom = ShapelyOps.union(sel)
                tM, w, h, poly, w_pts = ShapelyOps.min_bounding_rect(geom)
                self.object_location = (tM, w, h, poly, w_pts)
                self.startPoint = self._position_2d_from_coord(context, tM.translation)

                if self.action == 'rectangle':
                    self.complete(context)

        self._draw(context)

    def modal(self, context, event):
        if event.type in {'I', 'A', 'S', 'R', 'F', 'B'} and event.value == 'PRESS':

            self.keyboard(context, event)
        elif event.type in {'RIGHTMOUSE', 'ESC'}:
            if self.action == 'object':
                self._hide(context)
            else:
                self.complete(context)
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            return {'FINISHED'}
        elif event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            self.drag = True
            self.cursor_area.enable()
            self.cursor_fence.disable()
            if self.selectMode:
                self.startPoint = (event.mouse_region_x, event.mouse_region_y)
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            self.drag = False
            self.cursor_area.disable()
            self.cursor_fence.enable()
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
            if self.selectMode:
                self.select(context, self.startPoint, event)
            else:
                self.complete(context)
                if self.action == 'window':
                    self.feedback.instructions(context, "Select Polygons", "Click & Drag to select polygons in area", [
                        ('SHIFT', 'deselect'),
                        ('CTRL', 'contains'),
                        ('A', 'All'),
                        ('I', 'Inverse'),
                        ('B', 'Bigger than current'),
                        ('F', 'Create a window from selection'),
                        ('ESC or RIGHTMOUSE', 'exit tool when done')
                    ])
                elif self.action == 'door':
                    self.feedback.instructions(context, "Select Polygons", "Click & Drag to select polygons in area", [
                        ('SHIFT', 'deselect'),
                        ('CTRL', 'contains'),
                        ('A', 'All'),
                        ('I', 'Inverse'),
                        ('B', 'Bigger than current'),
                        ('F', 'Create a door from selection'),
                        ('ESC or RIGHTMOUSE', 'exit tool when done')
                    ])
            self.selectMode = True
        if event.type == 'MOUSEMOVE':
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
        return {'RUNNING_MODAL'}

    def _draw_2d_arc(self, context, c, p0, p1):
        """
            draw projection of 3d arc in 2d space
        """
        d0 = np.subtract(c, p0)
        d1 = np.subtract(p1, c)
        a0 = atan2(d0[1], d0[0])
        a1 = atan2(d1[1], d1[0])
        da = a1 - a0
        if da < pi:
            da += 2 * pi
        if da > pi:
            da -= 2 * pi
        da = da / 12
        r = np.linalg.norm(d1)
        pts = []
        for i in range(13):
            a = a0 + da * i
            p3d = c + Vector((cos(a) * r, sin(a) * r, 0))
            pts.append(self.coordsys.world * p3d)

        self.gl_arc.set_pos(pts)
        self.gl_arc.draw(context)
        self.gl_line.p = self.coordsys.world * c
        self.gl_line.v = pts[0] - self.gl_line.p
        self.gl_line.draw(context)

    def draw_callback(self, _self, context):
        """
            draw on screen feedback using gl.
        """
        self.feedback.draw(context)

        if self.selectMode:
            self.cursor_area.set_location(context, self.startPoint, self.endPoint)
            self.cursor_fence.set_location(context, self.endPoint)
            self.cursor_area.draw(context)
            self.cursor_fence.draw(context)
        else:
            if self.drag:
                x0, y0 = self.startPoint
                x1, y1 = self.endPoint
                # draw 2d line marker
                # self.gl.Line(x0, y0, x1, y1, self.gl.line_colour)

                # 2d line
                self.gl_side.p = Vector(self.startPoint)
                self.gl_side.v = Vector(self.endPoint) - Vector(self.startPoint)
                self.gl_side.draw(context)

                tM, w, h, l_pts, w_pts = self.object_location
                pt = self._position_3d_from_coord(context, self.endPoint)
                pt = tM.inverted() * Vector(pt)
                self.need_rotation = pt.y < 0
                if self.action == 'door':
                    # symbole porte
                    if pt.x > 0:
                        if pt.y > 0:
                            self.direction = 1
                            i_s, i_c, i_e = 3, 2, 1
                        else:
                            self.direction = 0
                            i_s, i_c, i_e = 2, 3, 0
                    else:
                        if pt.y > 0:
                            self.direction = 0
                            i_s, i_c, i_e = 0, 1, 2
                        else:
                            self.direction = 1
                            i_s, i_c, i_e = 1, 0, 3
                    self._draw_2d_arc(context, w_pts[i_c], w_pts[i_s], w_pts[i_e])
                elif self.action == 'window':
                    # symbole fenetre
                    if pt.y > 0:
                        i_s0, i_c0 = 0, 1
                        i_s1, i_c1 = 3, 2
                    else:
                        i_s0, i_c0 = 1, 0
                        i_s1, i_c1 = 2, 3
                    pc = w_pts[i_c0] + 0.5 * (w_pts[i_c1] - w_pts[i_c0])
                    self._draw_2d_arc(context, w_pts[i_c0], w_pts[i_s0], pc)
                    self._draw_2d_arc(context, w_pts[i_c1], w_pts[i_s1], pc)


class CoordSys(object):
    """
        reference coordsys
        world : matrix from local to world
        invert: matrix from world to local
        width, height: bonding region size
    """
    def __init__(self, objs):
        x = []
        y = []
        if len(objs) > 0:
            for obj in objs:
                pos = obj.location
                scale = obj.scale
                x.append(obj.bound_box[0][0] * scale.x + pos.x)
                x.append(obj.bound_box[6][0] * scale.x + pos.x)
                y.append(obj.bound_box[0][1] * scale.y + pos.y)
                y.append(obj.bound_box[6][1] * scale.y + pos.y)
        else:
            raise Exception("CoordSys require at least one object to initialize bounds")
        x0 = min(x)
        y0 = min(y)
        x1 = max(x)
        y1 = max(y)
        width, height = x1 - x0, y1 - y0
        midx, midy = x0 + width / 2.0, y0 + height / 2.0
        # reference coordsys bounding box center
        self.world = Matrix([
            [1, 0, 0, midx],
            [0, 1, 0, midy],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            ])
        self.invert = self.world.inverted()
        self.width = width
        self.height = height


class Prolongement():
    """ intersection of two segments outside segment (projection)
        c0 = point on current segment
        c1 = intersection point on oposite segment
        d = distance from ends to segment
        it = intersection validity
    """
    def __init__(self, c0, c1, d, it):
        self.c0 = c0
        self.c1 = c1
        self.d = d
        self.it = it

    def __str__(self):
        return "d:{0} valid:{1}".format(self.d, self.it[0].valid)


class Intersection():
    def __init__(self, isPoint=False):
        self.valid = True
        self.isPoint = isPoint


class Point(GeosPoint):
    """
     * A point class compatible with pygeos Coordinates
    """
    def __init__(self, coord, factory):
        GeosPoint.__init__(self, coord, factory)
        self.users = 0
        self.index = 0

    def distance(self, point):
        """ euclidian distance between points """
        return self.coord.distance(point.coord)

    def to_3d(self):
        return Vector((self.coord.x, self.coord.y, self.coord.z))

    def to_2d(self):
        return Vector((self.coord.x, self.coord.y))

    def add_user(self):
        self.users += 1


class Segment():

    def __init__(self, c0, c1, extend=EPSILON):

        # c0 c1 are Points
        self.c0 = c0
        self.c1 = c1

        self.p = c0.to_2d()
        self.v = c1.to_2d() - self.p

        self.shapeId = 0

        self._splits = []

        self.available = True
        # ensure uniqueness when merge

        self.opposite = False
        # this seg has an opposite

        self.original = False
        """
           NOTE : only original is stored into quadTree
           but quadtree is able to return _oposite when apply
        """
        self.envelope = Envelope(c0.coord, c1.coord)

    @property
    def length(self):
        """
            3d length
        """
        return self.v.length

    @property
    def cross_z(self):
        """
            2d Vector perpendicular on plane xy
            lie on the right side
            p1
            |--x
            p0
        """
        return Vector((self.v.y, -self.v.x))

    @property
    def splits(self):
        """
         splits occuring in the segment
         return points over segment
         sorted by parameter t
         filter valid ones
         remove dups
        """
        _s = sorted(self._splits, key=lambda s: s[0])
        _s = [s[1] for s in _s]  # if s[2][0].valid
        return [s for i, s in enumerate(_s) if i == 0 or _s[i - 1] is not s]

    @property
    def vect(self):
        """ vector p0-p1"""
        return self.v

    def lerp(self, t):
        """
            3d interpolation
        """
        return self.p + self.v * t

    def _point_sur_seg(self, point):
        dp = point.to_2d() - self.p
        dl = self.length
        t = (self.v * dp) / (dl * dl)
        return dp.length, t

    def _distance_point_seg(self, point):
        dp = point.to_2d() - self.p
        dl = self.length
        d = (self.v.x * dp.y - self.v.y * dp.x) / dl
        return d

    def _intersect_seg(self, segment):
        """ point_sur_segment return
            p: point d'intersection
            u: param t de l'intersection sur le segment courant
            v: param t de l'intersection sur le segment segment
            d: perpendicular distance of segment.p
        """
        c = segment.cross_z
        d = self.v * c
        if d == 0:
            d = self._distance_point_seg(segment.p)
            return False, 0, 0, 0, abs(d)
        dp = segment.p - self.p
        c2 = self.cross_z
        u = (c * dp) / d
        v = (c2 * dp) / d
        return True, self.lerp(u), u, v, 0

    def is_end(self, point):
        return point is self.c0 or point is self.c1

    def min_intersect_dist(self, t, point):
        """ distance intersection nearest end
            t: param t of intersection on segment
            point: intersection point
            return d: distance
        """
        if t > 0.5:
            return self.c1.distance(point)
        else:
            return self.c0.distance(point)

    def _append_splits(self, t, point, it):
        """
            append a unique split point
        """
        if point is self.c0 or point is self.c1:
            return
        self._splits.append((t, point, it))

    def slice(self, t, point, it):
        """
            t: param t on segment
            p: intersection point
            it:Intersection
            keep track of splits on segment
        """
        self._append_splits(t, point, it)

    def add_user(self):
        self.c0.add_user()
        self.c1.add_user()

    def consume(self):
        self.available = False

    def add_points(self, Q_segs):
        _splits = self.splits
        nSplits = len(_splits)
        if nSplits > 0:
            # kill silced segment
            self.consume()
            # build new segments
            Q_segs.newSegment(self.c0, _splits[0])
            for i in range(1, nSplits):
                Q_segs.newSegment(_splits[i - 1], _splits[i])
            Q_segs.newSegment(_splits[-1], self.c1)


class Io():

    def __init__(self, scene=None, coordsys=None, Q_segs=None, Q_points=None):
        self.scene = scene
        self.coordsys = coordsys
        self.Q_segs = Q_segs
        self.Q_points = Q_points

    @staticmethod
    def ensure_iterable(obj):
        try:
            iter(obj)
        except TypeError:
            obj = [obj]
            pass
        return obj

    # Input methods
    def _interpolate_bezier(self, pts: list, wM, p0, p1, resolution: int=12) -> None:
        # straight segment, worth testing here
        # since this can lower points count by a resolution factor
        # use normalized to handle non linear t
        if resolution == 0:
            pts.append(wM * p0.co.to_3d())
        else:
            v = (p1.co - p0.co).normalized()
            d1 = (p0.handle_right - p0.co).normalized()
            d2 = (p1.co - p1.handle_left).normalized()
            if d1 == v and d2 == v:
                pts.append(wM * p0.co.to_3d())
            else:
                seg = interpolate_bezier(wM * p0.co,
                    wM * p0.handle_right,
                    wM * p1.handle_left,
                    wM * p1.co,
                    resolution)
                for i in range(resolution - 1):
                    pts.append(seg[i].to_3d())

    def _coords_from_spline(self, wM, spline, resolution: int=12):
        pts = []
        if spline.type == 'POLY':
            pts = [wM * p.co.to_3d() for p in spline.points]
        elif spline.type == 'BEZIER':
            points = spline.bezier_points
            for i in range(1, len(points)):
                p0 = points[i - 1]
                p1 = points[i]
                self._interpolate_bezier(pts, wM, p0, p1, resolution)
            pts.append(wM * points[-1].co)
            if spline.use_cyclic_u:
                p0 = points[-1]
                p1 = points[0]
                self._interpolate_bezier(pts, wM, p0, p1, resolution)
        return pts

    def _add_curve(self, curve, resolution: int=12) -> None:
        """
         * Add a curve as segments and points in a tree
        """
        wM = self.coordsys.invert * curve.matrix_world
        for spline in curve.data.splines:
            pts = self._coords_from_spline(wM, spline, resolution)
            points = [self.Q_points.newPoint(pt) for pt in pts]
            # Ensure not unique
            if spline.use_cyclic_u:
                points.append(points[0])
            [self.Q_segs.newSegment(points[i], points[i + 1])
                for i in range(len(points) - 1)
                if points[i] is not points[i + 1]]

    def _curve_as_geom(self, gf, curve, resolution: int=12, geoms: list=[]) -> None:
        """
         * Input a curve as geom
        """
        wM = self.coordsys.invert * curve.matrix_world
        for spline in curve.data.splines:
            pts = self._coords_from_spline(wM, spline, resolution)
            # Ensure uniqueness of last point
            points = [self.Q_points.newPoint(pt).coord for pt in pts]
            if spline.use_cyclic_u:
                points.append(points[0].clone())
            geom = gf.createLineString(points)
            geoms.append(geom)

    @staticmethod
    def add_curves(Q_points, Q_segs, coordsys, curves: list, resolution: int=12) -> None:
        """
            @curves : blender curves collection
            Return coordsys for outputs
        """
        t = time.time()

        io = Io(Q_points=Q_points, Q_segs=Q_segs, coordsys=coordsys)
        for curve in curves:
            io._add_curve(curve, resolution)

        logger.debug("Io.add_curves() :%.2f seconds", time.time() - t)

    @staticmethod
    def getCoordsys(curves):
        return CoordSys(curves)

    @staticmethod
    def curves_to_geomcollection(curves, resolution: int=12, coordsys=None):
        """
         * Create lineStrings from curves
         * ensure points uniqueness using a Tree
        """
        t = time.time()
        gf = GeometryFactory()

        curves = Io.ensure_iterable(curves)

        if coordsys is None:
            coordsys = CoordSys(curves)

        Q_points = Qtree(coordsys, extend=0.5 * EPSILON)
        io = Io(Q_points=Q_points, coordsys=coordsys)

        geoms = []
        for curve in curves:
            io._curve_as_geom(gf, curve, resolution, geoms)

        # detect geometry type
        polys, dangles, cuts, invalids = PolygonizeOp.polygonize_full(geoms)

        # filter out nested touching polygon
        to_remove = []
        for poly in polys:
            for hole in poly.interiors:
                for i, other in enumerate(polys):
                    if (other is not poly and
                            hole.envelope.equals(other.envelope) and
                            CoordinateSequence.equals_unoriented(hole.coords, other.exterior.coords)):
                        to_remove.append(i)

        to_remove.sort()

        for i in reversed(to_remove):
            polys.pop(i)

        geom = gf.buildGeometry(polys + dangles + cuts + invalids)

        logger.debug("Io.curves_to_geomcollection() :%.2f seconds", time.time() - t)
        return geom

    @staticmethod
    def curves_to_geoms(curves, resolution: int=12, geoms: list=[], coordsys=None):
        """
         * Create lineStrings from curves
         * ensure points uniqueness using a Tree
        """
        t = time.time()
        gf = GeometryFactory()

        curves = Io.ensure_iterable(curves)
        if coordsys is None:
            coordsys = CoordSys(curves)

        Q_points = Qtree(coordsys, extend=0.5 * EPSILON)
        io = Io(Q_points=Q_points, coordsys=coordsys)
        for curve in curves:
            io._curve_as_geom(gf, curve, resolution, geoms)
        logger.debug("Io.curves_to_geoms() :%.2f seconds", time.time() - t)
        return coordsys

    # Output methods
    def _poly_to_wall(self, poly, height: float, name: str):

        curve = bpy.data.curves.new(name, type='CURVE')
        curve.dimensions = "2D"
        curve.fill_mode = 'BOTH'
        curve.extrude = height

        n_ext = len(poly.exterior.coords)
        n_int = len(poly.interiors)

        self._add_spline(curve, poly.exterior)

        for geom in poly.interiors:
            self._add_spline(curve, geom)

        curve_obj = bpy.data.objects.new(name, curve)
        curve_obj.matrix_world = self.coordsys.world

        self.scene.objects.link(curve_obj)

        curve_obj.select = True
        self.scene.objects.active = curve_obj

        return n_ext, n_int, curve_obj

    def wall_uv(self, me, bm):

        for face in bm.faces:
            face.select = face.material_index > 0

        bmesh.update_edit_mesh(me, True)
        bpy.ops.uv.cube_project(scale_to_bounds=False, correct_aspect=True)

        for face in bm.faces:
            face.select = face.material_index < 1

        bmesh.update_edit_mesh(me, True)
        bpy.ops.uv.smart_project(use_aspect=True, stretch_to_bounds=False)

    def _add_spline(self, curve, geometry):

        coords = list(geometry.coords)
        spline = curve.splines.new('POLY')
        spline.use_endpoint_u = False
        spline.use_cyclic_u = coords[0] == coords[-1]
        if spline.use_cyclic_u:
            coords.pop()
        spline.points.add(len(coords) - 1)
        for i, coord in enumerate(coords):
            x, y, z = coord.x, coord.y, coord.z
            spline.points[i].co = (x, y, z, 1)

    def _as_spline(self, curve, geometry):
        """
            add a spline into a blender curve
            @curve : blender curve
        """
        if geometry is None:
            logger.warning("Io._as_spline() Null geometry given")
            return
        if hasattr(geometry, 'exterior'):
            # Polygon
            self._add_spline(curve, geometry.exterior)
            for geom in geometry.interiors:
                self._add_spline(curve, geom)
        elif hasattr(geometry, 'geoms'):
            # Multi and Collections
            for geom in geometry.geoms:
                self._as_spline(curve, geom)
        else:
            # LinearRing, LineString and Shape
            self._add_spline(curve, geometry)

    def _to_curve(self, geoms, name: str, dimensions: str='3D'):
        geoms = Io.ensure_iterable(geoms)
        curve = bpy.data.curves.new(name, type='CURVE')
        curve.dimensions = dimensions
        for geom in geoms:
            self._as_spline(curve, geom)
        curve_obj = bpy.data.objects.new(name, curve)
        curve_obj.matrix_world = self.coordsys.world
        self.scene.objects.link(curve_obj)
        curve_obj.select = True
        return curve_obj

    @staticmethod
    def to_wall(scene, coordsys, geoms, height, name: str, walls: list=[]):
        """
            use curve extrude as it does respect vertices number and is not removing doubles
            so it is easy to set material index
            cap faces are tri, sides faces are quads
        """
        t = time.time()
        inside_mat = 0
        outside_mat = 1
        cut_mat = 2
        io = Io(scene=scene, coordsys=coordsys)
        bpy.ops.object.select_all(action='DESELECT')

        geoms = Io.ensure_iterable(geoms)

        for poly in geoms:
            if hasattr(poly, 'exterior'):

                half_height = height / 2.0

                n_ext, n_int, obj = io._poly_to_wall(poly, half_height, name)

                bpy.ops.object.convert(target="MESH")
                bpy.ops.object.mode_set(mode='EDIT')
                me = obj.data

                bm = bmesh.from_edit_mesh(me)
                bm.verts.ensure_lookup_table()
                bm.faces.ensure_lookup_table()

                for v in bm.verts:
                    v.co.z += half_height

                # model start from caps
                # then outside and end by inside
                # find first outside face index
                nfaces = 0
                for i, f in enumerate(bm.faces):
                    bm.faces[i].material_index = cut_mat
                    if len(f.verts) > 3:
                        nfaces = i
                        break

                # walls without interiors are inside
                if n_int > 0:
                    mat_index = outside_mat
                else:
                    mat_index = inside_mat

                # set outside if any
                for i in range(nfaces, nfaces + n_ext - 1):
                    bm.faces[i].material_index = mat_index

                # set inside
                for i in range(nfaces + n_ext - 1, len(bm.faces)):
                    bm.faces[i].material_index = inside_mat

                bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.003)
                bmesh.update_edit_mesh(me, True)

                io.wall_uv(me, bm)

                bpy.ops.mesh.dissolve_limited(angle_limit=0.00349066, delimit={'NORMAL'})
                bpy.ops.mesh.dissolve_degenerate()
                bpy.ops.object.mode_set(mode='OBJECT')
                bpy.ops.object.shade_flat()
                # MaterialUtils.add_wall_materials(obj)
                walls.append(obj)
            else:
                logger.debug("Io.to_wall() :skip %s", type(poly).__name__)

        logger.debug("Io.to_wall(%s) :%.2f seconds", len(walls), time.time() - t)

        return walls

    @staticmethod
    def to_curve(scene, coordsys, geoms, name: str="output", dimensions: str='3D'):
        t = time.time()
        io = Io(scene=scene, coordsys=coordsys)
        curve_obj = io._to_curve(geoms, name, dimensions)
        logger.debug("Io.to_curve() :%.2f seconds", time.time() - t)
        return curve_obj

    @staticmethod
    def to_curves(scene, coordsys, geoms, name: str="output", dimensions: str='3D'):
        t = time.time()
        io = Io(scene=scene, coordsys=coordsys)
        geoms = Io.ensure_iterable(geoms)
        curves = [io._to_curve(geom, name, dimensions) for geom in geoms]
        logger.debug("Io.to_curves() :%.2f seconds", time.time() - t)
        return curves

    def output(self, geoms, name: str="output", multiple: bool=True, dimensions: str='3D'):
        """
         * GeometryFactory.outputFactory use this method to
         * ouptut lines mostly for debug purposes
        """
        geoms = Io.ensure_iterable(geoms)

        if multiple:
            return [self._to_curve(geom, name, dimensions) for geom in geoms]

        if len(geoms) > 0:
            curve = self._to_curve(geoms, name, dimensions)
            return curve

    def outputCoord(self, coord, name: str="Error") -> None:
        target = bpy.data.objects.new(name, None)
        target.empty_draw_size = 1
        target.empty_draw_type = 'PLAIN_AXES'
        target.location = self.coordsys.world * Vector((coord.x, coord.y, 0))
        self.scene.objects.link(target)


class ShapelyOps():

    @staticmethod
    def min_bounding_rect(geom):
        """ min_bounding_rect
            minimum area oriented bounding rect
        """
        rect, transf_rect, inv_matrix = geom.computeMinimumRotatedRectangle()
        if rect is None:
            logger.debug("rect is None")
            return None

        w = transf_rect.envelope.width
        h = transf_rect.envelope.height
        centre = Coordinate()
        rect.envelope.centre(centre)
        ux, vx, uy, vy, px, py = inv_matrix

        if h > w:
            w, h = h, w
            ux, uy, vx, vy = vx, vy, -ux, -uy

        x = w / 2
        y = h / 2

        cs = geom._factory.coordinateSequenceFactory.create(
            [Coordinate(co[0], co[1]) for co in [(-x, -y), (-x, y), (x, y), (x, -y), (-x, -y)]])
            
        centered_rect = geom._factory.createLinearRing(cs)

        tM = Matrix([
            [ux, vx, 0, centre.x],
            [uy, vy, 0, centre.y],
            [0, 0, 1, 0],
            [0, 0, 0, 1]])

        return tM, w, h, centered_rect, rect.coords

    @staticmethod
    def detect_polygons(geoms):
        """ detect_polygons
        """
        print("Ops.detect_polygons()")
        t = time.time()
        result, dangles, cuts, invalids = PolygonizeOp.polygonize_full(geoms)
        print("Ops.detect_polygons() :%.2f seconds" % (time.time() - t))
        return result, dangles, cuts, invalids

    @staticmethod
    def optimize(geoms, tolerance=0.001, preserve_topology=False):
        """ optimize
        """
        t = time.time()
        geoms = Io.ensure_iterable(geoms)
        optimized = [geom.simplify(tolerance, preserve_topology) for geom in geoms]
        print("Ops.optimize() :%.2f seconds" % (time.time() - t))
        return optimized

    @staticmethod
    def union(geoms):
        """ fast union
            cascaded union - may require snap before use to fix precision issues
        """
        t = time.time()
        geoms = Io.ensure_iterable(geoms)
        union = PolygonsUnionOp.union(geoms)
        print("Ops.union() :%.2f seconds" % (time.time() - t))
        return union

    @staticmethod
    def boolean(a, b, opCode):
        if opCode == 1:
            return a.intersection(b)
        elif opCode == 2:
            return a.union(b)
        elif opCode == 3:
            return a.difference(b)
        elif opCode == 4:
            return a.symmetric_difference(b)


class Qtree(_QuadTree):
    """
        The top spatial index to be created by the user. Once created it can be
        populated with geographically placed members that can later be tested for
        intersection with a user inputted geographic bounding box.
    """
    def __init__(self, coordsys, extend=EPSILON, max_items=MAX_ITEMS, max_depth=MAX_DEPTH):
        """
            objs may be blender objects or shapely geoms
            extend: how much seek arround
        """
        self._extend = extend
        self._geoms = []

        # store input coordsys
        self.coordsys = coordsys

        super(Qtree, self).__init__(0, 0, coordsys.width, coordsys.height, max_items, max_depth)

        self._factory = GeometryFactory()

    @property
    def ngeoms(self):
        return len(self._geoms)

    def build(self, geoms):
        """
            Build a spacial index from shapely geoms
        """
        t = time.time()
        self._geoms = geoms
        for i, geom in enumerate(geoms):
            self._insert(i, self.getbounds(geom))
        logger.debug("Qtree.build() :%.2f seconds", time.time() - t)

    def insert(self, id, geom):
        self._geoms.append(geom)
        self._insert(id, self.getbounds(geom))

    def newPoint(self, co):

        point = Point(Coordinate(co.x, co.y), self._factory)
        count, found = self.intersects_ext(point, EPSILON)
        for id in found:
            return self._geoms[id]
        self.insert(self.ngeoms, point)
        return point

    def newSegment(self, c0, c1):
        """
            allow "opposite" segments,
            those segments are not found by intersects
            and not stored in self.geoms
        """
        new_seg = Segment(c0, c1, self._extend)
        count, found = self.intersects(new_seg)
        for id in found:
            old_seg = self._geoms[id]
            if (old_seg.c0 == c0 and old_seg.c1 == c1):
                return old_seg
            if (old_seg.c0 == c1 and old_seg.c1 == c0):
                if not old_seg.opposite:
                    old_seg.opposite = new_seg
                    new_seg.original = old_seg
                return old_seg.opposite
        self.insert(self.ngeoms, new_seg)
        return new_seg

    def getbounds(self, geom):
        if hasattr(geom, 'bounds'):
            return geom.bounds
        else:
            _b = geom.envelope
            return (_b.minx - EPSILON, _b.miny - EPSILON, _b.maxx + EPSILON, _b.maxy + EPSILON)

    def intersects_ext(self, geom, extend):
        bounds = self.getbounds(geom)
        _bounds = (bounds[0] - extend,
                    bounds[1] - extend,
                    bounds[2] + extend,
                    bounds[3] + extend)
        selection = list(self._intersect(_bounds))
        count = len(selection)
        return count, sorted(selection)

    def intersects(self, geom):
        bounds = self.getbounds(geom)
        selection = list(self._intersect(bounds))
        count = len(selection)
        return count, sorted(selection)


class Polygonizer():
    """
        Define collection of shapes as polylines and polygons
        detect polygons and classify boundary / interiors
        equivalent to Shapely one
    """
    def __init__(self, coordsys, extend=EPSILON):

        self.extend = extend

        self.coordsys = coordsys

        # Errors (shapes without left boundarys)
        self.err = []

        # points under walls
        self.inside_wall = []

    def _intersection_point(self, d, t, point, seg):
        if d > EPSILON:
            return point
        elif t > 0.5:
            return seg.c1
        else:
            return seg.c0

    def split(self, Q_points, Q_segs, extend=0.01, collinear=True):
        """ _split
            detect intersections between segments and create segments according
            is able to project segment ends on closest segment
            use point_tree and seg_tree
            assume merged segments on beforehand
        """
        t = time.time()

        # debug
        processed = 0
        process_collinear = 0
        process_point = 0

        segs = Q_segs._geoms
        nbsegs = Q_segs.ngeoms
        it_start = [None for x in range(nbsegs)]
        it_end = [None for x in range(nbsegs)]

        # points with 1 user are "extendable"
        for seg in Q_segs._geoms:
            seg.c0.add_user()
            seg.c1.add_user()

        for s, seg in enumerate(segs):

            # enlarge seek box for "extendable" segments
            if seg.c0.users < 2 or seg.c1.users < 2:
                count, idx = Q_segs.intersects_ext(seg, extend)
            else:
                count, idx = Q_segs.intersects(seg)

            for id in idx:

                if id > s:
                    # can't check for side determinant here
                    # as intersect may enlarge segment to nearest
                    # neighboor

                    _seg = segs[id]
                    intersect, co, u, v, d = seg._intersect_seg(_seg)

                    processed += 1

                    # POINT_INTERSECTION
                    if intersect:

                        process_point += 1

                        # store intersection state
                        # to be able to disable invalid ones
                        it = [Intersection()]

                        point = Q_points.newPoint(co)

                        # distance of nearest segment endpoint and intersection
                        du = seg.min_intersect_dist(u, point)
                        dv = _seg.min_intersect_dist(v, point)

                        # print("s:%s id:%s u:%7f v:%7f du:%7f dv:%7f" % (s, id, u, v, du, dv))

                        # does intersect realy occurs ?
                        if (((v <= 0 and _seg.c0.users < 2 and dv < extend) or
                                (0 <= v < 1) or
                                (v >= 1 and _seg.c1.users < 2 and dv < extend)) and
                                ((u <= 0 and seg.c0.users < 2 and du < extend) or
                                (0 <= u < 1) or
                                (u >= 1 and seg.c1.users < 2 and du < extend))):

                            # intersection point on segment id,
                            # segment end when distance is under precision

                            pt = self._intersection_point(dv, v, point, _seg)

                            # make last intersections invalid
                            # on both last and oposite segments
                            # prevent segment from being "extendable"
                            if pt is seg.c0:
                                if it_start[s] is not None:
                                    it_start[s].it[0].valid = False
                                    it_start[s].d = 0
                                # seg.c0.users = 2

                            elif pt is seg.c1:
                                if it_end[s] is not None:
                                    it_end[s].it[0].valid = False
                                    it_end[s].d = 0
                                # seg.c1.users = 2

                            elif u <= 0:
                                # enlarge segment s c0
                                last = it_start[s]
                                if last is None:
                                    it_start[s] = Prolongement(seg.c0, pt, du, it)
                                    # elif last.c1 is pt:
                                    #    it[0] = last.it[0]
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = seg.c0
                                    last.c1 = pt
                                    last.it = it
                            elif u < 1:
                                # intersection on segment s
                                seg.slice(u, pt, it)
                            else:
                                # enlarge segment s c1
                                last = it_end[s]
                                if last is None:
                                    it_end[s] = Prolongement(seg.c1, pt, du, it)
                                    # elif last.c1 is pt:
                                    # it[0] = last.it[0]
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = seg.c1
                                    last.c1 = pt
                                    last.it = it

                            # intersection point on segment seg,
                            # segment end when distance is under precision
                            pt = self._intersection_point(du, u, point, seg)

                            # make last intersections invalid
                            # on oposite segment
                            if pt is _seg.c0:
                                if it_start[id] is not None:
                                    it_start[id].it[0].valid = False
                                    it_start[id].d = 0
                                # _seg.c0.users = 2

                            elif pt is _seg.c1:
                                if it_end[id] is not None:
                                    it_end[id].it[0].valid = False
                                    it_end[id].d = 0
                                # _seg.c1.users = 2

                            elif v <= 0:
                                # enlarge segment id c0
                                last = it_start[id]
                                if last is None:
                                    it_start[id] = Prolongement(_seg.c0, pt, dv, it)
                                    # elif last.c1 is pt:
                                    # it[0] = last.it[0]
                                elif dv < last.d:
                                    last.it[0].valid = False
                                    last.d = dv
                                    last.c0 = _seg.c0
                                    last.c1 = pt
                                    last.it = it

                            elif v < 1:
                                # intersection on segment id
                                _seg.slice(v, pt, it)
                            else:
                                # enlarge segment id c1
                                last = it_end[id]
                                if last is None:
                                    it_end[id] = Prolongement(_seg.c1, pt, dv, it)
                                    # elif last.c1 is pt:
                                    # it[0] = last.it[0]
                                elif dv < last.d:
                                    last.it[0].valid = False
                                    last.d = dv
                                    last.c0 = _seg.c1
                                    last.c1 = pt
                                    last.it = it

                    # COLLINEAR_INTERSECTION
                    elif collinear and d < EPSILON:
                        # parallel segments, endpoint on segment
                        # skip NON_COLLINEAR aka when d > EPSILON
                        process_collinear += 1

                        # point _seg.c0 on segment seg
                        pt = _seg.c0
                        du, u = seg._point_sur_seg(pt)

                        if (((u <= 0 and seg.c0.users < 2 and du < extend) or
                                (0 <= u < 1) or
                                (u >= 1 and seg.c1.users < 2 and du < extend))):

                            it = [Intersection()]
                            if u <= 0:

                                # extend seg on c0 side
                                last = it_start[s]
                                if last is None:
                                    it_start[s] = Prolongement(seg.c0, pt, du, it)
                                    # elif last.c1 is pt:
                                    # it[0] = last.it[0]
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = seg.c0
                                    last.c1 = pt
                                    last.it = it

                            elif u < 1:
                                # occurs inside segment seg
                                seg.slice(u, pt, it)
                                # _seg.c0 must not extend
                                # _seg.c0.users = 2
                                if it_start[id] is not None:
                                    it_start[id].it[0].valid = False

                            else:
                                # extend seg on c1 side
                                last = it_end[s]
                                if last is None:
                                    it_end[s] = Prolongement(seg.c1, pt, du, it)
                                    # elif last.c1 is pt:
                                    # it[0] = last.it[0]
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = seg.c1
                                    last.c1 = pt
                                    last.it = it

                        # point _seg.c1 on segment seg
                        pt = _seg.c1
                        du, u = seg._point_sur_seg(pt)

                        if (((u <= 0 and seg.c0.users < 2 and du < extend) or
                                (0 <= u < 1) or
                                (u >= 1 and seg.c1.users < 2 and du < extend))):

                            it = [Intersection()]
                            if u <= 0:

                                # extend in c0 side
                                last = it_start[s]
                                if last is None:
                                    it_start[s] = Prolongement(seg.c0, pt, du, it)
                                    # elif last.c1 is pt:
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = seg.c0
                                    last.c1 = pt
                                    last.it = it

                            elif u < 1:
                                # occurs on segment seg
                                seg.slice(u, pt, it)
                                # _seg.c1.users = 2
                                if it_end[id] is not None:
                                    it_end[id].it[0].valid = False

                            else:
                                # extend
                                last = it_end[s]
                                if last is None:
                                    it_end[s] = Prolongement(seg.c1, pt, du, it)
                                    # elif last.c1 is pt:
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = seg.c1
                                    last.c1 = pt
                                    last.it = it

                        # point seg.c0 on segment _seg
                        pt = seg.c0
                        du, u = _seg._point_sur_seg(pt)

                        if (((u <= 0 and _seg.c0.users < 2 and du < extend) or
                                (0 <= u < 1) or
                                (u >= 1 and _seg.c1.users < 2 and du < extend))):

                            it = [Intersection()]
                            if u <= 0:

                                # extend in c0 side
                                last = it_start[id]
                                if last is None:
                                    it_start[id] = Prolongement(_seg.c0, pt, du, it)
                                    # elif last.c1 is pt:
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = _seg.c0
                                    last.c1 = pt
                                    last.it = it

                            elif u < 1:
                                # occurs on segment _seg
                                # always occurs so intersection
                                # doesent need to be "removable"
                                _seg.slice(u, pt, it)
                                # seg.c0.users = 2
                                if it_start[s] is not None:
                                    it_start[s].it[0].valid = False

                            else:
                                # extend
                                last = it_end[id]
                                if last is None:
                                    it_end[id] = Prolongement(_seg.c1, pt, du, it)
                                    # elif last.c1 is pt:
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = _seg.c1
                                    last.c1 = pt
                                    last.it = it

                        # point seg.c1 on segment _seg
                        pt = seg.c1
                        du, u = _seg._point_sur_seg(pt)

                        if (((u <= 0 and _seg.c0.users < 2 and du < extend) or
                                (0 <= u < 1) or
                                (u >= 1 and _seg.c1.users < 2 and du < extend))):

                            it = [Intersection()]
                            if u <= 0:

                                # extend _seg on c0 side
                                last = it_start[id]
                                if last is None:
                                    it_start[id] = Prolongement(_seg.c0, pt, du, it)
                                    # elif last.c1 is pt:
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = _seg.c0
                                    last.c1 = pt
                                    last.it = it

                            elif u < 1:
                                # occurs on segment _seg
                                _seg.slice(u, pt, it)
                                # seg.c1.users = 2
                                if it_end[s] is not None:
                                    it_end[s].it[0].valid = False

                            else:
                                # extend _seg on c1 side
                                last = it_end[id]
                                if last is None:
                                    it_end[id] = Prolongement(_seg.c1, pt, du, it)
                                    # elif last.c1 is pt:
                                elif du < last.d:
                                    last.it[0].valid = False
                                    last.d = du
                                    last.c0 = _seg.c1
                                    last.c1 = pt
                                    last.it = it

        for pro in it_start:
            if pro is not None and pro.it[0].valid:
                Q_segs.newSegment(pro.c0, pro.c1)

        for pro in it_end:
            if pro is not None and pro.it[0].valid:
                Q_segs.newSegment(pro.c0, pro.c1)

        logger.debug("Polygonizer.split() intersect (all:%s, points:%s, collinear:%s) :%.4f seconds",
            processed,
            process_point,
            process_collinear,
            (time.time() - t))

        t = time.time()
        _geoms = Q_segs._geoms
        for seg in _geoms:
            seg.add_points(Q_segs)

        logger.debug("Polygonizer.split() slice :%.4f seconds", (time.time() - t))

    @staticmethod
    def polygonize(context, curves, extend=0.0, resolution=12):

        t = time.time()
        curves = Io.ensure_iterable(curves)
        coordsys = CoordSys(curves)

        gf = GeometryFactory()
        gf.outputFactory = Io(scene=context.scene, coordsys=coordsys)

        op = Polygonizer(coordsys)
        # Ensure uniqueness of points and segments
        Q_segs = Qtree(coordsys, extend=extend)
        Q_points = Qtree(coordsys, extend=0.5 * EPSILON)

        Io.add_curves(Q_points, Q_segs, coordsys, curves, resolution)

        op.split(Q_points, Q_segs, extend)

        lines = gf.buildGeometry([gf.createLineString([seg.c0.coord, seg.c1.coord])
            for seg in Q_segs._geoms
            if seg.available and seg.c0 is not seg.c1])

        merged = lines.line_merge()
        # skip validity test for polygons
        polys, dangles, cuts, invalids = PolygonizeOp.polygonize_full(merged, skip_validity_check=True)
        vars_dict['select_polygons'] = SelectPolygons(polys, coordsys)
        vars_dict['select_lines'] = SelectLines(merged, coordsys)
        vars_dict['select_points'] = SelectPoints(Q_points._geoms, coordsys)

        logger.debug("Polygonizer.polygonize() :%.2f seconds polygons:%s invalids:%s",
            time.time() - t,
            len(polys),
            len(invalids))

        return coordsys, polys, dangles, cuts, invalids


class ARCHIPACK_OP_PolyLib_Pick2DPoints(Operator):
    bl_idname = "archipack.polylib_pick_2d_points"
    bl_label = "Pick lines"
    bl_description = "Pick lines"
    bl_options = {'REGISTER', 'UNDO'}
    pass_keys = ['NUMPAD_0', 'NUMPAD_1', 'NUMPAD_3', 'NUMPAD_4',
                 'NUMPAD_5', 'NUMPAD_6', 'NUMPAD_7', 'NUMPAD_8',
                 'NUMPAD_9', 'MIDDLEMOUSE', 'WHEELUPMOUSE', 'WHEELDOWNMOUSE']
    action = StringProperty(name="action", default="select")

    @classmethod
    def poll(self, context):
        global vars_dict
        return vars_dict['select_points'] is not None

    def modal(self, context, event):
        global vars_dict
        context.area.tag_redraw()
        if event.type in self.pass_keys:
            return {'PASS_THROUGH'}
        return vars_dict['select_points'].modal(context, event)

    def invoke(self, context, event):
        global vars_dict
        if vars_dict['select_points'] is None:
            self.report({'WARNING'}, "Use detect before")
            return {'CANCELLED'}
        elif context.space_data.type == 'VIEW_3D':
            vars_dict['select_points'].init(self, context, self.action)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}


class ARCHIPACK_OP_PolyLib_Pick2DLines(Operator):
    bl_idname = "archipack.polylib_pick_2d_lines"
    bl_label = "Pick lines"
    bl_description = "Pick lines"
    bl_options = {'REGISTER', 'UNDO'}
    pass_keys = ['NUMPAD_0', 'NUMPAD_1', 'NUMPAD_3', 'NUMPAD_4',
                 'NUMPAD_5', 'NUMPAD_6', 'NUMPAD_7', 'NUMPAD_8',
                 'NUMPAD_9', 'MIDDLEMOUSE', 'WHEELUPMOUSE', 'WHEELDOWNMOUSE']
    action = StringProperty(name="action", default="select")

    @classmethod
    def poll(self, context):
        global vars_dict
        return vars_dict['select_lines'] is not None

    def modal(self, context, event):
        global vars_dict
        context.area.tag_redraw()
        if event.type in self.pass_keys:
            return {'PASS_THROUGH'}
        return vars_dict['select_lines'].modal(context, event)

    def invoke(self, context, event):
        global vars_dict
        if vars_dict['select_lines'] is None:
            self.report({'WARNING'}, "Use detect before")
            return {'CANCELLED'}
        elif context.space_data.type == 'VIEW_3D':
            vars_dict['select_lines'].init(self, context, self.action)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}


class ARCHIPACK_OP_PolyLib_Pick2DPolygons(Operator):
    bl_idname = "archipack.polylib_pick_2d_polygons"
    bl_label = "Pick 2d"
    bl_description = "Pick polygons"
    bl_options = {'REGISTER', 'UNDO'}
    pass_keys = ['NUMPAD_0', 'NUMPAD_1', 'NUMPAD_3', 'NUMPAD_4',
                 'NUMPAD_5', 'NUMPAD_6', 'NUMPAD_7', 'NUMPAD_8',
                 'NUMPAD_9', 'MIDDLEMOUSE', 'WHEELUPMOUSE', 'WHEELDOWNMOUSE']
    action = StringProperty(name="action", default="select")

    @classmethod
    def poll(self, context):
        global vars_dict
        return vars_dict['select_polygons'] is not None

    def modal(self, context, event):
        global vars_dict
        context.area.tag_redraw()
        if event.type in self.pass_keys:
            return {'PASS_THROUGH'}
        return vars_dict['select_polygons'].modal(context, event)

    def invoke(self, context, event):
        global vars_dict
        if vars_dict['select_polygons'] is None:
            self.report({'WARNING'}, "Use detect before")
            return {'CANCELLED'}
        elif context.space_data.type == 'VIEW_3D':
            vars_dict['select_polygons'].init(self, context, self.action)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}


class ARCHIPACK_OP_PolyLib_Detect(Operator):
    bl_idname = "archipack.polylib_detect"
    bl_label = "Detect Polygons"
    bl_description = "Detect polygons from unordered splines"
    bl_options = {'REGISTER', 'UNDO'}
    extend = FloatProperty(name="extend", default=0.01, subtype='DISTANCE', unit='LENGTH', min=0)

    @classmethod
    def poll(self, context):
        return (len(context.selected_objects) > 0 and
            context.object is not None and
            context.object.type == 'CURVE')

    def execute(self, context):
        global vars_dict

        objs = [obj for obj in context.selected_objects if obj.type == 'CURVE']

        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}

        for obj in objs:
            obj.select = False

        try:
            coordsys, polys, dangles, cuts, invalids = Polygonizer.polygonize(
                context,
                objs,
                extend=self.extend,
                resolution=context.window_manager.archipack_polylib.resolution)
        except TopologyException as ex:
            self.report({'WARNING'}, "Topology error {}".format(ex))
            return {'CANCELLED'}
        except:
            self.report({'WARNING'}, "Unknown error")
            return {'CANCELLED'}

        if len(invalids) > 0:
            errs = Io.to_curve(context.scene, coordsys, invalids, "invalid_polygons")
            err_mat = vars_dict['select_polygons'].build_display_mat("Invalid_polygon", (1, 0, 0))
            # curve.data.bevel_depth = 0.02
            errs.color = (1, 0, 0, 1)
            if len(errs.data.materials) < 1:
                errs.data.materials.append(err_mat)
                errs.active_material = err_mat
            errs.select = True
            self.report({'WARNING'}, str(len(invalids)) + " invalid polygons detected")

        return {'FINISHED'}


class ARCHIPACK_OP_PolyLib_Offset(Operator):
    """Offset curves as separate objects
    """
    bl_idname = "archipack.polylib_offset"
    bl_label = "Offset"
    bl_description = "Offset lines"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return len(context.selected_objects) > 0 and context.object is not None and context.object.type == 'CURVE'

    def execute(self, context):
        t = time.time()
        wm = context.window_manager.archipack_polylib
        objs = list(obj for obj in context.selected_objects if obj.type == 'CURVE')

        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}

        for obj in objs:
            obj.select = False

        lines = []

        coordsys = Io.curves_to_geoms(objs, wm.resolution, lines)
        gf = GeometryFactory()
        gf.outputFactory = Io(scene=context.scene, coordsys=coordsys)
        offset = []

        distance = wm.offset_distance

        if wm.offset_side == 'right':
            distance = -distance
        """
        offset = coll.parallel_offset(distance,
                    resolution=wm.offset_resolution,
                    join_style=int(wm.offset_join_style),
                    mitre_limit=wm.offset_mitre_limit)
        """

        for line in lines:
            try:
                res = line.parallel_offset(distance, resolution=wm.offset_resolution,
                        join_style=int(wm.offset_join_style), mitre_limit=wm.offset_mitre_limit)
            except TopologyException as ex:
                self.report({'WARNING'}, "Topology error {}".format(ex))
                return {'CANCELLED'}
            except:
                self.report({'WARNING'}, "Unknown error")
                return {'CANCELLED'}
            offset.append(res)

        result = Io.to_curve(context.scene, coordsys, offset, 'offset')
        context.scene.objects.active = result
        logger.info("Offset :%.2f seconds", time.time() - t)
        
        return {'FINISHED'}


class ARCHIPACK_OP_PolyLib_Buffer(Operator):
    bl_idname = "archipack.polylib_buffer"
    bl_label = "Buffer"
    bl_description = "Buffer lines"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return len(context.selected_objects) > 0 and context.object is not None and context.object.type == 'CURVE'

    def execute(self, context):
        t = time.time()
        wm = context.window_manager.archipack_polylib
        objs = list(obj for obj in context.selected_objects if obj.type == 'CURVE')

        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}

        for obj in objs:
            obj.select = False

        lines = []

        coordsys = Io.curves_to_geoms(objs, wm.resolution, lines)
        gf = GeometryFactory()
        gf.outputFactory = Io(scene=context.scene, coordsys=coordsys)
        offset = []

        distance = wm.buffer_distance

        if wm.buffer_side == 'right':
            distance = -distance

        coll = gf.buildGeometry(lines)
        try:
            offset = coll.buffer(distance,
                        resolution=wm.buffer_resolution,
                        join_style=int(wm.buffer_join_style),
                        cap_style=int(wm.buffer_cap_style),
                        mitre_limit=wm.buffer_mitre_limit,
                        single_sided=wm.buffer_side != 'both'
                        )
        except TopologyException as ex:
            self.report({'WARNING'}, "Topology error {}".format(ex))
            return {'CANCELLED'}
        except:
            self.report({'WARNING'}, "Unknown error")
            return {'CANCELLED'}
        """
        for line in lines:
            res = line.buffer(distance, resolution=wm.buffer_resolution,
                        join_style=int(wm.buffer_join_style),
                        cap_style=int(wm.buffer_join_style),
                        mitre_limit=wm.buffer_mitre_limit,
                        single_sided=wm.buffer_side != 'both'
                        )
            offset.append(res)
        """

        result = Io.to_curve(context.scene, coordsys, offset, 'buffer')
        context.scene.objects.active = result
        logger.info("Buffer :%.2f seconds", time.time() - t)
        return {'FINISHED'}


class ARCHIPACK_OP_PolyLib_Boolean(Operator):
    bl_idname = "archipack.polylib_boolean"
    bl_label = "Boolean"
    bl_description = "Boolean operation"
    bl_options = {'REGISTER', 'UNDO'}
    opCode = IntProperty(default=0)

    @classmethod
    def poll(self, context):
        return (len(context.selected_objects) > 0 and
            context.active_object is not None and
            context.active_object.type == 'CURVE')

    def execute(self, context):
        t = time.time()
        a = context.active_object
        b = [obj for obj in context.selected_objects if obj.type == 'CURVE' and obj.name != a.name]
        wm = context.window_manager.archipack_polylib

        if len(b) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}

        for obj in b:
            obj.select = False

        a.select = False

        coordsys = Io.getCoordsys([a] + b)

        geom_a = Io.curves_to_geomcollection([a], wm.resolution, coordsys=coordsys)
        geom_b = Io.curves_to_geomcollection(b, wm.resolution, coordsys=coordsys)

        opCode = self.opCode
        if opCode > 10:
            opCode = opCode - 10
            geom_a, geom_b = geom_b, geom_a
            
        try:
            # Might throw TopologyException
            res = ShapelyOps.boolean(geom_a, geom_b, opCode)

        except TopologyException as ex:
            self.report({'WARNING'}, "Topology error {}".format(ex))
            return {'CANCELLED'}
        except:
            self.report({'WARNING'}, "Unknown error")
            return {'CANCELLED'}
        
        result = Io.to_curve(context.scene, coordsys, res, 'boolean')
        context.scene.objects.active = result
        logger.info("Boolean :%.2f seconds", time.time() - t)
        
        return {'FINISHED'}


class ARCHIPACK_OP_PolyLib_Simplify(Operator):
    bl_idname = "archipack.polylib_simplify"
    bl_label = "Simplify"
    bl_description = "Simplify lines"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return (len(context.selected_objects) > 0 and
            context.object is not None and
            context.object.type == 'CURVE')

    def execute(self, context):
        t = time.time()
        global vars_dict
        wm = context.window_manager.archipack_polylib
        objs = [obj for obj in context.selected_objects if obj.type == 'CURVE']
        if len(objs) < 1:
            self.report({'WARNING'}, "Select a curve object before")
            return {'CANCELLED'}
        for obj in objs:
            obj.select = False
        simple = []
        lines = []
        coordsys = Io.curves_to_geoms(objs, wm.resolution, lines)
        for line in lines:
            res = line.simplify(wm.simplify_tolerance, preserve_topology=wm.simplify_preserve_topology)
            simple.append(res)
        result = Io.to_curve(context.scene, coordsys, simple, 'simplify')
        context.scene.objects.active = result
        logger.info("Simplify :%.2f seconds", time.time() - t)
        return {'FINISHED'}


class ARCHIPACK_OP_PolyLib_OutputPolygons(Operator):
    bl_idname = "archipack.polylib_output_polygons"
    bl_label = "Output Polygons"
    bl_description = "Output all polygons"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        global vars_dict
        return vars_dict['select_polygons'] is not None

    def execute(self, context):
        global vars_dict
        result = Io.to_curve(context.scene, vars_dict['select_polygons'].coordsys,
                                vars_dict['select_polygons'].geoms, 'polygons')
        context.scene.objects.active = result
        return {'FINISHED'}


class ARCHIPACK_OP_PolyLib_OutputLines(Operator):
    bl_idname = "archipack.polylib_output_lines"
    bl_label = "Output lines"
    bl_description = "Output all lines"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        global vars_dict
        return vars_dict['select_lines'] is not None

    def execute(self, context):
        global vars_dict
        result = Io.to_curve(context.scene, vars_dict['select_lines'].coordsys,
                                vars_dict['select_lines'].geoms, 'lines')
        context.scene.objects.active = result
        return {'FINISHED'}


class ARCHIPACK_OP_PolyLib_Solidify(Operator):
    bl_idname = "archipack.polylib_solidify"
    bl_label = "Extrude"
    bl_description = "Extrude all polygons"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return (len(context.selected_objects) > 0 and
                context.object is not None and
                context.object.type == 'CURVE')

    def execute(self, context):
        wm = context.window_manager.archipack_polylib
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


class archipack_polylib(PropertyGroup):
    bl_idname = 'archipack.polylib_parameters'
    polygonize_expand = BoolProperty(default=False, description="Display polygonize tools")
    extend = FloatProperty(
            name="Extend",
            description="Extend to closest intersecting segment",
            default=0.01,
            subtype='DISTANCE', unit='LENGTH', min=0
            )
    resolution = IntProperty(
            name="Bezier resolution", min=0, default=12
            )

    offset_expand = BoolProperty(default=False, description="Display options")
    offset_distance = FloatProperty(
            name="Distance",
            default=0.05,
            subtype='DISTANCE', unit='LENGTH', min=0
            )
    offset_side = EnumProperty(
            name="Side", default='left',
            items=[('left', 'Left', 'Left'),
                ('right', 'Right', 'Right')]
            )
    offset_resolution = IntProperty(
            name="Resolution", default=16, min=0
            )
    offset_join_style = EnumProperty(
            name="Style", default='2',
            items=[('1', 'Round', 'Round'),
                    ('2', 'Mitre', 'Mitre'),
                    ('3', 'Bevel', 'Bevel')]
            )
    offset_mitre_limit = FloatProperty(
            name="Mitre limit",
            default=10.0,
            subtype='DISTANCE',
            unit='LENGTH', min=0
            )

    buffer_expand = BoolProperty(default=False, description="Display options")
    buffer_distance = FloatProperty(
            name="Distance",
            default=0.05,
            subtype='DISTANCE', unit='LENGTH', min=0
            )
    buffer_side = EnumProperty(
            name="Side", default='both',
            items=[('both', 'Both', 'Both'),
                ('left', 'Left', 'Left'),
                ('right', 'Right', 'Right')]
            )
    buffer_resolution = IntProperty(
            name="Resolution", default=16, min=0
            )
    buffer_join_style = EnumProperty(
            name="Join", default='2',
            items=[('1', 'Round', 'Round'),
                    ('2', 'Mitre', 'Mitre'),
                    ('3', 'Bevel', 'Bevel')]
            )
    buffer_cap_style = EnumProperty(
            name="Cap", default='3',
            items=[('1', 'Round', 'Round'),
                    ('2', 'Flat', 'Flat'),
                    ('3', 'Square', 'Square')]
            )
    buffer_mitre_limit = FloatProperty(
            name="Mitre limit",
            default=10.0,
            subtype='DISTANCE',
            unit='LENGTH', min=0
            )

    simplify_expand = BoolProperty(default=False, description="Display options")
    simplify_tolerance = FloatProperty(
            name="Tolerance",
            default=0.01,
            subtype='DISTANCE', unit='LENGTH', min=0
            )
    simplify_preserve_topology = BoolProperty(
            name="Preserve topology",
            description="Preserve topology (fast without, but may introduce self crossing)",
            default=False
            )
    solidify_thickness = FloatProperty(
            name="Thickness",
            default=2.7,
            subtype='DISTANCE', unit='LENGTH', min=0
            )
    boolean_expand = BoolProperty(default=False, description="Display 2d boolean tools")


@persistent
def load_handler(dummy):
    global vars_dict
    vars_dict['select_polygons'] = None
    vars_dict['select_lines'] = None
    vars_dict['select_points'] = None


def register():
    global vars_dict
    vars_dict = {
        # keep track of shapely geometry selection sets
        'select_polygons': None,
        'select_lines': None,
        'select_points': None
        }
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_Pick2DPolygons)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_Pick2DLines)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_Pick2DPoints)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_OutputPolygons)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_OutputLines)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_Offset)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_Buffer)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_Boolean)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_Simplify)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_Detect)
    bpy.utils.register_class(ARCHIPACK_OP_PolyLib_Solidify)
    bpy.utils.register_class(archipack_polylib)
    bpy.types.WindowManager.archipack_polylib = PointerProperty(type=archipack_polylib)
    bpy.app.handlers.load_post.append(load_handler)


def unregister():
    global vars_dict
    del vars_dict
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_Pick2DPolygons)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_Pick2DLines)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_Pick2DPoints)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_Detect)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_OutputPolygons)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_OutputLines)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_Offset)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_Buffer)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_Boolean)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_Simplify)
    bpy.utils.unregister_class(ARCHIPACK_OP_PolyLib_Solidify)
    bpy.utils.unregister_class(archipack_polylib)
    bpy.app.handlers.load_post.remove(load_handler)
    del bpy.types.WindowManager.archipack_polylib


"""
        # x, y = pt
        # loc = p.coordsys.world * Vector((x, y, 0))
        # bpy.ops.object.empty_add(type='PLAIN_AXES', radius=1,location=loc)


from archipack.archipack_polylines import Polygonizer
curves = C.selected_objects
coordsys, polys, dangles, cuts, invalids = Polygonizer.polygonize(C, curves, extend=0.0, resolution=12)


geoms = []
coords = Io.curves_to_geoms(C.selected_objects, 12, geoms)
ls = geoms[0]


b = ls.buffer(0.05, join_style=JOIN_STYLE.mitre, cap_style=CAP_STYLE.flat, single_sided=False)
Io.to_curve(C.scene, coords, b, name="buffer")




import time
import logging
logger = logging.getLogger("archipack")

from archipack.pygeos.op_polygonize import PolygonizeOp
from archipack.pygeos.geom import GeometryFactory
from archipack.pygeos.prepared import PreparedGeometryFactory
from archipack.pygeos.op_linemerge import LineMerger
from archipack.pygeos.op_overlay import OverlayOp, SnapOverlayOp
from archipack.pygeos.op_union import UnaryUnionOp
from archipack.pygeos.op_buffer import BufferOp
from archipack.archipack_polylines import Io, CoordSys, Polygonizer
from archipack.pygeos.shared import JOIN_STYLE, CAP_STYLE

geoms = []
coords = Io.curves_to_geoms(C.selected_objects, 12, geoms)

res = SnapOverlayOp.intersection(geoms[1], geoms[0])
Io.to_curve(C.scene, coords, res, name="res")

ls = geoms[0]
gf = GeometryFactory()
lr = gf.createLinearRing(ls.coords)
p = gf.createPolygon(lr)
b = p.buffer(0)

Io.to_curve(C.scene, coords, b, name="buffer")


walls = C.selected_objects
coords = CoordSys(walls)

gf = GeometryFactory()
gf.outputFactory = Io(scene=C.scene, coordsys=coords)

coordsys, poly1, dangles, cuts, invalids = Polygonizer.polygonize(C, walls[0], extend=0.0, resolution=12)
coordsys, poly2, dangles, cuts, invalids = Polygonizer.polygonize(C, walls[1], extend=0.0, resolution=12)

polys = poly1 + poly2

polys[0].intersects(polys[1])
polys[0].disjoint(polys[1])
polys[0].touches(polys[1])
polys[0].overlaps(polys[1])
polys[0].contains(polys[1])
polys[0].within(polys[1])
polys[0].crosses(polys[1])

prep = PreparedGeometryFactory.prepare(poly1[0])
prep.touches(poly2[0])

res = UnaryUnionOp.union(polys)
gf.output(res, name="cascaded_union", multiple=True)

opt = res.simplify(0.001, False)
gf.output(opt, name="simple", multiple=True)

buf = poly1[0].buffer(0.1, cap_style=1, join_style=2)
gf.output(buf, name="buffer", multiple=True)

hull = poly1[0].convex_hull
gf.output(hull, name="hull", multiple=True)

res = OverlayOp.overlayOp(poly1[0], poly2[0], OverlayOp.opDIFFERENCE)
gf.output(res, name="difference", multiple=True)

res = OverlayOp.overlayOp(poly1[0], poly2[0], OverlayOp.opUNION)
gf.output(res, name="union", multiple=True)

res = OverlayOp.overlayOp(poly1[0], poly2[0], OverlayOp.opINTERSECTION)
gf.output(res, name="intersection", multiple=True)

res = OverlayOp.overlayOp(poly1[0], poly2[0], OverlayOp.opSYMDIFFERENCE)
gf.output(res, name="sym_diff", multiple=True)





import time
import logging
logger = logging.getLogger("archipack")

walls = C.selected_objects
from archipack.archipack_polylines import Polygons
from archipack.pygeos.op_polygonize import PolygonizeOp
from archipack.pygeos.geom import GeometryFactory
from archipack.pygeos.op_linemerge import LineMerger
from archipack.pygeos.op_overlay import OverlayOp
from archipack.archipack_polylib import Io, OutputFactory

t = time.time()
curves = C.selected_objects
extend = 0.02

Polygons.polygonize(C, curves, extend=extend, resolution=12)


p = Polygons(coords)
p.from_curves(walls, 12)
p._split_segs(extend=extend)

t2 = time.time()
lines = [gf.createLineString([seg.c0, seg.c1]) for seg in p.Q_segs._geoms if seg.available and seg.c0 is not seg.c1]
logger.debug("Polygonize input lines :%.4f seconds" , (time.time() - t2))

t2 = time.time()
merged = LineMerger.merge(lines)
logger.debug("LineMerger :%.4f seconds" , (time.time() - t2))

t2 = time.time()
polys, dangles, cuts, invalids = PolygonizeOp.polygonize_full(merged)
logger.debug("PolygonizeOp.polygonize_full(lines) :%.4f seconds" , (time.time() - t2))


logger.debug("Polygonize :%.4f seconds" , (time.time() - t))
gf.output(polys, name="polys", multiple=True)



bpy.app.debug = True
C.object.data.show_extra_indices = True


p.Q_segs._geoms[1413]._intersect_seg(p.Q_segs._geoms[4455])




seg_index = 103
sel = p.Q_segs.intersects_ext(p.Q_segs._geoms[seg_index], extend)
sel = (0, sorted([i for i in sel[1] if p.Q_segs._geoms[i].available]))

def select_edges(context, segs):
    import bmesh
    obj = context.edit_object
    me = obj.data
    bm = bmesh.from_edit_mesh(me)
    for e in bm.edges:
        e.select = False
    for e in segs:
        bm.edges[e].select = True
    bmesh.update_edit_mesh(me, False)


select_edges(C, [11])




p.polygonize()

for i, shape in enumerate(p.shapes):
    gen = True
    for j, pt in enumerate(p.inside_wall):
        if shape.inside(pt):
            gen = False
            break
    if gen:
        curves = []
        shape.as_curves(C, curves, p.coordsys.world, True, True)
        boundary = curves.pop(0)
        boundary.select = True
        C.scene.objects.active = boundary
        bpy.ops.archipack.floor_from_curve(auto_manipulate=False)
        floor_name = C.active_object.name
        for hole in curves:
            hole.name = "Hole_{}_{}_{}_".format(i, shape.cw, shape.depth)
            bpy.ops.archipack.floor_cutter(parent=floor_name, curve=hole.name, auto_manipulate=False)
        bpy.ops.object.select_all(action='DESELECT')
        boundary.select = True
        C.scene.objects.active = boundary
        for hole in curves:
            hole.select = True
        bpy.ops.object.delete(use_global=False)


walls = C.selected_objects
from archipack.archipack_polylines import Polygons, CoordSys, Shape

extend = 0.01
coords = CoordSys(walls)
p = Polygons(coords, extend=extend)
p.from_curves(walls, 12)
p.cascaded_union(extend=extend)


p.Q_segs.intersects_ext(p.Q_segs._geoms[1413], 0.2)

p._reversed_copy()
p._neighboors()

p.polygonize()
p._relationship()


p.as_curves(C, shapes=p.err)

for i, shape in enumerate(p.shapes):
    c = shape.as_curve(C, p.coordsys.world, True, True)
    c.name = "Curve_{}_{}_{}_".format(i, shape.cw, shape.depth)


p.as_curves(C)

"""
