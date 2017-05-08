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
import bpy
import bgl
import blf
from math import sin, cos, atan2, pi
from mathutils import Vector, Matrix
from mathutils.geometry import intersect_line_plane, intersect_point_line, intersect_line_sphere
from bpy_extras import view3d_utils
from bpy.types import PropertyGroup
from bpy.props import FloatVectorProperty, StringProperty, CollectionProperty, BoolProperty
from bpy.app.handlers import persistent
from .archipack_utils import operator_exists

try:
    from np_station.np_point_move import snap_point
    HAS_NP_STATION = True
except:
    HAS_NP_STATION = False
    pass


# Arrow sizes (world units)
arrow_size = 0.1
# Handle area size (pixels)
handle_size = 10

# ------------------------------------------------------------------
# Define Gl Handle types
# ------------------------------------------------------------------


class Gl():
    """
        handle 3d -> 2d gl drawing
    """
    def __init__(self):
        self.width = 1
        self.pos_2d = Vector((0, 0))
        self.colour_active = (1.0, 0.0, 0.0, 1.0)
        self.colour_hover = (1.0, 1.0, 0.0, 1.0)
        self.colour_normal = (1.0, 1.0, 1.0, 1.0)
        self.colour_inactive = (0.0, 0.0, 0.0, 1.0)

    @property
    def colour(self):
        return self.colour_inactive

    def position_2d_from_coord(self, context, coord):
        """ coord given in local input coordsys
        """
        region = context.region
        rv3d = context.region_data
        loc = view3d_utils.location_3d_to_region_2d(region, rv3d, coord, self.pos_2d)
        return loc

    def _end(self):
        bgl.glEnd()
        bgl.glPopAttrib()
        bgl.glLineWidth(1)
        bgl.glDisable(bgl.GL_BLEND)
        bgl.glColor4f(0.0, 0.0, 0.0, 1.0)

    def _start_poly(self, colour):
        bgl.glPushAttrib(bgl.GL_ENABLE_BIT)
        bgl.glEnable(bgl.GL_BLEND)
        bgl.glColor4f(*colour)
        bgl.glBegin(bgl.GL_POLYGON)

    def _start_line(self, colour, width=1):
        bgl.glPushAttrib(bgl.GL_ENABLE_BIT)
        bgl.glEnable(bgl.GL_BLEND)
        bgl.glEnable(bgl.GL_LINE)
        bgl.glColor4f(*colour)
        bgl.glLineWidth(width)
        bgl.glBegin(bgl.GL_LINE_STRIP)

    def draw_text(self, text, x, y, angle, font_height, colour):
        # dirty fast assignment
        dpi, font_id = 72, 0
        bgl.glColor4f(*colour)
        blf.position(font_id, x, y, 0)
        blf.rotation(font_id, angle)
        blf.size(font_id, font_height, dpi)
        blf.draw(font_id, text)

    def draw(self, context):
        gl_type = type(self).__name__
        if 'Handle' in gl_type or gl_type in ['GlPolygon']:
            self._start_poly(self.colour)
        elif gl_type in ['GlLine', 'GlArc', 'GlPolyline']:
            self._start_line(self.colour, self.width)
        if gl_type == 'GlText':
            x, y = self.position_2d_from_coord(context, self.pts[0])
            self.draw_text(self.txt, x, y, self.angle, self.font_height, self.colour)
        else:
            for pt in self.pts:
                x, y = self.position_2d_from_coord(context, pt)
                bgl.glVertex2f(x, y)
            self._end()


class GlText(Gl):

    def __init__(self, round=2, label='', z_axis=Vector((0, 0, 1))):
        self.z_axis = z_axis
        self.value = 0
        self.round = round
        self.label = label
        self.font_height = 16
        Gl.__init__(self)

    @property
    def angle(self):
        return 0

    @property
    def pts(self):
        return [self.pos_3d]

    @property
    def txt(self):
        return self.label + str(round(self.value, self.round))

    def set_pos(self, context, value, pos_3d, direction, normal=Vector((0, 0, 1))):
        self.up_axis = direction.normalized()
        self.c_axis = self.up_axis.cross(normal)
        self.pos_3d = pos_3d
        self.value = value


class GlLine(Gl):

    def __init__(self, z_axis=Vector((0, 0, 1))):
        self.z_axis = z_axis
        self.p = Vector((0, 0, 0))
        self.v = Vector((0, 0, 0))
        Gl.__init__(self)

    @property
    def length(self):
        return self.v.length

    @property
    def angle(self):
        return atan2(self.v.y, self.v.x)

    @property
    def cross(self):
        return self.v.cross(self.z_axis)

    def normal(self, t=0):
        # perpendiculaire a droite du segment
        n = GlLine()
        n.p = self.lerp(t)
        n.v = self.cross
        return n

    def sized_normal(self, t, size):
        n = GlLine()
        n.p = self.lerp(t)
        n.v = size * self.cross.normalized()
        return n

    def lerp(self, t):
        return self.p + self.v * t

    def offset(self, offset):
        """
            offset > 0 on the right part
        """
        self.p += offset * self.cross.normalized()

    @property
    def pts(self):
        p0 = self.p
        p1 = self.p + self.v
        return [p0, p1]


class GlCircle(Gl):

    def __init__(self):
        self.r = 0
        self.c = Vector((0, 0, 0))
        Gl.__init__(self)


class GlArc(GlCircle):

    def __init__(self, z_axis=Vector((0, 0, 1))):
        """
            a0 and da arguments are in radians
            a0 = 0   on the right side
            a0 = pi on the left side
            da > 0 CCW contrary-clockwise
            da < 0 CW  clockwise
            stored internally as radians
        """
        GlCircle.__init__(self)
        if z_axis.z < 1:
            x_axis = z_axis.cross(Vector((0, 0, 1)))
            y_axis = x_axis.cross(z_axis)
        else:
            x_axis = Vector((1, 0, 0))
            y_axis = Vector((0, 1, 0))
        self.rM = Matrix([
            x_axis,
            y_axis,
            z_axis
        ])
        self.z_axis = z_axis
        self.a0 = 0
        self.da = 0

    @property
    def length(self):
        return self.r * abs(self.da)

    def normal(self, t=0):
        """
            always on the right side
        """
        n = GlLine(z_axis=self.z_axis)
        n.p = self.lerp(t)
        if self.da < 0:
            n.v = self.c - n.p
        else:
            n.v = n.p - self.c
        return n

    def sized_normal(self, t, size):
        n = GlLine(z_axis=self.z_axis)
        n.p = self.lerp(t)
        if self.da < 0:
            n.v = size * (self.c - n.p).normalized()
        else:
            n.v = size * (n.p - self.c).normalized()
        return n

    def lerp(self, t):
        a = self.a0 + t * self.da
        return self.c + self.rM * Vector((self.r * cos(a), self.r * sin(a), 0))

    def tangeant(self, t, length):
        a = self.a0 + t * self.da
        ca = cos(a)
        sa = sin(a)
        n = GlLine()
        n.p = self.c + self.rM * Vector((self.r * ca, self.r * sa, 0))
        n.v = self.rM * Vector((length * sa, -length * ca, 0))
        if self.da > 0:
            n.v = -n.v
        return n

    def offset(self, offset):
        """
            offset > 0 on the right part
        """
        if self.da > 0:
            radius = self.r + offset
        else:
            radius = self.r - offset
        return GlArc(self.c, radius, self.a0, self.da, z_axis=self.z_axis)

    @property
    def pts(self):
        n_pts = max(1, int(round(abs(self.da) / pi * 30, 0)))
        t_step = 1 / n_pts
        return [self.lerp(i * t_step) for i in range(n_pts + 1)]


class GlPolygon(Gl):
    def __init__(self, colour):
        self._colour = colour
        self.pts_3d = []
        Gl.__init__(self)

    def set_pos(self, pts_3d):
        self.pts_3d = pts_3d

    @property
    def colour(self):
        return self._colour

    @property
    def pts(self):
        return self.pts_3d


class GlPolyline(Gl):
    def __init__(self, colour):
        self._colour = colour
        self.pts_3d = []
        Gl.__init__(self)

    def set_pos(self, pts_3d):
        self.pts_3d = pts_3d
        # self.pts_3d.append(pts_3d[0])

    @property
    def colour(self):
        return self._colour

    @property
    def pts(self):
        return self.pts_3d


class GlHandle(Gl):

    def __init__(self, sensor_size, size, selectable=False):
        """
            sensor_size : 2d size in pixels of sensor area
            size : 3d size of handle
        """
        self.size = size
        self.sensor_size = sensor_size
        self.pos_3d = Vector((0, 0, 0))
        self.up_axis = Vector((0, 0, 0))
        self.c_axis = Vector((0, 0, 0))
        self.hover = False
        self.active = False
        self.selectable = selectable
        Gl.__init__(self)

    def set_pos(self, context, pos_3d, direction, normal=Vector((0, 0, 1))):
        self.up_axis = direction.normalized()
        self.c_axis = self.up_axis.cross(normal)
        self.pos_3d = pos_3d
        self.pos_2d = self.position_2d_from_coord(context, pos_3d)

    def check_hover(self, pos_2d):
        dp = pos_2d - self.pos_2d
        self.hover = abs(dp.x) < self.sensor_size and abs(dp.y) < self.sensor_size

    @property
    def pts(self):
        raise NotImplementedError

    @property
    def colour(self):
        if self.selectable:
            if self.active:
                return self.colour_active
            elif self.hover:
                return self.colour_hover
            return self.colour_normal
        else:
            return self.colour_inactive


class SquareHandle(GlHandle):

    def __init__(self, sensor_size, size, selectable=False):
        GlHandle.__init__(self, sensor_size, size, selectable)

    @property
    def pts(self):
        n = self.up_axis
        c = self.c_axis
        x = n * self.size / 2
        y = c * self.size / 2
        return [self.pos_3d - x - y, self.pos_3d + x - y, self.pos_3d + x + y, self.pos_3d - x + y]


class TriHandle(GlHandle):

    def __init__(self, sensor_size, size, selectable=False):
        GlHandle.__init__(self, sensor_size, size, selectable)

    @property
    def pts(self):
        n = self.up_axis
        c = self.c_axis
        x = n * self.size
        y = c * self.size / 2
        return [self.pos_3d - x + y, self.pos_3d - x - y, self.pos_3d]


# ------------------------------------------------------------------
# Define Manipulators
# ------------------------------------------------------------------


class Manipulator():

    def __init__(self, context, o, datablock, manipulator):
        """
            o : object to manipulate
            datablock : object data to manipulate
            manipulator: object archipack_manipulator datablock
        """
        self.pts_mode = 'SIZE'
        self.o = o
        self.datablock = datablock
        self.manipulator = manipulator
        self.origin = Vector((0, 0, 1))
        self.mouse_pos = Vector((0, 0))
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_callback, args, 'WINDOW', 'POST_PIXEL')

    @classmethod
    def poll(cls, context):
        return True

    def exit(self):
        # print("Manipulator.exit() %s" % (type(self).__name__))
        if self._handle is not None:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
        self._handle = None

    def press(self, context, event):
        """
            Manipulators must implement
            mouse press event handler
            return True to callback manipulable_manipulate
        """
        raise NotImplementedError

    def release(self, context, event):
        """
            Manipulators must implement
            mouse release event handler
            return False to callback manipulable_release
        """
        raise NotImplementedError

    def mouse_move(self, context, event):
        """
            Manipulators must implement
            mouse move event handler
            return True to callback manipulable_manipulate
        """
        raise NotImplementedError

    def timer(self, context, event):
        return False

    def modal(self, context, event):
        if event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            return self.press(context, event)
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            return self.release(context, event)
        elif event.type == 'MOUSEMOVE':
            return self.mouse_move(context, event)
        return False

    def mouse_position(self, event):
        self.mouse_pos.x, self.mouse_pos.y = event.mouse_region_x, event.mouse_region_y

    def get_pos3d(self, context):
        """
            convert mouse pos to 3d point over plane defined by origin and normal
        """
        region = context.region
        rv3d = context.region_data
        rM = context.active_object.matrix_world.to_3x3()
        view_vector_mouse = view3d_utils.region_2d_to_vector_3d(region, rv3d, self.mouse_pos)
        ray_origin_mouse = view3d_utils.region_2d_to_origin_3d(region, rv3d, self.mouse_pos)
        pt = intersect_line_plane(ray_origin_mouse, ray_origin_mouse + view_vector_mouse,
            self.origin, rM * self.manipulator.normal, False)
        # fix issue with parallel plane
        if pt is None:
            pt = intersect_line_plane(ray_origin_mouse, ray_origin_mouse + view_vector_mouse,
                self.origin, view_vector_mouse, False)
        return pt

    def get_value(self, data, attr, index=-1):
        try:
            if index > -1:
                return getattr(data, attr)[index]
            else:
                return getattr(data, attr)
        except:
            return 0

    def set_value(self, context, data, attr, value, index=-1):
        try:
            if self.get_value(data, attr, index) != value:
                # switch context so unselected object may be manipulable too
                old = context.active_object
                state = self.o.select
                self.o.select = True
                context.scene.objects.active = self.o
                if index > -1:
                    getattr(data, attr)[index] = value
                else:
                    setattr(data, attr, value)
                self.o.select = state
                old.select = True
                context.scene.objects.active = old
        except:
            pass

    def preTranslate(self, tM, vec):
        """
            return a preTranslated Matrix
            tM Matrix source
            vec Vector translation
        """
        return tM * Matrix([
        [1, 0, 0, vec.x],
        [0, 1, 0, vec.y],
        [0, 0, 1, vec.z],
        [0, 0, 0, 1]])

    def _move(self, o, axis, value):
        if axis == 'x':
            tM = self.preTranslate(o.matrix_world, Vector((value, 0, 0)))
        elif axis == 'y':
            tM = self.preTranslate(o.matrix_world, Vector((0, value, 0)))
        else:
            tM = self.preTranslate(o.matrix_world, Vector((0, 0, value)))
        o.matrix_world = tM

    def move_linked(self, context, axis, value):
        """
            Move an object along local axis
            takes care of linked too, fix issue #8
        """
        old = context.active_object
        bpy.ops.object.select_all(action='DESELECT')
        self.o.select = True
        context.scene.objects.active = self.o
        bpy.ops.object.select_linked(type='OBDATA')
        for o in context.selected_objects:
            if o != self.o:
                self._move(o, axis, value)
        bpy.ops.object.select_all(action='DESELECT')
        old.select = True
        context.scene.objects.active = old

    def move(self, context, axis, value):
        """
            Move an object along local axis
        """
        self._move(self.o, axis, value)


class SnapPointManipulator(Manipulator):
    """
        np_station based snap manipulator
        dosent update anything by itself.

        Use NP020PM.flag and NP020PM.helper.location
        to update in manipulable_manipulate

        Use prop2_name as string identifier
        Use p1 and p2 as number identifiers

    """
    def __init__(self, context, o, datablock, manipulator, handle_size):
        self.handle = SquareHandle(handle_size, 1.2 * arrow_size, selectable=True)
        Manipulator.__init__(self, context, o, datablock, manipulator)

    def check_hover(self):
        self.handle.check_hover(self.mouse_pos)

    def press(self, context, event):
        if self.handle.hover:
            # Start a np_station np_point_move snap tool
            self.handle.hover = False
            self.handle.active = True
            if HAS_NP_STATION:
                self.o.select = True
                takeloc = self.o.matrix_world * self.manipulator.p0
                print("Invoke np_point_move %s" % (takeloc))
                snap_point().invoke(takeloc=takeloc)
            # return True
        return False

    def release(self, context, event):
        # Release event should never be called
        # as events are captured by NP020PM modal
        self.check_hover()
        self.handle.active = False
        # Keep manipulator active until np_station exit
        # np_snap = snap_point()
        print("release should not occur")
        # self.origin = self.pos_3d
        # False to callback manipulable_release
        return False

    def update(self, context, event):
        # NOTE:
        # dosent set anything internally
        return
        # if MP020PM.flag is 'PLACE':
        # NP020PM.helper.location
        # pt = self.get_pos3d(context)
        # self.set_value(context, self.datablock, self.manipulator.prop1_name, pt)

    def mouse_move(self, context, event):
        """
            NOTE:
                mouse events are captured by NP020PM modal
                so we cant expect any direct response here
                while running, and thus we cant exepct a release event
                so wait for modal end, and on first move_event after
                disable handle and callback

        """
        self.mouse_position(event)
        if self.handle.active:
            np_snap = snap_point()
            print("mouse_move, np_snap.is_running:%s" % np_snap.is_running)
            # self.handle.active = np_snap.is_running
            # self.update(context)
            #
            # True here to callback manipulable_manipulate
            self.handle.active = np_snap.is_running
            return not self.handle.active
        else:
            self.check_hover()
        return False

    def draw_callback(self, _self, context):
        left, right, side, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.handle.set_pos(context, left, Vector((1, 0, 0)), normal=normal)
        self.handle.draw(context)


#
wall_pts = []


class WallSnapManipulator(Manipulator):
    """
        np_station based snap manipulator

        Use prop2_name as string identifier
        Use p1 and p2 as number identifiers

    """
    def __init__(self, context, o, datablock, manipulator, handle_size):
        self.np_snap = None
        self.wall_part1 = GlPolygon((0.5, 0, 0, 0.2))
        self.wall_line1 = GlPolyline((0.5, 0, 0, 0.8))
        self.wall_part2 = GlPolygon((0.5, 0, 0, 0.2))
        self.wall_line2 = GlPolyline((0.5, 0, 0, 0.8))
        # self.draw_placeholders = False
        # self._timer = None
        self.handle = SquareHandle(handle_size, 1.2 * arrow_size, selectable=True)
        Manipulator.__init__(self, context, o, datablock, manipulator)
        # disable own draw handler since np_snap handle it by it own
        # bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
        # self._handle = None

    @classmethod
    def poll(cls, context):
        return HAS_NP_STATION and operator_exists("OBJECT_OT_np020_point_move")

    def check_hover(self):
        self.handle.check_hover(self.mouse_pos)

    def press(self, context, event):
        global wall_pts
        if self.handle.hover:
            # Start a np_station np_point_move snap tool
            self.handle.hover = False
            self.handle.active = True
            if HAS_NP_STATION:
                self.o.select = True
                wall = self.o.data.archipack_wall2[0]
                wall_idx = int(self.manipulator.prop1_name)
                # init placeholders 3d absolute world positionned points
                # between 2 and 3 points, the 2nd being moved by np_snap.delta
                takeloc, p1, side, normal = self.manipulator.get_pts(self.o.matrix_world)
                wall_pts = [None, takeloc, p1]
                if wall_idx > 0:
                    p0, right, side, normal = wall.parts[wall_idx - 1].manipulators[2].get_pts(self.o.matrix_world)
                    wall_pts[0] = p0
                # enable placeholders drawing
                self.draw_placeholders = True
                print("wall_pts : %s" % wall_pts)
                # Invoke snap tool
                print("Invoke np_point_move %s" % (takeloc))
                self.np_snap = snap_point()
                # self._timer = context.window_manager.event_timer_add(0.1, context.window)
                self.np_snap.invoke(takeloc=takeloc, callback=self.np_callback, draw_callback=self.np_draw,
                    constrain=True, normal=Vector((0, 0, 1)))
            return True
        return False

    def release(self, context, event):
        # Release event should never be called
        # as events are captured by NP020PM modal
        self.check_hover()
        self.handle.active = False
        # Keep manipulator active until np_station exit
        # np_snap = snap_point()
        print("release should not occur")
        # self.origin = self.pos_3d
        # False to callback manipulable_release
        return False

    def np_callback(self, context, event, state):
        """
            np station callback on moving, place, or cancel
        """
        global wall_pts

        np_snap = snap_point()

        print("np_callback %s" % (state))
        if state == 'RUNNING':
            print("event.type: %s" % (event.type))
            # event here is mouse move
            # update gl placeholders location
            print("np_callback wall_pts : %s" % wall_pts)
            print("np_callback self: %s" % (type(self).__name__))

        else:
            # kill gl placeholder
            self.handle.active = False

            if state == 'SUCCESS':

                self.o.select = True
                # apply changes to wall
                wall = self.o.data.archipack_wall2[0]
                wall.auto_update = False

                g = wall.get_generator()

                # rotation relative to object
                rM = self.o.matrix_world.inverted().to_3x3()
                wall_idx = int(self.manipulator.prop1_name)

                # new point in object space
                pt = g.walls[wall_idx].lerp(0) + (rM * np_snap.delta).to_2d()
                da = 0

                # adjust size and rotation of segment before current
                if wall_idx > 0:
                    w = g.walls[wall_idx - 1]
                    part = wall.parts[wall_idx - 1]
                    dp = (pt - w.lerp(0))
                    part.length = dp.length
                    da = atan2(dp.y, dp.x) - w.straight(1).angle
                    a0 = part.a0 + da
                    if a0 > pi:
                        a0 -= 2 * pi
                    if a0 < -pi:
                        a0 += 2 * pi
                    print("a0:%.4f part.a0:%.4f da:%.4f" % (a0, part.a0, da))
                    part.a0 = a0

                # adjust length of current segment
                w = g.walls[wall_idx]
                part = wall.parts[wall_idx]
                p0 = w.lerp(1)
                dp = (p0 - pt)
                print("pt:%s delta:%s dp:%s idx:%s" % (pt, np_snap.delta, dp, wall_idx))
                part.length = dp.length

                da1 = atan2(dp.y, dp.x) - w.straight(1).angle
                a0 = part.a0 + da1 - da
                if a0 > pi:
                    a0 -= 2 * pi
                if a0 < -pi:
                    a0 += 2 * pi
                print("a0:%.4f part.a0:%.4f da:%.4f" % (a0, part.a0, da))
                part.a0 = a0

                # move object when point 0
                if wall_idx == 0:
                    self.o.location += np_snap.delta

                # adjust rotation on both sides of next segment
                if wall_idx + 1 < wall.n_parts:

                    part = wall.parts[wall_idx + 1]
                    a0 = part.a0 - da1
                    if a0 > pi:
                        a0 -= 2 * pi
                    if a0 < -pi:
                        a0 += 2 * pi
                    part.a0 = a0
                wall.auto_update = True
                wall.update(context)

        return

    def np_draw(self, context):
        # draw wall placeholders
        global wall_pts
        np_snap = snap_point()
        if self.o is None:
            return
        z = self.o.data.archipack_wall2[0].z
        p0 = wall_pts[1] + np_snap.delta
        p1 = wall_pts[2]
        self.wall_part1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.wall_line1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.wall_part1.draw(context)
        self.wall_line1.draw(context)
        if wall_pts[0] is not None:
            p0, p1 = wall_pts[0], p0
            self.wall_part2.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
            self.wall_line2.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
            self.wall_part2.draw(context)
            self.wall_line2.draw(context)

    def mouse_move(self, context, event):
        """
            NOTE:
                mouse events are captured by NP020PM modal
                so we cant expect any direct response here
                while running, and thus we cant exepct a release event
                so wait for modal end, and on first move_event after
                disable handle and callback

        """
        self.mouse_position(event)
        if self.handle.active:
            # False here to pass_through
            # print("i'm able to pick up mouse move event while transform running")
            return False
        else:
            self.check_hover()
        return False

    def draw_callback(self, _self, context):
        left, right, side, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.handle.set_pos(context, left, (left - right).normalized(), normal=normal)
        self.handle.draw(context)


#
fence_pts = []


class FenceSnapManipulator(Manipulator):
    """
        np_station based snap manipulator

        Use prop2_name as string identifier
        Use p1 and p2 as number identifiers

    """
    def __init__(self, context, o, datablock, manipulator, handle_size):
        self.np_snap = None
        self.fence_part1 = GlPolygon((0.5, 0, 0, 0.2))
        self.fence_line1 = GlPolyline((0.5, 0, 0, 0.8))
        self.fence_part2 = GlPolygon((0.5, 0, 0, 0.2))
        self.fence_line2 = GlPolyline((0.5, 0, 0, 0.8))
        # self.draw_placeholders = False
        # self._timer = None
        self.handle = SquareHandle(handle_size, 1.2 * arrow_size, selectable=True)
        Manipulator.__init__(self, context, o, datablock, manipulator)
        # disable own draw handler since np_snap handle it by it own
        # bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
        # self._handle = None

    @classmethod
    def poll(cls, context):
        return HAS_NP_STATION and operator_exists("OBJECT_OT_np020_point_move")

    def check_hover(self):
        self.handle.check_hover(self.mouse_pos)

    def press(self, context, event):
        global fence_pts
        if self.handle.hover:
            # Start a np_station np_point_move snap tool
            self.handle.hover = False
            self.handle.active = True
            if HAS_NP_STATION:
                self.o.select = True
                fence = self.o.data.archipack_fence[0]
                fence_idx = int(self.manipulator.prop1_name)
                # init placeholders 3d absolute world positionned points
                # between 2 and 3 points, the 2nd being moved by np_snap.delta
                takeloc, p1, side, normal = self.manipulator.get_pts(self.o.matrix_world)
                fence_pts = [None, takeloc, p1]
                if fence_idx > 0:
                    p0, right, side, normal = fence.parts[fence_idx - 1].manipulators[2].get_pts(self.o.matrix_world)
                    fence_pts[0] = p0
                # enable placeholders drawing
                self.draw_placeholders = True
                print("fence_pts : %s" % fence_pts)
                # Invoke snap tool
                print("Invoke np_point_move %s" % (takeloc))
                self.np_snap = snap_point()
                # self._timer = context.window_manager.event_timer_add(0.1, context.window)
                self.np_snap.invoke(takeloc=takeloc, callback=self.np_callback, draw_callback=self.np_draw,
                    constrain=True, normal=Vector((0, 0, 1)))
            return True
        return False

    def release(self, context, event):
        # Release event should never be called
        # as events are captured by NP020PM modal
        self.check_hover()
        self.handle.active = False
        # Keep manipulator active until np_station exit
        # np_snap = snap_point()
        print("release should not occur")
        # self.origin = self.pos_3d
        # False to callback manipulable_release
        return False

    def np_callback(self, context, event, state):
        """
            np station callback on moving, place, or cancel
        """
        global fence_pts

        np_snap = snap_point()

        print("np_callback %s" % (state))
        if state == 'RUNNING':
            print("event.type: %s" % (event.type))
            # event here is mouse move
            # update gl placeholders location
            print("np_callback fence_pts : %s" % fence_pts)
            print("np_callback self: %s" % (type(self).__name__))

        else:
            # kill gl placeholder
            self.handle.active = False

            if state == 'SUCCESS':

                self.o.select = True
                # apply changes to fence
                fence = self.o.data.archipack_fence[0]
                fence.auto_update = False

                g = fence.get_generator()

                # rotation relative to object
                rM = self.o.matrix_world.inverted().to_3x3()
                fence_idx = int(self.manipulator.prop1_name)

                # new point in object space
                pt = g.fences[fence_idx].lerp(0) + (rM * np_snap.delta).to_2d()
                da = 0

                # adjust size and rotation of segment before current
                if fence_idx > 0:
                    w = g.fences[fence_idx - 1]
                    part = fence.parts[fence_idx - 1]
                    dp = (pt - w.lerp(0))
                    part.length = dp.length
                    da = atan2(dp.y, dp.x) - w.straight(1).angle
                    a0 = part.a0 + da
                    if a0 > pi:
                        a0 -= 2 * pi
                    if a0 < -pi:
                        a0 += 2 * pi
                    print("a0:%.4f part.a0:%.4f da:%.4f" % (a0, part.a0, da))
                    part.a0 = a0

                # adjust length of current segment
                w = g.fences[fence_idx]
                part = fence.parts[fence_idx]
                p0 = w.lerp(1)
                dp = (p0 - pt)
                print("pt:%s delta:%s dp:%s idx:%s" % (pt, np_snap.delta, dp, fence_idx))
                part.length = dp.length

                da1 = atan2(dp.y, dp.x) - w.straight(1).angle
                a0 = part.a0 + da1 - da
                if a0 > pi:
                    a0 -= 2 * pi
                if a0 < -pi:
                    a0 += 2 * pi
                print("a0:%.4f part.a0:%.4f da:%.4f" % (a0, part.a0, da))
                part.a0 = a0

                # move object when point 0
                if fence_idx == 0:
                    self.o.location += np_snap.delta

                # adjust rotation on both sides of next segment
                if fence_idx + 1 < fence.n_parts:

                    part = fence.parts[fence_idx + 1]
                    a0 = part.a0 - da1
                    if a0 > pi:
                        a0 -= 2 * pi
                    if a0 < -pi:
                        a0 += 2 * pi
                    part.a0 = a0
                fence.auto_update = True
                fence.update(context)

        return

    def np_draw(self, context):
        # draw fence placeholders
        global fence_pts
        np_snap = snap_point()
        if self.o is None:
            return
        z = self.o.data.archipack_fence[0].post_z
        p0 = fence_pts[1] + np_snap.delta
        p1 = fence_pts[2]
        self.fence_part1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.fence_line1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.fence_part1.draw(context)
        self.fence_line1.draw(context)
        if fence_pts[0] is not None:
            p0, p1 = fence_pts[0], p0
            self.fence_part2.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
            self.fence_line2.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
            self.fence_part2.draw(context)
            self.fence_line2.draw(context)

    def mouse_move(self, context, event):
        """
            NOTE:
                mouse events are captured by NP020PM modal
                so we cant expect any direct response here
                while running, and thus we cant exepct a release event
                so wait for modal end, and on first move_event after
                disable handle and callback

        """
        self.mouse_position(event)
        if self.handle.active:
            # False here to pass_through
            # print("i'm able to pick up mouse move event while transform running")
            return False
        else:
            self.check_hover()
        return False

    def draw_callback(self, _self, context):
        left, right, side, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.handle.set_pos(context, left, (left - right).normalized(), normal=normal)
        self.handle.draw(context)


"""
class PointManipulator(Manipulator):
    def __init__(self, context, pos_3d, handle_size):
        self.handle = Handle(handle_size)

        self.pos_3d = pos_3d
        self.origin = pos_3d
        Manipulator.__init__(self, context)
    def check_hover(self):
        self.handle.check_hover(self.mouse_pos)
    def press(self, context, event):
        if self.handle.hover:
            self.handle.active = True
            return True
        return False
    def release(self, context, event):
        self.check_hover()
        self.handle.active = False
        self.origin = self.pos_3d
        return False
    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle.active:
            self.update_position_3d(context)
            return True
        else:
            self.check_hover()
        return False
    def update_position_3d(self, context):
        self.pos_3d = self.get_pos3d(context)
        self.delta = self.pos_3d-self.origin
    def draw_callback(self, _self, context):

        self.handle.draw(context, self.pos_3d, size, Vector((0,1,0)))
"""


class CounterManipulator(Manipulator):

    """
        increase or decrease an integer step by step on click
    """
    def __init__(self, context, o, datablock, manipulator, handle_size):
        self.handle_left = TriHandle(handle_size, arrow_size, selectable=True)
        self.handle_right = TriHandle(handle_size, arrow_size, selectable=True)
        self.line_0 = GlLine()
        self.label = GlText()
        Manipulator.__init__(self, context, o, datablock, manipulator)

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)
        self.handle_left.check_hover(self.mouse_pos)

    def press(self, context, event):
        if self.handle_right.hover:
            value = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.set_value(context, self.datablock, self.manipulator.prop1_name, value + 1)
            self.handle_right.active = True
            return True
        if self.handle_left.hover:
            value = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.set_value(context, self.datablock, self.manipulator.prop1_name, value - 1)
            self.handle_left.active = True
            return True
        return False

    def release(self, context, event):
        self.check_hover()
        self.handle_right.active = False
        self.handle_left.active = False
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active:
            return True
        if self.handle_left.active:
            return True
        else:
            self.check_hover()
        return False

    def draw_callback(self, _self, context):
        """
            draw on screen feedback using gl.
        """
        left, right, side, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.origin = left
        self.line_0.p = left
        self.line_0.v = right - left
        self.line_0.z_axis = normal
        self.label.z_axis = normal
        value = self.get_value(self.datablock, self.manipulator.prop1_name)
        self.handle_left.set_pos(context, self.line_0.p, -self.line_0.v, normal=normal)
        self.handle_right.set_pos(context, self.line_0.lerp(1), self.line_0.v, normal=normal)
        self.label.set_pos(context, value, self.line_0.lerp(0.5), self.line_0.v, normal=normal)
        self.label.draw(context)
        self.handle_left.draw(context)
        self.handle_right.draw(context)


class SizeManipulator(Manipulator):

    def __init__(self, context, o, datablock, manipulator, handle_size):
        self.handle_left = TriHandle(handle_size, arrow_size)
        self.handle_right = TriHandle(handle_size, arrow_size, selectable=True)
        self.line_0 = GlLine()
        self.line_1 = GlLine()
        self.line_2 = GlLine()
        self.label = GlText()
        Manipulator.__init__(self, context, o, datablock, manipulator)

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)

    def press(self, context, event):
        if self.handle_right.hover:
            self.handle_right.active = True
            return True
        return False

    def release(self, context, event):
        self.check_hover()
        self.handle_right.active = False
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active:
            self.update(context, event)
            return True
        else:
            self.check_hover()
        return False

    def update(self, context, event):
        # 0  1  2
        # |_____|
        #
        pt = self.get_pos3d(context)
        pt, t = intersect_point_line(pt, self.line_0.p, self.line_2.p)
        length = (self.line_0.p - pt).length
        if event.alt:
            length = round(length, 1)
        self.set_value(context, self.datablock, self.manipulator.prop1_name, length)

    def draw_callback(self, _self, context):
        """
            draw on screen feedback using gl.
        """
        left, right, side, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.origin = left
        self.line_1.p = left
        self.line_1.v = right - left
        self.line_0.z_axis = normal
        self.line_1.z_axis = normal
        self.line_2.z_axis = normal
        self.label.z_axis = normal
        self.line_0 = self.line_1.sized_normal(0, side.x * 1.1)
        self.line_2 = self.line_1.sized_normal(1, side.x * 1.1)
        self.line_1.offset(side.x * 1.0)
        self.handle_left.set_pos(context, self.line_1.p, -self.line_1.v, normal=normal)
        self.handle_right.set_pos(context, self.line_1.lerp(1), self.line_1.v, normal=normal)
        self.label.set_pos(context, self.line_1.length, self.line_1.lerp(0.5), self.line_1.v, normal=normal)
        self.label.draw(context)
        self.line_0.draw(context)
        self.line_1.draw(context)
        self.line_2.draw(context)
        self.handle_left.draw(context)
        self.handle_right.draw(context)


class SizeLocationManipulator(SizeManipulator):
    def __init__(self, context, o, datablock, manipulator, handle_size):
        SizeManipulator.__init__(self, context, o, datablock, manipulator, handle_size)
        self.handle_left.selectable = True
        self.dl = 0

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)
        self.handle_left.check_hover(self.mouse_pos)

    def press(self, context, event):
        self.dl = 0
        if self.handle_right.hover:
            self.handle_right.active = True
            return True
        if self.handle_left.hover:
            self.handle_left.active = True
            return True
        return False

    def release(self, context, event):
        self.check_hover()
        self.handle_right.active = False
        self.handle_left.active = False
        # self.move_linked(context, self.manipulator.prop2_name, self.dl)
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active or self.handle_left.active:
            self.update(context, event)
            return True
        else:
            self.check_hover()
        return False

    def update(self, context, event):
        # 0  1  2
        # |_____|
        #
        pt = self.get_pos3d(context)
        pt, t = intersect_point_line(pt, self.line_0.p, self.line_2.p)

        len_0 = (pt - self.line_0.p).length
        len_1 = (pt - self.line_2.p).length

        length = max(len_0, len_1)

        if event.alt:
            length = round(length, 1)

        dl = length - self.line_1.length

        if len_0 > len_1:
            dl = 0.5 * dl
        else:
            dl = -0.5 * dl
        self.dl += dl
        self.move(context, self.manipulator.prop2_name, dl)
        self.set_value(context, self.datablock, self.manipulator.prop1_name, length)
        self.move_linked(context, self.manipulator.prop2_name, dl)


class DeltaLocationManipulator(SizeManipulator):

    def __init__(self, context, o, datablock, manipulator, handle_size):
        SizeManipulator.__init__(self, context, o, datablock, manipulator, handle_size)

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)

    def press(self, context, event):
        if self.handle_right.hover:
            self.handle_right.active = True
            return True
        return False

    def release(self, context, event):
        self.check_hover()
        self.handle_right.active = False
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active:
            self.update(context, event)
            return True
        else:
            self.check_hover()
        return False

    def update(self, context, event):
        # 0  1  2
        # |_____|
        #
        p0 = self.line_1.p
        c = self.line_1.lerp(0.5)
        pt = self.get_pos3d(context)
        pt, t = intersect_point_line(pt, p0, c)
        len_0 = 0.5 * self.line_1.length
        len_1 = (pt - p0).length
        dl = (pt - c).length
        if event.alt:
            dl = round(dl, 1)
        if len_0 < len_1:
            dl = dl
        else:
            dl = -dl
        self.move(context, self.manipulator.prop1_name, dl)

    def draw_callback(self, _self, context):
        """
            draw on screen feedback using gl.
        """
        left, right, side, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.origin = left
        self.line_1.p = left
        self.line_1.v = right - left
        self.line_1.z_axis = normal
        self.handle_left.set_pos(context, self.line_1.lerp(0.5), -self.line_1.v, normal=normal)
        self.handle_right.set_pos(context, self.line_1.lerp(0.5), self.line_1.v, normal=normal)
        self.handle_left.draw(context)
        self.handle_right.draw(context)


class DumbSizeManipulator(SizeManipulator):
    """
        Show size while not being editable
    """

    def __init__(self, context, o, datablock, manipulator, handle_size):
        SizeManipulator.__init__(self, context, o, datablock, manipulator, handle_size)
        self.handle_right.selectable = False

    def mouse_move(self, context, event):
        return False


class AngleManipulator(Manipulator):
    """
        Manipulate angle between segments
        bound to [-pi, pi]
    """

    def __init__(self, context, o, datablock, manipulator, handle_size):

        # Angle
        self.handle_right = TriHandle(handle_size, arrow_size, selectable=True)
        self.handle_center = SquareHandle(handle_size, arrow_size)
        self.arc = GlArc()
        self.line_0 = GlLine()
        self.line_1 = GlLine()
        self.label_a = GlText()
        self.active_right = False
        Manipulator.__init__(self, context, o, datablock, manipulator)
        self.pts_mode = 'RADIUS'

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)

    def press(self, context, event):
        if self.handle_right.hover:
            self.handle_right.active = True
            return True
        return False

    def release(self, context, event):
        self.check_hover()
        self.handle_right.active = False
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active:
            self.update(context, event)
            return True
        else:
            self.check_hover()
        return False

    def update(self, context, event):
        pt = self.get_pos3d(context)
        c = self.arc.c
        v = 2 * self.arc.r * (pt - c).normalized()
        v0 = c - v
        v1 = c + v
        p0, p1 = intersect_line_sphere(v0, v1, c, self.arc.r)
        if p0 is not None and p1 is not None:
            if (p1 - pt).length < (p0 - pt).length:
                p0, p1 = p1, p0
            v = p0 - self.arc.c
            da = atan2(v.y, v.x) - self.line_0.angle
            if da > pi:
                da = da - 2 * pi
            if da < -pi:
                da = da + 2 * pi
            # print("a:%.4f da:%.4f a0:%.4f" % (atan2(v.y, v.x), da, self.line_0.angle))
            if da > pi:
                da = pi
            if da < -pi:
                da = -pi
            if event.alt:
                da = round(da / pi * 180, 0) / 180 * pi
            self.set_value(context, self.datablock, self.manipulator.prop1_name, da)

    def draw_callback(self, _self, context):
        c, left, right, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.line_0.z_axis = normal
        self.line_1.z_axis = normal
        self.arc.z_axis = normal
        self.label_a.z_axis = normal
        self.origin = c
        self.line_0.p = c
        self.line_1.p = c
        self.arc.c = c
        self.line_0.v = left
        self.line_0.v = -self.line_0.cross.normalized()
        self.line_1.v = right
        self.line_1.v = self.line_1.cross.normalized()
        self.arc.a0 = self.line_0.angle
        self.arc.da = self.get_value(self.datablock, self.manipulator.prop1_name)
        self.arc.r = 1.0
        self.handle_right.set_pos(context, self.line_1.lerp(1),
                                  self.line_1.sized_normal(1, -1 if self.arc.da > 0 else 1).v)
        self.handle_center.set_pos(context, self.arc.c, -self.line_0.v)
        self.label_a.set_pos(context, self.arc.da / pi * 180, self.arc.lerp(0.5), -self.line_0.v)
        self.arc.draw(context)
        self.line_0.draw(context)
        self.line_1.draw(context)
        self.handle_right.draw(context)
        self.handle_center.draw(context)
        self.label_a.draw(context)


class ArcAngleManipulator(Manipulator):
    """
        Manipulate angle of an arc
        when angle < 0 the arc center is on the left part of the circle
        when angle > 0 the arc center is on the right part of the circle
        bound to [-pi, pi]
    """

    def __init__(self, context, o, datablock, manipulator, handle_size):

        # Fixed
        self.handle_left = SquareHandle(handle_size, arrow_size)
        # Angle
        self.handle_right = TriHandle(handle_size, arrow_size, selectable=True)
        self.handle_center = SquareHandle(handle_size, arrow_size)
        self.arc = GlArc()
        self.line_0 = GlLine()
        self.line_1 = GlLine()
        self.label_a = GlText()
        self.label_r = GlText()
        self.active_right = False
        Manipulator.__init__(self, context, o, datablock, manipulator)
        self.pts_mode = 'RADIUS'

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)

    def press(self, context, event):
        if self.handle_right.hover:
            self.handle_right.active = True
            return True
        return False

    def release(self, context, event):
        self.check_hover()
        self.handle_right.active = False
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active:
            self.update(context, event)
            return True
        else:
            self.check_hover()
        return False

    def update(self, context, event):
        pt = self.get_pos3d(context)
        c = self.arc.c
        v = 2 * self.arc.r * (pt - c).normalized()
        v0 = c - v
        v1 = c + v
        p0, p1 = intersect_line_sphere(v0, v1, c, self.arc.r)
        if p0 is not None:
            if (p1 - pt).length < (p0 - pt).length:
                p0, p1 = p1, p0
            v = p0 - self.arc.c
            da = atan2(v.y, v.x) - self.line_0.angle
            if abs(da) > pi:
                da = pi
            # bottom and top points
            n = self.line_0.sized_normal(0, 1)
            top = pt - n.lerp(-1)
            bottom = pt - n.lerp(1)
            # left and right points
            n = self.arc.sized_normal(0, 1)
            right = pt - n.lerp(1)
            left = pt - n.lerp(-1)
            # we are on the right part
            if right.length < left.length:
                da = -abs(da)
            else:
                # on bottom part da = pi
                if bottom.length < top.length:
                    da = pi
                else:
                    da = abs(da)
            if event.alt:
                da = round(da / pi * 180, 0) / 180 * pi
            self.set_value(context, self.datablock, self.manipulator.prop1_name, da)

    def draw_callback(self, _self, context):
        # center : 3d points
        # left   : 3d vector pt-c
        # right  : 3d vector pt-c
        c, left, right, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.line_0.z_axis = normal
        self.line_1.z_axis = normal
        self.arc.z_axis = normal
        self.label_a.z_axis = normal
        self.label_r.z_axis = normal
        self.origin = c
        self.line_0.p = c
        self.line_1.p = c
        self.arc.c = c
        self.line_0.v = left
        self.line_1.v = right
        self.arc.a0 = self.line_0.angle
        self.arc.da = self.get_value(self.datablock, self.manipulator.prop1_name)
        self.arc.r = left.length
        self.handle_left.set_pos(context, self.line_0.lerp(1), self.line_0.v)
        self.handle_right.set_pos(context, self.line_1.lerp(1),
                                  self.line_1.sized_normal(1, -1 if self.arc.da > 0 else 1).v)
        self.handle_center.set_pos(context, self.arc.c, -self.line_0.v)
        self.label_a.set_pos(context, self.arc.da / pi * 180, self.arc.lerp(0.5), -self.line_0.v)
        self.label_r.set_pos(context, self.arc.r, self.line_0.lerp(0.5), self.line_0.v)
        self.arc.draw(context)
        self.line_0.draw(context)
        self.line_1.draw(context)
        self.handle_left.draw(context)
        self.handle_right.draw(context)
        self.handle_center.draw(context)
        self.label_r.draw(context)
        self.label_a.draw(context)


class ArcAngleRadiusManipulator(ArcAngleManipulator):

    def __init__(self, context, o, datablock, manipulator, handle_size):
        ArcAngleManipulator.__init__(self, context, o, datablock, manipulator, handle_size)
        self.handle_center = TriHandle(handle_size, arrow_size, selectable=True)
        self.active_center = False

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)
        self.handle_center.check_hover(self.mouse_pos)

    def press(self, context, event):
        if self.handle_right.hover:
            self.handle_right.active = True
            return True
        if self.handle_center.hover:
            self.handle_center.active = True
            return True
        return False

    def release(self, context, event):
        self.check_hover()
        self.handle_right.active = False
        self.handle_center.active = False
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active:
            self.update(context, event)
            return True
        elif self.handle_center.active:
            self.update_radius(context, event)
            return True
        else:
            self.check_hover()
        return False

    def update_radius(self, context, event):
        pt = self.get_pos3d(context)
        c = self.arc.c
        left = self.line_0.lerp(1)
        p, t = intersect_point_line(pt, c, left)
        radius = (left - p).length
        if event.alt:
            radius = round(radius, 1)
        self.set_value(context, self.datablock, self.manipulator.prop2_name, radius)


# ------------------------------------------------------------------
# Define a single Manipulator Properties to store on object
# ------------------------------------------------------------------


# Allow registering manipulators classes
manipulators_class_lookup = {}


def register_manipulator(type_key, manipulator_class):
    if type_key in manipulators_class_lookup.keys():
        raise RuntimeError("Manipulator of type {} allready exists, unable to override".format(type_key))
    manipulators_class_lookup[type_key] = manipulator_class


# Register default manipulators
register_manipulator('SIZE', SizeManipulator)
register_manipulator('SIZE_LOC', SizeLocationManipulator)
register_manipulator('ANGLE', AngleManipulator)
register_manipulator('ARC_ANGLE_RADIUS', ArcAngleRadiusManipulator)
register_manipulator('COUNTER', CounterManipulator)
register_manipulator('DUMB_SIZE', DumbSizeManipulator)
register_manipulator('DELTA_LOC', DeltaLocationManipulator)
# wall's np_station based snap
# register_manipulator('SNAP_POINT', SnapPointManipulator)
register_manipulator('WALL_SNAP', WallSnapManipulator)
register_manipulator('FENCE_SNAP', FenceSnapManipulator)


class archipack_manipulator(PropertyGroup):
    """
        A property group to add to manipulable objects
        type: type of manipulator
        prop1_name = the property name of object to modify
        prop2_name = another property name of object to modify (angle and radius)
        p0, p1, p2 3d Vectors as base points to represent manipulators on screen
        normal Vector normal of plane on with draw manipulator
    """
    type = StringProperty(default='SIZE')

    # How 3d points are stored in manipulators
    # SIZE = 2 absolute positionned and a scaling vector
    # RADIUS = 1 absolute positionned and 2 relatives
    pts_mode = StringProperty(default='SIZE')

    """
    EnumProperty(
        items=(
            ('SIZE', 'Size', 'Generic size manipulator', 0),
            ('SIZE_LOC', 'Size Location', 'Generic size from border manipulator', 1),
            ('ANGLE', 'Angle', 'Angle between two vectors', 2),
            ('ARC_ANGLE_RADIUS', 'Arc based angle', '', 3),
            ('COUNTER', 'Counter increase and decrease', '', 4),
            ('DUMB_SIZE', 'Dumb Size', 'Generic size not editable', 5),
            ('DELTA_LOC', 'Delta location', 'Move object on an axis', 6)
        ),
        default='SIZE'

    )
    """
    prop1_name = StringProperty()
    prop2_name = StringProperty()
    p0 = FloatVectorProperty(subtype='XYZ')
    p1 = FloatVectorProperty(subtype='XYZ')
    p2 = FloatVectorProperty(subtype='XYZ')
    normal = FloatVectorProperty(subtype='XYZ', default=(0, 0, 1))

    def set_pts(self, pts):
        self.p0, self.p1, self.p2 = pts

    def get_pts(self, tM):
        rM = tM.to_3x3()
        if self.pts_mode == 'SIZE':
            return tM * self.p0, tM * self.p1, self.p2, rM * self.normal
        else:
            return tM * self.p0, rM * self.p1, rM * self.p2, rM * self.normal

    def setup(self, context, o, datablock):
        """
            Factory return a manipulator object or None
            o:         object
            datablock: datablock to modify
        """
        global handle_size
        global manipulators_class_lookup

        if self.type not in manipulators_class_lookup.keys() or \
                not manipulators_class_lookup[self.type].poll(context):
            # RuntimeError is overkill here
            # silentely ignore allow skipping manipulators when deps as not meet
            # manip stack will simply be filled with None objects
            return None
            # raise RuntimeError("Manipulator of type {} not found".format(self.type))

        m = manipulators_class_lookup[self.type](context, o, datablock, self, handle_size)
        self.pts_mode = m.pts_mode
        return m


bpy.utils.register_class(archipack_manipulator)


# a global manipulator stack reference
# prevent Blender "ACCESS_VIOLATION" crashes
# use a dict to prevent potential
# collisions between many objects being in
# manipulate mode (at create time)
# use object names as loose keys
# NOTE : use app.drivers to reset all before file load
manip_stack = {}


@persistent
def empty_stack(dummy):
    """
        Empty manipulators stack on file load
    """
    global manip_stack
    for key in manip_stack.keys():
        for m in manip_stack[key]:
            if m is not None:
                m.exit()
        del manip_stack[key]


def exit_stack(key):
    """
        remove object from stack
    """
    global manip_stack
    if key in manip_stack.keys():
        for m in manip_stack[key]:
            if m is not None:
                m.exit()
        del manip_stack[key]


def clean_stack():
    """
        remove references to renamed or deleted objects from stack
    """
    global manip_stack
    for key in manip_stack.keys():
        if bpy.data.objects.find(key) is None:
            exit_stack(key)


def get_stack(key):
    """
        return reference to manipulator stack for given object
    """
    global manip_stack
    if key not in manip_stack.keys():
        manip_stack[key] = []
    return manip_stack[key]


bpy.app.handlers.load_pre.append(empty_stack)

# ------------------------------------------------------------------
# Define Manipulable to make a PropertyGroup manipulable
# ------------------------------------------------------------------


class Manipulable():
    """
        A class extending PropertyGroup to setup gl manipulators
        Beware : prevent crash calling manipulable_disable()
                 before changing manipulated data structure
    """
    manipulators = CollectionProperty(
            type=archipack_manipulator,
            description="store 3d points to draw gl manipulators"
            )
    manipulable_refresh = BoolProperty(
            default=False,
            options={'SKIP_SAVE'},
            description="Flag enable to rebuild manipulators when data model change"
            )
    manipulate_mode = BoolProperty(
            default=False,
            options={'SKIP_SAVE'},
            description="Flag to disable manipulators on click"
            )

    def manipulable_disable(self, context):
        """
            disable gl draw handlers
        """
        self.manipulate_mode = False

        clean_stack()

        o = context.active_object
        if o is not None:
            exit_stack(o.name)
            self.manip_stack = get_stack(o.name)

    def manipulable_setup(self, context):
        """
            TODO: Implement the setup part as per parent object basis
        """
        self.manipulable_disable(context)
        o = context.active_object
        for m in self.manipulators:
            self.manip_stack.append(m.setup(context, o, self))

    def manipulable_invoke(self, context):
        """
            call this in operator invoke()
        """
        # print("self.manipulate_mode:%s" % (self.manipulate_mode))
        if self.manipulate_mode:
            self.manipulable_disable(context)
            return False

        self.manip_stack = []
        self.manipulable_setup(context)
        self.manipulate_mode = True
        return True

    def manipulable_modal(self, context, event):
        """
            call in operator modal()
        """
        # setup again when manipulators type change
        if self.manipulable_refresh:
            self.manipulable_refresh = False
            self.manipulable_setup(context)
            self.manipulate_mode = True

        context.area.tag_redraw()

        # clean up manipulator on delete
        if event.type in {'X'}:
            self.manipulable_disable(context)
            bpy.ops.object.delete('INVOKE_DEFAULT', use_global=False)
            return {'FINISHED'}

        for m in self.manip_stack:
            # m should return false on left mouse release
            if m is not None and m.modal(context, event):
                self.manipulable_manipulate(context, event=event, manipulator=m)
                return {'RUNNING_MODAL'}

        # allow any action on release
        if event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            self.manipulable_release(context)

        if event.type in {'RIGHTMOUSE', 'ESC'}:
            self.manipulable_disable(context)
            self.manipulable_exit(context)
            return {'FINISHED'}

        return {'PASS_THROUGH'}

    # Callbacks
    def manipulable_release(self, context):
        """
            Override with action to do on mouse release
            eg: big update
        """
        return

    def manipulable_exit(self, context):
        """
            Override with action to do when modal exit
        """
        return

    def manipulable_manipulate(self, context, event=None, manipulator=None):
        """
            Override with action to do when a handle is active (pressed and mousemove)
        """
        return
