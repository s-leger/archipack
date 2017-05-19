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
from bpy_extras import view3d_utils, object_utils
from bpy.types import PropertyGroup
from bpy.props import FloatVectorProperty, StringProperty, CollectionProperty, BoolProperty
from bpy.app.handlers import persistent
from .archipack_snap import snap_point
from .archipack_keymaps import Keymaps

# ------------------------------------------------------------------
# Define Gl Handle types
# ------------------------------------------------------------------

# Arrow sizes (world units)
arrow_size = 0.1
# Handle area size (pixels)
handle_size = 10

# Font sizes and basic colour scheme
# kept outside of addon prefs until now
# as for a generic toolkit it is not appropriate
# we could provide a template for addon prefs
# matching those one
feedback_size_main = 16
feedback_size_title = 14
feedback_size_shortcut = 11
feedback_colour_main = (0.95, 0.95, 0.95, 1.0)
feedback_colour_key = (0.67, 0.67, 0.67, 1.0)
feedback_colour_shortcut = (0.51, 0.51, 0.51, 1.0)


# @TODO:
# 1 Make a clear separation of 2d (pixel position) and 3d (world position)
#   modes way to set gl coords
# 2 Unify methods to set points - currently set_pts, set_pos ...
# 3 Put all Gl part in a sub module as it may be used by other devs
#   as gl toolkit abstraction for screen feedback
# 4 Implement cursor badges (np_station sample)
# 5 Define a clear color scheme so it is easy to customize
# 6 Allow different arguments for each classes like
#   eg: for line p0 p1, p0 and vector (p1-p0)
#       raising exceptions when incomplete
# 7 Use correct words, normal is not realy a normal
#   but a perpendicular
# May be hard code more shapes ?
# Fine tuned text styles with shadows and surronding boxes / backgrounds
# Extending tests to hdr screens, ultra wide ones and so on
# Circular handle, handle styling (only border, filling ...)

# Keep point 3 in mind while doing this, to keep it simple and easy to use
# Take inspiration from other's feed back systems, talk to other devs
# and find who actually work on bgl future for 2.8 release


class Gl():
    """
        handle 3d -> 2d gl drawing
        d : dimensions
            3 to convert pos from 3d
            2 to keep pos as 2d absolute screen position
    """
    def __init__(self, d=3):
        # flag output over render image
        self.render = False
        # default line width
        self.width = 1
        # default line style
        self.style = bgl.GL_LINE
        # allow closed lines
        self.closed = False
        # nth dimensions of input coords 3=word coords 2=pixel screen coords
        self.d = d
        self.pos_2d = Vector((0, 0))
        self.colour_active = (1.0, 0.0, 0.0, 1.0)
        self.colour_hover = (1.0, 1.0, 0.0, 1.0)
        self.colour_normal = (1.0, 1.0, 1.0, 1.0)
        self.colour_inactive = (0.0, 0.0, 0.0, 1.0)

    @property
    def colour(self):
        return self.colour_inactive

    def position_2d_from_coord(self, context, coord, render=False):
        """ coord given in local input coordsys
        """
        if self.d == 2:
            return coord
        if render:
            return self.get_render_location(context, coord)
        region = context.region
        rv3d = context.region_data
        loc = view3d_utils.location_3d_to_region_2d(region, rv3d, coord, self.pos_2d)
        return loc

    def get_render_location(self, context, coord):
        scene = context.scene
        co_2d = object_utils.world_to_camera_view(scene, scene.camera, coord)
        # Get pixel coords
        render_scale = scene.render.resolution_percentage / 100
        render_size = (int(scene.render.resolution_x * render_scale),
                       int(scene.render.resolution_y * render_scale))
        return [round(co_2d.x * render_size[0]), round(co_2d.y * render_size[1])]

    def _end(self):
        bgl.glEnd()
        bgl.glPopAttrib()
        bgl.glLineWidth(1)
        bgl.glDisable(bgl.GL_BLEND)
        bgl.glColor4f(0.0, 0.0, 0.0, 1.0)

    def _start_poly(self):
        bgl.glPushAttrib(bgl.GL_ENABLE_BIT)
        bgl.glEnable(bgl.GL_BLEND)
        if self.render:
            # enable anti-alias on polygons
            bgl.glEnable(bgl.GL_POLYGON_SMOOTH)
        bgl.glColor4f(*self.colour)
        bgl.glBegin(bgl.GL_POLYGON)

    def _start_line(self):
        bgl.glPushAttrib(bgl.GL_ENABLE_BIT)
        if self.style == bgl.GL_LINE_STIPPLE:
            bgl.glLineStipple(1, 0x9999)
        bgl.glEnable(self.style)
        bgl.glEnable(bgl.GL_BLEND)
        if self.render:
            # enable anti-alias on lines
            bgl.glEnable(bgl.GL_LINE_SMOOTH)
        bgl.glColor4f(*self.colour)
        bgl.glLineWidth(self.width)
        if self.closed:
            bgl.glBegin(bgl.GL_LINE_LOOP)
        else:
            bgl.glBegin(bgl.GL_LINE_STRIP)

    def draw_text(self, x, y):
        # dirty fast assignment
        dpi, font_id = 72, 0
        bgl.glColor4f(*self.colour)
        if self.angle != 0:
            blf.enable(font_id, blf.ROTATION)
            blf.rotation(font_id, self.angle)
        blf.size(font_id, self.font_size, dpi)
        blf.position(font_id, x, y, 0)
        blf.draw(font_id, self.text)
        if self.angle != 0:
            blf.disable(font_id, blf.ROTATION)

    def draw(self, context, render=False):
        """
            render flag when rendering
        """
        self.render = render
        gl_type = type(self).__name__
        if 'Handle' in gl_type or gl_type in ['GlPolygon']:
            self._start_poly()
        elif gl_type in ['GlArc', 'GlLine', 'GlPolyline']:
            self._start_line()
        if 'Text' in gl_type:
            x, y = self.position_2d_from_coord(context, self.pts[0], render)
            self.draw_text(x, y)
        else:
            for pt in self.pts:
                x, y = self.position_2d_from_coord(context, pt, render)
                bgl.glVertex2f(x, y)
            self._end()


class GlText(Gl):

    def __init__(self, d=3, round=2, label='', font_size=16, colour=(1, 1, 1, 1), z_axis=Vector((0, 0, 1))):
        self.z_axis = z_axis
        self.value = None
        self.round = round
        self.label = label
        self.unit = ''
        self.font_size = font_size
        self.angle = 0
        Gl.__init__(self, d)
        self.colour_inactive = colour

    @property
    def text_size(self):
        dpi, font_id = 72, 0
        blf.enable(font_id, blf.ROTATION)
        blf.rotation(font_id, self.angle)
        blf.size(font_id, self.font_size, dpi)
        x, y = blf.dimensions(font_id, self.text)
        blf.disable(font_id, blf.ROTATION)
        return Vector((x, y))

    @property
    def pts(self):
        return [self.pos_3d]

    @property
    def text(self):
        if self.value is not None:
            return self.label + str(round(self.value, self.round)) + self.unit
        else:
            return self.label

    def set_pos(self, context, value, pos_3d, direction, angle=0, normal=Vector((0, 0, 1))):
        self.up_axis = direction.normalized()
        self.c_axis = self.up_axis.cross(normal)
        self.pos_3d = pos_3d
        self.value = value
        self.angle = angle


class GlLine(Gl):

    def __init__(self, d=3, z_axis=Vector((0, 0, 1))):
        self.z_axis = z_axis
        self.p = Vector((0, 0, 0))
        self.v = Vector((0, 0, 0))
        Gl.__init__(self, d)

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

    def __init__(self, d=3):
        self.r = 0
        self.c = Vector((0, 0, 0))
        Gl.__init__(self, d)


class GlArc(GlCircle):

    def __init__(self, d=3, z_axis=Vector((0, 0, 1))):
        """
            a0 and da arguments are in radians
            a0 = 0   on the right side
            a0 = pi on the left side
            da > 0 CCW contrary-clockwise
            da < 0 CW  clockwise
            stored internally as radians
        """
        GlCircle.__init__(self, d)
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
    def __init__(self, colour, d=3):
        self.pts_3d = []
        Gl.__init__(self, d)
        self.colour_inactive = colour

    def set_pos(self, pts_3d):
        self.pts_3d = pts_3d

    @property
    def pts(self):
        return self.pts_3d


class GlPolyline(Gl):
    def __init__(self, colour, d=3):
        self.pts_3d = []
        Gl.__init__(self, d=d)
        self.colour_inactive = colour

    def set_pos(self, pts_3d):
        self.pts_3d = pts_3d
        # self.pts_3d.append(pts_3d[0])

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
        self.sensor_width = sensor_size
        self.sensor_height = sensor_size
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
        self.pos_2d = self.position_2d_from_coord(context, self.sensor_center)

    def check_hover(self, pos_2d):
        if self.selectable:
            dp = pos_2d - self.pos_2d
            self.hover = abs(dp.x) < self.sensor_width and abs(dp.y) < self.sensor_height

    @property
    def sensor_center(self):
        pts = self.pts
        n = len(pts)
        x, y, z = 0, 0, 0
        for pt in pts:
            x += pt.x
            y += pt.y
            z += pt.z
        return Vector((x / n, y / n, z / n))

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


class EditableText(GlText, GlHandle):
    def __init__(self, sensor_size, size, selectable=False):
        GlHandle.__init__(self, sensor_size, size, selectable)
        GlText.__init__(self)

    def set_pos(self, context, value, pos_3d, direction, normal=Vector((0, 0, 1))):
        self.up_axis = direction.normalized()
        self.c_axis = self.up_axis.cross(normal)
        self.pos_3d = pos_3d
        self.value = value
        self.sensor_width, self.sensor_height = self.text_size
        self.pos_2d = self.position_2d_from_coord(context, pos_3d)

    @property
    def sensor_center(self):
        return self.pos_3d


class FeedbackPanel():
    """
        Feed-back panel
        inspired by np_station
    """
    def __init__(self, title='Archipack'):

        self.main_title = GlText(d=2, label=title, font_size=feedback_size_main, colour=feedback_colour_main)
        self.title = GlText(d=2, font_size=feedback_size_title, colour=feedback_colour_main)
        self.spacing = Vector((0.5 * feedback_size_shortcut, 0.5 * feedback_size_shortcut))
        self.margin = 50
        self.explanation = GlText(d=2, font_size=feedback_size_shortcut, colour=feedback_colour_main)
        self.shortcut_area = GlPolygon(colour=(0, 0.4, 0.6, 0.2), d=2)
        self.title_area = GlPolygon(colour=(0, 0.4, 0.6, 0.5), d=2)
        self.shortcuts = []
        self.on = False

    def disable(self):
        self.on = False

    def enable(self):
        self.on = True

    def instructions(self, context, title, explanation, shortcuts):
        """
            position from bottom to top
        """
        w = context.region.width
        # h = context.region.height
        # 0,0 = bottom left
        pos = Vector((self.margin + self.spacing.x, self.margin))
        self.shortcuts = []

        if len(shortcuts) > 0:
            pos.y += self.spacing.y

        add_y = 0

        for key, label in shortcuts:
            key = GlText(d=2, label=key + ' : ', font_size=feedback_size_shortcut, colour=feedback_colour_key)
            label = GlText(d=2, label=label, font_size=feedback_size_shortcut, colour=feedback_colour_shortcut)
            ks = key.text_size
            ls = label.text_size
            space = ks.x + ls.x + self.spacing.x
            add_y = ks.y + self.spacing.y
            if pos.x + space > w - 2 * self.margin:
                add_y = 0
                pos.y += ks.y + 2 * self.spacing.y
                pos.x = self.margin + self.spacing.x
            key.pos_3d = pos.copy()
            pos.x += ks.x
            label.pos_3d = pos.copy()
            pos.x += ls.x + 2 * self.spacing.x
            self.shortcuts.extend([key, label])

        if len(shortcuts) > 0:
            pos.y += add_y + 0.5 * self.spacing.y

        self.shortcut_area.pts_3d = [
            (self.margin, self.margin),
            (w - self.margin, self.margin),
            (w - self.margin, pos.y),
            (self.margin, pos.y)
            ]

        if len(shortcuts) > 0:
            pos.y += 0.5 * self.spacing.y

        self.title.label = ' : ' + title
        main_title_size = self.main_title.text_size
        self.title_area.pts_3d = [
            (self.margin, pos.y),
            (w - self.margin, pos.y),
            (w - self.margin, pos.y + main_title_size.y + 2 * self.spacing.y),
            (self.margin, pos.y + main_title_size.y + 2 * self.spacing.y)
            ]
        pos.y += self.spacing.y
        self.explanation.label = explanation
        self.explanation.pos_3d = (w - self.margin - self.spacing.x - self.explanation.text_size.x, pos.y)
        self.main_title.pos_3d = (self.margin + self.spacing.x, pos.y)
        self.title.pos_3d = (self.margin + self.spacing.x + main_title_size.x, pos.y)

    def draw(self, context, render=False):
        if self.on:
            """
                draw from bottom to top
                so we are able to always fit needs
            """
            self.shortcut_area.draw(context)
            self.title_area.draw(context)
            self.main_title.draw(context)
            self.title.draw(context)
            self.explanation.draw(context)
            for s in self.shortcuts:
                s.draw(context)


class GlCursorFence():
    """
        Cursor crossing Fence
    """
    def __init__(self, width=1, colour=(1.0, 1.0, 1.0, 0.5), style=bgl.GL_LINE_STIPPLE):
        self.line_x = GlLine(d=2)
        self.line_x.style = style
        self.line_x.width = width
        self.line_x.colour_inactive = colour
        self.line_y = GlLine(d=2)
        self.line_y.style = style
        self.line_y.width = width
        self.line_y.colour_inactive = colour
        self.on = True

    def set_location(self, context, location):
        w = context.region.width
        h = context.region.height
        x, y = location
        self.line_x.p = Vector((0, y))
        self.line_x.v = Vector((w, 0))
        self.line_y.p = Vector((x, 0))
        self.line_y.v = Vector((0, h))

    def enable(self):
        self.on = True

    def disable(self):
        self.on = False

    def draw(self, context, render=False):
        if self.on:
            self.line_x.draw(context)
            self.line_y.draw(context)


class GlCursorArea():
    def __init__(self,
                width=1,
                bordercolour=(1.0, 1.0, 1.0, 0.5),
                areacolour=(1.0, 1.0, 1.0, 0.1),
                style=bgl.GL_LINE_STIPPLE):

        self.border = GlPolyline(bordercolour, d=2)
        self.border.style = style
        self.border.width = width
        self.area = GlPolygon(areacolour, d=2)
        self.on = False

    def set_location(self, context, p0, p1):
        x0, y0 = p0
        x1, y1 = p1
        pos = [Vector((x0, y0)), Vector((x0, y1)), Vector((x1, y1)), Vector((x1, y0))]
        self.area.set_pos(pos)
        self.border.set_pos(pos)

    def enable(self):
        self.on = True

    def disable(self):
        self.on = False

    def draw(self, context, render=False):
        if self.on:
            self.area.draw(context)
            self.border.draw(context)


# ------------------------------------------------------------------
# Define Manipulators
# ------------------------------------------------------------------


class Manipulator():
    """
        Manipulator base class to derive other
        handle keyboard and modal events
        provide convenient funcs including getter and setter for datablock values
        store reference of base object, datablock and manipulator
    """
    keyboard_ascii = {
            ".", ",", "-", "+", "1", "2", "3",
            "4", "5", "6", "7", "8", "9", "0",
            "c", "m", "d", "k", "h", "a",
            " ", "/", "*", "'", "\""
            # "="
            }
    keyboard_type = {
            'BACK_SPACE', 'DEL',
            'LEFT_ARROW', 'RIGHT_ARROW'
            }

    def __init__(self, context, o, datablock, manipulator):
        """
            o : object to manipulate
            datablock : object data to manipulate
            manipulator: object archipack_manipulator datablock
        """
        self.keymap = Keymaps(context)
        self.feedback = FeedbackPanel()
        self.active = False
        # active text input value for manipulator
        self.keyboard_input_active = False
        self.label_value = 0
        # unit for keyboard input value
        self.value_type = 'LENGTH'
        self.pts_mode = 'SIZE'
        self.o = o
        self.datablock = datablock
        self.manipulator = manipulator
        self.origin = Vector((0, 0, 1))
        self.mouse_pos = Vector((0, 0))
        self.length_entered = ""
        self.line_pos = 0
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_callback, args, 'WINDOW', 'POST_PIXEL')

    @classmethod
    def poll(cls, context):
        """
            Allow manipulator enable/disable
            in given context
            handles will not show
        """
        return True

    def exit(self):
        """
            Modal exit, DONT EVEN TRY TO OVERRIDE
        """
        # print("Manipulator.exit() %s" % (type(self).__name__))
        if self._handle is not None:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
        self._handle = None

    # Mouse event handlers, MUST be overriden
    def mouse_press(self, context, event):
        """
            Manipulators must implement
            mouse press event handler
            return True to callback manipulable_manipulate
        """
        raise NotImplementedError

    def mouse_release(self, context, event):
        """
            Manipulators must implement
            mouse mouse_release event handler
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

    # Keyboard event handlers, MAY be overriden
    def keyboard_done(self, context, event, value):
        """
            Manipulators may implement
            keyboard value validated event handler
            value: changed by keyboard
            return True to callback manipulable_manipulate
        """
        return False

    def keyboard_editing(self, context, event, value):
        """
            Manipulators may implement
            keyboard value changed event handler
            value: string changed by keyboard
            allow realtime update of label
            return False to show edited value on window header
            return True when feedback show right on screen
        """
        if self.value_type == 'ROTATION':
            self.label_value = value / pi * 180
        else:
            self.label_value = value
        return True

    def keyboard_cancel(self, context, event):
        """
            Manipulators may implement
            keyboard entry cancelled
        """
        return

    def cancel(self, context, event):
        """
            Manipulators may implement
            cancelled event (ESC RIGHTCLICK)
        """
        self.active = False
        return

    def undo(self, context, event):
        """
            Manipulators may implement
            undo event (CTRL+Z)
        """
        return False

    # Internal, do not override unless you realy
    # realy realy deeply know what you are doing
    def keyboard_eval(self, context, event):
        """
            evaluate keyboard entry while typing
            do not override this one
        """
        c = event.ascii
        if c:
            if c == ",":
                c = "."
            self.length_entered = self.length_entered[:self.line_pos] + c + self.length_entered[self.line_pos:]
            self.line_pos += 1

        if self.length_entered:
            if event.type == 'BACK_SPACE':
                self.length_entered = self.length_entered[:self.line_pos - 1] + self.length_entered[self.line_pos:]
                self.line_pos -= 1

            elif event.type == 'DEL':
                self.length_entered = self.length_entered[:self.line_pos] + self.length_entered[self.line_pos + 1:]

            elif event.type == 'LEFT_ARROW':
                self.line_pos = (self.line_pos - 1) % (len(self.length_entered) + 1)

            elif event.type == 'RIGHT_ARROW':
                self.line_pos = (self.line_pos + 1) % (len(self.length_entered) + 1)

        try:
            value = bpy.utils.units.to_value(context.scene.unit_settings.system, self.value_type, self.length_entered)
            draw_on_header = self.keyboard_editing(context, event, value)
        except:  # ValueError:
            draw_on_header = True
            pass

        if draw_on_header:
            a = ""
            if self.length_entered:
                pos = self.line_pos
                a = self.length_entered[:pos] + '|' + self.length_entered[pos:]
            context.area.header_text_set("%s" % (a))

        # modal mode: do not let event bubble up
        return True

    def modal(self, context, event):
        """
            Modal handler
            handle mouse, and keyboard events
            enable and disable feedback
        """
        # print("Manipulator modal:%s %s" % (event.value, event.type))

        if event.type == 'MOUSEMOVE':
            return self.mouse_move(context, event)

        elif event.value == 'PRESS':

            if event.type == 'LEFTMOUSE':
                active = self.mouse_press(context, event)
                if active:
                    self.feedback.enable()
                return active
            
            elif event.shift:
                
                print("GOT a SHIFT")
                
            elif self.keymap.check(event, self.keymap.undo):
                if self.keyboard_input_active:
                    self.keyboard_input_active = False
                    self.keyboard_cancel(context, event)
                self.feedback.disable()
                # prevent undo CRASH
                return True

            elif self.keyboard_input_active and (
                    event.ascii in self.keyboard_ascii or
                    event.type in self.keyboard_type
                    ):
                # get keyboard input
                return self.keyboard_eval(context, event)

            elif event.type in {'ESC', 'RIGHTMOUSE'}:
                self.feedback.disable()
                if self.keyboard_input_active:
                    # allow keyboard exit without setting value
                    self.length_entered = ""
                    self.line_pos = 0
                    self.keyboard_input_active = False
                    self.keyboard_cancel(context, event)
                    return True
                elif self.active:
                    self.cancel(context, event)
                    return True
                return False

        elif event.value == 'RELEASE':

            if event.type == 'LEFTMOUSE':
                if not self.keyboard_input_active:
                    self.feedback.disable()
                return self.mouse_release(context, event)

            elif self.keyboard_input_active and event.type in {'RET', 'NUMPAD_ENTER'}:
                # validate keyboard input
                if self.length_entered != "":
                    try:
                        value = bpy.utils.units.to_value(
                            context.scene.unit_settings.system,
                            'LENGTH', self.length_entered)
                        self.length_entered = ""
                        ret = self.keyboard_done(context, event, value)
                    except:  # ValueError:
                        ret = False
                        self.keyboard_cancel(context, event)
                        self.report({'INFO'}, "Operation not supported yet")
                    context.area.header_text_set()
                    self.keyboard_input_active = False
                    self.feedback.disable()
                    return ret

        return False

    def mouse_position(self, event):
        """
            store mouse position in a 2d Vector
        """
        self.mouse_pos.x, self.mouse_pos.y = event.mouse_region_x, event.mouse_region_y

    def get_pos3d(self, context):
        """
            convert mouse pos to 3d point over plane defined by origin and normal
            pt is in world space
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
        """
            Datablock value getter with index support
        """
        try:
            if index > -1:
                return getattr(data, attr)[index]
            else:
                return getattr(data, attr)
        except:
            print("get_value of %s %s failed" % (data, attr))
            return 0

    def set_value(self, context, data, attr, value, index=-1):
        """
            Datablock value setter with index support
        """
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
            vec = Vector((value, 0, 0))
        elif axis == 'y':
            vec = Vector((0, value, 0))
        else:
            vec = Vector((0, 0, value))
        o.matrix_world = self.preTranslate(o.matrix_world, vec)

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


# OUT OF ORDER
class SnapPointManipulator(Manipulator):
    """
        np_station based snap manipulator
        dosent update anything by itself.
        NOTE : currently out of order
        and disabled in __init__
    """
    def __init__(self, context, o, datablock, manipulator, handle_size):

        raise NotImplementedError

        self.handle = SquareHandle(handle_size, 1.2 * arrow_size, selectable=True)
        Manipulator.__init__(self, context, o, datablock, manipulator)

    def check_hover(self):
        self.handle.check_hover(self.mouse_pos)

    def mouse_press(self, context, event):
        if self.handle.hover:
            self.handle.hover = False
            self.handle.active = True
            self.o.select = True
            # takeloc = self.o.matrix_world * self.manipulator.p0
            # print("Invoke sp_point_move %s" % (takeloc))
            # @TODO:
            # implement and add draw and callbacks
            # snap_point(takeloc, draw, callback)
            return True
        return False

    def mouse_release(self, context, event):
        self.check_hover()
        self.handle.active = False
        # False to callback manipulable_release
        return False

    def update(self, context, event):
        # NOTE:
        # dosent set anything internally
        return

    def mouse_move(self, context, event):
        """

        """
        self.mouse_position(event)
        if self.handle.active:
            # self.handle.active = np_snap.is_running
            # self.update(context)
            # True here to callback manipulable_manipulate
            return True
        else:
            self.check_hover()
        return False

    def draw_callback(self, _self, context, render=False):
        left, right, side, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.handle.set_pos(context, left, Vector((1, 0, 0)), normal=normal)
        self.handle.draw(context, render)


# Generic snap tool for line based archipack objects (fence, wall, maybe stair too)
gl_pts3d = []


class WallSnapManipulator(Manipulator):
    """
        np_station snap inspired manipulator
        Use prop1_name as string part index
        Use prop2_name as string identifier height property for placeholders

        Misnamed as it work for all line based archipack's
        primitives, currently wall and fences,
        but may also work with stairs (sharing same data structure)
    """
    def __init__(self, context, o, datablock, manipulator, handle_size):
        self.placeholder_part1 = GlPolygon((0.5, 0, 0, 0.2))
        self.placeholder_line1 = GlPolyline((0.5, 0, 0, 0.8))
        self.placeholder_part2 = GlPolygon((0.5, 0, 0, 0.2))
        self.placeholder_line2 = GlPolyline((0.5, 0, 0, 0.8))
        self.label = GlText()
        self.line = GlLine()
        self.handle = SquareHandle(handle_size, 1.2 * arrow_size, selectable=True)
        Manipulator.__init__(self, context, o, datablock, manipulator)

    def check_hover(self):
        self.handle.check_hover(self.mouse_pos)

    def mouse_press(self, context, event):
        global gl_pts3d
        if self.handle.hover:
            self.feedback.instructions(context, "Move / Snap", "Drag to move, use keyboard to input values", [
                ('CTRL', 'Snap'),
                ('MMBTN', 'Constraint to axis'),
                ('X Y', 'Constraint to axis'),
                ('SHIFT+Z', 'Constraint to xy plane'),
                ('RIGHTCLICK or ESC', 'exit without change')
                ])
            self.handle.hover = False
            self.handle.active = True
            self.o.select = True
            d = self.datablock
            idx = int(self.manipulator.prop1_name)
            # init placeholders 3d absolute world positionned points
            # between 2 and 3 points, the 2nd being moved by np_snap.delta
            takeloc, p1, side, normal = self.manipulator.get_pts(self.o.matrix_world)
            gl_pts3d = [None, takeloc, p1]
            if idx > 0:
                p0, right, side, normal = d.parts[idx - 1].manipulators[2].get_pts(self.o.matrix_world)
                gl_pts3d[0] = p0
            # enable placeholders drawing
            snap_point(takeloc, self.sp_draw, self.sp_callback,
                constraint_axis=(True, True, False))
            return True

        return False

    def mouse_release(self, context, event):
        self.check_hover()
        self.handle.active = False
        # False to callback manipulable_release
        return False

    def sp_callback(self, context, event, state, sp):
        """
            np station callback on moving, place, or cancel
        """
        global gl_pts3d

        # print("sp_callback %s" % (state))
        if state != 'RUNNING':

            # kill gl placeholder
            self.handle.active = False

            if state == 'SUCCESS':

                self.o.select = True
                # apply changes to wall
                d = self.datablock
                d.auto_update = False

                g = d.get_generator()

                # rotation relative to object
                rM = self.o.matrix_world.inverted().to_3x3()
                idx = int(self.manipulator.prop1_name)

                # new point in object space
                pt = g.segs[idx].lerp(0) + (rM * sp.delta).to_2d()
                da = 0

                # adjust size and rotation of segment before current
                if idx > 0:
                    w = g.segs[idx - 1]
                    part = d.parts[idx - 1]
                    dp = pt - w.p0

                    # adjust radius from distance between points..
                    # use p0-p1 distance as reference
                    if "C_" in part.type:
                        dw = (w.p1 - w.p0)
                        part.radius = part.radius / dw.length * dp.length
                        # angle pt - p0        - angle p0 p1
                        da = atan2(dp.y, dp.x) - atan2(dw.y, dw.x)
                    else:
                        part.length = dp.length
                        da = atan2(dp.y, dp.x) - w.straight(1).angle
                    a0 = part.a0 + da
                    if a0 > pi:
                        a0 -= 2 * pi
                    if a0 < -pi:
                        a0 += 2 * pi
                    # print("a0:%.4f part.a0:%.4f da:%.4f" % (a0, part.a0, da))
                    part.a0 = a0

                # adjust length of current segment
                w = g.segs[idx]
                part = d.parts[idx]
                dp = w.p1 - pt
                # adjust radius from distance between points..
                # use p0-p1 distance and angle as reference
                if "C_" in part.type:
                    dw = w.p1 - w.p0
                    part.radius = part.radius / dw.length * dp.length
                    da1 = atan2(dp.y, dp.x) - atan2(dw.y, dw.x)
                else:
                    part.length = dp.length
                    da1 = atan2(dp.y, dp.x) - w.straight(1).angle

                a0 = part.a0 + da1 - da
                if a0 > pi:
                    a0 -= 2 * pi
                if a0 < -pi:
                    a0 += 2 * pi
                # print("a0:%.4f part.a0:%.4f da:%.4f" % (a0, part.a0, da))
                part.a0 = a0

                # move object when point 0
                if idx == 0:
                    self.o.location += sp.delta

                # adjust rotation on both sides of next segment
                if idx + 1 < d.n_parts:

                    part = d.parts[idx + 1]
                    a0 = part.a0 - da1
                    if a0 > pi:
                        a0 -= 2 * pi
                    if a0 < -pi:
                        a0 += 2 * pi
                    part.a0 = a0
                d.auto_update = True

        return

    def sp_draw(self, sp, context):
        # draw wall placeholders
        global gl_pts3d
        if self.o is None:
            return
        z = self.get_value(self.datablock, self.manipulator.prop2_name)
        # print("z:%s, type:%s prop:%s" % (z, type(self.datablock).__name__, self.manipulator.prop2_name))
        p0 = gl_pts3d[1] + sp.delta
        p1 = gl_pts3d[2]
        self.placeholder_part1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.placeholder_line1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.placeholder_part1.draw(context, render=False)
        self.placeholder_line1.draw(context, render=False)
        if gl_pts3d[0] is not None:
            p0, p1 = gl_pts3d[0], p0
            self.placeholder_part2.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
            self.placeholder_line2.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
            self.placeholder_part2.draw(context, render=False)
            self.placeholder_line2.draw(context, render=False)

        self.line.p = gl_pts3d[1]
        self.line.v = sp.delta
        self.label.set_pos(context, self.line.length, self.line.lerp(0.5), self.line.v, normal=Vector((0, 0, 1)))
        self.line.draw(context, render=False)
        self.label.draw(context, render=False)

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle.active:
            # False here to pass_through
            # print("i'm able to pick up mouse move event while transform running")
            return False
        else:
            self.check_hover()
        return False

    def draw_callback(self, _self, context, render=False):
        left, right, side, normal = self.manipulator.get_pts(self.o.matrix_world)
        self.handle.set_pos(context, left, (left - right).normalized(), normal=normal)
        self.handle.draw(context, render)
        self.feedback.draw(context, render)


class CounterManipulator(Manipulator):
    """
        increase or decrease an integer step by step
        right on click to prevent misuse
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

    def mouse_press(self, context, event):
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

    def mouse_release(self, context, event):
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

    def draw_callback(self, _self, context, render=False):
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
        self.label.draw(context, render)
        self.handle_left.draw(context, render)
        self.handle_right.draw(context, render)


class SizeManipulator(Manipulator):

    def __init__(self, context, o, datablock, manipulator, handle_size):
        self.handle_left = TriHandle(handle_size, arrow_size)
        self.handle_right = TriHandle(handle_size, arrow_size, selectable=True)
        self.line_0 = GlLine()
        self.line_1 = GlLine()
        self.line_2 = GlLine()
        self.label = EditableText(handle_size, arrow_size, selectable=True)
        # self.label.label = 'S '
        Manipulator.__init__(self, context, o, datablock, manipulator)

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)
        self.label.check_hover(self.mouse_pos)

    def mouse_press(self, context, event):
        global gl_pts3d
        if self.handle_right.hover:
            self.active = True
            self.original_size = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.original_location = self.o.matrix_world.translation.copy()
            self.feedback.instructions(context, "Size", "Drag or Keyboard to modify size", [
                ('CTRL', 'Snap'), 
                ('SHIFT', 'Round'),
                ('RIGHTCLICK or ESC', 'cancel')
                ])
            left, right, side, dz = self.manipulator.get_pts(self.o.matrix_world)
            dx = (right - left).normalized()
            dy = dz.cross(dx)
            takemat = Matrix([
                [dx.x, dy.x, dz.x, right.x],
                [dx.y, dy.y, dz.y, right.y],
                [dx.z, dy.z, dz.z, right.z],
                [0, 0, 0, 1]
            ])
            t = context.scene.objects.find("Test")
            if t > -1:
                context.scene.objects[t].matrix_world = takemat
            gl_pts3d = [left, right]
            snap_point(takemat=takemat,
                draw=self.sp_draw,
                callback=self.sp_callback,
                constraint_axis=(True, False, False))
            self.handle_right.active = True
            """
            self.active = True
            self.original_size = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.feedback.instructions(context, "Size", "Drag to modify size", [
                ('ALT', 'Round value'), ('RIGHTCLICK or ESC', 'cancel')
                ])
            self.handle_right.active = True
            """
            return True
        if self.label.hover:
            self.feedback.instructions(context, "Size", "Use keyboard to modify size",
                [('ENTER', 'Validate'), ('RIGHTCLICK or ESC', 'cancel')])
            self.label.active = True
            self.keyboard_input_active = True
            return True
        return False

    def mouse_release(self, context, event):
        self.active = False
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

    def cancel(self, context, event):
        if self.active:
            self.mouse_release(context, event)
            self.set_value(context, self.datablock, self.manipulator.prop1_name, self.original_size)

    def keyboard_done(self, context, event, value):
        self.set_value(context, self.datablock, self.manipulator.prop1_name, value)
        self.label.active = False
        return True

    def keyboard_cancel(self, context, event):
        self.label.active = False
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

    def draw_callback(self, _self, context, render=False):
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
        if not self.keyboard_input_active:
            self.label_value = self.line_1.length
        self.label.set_pos(context, self.label_value, self.line_1.lerp(0.5), self.line_1.v, normal=normal)
        self.line_0.draw(context, render)
        self.line_1.draw(context, render)
        self.line_2.draw(context, render)
        self.handle_left.draw(context, render)
        self.handle_right.draw(context, render)
        self.label.draw(context, render)
        self.feedback.draw(context, render)

    def sp_draw(self, sp, context):
        global gl_pts3d
        if self.o is None:
            return
        p0 = gl_pts3d[0].copy()
        p1 = gl_pts3d[1].copy()
        p1 += sp.delta
        
        self.sp_update(context, p0, p1)
        """
        # manually update wall when manipulating wall child
        # as events are captured by snap
        if self.wall is not None:
            snap_helper = context.active_object
            context.scene.objects.active = self.wall
            self.wall.data.archipack_wall2[0].manipulable_manipulate(context, manipulator=self)
            context.scene.objects.active = snap_helper
        """
        return

    def sp_callback(self, context, event, state, sp):

        if state == 'SUCCESS':
            self.sp_draw(sp, context)
            self.mouse_release(context, event)

        if state == 'CANCEL':
            p0 = gl_pts3d[0].copy()
            p1 = gl_pts3d[1].copy()
            self.sp_update(context, p0, p1)
            self.mouse_release(context, event)

    def sp_update(self, context, p0, p1):
        length = (p0 - p1).length
        self.set_value(context, self.datablock, self.manipulator.prop1_name, length)
        
        

class SizeLocationManipulator(SizeManipulator):
    """
        Handle resizing by any of the boundaries
        of objects with centered pivots
        so when size change, object should move of the
        half of the change in the direction of change.

        Also take care of moving linked objects too
        Changing size is not necessary as link does
        allredy handle this and childs panels are
        updated by base object.
    """
    def __init__(self, context, o, datablock, manipulator, handle_size):
        SizeManipulator.__init__(self, context, o, datablock, manipulator, handle_size)
        self.handle_left.selectable = True

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)
        self.handle_left.check_hover(self.mouse_pos)
        self.label.check_hover(self.mouse_pos)

    def mouse_press(self, context, event):
        if self.handle_right.hover:
            self.active = True
            self.original_location = self.o.matrix_world.translation.copy()
            self.original_size = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.feedback.instructions(context, "Size", "Drag to modify size", [
                ('ALT', 'Round value'), ('RIGHTCLICK or ESC', 'cancel')
                ])
            self.handle_right.active = True
            return True
        if self.handle_left.hover:
            self.active = True
            self.original_location = self.o.matrix_world.translation.copy()
            self.original_size = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.feedback.instructions(context, "Size", "Drag to modify size", [
                ('ALT', 'Round value'), ('RIGHTCLICK or ESC', 'cancel')
                ])
            self.handle_left.active = True
            return True
        if self.label.hover:
            self.feedback.instructions(context, "Size", "Use keyboard to modify size",
                [('ENTER', 'Validate'), ('RIGHTCLICK or ESC', 'cancel')])
            self.label.active = True
            self.keyboard_input_active = True
            return True
        return False

    def mouse_release(self, context, event):
        self.active = False
        self.check_hover()
        self.handle_right.active = False
        self.handle_left.active = False
        if not self.keyboard_input_active:
            self.feedback.disable()
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active or self.handle_left.active:
            self.update(context, event)
            return True
        else:
            self.check_hover()
        return False

    def keyboard_done(self, context, event, value):
        """
        dl = value - self.get_value(self.datablock, self.manipulator.prop1_name)
        flip = self.get_value(self.datablock, 'flip')
        dl = 0.5 * dl
        if flip:
            dl = -dl
        self.move(context, self.manipulator.prop2_name, dl)
        """
        self.set_value(context, self.datablock, self.manipulator.prop1_name, value)
        # self.move_linked(context, self.manipulator.prop2_name, dl)
        self.label.active = False
        self.feedback.disable()
        return True

    def cancel(self, context, event):
        if self.active:
            self.mouse_release(context, event)
            # must move back to original location
            itM = self.o.matrix_world.inverted()
            dl = self.get_value(itM * self.original_location, self.manipulator.prop2_name)
            
            self.move(context, self.manipulator.prop2_name, dl)
            self.set_value(context, self.datablock, self.manipulator.prop1_name, self.original_size)
            self.move_linked(context, self.manipulator.prop2_name, dl)

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

        self.move(context, self.manipulator.prop2_name, dl)
        self.set_value(context, self.datablock, self.manipulator.prop1_name, length)
        self.move_linked(context, self.manipulator.prop2_name, dl)


class SnapSizeLocationManipulator(SizeLocationManipulator):
    """
        Snap aware extension of SizeLocationManipulator
        Handle resizing by any of the boundaries
        of objects with centered pivots
        so when size change, object should move of the
        half of the change in the direction of change.

        Also take care of moving linked objects too
        Changing size is not necessary as link does
        allredy handle this and childs panels are
        updated by base object.
        
        
    """
    def __init__(self, context, o, datablock, manipulator, handle_size):
        SizeLocationManipulator.__init__(self, context, o, datablock, manipulator, handle_size)
        if (context.active_object.data is not None and
                "archipack_wall2" in context.active_object.data):
            self.wall = context.active_object
        else:
            self.wall = None

    def mouse_press(self, context, event):
        global gl_pts3d
        if self.handle_right.hover:
            self.active = True
            self.original_size = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.original_location = self.o.matrix_world.translation.copy()
            self.feedback.instructions(context, "Size", "Drag or Keyboard to modify size", [
                ('CTRL', 'Snap'), 
                ('SHIFT', 'Round'),
                ('RIGHTCLICK or ESC', 'cancel')
                ])
            left, right, side, dz = self.manipulator.get_pts(self.o.matrix_world)
            dx = (right - left).normalized()
            dy = dz.cross(dx)
            takemat = Matrix([
                [dx.x, dy.x, dz.x, right.x],
                [dx.y, dy.y, dz.y, right.y],
                [dx.z, dy.z, dz.z, right.z],
                [0, 0, 0, 1]
            ])
            gl_pts3d = [left, right]
            snap_point(takemat=takemat,
            draw=self.sp_draw,
            callback=self.sp_callback,
            constraint_axis=(True, False, False))
        
            self.handle_right.active = True
            return True
        if self.handle_left.hover:
            self.active = True
            self.original_size = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.original_location = self.o.matrix_world.translation.copy()
            self.feedback.instructions(context, "Size", "Drag or Keyboard to modify size", [
                ('CTRL', 'Snap'), 
                ('SHIFT', 'Round'),
                ('RIGHTCLICK or ESC', 'cancel')
                ])
            left, right, side, dz = self.manipulator.get_pts(self.o.matrix_world)
            dx = (left - right).normalized()
            dy = dz.cross(dx)
            takemat = Matrix([
                [dx.x, dy.x, dz.x, right.x],
                [dx.y, dy.y, dz.y, right.y],
                [dx.z, dy.z, dz.z, right.z],
                [0, 0, 0, 1]
            ])
            gl_pts3d = [left, right]
            snap_point(takemat=takemat,
            draw=self.sp_draw,
            callback=self.sp_callback,
            constraint_axis=(True, False, False))
            self.handle_left.active = True
            return True
        if self.label.hover:
            self.feedback.instructions(context, "Size", "Use keyboard to modify size",
                [('ENTER', 'Validate'), ('RIGHTCLICK or ESC', 'cancel')])
            self.label.active = True
            self.keyboard_input_active = True
            return True
        return False

    def sp_draw(self, sp, context):
        global gl_pts3d
        if self.o is None:
            return
        p0 = gl_pts3d[0].copy()
        p1 = gl_pts3d[1].copy()
        if self.handle_right.active:
            p1 += sp.delta
        else:
            p0 += sp.delta
        self.sp_update(context, p0, p1)

        # manually update wall when manipulating wall child
        # as events are captured by snap
        if self.wall is not None:
            snap_helper = context.active_object
            context.scene.objects.active = self.wall
            self.wall.data.archipack_wall2[0].manipulable_manipulate(context, manipulator=self)
            context.scene.objects.active = snap_helper
        return

    def sp_callback(self, context, event, state, sp):

        if state == 'SUCCESS':
            self.sp_draw(sp, context)
            self.mouse_release(context, event)

        if state == 'CANCEL':
            p0 = gl_pts3d[0].copy()
            p1 = gl_pts3d[1].copy()
            self.sp_update(context, p0, p1)
            self.mouse_release(context, event)

    def sp_update(self, context, p0, p1):
        l0 = self.get_value(self.datablock, self.manipulator.prop1_name)
        length = (p0 - p1).length
        dp = length - l0
        if self.handle_left.active:
            dp = -dp
        dl = 0.5 * dp
        self.move(context, self.manipulator.prop2_name, dl)
        self.set_value(context, self.datablock, self.manipulator.prop1_name, length)
        self.move_linked(context, self.manipulator.prop2_name, dl)


class DeltaLocationManipulator(SizeManipulator):
    """
        Move a child window or door in wall segment
        not limited to this by the way
        
        @TODO:
            add snap feature on this one too
    """
    def __init__(self, context, o, datablock, manipulator, handle_size):
        SizeManipulator.__init__(self, context, o, datablock, manipulator, handle_size)
        self.label.label = ''
        self.feedback.instructions(context, "Move", "Drag to move", [
            ('ALT', 'Round value'), ('RIGHTCLICK or ESC', 'cancel')
            ])
        
    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)

    def mouse_press(self, context, event):
        if self.handle_right.hover:
            self.original_location = self.o.matrix_world.translation.copy()
            self.active = True
            self.feedback.enable()
            self.handle_right.active = True
            return True
        return False

    def mouse_release(self, context, event):
        self.check_hover()
        self.feedback.disable()
        self.active = False
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

    def cancel(self, context, event):
        if self.active:
            self.mouse_release(context, event)
            # must move back to original location
            itM = self.o.matrix_world.inverted()
            dl = self.get_value(itM * self.original_location, self.manipulator.prop1_name)
            self.move(context, self.manipulator.prop1_name, dl)
            
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

    def draw_callback(self, _self, context, render=False):
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
        self.handle_left.draw(context, render)
        self.handle_right.draw(context, render)
        self.feedback.draw(context)


class DumbSizeManipulator(SizeManipulator):
    """
        Show a size while not being editable
    """
    def __init__(self, context, o, datablock, manipulator, handle_size):
        SizeManipulator.__init__(self, context, o, datablock, manipulator, handle_size)
        self.handle_right.selectable = False
        self.label.selectable = False
        # self.label.label = 'Dumb '

    def mouse_move(self, context, event):
        return False


class AngleManipulator(Manipulator):
    """
        NOTE:
            There is a default shortcut to +5 and -5 on angles with left/right arrows
    
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
        self.label_a = EditableText(handle_size, arrow_size, selectable=True)
        self.label_a.unit = '°'
        Manipulator.__init__(self, context, o, datablock, manipulator)
        self.pts_mode = 'RADIUS'

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)
        self.label_a.check_hover(self.mouse_pos)

    def mouse_press(self, context, event):
        if self.handle_right.hover:
            self.active = True
            self.original_angle = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.feedback.instructions(context, "Angle", "Drag to modify angle", [
                ('ALT', 'Round value'), ('RIGHTCLICK or ESC', 'cancel')
                ])
            self.handle_right.active = True
            return True
        if self.label_a.hover:
            self.feedback.instructions(context, "Angle", "Use keyboard to modify angle",
                [('ENTER', 'validate'), ('RIGHTCLICK or ESC', 'cancel')])
            self.value_type = 'ROTATION'
            self.label_a.active = True
            self.keyboard_input_active = True
            return True
        return False

    def mouse_release(self, context, event):
        self.check_hover()
        self.handle_right.active = False
        self.active = False
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active:
            self.update(context, event)
            return True
        else:
            self.check_hover()
        return False

    def keyboard_done(self, context, event, value):
        self.set_value(context, self.datablock, self.manipulator.prop1_name, value / 180 * pi)
        self.label_a.active = False
        return True

    def keyboard_cancel(self, context, event):
        self.label_a.active = False
        return False

    def cancel(self, context, event):
        if self.active:
            self.mouse_release(context, event)
            self.set_value(context, self.datablock, self.manipulator.prop1_name, self.original_angle)

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
                da -= 2 * pi
            if da < -pi:
                da += 2 * pi
            # from there pi > da > -pi
            # print("a:%.4f da:%.4f a0:%.4f" % (atan2(v.y, v.x), da, self.line_0.angle))
            if da > pi:
                da = pi
            if da < -pi:
                da = -pi
            if event.alt:
                da = round(da / pi * 180, 0) / 180 * pi
            self.set_value(context, self.datablock, self.manipulator.prop1_name, da)

    def draw_callback(self, _self, context, render=False):
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
        label_value = self.arc.da / pi * 180
        if self.keyboard_input_active:
            label_value = self.label_value
        self.label_a.set_pos(context, label_value, self.arc.lerp(0.5), -self.line_0.v)
        self.arc.draw(context, render)
        self.line_0.draw(context, render)
        self.line_1.draw(context, render)
        self.handle_right.draw(context, render)
        self.handle_center.draw(context, render)
        self.label_a.draw(context, render)
        self.feedback.draw(context, render)


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
        self.label_a = EditableText(handle_size, arrow_size, selectable=True)
        self.label_r = EditableText(handle_size, arrow_size, selectable=False)
        self.label_a.unit = '°'
        Manipulator.__init__(self, context, o, datablock, manipulator)
        self.pts_mode = 'RADIUS'

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)
        self.label_a.check_hover(self.mouse_pos)

    def mouse_press(self, context, event):
        if self.handle_right.hover:
            self.active = True
            self.original_angle = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.feedback.instructions(context, "Angle", "Drag to modify angle", [
                ('ALT', 'Round value'), ('RIGHTCLICK or ESC', 'cancel')
                ])
            self.handle_right.active = True
            return True
        if self.label_a.hover:
            self.feedback.instructions(context, "Angle", "Use keyboard to modify angle",
                [('ENTER', 'validate'), ('RIGHTCLICK or ESC', 'cancel')])
            self.value_type = 'ROTATION'
            self.label_a.active = True
            self.keyboard_input_active = True
            return True
        if self.label_r.hover:
            self.feedback.instructions(context, "Radius", "Use keyboard to modify radius",
                [('ENTER', 'validate'), ('RIGHTCLICK or ESC', 'cancel')])
            self.value_type = 'LENGTH'
            self.label_r.active = True
            self.keyboard_input_active = True
            return True
        return False

    def mouse_release(self, context, event):
        self.check_hover()
        self.handle_right.active = False
        self.active = False
        return False

    def mouse_move(self, context, event):
        self.mouse_position(event)
        if self.handle_right.active:
            self.update(context, event)
            return True
        else:
            self.check_hover()
        return False

    def keyboard_done(self, context, event, value):
        self.set_value(context, self.datablock, self.manipulator.prop1_name, value / 180 * pi)
        self.label_a.active = False
        self.label_r.active = False
        return True

    def keyboard_cancel(self, context, event):
        self.label_a.active = False
        self.label_r.active = False
        return False

    def cancel(self, context, event):
        if self.active:
            self.mouse_release(context, event)
            self.set_value(context, self.datablock, self.manipulator.prop1_name, self.original_angle)

    def update(self, context, event):
        """
            NOTE:
                quadrant issue hidden somewhere
                mess with object's space and world space ?
        """
        # pt in world space here
        pt = self.get_pos3d(context)
        c = self.arc.c

        v = 2 * self.arc.r * (pt - c).normalized()
        v0 = c - v
        v1 = c + v
        p0, p1 = intersect_line_sphere(v0, v1, c, self.arc.r)

        if p0 is not None:
            # find nearest mouse intersection point
            if (p1 - pt).length < (p0 - pt).length:
                p0, p1 = p1, p0
            v = p0 - self.arc.c

            # a = atan2(v.y, v.x) - self.arc.sized_normal(0, 1).angle

            # a %= 2 * pi

            # print("a:%s >pi:%s <0:%s" % (a, a > pi, a < 0))

            # line0 from c to left (on radius side)
            da = atan2(v.y, v.x) - self.line_0.angle
            # da = atan2(v.y, v.x) - self.arc.sized_normal(0, 1).angle

            if abs(da) > pi:
                if da > 0:
                    da = pi
                else:
                    da = -pi

            # bottom and top points
            # s direction = "top"
            # take arc straight as reference as it wont never swap
            s = self.arc.tangeant(0, 1)
            top = s.lerp(1) - pt
            bottom = s.lerp(-1) - pt

            # left and right points
            # s direction = "right"
            n = s.sized_normal(0, 1)
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

    def draw_callback(self, _self, context, render=False):
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
        label_a_value = self.arc.da / pi * 180
        label_r_value = self.arc.r
        if self.keyboard_input_active:
            if self.value_type == 'LENGTH':
                label_r_value = self.label_value
            else:
                label_a_value = self.label_value
        self.label_a.set_pos(context, label_a_value, self.arc.lerp(0.5), -self.line_0.v)
        self.label_r.set_pos(context, label_r_value, self.line_0.lerp(0.5), self.line_0.v)
        self.arc.draw(context, render)
        self.line_0.draw(context, render)
        self.line_1.draw(context, render)
        self.handle_left.draw(context, render)
        self.handle_right.draw(context, render)
        self.handle_center.draw(context, render)
        self.label_r.draw(context, render)
        self.label_a.draw(context, render)
        self.feedback.draw(context, render)


class ArcAngleRadiusManipulator(ArcAngleManipulator):
    """
        Manipulate angle and radius of an arc
        when angle < 0 the arc center is on the left part of the circle
        when angle > 0 the arc center is on the right part of the circle
        bound to [-pi, pi]
    """

    def __init__(self, context, o, datablock, manipulator, handle_size):
        ArcAngleManipulator.__init__(self, context, o, datablock, manipulator, handle_size)
        self.handle_center = TriHandle(handle_size, arrow_size, selectable=True)
        self.label_r.selectable = True

    def check_hover(self):
        self.handle_right.check_hover(self.mouse_pos)
        self.handle_center.check_hover(self.mouse_pos)
        self.label_a.check_hover(self.mouse_pos)
        self.label_r.check_hover(self.mouse_pos)

    def mouse_press(self, context, event):
        if self.handle_right.hover:
            self.active = True
            self.original_angle = self.get_value(self.datablock, self.manipulator.prop1_name)
            self.feedback.instructions(context, "Angle", "Drag to modify angle", [
                ('ALT', 'Round value'), ('RIGHTCLICK or ESC', 'cancel')
                ])
            self.handle_right.active = True
            return True
        if self.handle_center.hover:
            self.active = True
            self.original_radius = self.get_value(self.datablock, self.manipulator.prop2_name)
            self.feedback.instructions(context, "Radius", "Drag to modify radius", [
                ('ALT', 'Round value'), ('RIGHTCLICK or ESC', 'cancel')
                ])
            self.handle_center.active = True
            return True
        if self.label_a.hover:
            self.feedback.instructions(context, "Angle", "Use keyboard to modify angle",
                [('ENTER', 'validate'), ('RIGHTCLICK or ESC', 'cancel')])
            self.value_type = 'ROTATION'
            self.label_a.active = True
            self.keyboard_input_active = True
            return True
        if self.label_r.hover:
            self.feedback.instructions(context, "Radius", "Use keyboard to modify radius",
                [('ENTER', 'validate'), ('RIGHTCLICK or ESC', 'cancel')])
            self.value_type = 'LENGTH'
            self.label_r.active = True
            self.keyboard_input_active = True
            return True
        return False

    def mouse_release(self, context, event):
        self.check_hover()
        self.active = False
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

    def keyboard_done(self, context, event, value):
        if self.value_type == 'LENGTH':
            self.set_value(context, self.datablock, self.manipulator.prop2_name, value)
            self.label_r.active = False
        else:
            self.set_value(context, self.datablock, self.manipulator.prop1_name, value / 180 * pi)
            self.label_a.active = False
        return True

    def update_radius(self, context, event):
        pt = self.get_pos3d(context)
        c = self.arc.c
        left = self.line_0.lerp(1)
        p, t = intersect_point_line(pt, c, left)
        radius = (left - p).length
        if event.alt:
            radius = round(radius, 1)
        self.set_value(context, self.datablock, self.manipulator.prop2_name, radius)

    def cancel(self, context, event):
        if self.handle_right.active:
            self.mouse_release(context, event)
            self.set_value(context, self.datablock, self.manipulator.prop1_name, self.original_angle)
        if self.handle_center.active:
            self.mouse_release(context, event)
            self.set_value(context, self.datablock, self.manipulator.prop2_name, self.original_radius)


# ------------------------------------------------------------------
# Define a single Manipulator Properties to store on object
# ------------------------------------------------------------------


# Allow registering manipulators classes
manipulators_class_lookup = {}


def register_manipulator(type_key, manipulator_class):
    if type_key in manipulators_class_lookup.keys():
        raise RuntimeError("Manipulator of type {} allready exists, unable to override".format(type_key))
    manipulators_class_lookup[type_key] = manipulator_class


class archipack_manipulator(PropertyGroup):
    """
        A property group to add to manipulable objects
        type_key: type of manipulator
        prop1_name = the property name of object to modify
        prop2_name = another property name of object to modify (eg: angle and radius)
        p0, p1, p2 3d Vectors as base points to represent manipulators on screen
        normal Vector normal of plane on with draw manipulator
    """
    type_key = StringProperty(default='SIZE')

    # How 3d points are stored in manipulators ?
    # SIZE = 2 absolute positionned and a scaling vector
    # RADIUS = 1 absolute positionned (center) and 2 relatives (sides)
    # POLYGON = 2 absolute positionned and a relative vector (for rect polygons)

    pts_mode = StringProperty(default='SIZE')
    prop1_name = StringProperty()
    prop2_name = StringProperty()
    p0 = FloatVectorProperty(subtype='XYZ')
    p1 = FloatVectorProperty(subtype='XYZ')
    p2 = FloatVectorProperty(subtype='XYZ')
    # allow orientation of manipulators by default on xy plane,
    # but may be used to constrain heights on local object space
    normal = FloatVectorProperty(subtype='XYZ', default=(0, 0, 1))

    def set_pts(self, pts, normal=None):
        """
            set 3d location of gl points (in object space)
            pts: array of 3 vectors 3d
            normal: optionnal vector 3d default to Z axis
        """
        pts = [Vector(p) for p in pts]
        self.p0, self.p1, self.p2 = pts
        if normal is not None:
            self.normal = Vector(normal)

    def get_pts(self, tM):
        """
            convert points from local to world absolute
            to draw them at the right place
            tM : object's world matrix
        """
        rM = tM.to_3x3()
        if self.pts_mode in ['SIZE', 'POLYGON']:
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

        if self.type_key not in manipulators_class_lookup.keys() or \
                not manipulators_class_lookup[self.type_key].poll(context):
            # RuntimeError is overkill but may be enabled for debug purposes
            # Silentely ignore allow skipping manipulators if / when deps as not meet
            # manip stack will simply be filled with None objects
            # raise RuntimeError("Manipulator of type {} not found".format(self.type_key))
            return None

        m = manipulators_class_lookup[self.type_key](context, o, datablock, self, handle_size)
        # points storage model as described upside
        self.pts_mode = m.pts_mode
        return m


# a global manipulator stack reference
# prevent Blender "ACCESS_VIOLATION" crashes
# use a dict to prevent potential
# collisions between many objects being in
# manipulate mode (at create time)
# use object names as loose keys
# NOTE : use app.drivers to reset before file load
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
    manip_stack = {}


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
    key_to_remove = []
    for key in manip_stack.keys():
        if bpy.data.objects.find(key) is None:
            key_to_remove.append(key)
    for key in key_to_remove:
        exit_stack(key)


def get_stack(key):
    """
        return reference to manipulator stack for given object
    """
    global manip_stack
    if key not in manip_stack.keys():
        manip_stack[key] = []
    return manip_stack[key]


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
            description="Flag manipulation state so we are able to toggle"
            )

    keymap = None

    def manipulable_getstack(self, context):
        o = context.active_object
        if o is not None:
            self.manip_stack = get_stack(o.name)

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
            should not be overriden
            as it provide all needed
            functionnality out of the box
        """
        # setup again when manipulators type change
        if self.manipulable_refresh:
            self.manipulable_refresh = False
            self.manipulable_setup(context)
            self.manipulate_mode = True

        context.area.tag_redraw()

        if self.keymap is None:
            self.keymap = Keymaps(context)

        if self.keymap.check(event, self.keymap.undo):
            # user feedback on undo by disabling manipulators
            self.manipulable_disable(context)
            return {'FINISHED'}

        # clean up manipulator on delete
        if self.keymap.check(event, self.keymap.delete):  # {'X'}:
            # @TODO:
            # for doors and windows, seek and destroy holes object if any
            # a dedicated delete method into those objects may be an option ?
            # A type check is required any way we choose
            #
            # Time for a generic archipack's datablock getter / filter into utils
            #
            # May also be implemented into nearly hidden "reference point"
            # to delete / duplicate / link duplicate / unlink of
            # a complete set of wall, doors and windows at once
            self.manipulable_disable(context)

            bpy.ops.object.delete('INVOKE_DEFAULT', use_global=False)
            return {'FINISHED'}

        for manipulator in self.manip_stack:
            # manipulator should return false on left mouse release
            # so proper release handler is called
            # and return true to call manipulate when required
            # print("manipulator:%s" % manipulator)
            if manipulator is not None and manipulator.modal(context, event):
                self.manipulable_manipulate(context, event, manipulator)
                return {'RUNNING_MODAL'}

        # allow any action on release
        if event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            self.manipulable_release(context)

        if event.type in {'RIGHTMOUSE', 'ESC'} and event.value == 'PRESS':
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

    def manipulable_manipulate(self, context, event, manipulator):
        """
            Override with action to do when a handle is active (pressed and mousemove)
        """
        return


def register():
    # Register default manipulators
    global manip_stack
    global manipulators_class_lookup
    manipulators_class_lookup = {}
    manip_stack = {}
    register_manipulator('SIZE', SizeManipulator)
    register_manipulator('SIZE_LOC', SizeLocationManipulator)
    register_manipulator('ANGLE', AngleManipulator)
    register_manipulator('ARC_ANGLE_RADIUS', ArcAngleRadiusManipulator)
    register_manipulator('COUNTER', CounterManipulator)
    register_manipulator('DUMB_SIZE', DumbSizeManipulator)
    register_manipulator('DELTA_LOC', DeltaLocationManipulator)
    # snap aware size loc
    register_manipulator('SNAP_SIZE_LOC', SnapSizeLocationManipulator)
    # register_manipulator('SNAP_POINT', SnapPointManipulator)
    # wall's line based object snap
    register_manipulator('WALL_SNAP', WallSnapManipulator)
    bpy.utils.register_class(archipack_manipulator)
    bpy.app.handlers.load_pre.append(empty_stack)


def unregister():
    global manip_stack
    global manipulators_class_lookup
    del manipulators_class_lookup
    del manip_stack
    bpy.utils.unregister_class(archipack_manipulator)
    bpy.app.handlers.load_pre.remove(empty_stack)
