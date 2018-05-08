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
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110- 1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

# ----------------------------------------------------------
# Author: Stephen Leger (s-leger)
#
# ----------------------------------------------------------
# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty, IntProperty, BoolProperty, BoolVectorProperty,
    CollectionProperty, FloatVectorProperty, EnumProperty, StringProperty
)
from mathutils import Vector, Matrix
from math import tan, sqrt, pi, sin, cos
from .bmesh_utils import BmeshEdit as bmed
from .panel import Panel as WindowPanel
from .archipack_handle import create_handle, window_handle_vertical_01, window_handle_vertical_02
# from .archipack_door_panel import ARCHIPACK_OT_select_parent
from .archipack_manipulator import Manipulable
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_gl import FeedbackPanel
from .archipack_object import (
    ArchipackObject, 
    ArchipackCreateTool, 
    ArchipackDrawTool,
    ArchipackObjectsManager
    )
from .archipack_segments import OpeningGenerator
from .archipack_keymaps import Keymaps
from .archipack_dimension import DimensionProvider


def update(self, context):
    self.update(context)


def update_childs(self, context):
    self.update(context, childs_only=True)


def update_portal(self, context):
    self.update_portal(context)


def set_cols(self, value):
    if self.n_cols != value:
        self.auto_update = False
        self._set_width(value)
        self.auto_update = True
        self.n_cols = value
    return None


def get_cols(self):
    return self.n_cols


class archipack_window_panelrow(PropertyGroup):
    width = FloatVectorProperty(
            name="Width",
            description="Panel width in percent of overall width",
            min=0.1,
            max=100.0,
            default=[
                50, 50, 50, 50, 50, 50, 50, 50,
                50, 50, 50, 50, 50, 50, 50, 50,
                50, 50, 50, 50, 50, 50, 50, 50,
                50, 50, 50, 50, 50, 50, 50
            ],
            size=31,
            update=update
            )
    fixed = BoolVectorProperty(
            name="Fixed",
            description="Fixed panel (generate a smallest frame)",
            default=[
                False, False, False, False, False, False, False, False,
                False, False, False, False, False, False, False, False,
                False, False, False, False, False, False, False, False,
                False, False, False, False, False, False, False, False
            ],
            size=32,
            update=update
            )
    cols = IntProperty(
            name="Panels",
            description="Number of panels on this row (up to 32)",
            min=1,
            max=32,
            default=2,
            get=get_cols, set=set_cols
            )
    n_cols = IntProperty(
            name="Panels",
            description="store number of panels, internal use only to avoid infinite recursion",
            min=1,
            max=32,
            default=2,
            update=update
            )
    height = FloatProperty(
            name="Height",
            min=0.1,
            default=1.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            name="auto_update",
            description="Disable auto update to avoid infinite recursion",
            default=True
            )

    def get_row(self, x, y):
        size = [Vector((x * self.width[w] / 100, y, 0)) for w in range(self.cols - 1)]
        sum_x = sum([s.x for s in size])
        size.append(Vector((x - sum_x, y, 0)))
        origin = []
        pivot = []
        ttl = 0
        xh = x / 2
        n_center = len(size) / 2
        for i, sx in enumerate(size):
            ttl += sx.x
            if i < n_center:
                # pivot left
                origin.append(Vector((ttl - xh - sx.x, 0)))
                pivot.append(1)
            else:
                # pivot right
                origin.append(Vector((ttl - xh, 0)))
                pivot.append(-1)
        return size, origin, pivot

    def _set_width(self, cols):
        width = 100 / cols
        for i in range(cols - 1):
            self.width[i] = width

    def find_datablock_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_window.datablock(o)
            if props:
                for row in props.rows:
                    if row == self:
                        return props
        return None

    def update(self, context):
        if self.auto_update:
            props = self.find_datablock_in_selection(context)
            if props is not None:
                props.update(context, childs_only=False)

    def draw(self, layout, context, last_row):
        # store parent at runtime to trigger update on parent
        row = layout.row()
        row.prop(self, "cols")
        row = layout.row()
        if not last_row:
            row.prop(self, "height")
        for i in range(self.cols - 1):
            row = layout.row()
            row.prop(self, "width", text="col " + str(i + 1), index=i)
            row.prop(self, "fixed", text="fixed", index=i)
        row = layout.row()
        row.label(text="col " + str(self.cols))
        row.prop(self, "fixed", text="fixed", index=(self.cols - 1))


class archipack_window_panel(ArchipackObject, PropertyGroup):
    # @TODO:
    # No need to store params here
    # pass those as update() params
    center = FloatVectorProperty(
            subtype='XYZ'
            )
    origin = FloatVectorProperty(
            subtype='XYZ'
            )
    size = FloatVectorProperty(
            subtype='XYZ'
            )
    radius = FloatVectorProperty(
            subtype='XYZ'
            )
    angle_y = FloatProperty(
            name='angle',
            unit='ROTATION',
            subtype='ANGLE',
            min=-1.5, max=1.5,
            default=0, precision=2,
            description='angle'
            )
    frame_y = FloatProperty(
            name='Depth',
            min=0,
            default=0.06, precision=2,
            unit='LENGTH', subtype='DISTANCE',
            description='frame depth'
            )
    frame_x = FloatProperty(
            name='Width',
            min=0,
            default=0.06, precision=2,
            unit='LENGTH', subtype='DISTANCE',
            description='frame width'
            )
    curve_steps = IntProperty(
            name="curve steps",
            min=1,
            max=128,
            default=1
            )
    shape = EnumProperty(
            name='Shape',
            items=(
                ('RECTANGLE', 'Rectangle', '', 0),
                ('ROUND', 'Top Round', '', 1),
                ('ELLIPSIS', 'Top elliptic', '', 2),
                ('QUADRI', 'Top oblique', '', 3),
                ('CIRCLE', 'Full circle', '', 4)
                ),
            default='RECTANGLE'
            )
    pivot = FloatProperty(
            name='pivot',
            min=-1, max=1,
            default=-1, precision=2,
            description='pivot'
            )
    side_material = IntProperty(
            name="side material",
            min=0,
            max=2,
            default=0
            )
    handle = EnumProperty(
            name='Shape',
            items=(
                ('NONE', 'No handle', '', 0),
                ('INSIDE', 'Inside', '', 1),
                ('BOTH', 'Inside and outside', '', 2)
                ),
            default='NONE'
            )
    handle_model = IntProperty(
            name="handle model",
            default=1,
            min=1,
            max=2
            )
    handle_altitude = FloatProperty(
            name='handle altitude',
            min=0,
            default=0.2, precision=2,
            unit='LENGTH', subtype='DISTANCE',
            description='handle altitude'
            )
    fixed = BoolProperty(
            name="Fixed",
            default=False
            )
    enable_glass = BoolProperty(
            name="Enable glass",
            default=True
            )

    @property
    def window(self):
        verre = 0.005
        chanfer = 0.004
        x0 = 0
        x1 = self.frame_x
        x2 = 0.75 * self.frame_x
        x3 = chanfer
        y0 = -self.frame_y
        y1 = 0
        y2 = -0.5 * self.frame_y
        y3 = -chanfer
        y4 = chanfer - self.frame_y

        if self.fixed:
            # profil carre avec support pour verre
            # p ______       y1
            # / |      y3
            # |       |___
            # x       |___   y2  verre
            # |       |      y4
            #  \______|      y0
            # x0 x3   x1
            #
            x1 = 0.5 * self.frame_x
            y1 = -0.45 * self.frame_y
            y3 = y1 - chanfer
            y4 = chanfer + y0
            y2 = (y0 + y2) / 2

            side_cap_front = -1
            side_cap_back = -1

            if self.enable_glass:
                side_cap_front = 6
                side_cap_back = 7

            return WindowPanel(
                True,  # closed
                [1, 0, 0, 0, 1, 2, 2, 2, 2],  # x index
                [x0, x3, x1],
                [y0, y4, y2, y3, y1, y1, y2 + verre, y2 - verre, y0],
                [0, 0, 1, 1, 1, 1, 0, 0, 0],  # materials
                side_cap_front=side_cap_front,
                side_cap_back=side_cap_back      # cap index
                )
        else:
            # profil avec chanfrein et joint et support pour verre
            # p ____         y1    inside
            # /     |_       y3
            # |       |___
            # x       |___   y2  verre
            # |      _|      y4
            #  \____|        y0
            # x0 x3 x2 x1          outside
            if self.side_material == 0:
                materials = [0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
            elif self.side_material == 1:
                # rail window exterior
                materials = [0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
            else:
                # rail window interior
                materials = [0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]

            side_cap_front = -1
            side_cap_back = -1

            if self.enable_glass:
                side_cap_front = 8
                side_cap_back = 9

            return WindowPanel(
                True,            # closed shape
                [1, 0, 0, 0, 1, 2, 2, 3, 3, 3, 3, 2, 2],     # x index
                [x0, x3, x2, x1],     # unique x positions
                [y0, y4, y2, y3, y1, y1, y3, y3, y2 + verre, y2 - verre, y4, y4, y0],
                materials,     # materials
                side_cap_front=side_cap_front,
                side_cap_back=side_cap_back      # cap index
                )

    @property
    def verts(self):
        offset = Vector((0, 0, 0))
        return self.window.vertices(self.curve_steps, offset, self.center, self.origin, self.size,
            self.radius, self.angle_y, self.pivot, shape_z=None, path_type=self.shape)

    @property
    def faces(self):
        return self.window.faces(self.curve_steps, path_type=self.shape)

    @property
    def matids(self):
        return self.window.mat(self.curve_steps, 2, 2, path_type=self.shape)

    @property
    def uvs(self):
        return self.window.uv(self.curve_steps, self.center, self.origin, self.size,
            self.radius, self.angle_y, self.pivot, 0, self.frame_x, path_type=self.shape)

    def find_handle(self, o):
        for child in o.children:
            if 'archipack_handle' in child:
                return child
        return None

    def update_handle(self, context, o):
        handle = self.find_handle(o)
        if handle is None:
            m = bpy.data.meshes.new("Handle")
            handle = create_handle(context, o, m)
            # MaterialUtils.add_handle_materials(handle)
        if self.handle_model == 1:
            verts, faces = window_handle_vertical_01(1)
        else:
            verts, faces = window_handle_vertical_02(1)
        handle.location = (self.pivot * (self.size.x - 0.4 * self.frame_x), 0, self.handle_altitude)
        bmed.buildmesh(context, handle, verts, faces)

    def remove_handle(self, context, o):
        handle = self.find_handle(o)
        if handle is not None:
            self.delete_object(context, handle)

    def update(self, context):

        o = self.find_in_selection(context)

        if o is None:
            return

        if self.handle == 'NONE':
            self.remove_handle(context, o)
        else:
            self.update_handle(context, o)

        bmed.buildmesh(context, o, self.verts, self.faces, self.matids, self.uvs)

        self.restore_context(context)


class archipack_window_shutter(ArchipackObject, PropertyGroup):
    center = FloatVectorProperty(
            subtype='XYZ'
            )
    origin = FloatVectorProperty(
            subtype='XYZ'
            )
    size = FloatVectorProperty(
            subtype='XYZ'
            )
    radius = FloatVectorProperty(
            subtype='XYZ'
            )
    angle_y = FloatProperty(
            name='angle',
            unit='ROTATION',
            subtype='ANGLE',
            min=-1.5, max=1.5,
            default=0, precision=2,
            description='angle'
            )
    depth = FloatProperty(
            name='Depth',
            min=0,
            default=0.06, precision=2,
            unit='LENGTH', subtype='DISTANCE',
            description='frame depth'
            )
    border = FloatProperty(
            name='Width',
            min=0,
            default=0.06, precision=2,
            unit='LENGTH', subtype='DISTANCE',
            description='frame width'
            )
    curve_steps = IntProperty(
            name="curve steps",
            min=1,
            max=128,
            default=1
            )
    shape = EnumProperty(
            name='Shape',
            items=(
                ('RECTANGLE', 'Rectangle', '', 0),
                ('ROUND', 'Top Round', '', 1),
                ('ELLIPSIS', 'Top elliptic', '', 2),
                ('QUADRI', 'Top oblique', '', 3),
                ('CIRCLE', 'Full circle', '', 4)
                ),
            default='RECTANGLE'
            )
    pivot = FloatProperty(
            name='pivot',
            min=-1, max=1,
            default=-1, precision=2,
            description='pivot'
            )
    offset = FloatProperty(
            name='offset',
            default=0, precision=2,
            description='x offset'
            )
    hinge_enable = BoolProperty(
            name="Hinge",
            default=False
            )
    hinge_count = IntProperty(
            name="#Hinge",
            min=0,
            max=4,
            default=2
            )
    hinge_space = FloatProperty(
            name='space',
            default=0, precision=2,
            description='Vertical space for hinges'
            )
    hinge_size = FloatProperty(
            name='size',
            default=0.03, min=0.001, precision=2,
            description='hinge vertical size'
            )

    @property
    def shutter(self):

        chanfer = 0.004
        border = self.border
        spacing = 0.75 * self.border
        x0 = 0
        x1 = border - 0.5 * spacing
        x3 = chanfer
        w = 0.5 * self.depth
        # offset pivot point on outside part
        y0 = 0
        y1 = y0 + w
        y2 = y1 - 0.5 * w
        y3 = y1 - chanfer
        y4 = y0 + chanfer

        # profil
        # p ______   y1
        #  /         y3
        # |
        # x          y2
        # |          y4
        #  \______   y0
        # x0 x3   x1
        #

        side = WindowPanel(
            False,  # closed
            [2, 1, 0, 0, 0, 1, 2],  # x index
            [x0, x3, x1],
            [y0, y0, y4, y2, y3, y1, y1],
            [5, 5, 5, 5, 5, 5, 5],  # materials
            closed_path=True,    #
            subdiv_x=0,
            subdiv_y=1
            )

        #     /   y2-y3
        #  __/    y1-y0
        #   x2 x3
        x2 = 0.5 * spacing
        x3 = x2 + chanfer
        y2 = y1 - chanfer
        y3 = y0 + chanfer

        face = WindowPanel(
            False,              # profil closed
            [0, 1, 2],          # x index
            [0, x2, x3],
            [y1, y1, y2],
            [5, 5, 5],          # material index
            side_cap_front=2,   # cap index
            closed_path=True
            )

        back = WindowPanel(
            False,              # profil closed
            [0, 1, 2],          # x index
            [x3, x2, 0],
            [y3, y0, y0],
            [5, 5, 5],          # material index
            side_cap_back=0,    # cap index
            closed_path=True
            )

        return side, face, back

    def hinge(self, altitude, verts):

        # panel chanfer
        chanfer = 0.004

        seg = 12
        deg = 2 * pi / seg
        radius = self.hinge_size / 6
        x = 0
        y = 0
        z = altitude

        d = (self.offset + self.pivot * chanfer) / radius
        tM = Matrix([
            [radius, 0, 0, x],
            [0, radius, 0, y],
            [0, 0, self.hinge_size, z - 0.5 * self.hinge_size],
            [0, 0, 0, 1]
        ])
        if self.pivot < 0:
            verts.extend([tM * Vector((sin(deg * a), cos(deg * a), 0)) for a in range(seg - 2)])
            verts.extend([
                tM * Vector((d, cos(deg * (seg - 3)), 0)),
                tM * Vector((d, 1, 0)),
            ])
        else:
            verts.extend([tM * Vector((sin(deg * a), cos(deg * a), 0)) for a in range(3, seg + 1)])
            verts.extend([
                tM * Vector((d, 1, 0)),
                tM * Vector((d, cos(deg * (seg - 3)), 0)),
            ])

        if self.pivot < 0:
            verts.extend([tM * Vector((sin(deg * a), cos(deg * a), 1)) for a in range(seg - 2)])
            verts.extend([
                tM * Vector((d, cos(deg * (seg - 3)), 1)),
                tM * Vector((d, 1, 1)),
            ])
        else:
            verts.extend([tM * Vector((sin(deg * a), cos(deg * a), 1)) for a in range(3, seg + 1)])
            verts.extend([
                tM * Vector((d, 1, 1)),
                tM * Vector((d, cos(deg * (seg - 3)), 1)),
            ])

    @property
    def verts(self):

        side, face, back = self.shutter
        border = self.border
        spacing = 0.75 * self.border

        x1 = border - 0.5 * spacing
        offset = Vector((self.offset, 0, 0))
        verts = side.vertices(self.curve_steps, offset, self.center, self.origin, self.size,
            self.radius, self.angle_y, self.pivot, shape_z=None, path_type=self.shape)

        p_radius = self.radius.copy()
        p_radius.x -= x1
        p_radius.y -= x1

        p_size = Vector((self.size.x - 2 * x1, (self.size.y - 2 * x1) / 2, 0))

        for j in range(2):
            if j < 1:
                shape = 'RECTANGLE'
            else:
                shape = self.shape
            offset = Vector((
                self.offset + self.pivot * x1,
                p_size.y * j + x1,
                0))

            origin = Vector((
                self.origin.x + self.pivot * x1,
                p_size.y * j + x1,
                0))

            verts += face.vertices(self.curve_steps, offset, self.center, origin,
                p_size, p_radius, self.angle_y, self.pivot, shape_z=None, path_type=shape)
            verts += back.vertices(self.curve_steps, offset, self.center, origin,
                p_size, p_radius, self.angle_y, self.pivot, shape_z=None, path_type=shape)

        if self.hinge_enable:
            z0 = 0.15
            dz = (self.hinge_space - 2 * z0) / (self.hinge_count - 1)
            for j in range(self.hinge_count):
                self.hinge(z0 + dz * j, verts)

        return verts

    @property
    def faces(self):

        side, face, back = self.shutter

        faces = side.faces(self.curve_steps, path_type=self.shape)
        faces_offset = side.n_verts(self.curve_steps, path_type=self.shape)

        for j in range(2):
            if j < 1:
                shape = 'RECTANGLE'
            else:
                shape = self.shape
            faces += face.faces(self.curve_steps, path_type=shape, offset=faces_offset)
            faces_offset += face.n_verts(self.curve_steps, path_type=shape)
            faces += back.faces(self.curve_steps, path_type=shape, offset=faces_offset)
            faces_offset += back.n_verts(self.curve_steps, path_type=shape)

        if self.hinge_enable:
            seg = 12
            for j in range(self.hinge_count):
                faces.append(tuple([faces_offset + i + seg for i in range(seg - 1, -1, -1)]))
                faces.append(tuple([faces_offset + i for i in range(seg)]))
                faces.extend([tuple([faces_offset + i + f for f in (1, 0, seg, seg + 1)]) for i in range(seg - 1)])
                faces.append((
                    faces_offset,
                    faces_offset + seg - 1,
                    faces_offset + 2 * seg - 1,
                    faces_offset + seg
                    ))

                faces_offset += 2 * seg

        return faces

    @property
    def uvs(self):

        side, face, back = self.shutter

        border = self.border
        spacing = 0.75 * self.border
        x1 = border - 0.5 * spacing

        uvs = side.uv(self.curve_steps,
            self.center,
            self.origin,
            self.size,
            self.radius,
            self.angle_y,
            self.pivot,
            self.border, 0,
            path_type=self.shape)

        p_radius = self.radius.copy()
        p_radius.x -= x1
        p_radius.y -= x1
        p_size = Vector((self.size.x - 2 * x1, (self.size.y - 2 * x1) / 2, 0))

        for j in range(2):
            if j < 1:
                shape = 'RECTANGLE'
            else:
                shape = self.shape
            origin = Vector((
                self.origin.x + self.pivot * x1,
                p_size.y * j + x1,
                0))
            uvs += face.uv(self.curve_steps, self.center, origin, p_size,
                p_radius, self.angle_y, self.pivot, 0, 0, path_type=shape)
            uvs += back.uv(self.curve_steps, self.center, origin, p_size,
                p_radius, self.angle_y, self.pivot, 0, 0, path_type=shape)

        if self.hinge_enable:
            seg = 12
            deg = 2 * pi / seg
            radius = 0.005
            x = 0
            y = 0
            z = 0
            size = 0.04
            tM = Matrix([
                [radius, 0, 0, x],
                [0, radius, 0, y],
                [0, 0, size, z],
                [0, 0, 0, 1]
            ])

            for j in range(self.hinge_count):
                uvs.append(tuple([(tM * Vector((sin(deg * a), cos(deg * a), 0))).to_2d() for a in range(seg)]))
                uvs.append(tuple([(tM * Vector((sin(deg * a), cos(deg * a), 0))).to_2d() for a in range(seg)]))
                uvs.extend([[(0, 0), (0, 1), (1, 1), (1, 0)] for i in range(seg)])

        return uvs

    @property
    def matids(self):

        side, face, back = self.shutter
        mat = side.mat(self.curve_steps, 5, 5, path_type=self.shape)
        for j in range(2):
            if j < 1:
                shape = 'RECTANGLE'
            else:
                shape = self.shape
            mat += face.mat(self.curve_steps, 5, 5, path_type=shape)
            mat += back.mat(self.curve_steps, 5, 5, path_type=shape)

        if self.hinge_enable:
            for j in range(self.hinge_count):
                seg = 12
                mat.extend([3, 3])
                mat.extend([3 for i in range(seg)])

        return mat

    def update(self, context):

        o = self.find_in_selection(context)

        if o is None:
            return

        bmed.buildmesh(context, o, self.verts, self.faces, self.matids, self.uvs)

        self.restore_context(context)


class archipack_window(ArchipackObject, Manipulable, DimensionProvider, PropertyGroup):
    x = FloatProperty(
            name='Width',
            min=0.1,
            default=100.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Width', update=update
            )
    y = FloatProperty(
            name='Depth',
            min=0.05,
            default=0.20, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Depth', update=update,
            )
    z = FloatProperty(
            name='Height',
            min=0.1,
            default=1.2, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='height', update=update,
            )
    angle_y = FloatProperty(
            name='Angle',
            unit='ROTATION',
            subtype='ANGLE',
            min=-1.5, max=1.5,
            default=0, precision=2,
            description='angle', update=update,
            )
    radius = FloatProperty(
            name='Radius',
            min=0.1,
            default=2.5, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='radius', update=update,
            )
    elipsis_b = FloatProperty(
            name='Ellipsis',
            min=0.1,
            default=0.5, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='ellipsis vertical size', update=update,
            )
    altitude = FloatProperty(
            name='Altitude',
            default=1.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='altitude', update=update,
            )
    offset = FloatProperty(
            name='Offset',
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='offset', update=update,
            )
    frame_y = FloatProperty(
            name='Depth',
            min=0,
            default=0.06, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='frame depth', update=update,
            )
    frame_x = FloatProperty(
            name='Width',
            min=0,
            default=0.06, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='frame width', update=update,
            )
    frame_overflow = FloatProperty(
            name='Overflow lateral',
            min=0,
            default=0.06, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='frame width', update=update,
            )
    overflow_out = FloatProperty(
            name='Finishing thickness',
            min=0,
            default=0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Thickness of outside wall finishing', update=update,
            )
    panel_x = FloatProperty(
            name='Width',
            min=0,
            default=0.06, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='panel width', update=update,
            )
    panel_y = FloatProperty(
            name='Depth',
            min=0,
            default=0.06, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='panel depth', update=update,
            )
    out_frame = BoolProperty(
            name="Out frame",
            default=False, update=update,
            )
    out_frame_y = FloatProperty(
            name='Side depth',
            min=0.001,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='frame side depth', update=update,
            )
    out_frame_y2 = FloatProperty(
            name='Front depth',
            min=0.001,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='frame front depth', update=update,
            )
    out_frame_x = FloatProperty(
            name='Front Width',
            min=0.0,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='frame width set to 0 disable front frame', update=update,
            )
    out_frame_offset = FloatProperty(
            name='Offset',
            min=0.0,
            default=0.0, precision=3, step=0.1,
            unit='LENGTH', subtype='DISTANCE',
            description='frame offset', update=update,
            )
    out_tablet_enable = BoolProperty(
            name="Sill out",
            default=True, update=update,
            )
    out_tablet_x = FloatProperty(
            name='Width',
            min=0.0,
            default=0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='tablet width', update=update,
            )
    out_tablet_y = FloatProperty(
            name='Depth',
            min=0.001,
            default=0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='tablet depth', update=update,
            )
    out_tablet_z = FloatProperty(
            name='Height',
            min=0.001,
            default=0.03, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='tablet height', update=update,
            )
    in_tablet_enable = BoolProperty(
            name="Sill in",
            default=True, update=update,
            )
    in_tablet_x = FloatProperty(
            name='Width',
            min=0.0,
            default=0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='tablet width', update=update,
            )
    in_tablet_y = FloatProperty(
            name='Depth',
            min=0.001,
            default=0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='tablet depth', update=update,
            )
    in_tablet_z = FloatProperty(
            name='Height',
            min=0.001,
            default=0.03, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='tablet height', update=update,
            )
    blind_inside = BoolProperty(
            default=False,
            name="Blind inside",
            description="Generate a blind inside",
            update=update
            )
    blind_outside = BoolProperty(
            default=False,
            name="Blind outside",
            description="Generate a blind outside",
            update=update
            )
    # internal blind, keep for compatibility
    blind_enable = BoolProperty(
            name="Blind",
            default=False, update=update,
            )
    blind_y = FloatProperty(
            name='Depth',
            min=0.001,
            default=0.002, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Store depth', update=update,
            )
    blind_z = FloatProperty(
            name='Height',
            min=0.001,
            default=0.03, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Store height', update=update,
            )
    blind_open = FloatProperty(
            name='Open',
            min=0.0, max=100,
            default=80, precision=1,
            subtype='PERCENTAGE',
            description="Blind open", 
            update=update,
            )
    rows = CollectionProperty(type=archipack_window_panelrow)
    n_rows = IntProperty(
            name="Number of rows",
            min=1,
            max=32,
            default=1, update=update,
            )
    curve_steps = IntProperty(
            name="Steps",
            min=6,
            max=128,
            default=16, update=update,
            )
    hole_outside_mat = IntProperty(
            name="Outside",
            description="Material index of wall for outside part of the hole",
            min=0,
            max=128,
            default=0, update=update,
            )
    hole_inside_mat = IntProperty(
            name="Inside",
            description="Material index of wall for inside part of the hole",
            min=0,
            max=128,
            default=1, update=update,
            )
    window_shape = EnumProperty(
            name='Shape',
            items=(
                ('RECTANGLE', 'Rectangle', '', 0),
                ('ROUND', 'Top Round', '', 1),
                ('ELLIPSIS', 'Top elliptic', '', 2),
                ('QUADRI', 'Top oblique', '', 3),
                ('CIRCLE', 'Full circle', '', 4)
                ),
            default='RECTANGLE', update=update,
            )
    window_type = EnumProperty(
            name='Type',
            items=(
                ('FLAT', 'Swing window', '', 0),
                ('RAIL', 'Rail window', '', 1)
                ),
            default='FLAT', update=update,
            )
    enable_glass = BoolProperty(
            name="Enable glass",
            default=True,
            update=update
            )
    warning = BoolProperty(
            options={'SKIP_SAVE'},
            name="Warning",
            default=False
            )
    handle_enable = BoolProperty(
            name='Handle',
            default=True, update=update_childs,
            )
    handle_altitude = FloatProperty(
            name="Altitude",
            min=0,
            default=1.4, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='handle altitude', update=update_childs,
            )
    shutter_enable = BoolProperty(
            name="Shutters",
            default=False, update=update,
            )
    shutter_left = IntProperty(
            name="#Left",
            default=1,
            update=update
            )
    shutter_right = IntProperty(
            name="#Right",
            default=1,
            update=update
            )
    shutter_border = FloatProperty(
            name='Border',
            min=0,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Shutter panels borders',
            update=update
            )
    shutter_depth = FloatProperty(
            name='Depth',
            min=0.01,
            default=0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Shutter panels depth',
            update=update
            )
    shutter_hinge = FloatProperty(
            name='Hinge',
            min=0.001,
            default=0.03, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Shutter hinge size',
            update=update
            )
    hole_margin = FloatProperty(
            name='Hole margin',
            min=0.0,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='how much hole surround wall'
            )
    flip = BoolProperty(
            default=False,
            update=update,
            description='flip outside and outside material of hole'
            )

    # layout related
    display_detail = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    display_panels = BoolProperty(
            options={'SKIP_SAVE'},
            default=True
            )
    display_materials = BoolProperty(
            options={'SKIP_SAVE'},
            default=True
            )

    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update
            )
    portal = BoolProperty(
            default=False,
            name="Portal",
            description="Makes this window a light portal, enabling better photon gathering during render",
            update=update
            )

    @property
    def shape(self):
        if self.window_type == 'RAIL':
            return 'RECTANGLE'
        else:
            return self.window_shape

    @property
    def window(self):
        # Flat window frame profil
        #  ___        y1
        # |   |__
        # |      |    y2
        # |______|    y0
        #
        x0 = 0
        x1 = -x0 - self.frame_x
        x2 = x0 + 0.5 * self.frame_x
        y0 = 0.5 * self.y - self.offset
        y2 = y0 + 0.5 * self.frame_y
        
        if self.window_type == 'FLAT':
            y1 = y0 + self.frame_y
            return WindowPanel(
                True,     # closed
                [0, 0, 1, 1, 2, 2],     # x index
                [x1, x0, x2],
                [y0, y1, y1, y2, y2, y0],
                [1, 1, 1, 1, 0, 0]  # material index
                )
        else:
            # Rail window frame profil
            #  ________       y1
            # |      __|      y5
            # |     |__       y4
            # |      __|      y3
            # |     |____
            # |          |    y2
            # |__________|    y0
            # -x1   0 x3 x2
            x2 = x0 + 0.5 * self.panel_x
            x3 = x0 + 0.2 * self.panel_x
            y1 = y0 + max(self.frame_y, 2.05 * self.panel_y)
            yc = y0 + 0.5 * (y1 - y0)
            y2 = yc - 0.995 * self.panel_y
            y3 = yc - 0.0175 * self.panel_y
            y4 = yc + 0.0175 * self.panel_y
            y5 = yc + 0.995 * self.panel_y

            return WindowPanel(
                True,     # closed
                [0, 0, 2, 2, 1, 1, 2, 2, 1, 1, 3, 3],     # x index
                [x1, x0, x3, x2],
                [y0, y1, y1, y5, y5, y4, y4, y3, y3, y2, y2, y0],
                [1, 1, 1, 3, 1, 3, 3, 3, 0, 3, 0, 0]  # material index
                )

    @property
    def hole(self):
        # profil percement                          ____
        #   _____  y_inside          vertical ___|      x1
        #  |
        #  |__     y0                outside   ___
        #     |___ y_outside                       |____ x1-shape_z     inside
        # -x1 x0
        y0 = 0.5 * self.y - self.offset
        x1 = self.frame_x     # sur-largeur percement interieur
        y_inside = 0.5 * self.y + self.hole_margin     # outside wall

        x0 = self._overflow - 0.001

        if self.out_frame:
            x0 -= min(self.frame_x, self.out_frame_y + self.out_frame_offset)

        outside_mat = self.hole_outside_mat
        inside_mat = self.hole_inside_mat
        # if self.flip:
        #    outside_mat, inside_mat = inside_mat, outside_mat

        y_outside = -y_inside           # inside wall
        y_outside -= self.overflow_out
        if self.frame_overflow > 0 or (
                (self.out_frame and self.out_frame_offset > 0) or
                self.out_tablet_enable):
            return WindowPanel(
                False,     # closed
                [1, 1, 0, 0],     # x index
                [-x1, x0],
                [y_outside, y0, y0, y_inside],
                [outside_mat, outside_mat, inside_mat],     # material index
                side_cap_front=3,     # cap index
                side_cap_back=0
                )
        else:
            # Hole without overflow
            return WindowPanel(
                False,     # closed
                [0, 0, 0],     # x index
                [x0],
                [y_outside, y0, y_inside],
                [outside_mat, inside_mat],     # material index
                side_cap_front=2,     # cap index
                side_cap_back=0
                )

    @property
    def inside_hole(self):
        # inside part of hole to setup line for floor
        # profil percement                          ____
        #   _____  y_inside          vertical ___|      x1
        #  |
        #  |__     y0                outside   ___
        #     |___ y_outside                       |____ x1-shape_z     inside
        # -x1 x0
        y0 = 0.5 * self.y - self.offset
        x1 = self.frame_x     # sur-largeur percement interieur
        y_inside = 0.5 * self.y + self.hole_margin     # outside wall
        if self.window_type == 'FLAT':
            y1 = y0 + self.frame_y
        else:
            y1 = y0 + max(self.frame_y, 2.05 * self.panel_y)
        inside_mat = self.hole_inside_mat
        # if self.flip:
        #    outside_mat, inside_mat = inside_mat, outside_mat

        return WindowPanel(
            False,     # closed
            [0, 0],     # x index
            [-x1],
            [y1, y_inside],
            [inside_mat, inside_mat],     # material index
            side_cap_front=1,     # cap index
            side_cap_back=0
            )

    @property
    def frame(self):
        # profil cadre
        #     ___     y0
        #  __|   |
        # |      |    y2
        # |______|    y1
        # x1 x2  x0
        y2 = -0.5 * self.y
        y0 = 0.5 * self.y - self.offset
        y1 = y2 - self.out_frame_y2 - self.overflow_out
        x0 = 0.001   # -min(self.frame_x - 0.001, self.out_frame_offset)
        x1 = x0 - self.out_frame_x
        x2 = x0 - self.out_frame_y
        # y = depth
        # x = width
        if self.out_frame_x <= self.out_frame_y:
            if self.out_frame_x == 0:
                pts_y = [y2, y0, y0, y2]
            else:
                pts_y = [y1, y0, y0, y1]
            return WindowPanel(
                True,     # closed profil
                [0, 0, 1, 1],     # x index
                [x2, x0],
                pts_y,
                [0, 0, 0, 0],     # material index
                closed_path=bool(self.shape == 'CIRCLE')  # closed path
                )
        else:
            return WindowPanel(
                True,     # closed profil
                [0, 0, 1, 1, 2, 2],     # x index
                [x1, x2, x0],
                [y1, y2, y2, y0, y0, y1],
                [0, 0, 0, 0, 0, 0],     # material index
                closed_path=bool(self.shape == 'CIRCLE')   # closed path
                )

    @property
    def out_tablet(self):
        # profil tablette
        #  __  y0
        # |  | y2
        # | / y3
        # |_| y1
        # x0 x2 x1
        y0 = 0.5 * self.y - self.offset
        y1 = -0.5 * self.y - self.out_tablet_y - self.overflow_out
        y2 = y0 - 0.01
        y3 = y2 - 0.04
        x2 = 0.001
        x0 = x2 - self.out_tablet_z
        x1 = x2 + 0.3 * self.frame_x
        # y = depth
        # x = width1
        return WindowPanel(
            True,     # closed profil
            [1, 1, 2, 2, 0, 0],     # x index
            [x0, x2, x1],
            [y1, y3, y2, y0, y0, y1],
            [4, 3, 3, 4, 4, 4],     # material index
            closed_path=False           # closed path
            )

    @property
    def in_tablet(self):
        # profil tablette
        #  __  y0
        # |  |
        # |  |
        # |__| y1
        # x0  x1
        y0 = 0.5 * self.y - self.offset + self.frame_y
        y1 = 0.5 * self.y + self.in_tablet_y
        if self.window_type == 'RAIL':
            y0 = 0.5 * self.y - self.offset + max(self.frame_y, 2.05 * self.panel_y)
        x0 = -self.frame_x
        x1 = min(x0 + self.in_tablet_z, x0 + self.frame_x - 0.001)
        # y = depth
        # x = width1
        return WindowPanel(
            True,     # closed profil
            [0, 0, 1, 1],     # x index
            [x0, x1],
            [y1, y0, y0, y1],
            [1, 1, 1, 1],     # material index
            closed_path=False           # closed path
            )

    @property
    def vertical_space(self):
        """
            avaliable space for hinges
        """
        center, origin, size, radius = self.get_radius(self.x, self.z)
        offset = Vector((0, self.altitude - self._overflow, 0))
        left, right = self.window.avaliable_vertical_space(self.curve_steps, offset, center, origin,
            size, radius, self.angle_y, 0, shape_z=None, path_type=self.shape)
        return left, right

    @property
    def verts(self):
        center, origin, size, radius = self.get_radius(self._x, self._z)
        is_not_circle = self.shape != 'CIRCLE'
        offset = Vector((0, self.altitude - self._overflow, 0))
        verts = self.window.vertices(self.curve_steps, offset, center, origin,
            size, radius, self.angle_y, 0, shape_z=None, path_type=self.shape)

        if self.out_frame:
            _size = Vector((self.x, self.z, 0))
            _offset = Vector((0, self.altitude, 0))
            _center = Vector((center.x, center.y - self._overflow, center.z))

            if self.shape == 'ELLIPSIS':
                _radius = Vector((
                    radius.x - self._overflow,
                    radius.y - self._overflow,
                    0))

            elif self.shape == 'QUADRI':
                _radius = Vector((self.x, 0, 0))

                if self.angle_y < 0:
                    _center.x = 0.5 * self.x
                else:
                    _center.x = -0.5 * self.x
                fx_z = self.z / self.x
                _center.y = min(
                    self.x / (self.x - self.frame_x) * self.z - self.frame_x * (1 + sqrt(1 + fx_z * fx_z)),
                    abs(tan(self.angle_y) * (self.x))
                    )
            else:
                _radius = Vector((radius.x - self._overflow, 0, 0))

            verts += self.frame.vertices(self.curve_steps, _offset, _center, origin,
                _size, _radius,
                self.angle_y, 0, shape_z=None, path_type=self.shape)

        if is_not_circle and self.out_tablet_enable:
            _offset = Vector((0, self.altitude, 0))
            _size = Vector((
                self.x + 2 * (self.out_tablet_x),
                size.y,
                size.z))
            verts += self.out_tablet.vertices(self.curve_steps, _offset, center, origin,
                _size, radius, self.angle_y, 0, shape_z=None, path_type='HORIZONTAL')

        if is_not_circle and self.in_tablet_enable:
            verts += self.in_tablet.vertices(self.curve_steps, offset, center, origin,
                Vector((size.x + 2 * (self.frame_x + self.in_tablet_x), size.y, size.z)),
                radius, self.angle_y, 0, shape_z=None, path_type='HORIZONTAL')

        return verts

    @property
    def faces(self):
        window = self.window
        faces = window.faces(self.curve_steps, path_type=self.shape)
        verts_offset = window.n_verts(self.curve_steps, path_type=self.shape)
        is_not_circle = self.shape != 'CIRCLE'
        if self.out_frame:
            frame = self.frame
            faces += frame.faces(self.curve_steps, path_type=self.shape, offset=verts_offset)
            verts_offset += frame.n_verts(self.curve_steps, path_type=self.shape)
        if is_not_circle and self.out_tablet_enable:
            tablet = self.out_tablet
            faces += tablet.faces(self.curve_steps, path_type='HORIZONTAL', offset=verts_offset)
            verts_offset += tablet.n_verts(self.curve_steps, path_type='HORIZONTAL')
        if is_not_circle and self.in_tablet_enable:
            tablet = self.in_tablet
            faces += tablet.faces(self.curve_steps, path_type='HORIZONTAL', offset=verts_offset)
            verts_offset += tablet.n_verts(self.curve_steps, path_type='HORIZONTAL')

        return faces

    @property
    def matids(self):
        mat = self.window.mat(self.curve_steps, 2, 2, path_type=self.shape)
        is_not_circle = self.shape != 'CIRCLE'
        if self.out_frame:
            mat += self.frame.mat(self.curve_steps, 0, 0, path_type=self.shape)
        if is_not_circle and self.out_tablet_enable:
            mat += self.out_tablet.mat(self.curve_steps, 0, 0, path_type='HORIZONTAL')
        if is_not_circle and self.in_tablet_enable:
            mat += self.in_tablet.mat(self.curve_steps, 0, 0, path_type='HORIZONTAL')
        return mat

    @property
    def uvs(self):
        center, origin, size, radius = self.get_radius(self._x, self._z)
        uvs = self.window.uv(self.curve_steps, center, origin, size, radius,
            self.angle_y, 0, 0, self.frame_x, path_type=self.shape)
        is_not_circle = self.shape != 'CIRCLE'
        if self.out_frame:
            uvs += self.frame.uv(self.curve_steps, center, origin, size, radius,
                self.angle_y, 0, 0, self.frame_x, path_type=self.shape)
        if is_not_circle and self.out_tablet_enable:
            uvs += self.out_tablet.uv(self.curve_steps, center, origin, size, radius,
                self.angle_y, 0, 0, self.frame_x, path_type='HORIZONTAL')
        if is_not_circle and self.in_tablet_enable:
            uvs += self.in_tablet.uv(self.curve_steps, center, origin, size, radius,
                self.angle_y, 0, 0, self.frame_x, path_type='HORIZONTAL')
        return uvs

    def find_blind(self, o, inside):
        for child in o.children:
            if child.type == 'MESH' and 'archipack_blind' in child.data:
                loc = o.matrix_world.inverted() * child.matrix_world.translation
                if inside:
                    if loc.y > 0:
                        return child
                elif loc.y < 0:
                    return child
        return None

    def update_blind(self, context, o, inside):

        blind = self.find_blind(o, inside)
        if inside:
            enabled = self.blind_inside
            overflow = 2 * self.frame_overflow
            style = 'VENITIAN'
            # half width + handle
            y = 0.02 + 0.08
        else:
            enabled = self.blind_outside
            overflow = 0
            style = 'SLAT'
            y = -0.5 * (self.y - self.offset)

        if enabled:

            x = self.x + overflow
            z = self.z + overflow
            a = self.altitude - 0.5 * overflow

            if blind is None:
                bpy.ops.archipack.blind(
                    x=x,
                    z=z,
                    offset_y=y,
                    altitude=a,
                    frame_enable=inside,
                    frame_depth=2 * y,
                    frame_height=0.04,
                    style=style,
                    randomize=True,
                    auto_manipulate=False
                    )
                blind = context.active_object
                blind.parent = o
                blind.select = False
            else:
                blind.select = True
                context.scene.objects.active = blind
                d = blind.data.archipack_blind[0]
                if (d.x != x or
                        d.z != z or
                        d.offset_y != y or
                        d.altitude != a):
                    d.auto_update = False
                    d.x = x
                    d.z = z
                    d.offset_y = y
                    d.altitude = a
                    d.auto_update = True
                blind.select = False

            if inside:
                tM = Matrix([
                    [-1, 0, 0, 0],
                    [0, -1, 0, 0.5 * self.y],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]
                ])
            else:
                tM = Matrix([
                    [1, 0, 0, 0],
                    [0, 1, 0, -0.5 * self.y],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]
                ])
            blind.matrix_world = o.matrix_world * tM

        else:
            self.delete_object(context, blind)

        context.scene.objects.active = o

    def find_portal(self, o):
        for child in o.children:
            if child.type == 'LAMP':
                return child
        return None

    def update_portal(self, context, o):

        lamp = self.find_portal(o)
        if self.portal:
            if lamp is None:
                bpy.ops.object.lamp_add(type='AREA')
                lamp = context.active_object
                lamp.name = "Portal"
                lamp.parent = o

            d = lamp.data
            d.cycles.is_portal = True
            d.shape = 'RECTANGLE'
            d.size = self.x
            d.size_y = self.z

            tM = Matrix([
                [1, 0, 0, 0],
                [0, 0, -1, -0.5 * self.y],
                [0, 1, 0, 0.5 * self.z + self.altitude],
                [0, 0, 0, 1]
            ])
            lamp.matrix_world = o.matrix_world * tM

        else:
            self.delete_object(context, lamp)

        context.scene.objects.active = o
    
    def get_generator(self, o=None):
        return OpeningGenerator(self, o)
        
    def setup_manipulators(self):
        if len(self.manipulators) == 4:
            return
        s = self.manipulators.add()
        s.prop1_name = "x"
        s.prop2_name = "x"
        s.type_key = "SNAP_SIZE_LOC"
        s = self.manipulators.add()
        s.prop1_name = "y"
        s.prop2_name = "y"
        s.type_key = "SNAP_SIZE_LOC"
        s = self.manipulators.add()
        s.prop1_name = "z"
        s.normal = Vector((0, 1, 0))
        s = self.manipulators.add()
        s.prop1_name = "altitude"
        s.normal = Vector((0, 1, 0))

    def remove_childs(self, context, o, to_remove):
        for child in o.children:
            if to_remove < 1:
                return
            if archipack_window_panel.filter(child):
                to_remove -= 1
                self.delete_object(context, child)

    def remove_shutters(self, context, childs, to_remove):
        for child in childs:
            if to_remove < 1:
                return
            to_remove -= 1
            self.delete_object(context, child)

    def update_rows(self, context, o):
        # remove rows
        for i in range(len(self.rows), self.n_rows, -1):
            self.rows.remove(i - 1)

        # add rows
        for i in range(len(self.rows), self.n_rows):
            self.rows.add()

        # wanted childs
        if self.shape == 'CIRCLE':
            w_childs = 1
        elif self.window_type == 'RAIL':
            w_childs = self.rows[0].cols
        else:
            w_childs = sum([row.cols for row in self.rows])

        # real childs
        childs = self.get_childs_panels(context, o)
        n_childs = len(childs)

        # remove child
        if n_childs > w_childs:
            self.remove_childs(context, o, n_childs - w_childs)

    def get_childs_panels(self, context, o):
        return [child for child in o.children if archipack_window_panel.filter(child)]

    def get_childs_shutters(self, context, o, left_side):
        if left_side:
            return [child for child in o.children
                if (archipack_window_shutter.filter(child) and
                    child.data.archipack_window_shutter[0].pivot > 0)
                ]
        else:
            return [child for child in o.children
                if (archipack_window_shutter.filter(child) and
                    child.data.archipack_window_shutter[0].pivot < 0)
                ]

    def adjust_size_and_origin(self, size, origin, pivot, materials):
        if len(size) > 1:
            size[0].x += 0.5 * self.frame_x
            size[-1].x += 0.5 * self.frame_x
        for i in range(1, len(size) - 1):
            size[i].x += 0.5 * self.frame_x
            origin[i].x += -0.25 * self.frame_x * pivot[i]
        for i, o in enumerate(origin):

            o.y = -(i % 2) * self.panel_y - min(0, 0.5 * self.frame_y - 1.025 * self.panel_y)
        for i, o in enumerate(origin):
            materials[i] = (1 - (i % 2)) + 1

    def find_handle(self, o):
        for handle in o.children:
            if 'archipack_handle' in handle:
                return handle
        return None

    def _synch_childs(self, context, o, linked, childs):
        """
            sub synch childs nodes of linked object
        """

        # remove childs not found on source
        l_childs = self.get_childs_panels(context, linked)
        c_names = [c.data.name for c in childs]
        for c in l_childs:
            try:
                id = c_names.index(c.data.name)
            except:
                self.delete_object(context, c)

        # children ordering may not be the same, so get the right l_childs order
        l_childs = self.get_childs_panels(context, linked)
        l_names = [c.data.name for c in l_childs]
        order = []
        for c in childs:
            try:
                id = l_names.index(c.data.name)
            except:
                id = -1
            order.append(id)

        # add missing childs and update other ones
        for i, child in enumerate(childs):
            if order[i] < 0:
                p = bpy.data.objects.new("Window Panel", child.data)
                # Link object into scene
                self.link_object_to_scene(context, p)
                p.show_transparent = True
                p.lock_location[1] = True
                p.lock_location[2] = True
                p.lock_rotation[1] = True
                p.lock_scale[0] = True
                p.lock_scale[1] = True
                p.lock_scale[2] = True
                p.parent = linked
                p.matrix_world = linked.matrix_world.copy()
                m = p.archipack_material.add()
                m.category = 'window'
                m.material = o.archipack_material[0].material
            else:
                p = l_childs[order[i]]

            self.synch_locks(p)

            # update handle
            handle = self.find_handle(child)
            h = self.find_handle(p)
            if handle is not None:
                if h is None:
                    h = create_handle(context, p, handle.data)
                h.location = handle.location.copy()
            elif h is not None:
                self.delete_object(context, h)

            p.location = child.location.copy()

        # restore context
        self.select_object(context, o, True)

    def _synch_shutters(self, context, o, linked, left_side):
        """
            sub synch childs nodes of linked object
        """
        childs = self.get_childs_shutters(context, o, left_side)
        # remove childs not found on source
        l_childs = self.get_childs_shutters(context, linked, left_side)
        c_names = [c.data.name for c in childs]
        for c in l_childs:
            try:
                id = c_names.index(c.data.name)
            except:
                self.delete_object(context, c)

        # children ordering may not be the same, so get the right l_childs order
        l_childs = self.get_childs_shutters(context, linked, left_side)
        l_names = [c.data.name for c in l_childs]
        order = []
        for c in childs:
            try:
                id = l_names.index(c.data.name)
            except:
                id = -1
            order.append(id)

        # add missing childs and update other ones
        for i, child in enumerate(childs):
            if order[i] < 0:
                p = bpy.data.objects.new("Shutter", child.data)
                # Link object into scene
                self.link_object_to_scene(context, p)
                p.lock_location[1] = True
                p.lock_location[2] = True
                p.lock_rotation[1] = True
                p.lock_scale[0] = True
                p.lock_scale[1] = True
                p.lock_scale[2] = True
                p.parent = linked
                p.matrix_world = linked.matrix_world.copy()
                m = p.archipack_material.add()
                m.category = 'window'
                m.material = o.archipack_material[0].material
            else:
                p = l_childs[order[i]]

            # self.synch_locks(p)

            p.location = child.location.copy()
            p.rotation_euler = child.rotation_euler.copy()

        # select and make active
        self.select_object(context, o, True)
        
    def _synch_hole(self, context, linked, hole):
        l_hole = self.find_hole(linked)
        if l_hole is None:
            l_hole = bpy.data.objects.new("hole", hole.data)
            l_hole['archipack_hole'] = True
            # Link object into scene
            self.link_object_to_scene(context, l_hole)
            l_hole.parent = linked
            l_hole.matrix_world = linked.matrix_world.copy()
            l_hole.location = hole.location.copy()
        else:
            l_hole.data = hole.data

    def synch_locks(self, p):
        p.lock_location[0] = self.window_type != 'RAIL'
        p.lock_rotation[0] = self.window_type == 'RAIL'
        p.lock_rotation[2] = self.window_type == 'RAIL'

    def synch_childs(self, context, o):
        """
            synch childs nodes of linked objects
        """
        bpy.ops.object.select_all(action='DESELECT')
        # select and make active
        self.select_object(context, o, True)
        childs = self.get_childs_panels(context, o)
        hole = self.find_hole(o)
        bpy.ops.object.select_linked(type='OBDATA')
        for linked in context.selected_objects:
            if linked != o:
                ld = archipack_window.datablock(linked)
                ld.update_portal(context, linked)
                ld.update_blind(context, linked, True)
                ld.update_blind(context, linked, False)
                self._synch_childs(context, o, linked, childs)
                self._synch_shutters(context, o, linked, True)
                self._synch_shutters(context, o, linked, False)
                if hole is not None:
                    self._synch_hole(context, linked, hole)

    def get_shutter_row(self, x, y, left_side):
        n_shutters = self.shutter_left + self.shutter_right
        size = Vector((x / n_shutters, y, 0))
        origin = []
        ttl = 0
        xh = x / 2
        # offset pivot
        if left_side:
            n_shutters = self.shutter_left
            ttl -= size.x
        else:
            ttl += self.shutter_left * size.x
            n_shutters = self.shutter_right

        for i in range(n_shutters):
            ttl += size.x
            origin.append(Vector((ttl - xh, 0)))
        return size, origin

    def update_shutter(self, context, o, left_side, hinge_space):
        # wanted childs
        if self.shutter_enable:
            if left_side:
                side = 0
                pivot = 1
                n_shutters = self.shutter_left
            else:
                pivot = -1
                n_shutters = self.shutter_right
                side = n_shutters - 1

        else:
            n_shutters = 0

        # real childs
        childs = self.get_childs_shutters(context, o, left_side)
        n_childs = len(childs)

        # remove child
        if n_childs > n_shutters:
            self.remove_shutters(context, childs, n_childs - n_shutters)

        if not self.shutter_enable or n_shutters == 0:
            return

        childs = self.get_childs_shutters(context, o, left_side)
        n_childs = len(childs)

        location_y = -0.5 * self.y - 0.25 * self.shutter_depth - self.overflow_out
        if self.out_frame:
            location_y -= self.out_frame_y2
        
        # Note: radius is slightly wrong: not taking overflow in account
        center, origin, size, radius = self.get_radius(self.x, self.z)
        offset = Vector((0.05, 0))
        size, origin = self.get_shutter_row(self.x, self.z, left_side)

        if hinge_space > 1.5:
            hinge_count = 3
        else:
            hinge_count = 2

        for panel in range(n_shutters):

            if panel >= n_childs:
                bpy.ops.archipack.window_shutter(
                    center=center,
                    origin=Vector((origin[panel].x, offset.y, 0)),
                    size=size,
                    radius=radius,
                    pivot=pivot,
                    shape=self.shape,
                    offset=pivot * offset.x,
                    curve_steps=self.curve_steps,
                    border=self.shutter_border,
                    depth=self.shutter_depth,
                    angle_y=self.angle_y,
                    hinge_enable=panel == side,
                    hinge_count=hinge_count,
                    hinge_space=hinge_space,
                    hinge_size=self.shutter_hinge,
                    material=o.archipack_material[0].material
                )
                child = context.active_object
                # parenting at 0, 0, 0 before set object matrix_world
                # so location remains local from frame
                child.parent = o
                child.matrix_world = o.matrix_world.copy()
                child.rotation_euler.z = pi
            else:
                child = childs[panel]
                # select and make active
                self.select_object(context, child, True)
                props = archipack_window_shutter.datablock(child)
                if props is not None:
                    props.origin = Vector((origin[panel].x, offset.y, 0))
                    props.center = center
                    props.radius = radius
                    props.size = size
                    props.pivot = pivot
                    props.shape = self.shape
                    props.offset = pivot * offset.x
                    props.curve_steps = self.curve_steps
                    props.border = self.shutter_border
                    props.depth = self.shutter_depth
                    props.angle_y = self.angle_y
                    props.hinge_enable = panel == side
                    props.hinge_count = hinge_count
                    props.hinge_space = hinge_space
                    props.hinge_size = self.shutter_hinge
                    props.update(context)

            # location y + frame width.
            child.location = Vector((
                origin[panel].x - pivot * offset.x + (side - panel) * size.x,
                origin[panel].y + location_y + (side - panel) * pivot * 0.5 * self.shutter_depth,
                self.altitude + offset.y
                ))

    def update_childs(self, context, o):
        """
            pass params to childrens
        """
        self.update_rows(context, o)
        childs = self.get_childs_panels(context, o)
        n_childs = len(childs)
        child_n = 0
        row_n = 0
        location_y = 0.5 * self.y - self.offset + 0.5 * self.frame_y
        center, origin, size, radius = self.get_radius(self._x, self._z)
        offset = Vector((0, 0))
        handle = 'NONE'
        if self.shape != 'CIRCLE':
            if self.handle_enable:
                if self._z > 1.8:
                    handle = 'BOTH'
                else:
                    handle = 'INSIDE'
            is_circle = False
        else:
            is_circle = True

        if self.window_type == 'RAIL':
            handle_model = 2
        else:
            handle_model = 1

        for row in self.rows:
            row_n += 1
            if row_n < self.n_rows and not is_circle and self.window_type != 'RAIL':
                z = row.height
                shape = 'RECTANGLE'
            else:
                z = max(2 * self.frame_x + 0.001, self._z - offset.y)
                shape = self.shape

            self.warning = bool(z > self._z - offset.y)
            if self.warning:
                break

            size, origin, pivot = row.get_row(self._x, z)

            # side materials
            materials = [0 for i in range(row.cols)]

            handle_altitude = min(
                max(4 * self.frame_x, self.handle_altitude + self._overflow - self.altitude),
                z - 4 * self.frame_x
                )

            if self.window_type == 'RAIL':
                self.adjust_size_and_origin(size, origin, pivot, materials)

            for panel in range(row.cols):
                child_n += 1

                if row.fixed[panel]:
                    enable_handle = 'NONE'
                else:
                    enable_handle = handle

                if child_n > n_childs:
                    bpy.ops.archipack.window_panel(
                        center=center,
                        origin=Vector((origin[panel].x, offset.y, 0)),
                        size=size[panel],
                        radius=radius,
                        pivot=pivot[panel],
                        shape=shape,
                        fixed=row.fixed[panel],
                        handle=enable_handle,
                        handle_model=handle_model,
                        handle_altitude=handle_altitude,
                        curve_steps=self.curve_steps,
                        side_material=materials[panel],
                        frame_x=self.panel_x,
                        frame_y=self.panel_y,
                        angle_y=self.angle_y,
                        enable_glass=self.enable_glass,
                        material=o.archipack_material[0].material
                    )
                    child = context.active_object
                    # parenting at 0, 0, 0 before set object matrix_world
                    # so location remains local from frame
                    child.parent = o
                    child.matrix_world = o.matrix_world.copy()

                else:
                    child = childs[child_n - 1]
                    # select and make active
                    self.select_object(context, child, True)
                    props = archipack_window_panel.datablock(child)
                    if props is not None:
                        props.origin = Vector((origin[panel].x, offset.y, 0))
                        props.center = center
                        props.radius = radius
                        props.size = size[panel]
                        props.pivot = pivot[panel]
                        props.shape = shape
                        props.fixed = row.fixed[panel]
                        props.handle = enable_handle
                        props.handle_model = handle_model
                        props.handle_altitude = handle_altitude
                        props.side_material = materials[panel]
                        props.curve_steps = self.curve_steps
                        props.frame_x = self.panel_x
                        props.frame_y = self.panel_y
                        props.angle_y = self.angle_y
                        props.enable_glass = self.enable_glass
                        props.update(context)
                # location y + frame width. frame depends on choosen profile (fixed or not)
                # update linked childs location too
                child.location = Vector((
                    origin[panel].x,
                    origin[panel].y + location_y + self.panel_y,
                    self.altitude - self._overflow + offset.y))

                self.synch_locks(child)

                if not row.fixed[panel]:
                    handle = 'NONE'

                # only one single panel allowed for circle
                if is_circle:
                    return

            # only one single row allowed for rail window
            if self.window_type == 'RAIL':
                return
            offset.y += row.height

    @property
    def _overflow(self):
        return min(0, self.frame_overflow - self.frame_x)

    @property
    def _x(self):
        return self.x + 2 * self._overflow

    @property
    def _z(self):
        return self.z + 2 * self._overflow

    def _get_tri_radius(self, _x, _z):
        return Vector((0, self.y, 0)), Vector((0, 0, 0)), \
            Vector((_x, _z, 0)), Vector((_x, 0, 0))

    def _get_quad_radius(self, _x, _z):
        fx_z = _z / _x
        center_y = min(_x / (_x - self.frame_x) * _z - self.frame_x * (1 + sqrt(1 + fx_z * fx_z)),
            abs(tan(self.angle_y) * _x))
        if self.angle_y < 0:
            center_x = 0.5 * _x
        else:
            center_x = -0.5 * _x
        return Vector((center_x, center_y, 0)), Vector((0, 0, 0)), \
            Vector((_x, _z, 0)), Vector((_x, 0, 0))

    def _get_round_radius(self, _x, _z):
        """
            bound radius to available space
            return center, origin, size, radius
        """
        x = 0.5 * _x - self.frame_x
        # minimum space available
        y = _z - sum([row.height for row in self.rows[:self.n_rows - 1]]) - 2 * self.frame_x
        y = min(y, x)
        # minimum radius inside
        r = y + x * (x - (y * y / x)) / (2 * y)
        radius = max(self.radius, 0.001 + self.frame_x + r)
        return Vector((0, _z - radius, 0)), Vector((0, 0, 0)), \
            Vector((_x, _z, 0)), Vector((radius, 0, 0))

    def _get_circle_radius(self, _x, _z):
        """
            return center, origin, size, radius
        """
        return Vector((0, 0.5 * _x, 0)), Vector((0, 0, 0)), \
            Vector((_x, _z, 0)), Vector((0.5 * _x, 0, 0))

    def _get_ellipsis_radius(self, _x, _z):
        """
            return center, origin, size, radius
        """
        y = self.z - sum([row.height for row in self.rows[:self.n_rows - 1]])
        radius_b = max(0, 0.001 - 2 * self.frame_x + min(y, self.elipsis_b))
        return Vector((0, _z - radius_b, 0)), Vector((0, 0, 0)), \
            Vector((_x, _z, 0)), Vector((_x / 2, radius_b, 0))

    def get_radius(self, _x, _z):
        """
            return center, origin, size, radius
        """
        if self.shape == 'ROUND':
            return self._get_round_radius(_x, _z)
        elif self.shape == 'ELLIPSIS':
            return self._get_ellipsis_radius(_x, _z)
        elif self.shape == 'CIRCLE':
            return self._get_circle_radius(_x, _z)
        elif self.shape == 'QUADRI':
            return self._get_quad_radius(_x, _z)
        elif self.shape in ['TRIANGLE', 'PENTAGON']:
            return self._get_tri_radius(_x, _z)
        else:
            return Vector((0, 0, 0)), Vector((0, 0, 0)), \
                Vector((_x, _z, 0)), Vector((0, 0, 0))

    def update(self, context, childs_only=False):
        # support for "copy to selected"
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        self.setup_manipulators()

        if childs_only is False:
            bmed.buildmesh(context, o, self.verts, self.faces, self.matids, self.uvs)

        self.update_portal(context, o)
        self.update_blind(context, o, True)
        self.update_blind(context, o, False)
        self.update_childs(context, o)

        left, right = self.vertical_space
        self.update_shutter(context, o, True, left)
        self.update_shutter(context, o, False, right)

        # update hole
        if childs_only is False and self.find_hole(o) is not None:
            self.interactive_hole(context, o)

        # store 3d points for gl manipulators
        x, y = 0.5 * self.x, 0.5 * self.y
        self.manipulators[0].set_pts([(-x, -y, 0), (x, -y, 0), (0.5, 0, 0)])
        self.manipulators[1].set_pts([(-x, -y, 0), (-x, y, 0), (-1, 0, 0)])
        self.manipulators[2].set_pts([(x, -y, self.altitude), (x, -y, self.altitude + self.z), (-1, 0, 0)])
        self.manipulators[3].set_pts([(x, -y, 0), (x, -y, self.altitude), (-1, 0, 0)])

        # self.update_dimension(context, o)
        self.add_dimension_point(0, Vector((-x, -y, 0)))
        self.add_dimension_point(1, Vector((x, -y, 0)))
        
        self.add_dimension_point(2, Vector((-0.5 * self._x - self.frame_x, y, 0)))
        self.add_dimension_point(3, Vector((0.5 * self._x + self.frame_x, y, 0)))
        
        # support for instances childs, update at object level
        self.synch_childs(context, o)

        # synch dimensions when apply
        if o.parent:
            self.update_dimensions(context, o)
            
        # restore context
        self.restore_context(context)

    def find_hole(self, o):
        for child in o.children:
            if 'archipack_hole' in child:
                return child
        return None

    def interactive_hole(self, context, o):
        hole_obj = self.find_hole(o)

        if hole_obj is None:
            m = bpy.data.meshes.new("hole")
            hole_obj = bpy.data.objects.new("hole", m)
            hole_obj['archipack_hole'] = True
            # Link object into scene
            self.link_object_to_scene(context, hole_obj)
            hole_obj.parent = o
            hole_obj.matrix_world = o.matrix_world.copy()
        else:
            # must ensure the hole layer is visible
            self.add_to_layer(context, hole_obj)
            
        hole = self.hole
        center, origin, size, radius = self.get_radius(self._x, self._z)
        x0 = 0

        if self.out_frame:
            x0 += min(self.frame_x + 0.001, self.out_frame_y + self.out_frame_offset)

        if self.out_tablet_enable:
            x0 -= self.out_tablet_z

        x0 = min(x0, -0.001)

        shape_z = [-0.001, x0]

        verts = hole.vertices(self.curve_steps,
            Vector((0, self.altitude - self._overflow, 0)),
            center, origin, size, radius,
            self.angle_y, 0, shape_z=shape_z, path_type=self.shape)

        faces = hole.faces(self.curve_steps, path_type=self.shape)

        matids = hole.mat(self.curve_steps, 2, 2, path_type=self.shape)

        uvs = hole.uv(self.curve_steps, center, origin, size, radius,
            self.angle_y, 0, 0, self.frame_x, path_type=self.shape)

        bmed.buildmesh(context, hole_obj, verts, faces, matids=matids, uvs=uvs, auto_smooth=False)
        return hole_obj

    def _add_spline(self, curve, coords):
        spline = curve.splines.new('POLY')
        spline.use_endpoint_u = False
        spline.use_cyclic_u = coords[-1] == coords[0]
        if coords[-1] == coords[0]:
            coords.pop()
        spline.points.add(len(coords) - 1)
        for i, coord in enumerate(coords):
            x, y, z = coord
            spline.points[i].co = (x, y, z, 1)

    def _to_curve(self, context, coords, name: str, dimensions: str='3D'):
        curve = bpy.data.curves.new(name, type='CURVE')
        curve.dimensions = dimensions
        for co in coords:
            self._add_spline(curve, co)
        curve_obj = bpy.data.objects.new(name, curve)
        # Link object into scene
        self.link_object_to_scene(context, curve_obj)
        # select and make active
        self.select_object(context, curve_obj, True)
        return curve_obj

    def as_2d(self, context, o):

        # frame
        center, origin, size, radius = self.get_radius(self._x, self._z)
        offset = Vector((0, 0, 0))

        # draw borders
        connect = 2
        if self.window_type == 'RAIL':
            connect = 3

        coords = self.window.as_2d(self.curve_steps, offset, center, origin,
            size, radius, 0, 0, connect=connect)

        # panels
        childs = self.get_childs_panels(context, o)
        
        row_n = 0
        location_y = 0.5 * self.y - self.offset + 0.5 * self.frame_y

        row = self.rows[0]
        row_n += 1

        materials = [0 for i in range(row.cols)]
        
        n_childs = len(childs)
        
        for panel in range(row.cols):
            if panel >= n_childs:
                break
            child = childs[panel]
            
            # location y + frame width. frame depends on choosen profile (fixed or not)
            # update linked childs location too
            d = archipack_window_panel.datablock(child) 
            if d.shape == 'CIRCLE':
                s = 2
            else: 
                s = 1
            location = Vector((
                d.origin.x,
                d.origin.y + location_y + self.panel_y,
                0))

            coords.extend(
                d.window.as_2d(
                    self.curve_steps, location, d.center, d.origin,
                    s * d.size, d.radius, 0, d.pivot, path_type=d.shape)
                )
            # arc
            if self.window_type == 'FLAT' and not d.fixed:
                x, y = location.x, location.y
                r = s * d.size.x
                steps = 8
                # 30 deg
                da = pi / (6 * steps)
                arc = [Vector((
                        x + d.pivot * r * cos(da * i),
                        y + r * sin(da * i),
                        0))
                    for i in range(steps + 1)]
                arc.append(location)
                coords.append(arc)

        curve = self._to_curve(context, coords, name="{}-2d".format(o.name), dimensions='2D')
        curve.matrix_world = o.matrix_world.copy()

        return curve

    def hole_2d(self, mode='SYMBOL'):
        """
          return coords of full / inside hole in 2d top ortho view
        """
        if mode == 'BOUND':
            x, y = 0.5 * self._x + self.frame_x, 0.5 * self.y + self.hole_margin
            return [(-x, -y, 0), (-x, y, 0), (x, y, 0), (x, -y, 0), (-x, -y, 0)]

        center, origin, size, radius = self.get_radius(self._x, self._z)
        if mode == 'FLOORS':
            coords = self.inside_hole.as_2d(self.curve_steps,
                Vector((0, 0, 0)),
                center, origin, size, radius,
                0, 0, shape_z=None, path_type=self.shape)
        else:
            coords = self.hole.as_2d(self.curve_steps,
                Vector((0, 0, 0)),
                center, origin, size, radius,
                0, 0, shape_z=None, path_type=self.shape)
        # Use only first curve
        return coords[0]


class ARCHIPACK_PT_window(Panel):
    bl_idname = "ARCHIPACK_PT_window"
    bl_label = "Window"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    # bl_context = 'object'
    bl_category = 'ArchiPack'

    # layout related
    display_detail = BoolProperty(
        default=False
    )
    display_panels = BoolProperty(
        default=True
    )

    @classmethod
    def poll(cls, context):
        return archipack_window.filter(context.active_object)

    def draw(self, context):
        o = context.active_object
        prop = archipack_window.datablock(o)
        if prop is None:
            return
        layout = self.layout
        layout.operator('archipack.manipulate', icon='HAND')
        row = layout.row(align=True)
        row.operator('archipack.window', text="Refresh", icon='FILE_REFRESH').mode = 'REFRESH'
        if o.data.users > 1:
            row.operator('archipack.window', text="Make unique ({})".format(o.data.users), icon='UNLINKED').mode = 'UNIQUE'
        row.operator('archipack.window', text="Delete", icon='ERROR').mode = 'DELETE'
        box = layout.box()
        # box.label(text="Styles")
        row = box.row(align=True)
        row.operator("archipack.window_preset_menu", text=bpy.types.ARCHIPACK_OT_window_preset_menu.bl_label)
        row.operator("archipack.window_preset", text="", icon='ZOOMIN')
        row.operator("archipack.window_preset", text="", icon='ZOOMOUT').remove_active = True
        box = layout.box()
        box.prop(prop, 'window_type')
        box.prop(prop, 'x')
        box.prop(prop, 'y')
        if prop.window_shape != 'CIRCLE':
            box.prop(prop, 'z')
            if prop.warning:
                box.label(text="Insufficient height", icon='ERROR')
        box.prop(prop, 'altitude')
        box.prop(prop, 'offset')
        box.prop(prop, 'overflow_out')
            
        if prop.window_type == 'FLAT':
            box = layout.box()
            box.prop(prop, 'window_shape')
            if prop.window_shape in ['ROUND', 'CIRCLE', 'ELLIPSIS']:
                box.prop(prop, 'curve_steps')
            if prop.window_shape in ['ROUND']:
                box.prop(prop, 'radius')
            elif prop.window_shape == 'ELLIPSIS':
                box.prop(prop, 'elipsis_b')
            elif prop.window_shape == 'QUADRI':
                box.prop(prop, 'angle_y')

        box = layout.box()
        icon = "TRIA_RIGHT"
        if prop.display_detail:
            icon = "TRIA_DOWN"
        
        box.prop(prop, 'display_detail', icon=icon, icon_only=True, text="Components", emboss=True)
        
        if prop.display_detail:
            box.prop(prop, 'enable_glass')
            box = layout.box()
            box.label("Window frame")
            box.prop(prop, 'frame_x')
            box.prop(prop, 'frame_y')
            box.prop(prop, 'frame_overflow')
            box = layout.box()
            box.label("Panel frames")
            box.prop(prop, 'panel_x')
            box.prop(prop, 'panel_y')
            if prop.window_shape != 'CIRCLE':
                box = layout.box()
                row = box.row(align=True)
                row.prop(prop, 'handle_enable')
                if prop.handle_enable:
                    box.prop(prop, 'handle_altitude')
            box = layout.box()
            row = box.row(align=True)
            row.prop(prop, 'out_frame')
            if prop.out_frame:
                box.prop(prop, 'out_frame_x')
                box.prop(prop, 'out_frame_y2')
                box.prop(prop, 'out_frame_y')
                box.prop(prop, 'out_frame_offset')
            if prop.window_shape != 'CIRCLE':
                box = layout.box()
                row = box.row(align=True)
                row.prop(prop, 'out_tablet_enable')
                if prop.out_tablet_enable:
                    box.prop(prop, 'out_tablet_x')
                    box.prop(prop, 'out_tablet_y')
                    box.prop(prop, 'out_tablet_z')
                box = layout.box()
                row = box.row(align=True)
                row.prop(prop, 'in_tablet_enable')
                if prop.in_tablet_enable:
                    box.prop(prop, 'in_tablet_x')
                    box.prop(prop, 'in_tablet_y')
                    box.prop(prop, 'in_tablet_z')
                box = layout.box()
                box.prop(prop, 'blind_inside')
                box.prop(prop, 'blind_outside')
                """
                row.prop(prop, 'blind_enable')
                if prop.blind_enable:
                    box.prop(prop, 'blind_open')
                    box.prop(prop, 'blind_y')
                    box.prop(prop, 'blind_z')
                """
                box = layout.box()
                row = box.row(align=True)
                row.prop(prop, 'shutter_enable')
                if prop.shutter_enable:
                    box.prop(prop, 'shutter_left')
                    box.prop(prop, 'shutter_right')
                    box.prop(prop, 'shutter_depth')
                    box.prop(prop, 'shutter_border')
                    box.prop(prop, 'shutter_hinge')

        if prop.window_shape != 'CIRCLE':
            box = layout.box()
            row = box.row(align=False)
            
            icon = "TRIA_RIGHT"
            if prop.display_panels:
                icon = "TRIA_DOWN"
            row.prop(prop, 'display_panels', icon=icon, icon_only=True, text="Rows", emboss=True)
            
            if prop.window_type != 'RAIL':
                row.prop(prop, 'n_rows')
            
            if prop.display_panels:
                if prop.window_type != 'RAIL':
                    last_row = prop.n_rows - 1
                    for i, row in enumerate(prop.rows):
                        box = layout.box()
                        box.label(text="Row " + str(i + 1))
                        row.draw(box, context, i == last_row)
                else:
                    box = layout.box()
                    row = prop.rows[0]
                    row.draw(box, context, True)
        box = layout.box()
        icon = "TRIA_RIGHT"
        if prop.display_materials:
            icon = "TRIA_DOWN"
        box.prop(prop, 'display_materials', icon=icon, icon_only=True, text="Materials", emboss=True)
        if prop.display_materials:
            box.label("Hole")
            box.prop(prop, 'hole_inside_mat')
            box.prop(prop, 'hole_outside_mat')

        layout.prop(prop, 'portal', icon="LAMP_AREA")


class ARCHIPACK_PT_window_panel(Panel):
    bl_idname = "ARCHIPACK_PT_window_panel"
    bl_label = "Window panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_window_panel.filter(context.active_object)

    def draw(self, context):
        layout = self.layout
        layout.operator("archipack.select_parent")


class ARCHIPACK_PT_window_shutter(Panel):
    bl_idname = "ARCHIPACK_PT_window_shutter"
    bl_label = "Shutter"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_window_shutter.filter(context.active_object)

    def draw(self, context):
        layout = self.layout
        layout.operator("archipack.select_parent")


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_window(ArchipackCreateTool, Operator):
    bl_idname = "archipack.window"
    bl_label = "Window"
    bl_description = "Window"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    x = FloatProperty(
        name='width',
        min=0.1, max=10000,
        default=2.0, precision=2,
        description='Width'
        )
    y = FloatProperty(
        name='depth',
        min=0.1, max=10000,
        default=0.20, precision=2,
        description='Depth'
        )
    z = FloatProperty(
        name='height',
        min=0.1, max=10000,
        default=1.2, precision=2,
        description='height'
        )
    altitude = FloatProperty(
        name='altitude',
        min=0.0, max=10000,
        default=1.0, precision=2,
        description='altitude'
        )
    mode = EnumProperty(
        items=(
        ('CREATE', 'Create', '', 0),
        ('DELETE', 'Delete', '', 1),
        ('REFRESH', 'Refresh', '', 2),
        ('UNIQUE', 'Make unique', '', 3),
        ),
        default='CREATE'
        )
    # auto_manipulate = BoolProperty(default=True)

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        m = bpy.data.meshes.new("Window")
        o = bpy.data.objects.new("Window", m)
        d = m.archipack_window.add()
        d.x = self.x
        d.y = self.y
        d.z = self.z
        d.altitude = self.altitude
        # Link object into scene
        self.link_object_to_scene(context, o)
        
        # select and make active
        self.select_object(context, o, True)
        self.add_material(o)
        self.load_preset(d)
        # select frame
        o.select = True
        context.scene.objects.active = o
        return o

    def delete(self, context):
        o = context.active_object
        if archipack_window.filter(o):
            bpy.ops.archipack.disable_manipulate()
            parent = o.parent
            self.delete_object(context, o)
            # synch wall dimensions when apply
            if parent:
                for c in parent.children:
                    if c.data and "archipack_wall2" in c.data:
                        c.data.archipack_wall2[0].synch_dimension(context, c)

    def update(self, context):
        o = context.active_object
        d = archipack_window.datablock(o)
        if d is not None:
            d.update(context)
            bpy.ops.object.select_linked(type='OBDATA')
            for linked in context.selected_objects:
                if linked != o:
                    archipack_window.datablock(linked).update(context)

        bpy.ops.object.select_all(action="DESELECT")
        # select and make active
        self.select_object(context, o, True)

    def unique(self, context):
        act = context.active_object
        sel = [o for o in context.selected_objects]
        bpy.ops.object.select_all(action="DESELECT")
        for o in sel:
            if archipack_window.filter(o):
                # select and make active
                self.select_object(context, o)
                for child in o.children:
                    d = child.data
                    if 'archipack_hole' in child or (
                            d is not None and (
                            'archipack_window_panel' in d or
                            'archipack_window_shutter' in d or
                            'archipack_dimension' in d
                            )):
                        # Not available in 2.8
                        child.hide_select = False
                        self.select_object(context, child)
                        
                        # select text and handles too
                        for c in child.children:
                            self.select_object(context, c)
                            
        if len(context.selected_objects) > 0:
            bpy.ops.object.make_single_user(type='SELECTED_OBJECTS', object=True,
                obdata=True, material=False, texture=False, animation=False)
            # Not available in 2.8
            for child in context.selected_objects:
                if 'archipack_hole' in child:
                    child.hide_select = True
        
        bpy.ops.object.select_all(action="DESELECT")
        self.select_object(context, act, True)
        for o in sel:
            self.select_object(context, o)
            
    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            if self.mode == 'CREATE':
                bpy.ops.object.select_all(action="DESELECT")
                o = self.create(context)
                o.location = bpy.context.scene.cursor_location
                # select and make active
                self.select_object(context, o, True)
                self.manipulate()
            elif self.mode == 'DELETE':
                self.delete(context)
            elif self.mode == 'REFRESH':
                self.update(context)
            elif self.mode == 'UNIQUE':
                self.unique(context)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_window_draw(ArchipackDrawTool, Operator):
    bl_idname = "archipack.window_draw"
    bl_label = "Draw Windows"
    bl_description = "Draw Windows over walls"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    filepath = StringProperty(default="")
    feedback = None
    stack = []
    object_name = ""

    @classmethod
    def poll(cls, context):
        return True

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def draw_callback(self, _self, context):
        self.feedback.draw(context)

    def add_object(self, context, event):

        o = context.active_object
        bpy.ops.object.select_all(action="DESELECT")

        if archipack_window.filter(o):

            # select and make active
            self.select_object(context, o, True)

            if event.shift:
                bpy.ops.archipack.window(mode="UNIQUE")

            # instance subs
            new_w = o.copy()
            new_w.data = o.data
            self.link_object_to_scene(context, new_w)
            for child in o.children:
                if "archipack_hole" not in child:
                    new_c = child.copy()
                    new_c.data = child.data
                    new_c.parent = new_w
                    self.link_object_to_scene(context, new_c)
                    # dup handle if any
                    for c in child.children:
                        new_h = c.copy()
                        new_h.data = c.data
                        new_h.parent = new_c
                        self.link_object_to_scene(context, new_h)

            o = new_w
            self.select_object(context, o, True)

        else:
            bpy.ops.archipack.window(auto_manipulate=False, filepath=self.filepath)
            o = context.active_object

        self.object_name = o.name

        bpy.ops.archipack.generate_hole('INVOKE_DEFAULT')
        self.select_object(context, o, True)
        
    def modal(self, context, event):

        context.area.tag_redraw()
        o = self.get_scene_object(context, self.object_name)

        if o is None:
            self.restore_walls(context)
            return {'FINISHED'}

        d = archipack_window.datablock(o)
        hole = None
        if d is not None:
            hole = d.find_hole(o)

        # hide hole from raycast
        to_hide = [o]
        to_hide.extend([child for child in o.children if archipack_window_panel.filter(child)])
        if hole is not None:
            to_hide.append(hole)

        for obj in to_hide:
            self.hide_object(hole)

        res, tM, wall, width, y, z_offset = self.mouse_hover_wall(context, event)

        for obj in to_hide:
            self.show_object(hole)
            
        if res and d is not None:
            o.matrix_world = tM.copy()
            if d.y != width:
                d.y = width
                
        if event.value == 'PRESS':

            if event.type in {'C'}:
                bpy.ops.archipack.window(mode='DELETE')
                self.feedback.disable()
                bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
                bpy.ops.archipack.window_preset_menu('INVOKE_DEFAULT', preset_operator="archipack.window_draw")
                self.restore_walls(context)
                return {'FINISHED'}

            if event.type in {'LEFTMOUSE', 'RET', 'NUMPAD_ENTER', 'SPACE'}:
                if wall is not None:

                    print(wall.name)
                    self.select_object(context, wall, True)
                    if bpy.ops.archipack.single_boolean.poll():
                        bpy.ops.archipack.single_boolean()
                    self.unselect_object(wall)
                    # o must be a window here
                    if d is not None:
                        context.scene.objects.active = o
                        self.stack.append(o)
                        self.add_object(context, event)
                        context.active_object.matrix_world = tM
                    return {'RUNNING_MODAL'}
            # prevent selection of other object
            if event.type in {'RIGHTMOUSE'}:
                return {'RUNNING_MODAL'}

        if self.keymap.check(event, self.keymap.undo) or (
                event.type in {'BACK_SPACE'} and event.value == 'RELEASE'
                ):
            if len(self.stack) > 0:
                last = self.stack.pop()
                self.select_object(context, last, True)
                bpy.ops.archipack.window(mode="DELETE")
                self.select_object(context, o, True)
            return {'RUNNING_MODAL'}

        if event.value == 'RELEASE':

            if event.type in {'ESC', 'RIGHTMOUSE'}:
                bpy.ops.archipack.window(mode='DELETE')
                self.feedback.disable()
                bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
                self.restore_walls(context)
                return {'FINISHED'}

        return {'PASS_THROUGH'}

    def invoke(self, context, event):

        if context.mode == "OBJECT":
            o = None
            self.stack = []
            self.keymap = Keymaps(context)
            # exit manipulate_mode if any
            bpy.ops.archipack.disable_manipulate()
            # invoke with shift pressed will use current object as basis for linked copy
            if self.filepath == '' and archipack_window.filter(context.active_object):
                o = context.active_object
            # hmm 2.8 ??
            context.scene.objects.active = None
            
            bpy.ops.object.select_all(action="DESELECT")
            self.select_object(context, o, True)
            self.add_object(context, event)
            self.feedback = FeedbackPanel()
            self.feedback.instructions(context, "Draw a window", "Click & Drag over a wall", [
                ('LEFTCLICK, RET, SPACE, ENTER', 'Create a window'),
                ('BACKSPACE, CTRL+Z', 'undo last'),
                ('C', 'Choose another window'),
                ('SHIFT', 'Make independant copy'),
                ('RIGHTCLICK or ESC', 'exit')
                ])
            self.feedback.enable()
            args = (self, context)

            self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_callback, args, 'WINDOW', 'POST_PIXEL')
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_window_portals(Operator):
    bl_idname = "archipack.window_portals"
    bl_label = "Portals"
    bl_description = "Create portal for each window"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return True

    def invoke(self, context, event):
        for o in context.scene.objects:
            d = archipack_window.datablock(o)
            if d is not None:
                d.update_portal(context)

        return {'FINISHED'}


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_window_panel(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.window_panel"
    bl_label = "Window panel"
    bl_description = "Window panel"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    center = FloatVectorProperty(
            subtype='XYZ'
            )
    origin = FloatVectorProperty(
            subtype='XYZ'
            )
    size = FloatVectorProperty(
            subtype='XYZ'
            )
    radius = FloatVectorProperty(
            subtype='XYZ'
            )
    angle_y = FloatProperty(
            name='angle',
            unit='ROTATION',
            subtype='ANGLE',
            min=-1.5, max=1.5,
            default=0, precision=2,
            description='angle'
            )
    frame_y = FloatProperty(
            name='Depth',
            min=0, max=100,
            default=0.06, precision=2,
            description='frame depth'
            )
    frame_x = FloatProperty(
            name='Width',
            min=0, max=100,
            default=0.06, precision=2,
            description='frame width'
            )
    curve_steps = IntProperty(
            name="curve steps",
            min=1,
            max=128,
            default=16
            )
    shape = EnumProperty(
            name='Shape',
            items=(
                ('RECTANGLE', 'Rectangle', '', 0),
                ('ROUND', 'Top Round', '', 1),
                ('ELLIPSIS', 'Top Elliptic', '', 2),
                ('QUADRI', 'Top oblique', '', 3),
                ('CIRCLE', 'Full circle', '', 4)
                ),
            default='RECTANGLE'
            )
    pivot = FloatProperty(
            name='pivot',
            min=-1, max=1,
            default=-1, precision=2,
            description='pivot'
            )
    side_material = IntProperty(
            name="side material",
            min=0,
            max=2,
            default=0
            )
    handle = EnumProperty(
            name='Handle',
            items=(
                ('NONE', 'No handle', '', 0),
                ('INSIDE', 'Inside', '', 1),
                ('BOTH', 'Inside and outside', '', 2)
                ),
            default='NONE'
            )
    handle_model = IntProperty(
            name="handle model",
            default=1,
            min=1,
            max=2
            )
    handle_altitude = FloatProperty(
            name='handle altitude',
            min=0, max=1000,
            default=0.2, precision=2,
            description='handle altitude'
            )
    fixed = BoolProperty(
            name="Fixed",
            default=False
            )
    material = StringProperty(
            name="material",
            default=""
            )
    enable_glass = BoolProperty(
            name="Enable glass",
            default=True
            )

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        m = bpy.data.meshes.new("Window Panel")
        o = bpy.data.objects.new("Window Panel", m)
        d = m.archipack_window_panel.add()
        d.center = self.center
        d.origin = self.origin
        d.size = self.size
        d.radius = self.radius
        d.frame_y = self.frame_y
        d.frame_x = self.frame_x
        d.curve_steps = self.curve_steps
        d.shape = self.shape
        d.fixed = self.fixed
        d.pivot = self.pivot
        d.angle_y = self.angle_y
        d.side_material = self.side_material
        d.handle = self.handle
        d.handle_model = self.handle_model
        d.handle_altitude = self.handle_altitude
        d.enable_glass = self.enable_glass
        # Link object into scene
        self.link_object_to_scene(context, o)
        
        # select and make active
        self.select_object(context, o, True)
        m = o.archipack_material.add()
        m.category = "window"
        m.material = self.material
        o.show_transparent = True
        o.lock_location = (False, True, True)
        o.lock_rotation = (False, True, False)
        o.lock_scale = (True, True, True)
        d.update(context)
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            o = self.create(context)
            # select and make active
            self.select_object(context, o, True)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_window_shutter(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.window_shutter"
    bl_label = "Window shutter"
    bl_description = "Window shutter"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    center = FloatVectorProperty(
            subtype='XYZ'
            )
    origin = FloatVectorProperty(
            subtype='XYZ'
            )
    size = FloatVectorProperty(
            subtype='XYZ'
            )
    radius = FloatVectorProperty(
            subtype='XYZ'
            )
    angle_y = FloatProperty(
            name='angle',
            unit='ROTATION',
            subtype='ANGLE',
            min=-1.5, max=1.5,
            default=0, precision=2,
            description='angle'
            )
    depth = FloatProperty(
            name='Depth',
            min=0, max=100,
            default=0.06, precision=2,
            description='frame depth'
            )
    border = FloatProperty(
            name='Width',
            min=0, max=100,
            default=0.06, precision=2,
            description='frame width'
            )
    curve_steps = IntProperty(
            name="curve steps",
            min=1,
            max=128,
            default=16
            )
    shape = EnumProperty(
            name='Shape',
            items=(
                ('RECTANGLE', 'Rectangle', '', 0),
                ('ROUND', 'Top Round', '', 1),
                ('ELLIPSIS', 'Top Elliptic', '', 2),
                ('QUADRI', 'Top oblique', '', 3),
                ('CIRCLE', 'Full circle', '', 4)
                ),
            default='RECTANGLE'
            )
    pivot = FloatProperty(
            name='pivot',
            min=-1, max=1,
            default=-1, precision=2,
            description='pivot'
            )
    material = StringProperty(
            name="material",
            default=""
            )
    offset = FloatProperty(
            name='offset',
            default=0, precision=2,
            description='x offset'
            )
    hinge_enable = BoolProperty(
            name="Hinge",
            default=False
            )
    hinge_count = IntProperty(
            name="#Hinge",
            min=0,
            max=4,
            default=2
            )
    hinge_space = FloatProperty(
            name='space',
            default=0, precision=2,
            description='Vertical space for hinges'
            )
    hinge_size = FloatProperty(
            name='size',
            default=0.03, precision=2,
            description='Vertical size of hinges'
            )

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        m = bpy.data.meshes.new("Shutter")
        o = bpy.data.objects.new("Shutter", m)
        d = m.archipack_window_shutter.add()
        d.center = self.center
        d.origin = self.origin
        d.size = self.size
        d.radius = self.radius
        d.depth = self.depth
        d.border = self.border
        d.curve_steps = self.curve_steps
        d.shape = self.shape
        d.pivot = self.pivot
        d.offset = self.offset
        d.angle_y = self.angle_y
        d.hinge_enable = self.hinge_enable
        d.hinge_count = self.hinge_count
        d.hinge_space = self.hinge_space
        d.hinge_size = self.hinge_size
        # Link object into scene
        self.link_object_to_scene(context, o)
        # select and make active
        self.select_object(context, o, True)
        m = o.archipack_material.add()
        m.category = "window"
        m.material = self.material
        o.lock_location = (False, True, True)
        o.lock_rotation = (False, True, False)
        o.lock_scale = (True, True, True)
        d.update(context)
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            o = self.create(context)
            # select and make active
            self.select_object(context, o, True)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to load / save presets
# ------------------------------------------------------------------

class ARCHIPACK_OT_window_preset_draw(PresetMenuOperator, Operator):
    bl_description = "Choose a preset and draw windows over wall"
    bl_idname = "archipack.window_preset_draw"
    bl_label = "Window Presets"
    preset_subdir = "archipack_window"
    
    
class ARCHIPACK_OT_window_preset_menu(PresetMenuOperator, Operator):
    bl_description = "Show Window Presets"
    bl_idname = "archipack.window_preset_menu"
    bl_label = "Window Presets"
    preset_subdir = "archipack_window"


class ARCHIPACK_OT_window_preset(ArchipackPreset, Operator):
    """Add a Window Preset"""
    bl_idname = "archipack.window_preset"
    bl_label = "Add Window Preset"
    preset_menu = "ARCHIPACK_OT_window_preset_menu"

    @property
    def blacklist(self):
        # 'x', 'y', 'z', 'altitude', 'window_shape'
        return ['manipulators']


def register():
    bpy.utils.register_class(archipack_window_panelrow)
    bpy.utils.register_class(archipack_window_panel)
    Mesh.archipack_window_panel = CollectionProperty(type=archipack_window_panel)
    bpy.utils.register_class(ARCHIPACK_PT_window_panel)
    bpy.utils.register_class(ARCHIPACK_OT_window_panel)

    bpy.utils.register_class(archipack_window_shutter)
    Mesh.archipack_window_shutter = CollectionProperty(type=archipack_window_shutter)
    bpy.utils.register_class(ARCHIPACK_PT_window_shutter)
    bpy.utils.register_class(ARCHIPACK_OT_window_shutter)

    bpy.utils.register_class(archipack_window)
    Mesh.archipack_window = CollectionProperty(type=archipack_window)
    bpy.utils.register_class(ARCHIPACK_OT_window_preset_menu)
    bpy.utils.register_class(ARCHIPACK_OT_window_preset_draw)
    bpy.utils.register_class(ARCHIPACK_PT_window)
    bpy.utils.register_class(ARCHIPACK_OT_window)
    bpy.utils.register_class(ARCHIPACK_OT_window_preset)
    bpy.utils.register_class(ARCHIPACK_OT_window_draw)
    bpy.utils.register_class(ARCHIPACK_OT_window_portals)


def unregister():
    bpy.utils.unregister_class(archipack_window_panelrow)
    bpy.utils.unregister_class(archipack_window_panel)
    bpy.utils.unregister_class(ARCHIPACK_PT_window_panel)
    del Mesh.archipack_window_panel
    bpy.utils.unregister_class(ARCHIPACK_OT_window_panel)

    bpy.utils.unregister_class(archipack_window_shutter)
    bpy.utils.unregister_class(ARCHIPACK_PT_window_shutter)
    del Mesh.archipack_window_shutter
    bpy.utils.unregister_class(ARCHIPACK_OT_window_shutter)

    bpy.utils.unregister_class(archipack_window)
    del Mesh.archipack_window
    bpy.utils.unregister_class(ARCHIPACK_OT_window_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_OT_window_preset_draw)
    bpy.utils.unregister_class(ARCHIPACK_PT_window)
    bpy.utils.unregister_class(ARCHIPACK_OT_window)
    bpy.utils.unregister_class(ARCHIPACK_OT_window_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_window_draw)
    bpy.utils.unregister_class(ARCHIPACK_OT_window_portals)
