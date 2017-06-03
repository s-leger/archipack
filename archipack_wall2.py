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
import bpy
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty, BoolProperty, IntProperty, StringProperty,
    FloatVectorProperty, CollectionProperty, EnumProperty
)
from .bmesh_utils import BmeshEdit as bmed
from .materialutils import MaterialUtils
from mathutils import Vector, Matrix
from mathutils.geometry import (
    intersect_line_plane,
    interpolate_bezier
    )
from bpy_extras.view3d_utils import (
    region_2d_to_origin_3d,
    region_2d_to_vector_3d
    )
from math import sin, cos, pi, atan2
from .archipack_manipulator import (
    Manipulable, archipack_manipulator,
    GlPolygon, GlPolyline,
    GlLine, GlText, FeedbackPanel
    )
from .archipack_2d import Line, Arc
from .archipack_snap import snap_point
from .archipack_keymaps import Keymaps


class Wall():
    def __init__(self, last, wall_z, z, t, flip):
        self.z = z
        self.wall_z = wall_z
        self.t = t
        self.flip = flip
        self.z_step = len(z)

    def get_z(self, t):
        t0 = self.t[0]
        z0 = self.z[0]
        for i in range(1, self.z_step):
            t1 = self.t[i]
            z1 = self.z[i]
            if t <= t1:
                return z0 + (t - t0) / (t1 - t0) * (z1 - z0)
            t0, z0 = t1, z1
        return self.z[-1]

    def make_faces(self, i, f, faces):
        if i < self.n_step:
            # 1 3   5 7
            # 0 2   4 6
            if self.flip:
                faces.append((f + 2, f, f + 1, f + 3))
            else:
                faces.append((f, f + 2, f + 3, f + 1))

    def p3d(self, verts, t):
        x, y = self.lerp(t)
        z = self.wall_z + self.get_z(t)
        verts.append((x, y, 0))
        verts.append((x, y, z))

    def make_wall(self, i, verts, faces):
        t = self.t_step[i]
        f = len(verts)
        self.p3d(verts, t)
        self.make_faces(i, f, faces)

    def straight_wall(self, length, da, wall_z, z, t):
        r = self.straight(length).rotate(da)
        return StraightWall(self, r.p, r.v, wall_z, z, t, self.flip)

    def curved_wall(self, a0, da, radius, wall_z, z, t):
        n = self.normal(1).rotate(a0).scale(radius)
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedWall(self, c, radius, a0, da, wall_z, z, t, self.flip)


class StraightWall(Wall, Line):
    def __init__(self, last, p, v, wall_z, z, t, flip):
        Line.__init__(self, p, v)
        Wall.__init__(self, last, wall_z, z, t, flip)

    def param_t(self, step_angle):
        self.t_step = self.t
        self.n_step = len(self.t) - 1


class CurvedWall(Wall, Arc):
    def __init__(self, last, c, radius, a0, da, wall_z, z, t, flip):
        Arc.__init__(self, c, radius, a0, da)
        Wall.__init__(self, last, wall_z, z, t, flip)

    def param_t(self, step_angle):
        t_step, n_step = self.steps_by_angle(step_angle)
        self.t_step = list(sorted([i * t_step for i in range(1, n_step)] + self.t))
        self.n_step = len(self.t_step) - 1


class WallGenerator():
    def __init__(self, parts):
        self.last_type = 'NONE'
        self.segs = []
        self.parts = parts
        self.faces_type = 'NONE'
        self.closed = False

    def add_part(self, type, radius, a0, da, length, wall_z, part_z, part_t, n_splits, flip):

        # TODO:
        # refactor this part (height manipulators)
        manip_index = []
        if len(self.segs) < 1:
            s = None
            z = [part_z[0]]
            manip_index.append(0)
        else:
            s = self.segs[-1]
            z = [s.z[-1]]

        t_cur = 0
        z_last = n_splits - 1
        t = [0]

        for i in range(n_splits):
            t_try = t[-1] + part_t[i]
            if t_try == t_cur:
                continue
            if t_try <= 1:
                t_cur = t_try
                t.append(t_cur)
                z.append(part_z[i])
                manip_index.append(i)
            else:
                z_last = i
                break

        if t_cur < 1:
            t.append(1)
            manip_index.append(z_last)
            z.append(part_z[z_last])

        # start a new wall
        if s is None:
            if type == 'S_WALL':
                p = Vector((0, 0))
                v = length * Vector((cos(a0), sin(a0)))
                s = StraightWall(s, p, v, wall_z, z, t, flip)
            elif type == 'C_WALL':
                c = -radius * Vector((cos(a0), sin(a0)))
                s = CurvedWall(s, c, radius, a0, da, wall_z, z, t, flip)
        else:
            if type == 'S_WALL':
                s = s.straight_wall(length, a0, wall_z, z, t)
            elif type == 'C_WALL':
                s = s.curved_wall(a0, da, radius, wall_z, z, t)

        self.segs.append(s)
        self.last_type = type

        return manip_index

    def close(self, closed):
        # Make last segment implicit closing one
        if closed:
            part = self.parts[-1]
            w = self.segs[-1]
            dp = self.segs[0].p0 - self.segs[-1].p0
            if "C_" in part.type:
                dw = (w.p1 - w.p0)
                w.r = part.radius / dw.length * dp.length
                # angle pt - p0        - angle p0 p1
                da = atan2(dp.y, dp.x) - atan2(dw.y, dw.x)
                a0 = w.a0 + da
                if a0 > pi:
                    a0 -= 2 * pi
                if a0 < -pi:
                    a0 += 2 * pi
                w.a0 = a0
            else:
                w.v = dp

    def make_wall(self, step_angle, flip, closed, verts, faces):

        # swap manipulators so they always face outside
        side = 1
        if flip:
            side = -1

        # Make last segment implicit closing one

        nb_segs = len(self.segs)
        if closed:
            nb_segs -= 1

        for i, wall in enumerate(self.segs):

            manipulators = self.parts[i].manipulators

            p0 = wall.p0.to_3d()
            p1 = wall.p1.to_3d()

            # angle from last to current segment
            if i > 0:

                if i < nb_segs:
                    manipulators[0].type_key = 'ANGLE'
                else:
                    manipulators[0].type_key = 'DUMB_ANGLE'

                v0 = self.segs[i - 1].straight(-side, 1).v.to_3d()
                v1 = wall.straight(side, 0).v.to_3d()
                manipulators[0].set_pts([p0, v0, v1])

            if type(wall).__name__ == "StraightWall":
                # segment length
                manipulators[1].type_key = 'SIZE'
                manipulators[1].prop1_name = "length"
                manipulators[1].set_pts([p0, p1, (side, 0, 0)])
            else:
                # segment radius + angle
                # scale to fix overlap with drag
                v0 = side * (wall.p0 - wall.c).to_3d()
                v1 = side * (wall.p1 - wall.c).to_3d()
                scale = 1.0 + (0.5 / v0.length)
                manipulators[1].type_key = 'ARC_ANGLE_RADIUS'
                manipulators[1].prop1_name = "da"
                manipulators[1].prop2_name = "radius"
                manipulators[1].set_pts([wall.c.to_3d(), scale * v0, scale * v1])

            # snap manipulator, dont change index !
            manipulators[2].set_pts([p0, p1, (1, 0, 0)])

            # dumb, segment index
            z = Vector((0, 0, 0.75 * wall.wall_z))
            manipulators[3].set_pts([p0 + z, p1 + z, (1, 0, 0)])

            wall.param_t(step_angle)
            if i < nb_segs:
                for j in range(wall.n_step + 1):
                    wall.make_wall(j, verts, faces)
            else:
                # last segment
                for j in range(wall.n_step):
                    print("%s" % (wall.n_step))
                    # wall.make_wall(j, verts, faces)

    def debug(self, verts):
        for wall in self.segs:
            for i in range(33):
                x, y = wall.lerp(i / 32)
                verts.append((x, y, 0))


def update(self, context):
    self.update(context)


def update_childs(self, context):
    self.update(context, update_childs=True, manipulable_refresh=True)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


def set_splits(self, value):
    if self.n_splits != value:
        self.auto_update = False
        self._set_t(value)
        self.auto_update = True
        self.n_splits = value
    return None


def get_splits(self):
    return self.n_splits


def update_type(self, context):

    d = self.find_in_selection(context)

    if d is not None and d.auto_update:

        d.auto_update = False
        idx = 0
        for i, part in enumerate(d.parts):
            if part == self:
                idx = i
                break
        a0 = 0
        if idx > 0:
            g = d.get_generator()
            w0 = g.segs[idx - 1]
            a0 = w0.straight(1).angle
            if "C_" in self.type:
                w = w0.straight_wall(self.length, self.a0, d.z, self.z, self.t)
            else:
                w = w0.curved_wall(self.a0, self.da, self.radius, d.z, self.z, self.t)
        else:
            g = WallGenerator(None)
            if "C_" in self.type:
                g.add_part("S_WALL", self.radius, self.a0, self.da,
                    self.length, d.z, self.z, self.t, self.n_splits, d.flip)
            else:
                g.add_part("C_WALL", self.radius, self.a0, self.da,
                    self.length, d.z, self.z, self.t, self.n_splits, d.flip)
            w = g.segs[0]
        # w0 - w - w1
        dp = w.p1 - w.p0
        if "C_" in self.type:
            self.radius = 0.5 * dp.length
            self.da = pi
            a0 = atan2(dp.y, dp.x) - pi / 2 - a0
        else:
            self.length = dp.length
            a0 = atan2(dp.y, dp.x) - a0

        if a0 > pi:
            a0 -= 2 * pi
        if a0 < -pi:
            a0 += 2 * pi
        self.a0 = a0

        if idx + 1 < d.n_parts:
            # adjust rotation of next part
            part1 = d.parts[idx + 1]
            if "C_" in self.type:
                a0 = part1.a0 - pi / 2
            else:
                a0 = part1.a0 + w.straight(1).angle - atan2(dp.y, dp.x)

            if a0 > pi:
                a0 -= 2 * pi
            if a0 < -pi:
                a0 += 2 * pi
            part1.a0 = a0

        d.auto_update = True


class archipack_wall2_part(PropertyGroup):
    type = EnumProperty(
            items=(
                ('S_WALL', 'Straight', '', 0),
                ('C_WALL', 'Curved', '', 1)
                ),
            default='S_WALL',
            update=update_type
            )
    length = FloatProperty(
            name="length",
            min=0.01,
            default=2.0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    radius = FloatProperty(
            name="radius",
            min=0.5,
            default=0.7,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    a0 = FloatProperty(
            name="angle depart",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    da = FloatProperty(
            name="angle",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    z = FloatVectorProperty(
            name="height",
            min=0,
            default=[
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0
            ],
            size=31,
            update=update
            )
    t = FloatVectorProperty(
            name="position",
            min=0,
            max=1,
            default=[
                1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1
            ],
            size=31,
            update=update
            )
    splits = IntProperty(
        name="splits",
        default=1,
        min=1,
        max=31,
        get=get_splits, set=set_splits
        )
    n_splits = IntProperty(
        name="splits",
        default=1,
        min=1,
        max=31,
        update=update
        )
    auto_update = BoolProperty(default=True)
    manipulators = CollectionProperty(type=archipack_manipulator)
    # ui related
    expand = BoolProperty(default=False)

    def _set_t(self, splits):
        t = 1 / splits
        for i in range(splits):
            self.t[i] = t

    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = ARCHIPACK_PT_wall2.params(o)
            if props:
                for part in props.parts:
                    if part == self:
                        return props
        return None

    def update(self, context, manipulable_refresh=False):
        if not self.auto_update:
            return
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context, manipulable_refresh)

    def draw(self, layout, context, index):

        row = layout.row(align=True)
        if self.expand:
            row.prop(self, 'expand', icon="TRIA_DOWN", icon_only=True, text="Part " + str(index + 1), emboss=False)
        else:
            row.prop(self, 'expand', icon="TRIA_RIGHT", icon_only=True, text="Part " + str(index + 1), emboss=False)

        row.prop(self, "type", text="")

        if self.expand:
            row = layout.row(align=True)
            row.operator("archipack.wall2_insert", text="Split").index = index
            row.operator("archipack.wall2_remove", text="Remove").index = index
            if self.type == 'C_WALL':
                row = layout.row()
                row.prop(self, "radius")
                row = layout.row()
                row.prop(self, "da")
            else:
                row = layout.row()
                row.prop(self, "length")
            row = layout.row()
            row.prop(self, "a0")
            row = layout.row()
            row.prop(self, "splits")
            for split in range(self.n_splits):
                row = layout.row()
                row.prop(self, "z", text="alt", index=split)
                row.prop(self, "t", text="pos", index=split)


class archipack_wall2_child(PropertyGroup):
    # Size  Loc
    # Delta Loc
    manipulators = CollectionProperty(type=archipack_manipulator)
    child_name = StringProperty()
    wall_idx = IntProperty()
    pos = FloatVectorProperty(subtype='XYZ')
    flip = BoolProperty(default=False)

    def get_child(self, context):
        d = None
        child = context.scene.objects.get(self.child_name)
        if child is not None and child.data is not None:
            if 'archipack_window' in child.data:
                d = child.data.archipack_window[0]
            elif 'archipack_door' in child.data:
                d = child.data.archipack_door[0]
        return child, d


class archipack_wall2(Manipulable, PropertyGroup):
    parts = CollectionProperty(type=archipack_wall2_part)
    n_parts = IntProperty(
            name="parts",
            min=1,
            max=32,
            default=1, update=update_manipulators
            )
    step_angle = FloatProperty(
            name="step angle",
            min=1 / 180 * pi,
            max=pi,
            default=6 / 180 * pi,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    width = FloatProperty(
            name="width",
            min=0.01,
            default=0.2,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    z = FloatProperty(
            name='height',
            min=0.1,
            default=2.7, precision=2,
            unit='LENGTH', subtype='DISTANCE',
            description='height', update=update,
            )
    x_offset = FloatProperty(
            name="x offset",
            min=-1, max=1,
            default=-1, precision=2, step=1,
            update=update
            )
    radius = FloatProperty(
            name="radius",
            min=0.5,
            default=0.7,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    da = FloatProperty(
            name="angle",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    flip = BoolProperty(default=False, update=update_childs)
    closed = BoolProperty(
            default=False,
            name="Close",
            update=update_manipulators
            )
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )
    realtime = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            name="RealTime",
            description="Relocate childs in realtime"
            )
    # dumb manipulators to show sizes between childs
    childs_manipulators = CollectionProperty(type=archipack_manipulator)
    childs = CollectionProperty(type=archipack_wall2_child)

    def insert_part(self, context, where):
        self.manipulable_disable(context)
        self.auto_update = False
        # the part we do split
        part_0 = self.parts[where]
        part_0.length /= 2
        part_0.da /= 2
        self.parts.add()
        part_1 = self.parts[len(self.parts) - 1]
        part_1.type = part_0.type
        part_1.length = part_0.length
        part_1.da = part_0.da
        part_1.a0 = 0
        # move after current one
        self.parts.move(len(self.parts) - 1, where + 1)
        self.n_parts += 1
        self.setup_parts_manipulators()
        self.auto_update = True

    def add_part(self, context, length):
        self.manipulable_disable(context)
        self.auto_update = False
        p = self.parts.add()
        p.length = length
        self.n_parts += 1
        self.setup_parts_manipulators()
        self.auto_update = True
        return p

    def remove_part(self, context, where):
        self.manipulable_disable(context)
        self.auto_update = False
        # preserve shape
        # using generator
        if where > 0:
            g = self.get_generator()
            w = g.segs[where - 1]
            dp = g.segs[where].p1 - w.p0
            if where + 1 < self.n_parts:
                a0 = g.segs[where + 1].straight(1).angle - atan2(dp.y, dp.x)
                part = self.parts[where + 1]
                if a0 > pi:
                    a0 -= 2 * pi
                if a0 < -pi:
                    a0 += 2 * pi
                part.a0 = a0
            part = self.parts[where - 1]
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

        self.parts.remove(where)
        self.n_parts -= 1
        # fix snap manipulators index
        self.setup_parts_manipulators()
        self.auto_update = True

    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        active = context.active_object
        selected = [o for o in context.selected_objects]
        for o in selected:
            if ARCHIPACK_PT_wall2.params(o) == self:
                return active, selected, o
        return active, selected, None

    def get_generator(self):
        # print("get_generator")
        g = WallGenerator(self.parts)
        for part in self.parts:
            g.add_part(part.type, part.radius, part.a0, part.da, part.length,
                self.z, part.z, part.t, part.n_splits, self.flip)
        g.close(self.closed)
        return g

    def update_parts(self, o, update_childs=False):
        # print("update_parts")
        # remove rows
        # NOTE:
        # n_parts+1
        # as last one is end point of last segment or closing one
        row_change = False
        for i in range(len(self.parts), self.n_parts, -1):
            row_change = True
            self.parts.remove(i - 1)

        # add rows
        for i in range(len(self.parts), self.n_parts):
            row_change = True
            self.parts.add()

        self.setup_parts_manipulators()

        g = self.get_generator()

        if o is not None and (row_change or update_childs):
            self.setup_childs(o, g)

        return g

    def setup_parts_manipulators(self):
        for i in range(self.n_parts):
            p = self.parts[i]
            n_manips = len(p.manipulators)
            if n_manips < 1:
                s = p.manipulators.add()
                s.type_key = "ANGLE"
                s.prop1_name = "a0"
            if n_manips < 2:
                s = p.manipulators.add()
                s.type_key = "SIZE"
                s.prop1_name = "length"
            if n_manips < 3:
                s = p.manipulators.add()
                s.type_key = 'WALL_SNAP'
                s.prop1_name = str(i)
                s.prop2_name = 'z'
            if n_manips < 4:
                s = p.manipulators.add()
                s.type_key = 'DUMB_STRING'
                s.prop1_name = str(i + 1)
            p.manipulators[2].prop1_name = str(i)
            p.manipulators[3].prop1_name = str(i + 1)

        self.manipulable_selectable = True

    def interpolate_bezier(self, pts, wM, p0, p1, resolution):
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
                    resolution + 1)
                for i in range(resolution):
                    pts.append(seg[i].to_3d())

    def from_spline(self, wM, resolution, spline):
        pts = []
        if spline.type == 'POLY':
            pts = [wM * p.co.to_3d() for p in spline.points]
            if spline.use_cyclic_u:
                pts.append(pts[0])
        elif spline.type == 'BEZIER':
            points = spline.bezier_points
            for i in range(1, len(points)):
                p0 = points[i - 1]
                p1 = points[i]
                self.interpolate_bezier(pts, wM, p0, p1, resolution)
            pts.append(wM * points[-1].co)
            if spline.use_cyclic_u:
                p0 = points[-1]
                p1 = points[0]
                self.interpolate_bezier(pts, wM, p0, p1, resolution)
                pts.append(pts[0])

        self.auto_update = False

        self.n_parts = len(pts) - 1
        self.update_parts(None)

        p0 = pts.pop(0)
        a0 = 0
        for i, p1 in enumerate(pts):
            dp = p1 - p0
            da = atan2(dp.y, dp.x) - a0
            if da > pi:
                da -= 2 * pi
            if da < -pi:
                da += 2 * pi
            p = self.parts[i]
            p.length = dp.to_2d().length
            p.dz = dp.z
            p.a0 = da
            a0 += da
            p0 = p1
        self.closed = True
        self.auto_update = True

    def update(self, context, manipulable_refresh=False, update_childs=False):
        # print("update manipulable_refresh:%s" % (manipulable_refresh))
        """

        """
        if not self.auto_update:
            return

        active, selected, o = self.find_in_selection(context)

        if o is None:
            return

        if manipulable_refresh:
            # prevent crash by removing all manipulators refs to datablock before changes
            self.manipulable_disable(context)

        verts = []
        faces = []

        g = self.update_parts(o, update_childs)
        # print("make_wall")
        g.make_wall(self.step_angle, self.flip, self.closed, verts, faces)

        if self.closed:
            f = len(verts)
            if self.flip:
                faces.append((0, f - 2, f - 1, 1))
            else:
                faces.append((f - 2, 0, 1, f - 1))

        # print("buildmesh")
        bmed.buildmesh(context, o, verts, faces, matids=None, uvs=None, weld=True)
        bpy.ops.object.shade_smooth()

        side = 1
        if self.flip:
            side = -1
        # Width
        offset = side * (0.5 * self.x_offset) * self.width
        self.manipulators[0].set_pts([
            g.segs[0].sized_normal(0, offset + 0.5 * side * self.width).v.to_3d(),
            g.segs[0].sized_normal(0, offset - 0.5 * side * self.width).v.to_3d(),
            (-side, 0, 0)
            ])

        # Parts COUNTER
        self.manipulators[1].set_pts([g.segs[-1].lerp(1.1).to_3d(),
            g.segs[-1].lerp(1.1 + 0.5 / g.segs[-1].length).to_3d(),
            (-side, 0, 0)
            ])

        # Height
        self.manipulators[2].set_pts([
            (0, 0, 0),
            (0, 0, self.z),
            (-1, 0, 0)
            ], normal=g.segs[0].straight(side, 0).v.to_3d())

        if self.realtime:
            # update child location and size
            self.relocate_childs(context, o, g)
            # store gl points
            self.update_childs(context, o, g)
        else:
            bpy.ops.archipack.wall2_throttle_update(name=o.name)

        modif = o.modifiers.get('Wall')
        if modif is None:
            modif = o.modifiers.new('Wall', 'SOLIDIFY')
            modif.use_quality_normals = True
            modif.use_even_offset = True
            modif.material_offset_rim = 2
            modif.material_offset = 1

        modif.thickness = self.width
        modif.offset = self.x_offset

        if manipulable_refresh:
            # print("manipulable_refresh=True")
            self.manipulable_refresh = True

        # restore context
        try:
            for o in selected:
                o.select = True
        except:
            pass

        active.select = True
        context.scene.objects.active = active

    # manipulable children objects like windows and doors
    def child_partition(self, array, begin, end):
        pivot = begin
        for i in range(begin + 1, end + 1):
            # wall idx
            if array[i][1] < array[begin][1]:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
            # param t on the wall
            elif array[i][1] == array[begin][1] and array[i][4] <= array[begin][4]:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
        array[pivot], array[begin] = array[begin], array[pivot]
        return pivot

    def sort_child(self, array, begin=0, end=None):
        # print("sort_child")
        if end is None:
            end = len(array) - 1

        def _quicksort(array, begin, end):
            if begin >= end:
                return
            pivot = self.child_partition(array, begin, end)
            _quicksort(array, begin, pivot - 1)
            _quicksort(array, pivot + 1, end)
        return _quicksort(array, begin, end)

    def add_child(self, name, wall_idx, pos, flip):
        # print("add_child %s %s" % (name, wall_idx))
        c = self.childs.add()
        c.child_name = name
        c.wall_idx = wall_idx
        c.pos = pos
        c.flip = flip
        m = c.manipulators.add()
        m.type_key = 'DELTA_LOC'
        m.prop1_name = "x"
        m = c.manipulators.add()
        m.type_key = 'SNAP_SIZE_LOC'
        m.prop1_name = "x"
        m.prop2_name = "x"

    def setup_childs(self, o, g):
        """
            Store childs
            create manipulators
            call after a boolean oop
        """
        # print("setup_childs")

        self.childs.clear()
        self.childs_manipulators.clear()
        if o.parent is None:
            return
        wall_with_childs = [0 for i in range(self.n_parts)]
        relocate = []
        dmax = 2 * self.width
        witM = o.matrix_world.inverted()
        itM = witM * o.parent.matrix_world
        rM = itM.to_3x3()
        for child in o.parent.children:
            if (child != o and
                    'archipack_robusthole' not in child and
                    'archipack_hole' not in child):
                tM = child.matrix_world.to_3x3()
                pt = (itM * child.location).to_2d()
                dir_y = (rM * tM * Vector((0, 1, 0))).to_2d()
                for wall_idx, wall in enumerate(g.segs):
                    # may be optimized with a bound check
                    res, d, t = wall.point_sur_segment(pt)
                    # outside is on the right side of the wall
                    #  p1
                    #  |-- x
                    #  p0
                    dir = wall.normal(t).v.normalized()
                    if res and t > 0 and t < 1 and abs(d) < dmax:
                        wall_with_childs[wall_idx] = 1
                        m = self.childs_manipulators.add()
                        m.type_key = 'DUMB_SIZE'
                        # always make window points outside
                        if child.data is not None and "archipack_window" in child.data:
                            flip = self.flip
                        else:
                            # let door point where user want
                            flip = (dir - dir_y).length > 0.5
                        # store z in wall space
                        relocate.append((
                            child.name,
                            wall_idx,
                            (t * wall.length, d, (itM * child.location).z),
                            flip,
                            t))

        self.sort_child(relocate)
        for child in relocate:
            name, wall_idx, pos, flip, t = child
            self.add_child(name, wall_idx, pos, flip)

        # add a dumb size from last child to end of wall segment
        for i in range(sum(wall_with_childs)):
            m = self.childs_manipulators.add()
            m.type_key = 'DUMB_SIZE'

    def relocate_childs(self, context, o, g):
        """
            Move and resize childs after wall edition
        """
        # print("relocate_childs")

        w = -self.x_offset * self.width
        if self.flip:
            w = -w
        tM = o.matrix_world
        for child in self.childs:
            c, d = child.get_child(context)
            if c is None:
                continue
            t = child.pos.x / g.segs[child.wall_idx].length
            n = g.segs[child.wall_idx].normal(t)
            rx, ry = -n.v.normalized()
            rx, ry = ry, -rx
            if child.flip:
                rx, ry = -rx, -ry

            if d is not None:
                c.select = True
                d.auto_update = False
                d.flip = child.flip
                d.y = self.width
                d.auto_update = True
                c.select = False
                x, y = n.p - (0.5 * w * n.v.normalized())
            else:
                x, y = n.p - (child.pos.y * n.v.normalized())

            context.scene.objects.active = o
            # preTranslate
            c.matrix_world = tM * Matrix([
                [rx, -ry, 0, x],
                [ry, rx, 0, y],
                [0, 0, 1, child.pos.z],
                [0, 0, 0, 1]
            ])

    def update_childs(self, context, o, g):
        """
            setup gl points for childs
        """
        # print("update_childs")

        if o.parent is None:
            return

        # swap manipulators so they always face outside
        manip_side = 1
        if self.flip:
            manip_side = -1

        itM = o.matrix_world.inverted() * o.parent.matrix_world
        m_idx = 0
        for wall_idx, wall in enumerate(g.segs):
            p0 = wall.lerp(0)
            wall_has_childs = False
            for child in self.childs:
                if child.wall_idx == wall_idx:
                    c, d = child.get_child(context)
                    if d is not None:
                        # child is either a window or a door
                        wall_has_childs = True
                        dt = 0.5 * d.x / wall.length
                        pt = (itM * c.location).to_2d()
                        res, y, t = wall.point_sur_segment(pt)
                        child.pos = (wall.length * t, y, child.pos.z)
                        p1 = wall.lerp(t - dt)
                        # dumb size between childs
                        self.childs_manipulators[m_idx].set_pts([
                            (p0.x, p0.y, 0),
                            (p1.x, p1.y, 0),
                            (manip_side * 0.5, 0, 0)])
                        m_idx += 1
                        x, y = 0.5 * d.x, -self.x_offset * 0.5 * d.y

                        if child.flip:
                            side = -manip_side
                        else:
                            side = manip_side

                        # delta loc
                        child.manipulators[0].set_pts([(-x, side * -y, 0), (x, side * -y, 0), (side, 0, 0)])
                        # loc size
                        child.manipulators[1].set_pts([
                            (-x, side * -y, 0),
                            (x, side * -y, 0),
                            (0.5 * side, 0, 0)])
                        p0 = wall.lerp(t + dt)
            p1 = wall.lerp(1)
            if wall_has_childs:
                # dub size after all childs
                self.childs_manipulators[m_idx].set_pts([
                    (p0.x, p0.y, 0),
                    (p1.x, p1.y, 0),
                    (manip_side * 0.5, 0, 0)])
                m_idx += 1

    def manipulate_childs(self, context):
        """
            setup child manipulators
        """
        # print("manipulate_childs")
        for wall_idx in range(self.n_parts):
            for child in self.childs:
                if child.wall_idx == wall_idx:
                    c, d = child.get_child(context)
                    if d is not None:
                        # delta loc
                        self.manip_stack.append(child.manipulators[0].setup(context, c, d, self.manipulate_callback))
                        # loc size
                        self.manip_stack.append(child.manipulators[1].setup(context, c, d, self.manipulate_callback))

    def manipulate_callback(self, context, o=None, manipulator=None):
        found = False
        if o.parent is not None:
            for c in o.parent.children:
                if (c.data is not None and
                        'archipack_wall2' in c.data and
                        c.data.archipack_wall2[0] == self):
                    context.scene.objects.active = c
                    found = True
                    break
        if found:
            self.manipulable_manipulate(context, manipulator=manipulator)

    def manipulable_manipulate(self, context, event=None, manipulator=None):
        type_name = type(manipulator).__name__
        # print("manipulable_manipulate %s" % (type_name))
        if type_name in [
                'DeltaLocationManipulator',
                'SizeLocationManipulator',
                'SnapSizeLocationManipulator'
                ]:
            # update manipulators pos of childs
            o = context.active_object
            if o.parent is None:
                return
            g = self.get_generator()
            itM = o.matrix_world.inverted() * o.parent.matrix_world
            for child in self.childs:
                c, d = child.get_child(context)
                if d is not None:
                    wall = g.segs[child.wall_idx]
                    pt = (itM * c.location).to_2d()
                    res, d, t = wall.point_sur_segment(pt)
                    child.pos = (t * wall.length, d, child.pos.z)
            # update childs manipulators
            self.update_childs(context, o, g)

    def manipulable_release(self, context):
        """
            Override with action to do on mouse release
            eg: big update
        """
        return

    def manipulable_setup(self, context):
        # print("manipulable_setup")
        self.manipulable_disable(context)
        o = context.active_object

        # setup childs manipulators
        self.manipulate_childs(context)
        nb_segs = self.n_parts
        if self.closed:
            nb_segs -= 1

        # update manipulators on version change
        self.setup_parts_manipulators()

        for i, part in enumerate(self.parts):

            if i < self.n_parts:
                if i > 0:
                    # start angle
                    self.manip_stack.append(part.manipulators[0].setup(context, o, part))

                # length / radius + angle
                self.manip_stack.append(part.manipulators[1].setup(context, o, part))

            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))

            # segment index
            self.manip_stack.append(part.manipulators[3].setup(context, o, self))

        # height as per segment will be here when done

        # TODO:
        # add last segment snap manipulator at end of the segment

        # width + counter
        for m in self.manipulators:
            self.manip_stack.append(m.setup(context, o, self))

        # dumb between childs
        for m in self.childs_manipulators:
            self.manip_stack.append(m.setup(context, o, self))

    def manipulable_exit(self, context):
        """
            Override with action to do when modal exit
        """
        return

    def manipulable_invoke(self, context):
        """
            call this in operator invoke()
        """
        # print("manipulable_invoke")
        if self.manipulate_mode:
            self.manipulable_disable(context)
            self.manipulate_mode = False
            return False

        self.manip_stack = []
        o = context.active_object
        g = self.get_generator()
        # setup childs manipulators
        self.setup_childs(o, g)
        # store gl points
        self.update_childs(context, o, g)
        # dont do anything ..
        # self.manipulable_release(context)
        self.manipulable_setup(context)
        self.manipulate_mode = True

        self._manipulable_invoke(context)

        return True


# Update throttle (smell hack here)
# use 2 globals to store a timer and state of update_action
update_timer = None
update_timer_updating = False


class ARCHIPACK_OT_wall2_throttle_update(Operator):
    bl_idname = "archipack.wall2_throttle_update"
    bl_label = "Update childs with a delay"

    name = StringProperty()

    def modal(self, context, event):
        global update_timer_updating
        if event.type == 'TIMER' and not update_timer_updating:
            update_timer_updating = True
            o = context.scene.objects.get(self.name)
            # print("delay update of %s" % (self.name))
            if o is not None:
                o.select = True
                context.scene.objects.active = o
                d = o.data.archipack_wall2[0]
                g = d.get_generator()
                # update child location and size
                d.relocate_childs(context, o, g)
                # store gl points
                d.update_childs(context, o, g)
                self.cancel(context)
        return {'PASS_THROUGH'}

    def execute(self, context):
        global update_timer
        global update_timer_updating
        if update_timer is not None:
            if update_timer_updating:
                return {'CANCELLED'}
            # reset update_timer so it only occurs once 0.1s after last action
            context.window_manager.event_timer_remove(update_timer)
            update_timer = context.window_manager.event_timer_add(0.1, context.window)
            return {'CANCELLED'}
        update_timer_updating = False
        context.window_manager.modal_handler_add(self)
        update_timer = context.window_manager.event_timer_add(0.1, context.window)
        return {'RUNNING_MODAL'}

    def cancel(self, context):
        global update_timer
        context.window_manager.event_timer_remove(update_timer)
        update_timer = None
        return {'CANCELLED'}


class ARCHIPACK_PT_wall2(Panel):
    bl_idname = "ARCHIPACK_PT_wall2"
    bl_label = "Wall"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    # bl_context = 'object'
    # bl_space_type = 'VIEW_3D'
    # bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    def draw(self, context):
        o = context.object
        prop = ARCHIPACK_PT_wall2.params(o)
        if prop is None:
            return
        layout = self.layout
        row = layout.row(align=True)
        row.operator("archipack.wall2_manipulate", icon='HAND')
        row = layout.row(align=True)
        row.prop(prop, 'realtime')
        box = layout.box()
        box.prop(prop, 'n_parts')
        box.prop(prop, 'step_angle')
        box.prop(prop, 'width')
        box.prop(prop, 'z')
        box.prop(prop, 'flip')
        box.prop(prop, 'x_offset')
        row = layout.row()
        row.prop(prop, "closed")
        for i, part in enumerate(prop.parts):
            box = layout.box()
            part.draw(box, context, i)

    @classmethod
    def params(cls, o):
        try:
            if 'archipack_wall2' not in o.data:
                return False
            else:
                return o.data.archipack_wall2[0]
        except:
            return False

    @classmethod
    def filter(cls, o):
        try:
            if 'archipack_wall2' not in o.data:
                return False
            else:
                return True
        except:
            return False

    @classmethod
    def poll(cls, context):
        o = context.object
        if o is None:
            return False
        return cls.filter(o)


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_wall2(Operator):
    bl_idname = "archipack.wall2"
    bl_label = "Wall"
    bl_description = "Create a Wall"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    auto_manipulate = BoolProperty(default=True)

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        m = bpy.data.meshes.new("Wall")
        o = bpy.data.objects.new("Wall", m)
        d = m.archipack_wall2.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        s = d.manipulators.add()
        s.prop1_name = "width"
        s = d.manipulators.add()
        s.prop1_name = "n_parts"
        s.type_key = 'COUNTER'
        s = d.manipulators.add()
        s.prop1_name = "z"
        s.normal = (0, 1, 0)
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        d.update(context)
        MaterialUtils.add_wall_materials(o)
        # around 12 degree
        m.auto_smooth_angle = 0.20944
        m.use_auto_smooth = True
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = bpy.context.scene.cursor_location
            o.select = True
            context.scene.objects.active = o
            if self.auto_manipulate:
                bpy.ops.archipack.wall2_manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_wall2_from_curve(Operator):
    bl_idname = "archipack.wall2_from_curve"
    bl_label = "Wall curve"
    bl_description = "Create a wall from a curve"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    auto_manipulate = BoolProperty(default=True)

    @classmethod
    def poll(self, context):
        return context.active_object is not None and context.active_object.type == 'CURVE'
    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        curve = context.active_object
        for spline in curve.data.splines:
            bpy.ops.archipack.wall2(auto_manipulate=False)
            o = context.scene.objects.active
            d = o.data.archipack_wall2[0]
            d.from_spline(curve.matrix_world, 12, spline)
            if spline.type == 'POLY':
                pt = spline.points[0].co
            elif spline.type == 'BEZIER':
                pt = spline.bezier_points[0].co
            else:
                pt = Vector((0, 0, 0))
            # pretranslate
            o.matrix_world = curve.matrix_world * Matrix([
                [1, 0, 0, pt.x],
                [0, 1, 0, pt.y],
                [0, 0, 1, pt.z],
                [0, 0, 0, 1]
                ])
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            if o is not None:
                o.select = True
                context.scene.objects.active = o
                if self.auto_manipulate:
                    bpy.ops.archipack.wall2_manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_wall2_from_slab(Operator):
    bl_idname = "archipack.wall2_from_slab"
    bl_label = "Slab -> Wall"
    bl_description = "Create a wall from a slab"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    auto_manipulate = BoolProperty(default=True)

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o is not None and o.data is not None and 'archipack_slab' in o.data
    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        slab = context.active_object
        wd = slab.data.archipack_slab[0]
        bpy.ops.archipack.wall2(auto_manipulate=False)
        o = context.scene.objects.active
        d = o.data.archipack_wall2[0]
        d.auto_update = False
        d.parts.clear()
        d.n_parts = wd.n_parts
        d.closed = True
        if not wd.closed:
            d.n_parts += 1
        for part in wd.parts:
            p = d.parts.add()
            if "S_" in part.type:
                p.type = "S_WALL"
            else:
                p.type = "C_WALL"
            p.length = part.length
            p.radius = part.radius
            p.da = part.da
            p.a0 = part.a0
        o.select = True
        context.scene.objects.active = o
        d.auto_update = True
        # pretranslate
        o.matrix_world = slab.matrix_world.copy()
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.select = True
            context.scene.objects.active = o
            if self.auto_manipulate:
                bpy.ops.archipack.wall2_manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to draw a wall
# ------------------------------------------------------------------


class ARCHIPACK_OT_wall2_draw(Operator):
    bl_idname = "archipack.wall2_draw"
    bl_label = "Draw a Wall"
    bl_description = "Draw a Wall"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    o = None
    state = 'RUNNING'
    flag_create = False
    flag_next = False
    wall_part1 = None
    wall_line1 = None
    line = None
    label = None
    feedback = None
    takeloc = Vector((0, 0, 0))

    def mouse_to_plane(self, context, event):
        """
            convert mouse pos to 3d point over plane defined by origin and normal
        """
        region = context.region
        rv3d = context.region_data
        co2d = (event.mouse_region_x, event.mouse_region_y)
        view_vector_mouse = region_2d_to_vector_3d(region, rv3d, co2d)
        ray_origin_mouse = region_2d_to_origin_3d(region, rv3d, co2d)
        pt = intersect_line_plane(ray_origin_mouse, ray_origin_mouse + view_vector_mouse,
            Vector((0, 0, 0)), Vector((0, 0, 1)), False)
        # fix issue with parallel plane
        if pt is None:
            pt = intersect_line_plane(ray_origin_mouse, ray_origin_mouse + view_vector_mouse,
                Vector((0, 0, 0)), view_vector_mouse, False)
        return pt

    @classmethod
    def poll(cls, context):
        return True

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def draw_callback(self, _self, context):
        self.feedback.draw(context)

    def sp_draw(self, sp, context):
        z = 2.7
        if self.state == 'CREATE':
            p0 = self.takeloc
        else:
            p0 = sp.takeloc

        p1 = sp.placeloc
        delta = p1 - p0
        # print("sp_draw state:%s delta:%s p0:%s p1:%s" % (self.state, delta.length, p0, p1))
        if delta.length == 0:
            return
        self.wall_part1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.wall_line1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.wall_part1.draw(context)
        self.wall_line1.draw(context)
        self.line.p = p0
        self.line.v = delta
        self.label.set_pos(context, self.line.length, self.line.lerp(0.5), self.line.v, normal=Vector((0, 0, 1)))
        self.label.draw(context)
        self.line.draw(context)

    def sp_callback(self, context, event, state, sp):
        # print("sp_callback event %s %s state:%s" % (event.type, event.value, state))

        if state == 'SUCCESS':

            if self.state == 'CREATE':
                takeloc = self.takeloc
                delta = sp.placeloc - self.takeloc
            else:
                takeloc = sp.takeloc
                delta = sp.delta

            old = context.active_object
            if self.o is None:
                bpy.ops.archipack.wall2(auto_manipulate=False)
                o = context.active_object
                o.location = takeloc
                self.o = o
                d = o.data.archipack_wall2[0]
                part = d.parts[0]
                part.length = delta.length
            else:
                o = self.o
                o.select = True
                context.scene.objects.active = o
                d = o.data.archipack_wall2[0]
                # Check for end close to start and close when applicable
                dp = sp.placeloc - o.location
                if dp.length < 0.01:
                    d.closed = True
                    state = 'CANCEL'
                    # return
                part = d.add_part(context, delta.length)

            # print("self.o :%s" % o.name)
            rM = o.matrix_world.inverted().to_3x3()
            g = d.get_generator()
            w = g.segs[-1]
            dp = rM * delta
            da = atan2(dp.y, dp.x) - w.straight(1).angle
            a0 = part.a0 + da
            if a0 > pi:
                a0 -= 2 * pi
            if a0 < -pi:
                a0 += 2 * pi
            part.a0 = a0

            context.scene.objects.active = old
            self.flag_next = True
            context.area.tag_redraw()
            # print("feedback.on:%s" % self.feedback.on)

        self.state = state

    def sp_init(self, context, event, state, sp):
        # print("sp_init event %s %s %s" % (event.type, event.value, state))
        if state == 'SUCCESS':
            self.state = 'RUNNING'
            self.takeloc = sp.placeloc.copy()
            # print("feedback.on:%s" % self.feedback.on)
        elif state == 'CANCEL':
            self.state = state
            return

    def modal(self, context, event):

        context.area.tag_redraw()
        # print("modal event %s %s" % (event.type, event.value))
        # if event.type == 'NONE':
        #    return {'PASS_THROUGH'}

        if self.state == 'STARTING':
            takeloc = self.mouse_to_plane(context, event)
            # print("STARTING")
            snap_point(takeloc=takeloc,
                callback=self.sp_init,
                constraint_axis=(True, True, False),
                release_confirm=True)
            return {'RUNNING_MODAL'}

        elif self.state == 'RUNNING':
            # print("RUNNING")
            self.state = 'CREATE'
            snap_point(takeloc=self.takeloc,
                draw=self.sp_draw,
                callback=self.sp_callback,
                constraint_axis=(True, True, False),
                release_confirm=True)
            return {'RUNNING_MODAL'}

        elif event.type in {'LEFTMOUSE', 'RET', 'NUMPAD_ENTER', 'SPACE'}:

            # print('LEFTMOUSE %s' % (event.value))
            self.feedback.instructions(context, "Draw a wall", "Click & Drag to add a segment", [
                ('CTRL', 'Snap'),
                ('MMBTN', 'Constraint to axis'),
                ('X Y', 'Constraint to axis'),
                ('BACK_SPACE', 'Remove part'),
                ('RIGHTCLICK or ESC', 'exit')
                ])

            if event.value == 'PRESS':

                if self.flag_next:
                    self.flag_next = False
                    o = self.o
                    o.select = True
                    context.scene.objects.active = o
                    d = o.data.archipack_wall2[0]
                    g = d.get_generator()
                    takeloc = o.matrix_world * g.segs[-1].p1.to_3d()
                    o.select = False
                else:
                    takeloc = self.mouse_to_plane(context, event)

                snap_point(takeloc=takeloc,
                    draw=self.sp_draw,
                    callback=self.sp_callback,
                    constraint_axis=(True, True, False),
                    release_confirm=True)

            return {'RUNNING_MODAL'}

        if self.keymap.check(event, self.keymap.undo) or (
                event.type in {'BACK_SPACE'} and event.value == 'RELEASE'
                ):
            if self.o is not None:
                o = self.o
                o.select = True
                context.scene.objects.active = o
                d = o.data.archipack_wall2[0]
                if d.n_parts > 1:
                    d.n_parts -= 1
            return {'RUNNING_MODAL'}

        if self.state == 'CANCEL' or (event.type in {'ESC', 'RIGHTMOUSE'} and
                event.value == 'RELEASE'):

            self.feedback.disable()
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')

            if self.o is not None:
                self.o.select = True
                context.scene.objects.active = self.o
                if bpy.ops.archipack.wall2_manipulate.poll():
                    bpy.ops.archipack.wall2_manipulate('INVOKE_DEFAULT')

            return {'FINISHED'}

        return {'PASS_THROUGH'}

    def invoke(self, context, event):

        if context.mode == "OBJECT":
            self.keymap = Keymaps(context)
            self.wall_part1 = GlPolygon((0.5, 0, 0, 0.2))
            self.wall_line1 = GlPolyline((0.5, 0, 0, 0.8))
            self.line = GlLine()
            self.label = GlText()
            self.feedback = FeedbackPanel()
            self.feedback.instructions(context, "Draw a wall", "Click & Drag to start", [
                ('CTRL', 'Snap'),
                ('MMBTN', 'Constraint to axis'),
                ('X Y', 'Constraint to axis'),
                ('RIGHTCLICK or ESC', 'exit without change')
                ])
            self.feedback.enable()
            args = (self, context)

            self.state = 'STARTING'

            self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_callback, args, 'WINDOW', 'POST_PIXEL')
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to manage parts
# ------------------------------------------------------------------


class ARCHIPACK_OT_wall2_insert(Operator):
    bl_idname = "archipack.wall2_insert"
    bl_label = "Insert"
    bl_description = "Insert part"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    index = IntProperty(default=0)

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            d = ARCHIPACK_PT_wall2.params(o)
            if d is None:
                return {'CANCELLED'}
            d.insert_part(context, self.index)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_wall2_remove(Operator):
    bl_idname = "archipack.wall2_remove"
    bl_label = "Remove"
    bl_description = "Remove part"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    index = IntProperty(default=0)

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            d = ARCHIPACK_PT_wall2.params(o)
            if d is None:
                return {'CANCELLED'}
            d.remove_part(context, self.index)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to manipulate object
# ------------------------------------------------------------------


class ARCHIPACK_OT_wall2_manipulate(Operator):
    bl_idname = "archipack.wall2_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}

    def __del__(self):
        print("ARCHIPACK_OT_wall2_manipulate End")

    @classmethod
    def poll(self, context):
        return ARCHIPACK_PT_wall2.filter(context.active_object)

    def invoke(self, context, event):
        o = context.active_object
        o.data.archipack_wall2[0].manipulable_invoke(context)
        return {'FINISHED'}


    def execute(self, context):
        """
            For use in boolean ops
        """
        if ARCHIPACK_PT_wall2.filter(context.active_object):
            o = context.active_object
            d = o.data.archipack_wall2[0]
            g = d.get_generator()
            d.setup_childs(o, g)
            d.update_childs(context, o, g)
            d.update(context)
            o.select = True
            context.scene.objects.active = o
        return {'FINISHED'}


def register():
    bpy.utils.register_class(archipack_wall2_part)
    bpy.utils.register_class(archipack_wall2_child)
    bpy.utils.register_class(archipack_wall2)
    Mesh.archipack_wall2 = CollectionProperty(type=archipack_wall2)
    bpy.utils.register_class(ARCHIPACK_PT_wall2)
    bpy.utils.register_class(ARCHIPACK_OT_wall2)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_draw)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_insert)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_remove)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_manipulate)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_from_slab)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_throttle_update)


def unregister():
    bpy.utils.unregister_class(archipack_wall2_part)
    bpy.utils.unregister_class(archipack_wall2_child)
    bpy.utils.unregister_class(archipack_wall2)
    del Mesh.archipack_wall2
    bpy.utils.unregister_class(ARCHIPACK_PT_wall2)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_draw)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_insert)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_remove)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_manipulate)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_from_slab)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_throttle_update)
