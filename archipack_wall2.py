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
from mathutils.geometry import intersect_line_plane
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

    def straight_wall(self, last, length, da, wall_z, z, t):
        r = last.straight(length).rotate(da)
        return StraightWall(last, r.p, r.v, wall_z, z, t, self.flip)

    def curved_wall(self, last, a0, da, radius, wall_z, z, t):
        n = last.normal(1).rotate(a0).scale(radius)
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedWall(last, c, radius, a0, da, wall_z, z, t, self.flip)


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

    def add_part(self, type, center, radius, a0, da, length, wall_z, part_z, part_t, n_splits, flip):

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
                if da < 0:
                    c = Vector((radius, 0))
                else:
                    c = Vector((-radius, 0))
                s = CurvedWall(s, c, radius, a0, da, wall_z, z, t, flip)

        else:
            if type == 'S_WALL':
                s = s.straight_wall(s, length, a0, wall_z, z, t)
            elif type == 'C_WALL':
                s = s.curved_wall(s, a0, da, radius, wall_z, z, t)

        self.segs.append(s)
        self.last_type = type
        return manip_index

    def make_wall(self, step_angle, verts, faces):
        for i, wall in enumerate(self.segs):
            manipulators = self.parts[i].manipulators
            p0 = wall.lerp(0).to_3d()
            p1 = wall.lerp(1).to_3d()
            # angle from last to current segment
            if i > 0:
                v0 = self.segs[i - 1].straight(-1, 1).v.to_3d()
                v1 = wall.straight(1, 0).v.to_3d()
                manipulators[0].set_pts([p0, v0, v1])

            if type(wall).__name__ == "StraightWall":
                # segment length
                manipulators[1].type_key = 'SIZE'
                manipulators[1].prop1_name = "length"
                manipulators[1].set_pts([p0, p1, (1, 0, 0)])
            else:
                # segment radius + angle
                v0 = (wall.lerp(0) - wall.c).to_3d()
                v1 = (wall.lerp(1) - wall.c).to_3d()
                manipulators[1].type_key = 'ARC_ANGLE_RADIUS'
                manipulators[1].prop1_name = "da"
                manipulators[1].prop2_name = "radius"
                manipulators[1].set_pts([wall.c.to_3d(), v0, v1])

            # snap manipulator, dont change index !
            manipulators[2].set_pts([p0, p1, (1, 0, 0)])

            wall.param_t(step_angle)
            for j in range(wall.n_step + 1):
                wall.make_wall(j, verts, faces)

    def debug(self, verts):
        for wall in self.segs:
            for i in range(33):
                x, y = wall.lerp(i / 32)
                verts.append((x, y, 0))


def update(self, context):
    self.update(context)


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


class archipack_wall2_part(PropertyGroup):
    type = EnumProperty(
            items=(
                ('S_WALL', 'Straight', '', 0),
                ('C_WALL', 'Curved', '', 1)
                ),
            default='S_WALL',
            update=update_manipulators
            )
    length = FloatProperty(
            name="length",
            min=0.01,
            max=100.0,
            default=2.0,
            update=update
            )
    radius = FloatProperty(
            name="radius",
            min=0.5,
            max=100.0,
            default=0.7,
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
            max=1000,
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
        if child.data is not None:
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
            max=100.0,
            default=0.2,
            update=update
            )
    z = FloatProperty(
            name='height',
            min=0.1, max=10000,
            default=2.7, precision=2,
            description='height', update=update,
            )
    x_offset = FloatProperty(
            name="x offset",
            min=-1, max=1,
            default=-1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    radius = FloatProperty(
            name="radius",
            min=0.5,
            max=100.0,
            default=0.7,
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
    flip = BoolProperty(default=False, update=update)
    closed = BoolProperty(
            default=False,
            name="Close",
            update=update
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
        part_0 = self.parts[where]
        part_0.length /= 2
        part_0.da /= 2
        p = self.parts.add()
        s = p.manipulators.add()
        s.type_key = "ANGLE"
        s.prop1_name = "a0"
        s = p.manipulators.add()
        s.prop1_name = "length"
        s = p.manipulators.add()
        # s.type_key = 'SNAP_POINT'
        s.type_key = 'WALL_SNAP'
        s.prop1_name = str(where + 1)
        s.prop2_name = 'z'
        part_1 = self.parts[len(self.parts) - 1]
        part_1.type = part_0.type
        # part_1.z_left = part_0.z_left
        # part_1.z_right = part_0.z_right
        part_1.length = part_0.length
        part_1.da = part_0.da
        part_1.a0 = part_0.a0
        part_0.a0 = 0
        self.parts.move(len(self.parts) - 1, where)
        self.n_parts += 1
        self.auto_update = True
        self.update(context, manipulable_refresh=True)

    def add_part(self, context, length):
        self.manipulable_disable(context)
        self.auto_update = False
        p = self.parts.add()
        s = p.manipulators.add()
        s.type_key = "ANGLE"
        s.prop1_name = "a0"
        s = p.manipulators.add()
        s.prop1_name = "length"
        s = p.manipulators.add()
        s.type_key = 'WALL_SNAP'
        s.prop1_name = str(self.n_parts)
        s.prop2_name = 'z'
        p.length = length
        self.n_parts += 1
        self.auto_update = True
        self.update(context, manipulable_refresh=True)
        return p

    def remove_part(self, context, where):
        self.manipulable_disable(context)
        self.auto_update = False
        self.parts.remove(where)
        self.n_parts -= 1
        self.auto_update = True
        self.update(context, manipulable_refresh=True)

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

    def from_lines(self, lines):
        """
            TODO:
                look at fence implementation for this part
        """
        raise NotImplementedError
        
        for i, line in enumerate(lines):
            if type(line).__name__ == 'Line':
                self.parts[i].a0 = line.angle
                self.parts[i].length = line.length
            else:
                self.parts[i].radius = line.r
                self.parts[i].a0 = line.a0
                self.parts[i].da = line.da

    def get_generator(self):
        # print("get_generator")
        center = Vector((0, 0))
        g = WallGenerator(self.parts)
        for part in self.parts:
            g.add_part(part.type, center, part.radius, part.a0, part.da, part.length,
                self.z, part.z, part.t, part.n_splits, self.flip)
        return g

    def update_parts(self, o):
        # print("update_parts")
        # remove rows
        row_change = False
        for i in range(len(self.parts), self.n_parts, -1):
            row_change = True
            self.parts.remove(i - 1)

        # add rows
        for i in range(len(self.parts), self.n_parts):
            row_change = True
            p = self.parts.add()
            s = p.manipulators.add()
            s.type_key = "ANGLE"
            s.prop1_name = "a0"
            s = p.manipulators.add()
            s.type_key = "SIZE"
            s.prop1_name = "length"
            s = p.manipulators.add()
            # s.type_key = 'SNAP_POINT'
            s.type_key = 'WALL_SNAP'
            s.prop1_name = str(i)
            s.prop2_name = 'z'

        g = self.get_generator()

        if row_change:
            self.setup_childs(o, g)

        return g

    def update(self, context, manipulable_refresh=False):
        # print("update manipulable_refresh:%s" % (manipulable_refresh))

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

        g = self.update_parts(o)
        # print("make_wall")
        g.make_wall(self.step_angle, verts, faces)

        if self.closed:
            f = len(verts)
            faces.append((f - 2, 0, 1, f - 1))

        # print("buildmesh")
        bmed.buildmesh(context, o, verts, faces, matids=None, uvs=None, weld=True)

        # Width
        self.manipulators[0].set_pts([(0, 0, 0), g.segs[0].sized_normal(0, -self.width).v.to_3d(), (-1, 0, 0)])

        # Parts COUNTER
        self.manipulators[1].set_pts([g.segs[-1].lerp(1.1).to_3d(),
            g.segs[-1].lerp(1.1 + 0.5 / g.segs[-1].length).to_3d(), (-1, 0, 0)])

        # Height
        self.manipulators[2].set_pts([(0, 0, 0), (0, 0, self.z), (-1, 0, 0)], normal=g.segs[0].straight(1, 0).v.to_3d())

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
            if array[i][0] < array[begin][0]:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
            # param t on the wall
            elif array[i][0] == array[begin][0] and array[i][3] <= array[begin][3]:
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
        m.type_key = 'SIZE_LOC'
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
                    dir = -wall.normal(t).v.normalized()
                    if res and t > 0 and t < 1 and abs(d) < dmax:
                        wall_with_childs[wall_idx] = 1
                        m = self.childs_manipulators.add()
                        m.type_key = 'DUMB_SIZE'
                        # store z in wall space
                        relocate.append((child.name, wall_idx, (t * wall.length, d, (itM * child.location).z),
                            (dir - dir_y).length > 0.5))

        self.sort_child(relocate)
        for child in relocate:
            name, wall_idx, pos, flip = child
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

        w = self.width
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
                d.y = w
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
                        self.childs_manipulators[m_idx].set_pts([(p0.x, p0.y, 0), (p1.x, p1.y, 0), (0.5, 0, 0)])
                        m_idx += 1
                        x, y = 0.5 * d.x, 0.5 * d.y
                        if child.flip:
                            side = -1
                        else:
                            side = 1
                        # delta loc
                        child.manipulators[0].set_pts([(-x, side * -y, 0), (x, side * -y, 0), (side, 0, 0)])
                        # loc size
                        child.manipulators[1].set_pts([(-x, side * -y, 0), (x, side * -y, 0), (0.5 * side, 0, 0)])
                        p0 = wall.lerp(t + dt)
            p1 = wall.lerp(1)
            if wall_has_childs:
                # dub size after all childs
                self.childs_manipulators[m_idx].set_pts([(p0.x, p0.y, 0), (p1.x, p1.y, 0), (0.5, 0, 0)])
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
                        self.manip_stack.append(child.manipulators[0].setup(context, c, d))
                        # loc size
                        self.manip_stack.append(child.manipulators[1].setup(context, c, d))

    def manipulable_manipulate(self, context, event=None, manipulator=None):
        type_name = type(manipulator).__name__
        # print("manipulable_manipulate %s" % (type_name))
        if type_name in ['DeltaLocationManipulator', 'SizeLocationManipulator']:
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
            if not self.realtime:
                self.update_childs(context, o, g)
            # trigger wall's update by hand as those manipulations
            # dosen't affect wall data 
            # NOTE: update manipulators should be sufficient here
            self.update(context)

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
        d = self

        # setup childs manipulators
        self.manipulate_childs(context)

        for i, part in enumerate(d.parts):
            if i >= d.n_parts:
                break
            if i > 0:
                # start angle
                self.manip_stack.append(part.manipulators[0].setup(context, o, part))

            # length / radius + angle
            self.manip_stack.append(part.manipulators[1].setup(context, o, part))

            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))
            
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
        # select frame
        o.select = True
        context.scene.objects.active = o
        if self.auto_manipulate:
            bpy.ops.archipack.wall2_manipulate('INVOKE_DEFAULT')
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = bpy.context.scene.cursor_location
            o.select = True
            context.scene.objects.active = o
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
        p0 = sp.takeloc
        p1 = sp.placeloc
        # self.feedback.draw(context)
        if sp.delta.length == 0:
            return
        self.wall_part1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.wall_line1.set_pos([p0, p1, Vector((p1.x, p1.y, p1.z + z)), Vector((p0.x, p0.y, p0.z + z))])
        self.wall_part1.draw(context)
        self.wall_line1.draw(context)
        self.line.p = p0
        self.line.v = sp.delta
        self.label.set_pos(context, self.line.length, self.line.lerp(0.5), self.line.v, normal=Vector((0, 0, 1)))
        self.label.draw(context)
        self.line.draw(context)

    def sp_callback(self, context, event, state, sp):

        self.state = state

        if state == 'SUCCESS':
            old = context.active_object
            if self.o is None:
                bpy.ops.archipack.wall2(auto_manipulate=False)
                o = context.active_object
                o.location = sp.takeloc
                self.o = o
                d = o.data.archipack_wall2[0]
                part = d.parts[0]
                part.length = sp.delta.length
            else:
                o = self.o
                o.select = True
                context.scene.objects.active = o
                d = o.data.archipack_wall2[0]
                part = d.add_part(context, sp.delta.length)

            # print("self.o :%s" % o.name)
            rM = o.matrix_world.inverted().to_3x3()
            g = d.get_generator()
            w = g.segs[-1]
            dp = rM * sp.delta
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
        elif state == 'CANCEL':
            return

    def modal(self, context, event):

        context.area.tag_redraw()

        if event.type in {'LEFTMOUSE', 'RET', 'NUMPAD_ENTER', 'SPACE'}:

            self.feedback.instructions(context, "Draw a wall", "Drag to move, click to add a segment", [
                ('CTRL', 'Snap'),
                ('MMBTN', 'Constraint to axis'),
                ('X Y', 'Constraint to axis'),
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
                    takeloc = o.matrix_world * g.segs[-1].lerp(1).to_3d()
                    o.select = False
                else:
                    takeloc = self.mouse_to_plane(context, event)

                snap_point(takeloc, self.sp_draw, self.sp_callback, constraint_axis=(True, True, False))
            return {'RUNNING_MODAL'}

        if self.state == 'CANCEL' or (event.type in {'ESC', 'RIGHTMOUSE'} and
                event.value == 'RELEASE'):

            self.feedback.disable()
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')

            if self.o is not None:
                self.o.select = True
                context.scene.objects.active = self.o
                bpy.ops.archipack.wall2_manipulate('INVOKE_DEFAULT')

            return {'FINISHED'}

        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        if context.mode == "OBJECT":
            self.wall_part1 = GlPolygon((0.5, 0, 0, 0.2))
            self.wall_line1 = GlPolyline((0.5, 0, 0, 0.8))
            self.line = GlLine()
            self.label = GlText()
            self.feedback = FeedbackPanel()
            self.feedback.instructions(context, "Draw a wall", "Click to start and drag to move", [
                ('CTRL', 'Snap'),
                ('MMBTN', 'Constraint to axis'),
                ('X Y', 'Constraint to axis'),
                ('RIGHTCLICK or ESC', 'exit without change')
                ])
            self.feedback.enable()
            args = (self, context)
            self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_callback, args, 'WINDOW', 'POST_PIXEL')
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to create object
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
            prop = ARCHIPACK_PT_wall2.params(o)
            if prop is None:
                return {'CANCELLED'}
            prop.insert_part(context, self.index)
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
            prop = ARCHIPACK_PT_wall2.params(o)
            if prop is None:
                return {'CANCELLED'}
            prop.remove_part(context, self.index)
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

    @classmethod
    def poll(self, context):
        return ARCHIPACK_PT_wall2.filter(context.active_object)

    def modal(self, context, event):
        return self.d.manipulable_modal(context, event)

    def invoke(self, context, event):
        if context.space_data.type == 'VIEW_3D':
            o = context.active_object
            self.d = o.data.archipack_wall2[0]
            if self.d.manipulable_invoke(context):
                context.window_manager.modal_handler_add(self)
                return {'RUNNING_MODAL'}
            else:
                return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}

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
            d.manipulable_release(context)
            d.update(context)
            o.select = True
            context.scene.objects.active = o
        return {'FINISHED'}


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
bpy.utils.register_class(ARCHIPACK_OT_wall2_throttle_update)
