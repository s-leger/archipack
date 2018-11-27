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
import time
import bpy
import bmesh
from math import sin, cos, pi, atan2, tan
from mathutils import Vector, Matrix
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty, BoolProperty, IntProperty, StringProperty,
    FloatVectorProperty, CollectionProperty, EnumProperty
)
from .bmesh_utils import BmeshEdit as bmed
from .archipack_object import (
    ArchipackObject, ArchipackCreateTool,
    ArchipackDrawTool, ArchipackObjectsManager
    )
from .archipack_snap import snap_point
from .archipack_keymaps import Keymaps
from .archipack_polylines import Io
from .archipack_dimension import DimensionProvider
from .archipack_curveman import ArchipackUserDefinedPath
from .archipack_segments import StraightSegment, CurvedSegment, Generator, ArchipackSegment
from .archipack_throttle import throttle
from .archipack_manipulator import (
    Manipulable, archipack_manipulator,
    GlPolygon, GlPolyline,
    GlLine, GlText, FeedbackPanel
    )
# use multi wall support to generate floors / moldings
FLOORS_2D_SUPPORT_MULTI = True
from .archipack_polylines import CoordSys, Qtree
import logging
logger = logging.getLogger("archipack")


class Q_tree(Qtree):
    """
     A quadtree to minimize relocate intersections
    """
    def getbounds(self, seg, extend=0):
        x0, y0 = seg.p0.x, seg.p0.y
        x1, y1 = seg.p1.x, seg.p1.y
        return (min(x0, x1),
            min(y0, y1),
            max(x0, x1),
            max(y0, y1))

    def intersects(self, pt, extend):
        x, y = pt.x, pt.y
        bounds = (x - extend,
            y - extend,
            x + extend,
            y + extend)
        selection = list(self._intersect(bounds))
        count = len(selection)
        return count, sorted(selection)

    def insert(self, seg):
        idx = self.ngeoms
        self._geoms.append(seg)
        self._insert(idx, self.getbounds(seg[2]))


class Wall():
    def __init__(self, wall_z, z, t, flip, matid):
        self.z = z
        self.wall_z = wall_z
        self.t = t
        self.flip = flip
        self.matid = matid
        self.z_step = len(z)

    def set_offset(self, offset, last=None):
        """
            Offset line and compute intersection point
            between segments
        """
        self.line = self.make_offset(offset, last)

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

    def make_faces(self, f, faces):
        faces.append((f, f + 2, f + 3, f + 1))

    def p3d(self, verts, t, z_min=0):
        p = self.line.lerp(t)
        z = z_min + self.wall_z + self.get_z(t)
        verts.append((p.x, p.y, z_min))
        verts.append((p.x, p.y, z))

    def make_wall(self, make_faces, j, z_min, verts, faces, matids):
        t = self.t_step[j]
        f = len(verts)
        self.p3d(verts, t, z_min)
        if make_faces:
            self.make_faces(f - 2, faces)

    def get_coord(self, i, verts, z0):
        t = self.t_step[i]
        p = self.line.lerp(t)
        verts.append((p.x, p.y, z0))


class StraightWall(Wall, StraightSegment):
    def __init__(self, p, v, wall_z, z, t, flip, matid):
        StraightSegment.__init__(self, p, v)
        Wall.__init__(self, wall_z, z, t, flip, matid)

    def param_t(self, step_angle):
        self.t_step = self.t
        self.n_step = len(self.t) - 1


class CurvedWall(Wall, CurvedSegment):
    def __init__(self, c, radius, a0, da, wall_z, z, t, flip, matid):
        CurvedSegment.__init__(self, c, radius, a0, da)
        Wall.__init__(self, wall_z, z, t, flip, matid)

    def param_t(self, step_angle):
        t_step, n_step = self.steps_by_angle(step_angle)
        self.t_step = list(sorted([i * t_step for i in range(1, n_step)] + self.t))
        self.n_step = len(self.t_step) - 1


class WallGenerator(Generator):
    def __init__(self, d, o=None):
        Generator.__init__(self, d, o)
        self.last_type = 'NONE'
        self.faces_type = 'NONE'

    def add_part(self, part):

        # TODO:
        # refactor this part (height manipulators)
        d = self.d

        manip_index = []
        if len(self.segs) < 1:
            s = None
            z = [part.z[0]]
            manip_index.append(0)
        else:
            s = self.segs[-1]
            z = [s.z[-1]]

        t_cur = 0
        z_last = part.n_splits - 1
        t = [0]

        for i in range(part.n_splits):
            t_try = t[-1] + part.t[i]
            if t_try == t_cur:
                continue
            if t_try <= 1:
                t_cur = t_try
                t.append(t_cur)
                z.append(part.z[i])
                manip_index.append(i)
            else:
                z_last = i
                break

        if t_cur < 1:
            t.append(1)
            manip_index.append(z_last)
            z.append(part.z[z_last])

        # start a new wall
        if 'C_' in part.type:
            if s is None:
                c = self.location - (self.rot * (part.radius * Vector((cos(part.a0), sin(part.a0), 0)))).to_2d()
                s = CurvedWall(c, part.radius, part.a0, part.da, d.z, z, t, d.flip, part.material_index - 1)
            else:
                n = s.normal(1).rotate(part.a0).scale(part.radius)
                if part.da < 0:
                    n.v = -n.v
                a0 = n.angle
                c = n.p - n.v
                s = CurvedWall(c, part.radius, a0, part.da, d.z, z, t, d.flip, part.material_index - 1)
        else:
            if s is None:
                p = self.location
                v = (self.rot * Vector((part.length, 0, 0))).to_2d()
                s = StraightWall(p, v, d.z, z, t, d.flip, part.material_index - 1).rotate(part.a0)
            else:
                r = s.straight(part.length).rotate(part.a0)
                s = StraightWall(r.p, r.v, d.z, z, t, d.flip, part.material_index - 1)

        self.segs.append(s)
        self.last_type = part.type

        return manip_index

    def make_wall(self, verts, faces, matids):

        d = self.d

        _segs = self.segs
        nb_segs = self.n_parts
        make_faces = False
        for i, wall in enumerate(_segs):
            wall.param_t(d.step_angle)

            if i < nb_segs:
                for j in range(wall.n_step):
                    wall.make_wall(make_faces, j, -d.z_offset, verts, faces, matids)
                    matids.append(wall.matid)
                    make_faces = True

        # when wall is open, add last section of vertices
        if not d.closed:
            wall = _segs[-2]
            wall.make_wall(True, wall.n_step, -d.z_offset, verts, faces, matids)

        # swap manipulators so they always face outside
        side = 1
        if d.flip:
            side = -1

        self.locate_manipulators(side)

    def debug(self, verts):
        for wall in self.segs:
            for i in range(33):
                x, y = wall.lerp(i / 32)
                verts.append((x, y, 0))

    def make_surface(self, o, verts, height):
        bm = bmesh.new()
        for v in verts:
            bm.verts.new(v)
        bm.verts.ensure_lookup_table()
        for i in range(1, len(verts)):
            bm.edges.new((bm.verts[i - 1], bm.verts[i]))
        bm.edges.new((bm.verts[-1], bm.verts[0]))
        bm.edges.ensure_lookup_table()
        bmesh.ops.contextual_create(bm, geom=bm.edges)
        geom = bm.faces[:]
        bmesh.ops.solidify(bm, geom=geom, thickness=height)
        bm.to_mesh(o.data)
        bm.free()

    def get_coords(self, offset, z, d, verts):
        self.set_offset(offset)
        self.close(offset)
        _segs = self.segs
        nb_segs = self.n_parts
        f = len(verts)

        for i, wall in enumerate(_segs):
            wall.param_t(d.step_angle)
            if i < nb_segs:
                for j in range(wall.n_step):
                    wall.get_coord(j, verts, z)

        if not d.closed:
            wall = _segs[-2]
            wall.get_coord(wall.n_step, verts, z)

        # fix precision issues by extending ends
        # on open walls
        if d.extend % 2 == 1:
            seg = self.segs[0].line
            p = seg.lerp(-0.02 / seg.length)
            verts[f] = (p.x, p.y, z)

        if d.extend > 1:
            seg = self.segs[-2].line
            p = seg.lerp(1 + 0.02 / seg.length)
            verts[-1] = (p.x, p.y, z)

    def get_ends(self, offset, o, itM, dist):
        """
         Return points outside wall at both ends
         for neighboor analysis
         only work on open walls
        """
        self.set_offset(offset)
        seg = self.segs[0].line
        p0 = itM * o.matrix_world * seg.lerp(-dist / seg.length).to_3d()
        seg = self.segs[-2].line
        p1 = itM * o.matrix_world * seg.lerp(1 + dist / seg.length).to_3d()
        return p0, p1


def update(self, context):
    self.update(context)
    return None


def update_childs(self, context):
    self.update(context, update_childs=True, manipulable_refresh=True)
    return None


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)
    return None


def update_relocate(self, context):
    self.update(context, relocate_childs=True)
    return None


def update_mat(self, context):
    self.update(context, relocate_childs=False)
    return None


def get_width(self):
    return self.width


def set_width(self, value):
    self.last_width = self.width
    self.width = value
    return None


def update_width(self, context):
    o = context.active_object
    if len(self.reloc_childs) == 0:
        self.setup_childs(context, o)
    # delta is on inside part (negative distance)
    delta = 0.5 * (self.last_width - self.width)
    for child in self.reloc_childs:
        if child.parent_name0 == o.name:
            # d0 is distance from axis
            if child.d0 > 0:
                child.d0 -= delta
            else:
                child.d0 += delta
        if child.parent_name1 == o.name:
            # d1 is distance from axis
            if child.d1 > 0:
                child.d1 -= delta
            else:
                child.d1 += delta
    self.update(context, relocate_childs=True)
    return None


def update_t_part(self, context):
    """
        Make this wall a T child of parent wall
        orient child so y points inside wall and x follow wall segment
        set child a0 according
    """
    o = self.find_in_selection(context)
    if o is not None:

        # w is parent wall
        w = self.get_scene_object(context, self.t_part)
        wd = archipack_wall2.datablock(w)

        if wd is not None:
            og = self.get_generator()
            self.setup_childs(context, o)

            bpy.ops.object.select_all(action="DESELECT")

            # 5 cases here:
            # 1 No parents at all
            # 2 o has parent
            # 3 w has parent
            # 4 o and w share same parent allready
            # 5 o and w dosent share parent
            link_to_parent = False

            # when both walls do have a reference point, we may delete one of them
            to_delete = None

            # select childs and make parent reference point active
            if w.parent is None:
                # Either link to o.parent or create new parent
                link_to_parent = True
                if o.parent is None:
                    # create a reference point and make it active
                    x, y, z = w.bound_box[0]
                    context.scene.cursor_location = w.matrix_world * Vector((x, y, z))
                    context.scene.cursor_location.z += self.z_offset
                    # fix issue #9
                    self.select_object(context, o, True)
                    bpy.ops.archipack.reference_point()
                    self.select_object(context, o)
                else:
                    self.select_object(context, o.parent, True)
                self.select_object(context, w)
            else:
                # w has parent
                if o.parent is not w.parent:
                    link_to_parent = True
                    self.select_object(context, w.parent, True)
                    self.select_object(context, o)
                    if o.parent is not None:
                        # store o.parent to delete it
                        to_delete = o.parent
                        for c in o.parent.children:
                            if c is not o:
                                c.hide_select = False
                                self.select_object(context, c)

            parent = context.object

            dmax = 2 * wd.width

            wg = wd.get_generator()

            otM = o.matrix_world
            orM = Matrix([
                otM[0].to_2d(),
                otM[1].to_2d()
                ])

            wtM = w.matrix_world
            wrM = Matrix([
                wtM[0].to_2d(),
                wtM[1].to_2d()
                ])

            # dir in absolute world coordsys
            dir = orM * og.segs[0].straight(1, 0).v

            # pt in w coordsys
            pos = otM.translation
            pt = (wtM.inverted() * pos).to_2d()

            for wall_idx, wall in enumerate(wg.segs):
                res, dist, t = wall.point_sur_segment(pt)
                # outside is on the right side of the wall
                #  p1
                #  |-- x
                #  p0

                # NOTE:
                # rotation here is wrong when w has no parent while o has parent
                t_bound = wall.length / self.width
                if res and t > -t_bound and t < 1 + t_bound and abs(dist) < dmax:
                    x = wrM * wall.straight(1, t).v
                    y = wrM * wall.normal(t).v.normalized()
                    self.parts[0].a0 = dir.angle_signed(x)
                    self.select_object(context, o, True)
                    self.update(context)
                    o.matrix_world = Matrix([
                        [x.x, -y.x, 0, pos.x],
                        [x.y, -y.y, 0, pos.y],
                        [0, 0, 1, pos.z],
                        [0, 0, 0, 1]
                    ])
                    break
            self.select_object(context, parent, True)

            if link_to_parent and bpy.ops.archipack.parent_to_reference.poll():
                bpy.ops.archipack.parent_to_reference('INVOKE_DEFAULT')

            # update generator to take new rotation in account
            # use this to relocate windows on wall after reparenting
            # g = self.get_generator()
            self.relocate_childs(context, o)

            # hide holes from select
            for c in parent.children:
                if "archipack_hybridhole" in c:
                    c.hide_select = True

            # delete unneeded reference point
            if to_delete is not None:
                bpy.ops.object.select_all(action="DESELECT")
                self.select_object(context, to_delete, True)
                if bpy.ops.object.delete.poll():
                    bpy.ops.object.delete(use_global=False)

        elif self.t_part != "":
            self.t_part = ""

    self.restore_context(context)
    return None


def set_splits(self, value):
    if self.n_splits != value:
        auto_update = self.auto_update
        self.auto_update = False
        self._set_t(value)
        self.auto_update = auto_update
        self.n_splits = value
    return None


def get_splits(self):
    return self.n_splits


class archipack_wall2_part(ArchipackSegment, PropertyGroup):

    z = FloatVectorProperty(
            name="Height",
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
            name="Position",
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
            name="Splits",
            default=1,
            min=1,
            max=31,
            get=get_splits, set=set_splits
            )
    n_splits = IntProperty(
            name="Splits",
            default=1,
            min=1,
            max=31,
            update=update
            )
    material_index = IntProperty(
            name="Material index",
            min=1,
            max=10,
            default=1,
            update=update_mat
            )
    manipulators = CollectionProperty(type=archipack_manipulator)

    # ui related
    expand = BoolProperty(default=False)

    def update(self, context, manipulable_refresh=False, relocate_childs=True):

        # Reset change side
        self.change_side = 'RIGHT'

        if self.auto_update:
            idx, o, d = self.find_datablock_in_selection(context)
            if d is not None:
                d.update(context, manipulable_refresh, relocate_childs=relocate_childs)

    def _set_t(self, splits):
        t = 1 / splits
        for i in range(splits):
            self.t[i] = t

    def get_datablock(self, o):
        return archipack_wall2.datablock(o)

    def draw(self, context, layout, index, draw_type=True):

        row = layout.row(align=True)
        icon = "TRIA_RIGHT"
        if self.expand:
            icon = "TRIA_DOWN"

        row.prop(self, 'expand', icon=icon, icon_only=True, text="Seg {}".format(index + 1), emboss=True)
        row.prop(self, "type_ui", text="")

        if self.expand:
            self.draw_insert(context, layout, index)
            if self.type == 'C_SEG':
                layout.prop(self, "r_ui")
                layout.prop(self, "da_ui")
            else:
                layout.prop(self, "l_ui")
            layout.prop(self, "a_ui")
            layout.prop(self, "material_index")
            layout.prop(self, "splits")
            for split in range(self.n_splits):
                row = layout.row()
                row.prop(self, "z", text="alt", index=split)
                row.prop(self, "t", text="pos", index=split)


class archipack_wall2_relocate_child(PropertyGroup):
    """
     Store child points relationship to parent segments
     as distance from 2 nearest segments
    """
    child_name = StringProperty(
        description="Name of child object"
        )
    child_idx = IntProperty(
        default=-1,
        description="Index of linked child part (startpoint) when -1 set object location"
        )
    seg0 = IntProperty(
        description="Index of parent segment 0 the point is bound to"
        )
    seg1 = IntProperty(
        description="Index of parent segment 1 the point is bound to"
        )
    parent_name0 = StringProperty(default="")
    parent_name1 = StringProperty(default="")
    d0 = FloatProperty(
        description="Distance from parent segment 0 the point is bound to"
        )
    d1 = FloatProperty(
        description="Distance from parent segment 1 the point is bound to"
        )
    t = FloatProperty(
        description="Distance from start of closest segment normalized"
        )

    def filter_child(self, o):
        d = None
        k = None
        if o and o.data:
            od = o.data
            for key in od.keys():
                if "archipack_" in key and key[10:] in (
                        "window",
                        "door",
                        "custom",
                        "wall2",
                        "slab",
                        "fence",
                        "floor",
                        "molding",
                        "slab_cutter",
                        "floor_cutter",
                        "roof_cutter"
                        ):
                    try:
                        d = getattr(od, key)[0]
                        k = key
                    except:
                        pass
                    break
        return d, k

    def get_child(self, context):
        d = None
        # 2.8 ???
        o = context.scene.objects.get(self.child_name)
        d, k = self.filter_child(o)
        return o, d

    def get_parent0(self, context):
        d = None
        # 2.8 ???
        o = context.scene.objects.get(self.parent_name0)
        d, k = self.filter_child(o)
        return o, d

    def get_parent1(self, context):
        d = None
        # 2.8 ???
        o = context.scene.objects.get(self.parent_name1)
        d, k = self.filter_child(o)
        return o, d


class archipack_wall2_child(PropertyGroup):
    manipulators = CollectionProperty(type=archipack_manipulator)
    child_name = StringProperty()
    child_idx = IntProperty(
        default=-1,
        description="Index of linked child part (startpoint)"
        )
    wall_idx = IntProperty(
        description="Index of wall part"
        )
    pos = FloatVectorProperty(subtype='XYZ')
    flip = BoolProperty(default=False)
    flippable = BoolProperty(default=False)

    def get_child(self, context):
        d = None
        # 2.8 ???
        child = context.scene.objects.get(self.child_name)
        if child is not None and child.data is not None:
            cd = child.data
            if 'archipack_window' in cd:
                d = cd.archipack_window[0]
            elif 'archipack_door' in cd:
                d = cd.archipack_door[0]
            elif 'archipack_custom' in cd:
                d = cd.archipack_custom[0]
        return child, d


def set_base_line(self, value):
    if value == 0:
        self.x_offset = -1
    elif value == 1:
        self.x_offset = 1
    else:
        self.x_offset = 0
    return None


def get_base_line(self):
    if self.x_offset == 0:
        return 2
    elif self.x_offset == 1:
        return 1
    else:
        return 0


class archipack_wall2_finish(PropertyGroup):
    location_index = IntProperty(
            name="Location",
            description="Add finish over segments with this material index",
            min=1,
            max=11,
            default=1,
            update=update
            )
    material_index = IntProperty(
            name="Material",
            description="Material index of finish",
            min=1,
            default=1,
            update=update
            )
    altitude = FloatProperty(
            name="Altitude",
            description="Altitude from bottom of wall",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    z = FloatProperty(
            name="Height",
            description="Height of finish",
            min=0.01,
            default=20,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    rotation = FloatProperty(
            name="Rotation",
            description="Rotate pattern",
            default=0,
            min=-pi/4,
            max=pi/4,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    pattern_type = EnumProperty(
        name="Pattern",
        items=(
            ('V_BOARD', "Vertical board", "Vertical board"),
            ('H_BOARD', "Horizontal board", "Horizontal board"),
            ('TILES', "Tiles", "Ceramic tiles"),
            ('LOGGING', "Logs", "Round logs"),
            ('USER_DEFINED', "User defined", "User defined object as finishing")
            ),
        default='H_BOARD',
        update=update
        )
    side = EnumProperty(
        name="Side",
        items=(
            ('INSIDE', "Inside", "Inside"),
            ('OUTSIDE', "Outside", "Outside")
            ),
        default='OUTSIDE',
        update=update
        )
    tile_x = FloatProperty(
            name="Width",
            description="Tile width",
            min=0.01,
            default=0.3,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_y = FloatProperty(
            name="Height",
            description="Tile height",
            default=0.6,
            min=0.01,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_z = FloatProperty(
            name="Thickness",
            description="Tile thickness",
            default=0.005,
            min=0.001,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_space = FloatProperty(
            name="Mortar",
            description="Mortar width",
            default=0.002,
            min=0.0001,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    board_x = FloatProperty(
            name="Width",
            description="Pattern width",
            min=0.01,
            default=0.15,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    board_z = FloatProperty(
            name="Thickness",
            description="Pattern thickness",
            default=0.02,
            min=0.001,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    expand = BoolProperty(
        default=False,
        options={'SKIP_SAVE'}
        )
    user_pattern = StringProperty(
        default="",
        update=update
        )

    def draw(self, context, layout, index):
        row = layout.row(align=True)
        icon = "TRIA_RIGHT"
        if self.expand:
            icon = "TRIA_DOWN"
        row.prop(self, 'expand', icon=icon, text="{}".format(index + 1), icon_only=True, emboss=True)
        row.prop(self, 'side', text="")
        row.prop(self, 'location_index')
        row.operator("archipack.wall2_remove_finish", icon="ZOOMOUT", text="").index = index

        if self.expand:
            layout.prop(self, 'altitude')
            layout.prop(self, 'z')
            # rotation still dosent work
            # layout.prop(self, 'rotation')
            row = layout.row(align=True)
            row.prop(self, 'pattern_type', text="")
            row.prop(self, 'material_index')
            if self.pattern_type == 'TILES':
                layout.prop(self, 'tile_x')
                layout.prop(self, 'tile_y')
                layout.prop(self, 'tile_z')
                layout.prop(self, 'tile_space')

            elif 'BOARD' in self.pattern_type or self.pattern_type == 'LOGGING':
                layout.prop(self, 'board_x')

            if 'BOARD' in self.pattern_type:
                layout.prop(self, 'board_z')

            if self.pattern_type == 'USER_DEFINED':
                row = layout.row()
                row.prop_search(self, "user_pattern", context.scene, "objects", text="")

    def find_datablock_in_selection(self, context):
        """
         Return selected object this instance belongs to
         provide support for "copy to selected"
        """
        selected = context.selected_objects[:]
        for o in selected:
            d = archipack_wall2.datablock(o)
            if d:
                for idx, part in enumerate(d.finish):
                    if part == self:
                        return o, d
        return None, None

    def update(self, context, manipulable_refresh=False, relocate_childs=False):

        o, d = self.find_datablock_in_selection(context)
        if d is not None:
            d.update(context, manipulable_refresh, relocate_childs=relocate_childs)


class archipack_wall2(ArchipackObject, ArchipackUserDefinedPath, Manipulable, DimensionProvider, PropertyGroup):
    parts = CollectionProperty(type=archipack_wall2_part)
    step_angle = FloatProperty(
            description="Curved parts segmentation",
            name="Step angle",
            min=1 / 180 * pi,
            max=pi,
            default=6 / 180 * pi,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )

    last_width = FloatProperty(
            options={'SKIP_SAVE'}
            )
    width_ui = FloatProperty(
            options={'SKIP_SAVE'},
            name="Width",
            min=0.01,
            default=0.2,
            unit='LENGTH', subtype='DISTANCE',
            get=get_width,
            set=set_width
            )
    width = FloatProperty(
            name="Width",
            min=0.01,
            default=0.2,
            unit='LENGTH', subtype='DISTANCE',
            update=update_width
            )
    z = FloatProperty(
            name='Height',
            min=0.1,
            default=2.7, precision=5,
            unit='LENGTH', subtype='DISTANCE',
            description='height', update=update,
            )
    base_line = EnumProperty(
            name="Base",
            description="Wall part on base line",
            items=(
                ('OUTSIDE', "Outside", "Outside"),
                ('INSIDE', "Inside", "Inside"),
                ('AXIS', "Axis", "Axis")),
            default='OUTSIDE',
            set=set_base_line,
            get=get_base_line,
            update=update_relocate
            )
    offset = FloatProperty(
            name="Offset",
            description="Lateral offset of wall",
            unit='LENGTH', subtype='DISTANCE',
            default=0, precision=5, step=1,
            update=update_relocate
            )
    x_offset = FloatProperty(
            name="Offset",
            unit='LENGTH', subtype='DISTANCE',
            min=-1, max=1,
            default=-1
            )
    z_offset = FloatProperty(
            name="Floor thickness",
            unit='LENGTH', subtype='DISTANCE',
            description='Thickness of floors',
            min=0,
            default=0.1, precision=5, step=1,
            update=update
            )
    radius = FloatProperty(
            name="Radius",
            min=0.5,
            default=0.7,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    da = FloatProperty(
            name="Angle",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    flip = BoolProperty(
            name="Flip",
            description="Flip inside and outside",
            default=False,
            update=update_relocate
            )

    closed = BoolProperty(
            default=False,
            name="Close",
            update=update_manipulators
            )
    always_closed = BoolProperty(
            default=False,
            name="Always closed geometry",
            description="Flag indicate whenever geometry path is always closed",
            options={'SKIP_SAVE'}
            )

    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )
    auto_synch = BoolProperty(
            default=True,
            name="Auto synchro",
            description="Synchronize objects"
            )
    dimensions = BoolProperty(
            default=False,
            name="Dimensions",
            description="Buid static dimensions",
            update=update_childs
            )
    # dumb manipulators to show sizes between childs
    childs_manipulators = CollectionProperty(type=archipack_manipulator)
    # store to manipulate windows and doors
    childs = CollectionProperty(type=archipack_wall2_child)
    # store childs with parts points locations to relocate
    reloc_childs = CollectionProperty(
        options={'SKIP_SAVE'},
        type=archipack_wall2_relocate_child
        )
    finish = CollectionProperty(type=archipack_wall2_finish)
    finish_expand = BoolProperty(
            name="Finishings",
            description="Display finishings options",
            default=False,
            options={'SKIP_SAVE'}
            )
    finish_enable = BoolProperty(
            name="Finishings",
            description="Enable finishings",
            default=True,
            update=update
            )
    t_part = StringProperty(
            name="Parent wall",
            description="This part will follow parent when set",
            default="",
            update=update_t_part
            )
    fit_roof = BoolProperty(
            name="Auto fit roof",
            description="Automatically fit surrounding roof",
            default=False,
            update=update
            )
    extend = IntProperty(
            description="Flag to extend wall 2d to prevent precision issues",
            default=0,
            options={'SKIP_SAVE'}
            )

    @property
    def side(self):
        side = 1
        if self.flip:
            side = -side
        return side

    @property
    def left_offset(self):
        """
         Overall offset of wall on left side of base line
         offset + [-1 | 0] * width
        """
        return self.offset + 0.5 * (self.x_offset * self.side - 1) * self.width

    @property
    def axis_offset(self):
        """
         Offset of axis line
         offset + [-1 | 0] * width
        """
        return self.offset + 0.5 * self.side * self.x_offset * self.width

    def before_insert_part(self, context, o, g):
        return

    def after_insert_part(self, context, o, where, distance):
        # get child location on split seg and
        # update so they are in where + 1 when needed
        self.setup_childs(context, o)

    def before_remove_part(self, context, o, g):
        self.setup_childs(context, o)

    def after_remove_part(self, context, o, where, distance):
        """
         Rebuild childs index and location
         When removing a part
        """
        # get child location on removed seg and
        # update so they are in where - 1
        for c in self.childs:
            # relocate childs of removed segment
            # on removed segment - 1
            if c.wall_idx == where:
                c.pos.x += distance

            # update child part index
            if c.wall_idx >= where:
                c.wall_idx -= 1

        for c in self.reloc_childs:
            if c.wall_id0 >= where:
                c.wall_id0 -= 1
            if c.wall_id1 >= where:
                c.wall_id1 -= 1

    def get_generator(self, o=None, axis=False):
        """
         o: Object or Matrix, make generator in coordsys
         axis: Generate offset at wall axis
        """
        # print("get_generator")
        g = WallGenerator(self, o)
        for part in self.parts:
            g.add_part(part)
        if axis:
            offset = self.axis_offset
        else:
            # left part of wall
            offset = self.left_offset
        g.set_offset(offset)
        g.close(offset)
        return g

    def setup_manipulators(self):

        if len(self.manipulators) == 0:
            # make manipulators selectable
            s = self.manipulators.add()
            s.prop1_name = "width"

            s = self.manipulators.add()
            s.prop1_name = "n_parts"
            s.type_key = 'COUNTER'

            s = self.manipulators.add()
            s.prop1_name = "z"
            s.normal = (0, 1, 0)

        if self.t_part != "" and len(self.manipulators) < 4:
            s = self.manipulators.add()
            s.prop1_name = "x"
            s.type_key = 'DELTA_LOC'

        self.setup_parts_manipulators('z')

    def from_spline(self, context, wM, resolution, spline):

        o = self.find_in_selection(context)

        if o is None:
            return

        pts = self.coords_from_spline(
            spline,
            wM,
            resolution,
            ccw=True
            )
        if len(pts) < 2:
            return

        # translation of target object
        o.matrix_world = Matrix.Translation(pts[0].copy())
        self.auto_update = False
        self.closed = spline.use_cyclic_u
        if self.closed:
            pts.pop()
        self.from_points(pts)
        self.auto_update = True

    def after_reverse(self, context, o):
        self.setup_childs(context, o)

        # flip does trigger relocate and keep childs orientation
        self.offset = -self.offset
        self.flip = not self.flip
        self.auto_update = True
        self.update(context, manipulable_refresh=True)

    def apply_modifier(self, o, modif):
        """
          move modifier on top of stack and apply
        """
        index = o.modifiers.find(modif.name)
        for i in range(index):
            bpy.ops.object.modifier_move_up(modifier=modif.name)
        bpy.ops.object.modifier_apply(apply_as='DATA', modifier=modif.name)

    def update(self,
            context,
            manipulable_refresh=False,
            update_childs=False,
            relocate_childs=False
            ):

        tim = time.time()
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        if manipulable_refresh:
            # prevent crash by removing all manipulators refs to datablock before changes
            self.manipulable_disable(context)

        roof = None
        rd = None
        if self.fit_roof:
            roof, rd = self.find_roof(o)
            if roof is not None:
                rg = rd.get_generator(roof)
                rg.make_roof(context)
                rg.make_wall_fit(
                    context,
                    roof,
                    o,
                    False,
                    False,
                    True)

        verts = []
        faces = []
        matids = []
        changed = self.update_parts()
        g = self.get_generator(axis=False)
        logger.debug("Wall2.update() 1 generator :%.2f seconds", time.time() - tim)

        if changed or update_childs:
            self.setup_childs(context, o)

        g.make_wall(verts, faces, matids)
        logger.debug("Wall2.update() 2 make_wall :%.2f seconds", time.time() - tim)

        if self.closed:
            f = len(verts)
            faces.append((f - 2, 0, 1, f - 1))

        matoffset = 1
        matinside = 0
        rim = 21

        if self.flip:
            matids = [2 * matid + 1 for matid in matids]
            matoffset = -matoffset
        else:
            matids = [2 * matid + matinside for matid in matids]

        bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=None, weld=False)
        logger.debug("Wall2.update() 3 buildmesh :%.2f seconds", time.time() - tim)

        side = self.side

        offset = self.offset + self.x_offset * side * 0.5 * self.width

        # Width
        self.manipulators[0].set_pts([
            g.segs[0].sized_normal(0, offset + 0.5 * side * self.width).v.to_3d(),
            g.segs[0].sized_normal(0, offset - 0.5 * side * self.width).v.to_3d(),
            (-0.5 * side, 0, 0)
            ])

        # Parts COUNTER
        self.manipulators[1].set_pts([g.segs[-2].lerp(1.1).to_3d(),
            g.segs[-2].lerp(1.1 + 0.5 / g.segs[-2].length).to_3d(),
            (-side, 0, 0)
            ])

        # Height
        self.manipulators[2].set_pts([
            Vector((0, 0, -self.z_offset)),
            Vector((0, 0, self.z - self.z_offset)),
            (-1, 0, 0)
            ], normal=g.segs[0].straight(side, 0).v.to_3d())

        if self.t_part != "":
            t = 0.3 / g.segs[0].length
            self.manipulators[3].set_pts([
                g.segs[0].sized_normal(t, 0.1).p1.to_3d(),
                g.segs[0].sized_normal(t, -0.1).p1.to_3d(),
                (1, 0, 0)
                ])

        logger.debug("Wall2.update() 4 manipulators :%.2f seconds", time.time() - tim)

        # if self.realtime:
        # update child location and size
        if relocate_childs and self.auto_synch:
            self.relocate_childs(context, o)

        # store gl points
        self.update_childs(context, o, g)

        logger.debug("Wall2.relocate() 5 relocate :%.2f seconds", time.time() - tim)

        # Add / remove providers to wall dimensions
        self.update_dimension(context, o, g)
        logger.debug("Wall2.update() 6 dimensions :%.2f seconds", time.time() - tim)

        # Update all dimensions
        if relocate_childs or update_childs:
            self.update_dimensions(context, o)

        logger.debug("Wall2.update() 7 dimensions :%.2f seconds", time.time() - tim)

        if self.fit_roof and roof is not None:
            # assign a 'top' vertex group
            vgroup = o.vertex_groups.get('Top')
            if vgroup is None:
                vgroup = o.vertex_groups.new()
                vgroup.name = 'Top'

            for i, v in enumerate(o.data.vertices):
                if i % 2 == 1:
                    vgroup.add([v.index], 1, 'REPLACE')
                else:
                    vgroup.remove([v.index])

        modif = o.modifiers.get('Wall')
        if modif is None:
            modif = o.modifiers.new('Wall', 'SOLIDIFY')
            modif.use_quality_normals = True
            modif.use_even_offset = True
            modif.offset = 1

        modif.material_offset_rim = rim
        modif.material_offset = matoffset
        modif.thickness = self.width

        logger.debug("Wall2.relocate() 8 modifier :%.2f seconds", time.time() - tim)
        # self.x_offset

        self.apply_modifier(o, modif)

        self.fit_roof_modifier(o, roof, rd, apply=True)

        if len(self.finish) > 0 and self.finish_enable:

            throttle.add(context, o, self, 1.0)

            if not throttle.is_active(o.name):

                # finish
                bm = bmed._start(context, o)
                bm.verts.index_update()
                bm.faces.index_update()

                patterns = []
                for finish in self.finish:
                    self.build_finish(context, o, finish, bm, patterns)

                bpy.ops.object.mode_set(mode='OBJECT')
                bm.free()

                bmed.bmesh_join(context, o, patterns, normal_update=True)

            """
            bm = bmed._start(context, o)
            bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.001)
            bmed._end(bm, o)
            """

        if manipulable_refresh:
            self.manipulable_refresh = True

        self.restore_context(context)

    # Finish
    def vertical_boards(self, profil, sx, sy, minx, maxx, miny, maxy, offset, profil_z, verts, faces):
        n_y = 1 + int(sx / profil_z)
        # start at 0,0
        verts.append(Vector((minx, miny, offset)))
        verts.append(Vector((minx, maxy, offset)))
        faces.append((0, 1, 3, 2))
        f = 2
        for i in range(n_y):
            x0 = minx + i * profil_z
            for co in profil:
                verts.append(Vector((x0 + co.z, miny, co.y)))
                verts.append(Vector((x0 + co.z, maxy, co.y)))
                faces.append((f, f + 1, f + 3, f + 2))
                f += 2

        x0 = minx + n_y * profil_z
        verts.append(Vector((x0, miny, offset)))
        verts.append(Vector((x0, maxy, offset)))
        faces.append((f, f + 1, 1, 0))

    def horizontal_boards(self, profil, sx, sy, minx, maxx, miny, maxy, offset, profil_z, verts, faces):
        n_z = 1 + int(sy / profil_z)
        # start at 0,0
        verts.append(Vector((minx, miny, offset)))
        verts.append(Vector((maxx, miny, offset)))
        faces.append((0, 2, 3, 1))
        f = 2
        for i in range(n_z):
            z0 = miny + i * profil_z
            for co in profil:
                z = min(maxy, z0 + co.z)
                verts.append(Vector((minx, z, co.y)))
                verts.append(Vector((maxx, z, co.y)))
                faces.append((f, f + 2, f + 3, f + 1))
                f += 2

        z0 = min(maxy, miny + n_z * profil_z)
        verts.append(Vector((minx, z0, offset)))
        verts.append(Vector((maxx, z0, offset)))
        faces.append((f, 0, 1, f + 1))

    def ceiling(self, sx, sy, minx, maxx, miny, maxy, offset, finish, verts, faces, matids):
        n_y = 1 + int(sy / finish.tile_y)
        n_x = 1 + int(sx / finish.tile_x)

        # start at 0,0
        #
        #
        #   x_x_______x_x   x_______x x
        #   x x_______x x   x_______x x
        #   | |       | |   |       | |
        #   x x_______x x   x_______x x
        #   x_x_______x_x
        #     x0     x1 x2
        #

        # half join
        h = 0.5 * finish.tile_space
        # ceiling
        x1 = finish.tile_x - 2 * h
        x2 = finish.tile_x - h
        y = 0
        z0 = finish.tile_z + offset

        y = miny
        dx = finish.tile_x
        dy = [h, finish.tile_y - h, finish.tile_y]

        verts.append(Vector((minx, y, z0)))
        for j in range(n_x):
            x0 = minx + h + j * dx
            verts.append(Vector((x0, y, z0)))
            verts.append(Vector((x0 + x1, y, z0)))
            verts.append(Vector((x0 + x2, y, z0)))

        f = 0
        fx = 1 + n_x * 3
        for i in range(n_y):
            y0 = miny + i * finish.tile_y
            for k in range(3):
                y = y0 + dy[k]
                # 1st col
                mat = 1
                if k == 1:
                    mat = 0
                verts.append(Vector((minx, y, z0)))
                if k > 0 and y > maxy:
                    if k == 1:
                        y = maxy - h
                    else:
                        y = maxy
                for j in range(n_x):
                    x0 = minx + h + j * dx
                    verts.append(Vector((x0, y, z0)))
                    verts.append(Vector((x0 + x1, y, z0)))
                    verts.append(Vector((x0 + x2, y, z0)))
                    # faces.append((f, f + 1, f + fx + 1, f + fx))
                    faces.append((f, f + fx, f + fx + 1, f + 1))
                    f += 1
                    faces.append((f, f + fx, f + fx + 1, f + 1))
                    f += 1
                    faces.append((f, f + fx, f + fx + 1, f + 1))
                    f += 1
                    matids.extend([1, mat, 1])
                f += 1

    def make_user_pattern(self,
                        sx, sy,
                        minx, maxx, miny, maxy, offset,
                        finish,
                        user_x, user_y, user_verts, user_faces, user_matids, user_uvs,
                        verts, faces, matids):

        n_y = 2 + int(sy / user_y)
        n_x = 2 + int(sx / user_x)

        for nx in range(n_x):
            x = minx + nx * user_x
            for ny in range(n_y):
                y = miny + ny * user_y
                tM = Matrix.Translation(Vector((x, y, offset)))
                f = len(verts)
                verts.extend([tM * p for p in user_verts])
                faces.extend([tuple([i + f for i in p]) for p in user_faces])
                matids.extend(user_matids)

    def build_finish(self, context, o, finish, bm, patterns):
        # offset from wall twice of remove double to prevent merge of faces
        offset = 0.002
        if finish.pattern_type in {'H_BOARD', 'V_BOARD'}:
            z = finish.board_z
            x = finish.board_x
            profil = [Vector(p) for p in [
                    (0, z + offset, 0), (0, z + offset, 0.7 * x),
                    (0, 0.6 * z + offset, 0.8 * x), (0, 0.6 * z + offset, x)
                    ]]

        elif finish.pattern_type == 'LOGGING':
            z = finish.board_x
            x = finish.board_x
            r = 0.5 * z
            a0 = -pi / 2
            segs = 16
            a = pi / segs
            profil = [Vector((0, r * cos(a0 + a * i) + 2 * offset, r + r * sin(a0 + a * i))) for i in range(segs)]

        elif finish.pattern_type == 'TILES':
            z = finish.tile_z

        elif finish.pattern_type == 'USER_DEFINED':
            pattern_obj = self.get_scene_object(context, finish.user_pattern)
            if pattern_obj is None:
                return
            user_x, user_y, z = pattern_obj.dimensions
            d = pattern_obj.data
            user_verts = [Vector((-v.co.x, v.co.y, v.co.z)) for v in d.vertices]
            user_faces = [tuple(p.vertices[:]) for p in d.polygons]
            if len(d.uv_layers) > 0:
                uv_layer = d.uv_layers[0].data
                user_uvs = [[uv_layer[v].uv for v in p.loop_indices] for p in d.polygons]
            else:
                user_uvs = [[user_verts[v].to_2d() for v in f] for f in user_faces]
            user_idmat = [p.material_index for p in d.polygons]

        z_bot = finish.altitude
        z_top = z_bot + finish.z

        # target segments material index
        location_index = 2 * (finish.location_index - 1)
        material_index = 2 * (finish.material_index - 1)
        if finish.side == 'INSIDE':
            location_index += 1
            material_index += 1

        # lateral profile in yz axis (last will be next first one)
        # so this is profile -last vertex
        z_axis = Vector((0, 0, 1))

        for face in bm.faces:

            if face.material_index == location_index and abs(face.normal.z) < 0.1:

                face_normal = face.normal
                x_face = face_normal.cross(z_axis)
                y_face = x_face.cross(face_normal)
                origin = face.calc_center_bounds()

                # Matrix of a plane at object's 0,0,0 and with z matching face normal
                # Prerotate to allow pattern rotations
                face_rot = Matrix([
                    [x_face.x, y_face.x, face_normal.x, 0],
                    [x_face.y, y_face.y, face_normal.y, 0],
                    [x_face.z, y_face.z, face_normal.z, 0],
                    [0, 0, 0, 1]
                    ]) * Matrix.Rotation(finish.rotation, 4, z_axis)

                # z of face center in face_rot coordsys
                delta_z = (face_rot.inverted() * origin).z

                # Pretranslate matrix along z to face
                face_mat = face_rot * Matrix.Translation(
                        Vector((0, 0, delta_z))
                        )

                face_inverse = face_mat.inverted()

                # loop in face coordsys
                bounds = [face_inverse * loop.vert.co for loop in face.loops]

                # in face_inverse space
                cx = [co.x for co in bounds]
                cy = [co.y for co in bounds]

                minx = min(cx)
                maxx = max(cx)
                miny = z_bot
                maxy = min(z_top, max(cy))

                print(minx, maxx, miny, maxy)

                # extend pattern for depth on edges (arbitrary * 5)
                dx = z * 5
                minx -= dx
                maxx += dx

                # make minx start at multiple of profile
                """
                dx = minx % profil_z
                if minx > 0:
                    minx -= dx
                else:
                    minx += dx - profil_z
                """
                sx = maxx - minx

                # enlarge for rotation
                dy = 2 * abs(tan(finish.rotation)) * sx
                miny -= dy
                maxy += dy

                # make miny start at multiple of profile
                """
                dy = miny % profil_z
                if miny > 0:
                    miny -= dy
                else:
                    miny += dy - profil_z
                """
                sy = maxy - miny

                if sx <= 0 or sy <= 0:
                    continue

                # build our pattern
                pattern = bmesh.new()

                # fill in pattern
                verts = []
                faces = []
                matids = []
                _verts = pattern.verts
                _faces = pattern.faces
                _new_vert = _verts.new
                _new_face = _faces.new

                if finish.pattern_type == 'H_BOARD':
                    self.horizontal_boards(profil, sx, sy, minx, maxx, miny, maxy, offset, finish.board_x, verts, faces)
                elif finish.pattern_type == 'V_BOARD':
                    self.vertical_boards(profil, sx, sy, minx, maxx, miny, maxy, offset, finish.board_x, verts, faces)
                elif finish.pattern_type == 'LOGGING':
                    self.horizontal_boards(profil, sx, sy, minx, maxx, miny, maxy, offset, finish.board_x, verts, faces)
                elif finish.pattern_type == 'TILES':
                    self.ceiling(sx, sy, minx, maxx, miny, maxy, offset, finish, verts, faces, matids)
                elif finish.pattern_type == 'USER_DEFINED':
                    self.make_user_pattern(
                        sx, sy,
                        minx, maxx, miny, maxy, offset,
                        finish,
                        user_x, user_y, user_verts, user_faces, user_idmat, user_uvs,
                        verts, faces, matids)

                for v in verts:
                    _new_vert(v)

                _verts.ensure_lookup_table()
                _verts.index_update()

                for f in faces:
                    _new_face([_verts[i] for i in f])

                _faces.ensure_lookup_table()
                _faces.index_update()

                if finish.pattern_type == 'TILES':
                    pattern.normal_update()
                    # 20 is for mortar
                    mats = [material_index, 20]
                    for i, mat in enumerate(matids):
                        _faces[i].material_index = mats[mat]

                    # inset
                    inset = 0.001
                    geom = [f for f in _faces if f.material_index == material_index]

                    if len(geom) > 0:
                        bmesh.ops.inset_individual(
                            pattern,
                            faces=geom,
                            thickness=inset,
                            depth=-inset,
                            use_even_offset=False,
                            use_interpolate=False,
                            use_relative_offset=False
                            )

                        _faces.ensure_lookup_table()

                    # extrude boundarys
                    geom = [ed for ed in pattern.edges[:] if ed.is_boundary]
                    if len(geom) > 0:
                        ret = bmesh.ops.extrude_edge_only(
                            pattern,
                            edges=geom)

                        geom = [v for v in ret["geom"] if isinstance(v, bmesh.types.BMVert)]
                        del ret
                        bmesh.ops.translate(
                            pattern,
                            verts=geom,
                            vec=Vector((0, 0, offset - finish.tile_z))
                            )

                elif finish.pattern_type == 'USER_DEFINED':

                    pattern.normal_update()

                    for i, mat in enumerate(matids):
                        _faces[i].material_index = mat

                    # weld
                    geom = [v for v in pattern.verts[:] if v.is_boundary]
                    if len(geom) > 0:
                        bmesh.ops.remove_doubles(
                            pattern,
                            verts=geom,
                            dist=0.0001
                            )

                    # extrude boundarys
                    geom = [ed for ed in pattern.edges[:] if ed.is_boundary]
                    if len(geom) > 0:
                        ret = bmesh.ops.extrude_edge_only(
                            pattern,
                            edges=geom)

                        geom = [v for v in ret["geom"] if isinstance(v, bmesh.types.BMVert)]
                        del ret
                        bmesh.ops.translate(
                            pattern,
                            verts=geom,
                            vec=Vector((0, 0, -0.5 * offset))
                            )

                else:
                    for i in range(len(_faces)):
                        _faces[i].material_index = material_index

                geom = [ed for ed in pattern.edges[:] if ed.is_boundary]
                if len(geom) > 0:
                    bmesh.ops.holes_fill(
                            pattern,
                            edges=geom,
                            sides=len(geom)
                            )

                # move our pattern in object's space
                pattern.transform(face_mat)

                # verts index in face
                f_index = [v.index for v in face.verts]

                # Slice using linked faces
                for edge in face.edges:

                    for link_face in edge.link_faces:
                        if link_face is not face:
                            slice_co = edge.verts[0].co
                            if abs(link_face.normal.z) < 0.1:
                                y = edge.verts[1].co - slice_co
                                # edge is not garantee to be in the right order
                                # ensure edge is in face order
                                last = f_index[-1]
                                for idx in f_index:
                                    if edge.verts[0].index == idx and edge.verts[1].index == last:
                                        break
                                    if edge.verts[1].index == idx and edge.verts[0].index == last:
                                        y = -y
                                        break
                                    last = idx
                                # intersection of faces using same material index
                                if link_face.material_index == location_index:
                                    # slice plane use both normals
                                    # perpendicular to edge, ensure always looking inside
                                    slice_no = (face_normal + link_face.normal).normalized().cross(y)
                                    slice_co -= 0.00005 * slice_no
                                else:
                                    # faces dosen't share material index
                                    slice_no = face_normal.cross(y)
                            else:
                                # either a top or bottom face
                                # slice plane use this face normal
                                slice_no = link_face.normal
                                # slice bottom part when pattern is under wall bottom
                                if link_face.normal.z < -0.5:
                                    slice_co = Vector((0, 0, finish.altitude))

                            # slice pattern using co and no
                            self.slice_pattern(pattern, slice_co, slice_no)

                # slice top
                slice_co = Vector((0, 0, finish.altitude + finish.z))
                slice_no = Vector((0, 0, 1))
                self.slice_pattern(pattern, slice_co, slice_no)

                # remove outside parts of loop
                patterns.append(pattern)

    def slice_pattern(self, pattern, slice_co, slice_no):
        geom = pattern.verts[:]
        geom.extend(pattern.edges[:])
        geom.extend(pattern.faces[:])
        bmesh.ops.bisect_plane(pattern,
            geom=geom,
            dist=0.001,
            plane_co=slice_co,
            plane_no=slice_no,
            use_snap_center=False,
            clear_outer=True,
            clear_inner=False
            )

        geom = [ed for ed in pattern.edges[:] if not ed.is_manifold]
        if len(geom) > 0:
            bmesh.ops.holes_fill(
                pattern,
                edges=geom,
                sides=len(geom)
                )

    def fit_roof_modifier(self, o, roof, rd, apply=False):
        if self.fit_roof and roof:
            target = rd.find_shrinkwrap(roof)
            modif = o.modifiers.new('Roof', 'SHRINKWRAP')
            modif.wrap_method = 'PROJECT'
            modif.vertex_group = 'Top'
            modif.use_project_z = True
            modif.use_negative_direction = True
            modif.target = target
            if apply:
                self.apply_modifier(o, modif)

    def add_child(self, name, wall_idx, pos, flip, flippable):
        # print("add_child %s %s" % (name, wall_idx))
        c = self.childs.add()
        c.child_name = name
        c.wall_idx = wall_idx
        c.pos = pos
        c.flip = flip
        c.flippable = flippable
        m = c.manipulators.add()
        m.type_key = 'DELTA_LOC'
        m.prop1_name = "x"
        m = c.manipulators.add()
        m.type_key = 'SNAP_SIZE_LOC'
        m.prop1_name = "x"
        m.prop2_name = "x"

    def _find_relocate_childs(self, o, childs):
        d, k = archipack_wall2_relocate_child.filter_child(self, o)
        if d:
            if k not in childs:
                childs[k] = {}
            childs[k][o.name] = (o, d)
        for c in o.children:
            self._find_relocate_childs(c, childs)

    def find_relocate_childs(self, o):
        # find childs and sort by kind
        childs = {}
        if o.parent:
            self._find_relocate_childs(o.parent, childs)
        return childs

    def add_relocate_child(
            self,
            tree,
            name,
            c_idx,
            pt,
            all_segs,
            maxdist,
            seg_idx=None,
            allow_single=False):
        """
         Find 2 closest segments
         Store distance to segments, segment indexes
         maxdist: distance of point from segment in front / back
         seg_idx: index of parent segment to project on closest
                  allow inner walls projection over outside
         allow_single: walls inside volume might have only isolated segment
                  this option enable link to isolated seg using a t and distance rule
        """
        closest = []

        count, selection = tree.intersects(pt, maxdist)

        # logger.debug("tree.intersects(%s) :%.4f seconds",  count, (time.time() - tim))

        for idx in selection:
            parent_name, parent_idx, seg = tree._geoms[idx]
            if parent_name != name and parent_name in all_segs:
                res, dist, t = seg.point_sur_segment(pt)
                abs_d = abs(dist)
                t_max = maxdist / seg.length

                # squared dist taking account of point location over segment
                d2 = dist ** 2
                if t > 1:
                    d2 += ((t - 1) * seg.length) ** 2
                elif t < 0:
                    d2 += (t * seg.length) ** 2

                if abs_d < maxdist and -t_max < t < 1 + t_max:
                    # logger.debug("%s %s %s %s", name, parent_name, d2, dist ** 2)
                    closest.append((d2, -dist, parent_name, parent_idx, t))

        # get 2 closest segments
        # childs walls use projection of themself through seg_idx

        if len(closest) > 1 or (seg_idx is not None and len(closest) > 0):
            closest.sort(key=lambda s: s[0])
            c = self.reloc_childs.add()
            c.child_name = name
            c.child_idx = c_idx
            if seg_idx is not None:
                # pick near segment with greateast angle (angle closest to 90)
                s0 = all_segs[name][seg_idx]
                i = 0
                a_min = pi
                for abs_d, d, p_name, p_idx, t in closest:
                    s1 = all_segs[p_name][p_idx]
                    a = abs(abs(s1.delta_angle(s0)) - pi / 2)
                    if a < a_min:
                        a_min = a
                        c.d0, c.parent_name0, c.seg0, c.t = d, p_name, p_idx, t
                # use child segment (eg for wall projection)
                c.d1, c.parent_name1, c.seg1 = 0, name, seg_idx
            else:
                d, d1, p1, i1, t1 = closest[0][0:5]
                # try to exclude colinear segments
                i = 1
                s0 = all_segs[p1][i1]
                d2, p2, i2 = closest[1][1:4]
                n_closest = len(closest) - 1
                a = abs(all_segs[p2][i2].delta_angle(s0))
                while i < n_closest and (a > 3.1405 or a < 0.001):
                    if closest[i][0] < d:
                        d, d1, p1, i1, t1 = closest[i][0:5]
                        s0 = all_segs[p1][i1]
                    i += 1
                    d2, p2, i2 = closest[i][1:4]
                    a = abs(all_segs[p2][i2].delta_angle(s0))

                c.d0, c.parent_name0, c.seg0, c.t = d1, p1, i1, t1
                c.d1, c.parent_name1, c.seg1 = d2, p2, i2

                # logger.debug("Multi {} c:{} p0:{} w0:{} p1:{} w1:{} d0:{} d1:{}".format(
                #   name, c_idx, c.parent_name0, c.seg0, c.parent_name1, c.seg1, c.d0, c.d1
                # ))

        elif allow_single and len(closest) > 0:
            # use wall's endpoint when isolated with a t + dist rule
            closest.sort(key=lambda s: s[0])
            c = self.reloc_childs.add()
            c.child_name = name
            c.child_idx = c_idx
            c.t = closest[0][4]
            c.d0, c.parent_name0, c.seg0 = closest[0][1:4]
            # logger.debug("Single {} c:{} p0:{} w0:{} d0:{} t:{}".format(
            #    name, c_idx, c.parent_name0, c.seg0, c.d0, c.t
            #    ))
        # else:
        #    logger.debug("Skip", name, len(closest), pt, seg_idx)

    def detect_openings(self, tree, name, pt, segs, maxdist, front, flippable=True):
        """
         Find closest segment
         Store distance to segments, segment indexes
        """
        closest = []
        count, selection = tree.intersects(pt, maxdist)
        # for parent_idx, seg in enumerate(segs):
        for idx in selection:
            parent_name, parent_idx, seg = tree._geoms[idx]
            res, dist, t = seg.point_sur_segment(pt)
            max_t = 2 * self.width / seg.length
            abs_d = abs(dist)
            if abs_d < maxdist and -max_t < t < 1 + max_t:
                closest.append((abs_d, -dist, parent_idx, t))

        # get closest segment
        if len(closest) > 0:
            closest.sort(key=lambda s: s[0])

            abs_d, dist, parent_idx, t = closest[0]

            seg = segs[parent_idx]

            # flippable objects (not always looking "outside")
            if flippable:
                dir = seg.sized_normal(t, 1).v
                flip = (front - dir).length > 0.5
            else:
                flip = self.flip

            m = self.childs_manipulators.add()
            m.type_key = 'DUMB_SIZE'

            return True, (
                name,
                parent_idx,
                (t * seg.length, dist, 0),
                flip,
                t,
                flippable)

        return False, ()

    def setup_childs(self, context, o, childs=None):
        """
            Store childs
            create manipulators
            call after a boolean oop

            Unlike other routines,
            generators use world coordsys
        """

        # logger.debug("Setup_childs %s", o.name)

        tim = time.time()
        self.childs.clear()
        self.reloc_childs.clear()
        self.childs_manipulators.clear()

        if o.parent is None:
            return 0

        # skip childs when not main wall
        skip_childs = True

        if childs is None:
            skip_childs = False
            childs = self.find_relocate_childs(o)

        logger.debug("Setup_childs %s  find_relocate_childs :%.4f seconds", o.name, (time.time() - tim))

        # retrieve objects to init quadtree
        objs = []
        for k, cat in childs.items():
            for name, child in cat.items():
                c, cd = child
                objs.append(c)

        # init a quadtree to minimize intersections tests on setup
        coordsys = CoordSys(objs)
        tree = Q_tree(coordsys)

        g = self.get_generator(o, axis=True)

        wall_with_childs = [0 for i in range(self.n_parts + 1)]
        relocate = []

        # g.segs contains explicit last segment
        # if not closed, remove last one
        if not self.closed:
            g.segs.pop()

        segs = [seg.line for seg in g.segs]

        # collect own walls segs
        all_segs = {}
        all_segs[o.name] = segs

        # Init Tree for openings
        for i, seg in enumerate(segs):
            tree.insert((o.name, i, seg))

        # walls_segs are generator ones
        walls_segs = {}
        walls_segs[o.name] = g.segs

        # walls segs left side including offset
        walls_offs = {}
        walls_offs[o.name] = segs

        logger.debug("Setup_childs %s  init tree :%.4f seconds", o.name, (time.time() - tim))

        # Order matter !
        # First update walls then other
        # as others depends on walls
        # so relocate_child order matter

        # detect openings on this wall
        # this envolve manipulators setup
        # each wall own her openings
        maxdist = self.width
        for k, cat in childs.items():
            if k == "archipack_window":
                for name, child in cat.items():
                    c, cd = child
                    p = c.matrix_world.translation.to_2d()
                    # outside Vector in y direction
                    tM = c.matrix_world.transposed()
                    front = -tM[1].to_2d()
                    res, reloc = self.detect_openings(tree, c.name, p, segs, maxdist, front, flippable=False)
                    if res:
                        relocate.append(reloc)
                        wall_with_childs[reloc[1]] = 1
                    # add 2 segs on both ends of openings
                    if cd.altitude == 0:
                        cg = cd.get_generator(c)
                        all_segs[name] = cg.segs

            elif k in ("archipack_door", "archipack_custom"):
                for name, child in cat.items():
                    c, cd = child
                    p = c.matrix_world.translation.to_2d()
                    # outside Vector in y direction
                    tM = c.matrix_world.transposed()
                    front = -tM[1].to_2d()
                    res, reloc = self.detect_openings(tree, c.name, p, segs, maxdist, front, flippable=True)
                    if res:
                        relocate.append(reloc)
                        wall_with_childs[reloc[1]] = 1
                        if k == "archipack_door":
                            # flip define material ids of hole
                            if cd.flip != reloc[3]:
                                self.select_object(context, c)
                                cd.flip = reloc[3]
                                self.unselect_object(c)
                    # add 2 segs on both ends of openings
                    if "archipack_door" in c.data or cd.altitude == 0:
                        cg = cd.get_generator(c)
                        all_segs[name] = cg.segs

        logger.debug("Setup_childs %s  process openings :%.4f seconds", o.name, (time.time() - tim))

        # sort openings along segments
        for child in sorted(relocate, key=lambda child: (child[1], child[4])):
            name, wall_idx, pos, flip, t, flippable = child
            # logger.debug("Setup_childs %s opening %s", o.name, name)
            self.add_child(name, wall_idx, pos, flip, flippable)

        logger.debug("Setup_childs %s  add childs :%.4f seconds", o.name, (time.time() - tim))

        # add a dumb size from last child to end of wall segment
        for i in range(sum(wall_with_childs)):
            m = self.childs_manipulators.add()
            m.type_key = 'DUMB_SIZE'

        # only run further on parent wall
        if skip_childs:
            logger.debug("Setup_childs %s  stop processing  :%.4f seconds\n", o.name, (time.time() - tim))
            # stop processing t childs
            return

        for name, segs in all_segs.items():
            if name != o.name:
                for i, seg in enumerate(segs):
                    tree.insert((name, i, seg))

        logger.debug("Setup_childs %s  build tree phase 2 :%.4f seconds", o.name, (time.time() - tim))

        # setup openings on other walls
        # get all walls segs so each wall is able to see all others
        if 'archipack_wall2' in childs:
            for name, child in childs['archipack_wall2'].items():
                c, cd = child
                if name != o.name:
                    # process only openings on other walls
                    cd.setup_childs(context, c, childs=childs)
                    # generator require valid manipulators
                    cd.setup_manipulators()
                    cg = cd.get_generator(c, axis=True)
                    if not cd.closed:
                        cg.segs.pop()
                    segs = [seg.line for seg in cg.segs]
                    walls_segs[name] = cg.segs
                    all_segs[name] = segs
                    # keep here to ensure stability
                    # against closest wall when editing inner wall
                    walls_offs[name] = segs

                    for i, seg in enumerate(segs):
                        tree.insert((name, i, seg))

            # Project t childs so they follow main without angle change
            # require a limited space check and use t child seg projection over parent
            # so find closest parent seg
            for name, child in childs['archipack_wall2'].items():
                if name != o.name:
                    c, cd = child
                    # check t childs ends segment only
                    segs = walls_segs[name]
                    s0 = segs[0]
                    s1 = segs[-1]
                    n_segs = len(segs)
                    self.add_relocate_child(
                            tree,
                            c.name,
                            0,
                            s0.p0,
                            all_segs,
                            maxdist,
                            seg_idx=0)

                    self.add_relocate_child(
                            tree,
                            c.name,
                            n_segs,
                            s1.p1,
                            all_segs,
                            maxdist,
                            seg_idx=n_segs - 1)
        logger.debug("Setup_childs %s  processing subs :%.4f seconds\n", o.name, (time.time() - tim))

        # curved segments are not supported
        for p in self.parts:
            if p.type == "C_SEG":
                return

        logger.debug("Setup_childs %s  processing soft :%.4f seconds\n", o.name, (time.time() - tim))

        # floor and moulding
        for k, cat in childs.items():
            if k in (
                    # "archipack_slab",
                    # "archipack_fence",
                    "archipack_molding",
                    "archipack_floor"):
                for name, child in cat.items():
                    c, cd = child
                    cd.setup_manipulators()
                    cg = cd.get_generator(c)
                    # allow single segments to support inside walls ends
                    for c_idx, seg in enumerate(cg.segs):
                        # point in world coordsys
                        p = seg.p0
                        self.add_relocate_child(
                            tree,
                            c.name,
                            c_idx,
                            p,
                            all_segs,
                            maxdist,
                            allow_single=True
                            )
                    if not cd.closed:
                        p = cg.segs[-1].p1
                        self.add_relocate_child(
                            tree,
                            c.name,
                            len(cg.segs),
                            p,
                            all_segs,
                            maxdist,
                            allow_single=True
                            )

        maxdist = 1e32
        if 'archipack_slab' in childs:
            # Explicit rule for slab, using wall segs only with larger dist
            # might use same rule for roof boundary cutters
            for name, child in childs['archipack_slab'].items():
                c, cd = child
                cd.setup_manipulators()
                cg = cd.get_generator(c)
                for c_idx, seg in enumerate(cg.segs):
                    # point in world coordsys
                    p = seg.p0
                    self.add_relocate_child(
                            tree,
                            c.name,
                            c_idx,
                            p,
                            walls_offs,
                            maxdist)

        # maxdist = 2 * self.width
        # NOTE:
        # cutters might move following parent
        # so we need to reloacate every points
        # in order to keep cutters location in synch
        # update cutters
        for k, cat in childs.items():
            if "cutter" in k:
                for name, child in cat.items():
                    c, cd = child
                    cd.setup_manipulators()
                    cg = cd.get_generator(c)
                    # allow single segments to support inside walls ends
                    for c_idx, seg in enumerate(cg.segs):
                        # point in world coordsys
                        p = seg.p0
                        self.add_relocate_child(
                            tree,
                            c.name,
                            c_idx,
                            p,
                            all_segs,
                            maxdist,
                            allow_single=True
                            )
                    if not cd.closed:
                        p = cg.segs[-1].p1
                        self.add_relocate_child(
                            tree,
                            c.name,
                            len(cg.segs),
                            p,
                            all_segs,
                            maxdist,
                            allow_single=True
                            )

        logger.debug("Setup childs end %s %.4f seconds\n", o.name, (time.time() - tim))

    def add_generator(self, name, c, d, generators, tM=None, force=False):
        if name != "" and c and (name not in generators or force):
            if tM is None:
                tM = c.matrix_world

            if archipack_wall2.filter(c):
                g = d.get_generator(tM, axis=True)
            else:
                g = d.get_generator(tM)

            generators[name] = [
                c,
                d,
                g,
                tM.inverted(),
                False
                ]

    def post_relocate(self, context):
        o = self.find_in_selection(context)
        if o is not None:
            logger.debug("post_relocate %s", o.name)
            self.relocate_childs(context, o)
            # update wall's finish
            if self.finish_enable and len(self.finish) > 0:
                self.update(context)

    def relocate_childs(self, context, o):
        """
            Move and resize childs after wall edition
            childs here are either doors or windows
            Unlike other routines,
            generators use world coordsys
            T childs walls only update doors and windows
        """
        tim = time.time()
        # throttle enable rules
        do_throttle_child = len(self.childs) > 10
        do_throttle_reloc = len(self.reloc_childs) > 5

        # throttle relocation of childs
        if do_throttle_child or do_throttle_reloc:
            throttle.add(context, o, self, 2.0, update_func="post_relocate")
            throttle.stack[o.name].update_func = "post_relocate"

        # only throttle openings when number of childs > 10
        if do_throttle_child and throttle.is_active(o.name):
            return

        logger.debug("Relocate_childs %s", o.name)
        g = self.get_generator(axis=True)

        side = -self.side

        tM = o.matrix_world
        loc, rot, scale = tM.decompose()

        # Generators: object, datablock, generator, mat world inverse, dirty flag
        generators = {}

        # relocate openings
        for child in self.childs:

            c, d = child.get_child(context)

            if c is None:
                continue

            # logger.debug("Relocate_opening %s %s", o.name, c.name)

            seg = g.segs[child.wall_idx].line
            t = child.pos.x / seg.length

            # y points inside
            # x is segment direction
            # when flip: -y -x
            n = seg.sized_normal(t, side)

            rx, ry = n.v
            rx, ry = ry, -rx
            if child.flippable:
                if child.flip:
                    rx, ry = -rx, -ry

            if d is not None:
                # print("change flip:%s width:%s" % (d.flip != child.flip, d.y != self.width))
                if d.y != self.width or d.flip != child.flip:
                    self.select_object(context, c, True)
                    d.auto_update = False
                    d.flip = child.flip
                    d.y = self.width
                    d.auto_update = True
                    self.unselect_object(c)
                x, y = n.p
            else:
                x, y = n.lerp(child.pos.y)

            self.select_object(context, o, True)
            # preTranslate
            wM = tM * Matrix([
                [rx, -ry, 0, x],
                [ry, rx, 0, y],
                [0, 0, 1, child.pos.z],
                [0, 0, 0, 1]
            ])
            c.matrix_world = wM
            if d is not None:
                self.add_generator(c.name, c, d, generators, tM=wM)
        logger.debug("Relocate_childs openings() :%.4f seconds", time.time() - tim)

        # Throttle relocation of reloc_childs
        if do_throttle_reloc and throttle.is_active(o.name):
            # adjust throttle timer
            # throttle.add(context, o, self, 2.0 * (time.time() - tim), update_func="post_relocate")
            logger.debug("Relocate_childs throttle(%s) :%.4f seconds", len(self.reloc_childs), time.time() - tim)
            return

        # build a dict for each child
        soft = {}

        # process childs in the right order
        child_names = []
        for child in self.reloc_childs:
            child_name = child.child_name
            if child_name not in soft:
                c, d = child.get_child(context)
                if c is None:
                    continue
                child_names.append(child_name)
                soft[child_name] = []
                self.add_generator(child_name, c, d, generators)

            soft[child_name].append(child)

            c, d = child.get_parent0(context)
            self.add_generator(child.parent_name0, c, d, generators)

            # parent_name1 could be empty
            c, d = child.get_parent1(context)
            self.add_generator(child.parent_name1, c, d, generators)
        logger.debug("Relocate_childs generators() :%.4f seconds", time.time() - tim)

        n_changes = 0

        for name in child_names:
            child = soft[name]
            # logger.debug("Relocate_childs start :%.2f seconds", time.time() - tim)
            # logger.debug("Relocate_childs %s child:%s", o.name, name)
            c, d, cg, itM, dirty = generators[name]
            if dirty:
                cg = d.get_generator(c)
                itM = c.matrix_world.inverted()
                generators[name] = [
                    c,
                    d,
                    cg,
                    itM,
                    False
                    ]
            # apply changes to generator
            n_segs = len(cg.segs)
            changed = False
            for cd in child:
                """
                if cd.parent_name1 != "":
                    logger.debug("Relocate Dual {} c:{} p0:{} w0:{} d0:{} p1:{} w1:{} d1:{}".format(
                        name,
                        cd.child_idx,
                        cd.parent_name0,
                        cd.seg0,
                        cd.d0,
                        cd.parent_name1,
                        cd.seg1,
                        cd.d1))
                else:
                    logger.debug("Relocate Single {} c:{} p0:{} w0:{} d0:{} t:{}".format(
                        name,
                        cd.child_idx,
                        cd.parent_name0,
                        cd.seg0,
                        cd.d0,
                        cd.t))
                """
                # find closest segments: use offset of main wall for t childs
                s0 = generators[cd.parent_name0][2].segs[cd.seg0].line.offset(cd.d0)

                if cd.parent_name1 == "":
                    # Single point rule (isolated inner wall end)
                    # use t and dist as rule to relocate
                    # p in world coordsysside *
                    p = s0.lerp(cd.t)
                else:
                    # Two segments intersection rule
                    # T childs use themself as 2nd segment without offset
                    if cd.parent_name1 != name:
                        s1 = generators[cd.parent_name1][2].segs[cd.seg1].line.offset(cd.d1)
                    else:
                        s1 = generators[cd.parent_name1][2].segs[cd.seg1].offset(cd.d1)
                    # p in world coordsys
                    res, p, u, v = s0.intersect_ext(s1)

                # if intersection fails p = 0
                if p != 0:
                    if cd.child_idx < n_segs:
                        if (cg.segs[cd.child_idx].p0 - p).length > 0.001:
                            n_changes += 1
                            changed = True
                        cg.segs[cd.child_idx].p0 = p
                        if cd.child_idx > 0:
                            cg.segs[cd.child_idx - 1].p1 = p
                        else:
                            d.move_object(c, p.to_3d())
                            # flag generator as dirty
                            generators[name][4] = True
                    else:
                        if (cg.segs[cd.child_idx - 1].p1 - p).length > 0.001:
                            n_changes += 1
                            changed = True
                        cg.segs[cd.child_idx - 1].p1 = p

            if changed:
                # update data from generator
                last = None
                for i, seg in enumerate(cg.segs):

                    p = d.parts[i]
                    if last is None:
                        # first a0 is absolute
                        p.a0 = cg.a0
                    else:
                        # a0 is delta between segs
                        p.a0 = seg.delta_angle(last)
                    if "C_" in p.type:
                        p.da = seg.da
                        p.radius = seg.r
                    else:
                        p.length = seg.length
                    last = seg

                # logger.debug("Relocate_childs change :%.2f seconds", time.time() - tim)
                self.select_object(context, c, True)
                d.update(context)
                # logger.debug("Relocate_childs update :%.2f seconds", time.time() - tim)
                self.unselect_object(c)

        logger.debug("Relocate_childs (%s of %s) :%.4f seconds", n_changes, len(self.reloc_childs), time.time() - tim)

        self.select_object(context, o, True)

    def update_childs(self, context, o, g):
        """
            setup gl points for childs
        """
        # print("update_childs")

        if o.parent is None:
            return

        # swap manipulators so they always face outside
        manip_side = self.side

        itM = o.matrix_world.inverted()
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
                        pt = (itM * c.matrix_world.translation).to_2d()
                        res, y, t = wall.point_sur_segment(pt)
                        # child.pos = (wall.length * t, y, child.pos.z)
                        p1 = wall.lerp(t - dt)
                        # dumb size between childs
                        self.childs_manipulators[m_idx].set_pts([
                            (p0.x, p0.y, 0),
                            (p1.x, p1.y, 0),
                            (0.5 * manip_side, 0, 0)])

                        m_idx += 1
                        x, y = 0.5 * d.x, 0.5 * d.y * self.x_offset + self.offset

                        # delta loc
                        child.manipulators[0].set_pts([
                            (-x, y, 0),
                            (x, y, 0),
                            (1, 0, 0)])

                        # loc size
                        child.manipulators[1].set_pts([
                            (-x, y, 0),
                            (x, y, 0),
                            (0.5, 0, 0)])

                        p0 = wall.lerp(t + dt)

            p1 = wall.lerp(1)
            if wall_has_childs:
                # dub size after all childs
                self.childs_manipulators[m_idx].set_pts([
                    (p0.x, p0.y, 0),
                    (p1.x, p1.y, 0),
                    (manip_side * 0.5, 0, 0)])
                m_idx += 1

    def remove_dimension(self, context, o):
        o.select = True
        context.scene.objects.active = o
        if bpy.ops.archipack.dimension_auto.poll():
            bpy.ops.archipack.dimension_auto(mode='DELETE')

    def update_dimension(self, context, o, g, synch_childs=True):

        # swap manipulators so they always face outside
        # dims = [child
        #    for child in o.children
        #    if child.data and "archipack_dimension_auto" in child.data
        #    ]

        dims = {child.data.archipack_dimension_auto[0].parent_uid: child
            for child in o.children
            if child.data and "archipack_dimension_auto" in child.data
            }

        if self.dimensions:

            # pts = [seg.p0.to_3d() for seg in g.segs]
            # is_cw = self.is_cw(pts)

            _parts = self.parts
            n_parts = len(_parts)
            # add right side measure points
            if not self.flip:
                offset = self.left_offset + self.width
                g.set_offset(offset)
                g.close(offset)
            _segs = g.segs
            for i, f in enumerate(_segs):
                part = _parts[i]
                if part.uid == 0:
                    self.create_uid(part, increment=2)
                # Dimensions points
                self.add_dimension_point(part.uid, f.line.p0.to_3d())
                if i == n_parts:
                    self.add_dimension_point(1 + part.uid, f.line.p1.to_3d())

            force_update = False
            for wall_idx, wall in enumerate(g.segs):

                parent_uid = self.parts[wall_idx].uid

                if parent_uid in dims:
                    dim = dims[parent_uid]
                    # remove last part dim when not closed
                    if self.closed or wall_idx < self.n_parts:
                        del dims[parent_uid]
                    self.select_object(context, dim, True)
                else:
                    # dont create missing dim for last segment
                    # when not closed
                    if not self.closed and wall_idx >= self.n_parts:
                        continue

                    bpy.ops.archipack.dimension_auto(
                        distance=1,
                        auto_parent=False,
                        flip_side=True,
                        auto_manipulate=False)

                    dim = context.active_object
                    dim.parent = o
                    dim.matrix_world = o.matrix_world.copy()

                dim.location = wall.line.p0.to_3d()

                rot_z = wall.angle

                if self.flip:
                    rot_z += pi

                dim.rotation_euler.z = rot_z

                # force dim matrix_world update
                if parent_uid not in dims:
                    force_update = True

                dim_d = dim.data.archipack_dimension_auto[0]
                dim_d.parent_uid = parent_uid

                dim_d.add_source(o, self.parts[wall_idx].uid, 'WALL2')

                # TODO:
                # when open last is wall_idx + 1
                # when closed last is parts[0]
                if self.closed and wall_idx >= self.n_parts:
                    dim_d.add_source(o, self.parts[0].uid, 'WALL2')
                elif wall_idx < self.n_parts:
                    dim_d.add_source(o, self.parts[wall_idx + 1].uid, 'WALL2')
                else:
                    dim_d.add_source(o, 1 + self.parts[wall_idx].uid, 'WALL2')

                for child in self.childs:
                    if child.wall_idx == wall_idx:
                        c, d = child.get_child(context)
                        if d is not None:
                            # child is either a window or a door
                            if "archipack_window" in c.data:
                                provider_type = 'WINDOW'
                            elif "archipack_door" in c.data:
                                provider_type = 'DOOR'
                            elif "archipack_custom" in c.data:
                                provider_type = 'CUSTOM'
                            dim_d.add_source(c, 0, provider_type)
                            dim_d.add_source(c, 1, provider_type)

                # dim_d.update(context)
                self.unselect_object(dim)

            if force_update:
                context.scene.update()

        for parent_uid in dims:
            if parent_uid != 0:
                dim = dims[parent_uid]
                self.delete_object(context, dim)

        self.select_object(context, o, True)

    def synch_dimension(self, context, o):
        """
          Synch dimension when manipulating window or doors
        """
        if self.dimensions:
            self.update_parts()
            g = self.get_generator()
            # self.setup_childs(o, g, maxiter=0)
            # self.update_childs(context, o, g)
            # prevent cyclic call
            self.update_dimension(context, o, g, synch_childs=False)

    def manipulate_childs(self, context):
        """
            setup child manipulators
        """
        # print("manipulate_childs")
        n_parts = self.n_parts
        if self.closed:
            n_parts += 1

        for wall_idx in range(n_parts):
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
                if (archipack_wall2.datablock(c) == self):
                    self.select_object(context, c, True)
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

            g = self.get_generator(axis=True)

            itM = o.matrix_world.inverted()
            for child in self.childs:
                c, d = child.get_child(context)
                if d is not None:
                    wall = g.segs[child.wall_idx].line
                    pt = (itM * c.matrix_world.translation).to_2d()
                    res, d, t = wall.point_sur_segment(pt)
                    child.pos = (t * wall.length, d, child.pos.z)

            # NOTE: object dosent update itself
            # explicit relocate
            self.relocate_childs(context, o)
            # update childs manipulators
            self.update_childs(context, o, g)
            self.update_dimensions(context, o)

    def manipulable_move_t_part(self, context, o=None, manipulator=None):
        """
            Callback for t_parts childs
        """
        type_name = type(manipulator).__name__
        # print("manipulable_manipulate %s" % (type_name))
        if type_name in [
                'DeltaLocationManipulator'
                ]:
            # update manipulators pos of childs
            if archipack_wall2.datablock(o) != self:
                return
            # update childs
            self.relocate_childs(context, o)

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
        n_parts = self.n_parts
        if self.closed:
            n_parts += 1

        # update manipulators on version change
        self.setup_manipulators()

        for i, part in enumerate(self.parts):

            if i < n_parts:
                if i > 0:
                    # start angle
                    self.manip_stack.append(part.manipulators[0].setup(context, o, part))

                # length / radius + angle
                self.manip_stack.append(part.manipulators[1].setup(context, o, part))
                # segment index
                self.manip_stack.append(part.manipulators[3].setup(context, o, self))

            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))

        # height as per segment will be here when done

        # width + counter
        for i, m in enumerate(self.manipulators):
            # skip counter on closed walls
            if i == 1 and self.closed:
                continue
            elif i == 3:
                # t child location along parent
                self.manip_stack.append(m.setup(context, o, self, self.manipulable_move_t_part))
            else:
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
            return False

        # self.manip_stack = []
        o = context.active_object
        # setup childs manipulators
        self.setup_childs(context, o)
        # if not res:
        #    return None
        # store gl points
        g = self.get_generator()
        self.update_childs(context, o, g)
        # dont do anything ..

        self.manipulable_setup(context)
        self.manipulate_mode = True

        self._manipulable_invoke(context)

        return True

    def find_roof(self, o):
        # TODO:
        # Handle multiple roofs
        p = self.get_topmost_parent(o)
        for c in p.children:
            if c.data is not None and "archipack_roof" in c.data:
                return c, c.data.archipack_roof[0]
        return None, None

    def get_childs_openings(self, context, wd, objs):
        # openings belong to this wall
        for child in wd.childs:
            c, d = child.get_child(context)
            if d is not None and "archipack_wall2" not in c.data:
                objs.append(c)

    def get_childs_geoms(self, context, o, wd, objs, t_childs):

        # rely on real childs
        childs = self.find_relocate_childs(o)
        # childs[k][c.name] = (c, d)
        if 'archipack_wall2' in childs:
            for name, child in childs['archipack_wall2'].items():
                if name != o.name:
                    c, d = child
                    # d.update_parts()
                    # update child so it find openings
                    # d.setup_childs(context, c)
                    t_childs.append((c, d))
                    objs.append(c)
                    self.get_childs_openings(context, d, objs)

    def determine_extend(self, context, coordsys, _factory, c, td, walls):
        """
         Perform neighboor analysis
         prevent precision issues
         by extending wall ends
         when points near ends are inside another wall
        """
        td.extend = 0
        if not td.closed:
            # extend a bit to prevent precision issue
            t_g = td.get_generator()
            side = 1
            if td.flip:
                side = -1
            # axis
            offset = side * 0.5 * td.width * td.x_offset
            # points outside wall ends by distance
            p0, p1 = t_g.get_ends(offset, c, coordsys.invert, 0.01)
            pt = _factory.createPoint(_factory.createCoordinate(p0))
            for poly in walls:
                if poly.contains(pt):
                    td.extend += 1
                    break
            pt = _factory.createPoint(_factory.createCoordinate(p1))
            for poly in walls:
                if poly.contains(pt):
                    td.extend += 2
                    break

    def as_geom(self, context, o, mode, inter, doors, windows, io=None):
        """
         Build 2d symbol of walls as pygeos entity for further processing
         cut windows and doors
         w, it = C.object.data.archipack_wall2[0].as_2d(C, C.object)
        """

        # objs contains wall and openings wich belongs to this wall
        objs = [o]

        self.update_parts()
        g = self.get_generator(axis=False)

        # t childs contains childs walls
        t_childs = []

        if mode in {
                'SYMBOL', 'MERGE', 'CONVERT',
                'FLOORS', 'FLOORS_CHILD',
                'FLOORS_MOLDINGS', 'FLOORS_MOLDINGS_CHILD'
                }:

            # MUST disable manipulators as setup_childs update
            # data structure and lead in crash
            bpy.ops.archipack.disable_manipulate()
            # setup all childs so we are able to see them
            # as this only occurs when manipulated
            if io is None:
                # setup childs is recursive
                self.setup_childs(context, o)

            # collect windows, doors for this wall
            self.get_childs_openings(context, self, objs)

            # collect other walls and openings only for main wall
            if io is None:
                self.get_childs_geoms(context, o, self, objs, t_childs)

        # t childs walls
        t_walls = []

        # Init coordsys and Io for main wall
        if io is None:
            # main wall
            coordsys = Io.getCoordsys(objs)
            io = Io(scene=context.scene, coordsys=coordsys)
            """
             Perform neighboor analysis
             prevent precision issues
             by extending wall ends
             when points near ends are inside another wall
            """
            walls = []

            io, wall, childs = self.as_geom(context, o, 'BOTH', [], [], [], io)
            walls.append(wall)
            for c, td in t_childs:
                io, wall, childs = td.as_geom(context, c, 'BOTH', [], [], [], io)
                walls.append(wall)

            self.determine_extend(context, coordsys, wall._factory, o, self, walls)
            for c, td in t_childs:
                self.determine_extend(context, coordsys, wall._factory, c, td, walls)

            # process t_childs
            if mode in {'SYMBOL', 'MERGE', 'CONVERT'}:
                # accumulate walls and intersections
                for c, td in t_childs:
                    io, wall, childs = td.as_geom(context, c, mode, inter, doors, windows, io)
                    t_walls.append(wall)
            elif mode in {'FLOORS', 'FLOORS_MOLDINGS'}:
                for c, td in t_childs:
                    io, wall, childs = td.as_geom(context, c, "{}_CHILD".format(mode), inter, doors, windows, io)
                    t_walls.append(wall)

        side = 1
        if self.flip:
            side = -1
        z = 0

        # coords of wall boundarys
        if mode in {
                'SYMBOL', 'MERGE', 'CONVERT',
                'OUTSIDE', 'BOTH',
                'FLOORS', 'FLOORS_CHILD',
                'FLOORS_MOLDINGS', 'FLOORS_MOLDINGS_CHILD'
                }:

            offset = self.offset + side * 0.5 * (1 + self.x_offset) * self.width
            outside = []
            g.get_coords(offset, z, self, outside)

        if mode in {
                'SYMBOL', 'MERGE', 'CONVERT',
                'INSIDE', 'BOTH',
                'FLOORS', 'FLOORS_CHILD',
                'FLOORS_MOLDINGS', 'FLOORS_MOLDINGS_CHILD'
                }:

            offset = self.offset - side * 0.5 * (1 - self.x_offset) * self.width
            inside = []
            g.get_coords(offset, z, self, inside)

        # reset extend flag
        self.extend = 0

        # create boundary line
        if self.closed:
            if (mode in {'SYMBOL', 'MERGE', 'CONVERT', 'BOTH', 'FLOORS_CHILD', 'FLOORS_MOLDINGS_CHILD'} or
                (FLOORS_2D_SUPPORT_MULTI and mode in {'FLOORS', 'FLOORS_MOLDINGS'})):
                wall = io.coords_to_polygon(o.matrix_world, outside, [inside], force_2d=True)
            elif mode in {'INSIDE', 'FLOORS', 'FLOORS_MOLDINGS'}:
                wall = io.coords_to_polygon(o.matrix_world, inside, force_2d=True)
            else:
                wall = io.coords_to_polygon(o.matrix_world, outside, force_2d=True)
        else:
            if mode in {
                    'SYMBOL', 'MERGE', 'CONVERT',
                    'BOTH',
                    'FLOORS', 'FLOORS_MOLDINGS',
                    'FLOORS_CHILD', 'FLOORS_MOLDINGS_CHILD'
                    }:
                outside.extend(list(reversed(inside)))
                wall = io.coords_to_polygon(o.matrix_world, outside, force_2d=True)
            elif mode == 'INSIDE':
                wall = io.coords_to_linestring(o.matrix_world, [inside], force_2d=True)
            else:
                wall = io.coords_to_linestring(o.matrix_world, [outside], force_2d=True)

        # merge t childs walls
        if mode in {'SYMBOL', 'MERGE'} or (FLOORS_2D_SUPPORT_MULTI and mode in {'FLOORS', 'FLOORS_MOLDINGS'}):
            # SYMBOL -> union t child walls
            for tw in t_walls:
                wall = wall.union(tw)
            # io.output(wall, name="wall_merge", multiple=False)

        elif mode in {'FLOORS', 'FLOORS_MOLDINGS'}:
            # FLOORS -> difference t child walls
            for tw in t_walls:
                wall = wall.difference(tw)
        elif mode == 'CONVERT':
            # merge walls into a single multipolygon
            t_childs.insert(0, (o, self))
            t_walls.insert(0, wall)
            wall = wall._factory.buildGeometry(t_walls)

        holes = []

        # get inside polygons
        if FLOORS_2D_SUPPORT_MULTI and mode in {'FLOORS', 'FLOORS_MOLDINGS'}:
            # wall is either a polygon or a multipolygon
            t_walls = []
            if wall.type_id == 6:
                polys = wall.geoms
            else:
                polys = [wall]

            # collect holes (inner walls not touching outside ones)
            holes = [
                poly
                for poly in polys
                if len(poly.interiors) < 1
                ]
            # print("len(holes): %s" % len(holes))

            for i, poly in enumerate(polys):
                # extract inside as polygons
                for interior in poly.interiors:

                    tmp = wall._factory.createPolygon(exterior=interior)
                    # find if any hole belong to this poly interior (wall inside but not touching)
                    contains = []
                    found = []
                    for j, hole in enumerate(holes):
                        if tmp.contains(hole):
                            contains.append(hole.exterior)
                            found.append(j)
                    for j in reversed(found):
                        holes.pop(j)
                    t_walls.append(
                        wall._factory.createPolygon(exterior=interior, interiors=contains)
                        )
                    # io.output(t_walls[-1], name="wall_inside {}".format(i), multiple=False)
            # add isolated walls
            t_walls.extend(holes)

            wall = wall._factory.buildGeometry(t_walls)
            # io.output(wall, name="wall_inside", multiple=False)

        if mode == 'FLOORS_MOLDINGS':
            # convert wall polygons to linestrings
            # so boolean apply on lines instead of areas
            # Keep track of wall polys so we are able to
            # know wich side is inside and revert lines according
            coords = []
            # type id 6: multipolygon
            if wall.type_id == 6:
                polys = wall.geoms
            else:
                polys = [wall]
            for poly in polys:
                coords.extend([p.coords for p in poly.interiors])
                coords.append(poly.exterior.coords)
            lines = [wall._factory.createLineString(coord) for coord in coords]
            wall = wall._factory.buildGeometry(lines)

        # openings
        if mode in {
                'SYMBOL',
                'FLOORS', 'FLOORS_MOLDINGS',
                'FLOORS_CHILD', 'FLOORS_MOLDINGS_CHILD'
                }:
            # collect child holes as polygons
            for child in self.childs:
                c, d = child.get_child(context)
                if d is not None:
                    # ignore windows when altitude > 0
                    if mode in {'FLOORS', 'FLOORS_CHILD'}:
                        if "archipack_window" in c.data:
                            if d.altitude > 0:
                                continue
                            h = d.hole_2d(mode)
                            hole = io.coords_to_polygon(c.matrix_world, h)
                            windows.append(hole)
                        # identify door side
                        if "archipack_door" in c.data:
                            # identify door side to get right
                            # hole for this wall
                            h = d.hole_2d('INSIDE')
                            hole = io.coords_to_polygon(c.matrix_world, h)
                            pt = hole._factory.createPoint(hole.coords[3])
                            doors.append((pt, hole))
                            h = d.hole_2d('OUTSIDE')
                            hole = io.coords_to_polygon(c.matrix_world, h)
                            pt = hole._factory.createPoint(hole.coords[0])
                            doors.append((pt, hole))
                    elif mode in {'FLOORS_MOLDINGS', 'FLOORS_MOLDINGS_CHILD'}:
                        # symbols and floor moldings require holes
                        if "archipack_window" in c.data and d.altitude > 0:
                            continue
                        h = d.hole_2d(mode)
                        hole = io.coords_to_polygon(c.matrix_world, h)
                        if "archipack_door" in c.data:
                            doors.append(hole)
                        else:
                            windows.append(hole)
                    else:
                        # symbols and floor moldings require holes
                        h = d.hole_2d(mode)
                        hole = io.coords_to_polygon(c.matrix_world, h)
                        if "archipack_door" in c.data:
                            doors.append(hole)
                        else:
                            windows.append(hole)

            # io._to_curve(wall, "Wall", '2D')
            # substract holes
            if len(windows) > 0:
                basis = windows[0]._factory.buildGeometry(windows)
                # io._to_curve(basis, "Holes", '2D')
                if mode == 'SYMBOL':
                    inter.append(wall.intersection(basis))
                    wall = wall.difference(basis)
                elif mode == 'FLOORS':
                    if FLOORS_2D_SUPPORT_MULTI:
                        # identify wich polygon the window belong to
                        # wall is either a multipolygon (type_id 6)
                        # or a polygon
                        if wall.type_id == 6:
                            polys = wall.geoms
                        else:
                            polys = [wall]
                        for i, poly in enumerate(polys):
                            for j, window in enumerate(windows):
                                if poly.intersects(window):
                                    # poly intersects
                                    if not poly.contains(window):
                                        # if window not wholy contained into poly
                                        # grow surface inside under window
                                        poly = poly.union(window)
                                        polys[i] = poly
                        if len(polys) == 1:
                            # single polygon, so wall is that one
                            wall = polys[0]
                    else:
                        wall = wall.union(basis)

                # elif mode == 'FLOORS_CHILD':
                #     wall = wall.difference(basis)
                elif mode == 'FLOORS_MOLDINGS':
                    wall = wall.difference(basis)

            # process doors only in inside openings mode
            if len(doors) > 0:
                # io._to_curve(basis, "Holes", '2D')
                if mode == 'SYMBOL':
                    basis = doors[0]._factory.buildGeometry(doors)
                    wall = wall.difference(basis)
                elif mode == 'FLOORS':
                    # identify wich polygon the door belong to
                    # wall is either a multipolygon (type_id 6)
                    # or a polygon
                    if wall.type_id == 6:
                        polys = wall.geoms
                    else:
                        polys = [wall]
                    for i, poly in enumerate(polys):
                        # @TODO:
                        # Doors inside on T child may touch 2 polygons (straight unclosed wall)
                        # in such case make an union of door inside and poly
                        for j, pt_door in enumerate(doors):
                            pt, door = pt_door
                            if poly.intersects(door):
                                # poly intersects
                                # if point inside poly then add to poly
                                if not poly.contains(door):
                                    if poly.contains(pt):
                                        poly = poly.union(door)
                                        polys[i] = poly
                                    else:
                                        # if door not wholy contained into poly
                                        # substract
                                        poly = poly.difference(door)
                                        polys[i] = poly
                    if len(polys) == 1:
                        # single polygon, so wall is that one
                        wall = polys[0]
                elif mode == 'FLOORS_MOLDINGS':
                    basis = doors[0]._factory.buildGeometry(doors)
                    wall = wall.difference(basis)

        # Boolean result might not always preserve wall direction
        # so ensure lines are in the right direction
        if mode == 'FLOORS_MOLDINGS':
            merged = wall.line_merge()
            # merged is an array of geoms
            # switch direction of line so molden are always outside of wall
            for line in merged:
                co = line.coords
                # a point in the middle of first segment at 0.01 * width on right side of line
                n = 0.01 * self.width * Vector((co[1].y - co[0].y, co[0].x - co[1].x, 0)).normalized()
                pos = Vector((0.5 * (co[0].x + co[1].x), 0.5 * (co[0].y + co[1].y), 0)) + n
                pt = wall._factory.createPoint(wall._factory.createCoordinate(pos))
                swap = True
                for poly in polys:
                    if poly.contains(pt):
                        swap = False
                        break
                if swap:
                    line.coords = wall._factory.coordinateSequenceFactory.create(line.coords[::-1])
            wall = wall._factory.buildGeometry(merged)

        # output only geometry from there
        # for further processing of t_childs
        return io, wall, t_childs


class ARCHIPACK_PT_wall2(Panel):
    bl_idname = "ARCHIPACK_PT_wall2"
    bl_label = "Wall"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    def draw(self, context):
        d = archipack_wall2.datablock(context.object)
        if d is None:
            return
        layout = self.layout
        row = layout.row(align=True)
        row.operator("archipack.manipulate", icon='HAND')
        row.operator("archipack.wall2", text="Delete", icon='ERROR').mode = 'DELETE'
        # row = layout.row(align=True)
        # row.d(d, 'realtime')
        # layout.label(text="Manip mode:{}".format(d.manipulate_mode))
        box = layout.box()
        d.template_user_path(context, box)
        box = layout.box()
        box.prop_search(d, "t_part", context.scene, "objects", text="T parent", icon='OBJECT_DATAMODE')

        box = layout.box()
        box.prop(d, 'auto_synch', icon="AUTO", emboss=True)
        box.prop(d, 'width_ui')
        box.prop(d, 'z')
        box.prop(d, 'z_offset')
        box.prop(d, 'offset')
        row = box.row()
        row.prop(d, 'base_line', text="")
        row.prop(d, 'flip')
        box.prop(d, 'step_angle')
        row = box.row()
        row.operator("archipack.wall2_fit_roof")
        row.prop(d, "fit_roof", text="Auto fit", icon="AUTO", emboss=True)
        # box = layout.box()
        #
        d.template_parts(context, layout, draw_type=True)

        box = layout.box()
        row = box.row()
        icon = "TRIA_RIGHT"
        if d.finish_expand:
            icon = "TRIA_DOWN"
        row.prop(d, 'finish_expand', icon=icon, text="Finishings ({})".format(len(d.finish)), icon_only=True, emboss=True)
        row.prop(d, 'finish_enable', text="")
        row.operator("archipack.wall2_add_finish", text="", icon="ZOOMIN")

        if d.finish_expand:
            for index, finish in enumerate(d.finish):
                box = layout.box()
                finish.draw(context, box, index)

        box = layout.box()
        box.prop(d, "dimensions")
        box.label(text="Create curves")
        row = box.row(align=True)
        row.operator("archipack.wall2_to_curve", text="Symbol").mode = 'SYMBOL'
        row.operator("archipack.wall2_to_curve", text="Floor").mode = 'FLOORS'
        row.operator("archipack.wall2_to_curve", text="Molding").mode = 'FLOORS_MOLDINGS'
        row = box.row(align=True)
        row.operator("archipack.wall2_to_curve", text="Bounds").mode = 'MERGE'
        row.operator("archipack.wall2_to_curve", text="In").mode = 'INSIDE'
        row.operator("archipack.wall2_to_curve", text="Out").mode = 'OUTSIDE'
        # row.operator("archipack.wall2_fit_roof", text="Inside").inside = True

    @classmethod
    def poll(cls, context):
        return archipack_wall2.poll(context.active_object)


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_wall2(ArchipackCreateTool, Operator):
    bl_idname = "archipack.wall2"
    bl_label = "Wall"
    bl_description = "Delete wall and childs"
    bl_category = 'Archipack'

    mode = EnumProperty(
        items=(
        ('CREATE', 'Create', '', 0),
        ('DELETE', 'Delete', '', 1)
        ),
        default='CREATE'
        )

    def delete(self, context):
        o = context.active_object
        if archipack_wall2.filter(o):
            bpy.ops.archipack.disable_manipulate()
            self.delete_object(context, o)

    def create(self, context):
        m = bpy.data.meshes.new("Wall")
        o = bpy.data.objects.new("Wall", m)
        d = m.archipack_wall2.add()
        d.manipulable_selectable = True
        # Link object into scene
        self.link_object_to_scene(context, o)

        # select and make active
        self.select_object(context, o, True)

        # around 12 degree
        m.auto_smooth_angle = 0.20944
        self.select_object(context, o, True)
        self.load_preset(d)
        self.add_material(o)
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            if self.mode == 'CREATE':
                bpy.ops.object.select_all(action="DESELECT")
                o = self.create(context)
                o.location = bpy.context.scene.cursor_location
                # select and make active
                self.select_object(context, o, True)
                self.manipulate()
            else:
                self.delete(context)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_wall2_from_curve(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.wall2_from_curve"
    bl_label = "Wall curve"
    bl_description = "Create wall(s) from a curve"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'INTERNAL', 'UNDO'}

    auto_manipulate = BoolProperty(default=False)

    @classmethod
    def poll(self, context):
        return context.active_object is not None and context.active_object.type == 'CURVE'

    def create(self, context):
        curve = context.active_object
        for spline in curve.data.splines:
            bpy.ops.archipack.wall2(auto_manipulate=self.auto_manipulate)
            o = context.active_object
            d = archipack_wall2.datablock(o)
            d.from_spline(context, curve.matrix_world, 12, spline)

        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            # select and make active
            self.select_object(context, o, True)

            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_wall2_from_slab(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.wall2_from_slab"
    bl_label = "->Wall"
    bl_description = "Create a wall from a slab"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'INTERNAL', 'UNDO'}

    auto_manipulate = BoolProperty(default=True)

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o is not None and o.data is not None and 'archipack_slab' in o.data

    def create(self, context):
        slab = context.active_object
        wd = slab.data.archipack_slab[0]
        bpy.ops.archipack.wall2(auto_manipulate=self.auto_manipulate)
        o = context.active_object
        d = archipack_wall2.datablock(o)
        d.auto_update = False
        d.parts.clear()
        d.n_parts = wd.n_parts - 1
        d.closed = True
        for part in wd.parts:
            p = d.parts.add()
            p.type = part.type
            p.length = part.length
            p.radius = part.radius
            p.da = part.da
            p.a0 = part.a0
        # select and make active
        self.select_object(context, o, True)

        d.auto_update = True

        # pretranslate
        o.matrix_world = Matrix([
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, d.z_offset],
                [0, 0, 0, 1],
                ]) * slab.matrix_world
        bpy.ops.object.select_all(action='DESELECT')
        self.select_object(context, slab)
        self.select_object(context, o, True)
        bpy.ops.archipack.add_reference_point()

        self.unselect_object(slab)
        self.unselect_object(o.parent)

        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)

            # select and make active
            self.select_object(context, o, True)

            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_wall2_fit_roof(Operator):
    bl_idname = "archipack.wall2_fit_roof"
    bl_label = "Fit roof"
    bl_description = "Fit roof"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'INTERNAL', 'UNDO'}

    inside = BoolProperty(default=False)
    skip_z = BoolProperty(default=False)
    auto_update = BoolProperty(default=True)

    @classmethod
    def poll(self, context):
        return archipack_wall2.poll(context.active_object)

    def execute(self, context):
        o = context.active_object
        d = archipack_wall2.datablock(o)
        r, rd = d.find_roof(o)
        if rd is not None:
            # d.setup_childs(o)
            rd.make_wall_fit(context,
                r,
                o,
                inside=self.inside,
                auto_update=self.auto_update,
                skip_z=self.skip_z)

        return {'FINISHED'}


class ARCHIPACK_OT_wall2_to_curve(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.wall2_to_curve"
    bl_label = "To curve"
    bl_description = "Create curve from wall"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'INTERNAL', 'UNDO'}
    mode = EnumProperty(
        items=(
            ('BOTH', 'Both sides', 'Wall boundary'),
            ('INSIDE', 'Inside', 'Wall inside'),
            ('FLOORS', 'Floors', 'Floors'),
            ('OUTSIDE', 'Outside', 'Wall outside'),
            ('SYMBOL', 'Symbol', 'Wall 2d Symbol'),
            ('MERGE', 'Merge', 'Merge 2d walls boundary'),
            ('CONVERT',
                'Convert',
                'Join 3d walls as a single mesh (BWARE !, use CARVE mode of boolean so might crash blender)'),
            ('FLOORS_MOLDINGS', 'Floor moldings', 'Floor moldings')
        ),
        default='SYMBOL'
        )

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            d = archipack_wall2.datablock(o)
            if d is None:
                return {'CANCELLED'}
            sel = []
            inter = []
            # doors holes
            doors = []
            windows = []

            io, wall, t_childs = d.as_geom(context, o, self.mode, inter, doors, windows)
            # windows openings for symbols
            if len(inter) > 0:
                inter = inter[0]._factory.buildGeometry(inter)
                res = io._to_curve(inter, "{}-w-2d".format(o.name), '2D')
                sel.append(res)

            if self.mode == 'CONVERT':
                # build walls as separated geometry and merge using a boolean
                roof, rd = d.find_roof(o)
                m = bpy.data.meshes.new("Wall")
                new_o = bpy.data.objects.new("Wall", m)
                new_o.matrix_world = io.coordsys.world.copy()
                # Link object into scene
                self.link_object_to_scene(context, new_o)

                if wall.geom_type == 'MultiPolygon':
                    polys = wall.geoms
                else:
                    polys = [wall]

                for i, poly in enumerate(polys):
                    # setup a random z diff for childs walls so boolean dosent fails
                    c, cd = t_childs[i]
                    dz = c.matrix_world.translation.z - o.matrix_world.translation.z - cd.z_offset

                    walls = Io.to_wall(context, io.coordsys, poly, cd.z, name="Wall", walls=[], clean=False)
                    for w in walls:
                        w.location.z += dz
                        # select and make active
                        self.select_object(context, w, True)
                        Io.assign_matindex_to_wall(w)

                    if cd.fit_roof:
                        for w in walls:
                            cd.fit_roof_modifier(w, roof, rd, apply=True)

                    self.select_object(context, new_o, True)

                    for w in walls:
                        modif = new_o.modifiers.new('AutoMerge', 'BOOLEAN')
                        modif.operation = 'UNION'
                        if hasattr(modif, 'solver'):
                            modif.solver = 'CARVE'
                        modif.object = w
                        d.apply_modifier(new_o, modif)
                        self.delete_object(context, w)

                self.select_object(context, new_o, True)
                bpy.ops.archipack.wall(z=d.z)

                self.select_object(context, new_o, True)
                bpy.ops.object.mode_set(mode='EDIT')
                bpy.ops.mesh.select_all(action='SELECT')
                bpy.ops.mesh.dissolve_limited(angle_limit=0.0174533, delimit={'MATERIAL'})
                bpy.ops.object.mode_set(mode='OBJECT')

                if bpy.ops.archipack.auto_boolean.poll():
                    bpy.ops.archipack.auto_boolean()

            else:
                res = io._to_curve(wall, "{}-2d".format(o.name), '2D')

                sel.append(res)
                bpy.ops.object.select_all(action="DESELECT")
                for o in sel:
                    self.select_object(context, o)
                self.select_object(context, res, True)

            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_wall2_add_finish(Operator):
    bl_idname = "archipack.wall2_add_finish"
    bl_label = "Add finishings"
    bl_description = "Add finishing"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'INTERNAL', 'UNDO'}

    @classmethod
    def poll(self, context):
        return archipack_wall2.poll(context.active_object)

    def execute(self, context):
        o = context.active_object
        d = archipack_wall2.datablock(o)
        d.finish.add()
        d.update(context)
        return {'FINISHED'}


class ARCHIPACK_OT_wall2_remove_finish(Operator):
    bl_idname = "archipack.wall2_remove_finish"
    bl_label = "Remove finishings"
    bl_description = "Remove finishing"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'INTERNAL', 'UNDO'}
    index = IntProperty(default=0)

    @classmethod
    def poll(self, context):
        return archipack_wall2.poll(context.active_object)

    def execute(self, context):
        o = context.active_object
        d = archipack_wall2.datablock(o)
        d.finish.remove(self.index)
        d.update(context)
        return {'FINISHED'}


# ------------------------------------------------------------------
# Define operator class to draw a wall
# ------------------------------------------------------------------


class ARCHIPACK_OT_wall2_draw(ArchipackDrawTool, Operator):
    bl_idname = "archipack.wall2_draw"
    bl_label = "Draw a Wall"
    bl_description = "Create a wall by drawing its baseline in 3D view"
    bl_category = 'Archipack'

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
    sel = []
    act = None

    # constraint to other wall and make a T child
    parent = None
    takemat = None

    @classmethod
    def poll(cls, context):
        return True

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
        logger.debug("ARCHIPACK_OT_wall2_draw.sp_callback event %s %s state:%s", event.type, event.value, state)

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
                d = archipack_wall2.datablock(o)
                part = d.parts[0]
                part.length = delta.length
            else:
                o = self.o
                # select and make active
                self.select_object(context, o, True)
                d = archipack_wall2.datablock(o)
                # Check for end close to start and close when applicable
                dp = sp.placeloc - o.location
                if dp.length < 0.01:
                    d.closed = True
                    self.state = 'CANCEL'
                    return

                part = d.add_part(context, delta.length)

            # print("self.o :%s" % o.name)
            rM = o.matrix_world.inverted().to_3x3()
            g = d.get_generator()
            w = g.segs[-2]
            dp = rM * delta
            da = atan2(dp.y, dp.x) - w.straight(1).angle
            a0 = part.a0 + da
            if a0 > pi:
                a0 -= 2 * pi
            if a0 < -pi:
                a0 += 2 * pi
            part.a0 = a0
            d.update(context)

            self.select_object(context, old, True)
            self.flag_next = True
            context.area.tag_redraw()
            # print("feedback.on:%s" % self.feedback.on)

        self.state = state

    def sp_init(self, context, event, state, sp):
        # print("sp_init event %s %s %s" % (event.type, event.value, state))
        if state == 'SUCCESS':
            # point placed, check if a wall was under mouse
            res, tM, wall, width, y, z_offset = self.mouse_hover_wall(context, event)
            if res:
                d = archipack_wall2.datablock(wall)
                if event.ctrl:
                    # user snap, use direction as constraint
                    tM.translation = sp.placeloc.copy()
                else:
                    # without snap, use wall's bottom
                    tM.translation -= y.normalized() * (0.5 * d.width)
                self.takeloc = tM.translation
                self.parent = wall.name
                self.takemat = tM
            else:
                self.takeloc = sp.placeloc.copy()

            self.state = 'RUNNING'
            # print("feedback.on:%s" % self.feedback.on)
        elif state == 'CANCEL':
            self.state = state
            return

    def ensure_ccw(self):
        """
            Wall to slab expect wall vertex order to be ccw
            so reverse order here when needed
        """
        d = archipack_wall2.datablock(self.o)
        g = d.get_generator(axis=False)
        pts = [seg.p0.to_3d() for seg in g.segs]

        if d.closed:
            pts.append(pts[0])

        if d.is_cw(pts):
            d.x_offset = 1
            pts = list(reversed(pts))
            self.o.location += pts[0] - pts[-1]

        d.from_points(pts, d.closed)

    def modal(self, context, event):

        context.area.tag_redraw()
        if event.type in {'NONE', 'TIMER', 'EVT_TWEAK_L', 'WINDOW_DEACTIVATE'}:
            return {'PASS_THROUGH'}

        if self.keymap.check(event, self.keymap.delete):
            self.feedback.disable()
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            self.o = None
            return {'FINISHED', 'PASS_THROUGH'}

        if self.state == 'STARTING' and event.type not in {'ESC', 'RIGHTMOUSE'}:
            # wait for takeloc being visible when button is over horizon
            takeloc = self.mouse_to_plane(context, event)
            if takeloc is not None:
                logger.debug("ARCHIPACK_OT_wall2_draw.modal(STARTING) location:%s", takeloc)
                snap_point(takeloc=takeloc,
                    callback=self.sp_init,
                    constraint_axis=(True, True, False),
                    release_confirm=True)
            return {'RUNNING_MODAL'}

        elif self.state == 'RUNNING':
            # print("RUNNING")
            logger.debug("ARCHIPACK_OT_wall2_draw.modal(RUNNING) location:%s", self.takeloc)
            self.state = 'CREATE'
            snap_point(takeloc=self.takeloc,
                draw=self.sp_draw,
                takemat=self.takemat,
                transform_orientation=context.space_data.transform_orientation,
                callback=self.sp_callback,
                constraint_axis=(True, True, False),
                release_confirm=self.max_style_draw_tool)
            return {'RUNNING_MODAL'}

        elif self.state != 'CANCEL' and event.type in {'C', 'c'}:

            logger.debug("ARCHIPACK_OT_wall2_draw.modal(%s) C pressed", self.state)
            self.feedback.disable()
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')

            o = self.o
            # select and make active
            self.select_object(context, o, True)

            d = archipack_wall2.datablock(o)
            d.closed = True

            if bpy.ops.archipack.manipulate.poll():
                bpy.ops.archipack.manipulate('INVOKE_DEFAULT')

            return {'FINISHED'}

        elif self.state != 'CANCEL' and event.type in {'LEFTMOUSE', 'RET', 'NUMPAD_ENTER', 'SPACE'}:

            # print('LEFTMOUSE %s' % (event.value))
            self.feedback.instructions(context, "Draw a wall", "Click & Drag to add a segment", [
                ('ENTER', 'Add part'),
                ('BACK_SPACE', 'Remove part'),
                ('CTRL', 'Snap'),
                ('C', 'Close wall and exit'),
                ('MMBTN', 'Constraint to axis'),
                ('X Y', 'Constraint to axis'),
                ('RIGHTCLICK or ESC', 'exit')
                ])

            # press with max mode release with blender mode
            if self.max_style_draw_tool:
                evt_value = 'PRESS'
            else:
                evt_value = 'RELEASE'

            if event.value == evt_value:

                if self.flag_next:
                    self.flag_next = False
                    o = self.o

                    # select and make active
                    self.select_object(context, o, True)

                    d = archipack_wall2.datablock(o)
                    g = d.get_generator()
                    p0 = g.segs[-2].p0
                    p1 = g.segs[-2].p1
                    dp = p1 - p0
                    takemat = o.matrix_world * Matrix([
                        [dp.x, dp.y, 0, p1.x],
                        [dp.y, -dp.x, 0, p1.y],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]
                    ])
                    takeloc = o.matrix_world * p1.to_3d()
                    self.unselect_object(o)
                else:
                    takemat = None
                    takeloc = self.mouse_to_plane(context, event)

                if takeloc is not None:

                    logger.debug("ARCHIPACK_OT_wall2_draw.modal(CREATE) location:%s", takeloc)

                    snap_point(takeloc=takeloc,
                        takemat=takemat,
                        draw=self.sp_draw,
                        callback=self.sp_callback,
                        constraint_axis=(True, True, False),
                        release_confirm=self.max_style_draw_tool)

            return {'RUNNING_MODAL'}

        if self.keymap.check(event, self.keymap.undo) or (
                event.type in {'BACK_SPACE'} and event.value == 'RELEASE'
                ):
            if self.o is not None:
                o = self.o

                # select and make active
                self.select_object(context, o, True)
                d = archipack_wall2.datablock(o)
                if d.n_parts > 1:
                    d.n_parts -= 1
            return {'RUNNING_MODAL'}

        if self.state == 'CANCEL' or (event.type in {'ESC', 'RIGHTMOUSE'} and
                event.value == 'RELEASE'):

            self.feedback.disable()
            bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
            logger.debug("ARCHIPACK_OT_wall2_draw.modal(CANCEL) %s", event.type)
            if self.o is None:

                # select and make active
                self.select_object(context, self.act, True)

                for o in self.sel:
                    self.select_object(context, o)
            else:
                self.select_object(context, self.o, True)

                # remove last segment with blender mode
                d = archipack_wall2.datablock(self.o)
                if not self.max_style_draw_tool:
                    if not d.closed and d.n_parts > 1:
                        d.n_parts -= 1
                self.select_object(context, self.o, True)
                # make T child
                if self.parent is not None:
                    d.t_part = self.parent

                if bpy.ops.archipack.manipulate.poll():
                    bpy.ops.archipack.manipulate('INVOKE_DEFAULT')

            return {'FINISHED'}

        return {'PASS_THROUGH'}

    def invoke(self, context, event):

        if context.mode == "OBJECT":
            prefs = context.user_preferences.addons[__name__.split('.')[0]].preferences
            self.max_style_draw_tool = prefs.max_style_draw_tool
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
                ('SHIFT+CTRL+TAB', 'Switch snap mode'),
                ('RIGHTCLICK or ESC', 'exit without change')
                ])
            self.feedback.enable()
            args = (self, context)

            self.sel = context.selected_objects[:]
            self.act = context.active_object
            bpy.ops.object.select_all(action="DESELECT")
            
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


def register():
    bpy.utils.register_class(archipack_wall2_part)
    bpy.utils.register_class(archipack_wall2_child)
    bpy.utils.register_class(archipack_wall2_finish)
    bpy.utils.register_class(archipack_wall2_relocate_child)
    bpy.utils.register_class(archipack_wall2)
    Mesh.archipack_wall2 = CollectionProperty(type=archipack_wall2)
    bpy.utils.register_class(ARCHIPACK_PT_wall2)
    bpy.utils.register_class(ARCHIPACK_OT_wall2)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_add_finish)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_remove_finish)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_draw)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_from_slab)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_fit_roof)
    bpy.utils.register_class(ARCHIPACK_OT_wall2_to_curve)


def unregister():
    bpy.utils.unregister_class(archipack_wall2_part)
    bpy.utils.unregister_class(archipack_wall2_child)
    bpy.utils.unregister_class(archipack_wall2_finish)
    bpy.utils.unregister_class(archipack_wall2_relocate_child)
    bpy.utils.unregister_class(archipack_wall2)
    del Mesh.archipack_wall2
    bpy.utils.unregister_class(ARCHIPACK_PT_wall2)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_add_finish)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_remove_finish)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_draw)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_from_slab)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_fit_roof)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall2_to_curve)
