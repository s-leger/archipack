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
import time
# noinspection PyUnresolvedReferences
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty, BoolProperty, IntProperty,
    StringProperty, EnumProperty,
    CollectionProperty
    )
from .bmesh_utils import BmeshEdit as bmed
from random import randint
import bmesh
from .materialutils import MaterialUtils
from mathutils import Vector, Matrix
from mathutils.geometry import interpolate_bezier
from math import sin, cos, pi, atan2, sqrt, tan
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_2d import Line, Arc
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_object import ArchipackCreateTool, ArchipackObject


class Roof():

    def __init__(self):
        # total distance from start
        self.dist = 0
        self.t_start = 0
        self.t_end = 0
        self.dz = 0
        self.z0 = 0
        self.angle_0 = 0
        self.v0_idx = 0
        self.v1_idx = 0
        self.constraint_type = None
        self.slope_left = 1
        self.slope_right = 1
        self.width_left = 1
        self.width_right = 1
        self.depends_left = 0
        self.depends_right = 0
        # force hip or valley
        self.enforce_part = 'AUTO'

    def copy(self):
        # if "Straight" in type(self).__name__:
        s = StraightRoof(self.p.copy(), self.v.copy())
        # else:
        # s = CurvedRoof(self.c.copy(), self.radius, self.a0, self.da)
        s.angle_0 = self.angle_0
        s.z0 = self.z0
        s.v0_idx = self.v0_idx
        s.v1_idx = self.v1_idx
        s.constraint_type = self.constraint_type
        s.slope_left = self.slope_left
        s.slope_right = self.slope_right
        s.width_left = self.width_left
        s.width_right = self.width_right
        s.depends_left = self.depends_left
        s.depends_right = self.depends_right
        s.enforce_part = self.enforce_part
        return s

    def straight(self, length, t=1):
        s = self.copy()
        s.p = self.lerp(t)
        s.v = self.v.normalized() * length
        return s

    def set_offset(self, offset, last=None):
        """
            Offset line and compute intersection point
            between segments
        """
        self.line = self.make_offset(offset, last)

    @property
    def t_diff(self):
        return self.t_end - self.t_start

    def straight_roof(self, a0, length):
        s = self.straight(length).rotate(a0)
        r = StraightRoof(s.p, s.v)
        r.angle_0 = a0
        return r

    def curved_roof(self, a0, da, radius):
        n = self.normal(1).rotate(a0).scale(radius)
        if da < 0:
            n.v = -n.v
        c = n.p - n.v
        r = CurvedRoof(c, radius, n.angle, da)
        r.angle_0 = a0
        return r


class StraightRoof(Roof, Line):
    def __str__(self):
        return "t_start:{} t_end:{} dist:{}".format(self.t_start, self.t_end, self.dist)

    def __init__(self, p, v):
        Line.__init__(self, p, v)
        Roof.__init__(self)


class CurvedRoof(Roof, Arc):
    def __str__(self):
        return "t_start:{} t_end:{} dist:{}".format(self.t_start, self.t_end, self.dist)

    def __init__(self, c, radius, a0, da):
        Arc.__init__(self, c, radius, a0, da)
        Roof.__init__(self)


"""
class RoofGenerator():

    def __init__(self, parts, origin=Vector((0, 0))):
        self.parts = parts
        self.segs = []
        self.length = 0
        self.origin = origin
        self.user_defined_post = None
        self.user_defined_uvs = None
        self.user_defined_mat = None
        self.z = 0
        self.slope = 0

    def add_part(self, part):

        if len(self.segs) < 1:
            s = None
        else:
            s = self.segs[-1]

        # start a new roof
        if s is None:
            if part.type == 'S_SEG':
                v = part.length * Vector((cos(part.a0), sin(part.a0)))
                s = StraightRoof(self.origin, v)
            elif part.type == 'C_SEG':
                c = self.origin - part.radius * Vector((cos(part.a0), sin(part.a0)))
                s = CurvedRoof(c, part.radius, part.a0, part.da)
        else:
            if part.type == 'S_SEG':
                s = s.straight_roof(part.a0, part.length)
            elif part.type == 'C_SEG':
                s = s.curved_roof(part.a0, part.da, part.radius)

        self.segs.append(s)
        self.last_type = part.type

    def set_offset(self):
        last = None
        for i, seg in enumerate(self.segs):
            seg.set_offset(self.parts[i].offset, last)
            last = seg.line

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

            if len(self.segs) > 1:
                w.line = w.make_offset(self.parts[-1].offset, self.segs[-2].line)

            p1 = self.segs[0].line.p1
            self.segs[0].line = self.segs[0].make_offset(self.parts[0].offset, w.line)
            self.segs[0].line.p1 = p1

        self.length = 0
        for i, seg in enumerate(self.segs):
            seg.line.dist = self.length
            self.length += seg.line.length

    def locate_manipulators(self):


        for i, f in enumerate(self.segs):

            manipulators = self.parts[i].manipulators
            p0 = f.p0.to_3d()
            p1 = f.p1.to_3d()
            # angle from last to current segment
            if i > 0:

                if i < len(self.segs) - 1:
                    manipulators[0].type_key = 'ANGLE'
                else:
                    manipulators[0].type_key = 'DUMB_ANGLE'

                v0 = self.segs[i - 1].straight(-1, 1).v.to_3d()
                v1 = f.straight(1, 0).v.to_3d()
                manipulators[0].set_pts([p0, v0, v1])

            if type(f).__name__ == "StraightRoof":
                # segment length
                manipulators[1].type_key = 'SIZE'
                manipulators[1].prop1_name = "length"
                manipulators[1].set_pts([p0, p1, (1.0, 0, 0)])
            else:
                # segment radius + angle
                v0 = (f.p0 - f.c).to_3d()
                v1 = (f.p1 - f.c).to_3d()
                manipulators[1].type_key = 'ARC_ANGLE_RADIUS'
                manipulators[1].prop1_name = "da"
                manipulators[1].prop2_name = "radius"
                manipulators[1].set_pts([f.c.to_3d(), v0, v1])

            # snap manipulator, dont change index !
            manipulators[2].set_pts([p0, p1, (1, 0, 0)])
            # dumb segment id
            manipulators[3].set_pts([p0, p1, (1, 0, 0)])

            # offset
            manipulators[4].set_pts([
                p0,
                p0 + f.sized_normal(0, max(0.0001, self.parts[i].offset)).v.to_3d(),
                (0.5, 0, 0)
            ])

    def make_profile(self, profile, idmat,
            x_offset, z_offset, extend, verts, faces, matids, uvs):

        # self.set_offset(x_offset)

        n_roofs = len(self.segs) - 1

        if n_roofs < 0:
            return

        sections = []

        f = self.segs[0]

        # first step
        if extend != 0:
            t = -extend / self.segs[0].line.length
            n = f.line.sized_normal(t, 1)
            # n.p = f.lerp(x_offset)
            sections.append((n, f.dz / f.length, f.z0 + f.dz * t))

        # add first section
        n = f.line.sized_normal(0, 1)
        # n.p = f.lerp(x_offset)
        sections.append((n, f.dz / f.length, f.z0))

        for s, f in enumerate(self.segs):
            n = f.line.sized_normal(1, 1)
            # n.p = f.lerp(x_offset)
            sections.append((n, f.dz / f.length, f.z0 + f.dz))

        if extend != 0:
            t = 1 + extend / self.segs[-1].line.length
            n = f.line.sized_normal(t, 1)
            # n.p = f.lerp(x_offset)
            sections.append((n, f.dz / f.length, f.z0 + f.dz * t))

        user_path_verts = len(sections)
        f = len(verts)
        if user_path_verts > 0:
            user_path_uv_v = []
            n, dz, z0 = sections[-1]
            sections[-1] = (n, dz, z0)
            n_sections = user_path_verts - 1
            n, dz, zl = sections[0]
            p0 = n.p
            v0 = n.v.normalized()
            for s, section in enumerate(sections):
                n, dz, zl = section
                p1 = n.p
                if s < n_sections:
                    v1 = sections[s + 1][0].v.normalized()
                dir = (v0 + v1).normalized()
                scale = 1 / cos(0.5 * acos(min(1, max(-1, v0 * v1))))
                for p in profile:
                    x, y = n.p + scale * (x_offset + p.x) * dir
                    z = zl + p.y + z_offset
                    verts.append((x, y, z))
                if s > 0:
                    user_path_uv_v.append((p1 - p0).length)
                p0 = p1
                v0 = v1

            # build faces using Panel
            lofter = Lofter(
                # closed_shape, index, x, y, idmat
                True,
                [i for i in range(len(profile))],
                [p.x for p in profile],
                [p.y for p in profile],
                [idmat for i in range(len(profile))],
                closed_path=False,
                user_path_uv_v=user_path_uv_v,
                user_path_verts=user_path_verts
                )
            faces += lofter.faces(16, offset=f, path_type='USER_DEFINED')
            matids += lofter.mat(16, idmat, idmat, path_type='USER_DEFINED')
            v = Vector((0, 0))
            uvs += lofter.uv(16, v, v, v, v, 0, v, 0, 0, path_type='USER_DEFINED')

    def debug(self, verts):
        for s, roof in enumerate(self.segs):
            for i in range(33):
                x, y = roof.line.lerp(i / 32)
                verts.append((x, y, 0))

    def find_z(self, axis, pt):
        d = 100000
        z = 0
        n_axis = len(axis.segs) - 1
        for i, s in enumerate(axis.segs):
            if i < n_axis:
                res, d0, t = s.point_sur_segment(pt)
                if res and abs(d0) < d:
                    d = abs(d0)
                print("res:%s i:%s d0:%s t:%s d:%s" % (res, i, d0, t, d))
                z = axis.z - d * axis.slope
        return z

    def intersection_partition(self, array, begin, end):
        pivot = begin
        for i in range(begin + 1, end + 1):
            # param t
            if array[i][1] < array[begin][1]:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
        array[pivot], array[begin] = array[begin], array[pivot]
        return pivot

    def sort_intersection(self, array, begin=0, end=None):
        # print("sort_child")
        if end is None:
            end = len(array) - 1

        def _quicksort(array, begin, end):
            if begin >= end:
                return
            pivot = self.intersection_partition(array, begin, end)
            _quicksort(array, begin, pivot - 1)
            _quicksort(array, pivot + 1, end)
        return _quicksort(array, begin, end)


    def get_verts(self, verts, edges, axis):




        if axis is None:
            return
        n_axis = len(axis.segs) - 1

        # boundary vertices without intersections
        n_boundary = len(self.segs) - 1

        # boundary vertices with intersections
        axis_id = n_boundary + 1
        print("n_boundary:%s n_axis:%s axis_id:%s" % (n_boundary, n_axis, axis_id))

        # add boundary verts

        for boundary in self.segs:
            p = boundary.line.p0
            x, y = p
            z = self.find_z(axis, p)
            verts.append((x, y, z))

        # intersections axis / boundary
        axis_it = []

        intersections = []

        # intersections for each boundary part
        for i, boundary in enumerate(self.segs):
            it = []
            bl = boundary.line
            for j, ax in enumerate(axis.segs):
                if j < n_axis:
                    res, p, u, v = bl.intersect_ext(ax)
                    # intersect with axis
                    if res:
                        axis_it.append([j, v, p])
                        if v > 0.5:
                            it.append([j + 1, u, p, False])
                        else:
                            it.append([j, u, p, False])

                    if j > 0 and len(axis_it) < 2:
                        # right of axis
                        n = ax.sized_normal(1000000, 0)

                        # demi-angle pour pente constante
                        # n.rotate(0.5 * ax.delta_angle(axis.segs[j - 1]))

                        res, p, u, v = bl.intersect_ext(n)
                        if res:
                            it.append([j, u, p, True])
                            axis_id += 1
                        # left of axis
                        n.v = -n.v
                        res, p, u, v = bl.intersect_ext(n)
                        if res:
                            it.append([j, u, p, True])
                            axis_id += 1

            # sort intersection for this boundary by param t
            self.sort_intersection(it)
            intersections.append(it)

        for i, it in enumerate(intersections):
            # add verts and edges for this boundary part
            i0 = i
            for j, t, p, do_edge in it:
                x, y = p
                z = self.find_z(axis, p)
                if do_edge:
                    i1 = len(verts)
                    verts.append((x, y, z))
                    # add edge from point to axis
                    edges.append([i1, axis_id + j])
                else:
                    # use axis vert
                    i1 = axis_id + j

                # append edge between it
                edges.append([i0, i1])
                i0 = i1

            # add last edge for boundary
            i1 = i + 1
            if i1 > n_boundary:
                i1 = 0
            edges.append([i0, i1])


        # add axis verts and edges
        start, end = axis_it

        if end[0] < start[0]:
            end, start = start, end
        elif end[0] == start[0]:
            if end[1] < start[1]:
                end, start = start, end

        x, y = start[2]
        z = self.z
        verts.append((x, y, z))

        x, y = end[2]
        verts.append((x, y, z))

        i0 = axis_id

        for j, s in enumerate(axis.segs):
            if j > 0 and j < end[0]:
                x, y = s.p0
                z = self.z
                i1 = len(verts)
                verts.append((x, y, z))
                if j > 1:
                    edges.append([i0, i0 + 1])
                i0 += 1

        edges.append([i0, i0 + 1])
"""


class RoofAxisSegment():
    def __init__(self, a0, idx, reversed):
        self.a0 = a0
        self.idx = idx
        self.reversed = reversed


class RoofAxisNode():
    def __init__(self):
        # axis segments
        self.segs = []
        self.root = None
        self.center = -1

    @property
    def count(self):
        return len(self.segs)

    @property
    def idx(self):
        """
            index of segments in this node
        """
        return [n.idx for n in self.segs]

    @property
    def last_idx(self):
        """
            index of last segments in this node
        """
        return self.segs[-1].idx

    @property
    def last(self):
        """
            last segments in this node
        """
        return self.segs[-1]

    def left_idx(self, index):
        if index + 1 > self.count:
            return self.segs[0].idx
        return self.segs[index + 1].idx

    def right_idx(self, index):
        return self.segs[index - 1].idx

    def add(self, a0, idx, reversed):
        s = RoofAxisSegment(a0, idx, reversed)
        if reversed:
            self.root = s
        self.segs.append(s)

    # sort tree segments by angle
    def partition(self, array, begin, end):
        pivot = begin
        for i in range(begin + 1, end + 1):
            if array[i].a0 < array[begin].a0:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
        array[pivot], array[begin] = array[begin], array[pivot]
        return pivot

    def sort(self):
        def _quicksort(array, begin, end):
            if begin >= end:
                return
            pivot = self.partition(array, begin, end)
            _quicksort(array, begin, pivot - 1)
            _quicksort(array, pivot + 1, end)

        end = len(self.segs) - 1
        _quicksort(self.segs, 0, end)

        for i, s in enumerate(self.segs):
            if s == self.root:
                self.center = i
                return



"""
import bmesh

m = C.object.data
[(round(v.co.x, 2), round(v.co.y, 2), round(v.co.z, 2)) for v in m.vertices]
[tuple(p.vertices) for p in m.polygons]

uvs = []
bpy.ops.object.mode_set(mode='EDIT')
bm = bmesh.from_edit_mesh(m)

layer = bm.loops.layers.uv.verify()
for i, face in enumerate(bm.faces):
    uv = []
    for j, loop in enumerate(face.loops):
        co = loop[layer].uv
        uv.append((round(co.x, 2), round(co.y, 2)))
    uvs.append(uv)

uvs
"""

class RoofGenerator():

    def __init__(self, d, origin=Vector((0, 0, 0))):
        self.parts = d.parts
        self.segs = []
        self.length = 0
        self.origin = origin.to_2d()
        self.z = origin.z
        self.width_right = d.width_right
        self.width_left = d.width_left
        self.slope_left = d.slope_left
        self.slope_right = d.slope_right
        self.user_defined_tile = None
        self.user_defined_uvs = None
        self.user_defined_mat = None

    def add_part(self, part):

        if len(self.segs) < 1 or part.bound_idx < 1:
            s = None
        else:
            s = self.segs[part.bound_idx - 1]

        # start a new roof
        if s is None:
            if part.type == 'S_SEG':
                v = part.length * Vector((cos(part.a0), sin(part.a0)))
                s = StraightRoof(self.origin, v)
            elif part.type == 'C_SEG':
                c = self.origin - part.radius * Vector((cos(part.a0), sin(part.a0)))
                s = CurvedRoof(c, part.radius, part.a0, part.da)
        else:
            if part.type == 'S_SEG':
                s = s.straight_roof(part.a0, part.length)
            elif part.type == 'C_SEG':
                s = s.curved_roof(part.a0, part.da, part.radius)
        # parent segment (root) index is  v0_idx - 1
        s.v0_idx = part.bound_idx
        s.constraint_type = part.constraint_type

        if part.constraint_type == 'SLOPE':
            s.enforce_part = part.enforce_part
        else:
            s.enforce_part = 'AUTO'

        s.angle_0 = part.a0
        s.take_precedence = part.take_precedence
        self.segs.append(s)

    def locate_manipulators(self):
        """

        """
        for i, f in enumerate(self.segs):

            manipulators = self.parts[i].manipulators
            p0 = f.p0.to_3d()
            p1 = f.p1.to_3d()
            # angle from last to current segment
            if i > 0:

                manipulators[0].type_key = 'ANGLE'
                v0 = self.segs[f.v0_idx - 1].straight(-1, 1).v.to_3d()
                v1 = f.straight(1, 0).v.to_3d()
                manipulators[0].set_pts([p0, v0, v1])

            if type(f).__name__ == "StraightRoof":
                # segment length
                manipulators[1].type_key = 'SIZE'
                manipulators[1].prop1_name = "length"
                manipulators[1].set_pts([p0, p1, (1.0, 0, 0)])
            else:
                # segment radius + angle
                v0 = (f.p0 - f.c).to_3d()
                v1 = (f.p1 - f.c).to_3d()
                manipulators[1].type_key = 'ARC_ANGLE_RADIUS'
                manipulators[1].prop1_name = "da"
                manipulators[1].prop2_name = "radius"
                manipulators[1].set_pts([f.c.to_3d(), v0, v1])

            # dumb segment id
            manipulators[2].set_pts([p0, p1, (1, 0, 0)])

    def debug(self, verts):
        for s, roof in enumerate(self.segs):
            for i in range(33):
                x, y = roof.lerp(i / 32)
                verts.append((x, y, 0))

    # sort tree segments by angle
    def seg_partition(self, array, begin, end):
        pivot = begin
        for i in range(begin + 1, end + 1):
            # wall idx
            if array[i].a0 < array[begin].a0:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
        array[pivot], array[begin] = array[begin], array[pivot]
        return pivot

    def sort_seg(self, array, begin=0, end=None):
        # print("sort_child")
        if end is None:
            end = len(array) - 1

        def _quicksort(array, begin, end):
            if begin >= end:
                return
            pivot = self.seg_partition(array, begin, end)
            _quicksort(array, begin, pivot - 1)
            _quicksort(array, pivot + 1, end)
        return _quicksort(array, begin, end)

    def intersect_chevron(self, line, side1, side2, slope, height, alt, segs, verts):
        res, p, u, v = line.intersect_ext(segs[side1])
        if res:
            x, y = p
        else:
            res, p, u, v = line.intersect_ext(segs[side2])
            if res:
                x, y = p
            else:
                x, y = line.p1
                u = 1
        z = self.z - slope * u + alt
        verts.append((x, y, z))
        verts.append((x, y, z - height))

    def make_roof(self, verts, edges):
        """
            Init data structure for possibly multi branched nodes

        """
        x, y = self.segs[0].p0
        z = self.z
        verts.append((x, y, z))

        # Axis verts
        for idx, s in enumerate(self.segs):
            s.v1_idx = idx + 1
            x, y = s.p1
            verts.append((x, y, z))
            s.z0 = z

        # node are connected segments
        # node
        # (segment idx)
        # (angle from root part > 0 right)
        # (reversed) denote a seg connected by p1
        #            "root" of node
        nodes = [RoofAxisNode() for s in range(len(self.segs) + 1)]
        new_idx = len(self.segs)

        for i, s in enumerate(self.segs):
            nodes[s.v0_idx].add(s.angle_0, i, False)
            nodes[s.v1_idx].add(-pi, i, True)

        # override first segment a0
        # so regular sort does work
        # nodes[0].segs[0].a0 = -pi
        nodes[0].root = nodes[0].segs[0]

        # Propagate slope and width
        # on node basis along axis
        # contigous -> same
        # T: and (x % 2 == 1)
        # First one take precedence over others
        # others inherit from side
        #
        #         l / l    l = left
        #          3       r = right
        #   l _1_ /
        #   r     \
        #          2
        #          r\ l
        #
        # X: rigth one r left one l (x % 2 == 0)
        # inherits from side
        #
        #    l 3 l          l = left
        # l__1_|_2_l        r = right
        # r    |   r
        #    r 4 r
        #

        # Init slope and width on seg 0
        seg = self.segs[0]
        seg.slope_right = self.slope_right
        seg.slope_left = self.slope_left
        seg.width_right = self.width_right
        seg.width_left = self.width_left

        for idx, node in enumerate(nodes):

            node.sort()

            """
            # special case: segment 0
            # there is no "root"
            if node.root is None:
                for i, n in enumerate(node.segs):
                    seg = self.segs[n.idx]
                    if seg.v1_idx == 1:
                        node.segs = node.segs[i:] + node.segs[:i]
            """

            nb_segs = node.count

            if nb_segs < 2:
                continue

            # get "root" slope and width
            s = node.root.idx
            seg = self.segs[s]
            slope_right = seg.slope_right
            slope_left = seg.slope_left
            width_right = seg.width_right
            width_left = seg.width_left
            depends_left = seg.depends_left
            depends_right = seg.depends_right

            # simple case: 2 contigous segments
            if nb_segs == 2:
                s = node.last_idx
                seg = self.segs[s]
                seg.slope_right = slope_right
                seg.slope_left = slope_left
                seg.width_right = width_right
                seg.width_left = width_left
                seg.depends_left = depends_left
                seg.depends_right = depends_right
                continue

            # More than 2 segments, uneven distribution
            if nb_segs % 2 == 1:
                # find wich child does take precedence
                # first one on rootline (arbitrary)
                center = nb_segs
                for s in node.idx:
                    seg = self.segs[s]
                    if seg.v1_idx < center:
                        center = seg.v1_idx
            else:
                # even distribution
                center = nb_segs / 2

            # user defined precedence if any
            for i, s in enumerate(node.idx):
                seg = self.segs[s]
                if seg.take_precedence:
                    center = i
                    break

            for i, s in enumerate(node.idx):
                # repartir les dependences
                # chaque segment depend du precedant
                # TODO:
                # chaque segment adjascent d'un meme cote d'un axe
                # depend du precedant
                # prendre depuis l'axe dans les 2 directions
                seg = self.segs[s]
                if i > 0:
                    if i < center:
                        seg.slope_left = slope_right
                        seg.slope_right = slope_right
                        seg.width_left = width_right
                        seg.width_right = width_right
                        seg.depends_right = depends_right
                        seg.depends_left = depends_right
                    elif i == center:
                        seg.slope_left = slope_left
                        seg.slope_right = slope_right
                        seg.width_left = width_left
                        seg.width_right = width_right
                        seg.depends_right = depends_right
                        seg.depends_left = depends_left
                    else:
                        seg.slope_left = slope_left
                        seg.slope_right = slope_left
                        seg.width_left = width_left
                        seg.width_right = width_left
                        seg.depends_right = depends_left
                        seg.depends_left = depends_left

        # vertices for slope between sections
        #
        #    2 slope            2 slope           2 slope
        #     |                  |                 |
        #     |______section_1___|___section_2_____|
        #     |                  |                 |
        #     |                  |                 |
        #
        segs = []
        for idx, node in enumerate(nodes):

            # count for 'HORIZONTAL' segments in node
            nb_segs = node.count
            # if node is empty:
            if nb_segs < 1:
                # print("empty node")

                break

            root = self.segs[node.root.idx]

            # check if there is a user defined slope
            # disable auto-slope when found
            slope_left = None
            slope_right = None
            for n in node.segs:
                seg = self.segs[n.idx]
                if not n.reversed and seg.constraint_type == 'SLOPE':
                    # print("node:%s seg:%s angle_0:%s" % (idx, n.idx, seg.angle_0))
                    nb_segs -= 1
                    if seg.angle_0 > 0:
                        slope_left = n.idx
                    else:
                        slope_right = n.idx

            # one 'HORIZONTAL' segment in node
            # with possibly 'SLOPE' segments
            # on roof extremity
            if nb_segs < 2:

                if node.root.reversed:
                    # node on end of roof segment
                    t = 1
                else:
                    # node on start of roof segment
                    segs.append(root)
                    t = 0

                # Ensure root is 'HORIZONTAL'
                # dosent make sense otherwhise
                if root.constraint_type == 'HORIZONTAL':
                    # left - auto slope
                    if slope_left is None:
                        length = root.width_left
                        new_idx += 1
                        new = root.straight(length, t).rotate(pi / 2)
                        new.constraint_type = 'SLOPE'
                        new.slope_right = root.slope_left
                        new.width_right = root.width_left
                        new.angle_0 = pi / 2
                        new.v0_idx = idx
                        new.v1_idx = new_idx
                        segs.append(new)
                        x, y = new.p1
                        z = root.z0 - length * root.slope_left
                        verts.append((x, y, z))

                    # user defined left slope
                    else:
                        cur = self.segs[slope_left]
                        # project side using user defined slope
                        # to fit avaliable width
                        da = cur.delta_angle(root)
                        # angle always ccw
                        if da < 0:
                            da = 2 * pi + da
                        if da == 0:
                            width = root.width_left
                        else:
                            width = min(3 * root.width_left, root.width_left / sin(da))
                        # print("length:%s da:%s" % (length, da))
                        cur = cur.copy()
                        cur.v = cur.v.normalized() * width
                        segs.append(cur)
                        x, y = cur.p1
                        z = root.z0 - root.width_left * root.slope_left
                        verts[cur.v1_idx] = (x, y, z)

                    # right - auto slope
                    if slope_right is None:
                        length = root.width_right
                        new_idx += 1
                        new = root.straight(length, t).rotate(-pi / 2)
                        new.slope_left = root.slope_right
                        new.width_left = root.width_right
                        new.constraint_type = 'SLOPE'
                        new.angle_0 = -pi / 2
                        new.v0_idx = idx
                        new.v1_idx = new_idx
                        segs.append(new)
                        x, y = new.p1
                        z = root.z0 - length * root.slope_right
                        verts.append((x, y, z))

                    # user defined right slope
                    else:

                        cur = self.segs[slope_right]
                        # project side using user defined slope
                        # to fit avaliable width
                        da = cur.delta_angle(root)
                        # angle always ccw
                        if da < 0:
                            da = 2 * pi + da
                        if da == 0:
                            width = root.width_right
                        else:
                            width = min(3 * root.width_right, root.width_right / sin(da))

                        # print("user def slope right length:%s da:%s" % (width, da))

                        cur = cur.copy()
                        cur.v = cur.v.normalized() * -width
                        segs.append(cur)
                        x, y = cur.p1
                        z = root.z0 - root.width_right * root.slope_right
                        verts[cur.v1_idx] = (x, y, z)

            # between segments
            else:

                # add slope segs between horizontal parts

                # last segment in ccw order
                last = self.segs[node.last_idx]

                # take inverse of node when reversed
                # so angle mesure remains constant
                if node.last.reversed:
                    last = last.copy()
                    last.p += last.v
                    last.v = -last.v

                # new_node.append(last)
                for n in node.segs:
                    cur = self.segs[n.idx]
                    width = cur.width_right
                    slope = cur.slope_right

                    # take inverse of node when reversed
                    # so angle mesure remains constant
                    if n.reversed:
                        width = cur.width_left
                        slope = cur.slope_left
                        cur = cur.copy()
                        cur.p += cur.v
                        cur.v = -cur.v

                    if cur.constraint_type == 'HORIZONTAL':
                        if last.constraint_type == cur.constraint_type:

                            # add slope between 2 'HORIZONTAL'
                            # sections

                            # cur will inherits from reproj as width
                            # and slope

                            # wich side here ?

                            da = cur.delta_angle(last)
                            # angle always ccw
                            if da < 0:
                                da = 2 * pi + da

                            reproj = min(3 * width, width / sin(0.5 * da))
                            new_idx += 1
                            new = last.straight(reproj, 0).rotate(0.5 * da)

                            # maybe store width here ?
                            new.slope_left = slope
                            new.slope_right = slope
                            new.width_left = width
                            new.width_right = width

                            new.constraint_type = 'SLOPE'
                            new.v0_idx = idx
                            new.v1_idx = new_idx
                            new.angle_0 = new.delta_angle(root)
                            segs.append(new)
                            x, y = new.p1
                            z = root.z0 - width * slope
                            verts.append((x, y, z))

                    # user defined 'SLOPE' segment
                    # reproject width and z of vertex
                    else:
                        da = cur.delta_angle(root)
                        # angle always ccw
                        if da < 0:
                            da = 2 * pi + da
                        reproj = min(3 * width, width / sin(da))
                        # print("length:%s da:%s" % (length, da))
                        cur = cur.copy()
                        cur.v = cur.v.normalized() * reproj
                        if cur.angle_0 < 0:
                            cur.v = -cur.v
                        x, y = cur.p1
                        z = root.z0 - width * slope
                        verts[cur.v1_idx] = (x, y, z)

                        if not n.reversed:
                            # update "rootline"
                            self.segs[n.idx].slope_right = slope * width / reproj
                            self.segs[n.idx].width_right = reproj
                            cur.v1_idx

                    # do not append "root" segments here
                    if not n.reversed:
                        segs.append(cur)

                    last = cur

        # axis and slope edges
        for s in segs:
            edges.append([s.v0_idx, s.v1_idx])

        # Store nodes including 'SLOPE' ones
        nodes = [RoofAxisNode() for s in range(new_idx + 1)]

        for i, s in enumerate(segs):
            nodes[s.v0_idx].add(s.angle_0, i, False)
            nodes[s.v1_idx].add(-pi, i, True)

        # override first segment a0
        # so regular sort does work
        # nodes[0].segs[0].a0 = -pi
        nodes[0].root = nodes[0].segs[0]

        # sort segs in nodes in ccw order
        for node in nodes:
            node.sort()

        self.nodes = nodes
        self.all_segs = segs

    def bissect(self, bm,
            plane_co,
            plane_no,
            dist=0.001,
            use_snap_center=False,
            clear_outer=True,
            clear_inner=False
            ):
        geom = bm.verts[:]
        geom.extend(bm.edges[:])
        geom.extend(bm.faces[:])

        bmesh.ops.bisect_plane(bm,
            geom=geom,
            dist=dist,
            plane_co=plane_co,
            plane_no=plane_no,
            use_snap_center=False,
            clear_outer=clear_outer,
            clear_inner=clear_inner
            )

    def lambris(self, d, verts, faces, edges, matids, uvs):

        idmat = 0
        segs = self.all_segs
        nodes = self.nodes

        for idx, node in enumerate(nodes):

            # find next node in segment
            for i, s in enumerate(node.segs):

                seg = segs[s.idx]

                if s.reversed:
                    continue

                next = nodes[seg.v1_idx]

                if next.count > 1:

                    center_0 = i
                    center_1 = next.center

                    # found 2nd node
                    if center_1 > -1:

                        i0 = seg.v0_idx
                        i1 = seg.v1_idx

                        ############
                        # right
                        ############

                        ir0 = next.right_idx(center_1)
                        ir1 = node.left_idx(center_0)

                        vi2 = segs[ir0].v1_idx
                        vi3 = segs[ir1].v1_idx

                        edges.append([vi2, vi3])
                        faces.append([i1, i0, vi3, vi2])

                        ############
                        # left
                        ############

                        il0 = node.right_idx(center_0)
                        il1 = next.left_idx(center_1)

                        vi2 = segs[il0].v1_idx
                        vi3 = segs[il1].v1_idx

                        edges.append([vi2, vi3])
                        faces.append([i0, i1, vi3, vi2])

                        matids.extend([idmat, idmat])
                        uvs.extend([
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)]
                        ])

    def couverture(self, context, o, d):

        idmat = 7
        rand = 3

        segs = self.all_segs
        nodes = self.nodes

        left, right = True, True
        angle_90 = round(pi / 2, 4)


        sx, sy, sz = d.tile_size_x, d.tile_size_y, d.tile_size_z

        """
        /* Bevel offset_type slot values */
        enum {
          BEVEL_AMT_OFFSET,
          BEVEL_AMT_WIDTH,
          BEVEL_AMT_DEPTH,
          BEVEL_AMT_PERCENT
        };
        """
        offset_type = 3

        if d.tile_offset > 0:
            offset = - d.tile_offset / 100
        else:
            offset = 0

        if d.tile_model == 'BRAAS2':
            t_pts = [Vector(p) for p in [
                (0.06, -1.0, 1.0), (0.19, -1.0, 0.5), (0.31, -1.0, 0.5), (0.44, -1.0, 1.0),
                (0.56, -1.0, 1.0), (0.69, -1.0, 0.5), (0.81, -1.0, 0.5), (0.94, -1.0, 1.0),
                (0.06, 0.0, 0.5), (0.19, 0.0, 0.0), (0.31, 0.0, 0.0), (0.44, 0.0, 0.5),
                (0.56, 0.0, 0.5), (0.69, 0.0, 0.0), (0.81, 0.0, 0.0), (0.94, 0.0, 0.5),
                (-0.0, -1.0, 1.0), (-0.0, 0.0, 0.5), (1.0, -1.0, 1.0), (1.0, 0.0, 0.5)]]
            t_faces = [
                (16, 0, 8, 17), (0, 1, 9, 8), (1, 2, 10, 9), (2, 3, 11, 10),
                (3, 4, 12, 11), (4, 5, 13, 12), (5, 6, 14, 13), (6, 7, 15, 14), (7, 18, 19, 15)]
        elif d.tile_model == 'BRAAS1':
            t_pts = [Vector(p) for p in [
                (0.1, -1.0, 1.0), (0.2, -1.0, 0.5), (0.6, -1.0, 0.5), (0.7, -1.0, 1.0),
                (0.1, 0.0, 0.5), (0.2, 0.0, 0.0), (0.6, 0.0, 0.0), (0.7, 0.0, 0.5),
                (-0.0, -1.0, 1.0), (-0.0, 0.0, 0.5), (1.0, -1.0, 1.0), (1.0, 0.0, 0.5)]]
            t_faces = [(8, 0, 4, 9), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 10, 11, 7)]
        elif d.tile_model == 'ETERNIT':
            t_pts = [Vector(p) for p in [
                (0.11, -1.0, 1.0), (0.9, -1.0, 1.0), (0.0, -0.79, 0.79),
                (1.0, -0.79, 0.79), (0.0, 2.0, -2.0), (1.0, 2.0, -2.0)]]
            t_faces = [(0, 1, 3, 5, 4, 2)]
        elif d.tile_model == 'ONDULEE':
            t_pts = [Vector(p) for p in [
                (0.0, -1.0, 0.1), (0.05, -1.0, 1.0), (0.1, -1.0, 0.1),
                (0.15, -1.0, 1.0), (0.2, -1.0, 0.1), (0.25, -1.0, 1.0),
                (0.3, -1.0, 0.1), (0.35, -1.0, 1.0), (0.4, -1.0, 0.1),
                (0.45, -1.0, 1.0), (0.5, -1.0, 0.1), (0.55, -1.0, 1.0),
                (0.6, -1.0, 0.1), (0.65, -1.0, 1.0), (0.7, -1.0, 0.1),
                (0.75, -1.0, 1.0), (0.8, -1.0, 0.1), (0.85, -1.0, 1.0),
                (0.9, -1.0, 0.1), (0.95, -1.0, 1.0), (1.0, -1.0, 0.1),
                (0.0, 0.0, 0.0), (0.05, 0.0, 0.9), (0.1, 0.0, 0.0),
                (0.15, 0.0, 0.9), (0.2, 0.0, 0.0), (0.25, 0.0, 0.9),
                (0.3, 0.0, 0.0), (0.35, 0.0, 0.9), (0.4, 0.0, 0.0),
                (0.45, 0.0, 0.9), (0.5, 0.0, 0.0), (0.55, 0.0, 0.9),
                (0.6, 0.0, 0.0), (0.65, 0.0, 0.9), (0.7, 0.0, 0.0),
                (0.75, 0.0, 0.9), (0.8, 0.0, 0.0), (0.85, 0.0, 0.9),
                (0.9, 0.0, 0.0), (0.95, 0.0, 0.9), (1.0, 0.0, 0.0)]]
            t_faces = [
                (0, 1, 22, 21), (1, 2, 23, 22), (2, 3, 24, 23),
                (3, 4, 25, 24), (4, 5, 26, 25), (5, 6, 27, 26),
                (6, 7, 28, 27), (7, 8, 29, 28), (8, 9, 30, 29),
                (9, 10, 31, 30), (10, 11, 32, 31), (11, 12, 33, 32),
                (12, 13, 34, 33), (13, 14, 35, 34), (14, 15, 36, 35),
                (15, 16, 37, 36), (16, 17, 38, 37), (17, 18, 39, 38),
                (18, 19, 40, 39), (19, 20, 41, 40)]
        elif d.tile_model == 'METAL':
            t_pts = [Vector(p) for p in [
                (0.0, -1.0, 0.0), (0.99, -1.0, 0.0), (1.0, -1.0, 0.0),
                (0.0, 0.0, 0.0), (0.99, 0.0, 0.0), (1.0, 0.0, 0.0),
                (0.99, -1.0, 1.0), (1.0, -1.0, 1.0), (1.0, 0.0, 1.0), (0.99, 0.0, 1.0)]]
            t_faces = [(0, 1, 4, 3), (7, 2, 5, 8), (1, 6, 9, 4), (6, 7, 8, 9)]
        elif d.tile_model == 'LAUZE':
            t_pts = [Vector(p) for p in [
                (0.75, -0.8, 0.8), (0.5, -1.0, 1.0), (0.25, -0.8, 0.8),
                (0.0, -0.5, 0.5), (1.0, -0.5, 0.5), (0.0, 0.5, -0.5), (1.0, 0.5, -0.5)]]
            t_faces = [(1, 0, 4, 6, 5, 3, 2)]
        elif d.tile_model == 'PLACEHOLDER':
            t_pts = [Vector(p) for p in [(0.0, -1.0, 1.0), (1.0, -1.0, 1.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]]
            t_faces = [(0, 1, 3, 2)]
        elif d.tile_model == 'ROMAN':
            t_pts = [Vector(p) for p in [
                (0.18, 0.0, 0.3), (0.24, 0.0, 0.58), (0.76, 0.0, 0.58), 
                (0.82, 0.0, 0.3), (0.05, -1.0, 0.5), (0.14, -1.0, 0.8), 
                (0.86, -1.0, 0.8), (0.95, -1.0, 0.5), (0.45, 0.0, 0.5), 
                (0.36, 0.0, 0.2), (-0.36, 0.0, 0.2), (-0.45, -0.0, 0.5), 
                (0.32, -1.0, 0.7), (0.26, -1.0, 0.42), (-0.26, -1.0, 0.42), 
                (-0.32, -1.0, 0.7), (0.5, 0.0, 0.74), (0.5, -1.0, 1.0), 
                (-0.0, -1.0, 0.26), (-0.0, 0.0, 0.0)]
            ]
            t_faces = [
                (0, 4, 5, 1), (16, 17, 6, 2), (2, 6, 7, 3), 
                (13, 12, 8, 9), (18, 13, 9, 19), (15, 14, 10, 11), 
                (14, 18, 19, 10), (1, 5, 17, 16)
            ]
        elif d.tile_model == 'ROUND':
            t_pts = [Vector(p) for p in [
                (0.0, -0.5, 0.5), (1.0, -0.5, 0.5), (0.0, 0.0, 0.0), 
                (1.0, 0.0, 0.0), (0.93, -0.71, 0.71), (0.78, -0.88, 0.88), 
                (0.39, -0.97, 0.97), (0.61, -0.97, 0.97), (0.07, -0.71, 0.71), 
                (0.22, -0.88, 0.88)]
            ]
            t_faces = [(6, 7, 5, 4, 1, 3, 2, 0, 8, 9)]
        else:
            return
        
        n_faces = len(t_faces)
        t_uvs = [[(t_pts[i].x, t_pts[i].y) for i in f] for f in t_faces]

        for idx, node in enumerate(nodes):
            # segment always between 2 nodes
            # create edge between rightmost of first node to leftmost of next node
            # same for other side

            # find next node in segment
            for i, s in enumerate(node.segs):

                seg = segs[s.idx]

                if s.reversed:
                    continue

                if nodes[seg.v1_idx].count > 1:

                    next = nodes[seg.v1_idx]
                    # segments sorted by angle from axis
                    # center are ids of axis segment on each node
                    # so center + 1 and center - 1
                    # are leftmost and rightmost slope segments

                    n_next = next.count
                    n_node = node.count

                    center_0 = i
                    center_1 = next.center

                    # found 2nd node
                    if center_1 > -1:

                        ############
                        # right
                        ############
                        if right:
                            dx, dy = d.tile_space_x, d.tile_space_y

                            ir0 = next.right_idx(center_1)
                            ir1 = node.left_idx(center_0)

                            s0 = segs[ir0]
                            s1 = segs[ir1]

                            # right part is larger than axis: compute t param in axis
                            res, d0, u = seg.point_sur_segment(s1.p1)
                            res, dr, v = seg.point_sur_segment(s0.p1)
                            trmin = min(0, u)
                            trmax = max(1, v)

                            #####################
                            # Right part
                            #####################

                            # compute base matrix top left of face
                            vx = -seg.v.normalized().to_3d()
                            vy = Vector((-vx.y, vx.x, seg.slope_left)).normalized()
                            vz = vx.cross(vy)

                            x0, y0 = seg.lerp(trmax)
                            z0 = self.z + d.tile_altitude

                            space_x = (trmax - trmin) * seg.length + 2 * d.tile_side
                            space_y = (d.tile_border + abs(dr)) * sqrt(1 + seg.slope_left * seg.slope_left)
                            n_x = 1 + int(space_x / dx)
                            n_y = 1 + int(space_y / dy)

                            if d.tile_fit_x:
                                dx = space_x / n_x

                            if d.tile_fit_y:
                                dy = space_y / n_y

                            if d.tile_alternate:
                                n_y += 1

                            tM = Matrix([
                                [vx.x, vy.x, vz.x, x0],
                                [vx.y, vy.y, vz.y, y0],
                                [vx.z, vy.z, vz.z, z0],
                                [0, 0, 0, 1]
                            ])

                            verts = []
                            faces = []
                            matids = []
                            uvs = []

                            for k in range(n_y):

                                y = k * dy

                                x0 = offset * dx - d.tile_side
                                nx = n_x

                                if d.tile_alternate and k % 2 == 1:
                                    x0 -= 0.5 * dx
                                    nx += 1

                                if d.tile_offset > 0:
                                    nx += 1

                                for j in range(nx):
                                    x = x0 + j * dx
                                    lM = tM * Matrix([
                                        [sx, 0, 0, x],
                                        [0, sy, 0, -y],
                                        [0, 0, sz, 0],
                                        [0, 0, 0, 1]
                                    ])

                                    v = len(verts)

                                    verts.extend([lM * p for p in t_pts])
                                    faces.extend([tuple(i + v for i in f) for f in t_faces])
                                    id = randint(idmat, idmat + rand)
                                    t_mats = [id for i in range(n_faces)]
                                    matids.extend(t_mats)
                                    uvs.extend(t_uvs)

                            # build temp bmesh and bissect
                            bm = bmed.buildmesh(
                                context, o, verts, faces, matids=matids, uvs=uvs,
                                weld=False, clean=False, auto_smooth=True, temporary=True)

                            # len(node/next) > 3  -> couloir ou faitiere
                            # len(node/next) < 3  -> terminaison

                            da0 = round(abs(seg.delta_angle(s1)), 4)
                            da1 = round(abs(seg.delta_angle(s0)), 4)

                            b0 = seg.p0.copy()
                            b1 = seg.p1.copy()
                            n0 = -s1.cross_z
                            n1 = s0.cross_z

                            if n_node > 3:
                                # angle from root > 90 = faitiere
                                # angle from root < 90 = couloir
                                if da0 < angle_90:
                                    # couloir
                                    b0 -= n0.normalized() * d.tile_couloir
                            else:
                                # bord
                                b0 += n0.normalized() * d.tile_side

                            if n_next > 3:

                                if da1 > angle_90:
                                    # couloir
                                    b1 -= n1.normalized() * d.tile_couloir
                            else:
                                # bord
                                b1 += n1.normalized() * d.tile_side

                            # -> couloir decouper vers l'interieur
                            # -> terminaison rallonger ou raccourcir
                            # -> faitiere: laisser tel quel
                            # origin and normal to bissect mesh

                            # node
                            self.bissect(bm, b0.to_3d(), n0.to_3d())
                            # next
                            self.bissect(bm, b1.to_3d(), n1.to_3d())
                            # top
                            self.bissect(bm, seg.p1.to_3d(), seg.cross_z.to_3d())
                            # bottom
                            x, y, z = seg.p0.to_3d() + space_y * -vy
                            self.bissect(bm, Vector((x, y, self.z + z)), -vy)

                            if d.tile_bevel:
                                geom = bm.verts[:]
                                geom.extend(bm.edges[:])
                                bmesh.ops.bevel(bm,
                                    geom=geom,
                                    offset=d.tile_bevel_amt,
                                    offset_type=offset_type,
                                    segments=d.tile_bevel_segs,
                                    profile=0.5,
                                    vertex_only=False,
                                    clamp_overlap=True,
                                    material=-1)

                            if d.tile_solidify:
                                geom = bm.faces[:]
                                verts = bm.verts[:]
                                bmesh.ops.solidify(bm, geom=geom, thickness=0.0001)
                                bmesh.ops.translate(bm, vec=vz * d.tile_height, space=o.matrix_world, verts=verts)

                            # merge with object
                            bmed.bmesh_join(context, o, [bm], normal_update=True)

                        ###################
                        # left
                        ###################
                        if left:

                            dx, dy = d.tile_space_x, d.tile_space_y

                            il0 = node.right_idx(center_0)
                            il1 = next.left_idx(center_1)

                            s0 = segs[il0]  # sur node
                            s1 = segs[il1]  # sur next

                            # left part is larger than axis: compute t param in axis
                            res, d0, u = seg.point_sur_segment(s0.p1)
                            res, dl, v = seg.point_sur_segment(s1.p1)
                            tlmin = min(0, u)
                            tlmax = max(1, v)

                            # compute matrix for face
                            vx = seg.v.normalized().to_3d()
                            vy = Vector((-vx.y, vx.x, seg.slope_right)).normalized()
                            vz = vx.cross(vy)

                            x0, y0 = seg.lerp(tlmin)
                            z0 = self.z + d.tile_altitude

                            space_x = (tlmax - tlmin) * seg.length + 2 * d.tile_side
                            space_y = (d.tile_border + abs(dl)) * sqrt(1 + seg.slope_right * seg.slope_right)
                            n_x = 1 + int(space_x / dx)
                            n_y = 1 + int(space_y / dy)


                            if d.tile_fit_x:
                                dx = space_x / n_x

                            if d.tile_fit_y:
                                dy = space_y / n_y

                            if d.tile_alternate:
                                n_y += 1

                            tM = Matrix([
                                [vx.x, vy.x, vz.x, x0],
                                [vx.y, vy.y, vz.y, y0],
                                [vx.z, vy.z, vz.z, z0],
                                [0, 0, 0, 1]
                            ])

                            verts = []
                            faces = []
                            matids = []
                            uvs = []

                            for k in range(n_y):

                                y = k * dy

                                x0 = offset * dx - d.tile_side
                                nx = n_x

                                if d.tile_alternate and k % 2 == 1:
                                    x0 -= 0.5 * dx
                                    nx += 1

                                if d.tile_offset > 0:
                                    nx += 1

                                for j in range(nx):
                                    x = x0 + j * dx
                                    lM = tM * Matrix([
                                        [sx, 0, 0, x],
                                        [0, sy, 0, -y],
                                        [0, 0, sz, 0],
                                        [0, 0, 0, 1]
                                    ])

                                    v = len(verts)

                                    verts.extend([lM * p for p in t_pts])
                                    faces.extend([tuple(i + v for i in f) for f in t_faces])
                                    id = randint(idmat, idmat + rand)
                                    t_mats = [id for i in range(n_faces)]
                                    matids.extend(t_mats)
                                    uvs.extend(t_uvs)

                            # build temp bmesh and bissect
                            bm = bmed.buildmesh(
                                context, o, verts, faces, matids=matids, uvs=uvs,
                                weld=False, clean=False, auto_smooth=True, temporary=True)

                            da0 = round(abs(seg.delta_angle(s0)), 4)
                            da1 = round(abs(seg.delta_angle(s1)), 4)

                            b0 = seg.p0.copy()
                            b1 = seg.p1.copy()
                            n0 = s0.cross_z
                            n1 = -s1.cross_z

                            if n_node > 3:
                                if da0 < angle_90:
                                    # couloir
                                    b0 -= n0.normalized() * d.tile_couloir
                            else:
                                # bord
                                b0 += n0.normalized() * d.tile_side

                            if n_next > 3:
                                # angle from root > 90 = faitiere
                                # angle from root < 90 = couloir
                                if da1 > angle_90:
                                    # couloir
                                    b1 -= n1.normalized() * d.tile_couloir
                            else:
                                # bord
                                b1 += n1.normalized() * d.tile_side
                            # node
                            self.bissect(bm, b0.to_3d(), n0.to_3d())
                            # next
                            self.bissect(bm, b1.to_3d(), n1.to_3d())
                            # top
                            self.bissect(bm, seg.p1.to_3d(), -seg.cross_z.to_3d())
                            # bottom
                            # w = seg.width_right + d.tile_border
                            # z = self.z + d.tile_altitude + d.tile_size_z - w * seg.slope_right
                            # x, y = w * -vy.to_2d().normalized()
                            x, y, z = seg.p0.to_3d() + space_y * -vy
                            self.bissect(bm, Vector((x, y, self.z + z)), -vy)

                            if d.tile_bevel:
                                geom = bm.verts[:]
                                geom.extend(bm.edges[:])
                                bmesh.ops.bevel(bm,
                                    geom=geom,
                                    offset=d.tile_bevel_amt,
                                    offset_type=offset_type,
                                    segments=d.tile_bevel_segs,
                                    profile=0.5,
                                    vertex_only=False,
                                    clamp_overlap=True,
                                    material=-1)

                            if d.tile_solidify:
                                geom = bm.faces[:]
                                verts = bm.verts[:]
                                bmesh.ops.solidify(bm, geom=geom, thickness=0.0001)
                                bmesh.ops.translate(bm, vec=vz * d.tile_height, space=o.matrix_world, verts=verts)

                            # merge with object
                            bmed.bmesh_join(context, o, [bm], normal_update=True)

        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.faces_shade_smooth()
        bpy.ops.object.mode_set(mode='OBJECT')

    def rake(self, d, verts, faces, edges, matids, uvs):

        idmat = 1

        segs = self.all_segs
        nodes = self.nodes

        for idx, node in enumerate(nodes):

            # find next node in segment
            for i, s in enumerate(node.segs):

                seg = segs[s.idx]

                if s.reversed:
                    continue

                next = nodes[seg.v1_idx]

                if next.count > 1:

                    center_0 = i
                    center_1 = next.center

                    # found 2nd node
                    if center_1 > -1:

                        ############
                        # right
                        ############

                        ir0 = next.right_idx(center_1)
                        ir1 = node.left_idx(center_0)

                        ############
                        # left
                        ############

                        il0 = node.right_idx(center_0)
                        il1 = next.left_idx(center_1)

                        #####################
                        # Vire-vents
                        #####################
                        if node.count == 3:

                            # right
                            s0 = segs[il0]
                            if s0.enforce_part is 'AUTO':
                                f = len(verts)
                                s1 = s0.offset(d.rake_offset)
                                res, d0, t = seg.point_sur_segment(s0.p1)
                                slope = abs(d0) * seg.slope_right
                                s0 = s0.offset(d.rake_offset - d.rake_width)
                                x0, y0 = s0.p0
                                x1, y1 = s1.p0
                                x2, y2 = s1.p1
                                x3, y3 = s0.p1
                                z0 = self.z + d.rake_altitude
                                z1 = z0 - slope
                                verts.extend([
                                    (x0, y0, z0),
                                    (x1, y1, z0),
                                    (x2, y2, z1),
                                    (x3, y3, z1),
                                ])
                                z0 -= d.rake_height
                                z1 -= d.rake_height
                                verts.extend([
                                    (x0, y0, z0),
                                    (x1, y1, z0),
                                    (x2, y2, z1),
                                    (x3, y3, z1),
                                ])

                                faces.extend([
                                    # top
                                    (f, f + 1, f + 2, f + 3),
                                    # sides
                                    (f, f + 4, f + 5, f + 1),
                                    (f + 1, f + 5, f + 6, f + 2),
                                    (f + 2, f + 6, f + 7, f + 3),
                                    (f + 3, f + 7, f + 4, f),
                                    # bottom
                                    (f + 4, f + 7, f + 6, f + 5)
                                ])
                                edges.append([f, f + 3])
                                edges.append([f + 1, f + 2])
                                edges.append([f + 4, f + 7])
                                edges.append([f + 5, f + 6])

                                matids.extend([idmat, idmat, idmat, idmat, idmat, idmat])
                                uvs.extend([
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)]
                                ])

                            # left
                            s1 = segs[ir1]
                            if s1.enforce_part is 'AUTO':
                                f = len(verts)
                                s0 = s1.offset(-d.rake_offset)
                                res, d0, t = seg.point_sur_segment(s1.p1)
                                slope = abs(d0) * seg.slope_left
                                s1 = s1.offset(d.rake_width - d.rake_offset )
                                x0, y0 = s0.p0
                                x1, y1 = s1.p0
                                x2, y2 = s1.p1
                                x3, y3 = s0.p1

                                z0 = self.z + d.rake_altitude
                                z1 = z0 - slope
                                verts.extend([
                                    (x0, y0, z0),
                                    (x1, y1, z0),
                                    (x2, y2, z1),
                                    (x3, y3, z1),
                                ])
                                z0 -= d.rake_height
                                z1 -= d.rake_height
                                verts.extend([
                                    (x0, y0, z0),
                                    (x1, y1, z0),
                                    (x2, y2, z1),
                                    (x3, y3, z1),
                                ])

                                faces.extend([
                                    # top
                                    (f, f + 1, f + 2, f + 3),
                                    # sides
                                    (f, f + 4, f + 5, f + 1),
                                    (f + 1, f + 5, f + 6, f + 2),
                                    (f + 2, f + 6, f + 7, f + 3),
                                    (f + 3, f + 7, f + 4, f),
                                    # bottom
                                    (f + 4, f + 7, f + 6, f + 5)
                                ])
                                edges.append([f, f + 3])
                                edges.append([f + 1, f + 2])
                                edges.append([f + 4, f + 7])
                                edges.append([f + 5, f + 6])
                                matids.extend([idmat, idmat, idmat, idmat, idmat, idmat])
                                uvs.extend([
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)]
                                ])

                        if next.count == 3:
                            # left
                            s0 = segs[ir0]
                            if s0.enforce_part is 'AUTO':
                                f = len(verts)
                                s1 = s0.offset(d.rake_offset )
                                res, d0, t = seg.point_sur_segment(s0.p1)
                                slope = abs(d0) * seg.slope_left
                                s0 = s0.offset(d.rake_offset - d.rake_width)
                                x0, y0 = s0.p0
                                x1, y1 = s1.p0
                                x2, y2 = s1.p1
                                x3, y3 = s0.p1

                                z0 = self.z + d.rake_altitude
                                z1 = z0 - slope
                                verts.extend([
                                    (x0, y0, z0),
                                    (x1, y1, z0),
                                    (x2, y2, z1),
                                    (x3, y3, z1),
                                ])
                                z0 -= d.rake_height
                                z1 -= d.rake_height
                                verts.extend([
                                    (x0, y0, z0),
                                    (x1, y1, z0),
                                    (x2, y2, z1),
                                    (x3, y3, z1),
                                ])

                                faces.extend([
                                    # top
                                    (f, f + 1, f + 2, f + 3),
                                    # sides
                                    (f, f + 4, f + 5, f + 1),
                                    (f + 1, f + 5, f + 6, f + 2),
                                    (f + 2, f + 6, f + 7, f + 3),
                                    (f + 3, f + 7, f + 4, f),
                                    # bottom
                                    (f + 4, f + 7, f + 6, f + 5)
                                ])
                                edges.append([f, f + 3])
                                edges.append([f + 1, f + 2])
                                edges.append([f + 4, f + 7])
                                edges.append([f + 5, f + 6])
                                matids.extend([idmat, idmat, idmat, idmat, idmat, idmat])
                                uvs.extend([
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)]
                                ])
                            # Right
                            s1 = segs[il1]
                            if s1.enforce_part is 'AUTO':
                                f = len(verts)

                                s0 = s1.offset(-d.rake_offset )
                                res, d0, t = seg.point_sur_segment(s0.p1)
                                slope = abs(d0) * seg.slope_right
                                s1 = s1.offset(d.rake_width - d.rake_offset)
                                x0, y0 = s0.p0
                                x1, y1 = s1.p0
                                x2, y2 = s1.p1
                                x3, y3 = s0.p1

                                z0 = self.z + d.rake_altitude
                                z1 = z0 - slope
                                verts.extend([
                                    (x0, y0, z0),
                                    (x1, y1, z0),
                                    (x2, y2, z1),
                                    (x3, y3, z1),
                                ])
                                z0 -= d.rake_height
                                z1 -= d.rake_height
                                verts.extend([
                                    (x0, y0, z0),
                                    (x1, y1, z0),
                                    (x2, y2, z1),
                                    (x3, y3, z1),
                                ])

                                faces.extend([
                                    # top
                                    (f, f + 1, f + 2, f + 3),
                                    # sides
                                    (f, f + 4, f + 5, f + 1),
                                    (f + 1, f + 5, f + 6, f + 2),
                                    (f + 2, f + 6, f + 7, f + 3),
                                    (f + 3, f + 7, f + 4, f),
                                    # bottom
                                    (f + 4, f + 7, f + 6, f + 5)
                                ])
                                edges.append([f, f + 3])
                                edges.append([f + 1, f + 2])
                                edges.append([f + 4, f + 7])
                                edges.append([f + 5, f + 6])
                                matids.extend([idmat, idmat, idmat, idmat, idmat, idmat])
                                uvs.extend([
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)]
                                ])

    def facia(self, d, verts, faces, edges, matids, uvs):

        idmat = 2

        segs = self.all_segs
        nodes = self.nodes

        for idx, node in enumerate(nodes):

            for i, s in enumerate(node.segs):

                seg = segs[s.idx]

                if s.reversed:
                    continue

                if nodes[seg.v1_idx].count > 1:
                    next = nodes[seg.v1_idx]
                    # segments sorted by angle from axis
                    # center are ids of axis segment on each node
                    # so center + 1 and center - 1
                    # are leftmost and rightmost slope segments
                    center_0 = i
                    center_1 = next.center

                    # found 2nd node
                    if center_1 > -1:

                        ############
                        # right
                        ############

                        ir0 = next.right_idx(center_1)
                        ir1 = node.left_idx(center_0)

                        ############
                        # left
                        ############

                        il0 = node.right_idx(center_0)
                        il1 = next.left_idx(center_1)

                        #####################
                        # Larmiers
                        #####################

                        # ! pas forcement perpendiculaire
                        # t  x/sin(a)
                        # angle axe / segment
                        """
                         1___________________2
                        0|___________________|3
                        """
                        f = len(verts)
                        s0 = segs[ir0]
                        s1 = segs[ir1]
                        res, d0, t = seg.point_sur_segment(s0.p1)
                        slope = abs(d0) * seg.slope_left
                        a0 = abs(seg.delta_angle(s0))
                        if a0 > 0:
                            t0 = 1 + min(3 * d.facia_width, d.facia_width / sin(a0)) / s0.length
                        else:
                            t0 = 1 + d.facia_width / s0.length
                        a1 = abs(seg.delta_angle(s1))
                        if a1 > 0:
                            t1 = 1 + min(3 * d.facia_width, d.facia_width / sin(a1)) / s1.length
                        else:
                            t1 = 1 + d.facia_width / s1.length

                        x0, y0 = s0.p1
                        x1, y1 = s0.lerp(t0)
                        x2, y2 = s1.lerp(t1)
                        x3, y3 = s1.p1
                        z = self.z + d.facia_altitude - slope
                        verts.extend([
                            (x0, y0, z),
                            (x1, y1, z),
                            (x2, y2, z),
                            (x3, y3, z),
                        ])
                        z -= d.facia_height
                        verts.extend([
                            (x0, y0, z),
                            (x1, y1, z),
                            (x2, y2, z),
                            (x3, y3, z),
                        ])

                        faces.extend([
                            # top
                            (f, f + 1, f + 2, f + 3),
                            # sides
                            (f, f + 4, f + 5, f + 1),
                            (f + 1, f + 5, f + 6, f + 2),
                            (f + 2, f + 6, f + 7, f + 3),
                            (f + 3, f + 7, f + 4, f),
                            # bottom
                            (f + 4, f + 7, f + 6, f + 5)
                        ])
                        edges.append([f, f + 3])
                        edges.append([f + 1, f + 2])
                        edges.append([f + 4, f + 7])
                        edges.append([f + 5, f + 6])
                        matids.extend([idmat, idmat, idmat, idmat, idmat, idmat])
                        uvs.extend([
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)]
                        ])

                        # right side

                        f = len(verts)
                        s0 = segs[il0]
                        s1 = segs[il1]
                        res, d0, t = seg.point_sur_segment(s0.p1)
                        slope = abs(d0) * seg.slope_right
                        a0 = abs(seg.delta_angle(s0))
                        if a0 > 0:
                            t0 = 1 + min(3 * d.facia_width, d.facia_width / sin(a0)) / s0.length
                        else:
                            t0 = 1 + d.facia_width / s0.length
                        a1 = abs(seg.delta_angle(s1))
                        if a1 > 0:
                            t1 = 1 + min(3 * d.facia_width, d.facia_width / sin(a1)) / s1.length
                        else:
                            t1 = 1 + d.facia_width / s1.length
                        x0, y0 = s0.p1
                        x1, y1 = s0.lerp(t0)
                        x2, y2 = s1.lerp(t1)
                        x3, y3 = s1.p1
                        z = self.z + d.facia_altitude - slope
                        verts.extend([
                            (x0, y0, z),
                            (x1, y1, z),
                            (x2, y2, z),
                            (x3, y3, z),
                        ])
                        z -= d.facia_height
                        verts.extend([
                            (x0, y0, z),
                            (x1, y1, z),
                            (x2, y2, z),
                            (x3, y3, z),
                        ])

                        faces.extend([
                            # top
                            (f, f + 1, f + 2, f + 3),
                            # sides
                            (f, f + 4, f + 5, f + 1),
                            (f + 1, f + 5, f + 6, f + 2),
                            (f + 2, f + 6, f + 7, f + 3),
                            (f + 3, f + 7, f + 4, f),
                            # bottom
                            (f + 4, f + 7, f + 6, f + 5)
                        ])
                        edges.append([f, f + 3])
                        edges.append([f + 1, f + 2])
                        edges.append([f + 4, f + 7])
                        edges.append([f + 5, f + 6])
                        matids.extend([idmat, idmat, idmat, idmat, idmat, idmat])
                        uvs.extend([
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)]
                        ])

    def gutter(self, d, verts, faces, edges, matids, uvs):

        idmat = 5

        left, right = True, True

        segs = self.all_segs
        nodes = self.nodes


        # caps at start and end
        if d.gutter_segs % 2 == 1:
            n_faces = int((d.gutter_segs - 1) / 2)
        else:
            n_faces = int((d.gutter_segs / 2) - 1)

        df = 2 * d.gutter_segs + 1

        for idx, node in enumerate(nodes):

            # find next node in segment
            for i, s in enumerate(node.segs):

                seg = segs[s.idx]

                if s.reversed:
                    continue

                next = nodes[seg.v1_idx]

                if next.count > 1:

                    # segments sorted by angle from axis
                    # center are ids of axis segment on each node
                    # so center + 1 and center - 1
                    # are leftmost and rightmost slope segments
                    center_0 = i
                    center_1 = next.center

                    # found 2nd node
                    if center_1 > -1:

                        ############
                        # right
                        ############

                        ir0 = next.right_idx(center_1)
                        ir1 = node.left_idx(center_0)

                        ############
                        # left
                        ############

                        il0 = node.right_idx(center_0)
                        il1 = next.left_idx(center_1)

                        #####################
                        # Chenaux
                        #####################
                        if left:

                            f = len(verts)
                            s0 = segs[ir0]
                            s1 = segs[ir1]
                            res, d0, t = seg.point_sur_segment(s0.p1)
                            slope = abs(d0) * seg.slope_left
                            a0 = abs(seg.delta_angle(s0))
                            if a0 > 0:
                                scale_0 = min(3, 1 / sin(a0)) / s0.length
                            else:
                                scale_0 = 1

                            a1 = abs(seg.delta_angle(s1))
                            if a1 > 0:
                                scale_1 = min(3, 1 / sin(a1)) / s1.length
                            else:
                                scale_1 = 1

                            zt = self.z + d.facia_altitude - slope
                            z0 = self.z - slope + d.gutter_alt
                            z1 = z0 - 0.5 * d.gutter_width
                            z2 = z1 - 0.5 * d.gutter_width
                            z3 = z1 - 0.5 * d.gutter_boudin
                            dz0 = z2 - z1
                            dz1 = z3 - z1

                            tt = 1 + scale_0 * d.facia_width
                            t0 = 1 + scale_0 * d.gutter_dist
                            t1 = t0 + scale_0 * (0.5 * d.gutter_width)
                            t2 = t1 + scale_0 * (0.5 * d.gutter_width)
                            t3 = t2 + scale_0 * (0.5 * d.gutter_boudin)

                            # bord tablette
                            xt, yt = s0.lerp(tt)

                            # bord
                            x0, y0 = s0.lerp(t0)
                            # axe chenaux
                            x1, y1 = s0.lerp(t1)
                            # bord boudin interieur
                            x2, y2 = s0.lerp(t2)
                            # axe boudin
                            x3, y3 = s0.lerp(t3)

                            dx = x0 - x1
                            dy = y0 - y1

                            verts.append((xt, yt, zt))
                            # chenaux
                            da = pi / d.gutter_segs
                            for i in range(d.gutter_segs):
                                sa = sin(i * da)
                                ca = cos(i * da)
                                verts.append((x1 + dx * ca, y1 + dy * ca, z1 + dz0 * sa))

                            dx = x2 - x3
                            dy = y2 - y3

                            # boudin
                            da = -pi / (0.75 * d.gutter_segs)
                            for i in range(d.gutter_segs):
                                sa = sin(i * da)
                                ca = cos(i * da)
                                verts.append((x3 + dx * ca, y3 + dy * ca, z1 + dz1 * sa))

                            tt = 1 + scale_1 * d.facia_width
                            t0 = 1 + scale_1 * d.gutter_dist
                            t1 = t0 + scale_1 * (0.5 * d.gutter_width)
                            t2 = t1 + scale_1 * (0.5 * d.gutter_width)
                            t3 = t2 + scale_1 * (0.5 * d.gutter_boudin)

                            # bord tablette
                            xt, yt = s1.lerp(tt)

                            # bord
                            x0, y0 = s1.lerp(t0)
                            # axe chenaux
                            x1, y1 = s1.lerp(t1)
                            # bord boudin interieur
                            x2, y2 = s1.lerp(t2)
                            # axe boudin
                            x3, y3 = s1.lerp(t3)

                            dx = x0 - x1
                            dy = y0 - y1

                            # tablette
                            verts.append((xt, yt, zt))
                            faces.append((f + df, f, f + 1, f + df + 1))
                            uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                            matids.append(idmat)

                            # chenaux
                            da = pi / d.gutter_segs
                            for i in range(d.gutter_segs):
                                sa = sin(i * da)
                                ca = cos(i * da)
                                verts.append((x1 + dx * ca, y1 + dy * ca, z1 + dz0 * sa))

                            dx = x2 - x3
                            dy = y2 - y3

                            # boudin
                            da = -pi / (0.75 * d.gutter_segs)
                            for i in range(d.gutter_segs):
                                sa = sin(i * da)
                                ca = cos(i * da)
                                verts.append((x3 + dx * ca, y3 + dy * ca, z1 + dz1 * sa))

                            df = 2 * d.gutter_segs + 1

                            for i in range(1, 2 * d.gutter_segs):
                                j = i + f
                                faces.append((j, j + df, j + df + 1, j + 1))
                                uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                                matids.append(idmat)

                            """
                                    segs = 6

                                    n_faces = segs / 2 - 1

                                        0           6
                                         1         5
                                           2     4
                                              3
                            """
                            # close start
                            if next.count == 3:

                                if d.gutter_segs % 2 == 0:
                                    faces.append((f + n_faces + 3, f + n_faces + 1, f + n_faces + 2))
                                    uvs.append([(0, 0), (1, 0), (0.5, -0.5)])
                                    matids.append(idmat)

                                for i in range(n_faces):

                                    j = i + f + 1
                                    k = f + d.gutter_segs - i
                                    faces.append((j + 1, k, k + 1, j))
                                    uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                                    matids.append(idmat)

                            # close end
                            if node.count == 3:

                                f += 2 * d.gutter_segs + 1

                                if d.gutter_segs % 2 == 0:
                                    faces.append((f + n_faces + 1, f + n_faces + 3, f + n_faces + 2))
                                    uvs.append([(0, 0), (1, 0), (0.5, -0.5)])
                                    matids.append(idmat)

                                for i in range(n_faces):

                                    j = i + f + 1
                                    k = f + d.gutter_segs - i
                                    faces.append((j, k + 1, k, j + 1))
                                    uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                                    matids.append(idmat)

                        if right:
                            f = len(verts)

                            s0 = segs[il0]
                            s1 = segs[il1]
                            res, d0, t = seg.point_sur_segment(s0.p1)
                            slope = abs(d0) * seg.slope_right
                            a0 = abs(seg.delta_angle(s0))
                            if a0 > 0:
                                scale_0 = min(3, 1 / sin(a0)) / s0.length
                            else:
                                scale_0 = 1

                            a1 = abs(seg.delta_angle(s1))
                            if a1 > 0:
                                scale_1 = min(3, 1 / sin(a1)) / s1.length
                            else:
                                scale_1 = 1

                            zt = self.z + d.facia_altitude - slope
                            z0 = self.z - slope + d.gutter_alt
                            z1 = z0 - 0.5 * d.gutter_width
                            z2 = z1 - 0.5 * d.gutter_width
                            z3 = z1 - 0.5 * d.gutter_boudin
                            dz0 = z2 - z1
                            dz1 = z3 - z1

                            tt = 1 + scale_0 * d.facia_width
                            t0 = 1 + scale_0 * d.gutter_dist
                            t1 = t0 + scale_0 * (0.5 * d.gutter_width)
                            t2 = t1 + scale_0 * (0.5 * d.gutter_width)
                            t3 = t2 + scale_0 * (0.5 * d.gutter_boudin)

                            # bord tablette
                            xt, yt = s0.lerp(tt)

                            # bord
                            x0, y0 = s0.lerp(t0)
                            # axe chenaux
                            x1, y1 = s0.lerp(t1)
                            # bord boudin interieur
                            x2, y2 = s0.lerp(t2)
                            # axe boudin
                            x3, y3 = s0.lerp(t3)

                            dx = x0 - x1
                            dy = y0 - y1

                            # tablette
                            verts.append((xt, yt, zt))

                            # chenaux
                            da = pi / d.gutter_segs
                            for i in range(d.gutter_segs):
                                sa = sin(i * da)
                                ca = cos(i * da)
                                verts.append((x1 + dx * ca, y1 + dy * ca, z1 + dz0 * sa))

                            dx = x2 - x3
                            dy = y2 - y3

                            # boudin
                            da = -pi / (0.75 * d.gutter_segs)
                            for i in range(d.gutter_segs):
                                sa = sin(i * da)
                                ca = cos(i * da)
                                verts.append((x3 + dx * ca, y3 + dy * ca, z1 + dz1 * sa))

                            tt = 1 + scale_1 * d.facia_width
                            t0 = 1 + scale_1 * d.gutter_dist
                            t1 = t0 + scale_1 * (0.5 * d.gutter_width)
                            t2 = t1 + scale_1 * (0.5 * d.gutter_width)
                            t3 = t2 + scale_1 * (0.5 * d.gutter_boudin)

                            # tablette
                            xt, yt = s1.lerp(tt)

                            # bord
                            x0, y0 = s1.lerp(t0)
                            # axe chenaux
                            x1, y1 = s1.lerp(t1)
                            # bord boudin interieur
                            x2, y2 = s1.lerp(t2)
                            # axe boudin
                            x3, y3 = s1.lerp(t3)

                            dx = x0 - x1
                            dy = y0 - y1

                            # tablette
                            verts.append((xt, yt, zt))
                            faces.append((f + df, f, f + 1, f + df + 1))
                            uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                            matids.append(idmat)

                            # chenaux
                            da = pi / d.gutter_segs
                            for i in range(d.gutter_segs):
                                sa = sin(i * da)
                                ca = cos(i * da)
                                verts.append((x1 + dx * ca, y1 + dy * ca, z1 + dz0 * sa))

                            dx = x2 - x3
                            dy = y2 - y3

                            # boudin
                            da = -pi / (0.75 * d.gutter_segs)
                            for i in range(d.gutter_segs):
                                sa = sin(i * da)
                                ca = cos(i * da)
                                verts.append((x3 + dx * ca, y3 + dy * ca, z1 + dz1 * sa))

                            for i in range(1, 2 * d.gutter_segs):
                                j = i + f
                                faces.append((j, j + df, j + df + 1, j + 1))
                                uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                                matids.append(idmat)

                            # close start
                            if node.count == 3:

                                if d.gutter_segs % 2 == 0:
                                    faces.append((f + n_faces + 3, f + n_faces + 1, f + n_faces + 2))
                                    uvs.append([(0, 0), (1, 0), (0.5, -0.5)])
                                    matids.append(idmat)

                                for i in range(n_faces):
                                    j = i + f + 1
                                    k = f + d.gutter_segs - i
                                    faces.append((j + 1, k, k + 1, j))
                                    uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                                    matids.append(idmat)

                            # close end
                            if next.count == 3:

                                f += 2 * d.gutter_segs + 1

                                if d.gutter_segs % 2 == 0:
                                    faces.append((f + n_faces + 1, f + n_faces + 3, f + n_faces + 2))
                                    uvs.append([(0, 0), (1, 0), (0.5, -0.5)])
                                    matids.append(idmat)

                                for i in range(n_faces):

                                    j = i + f + 1
                                    k = f + d.gutter_segs - i
                                    faces.append((j, k + 1, k, j + 1))
                                    uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                                    matids.append(idmat)

    def beam_primary(self, d, verts, faces, edges, matids, uvs):

        idmat = 3

        segs = self.all_segs
        nodes = self.nodes

        for idx, node in enumerate(nodes):
            # segment always between 2 nodes
            # create edge between rightmost of first node to leftmost of next node
            # same for other side

            # find next node in segment
            for i, s in enumerate(node.segs):

                seg = segs[s.idx]

                if s.reversed:
                    continue

                if nodes[seg.v1_idx].count > 1:
                    next = nodes[seg.v1_idx]
                    # segments sorted by angle from axis
                    # center are ids of axis segment on each node
                    # so center + 1 and center - 1
                    # are leftmost and rightmost slope segments
                    center_0 = i
                    center_1 = next.center

                    # found 2nd node
                    if center_1 > -1:

                        ############
                        # right
                        ############

                        ir0 = next.right_idx(center_1)
                        ir1 = node.left_idx(center_0)

                        ############
                        # left
                        ############

                        il0 = node.right_idx(center_0)
                        il1 = next.left_idx(center_1)

                        ####################
                        # Poutre Faitiere
                        ####################

                        """
                         1___________________2   left
                        0|___________________|3  axis
                         |___________________|   right
                         5                   4
                        """
                        f = len(verts)

                        left = seg.offset(0.5 * d.beam_primary_width)
                        right = seg.offset(-0.5 * d.beam_primary_width)

                        # offset from roof border
                        s0 = segs[il0]
                        s1 = segs[il1]
                        s2 = segs[ir0]
                        s3 = segs[ir1]
                        t0 = 0
                        t1 = 1
                        if node.count == 3:
                            s0 = s0.offset(d.beam_primary_offset)
                            s3 = s3.offset(-d.beam_primary_offset)
                            t0 = -d.beam_primary_offset / seg.length

                        if next.count == 3:
                            s1 = s1.offset(-d.beam_primary_offset)
                            s2 = s2.offset(d.beam_primary_offset)
                            t1 = 1 + d.beam_primary_offset / seg.length

                        res, p1, t = left.intersect(s0)
                        if not res:
                            continue
                        res, p2, t = left.intersect(s1)
                        if not res:
                            continue
                        res, p4, t = right.intersect(s2)
                        if not res:
                            continue
                        res, p5, t = right.intersect(s3)
                        if not res:
                            continue
                        x0, y0 = seg.lerp(t0)
                        x1, y1 = p1
                        x2, y2 = p2
                        x3, y3 = seg.lerp(t1)
                        x4, y4 = p4
                        x5, y5 = p5
                        z = self.z + d.beam_primary_alt
                        verts.extend([
                            (x0, y0, z),
                            (x1, y1, z),
                            (x2, y2, z),
                            (x3, y3, z),
                            (x4, y4, z),
                            (x5, y5, z)
                        ])
                        z -= d.beam_primary_height
                        verts.extend([
                            (x0, y0, z),
                            (x1, y1, z),
                            (x2, y2, z),
                            (x3, y3, z),
                            (x4, y4, z),
                            (x5, y5, z)
                        ])
                        faces.extend([
                            # top
                            (f, f + 1, f + 2, f + 3),
                            (f + 5, f, f + 3, f + 4),
                            # sides
                            (f, f + 6, f + 7, f + 1),
                            (f + 1, f + 7, f + 8, f + 2),
                            (f + 2, f + 8, f + 9, f + 3),
                            (f + 3, f + 9, f + 10, f + 4),
                            (f + 4, f + 10, f + 11, f + 5),
                            (f + 5, f + 11, f + 6, f),
                            # bottom
                            (f + 6, f + 9, f + 8, f + 7),
                            (f + 11, f + 10, f + 9, f + 6)
                        ])

                        edges.append([f + 1, f + 2])
                        edges.append([f + 5, f + 4])
                        edges.append([f + 7, f + 8])
                        edges.append([f + 11, f + 10])
                        matids.extend([
                            idmat, idmat, idmat, idmat, idmat,
                            idmat, idmat, idmat, idmat, idmat
                            ])
                        uvs.extend([
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)],
                            [(0, 0), (1, 0), (1, 1), (0, 1)]
                        ])

    def rafter(self, d, verts, faces, edges, matids, uvs):

        idmat = 4

        # Rafters / Chevrons
        start = max(0.001 + 0.5 * d.rafter_width, d.rafter_start)
        chevrons_alt = -0.001

        segs = self.all_segs
        nodes = self.nodes

        for idx, node in enumerate(nodes):
            # segment always between 2 nodes
            # create edge between rightmost of first node to leftmost of next node
            # same for other side

            # find next node in segment
            for i, s in enumerate(node.segs):

                seg = segs[s.idx]

                if s.reversed:
                    continue

                if nodes[seg.v1_idx].count > 1:
                    next = nodes[seg.v1_idx]
                    # segments sorted by angle from axis
                    # center are ids of axis segment on each node
                    # so center + 1 and center - 1
                    # are leftmost and rightmost slope segments
                    center_0 = i
                    center_1 = next.center

                    # found 2nd node
                    if center_1 > -1:

                        ############
                        # right
                        ############

                        ir0 = next.right_idx(center_1)
                        ir1 = node.left_idx(center_0)

                        ############
                        # left
                        ############

                        il0 = node.right_idx(center_0)
                        il1 = next.left_idx(center_1)

                        # right part is larger than axis: compute t param in axis
                        res, d0, u = seg.point_sur_segment(segs[ir1].p1)
                        res, dr, v = seg.point_sur_segment(segs[ir0].p1)
                        trmin = min(0, u)
                        trmax = max(1, v)

                        # left part is larger than axis: compute t param in axis
                        res, d0, u = seg.point_sur_segment(segs[il0].p1)
                        res, dl, v = seg.point_sur_segment(segs[il1].p1)
                        tlmin = min(0, u)
                        tlmax = max(1, v)

                        ######################
                        # Left part chevrons
                        ######################

                        f = len(verts)

                        t0 = trmin + (start - 0.5 * d.rafter_width) / seg.length
                        t1 = trmin + (start + 0.5 * d.rafter_width) / seg.length

                        tx = start / seg.length
                        dt = d.rafter_spacing / seg.length

                        n_items = max(1, round((trmax - trmin) / dt, 0))

                        dt = ((trmax - trmin) - 2 * tx) / n_items

                        for j in range(int(n_items) + 1):
                            #
                            n0 = seg.sized_normal(t1 + j * dt, - seg.width_left)
                            n1 = seg.sized_normal(t0 + j * dt, - seg.width_left)
                            slope = n0.length * seg.slope_left

                            # verts start from axis
                            right_before = t1 + j * dt > 0 and t1 + j * dt < 1
                            left_before = t0 + j * dt > 0 and t0 + j * dt < 1

                            if right_before:
                                z = self.z + d.rafter_alt
                                x, y = n0.p0
                                verts.append((x, y, z))
                                verts.append((x, y, z - d.rafter_height))

                            self.intersect_chevron(n0, ir0, ir1, slope, d.rafter_height, d.rafter_alt, segs, verts)

                            if not right_before:
                                z = self.z - slope + d.rafter_alt
                                x, y = n0.p1
                                verts.append((x, y, z))
                                verts.append((x, y, z - d.rafter_height))

                            if left_before:
                                z = self.z + d.rafter_alt
                                x, y = n1.p0
                                verts.append((x, y, z))
                                verts.append((x, y, z - d.rafter_height))

                            self.intersect_chevron(n1, ir0, ir1, slope, d.rafter_height, d.rafter_alt, segs, verts)

                            if not left_before:
                                z = self.z - slope + d.rafter_alt
                                x, y = n1.p1
                                verts.append((x, y, z))
                                verts.append((x, y, z - d.rafter_height))

                            edges.append([f, f + 2])
                            edges.append([f + 1, f + 3])
                            edges.append([f + 4, f + 6])
                            edges.append([f + 5, f + 7])
                            faces.extend([
                                (f, f + 4, f + 5, f + 1),
                                (f + 1, f + 5, f + 7, f + 3),
                                (f + 2, f + 3, f + 7, f + 6),
                                (f + 2, f + 6, f + 4, f),
                                (f, f + 1, f + 3, f + 2),
                                (f + 5, f + 4, f + 6, f + 7)
                            ])
                            matids.extend([idmat, idmat, idmat, idmat, idmat, idmat])
                            uvs.extend([
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)]
                            ])
                            f += 8

                        ######################
                        # Right part chevrons
                        ######################

                        f = len(verts)

                        t0 = tlmin + (start - 0.5 * d.rafter_width) / seg.length
                        t1 = tlmin + (start + 0.5 * d.rafter_width) / seg.length

                        tx = start / seg.length
                        dt = d.rafter_spacing / seg.length

                        n_items = max(1, round((tlmax - tlmin) / dt, 0))

                        dt = ((tlmax - tlmin) - 2 * tx) / n_items

                        for j in range(int(n_items) + 1):
                            n0 = seg.sized_normal(t0 + j * dt, seg.width_right)
                            n1 = seg.sized_normal(t1 + j * dt, seg.width_right)
                            slope = n0.length * seg.slope_right

                            # verts start from axis
                            right_before = t0 + j * dt > 0 and t0 + j * dt < 1
                            left_before = t1 + j * dt > 0 and t1 + j * dt < 1

                            if right_before:
                                z = self.z + d.rafter_alt
                                x, y = n0.p0
                                verts.append((x, y, z))
                                verts.append((x, y, z - d.rafter_height))

                            self.intersect_chevron(n0, il0, il1, slope, d.rafter_height, d.rafter_alt, segs, verts)

                            if not right_before:
                                z = self.z - slope + d.rafter_alt
                                x, y = n0.p1
                                verts.append((x, y, z))
                                verts.append((x, y, z - d.rafter_height))

                            if left_before:
                                z = self.z + d.rafter_alt
                                x, y = n1.p0
                                verts.append((x, y, z))
                                verts.append((x, y, z - d.rafter_height))

                            self.intersect_chevron(n1, il0, il1, slope, d.rafter_height, d.rafter_alt, segs, verts)

                            if not left_before:
                                z = self.z - slope + d.rafter_alt
                                x, y = n1.p1
                                verts.append((x, y, z))
                                verts.append((x, y, z - d.rafter_height))

                            edges.append([f, f + 2])
                            edges.append([f + 1, f + 3])
                            edges.append([f + 4, f + 6])
                            edges.append([f + 5, f + 7])

                            faces.extend([
                                (f, f + 4, f + 5, f + 1),
                                (f + 1, f + 5, f + 7, f + 3),
                                (f + 2, f + 3, f + 7, f + 6),
                                (f + 2, f + 6, f + 4, f),
                                (f, f + 1, f + 3, f + 2),
                                (f + 5, f + 4, f + 6, f + 7)
                            ])
                            matids.extend([idmat, idmat, idmat, idmat, idmat, idmat])
                            uvs.extend([
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)],
                                [(0, 0), (1, 0), (1, 1), (0, 1)]
                            ])
                            f += 8

    def hips(self, d, verts, faces, edges, matids, uvs):

        idmat_valley = 5
        idmat = 6
        idmat_poutre = 4

        sx, sy, sz = d.hip_size_x, d.hip_size_y, d.hip_size_z

        if d.hip_model == 'ROUND':

            # round hips
            t_pts = [Vector((sx * x, sy * y, sz * z)) for x, y, z in [
                (-0.5, 0.34, 0.08), (-0.5, 0.32, 0.19), (0.5, -0.4, -0.5),
                (0.5, 0.4, -0.5), (-0.5, 0.26, 0.28), (-0.5, 0.16, 0.34),
                (-0.5, 0.05, 0.37), (-0.5, -0.05, 0.37), (-0.5, -0.16, 0.34),
                (-0.5, -0.26, 0.28), (-0.5, -0.32, 0.19), (-0.5, -0.34, 0.08),
                (-0.5, -0.25, -0.5), (-0.5, 0.25, -0.5), (0.5, -0.08, 0.5),
                (0.5, -0.5, 0.08), (0.5, -0.24, 0.47), (0.5, -0.38, 0.38),
                (0.5, -0.47, 0.24), (0.5, 0.5, 0.08), (0.5, 0.08, 0.5),
                (0.5, 0.47, 0.24), (0.5, 0.38, 0.38), (0.5, 0.24, 0.47)
            ]]
            t_faces = [
                (23, 22, 4, 5), (3, 19, 21, 22, 23, 20, 14, 16, 17, 18, 15, 2), (14, 20, 6, 7),
                (18, 17, 9, 10), (15, 18, 10, 11), (21, 19, 0, 1), (17, 16, 8, 9),
                (13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 1, 0), (19, 3, 13, 0), (20, 23, 5, 6), (22, 21, 1, 4),
                (3, 2, 12, 13), (2, 15, 11, 12), (16, 14, 7, 8)
            ]
            t_uvs = [
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.5, 1.0), (0.75, 0.93), (0.93, 0.75),
                (1.0, 0.5), (0.93, 0.25), (0.75, 0.07),
                (0.5, 0.0), (0.25, 0.07), (0.07, 0.25),
                (0.0, 0.5), (0.07, 0.75), (0.25, 0.93)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.5, 1.0), (0.75, 0.93), (0.93, 0.75),
                (1.0, 0.5), (0.93, 0.25), (0.75, 0.07),
                (0.5, 0.0), (0.25, 0.07), (0.07, 0.25),
                (0.0, 0.5), (0.07, 0.75), (0.25, 0.93)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            ]
            # affect vertex with slope
            t_left = []
            t_right = []

        elif d.hip_model == 'ETERNIT':

            # square hips "eternit like"
            t_pts = [Vector((sx * x, sy * y, sz * z)) for x, y, z in [
                (0.5, 0.5, 0.0), (-0.5, 0.5, -0.5), (0.5, -0.5, 0.0),
                (-0.5, -0.5, -0.5), (0.5, 0.0, 0.0), (-0.5, -0.0, -0.5),
                (0.5, 0.0, 0.5), (0.5, -0.5, 0.5), (-0.5, -0.5, 0.0),
                (-0.5, -0.0, 0.0), (0.5, 0.5, 0.5), (-0.5, 0.5, 0.0)]
            ]
            t_faces = [
                (4, 2, 3, 5), (0, 4, 5, 1), (6, 9, 8, 7),
                (10, 11, 9, 6), (0, 10, 6, 4), (5, 9, 11, 1),
                (2, 7, 8, 3), (1, 11, 10, 0), (4, 6, 7, 2), (3, 8, 9, 5)
                ]
            t_uvs = [
                [(0.0, 0.5), (0.0, 1.0), (1.0, 1.0), (1.0, 0.5)], [(0.0, 0.0), (0.0, 0.5), (1.0, 0.5), (1.0, 0.0)],
                [(0.0, 0.5), (1.0, 0.5), (1.0, 1.0), (0.0, 1.0)], [(0.0, 0.0), (1.0, 0.0), (1.0, 0.5), (0.0, 0.5)],
                [(0.0, 0.5), (0.0, 1.0), (0.5, 1.0), (0.5, 0.5)], [(0.5, 0.5), (0.5, 1.0), (0.0, 1.0), (0.0, 0.5)],
                [(0.0, 0.5), (0.0, 1.0), (1.0, 1.0), (1.0, 0.5)], [(0.0, 0.5), (0.0, 1.0), (-1.0, 1.0), (-1.0, 0.5)],
                [(0.5, 0.5), (0.5, 1.0), (1.0, 1.0), (1.0, 0.5)], [(0.0, 0.5), (0.0, 1.0), (-0.5, 1.0), (-0.5, 0.5)]
            ]
            t_left = [2, 3, 7, 8]
            t_right = [0, 1, 10, 11]

        elif d.hip_model == 'FLAT':
            # square hips "eternit like"
            t_pts = [Vector((sx * x, sy * y, sz * z)) for x, y, z in [
                (-0.5, -0.4, 0.0), (-0.5, -0.4, 0.5), (-0.5, 0.4, 0.0),
                (-0.5, 0.4, 0.5), (0.5, -0.5, 0.5), (0.5, -0.5, 1.0),
                (0.5, 0.5, 0.5), (0.5, 0.5, 1.0), (-0.5, 0.33, 0.0),
                (-0.5, -0.33, 0.0), (0.5, -0.33, 0.5), (0.5, 0.33, 0.5),
                (-0.5, 0.33, -0.5), (-0.5, -0.33, -0.5), (0.5, -0.33, -0.5),
                (0.5, 0.33, -0.5)]
            ]
            t_faces = [
                (0, 1, 3, 2, 8, 9), (2, 3, 7, 6), (6, 7, 5, 4, 10, 11),
                (4, 5, 1, 0), (9, 10, 4, 0), (7, 3, 1, 5),
                (2, 6, 11, 8), (9, 8, 12, 13), (12, 15, 14, 13),
                (8, 11, 15, 12), (10, 9, 13, 14), (11, 10, 14, 15)]
            t_uvs = [
                [(0.5, 1.0), (0.93, 0.75), (0.93, 0.25), (0.5, 0.0), (0.07, 0.25), (0.07, 0.75)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.5, 1.0), (0.93, 0.75), (0.93, 0.25), (0.5, 0.0), (0.07, 0.25), (0.07, 0.75)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
                [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            ]
            t_left = []
            t_right = []

        t_idmats = [idmat for f in t_faces]

        angle_90 = round(pi / 2, 3)

        segs = self.all_segs
        nodes = self.nodes

        for idx, node in enumerate(nodes):

            n_node = node.count

            # faitieres sur parties en pente
            # seulement entre segments
            if n_node > 3:

                last = segs[node.last_idx]

                for i, s in enumerate(node.segs):
                    seg = segs[s.idx]

                    if s.reversed:
                        last = seg
                        continue

                    if seg.constraint_type == 'SLOPE':

                        # if n_node < 4 and seg.enforce_part == 'AUTO':
                        #     continue

                        da = round(seg.delta_angle(last), 3)
                        # t param at segment end + border
                        if da > 0:
                            tmax = 1 + min(3 * d.tile_border, d.tile_border / sin(da)) / seg.length
                        else:
                            tmax = 1

                        p1 = seg.lerp(tmax)

                        res, d0, t = last.point_sur_segment(p1)

                        # delta height at roof border
                        slope = abs(d0) * seg.slope_left

                        if abs(da) != angle_90:

                            # poutre inside
                            if d.beam_sec_enable:
                                f = len(verts)
                                s0 = seg.offset(0.5 * d.beam_sec_width)
                                s1 = seg.offset(-0.5 * d.beam_sec_width)
                                x0, y0 = s0.p0
                                x1, y1 = s0.p1
                                x2, y2 = s1.p0
                                x3, y3 = s1.p1
                                z0 = self.z + d.beam_sec_alt
                                z1 = z0 - d.beam_sec_height
                                z2 = z0 - slope
                                z3 = z2 - d.beam_sec_height
                                verts.extend([
                                    (x0, y0, z0),
                                    (x0, y0, z1),
                                    (x1, y1, z2),
                                    (x1, y1, z3),
                                    (x2, y2, z0),
                                    (x2, y2, z1),
                                    (x3, y3, z2),
                                    (x3, y3, z3),
                                ])
                                faces.extend([
                                    (f, f + 4, f + 5, f + 1),
                                    (f + 1, f + 5, f + 7, f + 3),
                                    (f + 2, f + 3, f + 7, f + 6),
                                    (f + 2, f + 6, f + 4, f),
                                    (f, f + 1, f + 3, f + 2),
                                    (f + 5, f + 4, f + 6, f + 7)
                                ])
                                matids.extend([
                                    idmat_poutre, idmat_poutre, idmat_poutre,
                                    idmat_poutre, idmat_poutre, idmat_poutre
                                ])
                                uvs.extend([
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)],
                                    [(0, 0), (1, 0), (1, 1), (0, 1)]
                                ])

                            # valley :
                            if t > 0 and t < 1:
                                if d.valley_enable:
                                    """
                                        Note:
                                        Enforce for valley
                                        only make sense with T childs
                                        and depends on parent
                                        enforce for valley seg.enforce_part == 'VALLEY' or
                                    """
                                    # decalage pour les 2 bords
                                    f = len(verts)
                                    if abs(da) > 0:
                                        op = abs(tan(angle_90 - da))
                                        tside = min(6 * d.tile_couloir, 2 * d.tile_couloir * op) / seg.length
                                    else:
                                        tside = 0

                                    s0 = seg.offset(2 * d.tile_couloir)
                                    s1 = seg.offset(-2 * d.tile_couloir)

                                    x0, y0 = seg.p0
                                    x1, y1 = s0.lerp(tside)
                                    x2, y2 = s0.lerp(tmax + tside)
                                    x3, y3 = p1
                                    x4, y4 = s1.lerp(tmax + tside)
                                    x5, y5 = s1.lerp(tside)
                                    z0 = self.z + d.valley_altitude
                                    z1 = z0
                                    z2 = z0 - slope
                                    z3 = z0 - slope

                                    verts.extend([
                                        (x0, y0, z0),
                                        (x1, y1, z1),
                                        (x2, y2, z2),
                                        (x3, y3, z3),
                                        (x4, y4, z2),
                                        (x5, y5, z1)
                                    ])
                                    faces.extend([
                                        # top
                                        (f, f + 1, f + 2, f + 3),
                                        (f + 5, f, f + 3, f + 4),
                                    ])

                                    matids.extend([
                                        idmat_valley, idmat_valley
                                        ])
                                    uvs.extend([
                                        [(0, 0), (1, 0), (1, 1), (0, 1)],
                                        [(0, 0), (1, 0), (1, 1), (0, 1)]
                                    ])

                            # hip seg.enforce_part == 'HIP' or
                            elif d.hip_enable and (t < 0 or t > 1):

                                vx = seg.v.to_3d()
                                vx.z = -seg.width_left * seg.slope_left
                                vx = vx.normalized()
                                vy = seg.cross_z.normalized().to_3d()
                                vz = vx.cross(-vy)
                                x0, y0 = seg.p0 + vx.to_2d().normalized() * 0.3 * d.hip_size_y
                                z0 = self.z + d.hip_alt
                                tM = Matrix([
                                    [vx.x, vy.x, vz.x, x0],
                                    [vx.y, vy.y, vz.y, y0],
                                    [vx.z, vy.z, vz.z, z0],
                                    [0, 0, 0, 1]
                                ])
                                space_2d = tmax * seg.length - 0.1 * d.hip_size_y
                                space_z = tmax * seg.width_left * seg.slope_left
                                space_x = sqrt(space_2d * space_2d + space_z * space_z)
                                n_x = 1 + int(space_x / d.hip_space_x)
                                dx = space_x / n_x
                                x0 = 0.5 * dx

                                t_verts = [p for p in t_pts]

                                # apply slope
                                dz = t_verts[i].y * cos(abs(da))

                                for i in t_left:
                                    t_verts[i] = t_verts[i].copy()
                                    t_verts[i].z -= dz
                                for i in t_right:
                                    t_verts[i] = t_verts[i].copy()
                                    t_verts[i].z -= dz

                                for k in range(n_x):
                                    lM = tM * Matrix([
                                        [1, 0, 0, x0 + k * dx],
                                        [0, -1, 0, 0],
                                        [0, 0, 1, 0],
                                        [0, 0, 0, 1]
                                    ])
                                    f = len(verts)

                                    verts.extend([lM * p for p in t_verts])
                                    faces.extend([tuple(i + f for i in p) for p in t_faces])
                                    matids.extend(t_idmats)
                                    uvs.extend(t_uvs)

                    last = seg

            if d.hip_enable:
                # hip horizontales
                for i, s in enumerate(node.segs):

                    seg = segs[s.idx]

                    if s.reversed:
                        continue

                    if nodes[seg.v1_idx].count > 1:

                        next = nodes[seg.v1_idx]

                        n_next = next.count

                        tmin = 0
                        tmax = 1

                        if n_node == 3:
                            tmin = 0 - d.tile_side / seg.length

                        if n_next == 3:
                            tmax = 1 + d.tile_side / seg.length

                        # print("tmin:%s tmax:%s" % (tmin, tmax))
                        ####################
                        # Faitiere
                        ####################

                        f = len(verts)
                        s_len = (tmax - tmin) * seg.length
                        n_obj = 1 + int(s_len / d.hip_space_x)
                        dx = s_len / n_obj
                        x0 = 0.5 * dx
                        v = seg.v.normalized()
                        p0 = seg.lerp(tmin)
                        tM = Matrix([
                            [v.x, v.y, 0, p0.x],
                            [v.y, -v.x, 0, p0.y],
                            [0, 0, 1, self.z + d.hip_alt],
                            [0, 0, 0, 1]
                        ])
                        t_verts = [p.copy() for p in t_pts]

                        # apply slope
                        for i in t_left:
                            t_verts[i].z += t_verts[i].y * (seg.slope_right - d.tile_size_z / d.tile_size_y)
                        for i in t_right:
                            t_verts[i].z -= t_verts[i].y * (seg.slope_left - d.tile_size_z / d.tile_size_y)

                        for k in range(n_obj):
                            lM = tM * Matrix([
                                [1, 0, 0, x0 + k * dx],
                                [0, -1, 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]
                            ])
                            v = len(verts)
                            verts.extend([lM * p for p in t_verts])
                            faces.extend([tuple(i + v for i in f) for f in t_faces])
                            matids.extend(t_idmats)
                            uvs.extend(t_uvs)

    def slope_for_angle(self, slope, angle):
        """
            slope by angle from axis
        """
        return sin(angle) * slope

    def width_for_angle(self, width, angle):
        """
            length of projection of width over axis

             / width perpendicular
            /|\
           / | \/
     angle/__|_|   axis
              l
        """
        return width / sin(angle)

    def make_hole(self, context, hole_obj, d):
        """
            Hole
            follow axis based on node roots
            take an offset from horiziontal parts
            build 2 verts

            top view                  front view
                   11-10 /
          5-4___s4___/  /                0-6      z0
       0-3  /__________/  / 6-9         /|\
            \_____s2_____/        5-11 / | \ 1-7  z1 - z2
         1-2            7-8           |__|__|
                                  4-10  3-9  2-8  z3
        """

        segs = self.all_segs
        node = self.nodes[0]
        root = segs[node.root.idx]
        il0 = node.left_idx(node.center)
        ir0 = node.right_idx(node.center)
        s0 = segs[il0]
        s1 = segs[ir0]
        s2 = root.offset(min(0.001 - root.width_left, d.hole_offset_left - root.width_left))
        s4 = root.offset(max(root.width_right - 0.001, root.width_right - d.hole_offset_right))
        x0, y0 = root.p0
        res, p, t = s0.intersect(s2)
        x1, y1 = p
        res, p, t = s1.intersect(s4)
        x4, y4 = p
        t = max(0.001, (root.length - d.hole_offset_front) /  root.length)
        x6, y6 = root.lerp(t)
        x7, y7 = s2.lerp(t)
        x10, y10 = s4.lerp(t)
        z0 = self.z + 1.0
        z1 = z0
        z2 = z0
        z3 = z0 - max(z1, z2) - 1.0

        verts = [
            (x0, y0, z0),
            # right
            (x1, y1, z1),
            (x1, y1, z3),
            # center
            (x0, y0, z3),
            # left
            (x4, y4, z3),
            (x4, y4, z2),

            (x6, y6, z0),
            # right
            (x7, y7, z1),
            (x7, y7, z3),
            # center
            (x6, y6, z3),
            # left
            (x10, y10, z3),
            (x10, y10, z2)
        ]
        faces = [
            (0, 1, 2, 3),
            (0, 3, 4, 5),

            (0, 6, 7, 1),
            (1, 7, 8, 2),
            (2, 8, 9, 3),
            (3, 9, 10, 4),
            (4, 10, 11, 5),
            (5, 11, 6, 0),

            (6, 9, 8, 7),
            (6, 11, 10, 9)
        ]
        matids = None
        uvs = None
        bmed.buildmesh(
            context, hole_obj, verts, faces, matids=matids, uvs=uvs,
            weld=False, clean=False, auto_smooth=False, temporary=False)


# bpy.app.debug = True

def update(self, context):
    self.update(context)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


def update_parent(self, context):

    # update part a0
    o = context.active_object
    p, d = self.find_parent(context)

    if d is not None:
        # add a reference point
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
        if p.parent is None:
            # Either link to o.parent or create new parent
            link_to_parent = True
            if o.parent is None:
                # create a reference point and make it active
                x, y, z = p.bound_box[0]
                context.scene.cursor_location = p.matrix_world * Vector((x, y, z))
                # fix issue #9
                context.scene.objects.active = o
                bpy.ops.archipack.reference_point()
                o.select = True
            else:
                context.scene.objects.active = o.parent
            p.select = True
        else:
            # w has parent
            if o.parent is not p.parent:
                link_to_parent = True
                context.scene.objects.active = p.parent
                o.select = True
                if o.parent is not None:
                    # store o.parent to delete it
                    to_delete = o.parent
                    for c in o.parent.children:
                        if c is not o:
                            c.hide_select = False
                            c.select = True

        parent = context.active_object

        if link_to_parent and bpy.ops.archipack.parent_to_reference.poll():
            bpy.ops.archipack.parent_to_reference('INVOKE_DEFAULT')

        # hide holes from select
        for c in parent.children:
            if "archipack_hybridhole" in c:
                c.hide_select = True

        # delete unneeded reference point
        if to_delete is not None:
            bpy.ops.object.select_all(action="DESELECT")
            to_delete.select = True
            context.scene.objects.active = to_delete
            if bpy.ops.object.delete.poll():
                bpy.ops.object.delete(use_global=False)


        bpy.ops.object.select_all(action="DESELECT")
        o.select = True
        context.scene.objects.active = o

        bpy.ops.archipack.generate_hole('INVOKE_DEFAULT')

        hole_obj = context.active_object
        hole_obj.select = False

        o.select = True
        context.scene.objects.active = p

        p.select = True
        if bpy.ops.archipack.single_boolean.poll():
            bpy.ops.archipack.single_boolean()

        hole_obj.matrix_world = o.matrix_world.copy()

        p.select = False

        context.scene.objects.active = o
        o.select = True
        # trigger object update
        self.parts[0].a0 = pi / 2


    elif self.t_parent != "":
        self.t_parent = ""


def update_location(self, context):
    self.update(context, relocate_only=True)


def update_childs(self, context):
    self.update(context, update_childs=True)


materials_enum = (
            ('0', 'Ceiling', '', 0),
            ('1', 'White', '', 1),
            ('2', 'Concrete', '', 2),
            ('3', 'Wood', '', 3),
            ('4', 'Metal', '', 4),
            ('5', 'Glass', '', 5)
            )


class archipack_roof_material(PropertyGroup):
    index = EnumProperty(
        items=materials_enum,
        default='4',
        update=update
        )

    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_roof.datablock(o)
            if props:
                for part in props.rail_mat:
                    if part == self:
                        return props
        return None

    def update(self, context):
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context)


class ArchipackSegment():
    type = EnumProperty(
            items=(
                ('S_SEG', 'Straight roof', '', 0),
                ('C_SEG', 'Curved roof', '', 1),
                ),
            default='S_SEG'
            )
    length = FloatProperty(
            name="length",
            min=0.01,
            max=1000.0,
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
    da = FloatProperty(
            name="angle",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    a0 = FloatProperty(
            name="angle",
            min=-2 * pi,
            max=2 * pi,
            default=0,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    dz = FloatProperty(
            name="delta z",
            default=0
            )
    offset = FloatProperty(
            name="offset",
            min=0,
            default=0,
            update=update
            )
    manipulators = CollectionProperty(type=archipack_manipulator)

    def find_in_selection(self, context):
        raise NotImplementedError

    def update(self, context, manipulable_refresh=False):
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context, manipulable_refresh)

    def draw(self, layout, context, index):
        box = layout.box()
        row = box.row()
        row.prop(self, "type", text="")
        if self.type in ['C_SEG']:
            row = box.row()
            row.prop(self, "radius")
            row = box.row()
            row.prop(self, "da")
        else:
            row = box.row()
            row.prop(self, "length")
        row = box.row()
        row.prop(self, "offset")
        row = box.row()
        row.prop(self, "a0")


class ArchipackLines():
    n_parts = IntProperty(
            name="parts",
            min=1,
            default=1, update=update_manipulators
            )
    # UI layout related
    parts_expand = BoolProperty(
            default=False
            )

    def update(self, context, manipulable_refresh=False):
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context, manipulable_refresh)

    def draw(self, layout, context):
        box = layout.box()
        row = box.row()
        if self.parts_expand:
            row.prop(self, 'parts_expand', icon="TRIA_DOWN", icon_only=True, text="Parts", emboss=False)
            box.prop(self, 'n_parts')
            for i, part in enumerate(self.parts):
                part.draw(layout, context, i)
        else:
            row.prop(self, 'parts_expand', icon="TRIA_RIGHT", icon_only=True, text="Parts", emboss=False)

    def update_parts(self):
        # print("update_parts")
        # remove rows
        # NOTE:
        # n_parts+1
        # as last one is end point of last segment or closing one
        for i in range(len(self.parts), self.n_parts + 1, -1):
            self.parts.remove(i - 1)

        # add rows
        for i in range(len(self.parts), self.n_parts + 1):
            self.parts.add()

        self.setup_manipulators()

    def setup_parts_manipulators(self):
        for i in range(self.n_parts + 1):
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
            if n_manips < 5:
                s = p.manipulators.add()
                s.type_key = "SIZE"
                s.prop1_name = "offset"
            p.manipulators[2].prop1_name = str(i)
            p.manipulators[3].prop1_name = str(i + 1)


class archipack_roof_segment(ArchipackSegment, PropertyGroup):

    bound_idx = IntProperty(
        default=0,
        min=0,
        update=update_manipulators
        )

    take_precedence = BoolProperty(
        name="Take precedence",
        description="On T segment take width precedence",
        default=False,
        update=update
        )

    constraint_type = EnumProperty(
        items=(
            ('HORIZONTAL', 'Horizontal', '', 0),
            ('SLOPE', 'Slope', '', 1)
            ),
        default='HORIZONTAL',
        update=update_manipulators
        )

    enforce_part = EnumProperty(
        name="Enforce part",
        items=(
            ('AUTO', 'Auto', '', 0),
            ('VALLEY', 'Valley', '', 1),
            ('HIP', 'Hip', '', 2)
            ),
        default='AUTO',
        update=update
        )

    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            d = archipack_roof.datablock(o)
            if d:
                for part in d.parts:
                    if part == self:
                        return d
        return None

    def draw(self, layout, context, index):
        box = layout.box()
        row = box.row()
        row.prop(self, "constraint_type", text=str(index + 1))
        if self.constraint_type == 'SLOPE':
            row = box.row()
            row.prop(self, "enforce_part", text="")

        if self.type in ['C_SEG']:
            row = box.row()
            row.prop(self, "radius")
            row = box.row()
            row.prop(self, "da")
        else:
            row = box.row()
            row.prop(self, "length")
        row = box.row()
        row.prop(self, "a0")
        row = box.row()
        row.prop(self, 'bound_idx')


class archipack_roof(ArchipackLines, ArchipackObject, Manipulable, PropertyGroup):
    parts = CollectionProperty(type=archipack_roof_segment)
    z = FloatProperty(
            name="z",
            default=3, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    slope_left = FloatProperty(
            name="L slope",
            default=0.5, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update_childs
            )
    slope_right = FloatProperty(
            name="R slope",
            default=0.5, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update_childs
            )
    width_left = FloatProperty(
            name="L width",
            default=3, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    width_right = FloatProperty(
            name="R width",
            default=3, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    do_faces = BoolProperty(
            name="Make faces",
            default=False,
            update=update
            )
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )
    quick_edit = BoolProperty(
            name="Quick Edit",
            default=True
            )
    force_update = BoolProperty(
            name="Throttle",
            default=True
            )

    tile_enable = BoolProperty(
            name="Enable",
            default=True,
            update=update
            )
    tile_solidify = BoolProperty(
            name="Solidify",
            default=True,
            update=update
            )
    tile_height = FloatProperty(
            name="Height",
            description="Amount for solidify",
            min=0,
            default=0.02,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_bevel = BoolProperty(
            name="Bevel",
            default=False,
            update=update
            )
    tile_bevel_amt = FloatProperty(
            name="Amount",
            description="Amount for bevel",
            min=0,
            default=0.02,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_bevel_segs = IntProperty(
            name="Segs",
            description="Bevel Segs",
            min=1,
            default=2,
            update=update
            )
    tile_alternate = BoolProperty(
            name="Alternate",
            default=False,
            update=update
            )
    tile_offset = FloatProperty(
            name="Offset",
            description="Offset from start",
            min=0,
            max=100,
            subtype="PERCENTAGE",
            update=update
            )
    tile_altitude = FloatProperty(
            name="Altitude",
            description="Altitude from roof",
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_size_x = FloatProperty(
            name="x",
            description="Size of tiles on x axis",
            min=0.01,
            default=0.2,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_size_y = FloatProperty(
            name="y",
            description="Size of tiles on y axis",
            min=0.01,
            default=0.3,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_size_z = FloatProperty(
            name="z",
            description="Size of tiles on z axis",
            min=0.0,
            default=0.02,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_space_x = FloatProperty(
            name="x",
            description="Space between tiles on x axis",
            min=0.01,
            default=0.2,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_space_y = FloatProperty(
            name="y",
            description="Space between tiles on y axis",
            min=0.01,
            default=0.3,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_fit_x = BoolProperty(
            name="Fit x",
            description="Fit roof on x axis",
            default=True,
            update=update
            )
    tile_fit_y = BoolProperty(
            name="Fit y",
            description="Fit roof on y axis",
            default=True,
            update=update
            )
    tile_expand = BoolProperty(
            name="Tiles",
            description="Expand tiles panel",
            default=False
            )
    tile_model = EnumProperty(
            name="model",
            items=(
                ('BRAAS1', 'Braas 1', '', 0),
                ('BRAAS2', 'Braas 2', '', 1),
                ('ETERNIT', 'Eternit', '', 2),
                ('LAUZE', 'Lauze', '', 3),
                ('ROMAN', 'Roman', '', 4),
                ('ROUND', 'Round', '', 5),
                ('PLACEHOLDER', 'Square', '', 6),
                ('ONDULEE', 'Ondule', '', 7),
                ('METAL', 'Metal', '', 8),
             #  ('USER', 'User defined', '', 7)
                ),
            default="BRAAS2",
            update=update
            )
    tile_side = FloatProperty(
            name="Side",
            description="Space on side",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_couloir = FloatProperty(
            name="Valley",
            description="Space between tiles on valley",
            min=0,
            default=0.05,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    tile_border = FloatProperty(
            name="Bottom",
            description="Tiles offset from bottom",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    gutter_expand = BoolProperty(
            name="Gutter",
            description="Expand gutter panel",
            default=False
            )
    gutter_enable = BoolProperty(
            name="Enable",
            default=True,
            update=update
            )
    gutter_alt = FloatProperty(
            name="Altitude",
            description="altitude",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    gutter_width = FloatProperty(
            name="Width",
            description="Width",
            min=0.01,
            default=0.15,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    gutter_dist = FloatProperty(
            name="Spacing",
            description="Spacing",
            min=0,
            default=0.05,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    gutter_boudin = FloatProperty(
            name="Small width",
            description="Small width",
            min=0,
            default=0.015,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    gutter_segs = IntProperty(
            default=6,
            min=1,
            name="Segs",
            update=update
            )

    beam_expand = BoolProperty(
            name="Beam",
            description="Expand beam panel",
            default=False
            )
    beam_primary_enable = BoolProperty(
            name="Primary",
            default=True,
            update=update
            )
    beam_primary_width = FloatProperty(
            name="Width",
            description="Width",
            min=0.01,
            default=0.15,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_primary_height = FloatProperty(
            name="Height",
            description="Height",
            min=0.01,
            default=0.15,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_primary_offset = FloatProperty(
            name="Offset",
            description="Distance from roof border",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_primary_alt = FloatProperty(
            name="Altitude",
            description="Altitude from roof",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_sec_enable = BoolProperty(
            name="Hip rafter",
            default=True,
            update=update
            )
    beam_sec_width = FloatProperty(
            name="Width",
            description="Width",
            min=0.01,
            default=0.15,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_sec_height = FloatProperty(
            name="Height",
            description="Height",
            min=0.01,
            default=0.15,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_sec_alt = FloatProperty(
            name="Altitude",
            description="Distance from roof",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    rafter_enable = BoolProperty(
            name="Rafter",
            default=True,
            update=update
            )
    rafter_width = FloatProperty(
            name="Width",
            description="Width",
            min=0.01,
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    rafter_height = FloatProperty(
            name="Height",
            description="Height",
            min=0.01,
            default=0.15,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    rafter_spacing = FloatProperty(
            name="Spacing",
            description="Spacing",
            min=0.1,
            default=0.7,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    rafter_start = FloatProperty(
            name="Offset",
            description="Spacing from roof border",
            min=0,
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    rafter_alt = FloatProperty(
            name="Altitude",
            description="Altitude from roof",
            max=-0.0001,
            default=-0.001,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    hip_enable = BoolProperty(
            name="Enable",
            default=True,
            update=update
            )
    hip_expand = BoolProperty(
            name="Hips",
            description="Expand hips panel",
            default=False
            )
    hip_alt = FloatProperty(
            name="Altitude",
            description="Hip altitude from roof",
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    hip_space_x = FloatProperty(
            name="Spacing",
            description="Space between hips",
            min=0.01,
            default=0.4,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    hip_size_x = FloatProperty(
            name="l",
            description="Length of hip",
            min=0.01,
            default=0.4,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    hip_size_y = FloatProperty(
            name="w",
            description="Width of hip",
            min=0.01,
            default=0.15,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    hip_size_z = FloatProperty(
            name="h",
            description="Height of hip",
            min=0.0,
            default=0.15,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    hip_model = EnumProperty(
            name="model",
            items=(
                ('ROUND', 'Round', '', 0),
                ('ETERNIT', 'Eternit', '', 1),
                ('FLAT', 'Flat', '', 2)
                ),
            default="ROUND",
            update=update
            )
    valley_altitude  = FloatProperty(
            name="x",
            description="Valley altitude from roof",
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    valley_enable = BoolProperty(
            name="Valley",
            default=True,
            update=update
            )

    facia_enable = BoolProperty(
            name="Enable",
            description="Enable Facia",
            default=True,
            update=update
            )
    facia_expand = BoolProperty(
            name="Facia",
            description="Expand facia panel",
            default=False
            )
    facia_height = FloatProperty(
            name="Height",
            description="Height",
            min=0.01,
            default=0.3,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    facia_width =FloatProperty(
            name="Width",
            description="Width",
            min=0.01,
            default=0.02,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    facia_altitude = FloatProperty(
            name="Altitude",
            description="Facia altitude from roof",
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    rake_enable = BoolProperty(
            name="Enable",
            description="Enable Rake",
            default=True,
            update=update
            )
    rake_expand = BoolProperty(
            name="Rake",
            description="Expand rake panel",
            default=False
            )
    rake_height = FloatProperty(
            name="Height",
            description="Height",
            min=0.01,
            default=0.3,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    rake_width =FloatProperty(
            name="Width",
            description="Width",
            min=0.01,
            default=0.02,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    rake_offset =FloatProperty(
            name="Offset",
            description="Offset from roof border",
            default=0.001,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    rake_altitude = FloatProperty(
            name="Altitude",
            description="Facia altitude from roof",
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    t_parent = StringProperty(
        name="Parent",
        default="",
        update=update_parent
        )
    t_part = IntProperty(
        name="Part",
        description="Parent part index",
        default=0,
        min=0,
        update=update_location
        )
    t_dist_x = FloatProperty(
        name="Dist x",
        description="Location on axis ",
        default=0,
        update=update_location
        )
    t_dist_y = FloatProperty(
        name="Dist y",
        description="Lateral distance from axis",
        min=0,
        default=0,
        update=update
        )

    hole_offset_left = FloatProperty(
        name="Left",
        description="Left distance from border",
        min=0,
        default=0,
        update=update
        )
    hole_offset_right = FloatProperty(
        name="Right",
        description="Right distance from border",
        min=0,
        default=0,
        update=update
        )
    hole_offset_front = FloatProperty(
        name="Front",
        description="Front distance from border",
        min=0,
        default=0,
        update=update
        )

    def update_parts(self):
        # print("update_parts")
        # remove rows
        # NOTE:
        # n_parts+1
        # as last one is end point of last segment or closing one
        for i in range(len(self.parts), self.n_parts, -1):
            self.parts.remove(i - 1)

        # add rows
        for i in range(len(self.parts), self.n_parts):
            bound_idx = len(self.parts)
            self.parts.add()
            self.parts[-1].bound_idx = bound_idx

        self.setup_manipulators()

    def setup_manipulators(self):
        if len(self.manipulators) < 1:
            s = self.manipulators.add()
            s.type_key = "SIZE"
            s.prop1_name = "z"
            s.normal = (0, 1, 0)
        if len(self.manipulators) < 2:
            s = self.manipulators.add()
            s.type_key = "SIZE"
            s.prop1_name = "width_left"
        if len(self.manipulators) < 3:
            s = self.manipulators.add()
            s.type_key = "SIZE"
            s.prop1_name = "width_right"

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
                s.type_key = 'DUMB_STRING'
                s.prop1_name = str(i + 1)
            p.manipulators[2].prop1_name = str(i + 1)

    def interpolate_bezier(self, pts, wM, p0, p1, resolution):
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
        self.update_parts()

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

        self.auto_update = True

    def update_path(self, context):
        user_def_path = context.scene.objects.get(self.user_defined_path)
        if user_def_path is not None and user_def_path.type == 'CURVE':
            self.from_spline(user_def_path.matrix_world, self.user_defined_resolution, user_def_path.data.splines[0])

    def get_generator(self, origin=Vector((0, 0, 0))):
        g = RoofGenerator(self, origin)

        # TODO: sort part by bound idx so deps always find parent

        for i, part in enumerate(self.parts):
            # skip part if bound_idx > parent
            # so deps always see parent
            if part.bound_idx <= i:
                g.add_part(part)
        g.locate_manipulators()
        return g

    def make_surface(self, o, verts, edges):
        bm = bmesh.new()
        for v in verts:
            bm.verts.new(v)
        bm.verts.ensure_lookup_table()
        for ed in edges:
            bm.edges.new((bm.verts[ed[0]], bm.verts[ed[1]]))
        bm.edges.ensure_lookup_table()
        # bmesh.ops.contextual_create(bm, geom=bm.edges)
        bm.to_mesh(o.data)
        bm.free()

    def find_parent(self, context):
        o = context.scene.objects.get(self.t_parent)
        return o, archipack_roof.datablock(o)

    def intersection_angle(self, t_slope, t_width, p_slope, angle):
        # 2d intersection angle between two roofs parts
        dy = abs(t_slope * t_width / p_slope)
        ca = cos(angle)
        ta = tan(angle)
        if ta == 0:
            w0 = 0
        else:
            w0 = dy * ta
        if ca == 0:
            w1 = 0
        else:
            w1 = t_width / ca
        dx = w1 - w0
        return atan2(dy, dx)

    def relocate_child(self, context, o, g, child):

        d = archipack_roof.datablock(child)

        if d is not None and d.t_part - 1 < len(g.segs):

            seg = g.segs[d.t_part]
            # adjust T part matrix_world from parent
            # T part origin located on parent axis
            # with y in parent direction
            t = (d.t_dist_x / seg.length)
            x, y, z = seg.lerp(t).to_3d()
            dy = seg.v.normalized()
            child.matrix_world = o.matrix_world * Matrix([
                [dy.x, dy.y, 0, x],
                [dy.y, -dy.x, 0, y],
                [0, 0, 1, z],
                [0, 0, 0, 1]
            ])

    def relocate_childs(self, context, o, g):
        if o.parent is not None:
            for child in o.parent.children:
                d = archipack_roof.datablock(child)
                if d is not None and d.t_parent == o.name:
                    self.relocate_child(context, o, g, child)

    def update_childs(self, context, o, g):
        if o.parent is not None:
            for child in o.parent.children:
                d = archipack_roof.datablock(child)
                if d is not None and d.t_parent == o.name:
                    child.select = True
                    context.scene.objects.active = child
                    d.update(context, force_update=True)
                    child.select = False
            o.select = True
            context.scene.objects.active = o

    def update(self,
                context,
                manipulable_refresh=False,
                relocate_only=False,
                update_childs=False,
                force_update=False):
        """
            relocate only: relocate childs
            update_childs: force childs update
            force_update: skip throttle
        """
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # clean up manipulators before any data model change
        if manipulable_refresh:
            self.manipulable_disable(context)

        self.update_parts()

        verts, edges, faces, matids, uvs = [], [], [], [], []

        y = 0
        z = self.z
        p, d = self.find_parent(context)
        g = None

        # t childs
        if d is not None:

            pg = d.get_generator()
            pg.make_roof([], [])

            if self.t_part - 1 < len(pg.segs):

                seg = pg.segs[self.t_part]

                d.relocate_child(context, p, pg, o)

                if relocate_only:
                    self.restore_context(context)
                    return

                a0 = self.parts[0].a0

                if a0 > 0:
                    # a_axis est mesure depuis la perpendiculaire à l'axe
                    a_axis = a0 - pi / 2
                    slope = seg.slope_right
                    y = self.t_dist_y
                else:
                    a_axis = a0 + pi / 2
                    slope = seg.slope_left
                    y = -self.t_dist_y

                if slope == 0:
                    slope = 0.0001

                z = d.z - self.t_dist_y * slope

                # a_right from axis cross z

                b_right = self.intersection_angle(
                    self.slope_left,
                    -self.width_left,
                    slope,
                    a_axis)

                a_right = b_right

                b_left = self.intersection_angle(
                    self.slope_right,
                    self.width_right,
                    slope,
                    a_axis)

                a_left = b_left - pi

                g = self.get_generator(origin=Vector((0, y, z)))

                # Add 'SLOPE' constraints for segment 0
                v = Vector((cos(a_left), sin(a_left)))
                s = StraightRoof(g.origin, v)
                s.v0_idx = 0
                s.constraint_type = 'SLOPE'
                # s.enforce_part = 'VALLEY'
                s.angle_0 = a_left
                s.take_precedence = False
                g.segs.append(s)

                v = Vector((cos(a_right), sin(a_right)))
                s = StraightRoof(g.origin, v)
                s.v0_idx = 0
                s.constraint_type = 'SLOPE'
                # s.enforce_part = 'VALLEY'
                s.angle_0 = a_right
                s.take_precedence = False
                g.segs.append(s)

        if g is None:
            g = self.get_generator(origin=Vector((0, y, z)))

        if len(g.segs) > 0:
            f = g.segs[0]
            # z
            n = f.straight(-1, 0).v.to_3d()
            self.manipulators[0].set_pts([(0, 0, 0), (0, 0, self.z), (1, 0, 0)], normal=n)
            # left width
            n = f.sized_normal(0, -self.width_left)
            self.manipulators[1].set_pts([n.p0.to_3d(), n.p1.to_3d(), (-1, 0, 0)])
            # right width
            n = f.sized_normal(0, self.width_right)
            self.manipulators[2].set_pts([n.p0.to_3d(), n.p1.to_3d(), (1, 0, 0)])

        g.make_roof(verts, edges)

        self.relocate_childs(context, o, g)

        if update_childs:
            self.update_childs(context, o, g)

        # print("%s" % (faces))
        g.lambris(self, verts, faces, edges, matids, uvs)
        n_verts = len(verts)

        if self.rake_enable:
            g.rake(self, verts, faces, edges, matids, uvs)

        if self.facia_enable:
            g.facia(self, verts, faces, edges, matids, uvs)

        if self.beam_primary_enable:
            g.beam_primary(self, verts, faces, edges, matids, uvs)

        if self.rafter_enable:
            g.rafter(self, verts, faces, edges, matids, uvs)

        g.hips(self, verts, faces, edges, matids, uvs)

        if self.gutter_enable:
            g.gutter(self, verts, faces, edges, matids, uvs)

        if force_update:
            self.force_update = True

        if self.do_faces or self.force_update:
                bmed.buildmesh(
                    context, o, verts, faces, matids=matids, uvs=uvs,
                    weld=False, clean=False, auto_smooth=True, temporary=False)
        else:
            self.make_surface(o, verts, edges)

        # quick edit disable
        # boolean and tiles
        # unless throttle force update
        if self.quick_edit and not force_update:
            update_hole = False
            bpy.ops.archipack.roof_throttle_update(name=o.name)

        if not self.quick_edit or self.force_update:

            update_hole = True

            if self.tile_enable:
                bpy.ops.object.mode_set(mode='EDIT')
                g.couverture(context, o, self)

            if self.quick_edit:
                self.force_update = False

            # throttle update hole too
            if d is not None:
                hole_obj = self.find_hole(o)
                if hole_obj is not None:
                    g.make_hole(context, hole_obj, self)

        # disable AutoBoolean on
        if o.parent is not None:
            for child in o.parent.children:
                m = child.modifiers.get('AutoMerge')
                if m is not None:
                    m.show_viewport = update_hole


        # enable manipulators rebuild
        if manipulable_refresh:
            self.manipulable_refresh = True


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
            context.scene.objects.link(hole_obj)
            hole_obj['archipack_hole'] = True
            hole_obj.parent = o
            hole_obj.matrix_world = o.matrix_world.copy()

        hole_obj.data.materials.clear()
        for mat in o.data.materials:
            hole_obj.data.materials.append(mat)
        g = self.get_generator()
        g.make_roof([], [])
        g.make_hole(context, hole_obj, self)
        return hole_obj

    def manipulable_setup(self, context):
        """
            NOTE:
            this one assume context.active_object is the instance this
            data belongs to, failing to do so will result in wrong
            manipulators set on active object
        """
        self.manipulable_disable(context)

        o = context.active_object

        self.setup_manipulators()

        for i, part in enumerate(self.parts):

            if i > 0:
                # start angle
                self.manip_stack.append(part.manipulators[0].setup(context, o, part))

            # length / radius + angle
            self.manip_stack.append(part.manipulators[1].setup(context, o, part))
            # index
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))


        for m in self.manipulators:
            self.manip_stack.append(m.setup(context, o, self))

    def draw(self, layout, context):
        box = layout.box()
        row = box.row()
        if self.parts_expand:
            row.prop(self, 'parts_expand', icon="TRIA_DOWN", icon_only=True, text="Parts", emboss=False)
            box.prop(self, 'n_parts')
            # box.prop(self, 'closed')
            for i, part in enumerate(self.parts):
                part.draw(layout, context, i)
        else:
            row.prop(self, 'parts_expand', icon="TRIA_RIGHT", icon_only=True, text="Parts", emboss=False)


"""
class archipack_roof_boundary(ArchipackSegment, PropertyGroup):

    def find_in_selection(self, context):
        selected = [o for o in context.selected_objects]
        for o in selected:
            d = archipack_roof.datablock(o)
            if d:
                for part in d.parts:
                    if part == self:
                        return d
        return None


class archipack_roof(ArchipackLines, ArchipackObject, Manipulable, PropertyGroup):
    # boundary
    parts = CollectionProperty(type=archipack_roof_boundary)
    # axis = PointerProperty(type=archipack_roof)

    user_defined_path = StringProperty(
            name="user defined",
            update=update_path
            )
    user_defined_resolution = IntProperty(
            name="resolution",
            min=1,
            max=128,
            default=12, update=update_path
            )
    x_offset = FloatProperty(
            name="x offset",
            min=-1000, max=1000,
            default=0.0, precision=2, step=1,
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
    angle_limit = FloatProperty(
            name="angle",
            min=0,
            max=2 * pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update_manipulators
            )

    z = FloatProperty(
            name="z",
            default=3, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    slope = FloatProperty(
            name="slope",
            default=1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    closed = BoolProperty(
            default=False,
            name="Close",
            update=update_manipulators
            )

    # Flag to prevent mesh update while making bulk changes over variables
    # use :
    # .auto_update = False
    # bulk changes
    # .auto_update = True
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )

    def setup_manipulators(self):
        if len(self.manipulators) < 2:
            s = self.manipulators.add()
            s.prop1_name = "width"
            s = self.manipulators.add()
            s.prop1_name = "height"
            s.normal = Vector((0, 1, 0))
        self.setup_parts_manipulators()

    def interpolate_bezier(self, pts, wM, p0, p1, resolution):
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
        self.update_parts()

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

        self.auto_update = True

    def update_path(self, context):
        user_def_path = context.scene.objects.get(self.user_defined_path)
        if user_def_path is not None and user_def_path.type == 'CURVE':
            self.from_spline(user_def_path.matrix_world, self.user_defined_resolution, user_def_path.data.splines[0])

    def make_surface(self, o, verts, edges):
        bm = bmesh.new()
        for v in verts:
            bm.verts.new(v)
        bm.verts.ensure_lookup_table()
        for ed in edges:
            bm.edges.new((bm.verts[ed[0]], bm.verts[ed[1]]))
        bm.edges.ensure_lookup_table()
        bmesh.ops.contextual_create(bm, geom=bm.edges)
        bm.to_mesh(o.data)
        bm.free()

    def update(self, context, manipulable_refresh=False):

        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # clean up manipulators before any data model change
        if manipulable_refresh:
            self.manipulable_disable(context)

        self.update_parts()

        verts = []
        edges = []
        faces = []
        matids = []
        uvs = []

        g = self.get_generator()

        ag = None

        for a in o.children:
            axis = archipack_roof.datablock(a)
            if axis is not None:
                ag = axis.get_generator(origin=a.matrix_world.translation - o.matrix_world.translation)
                ag.z = self.z
                ag.slope = self.slope

        # vertex index in order to build axis
        g.get_verts(verts, edges, ag)
        print("edges %s" % edges)

        if len(verts) > 2:
            self.make_surface(o, verts, edges)
        else:
            g.debug(verts)
            bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=uvs, weld=True, clean=False)

        # enable manipulators rebuild
        if manipulable_refresh:
            self.manipulable_refresh = True

        # restore context
        self.restore_context(context)

    def manipulable_setup(self, context):

        self.manipulable_disable(context)
        o = context.active_object

        n_parts = self.n_parts
        if self.closed:
            n_parts += 1

        self.setup_manipulators()

        for i, part in enumerate(self.parts):
            if i < n_parts:

                if i > 0:
                    # start angle
                    self.manip_stack.append(part.manipulators[0].setup(context, o, part))

                # length / radius + angle
                self.manip_stack.append(part.manipulators[1].setup(context, o, part))
                # index
                self.manip_stack.append(part.manipulators[3].setup(context, o, self))
                # offset
                self.manip_stack.append(part.manipulators[4].setup(context, o, part))

            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))


        for m in self.manipulators:
            self.manip_stack.append(m.setup(context, o, self))


class ARCHIPACK_PT_roof(Panel):
    bl_idname = "ARCHIPACK_PT_roof"
    bl_label = "Roof"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_roof.filter(context.active_object)

    def draw(self, context):
        prop = archipack_roof.datablock(context.active_object)
        if prop is None:
            return
        scene = context.scene
        layout = self.layout
        row = layout.row(align=True)
        row.operator('archipack.roof_manipulate')
        box = layout.box()
        row = box.row()
        row.operator("archipack.roof_axis")
        # box.label(text="Styles")
        row = box.row(align=True)
        row.operator("archipack.roof_preset_menu", text=bpy.types.ARCHIPACK_OT_roof_preset_menu.bl_label)
        row.operator("archipack.roof_preset", text="", icon='ZOOMIN')
        row.operator("archipack.roof_preset", text="", icon='ZOOMOUT').remove_active = True

        row = layout.row(align=True)
        row.prop_search(prop, "user_defined_path", scene, "objects", text="", icon='OUTLINER_OB_CURVE')
        box = layout.box()
        box.prop(prop, 'z')
        box.prop(prop, 'slope')
        box = layout.box()
        box.prop(prop, 'user_defined_resolution')
        box.prop(prop, 'x_offset')
        box.prop(prop, "closed")
        box.prop(prop, 'angle_limit')
        prop.draw(layout, context)
        # prop.axis.draw(layout, context)
"""


class ARCHIPACK_PT_roof(Panel):
    bl_idname = "ARCHIPACK_PT_roof"
    bl_label = "Roof"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_roof.filter(context.active_object)

    def draw(self, context):
        prop = archipack_roof.datablock(context.active_object)
        if prop is None:
            return
        scene = context.scene
        layout = self.layout
        row = layout.row(align=True)
        row.operator('archipack.roof_manipulate', icon='HAND')

        box = layout.box()
        row = box.row(align=True)
        row.operator("archipack.roof_preset_menu", text=bpy.types.ARCHIPACK_OT_roof_preset_menu.bl_label)
        row.operator("archipack.roof_preset", text="", icon='ZOOMIN')
        row.operator("archipack.roof_preset", text="", icon='ZOOMOUT').remove_active = True
        box = layout.box()
        box.prop_search(prop, "t_parent", scene, "objects", text="Parent", icon='OBJECT_DATA')
        p, d = prop.find_parent(context)
        if d is not None:
            box.prop(prop, 't_part')
            box.prop(prop, 't_dist_x')
            box.prop(prop, 't_dist_y')
            box.label(text="Hole")
            box.prop(prop, 'hole_offset_front')
            box.prop(prop, 'hole_offset_left')
            box.prop(prop, 'hole_offset_right')
        box = layout.box()
        box.prop(prop, 'quick_edit', icon="MOD_MULTIRES")
        box.prop(prop, 'do_faces')
        if d is None:
            box.prop(prop, 'z')
        box.prop(prop, 'slope_left')
        box.prop(prop, 'slope_right')
        box.prop(prop, 'width_left')
        box.prop(prop, 'width_right')
        # parts
        prop.draw(layout, context)
        # tiles
        box = layout.box()
        row = box.row(align=True)
        if prop.tile_expand:
            row.prop(prop, 'tile_expand', icon="TRIA_DOWN", text="Tiles", icon_only=True, emboss=False)
        else:
            row.prop(prop, 'tile_expand', icon="TRIA_RIGHT", text="Tiles", icon_only=True, emboss=False)
        row.prop(prop, 'tile_enable')
        if prop.tile_expand:
            box.prop(prop, 'tile_model', text="")

            box.prop(prop, 'tile_solidify', icon='MOD_SOLIDIFY')
            if prop.tile_solidify:
                box.prop(prop, 'tile_height')
                box.separator()
            box.prop(prop, 'tile_bevel', icon='MOD_BEVEL')
            if prop.tile_bevel:
                box.prop(prop, 'tile_bevel_amt')
                box.prop(prop, 'tile_bevel_segs')
                box.separator()
            box.label(text="Tile size")
            row = box.row(align=True)
            row.prop(prop, 'tile_size_x')
            row.prop(prop, 'tile_size_y')
            row.prop(prop, 'tile_size_z')
            box.prop(prop, 'tile_altitude')

            box.separator()
            box.label(text="Distribution")
            box.prop(prop, 'tile_alternate', icon='NLA')
            row = box.row(align=True)
            row.prop(prop, 'tile_fit_x', icon='ALIGN')
            row.prop(prop, 'tile_fit_y', icon='ALIGN')
            box.prop(prop, 'tile_offset')

            box.label(text="Spacing")
            row = box.row(align=True)
            row.prop(prop, 'tile_space_x')
            row.prop(prop, 'tile_space_y')

            box.separator()		# hip
            box.label(text="Borders")
            box.prop(prop, 'tile_side')
            box.prop(prop, 'tile_couloir')
            box.prop(prop, 'tile_border')

        box = layout.box()
        row = box.row(align=True)
        if prop.hip_expand:
            row.prop(prop, 'hip_expand', icon="TRIA_DOWN", text="Hip", icon_only=True, emboss=False)
        else:
            row.prop(prop, 'hip_expand', icon="TRIA_RIGHT", text="Hip", icon_only=True, emboss=False)
        row.prop(prop, 'hip_enable')
        if prop.hip_expand:
            box.prop(prop, 'hip_model', text="")

            box.label(text="Hip size")
            row = box.row(align=True)
            row.prop(prop, 'hip_size_x')
            row.prop(prop, 'hip_size_y')
            row.prop(prop, 'hip_size_z')
            box.prop(prop, 'hip_alt')
            box.prop(prop, 'hip_space_x')

        box = layout.box()
        row = box.row(align=True)

        if prop.beam_expand:
            row.prop(prop, 'beam_expand', icon="TRIA_DOWN", text="Beam", icon_only=True, emboss=False)
        else:
            row.prop(prop, 'beam_expand', icon="TRIA_RIGHT", text="Beam", icon_only=True, emboss=False)
        if prop.beam_expand:
            box.prop(prop, 'beam_primary_enable')
            if prop.beam_primary_enable:
                box.prop(prop, 'beam_primary_width')
                box.prop(prop, 'beam_primary_height')
                box.prop(prop, 'beam_primary_offset')
                box.prop(prop, 'beam_primary_alt')
            box.separator()
            box.prop(prop, 'beam_sec_enable')
            if prop.beam_sec_enable:
                box.prop(prop, 'beam_sec_width')
                box.prop(prop, 'beam_sec_height')
                box.prop(prop, 'beam_sec_alt')
            box.separator()
            box.prop(prop, 'rafter_enable')
            if prop.rafter_enable:
                box.prop(prop, 'rafter_height')
                box.prop(prop, 'rafter_width')
                box.prop(prop, 'rafter_spacing')
                box.prop(prop, 'rafter_start')
                box.prop(prop, 'rafter_alt')

        box = layout.box()
        row = box.row(align=True)
        if prop.gutter_expand:
            row.prop(prop, 'gutter_expand', icon="TRIA_DOWN", text="Gutter", icon_only=True, emboss=False)
        else:
            row.prop(prop, 'gutter_expand', icon="TRIA_RIGHT", text="Gutter", icon_only=True, emboss=False)
        row.prop(prop, 'gutter_enable')
        if prop.gutter_expand:
            box.prop(prop, 'gutter_alt')
            box.prop(prop, 'gutter_width')
            box.prop(prop, 'gutter_dist')
            box.prop(prop, 'gutter_boudin')
            box.prop(prop, 'gutter_segs')


        box = layout.box()
        row = box.row(align=True)
        if prop.facia_expand:
            row.prop(prop, 'facia_expand', icon="TRIA_DOWN", text="Facia", icon_only=True, emboss=False)
        else:
            row.prop(prop, 'facia_expand', icon="TRIA_RIGHT", text="Facia", icon_only=True, emboss=False)
        row.prop(prop, 'facia_enable')
        if prop.facia_expand:
            box.prop(prop, 'facia_altitude')
            box.prop(prop, 'facia_width')
            box.prop(prop, 'facia_height')

        box = layout.box()
        row = box.row(align=True)
        if prop.rake_expand:
            row.prop(prop, 'rake_expand', icon="TRIA_DOWN", text="Rake", icon_only=True, emboss=False)
        else:
            row.prop(prop, 'rake_expand', icon="TRIA_RIGHT", text="Rake", icon_only=True, emboss=False)
        row.prop(prop, 'rake_enable')
        if prop.rake_expand:
            box.prop(prop, 'rake_altitude')
            box.prop(prop, 'rake_width')
            box.prop(prop, 'rake_height')
            box.prop(prop, 'rake_offset')

        """
        box = layout.box()
        row.prop_search(prop, "user_defined_path", scene, "objects", text="", icon='OUTLINER_OB_CURVE')
        box.prop(prop, 'user_defined_resolution')
        box.prop(prop, 'angle_limit')
        """


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_roof(ArchipackCreateTool, Operator):
    bl_idname = "archipack.roof"
    bl_label = "Roof"
    bl_description = "Roof"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    def create(self, context):
        m = bpy.data.meshes.new("Roof")
        o = bpy.data.objects.new("Roof", m)
        d = m.archipack_roof.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        self.add_material(o)
        self.load_preset(d)
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = context.scene.cursor_location
            o.select = True
            context.scene.objects.active = o
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_roof_from_curve(Operator):
    bl_idname = "archipack.roof_from_curve"
    bl_label = "Roof curve"
    bl_description = "Create a roof from a curve"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    auto_manipulate = BoolProperty(default=True)

    @classmethod
    def poll(self, context):
        return context.active_object is not None and context.active_object.type == 'CURVE'

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        curve = context.active_object
        m = bpy.data.meshes.new("Roof")
        o = bpy.data.objects.new("Roof", m)
        d = m.archipack_roof.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        d.user_defined_path = curve.name
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        d.update_path(context)
        MaterialUtils.add_stair_materials(o)
        spline = curve.data.splines[0]
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
        o.select = True
        context.scene.objects.active = o
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            self.create(context)
            if self.auto_manipulate:
                bpy.ops.archipack.roof_manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to manipulate object
# ------------------------------------------------------------------


class ARCHIPACK_OT_roof_manipulate(Operator):
    bl_idname = "archipack.roof_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return archipack_roof.filter(context.active_object)

    def invoke(self, context, event):
        d = archipack_roof.datablock(context.active_object)
        d.manipulable_invoke(context)
        return {'FINISHED'}


# Update throttle (smell hack here)
# use 2 globals to store a timer and state of update_action
# TODO: dict for those per object name
# to
update_timer = None
update_timer_updating = False
throttle_delay = 1
throttle_start = 0


class ARCHIPACK_OT_roof_throttle_update(Operator):
    bl_idname = "archipack.roof_throttle_update"
    bl_label = "Update childs with a delay"

    name = StringProperty()

    def modal(self, context, event):
        global update_timer_updating
        global throttle_start
        global throttle_delay
        if event.type == 'TIMER' and not update_timer_updating:
            # cant rely on timer as another one may run
            if time.time() - throttle_start > throttle_delay:
                update_timer_updating = True
                o = context.scene.objects.get(self.name)
                # print("delay update of %s" % (self.name))
                if o is not None:
                    o.select = True
                    context.scene.objects.active = o
                    d = o.data.archipack_roof[0]
                    d.force_update = True
                    d.update(context)
                    return self.cancel(context)
        return {'PASS_THROUGH'}

    def execute(self, context):
        global update_timer
        global update_timer_updating
        if update_timer is not None:
            context.window_manager.event_timer_remove(update_timer)
            if update_timer_updating:
                return {'CANCELLED'}
            # reset update_timer so it only occurs once 0.1s after last action
            throttle_start = time.time()
            update_timer = context.window_manager.event_timer_add(throttle_delay, context.window)
            return {'CANCELLED'}
        throttle_start = time.time()
        update_timer_updating = False
        context.window_manager.modal_handler_add(self)
        update_timer = context.window_manager.event_timer_add(throttle_delay, context.window)
        return {'RUNNING_MODAL'}

    def cancel(self, context):
        global update_timer
        context.window_manager.event_timer_remove(update_timer)
        update_timer = None
        return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to load / save presets
# ------------------------------------------------------------------


class ARCHIPACK_OT_roof_preset_menu(PresetMenuOperator, Operator):
    bl_idname = "archipack.roof_preset_menu"
    bl_label = "Roof Styles"
    preset_subdir = "archipack_roof"


class ARCHIPACK_OT_roof_preset(ArchipackPreset, Operator):
    """Add a Roof Styles"""
    bl_idname = "archipack.roof_preset"
    bl_label = "Add Roof Style"
    preset_menu = "ARCHIPACK_OT_roof_preset_menu"

    @property
    def blacklist(self):
        return ['n_parts', 'parts', 'manipulators', 'user_defined_path']


def register():
    bpy.utils.register_class(archipack_roof_material)
    # bpy.utils.register_class(archipack_roof_boundary)
    bpy.utils.register_class(archipack_roof_segment)
    bpy.utils.register_class(archipack_roof)
    Mesh.archipack_roof = CollectionProperty(type=archipack_roof)
    bpy.utils.register_class(ARCHIPACK_OT_roof_preset_menu)
    bpy.utils.register_class(ARCHIPACK_PT_roof)
    bpy.utils.register_class(ARCHIPACK_OT_roof)
    bpy.utils.register_class(ARCHIPACK_OT_roof_preset)
    bpy.utils.register_class(ARCHIPACK_OT_roof_manipulate)
    bpy.utils.register_class(ARCHIPACK_OT_roof_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_roof_throttle_update)


def unregister():
    bpy.utils.unregister_class(archipack_roof_material)
    # bpy.utils.unregister_class(archipack_roof_boundary)
    bpy.utils.unregister_class(archipack_roof_segment)
    bpy.utils.unregister_class(archipack_roof)
    del Mesh.archipack_roof
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_PT_roof)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_manipulate)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_throttle_update)
