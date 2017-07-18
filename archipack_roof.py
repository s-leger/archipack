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
from math import sin, cos, pi, atan2, sqrt, tan, degrees
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_2d import Line, Arc
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_object import ArchipackCreateTool, ArchipackObject


class Roof():

    def __init__(self):
        self.angle_0 = 0
        self.v0_idx = 0
        self.v1_idx = 0
        self.constraint_type = None
        self.slope_left = 1
        self.slope_right = 1
        self.width_left = 1
        self.width_right = 1
        self.auto_left = 'AUTO'
        self.auto_right = 'AUTO'
        self.backward_left = 0
        self.backward_right = 0
        self.type = 'SIDE'
        # force hip or valley
        self.enforce_part = 'AUTO'

    @property
    def copy(self):
        s = StraightRoof(self.p.copy(), self.v.copy())
        s.angle_0 = self.angle_0
        s.v0_idx = self.v0_idx
        s.v1_idx = self.v1_idx
        s.constraint_type = self.constraint_type
        s.slope_left = self.slope_left
        s.slope_right = self.slope_right
        s.width_left = self.width_left
        s.width_right = self.width_right
        s.auto_left = self.auto_left
        s.auto_right = self.auto_right
        s.backward_left = self.backward_left
        s.backward_right = self.backward_right
        s.type = self.type
        s.enforce_part = self.enforce_part
        return s

    def straight(self, length, t=1):
        s = self.copy
        s.p = self.lerp(t)
        s.v = self.v.normalized() * length
        return s

    def set_offset(self, offset, last=None):
        """
            Offset line and compute intersection point
            between segments
        """
        self.line = self.make_offset(offset, last)

    def offset(self, offset):
        o = self.copy
        o.p += offset * self.cross_z.normalized()
        return o

    @property
    def oposite(self):
        o = self.copy
        o.p += o.v
        o.v = -o.v
        return o

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


class RoofBoundaryGenerator():

    def __init__(self, d, origin=Vector((0, 0))):
        self.parts = d.parts
        self.operation = d.operation
        self.segs = []
        self.length = 0
        self.origin = origin

    def add_part(self, part):

        if len(self.segs) < 1:
            s = None
        else:
            s = self.segs[-1]

        # start a new roof
        if s is None:
            v = part.length * Vector((cos(part.a0), sin(part.a0)))
            s = StraightRoof(self.origin, v)
        else:
            s = s.straight_roof(part.a0, part.length)

        s.type = part.type

        self.segs.append(s)

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
        if self.operation == 'DIFFERENCE':
            side = -1
        else:
            side = 1
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

                v0 = self.segs[i - 1].straight(-side, 1).v.to_3d()
                v1 = f.straight(side, 0).v.to_3d()
                manipulators[0].set_pts([p0, v0, v1])

            # segment length
            manipulators[1].type_key = 'SIZE'
            manipulators[1].prop1_name = "length"
            manipulators[1].set_pts([p0, p1, (side, 0, 0)])

            # snap manipulator, dont change index !
            manipulators[2].set_pts([p0, p1, (side, 0, 0)])
            # dumb segment id
            manipulators[3].set_pts([p0, p1, (side, 0, 0)])

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

    def change_coordsys(self, fromTM, toTM):
        """
            move shape fromTM into toTM coordsys
        """
        dp = (toTM.inverted() * fromTM.translation).to_2d()
        da = toTM.row[1].to_2d().angle_signed(fromTM.row[1].to_2d())
        ca = cos(da)
        sa = sin(da)
        rM = Matrix([
            [ca, -sa],
            [sa, ca]
            ])
        for s in self.segs:
            tp = (rM * s.p0) - s.p0 + dp
            s.rotate(da)
            s.translate(tp)

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

    def get_verts(self, verts, edges):

        # boundary vertices without intersections
        n_boundary = len(self.segs) - 1

        # add boundary verts

        for boundary in self.segs:
            x, y = boundary.line.p0
            verts.append((x, y, 0))

        for i in range(n_boundary):
            edges.append([i, i + 1])

    def get_index(self, index):
        n_segs = len(self.segs)
        if index >= n_segs:
            index -= n_segs
        return index

    def get_seg(self, index):
        idx = self.get_index(index)
        return self.segs[idx]


class RoofSegment():
    def __init__(self, seg, left, right):
        self.seg = seg
        self.left = left
        self.right = right
        self.a0 = 0
        self.reversed = False


class RoofAxisNode():
    def __init__(self):
        # axis segments
        self.segs = []
        self.root = None
        self.center = -1
        # store count of horizontal segs
        self.n_horizontal = 0
        # store count of slopes segs
        self.n_slope = 0

    @property
    def count(self):
        return len(self.segs)

    @property
    def last(self):
        """
            last segments in this node
        """
        return self.segs[-1]

    def left(self, index):
        if index + 1 >= self.count:
            return self.segs[0]
        return self.segs[index + 1]

    def right(self, index):
        return self.segs[index - 1]

    def add(self, a0, reversed, seg, left, right):

        if seg.constraint_type == 'HORIZONTAL':
            self.n_horizontal += 1
        elif seg.constraint_type == 'SLOPE':
            self.n_slope += 1

        s = RoofSegment(seg, left, right)
        s.a0 = a0
        s.reversed = reversed
        if reversed:
            self.root = s
        self.segs.append(s)

    def update_center(self):
        for i, s in enumerate(self.segs):
            if s is self.root:
                self.center = i
                return

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

        # index of root in segs array
        self.update_center()


class RoofPolygon():
    """
        ccw roof pitch boundary
        closed by explicit segment
        handle triangular shape with zero axis length

        mov  <_________________
            |                  /\
            |                  | rot
            |                  |   left     last <> next
            \/_____axis_______>|

     node   <_____axis________      next
            |                  /\
            |                  | rot
            |                  |   right    last <> next
        mov \/________________>|
                           side angle
    """
    def __init__(self, axis, side, fake_axis=None):
        """
            Create a default rectangle
            axis from node to next
            slope float -z for 1 in side direction
            side in ['LEFT', 'RIGHT'] in axis direction

            NOTE:
            when axis length is null (eg: triangular shape)
            use "fake_axis" with a 1 length to handle
            distance from segment
        """
        if side == 'LEFT':
            # slope
            self.slope = axis.slope_left
            # width
            self.width = axis.width_left
            # constraint width
            self.auto_mode = axis.auto_left
        else:
            # slope
            self.slope = axis.slope_right
            # width
            self.width = axis.width_right
            # constraint width
            self.auto_mode = axis.auto_right

        self.side = side
        # backward deps
        self.backward = False
        # pointers to neighboors along axis
        self.last = None
        self.next = None
        self.other_side = None

        # axis segment
        if side == 'RIGHT':
            self.axis = axis.oposite
            self.node_angle = pi / 2
            self.next_angle = -pi / 2
        else:
            self.axis = axis
            self.node_angle = -pi / 2
            self.next_angle = pi / 2

        self.fake_axis = None

        # _axis is either a fake one or real one
        # to prevent further check
        if fake_axis is None:
            self._axis = self.axis
        else:
            if side == 'RIGHT':
                self.fake_axis = fake_axis.oposite
            else:
                self.fake_axis = fake_axis
            self._axis = self.fake_axis

        self.convex = True
        # segments from axis end in ccw order
        # closed by explicit segment
        self.segs = []
        # holes
        self.holes = []
        #
        self.tmin = 0
        self.tmax = 1
        self.dt = 1
        self.ysize = 0
        self.xsize = 0

    def move_node(self, p):
        """
            Move slope point in node side
        """
        if self.side == 'LEFT':
            self.segs[-1].p0 = p
            self.segs[2].p1 = p
        else:
            self.segs[2].p0 = p
            self.segs[1].p1 = p

    def move_next(self, p):
        """
            Move slope point in next side
        """
        if self.side == 'LEFT':
            self.segs[2].p0 = p
            self.segs[1].p1 = p
        else:
            self.segs[-1].p0 = p
            self.segs[2].p1 = p

    def node_link(self):
        angle_90 = round(pi / 2, 4)
        if self.side == 'LEFT':
            idx = -1
        else:
            idx = 1
        da = abs(round(self.node_angle, 4))
        type = "LINK"
        if da < angle_90:
            type += "_VALLEY"
        elif da > angle_90:
            type += "_HIP"
        self.segs[idx].type = type

    def next_link(self):
        angle_90 = round(pi / 2, 4)
        if self.side == 'LEFT':
            idx = 1
        else:
            idx = -1
        da = abs(round(self.next_angle, 4))
        type = "LINK"
        if da < angle_90:
            type += "_VALLEY"
        elif da > angle_90:
            type += "_HIP"
        self.segs[idx].type = type

    def bind(self, last, ccw=False):
        """
            always in axis real direction
        """
        # backward dependancy relative to axis
        if last.backward:
            self.backward = self.side == last.side

        if self.backward:
            self.next = last
            last.last = self
        else:
            self.last = last
            last.next = self

        # width auto
        if self.auto_mode == 'AUTO':
            self.width = last.width
            self.slope = last.slope
        elif self.auto_mode == 'WIDTH':
            self.slope = last.slope * last.width / self.width
        elif self.auto_mode == 'SLOPE':
            self.width = last.width * last.slope / self.slope

        self.make_segments()
        last.make_segments()

        res, p, t = self.segs[2].intersect(last.segs[2])

        if res:
            # dont move anything when no intersection found
            # aka when delta angle == 0
            self.move_node(p)
            if self.side != last.side:
                last.move_node(p)
            else:
                last.move_next(p)

        # Free mode
        # move border
        # and find intersections
        # with sides
        if self.auto_mode == 'ALL':
            s0 = self._axis.offset(-self.width)
            res, p0, t = self.segs[1].intersect(s0)
            self.segs[2].p0 = p0
            self.segs[1].p1 = p0
            res, p1, t = self.segs[-1].intersect(s0)
            self.segs[2].p1 = p1
            self.segs[-1].p0 = p1

        #   /\
        #   |   angle
        #   |____>
        #
        # v1 node -> next
        if self.side == 'LEFT':
            v1 = self._axis.v
        else:
            v1 = -self._axis.v

        if last.side == self.side:
            # contigous, v0 node <- next

            # half angle between segments
            if self.side == 'LEFT':
                v0 = -last._axis.v
            else:
                v0 = last._axis.v
            da = v0.angle_signed(v1)
            if ccw:
                if da < 0:
                    da = 2 * pi + da
            elif da > 0:
                da = da - 2 * pi
            last.next_angle = 0.5 * da
            last.next_link()

        else:
            # alternate v0 node -> next
            # half angle between segments
            if last.side == 'LEFT':
                v0 = last._axis.v
            else:
                v0 = -last._axis.v
            da = v0.angle_signed(v1)
            # angle always ccw
            if ccw:
                if da < 0:
                    da = 2 * pi + da
            elif da > 0:
                da = da - 2 * pi
            last.node_angle = 0.5 * da
            last.node_link()

        self.node_angle = -0.5 * da
        self.node_link()
    
    def get_index(self, index):
        n_segs = len(self.segs)
        if index >= n_segs:
            index -= n_segs
        return index
    
    def next_seg(self, index):
        idx = self.get_index(index + 1)
        return self.segs[idx]
    
    def last_seg(self, index):
        return self.segs[index - 1]
    
    def make_segments(self):
        if len(self.segs) < 1:
            s0 = self._axis
            w = self.width
            s1 = s0.straight(w, 1).rotate(pi / 2)
            s1.type = 'SIDE'
            s3 = s0.straight(w, 0).rotate(pi / 2).oposite
            s3.type = 'SIDE'
            s2 = StraightRoof(s1.p1, s3.p0 - s1.p1)
            s2.type = 'BOTTOM'
            self.segs = [s0, s1, s2, s3]
    
    def move_side(self, pt):
        """
            offset side to point
        """
        s2 = self.segs[2]
        d0, t = self.distance(s2.p0)
        d1, t = self.distance(pt)
        self.slope = self.slope * d0 / d1
        self.segs[2] = s2.offset(d1 - d0)
    
    def propagate_backward(self, pt):
        """
            Propagate slope, keep 2d angle of slope
            Move first point and border
            keep border parallel
            adjust slope
            and next shape
        """
        # distance of p
        # offset side to point
        self.move_side(pt)

        # move verts on node side
        self.move_next(pt)

        if self.side == 'LEFT':
            # move verts on next side
            res, p, t = self.segs[-1].intersect(self.segs[2])
        else:
            # move verts on next side
            res, p, t = self.segs[1].intersect(self.segs[2])
        self.move_node(p)

        if self.next is not None:
            self.next.propagate_backward(p)
    
    def propagate_forward(self, pt):
        """
            Propagate slope, keep 2d angle of slope
            Move first point and border
            keep border parallel
            adjust slope
            and next shape
        """
        # offset side to point
        self.move_side(pt)

        # move verts on node side
        self.move_node(pt)
        if self.side == 'LEFT':
            # move verts on next side
            res, p, t = self.segs[1].intersect(self.segs[2])
        else:
            # move verts on next side
            res, p, t = self.segs[-1].intersect(self.segs[2])
        self.move_next(p)
        if self.next is not None:
            self.next.propagate_forward(p)
    
    def rotate_next_slope(self, a0):
        """
            Rotate next slope part
        """
        if self.side == 'LEFT':
            s0 = self.segs[1].rotate(a0)
            s1 = self.segs[2]
            res, p, t = s1.intersect(s0)
        else:
            s0 = self.segs[2]
            s1 = self.segs[-1]
            res, p, t = s1.oposite.rotate(-a0).intersect(s0)
        s1.p0 = p
        s0.p1 = p

        if self.next is not None:
            # if self.next.auto_mode == 'ALL':
            #    return
            if self.next.backward:
                self.next.propagate_backward(p)
            else:
                self.next.propagate_forward(p)
    
    def rotate_node_slope(self, a0):
        """
            Rotate node slope part
        """
        if self.side == 'LEFT':
            s0 = self.segs[2]
            s1 = self.segs[-1]
            res, p, t = s1.oposite.rotate(-a0).intersect(s0)
        else:
            s0 = self.segs[1].rotate(a0)
            s1 = self.segs[2]
            res, p, t = s1.intersect(s0)
        s1.p0 = p
        s0.p1 = p

        if self.next is not None:
            # if self.next.auto_mode == 'ALL':
            #    return
            if self.next.backward:
                self.next.propagate_backward(p)
            else:
                self.next.propagate_forward(p)
    
    def distance(self, pt):
        """
            distance from axis
        """
        res, d, t = self._axis.point_sur_segment(pt)
        return d, t
    
    def inside(self, pt, segs=None):
        """
            Point inside poly (raycast method)
            support concave polygons
            TODO:
            make s1 angle different than all othr segs
        """
        s1 = Line(pt, Vector((2 * self.xsize, 0.1)))
        counter = 0
        if segs is None:
            segs = self.segs
        for s in segs:
            res, p, t, u = s.intersect_ext(s1)
            if res:
                counter += 1
        return counter % 2 == 1
    
    def altitude(self, pt):
        d, t = self.distance(pt)
        return -d * self.slope
    
    def uv(self, pt):
        d, t = self.distance(pt)
        return ((t - self.tmin) * self.xsize, d)
    
    def intersect(self, seg):
        """
            compute intersections of a segment with boundaries
            segment must start on axis
            return segments inside
        """
        it = []
        for s in self.segs:
            res, p, t, u = seg.intersect_ext(s)
            if res:
                it.append((t, p))
        return it
    
    def is_convex(self):
        n_segs = len(self.segs)
        self.convex = True
        sign = False
        s0 = self.segs[-1]
        for i in range(n_segs):
            s1 = self.segs[i]
            c = s0.v.cross(s1.v)
            if i == 0:
                sign = (c > 0)
            elif sign != (c > 0):
                self.convex = False
                return
            s0 = s1
    
    def iter_segs(self, segs, index):
        n_segs = len(segs)
        if index >= n_segs:
            index -= n_segs
        return index
    
    def follow_until_cross(
            self, 
            f_segs, f_copy, f_start, 
            c_segs, c_copy, c_start, 
            end_point, 
            skip_last, 
            store):
        """
            f_ the line to follow
            c_ the crossing one
            end_point : the point we may reach to end the loop
            skip_last : when starting from intersection skip last from test 
                        to prevent redetection the same intersection
            store : array of segments found in the rootline
        """
        f_nsegs = len(f_segs)
        c_nsegs = len(c_segs)
        for f in range(f_nsegs):
            
            f_idx = self.iter_segs(f_segs, f_start + f)
            f_s = f_segs[f_idx]
            f_c = f_copy[f_idx]
            
            store.append(f_c)
            
            print("append r:%s f_c.p0:%s p1:%s" % (f_idx, f_c.p0, f_c.p1))
            
            # found starting point then exit
            if f_s.p1 == end_point:
                return True, f_start, c_start
            
            # when cross occur, modify crossing segments
            # so it fit with intersection
            for c in range(c_nsegs - skip_last):
                c_idx = self.iter_segs(c_segs, c_start + c)
                c_s = c_segs[c_idx]
                c_c = c_copy[c_idx]
                res, p, u, v = f_s.intersect_ext(c_s)
                if res:
                    f_c.p1 = p
                    c_c.p0 = p
                    c_start = c_idx
                    # a segment only cross another once
                    # so prevent detection of same intersection 
                    f_start = f_idx + 1
                    print("r:%s intersect b:%s at p:%s " % (f_idx, c_idx, p))
                    return False, f_start, c_start
                    
            # on next iterations 
            # check for the whole shape
            skip_last = 0
            
        return False, f_start, c_start
        
    def slice(self, boundary):
        """
            find a point of this pitch inside boundary
            follow pitch until boundary cross, storing each segment
            when cross occur, modify crossing segments to store crossing point
            then follow boundary storing segments until cross pitch again.
            repeat until pitch start point is reached
            
            4 cases:
            1 pitch has point in boundary -> start from this point
            2 boundary hea point in pitch -> start from this point
            3 no points inside -> find first crossing segment
            4 not points inside and no crossing segments 
        """
        
        # keep inside or cut inside
        # keep inside must be CCW
        # cut inside must be CW
        keep_inside = (boundary.operation == 'INTERSECTION')
        
        f_start = -1
        c_start = -1
        f_segs = self.segs
        c_segs = boundary.segs
        f_copy = [s.copy for s in f_segs]
        c_copy = [s.copy for s in c_segs]
        
        # find a point of this pitch inside boundary
        for i, s in enumerate(f_segs):
            res = self.inside(s.p0, c_segs)
            if res == keep_inside:
                f_start = i
                c_start = 0
                end_point = f_segs[f_start].p0
                print("pitch pt inside f_start:%s" % (f_start))
                break
            
        # no pitch point found inside boundary
        # seek for point of boundary inside pitch
        if c_start < 0:
            for i, s in enumerate(c_segs):
                res = self.inside(s.p0)
                if res == keep_inside:
                    c_start = i
                    f_start = 0
                    end_point = c_segs[c_start].p0
                    print("boundary pt inside c_start:%s" % (c_start))
                    # swap boundary / pitch so we start from boundary
                    c_start, f_start = f_start, c_start
                    c_segs, f_segs = f_segs, c_segs
                    c_copy, f_copy = f_copy, c_copy
                    break
            
        # no points found at all
        # seek for a crossing segment 
        if f_start < 0:
            print("no pt inside")
            return
                    
        store = []
        iter = 0
        
        print("slice end_point:%s" % (end_point))
        
        # on first iteration, use all segments
        skip_last = 0
            
        while iter < 10:    
            end, c_start, f_start = self.follow_until_cross(
                                        f_segs, f_copy, f_start, 
                                        c_segs, c_copy, c_start, 
                                        end_point, 
                                        skip_last, 
                                        store)                
            if end:
                break
            # after first cross, skip last segment
            # as it is the intersected one    
            skip_last = 1
            
            # follow other line
            # swap seek / cross
            c_segs, f_segs = f_segs, c_segs
            c_copy, f_copy = f_copy, c_copy
            iter += 1
        
        # points are found, but nothing intersects
        # so use boundary as hole
        if iter < 1:
            self.holes.append(boundary)
            return 
            
        # find "axis"
        a = 0
        for i, s in enumerate(store):
            if s.type == 'AXIS':
                a = i
                break
                
        self.segs = store[a:] + store[:a]
        self._axis = self.segs[0]

        print("types:%s" % ([s.type for s in self.segs]))
        self.is_convex()
        self.limits()

    def merge(self, other):

        raise NotImplementedErrot

    def draw(self, context, z, verts, edges):
        f = len(verts)
        #
        #   0_______1
        #   |_______|
        #   3       2
        verts.extend([(s.p0.x, s.p0.y, z + self.altitude(s.p0)) for s in self.segs])
        n_segs = len(self.segs) - 1
        edges.extend([[f + i, f + i + 1] for i in range(n_segs)])
        edges.append([f + n_segs, f])

        w = 1
        a0 = pi / 2
        if self.backward:
            a = -pi / 4
        else:
            a = pi / 4
        if self.side == 'RIGHT':
            a = -a
        s = self._axis.straight(w, 0.5).rotate(a0 + a)
        f = len(verts)
        p0 = s.p0.to_3d()
        p0.z = z
        p1 = s.p1.to_3d()
        p1.z = z
        verts.extend([p0, p1])
        edges.append([f, f + 1])
    
    def as_string(self):
        """
            Print strips relationships
        """
        if self.backward:
            a0 = degrees(self.next_angle)
            a1 = degrees(self.node_angle)
            dir = "/\\"
            print("%s next %s" % (dir, a0))
        else:
            a0 = degrees(self.node_angle)
            a1 = degrees(self.next_angle)
            dir = "\\/"
            print("%s node %s" % (dir, a0))
        print("%s %s" % (dir, self.side))
        if self.backward:
            print("%s node %s" % (dir, a1))
        else:
            print("%s next %s" % (dir, a1))
        if self.next:
            print("_________")
            self.next.as_string()
        else:
            print("#########")
    
    def limits(self):
        dist = []
        param_t = []
        for s in self.segs:
            res, d, t = self._axis.point_sur_segment(s.p0)
            param_t.append(t)
            dist.append(d)
        self.tmin = min(param_t)
        self.tmax = max(param_t)
        self.dt = self.tmax - self.tmin
        self.ysize = max(dist)
        self.xsize = self._axis.length * self.dt

        
"""
import bmesh
m = C.object.data
[(round(v.co.x, 3), round(v.co.y, 3), round(v.co.z, 3)) for v in m.vertices]
[tuple(p.vertices) for p in m.polygons]

uvs = []
bpy.ops.object.mode_set(mode='EDIT')
bm = bmesh.from_edit_mesh(m)
[tuple(i.index for i in edge.verts) for edge in bm.edges]

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
            v = part.length * Vector((cos(part.a0), sin(part.a0)))
            s = StraightRoof(self.origin, v)
        else:
            s = s.straight_roof(part.a0, part.length)
        # parent segment (root) index is  v0_idx - 1
        s.v0_idx = part.bound_idx
        s.constraint_type = part.constraint_type

        if part.constraint_type == 'SLOPE':
            s.enforce_part = part.enforce_part
        else:
            s.enforce_part = 'AUTO'

        s.angle_0 = part.a0
        s.take_precedence = part.take_precedence
        s.auto_right = part.auto_right
        s.auto_left = part.auto_left
        s.width_left = part.width_left
        s.width_right = part.width_right
        s.slope_left = part.slope_left
        s.slope_right = part.slope_right
        s.type = 'AXIS'
        self.segs.append(s)

    def locate_manipulators(self):
        """

        """
        for i, f in enumerate(self.segs):

            manipulators = self.parts[i].manipulators
            p0 = f.p0.to_3d()
            p0.z = self.z
            p1 = f.p1.to_3d()
            p1.z = self.z
            # angle from last to current segment
            if i > 0:

                manipulators[0].type_key = 'ANGLE'
                v0 = self.segs[f.v0_idx - 1].straight(-1, 1).v.to_3d()
                v1 = f.straight(1, 0).v.to_3d()
                manipulators[0].set_pts([p0, v0, v1])

            # segment length
            manipulators[1].type_key = 'SIZE'
            manipulators[1].prop1_name = "length"
            manipulators[1].set_pts([p0, p1, (1.0, 0, 0)])

            # dumb segment id
            manipulators[2].set_pts([p0, p1, (1, 0, 0)])

            p0 = f.lerp(0.5).to_3d()
            p0.z = self.z
            # size left
            p1 = f.sized_normal(0.5, -self.parts[i].width_left).p1.to_3d()
            p1.z = self.z
            manipulators[3].set_pts([p0, p1, (1, 0, 0)])

            # size right
            p1 = f.sized_normal(0.5, self.parts[i].width_right).p1.to_3d()
            p1.z = self.z
            manipulators[4].set_pts([p0, p1, (-1, 0, 0)])

            # slope left
            n0 = f.sized_normal(0.5, -1)
            p0 = n0.p1.to_3d()
            p0.z = self.z
            p1 = p0.copy()
            p1.z = self.z - self.parts[i].slope_left
            manipulators[5].set_pts([p0, p1, (-1, 0, 0)], normal=n0.v.to_3d())

            # slope right
            n0 = f.sized_normal(0.5, 1)
            p0 = n0.p1.to_3d()
            p0.z = self.z
            p1 = p0.copy()
            p1.z = self.z - self.parts[i].slope_right
            manipulators[6].set_pts([p0, p1, (1, 0, 0)], normal=n0.v.to_3d())

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

    
    def make_roof(self, context, verts, edges):
        """
            Init data structure for possibly multi branched nodes
        """
        z = self.z

        pans = []

        # node are connected segments
        # node
        # (segment idx)
        # (angle from root part > 0 right)
        # (reversed) denote a seg connected by p1
        #            "root" of node
        nodes = [RoofAxisNode() for s in range(len(self.segs) + 1)]

        # Init width on seg 0
        s0 = self.segs[0]
        if self.parts[0].auto_left in {'AUTO', 'SLOPE'}:
            s0.width_left = self.width_left
        if self.parts[0].auto_right in {'AUTO', 'SLOPE'}:
            s0.width_right = self.width_right
        if self.parts[0].auto_left in {'AUTO', 'WIDTH'}:
            s0.slope_left = self.slope_left
        if self.parts[0].auto_left in {'AUTO', 'WIDTH'}:
            s0.slope_right = self.slope_right

        # make nodes with HORIZONTAL constraints
        for idx, s in enumerate(self.segs):
            s.v1_idx = idx + 1
            if s.constraint_type == 'HORIZONTAL':
                left = RoofPolygon(s, 'LEFT')
                right = RoofPolygon(s, 'RIGHT')
                left.other_side = right
                right.other_side = left
                rs = RoofSegment(s, left, right)
                pans.append(rs)
                nodes[s.v0_idx].add(s.angle_0, False, s, left, right)
                nodes[s.v1_idx].add(-pi, True, s, left, right)

        # set first node root
        # so regular sort does work
        nodes[0].root = nodes[0].segs[0]

        # Propagate slope and width
        # on node basis along axis
        # bi-direction Radial around node
        # from left and right to center
        # contigous -> same
        # T: and (x % 2 == 1)
        # First one take precedence over others
        # others inherit from side
        #
        #         l / rb    l = left
        #          3        r = right
        #   l _1_ /         b = backward
        #   r     \
        #          2
        #          r\ l
        #
        # X: rigth one r left one l (x % 2 == 0)
        # inherits from side
        #
        #    l 3 lb         l = left
        # l__1_|_2_l        r = right
        # r    |   r        b = backward -> propagate in reverse axis direction
        #    r 4 rb
        #
        # for idx, node in enumerate(nodes):
        #    print("idx:%s node:%s" % (idx, node.root))

        for idx, node in enumerate(nodes):

            node.sort()

            nb_segs = node.count

            if node.root is None:
                continue

            left = node.root.left
            right = node.root.right

            # basic one single node
            if nb_segs < 2:
                left.make_segments()
                right.make_segments()
                continue

            # get "root" slope and width
            l_bind = left
            r_bind = right

            # simple case: 2 contigous segments
            if nb_segs == 2:
                s = node.last
                s.right.bind(r_bind, ccw=False)
                s.left.bind(l_bind, ccw=True)
                continue

            # More than 2 segments, uneven distribution
            if nb_segs % 2 == 1:
                # find wich child does take precedence
                # first one on rootline (arbitrary)
                center = (nb_segs - 1) / 2
            else:
                # even distribution
                center = nb_segs / 2

            # user defined precedence if any
            for i, s in enumerate(node.segs):
                if s.seg.take_precedence:
                    center = i
                    break

            # bind right side to center
            for i, s in enumerate(node.segs):
                # skip axis
                if i > 0:
                    if i < center:
                        # right contigous with last
                        s.right.bind(r_bind, ccw=False)

                        # next bind to left
                        r_bind = s.left

                        # left backward, not bound
                        # so setup width and slope
                        if s.left.auto_mode in {'AUTO', 'WIDTH'}:
                            s.left.slope = right.slope
                        if s.left.auto_mode in {'AUTO', 'SLOPE'}:
                            s.left.width = right.width
                        s.left.backward = True
                    else:
                        # right bound to last
                        s.right.bind(r_bind, ccw=False)
                        break

            # bind left side to center
            for i, s in enumerate(reversed(node.segs)):
                # skip axis
                if i < nb_segs - center - 1:
                    # left contigous with last
                    s.left.bind(l_bind, ccw=True)
                    # next bind to right
                    l_bind = s.right
                    # right backward, not bound
                    # so setup width and slope
                    if s.right.auto_mode in {'AUTO', 'WIDTH'}:
                        s.right.slope = left.slope
                    if s.right.auto_mode in {'AUTO', 'SLOPE'}:
                        s.right.width = left.width
                    s.right.backward = True
                else:
                    # right bound to last
                    s.left.bind(l_bind, ccw=True)
                    break

        # slope constraints allowed between segments
        # multiple (up to 2) on start and end
        # single between others
        #
        #    2 slope            2 slope           2 slope
        #     |                  |                 |
        #     |______section_1___|___section_2_____|
        #     |                  |                 |
        #     |                  |                 |
        #    multiple           single            multiple

        # add slopes constraints to nodes
        for i, s in enumerate(self.segs):
            if s.constraint_type == 'SLOPE':
                nodes[s.v0_idx].add(s.angle_0, False, s, None, None)

        # sort nodes, remove duplicate slopes between
        # horizontal, keeping only first one
        for idx, node in enumerate(nodes):
            to_remove = []
            node.sort()
            # remove dup between all
            # but start / end nodes
            if node.n_horizontal > 1:
                last = None
                for i, s in enumerate(node.segs):
                    if s.seg.constraint_type == last:
                        if s.seg.constraint_type == 'SLOPE':
                            to_remove.append(i)
                    last = s.seg.constraint_type
                for i in reversed(to_remove):
                    node.segs.pop(i)
                node.update_center()

        for idx, node in enumerate(nodes):

            root = node.root

            # a node may contain many slopes
            # 2 * (part starting from node - 1)
            #
            #        s0
            # root 0 |_______
            #        |
            #        s1
            #
            #               s1
            # root   _______|
            #               |
            #               s0
            #
            #       s3  3  s2
            #     l   \l|r/ l
            # root  ___\|/___ 2
            #     r    /|\  r
            #         /r|l\
            #       s0  1  s1
            #
            #        s2  s1=slope
            #        |r /
            #        | / l
            #        |/____s
            #
            # root to first child -> equal side
            # any other childs -> oposite sides

            if node.n_horizontal == 1:
                # slopes at start or end of segment
                if node.n_slope > 0:
                    # node has user def slope
                    s = node.root
                    s0 = node.left(node.center)
                    a0 = s0.seg.delta_angle(s.seg)
                    if node.root.reversed:
                        # slope at end of segment
                        # first one is right or left
                        if a0 < 0:
                            # right side
                            res, p, t = s0.seg.intersect(s.right.segs[2])
                            s.right.segs[-1].p0 = p
                            s.right.segs[2].p1 = p
                        else:
                            # left side
                            res, p, t = s0.seg.intersect(s.left.segs[2])
                            s.left.segs[1].p1 = p
                            s.left.segs[2].p0 = p
                        if node.n_slope > 1:
                            a1 = s1.seg.delta_angle(s.seg)
                            # last one must be left
                            s1 = node.right(node.center)
                            # both slopes on same side:
                            # skip this one
                            if a0 < 0 and a1 > 0:
                                # left side
                                res, p, t = s1.seg.intersect(s.left.segs[2])
                                s.left.segs[1].p1 = p
                                s.left.segs[2].p0 = p
                    else:
                        # slope at start of segment
                        if a0 < 0:
                            # right side
                            res, p, t = s0.seg.intersect(s.right.segs[2])
                            s.right.segs[1].p1 = p
                            s.right.segs[2].p0 = p
                        else:
                            # left side
                            res, p, t = s0.seg.intersect(s.left.segs[2])
                            s.left.segs[-1].p0 = p
                            s.left.segs[2].p1 = p
                        if node.n_slope > 1:
                            a1 = s1.seg.delta_angle(s.seg)
                            # last one must be right
                            s1 = node.right(node.center)
                            # both slopes on same side:
                            # skip this one
                            if a0 < 0 and a1 > 0:
                                # right side
                                res, p, t = s1.seg.intersect(s.right.segs[2])
                                s.right.segs[1].p1 = p
                                s.right.segs[2].p0 = p

            else:
                # slopes between segments
                for i, s0 in enumerate(node.segs):
                    s1 = node.left(i)
                    s2 = node.left(i + 1)

                    if s1.seg.constraint_type == 'SLOPE':

                        # 3 cases:
                        # s0 is root contigous -> sides are same
                        # s2 is root contigous -> sides are same
                        # back to back -> sides are not same

                        if s0.reversed:
                            # contigous right / right
                            # 2 cases
                            # right is backward
                            # right is forward
                            if s2.right.backward:
                                # s0 depends on s2
                                main = s2.right
                                v = main.segs[1].v
                            else:
                                # s2 depends on s0
                                main = s0.right
                                v = -main.segs[-1].v
                            res, p, t = s1.seg.intersect(main.segs[2])
                            if res:
                                # slope vector
                                dp = p - s1.seg.p0
                                a0 = dp.angle_signed(v)
                                if s2.right.backward:
                                    main.rotate_node_slope(a0)
                                else:
                                    main.rotate_next_slope(-a0)
                        elif s2.reversed:
                            # contigous left / left
                            # 2 cases
                            # left is backward
                            # left is forward
                            if s0.left.backward:
                                # s0 depends on s2
                                main = s0.left
                                v = -main.segs[-1].v
                            else:
                                # s2 depends on s0
                                main = s2.left
                                v = main.segs[1].v
                            res, p, t = s1.seg.intersect(main.segs[2])
                            if res:
                                # slope vector
                                dp = p - s1.seg.p0
                                a0 = dp.angle_signed(v)
                                if s0.left.backward:
                                    main.rotate_node_slope(-a0)
                                else:
                                    main.rotate_next_slope(a0)
                        else:
                            # back left / right
                            # 2 cases
                            # left is backward
                            # left is forward
                            if s0.left.backward:
                                # s2 depends on s0
                                main = s0.left
                                v = -main.segs[-1].v
                            else:
                                # s0 depends on s2
                                main = s2.right
                                v = main.segs[1].v

                            res, p, t = s1.seg.intersect(main.segs[2])
                            if res:
                                # slope vector
                                dp = p - s1.seg.p0
                                a0 = dp.angle_signed(v)
                                if s0.left.backward:
                                    main.rotate_node_slope(-a0)
                                else:
                                    main.rotate_node_slope(a0)

        # make flat array
        self.pans = []
        for pan in pans:
            self.pans.extend([pan.left, pan.right])

        # merge contigous with 0 angle diff
        to_remove = []
        for i, pan in enumerate(self.pans):
            next = pan.next
            if next is not None:
                # same side
                if next.side == pan.side:
                    if round(next._axis.delta_angle(pan._axis), 4) == 0:
                        next.last = pan.last
                        to_remove.append(i)
                        if pan.side == 'LEFT':
                            next.slope = pan.slope
                            if next.backward:
                                next.next_angle = pan.node_angle
                                next._axis.p1 = pan._axis.p1
                                next.segs[1] = pan.segs[1]
                                next.segs[2].p0 = pan.segs[2].p0
                            else:
                                next.node_angle = pan.node_angle
                                next._axis.p0 = pan._axis.p0
                                next.segs[-1] = pan.segs[-1]
                                next.segs[2].p1 = pan.segs[2].p1
                        else:
                            if next.backward:
                                next.next_angle = pan.node_angle
                                next._axis.p0 = pan._axis.p0
                                next.segs[-1] = pan.segs[-1]
                                next.segs[2].p1 = pan.segs[2].p1
                            else:
                                next.node_angle = pan.node_angle
                                next._axis.p1 = pan._axis.p1
                                next.segs[1] = pan.segs[1]
                                next.segs[2].p0 = pan.segs[2].p0

        for i in reversed(to_remove):
            self.pans.pop(i)

        # compute limits
        for pan in self.pans:
            pan.limits()

        """
        for pan in self.pans:
            if pan.last is None:
                pan.as_string()
        """
        for s in self.segs:
            if s.constraint_type == 'SLOPE':
                f = len(verts)
                p0 = s.p0.to_3d()
                p0.z = self.z
                p1 = s.p1.to_3d()
                p1.z = self.z
                verts.extend([p0, p1])
                edges.append([f, f + 1])

        #   Include user def boundarys if any
        #
        # node ________axis_______ next
        #   |                      |
        # s0|       part           | s1 borders segments
        #   |                      |
        #
        #       boundary segments
        #
        #
        #  Boundary may intersect borders or axis
        #  Boundary does cut part in 2 spaces
        #  - inside and outside
        #  - inside part is defined by boundary itself
        #
        #
        return

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

    def lambris(self, context, o, d):

        idmat = 0
        lambris_height = 0.02
        alt = self.z - lambris_height
        for pan in self.pans:
            
            verts = []
            faces = []
            matids = []
            uvs = []
            
            f = len(verts)
            verts.extend([(s.p0.x, s.p0.y, alt + pan.altitude(s.p0)) for s in pan.segs])
            uvs.append([pan.uv(s.p0) for s in pan.segs])
            n_segs = len(pan.segs)
            face = [f + i for i in range(n_segs)]
            faces.append(face)
            matids.append(idmat)
            
            bm = bmed.buildmesh(
                context, o, verts, faces, matids=matids, uvs=uvs,
                weld=False, clean=False, auto_smooth=True, temporary=True)
            
            for hole in pan.holes:
                for seg in hole.segs:
                    self.bissect(bm, seg.p0.to_3d(), seg.cross_z.to_3d(), clear_outer=False)
                # when hole segs are found clear parts inside hole
                f_geom = [f for f in bm.faces 
                    if pan.inside(
                        f.calc_center_median().to_2d(), 
                        segs=hole.segs)]
                bmesh.ops.delete(bm, geom=f_geom, context=5)
            
            geom = bm.faces[:]
            verts = bm.verts[:]
            bmesh.ops.solidify(bm, geom=geom, thickness=0.0001)
            bmesh.ops.translate(bm, vec=Vector((0, 0, lambris_height)), space=o.matrix_world, verts=verts)

            # merge with object
            bmed.bmesh_join(context, o, [bm], normal_update=True)
            
    def couverture(self, context, o, d):

        idmat = 7
        rand = 3

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

        dx, dy = d.tile_space_x, d.tile_space_y

        for pan in self.pans:

            seg = pan._axis
            # compute base matrix top left of face
            vx = -seg.v.normalized().to_3d()
            vy = Vector((-vx.y, vx.x, pan.slope)).normalized()
            vz = vx.cross(vy)

            x0, y0 = seg.lerp(pan.tmax)
            z0 = self.z + d.tile_altitude
            ysize_2d = (d.tile_border + pan.ysize)
            space_x = pan.xsize + 2 * d.tile_side
            space_y = ysize_2d * sqrt(1 + pan.slope * pan.slope)
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

            # clean outer on convex parts
            # pan.convex = False
            remove = pan.convex

            for seg in pan.segs:
                if seg.type  == 'AXIS':
                    self.bissect(bm, seg.p1.to_3d(), seg.cross_z.to_3d(), clear_outer=remove)
                elif seg.type == 'BOTTOM':
                    s0 = seg.offset(d.tile_border)
                    dz = pan.altitude(s0.p0)
                    vx = s0.v.to_3d()
                    vx.z = pan.altitude(s0.p1) - dz
                    vy = vz.cross(vx.normalized())
                    x, y = s0.p0
                    z = z0 + dz
                    self.bissect(bm, Vector((x, y, z)), -vy, clear_outer=remove)
                elif seg.type == 'SIDE':
                    p0 = seg.p0 + seg.cross_z.normalized() * d.tile_side
                    self.bissect(bm, p0.to_3d(), seg.cross_z.to_3d(), clear_outer=remove)
                elif seg.type == 'LINK_VALLEY':
                    p0 = seg.p0 - seg.cross_z.normalized() * d.tile_couloir
                    self.bissect(bm, p0.to_3d(), seg.cross_z.to_3d(), clear_outer=remove)
                elif seg.type in {'LINK_HIP', 'LINK'}:
                    self.bissect(bm, seg.p0.to_3d(), seg.cross_z.to_3d(), clear_outer=remove)
            
            for hole in pan.holes:
                for seg in hole.segs:
                    self.bissect(bm, seg.p0.to_3d(), seg.cross_z.to_3d(), clear_outer=False)
                # when hole segs are found clear parts inside hole
                f_geom = [f for f in bm.faces 
                    if pan.inside(
                        f.calc_center_median().to_2d(), 
                        segs=hole.segs)]
                bmesh.ops.delete(bm, geom=f_geom, context=5)
             
            # when not convex, select and remove outer parts
            if not pan.convex:
                """
                /* del "context" slot values, used for operator too */
                enum {
                    DEL_VERTS = 1,
                    DEL_EDGES,
                    DEL_ONLYFACES,
                    DEL_EDGESFACES,
                    DEL_FACES,
                    /* A version of 'DEL_FACES' that keeps edges on face boundaries,
                     * allowing the surrounding edge-loop to be kept from removed face regions. */
                    DEL_FACES_KEEP_BOUNDARY,
                    DEL_ONLYTAGGED
                };
                """
                f_geom = [f for f in bm.faces if not pan.inside(f.calc_center_median().to_2d())]
                bmesh.ops.delete(bm, geom=f_geom, context=5)

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
        bpy.ops.object.mode_set(mode='OBJECT')    
            
    def rake(self, d, verts, faces, edges, matids, uvs):

        #####################
        # Vire-vents
        #####################

        idmat = 1
        for pan in self.pans:
            for s0 in pan.segs:
                if s0.type == 'SIDE':
                    if s0.enforce_part is 'AUTO':
                        f = len(verts)
                        s1 = s0.offset(d.rake_offset)
                        s0 = s0.offset(d.rake_offset - d.rake_width)
                        x0, y0 = s0.p0
                        x1, y1 = s1.p0
                        x2, y2 = s1.p1
                        x3, y3 = s0.p1
                        z0 = self.z + d.rake_altitude + pan.altitude(s0.p0)
                        z1 = self.z + d.rake_altitude + pan.altitude(s0.p1)
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

        #####################
        # Larmiers
        #####################

        idmat = 2
        for pan in self.pans:
            for i, s in enumerate(pan.segs):
                if s.type == 'BOTTOM':
                    f = len(verts)
                    s0 = s.offset(d.facia_width)
                    s1 = pan.last_seg(i)
                    s2 = pan.next_seg(i)

                    if s1.type == 'BOTTOM':
                        s1 = s1.offset(d.facia_width)
                    res, p0, t = s1.intersect(s0)
                    if not res:
                        p0 = s0.p0

                    if s2.type == 'BOTTOM':
                        s2 = s2.offset(d.facia_width)
                    res, p1, t = s2.intersect(s0)
                    if not res:
                        p1 = s0.p1

                    x0, y0 = s.p0
                    x1, y1 = p0
                    x2, y2 = p1
                    x3, y3 = s.p1
                    z0 = self.z + d.facia_altitude + pan.altitude(s.p0)
                    z1 = self.z + d.facia_altitude + pan.altitude(s.p1)
                    verts.extend([
                        (x0, y0, z0),
                        (x1, y1, z0),
                        (x2, y2, z1),
                        (x3, y3, z1),
                    ])
                    z0 -= d.facia_height
                    z1 -= d.facia_height
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

    def gutter(self, d, verts, faces, edges, matids, uvs):

        #####################
        # Chenaux
        #####################

        idmat = 5

        # caps at start and end
        if d.gutter_segs % 2 == 1:
            n_faces = int((d.gutter_segs - 1) / 2)
        else:
            n_faces = int((d.gutter_segs / 2) - 1)

        df = 2 * d.gutter_segs + 1

        for pan in self.pans:
            for i, s in enumerate(pan.segs):
                if s.type == 'BOTTOM':
                    f = len(verts)
                    s2 = pan.last_seg(i)
                    s3 = pan.next_seg(i)
                    d2 = d.gutter_dist + 0.5 * d.gutter_width
                    s4 = s.offset(-d2)

                    a0 = pan.node_angle
                    a1 = pan.next_angle

                    if s2.type == 'BOTTOM':
                        s2 = s2.offset(-d2)
                        res, p, t = s2.intersect(s4)
                        if pan.side == 'LEFT':
                            a1 = -s2.v.angle_signed(p - s.p0)
                        else:
                            a0 = -s2.v.angle_signed(p - s.p0)

                    if s3.type == 'BOTTOM':
                        s3 = s3.offset(-d2)
                        res, p, t = s3.intersect(s4)
                        if pan.side == 'LEFT':
                            a0 = -s3.v.angle_signed(s.p1 - p)
                        else:
                            a1 = -s3.v.angle_signed(s.p1 - p)
                    
                    
                    if pan.side == 'LEFT':
                        s2, s3 = s3, s2
                        s0 = s.straight(-1, 1).rotate(-a0)
                        s1 = s.straight(1, 0).rotate(-a1)
                    else:
                        s0 = s.straight(1, 0).rotate(-a0)
                        s1 = s.straight(-1, 1).rotate(-a1)

                    if abs(a0) > 0:
                        scale_0 = min(3, 1 / sin(abs(a0)))
                    else:
                        scale_0 = 1


                    if abs(a1) > 0:
                        scale_1 = min(3, 1 / sin(abs(a1)))
                    else:
                        scale_1 = 1

                    # print("scale_0:%s scale_1:%s" % (scale_0, scale_1))

                    zt = self.z + d.facia_altitude + pan.altitude(s0.p0)
                    z0 = self.z + d.gutter_alt + pan.altitude(s0.p0)
                    z1 = z0 - 0.5 * d.gutter_width
                    z2 = z1 - 0.5 * d.gutter_width
                    z3 = z1 - 0.5 * d.gutter_boudin
                    dz0 = z2 - z1
                    dz1 = z3 - z1

                    tt = scale_0 * d.facia_width
                    t0 = scale_0 * d.gutter_dist
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

                    zt = self.z + d.facia_altitude + pan.altitude(s1.p0)
                    z0 = self.z + d.gutter_alt + pan.altitude(s1.p0)
                    z1 = z0 - 0.5 * d.gutter_width
                    z2 = z1 - 0.5 * d.gutter_width
                    z3 = z1 - 0.5 * d.gutter_boudin
                    dz0 = z2 - z1
                    dz1 = z3 - z1
                    tt = scale_1 * d.facia_width
                    t0 = scale_1 * d.gutter_dist
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
                    if s2.type == 'SIDE':

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
                    if s3.type == 'SIDE':

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

        for pan in self.pans:
            for i, s in enumerate(pan.segs):

                if s.type == 'AXIS':

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

                    s2 = s.offset(-0.5 * d.beam_width)

                    # offset from roof border
                    s0 = pan.last_seg(i)
                    s1 = pan.next_seg(i)
                    t0 = 0
                    t1 = 1

                    if s0.type == 'SIDE':
                        s0 = s0.offset(d.beam_offset)
                        t0 = -d.beam_offset / s.length

                    if s1.type == 'SIDE':
                        s1 = s1.offset(d.beam_offset)
                        t1 = 1 + d.beam_offset / s.length

                    res, p0, t = s2.intersect(s0)
                    if not res:
                        continue
                    res, p1, t = s2.intersect(s1)
                    if not res:
                        continue

                    x0, y0 = p0
                    x1, y1 = s.lerp(t0)
                    x2, y2 = p1
                    x3, y3 = s.lerp(t1)
                    z0 = self.z + d.beam_alt + pan.altitude(p0)
                    z1 = z0 - d.beam_height
                    z2 = self.z + d.beam_alt + pan.altitude(p1)
                    z3 = z2 - d.beam_height
                    verts.extend([
                        (x0, y0, z0),
                        (x1, y1, z0),
                        (x2, y2, z2),
                        (x3, y3, z2),
                        (x0, y0, z1),
                        (x1, y1, z1),
                        (x2, y2, z3),
                        (x3, y3, z3),
                    ])
                    if s0.type == 'SIDE':
                        faces.append((f, f + 4, f + 5, f + 1))
                        uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                        matids.append(idmat)
                    if s1.type == 'SIDE':
                        faces.append((f + 2, f + 3, f + 7, f + 6))
                        uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
                        matids.append(idmat)

                    faces.extend([
                        # internal side
                        # (f + 1, f + 5, f + 7, f + 3),
                        # external side
                        (f + 2, f + 6, f + 4, f),
                        # top
                        (f, f + 1, f + 3, f + 2),
                        # bottom
                        (f + 5, f + 4, f + 6, f + 7)
                    ])
                    matids.extend([
                        idmat, idmat, idmat
                    ])
                    uvs.extend([
                        [(0, 0), (1, 0), (1, 1), (0, 1)],
                        [(0, 0), (1, 0), (1, 1), (0, 1)],
                        [(0, 0), (1, 0), (1, 1), (0, 1)]
                    ])
    
    def intersect_chevron(self, line, pan, alt, verts):
        u0 = 2
        p1 = line.p1
        for i, seg in enumerate(pan.segs):
            if i == 0:
                continue
            res, p, u, v = line.intersect_ext(seg)
            if res and u < u0:
                p1 = p
                u0 = u
        
        z = self.z + alt + pan.altitude(p1)
        x, y = p1
        verts.append((x, y, z))
        
    def rafter(self, context, o, d):

        idmat = 4

        # Rafters / Chevrons
        start = max(0.001 + 0.5 * d.rafter_width, d.rafter_start)
        chevrons_alt = -0.001
        # build temp bmesh and bissect
             
        for pan in self.pans:
            tmin, tmax, ysize = pan.tmin, pan.tmax, pan.ysize
            
            f = 0
            
            verts = []
            faces = []
            matids = []
            uvs = []
            alt = d.rafter_alt - d.rafter_height
            seg = pan._axis

            t0 = tmin + (start - 0.5 * d.rafter_width) / seg.length
            t1 = tmin + (start + 0.5 * d.rafter_width) / seg.length

            tx = start / seg.length
            dt = d.rafter_spacing / seg.length

            n_items = max(1, round((tmax - tmin) / dt, 0))

            dt = ((tmax - tmin) - 2 * tx) / n_items

            for j in range(int(n_items) + 1):
                n0 = seg.sized_normal(t1 + j * dt, - ysize)
                n1 = seg.sized_normal(t0 + j * dt, - ysize)
                f = len(verts)
                # verts start from axis
                right_before = t1 + j * dt > 0 and t1 + j * dt < 1
                left_before = t0 + j * dt > 0 and t0 + j * dt < 1

                if right_before:
                    z = self.z + alt + pan.altitude(n0.p0)
                    x, y = n0.p0
                    verts.append((x, y, z))
                    
                self.intersect_chevron(n0, pan, alt, verts)

                if not right_before:
                    z = self.z + alt + pan.altitude(n0.p1)
                    x, y = n0.p1
                    verts.append((x, y, z))
                    
                if left_before:
                    z = self.z + alt + pan.altitude(n1.p0)
                    x, y = n1.p0
                    verts.append((x, y, z))
                    
                self.intersect_chevron(n1, pan, alt, verts)

                if not left_before:
                    z = self.z + alt + pan.altitude(n1.p1)
                    x, y = n1.p1
                    verts.append((x, y, z))
                    
                faces.append((f, f + 1, f + 3, f + 2))
                matids.append(idmat)
                uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
            
            bm = bmed.buildmesh(
                context, o, verts, faces, matids=matids, uvs=uvs,
                weld=False, clean=False, auto_smooth=True, temporary=True)
            
            for hole in pan.holes:
                for seg in hole.segs:
                    self.bissect(bm, seg.p0.to_3d(), seg.cross_z.to_3d(), clear_outer=False)
                # when hole segs are found clear parts inside hole
                f_geom = [f for f in bm.faces 
                    if pan.inside(
                        f.calc_center_median().to_2d(), 
                        segs=hole.segs)]
                bmesh.ops.delete(bm, geom=f_geom, context=5)
            
            geom = bm.faces[:]
            verts = bm.verts[:]
            bmesh.ops.solidify(bm, geom=geom, thickness=0.0001)
            bmesh.ops.translate(bm, vec=Vector((0, 0, d.rafter_height)), space=o.matrix_world, verts=verts)

            # merge with object
            bmed.bmesh_join(context, o, [bm], normal_update=True)
            
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

        for pan in self.pans:
            for i, s in enumerate(pan.segs):
                if 'LINK' in s.type and d.beam_sec_enable:
                    # poutre inside
                    f = len(verts)
                    s0 = s.offset(0.5 * d.beam_sec_width)
                    s1 = s.offset(-0.5 * d.beam_sec_width)
                    x0, y0 = s0.p0
                    x1, y1 = s0.p1
                    x2, y2 = s1.p0
                    x3, y3 = s1.p1
                    z0 = self.z + d.beam_sec_alt + pan.altitude(s.p0)
                    z1 = z0 - d.beam_sec_height
                    z2 = self.z + d.beam_sec_alt + pan.altitude(s.p1)
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
                if s.type == 'LINK_HIP':
                     if d.hip_enable:
                        # s0 is bottom seg
                        s0 = pan.last_seg(i)
                        if s0.type == 'AXIS':
                            s0 = pan.next_seg(i)
                            s2 = s
                        else:
                            continue

                        s1 = s0.offset(d.tile_border)
                        res, p1, tmax = s2.intersect(s1)

                        vx = s2.v.to_3d()
                        vx.z = -pan.ysize * pan.slope
                        vx = vx.normalized()
                        vy = s2.cross_z.normalized().to_3d()
                        vz = vx.cross(-vy)
                        x0, y0 = s2.p0 + vx.to_2d().normalized() * 0.3 * d.hip_size_y
                        z0 = self.z + d.hip_alt
                        tM = Matrix([
                            [vx.x, vy.x, vz.x, x0],
                            [vx.y, vy.y, vz.y, y0],
                            [vx.z, vy.z, vz.z, z0],
                            [0, 0, 0, 1]
                        ])

                        space_2d = tmax * pan.ysize - 0.1 * d.hip_size_y
                        space_z = pan.altitude(p1)
                        space_x = sqrt(space_2d * space_2d + space_z * space_z)
                        n_x = 1 + int(space_x / d.hip_space_x)
                        dx = space_x / n_x
                        x0 = 0.5 * dx

                        t_verts = [p for p in t_pts]

                        # apply slope
                        dz = 0
                        #t_verts[i].y * cos(abs(da))

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
                elif s.type == 'LINK_VALLEY':
                    if d.valley_enable:
                        f = len(verts)
                        s0 = s.offset(-2 * d.tile_couloir)
                        s1 = pan._axis
                        # s2 = border
                        s2 = pan.last_seg(i)
                        if s2.type == 'AXIS':
                            s3 = s
                            s2 = pan.next_seg(i)
                        else:
                            s3 = s.oposite
                        res, p0, t = s0.intersect(s1)
                        res, p1, t = s0.intersect(s2)

                        x0, y0 = s3.p0
                        x1, y1 = p0
                        x2, y2 = p1
                        x3, y3 = s3.p1
                        z0 = self.z + d.valley_altitude + pan.altitude(s3.p0)
                        z1 = self.z + d.valley_altitude + pan.altitude(s3.p1)

                        verts.extend([
                            (x0, y0, z0),
                            (x1, y1, z0),
                            (x2, y2, z1),
                            (x3, y3, z1),
                        ])
                        faces.extend([
                            (f, f + 1, f + 2, f + 3)
                        ])
                        matids.extend([
                            idmat_valley
                            ])
                        uvs.extend([
                            [(0, 0), (1, 0), (1, 1), (0, 1)]
                        ])
                elif s.type == 'AXIS' and d.hip_enable and pan.side == 'LEFT':

                    tmin = 0
                    tmax = 1
                    s0 = pan.last_seg(i)
                    if s0.type == 'SIDE':
                        tmin = 0 - d.tile_side / s.length
                    s1 = pan.next_seg(i)

                    if s1.type == 'SIDE':
                        tmax = 1 + d.tile_side / s.length

                    # print("tmin:%s tmax:%s" % (tmin, tmax))
                    ####################
                    # Faitiere
                    ####################

                    f = len(verts)
                    s_len = (tmax - tmin) * s.length
                    n_obj = 1 + int(s_len / d.hip_space_x)
                    dx = s_len / n_obj
                    x0 = 0.5 * dx
                    v = s.v.normalized()
                    p0 = s.lerp(tmin)
                    tM = Matrix([
                        [v.x, v.y, 0, p0.x],
                        [v.y, -v.x, 0, p0.y],
                        [0, 0, 1, self.z + d.hip_alt],
                        [0, 0, 0, 1]
                    ])
                    t_verts = [p.copy() for p in t_pts]

                    # apply slope
                    for i in t_left:
                        t_verts[i].z += t_verts[i].y * (pan.other_side.slope - d.tile_size_z / d.tile_size_y)
                    for i in t_right:
                        t_verts[i].z -= t_verts[i].y * (pan.slope - d.tile_size_z / d.tile_size_y)

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

    def change_coordsys(self, fromTM, toTM):
        """
            move shape fromTM into toTM coordsys
        """
        dp = (toTM.inverted() * fromTM.translation).to_2d()
        da = toTM.row[1].to_2d().angle_signed(fromTM.row[1].to_2d())
        ca = cos(da)
        sa = sin(da)
        rM = Matrix([
            [ca, -sa],
            [sa, ca]
            ])
        for s in self.segs:
            tp = (rM * s.p0) - s.p0 + dp
            s.rotate(da)
            s.translate(tp)

    def t_partition(self, array, begin, end):
        pivot = begin
        for i in range(begin + 1, end + 1):
            # wall idx
            if array[i][0] < array[begin][0]:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
        array[pivot], array[begin] = array[begin], array[pivot]
        return pivot

    def sort_t(self, array, begin=0, end=None):
        # print("sort_child")
        if end is None:
            end = len(array) - 1

        def _quicksort(array, begin, end):
            if begin >= end:
                return
            pivot = self.t_partition(array, begin, end)
            _quicksort(array, begin, pivot - 1)
            _quicksort(array, pivot + 1, end)
        return _quicksort(array, begin, end)

    def make_wall_fit(self, context, o, wall, inside):
        wd = wall.data.archipack_wall2[0]
        wg = wd.get_generator()
        z0 = self.z - wd.z

        segs = self.all_segs
        nodes = self.nodes

        # wg in roof coordsys
        wg.change_coordsys(wall.matrix_world, o.matrix_world)

        if inside:
            # fit inside
            offset = -0.5 * (1 - wd.x_offset) * wd.width
        else:
            # fit outside
            offset = 0

        wg.set_offset(offset)

        wall_t = [[] for w in wg.segs]

        for idx, node in enumerate(nodes):

            # find next node in segment
            for i, s in enumerate(node.segs):

                seg = segs[s.idx]

                if s.reversed:
                    continue

                if nodes[seg.v1_idx].count > 1:
                    next = nodes[seg.v1_idx]

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

                        for widx, wseg in enumerate(wg.segs):
                            # walls segment in roof coordsys
                            p0 = wseg.line.p0
                            p1 = wseg.line.p1
                            # side of wall segment vs axis
                            res, d0, t0 = seg.point_sur_segment(p0)
                            res, d1, t1 = seg.point_sur_segment(p1)

                            # print("w:%s s:%s d0:%s d1:%s" % (widx, idx, d0, d1))

                            # side of axis
                            if d0 > 0:
                                slope = seg.slope_left
                                width = seg.width_left
                                s0 = segs[ir0]
                                s1 = segs[ir1]

                            elif d0 <= 0:
                                slope = seg.slope_right
                                width = seg.width_right
                                s0 = segs[il0]
                                s1 = segs[il1]

                            # check if wall segment is between slopes
                            res, d2, t2 = s0.point_sur_segment(p0)
                            res, d3, t3 = s1.point_sur_segment(p0)

                            # p0 between slopes
                            if d2 > 0 and d3 < 0:
                                z = z0 - abs(d0) * slope
                                wall_t[widx].append((0, z, widx, idx, "p0"))

                            # slope segments intersection check
                            res, p, t, v = wseg.line.intersect_ext(s0)
                            if res and t > 0 and t < 1 and v > 0 and v < 1:
                                res, d0, t0 = seg.point_sur_segment(p)
                                z = z0 - abs(d0) * slope
                                wall_t[widx].append((t, z, widx, idx, "slope s0"))

                            res, p, t, v = wseg.line.intersect_ext(s1)
                            if res and t > 0 and t < 1 and v > 0 and v < 1:
                                res, d0, t0 = seg.point_sur_segment(p)
                                z = z0 - abs(d0) * slope
                                wall_t[widx].append((t, z, widx, idx, "slope s1"))

                            # axis intersection
                            res, p, t, v = wseg.line.intersect_ext(seg)
                            if res and t > 0 and t < 1 and v > 0 and v < 1:
                                wall_t[widx].append((t, z0, widx, idx, "axis"))

                            # TODO:
                            # Roof bottom border intersection
                            # if ever this make sense

        old = context.active_object
        old_sel = wall.select
        wall.select = True
        context.scene.objects.active = wall

        wd.auto_update = False
        for seg in wall_t:
            self.sort_t(seg)
            print("seg: %s" % seg)
            for s in seg:
                t, z, widx, idx, type = s
                wd.parts[widx].n_splits = len(seg) + 1
                wd.parts[widx].z[0] = 0
                wd.parts[widx].t[0] = 0
                break

        for seg in wall_t:
            t0 = 0
            last_t = -1
            id = 1
            for s in seg:
                t, z, widx, idx, type = s
                if t == 0:
                    # add at end of last segment
                    if widx > 0:
                        lid = wd.parts[widx - 1].n_splits - 1
                        wd.parts[widx - 1].z[lid] = z
                        wd.parts[widx - 1].t[lid] = 1
                    else:
                        wd.parts[widx].z[0] = z
                        wd.parts[widx].t[0] = t
                        id = 1
                else:
                    if t == last_t:
                        wd.parts[widx].n_splits -= 1
                        continue
                    wd.parts[widx].z[id] = z
                    wd.parts[widx].t[id] = t - t0
                    t0 += t
                    id += 1
                    last_t = t

        if wd.closed:
            last = wd.parts[wd.n_parts].n_splits - 1
            wd.parts[wd.n_parts].z[last] = wd.parts[0].z[0]
            wd.parts[wd.n_parts].t[last] = 1.0

        wd.auto_update = True
        """
        for s in self.segs:
            s.as_curve(context)
        for s in wg.segs:
            s.as_curve(context)
        """
        wall.select = old_sel
        context.scene.objects.active = old

    def boundary(self, context, o, verts, edges):
        for b in o.children:
            d = archipack_roof_boundary.datablock(b)
            if d is not None:
                g = d.get_generator()
                g.change_coordsys(b.matrix_world, o.matrix_world)
                for pan in self.pans:
                    pan.slice(g)


"""
wall = C.object
o = C.object
d = o.data.archipack_roof[0]
d.make_wall_fit(C, wall)
"""

# bpy.app.debug = True

def update(self, context):
    self.update(context)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


def update_path(self, context):
    self.update_path(context)


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
                bpy.ops.archipack.reference_point(symbol_type='ROOF')
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


class ArchipackSegment():
    length = FloatProperty(
            name="length",
            min=0.01,
            max=1000.0,
            default=2.0,
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
    manipulators = CollectionProperty(type=archipack_manipulator)

    def find_in_selection(self, context):
        raise NotImplementedError

    def update(self, context, manipulable_refresh=False):
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context, manipulable_refresh)


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
    width_left = FloatProperty(
            name="L Width",
            min=0.01,
            default=3.0,
            update=update
            )
    width_right =FloatProperty(
            name="R Width",
            min=0.01,
            default=3.0,
            update=update
            )
    slope_left = FloatProperty(
            name="L slope",
            min=0.0,
            default=0.3,
            update=update
            )
    slope_right = FloatProperty(
            name="R slope",
            min=0.0,
            default=0.3,
            update=update
            )
    auto_left = EnumProperty(
            description="Left mode",
            name="Left",
            items=(
                ('AUTO', 'Auto', '', 0),
                ('WIDTH', 'Width', '', 1),
                ('SLOPE', 'Slope', '', 2),
                ('ALL', 'All', '', 3),
                ),
            default="AUTO",
            update=update_manipulators
            )
    auto_right = EnumProperty(
            description="Right mode",
            name="Right",
            items=(
                ('AUTO', 'Auto', '', 0),
                ('WIDTH', 'Width', '', 1),
                ('SLOPE', 'Slope', '', 2),
                ('ALL', 'All', '', 3),
                ),
            default="AUTO",
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

        box.prop(self, "constraint_type", text=str(index + 1))
        if self.constraint_type == 'SLOPE':
            box.prop(self, "enforce_part", text="")

        box.prop(self, "length")
        box.prop(self, "a0")
        box.prop(self, 'bound_idx')

        if self.constraint_type == 'HORIZONTAL':
            row = box.row(align=True)
            row.prop(self, "auto_left", text="")
            row.prop(self, "auto_right", text="")
            if self.auto_left in {'ALL', 'WIDTH'}:
                box.prop(self, "width_left")
            if self.auto_left in {'ALL', 'SLOPE'}:
                box.prop(self, "slope_left")
            if self.auto_right in {'ALL', 'WIDTH'}:
                box.prop(self, "width_right")
            if self.auto_right in {'ALL', 'SLOPE'}:
                box.prop(self, "slope_right")


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
    draft = BoolProperty(
            options={'SKIP_SAVE'},
            name="Draft mode",
            default=False,
            update=update
            )
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )
    quick_edit = BoolProperty(
            options={'SKIP_SAVE'},
            name="Quick Edit",
            default=True
            )
    force_update = BoolProperty(
            options={'SKIP_SAVE'},
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
            options={'SKIP_SAVE'},
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
            options={'SKIP_SAVE'},
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
            options={'SKIP_SAVE'},
            name="Beam",
            description="Expand beam panel",
            default=False
            )
    beam_enable = BoolProperty(
            name="Primary",
            default=True,
            update=update
            )
    beam_width = FloatProperty(
            name="Width",
            description="Width",
            min=0.01,
            default=0.2,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_height = FloatProperty(
            name="Height",
            description="Height",
            min=0.01,
            default=0.35,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_offset = FloatProperty(
            name="Offset",
            description="Distance from roof border",
            default=0.02,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_alt = FloatProperty(
            name="Altitude",
            description="Altitude from roof",
            default=-0.15,
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
            default=0.2,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    beam_sec_alt = FloatProperty(
            name="Altitude",
            description="Distance from roof",
            default=-0.1,
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
            default=0.2,
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
            options={'SKIP_SAVE'},
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
            name="Altitude",
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
            options={'SKIP_SAVE'},
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
            options={'SKIP_SAVE'},
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


    def make_wall_fit(self, context, o, wall, inside=False):
        origin = Vector((0, 0, self.z))
        g = self.get_generator(origin)
        g.make_roof(context, [], [])
        g.make_wall_fit(context, o, wall, inside)

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
            if n_manips < 4:
                s = p.manipulators.add()
                s.type_key = 'SIZE'
                s.prop1_name = "width_left"
            if n_manips < 5:
                s = p.manipulators.add()
                s.type_key = 'SIZE'
                s.prop1_name = "width_right"
            if n_manips < 6:
                s = p.manipulators.add()
                s.type_key = 'SIZE'
                s.prop1_name = "slope_left"
            if n_manips < 7:
                s = p.manipulators.add()
                s.type_key = 'SIZE'
                s.prop1_name = "slope_right"

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
            pg.make_roof(context, [], [])

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

        g.make_roof(context, [], [])

        g.boundary(context, o, verts, edges)

        if self.draft:

            for pan in g.pans:
                pan.draw(context, self.z, verts, edges)

            self.make_surface(o, verts, edges)

        else:
            if self.rake_enable:
                g.rake(self, verts, faces, edges, matids, uvs)

            if self.facia_enable:
                g.facia(self, verts, faces, edges, matids, uvs)

            if self.beam_enable:
                g.beam_primary(self, verts, faces, edges, matids, uvs)

            g.hips(self, verts, faces, edges, matids, uvs)

            if self.gutter_enable:
                g.gutter(self, verts, faces, edges, matids, uvs)


            bmed.buildmesh(

                context, o, verts, faces, matids=matids, uvs=uvs,
                weld=False, clean=False, auto_smooth=True, temporary=False)
            
            bpy.ops.object.mode_set(mode='EDIT')
            g.lambris(context, o, self)

            if self.rafter_enable:
                bpy.ops.object.mode_set(mode='EDIT')
                g.rafter(context, o, self)

            if self.tile_enable:
                bpy.ops.object.mode_set(mode='EDIT')
                g.couverture(context, o, self)


        """

        self.relocate_childs(context, o, g)

        if update_childs:
            self.update_childs(context, o, g)

        # print("%s" % (faces))

        n_verts = len(verts)




        if force_update:
            self.force_update = True

        if self.draft or self.force_update:
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

        """

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

            if part.constraint_type == 'HORIZONTAL':
                # length / radius + angle
                self.manip_stack.append(part.manipulators[1].setup(context, o, part))

            # index
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))

            # size left
            if part.auto_left in {'WIDTH', 'ALL'}:
                self.manip_stack.append(part.manipulators[3].setup(context, o, part))
            # size right
            if part.auto_right in {'WIDTH', 'ALL'}:
                self.manip_stack.append(part.manipulators[4].setup(context, o, part))
            # slope left
            if part.auto_left in {'SLOPE', 'ALL'}:
                self.manip_stack.append(part.manipulators[5].setup(context, o, part))
            # slope right
            if part.auto_right in {'SLOPE', 'ALL'}:
                self.manip_stack.append(part.manipulators[6].setup(context, o, part))



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


class archipack_roof_boundary_segment(ArchipackSegment, PropertyGroup):

    offset = FloatProperty(
            name="offset",
            min=0,
            default=0,
            update=update
            )

    type = EnumProperty(
        name="Type",
        items=(
            ('SIDE', 'Side', 'Side with rake', 0),
            ('BOTTOM', 'Bottom', 'Bottom with gutter', 1),
            ('LINK', 'Side link', 'Side witout decoration', 2),
            ('AXIS', 'Top', 'Top part with hip and beam', 3)
            # ('LINK_VALLEY', 'Side valley', 'Side with valley', 3),
            # ('LINK_HIP', 'Side hip', 'Side with hip', 4)
            ),
        default='BOTTOM',
        update=update
        )

    def find_in_selection(self, context):
        selected = [o for o in context.selected_objects]
        for o in selected:
            d = archipack_roof_boundary.datablock(o)
            if d:
                for part in d.parts:
                    if part == self:
                        return d
        return None

    def draw(self, layout, context, index):
        box = layout.box()
        row = box.row()
        row.prop(self, "type", text=str(index + 1))
        row = box.row()
        row.prop(self, "length")
        row = box.row()
        row.prop(self, "offset")
        row = box.row()
        row.prop(self, "a0")

        
def update_operation(self, context):
    self.reverse(context, make_ccw=(self.operation == 'INTERSECTION'))
    
    
class archipack_roof_boundary(ArchipackLines, ArchipackObject, Manipulable, PropertyGroup):
    # boundary
    parts = CollectionProperty(type=archipack_roof_boundary_segment)
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )
    z = FloatProperty(
            name="dumb z",
            description="Dumb z for manipulator",
            default=0.1,
            options={'SKIP_SAVE'}
            )    
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

    angle_limit = FloatProperty(
            name="angle",
            min=0,
            max=2 * pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update_manipulators
            )

    operation = EnumProperty(
            items=(
                ('DIFFERENCE', 'Difference', 'Cut inside part', 0),
                ('INTERSECTION', 'Intersection', 'Keep inside part', 1)
                ),
            default='DIFFERENCE',
            update=update_operation
            )        
    
    closed = BoolProperty(
            default=True,
            name="Close",
            options={'SKIP_SAVE'}
            )

    def setup_manipulators(self):
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

    def get_generator(self, origin=Vector((0, 0, 0))):
        g = RoofBoundaryGenerator(self, origin)
        for i, part in enumerate(self.parts):
            g.add_part(part)
        g.set_offset()
        g.close(self.closed)
        return g

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
        
    def is_cw(self, pts):
        p0 = pts[0]
        d = 0
        for p in pts[1:]:
            d += (p.x * p0.y - p.y * p0.x)
            p0 = p
        return d > 0

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
            if spline.use_cyclic_u:
                p0 = points[-1]
                p1 = points[0]
                self.interpolate_bezier(pts, wM, p0, p1, resolution)
                pts.append(pts[0])
            else:
                pts.append(wM * points[-1].co)

        if self.is_cw(pts):
            pts = list(reversed(pts))

        self.auto_update = False
        self.from_points(pts, spline.use_cyclic_u)
        self.auto_update = True

    def from_points(self, pts, closed):

        self.n_parts = len(pts) - 1

        if closed:
            self.n_parts -= 1

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
            if i >= len(self.parts):
                print("Too many pts for parts")
                break
            p = self.parts[i]
            p.length = dp.to_2d().length
            p.dz = dp.z
            p.a0 = da
            a0 += da
            p0 = p1

        self.closed = closed

    def reverse(self, context, make_ccw=False):

        g = self.get_generator()
        
        pts = [seg.p0.to_3d() for seg in g.segs]
        
        if self.is_cw(pts) != make_ccw:
            return
            
        types = [p.type for p in self.parts]
        
        if self.closed:
            pts.append(pts[0])
        
        pts = list(reversed(pts))
        self.auto_update = False

        # location wont change for closed boundary
        if not self.closed:
            dp = pts[0] - pts[-1]
            # pre-translate as dp is in local coordsys
            o.matrix_world = o.matrix_world * Matrix([
                [1, 0, 0, dp.x],
                [0, 1, 0, dp.y],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
                ])
        
        self.from_points(pts, self.closed)
        
        for i, type in enumerate(reversed(types)):
            self.parts[i].type = type
        
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
        bm.edges.new((bm.verts[-1], bm.verts[0]))
        bm.edges.ensure_lookup_table()
        # bmesh.ops.contextual_create(bm, geom=bm.edges)
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
        g.locate_manipulators()
        
        # vertex index in order to build axis
        g.get_verts(verts, edges)

        if len(verts) > 2:
            self.make_surface(o, verts, edges)
        else:
            g.debug(verts)
            bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=uvs, weld=True, clean=False)

        # enable manipulators rebuild
        if manipulable_refresh:
            self.manipulable_refresh = True
        
        d = archipack_roof.datablock(o.parent)
        if d is not None:
            o.parent.select = True
            context.scene.objects.active = o.parent
            d.update(context)
            
        o.parent.select = False
        context.scene.objects.active = o
        
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

                # length
                self.manip_stack.append(part.manipulators[1].setup(context, o, part))
                # index
                self.manip_stack.append(part.manipulators[3].setup(context, o, self))
                # offset
                self.manip_stack.append(part.manipulators[4].setup(context, o, part))

            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))


class ARCHIPACK_PT_roof_boundary(Panel):
    bl_idname = "ARCHIPACK_PT_roof_boundary"
    bl_label = "Roof Boundary"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_roof_boundary.filter(context.active_object)

    def draw(self, context):
        prop = archipack_roof_boundary.datablock(context.active_object)
        if prop is None:
            return
        scene = context.scene
        layout = self.layout
        box = layout.box()
        box.operator('archipack.roof_boundary_manipulate')
        box.prop(prop, 'operation', text="")
        box = layout.box()
        box.prop_search(prop, "user_defined_path", scene, "objects", text="", icon='OUTLINER_OB_CURVE')
        box.prop(prop, 'user_defined_resolution')
        box.prop(prop, 'x_offset')
        box.prop(prop, "closed")
        box.prop(prop, 'angle_limit')
        prop.draw(layout, context)


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
        o = context.active_object
        prop = archipack_roof.datablock(o)
        if prop is None:
            return
        scene = context.scene
        layout = self.layout
        row = layout.row(align=True)
        row.operator('archipack.roof_manipulate', icon='HAND')
        layout.operator('archipack.roof_boundary').parent = o.name

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
        box.prop(prop, 'draft')
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

            box.separator()     # hip
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
            box.separator()
            box.prop(prop, 'valley_enable')
            box.prop(prop, 'valley_altitude')

        box = layout.box()
        row = box.row(align=True)

        if prop.beam_expand:
            row.prop(prop, 'beam_expand', icon="TRIA_DOWN", text="Beam", icon_only=True, emboss=False)
        else:
            row.prop(prop, 'beam_expand', icon="TRIA_RIGHT", text="Beam", icon_only=True, emboss=False)
        if prop.beam_expand:
            box.prop(prop, 'beam_enable')
            if prop.beam_enable:
                box.prop(prop, 'beam_width')
                box.prop(prop, 'beam_height')
                box.prop(prop, 'beam_offset')
                box.prop(prop, 'beam_alt')
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


class ARCHIPACK_OT_roof_boundary(ArchipackCreateTool, Operator):
    bl_idname = "archipack.roof_boundary"
    bl_label = "Roof Cutter"
    bl_description = "Roof Cutter"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    
    parent = StringProperty("")
    
    def create(self, context):
        m = bpy.data.meshes.new("Roof Boundary")
        o = bpy.data.objects.new("Roof Boundary", m)
        d = m.archipack_roof_boundary.add()
        parent = context.scene.objects.get(self.parent)
        if parent is not None:
            o.parent = parent
            bbox = parent.bound_box
            angle_90 = pi / 2
            x0, y0, z = bbox[0]
            x1, y1, z = bbox[6]
            x = 0.2 * (x1 - x0)
            y = 0.2 * (y1 - y0)
            o.matrix_world = parent.matrix_world * Matrix([
                [1, 0, 0, -3 * x],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
                ])
            p = d.parts.add()
            p.a0 = - angle_90
            p.length = y
            p = d.parts.add()
            p.a0 = angle_90
            p.length = x
            p = d.parts.add()
            p.a0 = angle_90
            p.length = y
            d.n_parts = 3
            d.close = True
            pd = archipack_roof.datablock(parent)
            pd.boundary = o.name
        else:
            o.location = context.scene.cursor_location
        # make manipulators selectable
        d.manipulable_selectable = True
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        self.add_material(o)
        self.load_preset(d)
        update_operation(d, context)
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


class ARCHIPACK_OT_roof_boundary_manipulate(Operator):
    bl_idname = "archipack.roof_boundary_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return archipack_roof_boundary.filter(context.active_object)

    def invoke(self, context, event):
        d = archipack_roof_boundary.datablock(context.active_object)
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
    # bpy.utils.register_class(archipack_roof_material)
    bpy.utils.register_class(archipack_roof_boundary_segment)
    bpy.utils.register_class(archipack_roof_boundary)
    bpy.utils.register_class(ARCHIPACK_PT_roof_boundary)
    bpy.utils.register_class(ARCHIPACK_OT_roof_boundary)
    bpy.utils.register_class(ARCHIPACK_OT_roof_boundary_manipulate)
    Mesh.archipack_roof_boundary = CollectionProperty(type=archipack_roof_boundary)
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
    # bpy.utils.unregister_class(archipack_roof_material)
    bpy.utils.unregister_class(archipack_roof_boundary_segment)
    bpy.utils.unregister_class(archipack_roof_boundary)
    bpy.utils.unregister_class(ARCHIPACK_PT_roof_boundary)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_boundary)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_boundary_manipulate)
    del Mesh.archipack_roof_boundary
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
