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
    FloatProperty, BoolProperty, IntProperty,
    StringProperty, EnumProperty,
    CollectionProperty, PointerProperty
    )
from .bmesh_utils import BmeshEdit as bmed
import bmesh
from .materialutils import MaterialUtils
from .panel import Panel as Lofter
from mathutils import Vector, Matrix
from mathutils.geometry import interpolate_bezier
from math import sin, cos, pi, acos, atan2
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
        """
            setup corners and roofs dz
            compute index of roofs wich belong to each group of roofs between corners
            compute t of each roof
        """
            
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
        
        """
           6 _________________ 5
            |         |       |
          __|_______b |_b'____|__  axe principal  
         0  |     4 a | a'    | 3
            |_________|_______|
         1                      2
         
         
        
        """
        
        
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
               

class RoofAxisNode():
    def __init__(self, a0, idx, reversed):
        self.a0 = a0
        self.idx = idx
        self.reversed = reversed

               
class RoofAxisGenerator():

    def __init__(self, parts, origin=Vector((0, 0))):
        self.parts = parts
        self.segs = []
        self.length = 0
        self.origin = origin
        self.user_defined_post = None
        self.user_defined_uvs = None
        self.user_defined_mat = None
        self.z = 0
        self.slope = 1
        self.width = 2.5
        
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
        s.v0_idx = part.bound_idx
        s.constraint_type = part.constraint_type
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
                v0 = self.segs[f.v0_idx].straight(-1, 1).v.to_3d()
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
 
    def intersect_chevron(self, line, side1, side2, slope, height, segs, verts):
        res, p, u, v = line.intersect_ext(segs[side1])
        if res :
            x, y = p
            z = self.z - slope * u
        else:
            res, p, u, v = line.intersect_ext(segs[side2])
            if res:
                x, y = p
                z = self.z - slope * u
            else:
                x, y = line.p1
                z = self.z - slope
        
        verts.append((x, y, z))
        verts.append((x, y, z - height))
 
    def make_roof(self, verts, edges, faces):
        
        x, y = self.segs[0].p0
        z = self.z
        verts.append((x, y, z))
        
        for idx, s in enumerate(self.segs):
            
            s.v1_idx = idx + 1
            x, y = s.p1
            verts.append((x, y, z))
            s.z0 = z
            
        # node are connected segments
        # node    
        # (segment idx) 
        # (angle from root part > 0 right) 
        # (reversed) connected by p1
        
        nodes = [[] for s in range(len(self.segs) + 1)]
        new_idx = len(self.segs)
        
        for i, s in enumerate(self.segs):
            nodes[s.v0_idx].append(RoofAxisNode(s.angle_0, i, False))
            nodes[s.v1_idx].append(RoofAxisNode(-pi, i, True))
        
        segs = []
        for idx, node in enumerate(nodes):
            nb_segs = len(node)
            # if node is empty:
            if nb_segs < 1:
                print("empty node")
                break
            
            # if node is alone:
            # build 2 slopes on each sides
            root = self.segs[node[0].idx]
            
            # check if there is a user defined slope
            # disable auto-slope when found
            slope_left = None
            slope_right = None
            for s in node:
                seg = self.segs[s.idx]
                if not s.reversed and seg.constraint_type == 'SLOPE':
                    nb_segs -= 1
                    if seg.angle_0 > 0:
                        slope_left = s
                    else:
                        slope_right = s
            
            # ends either start or end of roof parts
            # where segment is alone or with user defined slope
            if nb_segs < 2:
                
                if node[0].reversed:
                    t = 1
                else:
                    segs.append(root)
                    t = 0
                
                if root.constraint_type == 'HORIZONTAL':
                    
                    # left - auto slope 
                    if slope_left is None:
                        length = self.width
                        new_idx += 1
                        new = root.straight(length, t).rotate(pi / 2)
                        new.constraint_type = 'SLOPE'
                        new.angle_0 = pi / 2
                        new.v0_idx = idx
                        new.v1_idx = new_idx
                        segs.append(new)
                        x, y = new.p1
                        z = root.z0 - length * self.slope
                        verts.append((x, y, z))
                    
                    # user defined left slope
                    else: 
                        cur = self.segs[slope_left.idx]
                        # reproject user defined slope 
                        da = cur.delta_angle(root)
                        # angle always ccw
                        if da < 0:
                            da = 2 * pi + da
                        length = min(3 * self.width, self.width / sin(da))
                        # print("length:%s da:%s" % (length, da))
                        cur = cur.copy()
                        cur.v = cur.v.normalized() * length
                        segs.append(cur)
                        x, y = cur.p1
                        z = root.z0 - self.width * self.slope
                        verts[cur.v1_idx] = (x, y, z)
                        
                    # right - auto slope
                    if slope_right is None:    
                        length = self.width
                        new_idx += 1
                        new = root.straight(length, t).rotate(-pi / 2)
                        new.constraint_type = 'SLOPE'
                        new.angle_0 = -pi / 2
                        new.v0_idx = idx
                        new.v1_idx = new_idx
                        segs.append(new)
                        x, y = new.p1
                        z = root.z0 - length * self.slope
                        verts.append((x, y, z))
                    
                    # user defined right slope
                    else:
                        cur = self.segs[slope_right.idx]
                        # reproject user defined slope 
                        da = cur.delta_angle(root)
                        # angle always ccw
                        if da < 0:
                            da = 2 * pi + da
                        length = min(3 * self.width, self.width / sin(da))
                        # print("length:%s da:%s" % (length, da))
                        cur = cur.copy()
                        cur.v = cur.v.normalized() * -length
                        segs.append(cur)
                        x, y = cur.p1
                        z = root.z0 - self.width * self.slope
                        verts[cur.v1_idx] = (x, y, z)
            
            # middle parts of roof
            else:
                self.sort_seg(node)
                
                # add slope between horizontal parts
                last = self.segs[node[-1].idx]
                
                if node[-1].reversed:
                    last = last.copy()
                    last.p += last.v
                    last.v = -last.v
                
                # new_node.append(last)
                for n in node:
                    cur = self.segs[n.idx]
                    if n.reversed:
                        cur = cur.copy()
                        cur.p += cur.v
                        cur.v = -cur.v
                    
                    if cur.constraint_type == 'HORIZONTAL':
                        if last.constraint_type == cur.constraint_type:
                            # add slope in the middle
                            
                            da = cur.delta_angle(last)
                            # angle always ccw
                            if da < 0:
                                da = 2 * pi + da
                                
                            length = min(3 * self.width, self.width / sin(0.5 * da))
                            new_idx += 1
                            new = last.straight(length, 0).rotate(0.5 * da)
                            new.constraint_type = 'SLOPE'
                            new.v0_idx = idx
                            new.v1_idx = new_idx
                            new.angle_0 = new.delta_angle(root)
                            segs.append(new)
                            x, y = new.p1
                            z = root.z0 - self.width * self.slope
                            verts.append((x, y, z))
                    
                    # reproject user defined slope 
                    else:
                        da = cur.delta_angle(root)
                        # angle always ccw
                        if da < 0:
                            da = 2 * pi + da
                        length = min(3 * self.width, self.width / sin(da))
                        # print("length:%s da:%s" % (length, da))
                        cur = cur.copy()
                        cur.v = cur.v.normalized() * length
                        if cur.angle_0 < 0:
                            cur.v = -cur.v
                        x, y = cur.p1
                        z = root.z0 - self.width * self.slope
                        verts[cur.v1_idx] = (x, y, z)
                        
                    if not n.reversed:
                        segs.append(cur)
                    
                    last = cur
        
        # axis and slope edges
        for s in segs:
            edges.append([s.v0_idx, s.v1_idx])
        
        # boundary edges / faces
        nodes = [[] for s in range(new_idx + 1)]
        for i, s in enumerate(segs):
            nodes[s.v0_idx].append(RoofAxisNode(s.angle_0, i, False))
            nodes[s.v1_idx].append(RoofAxisNode(-pi, i, True))
        
        for idx, node in enumerate(nodes):
            self.sort_seg(node)
        
        # Chevrons
        spacing = 0.7
        start = 0.1
        depth = 0.1
        chevron_height = 0.1
        
        faitiere_width = 0.2
        faitiere_height = 0.3
        
        larmier_width = 0.02
        larmier_height = 0.25
        larmier_altitude = 0.1
        virevents_depth = 0.02
        virevents_height = 0.15
        virevents_altitude = 0.1
        
        for idx, node in enumerate(nodes):
            # segment always between 2 nodes
            # create edge between rightmost of first node to leftmost of next node
            # same for other side
            
            # find next node in segment
            for i, s in enumerate(node):
                
                seg = segs[s.idx]
                
                if s.reversed:
                    continue
                    
                if len(nodes[seg.v1_idx]) > 1:
                    next = nodes[seg.v1_idx]
                    # segments sorted by angle from axis 
                    # center are ids of axis segment on each node
                    # so center + 1 and center - 1
                    # are leftmost and rightmost slope segments
                    center_0 = i
                    center_1 = -1
                    for i, s in enumerate(next):
                        if s.reversed:
                            center_1 = i
                            break
                    # found 2nd node
                    if center_1 > -1:
                    
                        i0 = seg.v0_idx
                        i1 = seg.v1_idx
                        
                        ############
                        # right
                        ############
                        
                        ir0 = next[center_1 - 1].idx
                        if center_0 + 1 < len(node):
                            ir1 = node[center_0 + 1].idx
                        else:
                            ir1 = node[0].idx
                        
                        vi2 = segs[ir0].v1_idx
                        vi3 = segs[ir1].v1_idx
                        
                        edges.append([vi2, vi3])
                        faces.append([i0, i1, vi2, vi3])
 
                        
                        ############
                        # left
                        ############
                        
                        il0 = node[center_0 - 1].idx
                        if center_1 + 1 < len(next):
                            il1 = next[center_1 + 1].idx
                        else:
                            il1 = next[0].idx
                        
                        vi2 = segs[il0].v1_idx
                        vi3 = segs[il1].v1_idx
                        
                        edges.append([vi2, vi3])
                        faces.append([i1, i0, vi2, vi3])
                        
                        # right part is larger than axis: compute t param in axis
                        res, d, u = seg.point_sur_segment(segs[ir1].p1)
                        res, dr, v = seg.point_sur_segment(segs[ir0].p1)
                        rslope = abs(dr) * self.slope
                        trmin = min(0, u)
                        trmax = max(1, v)
                        
                        # left part is larger than axis: compute t param in axis
                        res, d, u = seg.point_sur_segment(segs[il0].p1)
                        res, dl, v = seg.point_sur_segment(segs[il1].p1)
                        lslope = abs(dl) * self.slope
                        tlmin = min(0, u)
                        tlmax = max(1, v)
                        
                        #####################
                        # Couverture
                        #####################
                        """
                            store verts, faces, and 4 matrix to cut borders
                                
                        
                        # compute matrix for face
                        vx = seg.v.normalized().to_3d()
                        vy = Vector((vx.y, -vx.x, self.slope))
                        vz = vx.cross(vy)
                        m0 = [seg.p0.to_3d(), segs[ir0].cross_z.to_3d()]
                        m1 = [seg.p1.to_3d(), segs[ir1].cross_z.to_3d()]
                        
                        for j in range(n_x):
                            x, y = seg.lerp(trmin + j * dt)
                            for k in range(n_y):
                                y += k * dy   
                                tM = Matrix([
                                    [vx.x, vy.x, vz.x, x],
                                    [vx.y, vy.y, vz.y, y],
                                    [vx.z, vy.z, vz.z, self.z],
                                    [0, 0, 0, 1]
                                ])
                                tuile_pts.extend([tM * p for p in t_pts]) 
                                tuile_fcs.extend([(i + v for v in f) for f in t_faces])
                        
                        tuiles.append([tuile_pts, tuile_fcs, m0, m1, m2, m3])   
                        """
                        
                        #####################
                        # Vire-vents
                        #####################
                        if len(node) == 3:
                            f = len(verts)
                            s0 = segs[il0]
                            s1 = s0.offset(virevents_depth)
                            res, d, t = seg.point_sur_segment(s0.p1)
                            slope = abs(d) * self.slope
                        
                            x0, y0 = s0.p0
                            x1, y1 = s1.p0
                            x2, y2 = s1.p1
                            x3, y3 = s1.p1
                            z0 = self.z + virevents_altitude
                            z1 = z0 - slope
                            verts.extend([
                                (x0, y0, z0),
                                (x1, y1, z0),
                                (x2, y2, z1),
                                (x3, y3, z1),
                            ])
                            z0 -= virevents_height
                            z1 -= virevents_height
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
                                (f + 3, f + 7, f + 4, f ),
                                # bottom
                                (f + 4, f + 7, f + 6, f + 5)
                            ])
                            edges.append([f, f + 3])
                            edges.append([f + 1, f + 2])
                            edges.append([f + 4, f + 7])
                            edges.append([f + 5, f + 6])
                                
                                
                            f = len(verts)
                            s1 = segs[ir1]
                            s0 = s1.offset(-virevents_depth)
                            res, d, t = seg.point_sur_segment(s0.p1)
                            slope = abs(d) * self.slope
                        
                            x0, y0 = s0.p0
                            x1, y1 = s1.p0
                            x2, y2 = s1.p1
                            x3, y3 = s1.p1
                            
                            z0 = self.z + virevents_altitude
                            z1 = z0 - slope
                            verts.extend([
                                (x0, y0, z0),
                                (x1, y1, z0),
                                (x2, y2, z1),
                                (x3, y3, z1),
                            ])
                            z0 -= virevents_height
                            z1 -= virevents_height
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
                                (f + 3, f + 7, f + 4, f ),
                                # bottom
                                (f + 4, f + 7, f + 6, f + 5)
                            ])
                            edges.append([f, f + 3])
                            edges.append([f + 1, f + 2])
                            edges.append([f + 4, f + 7])
                            edges.append([f + 5, f + 6])
                            
                        if len(next) == 3:
                            f = len(verts)
                            s0 = segs[ir0]
                            s1 = s0.offset(virevents_depth)
                            res, d, t = seg.point_sur_segment(s0.p1)
                            slope = abs(d) * self.slope
                        
                            x0, y0 = s0.p0
                            x1, y1 = s1.p0
                            x2, y2 = s1.p1
                            x3, y3 = s1.p1
                            
                            z0 = self.z + virevents_altitude
                            z1 = z0 - slope
                            verts.extend([
                                (x0, y0, z0),
                                (x1, y1, z0),
                                (x2, y2, z1),
                                (x3, y3, z1),
                            ])
                            z0 -= virevents_height
                            z1 -= virevents_height
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
                                (f + 3, f + 7, f + 4, f ),
                                # bottom
                                (f + 4, f + 7, f + 6, f + 5)
                            ])
                            edges.append([f, f + 3])
                            edges.append([f + 1, f + 2])
                            edges.append([f + 4, f + 7])
                            edges.append([f + 5, f + 6])
                            
                            
                            f = len(verts)
                            s1 = segs[il1]
                            s0 = s1.offset(-virevents_depth)
                            res, d, t = seg.point_sur_segment(s0.p1)
                            slope = abs(d) * self.slope
                        
                            x0, y0 = s0.p0
                            x1, y1 = s1.p0
                            x2, y2 = s1.p1
                            x3, y3 = s1.p1
                            
                            z0 = self.z + virevents_altitude
                            z1 = z0 - slope
                            verts.extend([
                                (x0, y0, z0),
                                (x1, y1, z0),
                                (x2, y2, z1),
                                (x3, y3, z1),
                            ])
                            z0 -= virevents_height
                            z1 -= virevents_height
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
                                (f + 3, f + 7, f + 4, f ),
                                # bottom
                                (f + 4, f + 7, f + 6, f + 5)
                            ])
                            edges.append([f, f + 3])
                            edges.append([f + 1, f + 2])
                            edges.append([f + 4, f + 7])
                            edges.append([f + 5, f + 6])
                            
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
                        res, d, t = seg.point_sur_segment(s0.p1)
                        slope = abs(d) * self.slope
                        a0 = abs(seg.delta_angle(s0))
                        t0 = 1 + min(3 * larmier_width, larmier_width / sin(a0)) / s0.length
                        a1 = abs(seg.delta_angle(s1))
                        t1 = 1 + min(3 * larmier_width, larmier_width / sin(a1)) / s1.length
                        x0, y0 = s0.p1
                        x1, y1 = s0.lerp(t0)
                        x2, y2 = s1.lerp(t1)
                        x3, y3 = s1.p1
                        z = self.z + larmier_altitude - slope
                        verts.extend([
                            (x0, y0, z),
                            (x1, y1, z),
                            (x2, y2, z),
                            (x3, y3, z),
                        ])
                        z -= larmier_height
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
                            (f + 3, f + 7, f + 4, f ),
                            # bottom
                            (f + 4, f + 7, f + 6, f + 5)
                        ])
                        edges.append([f, f + 3])
                        edges.append([f + 1, f + 2])
                        edges.append([f + 4, f + 7])
                        edges.append([f + 5, f + 6])
                        
                        f = len(verts)
                        s0 = segs[il0]
                        s1 = segs[il1]
                        res, d, t = seg.point_sur_segment(s0.p1)
                        slope = abs(d) * self.slope
                        a0 = abs(seg.delta_angle(s0))
                        t0 = 1 + min(3 * larmier_width, larmier_width / sin(a0)) / s0.length
                        a1 = abs(seg.delta_angle(s1))
                        t1 = 1 + min(3 * larmier_width, larmier_width / sin(a1)) / s1.length
                        x0, y0 = s0.p1
                        x1, y1 = s0.lerp(t0)
                        x2, y2 = s1.lerp(t1)
                        x3, y3 = s1.p1
                        z = self.z + larmier_altitude - slope
                        verts.extend([
                            (x0, y0, z),
                            (x1, y1, z),
                            (x2, y2, z),
                            (x3, y3, z),
                        ])
                        z -= larmier_height
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
                            (f + 3, f + 7, f + 4, f ),
                            # bottom
                            (f + 4, f + 7, f + 6, f + 5)
                        ])
                        edges.append([f, f + 3])
                        edges.append([f + 1, f + 2])
                        edges.append([f + 4, f + 7])
                        edges.append([f + 5, f + 6])
                        
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
                        
                        left = seg.offset(0.5 * faitiere_width)
                        right = seg.offset(-0.5 * faitiere_width)
                        res, p1, t = left.intersect(segs[il0])
                        res, p2, t = left.intersect(segs[il1])
                        res, p4, t = right.intersect(segs[ir0])
                        res, p5, t = right.intersect(segs[ir1])
                        x0, y0 = seg.p0
                        x1, y1 = p1
                        x2, y2 = p2
                        x3, y3 = seg.p1
                        x4, y4 = p4
                        x5, y5 = p5
                        z = self.z - chevron_height
                        verts.extend([
                            (x0, y0, z),
                            (x1, y1, z),
                            (x2, y2, z),
                            (x3, y3, z),
                            (x4, y4, z),
                            (x5, y5, z)
                        ])
                        z -= faitiere_height
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
                        
                        ######################
                        # right part chevrons
                        ######################
                        
                        f = len(verts)
                        
                        
                        t0 = trmin + (start - 0.5 * depth) / seg.length
                        t1 = trmin + (start + 0.5 * depth) / seg.length
                        
                        tx = start / seg.length
                        dt = spacing / seg.length
                        
                        n_items = max(1, round((trmax - trmin) / dt, 0))
                        
                        dt = ((trmax - trmin) - 2 * tx) / n_items
                        
                        for j in range(int(n_items) + 1):
                            # 
                            n0 = seg.sized_normal(t1 + j * dt, - self.width)
                            n1 = seg.sized_normal(t0 + j * dt, - self.width)
                            slope = n0.length * self.slope
                            
                            # verts start from axis
                            right_before = t1 + j * dt > 0 and t1 + j * dt < 1
                            left_before = t0 + j * dt > 0 and t0 + j * dt < 1
                            
                            if right_before:
                                z = self.z
                                x, y = n0.p0
                                verts.append((x, y, z))
                                verts.append((x, y, z - chevron_height))
                            
                            self.intersect_chevron(n0, ir0, ir1, slope, chevron_height, segs, verts)
                            
                            if not right_before:
                                z = self.z - slope
                                x, y = n0.p1
                                verts.append((x, y, z))
                                verts.append((x, y, z - chevron_height))
                            
                            
                            if left_before:
                                z = self.z
                                x, y = n1.p0
                                verts.append((x, y, z))
                                verts.append((x, y, z - chevron_height))
                            
                            self.intersect_chevron(n1, ir0, ir1, slope, chevron_height, segs, verts)
      
                            if not left_before:    
                                z = self.z - slope
                                x, y = n1.p1
                                verts.append((x, y, z))
                                verts.append((x, y, z - chevron_height))
                            
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
                            f += 8
                        
                        
                        
                        
                        
                        ######################
                        # left part chevrons
                        ######################
                        
                        f = len(verts)
                        
                        
                        
                        t0 = tlmin + (start - 0.5 * depth) / seg.length
                        t1 = tlmin + (start + 0.5 * depth) / seg.length
                        
                        tx = start / seg.length
                        dt = spacing / seg.length
                        
                        n_items = max(1, round((tlmax - tlmin) / dt, 0))
                        
                        dt = ((tlmax - tlmin) - 2 * tx) / n_items
                        
                        for j in range(int(n_items) + 1):
                            n0 = seg.sized_normal(t0 + j * dt, self.width)
                            n1 = seg.sized_normal(t1 + j * dt, self.width)
                            slope = n0.length * self.slope
                            
                            # verts start from axis
                            right_before = t0 + j * dt > 0 and t0 + j * dt < 1
                            left_before = t1 + j * dt > 0 and t1 + j * dt < 1
                            
                            if right_before:
                                z = self.z
                                x, y = n0.p0
                                verts.append((x, y, z))
                                verts.append((x, y, z - chevron_height))
                            
                            self.intersect_chevron(n0, il0, il1, slope, chevron_height, segs, verts)
                            
                            if not right_before:
                                z = self.z - slope
                                x, y = n0.p1
                                verts.append((x, y, z))
                                verts.append((x, y, z - chevron_height))
                            
                            if left_before:
                                z = self.z
                                x, y = n1.p0
                                verts.append((x, y, z))
                                verts.append((x, y, z - chevron_height))
                            
                            self.intersect_chevron(n1, il0, il1, slope, chevron_height, segs, verts)
                                
                            if not left_before:    
                                z = self.z - slope
                                x, y = n1.p1
                                verts.append((x, y, z))
                                verts.append((x, y, z - chevron_height))
                                
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
                            f += 8
    
# bpy.app.debug = True
                
def update(self, context):
    self.update(context)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


def update_path(self, context):
    self.update_path(context)


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
                w = w0.straight_roof(self.a0, self.length)
            else:
                w = w0.curved_roof(self.a0, self.da, self.radius)
        else:
            g = RoofGenerator(None)
            g.add_part(self)
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
            default='S_SEG',
            update=update_type
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
            box.prop(self, 'closed')
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
        row_change = False
        for i in range(len(self.parts), self.n_parts + 1, -1):
            row_change = True
            self.parts.remove(i - 1)

        # add rows
        for i in range(len(self.parts), self.n_parts + 1):
            row_change = True
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

    def get_generator(self, origin=Vector((0, 0))):
        g = RoofGenerator(self.parts, origin)
        for part in self.parts:
            g.add_part(part)
        g.set_offset()
        g.close(self.closed)
        g.locate_manipulators()
        return g


class archipack_roof_axis_segment(ArchipackSegment, PropertyGroup):
    
    bound_idx = IntProperty(
        default=0, 
        min=0,
        update=update_manipulators
        )
    constraint_type = EnumProperty(
        items=(
            ('HORIZONTAL', 'Horizontal', '', 0),
            ('SLOPE', 'Slope', '', 1)
            ),
        default='HORIZONTAL'    
        )
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            d = archipack_roof_axis.datablock(o)
            if d:
                for part in d.parts:
                    if part == self:
                        return d
        return None
    def draw(self, layout, context, index):
        box = layout.box()
        row = box.row()
        row.prop(self, "constraint_type", text=str(index + 1))
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
        

class archipack_roof_axis(ArchipackLines, ArchipackObject, Manipulable, PropertyGroup):
    auto_update = BoolProperty(default=True)
    parts = CollectionProperty(type=archipack_roof_axis_segment)
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
    width = FloatProperty(
            name="width",
            default=2.5, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    closed = BoolProperty(
            default=False,
            name="Close",
            update=update_manipulators
            )
    do_faces = BoolProperty(
            name="Make faces",
            default=False,
            update=update
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
    angle_limit = FloatProperty(
            name="angle",
            min=0,
            max=2 * pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update_manipulators
            )
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )
    
    def update_parts(self):
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
            bound_idx = len(self.parts)
            self.parts.add()
            self.parts[-1].bound_idx = bound_idx
            
        self.setup_manipulators()
    
    def setup_manipulators(self):
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
    
    def get_generator(self, origin=Vector((0, 0))):
        g = RoofAxisGenerator(self.parts, origin)
        g.slope = self.slope
        g.z = self.z
        g.width = self.width
        for part in self.parts:
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
        matids = None
        uvs = None

        g = self.get_generator()
        g.make_roof(verts, edges, faces)
        # print("%s" % (faces))
        
       
        if self.do_faces:
            bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=uvs, weld=False, clean=False, auto_smooth=False)
        else:
            self.make_surface(o, verts, edges)
        
        # enable manipulators rebuild
        if manipulable_refresh:
            self.manipulable_refresh = True

        # restore context
        self.restore_context(context)

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
            self.manip_stack.append(part.manipulators[3].setup(context, o, self))
        
            # snap point
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

      
class archipack_roof_boundary(ArchipackSegment, PropertyGroup):

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


class archipack_roof(ArchipackLines, ArchipackObject, Manipulable, PropertyGroup):
    # boundary
    parts = CollectionProperty(type=archipack_roof_boundary)
    # axis = PointerProperty(type=archipack_roof_axis)

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
            axis = archipack_roof_axis.datablock(a)
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
        """
            NOTE:
            this one assume context.active_object is the instance this
            data belongs to, failing to do so will result in wrong
            manipulators set on active object
        """
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

class ARCHIPACK_PT_roof_axis(Panel):
    bl_idname = "ARCHIPACK_PT_roof_axis"
    bl_label = "Roof axis"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_roof_axis.filter(context.active_object)

    def draw(self, context):
        prop = archipack_roof_axis.datablock(context.active_object)
        if prop is None:
            return
        scene = context.scene
        layout = self.layout
        row = layout.row(align=True)
        row.operator('archipack.roof_axis_manipulate')
        # box.label(text="Styles")
        box = layout.box()
        box.prop(prop, 'z')
        box.prop(prop, 'slope')
        box.prop(prop, 'width')
        box.prop(prop, 'do_faces')
        """
        box = layout.box()
        row.prop_search(prop, "user_defined_path", scene, "objects", text="", icon='OUTLINER_OB_CURVE')
        box.prop(prop, 'user_defined_resolution')
        box.prop(prop, 'angle_limit')
        """
        prop.draw(layout, context)

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


class ARCHIPACK_OT_roof_axis(ArchipackCreateTool, Operator):
    bl_idname = "archipack.roof_axis"
    bl_label = "Roof axis"
    bl_description = "Roof axis"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    def create(self, context):
        m = bpy.data.meshes.new("Roof axis")
        o = bpy.data.objects.new("Roof axis", m)
        d = m.archipack_roof_axis.add()
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
            parent = context.active_object
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = context.scene.cursor_location
            o.select = True
            context.scene.objects.active = o
            o.parent = parent
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
        
 
class ARCHIPACK_OT_roof_axis_manipulate(Operator):
    bl_idname = "archipack.roof_axis_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return archipack_roof_axis.filter(context.active_object)

    def invoke(self, context, event):
        d = archipack_roof_axis.datablock(context.active_object)
        d.manipulable_invoke(context)
        return {'FINISHED'}
        
 
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
    bpy.utils.register_class(archipack_roof_boundary)
    bpy.utils.register_class(archipack_roof_axis_segment)
    bpy.utils.register_class(archipack_roof_axis)
    bpy.utils.register_class(archipack_roof)
    Mesh.archipack_roof = CollectionProperty(type=archipack_roof)
    bpy.utils.register_class(ARCHIPACK_OT_roof_preset_menu)
    bpy.utils.register_class(ARCHIPACK_PT_roof)
    bpy.utils.register_class(ARCHIPACK_OT_roof)
    bpy.utils.register_class(ARCHIPACK_OT_roof_preset)
    bpy.utils.register_class(ARCHIPACK_OT_roof_manipulate)
    bpy.utils.register_class(ARCHIPACK_OT_roof_from_curve)
    Mesh.archipack_roof_axis = CollectionProperty(type=archipack_roof_axis)
    bpy.utils.register_class(ARCHIPACK_PT_roof_axis)
    bpy.utils.register_class(ARCHIPACK_OT_roof_axis)
    bpy.utils.register_class(ARCHIPACK_OT_roof_axis_manipulate)

def unregister():
    bpy.utils.unregister_class(archipack_roof_material)
    bpy.utils.unregister_class(archipack_roof_boundary)
    bpy.utils.unregister_class(archipack_roof_axis_segment)
    bpy.utils.unregister_class(archipack_roof_axis)
    bpy.utils.unregister_class(archipack_roof)
    del Mesh.archipack_roof
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_PT_roof)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_manipulate)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_PT_roof_axis)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_axis)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_axis_manipulate)
    del Mesh.archipack_roof_axis
    