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
    FloatProperty, BoolProperty, IntProperty, CollectionProperty,
    StringProperty, EnumProperty
    )
from .bmesh_utils import BmeshEdit as bmed
from .panel import Panel as Lofter
from mathutils import Vector, Matrix
from math import sin, cos, pi, acos
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_2d import Line, Arc
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_object import ArchipackCreateTool, ArchipackObject, ArchipackObjectsManager
from .archipack_dimension import DimensionProvider
from .archipack_curveman import ArchipackProfile, ArchipackUserDefinedPath
from .archipack_segments import ArchipackSegment, Segment, Generator


class Fence(Segment):

    def __init__(self):
        # total distance from start
        Segment.__init__(self)
        self.dist = 0
        self.t_start = 0
        self.t_end = 0
        self.dz = 0
        self.z0 = 0

    @property
    def t_diff(self):
        return self.t_end - self.t_start


class StraightSegment(Fence, Line):
    def __str__(self):
        return "t_start:{} t_end:{} dist:{}".format(self.t_start, self.t_end, self.dist)

    def __init__(self, p, v):
        Fence.__init__(self)
        Line.__init__(self, p, v)


class CurvedSegment(Fence, Arc):
    def __str__(self):
        return "t_start:{} t_end:{} dist:{}".format(self.t_start, self.t_end, self.dist)

    def __init__(self, c, radius, a0, da):
        Fence.__init__(self)
        Arc.__init__(self, c, radius, a0, da)


class FenceSegment():
    def __str__(self):
        return "t_start:{} t_end:{} n_step:{}  t_step:{} i_start:{} i_end:{}".format(
            self.t_start, self.t_end, self.n_step, self.t_step, self.i_start, self.i_end)

    def __init__(self, t_start, t_end, n_step, t_step, i_start, i_end):
        self.t_start = t_start
        self.t_end = t_end
        self.n_step = n_step
        self.t_step = t_step
        self.i_start = i_start
        self.i_end = i_end


class FenceGenerator(Generator):

    def __init__(self, d, o=None):
        Generator.__init__(self, d, o)
        self.length = 0
        self.user_defined_post = None
        self.user_defined_uvs = None
        self.user_defined_mat = None

    @property
    def n_parts(self):
        return len(self.parts)

    def add_part(self, part):

        last = None

        if len(self.segs) > 0:
            last = self.segs[-1]

        if 'S_' in part.type:
            if last is None:
                p = self.location
                v = (self.rot * Vector((part.length, 0, 0))).to_2d()
                s = StraightSegment(p, v).rotate(part.a0)
            else:
                seg = last.straight(part.length).rotate(part.a0)
                s = StraightSegment(seg.p, seg.v)
        else:
            if last is None:
                c = self.location - (self.rot * (part.radius * Vector((cos(part.a0), sin(part.a0), 0)))).to_2d()
                s = CurvedSegment(c, part.radius, part.a0, part.da)
            else:
                n = last.normal(1).rotate(part.a0).scale(part.radius)
                if part.da < 0:
                    n.v = -n.v
                a = n.angle
                c = n.p - n.v
                s = CurvedSegment(c, part.radius, a, part.da)

        self.segs.append(s)

        self.last_type = part.type

    def param_t(self, angle_limit, post_spacing):
        """
            setup corners and fences dz
            compute index of fences wich belong to each group of fences between corners
            compute t of each fence
        """
        # segments are group of parts separated by limit angle
        self.segments = []
        i_start = 0
        t_start = 0
        dist_0 = 0
        z = 0
        self.length = 0
        # n_parts = len(self.parts) - 1
        n_parts = len(self.segs) - 1
        for i, f in enumerate(self.segs):
            f.dist = self.length
            self.length += f.line.length

        vz0 = Vector((1, 0))
        angle_z = 0
        for i, f in enumerate(self.segs):
            dz = self.parts[i].dz
            if f.dist > 0:
                f.t_start = f.dist / self.length
            else:
                f.t_start = 0

            f.t_end = (f.dist + f.line.length) / self.length
            f.z0 = z
            f.dz = dz
            z += dz

            if i < n_parts:

                vz1 = Vector((self.segs[i + 1].length, self.parts[i + 1].dz))
                angle_z = abs(vz0.angle_signed(vz1))
                vz0 = vz1

                if (abs(self.parts[i + 1].a0) >= angle_limit or angle_z >= angle_limit):
                    l_seg = f.dist + f.line.length - dist_0
                    t_seg = f.t_end - t_start
                    n_fences = max(1, int(l_seg / post_spacing))
                    t_fence = t_seg / n_fences
                    segment = FenceSegment(t_start, f.t_end, n_fences, t_fence, i_start, i)
                    dist_0 = f.dist + f.line.length
                    t_start = f.t_end
                    i_start = i
                    self.segments.append(segment)

        f = self.segs[-1]
        l_seg = f.dist + f.line.length - dist_0
        t_seg = f.t_end - t_start
        n_fences = max(1, int(l_seg / post_spacing))
        t_fence = t_seg / n_fences
        segment = FenceSegment(t_start, f.t_end, n_fences, t_fence, i_start, len(self.segs) - 1)
        self.segments.append(segment)

    def setup_user_defined_post(self, o, post_x, post_y, post_z, post_rotation, use_matid, matid):
        self.user_defined_post = o
        x = o.bound_box[6][0] - o.bound_box[0][0]
        y = o.bound_box[6][1] - o.bound_box[0][1]
        z = o.bound_box[6][2] - o.bound_box[0][2]
        # Prevent 0 division error on objects with single vertex
        if x != 0:
            x = post_x / x
        if y != 0:
            y = post_y / y
        if z != 0:
            z = post_z / z
        self.user_defined_post_scale = Vector((x, -y, z))
        m = o.data
        # create vertex group lookup dictionary for names
        vgroup_names = {vgroup.index: vgroup.name for vgroup in o.vertex_groups}
        # create dictionary of vertex group assignments per vertex
        self.vertex_groups = [[vgroup_names[g.group] for g in v.groups] for v in m.vertices]
        # uvs
        uv_act = m.uv_layers.active
        if uv_act is not None:
            uv_layer = uv_act.data
            self.user_defined_uvs = [[uv_layer[li].uv for li in p.loop_indices] for p in m.polygons]
        else:
            self.user_defined_uvs = [[(0, 0) for i in p.vertices] for p in m.polygons]

        # material ids
        if use_matid:
            self.user_defined_mat = [p.material_index for p in m.polygons]
        else:
            self.user_defined_mat = [matid for p in m.polygons]

        ca = cos(post_rotation)
        sa = sin(post_rotation)
        self.user_rM = Matrix([
            [ca, -sa, 0, 0],
            [sa, ca, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])

    def get_user_defined_post(self, tM, z0, z1, z2, slope, post_z, verts, faces, matids, uvs):
        f = len(verts)
        m = self.user_defined_post.data

        for i, g in enumerate(self.vertex_groups):
            co = m.vertices[i].co.copy()
            co.x *= self.user_defined_post_scale.x
            co.y *= self.user_defined_post_scale.y
            co.z *= self.user_defined_post_scale.z
            co = self.user_rM * co
            if 'Slope' in g:
                co.z += co.y * slope
            verts.append(tM * co)
        matids.extend(self.user_defined_mat)
        faces.extend([tuple([i + f for i in p.vertices]) for p in m.polygons])
        uvs.extend(self.user_defined_uvs)

    def get_post(self, post, post_x, post_y, post_z, post_alt, sub_offset_x,
            id_mat, verts, faces, matids, uvs):

        n, dz, zl = post
        slope = dz * post_y

        if self.user_defined_post is not None:
            x, y = -n.v.normalized()
            p = n.p + sub_offset_x * n.v.normalized()
            tM = Matrix([
                [x, y, 0, p.x],
                [y, -x, 0, p.y],
                [0, 0, 1, zl + post_alt],
                [0, 0, 0, 1]
            ])
            self.get_user_defined_post(tM, zl, 0, 0, dz, post_z, verts, faces, matids, uvs)
            return

        z3 = zl + post_z + post_alt - slope
        z4 = zl + post_z + post_alt + slope
        z0 = zl + post_alt - slope
        z1 = zl + post_alt + slope
        vn = n.v.normalized()
        dx = post_x * vn
        dy = post_y * Vector((vn.y, -vn.x))
        oy = sub_offset_x * vn
        x0, y0 = n.p - dx + dy + oy
        x1, y1 = n.p - dx - dy + oy
        x2, y2 = n.p + dx - dy + oy
        x3, y3 = n.p + dx + dy + oy
        f = len(verts)
        verts.extend([(x0, y0, z0), (x0, y0, z3),
                    (x1, y1, z1), (x1, y1, z4),
                    (x2, y2, z1), (x2, y2, z4),
                    (x3, y3, z0), (x3, y3, z3)])
        faces.extend([(f, f + 1, f + 3, f + 2),
                    (f + 2, f + 3, f + 5, f + 4),
                    (f + 4, f + 5, f + 7, f + 6),
                    (f + 6, f + 7, f + 1, f),
                    (f, f + 2, f + 4, f + 6),
                    (f + 7, f + 5, f + 3, f + 1)])
        matids.extend([id_mat, id_mat, id_mat, id_mat, id_mat, id_mat])
        x = [(0, 0), (0, post_z), (post_x, post_z), (post_x, 0)]
        y = [(0, 0), (0, post_z), (post_y, post_z), (post_y, 0)]
        z = [(0, 0), (post_x, 0), (post_x, post_y), (0, post_y)]
        uvs.extend([x, y, x, y, z, z])

    def get_panel(self, subs, altitude, panel_x, panel_z, sub_offset_x, idmat, verts, faces, matids, uvs):
        n_subs = len(subs)
        if n_subs < 1:
            return
        f = len(verts)
        x0 = sub_offset_x - 0.5 * panel_x
        x1 = sub_offset_x + 0.5 * panel_x
        z0 = 0
        z1 = panel_z
        profile = [Vector((x0, z0)), Vector((x1, z0)), Vector((x1, z1)), Vector((x0, z1))]
        user_path_uv_v = []
        n_sections = n_subs - 1
        n, dz, zl = subs[0]
        p0 = n.p
        v0 = n.v.normalized()
        for s, section in enumerate(subs):
            n, dz, zl = section
            p1 = n.p
            if s < n_sections:
                v1 = subs[s + 1][0].v.normalized()
            dir = (v0 + v1).normalized()
            scale = 1 / cos(0.5 * acos(min(1, max(-1, v0 * v1))))
            for p in profile:
                x, y = n.p + scale * p.x * dir
                z = zl + p.y + altitude
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
            user_path_verts=n_subs
            )
        faces += lofter.faces(16, offset=f, path_type='USER_DEFINED')
        matids += lofter.mat(16, idmat, idmat, path_type='USER_DEFINED')
        v = Vector((0, 0))
        uvs += lofter.uv(16, v, v, v, v, 0, v, 0, 0, path_type='USER_DEFINED')

    def make_subs(self, x, y, z, post_y, altitude,
            sub_spacing, offset_x, sub_offset_x, mat, verts, faces, matids, uvs):

        t_post = (0.5 * post_y - y) / self.length
        t_spacing = (sub_spacing + y) / self.length

        for segment in self.segments:
            t_step = segment.t_step
            t_start = segment.t_start + t_post
            s = 0
            s_sub = t_step - 2 * t_post
            n_sub = int(s_sub / t_spacing)
            if n_sub > 0:
                t_sub = s_sub / n_sub
            else:
                t_sub = 1
            i = segment.i_start
            while s < segment.n_step:
                t_cur = t_start + s * t_step
                for j in range(1, n_sub):
                    t_s = t_cur + t_sub * j
                    while self.segs[i].t_end < t_s:
                        i += 1
                    f = self.segs[i]
                    t = (t_s - f.t_start) / f.t_diff
                    n = f.line.normal(t)
                    post = (n, f.dz / f.length, f.z0 + f.dz * t)
                    self.get_post(post, x, y, z, altitude, sub_offset_x, mat, verts, faces, matids, uvs)
                s += 1

    def make_post(self, x, y, z, altitude, x_offset, mat, verts, faces, matids, uvs):
        # skip first port when closed
        make_first_post = not self.closed

        for segment in self.segments:
            t_step = segment.t_step
            t_start = segment.t_start
            s = 0
            i = segment.i_start
            while s < segment.n_step:
                t_cur = t_start + s * t_step
                while self.segs[i].t_end < t_cur:
                    i += 1
                f = self.segs[i]
                t = (t_cur - f.t_start) / f.t_diff
                n = f.line.normal(t)
                post = (n, f.dz / f.line.length, f.z0 + f.dz * t)
                # self.get_post(post, x, y, z, altitude, x_offset, mat, verts, faces, matids, uvs)
                if make_first_post:
                    self.get_post(post, x, y, z, altitude, 0, mat, verts, faces, matids, uvs)
                make_first_post = True
                s += 1

            if segment.i_end + 1 == len(self.segs):
                f = self.segs[segment.i_end]
                n = f.line.normal(1)
                post = (n, f.dz / f.line.length, f.z0 + f.dz)
                # self.get_post(post, x, y, z, altitude, x_offset, mat, verts, faces, matids, uvs)
                self.get_post(post, x, y, z, altitude, 0, mat, verts, faces, matids, uvs)

    def make_panels(self, x, z, post_y, altitude, panel_dist,
            offset_x, sub_offset_x, idmat, verts, faces, matids, uvs):

        t_post = (0.5 * post_y + panel_dist) / self.length
        for segment in self.segments:
            t_step = segment.t_step
            t_start = segment.t_start
            s = 0
            i = segment.i_start
            while s < segment.n_step:
                subs = []
                t_cur = t_start + s * t_step + t_post
                t_end = t_start + (s + 1) * t_step - t_post
                # find first section
                while self.segs[i].t_end < t_cur and i < segment.i_end:
                    i += 1
                f = self.segs[i]
                # 1st section
                t = (t_cur - f.t_start) / f.t_diff
                n = f.line.normal(t)
                subs.append((n, f.dz / f.line.length, f.z0 + f.dz * t))
                # crossing sections -> new segment
                while i < segment.i_end:
                    f = self.segs[i]
                    if f.t_end < t_end:
                        if type(f).__name__ == 'CurvedSegment':
                            # cant end after segment
                            t0 = max(0, (t_cur - f.t_start) / f.t_diff)
                            t1 = min(1, (t_end - f.t_start) / f.t_diff)
                            n_s = int(max(1, abs(f.da) * (5) / pi - 1))
                            dt = (t1 - t0) / n_s
                            for j in range(1, n_s + 1):
                                t = t0 + dt * j
                                n = f.line.sized_normal(t, 1)
                                # n.p = f.lerp(x_offset)
                                subs.append((n, f.dz / f.line.length, f.z0 + f.dz * t))
                        else:
                            n = f.line.normal(1)
                            subs.append((n, f.dz / f.line.length, f.z0 + f.dz))
                    if f.t_end >= t_end:
                        break
                    elif f.t_start < t_end:
                        i += 1

                f = self.segs[i]
                # last section
                if type(f).__name__ == 'CurvedSegment':
                    # cant start before segment
                    t0 = max(0, (t_cur - f.t_start) / f.t_diff)
                    t1 = min(1, (t_end - f.t_start) / f.t_diff)
                    n_s = int(max(1, abs(f.da) * (5) / pi - 1))
                    dt = (t1 - t0) / n_s
                    for j in range(1, n_s + 1):
                        t = t0 + dt * j
                        n = f.line.sized_normal(t, 1)
                        # n.p = f.lerp(x_offset)
                        subs.append((n, f.dz / f.line.length, f.z0 + f.dz * t))
                else:
                    t = (t_end - f.t_start) / f.t_diff
                    n = f.line.normal(t)
                    subs.append((n, f.dz / f.line.length, f.z0 + f.dz * t))

                # self.get_panel(subs, altitude, x, z, 0, idmat, verts, faces, matids, uvs)
                self.get_panel(subs, altitude, x, z, sub_offset_x, idmat, verts, faces, matids, uvs)
                s += 1

    def make_profile(self, profile, idmat,
            x_offset, z_offset, extend, closed, verts, faces, matids, uvs):

        if self.closed:
            extend = 0

        self.set_offset(x_offset)
        self.close(x_offset)

        n_fences = len(self.segs) - 1

        if n_fences < 0:
            return

        sections = []

        f = self.segs[0]

        # first step
        if extend != 0 and f.line.length != 0:
            t = -extend / f.line.length
            n = f.line.sized_normal(t, 1)
            # n.p = f.lerp(x_offset)
            sections.append((n, f.dz / f.line.length, f.z0 + f.dz * t))

        # add first section
        if self.closed:
            n = self.segs[-1].line.sized_normal(1, 1)
        else:
            n = f.line.sized_normal(0, 1)

        sections.append((n, f.dz / f.line.length, f.z0))

        for s, f in enumerate(self.segs):
            if f.line.length == 0:
                continue
            if type(f).__name__ == 'CurvedSegment':
                n_s = int(max(1, abs(f.da) * 30 / pi - 1))
                for i in range(1, n_s + 1):
                    t = i / n_s
                    n = f.line.sized_normal(t, 1)
                    # n.p = f.lerp(x_offset)
                    sections.append((n, f.dz / f.line.length, f.z0 + f.dz * t))
            else:
                n = f.line.sized_normal(1, 1)
                # n.p = f.lerp(x_offset)
                sections.append((n, f.dz / f.line.length, f.z0 + f.dz))

        if extend != 0 and f.line.length != 0:
            t = 1 + extend / f.line.length
            n = f.line.sized_normal(t, 1)
            sections.append((n, f.dz / f.line.length, f.z0 + f.dz * t))

        user_path_verts = len(sections)

        offset = len(verts)
        if user_path_verts > 0:
            user_path_uv_v = []
            # n, dz, z0 = sections[-1]
            # sections[-1] = (n, dz, z0)
            n_sections = user_path_verts - 1

            n, dz, zl = sections[0]
            p0 = n.p
            v0 = n.v.normalized()
            for s, section in enumerate(sections):
                n, dz, zl = section
                p1 = n.p
                if s < n_sections:
                    v1 = sections[s + 1][0].v.normalized()

                if not self.closed or s < n_sections:
                    dir = (v0 + v1).normalized()
                    scale = min(10, 1 / cos(0.5 * acos(min(1, max(-1, v0 * v1)))))
                    for p in profile:
                        # x, y = n.p + scale * (x_offset + p.x) * dir
                        x, y = n.p + scale * p.x * dir
                        z = zl + p.y + z_offset
                        verts.append((x, y, z))

                if s > 0:
                    user_path_uv_v.append((p1 - p0).length)
                p0 = p1
                v0 = v1

            if self.closed:
                user_path_verts -= 1
                n, dz, zl = sections[0]
                p1 = n.p
                user_path_uv_v.append((p1 - p0).length)

            # build faces using Panel
            lofter = Lofter(
                # closed_shape, index, x, y, idmat
                closed,
                [i for i in range(len(profile))],
                [p.x for p in profile],
                [p.y for p in profile],
                [idmat for i in range(len(profile))],
                closed_path=self.closed,
                user_path_uv_v=user_path_uv_v,
                user_path_verts=user_path_verts
                )
            faces += lofter.faces(16, offset=offset, path_type='USER_DEFINED')
            matids += lofter.mat(16, idmat, idmat, path_type='USER_DEFINED')
            v = Vector((0, 0))
            uvs += lofter.uv(16, v, v, v, v, 0, v, 0, 0, path_type='USER_DEFINED')


def update(self, context):
    self.update(context)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


def update_path(self, context):
    self.update_path(context)


materials_enum = (
            ('0', 'Wood', '', 0),
            ('1', 'Metal', '', 1),
            ('2', 'Glass', '', 2)
            )


class archipack_fence_part(ArchipackSegment, PropertyGroup):
    dz = FloatProperty(
            name="delta z",
            default=0,
            unit='LENGTH', subtype='DISTANCE'
            )
    manipulators = CollectionProperty(type=archipack_manipulator)

    def get_datablock(self, o):
        return archipack_fence.datablock(o)


class archipack_fence_rail(ArchipackObjectsManager, ArchipackProfile, PropertyGroup):
    profil_x = FloatProperty(
            name="Width",
            min=0.001,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    profil_y = FloatProperty(
            name="Height",
            min=0.001,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    offset = FloatProperty(
            name="Offset",
            default=0,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            )
    extend = FloatProperty(
            name="Extend",
            default=0,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            )
    alt = FloatProperty(
            name="Altitude",
            default=1.0,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            )
    profil = EnumProperty(
            name="Profil",
            items=(
                ('SQUARE', 'Square', '', 0),
                ('CIRCLE', 'Circle', '', 1),
                ('SAFETY', 'Safety rail', '', 2),
                ('USER', 'User defined', '', 3)
                ),
            default='SQUARE',
            update=update
            )
    mat = EnumProperty(
        items=materials_enum,
        default='0',
        update=update
        )
    auto_update = BoolProperty(
        options={'SKIP_SAVE'},
        default=True
        )

    def refresh_profile_size(self, context, x, y):
        self.profil_x = x
        self.profil_y = y
        self.auto_update = True
        self.update(context)

    def find_datablock_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_fence.datablock(o)
            if props is not None:
                for rail in props.rails:
                    if rail == self:
                        return props
        return None

    def update(self, context, manipulable_refresh=False):
        if self.auto_update:
            props = self.find_datablock_in_selection(context)
            if props is not None:
                props.update(context, manipulable_refresh)

    def draw(self, context, layout):
        layout.prop(self, 'profil')
        if self.profil == 'USER':
            self.draw_user_profile(context, layout)
        layout.prop(self, 'profil_x')
        layout.prop(self, 'profil_y')
        layout.prop(self, 'alt')
        layout.prop(self, 'offset')
        layout.prop(self, 'extend')
        layout.prop(self, 'mat')


class archipack_fence(ArchipackObject, ArchipackUserDefinedPath, Manipulable, DimensionProvider, PropertyGroup):

    parts = CollectionProperty(type=archipack_fence_part)

    """
    n_parts = IntProperty(
            name="Parts",
            min=1,
            default=1, update=update_manipulators
            )
    """
    x_offset = FloatProperty(
            name="Offset",
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    radius = FloatProperty(
            name="Radius",
            min=0.01,
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
    angle_limit = FloatProperty(
            name="Angle",
            min=0,
            max=2 * pi,
            default=pi / 8,
            subtype='ANGLE', unit='ROTATION',
            update=update_manipulators
            )
    shape = EnumProperty(
            items=(
                ('RECTANGLE', 'Straight', '', 0),
                ('CIRCLE', 'Curved ', '', 1)
                ),
            default='RECTANGLE',
            update=update
            )
    post = BoolProperty(
            name='Enable',
            default=True,
            update=update
            )
    post_spacing = FloatProperty(
            name="Spacing",
            min=0.1,
            default=1.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    post_x = FloatProperty(
            name="Width",
            min=0.001,
            default=0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    post_y = FloatProperty(
            name="Length",
            min=0.001, max=1000,
            default=0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    post_z = FloatProperty(
            name="Height",
            min=0.001,
            default=1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    post_alt = FloatProperty(
            name="Altitude",
            default=0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    post_rotation = FloatProperty(
            name="Rotation",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    user_defined_post_enable = BoolProperty(
            name="User",
            update=update,
            default=True
            )
    user_defined_post = StringProperty(
            name="User defined",
            update=update
            )
    idmat_post = EnumProperty(
            name="Post",
            items=materials_enum,
            default='1',
            update=update
            )
    subs = BoolProperty(
            name='Enable',
            default=False,
            update=update
            )
    subs_spacing = FloatProperty(
            name="Spacing",
            min=0.05,
            default=0.10, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    subs_x = FloatProperty(
            name="Width",
            min=0.001,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    subs_y = FloatProperty(
            name="Length",
            min=0.001,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    subs_z = FloatProperty(
            name="Height",
            min=0.001,
            default=1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    subs_alt = FloatProperty(
            name="Altitude",
            default=0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    subs_offset_x = FloatProperty(
            name="Offset",
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    subs_bottom = EnumProperty(
            name="Bottom",
            items=(
                ('STEP', 'Follow step', '', 0),
                ('LINEAR', 'Linear', '', 1),
                ),
            default='STEP',
            update=update
            )
    subs_rotation = FloatProperty(
            name="Rotation",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    user_defined_subs_enable = BoolProperty(
            name="User",
            update=update,
            default=True
            )
    user_defined_subs = StringProperty(
            name="User defined",
            update=update
            )
    idmat_subs = EnumProperty(
            name="Subs",
            items=materials_enum,
            default='1',
            update=update
            )
    panel = BoolProperty(
            name='Enable',
            default=True,
            update=update
            )
    panel_alt = FloatProperty(
            name="Altitude",
            default=0.25, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    panel_x = FloatProperty(
            name="Width",
            min=0.001,
            default=0.01, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    panel_z = FloatProperty(
            name="Height",
            min=0.001,
            default=0.6, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    panel_dist = FloatProperty(
            name="Spacing",
            min=0.001,
            default=0.05, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    panel_offset_x = FloatProperty(
            name="Offset",
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    idmat_panel = EnumProperty(
            name="Panels",
            items=materials_enum,
            default='2',
            update=update
            )
    rail = BoolProperty(
            name="Enable",
            update=update,
            default=False
            )
    rail_n = IntProperty(
            name="#",
            default=1,
            min=0,
            max=31,
            update=update
            )
    rails = CollectionProperty(type=archipack_fence_rail)
    """
    # Rails v1.x might re-enable to make old presets loadable
    # but not certain it is a good thing
    rail_x = FloatVectorProperty(
            name="Width",
            default=[
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05
            ],
            size=31,
            min=0.001,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            )
    rail_z = FloatVectorProperty(
            name="Height",
            default=[
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05
            ],
            size=31,
            min=0.001,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            )
    rail_offset = FloatVectorProperty(
            name="Offset",
            default=[
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0
            ],
            size=31,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            )
    rail_alt = FloatVectorProperty(
            name="Altitude",
            default=[
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
            ],
            size=31,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            )
    rail_type = IntVectorProperty(
            name="Profil",
            default=[1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1
            ],
            size=31,
            update=update
            )
    rail_mat = CollectionProperty(type=archipack_fence_material)
    """
    handrail = BoolProperty(
            name="Enable",
            update=update,
            default=True
            )
    handrail_offset = FloatProperty(
            name="Offset",
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handrail_alt = FloatProperty(
            name="Altitude",
            default=1.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handrail_extend = FloatProperty(
            name="Extend",
            min=0,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handrail_slice = BoolProperty(
            name='Slice',
            default=True,
            update=update
            )
    handrail_slice_right = BoolProperty(
            name='Slice',
            default=True,
            update=update
            )
    handrail_profil = EnumProperty(
            name="Profil",
            items=(
                ('SQUARE', 'Square', '', 0),
                ('CIRCLE', 'Circle', '', 1),
                ('COMPLEX', 'Circle over square', '', 2)
                ),
            default='SQUARE',
            update=update
            )
    handrail_x = FloatProperty(
            name="Width",
            min=0.001,
            default=0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handrail_y = FloatProperty(
            name="Height",
            min=0.001,
            default=0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handrail_radius = FloatProperty(
            name="Radius",
            min=0.001,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    idmat_handrail = EnumProperty(
            name="Handrail",
            items=materials_enum,
            default='0',
            update=update
            )
    use_subs_material = BoolProperty(
            name="Use material",
            description="Use material indexes of object",
            update=update,
            default=False)

    use_post_material = BoolProperty(
            name="Use material",
            description="Use material indexes of object",
            update=update,
            default=False)

    # UI layout related
    parts_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    rail_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    idmats_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    handrail_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    post_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    panel_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    subs_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    closed = BoolProperty(
            default=False,
            update=update_manipulators
            )
    always_closed = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )

    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )

    def setup_manipulators(self):

        if len(self.manipulators) == 0:
            s = self.manipulators.add()
            s.prop1_name = "width"
            s = self.manipulators.add()
            s.prop1_name = "height"
            s.normal = Vector((0, 1, 0))
        self.setup_parts_manipulators('dumb_z')

    def update_parts(self):

        # remove rails materials
        for i in range(len(self.rails), self.rail_n, -1):
            self.rails.remove(i - 1)

        # add rails
        for i in range(len(self.rails), self.rail_n):
            self.rails.add()

        return ArchipackUserDefinedPath.update_parts(self)

    def from_spline(self, context, wM, resolution, spline):

        o = self.find_in_selection(context)

        if o is None:
            return

        pts = self.coords_from_spline(
                spline,
                wM,
                resolution
                )

        if len(pts) < 2:
            return

        if self.user_defined_reverse:
            pts = list(reversed(pts))

        o.matrix_world = Matrix.Translation(pts[0].copy())

        auto_update = self.auto_update
        self.auto_update = False
        self.closed = spline.use_cyclic_u
        if self.closed:
            pts.pop()
        self.from_points(pts)
        self.auto_update = auto_update

    def get_generator(self, o=None):
        g = FenceGenerator(self, o)
        for part in self.parts:
            # type, radius, da, length
            g.add_part(part)
        g.set_offset(self.x_offset)
        g.close(self.x_offset)
        return g

    def update(self, context, manipulable_refresh=False):

        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # clean up manipulators before any data model change
        if manipulable_refresh:
            self.manipulable_disable(context)

        self.update_parts()

        verts = []
        faces = []
        matids = []
        uvs = []

        g = self.get_generator()
        g.locate_manipulators()

        if not self.closed:
            g.segs.pop()

        g.param_t(self.angle_limit, self.post_spacing)

        # depth at bottom
        # self.manipulators[1].set_pts([(0, 0, 0), (0, 0, self.height), (1, 0, 0)])

        if self.user_defined_post_enable:
            # user defined posts
            user_def_post = context.scene.objects.get(self.user_defined_post)
            if user_def_post is not None and user_def_post.type == 'MESH':
                g.setup_user_defined_post(user_def_post,
                    self.post_x, self.post_y, self.post_z, self.post_rotation,
                    self.use_post_material, int(self.idmat_post))

        if self.post:
            g.make_post(0.5 * self.post_x, 0.5 * self.post_y, self.post_z,
                    self.post_alt, self.x_offset,
                    int(self.idmat_post), verts, faces, matids, uvs)

        # reset user def posts
        g.user_defined_post = None

        # user defined subs
        if self.user_defined_subs_enable:
            user_def_subs = context.scene.objects.get(self.user_defined_subs)
            if user_def_subs is not None and user_def_subs.type == 'MESH':
                g.setup_user_defined_post(user_def_subs,
                    self.subs_x, self.subs_y, self.subs_z, self.subs_rotation,
                    self.use_subs_material, int(self.idmat_subs))

        if self.subs:
            g.make_subs(0.5 * self.subs_x, 0.5 * self.subs_y, self.subs_z,
                    self.post_y, self.subs_alt, self.subs_spacing,
                    self.x_offset, self.subs_offset_x, int(self.idmat_subs), verts, faces, matids, uvs)

        g.user_defined_post = None

        if self.panel:
            g.make_panels(0.5 * self.panel_x, self.panel_z, self.post_y,
                    self.panel_alt, self.panel_dist, self.x_offset, self.panel_offset_x,
                    int(self.idmat_panel), verts, faces, matids, uvs)

        if self.rail:
            for i in range(self.rail_n):
                rd = self.rails[i]
                x = 0.5 * rd.profil_x
                y = rd.profil_y
                closed = True

                if rd.profil == 'SQUARE':
                    rail = [Vector((-x, y)), Vector((-x, 0)), Vector((x, 0)), Vector((x, y))]
                elif rd.profil == 'CIRCLE':
                    rail = [Vector((x * sin(0.1 * -a * pi), x * (0.5 + cos(0.1 * -a * pi)))) for a in range(0, 20)]

                elif rd.profil == 'SAFETY':
                    closed = False
                    rail = [Vector((i * x, j * y)) for i, j in [(0, -0.5),
                        (1, -0.35714),
                        (1, -0.21429),
                        (0, -0.07143),
                        (0, 0.07143),
                        (1, 0.21429),
                        (1, 0.35714),
                        (0, 0.5)]]
                elif rd.profil == 'USER':
                    curve = rd.update_profile(context)
                    if curve and curve.type == 'CURVE':
                        sx, sy = 1, 1
                        if rd.user_profile_dimension.x > 0:
                            sx = rd.profil_x / rd.user_profile_dimension.x
                        if rd.user_profile_dimension.y > 0:
                            sy = rd.profil_y / rd.user_profile_dimension.y
                        wM = Matrix([
                            [sx, 0, 0, 0],
                            [0, sy, 0, 0],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1]
                            ])
                        for spline in curve.data.splines:
                            rail = self.coords_from_spline(spline, wM, 12, ccw=True)
                            closed = rail[0] == rail[-1]
                            if closed:
                                rail.pop()
                            g.make_profile(rail, int(rd.mat), self.x_offset - rd.offset,
                                rd.alt, rd.extend, closed, verts, faces, matids, uvs)
                    else:
                        print("fence.update curve not found")

                    # dont call
                    continue

                g.make_profile(rail, int(rd.mat), self.x_offset - rd.offset,
                        rd.alt, rd.extend, closed, verts, faces, matids, uvs)

        if self.handrail:

            if self.handrail_profil == 'COMPLEX':
                sx = self.handrail_x
                sy = self.handrail_y
                handrail = [Vector((sx * x, sy * y)) for x, y in [
                (-0.28, 1.83), (-0.355, 1.77), (-0.415, 1.695), (-0.46, 1.605), (-0.49, 1.51), (-0.5, 1.415),
                (-0.49, 1.315), (-0.46, 1.225), (-0.415, 1.135), (-0.355, 1.06), (-0.28, 1.0), (-0.255, 0.925),
                (-0.33, 0.855), (-0.5, 0.855), (-0.5, 0.0), (0.5, 0.0), (0.5, 0.855), (0.33, 0.855), (0.255, 0.925),
                (0.28, 1.0), (0.355, 1.06), (0.415, 1.135), (0.46, 1.225), (0.49, 1.315), (0.5, 1.415),
                (0.49, 1.51), (0.46, 1.605), (0.415, 1.695), (0.355, 1.77), (0.28, 1.83), (0.19, 1.875),
                (0.1, 1.905), (0.0, 1.915), (-0.095, 1.905), (-0.19, 1.875)]]

            elif self.handrail_profil == 'SQUARE':
                x = 0.5 * self.handrail_x
                y = self.handrail_y
                handrail = [Vector((-x, y)), Vector((-x, 0)), Vector((x, 0)), Vector((x, y))]

            elif self.handrail_profil == 'CIRCLE':
                r = self.handrail_radius
                handrail = [Vector((r * sin(0.1 * -a * pi), r * (0.5 + cos(0.1 * -a * pi)))) for a in range(0, 20)]

            g.make_profile(handrail, int(self.idmat_handrail), self.x_offset - self.handrail_offset,
                self.handrail_alt, self.handrail_extend, True, verts, faces, matids, uvs)

        bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=uvs, weld=True, clean=False)

        self.update_dimensions(context, o)

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

        n_parts = self.n_parts
        if self.closed:
            n_parts += 1

        for i, part in enumerate(self.parts):
            if i < n_parts:
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


class ARCHIPACK_PT_fence(Panel):
    bl_idname = "ARCHIPACK_PT_fence"
    bl_label = "Fence"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_fence.poll(context.active_object)

    def draw(self, context):
        prop = archipack_fence.datablock(context.active_object)
        if prop is None:
            return
        scene = context.scene
        layout = self.layout
        row = layout.row(align=True)
        row.operator('archipack.manipulate', icon='HAND')
        box = layout.box()
        # box.label(text="Styles")
        row = box.row(align=True)
        row.operator("archipack.fence_preset_menu", text=bpy.types.ARCHIPACK_OT_fence_preset_menu.bl_label)
        row.operator("archipack.fence_preset", text="", icon='ZOOMIN')
        row.operator("archipack.fence_preset", text="", icon='ZOOMOUT').remove_active = True

        box = layout.box()
        expand = prop.template_user_path(context, box)
        if expand:
            if prop.user_defined_path is not "":
                box.prop(prop, 'user_defined_reverse')

            box.prop(prop, 'angle_limit')

        box = layout.box()
        box.prop(prop, 'x_offset')

        prop.template_parts(context, layout, draw_type=True)

        box = layout.box()
        row = box.row(align=True)
        icon = "TRIA_RIGHT"
        if prop.handrail_expand:
            icon = "TRIA_DOWN"

        row.prop(prop, 'handrail_expand', icon=icon, icon_only=True, text="Handrail", emboss=True)
        row.prop(prop, 'handrail')

        if prop.handrail_expand:
            box.prop(prop, 'handrail_alt')
            box.prop(prop, 'handrail_offset')
            box.prop(prop, 'handrail_extend')
            box.prop(prop, 'handrail_profil')
            if prop.handrail_profil != 'CIRCLE':
                box.prop(prop, 'handrail_x')
                box.prop(prop, 'handrail_y')
            else:
                box.prop(prop, 'handrail_radius')
            row = box.row(align=True)
            row.prop(prop, 'handrail_slice')

        box = layout.box()
        row = box.row(align=True)
        icon = "TRIA_RIGHT"
        if prop.post_expand:
            icon = "TRIA_DOWN"

        row.prop(prop, 'post_expand', icon=icon, icon_only=True, text="Post", emboss=True)
        row.prop(prop, 'post')

        if prop.post_expand:
            box.prop(prop, 'post_spacing')
            box.prop(prop, 'post_x')
            box.prop(prop, 'post_y')
            box.prop(prop, 'post_z')
            box.prop(prop, 'post_alt')
            row = box.row(align=True)
            row.prop(prop, 'user_defined_post_enable', text="")
            row.operator("archipack.fence_subpart_dimensions", text="", icon='BBOX').part = "POST"
            row.prop_search(prop, "user_defined_post", scene, "objects", text="")
            if prop.user_defined_post:
                box.prop(prop, 'post_rotation')
                box.prop(prop, 'use_post_material')

        box = layout.box()
        row = box.row(align=True)
        icon = "TRIA_RIGHT"
        if prop.subs_expand:
            icon = "TRIA_DOWN"

        row.prop(prop, 'subs_expand', icon=icon, icon_only=True, text="Subs", emboss=True)
        row.prop(prop, 'subs')

        if prop.subs_expand:
            box.prop(prop, 'subs_spacing')
            box.prop(prop, 'subs_x')
            box.prop(prop, 'subs_y')
            box.prop(prop, 'subs_z')
            box.prop(prop, 'subs_alt')
            box.prop(prop, 'subs_offset_x')
            row = box.row(align=True)
            row.prop(prop, 'user_defined_subs_enable', text="")
            row.operator("archipack.fence_subpart_dimensions", text="", icon='BBOX').part = "SUB"
            row.prop_search(prop, "user_defined_subs", scene, "objects", text="")
            if prop.user_defined_subs:
                box.prop(prop, 'subs_rotation')
                box.prop(prop, 'use_subs_material')

        box = layout.box()
        row = box.row(align=True)
        icon = "TRIA_RIGHT"
        if prop.panel_expand:
            icon = "TRIA_DOWN"

        row.prop(prop, 'panel_expand', icon=icon, icon_only=True, text="Panels", emboss=True)
        row.prop(prop, 'panel')

        if prop.panel_expand:
            box.prop(prop, 'panel_dist')
            box.prop(prop, 'panel_x')
            box.prop(prop, 'panel_z')
            box.prop(prop, 'panel_alt')
            box.prop(prop, 'panel_offset_x')

        box = layout.box()
        row = box.row(align=True)
        icon = "TRIA_RIGHT"
        if prop.rail_expand:
            icon = "TRIA_DOWN"

        row.prop(prop, 'rail_expand', icon=icon, icon_only=True, text="Rails", emboss=True)
        row.prop(prop, 'rail')

        if prop.rail_expand:
            box.prop(prop, 'rail_n')
            for i in range(prop.rail_n):
                box = layout.box()
                box.label(text="Rail " + str(i + 1))
                prop.rails[i].draw(context, box)

        box = layout.box()
        row = box.row()
        icon = "TRIA_RIGHT"
        if prop.idmats_expand:
            icon = "TRIA_DOWN"

        row.prop(prop, 'idmats_expand', icon=icon, icon_only=True, text="Material index", emboss=True)

        if prop.idmats_expand:
            box.prop(prop, 'idmat_handrail')
            box.prop(prop, 'idmat_panel')
            box.prop(prop, 'idmat_post')
            box.prop(prop, 'idmat_subs')


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_fence(ArchipackCreateTool, Operator):
    bl_idname = "archipack.fence"
    bl_label = "Fence"
    bl_description = "Fence"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    def create(self, context):
        m = bpy.data.meshes.new("Fence")
        o = bpy.data.objects.new("Fence", m)
        d = m.archipack_fence.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        # Link object into scene
        self.link_object_to_scene(context, o)

        # select and make active
        self.select_object(context, o, True)
        self.load_preset(d)
        self.add_material(o)
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = bpy.context.scene.cursor_location
            # select and make active
            self.select_object(context, o, True)
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


"""
 TODO:
 make fence child of curve
 use relationship to update (add/remove) childs
"""


class ARCHIPACK_OT_fence_from_curve(ArchipackCreateTool, Operator):
    bl_idname = "archipack.fence_from_curve"
    bl_label = "Fence curve"
    bl_description = "Create fence(s) from a curve"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o is not None and o.type == 'CURVE'

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create_one(self, context, curve, i):
        bpy.ops.archipack.fence('INVOKE_DEFAULT', auto_manipulate=False)
        o = context.active_object
        d = archipack_fence.datablock(o)
        d.auto_update = False
        d.user_defined_spline = i
        d.user_defined_path = curve.name
        d.user_defined_resolution = min(128, curve.data.resolution_u)
        d.auto_update = True
        return o

    def create(self, context):
        o = None
        curve = context.active_object
        sel = []
        for i, spline in enumerate(curve.data.splines):
            o = self.create_one(context, curve, i)
            sel.append(o)
        for obj in sel:
            self.select_object(context, o, True)
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)

            # select and make active
            self.select_object(context, o, True)
            # self.manipulate()
            return {'FINISHED'}

        elif context.mode == 'EDIT_CURVE':
            # Tynkatopi's contrib
            curve = context.active_object
            d = curve.data
            spline_index = 0

            if d.splines.active:
                spline_index = d.splines[:].index(d.splines.active)

            rot = curve.rotation_euler.copy()
            curve.rotation_euler = 0, 0, 0
            bpy.ops.object.mode_set(mode='OBJECT')
            o = self.create_one(context, curve, 0)
            if o:
                d = archipack_fence.datablock(o)
                o.rotation_euler = 0, 0, 0
                o.parent = curve
                # update here so fence location match spline vertex 1
                d.user_defined_spline = spline_index
                self.select_object(context, o)

            # select and make active
            self.select_object(context, curve, True)
            curve.rotation_euler = rot
            bpy.ops.object.mode_set(mode='EDIT')
            return {'FINISHED'}

        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object/Edit mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to load / save presets
# ------------------------------------------------------------------


class ARCHIPACK_OT_fence_preset_create(PresetMenuOperator, Operator):
    bl_description = "Show Fence presets and create object at cursor location"
    bl_idname = "archipack.fence_preset_create"
    bl_label = "Fence Styles"
    preset_subdir = "archipack_fence"


class ARCHIPACK_OT_fence_preset_menu(PresetMenuOperator, Operator):
    bl_description = "Show Fence Presets"
    bl_idname = "archipack.fence_preset_menu"
    bl_label = "Fence Styles"
    preset_subdir = "archipack_fence"


class ARCHIPACK_OT_fence_preset(ArchipackPreset, Operator):
    """Add a Fence Preset"""
    bl_idname = "archipack.fence_preset"
    bl_label = "Add Fence Style"
    preset_menu = "ARCHIPACK_OT_fence_preset_menu"

    @property
    def blacklist(self):
        return ['manipulators', 'n_parts', 'parts', 'user_defined_path', 'user_defined_spline']


class ARCHIPACK_OT_fence_subpart_dimensions(Operator):

    bl_idname = "archipack.fence_subpart_dimensions"
    bl_label = "Dimension"
    bl_description = "Use object's dimensions"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    part = StringProperty()

    @classmethod
    def poll(cls, context):
        return archipack_fence.filter(context.active_object)

    def execute(self, context):
        d = archipack_fence.datablock(context.active_object)

        if d is None:
            self.report({'WARNING'}, "Archipack: Operator only valid with fence")
            return {'CANCELLED'}

        if self.part == "SUB":
            part_obj = bpy.data.objects.get(d.user_defined_subs)
            if part_obj is None:
                self.report({'WARNING'}, "Archipack: User defined sub object not found")
                return {'CANCELLED'}
            d.subs_x, d.subs_y, d.subs_z = part_obj.dimensions.x, part_obj.dimensions.y, part_obj.dimensions.z
        else:
            part_obj = bpy.data.objects.get(d.user_defined_post)
            if part_obj is None:
                self.report({'WARNING'}, "Archipack: User defined post object not found")
                return {'CANCELLED'}
            d.post_x, d.post_y, d.post_z = part_obj.dimensions.x, part_obj.dimensions.y, part_obj.dimensions.z

        return {'FINISHED'}


def register():
    # bpy.utils.register_class(archipack_fence_material)
    bpy.utils.register_class(archipack_fence_rail)
    bpy.utils.register_class(archipack_fence_part)
    bpy.utils.register_class(archipack_fence)
    Mesh.archipack_fence = CollectionProperty(type=archipack_fence)
    bpy.utils.register_class(ARCHIPACK_OT_fence_preset_menu)
    bpy.utils.register_class(ARCHIPACK_PT_fence)
    bpy.utils.register_class(ARCHIPACK_OT_fence)
    bpy.utils.register_class(ARCHIPACK_OT_fence_preset)
    bpy.utils.register_class(ARCHIPACK_OT_fence_preset_create)
    bpy.utils.register_class(ARCHIPACK_OT_fence_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_fence_subpart_dimensions)


def unregister():
    # bpy.utils.unregister_class(archipack_fence_material)
    bpy.utils.unregister_class(archipack_fence_rail)
    bpy.utils.unregister_class(archipack_fence_part)
    bpy.utils.unregister_class(archipack_fence)
    del Mesh.archipack_fence
    bpy.utils.unregister_class(ARCHIPACK_OT_fence_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_PT_fence)
    bpy.utils.unregister_class(ARCHIPACK_OT_fence)
    bpy.utils.unregister_class(ARCHIPACK_OT_fence_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_fence_preset_create)
    bpy.utils.unregister_class(ARCHIPACK_OT_fence_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_fence_subpart_dimensions)
