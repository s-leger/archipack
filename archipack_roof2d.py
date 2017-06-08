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
from bpy.types import Operator, PropertyGroup, Mesh, Panel, Menu
from bpy.props import (
    FloatProperty, BoolProperty, IntProperty,
    StringProperty, EnumProperty, FloatVectorProperty,
    CollectionProperty, PointerProperty
    )
from .bmesh_utils import BmeshEdit as bmed
from .materialutils import MaterialUtils
from .panel import Panel as Lofter
from mathutils import Vector, Matrix
from mathutils.geometry import interpolate_bezier
from math import sin, cos, pi, acos, atan2
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_2d import Line, Arc
from .archipack_preset import ArchipackPreset
from .archipack_object import ArchipackCreateTool, ArchipackObject


class Roof():

    def __init__(self):
        # total distance from start
        self.dist = 0
        self.t_start = 0
        self.t_end = 0
        self.dz = 0
        self.z0 = 0
        self.a0 = 0

    def set_offset(self, offset, last=None):
        """
            Offset line and compute intersection point between
            straight segments
        """
        if last is not None:
            last = last.line
        self.line = self.make_offset(offset, last)
        
    @property
    def t_diff(self):
        return self.t_end - self.t_start

    def straight_roof(self, a0, length):
        s = self.straight(length).rotate(a0)
        return StraightRoof(s.p, s.v)

    def curved_roof(self, a0, da, radius):
        n = self.normal(1).rotate(a0).scale(radius)
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedRoof(c, radius, a0, da)


class StraightRoof(Roof, Line):
    def __str__(self):
        return "t_start:{} t_end:{} dist:{}".format(self.t_start, self.t_end, self.dist)

    def __init__(self, p, v):
        Roof.__init__(self)
        Line.__init__(self, p, v)


class CurvedRoof(Roof, Arc):
    def __str__(self):
        return "t_start:{} t_end:{} dist:{}".format(self.t_start, self.t_end, self.dist)

    def __init__(self, c, radius, a0, da):
        Roof.__init__(self)
        Arc.__init__(self, c, radius, a0, da)


class RoofGenerator():

    def __init__(self, parts):
        self.parts = parts
        self.segs = []
        self.length = 0
        self.origin = Vector((0, 0))
        self.user_defined_post = None
        self.user_defined_uvs = None
        self.user_defined_mat = None

    def add_part(self, part, origin):
        
        if len(self.segs) < 1:
            s = None
        else:
            s = self.segs[-1]
        last = s
        # start a new roof
        if s is None:
            self.origin = origin
            if part.type == 'S_SEG':
                v = part.length * Vector((cos(part.a0), sin(part.a0)))
                s = StraightRoof(origin, v)
            elif part.type == 'C_SEG':
                c = origin - part.radius * Vector((cos(part.a0), sin(part.a0)))
                s = CurvedRoof(c, part.radius, part.a0, part.da)
        else:
            if part.type == 'S_SEG':
                s = s.straight_roof(part.a0, part.length)
            elif part.type == 'C_SEG':
                s = s.curved_roof(part.a0, part.da, part.radius)
        
        self.segs.append(s)
        self.last_type = part.type
        if part == self.parts[-1] or part == self.parts[1]:
            s.set_offset(part.offset)
        else:
            s.set_offset(part.offset, last=last)
        
    def close(self, closed):
        # Make last segment implicit closing one
        if closed:
            part = self.parts[-1]
            w = self.segs[-1]
            # dp = self.segs[0].p0 - self.segs[-1].p0
            p1 = self.segs[0].line.p1
            self.segs[0].line = self.segs[0].make_offset(self.parts[0].offset, w.line)
            self.segs[0].line.p1 = p1
            
            return
            
            if "C_" in part.type:
                """
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
                """
            else:
                w.v = dp
                w.line = w.offset(part.offset)
            
            # wont work with arc
            if len(self.segs) > 2:
                i, p, t = self.segs[-2].line.intersect(w.line)
                if i:
                    if type(self.segs[-2]).__name__ == 'StraightRoof':
                        self.segs[-2].line.p1 = p
                    if type(w).__name__ == 'StraightRoof':
                        w.line.p0 = p
                    
            i, p, t = self.segs[0].line.intersect(w.line)
            if i:
                if type(self.segs[0]).__name__ == 'StraightRoof':
                    self.segs[0].line.p0 = p
                if type(w).__name__ == 'StraightRoof':
                    w.line.p1 = p
            
    def param_t(self):
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
            if s > 0 and s < len(self.segs) - 1:
                for i in range(33):
                    x, y = roof.line.lerp(i / 32)
                    verts.append((x, y, 0))


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
            g = WallGenerator(None)
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
            min=3,
            default=3, update=update_manipulators
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
        for i in range(len(self.parts), self.n_parts, -1):
            row_change = True
            self.parts.remove(i - 1)

        # add rows
        for i in range(len(self.parts), self.n_parts):
            self.parts.add()

        self.setup_parts_manipulators()
          
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
            if n_manips < 5:
                s = p.manipulators.add()
                s.type_key = "SIZE"
                s.prop1_name = "offset"
            p.manipulators[2].prop1_name = str(i)
            p.manipulators[3].prop1_name = str(i + 1)
    
    def get_generator(self, origin=Vector((0, 0))):
        g = RoofGenerator(self.parts)
        for part in self.parts:
            # type, radius, da, length
            g.add_part(part, origin)

        g.close(self.closed)    
        g.param_t()
        return g

    
class archipack_roof_boundary_segment(ArchipackSegment, PropertyGroup):
    
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_roof.datablock(o)
            if props:
                for part in props.boundary.parts:
                    if part == self:
                        return props.boundary
        return None


class archipack_roof_boundary(ArchipackLines, PropertyGroup):
    auto_update = BoolProperty(default=True)
    parts = CollectionProperty(type=archipack_roof_boundary_segment)
    closed = BoolProperty(
            default=False,
            name="Close",
            update=update_manipulators
            )
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_roof.datablock(o)
            if props and props.boundary == self:
                    return props
        return None


class archipack_roof_axis_segment(ArchipackSegment, PropertyGroup):
    
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_roof.datablock(o)
            if props:
                for part in props.axis.parts:
                    if part == self:
                        return props.axis
        return None


class archipack_roof_axis(ArchipackLines, PropertyGroup):
    auto_update = BoolProperty(default=True)
    parts = CollectionProperty(type=archipack_roof_axis_segment)
    closed = BoolProperty(
            default=False,
            name="Close",
            update=update_manipulators
            )
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_roof.datablock(o)
            if props and props.axis == self:
                    return props
        return None

        
class archipack_roof(ArchipackObject, Manipulable, PropertyGroup):
    # boundary
    boundary = PointerProperty(type=archipack_roof_boundary)
    axis = PointerProperty(type=archipack_roof_axis)
    
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
            default=0.01, precision=2, step=1,
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
        
    def update_parts(self):
        # print("update_parts")
        # remove rows
        # NOTE:
        # n_parts+1
        # as last one is end point of last segment or closing one
        self.boundary.update_parts()
        self.axis.update_parts()
        self.setup_manipulators()
        
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

    def update(self, context, manipulable_refresh=False):

        active, selected, o = self.find_in_selection(context)

        if o is None or not self.auto_update:
            return

        # clean up manipulators before any data model change
        if manipulable_refresh:
            self.manipulable_disable(context)

        self.update_parts()

        verts = []
        faces = []
        matids = []
        uvs = []
        
        # NOTE:
        # first segment is explicit origin
        # and is not ment to be a part of result
        # both axis and boundary do share origin
        
        g = self.boundary.get_generator()
        g.debug(verts)
        
        g = self.axis.get_generator()
        g.debug(verts)
        
        
        bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=uvs, weld=True, clean=False)

        # enable manipulators rebuild
        if manipulable_refresh:
            self.manipulable_refresh = True

        # restore context
        try:
            for o in selected:
                o.select = True
        except:
            pass

        active.select = True
        context.scene.objects.active = active

    def manipulable_setup(self, context):
        """
            TODO: Implement the setup part as per parent object basis

            self.manipulable_disable(context)
            o = context.active_object
            for m in self.manipulators:
                self.manip_stack.append(m.setup(context, o, self))

        """
        self.manipulable_disable(context)
        o = context.active_object
        d = self
        
        self.setup_manipulators()
        self.axis.setup_parts_manipulators()
        self.boundary.setup_parts_manipulators()

        for i, part in enumerate(d.boundary.parts):
            if i >= d.boundary.n_parts:
                break

            if i > 0 and i < d.boundary.n_parts - 1:
                # start angle
                self.manip_stack.append(part.manipulators[0].setup(context, o, part))

                # length / radius + angle
                self.manip_stack.append(part.manipulators[1].setup(context, o, part))
                # offset
                self.manip_stack.append(part.manipulators[4].setup(context, o, part))
                # index
                self.manip_stack.append(part.manipulators[3].setup(context, o, d.boundary))
                
            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, d.boundary))
            
        for i, part in enumerate(d.axis.parts):
            if i >= d.axis.n_parts:
                break

            if i > 0 and i < d.axis.n_parts - 1:
                # start angle
                self.manip_stack.append(part.manipulators[0].setup(context, o, part))
                
                # length / radius + angle
                self.manip_stack.append(part.manipulators[1].setup(context, o, part))
                # index
                self.manip_stack.append(part.manipulators[3].setup(context, o, d.axis))
                
            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, d.axis))
            # offset
            # self.manip_stack.append(part.manipulators[4].setup(context, o, part))

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
        # box.label(text="Styles")
        row = box.row(align=True)
        row.menu("ARCHIPACK_MT_roof_preset", text=bpy.types.ARCHIPACK_MT_roof_preset.bl_label)
        row.operator("archipack.roof_preset", text="", icon='ZOOMIN')
        row.operator("archipack.roof_preset", text="", icon='ZOOMOUT').remove_active = True
        row = layout.row(align=True)
        box = layout.box()
        row.prop_search(prop, "user_defined_path", scene, "objects", text="", icon='OUTLINER_OB_CURVE')
        box.prop(prop, 'user_defined_resolution')
        box.prop(prop, 'x_offset')
        box.prop(prop, "closed")
        box.prop(prop, 'angle_limit')
        prop.boundary.draw(layout, context)
        prop.axis.draw(layout, context)
       
    
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
            o.location = bpy.context.scene.cursor_location
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
        if context.space_data.type == 'VIEW_3D':
            d = archipack_roof.datablock(context.active_object)
            d.manipulable_invoke(context)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to load / save presets
# ------------------------------------------------------------------


class ARCHIPACK_MT_roof_preset(Menu):
    bl_label = "Roof Styles"
    preset_subdir = "archipack_roof"
    preset_operator = "script.execute_preset"
    draw = Menu.draw_preset


class ARCHIPACK_OT_roof_preset(ArchipackPreset, Operator):
    """Add a Roof Styles"""
    bl_idname = "archipack.roof_preset"
    bl_label = "Add Roof Style"
    preset_menu = "ARCHIPACK_MT_roof_preset"

    datablock_name = StringProperty(
        name="Datablock",
        default='archipack_roof',
        maxlen=64,
        options={'HIDDEN', 'SKIP_SAVE'},
        )

    @property
    def blacklist(self):
        return ['n_parts', 'parts', 'manipulators', 'user_defined_path']


def register():
    bpy.utils.register_class(archipack_roof_material)
    bpy.utils.register_class(archipack_roof_boundary_segment)
    bpy.utils.register_class(archipack_roof_boundary)
    bpy.utils.register_class(archipack_roof_axis_segment)
    bpy.utils.register_class(archipack_roof_axis)
    bpy.utils.register_class(archipack_roof)
    Mesh.archipack_roof = CollectionProperty(type=archipack_roof)
    bpy.utils.register_class(ARCHIPACK_MT_roof_preset)
    bpy.utils.register_class(ARCHIPACK_PT_roof)
    bpy.utils.register_class(ARCHIPACK_OT_roof)
    bpy.utils.register_class(ARCHIPACK_OT_roof_preset)
    bpy.utils.register_class(ARCHIPACK_OT_roof_manipulate)
    bpy.utils.register_class(ARCHIPACK_OT_roof_from_curve)


def unregister():
    bpy.utils.unregister_class(archipack_roof_material)
    bpy.utils.unregister_class(archipack_roof_boundary_segment)
    bpy.utils.unregister_class(archipack_roof_boundary)
    bpy.utils.unregister_class(archipack_roof_axis_segment)
    bpy.utils.unregister_class(archipack_roof_axis)
    bpy.utils.unregister_class(archipack_roof)
    del Mesh.archipack_roof
    bpy.utils.unregister_class(ARCHIPACK_MT_roof_preset)
    bpy.utils.unregister_class(ARCHIPACK_PT_roof)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_manipulate)
    bpy.utils.unregister_class(ARCHIPACK_OT_roof_from_curve)
