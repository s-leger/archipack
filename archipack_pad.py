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
import bmesh
from .bmesh_utils import BmeshEdit as bmed
from .materialutils import MaterialUtils
from .panel import Panel as Lofter
from mathutils import Vector, Matrix
from mathutils.geometry import interpolate_bezier
from math import sin, cos, pi, acos, atan2
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_2d import Line, Arc
from .archipack_preset import ArchipackPreset


class Pad():

    def __init__(self):
        pass
        
    def straight_pad(self, last, a0, length):
        s = last.straight(length).rotate(a0)
        return StraightPad(s.p, s.v)

    def curved_pad(self, last, a0, da, radius):
        n = last.normal(1).rotate(a0).scale(radius)
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedPad(c, radius, a0, da)


class StraightPad(Pad, Line):
    
    def __init__(self, p, v):
        Pad.__init__(self)
        Line.__init__(self, p, v)


class CurvedPad(Pad, Arc):
    
    def __init__(self, c, radius, a0, da):
        Pad.__init__(self)
        Arc.__init__(self, c, radius, a0, da)


class PadGenerator():

    def __init__(self, parts):
        self.parts = parts
        self.segs = []
        
    def add_part(self, type, radius, a0, da, length, offset):

        if len(self.segs) < 1:
            s = None
        else:
            s = self.segs[-1]
        
        # start a new pad
        if s is None:
            if type == 'S_SEG':
                p = Vector((0, 0))
                v = length * Vector((cos(a0), sin(a0)))
                s = StraightPad(p, v)
            elif type == 'C_SEG':
                c = -radius * Vector((cos(a0), sin(a0)))
                s = CurvedPad(c, radius, a0, da)
        else:
            if type == 'S_SEG':
                s = s.straight_pad(s, a0, length)
            elif type == 'C_SEG':
                s = s.curved_pad(s, a0, da, radius)
        
        self.segs.append(s)
        self.last_type = type
        
        
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
            
    def locate_manipulators(self):
        """
            setup manipulators
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

            if type(f).__name__ == "StraightPad":
                # segment length
                manipulators[1].type = 'SIZE'
                manipulators[1].prop1_name = "length"
                manipulators[1].set_pts([p0, p1, (1, 0, 0)])
            else:
                # segment radius + angle
                v0 = (f.p0 - f.c).to_3d()
                v1 = (f.p1 - f.c).to_3d()
                manipulators[1].type = 'ARC_ANGLE_RADIUS'
                manipulators[1].prop1_name = "da"
                manipulators[1].prop2_name = "radius"
                manipulators[1].set_pts([f.c.to_3d(), v0, v1])

            # snap manipulator, dont change index !
            manipulators[2].set_pts([p0, p1, (1, 0, 0)])
            # dumb segment id
            manipulators[3].set_pts([p0, p1, (1, 0, 0)])
    
    def get_verts(self, verts):
        for pad in self.segs:
            if "Curved" in type(pad).__name__:
                for i in range(16):
                    x, y = pad.lerp(i / 16)
                    verts.append((x, y, 0))
            else:
                x, y = pad.p0
                verts.append((x, y, 0))
            """
            for i in range(33):
                x, y = pad.line.lerp(i / 32)
                verts.append((x, y, 0))
            """


def update(self, context):
    self.update(context)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


def update_path(self, context):
    self.update_path(context)


materials_enum = (
            ('0', 'Ceiling', '', 0),
            ('1', 'White', '', 1),
            ('2', 'Concrete', '', 2),
            ('3', 'Wood', '', 3),
            ('4', 'Metal', '', 4),
            ('5', 'Glass', '', 5)
            )


class archipack_pad_material(PropertyGroup):
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
            props = ARCHIPACK_PT_pad.params(o)
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
    """
        A single manipulable segment
        either linear or arc based
    """
    type = EnumProperty(
            items=(
                ('S_SEG', 'Straight pad', '', 0),
                ('C_SEG', 'Curved pad', '', 1),
                ),
            default='S_SEG',
            update=update_manipulators
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
    
    def draw_insert(self, context, layout, index):
        """
            May implement draw for insert / remove operators
        """
        pass
        
    def draw(self, context, layout, index):
        box = layout.box()
        row = box.row()
        row.prop(self, "type", text=str(index + 1))
        self.draw_insert(context, box, index)
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



class archipack_pad_part(ArchipackSegment, PropertyGroup):
    
    def draw_insert(self, context, layout, index):
        row = layout.row(align=True)
        row.operator("archipack.pad_insert", text="Split").index = index
        row.operator("archipack.pad_remove", text="Remove").index = index
            
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = ARCHIPACK_PT_pad.params(o)
            if props:
                for part in props.parts:
                    if part == self:
                        return props
        return None


        
class archipack_pad(Manipulable, PropertyGroup):
    # boundary
    n_parts = IntProperty(
        name="parts",
            min=1,
            default=1, update=update_manipulators
            )
    parts = CollectionProperty(type=archipack_pad_part)
    closed = BoolProperty(
            default=False,
            name="Close",
            update=update_manipulators
            )
    # UI layout related
    parts_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
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
    
    z = FloatProperty(
            name="z",
            default=0.3, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
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

    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        active = context.active_object
        selected = [o for o in context.selected_objects]
        for o in selected:
            if ARCHIPACK_PT_pad.params(o) == self:
                return active, selected, o
        return active, selected, None

    def get_generator(self):
        g = PadGenerator(self.parts)
        for part in self.parts:
            # type, radius, da, length
            g.add_part(part.type, part.radius, part.a0, part.da, part.length, 0)

        g.close(self.closed)    
        g.locate_manipulators()
        return g    
    
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
        self.closed = True
        self.auto_update = True

    def update_path(self, context):
        user_def_path = context.scene.objects.get(self.user_defined_path)
        if user_def_path is not None and user_def_path.type == 'CURVE':
            self.from_spline(user_def_path.matrix_world, self.user_defined_resolution, user_def_path.data.splines[0])

    def make_surface(self, o, verts):
        bm = bmesh.new()
        for v in verts:
            bm.verts.new(v)
        bm.verts.ensure_lookup_table()
        for i in range(1, len(verts)):
            bm.edges.new((bm.verts[i - 1], bm.verts[i]))
        bm.edges.new((bm.verts[-1], bm.verts[0]))
        bm.edges.ensure_lookup_table()
        bmesh.ops.contextual_create(bm, geom=bm.edges)
        bm.to_mesh(o.data)
        bm.free()
        
        """
        import bmesh
        verts = [
            Vector((0,0,0)),
            Vector((1,0,0)),
            Vector((1,1,0)),
            Vector((0,1,0))
        ]
        bm = bmesh.new()
        for v in verts:
            bm.verts.new(v)
        bm.verts.ensure_lookup_table()
        for i in range(1, len(verts)):
            bm.edges.new((bm.verts[i - 1], bm.verts[i]))
        bm.edges.ensure_lookup_table()
        bmesh.ops.contextual_create(bm, geom=bm.edges)
        bm.to_mesh(o.data)
        bm.free()
        """
    
    def unwrap_uv(self, o):
        bm = bmesh.new()
        bm.from_mesh(o.data)
        for face in bm.faces:
            face.select = face.material_index > 0
        bm.to_mesh(o.data)
        bpy.ops.uv.cube_project(scale_to_bounds=False, correct_aspect=True)
        
        for face in bm.faces:
            face.select = face.material_index < 1
        bm.to_mesh(o.data)
        bpy.ops.uv.smart_project(use_aspect=True, stretch_to_bounds=False)
        bm.free()
             
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

        g = self.get_generator()
        g.get_verts(verts)
        self.make_surface(o, verts)
        
        modif = o.modifiers.get('Pad')
        if modif is None:
            modif = o.modifiers.new('Pad', 'SOLIDIFY')
            modif.use_quality_normals = True
            modif.use_even_offset = True
            modif.material_offset_rim = 2
            modif.material_offset = 1

        modif.thickness = self.z
        modif.offset = 1.0

        # Height
        self.manipulators[0].set_pts([
            (0, 0, 0),
            (0, 0, -self.z),
            (-1, 0, 0)
            ], normal=g.segs[0].straight(-1, 0).v.to_3d())

        # bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=uvs, weld=True, clean=False)

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

        for i, part in enumerate(d.parts):
            if i >= d.n_parts:
                break

            if i > 0:
                # start angle
                self.manip_stack.append(part.manipulators[0].setup(context, o, part))

            # length / radius + angle
            self.manip_stack.append(part.manipulators[1].setup(context, o, part))

            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, d))
            # index
            self.manip_stack.append(part.manipulators[3].setup(context, o, d))
            
       
        for m in self.manipulators:
            self.manip_stack.append(m.setup(context, o, self))


class ARCHIPACK_PT_pad(Panel):
    bl_idname = "ARCHIPACK_PT_pad"
    bl_label = "Pad"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    def draw(self, context):
        o = context.object
        scene = context.scene
        prop = ARCHIPACK_PT_pad.params(o)
        if prop is None:
            return
        layout = self.layout
        row = layout.row(align=True)
        row.operator('archipack.pad_manipulate')
        box = layout.box()
        # box.label(text="Styles")
        row = box.row(align=True)
        row.menu("ARCHIPACK_MT_pad_preset", text=bpy.types.ARCHIPACK_MT_pad_preset.bl_label)
        row.operator("archipack.pad_preset", text="", icon='ZOOMIN')
        row.operator("archipack.pad_preset", text="", icon='ZOOMOUT').remove_active = True
        row = layout.row(align=True)
        box = layout.box()
        row.prop_search(prop, "user_defined_path", scene, "objects", text="", icon='OUTLINER_OB_CURVE')
        box.prop(prop, 'user_defined_resolution')
        box.prop(prop, 'angle_limit')
        box.prop(prop, 'z')
        box = layout.box()
        row = box.row()
        if prop.parts_expand:
            row.prop(prop, 'parts_expand', icon="TRIA_DOWN", icon_only=True, text="Parts", emboss=False)
            box.prop(prop, 'n_parts')
            # box.prop(prop, 'closed')
            for i, part in enumerate(prop.parts):
                part.draw(context, layout, i)
        else:
            row.prop(prop, 'parts_expand', icon="TRIA_RIGHT", icon_only=True, text="Parts", emboss=False)
     
        
        
    @classmethod
    def params(cls, o):
        try:
            if 'archipack_pad' not in o.data:
                return False
            else:
                return o.data.archipack_pad[0]
        except:
            return False

    @classmethod
    def filter(cls, o):
        try:
            if 'archipack_pad' not in o.data:
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


class ARCHIPACK_OT_pad(Operator):
    bl_idname = "archipack.pad"
    bl_label = "Pad"
    bl_description = "Pad"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    auto_manipulate = BoolProperty(default=True)

    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        m = bpy.data.meshes.new("Pad")
        o = bpy.data.objects.new("Pad", m)
        d = m.archipack_pad.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        s = d.manipulators.add()
        s.prop1_name = "z"
        s.normal = Vector((0, 1, 0))
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        d.update(context)
        # MaterialUtils.add_pad_materials(o)
        o.location = bpy.context.scene.cursor_location
        o.select = True
        context.scene.objects.active = o
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = bpy.context.scene.cursor_location
            if self.auto_manipulate:
                bpy.ops.archipack.pad_manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}

           
class ARCHIPACK_OT_pad_insert(Operator):
    bl_idname = "archipack.pad_insert"
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
            d = ARCHIPACK_PT_pad.params(o)
            if d is None:
                return {'CANCELLED'}
            d.insert_part(context, self.index)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_pad_remove(Operator):
    bl_idname = "archipack.pad_remove"
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
            d = ARCHIPACK_PT_pad.params(o)
            if d is None:
                return {'CANCELLED'}
            d.remove_part(context, self.index)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}

# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_pad_from_curve(Operator):
    bl_idname = "archipack.pad_from_curve"
    bl_label = "Pad curve"
    bl_description = "Create a pad from a curve"
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
        m = bpy.data.meshes.new("Pad")
        o = bpy.data.objects.new("Pad", m)
        d = m.archipack_pad.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        s = d.manipulators.add()
        s.prop1_name = "z"
        s.normal = Vector((0, 1, 0))
        d.user_defined_path = curve.name
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        d.update_path(context)
        MaterialUtils.add_wall_materials(o)
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
                bpy.ops.archipack.pad_manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}

            
class ARCHIPACK_OT_pad_from_wall(Operator):
    bl_idname = "archipack.pad_from_wall"
    bl_label = "Wall -> Pad"
    bl_description = "Create a pad from a wall"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    auto_manipulate = BoolProperty(default=True)

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o is not None and o.data is not None and 'archipack_wall2' in o.data
    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        wall = context.active_object
        m = bpy.data.meshes.new("Pad")
        o = bpy.data.objects.new("Pad", m)
        d = m.archipack_pad.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        s = d.manipulators.add()
        s.prop1_name = "z"
        s.normal = Vector((0, 1, 0))
        wd = wall.data.archipack_wall2[0]
        d.auto_update = False
        d.closed = True
        d.parts.clear()
        d.n_parts = wd.n_parts
        if not wd.closed:
            d.n_parts += 1
        for part in wd.parts:
            p = d.parts.add()
            if "S_" in part.type:
                p.type = "S_SEG"
            else:
                p.type = "C_SEG"
            p.length = part.length
            p.radius = part.radius
            p.da = part.da
            p.a0 = part.a0
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        d.auto_update = True
        MaterialUtils.add_wall_materials(o)
        
        # pretranslate
        o.matrix_world = wall.matrix_world.copy()
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
                bpy.ops.archipack.pad_manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}

# ------------------------------------------------------------------
# Define operator class to manipulate object
# ------------------------------------------------------------------


class ARCHIPACK_OT_pad_manipulate(Operator):
    bl_idname = "archipack.pad_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return ARCHIPACK_PT_pad.filter(context.active_object)

    def modal(self, context, event):
        return self.d.manipulable_modal(context, event)

    def invoke(self, context, event):
        if context.space_data.type == 'VIEW_3D':
            o = context.active_object
            self.d = o.data.archipack_pad[0]
            if self.d.manipulable_invoke(context):
                context.window_manager.modal_handler_add(self)
                return {'RUNNING_MODAL'}
            else:
                return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to load / save presets
# ------------------------------------------------------------------


class ARCHIPACK_MT_pad_preset(Menu):
    bl_label = "Pad Styles"
    preset_subdir = "archipack_pad"
    preset_operator = "script.execute_preset"
    draw = Menu.draw_preset


class ARCHIPACK_OT_pad_preset(ArchipackPreset, Operator):
    """Add a Pad Styles"""
    bl_idname = "archipack.pad_preset"
    bl_label = "Add Pad Style"
    preset_menu = "ARCHIPACK_MT_pad_preset"

    datablock_name = StringProperty(
        name="Datablock",
        default='archipack_pad',
        maxlen=64,
        options={'HIDDEN', 'SKIP_SAVE'},
        )

    @property
    def blacklist(self):
        return ['n_parts', 'parts', 'manipulators', 'user_defined_path']


def register():
    bpy.utils.register_class(archipack_pad_material)
    bpy.utils.register_class(archipack_pad_part)
    bpy.utils.register_class(archipack_pad)
    Mesh.archipack_pad = CollectionProperty(type=archipack_pad)
    bpy.utils.register_class(ARCHIPACK_MT_pad_preset)
    bpy.utils.register_class(ARCHIPACK_PT_pad)
    bpy.utils.register_class(ARCHIPACK_OT_pad)
    bpy.utils.register_class(ARCHIPACK_OT_pad_insert)
    bpy.utils.register_class(ARCHIPACK_OT_pad_remove)
    bpy.utils.register_class(ARCHIPACK_OT_pad_preset)
    bpy.utils.register_class(ARCHIPACK_OT_pad_manipulate)
    bpy.utils.register_class(ARCHIPACK_OT_pad_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_pad_from_wall)

def unregister():
    bpy.utils.unregister_class(archipack_pad_material)
    bpy.utils.unregister_class(archipack_pad_part)
    bpy.utils.unregister_class(archipack_pad)
    del Mesh.archipack_pad
    bpy.utils.unregister_class(ARCHIPACK_MT_pad_preset)
    bpy.utils.unregister_class(ARCHIPACK_PT_pad)
    bpy.utils.unregister_class(ARCHIPACK_OT_pad)
    bpy.utils.unregister_class(ARCHIPACK_OT_pad_insert)
    bpy.utils.unregister_class(ARCHIPACK_OT_pad_remove)
    bpy.utils.unregister_class(ARCHIPACK_OT_pad_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_pad_manipulate)
    bpy.utils.unregister_class(ARCHIPACK_OT_pad_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_pad_from_wall)