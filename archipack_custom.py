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
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

# ----------------------------------------------------------
# Author: Stephen Leger (s-leger)
# ----------------------------------------------------------
# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty, FloatVectorProperty, BoolProperty,
    CollectionProperty, StringProperty
)
from bpy_extras import view3d_utils
from math import pi
from mathutils import Vector, Matrix
from mathutils.geometry import (
    intersect_line_plane, 
    intersect_point_line, 
    intersect_line_sphere
    )
from .archipack_manipulator import Manipulable
# from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_gl import FeedbackPanel, GlPolygon, SquareHandle, GlLine, GlText
from .archipack_object import ArchipackObject, ArchipackDrawTool
from .archipack_keymaps import Keymaps
from .archipack_dimension import DimensionProvider
xAxis = Vector((1, 0, 0))
yAxis = Vector((0, 1, 0))
zAxis = Vector((0, 0, 1))


def update(self, context):
    self.update(context)


class archipack_custom_part(ArchipackObject, PropertyGroup):
    """
      Defines and custom part
      Provide access to custom parent ui
      While parent resize from axis
      parts might have pivot not matching parent one
      so we must recompute pivot location
    """

    # pivot location relative to parent
    pivot_location = FloatVectorProperty(subtype='XYZ')

    def find_custom(self, context):
        o = self.find_in_selection(context)
        if o and o.parent:
            return o.parent


class archipack_custom(ArchipackObject, Manipulable, DimensionProvider, PropertyGroup):
    """
      This is the custom root
      use parent-child relationship to identify parts
    """
    last_size = FloatVectorProperty(subtype='XYZ')
    x = FloatProperty(
        name="width",
        min=0.01,
        update=update
        )
    y = FloatProperty(
        name="depth",
        min=0.01,
        update=update
        )
    z = FloatProperty(
        name="height",
        min=0.01,
        update=update
        )
    altitude = FloatProperty(
        name="altitude",
        update=update
        )
    flip = BoolProperty(
        description="Flag store object orientation in wall",
        default=False
        )   
    auto_update = BoolProperty(
        default=True,
        update=update,
        options={'SKIP_SAVE'}
        )

    def find_parts(self, o):
        return [
            c for c in o.children
            if archipack_custom_part.filter(c)
            ]

    def relocate(self, itM, o, dx, dy, alt):
        """
          apply delta to child location
        """
        d = archipack_custom_part.datablock(o)
        
        # child location relative to parent
        loc = d.pivot_location.copy()

        # delta location of pivot
        pivot = Vector()

        if loc.x > 0:
            pivot.x = dx
        else:
            pivot.x = -dx
        if loc.y > 0:
            pivot.y = dy
        else:
            pivot.y = -dy
        
        # delta for parent part
        loc += pivot

        # delta for child part
        # move verts so pivot remains constant in child
        cM = Matrix.Translation(-pivot)
        for v in o.data.vertices:
            v.co = cM * v.co
        # Move child so the pivot stay in same location relative to parent
        o.location = loc + Vector((0, 0, alt))
        # store location to child data
        d.pivot_location = loc.copy()

    def get_weights(self, o, group_index):
        for i, v in enumerate(o.data.vertices):
            for g in v.groups:
                if g.group == group_index:
                    yield (i, g.weight)
                    break

    def move_left(self, o, group_index):
        for i, v in enumerate(o.data.vertices):
            is_ingroup = False
            for g in v.groups:
                if g.group == group_index:
                    is_ingroup = True
                    break
            if is_ingroup:
                yield v.co.x <= 0
            else:
                yield v.co.x > 0

    def move_right(self, o, group_index):
        for i, v in enumerate(o.data.vertices):
            is_ingroup = False
            for g in v.groups:
                if g.group == group_index:
                    is_ingroup = True
                    break
            if is_ingroup:
                yield v.co.x >= 0
            else:
                yield v.co.x < 0

    def resize(self, context, o, group_name, delta, use_weights=True):
        # recompute locations in object coordsys

        m = o.data
        vgroups = {vgroup.name: vgroup.index for vgroup in o.vertex_groups}
        if group_name in vgroups:
            group_index = vgroups[group_name]
            weights = self.get_weights(o, group_index)
            if use_weights:
                for i, w in weights:
                    m.vertices[i].co += delta * w
            else:
                for i, w in weights:
                    m.vertices[i].co += delta

    def setup_manipulators(self):
        if len(self.manipulators) > 0:
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

    def manipulable_invoke(self, context):

        if self.manipulate_mode:
            self.manipulable_disable(context)
            return False

        # store current position before manipulating
        global last_pos

        # last_pos = context.active_object.matrix_world.translation.copy()

        self.manipulable_setup(context)
        self.manipulate_mode = True

        self._manipulable_invoke(context)

        return True

    def _synch_childs(self, o):
        childs = self.find_parts(o)
        # update linked child location relative to parent
        # vertices are updated through linked data
        for c in childs:
            d = archipack_custom_part.datablock(c)
            c.location = d.pivot_location.copy() + Vector((0, 0, self.altitude))

    def synch_childs(self, context, o):
        """
            synch childs nodes of linked objects
            update location relative to parent
        """
        bpy.ops.object.select_all(action='DESELECT')
        o.select = True
        context.scene.objects.active = o
        bpy.ops.object.select_linked(type='OBDATA')
        for linked in context.selected_objects:
            if linked != o:
                self._synch_childs(linked)

    def update(self, context):
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        self.setup_manipulators()

        parts = self.find_parts(o)

        last_x, last_y, last_z = self.last_size
        dx = 0.5 * (self.x - last_x)
        dy = 0.5 * (self.y - last_y)
        dz = self.z - last_z
        da = self.altitude - o.bound_box[0][2]
                
        new_size = Vector((self.x, self.y, self.z))
        parts.append(o)

        itM = o.matrix_world.inverted()
        
        if da != 0:
            for v in o.data.vertices:
                v.co.z += da
        
        for c in parts:
            if self.x != last_x:
                self.resize(context, c, "right", Vector((dx, 0, 0)))
                self.resize(context, c, "left", Vector((-dx, 0, 0)))
            elif self.y != last_y:
                self.resize(context, c, "front", Vector((0, dy, 0)))
                self.resize(context, c, "back", Vector((0, -dy, 0)))
            elif self.z != last_z:
                self.resize(context, c, "top", Vector((0, 0, dz)))
            
            # relocate child
            if o.name != c.name and (dx != 0 or dy != 0 or da != 0):
                self.relocate(itM, c, dx, dy, self.altitude)
            
        self.last_size = new_size

        x, y = 0.5 * self.x, 0.5 * self.y
        self.manipulators[0].set_pts([(-x, -y, 0), (x, -y, 0), (1, 0, 0)])
        self.manipulators[1].set_pts([(-x, -y, 0), (-x, y, 0), (-1, 0, 0)])
        self.manipulators[2].set_pts([(x, -y, self.altitude), (x, -y, self.altitude + self.z), (-1, 0, 0)])
        self.manipulators[3].set_pts([(x, -y, 0), (x, -y, self.altitude), (-1, 0, 0)])

        self.add_dimension_point(0, Vector((-x, -y, 0)))
        self.add_dimension_point(1, Vector((x, -y, 0)))
        self.update_dimensions(context, o)

        # synch linked childs location
        self.synch_childs(context, o)

        self.restore_context(context)


class ARCHIPACK_PT_custom(Panel):
    bl_idname = "ARCHIPACK_PT_custom"
    bl_label = "Custom object"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_custom.filter(context.active_object)

    def draw(self, context):
        layout = self.layout
        o = context.active_object
        d = archipack_custom.datablock(o)
        if d is None:
            return
        layout.operator('archipack.manipulate', icon='HAND')
        layout.prop(d, 'x')
        layout.prop(d, 'y')
        layout.prop(d, 'z')
        layout.prop(d, 'altitude')


class ARCHIPACK_PT_custom_part(Panel):
    bl_idname = "ARCHIPACK_PT_custom_part"
    bl_label = "Custom part"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_custom_part.filter(context.active_object)

    def draw(self, context):
        layout = self.layout
        layout.operator("archipack.select_parent")


class ARCHIPACK_OT_make_custom(Operator):
    bl_idname = "archipack.make_custom"
    bl_label = "Make Custom object"
    bl_description = "Add custom parametric ability to selection"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    x_min = FloatProperty(
        name="min"
        )
    x_max = FloatProperty(
        name="max"
        )
    x_soft = FloatProperty(
        name="soft",
        min=0
        )
    y_min = FloatProperty(
        name="min"
        )
    y_max = FloatProperty(
        name="max"
        )
    y_soft = FloatProperty(
        name="soft",
        min=0
        )
    z_max = FloatProperty(
        name="max",
        min=0
        )
    z_soft = FloatProperty(
        name="soft",
        min=0
        )

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def draw(self, context):
        layout = self.layout
        row = layout.row(align=True)
        row.label(text="min")
        row.label(text="max")
        row.label(text="soft")
        row = layout.row(align=True)
        row.prop(self, "x_min", text="x")
        row.prop(self, "x_max", text="")
        row.prop(self, "x_soft", text="")
        row = layout.row(align=True)
        row.prop(self, "y_min", text="y")
        row.prop(self, "y_max", text="")
        row.prop(self, "y_soft", text="")
        row = layout.row(align=True)
        row.label(text="z")
        row.prop(self, "z_max", text="")
        row.prop(self, "z_soft", text="")

    def weight(self, o, group, itM, soft, limit):
        """                 p.x
                    0          limit    p.x max
                          | soft |
                        start
        """
        for v in o.data.vertices:
            p = itM * v.co
            if soft > 0:
                start = max(0, limit - soft)
                w = min(1, max(0, (p.x - start) / soft))
            else:
                if p.x >= limit:
                    w = 1
                else:
                    w = 0
            if w > 0:
                group.add([v.index], w, 'REPLACE')
            else:
                group.remove([v.index])

    def add_vertex_groups(self, o, d, center):
        vgroups = {vgroup.name: vgroup.index for vgroup in o.vertex_groups}

        # matrix to retrieve vertex distance from center
        tM = Matrix.Translation(center)
        rM_top = Matrix.Rotation(-pi / 2, 4, Vector((0, 1, 0)))
        rM_front = Matrix.Rotation(pi / 2, 4, Vector((0, 0, 1)))
        rM_back = Matrix.Rotation(-pi / 2, 4, Vector((0, 0, 1)))
        rM_left = Matrix.Rotation(pi, 4, Vector((0, 0, 1)))
        rM_right = Matrix.Rotation(0, 4, Vector((0, 0, 1)))
        itM = [(tM * rM).inverted() * o.matrix_world for rM in [rM_left, rM_right, rM_back, rM_front, rM_top]]

        limit = [
            self.x_min, self.x_max,
            self.y_min, self.y_max,
            self.z_max,
            ]
        soft = [
            self.x_soft, self.x_soft,
            self.y_soft, self.y_soft,
            self.z_soft,
            ]

        for i, group_name in enumerate(['left', 'right', 'back', 'front', 'top']):
            if group_name not in vgroups:
                g = o.vertex_groups.new()
                g.name = group_name
            else:
                g = o.vertex_groups[group_name]
            self.weight(o, g, itM[i], soft[i], limit[i])
    
    def clear_parent_inverse(self, o):
        """
          Set matrix_parent_inverse to identity
          keeping visual transforms
          so .location of child is in parent coordsys
        """
        tM = o.parent.matrix_world
        itM = tM.inverted() * o.matrix_world
        o.matrix_world = tM.copy()
        o.matrix_parent_inverse = Matrix()
        o.matrix_local = itM
    
    def create(self, context):

        act = context.active_object
        d = archipack_custom.datablock(act)
        if not d:
            d = act.data.archipack_custom.add()
        if archipack_custom_part.filter(act):
            act.data.archipack_custom_part.remove(0)

        # 2d bound, bottom of geometry
        bound = Vector(act.bound_box[0]) + Vector(act.bound_box[7])
        center = act.matrix_world * (0.5 * bound)
        altitude = act.bound_box[0][2]
        itM = (Matrix.Translation(Vector((0, 0, altitude))) * act.matrix_world).inverted()
            
        d.auto_update = False
        d.last_size = act.dimensions.copy()
        d.x = act.dimensions.x
        d.y = act.dimensions.y
        d.z = act.dimensions.z
        d.altitude = altitude
        
        sel = [o for o in context.selected_objects if o.type == 'MESH']
        names = [o.name for o in sel]
        sel.extend([c for c in act.children if c.name not in names and c.type == 'MESH'])
        
        for o in sel:
                
            self.add_vertex_groups(o, d, center)
            
            if o.name != act.name:
                if archipack_custom.filter(o):
                    o.data.archipack_custom.remove(0)
                if not archipack_custom_part.filter(o):
                    o.data.archipack_custom_part.add()
                dc = archipack_custom_part.datablock(o)
                dc.pivot_location = itM * o.matrix_world.translation
                self.clear_parent_inverse(o)
            
        d.auto_update = True
        return act
       
    def check_hover(self):
        self.min.check_hover(self.mouse_pos)
        self.max.check_hover(self.mouse_pos)
        self.z.check_hover(self.mouse_pos)
    
    def mouse_release(self, context, event):
        self.active = False
        self.check_hover()
        self.min.active = False
        self.max.active = False
        self.z.active = False
        self.create(context)
        return False

    def mouse_move(self, context, event):
        self.mouse_pos = Vector((event.mouse_region_x, event.mouse_region_y))
        if self.active:
            # draggin mouse
            if self.min.active:
                pos_3d = self.get_pos3d(context)
                dp = pos_3d - self.origin
                
                self.x_min = -min(self.x_max, dp.x)
                self.y_min = -min(self.y_max, dp.y)
                
            if self.max.active:
                pos_3d = self.get_pos3d(context)
                dp = pos_3d - self.origin
                self.x_max = max(-self.x_min, dp.x)
                self.y_max = max(-self.y_min, dp.y)
                
            if self.z.active:
                pos_3d = self.get_pos3d(context, normal=Vector((1, 0, 0)))
                dp = pos_3d - self.origin
                self.z_max = dp.z
        else:
            self.check_hover()
    
    def mouse_press(self, context, event):
        if self.min.hover:
            self.active = True
            self.min.active = True
        elif self.max.hover:
            self.active = True
            self.max.active = True
        elif self.z.hover:
            self.active = True
            self.z.active = True  
                
    def get_pos3d(self, context, normal=Vector((0, 0, 1))):
        """
            convert mouse pos to 3d point over plane defined by origin and normal
            pt is in world space
        """
        region = context.region
        rv3d = context.region_data
        view_vector_mouse = view3d_utils.region_2d_to_vector_3d(region, rv3d, self.mouse_pos)
        ray_origin_mouse = view3d_utils.region_2d_to_origin_3d(region, rv3d, self.mouse_pos)
        pt = intersect_line_plane(ray_origin_mouse, ray_origin_mouse + view_vector_mouse,
            self.origin, normal, False)
        # fix issue with parallel plane
        if pt is None:
            pt = intersect_line_plane(ray_origin_mouse, ray_origin_mouse + view_vector_mouse,
                self.origin, view_vector_mouse, False)
        return pt
    
    def modal(self, context, event):
        
        context.area.tag_redraw()
        
        if event.type == 'MOUSEMOVE':
            self.mouse_move(context, event)

        elif event.value == 'PRESS':

            if event.type == 'LEFTMOUSE':
                self.mouse_press(context, event)
                
            elif event.type in {'ESC', 'RIGHTMOUSE'}:
                if self.active:
                    self.mouse_release(context, event)
                self.create(context)
                bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')        
                return {'FINISHED'}
        
        elif event.value == 'RELEASE':

            if event.type == 'LEFTMOUSE':
                self.mouse_release(context, event)
            
        return {'PASS_THROUGH'}
    
    def draw_handler(self, _self, context):
        p0 = self.origin + Vector((-self.x_min, -self.y_min, 0))
        p2 = self.origin + Vector((self.x_max, self.y_max, 0))
        pz = self.origin + Vector((0, 0, self.z_max))
        p1 = Vector((p2.x, p0.y, p0.z))
        p3 = Vector((p0.x, p2.y, p0.z))
        self.rect.set_pos([p0, p1, p2, p3])
        self.min.set_pos(context, p0, xAxis)
        self.max.set_pos(context, p2, xAxis)
        self.z.set_pos(context, pz, zAxis, normal=yAxis)
        
        self.line_xmin.p = Vector((self.pmin.x, self.origin.y - self.y_min, self.pmin.z))
        self.line_xmax.p = Vector((self.pmin.x, self.origin.y + self.y_max, self.pmin.z))
        self.line_ymin.p = Vector((self.origin.x - self.x_min, self.pmin.y, self.pmin.z))
        self.line_ymax.p = Vector((self.origin.x + self.x_max, self.pmin.y, self.pmin.z))
        self.line_z.p = Vector((self.pmin.x, self.origin.y, self.origin.z + self.z_max))
        
        cy = 0.5 * (self.y_max - self.y_min)
        cx = 0.5 * (self.x_max - self.x_min)
        self.text_xmin.set_pos(context, None, self.origin + Vector((-self.x_min - 0.1, cy, 0)), xAxis)
        self.text_xmax.set_pos(context, None, self.origin + Vector((self.x_max + 0.1, cy, 0)), xAxis)
        self.text_ymin.set_pos(context, None, self.origin + Vector((cx, -self.y_min - 0.1, 0)), xAxis)
        self.text_ymax.set_pos(context, None, self.origin + Vector((cx, self.y_max + 0.1, 0)), xAxis)
        self.text_z.set_pos(context, None, self.origin + Vector((0, 0, self.z_max + 0.1)), xAxis, normal=yAxis)       
        
        self.line_xmin.draw(context)
        self.line_xmax.draw(context)
        self.line_ymin.draw(context)
        self.line_ymax.draw(context)
        self.line_z.draw(context)
        
        self.text_xmin.draw(context)
        self.text_xmax.draw(context)
        self.text_ymin.draw(context)
        self.text_ymax.draw(context)
        self.text_z.draw(context)
        
        self.rect.draw(context)
        self.min.draw(context)
        self.max.draw(context)
        self.z.draw(context)
        
    def invoke(self, context, event):
        o = context.active_object
        p0 = Vector(o.bound_box[0])
        p1 = Vector(o.bound_box[7])
        bound = p0 + p1
        
        self.x_min = 0.5 * o.dimensions.x
        self.x_max = 0.5 * o.dimensions.x
        self.y_min = 0.5 * o.dimensions.y
        self.y_max = 0.5 * o.dimensions.y
        self.z_max = 0.5 * o.dimensions.z
        
        sensor_size = 10
        size = 0.05
        self.active = False
        self.pmin = o.matrix_world * p0
        self.origin = o.matrix_world * (0.5 * bound)
        self.mouse_pos = Vector((0, 0))
        self.rect = GlPolygon(d=3, colour=(1, 1, 1, 0.2))
        self.min = SquareHandle(sensor_size, size, draggable=True)
        self.max = SquareHandle(sensor_size, size, draggable=True)
        self.z = SquareHandle(sensor_size, size, draggable=True)
        
        p0 = o.matrix_world * p0
        p2 = o.matrix_world * p1
        p1 = Vector((p2.x, p0.y, p0.z))
        p3 = Vector((p0.x, p2.y, p0.z))
        pz = self.origin + Vector((0, 0, self.z_max))
        pz0 = self.origin + Vector((-self.x_min, 0, self.z_max))
        pz1 = self.origin + Vector((self.x_max, 0, self.z_max))
        
        self.line_xmin = GlLine(p0=p0, p1=p1)
        self.line_xmax = GlLine(p0=p3, p1=p2)
        self.line_ymin = GlLine(p0=p0, p1=p3)
        self.line_ymax = GlLine(p0=p1, p1=p2)
        self.line_z = GlLine(p0=pz0, p1=pz1)
        
        self.line_xmin.colour_inactive = (1, 0, 0, 1)
        self.line_xmax.colour_inactive = (1, 0, 0, 1)
        self.line_ymin.colour_inactive = (0, 1, 0, 1)
        self.line_ymax.colour_inactive = (0, 1, 0, 1)
        self.line_z.colour_inactive = (0, 0, 1, 1)
        
        self.text_xmin = GlText(label="Left", colour=(0, 1, 0, 1), font_size=16)
        self.text_xmax = GlText(label="Right", colour=(0, 1, 0, 1), font_size=16)
        self.text_ymin = GlText(label="Front", colour=(1, 0, 0, 1), font_size=16)
        self.text_ymax = GlText(label="Back", colour=(1, 0, 0, 1), font_size=16)
        self.text_z = GlText(label="Top", colour=(0, 0, 1, 1), font_size=16, z_axis=yAxis)
        
        cy = 0.5 * (self.y_max - self.y_min)
        cx = 0.5 * (self.x_max - self.x_min)
        self.text_xmin.set_pos(context, None, self.origin + Vector((-self.x_min - 0.1, cy, 0)), xAxis)
        self.text_xmax.set_pos(context, None, self.origin + Vector((self.x_max + 0.1, cy, 0)), xAxis)
        self.text_ymin.set_pos(context, None, self.origin + Vector((cx, -self.y_min - 0.1, 0)), xAxis)
        self.text_ymax.set_pos(context, None, self.origin + Vector((cx, self.y_max + 0.1, 0)), xAxis)
        self.text_z.set_pos(context, None, self.origin + Vector((0, 0, self.z_max + 0.1)), xAxis, normal=yAxis)
        
        self.min.set_pos(context, p0, xAxis)
        self.max.set_pos(context, p1, xAxis)
        self.z.set_pos(context, pz, zAxis, normal=yAxis)
        self.rect.set_pos([p0, p1, p2, p3])
        
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_handler,
                args, 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
        
    def execute(self, context):
        if context.mode == "OBJECT":
            o = self.create(context)
            o.select = True
            context.scene.objects.active = o
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}
          
            
class ARCHIPACK_OT_custom_draw(ArchipackDrawTool, Operator):
    bl_idname = "archipack.custom_draw"
    bl_label = "Draw Custom"
    bl_description = "Draw Custom object over walls"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    filepath = StringProperty(default="")
    feedback = None
    stack = []
    object_name = ""

    @classmethod
    def poll(cls, context):
        o = context.active_object
        return archipack_custom.filter(o)

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def draw_callback(self, _self, context):
        self.feedback.draw(context)

    def add_object(self, context, event):

        o = context.active_object
        bpy.ops.object.select_all(action="DESELECT")

        if archipack_custom.filter(o):

            o.select = True
            context.scene.objects.active = o

            # instance subs
            o.select = False

            # copy when shift pressed
            o = self.duplicate_object(context, o, not event.shift)

            o.select = True
            context.scene.objects.active = o

        self.object_name = o.name

        o.select = True
        context.scene.objects.active = o

    def modal(self, context, event):

        context.area.tag_redraw()
        o = context.scene.objects.get(self.object_name)

        if o is None:
            return {'FINISHED'}

        d = archipack_custom.datablock(o)

        # hide hole from raycast
        to_hide = [o]
        to_hide.extend([child for child in o.children if archipack_custom_part.filter(child)])

        for obj in to_hide:
            obj.hide = True

        res, tM, wall, width, y, z_offset = self.mouse_hover_wall(context, event)

        for obj in to_hide:
            obj.hide = False

        if res and d is not None:
            o.matrix_world = tM.copy()
            if d.y != width:
                d.y = width

        if event.value == 'PRESS':

            if event.type in {'LEFTMOUSE', 'RET', 'NUMPAD_ENTER', 'SPACE'}:
                if wall is not None:

                    wall.select = True
                    context.scene.objects.active = wall
                    if bpy.ops.archipack.single_boolean.poll():
                        bpy.ops.archipack.single_boolean()

                    wall.select = False
                    # o must be a custom here
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
                context.scene.objects.active = last
                self.delete_object(context, last)
                context.scene.objects.active = o
            return {'RUNNING_MODAL'}

        if event.value == 'RELEASE':

            if event.type in {'ESC', 'RIGHTMOUSE'}:
                self.delete_object(context, o)
                self.feedback.disable()
                bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
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
            o = context.active_object
            context.scene.objects.active = None
            bpy.ops.object.select_all(action="DESELECT")
            if o is not None:
                o.select = True
                context.scene.objects.active = o
            self.add_object(context, event)
            self.feedback = FeedbackPanel()
            self.feedback.instructions(context, "Draw a custom", "Click & Drag over a wall", [
                ('LEFTCLICK, RET, SPACE, ENTER', 'Create a custom'),
                ('BACKSPACE, CTRL+Z', 'undo last'),
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


def register():
    bpy.utils.register_class(archipack_custom_part)
    bpy.utils.register_class(archipack_custom)
    Mesh.archipack_custom = CollectionProperty(type=archipack_custom)
    Mesh.archipack_custom_part = CollectionProperty(type=archipack_custom_part)
    bpy.utils.register_class(ARCHIPACK_PT_custom)
    bpy.utils.register_class(ARCHIPACK_PT_custom_part)
    bpy.utils.register_class(ARCHIPACK_OT_custom_draw)
    bpy.utils.register_class(ARCHIPACK_OT_make_custom)


def unregister():
    bpy.utils.unregister_class(archipack_custom)
    bpy.utils.unregister_class(archipack_custom_part)
    bpy.utils.unregister_class(ARCHIPACK_PT_custom)
    bpy.utils.unregister_class(ARCHIPACK_PT_custom_part)
    bpy.utils.unregister_class(ARCHIPACK_OT_custom_draw)
    bpy.utils.unregister_class(ARCHIPACK_OT_make_custom)
    del Mesh.archipack_custom
    del Mesh.archipack_custom_part
