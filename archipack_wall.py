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
import bmesh
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import FloatProperty, CollectionProperty
from .archipack_object import ArchipackObject, ArchipackCreateTool


def update_wall(self, context):
    self.update(context)


class archipack_wall(ArchipackObject, PropertyGroup):
    z = FloatProperty(
            name='height',
            min=0.1, max=10000,
            default=2.7, precision=2,
            description='height', update=update_wall,
            )

    def update(self, context):
        # update height via bmesh to avoid loosing material ids
        # this should be the rule for other simple objects
        # as long as there is no topologic changes
        o = context.active_object
        if archipack_wall.datablock(o) != self:
            return
        # bpy.ops.object.mode_set(mode='EDIT')
        me = o.data
        # bm = bmesh.from_edit_mesh(me)
        # bm.verts.ensure_lookup_table()
        # bm.faces.ensure_lookup_table()
        new_z = self.z
        # use "Top" vertex group
        # create vertex group lookup dictionary for names
        vgroup_names = {vgroup.index: vgroup.name for vgroup in o.vertex_groups}
        # create dictionary of vertex group assignments per vertex
        vertex_groups = [[vgroup_names[g.group] for g in v.groups] for v in me.vertices]
        last_z = [me.vertices[i].co.z for i, g in enumerate(vertex_groups) if 'Top' in g]
        if len(last_z) > 0:
            max_z = max(last_z)
            for i, g in enumerate(vertex_groups):
                if 'Top' in g:
                    me.vertices[i].co.z = new_z
        # bpy.ops.object.mode_set(mode='OBJECT')


class ARCHIPACK_PT_wall(Panel):
    bl_idname = "ARCHIPACK_PT_wall"
    bl_label = "Wall"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_wall.filter(context.active_object)

    def draw(self, context):

        prop = archipack_wall.datablock(context.active_object)
        if prop is None:
            return        
        layout = self.layout
        layout.prop(prop, 'z')


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_wall(Operator, ArchipackCreateTool):
    bl_idname = "archipack.wall"
    bl_label = "Wall"
    bl_description = "Add wall parameters to active object to draw windows and doors and use autoboolean"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    z = FloatProperty(
        name="z",
        default=2.7
        )

    @classmethod
    def poll(cls, context):
        o = context.active_object
        return o is not None and o.type == 'MESH' and "archipack_wall2" not in o.data

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            if archipack_wall.filter(o):
                return {'CANCELLED'}
            params = o.data.archipack_wall.add()
            # @TODO: estimate z from mesh,
            # using "TOP" vertex group
            params.z = self.z
            self.add_material(o)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


def register():
    bpy.utils.register_class(archipack_wall)
    Mesh.archipack_wall = CollectionProperty(type=archipack_wall)
    bpy.utils.register_class(ARCHIPACK_PT_wall)
    bpy.utils.register_class(ARCHIPACK_OT_wall)


def unregister():
    bpy.utils.unregister_class(archipack_wall)
    del Mesh.archipack_wall
    bpy.utils.unregister_class(ARCHIPACK_PT_wall)
    bpy.utils.unregister_class(ARCHIPACK_OT_wall)
