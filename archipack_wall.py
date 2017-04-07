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

# ----------------------------------------------------------
# Author: Stephen Leger (s-leger)
#
# ----------------------------------------------------------
# noinspection PyUnresolvedReferences
import bpy
import bmesh
# noinspection PyUnresolvedReferences
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import FloatProperty, CollectionProperty, BoolProperty 
from mathutils import Vector, Matrix

def update_wall(self, context):
    self.update(context)
   
class WallProperty(PropertyGroup):
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
        if ARCHIPACK_PT_wall.params(o) != self:
            return 
        bpy.ops.object.mode_set(mode='EDIT')    
        me = o.data
        bm = bmesh.from_edit_mesh(me)
        bm.verts.ensure_lookup_table()
        bm.faces.ensure_lookup_table()
        new_z = self.z
        last_z = list(v.co.z for v in bm.verts)
        max_z = max(last_z)
        for v in bm.verts:
            if v.co.z == max_z:
                v.co.z = new_z
        bmesh.update_edit_mesh(me, True)
        bpy.ops.object.mode_set(mode='OBJECT')
        
class ARCHIPACK_PT_wall(Panel):
    bl_idname = "ARCHIPACK_PT_wall"
    bl_label = "Wall"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'
    
    def draw(self, context):
        o = context.object
        if not ARCHIPACK_PT_wall.filter(o):
            return 
        layout = self.layout
        prop = o.data.WallProperty[0]
        layout.prop(prop, 'z')
     
    @classmethod
    def params(cls, o):
        try:
            if 'WallProperty' not in o.data:
                return False
            else:
                return o.data.WallProperty[0]
        except:
            return False
    
    @classmethod
    def filter(cls, o):
        try:
            if 'WallProperty' not in o.data:
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
class ARCHIPACK_OT_wall(Operator):
    bl_idname = "archipack.wall"
    bl_label = "Wall"
    bl_description = "Add wall parameters to active object"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    z = FloatProperty(
        name="z",
        default=2.7
        )
        
    @classmethod
    def poll(cls, context):
        return context.active_object is not None
    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')
    
    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            if ARCHIPACK_PT_wall.filter(o):
                return {'CANCELLED'}
            params = o.data.WallProperty.add()
            params.z = self.z
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}

bpy.utils.register_class(WallProperty)
Mesh.WallProperty = CollectionProperty(type=WallProperty)
bpy.utils.register_class(ARCHIPACK_PT_wall)
bpy.utils.register_class(ARCHIPACK_OT_wall)

