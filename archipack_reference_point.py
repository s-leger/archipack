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
from bpy.types import Operator, PropertyGroup, Object, Panel
from bpy.props import FloatVectorProperty, CollectionProperty
from mathutils import Vector


class archipack_reference_point(PropertyGroup):
    location_2d = FloatVectorProperty(
        subtype='XYZ',
        name="position 2d",
        default=Vector((0, 0, 0))
        )
    location_3d = FloatVectorProperty(
        subtype='XYZ',
        name="position 3d",
        default=Vector((0, 0, 0))
        )

    @classmethod
    def filter(cls, o):
        """
            Filter object with this class in data
            return
            True when object contains this datablock
            False otherwhise
            usage:
            class_name.filter(object) from outside world
            self.__class__.filter(object) from instance
        """
        try:
            return cls.__name__ in o
        except:
            pass
        return False

    @classmethod
    def datablock(cls, o):
        """
            Retrieve datablock from base object
            return
                datablock when found
                None when not found
            usage:
                class_name.datablock(object) from outside world
                self.__class__.datablock(object) from instance
        """
        try:
            return getattr(o, cls.__name__)[0]
        except:
            pass
        return None


class ARCHIPACK_PT_reference_point(Panel):
    bl_idname = "ARCHIPACK_PT_reference_point"
    bl_label = "Reference point"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_reference_point.filter(context.active_object)

    def draw(self, context):
        o = context.active_object
        props = archipack_reference_point.datablock(o)
        if props is None:
            return
        layout = self.layout
        if (o.location - props.location_2d).length < 0.01:
            layout.operator('archipack.move_to_3d')
            layout.operator('archipack.move_2d_reference_to_cursor')
        else:
            layout.operator('archipack.move_to_2d')


class ARCHIPACK_OT_reference_point(Operator):
    """Add reference point"""
    bl_idname = "archipack.reference_point"
    bl_label = "Reference point"
    bl_description = "Add reference point"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    location_3d = FloatVectorProperty(
        subtype='XYZ',
        name="position 3d",
        default=Vector((0, 0, 0))
        )

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def execute(self, context):
        if context.mode == "OBJECT":
            x, y, z = context.scene.cursor_location
            bpy.ops.object.empty_add(type='ARROWS', radius=0.5, location=Vector((x, y, 0)))
            o = context.active_object
            o.name = "Reference"
            props = o.archipack_reference_point.add()
            props.location_2d = Vector((x, y, 0))
            props.location_3d = self.location_3d
            o.select = True
            context.scene.objects.active = o
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_move_to_3d(Operator):
    bl_idname = "archipack.move_to_3d"
    bl_label = "Move to 3d"
    bl_description = "Move point to 3d position"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return archipack_reference_point.filter(context.active_object)

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            props = archipack_reference_point.datablock(o)
            if props is None:
                return {'CANCELLED'}
            o.location = props.location_3d
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_move_to_2d(Operator):
    bl_idname = "archipack.move_to_2d"
    bl_label = "Move to 2d"
    bl_description = "Move point to 2d position"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return archipack_reference_point.filter(context.active_object)

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            props = archipack_reference_point.datablock(o)
            if props is None:
                return {'CANCELLED'}
            props.location_3d = o.location
            o.location = props.location_2d
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_store_2d_reference(Operator):
    bl_idname = "archipack.store_2d_reference"
    bl_label = "Set 2d"
    bl_description = "Set 2d reference position"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return archipack_reference_point.filter(context.active_object)

def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            props = archipack_reference_point.datablock(o)
            if props is None:
                return {'CANCELLED'}
            x, y, z = o.location
            props.location_2d = Vector((x, y, 0))
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_move_2d_reference_to_cursor(Operator):
    bl_idname = "archipack.move_2d_reference_to_cursor"
    bl_label = "Change 2d"
    bl_description = "Change 2d reference position to cursor location without moving childs"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return archipack_reference_point.filter(context.active_object)

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            props = archipack_reference_point.datablock(o)
            if props is None:
                return {'CANCELLED'}
            bpy.ops.object.select_all(action="DESELECT")
            bpy.ops.archipack.reference_point(location_3d=props.location_3d)
            for child in o.children:
                child.select = True
            bpy.ops.archipack.parent_to_reference()
            context.scene.objects.unlink(o)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_parent_to_reference(Operator):
    bl_idname = "archipack.parent_to_reference"
    bl_label = "Parent"
    bl_description = "Make selected object childs of parent reference point"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        return archipack_reference_point.filter(context.active_object)

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            props = archipack_reference_point.datablock(o)
            if props is None:
                return {'CANCELLED'}
            sel = [obj for obj in context.selected_objects if obj != o and obj.parent != o]
            itM = o.matrix_world.inverted()
            # print("parent_to_reference parenting:%s objects" % (len(sel)))
            for child in sel:
                rs = child.matrix_world.to_3x3().to_4x4()
                loc = itM * child.matrix_world.translation
                child.parent = None
                child.matrix_parent_inverse.identity()
                child.location = Vector((0, 0, 0))
                child.parent = o
                child.matrix_world = rs
                child.location = loc
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


def register():
    bpy.utils.register_class(archipack_reference_point)
    Object.archipack_reference_point = CollectionProperty(type=archipack_reference_point)
    bpy.utils.register_class(ARCHIPACK_PT_reference_point)
    bpy.utils.register_class(ARCHIPACK_OT_reference_point)
    bpy.utils.register_class(ARCHIPACK_OT_move_to_3d)
    bpy.utils.register_class(ARCHIPACK_OT_move_to_2d)
    bpy.utils.register_class(ARCHIPACK_OT_store_2d_reference)
    bpy.utils.register_class(ARCHIPACK_OT_move_2d_reference_to_cursor)
    bpy.utils.register_class(ARCHIPACK_OT_parent_to_reference)


def unregister():
    bpy.utils.unregister_class(archipack_reference_point)
    del Object.archipack_reference_point
    bpy.utils.unregister_class(ARCHIPACK_PT_reference_point)
    bpy.utils.unregister_class(ARCHIPACK_OT_reference_point)
    bpy.utils.unregister_class(ARCHIPACK_OT_move_to_3d)
    bpy.utils.unregister_class(ARCHIPACK_OT_move_to_2d)
    bpy.utils.unregister_class(ARCHIPACK_OT_store_2d_reference)
    bpy.utils.unregister_class(ARCHIPACK_OT_move_2d_reference_to_cursor)
    bpy.utils.unregister_class(ARCHIPACK_OT_parent_to_reference)
