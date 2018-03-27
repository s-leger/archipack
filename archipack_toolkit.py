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
from bpy.types import Operator, PropertyGroup, Mesh, Panel

# Minimal required property types
from bpy.props import (
    FloatProperty, BoolProperty, CollectionProperty
    )
from mathutils import Vector
# Bmesh made easy
from .bmesh_utils import BmeshEdit as bmed

# Manipulable
from .archipack_manipulator import Manipulable

# Preset system
from .archipack_preset import ArchipackPreset, PresetMenuOperator

# Base for Propertygroup and create tool Operator
from .archipack_object import ArchipackCreateTool, ArchipackObject


def update(self, context):
    self.update(context)


class archipack_myobject(ArchipackObject, Manipulable, PropertyGroup):
    """ Archipack toolkit sample"""
    x = FloatProperty(
            name="Width",
            default=2.0, min=0.01,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    y = FloatProperty(
            name="Depth",
            default=2.0, min=0.01,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    z = FloatProperty(
            name="Height",
            default=2.0, min=0.01,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    auto_update = BoolProperty(
            # Wont save auto_update state in any case
            options={'SKIP_SAVE'},
            default=True,
            update=update
            )

    @property
    def verts(self):
        """
            Object vertices coords
        """
        x = 0.5 * self.x
        y = 0.5 * self.y
        z = self.z
        return [
            (-x, y, 0),
            (-x, -y, 0),
            (x, -y, 0),
            (x, y, 0),
            (-x, y, z),
            (-x, -y, z),
            (x, -y, z),
            (x, y, z)
        ]

    @property
    def faces(self):
        """
            Object faces vertices index
        """
        return [
            (0, 1, 2, 3),
            (7, 6, 5, 4),
            (7, 4, 0, 3),
            (4, 5, 1, 0),
            (5, 6, 2, 1),
            (6, 7, 3, 2)
        ]

    @property
    def uvs(self):
        """
            Object faces uv coords
        """
        return [
            [(0, 0), (0, 1), (1, 1), (1, 0)],
            [(0, 0), (0, 1), (1, 1), (1, 0)],
            [(0, 0), (0, 1), (1, 1), (1, 0)],
            [(0, 0), (0, 1), (1, 1), (1, 0)],
            [(0, 0), (0, 1), (1, 1), (1, 0)],
            [(0, 0), (0, 1), (1, 1), (1, 0)]
        ]

    @property
    def matids(self):
        """
            Object material indexes for each face
        """
        return [0, 0, 0, 0, 0, 0]

    def setup_manipulators(self):
        if len(self.manipulators) < 1:
            # add manipulator for x property
            s = self.manipulators.add()
            s.prop1_name = "x"
            s.type_key = 'SIZE'

            # add manipulator for y property
            s = self.manipulators.add()
            s.prop1_name = "y"
            s.type_key = 'SIZE'

            # add manipulator for z property
            s = self.manipulators.add()
            s.prop1_name = "z"
            s.type_key = 'SIZE'
            # draw this one on xz plane
            s.normal = Vector((0, 1, 0))

    def update(self, context):

        # provide support for "copy to selected"
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # dynamically create manipulators when needed
        self.setup_manipulators()

        # update your mesh from parameters
        bmed.buildmesh(context,
                       o,
                       self.verts,
                       self.faces,
                       matids=self.matids,
                       uvs=self.uvs,
                       weld=False)

        # update manipulators location (3d location in object coordsystem)
        x, y = 0.5 * self.x, 0.5 * self.y
        self.manipulators[0].set_pts([(-x, -y, 0), (x, -y, 0), (1, 0, 0)])
        self.manipulators[1].set_pts([(-x, -y, 0), (-x, y, 0), (-1, 0, 0)])
        self.manipulators[2].set_pts([(x, -y, 0), (x, -y, self.z), (-1, 0, 0)])

        # always restore context
        self.restore_context()


class ARCHIPACK_PT_myobject(Panel):
    bl_idname = "ARCHIPACK_PT_myobject"
    bl_label = "MyObject"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        # ensure your object panel only show when active object is the right one
        return archipack_myobject.poll(context.active_object)

    def draw(self, context):
        o = context.active_object
        if not archipack_myobject.filter(o):
            return
        layout = self.layout

        # retrieve datablock of your object
        props = archipack_myobject.datablock(o)

        # Manipulate mode operator
        layout.operator('archipack.myobject_manipulate', icon='HAND')

        box = layout.box()
        row = box.row(align=True)

        # Presets operators
        row.operator("archipack.myobject_preset_menu",
                     text=bpy.types.ARCHIPACK_OT_myobject_preset_menu.bl_label)
        row.operator("archipack.myobject_preset",
                      text="",
                      icon='ZOOMIN')
        row.operator("archipack.myobject_preset",
                      text="",
                      icon='ZOOMOUT').remove_active = True

        row = layout.row()
        box = row.box()
        box.label(text="Size")
        box.prop(props, 'x')
        box.prop(props, 'y')
        box.prop(props, 'z')


class ARCHIPACK_OT_myobject(ArchipackCreateTool, Operator):
    bl_idname = "archipack.myobject"
    bl_label = "Myobject"
    bl_description = "Create Myobject"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    def create(self, context):

        # Create an empty mesh datablock
        m = bpy.data.meshes.new("Myobject")

        # Create an object using the mesh datablock
        o = bpy.data.objects.new("Myobject", m)

        # Add your properties on mesh datablock
        d = m.archipack_myobject.add()

        # Link object into scene
        context.scene.objects.link(o)

        # select and make active
        o.select = True
        context.scene.objects.active = o

        # Load preset into datablock
        self.load_preset(d)

        # add a material
        self.add_material(o)
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = bpy.context.scene.cursor_location
            o.select = True
            context.scene.objects.active = o

            # Start manipulate mode
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_myobject_preset_menu(PresetMenuOperator, Operator):
    bl_idname = "archipack.myobject_preset_menu"
    bl_label = "Myobject preset"
    preset_subdir = "archipack_myobject"


class ARCHIPACK_OT_myobject_preset(ArchipackPreset, Operator):
    """Add a Myobject Preset"""
    bl_idname = "archipack.myobject_preset"
    bl_label = "Add Myobject preset"
    preset_menu = "ARCHIPACK_OT_myobject_preset_menu"

    @property
    def blacklist(self):
        return ['manipulators']


class ARCHIPACK_OT_myobject_manipulate(Operator):
    bl_idname = "archipack.myobject_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return archipack_myobject.poll(context.active_object)

    def invoke(self, context, event):
        d = archipack_myobject.datablock(context.active_object)
        d.manipulable_invoke(context)
        return {'FINISHED'}


def register():
    bpy.utils.register_class(archipack_myobject)
    Mesh.archipack_myobject = CollectionProperty(type=archipack_myobject)
    bpy.utils.register_class(ARCHIPACK_PT_myobject)
    bpy.utils.register_class(ARCHIPACK_OT_myobject)
    bpy.utils.register_class(ARCHIPACK_OT_myobject_preset_menu)
    bpy.utils.register_class(ARCHIPACK_OT_myobject_preset)
    bpy.utils.register_class(ARCHIPACK_OT_myobject_manipulate)


def unregister():
    bpy.utils.unregister_class(archipack_myobject)
    del Mesh.archipack_myobject
    bpy.utils.unregister_class(ARCHIPACK_PT_myobject)
    bpy.utils.unregister_class(ARCHIPACK_OT_myobject)
    bpy.utils.unregister_class(ARCHIPACK_OT_myobject_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_OT_myobject_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_myobject_manipulate)
