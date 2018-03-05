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
from bpy.types import Operator, PropertyGroup, Curve, Panel
from bpy.props import (
    FloatProperty, EnumProperty, BoolProperty,
    CollectionProperty, IntProperty, StringProperty
    )
from mathutils import Vector
from .archipack_manipulator import Manipulable
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_object import ArchipackCreateTool, ArchipackObject


def update(self, context):
    self.update(context)


# Paper formats (mm)
# Note: must contain a 'USER_DEFINED' key
paper_formats = {
    'A0': [1189, 841],
    'A1': [841, 594],
    'A2': [594, 420],
    'A3': [420, 297],
    'A4': [297, 210],
    'A5': [210, 148],
    'A6': [148, 105],
    'B0': [1414, 1000],
    'B1': [1000, 707],
    'B2': [707, 500],
    'B3': [500, 353],
    'B4': [353, 250],
    'B5': [250, 176],
    'B6': [176, 125],
    'C0': [1297, 917],
    'C1': [917, 648],
    'C2': [648, 458],
    'C3': [458, 324],
    'C4': [324, 229],
    'C5': [229, 162],
    'C6': [162, 114],
    'D1': [771, 545],
    'D2': [545, 385],
    'D3': [385, 272],
    'D4': [272, 192],
    'D5': [192, 136],
    'D6': [136, 96],
    'E3': [560, 400],
    'E4': [400, 280],
    'E5': [280, 200],
    'E6': [200, 140],
    'Arch A': [304.8, 228.6],
    'Arch B': [457.2, 304.8],
    'Arch C': [609.6, 457.2],
    'Arch D': [914.4, 609.6],
    'Arch E': [1219.2, 914.4],
    'Arch E1': [1066.8, 762],
    'USER_DEFINED': [297, 210]
    }


# Enumerator for paper chooser
paper_keys = list(paper_formats.keys())
paper_keys.sort()
paper_enum = tuple([(
        k,
        "{} ({}/{} mm)".format(k, paper_formats[k][0], paper_formats[k][1]),
        "{} ({}/{} mm)".format(k, paper_formats[k][0], paper_formats[k][1])
        ) for k in paper_keys])


class archipack_layout(ArchipackObject, Manipulable, PropertyGroup):
    """ Archipack Layout frame"""

    detail_scale = EnumProperty(
        name="Scale",
        items=(
            ('1000', '1/1000', '1/1000'),
            ('500', '1/500', '1/500'),
            ('200', '1/200', '1/200'),
            ('100', '1/100', '1/100'),
            ('50', '1/50', '1/50'),
            ('20', '1/20', '1/20'),
            ('10', '1/10', '1/10'),
            ('5', '1/5', '1/5'),
            ('2', '1/2', '1/2'),
            ('1', '1/1', '1/1'),
            ('USER_DEFINED', 'User defined', 'User defined')
        ),
        default='50',
        update=update
        )
    detail_user = FloatProperty(
        name="User defined scale",
        default=0.05,
        precision=4,
        update=update
        )
    paper_format = EnumProperty(
        name="Paper format",
        items=paper_enum,
        default='A4',
        update=update
        )
    paper_w = FloatProperty(
        name="Width (mm)",
        default=297,
        precision=4,
        update=update
        )
    paper_h = FloatProperty(
        name="Height (mm)",
        default=210,
        precision=4,
        update=update
        )
    paper_orientation = EnumProperty(
        name='Orientation',
        items=(
            ('PORTRAIT', 'Portrait', 'Portrait'),
            ('LANDSCAPE', 'Landscape', 'Landscape')
        ),
        default='LANDSCAPE',
        update=update
        )

    resolution = IntProperty(
        name="resolution (dpi)",
        default=90
        )

    scale_prefix = StringProperty(
            name="Scale prefix",
            default="Scale",
            update=update
            )
    text_size = FloatProperty(
            name="Text size",
            default=0.1,
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
    def scale(self):
        if self.detail_scale == 'USER_DEFINED':
            s = self.detail_user
        else:
            s = float(self.detail_scale)
        return s

    def canvas_size(self, context):
        """
          Canvas size in world coordinates (m)
        """
        w, h = self.paper_size
        # Unit scale: 1 Blender unit = 1 m * this factor
        scale = context.scene.unit_settings.scale_length * self.scale
        return w / 1000 * scale, h / 1000 * scale

    def canvas_scale(self, context):
        """
         * scale factor from world to paper unit (mm)
        """
        # Unit scale: 1 Blender unit = 1 m * this factor
        return context.scene.unit_settings.scale_length * self.scale / 1000

    @property
    def pixel_size(self):
        # pixels / mm
        return self.resolution / 25.4

    @property
    def paper_size(self):
        """
          Paper size in mm
        """
        paper_format = self.paper_format
        w, h = paper_formats[paper_format]
        if paper_format == 'USER_DEFINED':
            w, h = self.paper_w, self.paper_h
        elif self.paper_orientation == 'PORTRAIT':
            w, h = h, w
        return w, h

    def setup_manipulators(self):
        if len(self.manipulators) < 1:
            # add manipulator for x property
            s = self.manipulators.add()
            s.prop1_name = "paper_w"
            s.type_key = 'SIZE'

    def _add_spline(self, curve, closed, coords):
        spline = curve.splines.new('POLY')
        spline.use_endpoint_u = False
        spline.use_cyclic_u = closed
        spline.points.add(len(coords) - 1)
        for i, coord in enumerate(coords):
            x, y, z = coord
            spline.points[i].co = (x, y, z, 1)

    def update_child(self, context, o, location):

        t_o = None
        text = None
        # find text if any
        for child in o.children:
            if child.type == 'FONT':
                t_o = child
                text = t_o.data
                break

        if t_o is None:
            name = "Scale"
            if self.scale_prefix != "":
                name = self.scale_prefix
            text = bpy.data.curves.new(name, type='FONT')
            text.dimensions = '2D'
            t_o = bpy.data.objects.new(name, text)
            context.scene.objects.link(t_o)
            t_o.parent = o

        t_o.location = location
        if self.detail_scale == 'USER_DEFINED':
            text.body = "{}: 1/{}".format(self.scale_prefix, self.user_scale)
        else:
            text.body = "{}: 1/{}".format(self.scale_prefix, self.detail_scale)

        text.size = self.text_size
        text.align_x = 'CENTER'
        return t_o

    def update(self, context):

        # provide support for "copy to selected"
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # dynamically create manipulators when needed
        self.setup_manipulators()

        w, h = self.canvas_size(context)
        p0 = Vector((0, 0, 0))
        p1 = Vector((w, 0, 0))
        p2 = Vector((w, h, 0))
        p3 = Vector((0, h, 0))
        # sides
        curve = o.data
        curve.splines.clear()
        self._add_spline(curve, True, [p0, p1, p2, p3])

        # update manipulators location (3d location in object coordsystem)
        self.manipulators[0].set_pts([p0, p1, (-1, 0, 0)])

        self.update_child(context, o, Vector((0.5 * w, 0.5 * self.text_size, 0)))

        # always restore context
        self.restore_context(context)


class ARCHIPACK_PT_layout(Panel):
    bl_idname = "ARCHIPACK_PT_layout"
    bl_label = "Layout"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        # ensure your object panel only show when active object is the right one
        return archipack_layout.filter(context.active_object)

    def draw(self, context):
        o = context.active_object
        if not archipack_layout.filter(o):
            return
        layout = self.layout

        # retrieve datablock of your object
        props = archipack_layout.datablock(o)

        # Manipulate mode operator
        layout.operator('archipack.manipulate', icon='HAND')

        box = layout.box()
        row = box.row(align=True)

        # Presets operators
        row.operator("archipack.layout_preset_menu",
                     text=bpy.types.ARCHIPACK_OT_layout_preset_menu.bl_label)
        row.operator("archipack.layout_preset",
                      text="",
                      icon='ZOOMIN')
        row.operator("archipack.layout_preset",
                      text="",
                      icon='ZOOMOUT').remove_active = True

        box = layout.box()
        box.prop(props, 'detail_scale')
        if props.detail_scale == 'USER_DEFINED':
            box.prop(props, 'detail_user')
        box.prop(props, 'paper_format')
        if props.paper_format == 'USER_DEFINED':
            box.prop(props, 'paper_w')
            box.prop(props, 'paper_h')
        box.prop(props, 'paper_orientation')
        box = layout.box()
        box.prop(props, 'scale_prefix')
        box.prop(props, 'text_size')
        layout.operator("archipack.export_svg", text="Export .svg")


class ARCHIPACK_OT_layout(ArchipackCreateTool, Operator):
    bl_idname = "archipack.layout"
    bl_label = "Layout"
    bl_description = "Create Layout"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    def create(self, context):

        # Create an empty curve datablock
        c = bpy.data.curves.new("Layout", type='CURVE')
        c.dimensions = '3D'
        o = bpy.data.objects.new("Layout", c)
        o.lock_scale = (True, True, True)
        # Add your properties on curve datablock
        d = c.archipack_layout.add()

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
            x, y, z = bpy.context.scene.cursor_location
            o.location = Vector((x, y, 0))
            o.select = True
            context.scene.objects.active = o

            # Start manipulate mode
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_layout_preset_menu(PresetMenuOperator, Operator):
    bl_idname = "archipack.layout_preset_menu"
    bl_label = "Layout preset"
    preset_subdir = "archipack_layout"


class ARCHIPACK_OT_layout_preset(ArchipackPreset, Operator):
    """Add a Layout Preset"""
    bl_idname = "archipack.layout_preset"
    bl_label = "Add Layout preset"
    preset_menu = "ARCHIPACK_OT_layout_preset_menu"

    @property
    def blacklist(self):
        return ['manipulators']


def register():
    bpy.utils.register_class(archipack_layout)
    Curve.archipack_layout = CollectionProperty(type=archipack_layout)
    bpy.utils.register_class(ARCHIPACK_PT_layout)
    bpy.utils.register_class(ARCHIPACK_OT_layout)
    bpy.utils.register_class(ARCHIPACK_OT_layout_preset_menu)
    bpy.utils.register_class(ARCHIPACK_OT_layout_preset)


def unregister():
    bpy.utils.unregister_class(archipack_layout)
    del Curve.archipack_layout
    bpy.utils.unregister_class(ARCHIPACK_PT_layout)
    bpy.utils.unregister_class(ARCHIPACK_OT_layout)
    bpy.utils.unregister_class(ARCHIPACK_OT_layout_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_OT_layout_preset)
