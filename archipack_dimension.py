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
    FloatProperty, EnumProperty, BoolProperty, CollectionProperty
    )
from mathutils import Vector, Matrix
from math import cos, sin, pi
from .archipack_manipulator import Manipulable
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_object import ArchipackCreateTool, ArchipackObject
from .archipack_gl import GlText


def update(self, context):
    self.update(context)


class archipack_dimension(ArchipackObject, Manipulable, PropertyGroup):
    """ Archipack Dimension curve"""
    distance = FloatProperty(
            description="Symbol lateral distance",
            name="Distance",
            default=1.0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    size = FloatProperty(
            description="Mesured dimension",
            name="Size",
            default=1.0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    origin = FloatProperty(
            description="Altitude offset",
            name="Origin",
            default=0.0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    angle = FloatProperty(
            description="Mesured dimension",
            name="Angle",
            default=0.0,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    area = FloatProperty(
            description="Mesured dimension",
            name="Area",
            default=1.0,
            unit='AREA', subtype='DISTANCE',
            update=update
            )
    volume = FloatProperty(
            description="Mesured dimension",
            name="Volume",
            default=1.0,
            unit='VOLUME', subtype='DISTANCE',
            update=update
            )
    type = EnumProperty(
            name="Type",
            items=(
                ("SIZE", "Size", "Size"),
                ("ANGLE", "Angle", "Angle"),
                ("ALTITUDE", "Altitude", "Altitude"),
                ("AREA", "Area", "Area"),
                ("VOLUME", "Volume", "Volume")
                ),
            default="SIZE",
            update=update
            )
    symbol_size = FloatProperty(
            name="Symbol Size",
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    symbol_thickness = FloatProperty(
            name="Thickness",
            description="Arrow thickness factor",
            default=0.25, min=0.01,
            update=update
            )

    symbol_type = EnumProperty(
            name="Symbol shape",
            items=(
                ("ARROW_INSIDE", "Arrow inside", "Arrow inside"),
                ("ARROW_OUTSIDE", "Arrow outside", "Arrow outside"),
                ("TICK", "Architectural tick", "Architectural tick"),
                ("DOT", "Dot", "Dot")
                ),
            default="TICK",
            update=update
            )
    text_size = FloatProperty(
            name="Text Size",
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    text_angle = FloatProperty(
            name="Text Angle",
            default=0.0,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    text_location = EnumProperty(
            items=(
                ("TOP_CENTER", "Top center", "Top center"),
                ("BOTTOM_CENTER", "Bottom center", "Bottom center")
                ),
            default="TOP_CENTER",
            update=update
            )

    auto_update = BoolProperty(
            # Wont save auto_update state in any case
            options={'SKIP_SAVE'},
            default=True,
            update=update
            )

    def text(self, context, value, precision=2):

        dimension = 1

        if self.type == 'AREA':
            dimension = 2
        elif self.type == 'VOLUME':
            dimension = 3

        # Either 'ANGLE' or 'SIZE'
        unit_type = self.type
        unit_mode = 'METER'

        if self.type in {'AREA', 'VOLUME', 'ALTITUDE'}:
            unit_type = 'SIZE'

        if self.type == 'ANGLE':
            unit_mode = 'AUTO'

        label = GlText(
            label="",
            value=value,
            precision=precision,
            unit_mode=unit_mode,
            unit_type=unit_type,
            dimension=dimension
            )
        return label.add_units(context)

    def setup_manipulators(self):
        if len(self.manipulators) < 1:
            # add manipulator for x property
            s = self.manipulators.add()
            s.prop1_name = "size"
            s.type_key = 'SIZE'

    def _add_spline(self, curve, closed, coords):
        spline = curve.splines.new('POLY')
        spline.use_endpoint_u = False
        spline.use_cyclic_u = closed
        spline.points.add(len(coords) - 1)
        for i, coord in enumerate(coords):
            x, y, z = coord
            spline.points[i].co = (x, y, z, 1)

    def text_frame(self, context, curve, p0, t_o):
        # need scene update to retrieve
        # current text dimension
        context.scene.update()
        d = t_o.dimensions
        a = -self.text_angle
        rM = Matrix([
            [cos(a), sin(a), 0, 0],
            [-sin(a), cos(a), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        p0 += rM * Vector((0, 0.5 * d.y, 0))
        x = rM * Vector((0.5 * (d.x + d.y), 0, 0))
        y = rM * Vector((0, d.y, 0))
        self._add_spline(curve, True, [p0 - x - y, p0 + x - y, p0 + x + y, p0 - x + y])

    def update_child(self, context, o, center, value):
        t_o = None
        text = None
        # find text if any
        for child in o.children:
            if child.type == 'FONT':
                t_o = child
                text = t_o.data
                break

        if t_o is None:
            text = bpy.data.curves.new(o.name, type='FONT')
            text.dimensions = '2D'
            t_o = bpy.data.objects.new(o.name, text)
            context.scene.objects.link(t_o)
            t_o.parent = o

        t_o.location = center
        t_o.rotation_euler.z = self.text_angle

        text.body = self.text(context, value)
        text.size = self.text_size
        text.align_x = 'CENTER'
        return t_o

    def symbol(self, curve, p0, side=1, orientation='X'):
        """
         p0: location of tip
         side: side factor in [-1, 1]
        """
        symbol_type = self.symbol_type

        s = side * self.symbol_size

        if 'ARROW' in symbol_type:

            if orientation == 'X':
                x = Vector((s, 0, 0))
                y = Vector((0, s * self.symbol_thickness, 0))
            else:
                x = Vector((0, s, 0))
                y = Vector((s * self.symbol_thickness, 0, 0))

            if 'ARROW_OUTSIDE' in symbol_type:
                x = -x

            self._add_spline(curve, True, [p0 + x + y, p0 + x - y, p0])

        elif symbol_type == 'DOT':
            seg = 16
            da = 2 * pi / seg
            r = 0.5 * s
            self._add_spline(curve, True, [
                p0 + r * Vector((cos(da * i), sin(da * i), 0))
                for i in range(seg)
                ])

        elif symbol_type == 'TICK':
            r = s * 0.25 * (2 ** 0.5)
            x = Vector((r, r, 0))
            self._add_spline(curve, False, [p0 - x, p0 + x])

        elif symbol_type == 'CROSS':
            r = s * 0.25 * (2 ** 0.5)
            x = Vector((r, r, 0))
            x2 = Vector((r, -r, 0))
            self._add_spline(curve, False, [p0 - x, p0 + x])
            self._add_spline(curve, False, [p0 - x2, p0 + x2])

    def update(self, context):

        # provide support for "copy to selected"
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # dynamically create manipulators when needed
        self.setup_manipulators()

        curve = o.data
        curve.splines.clear()

        s = self.symbol_size
        dist = self.distance
        if dist > 0:
            side = 1
        else:
            side = -1
        p0 = Vector((0, 0, 0))

        # update your mesh from parameters
        if self.type == 'SIZE':
            p1 = Vector((self.size, 0, 0))
            p2 = Vector((0, dist, 0))
            p3 = Vector((self.size, dist, 0))
            s = self.symbol_size
            x = Vector((s, 0, 0))
            y = Vector((0, s, 0))
            # sides
            self._add_spline(curve, False, [p0, p2 + side * y])
            self._add_spline(curve, False, [p1, p3 + side * y])
            # bar
            if self.symbol_type == 'ARROW_INSIDE':
                self._add_spline(curve, False, [p2 + x, p3 - x])

            elif self.symbol_type == 'ARROW_OUTSIDE':
                self._add_spline(curve, False, [p2, p3])

            else:
                self._add_spline(curve, False, [p2 - x, p3 + x])

            # arrows
            self.symbol(curve, p2, 1)
            self.symbol(curve, p3, -1)
            center = p2 + 0.5 * Vector((self.size, self.text_size, 0))
            t_o = self.update_child(context, o, center, self.size)
            self.manipulators[0].set_pts([p0, p1, (-dist, 0, 0)])

        elif self.type == 'ALTITUDE':
            p1 = Vector((0, self.size, 0))
            p2 = Vector((dist, 0, 0))
            p3 = Vector((dist, self.size, 0))
            s = self.symbol_size
            y = Vector((side * s, 0, 0))
            # sides
            self._add_spline(curve, False, [p1, p3 + y])

            # arrows
            self.symbol(curve, p3, -1, 'Y')
            center = p3 + 0.5 * Vector((0, self.text_size, 0))
            t_o = self.update_child(context, o, center, self.size + self.origin)
            self.manipulators[0].set_pts([p0, p1, (dist, 0, 0)], normal=Vector((0, 0, 1)))

        elif self.type == 'ANGLE':
            t_o = self.update_child(context, o, p0, self.angle)

        elif self.type == 'AREA':
            t_o = self.update_child(context, o, p0, self.area)
            self.text_frame(context, curve, p0, t_o)

        elif self.type == 'VOLUME':
            t_o = self.update_child(context, o, p0, self.volume)
            self.text_frame(context, curve, p0, t_o)

        # always restore context
        self.restore_context(context)


class ARCHIPACK_PT_dimension(Panel):
    bl_idname = "ARCHIPACK_PT_dimension"
    bl_label = "Dimension"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        # ensure your object panel only show when active object is the right one
        return archipack_dimension.filter(context.active_object)

    def draw(self, context):
        o = context.active_object
        if not archipack_dimension.filter(o):
            return
        layout = self.layout

        # retrieve datablock of your object
        props = archipack_dimension.datablock(o)

        # Manipulate mode operator
        row = layout.row(align=True)
        row.operator('archipack.manipulate', icon='HAND')
        row.operator('archipack.dimension', text="Delete", icon='ERROR').mode = 'DELETE'

        box = layout.box()
        row = box.row(align=True)

        # Presets operators
        row.operator("archipack.dimension_preset_menu",
                     text=bpy.types.ARCHIPACK_OT_dimension_preset_menu.bl_label)
        row.operator("archipack.dimension_preset",
                      text="",
                      icon='ZOOMIN')
        row.operator("archipack.dimension_preset",
                      text="",
                      icon='ZOOMOUT').remove_active = True

        box = layout.box()
        box.prop(props, 'type')
        box.prop(props, 'distance')
        if props.type == 'SIZE':
            box.prop(props, 'size')
        elif props.type == 'ANGLE':
            box.prop(props, 'angle')
        elif props.type == 'ALTITUDE':
            box.prop(props, 'size', text='Altitude')
            box.prop(props, 'origin')
        elif props.type == 'AREA':
            box.prop(props, 'area')
        elif props.type == 'VOLUME':
            box.prop(props, 'volume')
        box = layout.box()
        box.prop(props, 'symbol_type')
        box.prop(props, 'symbol_size')
        if "ARROW" in props.symbol_type:
            box.prop(props, 'symbol_thickness')
        box.prop(props, 'text_size')
        box.prop(props, 'text_angle')


class ARCHIPACK_OT_dimension(ArchipackCreateTool, Operator):
    bl_idname = "archipack.dimension"
    bl_label = "Dimension"
    bl_description = "Create Dimension"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    mode = EnumProperty(
            items=(
            ('CREATE', 'Create', '', 0),
            ('DELETE', 'Delete', '', 1)
            ),
            default='CREATE'
            )

    def create(self, context):

        # Create an empty curve datablock
        c = bpy.data.curves.new("Dimension", type='CURVE')
        c.dimensions = '2D'
        o = bpy.data.objects.new("Dimension", c)

        # Add your properties on curve datablock
        d = c.archipack_dimension.add()

        # Link object into scene
        context.scene.objects.link(o)

        # select and make active
        o.select = True
        context.scene.objects.active = o

        # Load preset into datablock
        self.load_preset(d)

        # add a material
        self.add_material(o)

        for child in o.children:
            m = child.archipack_material.add()
            m.category = "dimension"
            m.material = o.archipack_material[0].material
        return o

    def delete(self, context, o):
        if archipack_dimension.filter(o):
            self.delete_object(context, o)

    def execute(self, context):
        if context.mode == "OBJECT":
            if self.mode == 'DELETE':
                o = context.active_object
                bpy.ops.archipack.disable_manipulate()
                self.delete(context, o)
            else:
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


class ARCHIPACK_OT_dimension_preset_menu(PresetMenuOperator, Operator):
    bl_idname = "archipack.dimension_preset_menu"
    bl_label = "Dimension preset"
    preset_subdir = "archipack_dimension"


class ARCHIPACK_OT_dimension_preset(ArchipackPreset, Operator):
    """Add a Dimension Preset"""
    bl_idname = "archipack.dimension_preset"
    bl_label = "Add Dimension preset"
    preset_menu = "ARCHIPACK_OT_dimension_preset_menu"

    @property
    def blacklist(self):
        return ['manipulators']


def register():
    bpy.utils.register_class(archipack_dimension)
    Curve.archipack_dimension = CollectionProperty(type=archipack_dimension)
    bpy.utils.register_class(ARCHIPACK_PT_dimension)
    bpy.utils.register_class(ARCHIPACK_OT_dimension)
    bpy.utils.register_class(ARCHIPACK_OT_dimension_preset_menu)
    bpy.utils.register_class(ARCHIPACK_OT_dimension_preset)


def unregister():
    bpy.utils.unregister_class(archipack_dimension)
    del Curve.archipack_dimension
    bpy.utils.unregister_class(ARCHIPACK_PT_dimension)
    bpy.utils.unregister_class(ARCHIPACK_OT_dimension)
    bpy.utils.unregister_class(ARCHIPACK_OT_dimension_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_OT_dimension_preset)
