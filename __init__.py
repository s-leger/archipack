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
#
# ----------------------------------------------------------

bl_info = {
    'name': 'Archipack',
    'description': 'Architectural objects and 2d polygons detection from unordered splines',
    'author': 's-leger',
    'license': 'GPL',
    'deps': '',
    'version': (1, 4, 1),
    'blender': (2, 7, 8),
    'location': 'View3D > Tools > Create > Archipack',
    'warning': '',
    'wiki_url': 'https://s-leger.github.io/archipack/index.html',
    'tracker_url': 'https://github.com/s-leger/archipack/issues',
    'link': 'https://github.com/s-leger/archipack',
    'support': 'COMMUNITY',
    'category': 'Add Mesh'
    }


__version__ = ".".join([str(i) for i in bl_info['version']])


import os
import importlib as imp


submodules = (
    'iconmanager',
    'progressbar',
    'material',
    'snap',
    'manipulator',
    'curveman',
    'throttle',
    'segments',
    'dimension',
    'reference_point',
    'autoboolean',
    'door',
    'window',
    'stair',
    'wall',
    'wall2',
    'slab',
    'roof',
    'fence',
    'truss',
    'custom',
    # toolkit
    'floor',
    'floor_heating',
    'blind',
    'kitchen',
    'molding',
    'rendering',
    'section',
    'animation',
    # envi
    'io',
    '2d_layout',
    'io_export_svg',
    'polylines',
)


if "bpy" in locals():
    glob = globals()
    for sub in submodules:
        imp.reload(glob["archipack_{}".format(sub)])
    imp.reload(addon_updater_ops)
    print("archipack: reload ready")

else:
    glob = globals()
    for sub in submodules:
        glob["archipack_{}".format(sub)] = imp.import_module(".archipack_{}".format(sub), __package__)
    from . import addon_updater_ops
    print("archipack: ready")


# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import (
    Panel, AddonPreferences, Menu
    )
from bpy.props import (
    StringProperty, BoolProperty,
    IntProperty, FloatProperty,
    FloatVectorProperty
    )


icon_man = None


# ----------------------------------------------------
# Addon preferences
# ----------------------------------------------------

def update_panel(self, context):
    try:
        bpy.utils.unregister_class(TOOLS_PT_Archipack_PolyLib)
        bpy.utils.unregister_class(TOOLS_PT_Archipack_Tools)
        bpy.utils.unregister_class(TOOLS_PT_Archipack_Create)
        bpy.utils.unregister_class(TOOLS_PT_Archipack_Create_2d)
        bpy.utils.unregister_class(TOOLS_PT_Archipack_Create_Custom)
        bpy.utils.unregister_class(TOOLS_PT_Archipack_About)
    except:
        pass
    prefs = context.user_preferences.addons[__name__].preferences
    TOOLS_PT_Archipack_Create.bl_category = prefs.create_category
    bpy.utils.register_class(TOOLS_PT_Archipack_Create)
    TOOLS_PT_Archipack_Create_2d.bl_category = prefs.create_category
    bpy.utils.register_class(TOOLS_PT_Archipack_Create_2d)
    TOOLS_PT_Archipack_Create_Custom.bl_category = prefs.create_category
    bpy.utils.register_class(TOOLS_PT_Archipack_Create_Custom)
    TOOLS_PT_Archipack_Tools.bl_category = prefs.tools_category
    bpy.utils.register_class(TOOLS_PT_Archipack_Tools)
    TOOLS_PT_Archipack_PolyLib.bl_category = prefs.tools_category
    bpy.utils.register_class(TOOLS_PT_Archipack_PolyLib)
    TOOLS_PT_Archipack_About.bl_category = prefs.create_category
    bpy.utils.register_class(TOOLS_PT_Archipack_About)


class Archipack_Pref(AddonPreferences):
    bl_idname = __name__
    tools_category = StringProperty(
            name="Tools",
            description="Choose a name for the category of the Tools panel",
            default="Archipack",
            update=update_panel
            )
    create_category = StringProperty(
            name="Create",
            description="Choose a name for the category of the Create panel",
            default="Archipack",
            update=update_panel
            )
    create_submenu = BoolProperty(
            name="Use Sub-menu",
            description="Put Achipack's object into a sub menu (shift+a)",
            default=True
            )
    max_style_draw_tool = BoolProperty(
            name="Draw a wall use 3dsmax style",
            description="Reverse clic / release & drag cycle for Draw a wall",
            default=True
            )
    throttle_enable = BoolProperty(
            name="Quick edit",
            description="When enabled, prevent complex objects to update in real time",
            default=True
            )
    throttle_delay = IntProperty(
            name="Delay",
            description="Quick edit, how much time to wait before updating complex objects (seconds)",
            default=5
            )
    # Arrow sizes (world units)
    arrow_size = FloatProperty(
            name="Arrow",
            description="Manipulators arrow size (blender units)",
            default=0.05
            )
    # Handle area size (pixels)
    handle_size = IntProperty(
            name="Handle",
            description="Manipulators handle sensitive area size (pixels)",
            min=2,
            default=10
            )
    constant_handle_size = BoolProperty(
            name="Constant handle size",
            description="When checked, handle size on scree remains constant (handle size pixels)",
            default=False
            )
    handle_colour_selected = FloatVectorProperty(
            name="Selected handle colour",
            description="Handle color when selected",
            subtype='COLOR_GAMMA',
            default=(0.0, 0.0, 0.7, 1.0),
            size=4,
            min=0, max=1
            )
    handle_colour_inactive = FloatVectorProperty(
            name="Inactive handle colour",
            description="Handle color when disabled",
            subtype='COLOR_GAMMA',
            default=(0.3, 0.3, 0.3, 1.0),
            size=4,
            min=0, max=1
            )
    handle_colour_normal = FloatVectorProperty(
            name="Handle colour normal",
            description="Base handle color when not selected",
            subtype='COLOR_GAMMA',
            default=(1.0, 1.0, 1.0, 1.0),
            size=4,
            min=0, max=1
            )
    handle_colour_hover = FloatVectorProperty(
            name="Handle colour hover",
            description="Handle color when mouse hover",
            subtype='COLOR_GAMMA',
            default=(1.0, 1.0, 0.0, 1.0),
            size=4,
            min=0, max=1
            )
    handle_colour_active = FloatVectorProperty(
            name="Handle colour active",
            description="Handle colour when moving",
            subtype='COLOR_GAMMA',
            default=(1.0, 0.0, 0.0, 1.0),
            size=4,
            min=0, max=1
            )
    matlib_path = StringProperty(
            name="Folder path",
            subtype="DIR_PATH",
            description="Absolute path to material library folder",
            default=""
            )
    # Font sizes and basic colour scheme
    # kept outside of addon prefs until now
    # as for a generic toolkit it is not appropriate
    # we could provide a template for addon prefs
    # matching those one
    feedback_size_main = IntProperty(
            name="Main",
            description="Main title font size (pixels)",
            min=2,
            default=16
            )
    feedback_size_title = IntProperty(
            name="Title",
            description="Tool name font size (pixels)",
            min=2,
            default=14
            )
    feedback_size_shortcut = IntProperty(
            name="Shortcut",
            description="Shortcuts font size (pixels)",
            min=2,
            default=11
            )
    feedback_shortcut_area = FloatVectorProperty(
            name="Background Shortcut",
            description="Shortcut area background color",
            subtype='COLOR_GAMMA',
            default=(0, 0.4, 0.6, 0.2),
            size=4,
            min=0, max=1
            )
    feedback_title_area = FloatVectorProperty(
            name="Background Main",
            description="Title area background color",
            subtype='COLOR_GAMMA',
            default=(0, 0.4, 0.6, 0.5),
            size=4,
            min=0, max=1
            )
    feedback_colour_main = FloatVectorProperty(
            name="Font Main",
            description="Title color",
            subtype='COLOR_GAMMA',
            default=(0.95, 0.95, 0.95, 1.0),
            size=4,
            min=0, max=1
            )
    feedback_colour_key = FloatVectorProperty(
            name="Font Shortcut key",
            description="KEY label color",
            subtype='COLOR_GAMMA',
            default=(0.67, 0.67, 0.67, 1.0),
            size=4,
            min=0, max=1
            )
    feedback_colour_shortcut = FloatVectorProperty(
            name="Font Shortcut hint",
            description="Shortcuts text color",
            subtype='COLOR_GAMMA',
            default=(0.51, 0.51, 0.51, 1.0),
            size=4,
            min=0, max=1
            )
    experimental_features = BoolProperty(
            name="Experimental features",
            description="Enable experimental features (may be unstable)",
            default=False
            )

    # addon updater preferences
    auto_check_update = BoolProperty(
            name="Auto-check for Update",
            description="If enabled, auto-check for updates using an interval",
            default=False,
            )

    updater_intrval_months = IntProperty(
            name='Months',
            description="Number of months between checking for updates",
            default=0,
            min=0
            )
    updater_intrval_days = IntProperty(
            name='Days',
            description="Number of days between checking for updates",
            default=7,
            min=0,
            )
    updater_intrval_hours = IntProperty(
            name='Hours',
            description="Number of hours between checking for updates",
            default=0,
            min=0,
            max=23
            )
    updater_intrval_minutes = IntProperty(
            name='Minutes',
            description="Number of minutes between checking for updates",
            default=0,
            min=0,
            max=59
            )

    def draw(self, context):
        layout = self.layout
        box = layout.box()
        box.label("Setup actions (see Documentation)")
        row = box.row()
        col = row.column()
        col.label(text="Presets (may take up to 3 min):")
        col.operator("archipack.render_thumbs", icon="RENDER_RESULT")

        col = row.column()
        col.label(text="Material library:")
        col.prop(self, "matlib_path")

        box = layout.box()
        row = box.row()
        col = row.column()
        col.label(text="Tab Category:")
        col.prop(self, "tools_category")
        col.prop(self, "create_category")
        col.prop(self, "create_submenu")

        box = layout.box()
        box.label("Features")
        row = box.row()
        row.prop(self, "throttle_enable")
        row.prop(self, "throttle_delay")
        box.prop(self, "max_style_draw_tool")
        box.prop(self, "experimental_features")
        box = layout.box()
        row = box.row()
        split = row.split(percentage=0.5)
        col = split.column()
        col.label(text="Colors:")
        row = col.row(align=True)
        row.prop(self, "feedback_title_area")
        row = col.row(align=True)
        row.prop(self, "feedback_shortcut_area")
        row = col.row(align=True)
        row.prop(self, "feedback_colour_main")
        row = col.row(align=True)
        row.prop(self, "feedback_colour_key")
        row = col.row(align=True)
        row.prop(self, "feedback_colour_shortcut")

        row = col.row(align=True)
        row.prop(self, "handle_colour_normal")
        row = col.row(align=True)
        row.prop(self, "handle_colour_hover")
        row = col.row(align=True)
        row.prop(self, "handle_colour_active")
        row = col.row(align=True)
        row.prop(self, "handle_colour_selected")
        row = col.row(align=True)
        row.prop(self, "handle_colour_inactive")

        col = split.column()
        col.label(text="Font size:")
        col.prop(self, "feedback_size_main")
        col.prop(self, "feedback_size_title")
        col.prop(self, "feedback_size_shortcut")
        col.label(text="Manipulators:")
        col.prop(self, "arrow_size")
        col.prop(self, "handle_size")
        col.prop(self, "constant_handle_size")
        addon_updater_ops.update_settings_ui(self, context)


# ----------------------------------------------------
# Archipack panels
# ----------------------------------------------------

class TOOLS_PT_Archipack_PolyLib(Panel):
    bl_label = "2d to 3d"
    bl_idname = "TOOLS_PT_Archipack_PolyLib"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Tools"
    bl_context = "objectmode"

    @classmethod
    def poll(self, context):

        global archipack_polylines
        return ((archipack_polylines.vars_dict['select_polygons'] is not None) or
                (context.object is not None and context.object.type == 'CURVE'))

    def draw(self, context):
        global icon_man
        icons = icon_man["main"]
        layout = self.layout
        params = context.window_manager.archipack_polylib
        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        if params.polygonize_expand:
                row.prop(params, "polygonize_expand", icon='TRIA_DOWN', text="")
        else:
            row.prop(params, "polygonize_expand", icon='TRIA_RIGHT', text="")

        row.operator(
            "archipack.polylib_polygonize",
            icon_value=icons["detect"].icon_id,
            text='Detect'
            )

        if params.polygonize_expand:
            box.prop(params, "polygonize_bezier_resolution")
            box.prop(params, "polygonize_extend")
            box.prop(params, "polygonize_all_segs")

            box.operator(
                "archipack.polylib_pick_2d_polygons",
                icon_value=icons["selection"].icon_id,
                text='Polygons'
                )

            box.operator(
                "archipack.polylib_pick_2d_lines",
                icon_value=icons["selection"].icon_id,
                text='Lines'
                )

            box.operator(
                "archipack.polylib_pick_2d_points",
                icon_value=icons["selection"].icon_id,
                text='Points'
                )

            box.label(text="Walls")
            box.prop(params, "polygonize_thickness")

        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        if params.simplify_expand:
            row.prop(params, "simplify_expand", icon='TRIA_DOWN', text="")
        else:
            row.prop(params, "simplify_expand", icon='TRIA_RIGHT', text="")
        row.operator("archipack.polylib_simplify")
        if params.simplify_expand:
            box.prop(params, "simplify_bezier_resolution")
            box.prop(params, "simplify_tolerance")
            box.prop(params, "simplify_preserve_topology")

        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        if params.offset_expand:
            row.prop(params, "offset_expand", icon='TRIA_DOWN', text="")
        else:
            row.prop(params, "offset_expand", icon='TRIA_RIGHT', text="")
        row.operator("archipack.polylib_offset")
        if params.offset_expand:
            box.prop(params, "offset_bezier_resolution")
            box.prop(params, "offset_distance")
            box.prop(params, "offset_side")
            box.prop(params, "offset_resolution")
            box.prop(params, "offset_join_style")
            box.prop(params, "offset_mitre_limit")

        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        if params.buffer_expand:
            row.prop(params, "buffer_expand", icon='TRIA_DOWN', text="")
        else:
            row.prop(params, "buffer_expand", icon='TRIA_RIGHT', text="")
        row.operator("archipack.polylib_buffer")
        if params.buffer_expand:
            box.prop(params, "buffer_bezier_resolution")
            box.prop(params, "buffer_distance")
            box.prop(params, "buffer_side")
            box.prop(params, "buffer_resolution")
            box.prop(params, "buffer_join_style")
            box.prop(params, "buffer_cap_style")
            box.prop(params, "buffer_mitre_limit")

        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        if params.boolean_expand:
            row.prop(params, "boolean_expand", icon='TRIA_DOWN', text="")
        else:
            row.prop(params, "boolean_expand", icon='TRIA_RIGHT', text="")
        row.label(text="2d Boolean")
        if params.boolean_expand:
            box.operator("archipack.polylib_boolean", text="Active - Selected").opCode = 'DIFFERENCE'
            box.operator("archipack.polylib_boolean", text="Selected - Active").opCode = 'REVDIFFERENCE'
            box.operator("archipack.polylib_boolean", text="Intersection").opCode = 'INTERSECTION'
            box.operator("archipack.polylib_boolean", text="Union").opCode = 'UNION'
            box.operator("archipack.polylib_boolean", text="SymDifference").opCode = 'SYMDIFFERENCE'
            box.prop(params, "boolean_bezier_resolution")


class TOOLS_PT_Archipack_Tools(Panel):
    bl_label = "Tools"
    bl_idname = "TOOLS_PT_Archipack_Tools"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Tools"
    bl_context = "objectmode"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(self, context):
        return True

    def draw(self, context):
        prefs = context.user_preferences.addons[__name__].preferences
        layout = self.layout

        box = layout.box()
        box.label("Auto boolean")
        box.operator("archipack.auto_boolean", text="AutoBoolean", icon='AUTO')
        box.label("Apply holes")
        row = box.row(align=True)
        row.operator("archipack.apply_holes", text="selected").selected_only = True
        row.operator("archipack.apply_holes", text="all").selected_only = False

        box = layout.box()
        box.label("Kill parameters")
        row = box.row(align=True)
        row.operator("archipack.kill_archipack", text="selected", icon="ERROR").selected_only = True
        row.operator("archipack.kill_archipack", text="all", icon="ERROR").selected_only = False

        if prefs.experimental_features:
            box = layout.box()
            box.label("Animation support")
            row = box.row(align=True)
            row.operator("archipack.animation", text="Add", icon='ZOOMIN').mode = 'ENABLE'
            row.operator("archipack.animation", text="Remove", icon='ZOOMOUT').mode = 'DISABLE'
            row.operator("archipack.animation", text="Clear", icon='X').mode = 'CLEAR'

        """
        box = layout.box()
        box.label("Generate preset thumbs")
        box.operator("archipack.render_thumbs", icon="RENDER_RESULT")
        """


class TOOLS_PT_Archipack_About(Panel):
    bl_label = "About"
    bl_idname = "TOOLS_PT_Archipack_About"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Create"
    bl_context = "objectmode"

    def draw(self, context):
        layout = self.layout
        box = layout
        box.label(text="Archipack {}".format(__version__))
        box.label(text="by Stephen Leger")


class TOOLS_PT_Archipack_Create(Panel):
    bl_label = "Add Objects"
    bl_idname = "TOOLS_PT_Archipack_Create"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Create"
    bl_context = "objectmode"

    @classmethod
    def poll(self, context):
        return True

    def draw(self, context):
        global icon_man
        prefs = context.user_preferences.addons[__name__].preferences
        addon_updater_ops.check_for_update_background(context)

        icons = icon_man["main"]
        layout = self.layout
        box = layout.box()
        box.prop(prefs, "throttle_enable", icon="MOD_MULTIRES")
        box.prop(prefs, "throttle_delay")

        box = layout
        row = box.row(align=True)
        row.operator("archipack.wall2_draw",
                    text="Wall",
                    icon_value=icons["wall"].icon_id
                    )
        row.operator("archipack.wall2_from_curve", text="", icon='CURVE_DATA')
        row = box.row(align=True)
        # col = row.column()
        # subrow = col.row(align=True)
        row.operator("archipack.window_preset_draw",
                    text="Window",
                    icon_value=icons["window"].icon_id
                    ).preset_operator = "archipack.window_draw"
        # col = row.column()
        # subrow = col.row(align=True)
        row = box.row(align=True)
        row.operator("archipack.door_preset_draw",
                    text="Door",
                    icon_value=icons["door"].icon_id
                    ).preset_operator = "archipack.door_draw"
        box.operator("archipack.stair_preset_create",
                    text="Stairs",
                    icon_value=icons["stair"].icon_id
                    ).preset_operator = "archipack.stair"
        row = box.row(align=True)
        row.operator("archipack.fence_preset_create",
                    text="Fence",
                    icon_value=icons["fence"].icon_id
                    ).preset_operator = "archipack.fence"
        row.operator("archipack.fence_from_curve", text="", icon='CURVE_DATA')
        box.operator("archipack.wall2_from_slab",
                    icon_value=icons["wall_from_slab"].icon_id)

        row = box.row(align=True)
        row.operator("archipack.slab_from_wall",
                    icon_value=icons["slab_from_wall"].icon_id
                    ).ceiling = False
        row.operator("archipack.slab_from_wall",
                    text="->Ceiling",
                    icon_value=icons["ceiling_from_wall"].icon_id
                    ).ceiling = True
        row.operator("archipack.slab_from_curve",
                    text="", icon='CURVE_DATA'
                    )

        row = box.row(align=True)
        row.operator("archipack.floor_preset_create",
                    text="Floor",
                    icon_value=icons["floor"].icon_id
                    ).preset_operator = "archipack.floor"
        row.operator("archipack.floor_preset_from_wall",
                    text="->Floor",
                    icon_value=icons["floor_from_wall"].icon_id
                    )
        row.operator("archipack.floor_preset_from_curve",
                    text="",
                    icon='CURVE_DATA')

        row = box.row(align=True)
        row.operator("archipack.molding_preset_create",
                    text="Molding",
                    icon_value=icons["molding"].icon_id
                    ).preset_operator = "archipack.molding"
        row.operator("archipack.molding_preset_from_wall",
                    text="->Molding",
                    icon_value=icons["molding_from_wall"].icon_id
                    )
        row.operator("archipack.molding_preset_from_curve",
                    text="",
                    icon='CURVE_DATA'
                    )

        row = box.row(align=True)
        row.operator("archipack.roof_preset_create",
                    text="Roof",
                    icon_value=icons["roof"].icon_id
                    ).preset_operator = "archipack.roof"
        row.operator("archipack.roof_from_wall",
                    text="->Roof",
                    icon_value=icons["roof_from_wall"].icon_id
                    )
        # toolkit
        # row = box.row(align=True)
        # row.operator("archipack.myobject")
        box.operator("archipack.kitchen_preset_create",
                    text="Kitchen",
                    icon_value=icons["kitchen"].icon_id
                    ).preset_operator = "archipack.kitchen"

        box.operator("archipack.blind_preset_create",
                    text="Blind",
                    icon_value=icons["blind"].icon_id
                    ).preset_operator = "archipack.blind"
        box.operator("archipack.truss",
                    icon_value=icons["truss"].icon_id
                    )

        addon_updater_ops.update_notice_box_ui(self, context)


class TOOLS_PT_Archipack_Create_2d(Panel):
    bl_label = "Add 2d Objects"
    bl_idname = "TOOLS_PT_Archipack_Create_2d"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Create"
    bl_context = "objectmode"

    @classmethod
    def poll(self, context):
        return True

    def draw(self, context):
        global icon_man
        prefs = context.user_preferences.addons[__name__].preferences

        icons = icon_man["main"]
        layout = self.layout

        box = layout
        box.operator("archipack.dimension_auto",
                    icon_value=icons["dimension_auto"].icon_id
                    ).mode = 'CREATE'
        box.operator("archipack.layout",
                    icon_value=icons["layout"].icon_id
                    )
        box.operator("archipack.section",
                    icon_value=icons["section"].icon_id
                    ).mode = 'CREATE'
        box.operator("archipack.section_camera",
                    text="Section cam",
                    icon='CAMERA_DATA'
                    ).mode = 'CREATE'

        if prefs.experimental_features:
            box = layout.box()
            box.label(text="Experimental features")
            box.operator("archipack.floor_heating")
            box.operator("archipack.dimension")


class TOOLS_PT_Archipack_Create_Custom(Panel):
    bl_label = "Custom Objects"
    bl_idname = "TOOLS_PT_Archipack_Create_Custom"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Create"
    bl_context = "objectmode"

    @classmethod
    def poll(self, context):
        return True

    def draw(self, context):
        global icon_man
        prefs = context.user_preferences.addons[__name__].preferences

        # icons = icon_man["main"]
        layout = self.layout
        box = layout
        row = box.row(align=True)
        row.operator("archipack.custom_wall", text="Custom wall")
        row.operator("archipack.custom_wall_remove", text="", icon='X')

        row = box.row(align=True)
        row.operator("archipack.custom_hole", text="Custom hole")
        row.operator("archipack.custom_hole_remove", text="", icon='X')

        if prefs.experimental_features:
            box = layout
            row = box.row(align=True)
            row.operator("archipack.make_custom", text="Custom", icon='MONKEY')
            row.operator("archipack.custom_draw", text="Draw", icon='GREASEPENCIL')


# ----------------------------------------------------
# ALT + A menu
# ----------------------------------------------------


def draw_menu(self, context):
    global icon_man
    icons = icon_man["main"]
    layout = self.layout
    layout.operator_context = 'INVOKE_REGION_WIN'

    layout.operator("archipack.wall2",
                    text="Wall",
                    icon_value=icons["wall"].icon_id
                    ).mode = 'CREATE'
    layout.operator("archipack.window_preset_menu",
                    text="Window",
                    icon_value=icons["window"].icon_id
                    ).preset_operator = "archipack.window"
    layout.operator("archipack.door_preset_menu",
                    text="Door",
                    icon_value=icons["door"].icon_id
                    ).preset_operator = "archipack.door"
    layout.operator("archipack.stair_preset_menu",
                    text="Stair",
                    icon_value=icons["stair"].icon_id
                    ).preset_operator = "archipack.stair"
    layout.operator("archipack.fence_preset_menu",
                    text="Fence",
                    icon_value=icons["fence"].icon_id
                    ).preset_operator = "archipack.fence"
    layout.operator("archipack.floor_preset_menu",
                    text="Floor",
                    icon_value=icons["floor"].icon_id
                    ).preset_operator = "archipack.floor"
    layout.operator("archipack.molding_preset_menu",
                    text="Molding",
                    icon_value=icons["molding"].icon_id
                    ).preset_operator = "archipack.molding"
    layout.operator("archipack.roof_preset_menu",
                    text="Roof",
                    icon_value=icons["roof"].icon_id
                    ).preset_operator = "archipack.roof"
    layout.operator("archipack.kitchen_preset_menu",
                    text="Kitchen",
                    icon_value=icons["kitchen"].icon_id
                    ).preset_operator = "archipack.kitchen"
    layout.operator("archipack.blind_preset_menu",
                    text="Blind",
                    icon_value=icons["blind"].icon_id
                    ).preset_operator = "archipack.blind"
    layout.operator("archipack.truss",
                    text="Truss",
                    icon_value=icons["truss"].icon_id
                    )
    layout.operator("archipack.dimension_auto",
                    icon_value=icons["dimension_auto"].icon_id
                    ).mode = 'CREATE'
    layout.operator("archipack.layout",
                    icon_value=icons["layout"].icon_id
                    )
    layout.operator("archipack.section",
                    icon_value=icons["section"].icon_id
                    ).mode = 'CREATE'
    layout.operator("archipack.section_camera",
                    text="Section cam",
                    icon='CAMERA_DATA'
                    ).mode = 'CREATE'


class ARCHIPACK_MT_create_menu(Menu):
    bl_label = 'Archipack'
    bl_idname = 'ARCHIPACK_MT_create_menu'

    def draw(self, context):
        draw_menu(self, context)


def menu_func(self, context):
    layout = self.layout
    layout.separator()
    global icon_man
    icons = icon_man["main"]

    # either draw sub menu or right at end of this one
    if context.user_preferences.addons[__name__].preferences.create_submenu:
        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.menu("ARCHIPACK_MT_create_menu", icon_value=icons["archipack"].icon_id)
    else:
        draw_menu(self, context)


def register():
    global icon_man

    glob = globals()
    for sub in submodules:
        module = glob["archipack_{}".format(sub)]
        if hasattr(module, "register"):
            module.register()

    # load icons, icon_man automatically free on unregister
    icon_man = glob["archipack_iconmanager"].icons
    icon_man.load(os.path.join(os.path.dirname(__file__), "icons"), 'main')

    bpy.utils.register_class(Archipack_Pref)
    update_panel(None, bpy.context)

    bpy.utils.register_class(ARCHIPACK_MT_create_menu)
    bpy.types.INFO_MT_mesh_add.append(menu_func)

    addon_updater_ops.register(bl_info)

    # alt+M
    try:
        wm = bpy.context.window_manager
        kc = wm.keyconfigs.addon
        if kc:
            km = kc.keymaps.new(name='3D View', space_type='VIEW_3D')
            # Change here the Type of Event for your own shortcut, default ALT+M
            km.keymap_items.new('archipack.manipulate', 'M', 'PRESS', alt=True)
            # ctrl = True, shift = True)
    except:
        pass


def unregister():
    bpy.types.INFO_MT_mesh_add.remove(menu_func)
    bpy.utils.unregister_class(ARCHIPACK_MT_create_menu)
    bpy.utils.unregister_class(TOOLS_PT_Archipack_PolyLib)
    bpy.utils.unregister_class(TOOLS_PT_Archipack_Tools)
    bpy.utils.unregister_class(TOOLS_PT_Archipack_Create)
    bpy.utils.unregister_class(TOOLS_PT_Archipack_Create_2d)
    bpy.utils.unregister_class(TOOLS_PT_Archipack_Create_Custom)
    bpy.utils.unregister_class(TOOLS_PT_Archipack_About)
    bpy.utils.unregister_class(Archipack_Pref)

    # unregister submodules
    glob = globals()
    for sub in reversed(submodules):
        module = glob["archipack_{}".format(sub)]
        if hasattr(module, "unregister"):
            module.unregister()

    addon_updater_ops.unregister()

    # alt+M
    try:
        wm = bpy.context.window_manager
        kc = wm.keyconfigs.addon
        if kc:
            km = kc.keymaps['3D View']
            for kmi in km.keymap_items:
                if kmi.idname == 'archipack.manipulate':
                    km.keymap_items.remove(kmi)
                    break
    except:
        pass


if __name__ == "__main__":
    register()
