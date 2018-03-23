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
    'version': (1, 3, 8),
    'blender': (2, 7, 8),
    'location': 'View3D > Tools > Create > Archipack',
    'warning': '',
    'wiki_url': 'https://s-leger.github.io/archipack/index.html',
    'tracker_url': 'https://github.com/s-leger/archipack/issues',
    'link': 'https://github.com/s-leger/archipack',
    'support': 'COMMUNITY',
    'category': 'Add Mesh'
    }

import os

if "bpy" in locals():
    import importlib as imp
    imp.reload(archipack_progressbar)
    imp.reload(archipack_material)
    imp.reload(archipack_snap)
    imp.reload(archipack_manipulator)
    imp.reload(archipack_curveman)
    imp.reload(archipack_dimension)
    imp.reload(archipack_reference_point)
    imp.reload(archipack_autoboolean)
    imp.reload(archipack_door)
    imp.reload(archipack_window)
    imp.reload(archipack_stair)
    imp.reload(archipack_wall)
    imp.reload(archipack_wall2)
    imp.reload(archipack_slab)
    imp.reload(archipack_roof)
    imp.reload(archipack_fence)
    imp.reload(archipack_truss)
    # imp.reload(archipack_toolkit)
    imp.reload(archipack_custom)
    imp.reload(archipack_floor)
    imp.reload(archipack_floor_heating)
    imp.reload(archipack_blind)
    imp.reload(archipack_kitchen)
    imp.reload(archipack_molding)
    imp.reload(archipack_rendering)
    imp.reload(archipack_section)
    imp.reload(archipack_animation)
    # imp.reload(archipack_envi)
    imp.reload(archipack_io)
    imp.reload(archipack_2d_layout)
    imp.reload(archipack_io_export_svg)
    imp.reload(archipack_polylines)
    imp.reload(addon_updater_ops)
    # imp.reload(archipack_i18n)

    print("archipack: reload ready")
else:
    from . import archipack_progressbar
    from . import archipack_material
    from . import archipack_snap
    from . import archipack_manipulator
    from . import archipack_curveman
    from . import archipack_dimension
    from . import archipack_reference_point
    from . import archipack_autoboolean
    from . import archipack_door
    from . import archipack_window
    from . import archipack_stair
    from . import archipack_wall
    from . import archipack_wall2
    from . import archipack_slab
    from . import archipack_roof
    from . import archipack_fence
    from . import archipack_truss
    from . import archipack_custom
    # from . import archipack_toolkit
    from . import archipack_floor
    from . import archipack_floor_heating
    from . import archipack_blind
    from . import archipack_kitchen
    from . import archipack_molding
    from . import archipack_rendering
    from . import archipack_section
    from . import archipack_animation
    # from . import archipack_envi
    from . import archipack_io
    from . import archipack_2d_layout
    from . import archipack_io_export_svg
    from . import archipack_polylines
    from . import addon_updater_ops
    # from . import archipack_i18n

    print("archipack: ready")

# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import (
    Panel, WindowManager, PropertyGroup,
    AddonPreferences, Menu
    )
from bpy.props import (
    EnumProperty, PointerProperty,
    StringProperty, BoolProperty,
    IntProperty, FloatProperty, FloatVectorProperty
    )

from bpy.utils import previews
icons_collection = {}


# ----------------------------------------------------
# Addon preferences
# ----------------------------------------------------

def update_panel(self, context):
    try:
        bpy.utils.unregister_class(TOOLS_PT_Archipack_PolyLib)
        bpy.utils.unregister_class(TOOLS_PT_Archipack_Tools)
        bpy.utils.unregister_class(TOOLS_PT_Archipack_Create)
    except:
        pass
    prefs = context.user_preferences.addons[__name__].preferences
    TOOLS_PT_Archipack_PolyLib.bl_category = prefs.tools_category
    bpy.utils.register_class(TOOLS_PT_Archipack_PolyLib)
    TOOLS_PT_Archipack_Tools.bl_category = prefs.tools_category
    bpy.utils.register_class(TOOLS_PT_Archipack_Tools)
    TOOLS_PT_Archipack_Create.bl_category = prefs.create_category
    bpy.utils.register_class(TOOLS_PT_Archipack_Create)


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
    matlib_path = StringProperty(
            name="Folder path",
            subtype="DIR_PATH",
            description="absolute path to material library folder",
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
        col = split.column()
        col.label(text="Font size:")
        col.prop(self, "feedback_size_main")
        col.prop(self, "feedback_size_title")
        col.prop(self, "feedback_size_shortcut")
        col.label(text="Manipulators:")
        col.prop(self, "arrow_size")
        col.prop(self, "handle_size")
        addon_updater_ops.update_settings_ui(self, context)


# ----------------------------------------------------
# Archipack panels
# ----------------------------------------------------


class TOOLS_PT_Archipack_PolyLib(Panel):
    bl_label = "Archipack 2d to 3d"
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
        global icons_collection
        icons = icons_collection["main"]
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
    bl_label = "Archipack Tools"
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
        wm = context.window_manager
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

        box = layout.box()
        box.label("Rendering")
        box.prop(wm.archipack, 'render_type', text="")
        box.operator("archipack.render", icon='RENDER_STILL')
        
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


class TOOLS_PT_Archipack_Create(Panel):
    bl_label = "Add Archipack"
    bl_idname = "TOOLS_PT_Archipack_Create"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Create"
    bl_context = "objectmode"

    @classmethod
    def poll(self, context):
        return True

    def draw(self, context):
        global icons_collection
        prefs = context.user_preferences.addons[__name__].preferences
        addon_updater_ops.check_for_update_background(context)

        icons = icons_collection["main"]
        layout = self.layout
        box = layout.box()
        box.label("Objects")
        row = box.row(align=True)
        row.operator("archipack.wall2",
                    icon_value=icons["wall"].icon_id
                    ).mode = 'CREATE'
        row.operator("archipack.wall2_draw", text="Draw", icon='GREASEPENCIL')
        row.operator("archipack.wall2_from_curve", text="", icon='CURVE_DATA')
        row = box.row(align=True)
        # col = row.column()
        # subrow = col.row(align=True)
        row.operator("archipack.window_preset_menu",
                    text="Window",
                    icon_value=icons["window"].icon_id
                    ).preset_operator = "archipack.window"
        row.operator("archipack.window_preset_menu",
                    text="Draw",
                    icon='GREASEPENCIL'
                    ).preset_operator = "archipack.window_draw"
        # col = row.column()
        # subrow = col.row(align=True)
        row = box.row(align=True)
        row.operator("archipack.door_preset_menu",
                    text="Door",
                    icon_value=icons["door"].icon_id
                    ).preset_operator = "archipack.door"
        row.operator("archipack.door_preset_menu",
                    text="Draw",
                    icon='GREASEPENCIL'
                    ).preset_operator = "archipack.door_draw"
        box.operator("archipack.stair_preset_menu",
                    text="Stair",
                    icon_value=icons["stair"].icon_id
                    ).preset_operator = "archipack.stair"
        row = box.row(align=True)
        row.operator("archipack.fence_preset_menu",
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
        row.operator("archipack.floor_preset_menu",
                    text="Floor",
                    icon_value=icons["floor"].icon_id
                    ).preset_operator = "archipack.floor"
        row.operator("archipack.floor_preset_from_wall",
                    text="->Floor",
                    icon_value=icons["floor_from_wall"].icon_id
                    )
                    # .preset_operator = "archipack.floor_from_wall"
        row.operator("archipack.floor_preset_from_curve",
                    text="",
                    icon='CURVE_DATA')
                    # .preset_operator = "archipack.floor_from_curve"
        row = box.row(align=True)
        row.operator("archipack.molding_preset_menu",
                    text="Molding",
                    icon_value=icons["molding"].icon_id
                    ).preset_operator = "archipack.molding"
        row.operator("archipack.molding_preset_from_wall",
                    text="->Molding",
                    icon_value=icons["molding_from_wall"].icon_id
                    )
                    # .preset_operator = "archipack.molding_from_wall"
        row.operator("archipack.molding_preset_from_curve",
                    text="",
                    icon='CURVE_DATA'
                    )
                    # .preset_operator = "archipack.molding_from_curve"

        row = box.row(align=True)
        row.operator("archipack.roof_preset_menu",
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
        box.operator("archipack.kitchen_preset_menu",
                    text="Kitchen",
                    icon_value=icons["kitchen"].icon_id
                    ).preset_operator = "archipack.kitchen"

        box.operator("archipack.blind_preset_menu",
                    text="Blind",
                    icon_value=icons["blind"].icon_id
                    ).preset_operator = "archipack.blind"
        box.operator("archipack.truss",
                    icon_value=icons["truss"].icon_id
                    )

        box = layout.box()
        box.label("2d Objects")
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
            row = box.row(align=True)
            row.operator("archipack.make_custom", text="Custom", icon='MONKEY')
            row.operator("archipack.custom_draw", text="Draw", icon='GREASEPENCIL')
                
        box = layout.box()
        box.label(text="Custom objects")
        box.operator("archipack.wall", text="Custom wall")

        row = box.row(align=True)
        row.operator("archipack.custom_hole", text="Custom hole").remove = False
        row.operator("archipack.custom_hole", text="", icon='X').remove = True

        addon_updater_ops.update_notice_box_ui(self, context)


# ----------------------------------------------------
# ALT + A menu
# ----------------------------------------------------


def draw_menu(self, context):
    global icons_collection
    icons = icons_collection["main"]
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
    global icons_collection
    icons = icons_collection["main"]

    # either draw sub menu or right at end of this one
    if context.user_preferences.addons[__name__].preferences.create_submenu:
        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.menu("ARCHIPACK_MT_create_menu", icon_value=icons["archipack"].icon_id)
    else:
        draw_menu(self, context)


# ----------------------------------------------------
# Datablock to store global addon variables
# ----------------------------------------------------


class archipack_data(PropertyGroup):
    render_type = EnumProperty(
        items=(
            ('1', "Draw over", "Draw over last rendered image"),
            ('2', "OpenGL", ""),
            ('3', "Animation OpenGL", ""),
            ('4', "Image", "Render image and draw over"),
            ('5', "Animation", "Draw on each frame")
            ),
        name="Render type",
        description="Render method"
        )


def register():
    global icons_collection
    icons = previews.new()
    icons_dir = os.path.join(os.path.dirname(__file__), "icons")
    for icon in os.listdir(icons_dir):
        name, ext = os.path.splitext(icon)
        icons.load(name, os.path.join(icons_dir, icon), 'IMAGE')
    icons_collection["main"] = icons

    archipack_progressbar.register()
    archipack_material.register()
    archipack_snap.register()
    archipack_manipulator.register()
    archipack_curveman.register()
    archipack_dimension.register()
    archipack_reference_point.register()
    archipack_autoboolean.register()
    archipack_door.register()
    archipack_window.register()
    archipack_stair.register()
    archipack_wall.register()
    archipack_wall2.register()
    archipack_roof.register()
    archipack_slab.register()
    archipack_fence.register()
    archipack_truss.register()
    archipack_blind.register()
    archipack_kitchen.register()
    archipack_molding.register()
    archipack_custom.register()
    # archipack_toolkit.register()
    archipack_floor.register()
    archipack_floor_heating.register()
    archipack_rendering.register()
    archipack_section.register()
    archipack_animation.register()
    archipack_io.register()
    archipack_2d_layout.register()
    archipack_io_export_svg.register()
    archipack_polylines.register()

    bpy.utils.register_class(archipack_data)
    WindowManager.archipack = PointerProperty(type=archipack_data)
    bpy.utils.register_class(Archipack_Pref)
    update_panel(None, bpy.context)
    bpy.utils.register_class(ARCHIPACK_MT_create_menu)
    bpy.types.INFO_MT_mesh_add.append(menu_func)

    addon_updater_ops.register(bl_info)


def unregister():
    global icons_collection
    bpy.types.INFO_MT_mesh_add.remove(menu_func)
    bpy.utils.unregister_class(ARCHIPACK_MT_create_menu)

    bpy.utils.unregister_class(TOOLS_PT_Archipack_PolyLib)
    bpy.utils.unregister_class(TOOLS_PT_Archipack_Tools)
    bpy.utils.unregister_class(TOOLS_PT_Archipack_Create)
    bpy.utils.unregister_class(Archipack_Pref)
    # unregister subs
    archipack_progressbar.unregister()
    archipack_material.unregister()
    archipack_snap.unregister()
    archipack_autoboolean.unregister()
    archipack_reference_point.unregister()
    archipack_door.unregister()
    archipack_window.unregister()
    archipack_stair.unregister()
    archipack_wall.unregister()
    archipack_wall2.unregister()
    archipack_roof.unregister()
    archipack_slab.unregister()
    archipack_fence.unregister()
    archipack_truss.unregister()
    archipack_blind.unregister()
    archipack_kitchen.unregister()
    archipack_molding.unregister()
    # archipack_toolkit.unregister()
    archipack_custom.unregister()
    archipack_floor.unregister()
    archipack_floor_heating.unregister()
    archipack_rendering.unregister()
    archipack_section.unregister()
    archipack_animation.unregister()
    archipack_io.unregister()
    archipack_2d_layout.unregister()
    archipack_io_export_svg.unregister()
    archipack_polylines.unregister()
    archipack_manipulator.unregister()
    archipack_dimension.unregister()
    archipack_curveman.unregister()
    bpy.utils.unregister_class(archipack_data)
    del WindowManager.archipack
    # print("Archipack collection.remove start")
    for icons in icons_collection.values():
        previews.remove(icons)
    icons_collection.clear()
    # print("Archipack collection.remove success")
    addon_updater_ops.unregister()


if __name__ == "__main__":
    register()
