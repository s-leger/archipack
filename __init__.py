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
    'deps': 'shapely',
    'version': (1, 2, 2),
    'blender': (2, 7, 8),
    'location': 'View3D > Tools > Create > Archipack',
    'warning': '',
    'wiki_url': 'https://github.com/s-leger/blenderPolygons/wiki',
    'tracker_url': 'https://github.com/s-leger/blenderPolygons/issues',
    'link': 'https://github.com/s-leger/blenderPolygons',
    'support': 'COMMUNITY',
    'category': '3D View'
    }

import os

if "bpy" in locals():
    import importlib as imp
    imp.reload(archipack_autoboolean)
    imp.reload(archipack_door)
    imp.reload(archipack_window)
    imp.reload(archipack_stair)
    imp.reload(archipack_wall2)
    imp.reload(archipack_fence)
    imp.reload(archipack_wall)
    imp.reload(archipack_rendering)
    try:
        imp.reload(archipack_polylib)
        HAS_POLYLIB = True
    except:
        HAS_POLYLIB = False
        pass

    print("archipack: reload ready")
else:
    from . import archipack_autoboolean
    from . import archipack_door
    from . import archipack_window
    from . import archipack_stair
    from . import archipack_wall
    from . import archipack_wall2
    from . import archipack_fence
    from . import archipack_rendering
    try:
        """
            polylib depends on shapely
            raise ImportError when not meet
        """
        from . import archipack_polylib
        HAS_POLYLIB = True
    except:
        print("archipack: Polylib failed to load, missing shapely ?")
        HAS_POLYLIB = False
        pass

    # from . import archipack_polylib

    print("archipack: ready")

# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import Panel, WindowManager, PropertyGroup
from bpy.props import EnumProperty, PointerProperty
from bpy.utils import previews as iconsLib
icons_dict = {}


class TOOLS_PT_Archipack_PolyLib(Panel):
    bl_label = "Archipack 2d to 3d"
    bl_idname = "TOOLS_PT_Archipack_PolyLib"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Tools"

    @classmethod
    def poll(self, context):

        global archipack_polylib
        return HAS_POLYLIB and ((archipack_polylib.vars_dict['select_polygons'] is not None) or
                (context.object is not None and context.object.type == 'CURVE'))

    def draw(self, context):
        layout = self.layout
        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        row.operator(
            "archipack.polylib_detect",
            icon_value=icons_dict["detect"].icon_id,
            text='Detect'
            ).extend = context.window_manager.archipack_polylib.extend
        row.prop(context.window_manager.archipack_polylib, "extend")
        row = box.row(align=True)
        row.prop(context.window_manager.archipack_polylib, "resolution")
        row = box.row(align=True)
        row.label(text="Polygons")
        row = box.row(align=True)
        row.operator(
            "archipack.polylib_pick_2d_polygons",
            icon_value=icons_dict["selection"].icon_id,
            text='Select'
            ).action = 'select'
        row.operator(
            "archipack.polylib_pick_2d_polygons",
            icon_value=icons_dict["union"].icon_id,
            text='Union'
            ).action = 'union'
        row.operator(
            "archipack.polylib_output_polygons",
            icon_value=icons_dict["polygons"].icon_id,
            text='All')
        row = box.row(align=True)
        row.operator(
            "archipack.polylib_pick_2d_polygons",
            text='Wall',
            icon_value=icons_dict["wall"].icon_id).action = 'wall'
        row.prop(context.window_manager.archipack_polylib, "solidify_thickness")
        row = box.row(align=True)
        row.operator("archipack.polylib_pick_2d_polygons",
            text='Window',
            icon_value=icons_dict["window"].icon_id).action = 'window'
        row.operator("archipack.polylib_pick_2d_polygons",
            text='Door',
            icon_value=icons_dict["door"].icon_id).action = 'door'
        row.operator("archipack.polylib_pick_2d_polygons", text='Rectangle').action = 'rectangle'
        row = box.row(align=True)
        row.label(text="Lines")
        row = box.row(align=True)
        row.operator(
            "archipack.polylib_pick_2d_lines",
            icon_value=icons_dict["selection"].icon_id,
            text='Lines').action = 'select'
        row.operator(
            "archipack.polylib_pick_2d_lines",
            icon_value=icons_dict["union"].icon_id,
            text='Union').action = 'union'
        row.operator(
            "archipack.polylib_output_lines",
            icon_value=icons_dict["polygons"].icon_id,
            text='All')
        # row = layout.row(align=True)
        # box = row.box()
        # row = box.row(align=True)
        # row.operator("archipack.polylib_solidify")
        row = box.row(align=True)
        row.label(text="Points")
        row = box.row(align=True)
        row.operator(
            "archipack.polylib_pick_2d_points",
            icon_value=icons_dict["selection"].icon_id,
            text='Points').action = 'select'
        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        row.operator("archipack.polylib_simplify")
        row.prop(context.window_manager.archipack_polylib, "simplify_tolerance")
        row = box.row(align=True)
        row.prop(context.window_manager.archipack_polylib, "simplify_preserve_topology")
        row = layout.row(align=True)
        box = row.box()
        row = box.row(align=True)
        row.operator("archipack.polylib_offset")
        row = box.row(align=True)
        row.prop(context.window_manager.archipack_polylib, "offset_distance")
        row = box.row(align=True)
        row.prop(context.window_manager.archipack_polylib, "offset_side")
        row = box.row(align=True)
        row.prop(context.window_manager.archipack_polylib, "offset_resolution")
        row = box.row(align=True)
        row.prop(context.window_manager.archipack_polylib, "offset_join_style")
        row = box.row(align=True)
        row.prop(context.window_manager.archipack_polylib, "offset_mitre_limit")


class TOOLS_PT_Archipack_Tools(Panel):
    bl_label = "Archipack"
    bl_idname = "TOOLS_PT_Archipack_Tools"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Tools"

    @classmethod
    def poll(self, context):
        return True

    def draw(self, context):
        wm = context.window_manager
        layout = self.layout
        row = layout.row(align=True)
        box = row.box()
        box.label("Auto boolean")
        row = box.row(align=True)
        row.operator("archipack.auto_boolean", text="Robust", icon='HAND').interactive = False
        row.operator("archipack.auto_boolean", text="Interactive", icon='AUTO').interactive = True
        row = layout.row(align=True)
        box = row.box()
        box.label("Rendering")
        row = box.row(align=True)
        row.prop(wm.archipack, 'render_type', text="")
        row = box.row(align=True)
        row.operator("archipack.render", icon='RENDER_STILL')


class TOOLS_PT_Archipack_Create(Panel):
    bl_label = "Archipack"
    bl_idname = "TOOLS_PT_Archipack_Create"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "Create"

    @classmethod
    def poll(self, context):
        return True

    def draw(self, context):
        layout = self.layout
        row = layout.row(align=True)
        box = row.box()
        box.label("Objects")
        row = box.row(align=True)
        row.operator("archipack.window", icon_value=icons_dict["window"].icon_id).mode = 'CREATE'
        row.operator("archipack.door", icon_value=icons_dict["door"].icon_id).mode = 'CREATE'
        row = box.row(align=True)
        row.operator("archipack.stair", icon_value=icons_dict["stair"].icon_id)
        row = box.row(align=True)
        row.operator("archipack.wall2", icon_value=icons_dict["wall"].icon_id)
        row.operator("archipack.wall2_draw", icon='GREASEPENCIL')
        row = box.row(align=True)
        row.operator("archipack.fence")
        row.operator("archipack.fence_from_curve", icon='CURVE_DATA')


# ----------------------------------------------------
# ALT + A menu
# ----------------------------------------------------

# Define "Archipack" menu
def menu_func(self, context):
    layout = self.layout
    layout.separator()
    layout.operator_context = 'INVOKE_REGION_WIN'
    layout.operator("archipack.wall2",
                    text="Wall",
                    icon_value=icons_dict["wall"].icon_id
                    )
    layout.operator("archipack.window",
                    text="Window",
                    icon_value=icons_dict["window"].icon_id
                    ).mode = 'CREATE'
    layout.operator("archipack.door",
                    text="Door",
                    icon_value=icons_dict["door"].icon_id
                    ).mode = 'CREATE'
    layout.operator("archipack.fence",
                    text="Fence")


# ----------------------------------------------------
# Datablock to store global addon variables
# ----------------------------------------------------


class archipack_data(PropertyGroup):
    render_type = EnumProperty (
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
    global icons_dict
    icons_dict = iconsLib.new()
    icons_dir = os.path.join(os.path.dirname(__file__), "icons")
    for icon in os.listdir(icons_dir):
        name, ext = os.path.splitext(icon)
        icons_dict.load(name, os.path.join(icons_dir, icon), 'IMAGE')
    bpy.types.INFO_MT_mesh_add.append(menu_func)
    bpy.utils.register_class(archipack_data)
    WindowManager.archipack = PointerProperty(type=archipack_data)

    bpy.utils.register_module(__name__)


def unregister():
    bpy.types.INFO_MT_mesh_add.remove(menu_func)
    bpy.utils.unregister_class(archipack_data)
    del WindowManager.archipack
    global icons_dict
    iconsLib.remove(icons_dict)
    bpy.utils.unregister_module(__name__)


if __name__ == "__main__":
    register()
