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
import os
import bpy
import zipfile
from bpy.types import Operator
from bpy.props import StringProperty
from bpy_extras.io_utils import ImportHelper, ExportHelper


class ARCHIPACK_IO_Import_Preset(Operator, ImportHelper):
    """Import archipack preset file in the apk format (.apk)"""
    bl_idname = "archipack.import_presets"
    bl_label = "Import archipack presets"
    bl_options = {"UNDO"}

    # ImportHelper mixin class uses this
    filename_ext = ".apk"

    filter_glob = StringProperty(
        default="*.apk",
        options={"HIDDEN"},
    )

    def execute(self, context):
        outdir = bpy.utils.user_resource('SCRIPTS',
                          "presets",
                          create=True)
        with zipfile.ZipFile(self.filepath, "r") as z:
            z.extractall(outdir)

        return {"FINISHED"}


class ARCHIPACK_IO_Export_Preset(Operator, ExportHelper):
    """Export archipack preset file in the apk format (.apk)"""
    bl_idname = "archipack.export_presets"
    bl_label = "Export archipack presets"
    bl_options = {"UNDO"}

    # ImportHelper mixin class uses this
    filename_ext = ".apk"

    filter_glob = StringProperty(
        default="*.apk",
        options={"HIDDEN"},
    )

    def make_zipfile(self, output_filename, source_dir):
        relroot = os.path.abspath(source_dir)
        with zipfile.ZipFile(output_filename, "w") as zip:
            for root, dirs, files in os.walk(source_dir):
                if "archipack" in root:
                    zip.write(root, os.path.relpath(root, relroot))
                    for file in files:
                        filename = os.path.join(root, file)
                        if os.path.isfile(filename):
                            arcname = os.path.join(os.path.relpath(root, relroot), file)
                            zip.write(filename, arcname)

    def execute(self, context):
        srcdir = bpy.utils.user_resource('SCRIPTS',
                          "presets",
                          create=False)

        self.make_zipfile(self.filepath, srcdir)

        return {"FINISHED"}


def menu_func_export(self, context):
    self.layout.operator(ARCHIPACK_IO_Export_Preset.bl_idname, text="Archipack preset APK (.apk)")


def menu_func_import(self, context):
    self.layout.operator(ARCHIPACK_IO_Import_Preset.bl_idname, text="Archipack preset APK (.apk)")


def register():
    bpy.utils.register_class(ARCHIPACK_IO_Import_Preset)
    bpy.types.INFO_MT_file_import.append(menu_func_import)
    bpy.utils.register_class(ARCHIPACK_IO_Export_Preset)
    bpy.types.INFO_MT_file_export.append(menu_func_export)


def unregister():
    bpy.utils.unregister_class(ARCHIPACK_IO_Import_Preset)
    bpy.types.INFO_MT_file_import.remove(menu_func_import)
    bpy.utils.unregister_class(ARCHIPACK_IO_Export_Preset)
    bpy.types.INFO_MT_file_import.remove(menu_func_export)
