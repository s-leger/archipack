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
from bl_operators.presets import AddPresetBase


class ArchipackPreset(AddPresetBase):

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o is not None and \
            o.data is not None

    @property
    def preset_subdir(self):
        return self.datablock_name

    @property
    def blacklist(self):
        """
            properties black list for presets
            may override on addon basis
        """
        return []

    @property
    def preset_values(self):
        blacklist = self.blacklist
        blacklist.extend(bpy.types.Mesh.bl_rna.properties.keys())
        d = getattr(bpy.context.active_object.data, self.datablock_name)[0]
        props = d.rna_type.bl_rna.properties.items()
        ret = []
        for prop_id, prop in props:
            if prop_id not in blacklist:
                if not (prop.is_hidden or prop.is_skip_save):
                    ret.append("d.%s" % prop_id)
        return ret

    @property
    def preset_defines(self):
        return ["d = bpy.context.active_object.data." + self.datablock_name + "[0]"]

    def pre_cb(self, context):
        getattr(bpy.context.active_object.data, self.datablock_name)[0].auto_update = False

    def post_cb(self, context):
        getattr(bpy.context.active_object.data, self.datablock_name)[0].auto_update = True
