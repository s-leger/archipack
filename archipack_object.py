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
# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.props import BoolProperty, StringProperty
from .materialutils import MaterialUtils


class ArchipackObject():
    """
        Shared property of archipack's objects PropertyGroup
        provide basic support for copy to selected
        and datablock access / filtering by object
    """

    def iskindof(self, o, typ):
        """
            return true if object contains databloc of typ name
        """
        return o.data is not None and typ in o.data

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
            return cls.__name__ in o.data
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
            return getattr(o.data, cls.__name__)[0]
        except:
            pass
        return None

    def find_in_selection(self, context, auto_update=True):
        """
            find witch selected object this datablock instance belongs to
            store context to be able to restore after oops
            provide support for "copy to selected"
            return
            object or None when instance not found in selected objects
        """
        if not auto_update:
            return None
            
        active = context.active_object
        selected = [o for o in context.selected_objects]
        
        for o in selected:
            if self.__class__.datablock(o) == self:
                self.previously_selected = selected
                self.previously_active = active
                return o
        
        return None
    
    def restore_context(self, context):
        # restore context
        bpy.ops.object.select_all(action="DESELECT")
        
        try:
            for o in self.previously_selected:
                o.select = True
        except:
            pass

        self.previously_active.select = True
        context.scene.objects.active = self.previously_active
        self.previously_selected = None
        self.previously_active = None
        
        
class ArchipackCreateTool():
    """
        Shared property of archipack's create tool Operator
    """
    auto_manipulate = BoolProperty(
            name="Auto manipulate",
            description="Enable object's manipulators after create",
            options={'SKIP_SAVE'},
            default=True
            )
    filepath = StringProperty(
            options={'SKIP_SAVE'},
            name="Preset",
            description="Full filename of python preset to load at create time",
            default=""
            )

    @property
    def archipack_category(self):
        """
            return target object name from ARCHIPACK_OT_object_name
        """
        return self.bl_idname[13:]

    def load_preset(self, d):
        """
            Load python preset
            preset: full filename.py with path
        """
        d.auto_update = False
        if self.filepath != "":
            try:
                bpy.ops.script.python_file_run(filepath=self.filepath)
            except:
                print("Archipack unable to load preset file : %s" % (self.filepath))
                pass
        d.auto_update = True

    def add_material(self, o):
        try:
            getattr(MaterialUtils, "add_" + self.archipack_category + "_materials")(o)
        except:
            print("Archipack MaterialUtils.add_%s_materials not found" % (self.archipack_category))
            pass

    def manipulate(self):
        if self.auto_manipulate:
            try:
                op = getattr(bpy.ops.archipack, self.archipack_category + "_manipulate")
                if op.poll():
                    op('INVOKE_DEFAULT')
            except:
                print("Archipack bpy.ops.archipack.%s_manipulate not found" % (self.archipack_category))
                pass

              
"""
d = archipack_window.datablock(o)
archipack_window.filter(o)
"""
