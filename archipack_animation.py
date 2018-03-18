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
from bpy.app.handlers import persistent
from bpy.types import (
    Object, Operator
    )
from bpy.props import (
    EnumProperty, BoolProperty
    )


# Store animated objects for current session
# key: object name, value: datablock name
animated = {}


@persistent
def archipack_animation_onload(dummy):
    """
      Fill in animated dict on file load
    """
    global animated
    for o in bpy.data.objects:
        if "archipack_animation" in o and o.archipack_animation:
            for key in o.data.keys():
                if "archipack_" in key:
                    animated[o.name] = key
    if len(animated) > 0:
        if archipack_animation_updater not in bpy.app.handlers.frame_change_pre:
            bpy.app.handlers.frame_change_pre.append(archipack_animation_updater)


@persistent
def archipack_animation_onunload(dummy):
    """
        Cleanup animated dict on file unload
    """
    global animated
    animated.clear()
    if archipack_animation_updater in bpy.app.handlers.frame_change_pre:
        bpy.app.handlers.frame_change_pre.remove(archipack_animation_updater)


def archipack_animation_updater(dummy):
    global animated
    if len(animated) > 0:
        context = bpy.context
        scene = context.scene
        act = scene.objects.active
        sel = context.selected_objects[:]
        for name in animated:
            o = scene.objects.get(name)
            if o and o.archipack_animation:
                d = getattr(o.data, animated[name])[0]
                o.select = True
                scene.objects.active = o
                d.update(context)
                o.select = False
        for o in sel:
            o.select = True
        scene.objects.active = act


class ARCHIPACK_OT_animation(Operator):
    bl_idname = "archipack.animation"
    bl_label = "Animation"
    bl_description = "Manage animation support for object parameters (Add / Remove / Clear)"
    bl_options = {'REGISTER', 'UNDO'}
    mode = EnumProperty(
        items=(
            ('ENABLE', 'Enable', 'Enable animation support for selected objects'),
            ('DISABLE', 'Disable', 'Disable animation support for selected objects'),
            ('CLEAR', 'Cleanup', 'Remove animation support for all objects')
            )
        )

    @classmethod
    def poll(cls, context):
        return True

    def enable_animation(self, sel):
        """
          Add animation support on selected objects
        """
        global animated
        for o in sel:
            if o.data:
                for key in o.data.keys():
                    if "archipack_" in key:
                        d = getattr(o.data, key)[0]
                        if hasattr(d, "update"):
                            o.archipack_animation = True
                            animated[o.name] = key

        if len(animated) > 0:
            if archipack_animation_updater not in bpy.app.handlers.frame_change_pre:
                bpy.app.handlers.frame_change_pre.append(archipack_animation_updater)

    def disable_animation(self, sel):
        """
          Remove animation support on selected objects
        """
        global animated
        for o in sel:
            if "archipack_animation" in o:
                o.archipack_animation = False
                if o.name in animated:
                    del animated[o.name]

    def execute(self, context):

        if self.mode == 'ENABLE':
            sel = context.selected_objects
            self.enable_animation(sel)

        elif self.mode == 'DISABLE':
            sel = context.selected_objects
            self.disable_animation(sel)

        elif self.mode == 'CLEAR':
            sel = context.scene.objects
            self.disable_animation(sel)
            archipack_animation_onunload(None)

        return {'FINISHED'}


def register():
    global animated
    animated = {}
    Object.archipack_animation = BoolProperty(name="Archipack animation", default=False)
    bpy.app.handlers.load_pre.append(archipack_animation_onunload)
    bpy.app.handlers.load_post.append(archipack_animation_onload)
    bpy.utils.register_class(ARCHIPACK_OT_animation)


def unregister():
    global animated
    animated.clear()
    bpy.app.handlers.load_pre.remove(archipack_animation_onunload)
    bpy.app.handlers.load_post.remove(archipack_animation_onload)
    bpy.utils.unregister_class(ARCHIPACK_OT_animation)
    del Object.archipack_animation
