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
# ----------------------------------------------------------
import bpy
import time
from bpy.app.handlers import persistent
from bpy.types import Operator
from mathutils import Vector
from .archipack_object import ArchipackObjectsManager
from .archipack_gl import GlText


class ThrottleObject(ArchipackObjectsManager):
    """
     Represent an object waiting for update
    """
    def __init__(self, datablock, duration, update_func):
        d = None
        try:
            d = datablock.__class__.datablock
        except:
            pass
        self.update_func = update_func
        self.datablock = d
        self.timeout = time.time() + duration
        self.waiting = True

    def update(self, context, name, now):
        if now < self.timeout:
            return False

        o = self.get_scene_object(context, name)
        act = context.active_object
        try:
            selected = self.is_selected(o)
            d = self.datablock(o)
            if d:
                self.select_object(context, o, True)
                self.waiting = False
                getattr(d, self.update_func)(context)
                if not selected:
                    self.unselect_object(o)

        except:
            pass
        self.select_object(context, act, True)
        return True


class ThrottleHandler():
    """
     Throtteling core
     Provide throtteling update ability for complex objects
     Run throttle modal and update objects when needed

     @TODO:
     allow different delay on per object basis
     exit modal when stack is empty
     modal MUST prevent undo

     Objects throttling implementation:

     from .archipack_throttle import throttle

     - must call throttle.add(context, o, self)
        on top of .update

     - must use throttle.active state
        to provide degenerated geometry on .update
    """
    def __init__(self):
        self.stack = {}
        self._active = False

    def add(self, context, o, d, duration=-1, update_func="update"):

        if o is None or d is None:
            return

        addon_name = __name__.split('.')[0]
        prefs = context.user_preferences.addons[addon_name].preferences

        if not prefs.throttle_enable:
            self.clear()
            return

        if duration < 0:
            duration = prefs.throttle_delay

        start = False

        if len(self.stack) < 1 and not self._active:
            start = True
            self._active = True

        if o.name not in self.stack:
            self.stack[o.name] = ThrottleObject(d, duration, update_func)
        else:
            tim = time.time() + duration
            if tim > self.stack[o.name].timeout:
                self.stack[o.name].timeout = tim
            
        if start:
            bpy.ops.archipack.throttle_update('INVOKE_DEFAULT')

    def is_active(self, name):
        return (self._active and
                name in self.stack and
                self.stack[name].waiting)

    def update(self, context):

        now = time.time()

        stack = list(self.stack.items())
        for name, u in stack:
            if u.update(context, name, now):
                if name in self.stack:
                    del self.stack[name]

        if len(self.stack) < 1:
            self._active = False
            return True

        return False

    def clear(self):
        self.active = False
        self.stack.clear()


throttle = ThrottleHandler()


"""
from archipack import archipack_throttle
archipack_throttle.register()
throttle = archipack_throttle.throttle
throttle.add(C, C.object, None)
"""


@persistent
def cleanup(dummy):
    global throttle
    throttle.clear()


class ARCHIPACK_OT_throttle_update(Operator):
    bl_idname = "archipack.throttle_update"
    bl_label = "Update objects with a delay"

    def draw_handler(self, _self, context):
        global throttle
        label = "Archipack {} objects waiting for update".format(len(throttle.stack))
        self.text.label = label
        self.text.draw(context)

    def modal(self, context, event):
        global throttle
        context.area.tag_redraw()
        if event.type == 'TIMER' and throttle.update(context):
            if self._timer is not None:
                context.window_manager.event_timer_remove(self._timer)
                self._timer = None
            if self._handle is not None:
                bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
                self._handle = None
            return {'FINISHED'}
        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        self.text = GlText(label="Throttle", colour=(0, 1, 0, 1), d=2, font_size=12)
        self.text.set_pos(context, None, Vector((20, 80)), Vector((1, 0, 0)))
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(
                self.draw_handler,
                args,
                'WINDOW',
                'POST_PIXEL')
        self._timer = context.window_manager.event_timer_add(0.1, context.window)
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}


def register():
    bpy.utils.register_class(ARCHIPACK_OT_throttle_update)
    bpy.app.handlers.load_pre.append(cleanup)


def unregister():
    global throttle
    throttle.clear()
    bpy.utils.unregister_class(ARCHIPACK_OT_throttle_update)
    bpy.app.handlers.load_pre.remove(cleanup)
