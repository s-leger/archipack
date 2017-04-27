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
import bpy
from bpy.types import Operator
from bpy.props import BoolProperty
from mathutils import Vector
# noinspection PyUnresolvedReferences
from . import archipack_reference_point
# noinspection PyUnresolvedReferences


class ARCHIPACK_OT_auto_boolean(Operator):
    bl_idname = "archipack.auto_boolean"
    bl_label = "AutoBoolean"
    bl_description = "Automatic boolean for doors and windows"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    interactive = BoolProperty(
        name="interactive",
        default=True,
        description="Make one boolean for each objects"
    )

    # internal variables
    itM = None
    min_x = 0
    min_y = 0
    min_z = 0
    max_x = 0
    max_y = 0
    max_z = 0

    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.prop(self, 'interactive')

    def autoboolean(self, context, wall):
        bpy.ops.object.select_all(action='DESELECT')
        context.scene.objects.active = None
        childs = []
        # get wall bounds to find what's inside
        self._get_bounding_box(wall)
        # generate holes for crossing window and doors
        # hole(s) are selected and active after this one
        self._generate_holes(context, childs)
        # update / remove / add  boolean modifier
        if self.interactive:
            self._update_interactive_boolean_modif(context, wall, childs)
        else:
            self._update_boolean_modif(context, wall, childs)
        # parenting childs to wall reference point
        if wall.parent is None:
            x, y, z = wall.bound_box[0]
            context.scene.cursor_location = wall.matrix_world * Vector((x, y, z))
            bpy.ops.archipack.reference_point()
        else:
            context.scene.objects.active = wall.parent
        bpy.ops.object.select_all(action='DESELECT')
        wall.select = True
        for o in childs:
            o.select = True
        bpy.ops.archipack.parent_to_reference()

    def _get_bounding_box(self, wall):
        self.itM = wall.matrix_world.inverted()
        x, y, z = wall.bound_box[0]
        self.min_x = x
        self.min_y = y
        self.min_z = z
        x, y, z = wall.bound_box[6]
        self.max_x = x
        self.max_y = y
        self.max_z = z
        self.center = Vector((
            self.min_x + 0.5 * (self.max_x - self.min_x),
            self.min_y + 0.5 * (self.max_y - self.min_y),
            self.min_z + 0.5 * (self.max_z - self.min_z)))

    def _contains(self, pt):
        p = self.itM * pt
        return (p.x >= self.min_x and p.x <= self.max_x and
            p.y >= self.min_y and p.y <= self.max_y and
            p.z >= self.min_z and p.z <= self.max_z)

    def _remove_boolean_modif_object(self, context, old):
        if old is not None:
            context.scene.objects.unlink(old)
            bpy.data.objects.remove(old, do_unlink=True)

    def _prepare_hole(self, hole):
        hole.lock_location = (True, True, True)
        hole.lock_rotation = (True, True, True)
        hole.lock_scale = (True, True, True)
        hole.draw_type = 'WIRE'
        hole.hide_render = True
        hole.select = True

    def _generate_holes(self, context, childs):
        # generate holes from archipack primitives
        for o in context.scene.objects:
            if o.data is not None:
                if ('archipack_window' in o.data and
                    (self._contains(o.location) or
                     self._contains(o.matrix_world * Vector((0, 0, 0.5 * o.data.archipack_window[0].z))))):
                    if self.interactive:
                        hole = o.data.archipack_window[0].interactive_hole(context, o)
                    else:
                        hole = o.data.archipack_window[0].robust_hole(context, o.matrix_world)
                    self._prepare_hole(hole)
                    childs.append(o)
                elif ('archipack_door' in o.data and
                    (self._contains(o.location) or
                     self._contains(o.matrix_world * Vector((0, 0, 0.5 * o.data.archipack_door[0].z))))):
                    if self.interactive:
                        hole = o.data.archipack_door[0].interactive_hole(context, o)
                    else:
                        hole = o.data.archipack_door[0].robust_hole(context, o.matrix_world)
                    self._prepare_hole(hole)
                    childs.append(o)

    def _remove_boolean_modif(self, context, obj, modif):
        old = modif.object
        obj.modifiers.remove(modif)
        self._remove_boolean_modif_object(context, old)

    def partition(self, array, begin, end):
        pivot = begin
        for i in range(begin + 1, end + 1):
            if array[i][1] <= array[begin][1]:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
        array[pivot], array[begin] = array[begin], array[pivot]
        return pivot

    def quicksort(self, array, begin=0, end=None):
        if end is None:
            end = len(array) - 1

        def _quicksort(array, begin, end):
            if begin >= end:
                return
            pivot = self.partition(array, begin, end)
            _quicksort(array, begin, pivot - 1)
            _quicksort(array, pivot + 1, end)
        return _quicksort(array, begin, end)

    def _update_interactive_boolean_modif(self, context, wall, childs):
        modifs = [wall.modifiers.get(modif.name) for modif in wall.modifiers if modif.type == 'BOOLEAN']
        n_modifs = len(modifs)
        # sort holes by distance from wall center from closest to farthest
        # as multiple boolean may be a bit more robust like that
        center = wall.matrix_world * self.center
        holes = [(o, (o.matrix_world.translation - center).length) for o in context.selected_objects]
        n_holes = len(context.selected_objects)
        self.quicksort(holes)
        holes = [o[0] for o in holes]
        # remove modifiers
        for i in range(n_modifs, n_holes, -1):
            self._remove_boolean_modif(context, wall, modifs[i - 1])
        # add modifiers
        for i in range(n_holes):
            if i < n_modifs:
                modif = modifs[i]
                if modif.object not in holes:
                    self._remove_boolean_modif_object(context, modif.object)
            else:
                modif = wall.modifiers.new('AutoBoolean', 'BOOLEAN')
                modif.operation = 'DIFFERENCE'
            modif.object = holes[i]

    def _update_boolean_modif(self, context, wall, childs):
        modifs = [wall.modifiers.get(modif.name) for modif in wall.modifiers if modif.type == 'BOOLEAN']
        modif = wall.modifiers.get('AutoBoolean')
        for m in modifs:
            if m != modif:
                self._remove_boolean_modif(context, wall, m)
        if bool(len(context.selected_objects) > 0):
            # more than one hole : join, result becomes context.object
            if len(context.selected_objects) > 1:
                bpy.ops.object.join()
                context.object['archipack_robusthole'] = True
            hole = context.object
            childs.append(hole)
            old = None
            if modif is None:
                modif = wall.modifiers.new('AutoBoolean', 'BOOLEAN')
                modif.operation = 'DIFFERENCE'
            else:
                old = modif.object
                self._remove_boolean_modif_object(context, old)
            modif.object = hole
        elif modif is not None:
            self._remove_boolean_modif(context, wall, modif)

    def filter_wall(self, wall):
        d = wall.data
        return (d is None or
               'archipack_window' in d or
               'archipack_window_panel' in d or
               'archipack_door' in d or
               'archipack_doorpanel' in d or
               'archipack_hole' in wall or
               'archipack_robusthole' in wall or
               'archipack_handle' in wall)

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            active = context.scene.objects.active
            walls = [wall for wall in context.selected_objects if not self.filter_wall(wall)]
            bpy.ops.object.select_all(action='DESELECT')
            for wall in walls:
                self.autoboolean(context, wall)
                wall.select = True
                context.scene.objects.active = wall
                if wall.data is not None and 'archipack_wall2' in wall.data:
                    bpy.ops.archipack.wall2_manipulate('EXEC_DEFAULT')
            # reselect walls
            bpy.ops.object.select_all(action='DESELECT')
            for wall in walls:
                wall.select = True
            context.scene.objects.active = active
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


bpy.utils.register_class(ARCHIPACK_OT_auto_boolean)
