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
import numpy as np
from bpy.types import Operator
from bpy.props import BoolProperty
from mathutils import Vector
from .archipack_object import ArchipackObjectsManager


class ArchipackBoolManager(ArchipackObjectsManager):
    """
        Handle hybrid methods for booleans
        merge holes with boolean and use result on wall
    """
    def __init__(self):
        """
            mode in 'ROBUST', 'INTERACTIVE', 'HYBRID'
        """
        # internal variables
        self.itM = None
        self.min_x = 0
        self.min_y = 0
        self.min_z = 0
        self.max_x = 0
        self.max_y = 0
        self.max_z = 0

    def _world_bounding_box(self, o):
        tM = o.matrix_world
        bpts = np.array([list(tM * Vector(b)) for b in o.bound_box]).transpose()
        minx, miny, minz = bpts.min(axis=1)
        maxx, maxy, maxz = bpts.max(axis=1)
        return minx, miny, minz, maxx, maxy, maxz

    def _init_bounding_box(self, wall):
        self.minx, self.miny, self.minz, \
            self.maxx, self.maxy, self.maxz = self._world_bounding_box(wall)

        self.center = Vector((
            0.5 * (self.maxx + self.minx),
            0.5 * (self.maxy + self.miny),
            0.5 * (self.maxz + self.minz)))

    def _contains(self, o):
        # check for bounding boxes intersections
        minx, miny, minz, maxx, maxy, maxz = self._world_bounding_box(o)
        if (maxx < self.minx or minx > self.maxx or
                maxy < self.miny or miny > self.maxy or
                maxz < self.minz or minz > self.maxz):
            return False

        return True

    def filter_wall(self, wall):
        d = wall.data
        return d is not None and (
               'archipack_wall2' in d or 'archipack_wall' in d)

    def datablock(self, o):
        """
            get datablock from windows and doors
            return
                datablock if found
                None when not found
        """
        d = None
        if o.data is not None:
            if "archipack_window" in o.data:
                d = o.data.archipack_window[0]
            elif "archipack_door" in o.data:
                d = o.data.archipack_door[0]
        return d

    def prepare_hole(self, context, hole):
        hole.draw_type = 'WIRE'
        hole.hide_render = True
        self.select_object(context, hole)
        if "archipack_custom_hole" not in hole:
            hole.lock_location = (True, True, True)
            hole.lock_rotation = (True, True, True)
            hole.lock_scale = (True, True, True)
            hole.hide_select = True

        # hide hole from cycles
        hole.cycles_visibility.camera = False
        hole.cycles_visibility.diffuse = False
        hole.cycles_visibility.glossy = False
        hole.cycles_visibility.shadow = False
        hole.cycles_visibility.scatter = False
        hole.cycles_visibility.transmission = False

    def get_child_hole(self, o):

        # Handle custom holes : objects tagged with "archipack_custom_hole"
        if "archipack_custom_hole" in o:
            return o

        for hole in o.children:
            if "archipack_hole" in hole or "archipack_custom_hole" in hole:
                return hole
        return None

    def _generate_hole(self, context, o):
        # use existing one
        hole = self.get_child_hole(o)
        if hole is not None:
            return hole

        # generate single hole from archipack primitives
        # regardless of location to allow draw tools to show holes
        d = self.datablock(o)
        if d is not None:
            hole = d.interactive_hole(context, o)

        return hole

    def sort_holes(self, wall, holes):
        """
            sort hole from center to borders by distance from center
            may improve nested booleans
        """
        holes = sorted(
            [(o, (o.matrix_world.translation - self.center).length) for o in holes],
            key=lambda x: x[1])

        return [o[0] for o in holes]

    def difference(self, basis, hole):
        m = basis.modifiers.new('AutoBoolean', 'BOOLEAN')
        m.operation = 'DIFFERENCE'
        m.object = hole

    def union(self, basis, hole):
        m = basis.modifiers.new('AutoMerge', 'BOOLEAN')
        m.operation = 'UNION'
        m.object = hole

    def remove_modif_and_object(self, context, o, to_delete):
        for m, h in to_delete:
            if m is not None:
                if m.object is not None:
                    m.object = None
                o.modifiers.remove(m)
            if h is not None:
                if "archipack_custom_hole" not in h:
                    self.delete_object(context, h)

    def create_merge_basis(self, context, wall):
        """
            Create object to merge all holes using booleans
        """
        h = bpy.data.meshes.new("AutoBoolean")
        hole_obj = bpy.data.objects.new("AutoBoolean", h)
        hole_obj['archipack_hybridhole'] = True
        
        self.link_object_to_scene(context, hole_obj)

        if wall.parent is not None:
            hole_obj.parent = wall.parent
        hole_obj.matrix_world = wall.matrix_world.copy()
        for mat in wall.data.materials:
            hole_obj.data.materials.append(mat)
        return hole_obj

    def update_hybrid(self, context, wall, childs, holes):
        """
            Update all holes modifiers
            remove holes not found in childs

            robust -> mixed:
                there is only one object taged with "archipack_robusthole"
            interactive -> mixed:
                many modifisers on wall taged with "archipack_hole"
                keep objects
        """
        existing = []
        to_delete = []

        # remove modifier and holes not found in new list
        self.remove_modif_and_object(context, wall, to_delete)

        m = wall.modifiers.get("AutoMixedBoolean")
        if m is None:
            m = wall.modifiers.new('AutoMixedBoolean', 'BOOLEAN')
            m.operation = 'DIFFERENCE'

        if m.object is None:
            hole_obj = self.create_merge_basis(context, wall)
            m.object = hole_obj
        else:
            hole_obj = m.object

        self.prepare_hole(context, hole_obj)

        to_delete = []

        # mixed-> mixed
        for m in hole_obj.modifiers:
            h = m.object
            if h in holes:
                existing.append(h)
            else:
                to_delete.append([m, h])

        # remove modifier and holes not found in new list
        self.remove_modif_and_object(context, hole_obj, to_delete)

        # add modifier and holes not found in existing
        for h in holes:
            if h not in existing:
                self.union(hole_obj, h)

        # AutoBoolean will be child of reference point
        childs.append(hole_obj)

    def autoboolean(self, context, wall):
        """
            Entry point for multi-boolean operations like
            in T panel autoBoolean
        """
        io = None
        wall2d = None
        if wall.data is not None and "archipack_wall2" in wall.data:
            # ensure wall modifier is there before any boolean
            # to support "revival" of applied modifiers
            m = wall.modifiers.get("Wall")
            wd = wall.data.archipack_wall2[0]
            io, wall2d, childs = wd.as_geom(context, wall, 'BOTH', [], [], [])
            if m is None:
                self.select_object(context, wall, True)
                wd.update(context)

        bpy.ops.object.select_all(action='DESELECT')
        context.scene.objects.active = None
        childs = []
        holes = []
        # get wall bounds to find what's inside
        self._init_bounding_box(wall)

        # either generate hole or get existing one


        for o in self.get_scene_objects(context):
            # filter holes found in wall bounding box
            if o.type == 'MESH' and self._contains(o):

                d = self.datablock(o)
                intersect = True

                # deep check to ensure neighboors walls
                # dosent interfer with this one
                # using a pygeos based 2d check
                if d is not None and wall2d is not None:
                    coords = d.hole_2d('BOUND')
                    hole2d = io.coords_to_polygon(o.matrix_world, coords)
                    intersect = wall2d.intersects(hole2d)

                if intersect:
                    h = self._generate_hole(context, o)
                    if h is not None:
                        holes.append(h)
                        childs.append(o)

        # sort from center to border
        self.sort_holes(wall, holes)

        # hole(s) are selected and active after this one
        for hole in holes:
            # copy wall material to hole
            hole.data.materials.clear()
            for mat in wall.data.materials:
                hole.data.materials.append(mat)

            self.prepare_hole(context, hole)

        # update / remove / add  boolean modifier
        self.update_hybrid(context, wall, childs, holes)

        bpy.ops.object.select_all(action='DESELECT')
        # parenting childs to wall reference point
        """
        if wall.parent is None:
            x, y, z = wall.bound_box[0]
            context.scene.cursor_location = wall.matrix_world * Vector((x, y, z))
            # fix issue #9
            self.select_object(context, wall, True)
            bpy.ops.archipack.reference_point()
        else:
            self.show_object(wall.parent)
            wall.parent.hide_select = False
            self.select_object(context, wall.parent, True)
        """
        self.select_object(context, wall, True)
        for o in childs:
            # parent archipack_custom
            if o.parent and o.data and "archipack_custom_part" in o.data:
                o.parent.hide_select = False
                self.select_object(context, o.parent)
            else:
                o.hide_select = False
                self.select_object(context, o)
        
        bpy.ops.archipack.add_reference_point()
        
        # if bpy.ops.archipack.parent_to_reference.poll():
        #    bpy.ops.archipack.parent_to_reference()

        for o in childs:
            if "archipack_hole" in o or "archipack_hybridhole" in o:
                o.hide_select = True

    def singleboolean(self, context, wall, o):
        """
            Entry point for single boolean operations
            in use in draw door and windows over wall
            o is either a window or a door
        """

        # generate holes for crossing window and doors
        hole = self._generate_hole(context, o)
        hole_obj = None

        if hole is None:
            return

        hole.data.materials.clear()
        for mat in wall.data.materials:
            hole.data.materials.append(mat)

        self.prepare_hole(context, hole)

        # find or add merge basis to wall
        m = wall.modifiers.get('AutoMixedBoolean')

        if m is None:
            m = wall.modifiers.new('AutoMixedBoolean', 'BOOLEAN')
            m.operation = 'DIFFERENCE'

        if m.object is None:
            hole_obj = self.create_merge_basis(context, wall)
            m.object = hole_obj
        else:
            hole_obj = m.object

        # add hole to merge basis
        self.union(hole_obj, hole)

        bpy.ops.object.select_all(action='DESELECT')
        
        self.select_object(context, wall, True)
        self.select_object(context, hole_obj)
        self.select_object(context, o)
        bpy.ops.archipack.add_reference_point()
        
        self.select_object(context, wall, True)

        if "archipack_wall2" in wall.data:
            d = wall.data.archipack_wall2[0]
            d.setup_childs(context, wall)
            if d.dimensions:
                # update manipulators
                d.update(context, manipulable_refresh=True)
            else:
                d.relocate_childs(context, wall)

        if hole_obj is not None:
            self.prepare_hole(context, hole_obj)


class ARCHIPACK_OT_single_boolean(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.single_boolean"
    bl_label = "SingleBoolean"
    bl_description = "Add single boolean for doors and windows"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    """
        Wall must be active object
        window or door must be selected
    """

    @classmethod
    def poll(cls, context):
        w = context.active_object
        return (w.data is not None and
            ("archipack_wall2" in w.data or
            "archipack_wall" in w.data) and
            len(context.selected_objects) == 2
            )

    def draw(self, context):
        pass

    def execute(self, context):
        if context.mode == "OBJECT":
            wall = context.active_object
            manager = ArchipackBoolManager()
            for o in context.selected_objects:
                if o != wall:
                    manager.singleboolean(context, wall, o)
                    break
            self.unselect_object(o)
            self.select_object(context, wall, True)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_auto_boolean(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.auto_boolean"
    bl_label = "AutoBoolean"
    bl_description = "Automatic boolean for doors and windows. Select your wall(s) then push"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.archipack.disable_manipulate()
            manager = ArchipackBoolManager()
            active = context.active_object
            walls = [wall for wall in context.selected_objects if manager.filter_wall(wall)]
            bpy.ops.object.select_all(action='DESELECT')
            for wall in walls:
                manager.autoboolean(context, wall)
                bpy.ops.object.select_all(action='DESELECT')
                self.select_object(context, wall, True)
                
            # reselect walls
            bpy.ops.object.select_all(action='DESELECT')
            for wall in walls:
                self.select_object(context, wall)
            self.select_object(context, active, True)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_generate_hole(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.generate_hole"
    bl_label = "Generate hole"
    bl_description = "Generate interactive hole for doors and windows"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        if context.mode == "OBJECT":
            manager = ArchipackBoolManager()
            o = context.active_object

            d = manager.datablock(o)
            if d is None:
                self.report({'WARNING'}, "Archipack: active object must be a door or a window")
                return {'CANCELLED'}
            bpy.ops.object.select_all(action='DESELECT')
            self.select_object(context, o, True)
            hole = manager._generate_hole(context, o)
            manager.prepare_hole(context, hole)
            self.unselect_object(hole)
            self.select_object(context, o, True)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_custom_hole(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.custom_hole"
    bl_label = "Custom hole"
    bl_description = "Make active object a hole for autoboolean"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    remove = BoolProperty(name="Clear hole parameter", default=False)

    @classmethod
    def poll(self, context):
        o = context.active_object
        return (o is not None and
            o.type == 'MESH' and
            not ArchipackBoolManager.filter_wall(None, o))

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            if self.remove:
                if "archipack_custom_hole" in o:
                    del o["archipack_custom_hole"]
                    o.draw_type = 'TEXTURED'
                    for hole in self.get_scene_objects(context):
                        if "archipack_hybridhole" in hole:
                            for m in hole.modifiers:
                                if m.object == o:
                                    m.object = None
                                    hole.modifiers.remove(m)
            else:
                o["archipack_custom_hole"] = 1
                o.draw_type = 'WIRE'
                manager = ArchipackBoolManager()
                manager._init_bounding_box(o)
                walls = [wall for wall in self.get_scene_objects(context)
                    if manager.filter_wall(wall) and manager._contains(wall)]
                for wall in walls:
                    self.select_object(context, wall, True)
                    bpy.ops.archipack.auto_boolean()
                    self.unselect_object(wall)
                self.select_object(context, o, True)

            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_apply_holes(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.apply_holes"
    bl_label = "Apply holes"
    bl_description = "Apply modifiers and remove holes from scene"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    selected_only = BoolProperty(default=False)

    @classmethod
    def poll(cls, context):
        return context.mode == "OBJECT"

    def modifiers_apply(self, context, o):
        ctx = bpy.context.copy()
        ctx['object'] = o
        for mod in o.modifiers[:]:
            ctx['modifier'] = mod
            try:
                bpy.ops.object.modifier_apply(ctx, apply_as='DATA',
                                              modifier=ctx['modifier'].name)
            except:
                pass

    def get_boolobjects(self, o, to_remove):
        modifiers = [m for m in o.modifiers if m.type == 'BOOLEAN']
        for m in modifiers:
            if m.object is None:
                o.modifiers.remove(m)
            else:
                if 'archipack_custom_hole' not in m.object:
                    to_remove.append(m.object)
                if 'archipack_hybridhole' in m.object:
                    self.get_boolobjects(m.object, to_remove)

    def apply(self, context, objects):
        to_remove = []
        for o in objects:
            self.get_boolobjects(o, to_remove)

        for o in objects:
            if o.data is not None and ("archipack_wall2" in o.data or "archipack_wall" in o.data):
                self.modifiers_apply(context, o)

        bpy.ops.object.select_all(action="DESELECT")
        for r in to_remove:
            r.hide_select = False
            self.select_object(context, r, True)
        bpy.ops.object.delete(use_global=False)

    def execute(self, context):
        if context.mode == "OBJECT":

            if self.selected_only:
                objects = context.selected_objects[:]
            else:
                objects = self.get_scene_objects(context)

            self.apply(context, [o for o in objects if o.type == 'MESH'])

            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


def register():
    bpy.utils.register_class(ARCHIPACK_OT_generate_hole)
    bpy.utils.register_class(ARCHIPACK_OT_single_boolean)
    bpy.utils.register_class(ARCHIPACK_OT_auto_boolean)
    bpy.utils.register_class(ARCHIPACK_OT_custom_hole)
    bpy.utils.register_class(ARCHIPACK_OT_apply_holes)


def unregister():
    bpy.utils.unregister_class(ARCHIPACK_OT_generate_hole)
    bpy.utils.unregister_class(ARCHIPACK_OT_single_boolean)
    bpy.utils.unregister_class(ARCHIPACK_OT_auto_boolean)
    bpy.utils.unregister_class(ARCHIPACK_OT_custom_hole)
    bpy.utils.unregister_class(ARCHIPACK_OT_apply_holes)
