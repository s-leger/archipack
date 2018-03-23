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
from mathutils import Vector, Matrix, Quaternion


class ViewManager():
    """
      A class to set view so it is able to perform knife project
    """
    def __init__(self, context):
        self.view3d = context.space_data
        self.scene = context.scene
        self.distance = 0
        self.location = Vector()
        self.rotation = Quaternion()
        self.perspective = ""
        self.lens = 0
        self.clip_start = 0
        self.clip_end = 0
        self.camera_name = ""
        self.camera_type = ""
        self.lock_object_name = ""

    def save(self):
        view3d = self.view3d
        region3d = view3d.region_3d
        self.distance = region3d.view_distance
        self.location = region3d.view_location.copy()
        self.rotation = region3d.view_rotation.copy()
        self.perspective = region3d.view_perspective
        self.lens = view3d.lens
        self.clip_start = view3d.clip_start
        self.clip_end = view3d.clip_end
        if region3d.view_perspective == 'CAMERA':
            self.camera_type = view3d.camera.type
            self.camera_name = view3d.camera.name
        if view3d.lock_object is not None:
            self.lock_object_name = view3d.lock_object.name
        else:
            self.lock_object_name = ""
    
    def safe_ortho_view(self, radius, location):
        """
         Set view projection so it is able to handle knife project
         Knife project work best in ortho oblique view
         Center view on target, set distance and
         optimize view clip so whole target is visible
        """
        distance = 2 * radius
        bound = Vector((radius, -radius, 0))
        v = Vector((1, -1, 1)).normalized()
        z = -v
        x = z.cross(Vector((0, 0, 1)))
        y = x.cross(z)
        loc = location + distance * v
        # Inverse view matrix
        itM = Matrix([
            [x.x, y.x, z.x, loc.x],
            [x.y, y.y, z.y, loc.y],
            [x.z, y.z, z.z, loc.z],
            [0, 0, 0, 1]]).inverted()
        view3d = self.view3d
        view3d.lens = 35
        start = abs((itM * -bound).z)
        end = abs((itM * bound).z)
        if start > end:
            end, start = start, end
        view3d.clip_start = 0.1
        view3d.clip_end = 20 * end
        region3d = view3d.region_3d
        region3d.view_perspective = 'ORTHO'
        region3d.view_rotation = Quaternion((
            0.8204732537269592,
            0.42470818758010864,
            0.17591983079910278,
            0.33985117077827454
            ))
        region3d.view_distance = distance
        region3d.view_location = location
        region3d.update()
    
    def safe_knife_project(self, radius, location):
        """
         Set view projection so it is able to handle knife project
         Knife project work best in ortho oblique view
         Center view on target, set distance and
         optimize view clip so whole target is visible
        """
        distance = 2 * radius
        bound = Vector((radius, -radius, 0))
        v = Vector((1, -1, 1)).normalized()
        z = -v
        x = z.cross(Vector((0, 0, 1)))
        y = x.cross(z)
        loc = location + distance * v
        # Inverse view matrix
        itM = Matrix([
            [x.x, y.x, z.x, loc.x],
            [x.y, y.y, z.y, loc.y],
            [x.z, y.z, z.z, loc.z],
            [0, 0, 0, 1]]).inverted()
        view3d = self.view3d
        view3d.lens = 35
        start = abs((itM * -bound).z)
        end = abs((itM * bound).z)
        if start > end:
            end, start = start, end
        view3d.clip_start = max(0, start - 1)
        view3d.clip_end = 1 + end
        region3d = view3d.region_3d
        region3d.view_perspective = 'ORTHO'
        region3d.view_rotation = Quaternion((
            0.8204732537269592,
            0.42470818758010864,
            0.17591983079910278,
            0.33985117077827454
            ))
        region3d.view_distance = distance
        region3d.view_location = location
        region3d.update()

    def restore(self):
        view3d = self.view3d
        region3d = view3d.region_3d
        region3d.view_distance = self.distance
        region3d.view_location = self.location
        region3d.view_rotation = self.rotation
        region3d.view_perspective = self.perspective
        view3d.lens = self.lens
        view3d.clip_start = self.clip_start
        view3d.clip_end = self.clip_end
        if self.perspective == "CAMERA":
            lock_obj = self._get_object(self.lock_object_name)
            if lock_obj:
                view3d.lock_object = lock_obj
            else:
                cam = self._get_object(self.camera_name)
                if cam:
                    view3d.camera = cam
        region3d.update()

    def _get_object(self, name):
        return bpy.data.objects.get(name)
