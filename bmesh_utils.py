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
import bmesh


class BmeshEdit():
    @staticmethod
    def _start(context, o):
        """
            private, start bmesh editing of active object
        """
        o.select = True
        context.scene.objects.active = o
        bpy.ops.object.mode_set(mode='EDIT')
        bm = bmesh.from_edit_mesh(o.data)
        bm.verts.ensure_lookup_table()
        bm.faces.ensure_lookup_table()
        return bm

    @staticmethod
    def _end(bm, o):
        """
            private, end bmesh editing of active object
        """
        bmesh.update_edit_mesh(o.data, True)
        bpy.ops.object.mode_set(mode='OBJECT')
        bm.free()

    @staticmethod
    def _matids(bm, matids):
        for i, matid in enumerate(matids):
            bm.faces[i].material_index = matid

    @staticmethod
    def _uvs(bm, uvs):
        layer = bm.loops.layers.uv.verify()
        l_i = len(uvs)
        for i, face in enumerate(bm.faces):
            if i > l_i:
                raise RuntimeError("Missing uvs for face {}".format(i))
            l_j = len(uvs[i])
            for j, loop in enumerate(face.loops):
                if j > l_j:
                    raise RuntimeError("Missing uv {} for face {}".format(j, i))
                loop[layer].uv = uvs[i][j]

    @staticmethod
    def _verts(bm, verts):
        for i, v in enumerate(verts):
            bm.verts[i].co = v

    @staticmethod
    def buildmesh(context, o, verts, faces, matids=None, uvs=None, weld=False, clean=False):
        bm = BmeshEdit._start(context, o)
        bm.clear()
        for v in verts:
            bm.verts.new(v)
        bm.verts.ensure_lookup_table()
        for f in faces:
            bm.faces.new([bm.verts[i] for i in f])
        bm.faces.ensure_lookup_table()
        if matids is not None:
            BmeshEdit._matids(bm, matids)
        if uvs is not None:
            BmeshEdit._uvs(bm, uvs)
        if weld:
            bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.001)
        BmeshEdit._end(bm, o)
        if clean:
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_all(action='SELECT')
            bpy.ops.mesh.delete_loose()
            bpy.ops.object.mode_set(mode='OBJECT')

    @staticmethod
    def verts(context, o, verts):
        """
            update vertex position of active object
        """
        bm = BmeshEdit._start(context, o)
        BmeshEdit._verts(bm, verts)
        BmeshEdit._end(bm, o)

    @staticmethod
    def aspect(context, o, matids, uvs):
        """
            update material id and uvmap of active object
        """
        bm = BmeshEdit._start(context, o)
        BmeshEdit._matids(bm, matids)
        BmeshEdit._uvs(bm, uvs)
        BmeshEdit._end(bm, o)
