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
import time
import logging
logger = logging.getLogger("archipack_bmesh")
import bpy
import bmesh
from .archipack_object import objman


class BmeshEdit():
    
    @staticmethod
    def ensure_bmesh(bm):
        bm.verts.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        bm.faces.ensure_lookup_table()
        
    @staticmethod
    def _start(context, o):
        """
            private, start bmesh editing of active object
        """
        objman.select_object(context, o, True)
        # bpy.ops.object.mode_set(mode='EDIT')
        # bm = bmesh.from_edit_mesh(o.data)
        # BmeshEdit.ensure_bmesh(bm)
        bm = bmesh.new()
        bm.from_mesh(o.data)
        return bm
    
    @staticmethod
    def bmesh_join(context, o, list_of_bmeshes, normal_update=False):
        """
            takes as input a list of bm references and outputs a single merged bmesh
            allows an additional 'normal_update=True' to force _normal_ calculations.
        """
        bm = BmeshEdit._start(context, o)
        
        add_vert = bm.verts.new
        add_face = bm.faces.new
        add_edge = bm.edges.new

        for bm_to_add in list_of_bmeshes:
            offset = len(bm.verts)
            for v in bm_to_add.verts:
                add_vert(v.co)

            bm.verts.index_update()
            BmeshEdit.ensure_bmesh(bm)
            
            if bm_to_add.faces:
                layer = bm_to_add.loops.layers.uv.verify()
                dest = bm.loops.layers.uv.verify()
                for face in bm_to_add.faces:
                    f = add_face(tuple(bm.verts[i.index + offset] for i in face.verts))
                    f.material_index = face.material_index
                    for j, loop in enumerate(face.loops):
                        f.loops[j][dest].uv = loop[layer].uv
                bm.faces.index_update()

            if bm_to_add.edges:
                for edge in bm_to_add.edges:
                    edge_seq = tuple(bm.verts[i.index + offset] for i in edge.verts)
                    try:
                        add_edge(edge_seq)
                    except ValueError:
                        # edge exists!
                        pass
                bm.edges.index_update()

        # cleanup
        for old_bm in list_of_bmeshes:
            old_bm.free()

        if normal_update:
            bm.normal_update()

        BmeshEdit._end(bm, o)
    
    @staticmethod
    def _end(bm, o):
        """
            private, end bmesh editing of active object
        """
        bm.normal_update()
        bm.to_mesh(o.data)
        # bmesh.update_edit_mesh(o.data, True)
        # bpy.ops.object.mode_set(mode='OBJECT')
        # bm.free()

    @staticmethod
    def _matids(bm, matids):
        _faces = bm.faces
        for i, matid in enumerate(matids):
            _faces[i].material_index = matid

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
    def emptymesh(context, o):
        bm = BmeshEdit._start(context, o)
        bm.clear()
        BmeshEdit.ensure_bmesh(bm)
        BmeshEdit._end(bm, o)
        
    @staticmethod
    def buildmesh(context, o, verts, faces,
            matids=None, uvs=None, weld=False,
            clean=False, auto_smooth=True, temporary=False):
        
        tim = time.time()
        
        if o is not None:
            # ensure object is visible 
            # otherwhise it is not editable
            vis_state = objman.is_visible(o)
            objman.show_object(o)
            
        if temporary:
            bm = bmesh.new()
        else:
            bm = BmeshEdit._start(context, o)
            bm.clear()
            
        logger.debug("BmeshEdit.buildmesh() %s :%.2f seconds", o.name, time.time() - tim)
        
        # BmeshEdit.ensure_bmesh(bm)
        _verts = bm.verts
        _faces = bm.faces
        _new_vert = _verts.new
        _new_face = _faces.new
        
        for v in verts:
            _new_vert(v)
        
        _verts.ensure_lookup_table()
        _verts.index_update()
        logger.debug("BmeshEdit.buildmesh() verts :%.2f seconds", time.time() - tim)
        
        for f in faces:
            _new_face([_verts[i] for i in f])
        
        # bm.edges.ensure_lookup_table()
        _faces.ensure_lookup_table()
        _faces.index_update()
        logger.debug("BmeshEdit.buildmesh() faces :%.2f seconds", time.time() - tim)
        
        if matids is not None:
            BmeshEdit._matids(bm, matids)
            logger.debug("BmeshEdit.buildmesh() _matids :%.2f seconds", time.time() - tim)
        
        if uvs is not None:
            BmeshEdit._uvs(bm, uvs)

        if temporary:
            return bm

        if weld:
            bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.001)
        
        BmeshEdit._end(bm, o)
        logger.debug("BmeshEdit.buildmesh() _end :%.2f seconds", time.time() - tim)
        
        if auto_smooth:
            # bpy.ops.mesh.faces_shade_smooth()
            bpy.ops.object.shade_smooth()
            o.data.use_auto_smooth = True
        else:
            bpy.ops.object.shade_flat()
            # bpy.ops.mesh.faces_shade_flat()
            
        if clean:
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_all(action='SELECT')
        
            bpy.ops.mesh.delete_loose()
        
            bpy.ops.object.mode_set(mode='OBJECT')
        
        logger.debug("BmeshEdit.buildmesh() mesh ops :%.2f seconds", time.time() - tim)
        
        if o is not None and not vis_state:
            objman.hide_object(o)
               
    @staticmethod
    def addmesh(context, o, verts, faces, matids=None, uvs=None, weld=False, clean=False, auto_smooth=True):
        bm = BmeshEdit._start(context, o)
        nv = len(bm.verts)
        nf = len(bm.faces)

        for v in verts:
            bm.verts.new(v)

        bm.verts.ensure_lookup_table()

        for f in faces:
            bm.faces.new([bm.verts[nv + i] for i in f])

        bm.faces.ensure_lookup_table()

        if matids is not None:
            for i, matid in enumerate(matids):
                bm.faces[nf + i].material_index = matid

        if uvs is not None:
            layer = bm.loops.layers.uv.verify()
            for i, face in enumerate(bm.faces[nf:]):
                for j, loop in enumerate(face.loops):
                    loop[layer].uv = uvs[i][j]

        if weld:
            bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.001)
        BmeshEdit._end(bm, o)
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        if auto_smooth:
            bpy.ops.mesh.faces_shade_smooth()
            o.data.use_auto_smooth = True
        else:
            bpy.ops.mesh.faces_shade_flat()
        if clean:
            bpy.ops.mesh.delete_loose()
        bpy.ops.object.mode_set(mode='OBJECT')

    @staticmethod
    def bevel(context, o,
            offset,
            offset_type=0,
            segments=1,
            profile=0.5,
            vertex_only=False,
            clamp_overlap=True,
            material=-1,
            use_selection=True):
        """
        /* Bevel offset_type slot values */
        enum {
          BEVEL_AMT_OFFSET,
          BEVEL_AMT_WIDTH,
          BEVEL_AMT_DEPTH,
          BEVEL_AMT_PERCENT
        };
        """
        bm = bmesh.new()
        bm.from_mesh(o.data)
        BmeshEdit.ensure_bmesh(bm)
        if use_selection:
            geom = [v for v in bm.verts if v.select]
            geom.extend([ed for ed in bm.edges if ed.select])
        else:
            geom = bm.verts[:]
            geom.extend(bm.edges[:])

        bmesh.ops.bevel(bm,
            geom=geom,
            offset=offset,
            offset_type=offset_type,
            segments=segments,
            profile=profile,
            vertex_only=vertex_only,
            clamp_overlap=clamp_overlap,
            material=material)

        bm.to_mesh(o.data)
        bm.free()

    @staticmethod
    def bissect(context, o,
            plane_co,
            plane_no,
            dist=0.001,
            use_snap_center=False,
            clear_outer=True,
            clear_inner=False,
            in_place=True
            ):

        bm = bmesh.new()
        bm.from_mesh(o.data)
        BmeshEdit.ensure_bmesh(bm)
        geom = bm.verts[:]
        geom.extend(bm.edges[:])
        geom.extend(bm.faces[:])

        bmesh.ops.bisect_plane(bm,
            geom=geom,
            dist=dist,
            plane_co=plane_co,
            plane_no=plane_no,
            use_snap_center=False,
            clear_outer=clear_outer,
            clear_inner=clear_inner
            )
        if in_place:
            bm.to_mesh(o.data)
            bm.free()
            bm = None
        return bm
        
    @staticmethod
    def solidify(context, o, amt, floor_bottom=False, altitude=0):
        bm = bmesh.new()
        bm.from_mesh(o.data)
        BmeshEdit.ensure_bmesh(bm)
        
        geom = bm.faces[:]
        bmesh.ops.solidify(bm, geom=geom, thickness=amt)
        if floor_bottom:
            for v in bm.verts:
                if not v.select:
                    v.co.z = altitude
        bm.to_mesh(o.data)
        bm.free()

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
