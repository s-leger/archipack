



import bpy
import bmesh

from .bmesh_utils import BmeshEdit as bmed
from .archipack_autoboolean import ArchipackBoolManager

class archipack_envi_zone():
    """
        Envi zone mesh generator
    """
    def __init__(self):
        self.walls = []
        self.floors = []
        self.slabs = []
        self.roofs = []
        self.ceils = []
        self.windows = []
        self.doors = []

    def add_zone(self, wall):
        """
            Zone
            @TODO:
                use a zone line
        """
        self.walls.append(wall)
        z = wall.matrix_world.translation.z
        if wall.parent is not None:
            for child in wall.parent.children:
                if child.data:
                    if "archipack_window" in child.data:
                        self.windows.append(child)
                    elif "archipack_door" in child.data:
                        self.doors.append(child)
                    elif "archipack_floor" in child.data:
                        self.floors.append(child)
                    elif "archipack_roof" in child.data:
                        if abs(child.matrix_world.translation.z - z) < 0.01:
                            self.roofs.append(child)
                    elif "archipack_slab" in child.data:
                        if child.matrix_world.translation.z - z > 0.1:
                            self.ceils.append(child)
                        else:
                            self.slabs.append(child)

    def p3d(self, wall, verts, t, z_min=0):
        x, y = wall.line.lerp(t)
        z = wall.wall_z + wall.get_z(t)
        verts.append((x, y, z_min))
        verts.append((x, y, z))                            
                            
    def make_wall_mesh(self, wall, z_min, i, verts, faces):
        t = wall.t_step[i]
        f = len(verts)
        self.p3d(wall, verts, t, z_min)
        wall.make_faces(i, f, faces)

    def make_wall(self, d, g, z_min, verts, faces):

        offset = -0.5 * (1 - d.x_offset) * d.width

        g.set_offset(offset)

        nb_segs = len(g.segs) - 1
        if d.closed:
            nb_segs += 1

        for i, wall in enumerate(g.segs):
            wall.param_t(d.step_angle)
            if i < nb_segs:
                for j in range(wall.n_step + 1):
                    self.make_wall_mesh(wall, z_min, j, verts, faces)

    def make_surface(self, verts, idmat):
        #
        bm = bmesh.new()
        for v in verts:
            bm.verts.new(v)
        bm.verts.ensure_lookup_table()
        for i in range(1, len(verts)):
            bm.edges.new((bm.verts[i - 1], bm.verts[i]))
        bm.edges.new((bm.verts[-1], bm.verts[0]))
        bm.edges.ensure_lookup_table()
        bmesh.ops.contextual_create(bm, geom=bm.edges)
        bm.faces.ensure_lookup_table()
        for f in bm.faces:
            f.material_index = idmat
        return bm
    
    def make_hole(self, context, o, d):
        bpy.ops.object.select_all(action='DESELECT')
        o.select = True
        context.scene.objects.active = o
        hole = d.robust_hole(context, o.matrix_world)
        hole.select = True
        context.scene.objects.active = hole
        bpy.ops.object.mode_set(mode='EDIT')
        bm = bmesh.from_edit_mesh(hole.data)
        for f in bm.faces:
            f.select = True
        bmesh.update_edit_mesh(hole.data, True)
        bm.free()
        bpy.ops.object.mode_set(mode='OBJECT')
        return hole
        
    def make_geometry(self, context):

        # floor:
        # either we have a floor object
        # or use a slab one
        z_min = -1000000
        if len(self.floors) > 0:
            for floor in self.floors:
                z = floor.matrix_world.translation.z + floor.data.archipack_floor[0].thickness
                if z > z_min:
                    z_min = z
        else:
            for slab in self.slabs:
                z = slab.matrix_world.translation.z
                if z > z_min:
                    z_min = z
        # ceil:
        # either we have a roof
        # or a ceil
        z_max = -1000000
        if len(self.roofs) > 0:
            has_roof = True
            roof = self.roofs[0]
            z_max = roof.matrix_world.translation.z + roof.data.archipack_roof[0].z
        else:
            has_roof = False
            for ceil in self.ceils:
                z = ceil.matrix_world.translation.z - ceil.data.archipack_slab[0].z
                if z > z_max:
                    z_max = z
        
        hole_mat = bpy.data.materials.new("hole_mat")
        wall_mat = bpy.data.materials.new("wall_mat")
        ceil_mat = bpy.data.materials.new("ceil_mat")
        floor_mat = bpy.data.materials.new("floor_mat")
        win_mat = bpy.data.materials.new("win_mat")
        door_mat = bpy.data.materials.new("door_mat")
        
        holes = []
        for door in self.doors:
            hole = self.make_hole(context, door, door.data.archipack_door[0])
            holes.append((door_mat, hole))
        
        for window in self.windows:
            hole = self.make_hole(context, window, window.data.archipack_window[0])
            holes.append((win_mat, hole))

        # now we know z_min and max
        # build inside geometry for walls
        for wall in self.walls:

            m = bpy.data.meshes.new("envi")
            wall_obj = bpy.data.objects.new("envi", m)
            context.scene.objects.link(wall_obj)
            wall_obj.matrix_world = wall.matrix_world.copy()
            wall_obj.data.materials.append(wall_mat)
            wall_obj.data.materials.append(floor_mat)
            wall_obj.data.materials.append(ceil_mat)
            verts = []
            faces = []
            d = wall.data.archipack_wall2[0]
            g = d.get_generator()
            self.make_wall(d, g, z_min, verts, faces)

            # ground
            g_verts = []
            c_verts = []
            for i, v in enumerate(verts):
                if i % 2 == 0:
                    g_verts.append(v)
                else:
                    c_verts.append(v)
                    
            g_bm = self.make_surface(g_verts, 1)
            
            # ceiling
            c_bm = self.make_surface(c_verts, 2)

            
            # build wall surface
            bmed.buildmesh(context, wall_obj, verts, faces, matids=[0 for f in faces])
            bmed.bmesh_join(context, wall_obj, [g_bm, c_bm], normal_update=True)
            bpy.ops.object.mode_set(mode='EDIT')
            bm = bmesh.from_edit_mesh(wall_obj.data)
            for f in bm.faces:
                f.select = False
            bmesh.update_edit_mesh(wall_obj.data, True)
            bm.free()
            
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.select_all(action="DESELECT")
            
            # create windows and doors surfaces
            to_merge = []
            for mat, hole in holes:
                wc = wall_obj.copy()
                wc.data = wall_obj.data.copy()
                context.scene.objects.link(wc)
                m = wc.modifiers.new('AutoBoolean', 'BOOLEAN')
                m.operation = 'INTERSECT'
                m.object = hole
                wc.select = True
                context.scene.objects.active = wc
                
                bpy.ops.object.modifier_apply(apply_as='DATA', modifier="AutoBoolean")
                wc.data.materials.clear()
                wc.data.materials.append(mat)
                bpy.ops.object.mode_set(mode='EDIT')
                bpy.ops.mesh.delete(type='FACE')
                bm = bmesh.from_edit_mesh(wc.data)
                for f in bm.faces:
                    f.material_index = 0
                bmesh.update_edit_mesh(wc.data, True)
                bm.free()
                bpy.ops.object.mode_set(mode='OBJECT')
                to_merge.append(wc)
            
            # cut wall
            bpy.ops.object.select_all(action="DESELECT")
            hc = None
            for mat, hole in holes:
                hole.select = True
                context.scene.objects.active = hole
                hc = hole
                
            bpy.ops.object.join()
            hc.select = True
            
            bpy.ops.object.mode_set(mode='EDIT')
            bm = bmesh.from_edit_mesh(hc.data)
            for f in bm.faces:
                f.select = True
            bmesh.update_edit_mesh(hc.data, True)
            bm.free()
            bpy.ops.object.mode_set(mode='OBJECT')
            
            m = wall_obj.modifiers.new('AutoBoolean', 'BOOLEAN')
            m.operation = 'DIFFERENCE'
            m.object = hc
            
            bpy.ops.object.select_all(action="DESELECT")
            wall_obj.select = True
            context.scene.objects.active = wall_obj
            bpy.ops.object.modifier_apply(apply_as='DATA', modifier="AutoBoolean")
            
            bpy.ops.object.select_all(action="DESELECT")
            hc.select = True
            context.scene.objects.active = hc
            bpy.ops.object.delete(use_global=False)
            
            bpy.ops.object.select_all(action="DESELECT")
            wall_obj.select = True
            context.scene.objects.active = wall_obj
            
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.delete(type='FACE')
            bpy.ops.object.mode_set(mode='OBJECT')
            
            for hole in to_merge:
                hole.select = True
                
            bpy.ops.object.join()
            
            bpy.ops.object.mode_set(mode='EDIT')
            bm = bmesh.from_edit_mesh(wall_obj.data)
            bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.001)
            bmesh.update_edit_mesh(wall_obj.data, True)
            bm.free()
            
            bpy.ops.mesh.select_all(action='SELECT')
            bpy.ops.mesh.delete_loose()
            
            bpy.ops.mesh.select_all(action='SELECT')
            bpy.ops.mesh.fill_holes(sides=0)
            bpy.ops.mesh.select_all(action='SELECT')
            bpy.ops.mesh.dissolve_degenerate()
            bpy.ops.mesh.select_all(action='SELECT')
            bpy.ops.mesh.faces_shade_flat()
            bpy.ops.object.mode_set(mode='OBJECT')
            
"""
wall = C.object
from archipack.archipack_envi import archipack_envi_zone as envi
vi = envi()
vi.add_zone(wall)
vi.make_geometry(C)
"""