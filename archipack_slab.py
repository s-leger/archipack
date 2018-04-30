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
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty, BoolProperty, IntProperty,
    StringProperty, EnumProperty,
    CollectionProperty
    )
import bmesh
from mathutils import Vector, Matrix
from math import pi
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_object import (
    ArchipackCreateTool,
    ArchipackObject,
    ArchipackObjectsManager
    )
from .archipack_cutter import (
    CutAblePolygon, CutAbleGenerator,
    ArchipackCutter,
    ArchipackCutterPart,
    update_operation
    )
from .archipack_dimension import DimensionProvider
from .archipack_curveman import ArchipackUserDefinedPath
from .archipack_segments import ArchipackSegment


class SlabGenerator(CutAblePolygon, CutAbleGenerator):

    def __init__(self, d, o=None):
        CutAbleGenerator.__init__(self, d, o)

    def get_verts(self, verts):
        verts.extend([s.p0.to_3d() for s in self.segs])

    def cut(self, context, o):
        """
            either external or holes cuts
        """
        # use offset segs as base
        self.as_lines(step_angle=0.0502)
        self.limits()

        for b in o.children:
            d = archipack_slab_cutter.datablock(b)
            if d is not None:
                tM = o.matrix_world.inverted() * b.matrix_world
                g = d.ensure_direction(tM)
                # g.change_coordsys(b.matrix_world, o.matrix_world)
                self.slice(g)

    def slab(self, context, o, d):

        verts = []
        self.get_verts(verts)
        if len(verts) > 2:

            # ensure verts are CCW
            if d.is_cw(verts):
                verts = list(reversed(verts))

            bm = bmesh.new()

            for v in verts:
                bm.verts.new(v)
            bm.verts.ensure_lookup_table()
            for i in range(1, len(verts)):
                bm.edges.new((bm.verts[i - 1], bm.verts[i]))
            bm.edges.new((bm.verts[-1], bm.verts[0]))
            bm.edges.ensure_lookup_table()
            bmesh.ops.contextual_create(bm, geom=bm.edges)

            self.cut_holes(bm, self)

            bmesh.ops.dissolve_limit(bm,
                        angle_limit=0.01,
                        use_dissolve_boundaries=False,
                        verts=bm.verts,
                        edges=bm.edges,
                        delimit=1)

            bm.to_mesh(o.data)
            bm.free()
            # geom = bm.faces[:]
            # verts = bm.verts[:]
            # bmesh.ops.solidify(bm, geom=geom, thickness=d.z)

            # merge with object
            # bmed.bmesh_join(context, o, [bm], normal_update=True)

            bpy.ops.object.mode_set(mode='OBJECT')


def update(self, context):
    self.update(context)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


def update_path(self, context):
    self.update_path(context)


materials_enum = (
            ('0', 'Ceiling', '', 0),
            ('1', 'White', '', 1),
            ('2', 'Concrete', '', 2),
            ('3', 'Wood', '', 3),
            ('4', 'Metal', '', 4),
            ('5', 'Glass', '', 5)
            )


class archipack_slab_material(PropertyGroup):
    index = EnumProperty(
        items=materials_enum,
        default='4',
        update=update
        )

    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            d = archipack_slab.datablock(o)
            if d:
                for part in d.rail_mat:
                    if part == self:
                        return d
        return None

    def update(self, context):
        d = self.find_in_selection(context)
        if d is not None:
            d.update(context)


class archipack_slab_relocate_child(PropertyGroup):
    child_name = StringProperty()
    child_idx = IntProperty(
        default=-1,
        description="Index of linked child part (startpoint)"
        )
    wall_id0 = IntProperty(
        description="Index of wall part"
        )
    wall_id1 = IntProperty(
        description="Index of wall part"
        )
    d0 = FloatProperty()
    d1 = FloatProperty()

    def get_child(self, context):
        d = None
        o = context.scene.objects.get(self.child_name)
        if o is not None and o.data is not None:
            cd = o.data
            if 'archipack_fence' in cd:
                d = cd.archipack_fence[0]
        return o, d


class archipack_slab_part(ArchipackSegment, PropertyGroup):
    manipulators = CollectionProperty(type=archipack_manipulator)

    linked_idx = IntProperty(default=-1)

    def get_datablock(self, o):
        return archipack_slab.datablock(o)

    def draw_insert(self, context, layout, index):
        row = layout.row(align=True)
        row.operator("archipack.segment_insert", text="Split").index = index
        row.operator("archipack.slab_balcony", text="Balcony").index = index
        if index > 0:
            row.operator("archipack.segment_remove", text="Remove").index = index


class archipack_slab(ArchipackObject, ArchipackUserDefinedPath, Manipulable, DimensionProvider, PropertyGroup):
    # boundary
    parts = CollectionProperty(type=archipack_slab_part)
    closed = BoolProperty(
            default=True,
            name="Close",
            options={'SKIP_SAVE'}
            )
    always_closed = BoolProperty(
            default=True,
            name="Always closed geometry",
            description="Flag indicate whenever geometry path is always closed",
            options={'SKIP_SAVE'}
            )

    # UI layout related
    x_offset = FloatProperty(
            name="Offset",
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    z = FloatProperty(
            name="Thickness",
            default=0.3, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    # store childs with parts points locations to relocate
    reloc_childs = CollectionProperty(type=archipack_slab_relocate_child)

    # Flag to prevent mesh update while making bulk changes over variables
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )

    def get_generator(self, o=None):
        g = SlabGenerator(self, o)
        for part in self.parts:
            # type, radius, da, length
            g.add_part(part)
        g.set_offset(self.x_offset)
        g.close(self.x_offset)
        g.locate_manipulators()
        # here segs are without offset
        # seg.line contains offset
        return g

    def new_part(self, where, type, length, radius, offset, a0, da):
        idx = len(self.parts)
        p = self.parts.add()
        p.type = type
        p.length = length
        p.offset = offset
        p.da = da
        p.a0 = a0
        self.n_parts += 1
        self.parts.move(idx, where + 1)
        return p

    def after_insert_part(self, context, o, where, distance):
        for c in self.reloc_childs:
            if c.wall_id0 > where:
                c.wall_id0 += 1
            if c.wall_id1 > where:
                c.wall_id1 += 1

    def insert_balcony(self, context, where):

        self.manipulable_disable(context)
        self.auto_update = False

        # the part we do split
        part_0 = self.parts[where]
        part_0.length /= 3
        part_0.da /= 3

        # store here to keep a valid ref
        type = part_0.type
        length = part_0.length
        radius = part_0.radius
        offset = part_0.offset
        da = part_0.da

        # 1st part 90deg
        self.new_part(where, "S_SEG", 1.5, 0, 0, -pi / 2, da)

        # 2nd part -90deg
        self.new_part(where + 1, type, length, radius + 1.5, 0, pi / 2, da)

        # 3th part -90deg
        self.new_part(where + 2, "S_SEG", 1.5, 0, 0, pi / 2, da)

        # 4th part -90deg
        self.new_part(where + 3, type, length, radius, offset, -pi / 2, da)

        self.setup_manipulators()

        for c in self.reloc_childs:
            if c.wall_id0 > where:
                c.wall_id0 += 4
            if c.wall_id1 > where:
                c.wall_id1 += 4

        self.auto_update = True

        g = self.get_generator()
        o = context.active_object
        bpy.ops.archipack.fence(auto_manipulate=False)
        c = context.active_object
        c.location.z = o.matrix_world.translation.z

        self.select_object(context, c)

        d = c.data.archipack_fence[0]
        d.n_parts = 3
        self.unselect_object(c)
        # link to o
        c.location = Vector((0, 0, 0))
        c.parent = o
        c.location = g.segs[where + 1].p0.to_3d()

        self.add_relocate_child(c.name, 0, g.segs[where + 1].p0, g.segs)
        self.add_relocate_child(c.name, 1, g.segs[where + 2].p0, g.segs)
        self.add_relocate_child(c.name, 2, g.segs[where + 3].p0, g.segs)
        self.add_relocate_child(c.name, 3, g.segs[where + 4].p0, g.segs)
        # c.matrix_world.translation = g.segs[where + 1].p1.to_3d()
        self.select_object(context, o, True)
        self.relocate_childs(context, o, g)

    def add_relocate_child(self, name, c_idx, pt, segs):
        """
         Find 2 closest segments
         Store distance to segments, segment indexes
        """
        closest = []
        for idx, seg in enumerate(segs):
            res, dist, t = seg.point_sur_segment(pt)
            abs_d = abs(dist)
            closest.append((abs_d, -dist, idx))

        if len(closest) > 1:
            closest.sort(key=lambda s: s[0])
            # filter out colinear segment delta angle < 0.057 deg
            s0 = segs[closest[0][2]]
            i = 1
            i2 = closest[1][2]
            while i < len(closest) and abs(segs[i2].delta_angle(s0)) < 0.001:
                i += 1
                i2 = closest[i][2]
            c = self.reloc_childs.add()
            c.child_name = name
            c.child_idx = c_idx
            c.d0, c.wall_id0 = closest[0][1:3]
            c.d1, c.wall_id1 = closest[i][1:3]
            print("{} c:{} w0:{} w1:{} d0:{} d1:{}".format(name, c_idx, c.wall_id0, c.wall_id1, c.d0, c.d1))
        else:
            print("skip", name, len(closest))

    def setup_childs(self, context, o):
        """
            Store childs
        """
        self.reloc_childs.clear()
        g = self.get_generator(o)

        for c in o.children:
            cd = c.data
            if cd and "archipack_fence" in cd:
                d = cd.archipack_fence[0]
                if d is not None:
                    cg = d.get_generator(c)

                    for c_idx, seg in enumerate(cg.segs):
                        p = seg.p0
                        self.add_relocate_child(c.name, c_idx, p, g.segs)

                    if not d.closed:
                        p = cg.segs[-1].p1
                        self.add_relocate_child(c.name, len(cg.segs), p, g.segs)

    def relocate_childs(self, context, o, g):
        """
            Move and resize childs after edition
        """
        tM = o.matrix_world
        # build a dict for each child
        soft = {}
        for c in self.reloc_childs:
            if c.child_name not in soft:
                soft[c.child_name] = []
            soft[c.child_name].append(c)

        for name, child in soft.items():
            c, d = child[0].get_child(context)
            if c is None:
                print("skip", name)
                continue

            itM = tM.inverted() * c.matrix_world
            # apply changes to generator
            cg = d.get_generator(itM)
            n_segs = len(cg.segs)
            for cd in child:
                print("relocate {} c:{} w0:{} w1:{} d0:{} d1:{}".format(
                    name,
                    cd.child_idx,
                    cd.wall_id0,
                    cd.wall_id1,
                    cd.d0,
                    cd.d1))
                s0 = g.segs[cd.wall_id0].offset(cd.d0)
                s1 = g.segs[cd.wall_id1].offset(cd.d1)
                
                # p is in o coordsys
                if "C_" in self.parts[cd.wall_id0].type:
                    res, p, u, v = s0.intersect_arc(s1)
                elif "C_" in self.parts[cd.wall_id1].type:
                    res, p, v, u = s1.intersect_arc(s0)
                else:
                    res, p, u, v = s0.intersect_ext(s1)
                    
                if p != 0:
                    if cd.child_idx < n_segs:
                        cg.segs[cd.child_idx].p0 = p
                        if cd.child_idx > 0:
                            cg.segs[cd.child_idx - 1].p1 = p
                        else:
                            d.move_object(c, tM * p.to_3d())
                    else:
                        cg.segs[cd.child_idx - 1].p1 = p

            # update data from generator
            last = None
            for i, seg in enumerate(cg.segs):
                p = d.parts[i]
                if last is None:
                    # first a0 is absolute
                    p.a0 = cg.a0
                else:
                    # a0 is delta between segs
                    p.a0 = seg.delta_angle(last)
                if "C_" in p.type:
                    p.da = seg.da
                    p.radius = seg.r
                else:
                    p.length = seg.length
                last = seg

            self.select_object(context, c, True)
            d.update(context)
            self.unselect_object(c)

        self.select_object(context, o, True)

    def after_remove_part(self, context, o, where, distance):
        """
         Rebuild childs index and location
         When removing a part
        """
        for c in self.reloc_childs:
            if c.wall_id0 >= where:
                c.wall_id0 -= 1
            if c.wall_id1 >= where:
                c.wall_id1 -= 1

    def setup_manipulators(self):

        if len(self.manipulators) < 1:
            s = self.manipulators.add()
            s.type_key = "SIZE"
            s.prop1_name = "z"
            s.normal = Vector((0, 1, 0))

        self.setup_parts_manipulators('z')

    def from_spline(self, context, wM, resolution, spline):
        o = self.find_in_selection(context)

        if o is None:
            return

        pts = self.coords_from_spline(
            spline,
            wM,
            resolution,
            ccw=True)

        if len(pts) < 3:
            return

        o.matrix_world = Matrix.Translation(pts[0].copy())
        auto_update = self.auto_update
        self.auto_update = False
        self.from_points(pts)
        self.auto_update = auto_update

    def make_surface(self, o, verts):
        bm = bmesh.new()
        for v in verts:
            bm.verts.new(v)
        bm.verts.ensure_lookup_table()
        for i in range(1, len(verts)):
            bm.edges.new((bm.verts[i - 1], bm.verts[i]))
        bm.edges.new((bm.verts[-1], bm.verts[0]))
        bm.edges.ensure_lookup_table()
        bmesh.ops.contextual_create(bm, geom=bm.edges)
        bm.to_mesh(o.data)
        bm.free()

    def unwrap_uv(self, o):
        bm = bmesh.new()
        bm.from_mesh(o.data)
        for face in bm.faces:
            face.select = face.material_index > 0
        bm.to_mesh(o.data)
        bpy.ops.uv.cube_project(scale_to_bounds=False, correct_aspect=True)

        for face in bm.faces:
            face.select = face.material_index < 1
        bm.to_mesh(o.data)
        bpy.ops.uv.smart_project(use_aspect=True, stretch_to_bounds=False)
        bm.free()

    def update(self, context, manipulable_refresh=False, update_childs=False):

        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # clean up manipulators before any data model change
        if manipulable_refresh:
            self.manipulable_disable(context)

        changed = self.update_parts()
        g = self.get_generator()

        if changed or update_childs:
            self.setup_childs(context, o)

        # relocate before cutting segs
        self.relocate_childs(context, o, g)

        self.select_object(context, o, True)

        # cut transfer .line with offset into g.segs
        g.cut(context, o)
        g.slab(context, o, self)

        modif = o.modifiers.get('Slab')
        if modif is None:
            modif = o.modifiers.new('Slab', 'SOLIDIFY')
            modif.use_quality_normals = True
            modif.use_even_offset = True
            modif.material_offset_rim = 2
            modif.material_offset = 1

        modif.thickness = self.z
        modif.offset = 1.0
        o.data.use_auto_smooth = True
        bpy.ops.object.shade_smooth()

        # Height
        self.manipulators[0].set_pts([
            self.origin,
            self.origin + Vector((0, 0, -self.z)),
            (-1, 0, 0)
            ], normal=g.segs[0].straight(-1, 0).v.to_3d())

        self.update_dimensions(context, o)

        # enable manipulators rebuild
        if manipulable_refresh:
            self.manipulable_refresh = True

        # restore context
        self.restore_context(context)

    def manipulable_setup(self, context):
        """
            NOTE:
            this one assume context.active_object is the instance this
            data belongs to, failing to do so will result in wrong
            manipulators set on active object
        """
        self.manipulable_disable(context)

        o = context.active_object

        self.setup_manipulators()

        for i, part in enumerate(self.parts):
            if i >= self.n_parts:
                break

            if i > 0:
                # start angle
                self.manip_stack.append(part.manipulators[0].setup(context, o, part))

            # length / radius + angle
            self.manip_stack.append(part.manipulators[1].setup(context, o, part))

            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))
            # index
            self.manip_stack.append(part.manipulators[3].setup(context, o, self))

        for m in self.manipulators:
            self.manip_stack.append(m.setup(context, o, self))

    def manipulable_invoke(self, context):
        """
            call this in operator invoke()
        """
        # print("manipulable_invoke")
        if self.manipulate_mode:
            self.manipulable_disable(context)
            return False

        o = context.active_object
        #
        # setup childs manipulators
        self.setup_childs(context, o)
        # generator does update manipulators location
        self.get_generator()
        self.manipulable_setup(context)
        self.manipulate_mode = True

        self._manipulable_invoke(context)

        return True


def update_hole(self, context):
    # update parent only when manipulated
    self.update(context, update_parent=True)


class archipack_slab_cutter_segment(ArchipackCutterPart, PropertyGroup):
    manipulators = CollectionProperty(type=archipack_manipulator)

    def get_datablock(self, o):
        return archipack_slab_cutter.datablock(o)


class archipack_slab_cutter(ArchipackCutter, ArchipackObject, Manipulable, DimensionProvider, PropertyGroup):
    parts = CollectionProperty(type=archipack_slab_cutter_segment)

    def update_points(self, context, o, pts, update_parent=False):
        self.auto_update = False
        self.from_points(pts)
        self.auto_update = True
        if update_parent:
            self.update_parent(context, o)

    def update_parent(self, context, o):
        if o is not None:
            d = archipack_slab.datablock(o.parent)
            if d is not None:
                self.select_object(context, o.parent, True)
                d.update(context)
                self.unselect_object(o.parent)
            self.select_object(context, o, True)


class ARCHIPACK_PT_slab(Panel):
    """Archipack Slab"""
    bl_idname = "ARCHIPACK_PT_slab"
    bl_label = "Slab"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    # bl_context = 'object'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_slab.poll(context.active_object)

    def draw(self, context):
        o = context.active_object
        d = archipack_slab.datablock(o)
        if d is None:
            return
        layout = self.layout
        row = layout.row(align=True)
        row.operator('archipack.manipulate', icon='HAND')
        row.operator("archipack.slab", text="Delete", icon='ERROR').mode = 'DELETE'
        box = layout.box()
        box.operator('archipack.slab_cutter').parent = o.name
        box = layout.box()
        d.template_user_path(context, box)
        box = layout.box()
        box.prop(d, 'z')
        box.prop(d, 'x_offset')
        # box = layout.box()
        # box.prop(d, 'auto_synch')
        d.template_parts(context, layout, draw_type=True)


class ARCHIPACK_PT_slab_cutter(Panel):
    bl_idname = "ARCHIPACK_PT_slab_cutter"
    bl_label = "Slab Cutter"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_slab_cutter.poll(context.active_object)

    def draw(self, context):
        d = archipack_slab_cutter.datablock(context.active_object)
        if d is None:
            return
        layout = self.layout
        d.draw(context, layout)


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_slab_balcony(Operator):
    bl_idname = "archipack.slab_balcony"
    bl_label = "Insert"
    bl_description = "Insert part"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    index = IntProperty(default=0)

    def execute(self, context):
        if context.mode == "OBJECT":
            d = archipack_slab.datablock(context.active_object)
            if d is None:
                return {'CANCELLED'}
            d.insert_balcony(context, self.index)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_slab(ArchipackCreateTool, Operator):
    bl_idname = "archipack.slab"
    bl_label = "Slab"
    bl_description = "Slab"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    mode = EnumProperty(
            items=(
            ('CREATE', 'Create', '', 0),
            ('DELETE', 'Delete', '', 1)
            ),
            default='CREATE'
            )

    def create(self, context):
        m = bpy.data.meshes.new("Slab")
        o = bpy.data.objects.new("Slab", m)
        d = m.archipack_slab.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        # Link object into scene
        self.link_object_to_scene(context, o)

        # select and make active
        self.select_object(context, o, True)

        self.load_preset(d)
        self.add_material(o)
        return o

    def delete(self, context):
        o = context.active_object
        if archipack_slab.filter(o):
            bpy.ops.archipack.disable_manipulate()
            self.delete_object(context, o)

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            if self.mode == 'CREATE':
                bpy.ops.object.select_all(action="DESELECT")
                o = self.create(context)
                o.location = bpy.context.scene.cursor_location
                # select and make active
                self.select_object(context, o, True)
                self.manipulate()
            else:
                self.delete(context)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_slab_from_curve(Operator):
    bl_idname = "archipack.slab_from_curve"
    bl_label = "Slab curve"
    bl_description = "Create a slab from a curve"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    auto_manipulate = BoolProperty(default=True)

    @classmethod
    def poll(self, context):
        return context.active_object is not None and context.active_object.type == 'CURVE'
    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        curve = context.active_object
        bpy.ops.archipack.slab(auto_manipulate=self.auto_manipulate)
        o = context.active_object
        d = archipack_slab.datablock(o)
        spline = curve.data.splines[0]
        d.from_spline(context, curve.matrix_world, 12, spline)

        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            self.create(context)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_slab_from_wall(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.slab_from_wall"
    bl_label = "->Slab"
    bl_description = "Create a slab from a wall"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    auto_manipulate = BoolProperty(default=True)
    ceiling = BoolProperty(default=False)

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o is not None and o.data is not None and 'archipack_wall2' in o.data

    def create(self, context):
        wall = context.active_object
        wd = wall.data.archipack_wall2[0]
        bpy.ops.archipack.slab(auto_manipulate=False)
        o = context.active_object
        d = archipack_slab.datablock(o)
        d.auto_update = False
        # use wall generator
        # so slab match wall's outside
        g = wd.get_generator()
        pts = [s.p0 for s in g.segs]
        if wd.is_cw(pts):
            offset = wd.left_offset
        else:
            offset = wd.left_offset + wd.width
        g.set_offset(offset)
        g.close(offset)
        d.parts.clear()
        d.n_parts = wd.n_parts + 1
        a0 = g.a0
        last = None
        for i, part in enumerate(wd.parts):
            seg = g.segs[i].line
            p = d.parts.add()
            p.type = part.type
            if last is not None:
                a0 = seg.delta_angle(last)
            last = seg
            p.a0 = a0
            if "S_" in part.type:
                p.length = seg.length
            else:
                p.radius = seg.r
                p.da = seg.da

        d.auto_update = True
        x, y = g.segs[0].line.p0
        # pretranslate
        if self.ceiling:
            o.matrix_world = Matrix.Translation(
                Vector((x, y, wd.z - wd.z_offset + d.z))
                ) * wall.matrix_world
        else:
            o.matrix_world = Matrix.Translation(
                Vector((x, y, -wd.z_offset))
                ) * wall.matrix_world
            bpy.ops.object.select_all(action='DESELECT')
            # parenting childs to wall reference point
            self.select_object(context, o)
            self.select_object(context, wall, True)
            bpy.ops.archipack.add_reference_point()
            self.unselect_object(wall)
            self.unselect_object(wall.parent)
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            self.select_object(context, o, True)
            if self.auto_manipulate:
                bpy.ops.archipack.manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_slab_cutter(ArchipackCreateTool, Operator):
    bl_idname = "archipack.slab_cutter"
    bl_label = "Slab Cutter"
    bl_description = "Slab Cutter"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    parent = StringProperty("")
    curve = StringProperty("")

    def create(self, context):
        m = bpy.data.meshes.new("Slab Cutter")
        o = bpy.data.objects.new("Slab Cutter", m)
        d = m.archipack_slab_cutter.add()
        parent = context.scene.objects.get(self.parent)
        curve = context.scene.objects.get(self.curve)

        if parent is not None:
            o.parent = parent
            bbox = parent.bound_box
            angle_90 = pi / 2
            x0, y0, z = bbox[0]
            x1, y1, z = bbox[6]
            x = 0.2 * (x1 - x0)
            y = 0.2 * (y1 - y0)
            o.matrix_world = parent.matrix_world * Matrix.Translation(Vector((-3 * x, 0, 0)))
            # d.auto_update = False
            p = d.parts.add()
            p.a0 = - angle_90
            p.length = y
            p = d.parts.add()
            p.a0 = angle_90
            p.length = x
            p = d.parts.add()
            p.a0 = angle_90
            p.length = y
            p = d.parts.add()
            p.a0 = angle_90
            p.length = x
            d.n_parts = 4
            # d.auto_update = True

        else:
            o.location = context.scene.cursor_location

        # make manipulators selectable
        d.manipulable_selectable = True
        self.link_object_to_scene(context, o)
        self.select_object(context, o, True)
        self.load_preset(d)
        update_operation(d, context)

        if curve is not None:
            d.user_defined_path = curve.name

        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            self.select_object(context, o, True)
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


def register():
    bpy.utils.register_class(archipack_slab_cutter_segment)
    bpy.utils.register_class(archipack_slab_cutter)
    Mesh.archipack_slab_cutter = CollectionProperty(type=archipack_slab_cutter)
    bpy.utils.register_class(ARCHIPACK_OT_slab_cutter)
    bpy.utils.register_class(ARCHIPACK_PT_slab_cutter)

    bpy.utils.register_class(archipack_slab_material)
    bpy.utils.register_class(archipack_slab_relocate_child)
    bpy.utils.register_class(archipack_slab_part)
    bpy.utils.register_class(archipack_slab)
    Mesh.archipack_slab = CollectionProperty(type=archipack_slab)
    bpy.utils.register_class(ARCHIPACK_PT_slab)
    bpy.utils.register_class(ARCHIPACK_OT_slab)
    bpy.utils.register_class(ARCHIPACK_OT_slab_balcony)
    bpy.utils.register_class(ARCHIPACK_OT_slab_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_slab_from_wall)


def unregister():
    bpy.utils.unregister_class(archipack_slab_material)
    bpy.utils.unregister_class(archipack_slab_relocate_child)
    bpy.utils.unregister_class(archipack_slab_part)
    bpy.utils.unregister_class(archipack_slab)
    del Mesh.archipack_slab
    bpy.utils.unregister_class(ARCHIPACK_PT_slab)
    bpy.utils.unregister_class(ARCHIPACK_OT_slab)
    bpy.utils.unregister_class(ARCHIPACK_OT_slab_balcony)
    bpy.utils.unregister_class(ARCHIPACK_OT_slab_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_slab_from_wall)
    del Mesh.archipack_slab_cutter
    bpy.utils.unregister_class(archipack_slab_cutter_segment)
    bpy.utils.unregister_class(archipack_slab_cutter)
    bpy.utils.unregister_class(ARCHIPACK_OT_slab_cutter)
    bpy.utils.unregister_class(ARCHIPACK_PT_slab_cutter)
