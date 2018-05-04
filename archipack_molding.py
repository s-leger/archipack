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
# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty, BoolProperty, CollectionProperty,
    StringProperty, EnumProperty
    )
from .bmesh_utils import BmeshEdit as bmed
from .panel import Panel as Lofter
from mathutils import Vector, Matrix
from math import sin, cos, pi, acos
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_object import ArchipackCreateTool, ArchipackObject
from .archipack_dimension import DimensionProvider
from .archipack_curveman import ArchipackProfile, ArchipackUserDefinedPath
from .archipack_segments import ArchipackSegment, Generator
import logging
logger = logging.getLogger("archipack")


class MoldingGenerator(Generator):

    def __init__(self, d, o=None):
        Generator.__init__(self, d, o)

    def make_profile(self, profile, idmat,
            x_offset, z_offset, extend, closed, verts, faces, matids, uvs):

        if not self.closed:
            self.segs.pop()

        n_moldings = len(self.segs) - 1

        if n_moldings < 0:
            return

        sections = []

        f = self.segs[0]
        f_last = self.segs[-1].line
        closed_path = (f.line.p - (f_last.p + f_last.v)).length < 0.0001

        # add first section
        if closed_path:
            n = f_last.sized_normal(1, 1)
        else:
            n = f.line.sized_normal(0, 1)

        # n.p = f.lerp(x_offset)
        sections.append(n)

        for s, f in enumerate(self.segs):
            if f.line.length == 0:
                continue
            n = f.line.sized_normal(1, 1)
            sections.append(n)

        user_path_verts = len(sections)
        offset = len(verts)
        if user_path_verts > 0:
            user_path_uv_v = []
            n_sections = user_path_verts - 1
            n = sections[0]
            p0 = n.p
            v0 = n.v.normalized()
            for s, n in enumerate(sections):
                p1 = n.p
                if s > 0:
                    user_path_uv_v.append((p1 - p0).length)
                if s < n_sections:
                    v1 = sections[s + 1].v.normalized()
                elif closed_path:
                    break
                dir = (v0 + v1).normalized()
                scale = min(10, 1 / cos(0.5 * acos(min(1, max(-1, v0 * v1)))))
                for p in profile:
                    # x, y = n.p + scale * (x_offset + p.x) * dir
                    x, y = n.p + scale * p.x * dir
                    z = p.y + z_offset
                    verts.append((x, y, z))
                p0 = p1
                v0 = v1
            if closed_path:
                user_path_verts -= 1
            # build faces using Panel
            lofter = Lofter(
                # closed_shape, index, x, y, idmat
                closed,
                [i for i in range(len(profile))],
                [p.x for p in profile],
                [p.y for p in profile],
                [idmat for i in range(len(profile))],
                closed_path=closed_path,
                user_path_uv_v=user_path_uv_v,
                user_path_verts=user_path_verts
                )
            faces += lofter.faces(1, offset=offset, path_type='USER_DEFINED')
            matids += lofter.mat(1, idmat, idmat, path_type='USER_DEFINED')
            v = Vector((0, 0))
            uvs += lofter.uv(1, v, v, v, v, 0, v, 0, 0, path_type='USER_DEFINED')


def update(self, context):
    if self.auto_update:
        self.update(context)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


class archipack_molding_part(ArchipackSegment, PropertyGroup):
    manipulators = CollectionProperty(type=archipack_manipulator)

    def get_datablock(self, o):
        return archipack_molding.datablock(o)


class archipack_molding(
        ArchipackObject,
        ArchipackProfile,
        ArchipackUserDefinedPath,
        Manipulable,
        DimensionProvider,
        PropertyGroup):

    parts = CollectionProperty(type=archipack_molding_part)

    x_offset = FloatProperty(
            name="Offset",
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    profil = EnumProperty(
            name="Profil",
            items=(
                ('SQUARE', 'Square', '', 0),
                ('CIRCLE', 'Circle', '', 1),
                ('USER', 'User defined', '', 2)
                ),
            default='SQUARE',
            update=update
            )
    profil_x = FloatProperty(
            name="Width",
            min=0.001,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    profil_y = FloatProperty(
            name="Height",
            min=0.001,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    profil_radius = FloatProperty(
            name="Radius",
            min=0.001,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    closed = BoolProperty(
            default=False,
            update=update_manipulators
            )
    always_closed = BoolProperty(
            default=False,
            options={'SKIP_SAVE'}
            )

    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )

    def setup_manipulators(self):

        if len(self.manipulators) == 0:
            s = self.manipulators.add()
            s.prop1_name = "width"
            s = self.manipulators.add()
            s.prop1_name = "height"
            s.normal = Vector((0, 1, 0))

        self.setup_parts_manipulators('profil_y')

    def from_spline(self, context, wM, resolution, spline):

        o = self.find_in_selection(context)

        if o is None:
            return

        pts = self.coords_from_spline(
            spline,
            wM,
            resolution
            )

        if len(pts) < 2:
            return

        if self.user_defined_reverse:
            pts = list(reversed(pts))

        o.matrix_world = Matrix.Translation(pts[0].copy())

        auto_update = self.auto_update
        self.auto_update = False

        self.closed = spline.use_cyclic_u
        if self.closed:
            pts.pop()

        self.from_points(pts)
        self.auto_update = auto_update

    def refresh_profile_size(self, context, x, y):
        self.profil_x = x
        self.profil_y = y

    def get_generator(self, o=None):
        g = MoldingGenerator(self, o)
        for part in self.parts:
            g.add_part(part)
        g.set_offset(self.x_offset)
        g.close(self.x_offset)
        return g

    def update(self, context, manipulable_refresh=False):

        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # clean up manipulators before any data model change
        if manipulable_refresh:
            self.manipulable_disable(context)

        self.update_parts()

        verts = []
        faces = []
        matids = []
        uvs = []

        g = self.get_generator()
        g.locate_manipulators()

        # reset user def posts

        if self.profil == 'USER':
            curve = self.update_profile(context)
            if curve and curve.type == 'CURVE':
                sx, sy = 1, 1
                if self.user_profile_dimension.x > 0:
                    sx = self.profil_x / self.user_profile_dimension.x
                if self.user_profile_dimension.y > 0:
                    sy = self.profil_y / self.user_profile_dimension.y

                wM = Matrix([
                    [sx, 0, 0, 0],
                    [0, sy, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]
                    ])

                for spline in curve.data.splines:

                    molding = self.coords_from_spline(spline, wM, 12, ccw=True)
                    closed = molding[0] == molding[-1]
                    if closed:
                        molding.pop()
                    g.make_profile(molding, 0, self.x_offset,
                        0, 0, True, verts, faces, matids, uvs)
            else:
                x = self.profil_x
                y = self.profil_y
                molding = [Vector((0, y)), Vector((0, 0)), Vector((x, 0)), Vector((x, y))]
                g.make_profile(molding, 0, self.x_offset,
                    0, 0, True, verts, faces, matids, uvs)
        else:
            if self.profil == 'SQUARE':
                x = self.profil_x
                y = self.profil_y
                molding = [Vector((0, y)), Vector((0, 0)), Vector((x, 0)), Vector((x, y))]

            elif self.profil == 'CIRCLE':
                x = self.profil_x
                y = self.profil_y
                r = min(self.profil_radius, x, y)
                segs = 6
                da = pi / (2 * segs)
                molding = [Vector((0, y)), Vector((0, 0)), Vector((x, 0))]
                molding.extend([
                    Vector((x + r * (cos(a * da) - 1), y + r * (sin(a * da) - 1)))
                    for a in range(segs + 1)
                    ])

            g.make_profile(molding, 0, self.x_offset,
                0, 0, True, verts, faces, matids, uvs)

        bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=uvs, weld=False, clean=False)

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

        n_parts = self.n_parts
        if self.closed:
            n_parts += 1

        for i, part in enumerate(self.parts):

            if i < n_parts:

                if i > 0:
                    # start angle
                    self.manip_stack.append(part.manipulators[0].setup(context, o, part))

                # length / radius + angle
                self.manip_stack.append(part.manipulators[1].setup(context, o, part))

                # index
                self.manip_stack.append(part.manipulators[3].setup(context, o, self))

            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))

        for m in self.manipulators:
            self.manip_stack.append(m.setup(context, o, self))


class ARCHIPACK_PT_molding(Panel):
    bl_idname = "ARCHIPACK_PT_molding"
    bl_label = "Molding"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_molding.poll(context.active_object)

    def draw(self, context):
        prop = archipack_molding.datablock(context.active_object)
        if prop is None:
            return
        layout = self.layout

        # template with icon right on object
        # layout.template_icon_view(prop, "preset", show_labels=True, scale=10)

        row = layout.row(align=True)
        row.operator('archipack.manipulate', icon='HAND')
        box = layout.box()
        # box.label(text="Styles")
        row = box.row(align=True)
        row.operator("archipack.molding_preset_menu", text=bpy.types.ARCHIPACK_OT_molding_preset_menu.bl_label)
        row.operator("archipack.molding_preset", text="", icon='ZOOMIN')
        row.operator("archipack.molding_preset", text="", icon='ZOOMOUT').remove_active = True
        box = layout.box()
        expand = prop.template_user_path(context, box)
        if expand:
            if prop.user_defined_path is not "":
                box.prop(prop, 'user_defined_reverse')

        box = layout.box()
        box.prop(prop, 'x_offset')

        prop.template_parts(context, layout, draw_type=False)

        box = layout.box()
        row = box.row(align=False)
        icon = "TRIA_RIGHT"
        if prop.profile_expand:
            icon = "TRIA_DOWN"

        row.prop(prop, 'profile_expand', icon=icon, icon_only=True, text="Profil", emboss=True)
        row.prop(prop, 'profil', text="")

        if prop.profile_expand:
            box.prop(prop, 'profil_x')
            box.prop(prop, 'profil_y')
            if prop.profil == 'CIRCLE':
                box.prop(prop, 'profil_radius')
            if prop.profil == 'USER':
                prop.draw_user_profile(context, box)


# ------------------------------------------------------------------
# Define operator class to create object
# TODO: turn into internal operator, allow molding from wall, from curve and regular
# ------------------------------------------------------------------


class ARCHIPACK_OT_molding(ArchipackCreateTool, Operator):
    bl_idname = "archipack.molding"
    bl_label = "Molding"
    bl_description = "Molding"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    def create(self, context):
        m = bpy.data.meshes.new("Molding")
        o = bpy.data.objects.new("Molding", m)
        d = m.archipack_molding.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        self.link_object_to_scene(context, o)
        self.select_object(context, o, True)
        self.load_preset(d)
        self.add_material(o)
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = context.scene.cursor_location
            self.select_object(context, o, True)
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_molding_from_curve(ArchipackCreateTool, Operator):
    bl_idname = "archipack.molding_from_curve"
    bl_label = "Molding curve"
    bl_description = "Create molding from curve"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return context.active_object is not None and context.active_object.type == 'CURVE'

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        o = None
        curve = context.active_object
        sel = []
        for i, spline in enumerate(curve.data.splines):
            bpy.ops.archipack.molding('INVOKE_DEFAULT', auto_manipulate=False)
            o = context.active_object
            d = archipack_molding.datablock(o)
            d.auto_update = False
            d.user_defined_spline = i
            d.user_defined_path = curve.name
            d.auto_update = True
            sel.append(o)

        for obj in sel:
            self.select_object(context, obj)

        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            self.select_object(context, o, True)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_molding_from_wall(ArchipackCreateTool, Operator):
    bl_idname = "archipack.molding_from_wall"
    bl_label = "->Floor"
    bl_description = "Create molding from a wall"
    bl_category = 'Archipack'
    bl_options = {'REGISTER'}

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o is not None and o.data is not None and 'archipack_wall2' in o.data

    def create(self, context):
        tim = time.time()
        m = bpy.data.meshes.new("Molding")
        o = bpy.data.objects.new("Molding", m)
        d = m.archipack_molding.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        d.auto_update = False
        self.link_object_to_scene(context, o)
        logger.debug("create link_object_to_scene() :%.4f seconds", time.time() - tim)
        self.select_object(context, o, True)
        bpy.ops.script.python_file_run(filepath=self.filepath)
        logger.debug("create python_file_run() :%.4f seconds", time.time() - tim)
        return o

    def molding_from_wall(self, context, w, wd):
        """
         Create flooring from surrounding wall
         Use slab cutters, windows and doors, T childs walls
        """
        tim = time.time()
        # wall is either a single or collection of polygons
        io, wall, childs = wd.as_geom(context, w, 'FLOOR_MOLDINGS', [], [], [])
        logger.debug("molding_from_wall wd.as_geom :%.4f seconds", time.time() - tim)

        # find slab holes if any
        o = None
        sel = []
        # MultiLineString
        if wall.type_id == 5:
            lines = wall.geoms
        else:
            lines = [wall]

        for i, line in enumerate(lines):
            boundary = io._to_curve(line, "{}-boundary".format(w.name), '2D')
            boundary.location.z = w.matrix_world.translation.z
            logger.debug("molding_from_wall boundary :%.4f seconds", time.time() - tim)
            o = self.create(context)
            logger.debug("molding_from_wall create :%.4f seconds", time.time() - tim)
            o.matrix_world = w.matrix_world.copy()
            d = archipack_molding.datablock(o)
            d.user_defined_path = boundary.name
            d.user_defined_spline = i
            d.auto_update = True
            logger.debug("molding_from_wall update :%.4f seconds", time.time() - tim)
            self.delete_object(context, boundary)
            d.user_defined_path = ""
            sel.append(o)
            self.unselect_object(o)

        self.select_object(context, w, True)
        for obj in sel:
            self.select_object(context, obj)
        bpy.ops.archipack.add_reference_point()
        self.unselect_object(w)
        logger.debug("molding_from_wall() :%.4f seconds", time.time() - tim)
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            wall = context.active_object
            wd = wall.data.archipack_wall2[0]
            o = self.molding_from_wall(context, wall, wd)
            self.select_object(context, o, True)
            # if self.auto_manipulate:
            #    bpy.ops.archipack.manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to load / save presets
# ------------------------------------------------------------------


class ARCHIPACK_OT_molding_preset_from_wall(PresetMenuOperator, Operator):
    bl_description = "Create molding from wall"
    bl_idname = "archipack.molding_preset_from_wall"
    bl_label = "-> Molding"
    preset_subdir = "archipack_molding"
    preset_operator = StringProperty(
        options={'SKIP_SAVE'},
        default="archipack.molding_from_wall"
    )

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o and o.data and "archipack_wall2" in o.data


class ARCHIPACK_OT_molding_preset_from_curve(PresetMenuOperator, Operator):
    bl_description = "Create molding(s) from a curve"
    bl_idname = "archipack.molding_preset_from_curve"
    bl_label = "-> Molding"
    preset_subdir = "archipack_molding"
    preset_operator = StringProperty(
        options={'SKIP_SAVE'},
        default="archipack.molding_from_curve"
    )

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o and o.type == 'CURVE'


class ARCHIPACK_OT_molding_preset_create(PresetMenuOperator, Operator):
    bl_description = "Show Molding presets and create object at cursor location"
    bl_idname = "archipack.molding_preset_create"
    bl_label = "Molding Styles"
    preset_subdir = "archipack_molding"


class ARCHIPACK_OT_molding_preset_menu(PresetMenuOperator, Operator):
    bl_description = "Display Molding presets"
    bl_idname = "archipack.molding_preset_menu"
    bl_label = "Molding Styles"
    preset_subdir = "archipack_molding"


class ARCHIPACK_OT_molding_preset(ArchipackPreset, Operator):
    """Add a Molding Preset"""
    bl_idname = "archipack.molding_preset"
    bl_label = "Add Molding Style"
    preset_menu = "ARCHIPACK_OT_molding_preset_menu"

    @property
    def blacklist(self):
        return ['manipulators', 'n_parts', 'parts', 'user_defined_path', 'user_defined_spline']


def register():
    bpy.utils.register_class(archipack_molding_part)
    bpy.utils.register_class(archipack_molding)
    Mesh.archipack_molding = CollectionProperty(type=archipack_molding)
    bpy.utils.register_class(ARCHIPACK_OT_molding_preset_menu)
    bpy.utils.register_class(ARCHIPACK_OT_molding_preset_create)
    bpy.utils.register_class(ARCHIPACK_OT_molding_preset_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_molding_preset_from_wall)
    bpy.utils.register_class(ARCHIPACK_PT_molding)
    bpy.utils.register_class(ARCHIPACK_OT_molding)
    bpy.utils.register_class(ARCHIPACK_OT_molding_preset)
    bpy.utils.register_class(ARCHIPACK_OT_molding_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_molding_from_wall)


def unregister():
    bpy.utils.unregister_class(archipack_molding_part)
    bpy.utils.unregister_class(archipack_molding)
    del Mesh.archipack_molding
    bpy.utils.unregister_class(ARCHIPACK_OT_molding_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_OT_molding_preset_create)
    bpy.utils.unregister_class(ARCHIPACK_OT_molding_preset_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_molding_preset_from_wall)
    bpy.utils.unregister_class(ARCHIPACK_PT_molding)
    bpy.utils.unregister_class(ARCHIPACK_OT_molding)
    bpy.utils.unregister_class(ARCHIPACK_OT_molding_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_molding_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_molding_from_wall)
