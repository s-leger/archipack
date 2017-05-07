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

# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import Operator, PropertyGroup, Mesh, Panel, Menu
from bpy.props import (
    FloatProperty, IntProperty, CollectionProperty, 
    EnumProperty, BoolProperty, StringProperty
    )
from mathutils import Vector

# door component objects (panels, handles ..)
from .bmesh_utils import BmeshEdit as bmed
from .panel import Panel as DoorPanel
from .materialutils import MaterialUtils
from .archipack_handle import create_handle
from .archipack_door_panel import ARCHIPACK_PT_door_panel, ARCHIPACK_OT_door_panel
from .archipack_manipulator import Manipulable
from .archipack_preset import ArchipackPreset

SPACING = 0.005
BATTUE = 0.01
BOTTOM_HOLE_MARGIN = 0.001
FRONT_HOLE_MARGIN = 0.1


def update(self, context):
    self.update(context)


def update_childs(self, context):
    self.update(context, childs_only=True)


class archipack_door(Manipulable, PropertyGroup):
    """
        The frame is the door main object
        parent parametric object
        create/remove/update her own childs
    """
    x = FloatProperty(
            name='width',
            min=0.25, max=10000,
            default=100.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Width', update=update,
            )
    y = FloatProperty(
            name='depth',
            min=0.1, max=10000,
            default=0.20, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='Depth', update=update,
            )
    z = FloatProperty(
            name='height',
            min=0.1, max=10000,
            default=2.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='height', update=update,
            )
    frame_x = FloatProperty(
            name='Width',
            min=0, max=100,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='frame width', update=update,
            )
    frame_y = FloatProperty(
            name='Depth',
            min=-100, max=100,
            default=0.03, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='frame depth', update=update,
            )
    direction = IntProperty(
            name="Direction",
            min=0,
            max=1,
            description="open direction", update=update,
            )
    door_y = FloatProperty(
            name='Depth',
            min=0.001, max=10000,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='depth', update=update,
            )
    door_offset = FloatProperty(
            name='Offset',
            min=0, max=100,
            default=0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='offset', update=update,
            )
    model = IntProperty(
            name="Model",
            min=0,
            max=3,
            default=0,
            description="Model", update=update,
            )
    n_panels = IntProperty(
            name="Panels",
            min=1,
            max=2,
            default=1,
            description="number of panels", update=update
            )
    chanfer = FloatProperty(
            name='chanfer',
            min=0.001, max=10000,
            default=0.005, precision=3, step=0.01,
            unit='LENGTH', subtype='DISTANCE',
            description='chanfer', update=update_childs,
            )
    panel_spacing = FloatProperty(
            name='spacing',
            min=0.001, max=10000,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='distance between panels', update=update_childs,
            )
    panel_bottom = FloatProperty(
            name='bottom',
            min=0.0, max=10000,
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='distance from bottom', update=update_childs,
            )
    panel_border = FloatProperty(
            name='border',
            min=0.001, max=10000,
            default=0.2, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='distance from border', update=update_childs,
            )
    panels_x = IntProperty(
            name="panels h",
            min=1,
            max=50,
            default=1,
            description="panels h", update=update_childs,
            )
    panels_y = IntProperty(
            name="panels v",
            min=1,
            max=50,
            default=1,
            description="panels v", update=update_childs,
            )
    panels_distrib = EnumProperty(
            name='distribution',
            items=(
                ('REGULAR', 'Regular', '', 0),
                ('ONE_THIRD', '1/3 2/3', '', 1)
                ),
            default='REGULAR', update=update_childs,
            )
    handle = EnumProperty(
            name='Handle',
            items=(
                ('NONE', 'No handle', '', 0),
                ('BOTH', 'Inside and outside', '', 1)
                ),
            default='BOTH', update=update_childs,
            )
    hole_margin = FloatProperty(
            name='hole margin',
            min=0.0, max=10000,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            description='how much hole surround wall'
            )
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update
            )
            
    @property
    def frame(self):

        #
        #    _____        y0
        #   |     |___    y1
        #   x         |   y3
        #   |         |
        #   |_________|   y2
        #
        #   x2    x1  x0
        x0 = 0
        x1 = -BATTUE
        x2 = -self.frame_x
        y0 = max(0.25 * self.door_y + 0.0005, self.y / 2 + self.frame_y)
        y1 = max(y0 - 0.5 * self.door_y - self.door_offset, -y0 + 0.001)
        y2 = -y0
        y3 = 0
        return DoorPanel(
            True,           # closed
            [0, 0, 0, 1, 1, 2, 2],  # x index
            [x2, x1, x0],
            [y2, y3, y0, y0, y1, y1, y2],
            [0, 1, 1, 1, 1, 0, 0],  # material index
            closed_path=False
            )

    @property
    def hole(self):
        #
        #    _____   y0
        #   |
        #   x        y2
        #   |
        #   |_____   y1
        #
        #   x0
        x0 = 0
        y0 = self.y / 2 + self.hole_margin
        y1 = -y0
        y2 = 0
        return DoorPanel(
            False,       # closed
            [0, 0, 0],  # x index
            [x0],
            [y1, y2, y0],
            [0, 1, 1],  # material index
            closed_path=True,
            side_cap_front=2,
            side_cap_back=0     # cap index
            )

    @property
    def verts(self):
        # door inner space
        v = Vector((0, 0, 0))
        size = Vector((self.x, self.z, self.y))
        return self.frame.vertices(16, v, v, v, size, v, 0, 0, shape_z=None, path_type='RECTANGLE')

    @property
    def faces(self):
        return self.frame.faces(16, path_type='RECTANGLE')

    @property
    def matids(self):
        return self.frame.mat(16, 0, 0, path_type='RECTANGLE')

    @property
    def uvs(self):
        v = Vector((0, 0, 0))
        size = Vector((self.x, self.z, self.y))
        return self.frame.uv(16, v, v, size, v, 0, 0, 0, 0, path_type='RECTANGLE')

    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        active = context.active_object
        selected = [o for o in context.selected_objects]
        for o in selected:
            if ARCHIPACK_PT_door.params(o) == self:
                return active, selected, o
        return active, selected, None

    def remove_childs(self, context, o, to_remove):
        for child in o.children:
            if to_remove < 1:
                return
            if ARCHIPACK_PT_door_panel.filter(child):
                self.remove_handle(context, child)
                to_remove -= 1
                context.scene.objects.unlink(child)
                bpy.data.objects.remove(child, do_unlink=True)

    def remove_handle(self, context, o):
        handle = self.find_handle(o)
        if handle is not None:
            context.scene.objects.unlink(handle)
            bpy.data.objects.remove(handle, do_unlink=True)

    def create_childs(self, context, o):

        n_childs = 0
        for child in o.children:
            if ARCHIPACK_PT_door_panel.filter(child):
                n_childs += 1

        # remove child
        if n_childs > self.n_panels:
            self.remove_childs(context, o, n_childs - self.n_panels)

        if n_childs < 1:
            # create one door panel
            bpy.ops.archipack.door_panel(x=self.x, z=self.z, door_y=self.door_y,
                    n_panels=self.n_panels, direction=self.direction)
            child = context.object
            child.parent = o
            child.matrix_world = o.matrix_world.copy()
            location = self.x / 2 + BATTUE - SPACING
            if self.direction == 0:
                location = -location
            child.location.x = location
            child.location.y = self.door_y

        if self.n_panels == 2 and n_childs < 2:
            # create 2nth door panel
            bpy.ops.archipack.door_panel(x=self.x, z=self.z, door_y=self.door_y,
                    n_panels=self.n_panels, direction=1 - self.direction)
            child = context.object
            child.parent = o
            child.matrix_world = o.matrix_world.copy()
            location = self.x / 2 + BATTUE - SPACING
            if self.direction == 1:
                location = -location
            child.location.x = location
            child.location.y = self.door_y

    def find_handle(self, o):
        for handle in o.children:
            if 'archipack_handle' in handle:
                return handle
        return None

    def get_childs_panels(self, context, o):
        childs = []
        for child in o.children:
            c, props = ARCHIPACK_PT_door_panel.params(child)
            if props is not None:
                childs.append(c)
        return childs

    def _synch_childs(self, context, o, linked, childs):
        """
            sub synch childs nodes of linked object
        """
        # remove childs not found on source
        l_childs = self.get_childs_panels(context, linked)
        c_names = [c.data.name for c in childs]
        for c in l_childs:
            try:
                id = c_names.index(c.data.name)
            except:
                self.remove_handle(context, c)
                context.scene.objects.unlink(c)
                bpy.data.objects.remove(c, do_unlink=True)

        # children ordering may not be the same, so get the right l_childs order
        l_childs = self.get_childs_panels(context, linked)
        l_names = [c.data.name for c in l_childs]
        order = []
        for c in childs:
            try:
                id = l_names.index(c.data.name)
            except:
                id = -1
            order.append(id)

        # add missing childs and update other ones
        for i, child in enumerate(childs):
            if order[i] < 0:
                p = bpy.data.objects.new("DoorPanel", child.data)
                context.scene.objects.link(p)
                p.lock_location[0] = True
                p.lock_location[1] = True
                p.lock_location[2] = True
                p.lock_rotation[0] = True
                p.lock_rotation[1] = True
                p.lock_scale[0] = True
                p.lock_scale[1] = True
                p.lock_scale[2] = True
                p.parent = linked
                p.matrix_world = linked.matrix_world.copy()
                p.location = child.location.copy()
            else:
                p = l_childs[order[i]]

            p.location = child.location.copy()

            # update handle
            handle = self.find_handle(child)
            h = self.find_handle(p)
            if handle is not None:
                if h is None:
                    h = create_handle(context, p, handle.data)
                    MaterialUtils.add_handle_materials(h)
                h.location = handle.location.copy()
            elif h is not None:
                context.scene.objects.unlink(h)
                bpy.data.objects.remove(h, do_unlink=True)

    def _synch_hole(self, context, linked, hole):
        l_hole = self.find_hole(linked)
        if l_hole is None:
            l_hole = bpy.data.objects.new("hole", hole.data)
            l_hole['archipack_hole'] = True
            context.scene.objects.link(l_hole)
            l_hole.parent = linked
            l_hole.matrix_world = linked.matrix_world.copy()
            l_hole.location = hole.location.copy()
        else:
            l_hole.data = hole.data

    def synch_childs(self, context, o):
        """
            synch childs nodes of linked objects
        """
        bpy.ops.object.select_all(action='DESELECT')
        o.select = True
        context.scene.objects.active = o
        childs = self.get_childs_panels(context, o)
        hole = self.find_hole(o)
        bpy.ops.object.select_linked(type='OBDATA')
        for linked in context.selected_objects:
            if linked != o:
                self._synch_childs(context, o, linked, childs)
                if hole is not None:
                    self._synch_hole(context, linked, hole)

    def update_childs(self, context, o):
        """
            pass params to childrens
        """
        childs = self.get_childs_panels(context, o)
        n_childs = len(childs)
        self.remove_childs(context, o, n_childs - self.n_panels)

        childs = self.get_childs_panels(context, o)
        n_childs = len(childs)
        child_n = 0

        # location_y = self.y / 2 + self.frame_y - SPACING
        # location_y = min(max(self.door_offset, - location_y), location_y) + self.door_y

        location_y = max(0.25 * self.door_y + 0.0005, self.y / 2 + self.frame_y)
        location_y = max(location_y - self.door_offset + 0.5 * self.door_y, -location_y + self.door_y + 0.001)

        x = self.x / self.n_panels + (3 - self.n_panels) * (BATTUE - SPACING)
        y = self.door_y
        z = self.z + BATTUE - SPACING

        if self.n_panels < 2:
            direction = self.direction
        else:
            direction = 0

        for panel in range(self.n_panels):
            child_n += 1

            if child_n == 1:
                handle = self.handle
            else:
                handle = 'NONE'

            if child_n > 1:
                direction = 1 - direction

            location_x = (2 * direction - 1) * (self.x / 2 + BATTUE - SPACING)

            if child_n > n_childs:
                bpy.ops.archipack.door_panel(
                    x=x,
                    y=y,
                    z=z,
                    model=self.model,
                    direction=direction,
                    chanfer=self.chanfer,
                    panel_border=self.panel_border,
                    panel_bottom=self.panel_bottom,
                    panel_spacing=self.panel_spacing,
                    panels_distrib=self.panels_distrib,
                    panels_x=self.panels_x,
                    panels_y=self.panels_y,
                    handle=handle
                    )
                child = context.object
                # parenting at 0, 0, 0 before set object matrix_world
                # so location remains local from frame
                child.parent = o
                child.matrix_world = o.matrix_world.copy()
            else:
                child = childs[child_n - 1]
                child.select = True
                context.scene.objects.active = child
                c, props = ARCHIPACK_PT_door_panel.params(child)
                if props is not None:
                    props.x = x
                    props.y = y
                    props.z = z
                    props.model = self.model
                    props.direction = direction
                    props.chanfer = self.chanfer
                    props.panel_border = self.panel_border
                    props.panel_bottom = self.panel_bottom
                    props.panel_spacing = self.panel_spacing
                    props.panels_distrib = self.panels_distrib
                    props.panels_x = self.panels_x
                    props.panels_y = self.panels_y
                    props.handle = handle
                    props.update(context)
            child.location = Vector((location_x, location_y, 0))

    def update(self, context, childs_only=False):

        # support for "copy to selected"
        active, selected, o = self.find_in_selection(context)

        if o is None or not self.auto_update:
            return

        if childs_only is False:
            bmed.buildmesh(context, o, self.verts, self.faces, self.matids, self.uvs)

        self.update_childs(context, o)

        if childs_only is False and self.find_hole(o) is not None:
            self.interactive_hole(context, o)

        # support for instances childs, update at object level
        self.synch_childs(context, o)

        # setup 3d points for gl manipulators
        x, y = 0.5 * self.x, 0.5 * self.y
        self.manipulators[0].set_pts([(-x, -y, 0), (x, -y, 0), (1, 0, 0)])
        self.manipulators[1].set_pts([(-x, -y, 0), (-x, y, 0), (-1, 0, 0)])
        self.manipulators[2].set_pts([(x, -y, 0), (x, -y, self.z), (-1, 0, 0)])

        # restore context
        bpy.ops.object.select_all(action="DESELECT")

        try:
            for o in selected:
                o.select = True
        except:
            pass

        active.select = True
        context.scene.objects.active = active

    def find_hole(self, o):
        for child in o.children:
            if 'archipack_hole' in child:
                return child
                return None

    def interactive_hole(self, context, o):
        hole_obj = self.find_hole(o)
        if hole_obj is None:
            m = bpy.data.meshes.new("hole")
            hole_obj = bpy.data.objects.new("hole", m)
            context.scene.objects.link(hole_obj)
            hole_obj['archipack_hole'] = True
            hole_obj.parent = o
            hole_obj.matrix_world = o.matrix_world.copy()
            MaterialUtils.add_wall_materials(hole_obj)
        hole = self.hole
        v = Vector((0, 0, 0))
        offset = Vector((0, -0.001, 0))
        size = Vector((self.x + 2 * self.frame_x, self.z + self.frame_x + 0.001, self.y))
        verts = hole.vertices(16, offset, v, v, size, v, 0, 0, shape_z=None, path_type='RECTANGLE')
        faces = hole.faces(16, path_type='RECTANGLE')
        matids = hole.mat(16, 0, 1, path_type='RECTANGLE')
        uvs = hole.uv(16, v, v, size, v, 0, 0, 0, 0, path_type='RECTANGLE')
        bmed.buildmesh(context, hole_obj, verts, faces, matids=matids, uvs=uvs)
        return hole_obj

    def robust_hole(self, context, tM):
        hole = self.hole
        m = bpy.data.meshes.new("hole")
        o = bpy.data.objects.new("hole", m)
        o['archipack_robusthole'] = True
        context.scene.objects.link(o)
        v = Vector((0, 0, 0))
        offset = Vector((0, -0.001, 0))
        size = Vector((self.x + 2 * self.frame_x, self.z + self.frame_x + 0.001, self.y))
        verts = hole.vertices(16, offset, v, v, size, v, 0, 0, shape_z=None, path_type='RECTANGLE')
        verts = [tM * Vector(v) for v in verts]
        faces = hole.faces(16, path_type='RECTANGLE')
        matids = hole.mat(16, 0, 1, path_type='RECTANGLE')
        uvs = hole.uv(16, v, v, size, v, 0, 0, 0, 0, path_type='RECTANGLE')
        bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=uvs)
        MaterialUtils.add_wall_materials(o)
        o.select = True
        context.scene.objects.active = o
        return o


class ARCHIPACK_PT_door(Panel):
    bl_idname = "ARCHIPACK_PT_door"
    bl_label = "Door"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    def draw(self, context):
        o = context.object
        if not ARCHIPACK_PT_door.filter(o):
            return
        layout = self.layout
        layout.operator('archipack.door_manipulate', icon='HAND')
        props = o.data.archipack_door[0]
        row = layout.row(align=True)
        row.operator('archipack.door', text="Refresh", icon='FILE_REFRESH').mode = 'REFRESH'
        if o.data.users > 1:
            row.operator('archipack.door', text="Make unique", icon='UNLINKED').mode = 'UNIQUE'
        row.operator('archipack.door', text="Delete", icon='ERROR').mode = 'DELETE'
        box = layout.box()
        # box.label(text="Styles")
        row = box.row(align=True)
        row.menu("ARCHIPACK_MT_door_preset", text=bpy.types.ARCHIPACK_MT_door_preset.bl_label)
        row.operator("archipack.door_preset", text="", icon='ZOOMIN')
        row.operator("archipack.door_preset", text="", icon='ZOOMOUT').remove_active = True
        row = layout.row()
        box = row.box()
        box.label(text="Size")
        box.prop(props, 'x')
        box.prop(props, 'y')
        box.prop(props, 'z')
        box.prop(props, 'door_offset')
        row = layout.row()
        box = row.box()
        row = box.row()
        row.label(text="Door")
        box.prop(props, 'direction')
        box.prop(props, 'n_panels')
        box.prop(props, 'door_y')
        box.prop(props, 'handle')
        row = layout.row()
        box = row.box()
        row = box.row()
        row.label(text="Frame")
        row = box.row(align=True)
        row.prop(props, 'frame_x')
        row.prop(props, 'frame_y')
        row = layout.row()
        box = row.box()
        row = box.row()
        row.label(text="Panels")
        box.prop(props, 'model')
        if props.model > 0:
            box.prop(props, 'panels_distrib', text="")
            row = box.row(align=True)
            row.prop(props, 'panels_x')
            if props.panels_distrib == 'REGULAR':
                row.prop(props, 'panels_y')
            box.prop(props, 'panel_bottom')
            box.prop(props, 'panel_spacing')
            box.prop(props, 'panel_border')
            box.prop(props, 'chanfer')

    @classmethod
    def params(cls, o):
        if cls.filter(o):
            return o.data.archipack_door[0]
        return None

    @classmethod
    def filter(cls, o):
        try:
            return bool('archipack_door' in o.data)
        except:
            return False

    @classmethod
    def poll(cls, context):
        o = context.object
        if o is None:
            return False
        return cls.filter(o)

# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_door(Operator):
    bl_idname = "archipack.door"
    bl_label = "Door"
    bl_description = "Door"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    x = FloatProperty(
            name='width',
            min=0.1, max=10000,
            default=0.80, precision=2,
            description='Width'
            )
    y = FloatProperty(
            name='depth',
            min=0.1, max=10000,
            default=0.20, precision=2,
            description='Depth'
            )
    z = FloatProperty(
            name='height',
            min=0.1, max=10000,
            default=2.0, precision=2,
            description='height'
            )
    direction = IntProperty(
            name="direction",
            min=0,
            max=1,
            description="open direction"
            )
    n_panels = IntProperty(
            name="panels",
            min=1,
            max=2,
            default=1,
            description="number of panels"
            )
    chanfer = FloatProperty(
            name='chanfer',
            min=0.001, max=10000,
            default=0.005, precision=3,
            description='chanfer'
            )
    panel_spacing = FloatProperty(
            name='spacing',
            min=0.001, max=10000,
            default=0.1, precision=2,
            description='distance between panels'
            )
    panel_bottom = FloatProperty(
            name='bottom',
            min=0.0, max=10000,
            default=0.0, precision=2,
            description='distance from bottom'
            )
    panel_border = FloatProperty(
            name='border',
            min=0.001, max=10000,
            default=0.2, precision=2,
            description='distance from border'
            )
    panels_x = IntProperty(
            name="panels h",
            min=1,
            max=50,
            default=1,
            description="panels h"
            )
    panels_y = IntProperty(
            name="panels v",
            min=1,
            max=50,
            default=1,
            description="panels v"
            )
    panels_distrib = EnumProperty(
            name='distribution',
            items=(
                ('REGULAR', 'Regular', '', 0),
                ('ONE_THIRD', '1/3 2/3', '', 1)
                ),
            default='REGULAR'
            )
    handle = EnumProperty(
            name='Shape',
            items=(
                ('NONE', 'No handle', '', 0),
                ('BOTH', 'Inside and outside', '', 1)
                ),
            default='BOTH'
            )
    mode = EnumProperty(
            items=(
            ('CREATE', 'Create', '', 0),
            ('DELETE', 'Delete', '', 1),
            ('REFRESH', 'Refresh', '', 2),
            ('UNIQUE', 'Make unique', '', 3),
            ),
            default='CREATE'
            )
    auto_manipulate = BoolProperty(
            default=True
            )

    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def create(self, context):
        """
            expose only basic params in operator
            use object property for other params
        """
        m = bpy.data.meshes.new("Frame")
        o = bpy.data.objects.new("Frame", m)
        d = m.archipack_door.add()
        d.x = self.x
        d.y = self.y
        d.z = self.z
        d.direction = self.direction
        d.n_panels = self.n_panels
        d.chanfer = self.chanfer
        d.panel_border = self.panel_border
        d.panel_bottom = self.panel_bottom
        d.panel_spacing = self.panel_spacing
        d.panels_distrib = self.panels_distrib
        d.panels_x = self.panels_x
        d.panels_y = self.panels_y
        d.handle = self.handle
        s = d.manipulators.add()
        s.prop1_name = "x"
        s.prop2_name = "x"
        s.type = "SIZE_LOC"
        s = d.manipulators.add()
        s.prop1_name = "y"
        s.prop2_name = "y"
        s.type = "SIZE_LOC"
        s = d.manipulators.add()
        s.prop1_name = "z"
        s.normal = Vector((0, 1, 0))
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        d.update(context)
        MaterialUtils.add_door_materials(o)
        o.select = True
        context.scene.objects.active = o
        return o

    def delete(self, context):
        o = context.active_object
        if o.data is not None and 'archipack_door' in o.data:
            for child in o.children:
                if 'archipack_hole' in child:
                    context.scene.objects.unlink(child)
                    bpy.data.objects.remove(child, do_unlink=True)
                elif child.data is not None and 'archipack_doorpanel' in child.data:
                    for handle in child.children:
                        if 'archipack_handle' in handle:
                            context.scene.objects.unlink(handle)
                            bpy.data.objects.remove(handle, do_unlink=True)
                    context.scene.objects.unlink(child)
                    bpy.data.objects.remove(child, do_unlink=True)
            context.scene.objects.unlink(o)
            bpy.data.objects.remove(o, do_unlink=True)

    def update(self, context):
        o = context.active_object
        if o.data is not None and 'archipack_door' in o.data:
            o.data.archipack_door[0].update(context)

    def unique(self, context):
        sel = [o for o in context.selected_objects]
        bpy.ops.object.select_all(action="DESELECT")
        for o in sel:
            if o.data is not None and 'archipack_door' in o.data:
                o.select = True
                for child in o.children:
                    if 'archipack_hole' in child or (child.data is not None and
                       'archipack_doorpanel' in child.data):
                        child.select = True
        if len(context.selected_objects) > 0:
            bpy.ops.object.make_single_user(type='SELECTED_OBJECTS', object=True,
                obdata=True, material=False, texture=False, animation=False)

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            if self.mode == 'CREATE':
                bpy.ops.object.select_all(action="DESELECT")
                o = self.create(context)
                o.location = bpy.context.scene.cursor_location
                if self.auto_manipulate:
                    bpy.ops.archipack.door_manipulate('INVOKE_DEFAULT')
            elif self.mode == 'DELETE':
                self.delete(context)
            elif self.mode == 'REFRESH':
                self.update(context)
            elif self.mode == 'UNIQUE':
                self.unique(context)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}

# ------------------------------------------------------------------
# Define operator class to manipulate object
# ------------------------------------------------------------------


class ARCHIPACK_OT_door_manipulate(Operator):
    bl_idname = "archipack.door_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return ARCHIPACK_PT_door.filter(context.active_object)

    def modal(self, context, event):
        return self.d.manipulable_modal(context, event)

    def invoke(self, context, event):
        if context.space_data.type == 'VIEW_3D':
            o = context.active_object
            self.d = o.data.archipack_door[0]
            if self.d.manipulable_invoke(context):
                context.window_manager.modal_handler_add(self)
                return {'RUNNING_MODAL'}
            else:
                return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to load / save presets
# ------------------------------------------------------------------


class ARCHIPACK_MT_door_preset(Menu):
    bl_label = "Door Styles"
    preset_subdir = "archipack_door" # you might wanna change this
    preset_operator = "script.execute_preset" # but not this
    draw = Menu.draw_preset # or that


class ARCHIPACK_OT_door_preset(ArchipackPreset, Operator):
    """Add a Door Styles"""
    bl_idname = "archipack.door_preset"
    bl_label = "Add Door Style"
    preset_menu = "ARCHIPACK_MT_door_preset"

    datablock_name = StringProperty(
        name="Datablock",
        default='archipack_door',
        maxlen=64,
        options={'HIDDEN', 'SKIP_SAVE'},
        )

    @property
    def blacklist(self):
        return ['x', 'y', 'z', 'direction', 'manipulators']


bpy.utils.register_class(archipack_door)
Mesh.archipack_door = CollectionProperty(type=archipack_door)
bpy.utils.register_class(ARCHIPACK_MT_door_preset)
bpy.utils.register_class(ARCHIPACK_PT_door)
bpy.utils.register_class(ARCHIPACK_OT_door)
bpy.utils.register_class(ARCHIPACK_OT_door_preset)
bpy.utils.register_class(ARCHIPACK_OT_door_manipulate)
