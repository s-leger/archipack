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
import os
import json
import subprocess
from math import atan2, pi
from mathutils.geometry import interpolate_bezier
from mathutils import Matrix, Vector
from bpy.types import Operator
from bpy.props import (
    StringProperty, EnumProperty, BoolProperty,
    FloatVectorProperty, IntProperty
    )
from .archipack_viewmanager import ViewManager
from .archipack_iconmanager import icons
from .archipack_object import ArchipackGenericOperator, ArchipackObjectsManager


precision = 6


"""Generate presets from selected curves in a files

sel = C.selected_objects
for o in sel:
    bpy.ops.archipack.profile_save(
        'INVOKE_DEFAULT',
        curve_name=o.name,
        preset_name="Molding_{}".format(o.name.capitalize())
        )
"""


def get_enum(self, context):
    return icons.enum(context, "archipack_curves")


def preset_operator(self, context):
    if self.auto_update:
        act = context.active_object
        sel = context.selected_objects
        print("preset_operator:", act.name)

        self.user_profile_savename = bpy.path.display_name(self.user_profile_filename)
        self.load_profile(context, update=True)

        cls = self.__class__.__name__
        x, y, z = self.user_profile_dimension
        for o in sel:
            if o.data and o.name != act.name:
                if cls in o.data:
                    ArchipackObjectsManager.select_object(self, context, o, True)
                    d = getattr(o.data, cls)[0]
                    d.auto_update = False
                    # enable user profile
                    d.profil = self.profil
                    d.user_profile_dimension = self.user_profile_dimension.copy()
                    d.user_profile_savename = self.user_profile_savename
                    d.user_profile_filename = self.user_profile_filename
                    # profile curve object name
                    d.user_profile = self.user_profile
                    d.refresh_profile_size(context, x, y)
                    d.auto_update = True

        ArchipackObjectsManager.select_object(self, context, act, True)

    return None


class ArchipackCurveManager():
    """
     Provide utility to manage curve input
    """
    def is_cw(self, pts):
        p0 = pts[0]
        d = 0
        for p in pts[1:]:
            d += (p.x * p0.y - p.y * p0.x)
            p0 = p
        return d > 0

    def interpolate_bezier(self, pts, wM, p0, p1, resolution):
        """
         Bezier curve approximation
        """
        if (resolution == 0 or
                (p0.handle_right_type == 'VECTOR' and
                p1.handle_left_type == 'VECTOR')):
            pts.append(wM * p0.co.to_3d())
        else:
            v = (p1.co - p0.co).normalized()
            d1 = (p0.handle_right - p0.co).normalized()
            d2 = (p1.co - p1.handle_left).normalized()
            if d1 == v and d2 == v:
                pts.append(wM * p0.co.to_3d())
            else:
                seg = interpolate_bezier(wM * p0.co,
                    wM * p0.handle_right,
                    wM * p1.handle_left,
                    wM * p1.co,
                    resolution + 1)
                for i in range(resolution):
                    pts.append(seg[i].to_3d())

    def coords_from_spline(self, spline, wM, resolution, ccw=False, cw=False, close=False):
        """
            Return coords from spline
            Explicitely closed: first coord = last coord
            wM: matrix to make points absolute world
            resolution: bezier resolution
            ccw: return points in ccw order
            cw: return points in cw order
            close: force closed spline
        """
        pts = []
        if spline.type == 'POLY':

            pts = [wM * p.co.to_3d() for p in spline.points]

            if close or spline.use_cyclic_u:
                pts.append(pts[0])

        elif spline.type == 'BEZIER':

            points = spline.bezier_points

            for i in range(1, len(points)):
                p0 = points[i - 1]
                p1 = points[i]
                self.interpolate_bezier(pts, wM, p0, p1, resolution)

            if close or spline.use_cyclic_u:
                p0 = points[-1]
                p1 = points[0]
                self.interpolate_bezier(pts, wM, p0, p1, resolution)
                pts.append(pts[0])

            else:
                pts.append(wM * points[-1].co)

        if ccw or cw:
            is_cw = self.is_cw(pts)
            if (ccw and is_cw) or (cw and not is_cw):
                pts = list(reversed(pts))

        if len(pts) < 1:
            return []

        # skip duplicate points using arbitrary set distance threshold
        # to handle precision issues
        new_pts = [pts[0]]
        for p in pts:
            if (p - new_pts[-1]).length < 0.0001:
                continue
            new_pts.append(p)

        return new_pts


def update_path(self, context):
    self.update_path(context)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


class ArchipackUserDefinedPath(ArchipackCurveManager):
    """
      Archipack path based base class
      base for wall, slab, floor, fence, molding and cutters
      provide "from curve" and parts management (insert/remove)
    """

    """
     Origin allow to move first point of an object
     without actually changing the pivot location
     currently for tests only, might break T childs walls and roofs
    """
    origin = FloatVectorProperty(
            name="Origin",
            description="Object origin in local space (offset from pivot)",
            subtype="XYZ"
            )

    user_defined_path = StringProperty(
            description="Use a curve to define shape",
            name="User defined",
            update=update_path
            )

    user_defined_reverse = BoolProperty(
            name="Reverse",
            default=False,
            update=update_path
            )

    user_defined_spline = IntProperty(
            name="Spline index",
            min=0,
            default=0,
            update=update_path
            )

    user_defined_resolution = IntProperty(
            name="Resolution",
            min=1,
            max=128,
            default=12,
            update=update_path
            )

    n_parts = IntProperty(
            name="Parts",
            min=1,
            default=1,
            update=update_manipulators
            )

    parts_expand = BoolProperty(
            description="Expand Segments panel",
            options={'SKIP_SAVE'},
            default=False
            )
    user_defined_path_expand = BoolProperty(
            description="Expand from curve panel",
            name="Expand",
            options={'SKIP_SAVE'},
            default=False
            )

    def move_object(self, o, p):
        """
         When firstpoint is moving we must move object according
         p is new x, y location in world coordsys
        """
        p = Vector((p.x, p.y, o.matrix_world.translation.z))
        # p is in o coordsys
        if o.parent:
            o.location = p * o.parent.matrix_world.inverted()
            o.matrix_world.translation = p
            """
            # Update linked objects location too
            act = context.active_object
            sel = context.selected_objects[:]
            bpy.ops.object.select_all(action="DESELECT")
            self.select_object(context, o, True)
            bpy.ops.object.select_linked()
            self.unselect_object(o)
            for c in context.selected_objects:
                c.location = p * c.parent.matrix_world.inverted()
                c.matrix_world.translation = p
            bpy.ops.object.select_all(action="DESELECT")

            # restore selection and active
            for o in sel:
                self.select_object(context, o)
            self.select_object(context, act, True)
            """
        else:
            o.location = p
            o.matrix_world.translation = p

    def from_points(self, pts):
        """
         Make parts from points
        """
        self.n_parts = len(pts) - 1

        self.update_parts()
        p0 = pts.pop(0)
        a0 = 0
        for i, p1 in enumerate(pts):
            if i >= len(self.parts):
                break
            dp = p1 - p0
            da = atan2(dp.y, dp.x) - a0
            # keep da in range -180 | +180 degree
            if da > pi:
                da -= 2 * pi
            if da < -pi:
                da += 2 * pi
            p = self.parts[i]
            p.type = "S_SEG"
            p.length = dp.to_2d().length
            p.dz = dp.z
            p.a0 = da
            a0 += da
            # keep a0 in range -180 | +180 degree
            if a0 > pi:
                a0 -= 2 * pi
            if a0 < -pi:
                a0 += 2 * pi
            p0 = p1

    def relocate_childs(self, context, o):
        """
         Override this method to synch childs when changes are done
        """
        return

    def after_reverse(self, context, o):
        self.auto_update = True

    def reverse(self, context, o):

        g = self.get_generator(o)

        # 2nd check for fences
        if not self.closed and len(g.segs) > self.n_parts:
            g.segs.pop()

        s_type = list(reversed([p.type for p in self.parts]))
        g_segs = [s.oposite for s in reversed(g.segs)]

        self.auto_update = False

        last = None
        for i, s in enumerate(g_segs):
            p = self.parts[i]
            p.type = s_type[i]
            if "C_" in s_type[i]:
                p.radius = s.r
                p.da = s.da
            else:
                p.length = s.length
            if last is None:
                p.a0 = s.a0
            else:
                p.a0 = s.delta_angle(last)
            last = s

        # location wont change when closed
        if not self.closed:
            self.move_object(o, g_segs[0].p0)

        self.after_reverse(context, o)

    def update_path(self, context):
        """
         Handle curve Io
         including multiple splines
        """
        path = self.get_scene_object(context, self.user_defined_path)
        if path is not None and path.type == 'CURVE':
            splines = path.data.splines
            if len(splines) > self.user_defined_spline:
                self.from_spline(
                    context,
                    path.matrix_world,
                    self.user_defined_resolution,
                    splines[self.user_defined_spline])

    def setup_parts_manipulators(self, z_prop='z'):
        """
         z_prop: prop name of dumb or real z size for placeholder
                 when drawing wall snap manipulator
        """
        for i, p in enumerate(self.parts):
            n_manips = len(p.manipulators)
            if n_manips < 1:
                s = p.manipulators.add()
                s.type_key = "ANGLE"
                s.prop1_name = "a_ui"
            if n_manips < 2:
                s = p.manipulators.add()
                s.type_key = "SIZE"
                s.prop1_name = "l_ui"
            if n_manips < 3:
                s = p.manipulators.add()
                s.type_key = 'WALL_SNAP'
                s.prop1_name = str(i)
                s.prop2_name = z_prop
            if n_manips < 4:
                s = p.manipulators.add()
                s.type_key = 'DUMB_STRING'
                s.prop1_name = str(i + 1)

            p.manipulators[2].prop1_name = str(i)
            p.manipulators[2].type_key = 'WALL_SNAP'
            p.manipulators[3].prop1_name = str(i + 1)

    def add_part(self, context, length):
        self.manipulable_disable(context)
        self.auto_update = False
        p = self.parts.add()
        p.length = length
        if not self.always_closed:
            self.parts.move(len(self.parts) - 1, self.n_parts)
            p = self.parts[self.n_parts]
        self.n_parts += 1
        self.setup_manipulators()
        self.auto_update = True
        return p

    def before_insert_part(self, context, o, g):
        return

    def after_insert_part(self, context, o, where, distance):
        """
         Rebuild childs index and location
         When removing a part
         Override in objects synch childs
        """
        return

    def insert_part(self, context, o, where):

        # disable manipulate as we are changing structure
        self.manipulable_disable(context)
        self.auto_update = False

        g = self.get_generator()
        length = g.segs[where].length / 2

        # detect childs location
        self.before_insert_part(context, o, g)

        # the part we do split
        part_0 = self.parts[where]
        typ = part_0.type
        # length = part_0.length / 2
        da = part_0.da / 2

        part_0.length = length
        part_0.da = da

        part_1 = self.parts.add()
        part_1.type = typ
        part_1.length = length
        part_1.da = da
        part_1.a0 = 0

        # move after current one
        self.n_parts += 1
        self.parts.move(len(self.parts) - 1, where + 1)

        self.after_insert_part(context, o, where, length)

        self.setup_manipulators()
        self.auto_update = True

    def before_remove_part(self, context, o, g):
        return

    def after_remove_part(self, context, o, where, distance):
        return

    def remove_part(self, context, o, where):

        if self.n_parts < 2 or where < 1:
            return

        # disable manipulate as we are changing structure
        self.manipulable_disable(context)
        self.auto_update = False

        g = self.get_generator()

        # detect childs location
        self.before_remove_part(context, o, g)

        # preserve shape
        # using generator
        w = g.segs[where - 1]
        part = self.parts[where - 1]

        # length of merged segment
        # to be added to childs of [where]
        length = w.length

        w.p1 = g.segs[where].p1

        if where + 1 < self.n_parts:
            self.parts[where + 1].a0 = g.segs[where + 1].delta_angle(w)

        if "C_" in part.type:
            part.radius = w.r
        else:
            part.length = w.length

        if where > 1:
            part.a0 = w.delta_angle(g.segs[where - 2])
        else:
            part.a0 = w.straight(1, 0).angle

        self.after_remove_part(context, o, where, length)

        self.parts.remove(where)
        self.n_parts -= 1

        # fix snap manipulators index
        self.setup_manipulators()
        self.auto_update = True

    def update_parts(self):

        # For some reason floors dosent follow the + 1 rule ??

        # flag indicate parts change
        changed = False

        # when object is .always_closed
        # n_parts match real parts count
        n_parts = self.n_parts
        if not self.always_closed:
            # otherwhise n_parts is 1 part less than real parts count
            n_parts += 1

        # remove parts
        for i in range(len(self.parts), n_parts, -1):
            changed = True
            self.parts.remove(i - 1)

        # add missing parts
        for i in range(len(self.parts), n_parts):
            changed = True
            self.parts.add()

        for p in self.parts:
            if p.uid == 0:
                self.create_uid(p)

        self.setup_manipulators()

        return changed

    def template_user_path(self, context, layout):
        """
         Draw from curve in ui
        """
        icon = "TRIA_RIGHT"
        if self.user_defined_path_expand:
            icon = "TRIA_DOWN"

        row = layout.row(align=True)
        row.prop(self, 'user_defined_path_expand', icon=icon, text="From curve", icon_only=True, emboss=True)

        if self.user_defined_path_expand:
            row = layout.row(align=True)
            if self.user_defined_path != "":
                op = row.operator("archipack.curve_edit", text="", icon="EDITMODE_HLT")
                op.curve_name = self.user_defined_path
                op.update_func_name = "update_path"
                row.operator(
                    "archipack.object_update",
                    text="",
                    icon='FILE_REFRESH'
                    ).update_func_name = "update_path"
            row.prop_search(self, "user_defined_path", context.scene, "objects", text="", icon='OUTLINER_OB_CURVE')

            if self.user_defined_path != "":
                layout.prop(self, 'user_defined_spline')
                layout.prop(self, 'user_defined_resolution')

        return self.user_defined_path_expand

    def template_parts(self, context, layout, draw_type=False):
        """
         Draw parts in ui
        """
        box = layout.box()
        row = box.row(align=False)

        icon = "TRIA_RIGHT"
        if self.parts_expand:
            icon = "TRIA_DOWN"

        row.prop(self, 'parts_expand', icon=icon, text="Segs", icon_only=True, emboss=True)
        row.prop(self, "n_parts", text="")

        if self.parts_expand:
            if not self.always_closed:
                box.prop(self, "closed")
                box.operator("archipack.path_reverse", icon='FILE_REFRESH')

            n_parts = self.n_parts
            if self.closed or self.always_closed:
                n_parts += 1

            for i, part in enumerate(self.parts):
                if i < n_parts:
                    box = layout.box()
                    part.draw(context, box, i, draw_type=draw_type)


class ARCHIPACK_OT_path_reverse(ArchipackGenericOperator, Operator):
    bl_idname = "archipack.path_reverse"
    bl_label = "Reverse"
    bl_description = "Reverse parts order"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'INTERNAL', 'UNDO'}

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            d = self.datablock(o)
            if d is None:
                return {'CANCELLED'}
            d.reverse(context, o)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ArchipackProfileManager(ArchipackCurveManager):
    """
     IO thumbs and .json files for user defined profiles
    """
    def _round_co(self, co, prec):
        return (round(co.x, prec), round(co.y, prec), round(co.z, prec))

    def _from_json(self, curve, json_str):
        """
          Adds a spline to a curve
        """
        js = json.loads(json_str)
        curve.resolution_u = js['resolution_u']
        curve.fill_mode = js['fill_mode']

        for spl in js['splines']:

            spline = curve.splines.new(spl['type'])
            spline.use_cyclic_u = spl['use_cyclic_u']
            handle_left = spl['handle_left']
            handle_right = spl['handle_right']
            pts = spl['points']

            if spl['type'] == 'BEZIER':
                spline.bezier_points.add(len(pts) - 1)
                for i, p in enumerate(pts):
                    spline.bezier_points[i].co = p
                    spline.bezier_points[i].handle_left = handle_left[i]
                    spline.bezier_points[i].handle_right = handle_right[i]
            elif spl['type'] == 'POLY':
                spline.points.add(len(pts) - 1)
                for i, p in enumerate(pts):
                    spline.bezier_points[i].co = p

        return js['x'], js['y']

    def _to_json(self, curve):
        curve.dimensions = '2D'
        js = {
            'resolution_u': curve.resolution_u,
            'fill_mode': curve.fill_mode,
            'splines': []
            }

        for spline in curve.splines:
            spl = {
                'type': spline.type,
                'use_cyclic_u': spline.use_cyclic_u,
            }
            if spline.type == 'BEZIER':
                pts = spline.bezier_points
                spl['handle_left'] = [self._round_co(p.handle_left, precision) for p in pts]
                spl['handle_right'] = [self._round_co(p.handle_right, precision) for p in pts]

            elif spline.type == 'POLY':
                pts = spline.points

            spl['points'] = [self._round_co(p.co, precision) for p in pts]
            js['splines'].append(spl)

        js['x'], js['y'] = self._curve_bound_box(curve)

        json_str = json.dumps(js, sort_keys=True)
        return json_str

    def _find_preset_full_path(self, preset_name):
        """
         @preset_name: preset clean file name without extension
         Return preset full path with .json extension
        """
        filename = "{}.json".format(preset_name)
        folders = bpy.utils.script_paths("presets")
        # append after user script so user override factory
        folders.append(os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "presets"))
        for folder in folders:
            user_preset = os.path.join(
                folder,
                "archipack_curves",
                filename)
            if os.path.isfile(user_preset):
                return user_preset
        return None

    def _user_preset_full_path(self, preset_name):
        """
         @preset_name: preset clean file name without extension
         Return absolute full path of preset in user presets folder
        """
        preset_file = "{}.json".format(preset_name)
        presets_path = bpy.utils.user_resource('SCRIPTS',
                          os.path.join("presets", "archipack_curves"),
                          create=True)
        return os.path.join(presets_path, preset_file)

    def load_curve(self, preset_name):
        """
          Load a json spline definition
          return curve data not linked to scene
        """
        curve = None
        x, y = 0, 0
        # load from current directory
        full_path = self._find_preset_full_path(preset_name)

        if full_path is None:
            print("Archipack: Error preset file not found")
            return curve, x, y

        try:
            with open(full_path) as fh:
                json_str = ''.join(fh.readlines())
            curve = bpy.data.curves.new(name=preset_name, type='CURVE')
            curve.dimensions = '2D'
            x, y = self._from_json(curve, json_str)
        except:
            print("Archipack: Error while loading {}".format(preset_name))
            pass
        return curve, x, y

    def _curve_bound_box(self, curve):
        # estimate curve size
        pts = []
        for spline in curve.splines:
            pts.extend(self.coords_from_spline(spline, Matrix(), 12))
        x = [co.x for co in pts]
        y = [co.y for co in pts]
        sx = round(max(0.0001, max(x) - min(x)), precision)
        sy = round(max(0.0001, max(y) - min(y)), precision)
        return sx, sy

    def save_curve(self, o, preset_name):
        """
         Save a json curve spline definition
        """
        curve = o.data
        if curve:
            file_name = bpy.path.clean_name(preset_name)
            # always save in user presets
            full_path = self._user_preset_full_path(file_name)
            json_str = self._to_json(curve)
            with open(full_path, 'w') as fh:
                fh.write(json_str)
            # clear file when present
            thumb_path = full_path[:-5] + '.png'
            if os.path.isfile(thumb_path):
                os.remove(thumb_path)
            # update profiles enum
            icons.add_preset("archipack_curves", file_name)
            self.background_render(thumb_path)

    def background_render(self, preset):
        generator = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + "archipack_thumbs.py"
        addon_name = __name__.split('.')[0]
        # Run external instance of blender like the original thumbnail generator.
        cmd = [
            bpy.app.binary_path,
            "--background",
            "--factory-startup",
            "-noaudio",
            "--python", generator,
            "--",
            "addon:" + addon_name,
            "matlib:none",
            "cls:archipack_curves",
            "preset:" + preset
            ]
        subprocess.Popen(cmd)


def update(self, context):
    if self.auto_update:
        o = ArchipackObjectsManager.get_scene_object(self, context, self.user_profile)
        if o and o.type == 'CURVE':
            self.auto_update = False
            sx, sy = self._curve_bound_box(o.data)
            self.user_profile_dimension.x = sx
            self.user_profile_dimension.y = sy
            self.refresh_profile_size(sx, sy)
            self.auto_update = True
            self.update(context)
        else:
            self.user_profile = ""


class ArchipackProfile(ArchipackProfileManager):
    """
      Propertygroup to manage user profiles
      Must inherit ArchipackObjectsManager
    """
    user_profile_filename = EnumProperty(
            name="Profiles",
            description="Available profiles presets",
            default=None,
            items=get_enum,
            update=preset_operator
            )
    user_profile_savename = StringProperty(
            options={'SKIP_SAVE'},
            default="",
            description="Profile filename without extension",
            name="Save"
            )
    user_profile = StringProperty(
            name="Profile",
            description="Use curve as profile",
            default="",
            update=update
            )
    user_profile_dimension = FloatVectorProperty(subtype='XYZ')
    profile_expand = BoolProperty(
        options={'SKIP_SAVE'},
        default=False,
        description="Expand panel display",
        name="Expand"
        )

    def refresh_profile_size(self, context, x, y):
        """
         Todo: override this method
         to update profile size ui
         eg:
         self.x = x
         self.z = y

        """
        return

    def load_profile(self, context, update=False):
        if self.user_profile_filename is None:
            return None
        o = self.get_scene_object(context, self.user_profile)
        curve, x, y = self.load_curve(self.user_profile_filename)
        if curve:
            if o is None:
                o = bpy.data.objects.new(curve.name, curve)
                o.show_name = True
                self.link_object_to_scene(context, o, layer_name="profile")
                self.user_profile_savename = o.name
            else:
                d = o.data
                o.name = curve.name
                o.data = curve
                if d.users < 1:
                    bpy.data.curves.remove(d)

            self.auto_update = False
            self.user_profile_dimension = Vector((x, y, 0))
            self.user_profile = o.name
            # Sub objects without auto_update must update parent using this call
            self.refresh_profile_size(context, x, y)
            self.auto_update = True
        return o

    def update_profile(self, context):
        if self.user_profile == "":
            return
        o = self.get_scene_object(context, self.user_profile)
        if o is None:
            o = self.load_profile(context)
        return o

    def draw_user_profile(self, context, layout, update_func_name="update"):
        # layout.prop(self, 'user_profile_filename')

        layout.template_icon_view(self, "user_profile_filename", show_labels=True, scale=5)
        row = layout.row(align=True)
        if self.user_profile != "":
            op = row.operator("archipack.curve_edit", text="", icon="EDITMODE_HLT")
            op.curve_name = self.user_profile
            op.update_func_name = update_func_name
            row.operator("archipack.object_update", text="", icon="FILE_REFRESH")
        row.prop_search(self, "user_profile", context.scene, "objects", text="", icon='OUTLINER_OB_CURVE')

        if self.user_profile != "":
            row = layout.row(align=True)
            row.prop(self, 'user_profile_savename')
            op = row.operator("archipack.profile_save", text="", icon="SAVE_COPY")
            op.curve_name = self.user_profile
            op.preset_name = self.user_profile_savename


class ARCHIPACK_OT_object_update(ArchipackGenericOperator, Operator):
    bl_idname = "archipack.object_update"
    bl_label = "Update"
    bl_description = "Update object"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    update_func_name = StringProperty(
        description="Update function name",
        default="update"
        )

    def update(self, context, o):
        last_state = self.is_selected(o)
        self.select_object(context, o, True)
        d = self.datablock(o)
        try:
            getattr(d, self.update_func_name)(context)
        except:
            pass
        if last_state:
            self.select_object(context, o)
        else:
            self.unselect_object(o)

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            sel = context.selected_objects
            for c in sel:
                self.update(context, c)
            self.select_object(context, o, True)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_curve_edit(ArchipackGenericOperator, Operator):
    bl_idname = "archipack.curve_edit"
    bl_label = "Edit"
    bl_description = "Edit a curve"
    bl_options = {'REGISTER', 'UNDO', 'INTERNAL'}

    update_func_name = StringProperty(
        description="Update function name",
        default="update"
        )
    curve_name = StringProperty(
        description="Curve name",
        default=""
        )

    object_name = None
    wm = None

    def modal(self, context, event):

        if event.type == 'ESC':
            bpy.ops.object.mode_set(mode="OBJECT")

        if context.mode != 'EDIT_CURVE':
            o = self.get_scene_object(context, self.object_name)
            if o and self.update_func_name != "":
                d = self.datablock(o)
                if d:
                    try:
                        self.select_object(context, o, True)
                        getattr(d, self.update_func_name)(context)
                    except:
                        pass
            self.wm.restore()
            return {'FINISHED'}

        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        if context.space_data is not None and context.space_data.type == 'VIEW_3D':
            self.object_name = context.active_object.name
            o = self.get_scene_object(context, self.curve_name)
            if o and o.type == 'CURVE':
                self.wm = ViewManager(context)
                self.wm.save()
                self.wm.safe_ortho_view(max(o.dimensions.x, o.dimensions.y), o.matrix_world.translation)
                bpy.ops.object.select_all(action="DESELECT")
                self.select_object(context, o, True)
                bpy.ops.object.mode_set(mode="EDIT")
                context.window_manager.modal_handler_add(self)
                return {'RUNNING_MODAL'}
            else:
                self.report({'WARNING'}, "Profile not found")
                return {'CANCELLED'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}


class ARCHIPACK_OT_profile_save(ArchipackObjectsManager, Operator):
    bl_idname = "archipack.profile_save"
    bl_label = "Save"
    bl_description = "Save a profile curve"
    bl_options = {'REGISTER', 'UNDO', 'INTERNAL'}

    curve_name = StringProperty(
        description="Curve object name",
        default=""
        )
    preset_name = StringProperty(
        description="Preset file name without extension",
        default=""
        )

    @classmethod
    def poll(self, context):
        return True

    def invoke(self, context, event):
        # if context.space_data is not None and context.space_data.type == 'VIEW_3D':
        o = self.get_scene_object(context, self.curve_name)
        if o and o.type == 'CURVE':
            manager = ArchipackProfileManager()
            manager.save_curve(o, self.preset_name)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Invalid curve")
            return {'CANCELLED'}


def register():
    bpy.utils.register_class(ARCHIPACK_OT_path_reverse)
    bpy.utils.register_class(ARCHIPACK_OT_curve_edit)
    bpy.utils.register_class(ARCHIPACK_OT_object_update)
    bpy.utils.register_class(ARCHIPACK_OT_profile_save)


def unregister():
    bpy.utils.unregister_class(ARCHIPACK_OT_path_reverse)
    bpy.utils.unregister_class(ARCHIPACK_OT_curve_edit)
    bpy.utils.unregister_class(ARCHIPACK_OT_object_update)
    bpy.utils.unregister_class(ARCHIPACK_OT_profile_save)


"""Test load/save
from archipack.archipack_curveman import ArchipackCurveManager, ArchipackProfileManager

cman = ArchipackCurveManager()

curve = C.object.data
cman.save_curve(curve, "test")

curve = cman.load_curve("test")
o = bpy.data.objects.new(curve.name, curve)
C.scene.objects.link(o)

"""
