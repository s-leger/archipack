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
from mathutils.geometry import interpolate_bezier
from mathutils import Matrix, Vector
from bpy.types import Operator
from bpy.props import (
    StringProperty, EnumProperty,
    FloatVectorProperty, IntProperty
    )
from .archipack_viewmanager import ViewManager
from .archipack_iconmanager import icons


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

        self.user_profile_savename = bpy.path.display_name(self.user_profile_filename)
        self.load_profile(context, update=True)

        act = context.active_object
        sel = context.selected_objects
        cls = self.__class__.__name__
        x, y, z = self.user_profile_dimension
        for o in sel:
            if o.data and o.name != act.name:
                if cls in o.data:
                    context.scene.objects.active = o
                    d = getattr(o.data, cls)[0]
                    d.auto_update = False
                    # enable user profile
                    d.profil = self.profil
                    d.user_profile_dimension = self.user_profile_dimension.copy()
                    d.user_profile_savename = self.user_profile_savename
                    d.user_profile_filename = self.user_profile_filename
                    # profile curve object name
                    d.user_profile = self.user_profile
                    d.refresh_profile_size(x, y)
                    d.auto_update = True
                    d.update(context)
        context.scene.objects.active = act

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

        return pts


def update_path(self, context):
    self.update_path(context)


class ArchipackUserDefinedPath(ArchipackCurveManager):
    """
      User defined path shared properties
      and operators
    """
    user_defined_path = StringProperty(
            description="Use a curve to define shape",
            name="User defined",
            update=update_path
            )
    user_defined_resolution = IntProperty(
            name="Resolution",
            min=1,
            max=128,
            default=12,
            update=update_path
            )

    def draw_user_path(self, layout, context):
        layout.label(text="From curve")
        row = layout.row(align=True)
        if self.user_defined_path != "":
            op = row.operator("archipack.profile_edit", text="", icon="EDITMODE_HLT")
            op.curve_name = self.user_defined_path
            op.update_func_name = "update_path"
            row.operator(
                "archipack.object_update",
                text="",
                icon='FILE_REFRESH'
                ).update_func_name = "update_path"
        row.prop_search(self, "user_defined_path", context.scene, "objects", text="", icon='OUTLINER_OB_CURVE')

        if self.user_defined_path is not "":
            layout.prop(self, 'user_defined_resolution')


class ArchipackProfileManager(ArchipackCurveManager):
    """
     IO for user defined profiles
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
    o = context.scene.objects.get(self.user_profile)
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
    """
    user_profile_filename = EnumProperty(
        name="Profiles",
        description="Available profiles presets",
        default=None,
        items=get_enum,
        update=preset_operator
        )
    user_profile = StringProperty(
        default="",
        description="Use curve as profile",
        name="Profile",
        update=update
        )
    user_profile_dimension = FloatVectorProperty(subtype='XYZ')

    user_profile_savename = StringProperty(
        options={'SKIP_SAVE'},
        default="",
        description="Profile filename without extension",
        name="Save"
        )

    def refresh_profile_size(self, x, y):
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
        o = context.scene.objects.get(self.user_profile)
        curve, x, y = self.load_curve(self.user_profile_filename)
        if curve:
            if o is None:
                o = bpy.data.objects.new(curve.name, curve)
                o.show_name = True
                context.scene.objects.link(o)
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
            self.refresh_profile_size(x, y)
            self.auto_update = True
            if update:
                self.update(context)
        return o

    def update_profile(self, context):
        if self.user_profile == "":
            return
        o = context.scene.objects.get(self.user_profile)
        if o is None:
            o = self.load_profile(context)
        return o

    def draw_user_profile(self, context, layout):
        # layout.prop(self, 'user_profile_filename')
        layout.template_icon_view(self, "user_profile_filename", show_labels=True, scale=5)
        row = layout.row(align=True)
        if self.user_profile != "":
            row.operator("archipack.profile_edit", text="", icon="EDITMODE_HLT").curve_name = self.user_profile
            row.operator("archipack.object_update", text="", icon="FILE_REFRESH")
        row.prop_search(self, "user_profile", context.scene, "objects", text="", icon='OUTLINER_OB_CURVE')

        if self.user_profile != "":
            row = layout.row(align=True)
            row.prop(self, 'user_profile_savename')
            op = row.operator("archipack.profile_save", text="", icon="SAVE_COPY")
            op.curve_name = self.user_profile
            op.preset_name = self.user_profile_savename


class ARCHIPACK_OT_object_update(Operator):
    bl_idname = "archipack.object_update"
    bl_label = "Update"
    bl_description = "Update object"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    update_func_name = StringProperty(
        description="Update function name",
        default="update"
        )

    @classmethod
    def poll(self, context):
        return context.active_object is not None

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def update(self, context, o, key):
        if o.data and key in o.data:
            context.scene.objects.active = o
            try:
                d = getattr(o.data, key)[0]
                getattr(d, self.update_func_name)(context)
            except:
                pass

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            sel = context.selected_objects
            if o:
                for key in o.data.keys():
                    if "archipack_" in key:
                        self.update(context, o, key)
                        for c in sel:
                            self.update(context, c, key)
                o.select = True
                context.scene.objects.active = o
            else:
                self.report({'WARNING'}, "Object not found")
                return {'CANCELLED'}
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_profile_edit(Operator):
    bl_idname = "archipack.profile_edit"
    bl_label = "Edit"
    bl_description = "Edit a profile curve"
    bl_options = {'REGISTER', 'UNDO', 'INTERNAL'}
    update_func_name = StringProperty(
        description="Update function name",
        default="update"
        )
    curve_name = StringProperty(
        description="Profile curve name",
        default=""
        )

    object_name = None
    wm = None

    @classmethod
    def poll(self, context):
        return context.active_object is not None

    def modal(self, context, event):

        if event.type == 'ESC':
            bpy.ops.object.mode_set(mode="OBJECT")

        if context.mode != 'EDIT_CURVE':
            o = context.scene.objects.get(self.object_name)
            if o and self.update_func_name != "":
                for key in o.data.keys():
                    if "archipack_" in key:
                        bpy.ops.object.select_all(action="DESELECT")
                        o.select = True
                        context.scene.objects.active = o
                        try:
                            d = getattr(o.data, key)[0]
                            getattr(d, self.update_func_name)(context)
                        except:
                            pass
            self.wm.restore()
            return {'FINISHED'}

        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        if context.space_data is not None and context.space_data.type == 'VIEW_3D':
            self.object_name = context.active_object.name
            o = context.scene.objects.get(self.curve_name)
            if o:
                self.wm = ViewManager(context)
                self.wm.save()
                self.wm.safe_ortho_view(max(o.dimensions.x, o.dimensions.y), o.matrix_world.translation)
                bpy.ops.object.select_all(action="DESELECT")
                o.select = True
                context.scene.objects.active = o
                bpy.ops.object.mode_set(mode="EDIT")
                context.window_manager.modal_handler_add(self)
                return {'RUNNING_MODAL'}
            else:
                self.report({'WARNING'}, "Profile not found")
                return {'CANCELLED'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}


class ARCHIPACK_OT_profile_save(Operator):
    bl_idname = "archipack.profile_save"
    bl_label = "Save"
    bl_description = "Save a profile curve"
    bl_options = {'REGISTER', 'UNDO', 'INTERNAL'}

    curve_name = StringProperty(
        description="Profile curve name",
        default=""
        )
    preset_name = StringProperty(
        description="Preset name",
        default=""
        )

    @classmethod
    def poll(self, context):
        return context.active_object is not None

    def invoke(self, context, event):
        # if context.space_data is not None and context.space_data.type == 'VIEW_3D':
        o = context.scene.objects.get(self.curve_name)
        if o and o.type == 'CURVE':
            manager = ArchipackProfileManager()
            manager.save_curve(o, self.preset_name)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'}


def register():
    bpy.utils.register_class(ARCHIPACK_OT_profile_edit)
    bpy.utils.register_class(ARCHIPACK_OT_object_update)
    bpy.utils.register_class(ARCHIPACK_OT_profile_save)


def unregister():
    bpy.utils.unregister_class(ARCHIPACK_OT_profile_edit)
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
