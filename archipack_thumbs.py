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
# Inspired by Asset-Flinguer
# ----------------------------------------------------------
import sys
from mathutils import Vector
from mathutils.geometry import interpolate_bezier
import bpy


class ArchipackProfileLoader():
    """
     IO for user defined profiles
    """
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

        return pts

    def _from_json(self, curve, json_str):
        """
          Adds a spline to a curve
        """
        import json
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

    def curve_center(self, curve):
        # estimate curve size
        pts = []
        for spline in curve.data.splines:
            pts.extend(self.coords_from_spline(spline, curve.matrix_world, 12))
        x = [co.x for co in pts]
        y = [co.y for co in pts]
        cx = 0.5 * (max(x) + min(x))
        cy = 0.5 * (max(y) + min(y))
        return cx, cy

    def load_curve(self, preset_name):
        """
          Load a json spline definition
          return curve data not linked to scene
        """
        curve = None
        x, y = 0, 0
        # load from current directory
        full_path = "{}.json".format(preset_name[:-4])
        print("load_curve", full_path)
        try:
            with open(full_path) as fh:
                json_str = ''.join(fh.readlines())
            curve = bpy.data.curves.new(name=bpy.path.display_name_from_filepath(full_path), type='CURVE')
            curve.dimensions = '2D'
            curve.fill_mode = 'BOTH'
            x, y = self._from_json(curve, json_str)
        except:
            print("Archipack: Error while loading {}".format(preset_name))
            pass
        return curve, x, y


def log(s):
    print("[log]" + s)


def create_lamp(context, loc):
    bpy.ops.object.lamp_add(
        type='POINT',
        radius=1,
        view_align=False,
        location=loc)
    lamp = context.active_object
    lamp.data.use_nodes = True
    tree = lamp.data.node_tree
    return tree, tree.nodes


def create_camera(context, loc, rot):
    bpy.ops.object.camera_add(
        view_align=True,
        enter_editmode=False,
        location=loc,
        rotation=rot)
    cam = context.active_object
    context.scene.camera = cam
    return cam


def get_center(o):
    x, y, z = o.bound_box[0]
    min_x = x
    min_y = y
    min_z = z
    x, y, z = o.bound_box[6]
    max_x = x
    max_y = y
    max_z = z
    return Vector((
        min_x + 0.5 * (max_x - min_x),
        min_y + 0.5 * (max_y - min_y),
        min_z + 0.5 * (max_z - min_z)))


def apply_simple_material(o, name, color):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    m.node_tree.nodes[1].inputs[0].default_value = color
    o.data.materials.append(m)


def generateThumb(context, cls, preset):
    log("### RENDER THUMB ############################")

    # Cleanup scene
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()

    log("Start generating: " + cls)

    # setup render
    context.scene.render.engine = 'CYCLES'
    render = context.scene.cycles
    render.progressive = 'PATH'
    render.samples = 24
    try:
        render.use_square_samples = True
    except:
        pass
    render.preview_samples = 24
    render.aa_samples = 24
    render.transparent_max_bounces = 8
    render.transparent_min_bounces = 8
    render.transmission_bounces = 8
    render.max_bounces = 8
    render.min_bounces = 6
    render.caustics_refractive = False
    render.caustics_reflective = False
    render.use_transparent_shadows = True
    render.diffuse_bounces = 1
    render.glossy_bounces = 4
    render = context.scene.render

    # curve thumbs
    if cls == 'archipack_curves':
        render.resolution_x = 128
        render.resolution_y = 128
        render.filepath = preset
        render.alpha_mode = 'TRANSPARENT'

        context.scene.world.horizon_color = (0.7, 0.7, 0.7)

        man = ArchipackProfileLoader()
        # preset is image name
        d, x, y = man.load_curve(preset)
        curve = bpy.data.objects.new("Curve", d)
        try:
            # 2.7
            context.scene.objects.link(curve)
        except:
            pass
        try:
            # 2.8
            bpy.context.scene_collection.objects.link(curve)
        except:
            pass
            
        # bevel open splines
        for spline in curve.data.splines:
            if not spline.use_cyclic_u:
                d.bevel_resolution = 2
                d.bevel_depth = 0.0125 * max(x, y)
                d.fill_mode = 'NONE'

        apply_simple_material(curve, "Curve", (0, 0, 0, 1))

        # setup ortho camera on top view
        log("Prepare camera")
        cx, cy = man.curve_center(curve)
        cam = create_camera(context, (cx, cy, 1), (0, 0, 0))
        cam.data.type = 'ORTHO'
        cam.data.ortho_scale = max(x, y)
        cam.data.clip_start = 0.001
        cam.data.clip_end = 100

    else:
        # engine settings
        render.resolution_x = 150
        render.resolution_y = 100
        render.filepath = preset[:-3] + ".png"

        # create object, loading preset
        getattr(bpy.ops.archipack, cls)('INVOKE_DEFAULT', filepath=preset, auto_manipulate=False)

        o = context.active_object
        size = o.dimensions
        center = get_center(o)

        # oposite / tan (0.5 * fov)  where fov is 49.134 deg
        dist = max(size) / 0.32
        loc = center + dist * Vector((-0.5, -1, 0.5)).normalized()

        log("Prepare camera")
        cam = create_camera(context, loc, (1.150952, 0.0, -0.462509))
        cam.data.lens = 50

        bpy.ops.object.select_all(action="DESELECT")
        try:
            o.select = True
        except:
            pass
        try:
            o.select_set(action="SELECT")
        except:
            pass
        bpy.ops.view3d.camera_to_view_selected()

        log("Prepare scene")
        # add plane
        bpy.ops.mesh.primitive_plane_add(
            radius=1000,
            view_align=False,
            enter_editmode=False,
            location=(0, 0, 0)
            )
        p = context.active_object
        apply_simple_material(p, "Plane", (1, 1, 1, 1))

        # add 3 lights
        tree, nodes = create_lamp(context, (3.69736, -7, 6.0))
        emit = nodes["Emission"]
        emit.inputs[1].default_value = 2000.0

        tree, nodes = create_lamp(context, (9.414563179016113, 5.446230888366699, 5.903861999511719))
        emit = nodes["Emission"]
        falloff = nodes.new(type="ShaderNodeLightFalloff")
        falloff.inputs[0].default_value = 5
        tree.links.new(falloff.outputs[2], emit.inputs[1])

        tree, nodes = create_lamp(context, (-7.847615718841553, 1.03135085105896, 5.903861999511719))
        emit = nodes["Emission"]
        falloff = nodes.new(type="ShaderNodeLightFalloff")
        falloff.inputs[0].default_value = 5
        tree.links.new(falloff.outputs[2], emit.inputs[1])

    # Set output filename.
    render.use_file_extension = True
    render.use_overwrite = True
    render.use_compositing = False
    render.use_sequencer = False
    render.resolution_percentage = 100
    # render.image_settings.file_format = 'PNG'
    # render.image_settings.color_mode = 'RGBA'
    # render.image_settings.color_depth = '8'

    # Configure output size.
    log("Render")

    # Render thumbnail
    bpy.ops.render.render(write_still=True)

    log("### COMPLETED ############################")


if __name__ == "__main__":
    preset = ""

    for arg in sys.argv:
        if arg.startswith("cls:"):
            cls = arg[4:]
        if arg.startswith("preset:"):
            preset = arg[7:]
        if arg.startswith("matlib:"):
            matlib = arg[7:]
        if arg.startswith("addon:"):
            module = arg[6:]
    try:
        bpy.ops.wm.addon_enable(module=module)
        bpy.context.user_preferences.addons[module].preferences.matlib_path = matlib
    except:
        raise RuntimeError("module name not found")
    generateThumb(bpy.context, cls, preset)
