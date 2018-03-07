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
# ----------------------------------------------------------

import bpy
from mathutils import Vector, Matrix
from bpy_extras.io_utils import ExportHelper
from bpy.types import BezierSplinePoint, Operator
from .archipack_autoboolean import ArchipackBoolManager


XML_HEADER = """<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 20010904//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg id="document" width="{0}px" height="{1}px"  viewBox="0 0 {0} {1}" xmlns="http://www.w3.org/2000/svg">
    <title id="title">{2}</title>
    <desc id="desc">Archipack svg exporter by S.Leger</desc>
"""
XML_PATH = """<path id="{}" style="
    fill:{};
    stroke:{};
    stroke-width:{}px;
    stroke-linecap:butt;
    stroke-linejoin:round;
    stroke-opacity:1"
    d="{}"/>
"""
XML_TEXT = """<text id="{}" style="
    font-style:normal;
    font-weight:normal;
    font-size:{}px;
    line-height:125%;
    font-family:sans-serif;
    letter-spacing:0px;
    word-spacing:0px;
    fill:{};
    fill-opacity:1;
    stroke:none;
    stroke-width:1px;
    stroke-linecap:butt;
    stroke-linejoin:miter;
    stroke-opacity:1;
    text-anchor:middle;
    text-align:center;"
    x="0" y="0" transform="{}">{}</text>
"""
XML_GROUP = """<g id="{}">
"""
XML_GROUP_END = """</g>
"""

XML_END = "</svg>"
MOVE_COMMAND = 'M{},{} '
LINE_COMMAND = 'L{},{} '
CURVE_COMMAND = 'C{},{} {},{} {},{} '
JOIN_COMMAND = 'Z '


class SvgGroup:
    def __init__(self, index, components):
        self.index = index
        self.components = components

    def output(self, f):
        f.write(XML_GROUP.format(self.index))
        for c in self.components:
            c.output(f)
        f.write(XML_GROUP_END)


class SvgPath:
    """
      path might contains multiple curves
      create with shape 0 and
      call add_shape to add shapes > 0
    """
    def __init__(self, index, tM, spline, width, fill_color, stroke_color):
        self.index = index
        self.width = width
        self.stroke = stroke_color
        self.path = []
        self.add_shape(tM, spline, fill_color)

    def add_shape(self, tM, spline, fill_color):

        self.fill = 'none'

        if spline.type == 'BEZIER':
            first_curve_point = None
            previous = None
            for n, point in enumerate(spline.bezier_points):
                co = tM * point.co.to_3d()
                if n == 0:
                    first_curve_point = point
                    self.path.append(MOVE_COMMAND.format(co.x, co.y))
                else:
                    self.path.append(self.make_curve_command(previous, tM, point))
                previous = point

            if spline.use_cyclic_u:
                self.path.append(self.make_curve_command(previous, tM, first_curve_point))

        elif spline.type == 'POLY':
            for n, point in enumerate(spline.points):
                co = tM * point.co.to_3d()
                if n == 0:
                    self.path.append(MOVE_COMMAND.format(co.x, co.y))
                else:
                    self.path.append(LINE_COMMAND.format(co.x, co.y))

        if spline.use_cyclic_u:
            self.path.append(JOIN_COMMAND)
            self.fill = fill_color

    def make_curve_command(self, previous, tM, point):
        co = tM * point.co.to_3d()
        left = tM * point.handle_left.to_3d()
        right = tM * previous.handle_right.to_3d()
        return CURVE_COMMAND.format(right.x, right.y, left.x, left.y, co.x, co.y)

    def output(self, f):
        f.write(XML_PATH.format(self.index, self.fill, self.stroke, self.width, "".join(self.path)))


class SvgText:
    def __init__(self, tM, curve, scale, stroke_color):
        x, y, z = tM * Vector((0, 0, 0))
        rM = tM.normalized().to_3x3()
        rxx, rxy = rM[0][0:2]
        ryx, ryy = rM[1][0:2]
        self.index = curve.name
        self.size = curve.data.size / scale
        self.stroke = stroke_color
        self.matrix = "matrix({},{},{},{},{},{})".format(rxx, rxy, -ryx, -ryy, x, y)
        self.body = curve.data.body

    def output(self, f):
        f.write(XML_TEXT.format(self.index, self.size, self.stroke, self.matrix, self.body))


class SVGStyle:
    def __init__(self, width, stroke, fill):
        self.width = width
        self.stroke = stroke
        self.fill = fill


def remove_curve(context, curve):
    d = curve.data
    context.scene.objects.unlink(curve)
    bpy.data.objects.remove(curve, do_unlink=True)
    bpy.data.curves.remove(d)


def SVG_blenderCurve(itM, curve, style, components, as_group=True):
    """
      when as_group is true, add components and
      append svg_path for each shape of curve into components and return a SvgGroup
      when as_group is false, add all shapes into same SvgPath entity
      append svg_path to components and return svg_path
    """
    name = curve.name
    tM = itM * curve.matrix_world
    if as_group:
        for i, spline in enumerate(curve.data.splines):
            components.append(SvgPath("{}-{}".format(name, i), tM, spline, style.width, style.fill, style.stroke))
        return SvgGroup(name, components)
    else:
        for i, spline in enumerate(curve.data.splines):
            if i == 0:
                svg_path = SvgPath("{}-0".format(name), tM, spline, style.width, style.fill, style.stroke)
            else:
                svg_path.add_shape(tM, spline, style.fill)
        components.append(svg_path)
        return svg_path


def SVG_dimension(itM, dimension, scale, style):
    components = []
    for txt in dimension.children:
        components.append(SvgText(itM * txt.matrix_world, txt, scale, style.fill))
    # add the curve into components
    SVG_blenderCurve(itM, dimension, style, components, as_group=False)
    return SvgGroup("{}".format(dimension.name), components)


def SVG_wall_childs(context, itM, wall, scale, styles, openings, dimensions):

    wd = wall.data.archipack_wall2[0]
    # windows / doors
    for child in wd.childs:
        c, d = child.get_child(context)
        if d is None:
            if "archipack_wall2" in c.data:
                SVG_wall_childs(context, itM, c, scale, styles, openings, dimensions)
        else:
            symbol = d.as_2d(context, c)
            # add window in her own group, let separate curves so panels fill override frame
            openings.append(SVG_blenderCurve(itM, symbol, styles['openings'], [], as_group=True))
            remove_curve(context, symbol)
            # window / door dimensions
            for child in c.children:
                d = child.data
                if d:
                    if "archipack_dimension_auto" in d:
                        dimensions.append(SVG_dimension(itM, child, scale, styles['dim']))
    # wall dimensions
    for child in wall.children:
        d = child.data
        if d:
            if "archipack_dimension_auto" in d:
                dimensions.append(SVG_dimension(itM, child, scale, styles['dim']))


def SVG_wall(context, itM, wall, scale, styles):
    """
     Export wall with dimension and openings symbols
    """
    openings = []
    dimensions = []
    parts = []
    bpy.ops.object.select_all(action="DESELECT")
    wall.select = True
    context.scene.objects.active = wall
    bpy.ops.archipack.wall2_to_curve(mode='SYMBOL')
    w = context.active_object
    f = [o for o in context.selected_objects if o.name != w.name]

    # wall plain parts
    SVG_blenderCurve(itM, w, styles['wall'], parts, as_group=False)
    remove_curve(context, w)

    # wall fill under windows
    if len(f) > 0:
        f = f[0]
        SVG_blenderCurve(itM, f, styles['hole'], parts, as_group=False)
        remove_curve(context, f)

    # windows / doors
    SVG_wall_childs(context, itM, wall, scale, styles, openings, dimensions)

    parts.append(SvgGroup("{}-openings".format(wall.name), openings))
    parts.append(SvgGroup("{}-dimensions".format(wall.name), dimensions))

    return SvgGroup(wall.name, parts)


def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


def getZ(curve):
    return curve.matrix_world.translation.z


def point_str(itM, point):
    co = itM * point.co
    if isinstance(point, BezierSplinePoint):
        left = itM * point.handle_left
        right = itM * point.handle_right
        return "[x: {}, y: {}, z: {},\n" \
               " x: {}, y: {}, z: {},\n" \
               " x: {}, y: {}, z: {}]".format(left.x, left.y, left.z,
                                              co.x, co.y, co.z,
                                              right.x, right.y, right.z)
    else:
        return "[x: {}, y: {}, z: {}]".format(co.x, co.y, co.z)


class ARCHIPACK_OP_ExportSvg(Operator, ExportHelper):
    bl_idname = "archipack.export_svg"
    bl_label = "Archipack Inkscape SVG Exporter (2d curves only)"
    bl_options = {'PRESET'}
    filename_ext = ".svg"
    # ExportHelper class properties
    filter_glob = bpy.props.StringProperty(
            default="*.svg",
            options={'HIDDEN'},
            )
            
    def get_topmost_parent(self, o):
        if o.parent:
            return self.get_topmost_parent(o.parent)
        else:
            return o
    
    def get_dimensions(self, o, dims):
        if o.data and "archipack_dimension_auto" in o.data:
            dims.append(o)
        for c in o.children:
            self.get_dimensions(c, dims)
    
    def execute(self, context):
        bpy.ops.object.mode_set(mode='OBJECT')
        result = {'FINISHED'}
        scene = context.scene

        layout = context.active_object
        d = layout.data.archipack_layout[0]
        curves = [o for o in scene.objects if o.type in {'CURVE', 'FONT'}]

        if len(curves) < 1:
            self.report({'ERROR'}, "Nothing to export, 0 curve(s) selected")
            return {'CANCELLED'}

        # use manager to find curves in bounding box
        manager = ArchipackBoolManager()
        manager._init_bounding_box(layout)

        # exclude z filtering
        manager.minz = -1e32
        manager.maxz = 1e32
        manager.center.z = 0

        # retrieve archipack's entity
        # dimension are childs of wall
        walls = [
            o for o in scene.objects
            if o.data and
            "archipack_wall2" in o.data and
            manager._contains(o)
            ]
        stairs = [
            o for o in scene.objects
            if o.data and
            "archipack_stair" in o.data and
            manager._contains(o)
            ]

        # document size according layout (mm)
        width, height = d.paper_size
        # pixels / mm
        pixel_size = d.pixel_size
        w, h = d.canvas_size(context)

        # a matrix wich convert world coords
        # into paper coords pixels
        s = d.canvas_scale(context) / pixel_size
        tM = layout.matrix_world.copy()
        # tM[1][3] += Vector((0, h, 0))
        itM = (tM * Matrix([
            [s, 0, 0, 0],
            [0, -s, 0, h],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
            ])).inverted()

        # filter curves found in layout
        curves = {o.name: o for o in curves if manager._contains(o) and o.name != layout.name}
        
        # remove dimension
        parents = {}
        for wall in walls:
            p = self.get_topmost_parent(wall)
            if p.name not in parents:
                parents[p.name] = p
        
        dims = []
        for p in parents.values():
            self.get_dimensions(p, dims)
        
        for dim in dims:
            if dim.name in curves:
                del curves[dim.name]
                for child in dim.children:
                    if child.name in curves:
                        del curves[child.name]
                    
        # line width in mm
        line_width = 0.2 * pixel_size
        styles = {
            #               line_width        stroke    fill
            'wall': SVGStyle(line_width, "#000000", "#a0a0a0"),
            'hole': SVGStyle(line_width, "#000000", "#c0c0c0"),
            'openings': SVGStyle(line_width, "#000000", "#FFFFFF"),
            'dim': SVGStyle(0.5 * line_width, "#000000", "#000000"),
            'curves': SVGStyle(line_width, "#000000", "#000000"),
            }

        # walls with symbols and dimensions
        svg_walls = [
            SVG_wall(
                context, itM, w, s, styles
                ) for w in walls if w.data.archipack_wall2[0].t_part == ""]

        # lines / text not part of walls
        svg_lines = []
        for curve in curves.values():
            tM = itM * curve.matrix_world
            # @TODO:
            # define color override as per curve object
            # with global default
            style = styles['curves']
            if len(curve.data.materials) > 0:
                style.stroke = rgb_to_hex((
                    int(curve.data.materials[0].diffuse_color.r * 255),
                    int(curve.data.materials[0].diffuse_color.g * 255),
                    int(curve.data.materials[0].diffuse_color.b * 255)))
                style.fill = rgb_to_hex((
                    int(curve.data.materials[0].specular_color.r * 255),
                    int(curve.data.materials[0].specular_color.g * 255),
                    int(curve.data.materials[0].specular_color.b * 255)))
            else:
                style.stroke = '#000000'
                style.fill = '#808080'

            if curve.type == 'CURVE':
                svg_lines.append(SVG_blenderCurve(itM, curve, style, []))

            elif curve.type == 'FONT':
                # Export text
                svg_lines.append(SvgText(tM, curve, s, style.stroke))

        for stair in stairs:
            stair.select = True
            context.scene.objects.active = stair
            bpy.ops.archipack.stair_to_curve(mode="SYMBOL")
            curve = context.active_object
            svg_lines.append(SVG_blenderCurve(itM, curve, styles['openings'], [], as_group=False))
            remove_curve(context, curve)
            stair.select = False

        # Open the file for writing
        with open(self.filepath, 'w') as f:
            f.write(XML_HEADER.format(width * pixel_size, height * pixel_size, layout.name))
            for line in svg_lines:
                line.output(f)
            for wall in svg_walls:
                wall.output(f)
            f.write(XML_END)
        layout.select = True
        context.scene.objects.active = layout
        return result


# Register in File > Import menu
def menu_func_export(self, context):
    self.layout.operator('archipack.export_svg')


def register():
    bpy.utils.register_class(ARCHIPACK_OP_ExportSvg)
    # bpy.types.INFO_MT_file_export.append(menu_func_export)


def unregister():
    import bpy
    # bpy.types.INFO_MT_file_export.remove(menu_func_export)
    bpy.utils.unregister_class(ARCHIPACK_OP_ExportSvg)
