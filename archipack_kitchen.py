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
# Inspired by archimesh
# Automatic generation of kitchen cabinet
# Author: Antonio Vazquez (antonioya)
#
# ----------------------------------------------------------
# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty, IntProperty, BoolProperty,
    CollectionProperty, EnumProperty
)
from mathutils import Vector, Matrix
from math import pi, sin, cos
from .bmesh_utils import BmeshEdit as bmed
from .panel import Panel as Lofter
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_preset import ArchipackPreset, PresetMenuOperator
# from .archipack_gl import FeedbackPanel
from .archipack_object import ArchipackObject, ArchipackCreateTool
from .archipack_dimension import DimensionProvider
# from .archipack_keymaps import Keymaps
tan22_5 = (2 ** 0.5 - 1)
sq2 = 0.5 * 2 ** 0.5


# ----------------------------------------------------------
#  Rotation types
# ----------------------------------------------------------
RotationType_Default = 0
RotationType_R90CW = 1
RotationType_R90CCW = 2
RotationType_R180 = 3
RotationType_User = 4

# ----------------------------------------------------------
# Material index
# ----------------------------------------------------------
mat_inside = 0
mat_base = 1
mat_wall = 2
mat_full = 3
mat_border = 4
mat_counter = 5
mat_baseboard = 6
mat_handle = 7
mat_glass = 8
mat_metal = 9
mat_owen_handle = 10
mat_owen_inside = 11
mat_owen_metal = 9
mat_owen_glass = 8
mat_owen_panel = 12
mat_dishwasher_panel = 9
mat_dishwasher_face = 1
mat_dishwasher_handle = 12
mat_cooker_side = 9
mat_cooker_top = 12
mat_sink = 9
mat_range_body = 9
mat_range_filter = 12


def update(self, context):
    self.update(context)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


def update_parent(self, context):
    self.update_parent(context)


def chanfer_square(x, y, z, r, sx, sy):
    """
        x, y, z of center
        sx, sy size
        r radius
        back
           x0 x1
        y0 2   3
        y1 1   4
        front
    """

    verts = []
    dx = sx - 2 * r
    dy = sy - 2 * r
    seg = 12
    q = int(1 + seg / 4)
    da = 2 * pi / seg
    xc = x + r - 0.5 * sx
    yc = y + r - 0.5 * sy
    a = 9 * da
    for j in range(2):
        for k in range(2):
            for n in range(q):
                verts.append((round(xc + r * cos(a), 3), round(yc + r * sin(a), 3), z))
                a -= da
            yc += dy
            a += da
        yc = y - r + 0.5 * sy
        dy = -dy
        xc += dx
    return verts


def make_box(tM, x0, x1, y0, y1, z0, z1, r, idmat, chanfer, verts, faces, matids, uvs):
    f = len(verts)
    if chanfer == 3 and r != 0:
        # Chanfer right front
        s = 5
        dx, dy = x1 - x0, y1 - y0
        r = min(abs(r), abs(dx), abs(dy))
        ry = r
        if dy < 0:
            ry = -ry
        rx = r
        if dx < 0:
            rx = -rx
        y2 = y1 - ry
        x2 = x0 + rx
        x3 = x1 - rx
        verts.extend([tM * Vector((v)) for v in [
            (x1, y0, z0), (x1, y1, z0), (x2, y1, z0), (x0, y2, z0), (x0, y0, z0),
            (x1, y0, z1), (x1, y1, z1), (x2, y1, z1), (x0, y2, z1), (x0, y0, z1)
            ]])
        uvs.extend([[(0, 0), (1, 0), (1, 0.9), (0.9, 1), (0, 1)] for i in range(2)])

    elif chanfer == 2 and r != 0:
        s = 6
        # chanfer front side
        dx, dy = x1 - x0, y1 - y0
        r = min(abs(r), 0.3 * abs(dx), 0.3 * abs(dy))
        ry = r
        if dy < 0:
            ry = -ry
        rx = r
        if dx < 0:
            rx = -rx
        x2 = x0 + rx
        y2 = y1 - ry
        y3 = y0 + ry
        verts.extend([tM * Vector((v)) for v in [
            (x1, y0, z0), (x1, y1, z0), (x2, y1, z0), (x0, y2, z0), (x0, y3, z0), (x2, y0, z0),
            (x1, y0, z1), (x1, y1, z1), (x2, y1, z1), (x0, y2, z1), (x0, y3, z1), (x2, y0, z1)
            ]])
        uvs.extend([[(0, 0), (1, 0), (1, 0.9), (0.9, 1), (0.1, 1), (0, 0.9)] for i in range(2)])

    elif chanfer == 1 and r != 0:
        # chanfer right side
        s = 6
        dx, dy = x1 - x0, y1 - y0
        r = min(abs(r), 0.3 * abs(dx), 0.3 * abs(dy))
        ry = r
        if dy < 0:
            ry = -ry
        rx = r
        if dx < 0:
            rx = -rx
        y2 = y1 - ry
        x2 = x0 + rx
        x3 = x1 - rx
        verts.extend([tM * Vector((v)) for v in [
            (x1, y0, z0), (x1, y2, z0), (x3, y1, z0), (x2, y1, z0), (x0, y2, z0), (x0, y0, z0),
            (x1, y0, z1), (x1, y2, z1), (x3, y1, z1), (x2, y1, z1), (x0, y2, z1), (x0, y0, z1)
            ]])
        uvs.extend([[(0, 0), (1, 0), (1, 0.9), (0.9, 1), (0.1, 1), (0, 0.9)] for i in range(2)])

    else:
        s = 4
        verts.extend([tM * Vector((v)) for v in [
            (x0, y1, z0), (x0, y0, z0), (x1, y0, z0), (x1, y1, z0),
            (x0, y1, z1), (x0, y0, z1), (x1, y0, z1), (x1, y1, z1)
            ]])
        uvs.extend([[(0, 0), (1, 0), (1, 1), (0, 1)] for i in range(2)])
    # sides
    faces.extend([(f + i + 1, f + i, f + s + i, f + s + i + 1) for i in range(s - 1)])
    # back
    faces.append((f, f + s - 1, f + 2 * s - 1, f + s))
    # bottom
    faces.append(tuple([f + i for i in range(s)]))
    # top
    faces.append(tuple([f + 2 * s - 1 - i for i in range(s)]))
    matids.extend([idmat for i in range(s + 2)])
    uvs.extend([[(0, 0), (1, 0), (1, 1), (0, 1)] for i in range(s)])


def dishwasher_door(size, verts, faces, matids, uvs):
    f = len(verts)
    cm = size.x / 60
    # x from left = 0 to right
    x0 = 0
    x1 = x0 + 0.5 * size.x - 5 * cm
    x2 = x0 + 0.5 * size.x + 5 * cm
    x3 = size.x

    # y from inside = th to outside
    y0 = size.z - 0.1 * cm
    y1 = y0 - 0.1 * cm
    y2 = 0

    # z from bottom = 0 to top
    z0 = 0.5 * cm
    z4 = size.y
    z5 = z4 - 2 * cm
    z3 = z5 - 3.35 * cm
    z2 = z3 - 4 * cm
    z1 = z2 - 3.35 * cm
    verts.extend([Vector(v) for v in [
        (x0, y0, z5), (x0, y0, z1), (x3, y0, z1),
        (x3, y0, z5), (x2, y2, z3), (x2, y2, z2),
        (x1, y2, z2), (x1, y2, z3), (x3, y0, z4),
        (x0, y0, z4), (x0, y2, z4), (x3, y2, z4),
        (x0, y0, z0), (x0, y2, z0), (x3, y0, z0),
        (x3, y2, z0), (x3, y2, z1), (x0, y2, z1),
        (x3, y2, z5), (x0, y2, z5), (x2, y1, z2),
        (x1, y1, z2), (x2, y1, z3), (x1, y1, z3)
        ]])
    faces.extend([tuple([f + i for i in v]) for v in [
        (19, 0, 1, 17), (16, 2, 3, 18), (11, 10, 19, 18),
        (9, 0, 19, 10), (17, 1, 12, 13), (3, 8, 11, 18),
        (10, 11, 8, 9), (0, 9, 8, 3, 2, 1), (4, 7, 23, 22),
        (15, 13, 12, 14), (1, 2, 14, 12), (2, 16, 15, 14),
        (16, 17, 13, 15), (5, 6, 17, 16), (4, 5, 16, 18),
        (7, 4, 18, 19), (6, 7, 19, 17), (20, 22, 23, 21),
        (6, 5, 20, 21), (7, 6, 21, 23), (5, 4, 22, 20)
        ]])
    matids.extend([
        mat_dishwasher_panel, mat_dishwasher_panel, mat_dishwasher_face,
        mat_dishwasher_face, mat_dishwasher_face, mat_dishwasher_face,
        mat_dishwasher_face, mat_dishwasher_panel, mat_dishwasher_panel,
        mat_dishwasher_handle, mat_dishwasher_panel, mat_dishwasher_face,
        mat_dishwasher_face, mat_dishwasher_panel, mat_dishwasher_panel,
        mat_dishwasher_panel, mat_dishwasher_panel, mat_dishwasher_panel,
        mat_dishwasher_panel, mat_dishwasher_panel, mat_dishwasher_panel
        ])
    uvs.extend([
        [(0.0, 0.97), (0.041, 0.97), (0.041, 0.811), (0.0, 0.811)],
        [(0.0, 0.811), (0.041, 0.811), (0.041, 0.97), (0.0, 0.97)],
        [(1.0, 1.0), (0.042, 1.0), (0.042, 0.97), (1.0, 0.97)],
        [(0.041, 1.0), (0.041, 0.97), (0.0, 0.97), (0.0, 1.0)],
        [(0.0, 0.811), (0.041, 0.811), (0.041, 0.811), (0.0, 0.811)],
        [(0.041, 0.97), (0.041, 1.0), (0.0, 1.0), (0.0, 0.97)],
        [(0.042, 0.0), (1.0, 0.0), (1.0, 0.039), (0.042, 0.039)],
        [(0.042, 0.97), (0.042, 1.0), (1.0, 1.0), (1.0, 0.97), (1.0, 0.811), (0.042, 0.811)],
        [(0.942, 0.914), (0.101, 0.914), (0.101, 0.914), (0.942, 0.914)],
        [(1.0, 0.0), (0.042, 0.0), (0.042, 0.039), (1.0, 0.039)],
        [(0.042, 0.811), (1.0, 0.811), (1.0, 0.811), (0.042, 0.811)],
        [(0.041, 0.811), (0.0, 0.811), (0.0, 0.811), (0.041, 0.811)],
        [(1.0, 0.811), (0.042, 0.811), (0.042, 0.811), (1.0, 0.811)],
        [(0.942, 0.866), (0.101, 0.866), (0.042, 0.811), (1.0, 0.811)],
        [(0.942, 0.914), (0.942, 0.866), (1.0, 0.811), (1.0, 0.97)],
        [(0.101, 0.914), (0.942, 0.914), (1.0, 0.97), (0.042, 0.97)],
        [(0.101, 0.866), (0.101, 0.914), (0.042, 0.97), (0.042, 0.811)],
        [(0.942, 0.866), (0.942, 0.914), (0.101, 0.914), (0.101, 0.866)],
        [(0.101, 0.866), (0.942, 0.866), (0.942, 0.866), (0.101, 0.866)],
        [(0.101, 0.914), (0.101, 0.866), (0.101, 0.866), (0.101, 0.914)],
        [(0.942, 0.866), (0.942, 0.914), (0.942, 0.914), (0.942, 0.866)]
        ])


def dishwasher_cab(tM, x, y, z, z0, th, dy, verts, faces, matids, uvs):
    f = len(verts)
    cm = x / 60
    # x from left to right
    x0 = th
    x1 = x0 + 4 * cm
    x3 = x - th
    x2 = x3 - 4 * cm
    # y from back to front
    y0 = -th
    y1 = y0 - 4 * cm
    y2 = -y
    # z from bottom to top
    z0 += th
    z3 = z0 + z - 2 * th
    z1 = z0 + 3 * cm
    z2 = z3 - 3 * cm

    verts.extend([tM * Vector(v) for v in [
        (x3, y2, z3), (x3, y2, z0), (x0, y2, z0),
        (x0, y2, z3), (x0, y0, z3), (x0, y0, z0),
        (x3, y0, z0), (x3, y0, z3), (x1, y2, z2),
        (x1, y2, z1), (x2, y2, z1), (x2, y2, z2),
        (x2, y1, z2), (x2, y1, z1), (x1, y1, z1),
        (x1, y1, z2)
        ]])
    faces.extend([tuple([f + i for i in v]) for v in [
        (2, 3, 4, 5), (3, 0, 7, 4), (0, 1, 6, 7),
        (1, 2, 5, 6), (5, 4, 7, 6), (3, 2, 9, 8),
        (0, 3, 8, 11), (1, 0, 11, 10), (2, 1, 10, 9),
        (9, 10, 13, 14), (10, 11, 12, 13), (11, 8, 15, 12),
        (8, 9, 14, 15), (12, 15, 14, 13)
        ]])
    matids.extend([mat_metal for i in range(14)])
    uvs.extend([
        [(0.042, 0.04), (0.042, 1.0), (0.938, 1.0), (0.938, 0.04)],
        [(0.074, 0.04), (0.968, 0.04), (0.968, 0.889), (0.074, 0.889)],
        [(0.042, 1.0), (0.042, 0.04), (0.938, 0.04), (0.938, 1.0)],
        [(0.968, 0.04), (0.074, 0.04), (0.074, 0.889), (0.968, 0.889)],
        [(0.074, 0.04), (0.074, 1.0), (0.968, 1.0), (0.968, 0.04)],
        [(0.074, 1.0), (0.074, 0.04), (0.106, 0.101), (0.106, 0.758)],
        [(0.968, 1.0), (0.074, 1.0), (0.106, 0.758), (0.936, 0.758)],
        [(0.968, 0.04), (0.968, 1.0), (0.936, 0.758), (0.936, 0.101)],
        [(0.074, 0.04), (0.968, 0.04), (0.936, 0.101), (0.106, 0.101)],
        [(0.106, 0.04), (0.936, 0.04), (0.936, 0.796), (0.106, 0.796)],
        [(0.042, 0.101), (0.042, 0.758), (0.84, 0.758), (0.84, 0.101)],
        [(0.936, 0.04), (0.106, 0.04), (0.106, 0.796), (0.936, 0.796)],
        [(0.042, 0.758), (0.042, 0.101), (0.84, 0.101), (0.84, 0.758)],
        [(0.936, 0.758), (0.106, 0.758), (0.106, 0.101), (0.936, 0.101)]
        ])


def oven_door(size, verts, faces, matids, uvs):
    cm = size.x / 60
    # x from left = 0 to right
    x0 = 0
    x5 = size.x
    x1 = x0 + 6.5 * cm
    x2 = x1 + 3 * cm
    x4 = x5 - 6.5 * cm
    x3 = x4 - 3 * cm
    # y from inside = th to outside
    y0 = size.z - 0.1 * cm
    y1 = 0
    y2 = y1 - 2 * cm
    y3 = y2 - 1 * cm
    # z from bottom = 0 to top
    z0 = 0.5 * cm
    z1 = z0 + 6.5 * cm
    z5 = size.y - 12.7 * cm
    z4 = z5 - 2.75 * cm
    z3 = z4 - 1 * cm
    z2 = z5 - 6.5 * cm
    f = len(verts)
    verts.extend([Vector(v) for v in [
        (x5, y1, z5), (x5, y1, z0), (x0, y1, z0),
        (x0, y1, z5), (x0, y0, z5), (x0, y0, z0),
        (x5, y0, z0), (x5, y0, z5), (x0, y1, z2),
        (x0, y1, z1), (x5, y1, z1), (x5, y1, z2),
        (x5, y0, z2), (x5, y0, z1), (x0, y0, z1),
        (x0, y0, z2), (x2, y2, z3), (x3, y2, z3),
        (x2, y2, z4), (x3, y2, z4), (x1, y2, z3),
        (x1, y2, z4), (x1, y3, z4), (x1, y3, z3),
        (x4, y2, z4), (x4, y2, z3), (x4, y3, z3),
        (x4, y3, z4), (x3, y1, z4), (x3, y1, z3),
        (x2, y1, z3), (x2, y1, z4), (x1, y1, z4),
        (x1, y1, z3), (x4, y1, z3), (x4, y1, z4),
        (x1, y1, z2), (x4, y1, z2), (x1, y1, z1),
        (x4, y1, z1), (x4, y0, z2), (x1, y0, z2),
        (x4, y0, z1), (x1, y0, z1)
        ]])
    faces.extend([tuple([f + i for i in v]) for v in [
        (3, 0, 7, 4), (1, 2, 5, 6), (7, 12, 40, 41, 15, 4),
        (10, 1, 6, 13), (0, 3, 8, 36, 37, 11), (0, 11, 12, 7),
        (2, 1, 10, 39, 38, 9), (39, 10, 13, 42), (10, 13, 12, 11),
        (36, 8, 15, 41), (8, 15, 14, 9), (8, 3, 4, 15),
        (37, 39, 10, 11), (2, 9, 14, 5), (5, 14, 43, 42, 13, 6),
        (40, 12, 13, 42), (17, 16, 18, 19), (27, 22, 23, 26),
        (27, 24, 19, 18, 21, 22), (20, 23, 22, 21), (24, 19, 28, 35),
        (24, 27, 26, 25), (18, 21, 32, 31), (16, 17, 25, 26, 23, 20),
        (16, 30, 31, 18), (21, 32, 33, 20), (20, 16, 30, 33),
        (17, 25, 34, 29), (19, 28, 29, 17), (25, 34, 35, 24),
        (15, 41, 43, 14), (41, 40, 42, 43), (8, 9, 38, 36),
        (36, 38, 39, 37), (11, 37, 40, 12), (37, 36, 41, 40),
        (9, 38, 43, 14), (38, 39, 42, 43)
        ]])
    matids.extend([
        mat_owen_metal, mat_owen_metal, mat_owen_metal,
        mat_owen_metal, mat_owen_metal, mat_owen_metal,
        mat_owen_metal, mat_owen_metal, mat_owen_glass,
        mat_owen_metal, mat_owen_glass, mat_owen_metal,
        mat_owen_glass, mat_owen_metal, mat_owen_metal,
        mat_owen_panel, mat_owen_metal, mat_owen_metal,
        mat_owen_metal, mat_owen_metal, mat_owen_handle,
        mat_owen_metal, mat_owen_handle, mat_owen_metal,
        mat_owen_handle, mat_owen_handle, mat_owen_handle,
        mat_owen_handle, mat_owen_handle, mat_owen_handle,
        mat_owen_panel, mat_owen_glass, mat_owen_glass,
        mat_owen_glass, mat_owen_metal, mat_owen_metal,
        mat_owen_metal, mat_owen_metal
        ])
    uvs.extend([
        [(0.054, 0.62), (0.451, 0.604), (0.453, 0.623), (0.056, 0.639)],
        [(0.451, 0.477), (0.056, 0.494), (0.054, 0.475), (0.451, 0.458)],
        [(0.453, 0.623), (0.507, 0.628), (0.452, 0.693), (0.067, 0.693), (0.001, 0.65), (0.056, 0.639)],
        [(0.509, 0.48), (0.451, 0.477), (0.451, 0.458), (0.502, 0.444)],
        [(0.451, 0.604), (0.054, 0.62), (0.0, 0.614), (0.055, 0.549), (0.44, 0.549), (0.506, 0.592)],
        [(0.451, 0.604), (0.506, 0.592), (0.507, 0.628), (0.453, 0.623)],
        [(0.056, 0.494), (0.451, 0.477), (0.509, 0.48), (0.453, 0.549), (0.067, 0.548), (0.001, 0.505)],
        [(0.444, 0.741), (0.498, 0.741), (0.498, 0.763), (0.444, 0.763)],
        [(0.604, 0.246), (0.609, 0.268), (1.0, 0.168), (0.863, 0.28)],
        [(0.054, 0.784), (0.0, 0.784), (0.0, 0.763), (0.054, 0.763)],
        [(0.243, 0.274), (0.243, 0.124), (0.0, 0.246), (0.073, 0.261)],
        [(0.0, 0.614), (0.054, 0.62), (0.056, 0.639), (0.001, 0.65)],
        [(0.702, 0.386), (0.593, 0.227), (0.604, 0.246), (0.863, 0.28)],
        [(0.056, 0.494), (0.001, 0.505), (0.0, 0.469), (0.054, 0.475)],
        [(0.054, 0.475), (0.0, 0.469), (0.055, 0.404), (0.44, 0.404), (0.502, 0.444), (0.451, 0.458)],
        [(0.793, 0.06), (1.0, 0.168), (0.609, 0.268), (0.603, 0.274)],
        [(0.753, 0.694), (0.045, 0.729), (0.045, 0.741), (0.754, 0.705)],
        [(0.8, 0.698), (0.0, 0.739), (0.002, 0.734), (0.798, 0.695)],
        [(0.8, 0.698), (0.79, 0.699), (0.754, 0.705), (0.045, 0.741), (0.011, 0.74), (0.0, 0.739)],
        [(0.012, 0.733), (0.002, 0.734), (0.0, 0.739), (0.011, 0.74)],
        [(0.79, 0.699), (0.754, 0.705), (0.757, 0.696), (0.786, 0.694)],
        [(0.79, 0.699), (0.8, 0.698), (0.798, 0.695), (0.789, 0.696)],
        [(0.045, 0.741), (0.011, 0.74), (0.019, 0.737), (0.037, 0.737)],
        [(0.045, 0.729), (0.753, 0.694), (0.789, 0.696), (0.798, 0.695), (0.002, 0.734), (0.012, 0.733)],
        [(0.045, 0.729), (0.037, 0.734), (0.037, 0.737), (0.045, 0.741)],
        [(0.011, 0.74), (0.019, 0.737), (0.019, 0.735), (0.012, 0.733)],
        [(0.012, 0.733), (0.045, 0.729), (0.037, 0.734), (0.019, 0.735)],
        [(0.753, 0.694), (0.789, 0.696), (0.786, 0.702), (0.757, 0.702)],
        [(0.754, 0.705), (0.757, 0.696), (0.757, 0.702), (0.753, 0.694)],
        [(0.789, 0.696), (0.786, 0.702), (0.786, 0.694), (0.79, 0.699)],
        [(0.243, 0.124), (0.295, 0.0), (0.172, 0.232), (0.0, 0.246)],
        [(0.295, 0.0), (0.793, 0.06), (0.603, 0.274), (0.172, 0.232)],
        [(0.243, 0.274), (0.073, 0.261), (0.191, 0.235), (0.283, 0.404)],
        [(0.283, 0.404), (0.191, 0.235), (0.593, 0.227), (0.702, 0.386)],
        [(0.498, 0.784), (0.444, 0.784), (0.444, 0.763), (0.498, 0.763)],
        [(0.444, 0.784), (0.054, 0.784), (0.054, 0.763), (0.444, 0.763)],
        [(0.0, 0.741), (0.054, 0.741), (0.054, 0.763), (0.0, 0.763)],
        [(0.054, 0.741), (0.444, 0.741), (0.444, 0.763), (0.054, 0.763)]
        ])


def oven_cab(tM, x, y, z, z0, th, dy, verts, faces, matids, uvs):
    f = len(verts)
    cm = x / 60
    # x from left to right
    x0 = 0
    x5 = x
    x1 = x0 + th
    x2 = x1 + 4 * cm
    x4 = x5 - th
    x3 = x4 - 4 * cm
    # y from back to front
    y0 = -th
    y1 = y0 - 4 * cm
    y2 = -y
    y3 = y2 - 0.003 * cm
    y4 = y2 - dy
    # z from bottom to top

    z5 = z0 + z
    z6 = z5 - th
    z0 += th
    z1 = z0 + 4 * cm
    z3 = z5 - 12.5 * cm
    z2 = z3 - 4 * cm - th
    z4 = z5 - 2 * cm
    verts.extend([tM * Vector(v) for v in [
        (x4, y2, z6), (x4, y2, z0), (x1, y2, z0),
        (x1, y2, z6), (x1, y0, z6), (x1, y0, z0),
        (x4, y0, z0), (x4, y0, z6), (x2, y2, z2),
        (x2, y2, z1), (x3, y2, z1), (x3, y2, z2),
        (x3, y1, z2), (x3, y1, z1), (x2, y1, z1),
        (x2, y1, z2), (x0, y3, z4), (x0, y3, z3),
        (x5, y3, z3), (x5, y3, z4), (x5, y4, z4),
        (x5, y4, z3), (x0, y4, z3), (x0, y4, z4),
        (x5, y3, z5), (x0, y3, z5), (x0, y4, z5),
        (x5, y4, z5)
        ]])
    faces.extend([tuple([f + i for i in v]) for v in [
        (2, 3, 4, 5), (3, 0, 7, 4), (0, 1, 6, 7),
        (1, 2, 5, 6), (5, 4, 7, 6), (3, 2, 9, 8),
        (0, 3, 8, 11), (1, 0, 11, 10), (2, 1, 10, 9),
        (9, 10, 13, 14), (10, 11, 12, 13), (11, 8, 15, 12),
        (8, 9, 14, 15), (12, 15, 14, 13), (23, 16, 17, 22),
        (21, 18, 19, 20), (27, 26, 23, 20), (25, 16, 23, 26),
        (21, 22, 17, 18), (19, 24, 27, 20), (26, 27, 24, 25),
        (16, 25, 24, 19, 18, 17), (21, 20, 23, 22)
        ]])
    matids.extend([
        mat_owen_metal, mat_owen_metal, mat_owen_metal,
        mat_owen_metal, mat_owen_metal, mat_owen_inside,
        mat_owen_inside, mat_owen_inside, mat_owen_inside,
        mat_owen_inside, mat_owen_inside, mat_owen_inside,
        mat_owen_inside, mat_owen_inside, mat_owen_glass,
        mat_owen_glass, mat_owen_metal, mat_owen_metal,
        mat_owen_glass, mat_owen_metal, mat_owen_metal,
        mat_owen_panel, mat_owen_glass
        ])
    uvs.extend([
        [(0.042, 0.04), (0.042, 1.0), (0.938, 1.0), (0.938, 0.04)],
        [(0.074, 0.04), (0.968, 0.04), (0.968, 0.889), (0.074, 0.889)],
        [(0.042, 1.0), (0.042, 0.04), (0.938, 0.04), (0.938, 1.0)],
        [(0.968, 0.04), (0.074, 0.04), (0.074, 0.889), (0.968, 0.889)],
        [(0.074, 0.04), (0.074, 1.0), (0.968, 1.0), (0.968, 0.04)],
        [(0.074, 1.0), (0.074, 0.04), (0.106, 0.101), (0.106, 0.758)],
        [(0.968, 1.0), (0.074, 1.0), (0.106, 0.758), (0.936, 0.758)],
        [(0.968, 0.04), (0.968, 1.0), (0.936, 0.758), (0.936, 0.101)],
        [(0.074, 0.04), (0.968, 0.04), (0.936, 0.101), (0.106, 0.101)],
        [(0.106, 0.04), (0.936, 0.04), (0.936, 0.796), (0.106, 0.796)],
        [(0.042, 0.101), (0.042, 0.758), (0.84, 0.758), (0.84, 0.101)],
        [(0.936, 0.04), (0.106, 0.04), (0.106, 0.796), (0.936, 0.796)],
        [(0.042, 0.758), (0.042, 0.101), (0.84, 0.101), (0.84, 0.758)],
        [(0.936, 0.758), (0.106, 0.758), (0.106, 0.101), (0.936, 0.101)],
        [(0.0, 0.97), (0.041, 0.97), (0.041, 0.811), (0.0, 0.811)],
        [(0.0, 0.811), (0.041, 0.811), (0.041, 0.97), (0.0, 0.97)],
        [(1.0, 1.0), (0.042, 1.0), (0.042, 0.97), (1.0, 0.97)],
        [(0.041, 1.0), (0.041, 0.97), (0.0, 0.97), (0.0, 1.0)],
        [(1.0, 0.0), (0.042, 0.0), (0.042, 0.039), (1.0, 0.039)],
        [(0.041, 0.97), (0.041, 1.0), (0.0, 1.0), (0.0, 0.97)],
        [(0.042, 0.0), (1.0, 0.0), (1.0, 0.039), (0.042, 0.039)],
        [(0.042, 0.97), (0.042, 1.0), (1.0, 1.0), (1.0, 0.97), (1.0, 0.811), (0.042, 0.811)],
        [(1.0, 0.811), (1.0, 0.97), (0.042, 0.97), (0.042, 0.811)]
        ])


def sink(tM, x, y, z, r, sx, sy, verts, faces, matids, uvs):
    r2 = 2 * r + 0.002
    f = len(verts)
    verts.extend([tM * Vector(v) for v in chanfer_square(x, y, z, r, sx + 0.02, sy + 0.02)])
    verts.extend([tM * Vector(v) for v in chanfer_square(x, y, z + 0.002, r, sx + 0.02, sy + 0.02)])
    verts.extend([tM * Vector(v) for v in chanfer_square(x, y, z + 0.002, r, sx - 0.004, sy - 0.004)])
    verts.extend([tM * Vector(v) for v in chanfer_square(x, y, z - 0.2, r, sx - 0.02, sy - 0.02)])
    verts.extend([tM * Vector(v) for v in chanfer_square(x, y + 0.35 * sy, z - 0.2, r, r2, r2)])
    verts.extend([tM * Vector(v) for v in chanfer_square(x, y + 0.35 * sy, z - 0.25, r, r2, r2)])
    s = 16
    faces.extend([(
        f + j * s + i + 1,
        f + j * s + i,
        f + (j + 1) * s + i,
        f + (j + 1) * s + i + 1
        ) for i in range(s - 1) for j in range(5)])
    faces.extend([(
        f + j * s,
        f + (j + 1) * s - 1,
        f + (j + 2) * s - 1,
        f + (j + 1) * s
        ) for j in range(5)])
    matids.extend([mat_sink for i in range(s * 5)])


def cook_top(tM, x, y, z, r, sx, sy, verts, faces, matids, uvs):
    f = len(verts)
    verts.extend([tM * Vector(v) for v in chanfer_square(x, y, z, r, sx + 0.02, sy + 0.02)])
    verts.extend([tM * Vector(v) for v in chanfer_square(x, y, z + 0.002, r, sx + 0.02, sy + 0.02)])
    verts.extend([tM * Vector(v) for v in chanfer_square(x, y, z + 0.002, r, sx, sy)])
    s = 16
    faces.extend([(
        f + j * s + i + 1,
        f + j * s + i,
        f + (j + 1) * s + i,
        f + (j + 1) * s + i + 1
        ) for i in range(s - 1) for j in range(2)])
    faces.extend([(
        f + j * s,
        f + (j + 1) * s - 1,
        f + (j + 2) * s - 1,
        f + (j + 1) * s
        ) for j in range(2)])
    f += 2 * s
    faces.append(tuple([f + i for i in range(s)]))
    matids.extend([mat_cooker_side for i in range(s * 2)])
    matids.append(mat_cooker_top)


def rangehood(tM, x, y, z, z0, th, door_y, verts, faces, matids, uvs):
    cm = x / 60
    x0 = 0
    x1 = th
    x2 = 0.5 * (x - th)
    x3 = 0.5 * (x + th)
    x4 = x - th
    x5 = x

    y0 = 0
    y1 = -th
    y2 = -th - 10 * cm
    y4 = -y - 15 * cm
    y3 = y4 + th

    z1 = z0 + th
    z2 = z0 + z
    f = len(verts)
    verts.extend([tM * Vector(v) for v in [
        (x0, y0, z2), (x0, y4, z2), (x5, y0, z2),
        (x5, y4, z2), (x1, y1, z0), (x1, y3, z0),
        (x4, y3, z0), (x4, y1, z0), (x0, y0, z0),
        (x0, y4, z0), (x5, y4, z0), (x5, y0, z0),
        (x1, y1, z1), (x1, y3, z1), (x4, y3, z1),
        (x4, y1, z1), (x1, y2, z1), (x4, y2, z1),
        (x2, y3, z1), (x2, y2, z1), (x3, y3, z1),
        (x3, y2, z1)
        ]])
    faces.extend([tuple([f + i for i in v]) for v in [
        (0, 1, 3, 2), (6, 5, 13, 18, 20, 14), (0, 2, 11, 8),
        (1, 0, 8, 9), (3, 1, 9, 10), (2, 3, 10, 11),
        (4, 5, 9, 8), (6, 7, 11, 10), (7, 4, 8, 11),
        (5, 6, 10, 9), (21, 17, 14, 20), (7, 6, 14, 17, 15),
        (4, 7, 15, 12), (5, 4, 12, 16, 13), (12, 15, 17, 21, 19, 16),
        (16, 19, 18, 13), (19, 21, 20, 18)
        ]])
    matids.extend([
        mat_range_body, mat_range_body, mat_range_body,
        mat_range_body, mat_range_body, mat_range_body,
        mat_range_body, mat_range_body, mat_range_body,
        mat_range_body, mat_range_filter, mat_range_body,
        mat_range_body, mat_range_body, mat_range_body,
        mat_range_filter, mat_range_body
        ])
    uvs.extend([
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.5, 1.0), (0.933, 0.75), (0.933, 0.25), (0.5, 0.0), (0.067, 0.25), (0.067, 0.75)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.5, 1.0), (0.976, 0.655), (0.794, 0.095), (0.206, 0.095), (0.024, 0.655)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.5, 1.0), (0.976, 0.655), (0.794, 0.095), (0.206, 0.095), (0.024, 0.655)],
        [(0.5, 1.0), (0.933, 0.75), (0.933, 0.25), (0.5, 0.0), (0.067, 0.25), (0.067, 0.75)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        ])


def handle_02(tM, x, y, z, verts, faces, matids, uvs):
    # Rail
    f = len(verts)
    hx = 0.5 * x
    hz = 0.5 * z
    z1 = 0.4 * z
    verts.extend([tM * Vector(v) for v in [
        (hx, 0.0, hz), (hx, 0.0, z1), (hx, -0.9 * y, z1),
        (hx, -y, hz), (hx, -0.9 * y, -z1), (hx, -y, -hz),
        (hx, -0.8 * y, -z1), (hx, -0.8 * y, -hz), (-hx, 0.0, hz),
        (-hx, 0.0, z1), (-hx, -0.9 * y, z1), (-hx, -y, hz),
        (-hx, -0.9 * y, -z1), (-hx, -y, -hz), (-hx, -0.8 * y, -z1),
        (-hx, -0.8 * y, -hz)
        ]])
    faces.extend([tuple([f + i for i in v]) for v in [
        (7, 15, 14, 6), (9, 8, 0, 1), (10, 9, 1, 2),
        (12, 10, 2, 4), (3, 5, 7, 6, 4, 2, 1, 0), (14, 12, 4, 6),
        (15, 7, 5, 13), (3, 11, 13, 5), (8, 11, 3, 0),
        (11, 8, 9, 10, 12, 14, 15, 13)
        ]])
    matids.extend([mat_handle for i in range(10)])
    uvs.extend([
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.5, 1.0), (0.9, 0.9), (1.0, 0.5), (0.9, 0.1),
        (0.5, 0.0), (0.1, 0.1), (0.0, 0.5), (0.1, 0.9)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.5, 1.0), (0.9, 0.9), (1.0, 0.5), (0.9, 0.1),
        (0.5, 0.0), (0.1, 0.1), (0.0, 0.5), (0.1, 0.9)]
    ])


def handle_08(tM, dy, verts, faces, matids, uvs):
    # normalized gap rail inside door
    f = len(verts)
    verts.extend([tM * Vector(v) for v in [
        (-0.5, dy, 0.5), (-0.5, 0.9 * dy, 0.5), (-0.5, 0.9 * dy, -0.389),
        (-0.5, dy, -0.5), (-0.5, 0.1 * dy, -0.389), (-0.5, 0.0, -0.5),
        (-0.5, 0.1 * dy, -0.278), (-0.5, 0.0, -0.278), (0.5, dy, 0.5),
        (0.5, 0.9 * dy, 0.5), (0.5, 0.9 * dy, -0.389), (0.5, dy, -0.5),
        (0.5, 0.1 * dy, -0.389), (0.5, 0.0, -0.5), (0.5, 0.1 * dy, -0.278),
        (0.5, 0.0, -0.278)
        ]])
    faces.extend([tuple([f + i for i in v]) for v in [
        (7, 15, 14, 6), (9, 8, 0, 1), (10, 9, 1, 2),
        (12, 10, 2, 4), (3, 5, 7, 6, 4, 2, 1, 0), (14, 12, 4, 6),
        (15, 7, 5, 13), (3, 11, 13, 5), (8, 11, 3, 0),
        (11, 8, 9, 10, 12, 14, 15, 13)
        ]])
    matids.extend([mat_handle for i in range(10)])
    uvs.extend([
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.5, 1.0), (0.854, 0.854), (1.0, 0.5), (0.854, 0.146),
        (0.5, 0.0), (0.146, 0.146), (0.0, 0.5), (0.146, 0.854)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.5, 1.0), (0.854, 0.854), (1.0, 0.5), (0.854, 0.146),
        (0.5, 0.0), (0.146, 0.146), (0.0, 0.5), (0.146, 0.854)]
        ])


def handle_09(tM, dy, y, verts, faces, matids, uvs):
    # normalized rail on top of door
    # dy: door depth
    # y: depth from door front
    # use matrix to scale in x/z
    # y >= 0 inside door
    f = len(verts)
    verts.extend([tM * Vector(v) for v in [
        (0.5, dy, -0.5), (0.5, 0.0, -0.5), (0.5, 0.0, -0.4),
        (0.5, dy, 0.5),
        (0.5, -0.515 * y, -0.4), (0.5, -0.515 * y, -0.5),
        (-0.5, dy, -0.5), (-0.5, 0.0, -0.5), (-0.5, 0.0, -0.4),
        (-0.5, dy, 0.5),
        (-0.5, -0.515 * y, -0.4), (-0.5, -0.515 * y, -0.5),
        (0.5, -0.7 * y, 0.5), (0.5, -y, 0.2), (0.5, -y, -0.2),
        (0.5, -0.7 * y, -0.5), (-0.5, -y, 0.2), (-0.5, -0.7 * y, 0.5),
        (-0.5, -0.7 * y, -0.5), (-0.5, -y, -0.2)
        ]])
    faces.extend([tuple([f + i for i in v]) for v in [
        (9, 6, 7, 8, 10, 11, 18, 19, 16, 17), (7, 6, 0, 1), (8, 7, 1, 2),
        (10, 8, 2, 4), (19, 14, 13, 16), (11, 10, 4, 5),
        (15, 18, 11, 5), (3, 12, 13, 14, 15, 5, 4, 2, 1, 0), (6, 9, 3, 0),
        (14, 19, 18, 15), (12, 17, 16, 13), (3, 9, 17, 12)
        ]])
    matids.extend([mat_handle for i in range(10)])
    uvs.extend([
        [(0.5, 1.0), (0.794, 0.905), (0.976, 0.655), (0.976, 0.345), (0.794, 0.095),
        (0.5, 0.0), (0.206, 0.095), (0.024, 0.345), (0.024, 0.655), (0.206, 0.905)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.5, 1.0), (0.794, 0.905), (0.976, 0.655), (0.976, 0.345), (0.794, 0.095),
        (0.5, 0.0), (0.206, 0.095), (0.024, 0.345), (0.024, 0.655), (0.206, 0.905)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        ])


def handle_12(tM, x, y, px, pz, rx, ry, sx, sz, s0, s1, start_angle, angle_profil, rosace, verts, faces, matids, uvs):
    """
      Full parametric Rounded pipes
        x: width
        y: depth
        prx, prz: pipe radius
        rosace: boolean base
        rx: radius in x
        ry: radius in y
        sx: scale along x from start to center
        sz: scale along z from start to center
        s0: segments along half path, minimum 1
        s1: profile section minimum 3
        a0: start angle for sections
        start_angle: start / end angle [0 - pi/2]
    """
    # rounded pipe |_|
    prx = 0.5 * px
    prz = 0.5 * pz

    rx = min(0.4999 * x, rx)
    ry = min(y, ry)

    f = len(verts)
    dx = x - 2 * rx
    a0 = pi + start_angle
    da = (pi - 2 * start_angle) / (2 * s0)
    a1 = angle_profil
    da1 = 2 * pi / s1
    xc = rx - 0.5 * x
    yc = ry - y
    if rosace:
        dy = 0
        for j in range(2):
            for i in range(s1):
                verts.append(tM *
                Vector((xc - rx - cos(a1) * prx * (0.5 + sx), dy, sin(a1) * prz * (0.5 + sz)))
                )
                a1 += da1
            dy -= 0.25 * prz
        for i in range(s1):
            verts.append(tM *
            Vector((xc - rx - cos(a1) * prx * sx, dy, sin(a1) * prz * sz))
            )
            a1 += da1
    else:
        for i in range(s1):
            verts.append(tM *
            Vector((xc - rx - cos(a1) * prx * sx, 0, sin(a1) * prz * sz))
            )
            a1 += da1
    scx = sx
    scz = sz
    dsx = (sx - 1) / (s0 - 1)
    dsz = (sz - 1) / (s0 - 1)
    if s0 > 2:
        s0 += 1
    for j in range(2):
        for n in range(s0):
            ca = cos(a0)
            sa = sin(a0)
            for i in range(s1):
                vxy = scx * cos(a1) * prx
                vz = scz * sin(a1) * prz
                verts.append(tM *
                Vector((xc + ca * (rx + vxy), yc + sa * (ry + vxy), vz))
                )
                a1 += da1
            scx -= dsx
            scz -= dsz
            a0 += da
        scx = sx - (s0 - 1) * dsx
        scz = sz - (s0 - 1) * dsz
        a0 = 2 * pi - start_angle - (s0 - 1) * da
        dsx = -dsx
        dsz = -dsz
        xc += dx
    xc -= dx
    if rosace:
        for i in range(s1):
            verts.append(tM *
            Vector((xc + rx + cos(a1) * prx * sx, dy, sin(a1) * prz * sz))
            )
            a1 += da1
        for j in range(2):
            dy += 0.25 * prz
            for i in range(s1):
                verts.append(tM *
                Vector((xc + rx + cos(a1) * prx * (0.5 + sx), dy, sin(a1) * prz * (0.5 + sz)))
                )
                a1 += da1
        s = 2 * s0 + 5
    else:
        s = 2 * s0 + 1
        for i in range(s1):
            verts.append(tM *
            Vector((xc + rx + cos(a1) * prx * sx, 0, sin(a1) * prz * sz))
            )
            a1 += da1
    faces.extend([(
        f + j * s1 + i + 1,
        f + j * s1 + i,
        f + (j + 1) * s1 + i,
        f + (j + 1) * s1 + i + 1
        ) for i in range(s1 - 1) for j in range(s)])
    faces.extend([(
        f + j * s1,
        f + (j + 1) * s1 - 1,
        f + (j + 2) * s1 - 1,
        f + (j + 1) * s1
        ) for j in range(s)])
    matids.extend([mat_handle for i in range(s * s1)])
    uvs.extend([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] for i in range(s * s1)])


def handle_14(tM, x, y, r, verts, faces, matids, uvs):
    # modern T pipes
    f = len(verts)
    r1 = 0.75 * r
    s = 18
    da = 2 * pi / s
    a = 0
    for i in range(2):
        dx = (1 - 2 * i) * 0.5 * (x - y)
        for j in range(2):
            for k in range(s):
                a -= da
                verts.append(tM * Vector((dx + cos(a) * r1, (j - 1) * y, sin(a) * r1)))
    for i in range(2):
        dx = 0.5 * (1 - 2 * i) * x
        for k in range(s):
            a -= da
            verts.append(tM * Vector((dx, -y + cos(a) * r, sin(a) * r)))
    for j in range(3):
        faces.extend([tuple([f + i + v for v in (0, 1, s + 1, s)]) for i in range(s - 1)])
        faces.append((f + s - 1, f, f + s, f + 2 * s - 1))
        faces.append(tuple([f + s - 1 - i for i in range(s)]))
        faces.append(tuple([f + s + i for i in range(s)]))
        f += 2 * s
        uvs.extend([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] for i in range(s)])
        for k in range(2):
            uvs.append([(cos(da * i), sin(da * i)) for i in range(s)])
        matids.extend([mat_handle for i in range(s + 2)])


def handle_15(tM, profil, s, verts, faces, matids, uvs):
    # Buttons
    f = len(verts)
    s1 = len(profil) - 1
    da = 2 * pi / s
    for x, y in profil:
        a = pi / 4
        for k in range(s):
            a -= da
            verts.append(tM * Vector((cos(a) * x, -y, sin(a) * x)))
    faces.extend([(
        f + j * s + i + 1,
        f + j * s + i,
        f + (j + 1) * s + i,
        f + (j + 1) * s + i + 1
        ) for i in range(s - 1) for j in range(s1)])
    faces.extend([(
        f + j * s,
        f + (j + 1) * s - 1,
        f + (j + 2) * s - 1,
        f + (j + 1) * s
        ) for j in range(s1)])
    f += s * (s1 + 1) - 1
    faces.append(tuple([f - i for i in range(s)]))
    uvs.extend([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] for i in range(s * s1)])
    uvs.append([(cos(da * i), sin(da * i)) for i in range(s)])
    matids.extend([mat_handle for i in range(s * s1 + 1)])


def handle_05(tM, verts, faces, matids, uvs):
    f = len(verts)
    verts.extend([tM * Vector(v) for v in [

        ]])
    faces.extend([tuple([f + i for i in v]) for v in [
        ]])
    matids.extend([2 for i in range(9)])
    uvs.extend([
        ])


# ------------------------------------------------------------------
# Define property group class for cabinet doors
# This is managed as an array of objects
# ------------------------------------------------------------------
class archipack_kitchen_module(ArchipackObject, PropertyGroup):
    """
     Vertical module inside kitchen
    """
    type = EnumProperty(
            name="Type",
            items=(
                ("1", "Drawer", "Drawer"),
                ("2", "Door L", "Door Left"),
                ("3", "Door R", "Door Right"),
                ("4", "Door T", "Door Top"),
                # ("6", "Door B", "Door Bottom"),
                ("5", "Double door", "Double Door"),
                ("6", "Board only", "Board only"),
                ("12", "Door L w glass", "Door Left with glass"),
                ("13", "Door R w glass", "Door Right with glass"),
                ("14", "Door T w glass", "Door Top with glass"),
                ("15", "Double w glass", "Double with glass"),
                ("54", "Oven", "Oven"),
                ("55", "Range hood", "Range hood"),
                # ("56", "Microwave", "Microwave"),
                # ("58", "Freezer", "Freezer"),
                ("59", "Dishwasher", "Dishwasher"),
                ("0", "None", "None"),
                ),
            default="1",
            update=update_parent
            )
    modules = IntProperty(
            name="Height",
            description="#modules 0=user defined 1=12.7 cm - 5",
            min=0, default=6,
            update=update_parent
            )
    z = FloatProperty(
            name="Height",
            min=0.01, default=0.762,
            unit='LENGTH', subtype='DISTANCE',
            update=update_parent
            )
    shelves = IntProperty(
            name='Shelves', min=0, default=1,
            description='Number total of shelves',
            update=update_parent,
            )
    handle = BoolProperty(
            name="Handle",
            description="Create a handle", default=True,
            update=update_parent,
            )
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            name="auto_update",
            description="disable auto update to avoid infinite recursion",
            default=True
            )

    @property
    def n_panels(self):
        if self.type in {'5', '15'}:
            return 2
        elif self.type == '0':
            return 0
        return 1

    def draw(self, layout, num, z_mode):
        module_type = int(self.type)

        row = layout.row(align=True)
        row.prop(self, 'type', text="Door " + str(num))
        if module_type < 50:
            row.prop(self, 'handle', text="")

        row = layout.row(align=True)
        if z_mode == 1:
            row.prop(self, 'modules')

        if z_mode > 1 or self.modules == 0:
            row.prop(self, 'z')

        if module_type < 50 and module_type not in {1}:
            row.prop(self, 'shelves')

    def door(self, size, kitchen, style, mat, mat_ext):
        border, w = kitchen.door_x, kitchen.door_y
        chanfer, board_chanfer = kitchen.door_chanfer, 0.5 * kitchen.door_board_chanfer
        x0 = 0
        x1 = border
        x2 = x1 - chanfer
        x3 = chanfer
        x4 = x1 + board_chanfer
        x5 = x4 + board_chanfer
        # offset pivot point on outside part
        y0 = 0
        y1 = y0 + w
        y2 = y1 - 0.5 * w
        y3 = y1 - chanfer
        y4 = y0 + chanfer
        y5 = y2 - tan22_5 * board_chanfer
        verre = 0.002

        # profil carre avec support pour verre
        # p ______       y1
        # / |      y3
        # |       |___
        # x       |___   y2  verre
        # |       |      y4
        #  \______|      y0
        # x0 x3 x2 x1
        # outside

        # flat board
        if style == 1:
            return Lofter(
                False,  # closed
                [0, 0],  # x index
                [x0],
                [y0, y1],
                [mat, mat],  # materials
                side_cap_front=1,
                side_cap_back=0      # cap index
                )

        # board with frame
        elif style == 2:
            return Lofter(
                False,  # closed
                [1, 1, 0, 0],  # x index
                [x0, x1],
                [y2, y0, y0, y1],
                [mat_border for i in range(4)],  # materials
                side_cap_front=3,
                side_cap_back=0      # cap index
                )

        # board with chanfer
        elif style == 11:
            return Lofter(
                False,  # closed
                [1, 0, 0, 1],  # x index
                [x0, x3],
                [y0, y4, y3, y1],
                [mat, mat, mat, mat],  # materials
                side_cap_front=3,
                side_cap_back=0      # cap index
                )

        # board with chanfer and frame
        elif style == 12:
            # fallback to panel without border
            # when border > size
            return Lofter(
                False,  # closed
                [3, 3, 2, 1, 0, 0],  # x index
                [x0, x3, x2, x1],
                [y2, y4, y0, y0, y4, y1],
                [mat_border for i in range(6)],  # materials
                side_cap_front=5,
                side_cap_back=0      # cap index
                )
        # board with chanfer and chanfered frame
        elif style == 13:
            # fallback to panel without border
            # when border > size
            return Lofter(
                False,  # closed
                [5, 4, 3, 3, 2, 1, 0, 0],  # x index
                [x0, x3, x2, x1, x4, x5],
                [y5, y2, y2, y4, y0, y0, y4, y1],
                [mat_ext, mat_ext,
                    mat_border, mat_border, mat_border,
                    mat_border, mat_border, mat_border],  # materials
                side_cap_front=7,
                side_cap_back=0      # cap index
                )
        # board with glass
        elif style == 21:
            return Lofter(
                True,  # closed
                [0, 0, 1, 1, 1, 1],  # x index
                [x0, x1],
                [y0, y1, y1, y2 + verre, y2 - verre, y0],
                [mat for i in range(6)],  # materials
                side_cap_front=3,
                side_cap_back=4      # cap index
                )

        # board with chanfer and glass
        elif style == 22:
            # fallback to panel without border
            # when border > size
            return Lofter(
                True,  # closed
                [1, 0, 0, 2, 2, 2, 2, 3],  # x index
                [x0, x3, x1, x2],
                [y0, y4, y1, y1, y2 + verre, y2 - verre, y4, y0],
                [mat for i in range(8)],  # materials
                side_cap_front=4,
                side_cap_back=5      # cap index
                )

    def find_datablock_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_kitchen.datablock(o)
            if props:
                for cabinet in props.cabinets:
                    for module in cabinet.modules:
                        if module == self:
                            return props
        return None

    def from_module(self, module):
        self.auto_update = False
        self.type = module.type
        self.modules = module.modules
        self.z = module.z
        self.shelves = module.shelves
        self.handle = module.handle
        self.auto_update = True

    def update_parent(self, context):

        if self.auto_update:
            props = self.find_datablock_in_selection(context)
            if props is not None:
                props.update(context)

    def update(self, context, size, pivot, kitchen, cab, module, handle_location):
        """
          kitchen_module may produce more than one child object
          so not every object's kitchen_module belong to cabs
          this update func dosent rely on this particular datablock
        """
        o = self.find_in_selection(context)

        if o is None:
            return

        module_type = int(module.type)
        style = int(kitchen.door_style)
        border = kitchen.door_x
        if style == 13:
            border += kitchen.door_board_chanfer
        handle_style = int(kitchen.handle)
        cab_location = cab.location
        th = kitchen.thickness

        if 10 < module_type < 20:
            # module with glass
            style = 20 + style % 10

        mat_side, mat_ext, mat_int = cab_location, cab_location, mat_inside
        border_overflow = 2 * border > size.y or 2 * border > size.x

        if border_overflow:
            # style with frame
            if style % 10 > 1:
                # face is of border mat
                mat_ext = mat_border
                mat_side = mat_border
            # fallback to panel without border
            if style in {2, 21}:
                style = 1
            elif style in {12, 13, 22}:
                style = 11

        # glass and no overflow
        if style > 20:
            mat_ext, mat_int = 8, 8

        if module_type > 50:
            verts, faces, matids, uvs = [], [], [], []

            if module_type == 54:
                oven_door(size, verts, faces, matids, uvs)
            elif module_type == 59:
                # might use resized panel and add part on top
                dishwasher_door(size, verts, faces, matids, uvs)

        else:

            offset = Vector((0, 0, 0))
            center = Vector((0, 0, 0))
            origin = Vector((0, 0, 0))
            radius = Vector((0, 0, 0))
            board_size = size.copy()

            # negative handle follow special rule
            # when door size > bottom size handle goes on panel side
            # otherwise handle on panel top / bottom
            handle_on_side = module_type not in {1, 4}
            if module.handle and handle_style in {180, 190}:
                # handles inside top/side of doors
                if size.y <= kitchen.height_default:
                    handle_on_side = False
                    board_size.y -= kitchen.handle_z
                    # offset door model for bottom handle
                    if handle_location < 0:
                        offset.y = kitchen.handle_z
                else:
                    handle_on_side = True
                    board_size.x -= kitchen.handle_z

            door = self.door(board_size, kitchen, style, mat_side, mat_ext)

            verts = door.vertices(1, offset, center, origin, board_size,
                radius, 0, pivot, shape_z=None, path_type='RECTANGLE')

            faces = door.faces(1, path_type='RECTANGLE')

            uvs = door.uv(1, center, origin, board_size,
                radius, 0, pivot, border, 0, path_type='RECTANGLE')

            matids = door.mat(1, mat_int, mat_ext, path_type='RECTANGLE')

            # -------------
            # Drawer
            # -------------
            if module_type == 1:
                # drawer
                sy = kitchen.y + cab.dy - th
                if cab_location == 2:
                    sy = kitchen.yw + cab.dy - th
                x0 = 0.5 * size.x * pivot
                y0 = kitchen.door_y
                y1 = y0 + sy - th
                sx = 0.5 * size.x - th
                z0 = 2 * th
                z1 = 0.5 * size.y
                tM = Matrix()
                # sides
                make_box(tM, x0 + sx - th, x0 + sx, y1, y0, z0, z1, 0,
                    mat_inside, False, verts, faces, matids, uvs)
                make_box(tM, x0 - sx, x0 + th - sx, y1, y0, z0, z1, 0,
                    mat_inside, False, verts, faces, matids, uvs)
                # back
                make_box(tM, x0 - sx, x0 + sx, y1 + th, sy, z0, z1, 0,
                    mat_inside, False, verts, faces, matids, uvs)
                # bottom
                make_box(tM, x0 - sx, x0 + sx, y1 + th, y0, z0 - th, z0, 0,
                    mat_inside, False, verts, faces, matids, uvs)

            # -------------
            # Handles
            # -------------
            if module_type != 6 and module_type < 50 and module.handle:

                # offset from borders
                x, y, z = kitchen.handle_x, kitchen.handle_y, kitchen.handle_z
                r = 0.5 * kitchen.handle_r

                # limit handles offset
                h = r
                if handle_style < 150:
                    # bars
                    r = 0.5 * x
                    h = 0.5 * z

                offset_z = max(kitchen.handle_dz, r)
                offset_x = max(kitchen.handle_dx, h)

                if handle_on_side:
                    handle_offset_x = min(offset_x, size.x - h)
                    handle_offset_z = min(offset_z, size.y - r)
                else:
                    handle_offset_x = min(offset_x, size.y - h)
                    handle_offset_z = min(offset_z, size.x - r)

                s = Vector((1, 1, 1))

                # [-1: bottom, 0:center, 1:top]
                zs = handle_location

                if module_type in {1, 4}:
                    # drawer / door top - bottom
                    handle_offset = zs * handle_offset_x
                else:
                    handle_offset = zs * handle_offset_z

                dz = 0.5 * (1 + zs) * size.y - handle_offset
                # handles inside doors top/bottom/side
                if handle_style in {180, 190}:
                    # offset is half handle z size
                    handle_offset_x = 0.5 * min(kitchen.handle_z, size.x, size.y)
                    if zs == 0:
                        zs = 1
                    if handle_on_side:
                        # z offset
                        dz = 0.5 * size.y
                        # force side / top-bottom rules
                        module_type = 0
                        # scale to fit door size, use raw door_y
                        s = Vector((zs * 2 * handle_offset_x, 1, zs * size.y))
                    else:
                        # z offset
                        dz = 0.5 * (1 + zs) * size.y - zs * handle_offset_x
                        # force side / top-bottom rules
                        module_type = 1
                        # scale to fit door size, use raw door_y
                        s = Vector((zs * size.x, 1, zs * 2 * handle_offset_x))

                # make handle fit
                # scale to fit door size
                elif kitchen.handle_fit and handle_style < 150:
                    if zs == 0:
                        zs = 1
                    if module_type in {1, 4}:
                        dz = 0.5 * (1 + zs) * size.y - zs * handle_offset_x
                        x = max(x, size.x - kitchen.handle_space)
                    else:
                        dz = 0.5 * size.y
                        x = max(x, size.y - kitchen.handle_space)
                else:
                    if handle_on_side:
                        x = min(x, size.y - h)
                    else:
                        x = min(x, size.x - h)

                if module_type in {1, 4}:
                    # drawer / door top - bottom
                    # Drawer - horizontal handle
                    tM = Matrix([
                        [s.x, 0, 0, pivot * (0.5 * size.x)],
                        [0, s.y, 0, 0],
                        [0, 0, s.z, dz],
                        [0, 0, 0, 1]
                    ])
                else:

                    # Door - vertical handle
                    tM = Matrix([
                        [0, 0, pivot * s.x, pivot * (size.x - handle_offset_x)],
                        [0, s.y, 0, 0],
                        [-pivot * s.z, 0, 0, dz],
                        [0, 0, 0, 1]
                    ])

                if handle_style == 2:
                    handle_02(tM, x, y, z, verts, faces, matids, uvs)

                # parametric rails
                elif handle_style == 180:
                    handle_08(tM, kitchen.door_y, verts, faces, matids, uvs)

                elif handle_style == 190:
                    handle_09(tM, kitchen.door_y, y, verts, faces, matids, uvs)

                # Loft based parametric bars
                elif handle_style == 120:
                    # 120 Modern tri bar Trida
                    handle_12(tM, x, y, z, z, 0, z, 0.707, 1, 2, 3, 0, pi, False,
                        verts, faces, matids, uvs)

                elif handle_style == 121:
                    # 121 Rounded square bar Borghamn
                    handle_12(tM, x, y, z, z, z, z, 1, 1, 9, 4, 0, pi / 4, False,
                        verts, faces, matids, uvs)

                elif handle_style == 122:
                    # 122 Old school rounded pipe bar Mollarp
                    handle_12(tM, x, y, z, z, 4 * z, 4 * z, 0.7, 0.7, 9, 18, 0, 0, True,
                        verts, faces, matids, uvs)

                elif handle_style == 123:
                    # 123 Modern Pipe bar Bagganas
                    handle_12(tM, x, y, z, z, 1.5 * z, 1.5 * z, 1, 1, 9, 18, 0, 0, True,
                        verts, faces, matids, uvs)

                elif handle_style == 124:
                    # 124 Modern Square Bar
                    handle_12(tM, x, y, z, z, 0, z, 0.707, 1, 2, 4, 0, pi / 4, False,
                        verts, faces, matids, uvs)

                elif handle_style == 125:
                    # 125 Old school elliptic
                    handle_12(tM, x, y, z, 2 * z, 0.45 * x, y, 1, 0.5, 9, 18, 0, pi / 4, False,
                        verts, faces, matids, uvs)

                elif handle_style == 126:
                    # 126 Elliptic bar
                    handle_12(tM, x, y, z, z, 0.45 * x, 0.3 * y, 1.5, 1, 9, 4, 0, pi / 4, False,
                        verts, faces, matids, uvs)

                elif handle_style == 127:
                    # 127 Old school rounded pipe bar 2 Eneryda
                    handle_12(tM, x, y, z * 0.666, z, 0.3 * x, y - 2 * z, 1, 0.5, 9, 18, 0, 0, True,
                        verts, faces, matids, uvs)

                # parametric t pipe bars
                elif handle_style == 140:
                    # 140 T pipe bar
                    r = 0.5 * z
                    handle_14(tM, x, y, r, verts, faces, matids, uvs)

                # Lathe based parametric handles
                elif handle_style == 150:
                    # 150 Circle handle
                    profil = [
                        (r * 0.3, 0),
                        (r * 0.3, 0.9 * y),
                        (r, 0.9 * y),
                        (r, y)
                        ]
                    handle_15(tM, profil, 18, verts, faces, matids, uvs)

                elif handle_style == 151:
                    # 151 Square handle
                    profil = [
                        (r * 0.3, 0),
                        (r * 0.3, 0.9 * y),
                        (r, 0.9 * y),
                        (r, y)
                        ]
                    handle_15(tM, profil, 4, verts, faces, matids, uvs)

                elif handle_style == 152:
                    # 152 Spherical handle
                    s = 18
                    y = max(2 * r, y)
                    da = pi / s
                    a = -pi / 2
                    r0 = r * cos(a + da)
                    profil = [(3.5 * r0, 0), (3.5 * r0, 0.03 * y), (2.5 * r0, 0.06 * y)]
                    a += da
                    for i in range(s - 2):
                        a += da
                        profil.append((cos(a) * r, y - z + sin(a) * r))
                    handle_15(tM, profil, s, verts, faces, matids, uvs)

                elif handle_style == 153:
                    # 153 Elliptic handle
                    s = 18
                    y = max(1.8 * r, y)
                    da = pi / s
                    a = -pi / 2
                    r1 = 0.2 * r
                    profil = [(0.5 * r + cos(a - 2 * da) * r1, 0)]

                    cy = y - sin(7 * da) * r - 2 * r1
                    a -= 2 * da
                    for i in range(2, s):
                        profil.append((
                            0.5 * r + cos(a) * r1,
                            cy + sin(a) * 2 * r1
                            ))
                        a -= da

                    a += 4 * da
                    cy = y - 0.5 * r
                    for i in range(s - 3):
                        profil.append((
                            -cos(a) * r,
                            cy - 0.5 * sin(a) * r
                            ))
                        a += da

                    handle_15(tM, profil, s, verts, faces, matids, uvs)

                elif handle_style == 154:
                    # 154 Mushroom like handle
                    profil = [
                        (0.5 * r, 0),
                        (0.5 * r, 0.5 * y),
                        (r, 0.8 * y),
                        (r, y)
                        ]
                    handle_15(tM, profil, 18, verts, faces, matids, uvs)

                elif handle_style == 155:
                    # 155 Mushroom like handle 2
                    profil = [
                        (0.55 * r, 0),
                        (0.55 * r, 0.2 * y),
                        (0.4 * r, 0.22 * y),
                        (0.4 * r, 0.7 * y),
                        (r, 0.8 * y),
                        (r, y)
                        ]
                    handle_15(tM, profil, 18, verts, faces, matids, uvs)

        bmed.buildmesh(context, o, verts, faces, matids, uvs)

        self.restore_context(context)


class archipack_kitchen_cabinet(ArchipackObject, PropertyGroup):
    # Define properties
    type = EnumProperty(
            # 0 for regular, 1 corner L, 2 corner R, 3 corner L+R
            items=(
                ('0', "Cabinet", ""),
                ('1', "Corner L", ""),
                ('2', "Corner R", ""),
                ('3', "Corner 45 degree", ""),
                ('4', "Cabinet without left side", ""),
                ('5', "Cabinet without right side", ""),
                ('6', "Cabinet without sides", ""),
                ('7', "Corner 45 degree without sides", ""),
                ('8', "Empty space", "Empty space (with countertop)"),
                ),
            name="Type",
            default="0",
            description="Type of cabinet",
            update=update_manipulators
            )
    cab_location = EnumProperty(
            items=(
                ('1', "Floor", ""),
                ('2', "Wall", ""),
                ('3', "Full", ""),
                ),
            name="Location",
            default="1",
            description="Location of cabinet",
            update=update
            )
    # Cabinet width
    x = FloatProperty(
            name='Width', min=0.001, default=0.60, precision=3,
            description='Cabinet width',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    dy = FloatProperty(
            name='y', default=0, precision=3,
            description='Modify depth',
            unit='LENGTH', subtype='DISTANCE',
            update=update,
            )
    dz = FloatProperty(
            name='z', default=0, precision=3,
            description='Modify height',
            unit='LENGTH', subtype='DISTANCE',
            update=update,
            )

    # Cabinet position shift
    reset_location = BoolProperty(
            name="Reset",
            description="Reset location", default=False,
            update=update,
            )
    px = FloatProperty(
            name='x', default=0, precision=3,
            description='Position x shift',
            unit='LENGTH', subtype='DISTANCE',
            update=update,
            )
    py = FloatProperty(
            name='y', default=0, precision=3,
            description='Position y shift',
            unit='LENGTH', subtype='DISTANCE',
            update=update,
            )
    pz = FloatProperty(
            name='z', default=0, precision=3,
            description='Position z shift',
            unit='LENGTH', subtype='DISTANCE',
            update=update,
            )
    lock_p = BoolProperty(
        name="Lock",
        description="Move next cabs too",
        default=True,
        update=update
        )

    # Doors
    n_modules = IntProperty(
        name="Doors",
        min=0, default=1,
        update=update
        )
    modules = CollectionProperty(type=archipack_kitchen_module)
    expand = BoolProperty(
            options={'SKIP_SAVE'},
            name="expand",
            description="Expand cabinet panel",
            default=False
            )

    # baseboard
    baseboard = BoolProperty(
            name="Baseboard",
            description="Create a baseboard automatically",
            default=True,
            update=update
            )
    base_left = BoolProperty(
            name="Left",
            description="Create a left baseboard", default=False,
            update=update,
            )
    base_right = BoolProperty(
            name="Right",
            description="Create a right baseboard", default=False,
            update=update,
            )
    base_sink = FloatProperty(
        name='sink', min=0,
        default=0.05, precision=3,
        description='Baseboard side sink',
        unit='LENGTH', subtype='DISTANCE',
        update=update
        )
    base_front = FloatProperty(
        name='front',
        default=0.0, precision=3,
        description='Baseboard front delta sink',
        unit='LENGTH', subtype='DISTANCE',
        update=update
        )

    # side boards
    panel_left_width = FloatProperty(
            name="L board",
            description="Width of left side board",
            min=0.001, default=0.018,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    panel_right_width = FloatProperty(
            name="R board",
            description="Width of right side board",
            min=0.001, default=0.018,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    panel_left = BoolProperty(
            name="Left",
            description="Create a left side board", default=False,
            update=update,
            )
    panel_right = BoolProperty(
            name="Right",
            description="Create a right side board", default=False,
            update=update,
            )

    # Countertop
    counter = EnumProperty(
            name="Countertop",
            items=(
                ('0', 'No Countertop', 'Disable countertop'),
                ('1', 'Countertop', 'Countertop regular'),
                ('2', 'Cook top', 'Countertop with cook top'),
                # 10+ create hole in countertop
                ('11', 'Sink', 'Countertop with sink'),
                ('12', 'Countertop hole', 'Countertop with hole'),
                ),
            default='1',
            update=update,
            description="Create counter"
            )
    counter_x = FloatProperty(
            name="Width",
            description="Width",
            min=0.001, default=0.5,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    counter_y = FloatProperty(
            name="Length",
            description="Length",
            min=0.001, default=0.4,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    # cabinet rotation
    rotate = EnumProperty(
            items=(
                ('0', "No rotation", ""),
                ('1', "90 CW", ""),
                ('2', "90 CCW", ""),
                ('3', "180", ""),
                ('4', "User defined", "")
                ),
            name="Rot",
            description="Rotate cabinet relative to previous one",
            update=update,
            )
    angle = FloatProperty(
            name="Angle",
            description="User defined rotation",
            default=0,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )

    manipulators = CollectionProperty(type=archipack_manipulator)
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            name="auto_update",
            description="disable auto update to avoid infinite recursion",
            default=True
            )
            
    # DimensionProvider        
    uid = IntProperty(default=0)
    
    @property
    def location(self):
        """
         also define material id for doors and sides
         1 for floor, 2 for wall, 3 for full
        """
        return int(self.cab_location)

    @property
    def cab_type(self):
        """
         0 for regular, 1 corner L, 2 corner R, 3 corner L+R
        """
        return int(self.type)

    @property
    def countertop_hole(self):
        """ Test if hole is needed in countertop
        """
        return int(self.counter) > 10 and self.location == 1

    @property
    def board_left(self):
        """ Space in use for left board and """
        if self.panel_left and self.cab_type != 2:
            return self.panel_left_width
        else:
            return 0

    @property
    def board_right(self):
        """ Space in use for right board"""
        if self.panel_right and self.cab_type != 1:
            return self.panel_right_width
        else:
            return 0

    def from_cabinet(self, cab):
        self.auto_update = False
        self.type = cab.type
        self.cab_location = cab.cab_location
        self.x = cab.x
        self.dy = cab.dy
        self.dz = cab.dz
        self.baseboard = cab.baseboard
        self.counter = cab.counter
        self.counter_x = cab.counter_x
        self.counter_y = cab.counter_y
        self.panel_left = cab.panel_left
        self.panel_left_width = cab.panel_left_width
        self.panel_right = cab.panel_right
        self.panel_right_width = cab.panel_right_width
        self.n_modules = cab.n_modules
        for module in cab.modules:
            m = self.modules.add()
            m.from_module(module)
        self.auto_update = True

    def find_datablock_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_kitchen.datablock(o)
            if props:
                for cab in props.cabinets:
                    if cab == self:
                        return props
        return None

    def update_parts(self):
        # remove modules
        for i in range(len(self.modules), self.n_modules, -1):
            self.modules.remove(i - 1)

        # add modules
        for i in range(len(self.modules), self.n_modules):
            self.modules.add()

    def update(self, context, manipulable_refresh=False):
        if self.auto_update:
            props = self.find_datablock_in_selection(context)
            if props is not None:
                props.update(context, manipulable_refresh)

    def setup_manipulators(self, i):
        n_manips = len(self.manipulators)

        if n_manips < 1:
            s = self.manipulators.add()
            s.type_key = 'DUMB_STRING'
        self.manipulators[0].prop1_name = str(i + 1)

        if n_manips < 2:
            s = self.manipulators.add()
            s.type_key = 'SIZE'
            s.prop1_name = 'x'

        if n_manips < 3:
            # dumb oposite side
            s = self.manipulators.add()
            s.type_key = 'DUMB_SIZE'

        if n_manips < 4:
            # dumb door size
            s = self.manipulators.add()
            s.type_key = 'DUMB_SIZE'

    def draw(self, box, num, prop):

        row = box.row(align=True)

        if self.expand:
            row.prop(self, "expand", icon="TRIA_DOWN", icon_only=True, text="Cab " + str(num + 1), emboss=False)
        else:
            row.prop(self, "expand", icon="TRIA_RIGHT", icon_only=True, text="Cab " + str(num + 1), emboss=False)
        row.prop(self, 'cab_location', text="")
        row.prop(self, 'type', text="")
        row.operator("archipack.kitchen_insert", icon="ZOOMIN", text="").index = num
        row.operator("archipack.kitchen_remove", icon="ZOOMOUT", text="").index = num

        if self.expand:
            cab_type = self.cab_type
            row = box.row(align=True)
            row.prop(self, 'x')
            row.prop(self, 'dy')
            row.prop(self, 'dz')

            row = box.row(align=True)
            row.prop(self, 'px')
            row.prop(self, 'py')
            row.prop(self, 'pz')
            if self.lock_p:
                row.prop(self, 'lock_p', icon="LOCKED", text="")
            else:
                row.prop(self, 'lock_p', icon="UNLOCKED", text="")

            row = box.row(align=True)
            row.prop(self, 'reset_location')
            row.prop(self, 'rotate', text="")
            if self.rotate == '4':
                row.prop(self, 'angle', text="")

            if self.location == 1:
                if prop.counter:
                    row = box.row(align=True)
                    row.prop(self, 'counter', text="")
                    if int(self.counter) > 1:
                        row = box.row(align=True)
                        row.prop(self, 'counter_x')
                        row.prop(self, 'counter_y')

            if self.location != 2 and cab_type != 8:
                if prop.baseboard:
                    row = box.row(align=True)
                    row.prop(self, 'baseboard', text="Baseboard")
                    if self.baseboard:
                        row.prop(self, 'base_front', text="")
                        row.prop(self, 'base_sink', text="")
                        row.prop(self, 'base_left', text="")
                        row.prop(self, 'base_right', text="")

            if cab_type != 8:
                row = box.row(align=True)
                row.prop(self, 'panel_left_width')
                row.prop(self, 'panel_left', text="")
                row.prop(self, 'panel_right_width')
                row.prop(self, 'panel_right', text="")

                row = box.row()
                row.prop(self, 'n_modules')
                for i, module in enumerate(reversed(self.modules)):
                    module.draw(box, self.n_modules - i, int(prop.z_mode))


class archipack_kitchen(ArchipackObject, Manipulable, DimensionProvider, PropertyGroup):

    door_style = EnumProperty(
            name='Doors',
            items=(
                ('1', 'Flat', 'Simple board without chanfer'),
                ('2', 'Flat with frame', 'Board with frame'),
                ('11', 'Chanfer', 'Simple board with chanfer'),
                ('12', 'Frame with chanfer, flat board', 'Board with frame and chanfer'),
                ('13', 'Frame and board with chanfer', 'Board with frame and chanfer'),
                ),
            default='11',
            update=update
            )
    door_x = FloatProperty(
            name='Border', min=0.001, default=0.08, precision=3,
            description='Door border width',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    door_y = FloatProperty(
            name='Thick', min=0.001, default=0.018, precision=3,
            description='Door thickness',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    door_chanfer = FloatProperty(
            name='Chanfer', min=0, default=0.001, precision=3,
            description='Door chanfer',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    door_gap = FloatProperty(
            name='Gap', min=0, default=0.002, precision=3,
            description='Gap between doors',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    door_board_chanfer = FloatProperty(
            name='Board chanfer', min=0, default=0.02, precision=3,
            description='Board chanfer width',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    thickness = FloatProperty(
            name='Thickness', min=0.001, default=0.018, precision=3,
            description='Board thickness',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    y = FloatProperty(
            name='Depth', min=0.001, default=0.59, precision=3,
            description='Default cabinet depth',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    yw = FloatProperty(
            name='Wall Depth', min=0.001, default=0.39, precision=3,
            description='Default wall cabinet depth',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    z_mode = EnumProperty(
            name='Z mode',
            items=(
                ('1', 'Height Multiplier', 'Height are multiples of base height, default 6/6/16 modules'),
                ('2', 'Height Absolute', 'Height are absolute values'),
                ),
            default='1',
            update=update
            )
    module_size = FloatProperty(
            name='Height', min=0.001, default=0.127, precision=3,
            description='Base height for multiplier, default (12.7cm / 5")',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    modules_default = IntProperty(
            name="Default",
            description='#Modules for default cabinet (default = 6)',
            min=1, default=6,
            update=update
            )
    modules_wall = IntProperty(
            name="Wall",
            description='#Modules for wall cabinet (default = 6)',
            min=1, default=6,
            update=update)
    modules_full = IntProperty(
            name="Full",
            description='#Modules for full height cabinet (default = 16)',
            min=1, default=16,
            update=update
            )
    z_default = FloatProperty(
            name="Default",
            unit='LENGTH', subtype='DISTANCE',
            description='Height for default cabinet',
            min=0.001, default=0.762, precision=3,
            update=update
            )
    z_wall = FloatProperty(
            name="Wall",
            unit='LENGTH', subtype='DISTANCE',
            description='Height for wall cabinet',
            min=0.001, default=0.762, precision=3,
            update=update)
    z_full = FloatProperty(
            name="Full",
            unit='LENGTH', subtype='DISTANCE',
            description='Height for full height cabinet',
            min=0.001, default=2.032, precision=3,
            update=update
            )

    handle = EnumProperty(
            items=(
                ('2', "Rail", "Rail"),
                ('180', "Gap", "Gap on door with rail"),
                ('190', "Small rail on border", "Gap on door with small rail"),

                # Parametric T pipe
                ('140', "T pipe bar", "T pipe bar"),

                # Parametric handles
                ('150', "Circle handle", "Circle handle"),
                ('151', "Square handle", "Square handle"),
                ('152', "Spherical handle", "Spherical handle"),
                ('153', "Elliptic handle", "Elliptic handle"),
                ('154', "Mushroom like handle", "Mushroom like handle"),
                ('155', "Mushroom like handle 2", "Mushroom like handle 2"),

                # parametric bars
                ('120', "Modern tri bar", "Modern tri bar"),
                ('121', "Rounded square bar", "Rounded square bar"),
                ('122', "Old school rounded pipe bar", "Old school rounded pipe bar"),
                ('127', "Old school rounded pipe bar 2", "Old school rounded pipe bar 2"),
                ('123', "Modern Pipe bar", "Modern Pipe bar"),
                ('124', "Modern Square Bar", "Modern Square Bar"),
                ('125', "Old school elliptic bar", "Old school elliptic bar"),
                ('126', "Elliptic bar", "Elliptic bar"),

                ('0', "None", ""),

                ),
            name="Handle",
            default='140',
            description="Type of handle",
            update=update
            )
    handle_x = FloatProperty(
            name='Width', min=0.001,
            default=0.12, precision=3,
            description='Displacement in X relative position (limited to door size)',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handle_y = FloatProperty(
            name='Depth', min=0.001,
            default=0.035, precision=3,
            description='Depth',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handle_z = FloatProperty(
            name='Height', min=0.001,
            default=0.01, precision=3,
            description='Height',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handle_r = FloatProperty(
            name='Diameter', min=0.001,
            default=0.04, precision=3,
            description='Handle diameter',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handle_dx = FloatProperty(
            name='Border', min=0.001,
            default=0.05, precision=3,
            description='Distance from border',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handle_dz = FloatProperty(
            name='Altitude', min=0.001,
            default=0.15, precision=3,
            description='Displacement in Z relative position (limited to door size)',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    handle_fit = BoolProperty(
            name="Fit width / height",
            description="Fit door width / height",
            default=True,
            update=update
            )
    handle_space = FloatProperty(
            name='Space', min=0.0, default=0.05,
            precision=3,
            description='Handle space from border',
            update=update
            )

    baseboard = BoolProperty(
            name="Baseboard",
            description="Create a baseboard automatically",
            default=True,
            update=update
            )
    base_height = FloatProperty(
            name='Height', min=0.001,
            default=0.16, precision=3,
            description='Baseboard height',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    base_sink = FloatProperty(
            name='Sink', min=0,
            default=0.05, precision=3,
            description='Baseboard sink',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    counter = BoolProperty(
            name="Countertop",
            description="Create a countertop automatically (only default cabinet)",
            default=True,
            update=update
            )
    counter_z = FloatProperty(
            name='Height', min=0.001,
            default=0.02, precision=3,
            description='Countertop height',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    counter_y = FloatProperty(
            name='Extend', min=0.001,
            default=0.03, precision=3,
            description='Countertop extent',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    counter_chanfer = FloatProperty(
            name='Chanfer', min=0.0001,
            default=0.001, precision=3,
            description='Countertop chanfer',
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    cabinet_num = IntProperty(
            name='#Cabinets', min=1,
            default=1,
            description='Number total of cabinets in the Kitchen',
            update=update_manipulators
            )
    cabinets = CollectionProperty(type=archipack_kitchen_cabinet)
    expand = BoolProperty(
            options={'SKIP_SAVE'},
            name="expand",
            description="Expand kitchen properties",
            default=False
            )

    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )

    @property
    def height_default(self):
        if self.z_mode == '1':
            return self.module_size * self.modules_default
        else:
            return self.z_default

    @property
    def height_wall(self):
        if self.z_mode == '1':
            return self.module_size * self.modules_wall
        else:
            return self.z_wall

    @property
    def height_full(self):
        if self.z_mode == '1':
            return self.module_size * self.modules_full
        else:
            return self.z_full

    @property
    def altitude_wall(self):
        if self.z_mode == '1':
            return self.module_size * (self.modules_full - self.modules_wall)
        else:
            return self.z_full - self.z_wall

    def cabinet_depth(self, cab):
        if cab.location == 2:
            y = self.yw
        else:
            y = self.y
        return y

    def cabinet_height(self, cab):
        cab_location = cab.location
        if cab_location == 1:
            return self.height_default + cab.dz
        elif cab_location == 2:
            return self.height_wall + cab.dz
        else:
            return self.height_full + cab.dz

    def door_height(self, module):
        if module.modules == 0 or self.z_mode == '2':
            # user defined height / absolute z mode
            return module.z
        else:
            return self.module_size * module.modules

    def insert_part(self, context, where):
        self.manipulable_disable(context)
        self.auto_update = False
        ref = self.cabinets[where]
        c = self.cabinets.add()
        c.from_cabinet(ref)
        # move after current one
        self.cabinets.move(len(self.cabinets) - 1, where + 1)
        self.cabinet_num += 1
        self.setup_manipulators()
        self.auto_update = True

    def remove_part(self, context, where):
        self.manipulable_disable(context)
        self.auto_update = False
        self.cabinets.remove(where)
        self.cabinet_num -= 1
        self.setup_manipulators()
        self.auto_update = True

    def create_cabinet(self, cab, tM, verts, faces, matids, uvs):
        cab_depth = self.cabinet_depth(cab)
        th = self.thickness
        sx, sy, sz = cab.x, max(th + 0.001, cab_depth + cab.dy), self.cabinet_height(cab)
        door_y = self.door_y
        """     Cab                   Corner R                 Corner  L                Corner 45
                 x0 x1   x2 x3          x0 x1  x4 x5   x2 x3   x0 x1  x4 x5  x2 x3        x0 x1 x4 x2 x3
           y0    _____________         ___________________      __________________       ______________
           y1   | |____B____| |       | |_______B_________|    |________B_______| |     | |_____B______|
                | | |     | | |       | | |     | |     | |    | |    | |     | | |     | | |   |    | |
                |S|L|     |R|S|       |S|L|     |S|     |R|    |R|    |S|     |L|S|     |S|R|   |    | |
           y2   | |_|_____|_| |       | |_|_____|_|_____|_|    |_|____|_|_____|_| |     | |_|___|____| |
           y3   |_|         |_|       |_|____F____|                   |____F____|_|     |_|\    |    |L|
                                                                                            \   |    | |
                                                                                             \  |    | |
           y4                                                                                 \ |____| |
           y5                                                                                  \|__F_|_|
                                                                                              |____S___|

        """

        # also define material id
        cab_location = cab.location

        cab_type = cab.cab_type
        countertop_hole = cab.countertop_hole
        door_style = int(self.door_style)

        y = sy
        radius = 0
        shelve_radius = 0
        chanfer = 0

        if cab_type == 3:
            chanfer = 3
            # Corner 45
            x = max(sy, sx)
            radius = max(0, x - sy)
            shelve_radius = radius - th
            y = x
            x4 = x - sy

        elif cab_type == 7:
            chanfer = 3
            # Corner 45
            x = max(sx, sy)
            radius = max(0, x - sy)
            shelve_radius = radius
            y = x
            x4 = x - sy

        elif cab_type == 1:
            # Corner L
            s_left = door_y + cab.panel_right_width
            x = max(cab_depth + s_left, sx)
            x4 = x - cab_depth - s_left

        elif cab_type == 2:
            # Corner R
            s_right = door_y + cab.panel_left_width
            x = max(cab_depth + s_right, sx)
            x4 = cab_depth + s_right
        else:
            x = sx
            x4 = 0

        # cabinet box size
        x0 = 0
        x1 = th
        x2 = x - th
        x3 = x

        y0 = 0
        y1 = -th
        y2 = -sy
        y3 = y2 - door_y
        y4 = th - y
        y5 = -y

        z0 = 0

        # wihtout left side
        if cab_type in {4, 6, 7}:
            x1 = x0

        # wihtout right side
        if cab_type in {5, 6}:
            x2 = x3
        elif cab_type == 7:
            y4 = y5

        # wall cabinet
        if cab_location == 2:
            # range hood
            if len(cab.modules) > 0:
                module = cab.modules[0]
                if module.type == '55':
                    if module.modules == 0 or self.z_mode == '2':
                        # user defined height / absolute z mode
                        zd = module.z
                    else:
                        zd = self.module_size * module.modules
                    z0 += zd

        z1 = z0 + th
        # z top inside variable for cabinets hole on top
        z2 = sz - th
        z3 = sz

        if countertop_hole:
            # cabinet with hole on top
            z2 = z3

        # Side boards on corners to prevent door opening conflict
        if cab_type == 1:
            # Corner L
            left_chanfer = 0
            if door_style > 10:
                left_chanfer = 2
            make_box(tM, x4, x3, y2, y3, z0, z3, self.door_chanfer,
                cab_location, left_chanfer, verts, faces, matids, uvs)

        elif cab_type == 2:
            # Corner R
            right_chanfer = 0
            if door_style > 10:
                right_chanfer = 2
            make_box(tM, x4, x0, y2, y3, z3, z0, self.door_chanfer,
                cab_location, right_chanfer, verts, faces, matids, uvs)

        # ------------
        # manipulators
        # ------------
        m_dist = (1 + (cab.location % 2)) * 0.25
        # gl location of cab number
        m_pos = tM * Vector((0.5 * x, -0.5 * y, sz))
        cab.manipulators[0].set_pts([
                m_pos,
                m_pos,
                (-3, 0, 0)
                ])

        # cabinet width
        m_left = tM * Vector((0, 0, z0))
        m_right = tM * Vector((x, 0, z0))
        cab.manipulators[1].set_pts([
                m_left,
                m_right,
                (-m_dist, 0, 0)
                ])
                
        # Dimension points of cab
        self.add_dimension_point(cab.uid, m_left)
        self.add_dimension_point(cab.uid + 1, m_right)
        
        if cab_type == 1:
            # Corner L side
            m_left = tM * Vector((x, - door_y - y, z0))

        elif cab_type == 2:
            # Corner R side
            m_right = tM * Vector((0, - door_y - y, z0))

        elif cab_type in {3, 7}:
            # Corner 45 side
            m_left = tM * Vector((x, -x, z0))

        # Corner dumb side
        cab.manipulators[2].set_pts([
                m_left,
                m_right,
                (m_dist, 0, 0)
                ])

        if cab_type != 8:
            # bottom
            make_box(tM, x0, x3, y0, y5, z0, z1, radius,
                cab_location, chanfer, verts, faces, matids, uvs)
            matids[-1] = mat_inside

            # back
            make_box(tM, x0, x3, y0, y1, z1, z3, 0,
                cab_location, False, verts, faces, matids, uvs)
            matids[-3] = mat_inside

            # L side
            if cab_type not in {4, 6, 7}:
                make_box(tM, x0, x1, y1, y2, z1, z2, 0,
                    cab_location, False, verts, faces, matids, uvs)
                matids[-4] = mat_inside

            # R side (is back part 2 of 45 degree)
            if cab_type not in {5, 6}:
                make_box(tM, x2, x3, y1, y5, z1, z2, 0,
                    cab_location, False, verts, faces, matids, uvs)
                matids[-6] = mat_inside

            # top
            if not countertop_hole:
                make_box(tM, x0, x3, y1, y5, z2, z3, radius,
                    cab_location, chanfer, verts, faces, matids, uvs)
                matids[-2] = mat_inside

            # corner 45 "Right" side of in front
            if cab_type == 3:
                make_box(tM, x4, x2, y4, y5, z1, z2, 0,
                    cab_location, False, verts, faces, matids, uvs)
                matids[-5] = mat_inside

            # -----------------
            # side boards
            # -----------------
            if cab.panel_left and cab_type != 2:
                left_chanfer = 0
                if door_style > 10:
                    left_chanfer = 1
                make_box(tM, x0 - cab.board_left, x0, y0, y3, z0, z3, self.door_chanfer,
                    cab_location, left_chanfer, verts, faces, matids, uvs)

            if cab.panel_right and cab_type != 1:
                right_chanfer = 0
                # corner 45 "Right" board in front
                if cab_type in {3, 7}:
                    if door_style > 10:
                        right_chanfer = 2
                    make_box(tM, x4 - door_y, x3, y5, y5 - cab.board_right, z0, z3, self.door_chanfer,
                        cab_location, right_chanfer, verts, faces, matids, uvs)
                else:
                    if door_style > 10:
                        right_chanfer = 1
                    make_box(tM, x3, x3 + cab.board_right, y0, y3, z0, z3, self.door_chanfer,
                        cab_location, right_chanfer, verts, faces, matids, uvs)

        # -----------------
        # sink / cook top
        # -----------------
        if cab_location == 1:
            counter_type = int(cab.counter)
            hx = cab.counter_x
            hy = cab.counter_y
            tM2 = tM.copy()
            cx = 0.5 * cab.x
            if cab_type == 2:
                cx += self.y
            cy = -0.5 * max(0, self.y + cab.dy)
            r = 0.02

            # corner 45 rotate 45 deg
            if cab_type in {3, 7}:
                x = max(sy, sx)
                radius = max(0, x - sy)
                tM2 = tM * Matrix([
                    [sq2, sq2, 0, 0.5 * radius],
                    [-sq2, sq2, 0, 0.5 * radius - x],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]
                    ])
                cx = 0
                cy = 0.5 * hy + 2 * self.thickness

            if counter_type == 11:
                sink(tM2, cx, cy, sz + self.counter_z, r, hx, hy, verts, faces, matids, uvs)
            elif counter_type == 2:
                cook_top(tM2, cx, cy, sz + self.counter_z, r, hx, hy, verts, faces, matids, uvs)

        # -----------------
        # shelves
        # -----------------
        tM2 = tM.copy()
        if cab_type != 8:
            if cab_type in {3, 7}:
                a = pi / 4
                ca = cos(a)
                sa = sin(a)
                sx = (2 ** 0.5) * radius - 2 * sq2 * door_y
                dx = -0.5 * sx
                dy = sy + (1 - sq2) * door_y
                tM2 = tM * Matrix([
                    [sq2, sq2, 0, 0.5 * radius + dx * ca + dy * sa],
                    [-sq2, sq2, 0, 0.5 * radius - x + dx * -sa + dy * ca],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]
                    ])
            else:
                y4 = y2

            for i, module in enumerate(cab.modules):

                module_type = int(module.type)

                # Door height
                zd = self.door_height(module)
                zd = min(zd, sz - z0)

                if zd == 0:
                    continue

                # Drawer and raw face board dosent allow shelves
                if module_type in {1, 6}:
                    z0 += zd
                    continue

                if module_type == 54:
                    oven_cab(tM2, sx, sy, zd, z0, th, door_y, verts, faces, matids, uvs)
                elif module_type == 55:
                    rangehood(tM2, sx, sy, zd, z0 - zd, th, door_y, verts, faces, matids, uvs)
                    zd = 0
                elif module_type == 59:
                    dishwasher_cab(tM2, sx, sy, zd, z0, th, door_y, verts, faces, matids, uvs)
                else:

                    # calculate separation
                    space = zd / (module.shelves + 1)
                    z1 = z0

                    # Ground not needed on 1 element
                    if i == 0:
                        start = 1
                        z1 += space
                    else:
                        start = 0

                    # shelves
                    for x in range(start, module.shelves + 1):
                        make_box(tM, x1, x2, y1, y4, z1, z1 + th, shelve_radius,
                            mat_inside, chanfer, verts, faces, matids, uvs)

                        z1 += space

                z0 += zd

    def update_cabinets(self, context, verts, faces, matids, uvs):

        # Vector for axis location
        loc = Vector((0, 0, 0))
        a = 0
        ca = 1
        sa = 0
        door_y = self.door_y

        for cab in self.cabinets:

            # side boards size
            board_left = cab.board_left
            board_right = cab.board_right
            cab_depth = self.cabinet_depth(cab)

            z = cab.pz
            cab_type = cab.cab_type
            cab_location = cab.location

            if self.baseboard:
                z += self.base_height

            # Wall
            if cab_location == 2:
                z += self.altitude_wall - cab.dz

            if cab.reset_location:
                loc = Vector((0, 0, 0))
                a = 0

            # rotate R cab
            if cab_type == 2:
                dx = cab_depth + door_y
                a += pi / 2
                # use last rotation and move basis location
                loc += Vector((dx * ca, dx * -sa, 0))

            else:
                rot_type = int(cab.rotate)
                if rot_type == RotationType_R90CW:
                    a += pi / 2
                elif rot_type == RotationType_R90CCW:
                    a -= pi / 2
                elif rot_type == RotationType_R180:
                    a += pi
                elif rot_type == RotationType_User:
                    a += cab.angle

            # add translation between last and current cab
            ca = cos(a)
            sa = sin(a)
            loc += Vector((
                (board_left + cab.px) * ca + cab.py * sa,
                (board_left + cab.px) * -sa + cab.py * ca,
                z))

            # matrix at cabinet border with z at bottom
            tM = Matrix([
                [ca, sa, 0, loc.x],
                [-sa, ca, 0, loc.y],
                [0, 0, 1, loc.z],
                [0, 0, 0, 1]
                ])

            self.create_cabinet(cab, tM, verts, faces, matids, uvs)

            if cab.lock_p:
                dx = 0
                dy = 0
                dz = cab.pz - z
            else:
                dx = -cab.px
                dy = -cab.py
                dz = -z

            # add translation for current cab to border
            if cab_type == 1:
                # Corner L
                a += pi / 2
                s_right = door_y + cab.panel_right_width
                dx += max(cab_depth + s_right, cab.x)
                loc = tM * Vector((dx, dy - cab_depth - door_y, dz))
            elif cab_type == 2:
                # Corner R
                s_left = door_y + cab.panel_left_width
                dx += max(cab_depth + s_left, cab.x) + board_right
                loc = tM * Vector((dx, 0, dz))
            elif cab_type in {3, 7}:
                # Corner 45
                a += pi / 2
                dx += max(cab_depth + max(0, cab.dy), cab.x)
                dy -= max(cab_depth + max(0, cab.dy), cab.x) + board_right
                loc = tM * Vector((dx, dy, dz))
            else:
                # add translation for current cab to border
                loc = tM * Vector((dx + cab.x + board_right, dy, dz))

    def create_baseboard(self, cab, dz, tM, verts, faces, matids, uvs):
        cab_depth = self.cabinet_depth(cab)
        sx, sy, sz = cab.x, self.y + cab.dy, self.base_height + dz + cab.pz
        th = self.thickness
        door_y = self.door_y
        if sz <= 0:
            return

        board_left = cab.board_left
        board_right = cab.board_right
        cab_type = cab.cab_type
        base_sink = self.base_sink + cab.base_front

        y = sy

        if cab_type in {3, 7}:
            # Corner 45
            x = max(cab_depth, sx)
            y = x
            x4 = max(
                base_sink,
                min(x - th, x - cab_depth - cab.dy + base_sink)
                )

        elif cab_type == 1:
            # Corner L
            s_left = door_y + cab.panel_right_width
            x = max(cab_depth + s_left, sx)
            x4 = x - cab_depth + base_sink

        elif cab_type == 2:
            # Corner R
            s_right = door_y + cab.panel_left_width
            x = max(cab_depth + s_right, sx)
            x4 = cab_depth - th - base_sink
        else:
            x = sx
            x4 = 0

        # cabinet box size
        x0 = 0
        x1 = x0 + th
        x2 = x - th
        x3 = x
        x5 = x4 + th
        y0 = 0
        y2 = min(-th, base_sink - sy)
        y1 = y2 + th
        y5 = -y

        z0 = 0
        z1 = sz

        # side sink
        cab_sink = min(0.5 * x - th - 0.0001, cab.base_sink)

        """     Cab                   Corner R                 Corner  L                Corner 45
                 x0 x1   x2 x3        x0 x1  x4 x5   x2 x3   x0 x1  x4 x5  x2 x3      x0 x1  x4  x2 x3
           y0    _____________         _________________      ________________        ______________
           y1   | |____B____| |       |_______B_________|    |________B_______|      | |_____B______|
                | | |     | | |       | |     | |     | |    | |    | |     | |      | | |   |    | |
                |S|L|     |R|S|       |L|     |S|     |R|    |R|    |S|     |L|      |S|R|   |    | |
           y2   | |_|_____|_| |       |_|_____|_|_____|_|    |_|____|_|_____|_|      | |_|___|____| |
           y3   |_|         |_|       |____D____|                   |____D____|      |_|\ \  |    |L|
                                      |____S____|                   |____S____|          \ \ |    | |
                                                                                          \ \|    | |
           y4                                                                              \ |____| |
           y5                                                                               \|__F_|_|
                                                                                           |____S___|

        """

        if cab_type in {3, 7}:

            sx += self.cabinet_depth(cab)

            r = base_sink * tan22_5
            x1 = x0 + r
            x4 = x + y2
            x5 = x4 + th
            y4 = x1 - x
            x0 -= board_left
            y5 -= board_right
            # offset sides
            if cab.base_left:
                # left side board
                x0 += min(board_left + r - 0.001, cab.base_sink)
                make_box(tM, x0, x0 + th, y0, y1, z0, z1, 0, mat_baseboard, False, verts, faces, matids, uvs)

            if cab.base_right:
                # right side board
                y5 += min(board_right + r - 0.001, cab.base_sink)
                make_box(tM, x5, x, y5 + th, y5, z0, z1, 0, mat_baseboard, False, verts, faces, matids, uvs)

            f = len(verts)
            # left board
            rt = th * tan22_5

            # center board
            verts.extend([tM * Vector(v) for v in [
                (x0, y1, z0), (x0, y2, z0), (x0, y2, z1), (x0, y1, z1),
                (x1 + rt, y1, z0), (x1, y2, z0), (x1, y2, z1), (x1 + rt, y1, z1),
                (x5, y4 + rt, z0), (x4, y4, z0), (x4, y4, z1), (x5, y4 + rt, z1),
                (x5, y5, z0), (x4, y5, z0), (x4, y5, z1), (x5, y5, z1)
                ]])
            s = 4
            faces.extend([(
                f + j * s + i + 1,
                f + j * s + i,
                f + (j + 1) * s + i,
                f + (j + 1) * s + i + 1
                ) for i in range(s - 1) for j in range(3)])
            faces.extend([(
                f + j * s,
                f + (j + 1) * s - 1,
                f + (j + 2) * s - 1,
                f + (j + 1) * s
                ) for j in range(3)])
            faces.append(tuple([f + i for i in range(s)]))
            faces.append(tuple([f + 4 * s - 1 - i for i in range(s)]))
            matids.extend([mat_baseboard for i in range(s * 3 + 2)])

        else:

            # offset sides
            if cab_type == 1:
                # Left corner
                y5 -= door_y
                x2 = min(x4, x2)
                x3 = x2 + th
                x0 -= board_left
                x1 -= board_left

            elif cab_type == 2:
                # Right corner
                y5 -= door_y
                x0 = max(x4, x0)
                x1 = x0 + th
                x2 += board_right
                x3 += board_right

            elif cab_type in {0, 4, 5, 6}:
                # straight cabinets
                x0 -= board_left
                x1 -= board_left
                x2 += board_right
                x3 += board_right

            if cab.base_left and cab_type != 2:
                x0 += cab_sink
                x1 += cab_sink

            if cab.base_right and cab_type != 1:
                x2 -= cab_sink
                x3 -= cab_sink

            # external faces
            make_box(tM, x0, x3, y1, y2, z0, z1, 0, mat_baseboard, False, verts, faces, matids, uvs)

            # left side
            if cab.base_left and cab_type != 2:
                make_box(tM, x0, x1, y0, y1, z0, z1, 0, mat_baseboard, False, verts, faces, matids, uvs)

            # right side
            if cab.base_right and cab_type != 1:
                make_box(tM, x2, x3, y0, y1, z0, z1, 0, mat_baseboard, False, verts, faces, matids, uvs)

            # Corners
            if cab_type in {1, 2}:
                make_box(tM, x4, x5, y2, y5, z0, z1, 0, mat_baseboard, False, verts, faces, matids, uvs)

    def update_baseboard(self, verts, faces, matids, uvs):
        loc = Vector((0, 0, 0))
        a = 0
        ca = 1
        sa = 0
        door_y = self.door_y

        # accumulate pz for lock mode
        z_accum = 0

        for cab in self.cabinets:
            cab_depth = self.cabinet_depth(cab)
            # side boards size
            board_left = cab.board_left
            board_right = cab.board_right

            if cab.reset_location:
                loc = Vector((0, 0, 0))
                a = 0
                z_accum = 0

            cab_location = cab.location
            cab_type = cab.cab_type

            # rotate R cab
            if cab_type == 2:
                dx = cab_depth + door_y
                a += pi / 2
                loc += Vector((dx * ca, dx * -sa, 0))
            else:

                rot_type = int(cab.rotate)

                if rot_type == RotationType_R90CW:
                    a += pi / 2
                elif rot_type == RotationType_R90CCW:
                    a -= pi / 2
                elif rot_type == RotationType_R180:
                    a += pi
                elif rot_type == RotationType_User:
                    a += cab.angle

            a = a % (2 * pi)

            # add translation from last to current cab
            ca = cos(a)
            sa = sin(a)
            loc += Vector((
                (board_left + cab.px) * ca + cab.py * sa,
                (board_left + cab.px) * -sa + cab.py * ca,
                0))

            # matrix at cabinet axis with z at bottom
            tM = Matrix([
                [ca, sa, 0, loc.x],
                [-sa, ca, 0, loc.y],
                [0, 0, 1, loc.z],
                [0, 0, 0, 1]
                ])

            if cab_location != 2 and cab.baseboard and cab_type != 8:
                self.create_baseboard(cab, z_accum, tM, verts, faces, matids, uvs)

            if cab.lock_p:
                dx = 0
                dy = 0
                dz = 0
                z_accum += cab.pz
            else:
                dx = -cab.px
                dy = -cab.py
                dz = 0

            if cab_type == 1:
                a += pi / 2
                # add translation for current cab to border
                s_right = door_y + cab.panel_right_width
                dx += max(cab_depth + s_right, cab.x)
                loc = tM * Vector((dx, dy - cab_depth - door_y, dz))
            elif cab_type == 2:
                s_left = door_y + cab.panel_left_width
                dx += max(cab_depth + s_left, cab.x) + board_right
                loc = tM * Vector((dx, 0, dz))
            elif cab_type in {3, 7}:
                # Corner 45
                a += pi / 2
                dx += max(cab_depth + max(0, cab.dy), cab.x)
                dy -= max(cab_depth + max(0, cab.dy), cab.x) + board_right
                loc = tM * Vector((dx, dy, dz))
            else:
                # add translation for current cab to border
                loc = tM * Vector((dx + cab.x + board_right, dy, dz))

    def create_counter(self, cab, tM, verts, faces, matids, uvs, start_new, need_endsection, close, last_countertop):
        """
          start new: generate a section at start point
          with L and R corners, the corner section is added
          with 45 corners, 2 sections are added
          need_endsection: generate faces and vertices for this section
                with regular cabs generate faces from start to end
                with L and R corner, generate faces from start to edge
                with 45 corner, generate faces from start to corner 1
        """
        sx, sy, sz = cab.x, max(0, self.y + cab.dy), self.counter_z
        over = self.counter_y
        cab_depth = self.cabinet_depth(cab)
        chanfer = self.counter_chanfer
        board_left = cab.board_left
        board_right = cab.board_right
        door_y = self.door_y

        chanfer_225 = 0
        over_255 = 0
        z0 = cab.dz + cab.pz
        z2 = z0 + sz
        z1 = z2 - chanfer

        # if corner the size is less
        ts = -board_left

        countertop_hole = cab.countertop_hole
        cab_type = cab.cab_type

        # use regular faces on last section
        last_section_faces = True
        last_section_verts = True

        y0 = 0
        y2 = -(sy + over)
        y1 = y2 + chanfer

        # first section
        if start_new:

            # corner 45 add special section when left_board is not there
            if board_left == 0 and cab_type in {3, 7}:
                chanfer_225 = tan22_5 * chanfer
                over_255 = tan22_5 * over

            f = len(verts)

            verts.extend([tM * Vector(v) for v in [
                (ts, y0, z0), (ts, y2 - over_255, z0), (ts, y2 - over_255, z1),
                (ts, y1 - over_255 + chanfer_225, z2), (ts, y0, z2)
                ]])

            # start face
            matids.append(mat_counter)
            faces.extend([tuple([f + i for i in v]) for v in [
                (0, 1, 2, 3, 4)
                ]])

        f = len(verts) - 5

        # 45 corners first section
        # add special corner section after board_left when present
        # and when cab is not starting one
        if cab_type in {3, 7} and (not start_new or board_left > 0):

            chanfer_225 = tan22_5 * chanfer
            over_255 = tan22_5 * over
            ts += board_left

            verts.extend([tM * Vector(v) for v in [
                (ts, y0, z0), (ts - over_255, y2, z0), (ts - over_255, y2, z1),
                (ts - over_255 + chanfer_225, y1, z2), (ts, y0, z2)
                ]])

            if not last_countertop:
                matids.extend([mat_counter for i in range(5)])
                faces.extend([tuple([f + i for i in v]) for v in [
                    # bottom,      back,        up,           chanfer,      front
                    (0, 4, 9, 5), (4, 3, 8, 9), (3, 2, 7, 8), (2, 1, 6, 7), (1, 0, 5, 6)
                    ]])

            f = len(verts) - 5

        # R corner section before hole
        # add faces here
        if cab_type == 2:
            f = len(verts) - 5
            # Corner R
            x0 = -cab.px
            x2 = -y1 - cab.px
            x1 = x2 - self.counter_chanfer
            # add corner section
            verts.extend([tM * Vector(v) for v in [
                (x0, y0, z0), (x2, y2, z0), (x2, y2, z1), (x1, y1, z2), (x0, y0, z2)
                ]])

            if not last_countertop:
                matids.extend([mat_counter for i in range(5)])
                faces.extend([tuple([f + i for i in v]) for v in [
                    # bottom,      back,        up,           chanfer,      front
                    (0, 4, 9, 5), (4, 3, 8, 9), (3, 2, 7, 8), (2, 1, 6, 7), (1, 0, 5, 6)
                    ]])

        # at this point f is first vert of last section
        if countertop_hole:

            hx = cab.counter_x
            hy = cab.counter_y
            cx = 0.5 * sx
            if cab_type == 2:
                cx += self.y
            cy = 0.5 * sy
            r = 0.02
            tM2 = tM.copy()

            # corner 45 rotate 45 deg
            if cab_type in {3, 7}:

                x = max(sy, sx)
                radius = max(0, x - sy)
                tM2 = tM * Matrix([
                    [sq2, sq2, 0, 0.5 * radius],
                    [-sq2, sq2, 0, 0.5 * radius - x],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]
                    ])
                cx = 0
                cy = -0.5 * hy - 2 * self.thickness

            # add verts for hole
            verts.extend([tM2 * Vector(v) for v in chanfer_square(cx, -cy, z0, r, hx, hy)])
            verts.extend([tM2 * Vector(v) for v in chanfer_square(cx, -cy, z2, r, hx, hy)])

            # faces for regular and L R corners
            if cab_type not in {3, 7}:
                matids.extend([mat_counter for i in range(27)])
                faces.extend([tuple([f + i for i in v]) for v in [
                    (17, 33, 32, 16), (10, 26, 25, 9), (18, 34, 33, 17),
                    (11, 27, 26, 10), (19, 35, 34, 18), (12, 28, 27, 11),
                    (20, 36, 35, 19), (13, 29, 28, 12), (6, 22, 21, 5),
                    (5, 21, 36, 20), (14, 30, 29, 13), (7, 23, 22, 6),
                    (15, 31, 30, 14), (8, 24, 23, 7), (16, 32, 31, 15),
                    (9, 25, 24, 8),
                    (12, 11, 10, 9, 8, 7, 6, 5, 1, 0),
                    (13, 37, 38, 20, 19, 18, 17, 16, 15, 14),
                    (13, 12, 0, 37), (20, 38, 1, 5), (0, 4, 41, 37),
                    (4, 3, 21, 22, 23, 24, 25, 26, 27, 28),
                    (29, 30, 31, 32, 33, 34, 35, 36, 40, 41),
                    (36, 21, 3, 40),
                    (29, 41, 4, 28), (3, 2, 39, 40), (2, 1, 38, 39)
                    ]])

        # Corner L section after hole
        if cab_type == 1:
            # Corner L
            s_left = door_y + cab.panel_right_width
            x = max(cab_depth + s_left, sx)
            x0 = x - cab.px
            x2 = max(
                over,
                min(x, x - cab_depth - over)
                )
            x1 = x2 + chanfer

            # when need_endsection
            # add corner here and use last
            # when not, use last as corner
            # with hole faces are allready there
            # without need faces
            if need_endsection:

                f = len(verts) - 5

                verts.extend([tM * Vector(v) for v in [
                    (x0, y0, z0), (x2, y2, z0), (x2, y2, z1), (x1, y1, z2), (x0, y0, z2)
                    ]])

                if not countertop_hole:
                    matids.extend([mat_counter for i in range(5)])
                    faces.extend([tuple([f + i for i in v]) for v in [
                        # bottom,      back,        up,           chanfer,      front
                        (0, 4, 9, 5), (4, 3, 8, 9), (3, 2, 7, 8), (2, 1, 6, 7), (1, 0, 5, 6)
                        ]])

                y0 = -(self.y + door_y)
                y1 = y0
                y2 = y0

        # Corner R section required if there is a hole
        # or need_endsection
        elif cab_type == 2:

            last_section_faces = False
            last_section_verts = False

            if countertop_hole or need_endsection:
                last_section_faces = True
                last_section_verts = True

            s_right = door_y + cab.panel_left_width
            x0 = max(sx, cab_depth + s_right) + board_right
            x1 = x0
            x2 = x0
            y0 = 0
            y2 = -(sy + over)
            y1 = y2 + self.counter_chanfer

        # Corner 45
        # corner section faces and verts
        elif cab_type in {3, 7}:
            # Corner 45
            chanfer_225 = tan22_5 * chanfer
            over_255 = tan22_5 * over
            x = max(cab_depth, sx)
            x0 = x
            x2 = max(
                over,
                min(x, x - cab_depth - cab.py - cab.dy - over)
                )
            x1 = x2 + chanfer
            y1 = -x
            # corner
            verts.extend([tM * Vector(v) for v in [
                (x0, 0, z0), (x0, 0, z2)
                ]])

            last_section_faces = False

            # dont enlarge when element is last one, unless board_right > 0
            if need_endsection and board_right == 0:
                # single end section dosent use additional faces
                y0 = y1
                y2 = y1
                x2 -= over_255
                x1 += chanfer_225 - over_255
            else:
                # this section goes over next element
                y0 = y1
                y2 = y1 - over_255
                y1 += chanfer_225 - over_255
                # if endsection or board right > 0
                # we need a corner section before end one
                if need_endsection or board_right > 0:
                    verts.extend([tM * Vector(v) for v in [
                        (x0, y0, z0), (x2, y2, z0), (x2, y2, z1), (x1, y1, z2), (x0, y0, z2)
                        ]])
                    last_section_faces = True
                    y0 -= board_right
                    y1 = y0
                    y2 = y0

            # faces of center section
            if countertop_hole:
                matids.extend([mat_counter for i in range(28)])
                faces.extend([tuple([f + i for i in v]) for v in [
                    (17, 33, 32, 16), (10, 26, 25, 9),
                    (18, 34, 33, 17), (11, 27, 26, 10), (19, 35, 34, 18),
                    (12, 28, 27, 11), (20, 36, 35, 19), (13, 29, 28, 12),
                    (6, 22, 21, 5), (5, 21, 36, 20), (14, 30, 29, 13),
                    (7, 23, 22, 6), (15, 31, 30, 14), (8, 24, 23, 7),
                    (16, 32, 31, 15), (9, 25, 24, 8), (12, 11, 10, 9, 8, 7, 6, 5, 1, 0),
                    (4, 3, 21, 22, 23, 24, 25, 26, 27, 28), (0, 4, 38, 37),
                    (37, 38, 43, 39), (2, 1, 40, 41), (3, 2, 41, 42),
                    (36, 21, 3, 42), (20, 40, 1, 5), (32, 33, 34, 35, 36, 42, 43, 29, 30, 31),
                    (29, 43, 38, 4, 28), (13, 12, 0, 37, 39), (17, 16, 15, 14, 13, 39, 40, 20, 19, 18)
                    ]])
            else:
                matids.extend([mat_counter for i in range(6)])
                faces.extend([tuple([f + i for i in v]) for v in [
                    (1, 0, 5, 7, 8), (3, 10, 11, 6, 4),
                    (0, 4, 6, 5), (5, 6, 11, 7), (2, 1, 8, 9),
                    (3, 2, 9, 10)
                    ]])
        else:

            if not need_endsection:
                last_section_faces = False
                last_section_verts = False

            elif countertop_hole:
                last_section_faces = False

            tx = sx + board_right
            x0, x1, x2 = tx, tx, tx

        """
        xy2 xy1               xy0
           __________________  z2
          /                  | z1
         |___________________| z0

        """
        if last_section_verts:
            # add last section verts
            verts.extend([tM * Vector(v) for v in [
                (x0, y0, z0), (x2, y2, z0), (x2, y2, z1), (x1, y1, z2), (x0, y0, z2)
                ]])

        if last_section_faces:
            # last section faces
            f = len(verts) - 10
            matids.extend([mat_counter for i in range(5)])
            faces.extend([tuple([f + i for i in v]) for v in [
                # bottom,      back,        up,           chanfer,      front
                (0, 4, 9, 5), (4, 3, 8, 9), (3, 2, 7, 8), (2, 1, 6, 7), (1, 0, 5, 6)
                ]])

        # close
        if close:
            # end side faces
            f = len(verts) - 5
            matids.append(mat_counter)
            faces.extend([tuple([f + i for i in v]) for v in [
                (4, 3, 2, 1, 0)
                ]])

    def update_counter(self, verts, faces, matids, uvs):
        loc = Vector((0, 0, 0))
        a = 0
        ca = 1
        sa = 0
        z = self.height_default
        door_y = self.door_y

        # start new countertop
        start_new = True

        # next one should not add faces on first section
        # when section is countertop with hole
        last_countertop = False

        if self.baseboard:
            z += self.base_height

        n_cabs = self.cabinet_num
        next_type = -1
        for i, cab in enumerate(self.cabinets):

            cab_depth = self.cabinet_depth(cab)
            # side boards size
            board_left = cab.board_left
            board_right = cab.board_right

            if cab.reset_location:
                loc = Vector((0, 0, 0))
                a = 0

            cab_type = cab.cab_type
            cab_location = cab.location

            # rotate R cab
            if cab_type == 2:
                dx = cab_depth + door_y
                a += pi / 2
                loc += Vector((dx * ca, dx * -sa, 0))
            else:

                rot_type = int(cab.rotate)

                if rot_type == RotationType_R90CW:
                    a += pi / 2
                elif rot_type == RotationType_R90CCW:
                    a -= pi / 2
                elif rot_type == RotationType_R180:
                    a += pi
                elif rot_type == RotationType_User:
                    a += cab.angle

            a = a % (2 * pi)

            # add translation from last to current cab
            ca = cos(a)
            sa = sin(a)
            loc += Vector((
                (board_left + cab.px) * ca + cab.py * sa,
                (board_left + cab.px) * -sa + cab.py * ca,
                z))

            # matrix at cabinet axis with z at bottom
            tM = Matrix([
                [ca, sa, 0, loc.x],
                [-sa, ca, 0, loc.y],
                [0, 0, 1, loc.z],
                [0, 0, 0, 1]
                ])

            # next cab restart
            restart = False

            # current cab need end section
            need_endsection = False

            # add ending faces
            close = True

            # Compare next to current to know if
            # need start / end section
            # need full restart
            # By default next one use current last vertex section
            # we may remove last section unless there is a hole
            # even with hole, corner 45 and L corner use own section

            if cab_location == 1 and cab.counter != "0":
                if i + 1 < n_cabs:
                    next_cab = self.cabinets[i + 1]
                    next_type = next_cab.cab_type

                    # restart is needed
                    # when cab not on same altitude
                    # seek for diff between current and next
                    # change when any != 0 and not lock_p
                    delta_z = next_cab.dz + next_cab.pz - cab.dz

                    if not cab.lock_p:
                        delta_z -= cab.pz

                    if (next_cab.location != 1 or
                            next_cab.px != 0 or
                            next_cab.py != 0 or
                            next_cab.dy != cab.dy or
                            abs(delta_z) > 0.0001 or
                            next_cab.counter == "0" or
                            next_cab.reset_location or
                            next_cab.rotate != "0"
                            ):
                        restart = True
                        need_endsection = True
                    else:
                        # we may remove end section unless there is a hole
                        # even with hole, corner 45 and R corner use own section
                        close = False
                        if cab.countertop_hole:
                            if next_type not in {2, 3, 7}:
                                need_endsection = True
                else:
                    restart = True
                    need_endsection = True

                # only build when section is needed and for cabs adding sections
                if start_new or need_endsection or cab_type in {1, 2, 3, 7}:
                    self.create_counter(cab, tM, verts, faces, matids, uvs,
                        start_new, need_endsection, close, last_countertop)

                # next one should not add regular faces on section before corner
                last_countertop = cab.countertop_hole and next_type in {2, 3, 7}
                start_new = restart

            else:
                last_countertop = False
                start_new = True

            if cab.lock_p:
                dx = 0
                dy = 0
                dz = cab.pz - z
            else:
                dx = -cab.px
                dy = -cab.py
                dz = -z

            if cab_type == 1:
                a += pi / 2
                # add translation for current cab to border
                s_right = door_y + cab.panel_right_width
                dx += max(cab_depth + s_right, cab.x)
                loc = tM * Vector((dx, dy - cab_depth - door_y, dz))
            elif cab_type == 2:
                s_left = door_y + cab.panel_left_width
                dx += max(cab_depth + s_left, cab.x) + board_right
                loc = tM * Vector((dx, 0, dz))
            elif cab_type in {3, 7}:
                # Corner 45
                a += pi / 2
                dx += max(cab_depth + max(0, cab.dy), cab.x)
                dy -= max(cab_depth + max(0, cab.dy), cab.x) + board_right
                loc = tM * Vector((dx, dy, dz))
            else:
                # add translation for current cab to border
                loc = tM * Vector((dx + cab.x + board_right, dy, dz))

    def get_childs_modules(self, o):
        return [child for child in o.children if archipack_kitchen_module.filter(child)]

    def remove_modules(self, context, childs, to_remove):
        for child in childs:
            if to_remove < 1:
                return
            to_remove -= 1
            self.delete_object(context, child)

    def update_modules(self, context, o):

        # real childs
        childs = self.get_childs_modules(o)

        # wanted childs
        w_childs = 0
        for cab in self.cabinets:

            z0 = 0
            zmax = self.cabinet_height(cab)

            cab.update_parts()

            if cab.cab_type == 8:
                continue

            for module in cab.modules:

                n_panels = module.n_panels

                # Door height
                zd = self.door_height(module)

                zd = min(zd, zmax - z0)
                if zd == 0:
                    continue

                # range hood, no door
                if module.type == '55':
                    n_panels = 0

                # No door, accumulate spacing
                w_childs += n_panels
                z0 += zd

        # real childs
        n_childs = len(childs)

        # remove child
        if n_childs > w_childs:
            self.remove_modules(context, childs, n_childs - w_childs)

        if w_childs == 0:
            return

        childs = self.get_childs_modules(o)
        n_childs = len(childs)

        loc = Vector((0, 0, 0))
        tM = Matrix()
        a = 0
        panel = 0
        ca = 1
        sa = 0
        door_y = self.door_y

        for cab in self.cabinets:

            cab_depth = self.cabinet_depth(cab)

            # side boards size
            board_left = cab.board_left
            board_right = cab.board_right
            cab_type = cab.cab_type
            cab_location = cab.location

            z = cab.pz

            if self.baseboard:
                z += self.base_height

            sx = cab.x
            sy = cab_depth + cab.dy
            dx = 0

            if cab_type == 1:
                # Corner L
                s_right = door_y + cab.panel_right_width
                sx = cab.x - min(cab_depth + s_right, cab.x)
            elif cab_type == 2:
                # Corner R
                s_left = door_y + cab.panel_left_width
                sx = cab.x - min(cab_depth + s_left, cab.x)
            elif cab_type in {3, 7}:
                # Corner 45
                x = max(sy, sx)
                radius = max(0, x - sy)
                sx = (2 ** 0.5) * radius

            sy += door_y

            # Wall
            if cab_location == 2:
                z += self.altitude_wall - cab.dz

            if cab.reset_location:
                loc = Vector((0, 0, 0))
                a = 0

            if cab_type == 2:
                # rotate R cab
                board_left = cab.panel_left_width
                dx = cab_depth + door_y
                a += pi / 2
                loc += Vector((dx * ca, dx * -sa, 0))

            else:

                rot_type = int(cab.rotate)

                if rot_type == RotationType_R90CW:
                    a += pi / 2
                elif rot_type == RotationType_R90CCW:
                    a -= pi / 2
                elif rot_type == RotationType_R180:
                    a += pi
                elif rot_type == RotationType_User:
                    a += cab.angle

            a = a % (2 * pi)

            # add translation from last to current cab
            ca = cos(a)
            sa = sin(a)
            loc += Vector((
                (board_left + cab.px) * ca + cab.py * sa,
                (board_left + cab.px) * -sa + cab.py * ca,
                0))
            # matrix at cabinet axis with z at bottom
            tM = Matrix([
                [ca, sa, 0, loc.x],
                [-sa, ca, 0, loc.y],
                [0, 0, 1, loc.z],
                [0, 0, 0, 1]
                ])

            # corner 45 relocate doors
            if cab_type in {3, 7}:
                tM2 = tM.copy()
                a += pi / 4
                sy = sq2 * door_y
                sx -= 2 * sy
                dx = -0.5 * sx
                tM = tM * Matrix([
                    [sq2, sq2, 0, radius * 0.5],
                    [-sq2, sq2, 0, radius * 0.5 - x],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]
                    ])
            # ------------
            # manipulators
            # ------------
            # dumb manipulator for door size
            m_dist = (1 + (cab.location % 2)) * 0.25
            m_left = tM * Vector((dx, -sy, z))
            m_right = tM * Vector((dx + sx, -sy, z))
            cab.manipulators[3].set_pts([
                    m_left,
                    m_right,
                    (m_dist, 0, 0)
                    ])

            z0 = 0
            zmax = self.cabinet_height(cab)

            if cab_type != 8:
                for module in cab.modules:

                    n_panels = module.n_panels

                    pivot = 1
                    if module.type in {'2', '12'}:
                        pivot = -1

                    # Door height
                    zd = self.door_height(module)
                    zd = min(zd, zmax - z0)
                    if zd == 0:
                        continue

                    # No door, accumulate spacing
                    if n_panels == 0:
                        z0 += zd
                        continue

                    # range hood
                    # No door, accumulate spacing
                    if module.type == '55':
                        z0 += zd
                        continue

                    size = Vector((sx / n_panels - self.door_gap, zd - self.door_gap, door_y))

                    # handle on top or bottom of door panel
                    # choose location nearest half height
                    handle_location = -1
                    altitude_ref = 0.5 * self.height_full - (z + z0)
                    dist = abs(altitude_ref)
                    if abs(altitude_ref - 0.5 * zd) < dist:
                        handle_location = 0
                    if abs(altitude_ref - zd) < dist:
                        handle_location = 1

                    for i in range(n_panels):

                        if i > 0:
                            pivot = -pivot

                        if panel >= n_childs:
                            name = 'Cabinet door'
                            # Create a panel mesh
                            m = bpy.data.meshes.new(name)
                            # 22.4 deg
                            m.auto_smooth_angle = 0.390953
                            child = bpy.data.objects.new(name, m)
                            props = m.archipack_kitchen_module.add()
                            context.scene.objects.link(child)
                            child.select = True
                            context.scene.objects.active = child
                            m = child.archipack_material.add()
                            m.category = "kitchen"
                            m.material = o.archipack_material[0].material
                            child.lock_location = (False, False, True)
                            child.lock_rotation = (False, True, False)
                            child.lock_scale = (True, True, True)
                            child.show_transparent = True
                            # parenting at 0, 0, 0 before set object matrix_world
                            # so location remains local from frame
                            child.parent = o
                            child.matrix_world = o.matrix_world.copy()

                        else:
                            child = childs[panel]
                            child.select = True
                            context.scene.objects.active = child
                            props = archipack_kitchen_module.datablock(child)

                        if props is not None:
                            props.update(context, size, pivot, self, cab, module, handle_location)

                        child.select = False
                        half_gap = 0.5 * self.door_gap
                        # location y + frame width.
                        child.location = tM * Vector((
                            dx + half_gap + 0.5 * (1 - pivot) * (sx - 2 * half_gap),
                            -sy,
                            z + z0 + half_gap))

                        child.rotation_euler.z = -a
                        panel += 1

                    z0 += zd

            if cab.lock_p:
                dx = 0
                dy = 0
                dz = cab.pz
            else:
                dx = -cab.px
                dy = -cab.py
                dz = 0

            if cab_type == 1:
                # Corner L
                a += pi / 2
                # add translation for current cab to border
                s_right = door_y + cab.panel_right_width
                dx += max(cab_depth + s_right, cab.x)
                loc = tM * Vector((dx, dy - cab_depth - door_y, dz))
            elif cab_type == 2:
                # Corner R
                s_left = door_y + cab.panel_left_width
                dx += max(cab_depth + s_left, cab.x) + board_right - board_left
                loc = tM * Vector((dx, 0, dz))
            elif cab_type in {3, 7}:
                # Corner 45
                a += pi / 4
                dx += max(cab_depth + max(0, cab.dy), cab.x)
                dy -= max(cab_depth + max(0, cab.dy), cab.x) + board_right
                loc = tM2 * Vector((dx, dy, dz))
            else:
                # add translation for current cab to border
                loc = tM * Vector((dx + cab.x + board_right, dy, dz))

    def setup_manipulators(self):
        # dumb base height
        n_manips = len(self.manipulators)

        if n_manips < 1:
            s = self.manipulators.add()
            # dumb oposite side
            s.type_key = 'DUMB_SIZE'

        for i, cab in enumerate(self.cabinets):
            cab.setup_manipulators(i)

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

        # dumb height of base cabs
        self.manip_stack.append(self.manipulators[0].setup(context, o, self))

        for i, part in enumerate(self.cabinets):
            # dumb number
            self.manip_stack.append(part.manipulators[0].setup(context, o, part))
            # width
            self.manip_stack.append(part.manipulators[1].setup(context, o, part))
            if part.cab_type in {1, 2, 3, 7}:
                # dumb size side 2
                self.manip_stack.append(part.manipulators[2].setup(context, o, part))
                # dumb size door
                self.manip_stack.append(part.manipulators[3].setup(context, o, part))

    def manipulable_invoke(self, context):
        """
            call this in operator invoke()
        """
        # print("manipulable_invoke")
        if self.manipulate_mode:
            self.manipulable_disable(context)
            return False

        self.manipulable_setup(context)
        self.manipulate_mode = True

        self._manipulable_invoke(context)

        return True

    def update_parts(self):

        # remove cabinets
        for i in range(len(self.cabinets), self.cabinet_num, -1):
            self.cabinets.remove(i - 1)

        # add cabinets
        for i in range(len(self.cabinets), self.cabinet_num):
            ref = self.cabinets[-1]
            c = self.cabinets.add()
            c.from_cabinet(ref)

        # here to prevent cab update
        # when loading presets
        for cab in self.cabinets:
            cab.update_parts()
            if cab.uid == 0:
                self.create_uid(cab, increment=4)
        
        self.setup_manipulators()

    def update(self, context, manipulable_refresh=False):
        # support for "copy to selected"
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        if manipulable_refresh:
            self.manipulable_disable(context)

        self.update_parts()

        verts = []
        faces = []
        matids = []
        uvs = []

        z = self.height_default
        if self.counter:
            z += self.counter_z
        if self.baseboard:
            z += self.base_height

        self.manipulators[0].set_pts([(0, 0, 0), (0, 0, z), (0.5, 0, 0)], normal=o.matrix_world[1].to_3d())

        if self.counter:
            self.update_counter(verts, faces, matids, uvs)

        self.update_cabinets(context, verts, faces, matids, uvs)

        if self.baseboard:
            self.update_baseboard(verts, faces, matids, uvs)

        bmed.buildmesh(context, o, verts, faces, matids)
        # unwrap mesh
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.uv.cube_project()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        #
        self.update_modules(context, o)

        if manipulable_refresh:
            self.manipulable_refresh = True
        
        self.update_dimensions(context, o)
        
        # restore context
        self.restore_context(context)


class ARCHIPACK_PT_kitchen(Panel):
    bl_idname = "ARCHIPACK_PT_kitchen"
    bl_label = "Kitchen"
    bl_description = "Cabinet Generator"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    # bl_context = 'object'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_kitchen.filter(context.active_object)

    def draw(self, context):
        o = context.active_object
        prop = archipack_kitchen.datablock(o)

        if prop is None:
            return

        layout = self.layout
        layout.operator('archipack.manipulate', icon='HAND')
        """
        row = layout.row(align=True)
        row.operator('archipack.kitchen', text="Refresh", icon='FILE_REFRESH').mode = 'REFRESH'
        if o.data.users > 1:
            row.operator('archipack.kitchen', text="Make unique", icon='UNLINKED').mode = 'UNIQUE'
        """
        layout.operator('archipack.kitchen', text="Delete", icon='ERROR').mode = 'DELETE'
        box = layout.box()
        row = box.row(align=True)
        row.operator("archipack.kitchen_preset_menu", text=bpy.types.ARCHIPACK_OT_kitchen_preset_menu.bl_label)
        row.operator("archipack.kitchen_preset", text="", icon='ZOOMIN')
        row.operator("archipack.kitchen_preset", text="", icon='ZOOMOUT').remove_active = True

        box = layout.box()
        if prop.expand:
            box.prop(prop, "expand", icon="TRIA_DOWN", icon_only=True, text="Kitchen", emboss=False)

            row = box.row()
            row.prop(prop, 'z_mode', text="")
            if prop.z_mode == '1':
                row.prop(prop, 'module_size')
                row = box.row()
                row.prop(prop, 'modules_default')
                row.prop(prop, 'y')
                row = box.row()
                row.prop(prop, 'modules_wall')
                row.prop(prop, 'yw')
                box.prop(prop, 'modules_full')
            else:
                row = box.row()
                row.prop(prop, 'z_default')
                row.prop(prop, 'y')
                row = box.row()
                row.prop(prop, 'z_wall')
                row.prop(prop, 'yw')
                box.prop(prop, 'z_full')

            row = box.row()
            row.prop(prop, 'thickness')

            box = layout.box()
            box.prop(prop, 'door_style')
            row = box.row()
            row.prop(prop, 'door_y')
            row.prop(prop, 'door_gap')

            door_style = int(prop.door_style)
            if door_style > 10:
                row = box.row()
                row.prop(prop, 'door_chanfer')

            if door_style in {2, 12, 13}:
                row.prop(prop, 'door_x')

            if door_style == 13:
                box.prop(prop, 'door_board_chanfer')

            box = layout.box()
            box.prop(prop, 'handle', text="Handle")
            if prop.handle != "0":
                handle_style = int(prop.handle)
                handle_fit = handle_style < 150
                if handle_fit:
                    box.prop(prop, 'handle_fit')
                row = box.row()
                if handle_style not in {180}:
                    row.prop(prop, 'handle_y')
                if handle_style > 149 and handle_style < 160:
                    # Handles
                    row.prop(prop, 'handle_r')
                else:
                    # Bars
                    if handle_style not in {180, 190}:
                        row.prop(prop, 'handle_x')
                    row.prop(prop, 'handle_z')

                if handle_style not in {180, 190}:
                    row = box.row()
                    row.prop(prop, 'handle_dx')
                    if handle_fit and prop.handle_fit:
                        row.prop(prop, 'handle_space')
                    else:
                        row.prop(prop, 'handle_dz')

            box = layout.box()
            box.prop(prop, "counter")
            if prop.counter:
                row = box.row()
                row.prop(prop, "counter_z")
                row.prop(prop, "counter_y")
                row.prop(prop, "counter_chanfer")
            box = layout.box()
            box.prop(prop, 'baseboard')
            if prop.baseboard:
                row = box.row()
                row.prop(prop, 'base_height')
                row.prop(prop, 'base_sink')

            box = layout.box()
            row = box.row()

        else:
            row = box.row()
            row.prop(prop, "expand", icon="TRIA_RIGHT", icon_only=True, text="Kitchen", emboss=False)

        # Cabinet number
        row.prop(prop, 'cabinet_num')
        # Add menu for cabinets
        for i, cab in enumerate(prop.cabinets):
            box = layout.box()
            cab.draw(box, i, prop)


class ARCHIPACK_PT_kitchen_module(Panel):
    bl_idname = "ARCHIPACK_PT_kitchen_module"
    bl_label = "Kitchen module"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_kitchen_module.filter(context.active_object)

    def draw(self, context):
        layout = self.layout
        layout.operator("archipack.select_parent")


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_kitchen(ArchipackCreateTool, Operator):
    bl_idname = "archipack.kitchen"
    bl_label = "Kitchen"
    bl_description = "Kitchen"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    thickness = FloatProperty(
            name='Thickness', min=0.001, default=0.018, precision=3,
            description='Board thickness',
            )
    y = FloatProperty(
            name='Depth', min=0.001, default=0.59, precision=3,
            description='Default cabinet depth',
            )
    z = FloatProperty(
            name='Height', min=0.001, default=0.70, precision=3,
            description='Default cabinet height',
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

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')

    def delete(self, context, o):
        if archipack_kitchen.filter(o):
            self.delete_object(context, o)

    def create(self, context):
        m = bpy.data.meshes.new("Kitchen")
        o = bpy.data.objects.new("Kitchen", m)
        d = m.archipack_kitchen.add()
        d.thickness = self.thickness
        d.y = self.y
        d.z = self.z
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        self.add_material(o)
        self.load_preset(d)
        # select frame
        o.select = True
        context.scene.objects.active = o
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            if self.mode == 'CREATE':
                o = self.create(context)
                o.location = bpy.context.scene.cursor_location
                o.select = True
                context.scene.objects.active = o
                self.manipulate()
            else:
                o = context.active_object
                bpy.ops.archipack.disable_manipulate()
                self.delete(context, o)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to add / remove parts
# ------------------------------------------------------------------


class ARCHIPACK_OT_kitchen_insert(Operator):
    bl_idname = "archipack.kitchen_insert"
    bl_label = "Insert"
    bl_description = "Insert part"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    index = IntProperty(default=0)

    def execute(self, context):
        if context.mode == "OBJECT":
            d = archipack_kitchen.datablock(context.active_object)
            if d is None:
                return {'CANCELLED'}
            d.insert_part(context, self.index)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_kitchen_remove(Operator):
    bl_idname = "archipack.kitchen_remove"
    bl_label = "Remove"
    bl_description = "Remove part"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    index = IntProperty(default=0)

    def execute(self, context):
        if context.mode == "OBJECT":
            d = archipack_kitchen.datablock(context.active_object)
            if d is None:
                return {'CANCELLED'}
            d.remove_part(context, self.index)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to load / save presets
# ------------------------------------------------------------------


class ARCHIPACK_OT_kitchen_preset_menu(PresetMenuOperator, Operator):
    bl_description = "Show Kitchen Presets"
    bl_idname = "archipack.kitchen_preset_menu"
    bl_label = "Kitchen Presets"
    preset_subdir = "archipack_kitchen"


class ARCHIPACK_OT_kitchen_preset(ArchipackPreset, Operator):
    """Add a Kitchen Preset"""
    bl_idname = "archipack.kitchen_preset"
    bl_label = "Add Kitchen Preset"
    preset_menu = "ARCHIPACK_OT_kitchen_preset_menu"

    @property
    def blacklist(self):
        return ['manipulators']


def register():
    bpy.utils.register_class(archipack_kitchen_module)
    Mesh.archipack_kitchen_module = CollectionProperty(type=archipack_kitchen_module)
    bpy.utils.register_class(ARCHIPACK_PT_kitchen_module)
    bpy.utils.register_class(archipack_kitchen_cabinet)
    bpy.utils.register_class(archipack_kitchen)
    Mesh.archipack_kitchen = CollectionProperty(type=archipack_kitchen)
    bpy.utils.register_class(ARCHIPACK_OT_kitchen_preset_menu)
    bpy.utils.register_class(ARCHIPACK_PT_kitchen)
    bpy.utils.register_class(ARCHIPACK_OT_kitchen)
    bpy.utils.register_class(ARCHIPACK_OT_kitchen_remove)
    bpy.utils.register_class(ARCHIPACK_OT_kitchen_insert)
    bpy.utils.register_class(ARCHIPACK_OT_kitchen_preset)


def unregister():
    bpy.utils.unregister_class(archipack_kitchen_cabinet)
    bpy.utils.unregister_class(archipack_kitchen_module)
    bpy.utils.unregister_class(ARCHIPACK_PT_kitchen_module)
    del Mesh.archipack_kitchen_module
    bpy.utils.unregister_class(archipack_kitchen)
    del Mesh.archipack_kitchen
    bpy.utils.unregister_class(ARCHIPACK_OT_kitchen_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_PT_kitchen)
    bpy.utils.unregister_class(ARCHIPACK_OT_kitchen)
    bpy.utils.unregister_class(ARCHIPACK_OT_kitchen_remove)
    bpy.utils.unregister_class(ARCHIPACK_OT_kitchen_insert)
    bpy.utils.unregister_class(ARCHIPACK_OT_kitchen_preset)
