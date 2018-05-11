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
# Author: Jacob Morris - Stephen Leger (s-leger)
# ----------------------------------------------------------
import time
import bpy
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty, CollectionProperty, StringProperty,
    BoolProperty, IntProperty, EnumProperty
    )
from mathutils import Vector, Matrix
from random import uniform
from math import radians, cos, sin, pi, sqrt
import bmesh
from .bmesh_utils import BmeshEdit as bmed
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_object import ArchipackCreateTool, ArchipackObject
from .archipack_cutter import (
    CutAblePolygon, CutAbleGenerator,
    ArchipackCutter,
    ArchipackCutterPart,
    update_operation
    )
from .archipack_dimension import DimensionProvider
from .archipack_curveman import ArchipackUserDefinedPath
from .archipack_segments import ArchipackSegment
from .archipack_throttle import throttle
import logging
logger = logging.getLogger("archipack_floor")

# ------------------------------------------------------------------
# Define property class to store object parameters and update mesh
# ------------------------------------------------------------------


class FloorGenerator(CutAblePolygon, CutAbleGenerator):

    def __init__(self, d, o=None):
        CutAbleGenerator.__init__(self, d, o)
        self.tM = Matrix.Rotation(d.rotation, 4, Vector((0, 0, 1)))

    def limits(self):
        itM = self.tM.inverted()
        pts = [itM * s.p0.to_3d() for s in self.segs]
        x_size = [p.x for p in pts]
        y_size = [p.y for p in pts]
        for s in self.segs:
            if "Curved" in type(s).__name__:
                c = itM * s.c.to_3d()
                x_size.append(c.x + s.r)
                x_size.append(c.x - s.r)
                y_size.append(c.y + s.r)
                y_size.append(c.y - s.r)
        self.xmin = min(x_size)
        self.xmax = max(x_size)
        self.xsize = self.xmax - self.xmin
        self.ymin = min(y_size)
        self.ymax = max(y_size)
        self.ysize = self.ymax - self.ymin

    def cut(self, context, o, d):
        """
            either external or holes cuts
        """
        self.as_lines()
        self.limits()
        self.is_convex()
        for b in o.children:
            d = archipack_floor_cutter.datablock(b)
            if d is not None:
                tM = o.matrix_world.inverted() * b.matrix_world
                g = d.ensure_direction(tM)
                # g.change_coordsys(b.matrix_world, o.matrix_world)
                self.slice(g)

    def floor(self, context, o, d, active):
        """
         active: throttle mode enabled
        """
        verts, faces, matids, uvs = [], [], [], []

        if d.bevel:
            bevel = d.bevel_amount
        else:
            bevel = 0

        # mesh for realtime
        realtime = active or d.pattern == "realtime"

        thickness = d.thickness
        bottom = 0

        if d.add_grout:
            thickness = min(d.thickness - d.mortar_depth, d.thickness - 0.0001)
            bottom = min(d.thickness - (d.mortar_depth + bevel), d.thickness - 0.0001)

        self.top = d.thickness

        if not realtime:
            self.generate_pattern(d, verts, faces, matids, uvs)
            bm = bmed.buildmesh(
                    context, o, verts, faces, matids=matids, uvs=uvs,
                    weld=False, clean=False, auto_smooth=True, temporary=True)

            self.cut_holes(bm, self)
            self.cut_boundary(bm, self)

            bmesh.ops.dissolve_limit(bm,
                        angle_limit=0.01,
                        use_dissolve_boundaries=False,
                        verts=bm.verts,
                        edges=bm.edges,
                        delimit=1)

            bm.verts.ensure_lookup_table()

            if d.solidify:
                # solidify and floor bottom
                geom = bm.faces[:]
                verts = bm.verts[:]
                edges = bm.edges[:]
                bmesh.ops.solidify(bm, geom=geom, thickness=0.0001)
                for v in verts:
                    v.co.z = bottom

                # bevel
                if d.bevel:
                    for v in bm.verts:
                        v.select = True
                    for v in verts:
                        v.select = False
                    for v in bm.edges:
                        v.select = True
                    for v in edges:
                        v.select = False
                    geom = [v for v in bm.verts if v.select]
                    geom.extend([v for v in bm.edges if v.select])
                    bmesh.ops.bevel(bm,
                        geom=geom,
                        offset=d.bevel_amount,
                        offset_type=0,
                        segments=1,     # d.bevel_res
                        profile=0.5,
                        vertex_only=False,
                        clamp_overlap=False,
                        material=-1)

            bm.to_mesh(o.data)
            bm.free()
        else:
            bmed.emptymesh(context, o)

        # Grout
        if realtime or d.add_grout:
            verts = []
            self.get_verts(verts)

            if realtime:
                for v in verts:
                    v.z = self.top

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
            self.cut_boundary(bm, self)

            bmesh.ops.dissolve_limit(bm,
                        angle_limit=0.01,
                        use_dissolve_boundaries=False,
                        verts=bm.verts,
                        edges=bm.edges,
                        delimit=1)

            bm.verts.ensure_lookup_table()
            bm.faces.ensure_lookup_table()

            geom = bm.faces[:]
            # realtime might use extrude instead ?
            if not realtime:
                bmesh.ops.solidify(bm, geom=geom, thickness=thickness)
            bmed.bmesh_join(context, o, [bm], normal_update=True)

        bpy.ops.object.mode_set(mode='OBJECT')

    # ---------------------------------------------------
    # Patterns
    # ---------------------------------------------------

    def regular_tile(self, d, verts, faces, matids, uvs):
        """
         ____  ____  ____
        |    ||    ||    | Regular tile, rows can be offset, either manually or randomly
        |____||____||____|
           ____  ____  ____
          |    ||    ||    |
          |____||____||____|
        """
        off = False
        o = 1
        if d.offset != 0:
            o = 1 / (100 / d.offset)
        y = self.ymin

        while y < self.ymax:
            x = self.xmin
            tl2 = d.tile_length
            if y < self.ymax < y + d.tile_length:
                tl2 = self.ymax - y

            while x < self.xmax:
                tw2 = d.tile_width

                if x < self.xmax < x + d.tile_width:
                    tw2 = self.xmax - x
                elif x == self.xmin and off and not d.random_offset:
                    tw2 = d.tile_width * o
                elif x == self.xmin and d.random_offset:
                    v = d.tile_width * d.offset_variance * 0.0049
                    tw2 = (d.tile_width / 2) + uniform(-v, v)

                self.add_plane(d, verts, faces, matids, uvs, x, y, tw2, tl2)
                x += tw2 + d.spacing

            y += tl2 + d.spacing
            off = not off

    def hopscotch(self, d, verts, faces, matids, uvs):
        """
         ____  _  Large tile, plus small one on top right corner
        |    ||_|
        |____| ____  _  But shifted up so next large one is right below previous small one
              |    ||_|
              |____|
        """
        sp = d.spacing

        # movement variables
        row = 0

        tw = d.tile_width
        tl = d.tile_length
        s_tw = (tw - sp) / 2  # small tile width
        s_tl = (tl - sp) / 2  # small tile length
        y = self.ymin - s_tl

        pre_y = y
        while y < self.ymax + s_tl or (row == 2 and y - sp < self.ymax):
            x = self.xmin
            step_back = True

            if row == 1:  # row start indented slightly
                x = self.xmin + s_tw + sp

            while x < self.xmax:
                if row == 0 or row == 1:
                    # adjust for if there is a need to cut off the bottom of the tile
                    if y < self.ymin - s_tl:
                        self.add_plane(d, verts, faces, matids, uvs, x, y, tw, tl + y - self.ymin)  # large one
                    else:
                        self.add_plane(d, verts, faces, matids, uvs, x, y, tw, tl)  # large one

                    self.add_plane(d, verts, faces, matids, uvs, x + tw + sp, y + s_tl + sp, s_tw, s_tl)  # small one

                    if step_back:
                        x += tw + sp
                        y -= s_tl + sp
                    else:
                        x += tw + s_tw + 2 * sp
                        y += s_tl + sp

                    step_back = not step_back
                else:
                    if x == self.xmin:  # half width for starting position
                        self.add_plane(d, verts, faces, matids, uvs, x, y, s_tw, tl)  # large one
                        # small one on right
                        self.add_plane(d, verts, faces, matids, uvs, x + s_tw + sp, y + s_tl + sp, s_tw, s_tl)
                        # small one on bottom
                        self.add_plane(d, verts, faces, matids, uvs, x, y - sp - s_tl, s_tw, s_tl)
                        x += (2 * s_tw) + tw + (3 * sp)
                    else:
                        self.add_plane(d, verts, faces, matids, uvs, x, y, tw, tl)  # large one
                        # small one on right
                        self.add_plane(d, verts, faces, matids, uvs, x + tw + sp, y + s_tl + sp, s_tw, s_tl)
                        x += (2 * tw) + (3 * sp) + s_tw

            if row == 0 or row == 2:
                y = pre_y + tl + sp
            else:
                y = pre_y + s_tl + sp
            pre_y = y

            row = (row + 1) % 3  # keep wrapping rows

    def stepping_stone(self, d, verts, faces, matids, uvs):
        """
         ____  __  ____
        |    ||__||    | Row of large one, then two small ones stacked beside it
        |    | __ |    |
        |____||__||____|
         __  __  __  __
        |__||__||__||__| Row of smalls
        """
        sp = d.spacing
        y = self.ymin
        row = 0

        tw = d.tile_width
        tl = d.tile_length
        s_tw = (tw - sp) / 2
        s_tl = (tl - sp) / 2

        while y < self.ymax:
            x = self.xmin

            while x < self.xmax:
                if row == 0:  # large one then two small ones stacked beside it
                    self.add_plane(d, verts, faces, matids, uvs, x, y, tw, tl)
                    self.add_plane(d, verts, faces, matids, uvs, x + tw + sp, y, s_tw, s_tl,)
                    self.add_plane(d, verts, faces, matids, uvs, x + tw + sp, y + s_tl + sp, s_tw, s_tl)
                    x += tw + s_tw + (2 * sp)
                else:  # row of small ones
                    self.add_plane(d, verts, faces, matids, uvs, x, y, s_tw, s_tl)
                    self.add_plane(d, verts, faces, matids, uvs, x + s_tw + sp, y, s_tw, s_tl)
                    x += tw + sp

            if row == 0:
                y += tl + sp
            else:
                y += s_tl + sp

            row = (row + 1) % 2

    def hexagon(self, d, verts, faces, matids, uvs):
        """
          __  Hexagon tiles
        /   \
        \___/
        """
        sp = d.spacing
        width = d.tile_width
        dia = (width / 2) / cos(radians(30))
        #               top of current, half way up next,    vertical spacing component
        vertical_spacing = dia * (1 + sin(radians(30))) + (sp * sin(radians(60)))  # center of one row to next row
        da = pi / 3
        base_points = [(sin(i * da), cos(i * da)) for i in range(6)]

        y = self.ymin
        offset = False
        while y - width / 2 < self.ymax:  # place tile as long as bottom is still within bounds
            if offset:
                x = self.xmin + width / 2
            else:
                x = self.xmin - sp / 2

            while x - width / 2 < self.xmax:  # place tile as long as left is still within bounds
                f = len(verts)

                if d.vary_thickness and d.thickness_variance > 0:
                    v = d.thickness / 100 * d.thickness_variance
                    z = uniform(self.top, self.top + v)
                else:
                    z = self.top

                for pt in base_points:
                    verts.append(self.tM * Vector((dia * pt[0] + x, dia * pt[1] + y, z)))

                faces.append([f] + [i for i in range(f + 1, len(verts))])
                uvs.append(base_points)
                self.add_matid(d, matids)

                x += width + sp

            y += vertical_spacing
            offset = not offset

    def windmill(self, d, verts, faces, matids, uvs):
        """
         __  ____
        |  ||____| This also has a square one in the middle, totaling 5 tiles per pattern
        |__|   __
         ____ |  |
        |____||__|
        """
        sp = d.spacing

        tw = d.tile_width
        tl = d.tile_length
        s_tw = (tw - sp) / 2
        s_tl = (tl - sp) / 2

        y = self.ymin
        while y < self.ymax:
            x = self.xmin

            while x < self.xmax:
                self.add_plane(d, verts, faces, matids, uvs, x, y, tw, s_tl)  # bottom
                self.add_plane(d, verts, faces, matids, uvs, x + tw + sp, y, s_tw, tl, rotate_uv=True)  # right
                self.add_plane(d, verts, faces, matids, uvs, x + s_tw + sp, y + tl + sp, tw, s_tl)  # top
                self.add_plane(d, verts, faces, matids, uvs, x, y + s_tl + sp, s_tw, tl, rotate_uv=True)  # left
                self.add_plane(d, verts, faces, matids, uvs, x + s_tw + sp, y + s_tl + sp, s_tw, s_tl)  # center

                x += tw + s_tw + (2 * sp)
            y += tl + s_tl + (2 * sp)

    def boards(self, d, verts, faces, matids, uvs):
        """
        ||| Typical wood boards
        |||
        """
        x = self.xmin
        bw, bl = d.board_width, d.board_length
        off = False
        o = 1 / (100 / d.offset) if d.offset != 0 else 0

        while x < self.xmax:
            if d.vary_width:
                v = bw * (d.width_variance / 100) * 0.99
                bw2 = bw + uniform(-v, v)
            else:
                bw2 = bw

            if bw2 + x > self.xmax:
                bw2 = self.xmax - x
            y = self.ymin

            counter = 1
            while y < self.ymax:
                bl2 = bl
                if d.vary_length:
                    v = bl * (d.length_variance / 100) * 0.99
                    bl2 = bl + uniform(-v, v)
                elif y == self.ymin and off and not d.random_offset:
                    bl2 = bl * o
                elif y == self.ymin and d.random_offset:
                    v = bl * d.offset_variance * 0.0049
                    bl2 = (bl / 2) + uniform(-v, v)

                if (counter >= d.max_boards and d.vary_length) or y + bl2 > self.ymax:
                    bl2 = self.ymax - y

                self.add_plane(d, verts, faces, matids, uvs, x, y, bw2, bl2, rotate_uv=True)
                y += bl2 + d.length_spacing
                counter += 1
            off = not off
            x += bw2 + d.width_spacing

    def square_parquet(self, d, verts, faces, matids, uvs):
        """
        ||--||-- Alternating groups oriented either horizontally, or forwards and backwards.
        ||--||-- self.spacing is used because it is the same spacing for width and length
        --||--|| Board width is calculated using number of boards and the length.
        --||--||
        """
        x = self.xmin
        start_orient_length = True

        # figure board width
        bl = d.short_board_length
        bw = (bl - (d.boards_in_group - 1) * d.spacing) / d.boards_in_group
        while x < self.xmax:
            y = self.ymin
            orient_length = start_orient_length
            while y < self.ymax:

                if orient_length:
                    start_x = x

                    for i in range(d.boards_in_group):
                        if x < self.xmax and y < self.ymax:
                            self.add_plane(d, verts, faces, matids, uvs, x, y, bw, bl, rotate_uv=True)
                            x += bw + d.spacing

                    x = start_x
                    y += bl + d.spacing

                else:
                    for i in range(d.boards_in_group):
                        if x < self.xmax and y < self.ymax:
                            self.add_plane(d, verts, faces, matids, uvs, x, y, bl, bw)
                            y += bw + d.spacing

                orient_length = not orient_length

            start_orient_length = not start_orient_length
            x += bl + d.spacing

    def herringbone(self, d, verts, faces, matids, uvs):
        """
        Boards are at 45 degree angle, in chevron pattern, ends are angled
        """
        width_dif = d.board_width / cos(radians(45))
        x_dif = d.short_board_length * cos(radians(45))
        y_dif = d.short_board_length * sin(radians(45))
        total_y_dif = width_dif + y_dif
        sp_dif = d.spacing / cos(radians(45))

        y = self.ymin - y_dif
        while y < self.ymax:
            x = self.xmin

            while x < self.xmax:
                # left side

                self.add_face(d, verts, faces, matids, uvs,
                    (x, y, 0), (x + x_dif, y + y_dif, 0),
                    (x + x_dif, y + total_y_dif, 0), (x, y + width_dif, 0))

                x += x_dif + d.spacing

                # right side
                if x < self.xmax:
                    self.add_face(d, verts, faces, matids, uvs,
                        (x, y + y_dif, 0), (x + x_dif, y, 0),
                        (x + x_dif, y + width_dif, 0), (x, y + total_y_dif, 0))
                    x += x_dif + d.spacing

            y += width_dif + sp_dif  # adjust spacing amount for 45 degree angle

    def herringbone_parquet(self, d, verts, faces, matids, uvs):
        """
        Boards are at 45 degree angle, in chevron pattern, ends are square, not angled
        """

        an_45 = 0.5 * sqrt(2)

        x_dif = d.short_board_length * an_45
        y_dif = d.short_board_length * an_45
        y_dif_45 = d.board_width * an_45
        x_dif_45 = d.board_width * an_45
        total_y_dif = y_dif + y_dif_45

        sp_dif = (d.spacing / an_45) / 2  # divide by two since it is used for both x and y
        width_dif = d.board_width / an_45

        y = self.ymin - y_dif
        while y - y_dif_45 < self.ymax:  # continue as long as bottom left corner is still good
            x = self.xmin

            while x - x_dif_45 < self.xmax:  # continue as long as top left corner is still good
                # left side

                self.add_face(d, verts, faces, matids, uvs,
                    (x, y, 0),
                    (x + x_dif, y + y_dif, 0),
                    (x + x_dif - x_dif_45, y + total_y_dif, 0),
                    (x - x_dif_45, y + y_dif_45, 0))

                x += x_dif - x_dif_45 + sp_dif
                y0 = y + y_dif - y_dif_45 - sp_dif

                if x < self.xmax:
                    self.add_face(d, verts, faces, matids, uvs,
                        (x, y0, 0),
                        (x + x_dif, y0 - y_dif, 0),
                        (x + x_dif + x_dif_45, y0 - y_dif + y_dif_45, 0),
                        (x + x_dif_45, y0 + y_dif_45, 0))

                    x += x_dif + x_dif_45 + sp_dif

                else:  # we didn't place the right board, so step ahead far enough the the while loop for x breaks
                    break

            y += width_dif + (2 * sp_dif)

    def add_matid(self, d, matids):
        if d.vary_materials:
            matid = uniform(1, d.matid)
        else:
            matid = d.matid
        matids.append(matid)

    def add_plane(self, d, verts, faces, matids, uvs, x, y, w, l, rotate_uv=False):
        """
        Adds vertices and faces for a place, clip to outer boundaries if clip is True
        :param x: start x position
        :param y: start y position
        :param w: width (in x direction)
        :param l: length (in y direction)
        """

        x1 = x + w
        y1 = y + l

        if d.vary_thickness and d.thickness_variance > 0:
            v = d.thickness / 100 * d.thickness_variance
            z = uniform(self.top, self.top + v)
        else:
            z = self.top

        p = len(verts)
        verts.extend([self.tM * Vector(p) for p in [(x, y, z), (x1, y, z), (x1, y1, z), (x, y1, z)]])
        faces.append([p + 3, p + 2, p + 1, p])
        if rotate_uv:
            uvs.append([(0, 0), (0, 1), (1, 1), (1, 0)])
        else:
            uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])
        self.add_matid(d, matids)

    def add_face(self, d, verts, faces, matids, uvs, p0, p1, p2, p3):
        """
        Adds vertices and faces for a place, clip to outer boundaries if clip is True
        :param x: start x position
        :param y: start y position
        :param w: width (in x direction)
        :param l: length (in y direction)
        """

        if d.vary_thickness and d.thickness_variance > 0:
            v = d.thickness / 100 * d.thickness_variance
            z = uniform(self.top, self.top + v)
        else:
            z = self.top

        p = len(verts)
        verts.extend([self.tM * Vector((v[0], v[1], z)) for v in [p0, p1, p2, p3]])
        faces.append([p + 3, p + 2, p + 1, p])

        uvs.append([(0, 0), (1, 0), (1, 1), (0, 1)])

        self.add_matid(d, matids)

    def add_manipulator(self, name, pt1, pt2, pt3):
        m = self.manipulators.add()
        m.prop1_name = name
        m.set_pts([pt1, pt2, pt3])

    def generate_pattern(self, d, verts, faces, matids, uvs):

        if d.pattern == "boards":
            self.boards(d, verts, faces, matids, uvs)
        elif d.pattern == "square_parquet":
            self.square_parquet(d, verts, faces, matids, uvs)
        elif d.pattern == "herringbone":
            self.herringbone(d, verts, faces, matids, uvs)
        elif d.pattern == "herringbone_parquet":
            self.herringbone_parquet(d, verts, faces, matids, uvs)
        elif d.pattern == "regular_tile":
            self.regular_tile(d, verts, faces, matids, uvs)
        elif d.pattern == "hopscotch":
            self.hopscotch(d, verts, faces, matids, uvs)
        elif d.pattern == "stepping_stone":
            self.stepping_stone(d, verts, faces, matids, uvs)
        elif d.pattern == "hexagon":
            self.hexagon(d, verts, faces, matids, uvs)
        elif d.pattern == "windmill":
            self.windmill(d, verts, faces, matids, uvs)


def update(self, context):
    if self.auto_update:
        self.update(context)


def update_manipulators(self, context):
    if self.auto_update:
        self.update(context, manipulable_refresh=True)


def update_path(self, context):
    self.update_path(context)


class archipack_floor_part(ArchipackSegment, PropertyGroup):
    manipulators = CollectionProperty(type=archipack_manipulator)

    def get_datablock(self, o):
        return archipack_floor.datablock(o)


class archipack_floor(ArchipackObject, ArchipackUserDefinedPath, Manipulable, DimensionProvider, PropertyGroup):
    """
    n_parts = IntProperty(
            name="Parts",
            min=1,
            default=1, update=update_manipulators
            )
    """
    parts = CollectionProperty(type=archipack_floor_part)

    # UI layout related
    """
    parts_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    """
    pattern = EnumProperty(
            name='Floor Pattern',
            items=(("boards", "Boards", ""),
                    ("square_parquet", "Square Parquet", ""),
                    ("herringbone_parquet", "Herringbone Parquet", ""),
                    ("herringbone", "Herringbone", ""),
                    ("regular_tile", "Regular Tile", ""),
                    ("hopscotch", "Hopscotch", ""),
                    ("stepping_stone", "Stepping Stone", ""),
                    ("hexagon", "Hexagon", ""),
                    ("windmill", "Windmill", ""),
                    ("realtime", "Realtime", "Realtime optimized simple mesh")),
            default="boards",
            update=update
            )
    rotation = FloatProperty(
            name='Rotation',
            subtype='ANGLE', unit='ROTATION',
            min=-pi,
            max=pi,
            update=update
            )
    spacing = FloatProperty(
            name='Spacing',
            description='The amount of space between boards or tiles in both directions',
            unit='LENGTH', subtype='DISTANCE',
            min=0,
            default=0.005,
            precision=2,
            update=update
            )
    thickness = FloatProperty(
            name='Thickness',
            description='Thickness',
            unit='LENGTH', subtype='DISTANCE',
            min=0.0,
            default=0.005,
            precision=2,
            update=update
            )
    vary_thickness = BoolProperty(
            name='Random Thickness',
            description='Vary thickness',
            default=False,
            update=update
            )
    thickness_variance = FloatProperty(
            name='Variance',
            description='How much vary by',
            min=0, max=100,
            default=25,
            precision=2,
            subtype='PERCENTAGE',
            update=update
            )

    board_width = FloatProperty(
            name='Width',
            description='The width',
            unit='LENGTH', subtype='DISTANCE',
            min=0.02,
            default=0.2,
            precision=2,
            update=update
            )
    vary_width = BoolProperty(
            name='Random Width',
            description='Vary width',
            default=False,
            update=update
            )
    width_variance = FloatProperty(
            name='Variance',
            description='How much vary by',
            subtype='PERCENTAGE',
            min=1, max=100, default=50,
            precision=2,
            update=update
            )
    width_spacing = FloatProperty(
            name='Width Spacing',
            description='The amount of space between boards in the width direction',
            unit='LENGTH', subtype='DISTANCE',
            min=0,
            default=0.002,
            precision=2,
            update=update
            )

    board_length = FloatProperty(
            name='Length',
            description='The length of the boards',
            unit='LENGTH', subtype='DISTANCE',
            precision=2,
            min=0.02,
            default=2,
            update=update
            )
    short_board_length = FloatProperty(
            name='Length',
            description='The length of the boards',
            unit='LENGTH', subtype='DISTANCE',
            precision=2,
            min=0.02,
            default=2,
            update=update
            )
    vary_length = BoolProperty(
            name='Random Length',
            description='Vary board length',
            default=False,
            update=update
            )
    length_variance = FloatProperty(
            name='Variance',
            description='How much board length can vary by',
            subtype='PERCENTAGE',
            min=1, max=100, default=50,
            precision=2, update=update
            )
    max_boards = IntProperty(
            name='Max Boards',
            description='Max number of boards in one row',
            min=1,
            default=20,
            update=update
            )
    length_spacing = FloatProperty(
            name='Length Spacing',
            description='The amount of space between boards in the length direction',
            unit='LENGTH', subtype='DISTANCE',
            min=0,
            default=0.002,
            precision=2,
            update=update
            )

    # parquet specific
    boards_in_group = IntProperty(
            name='Boards in Group',
            description='Number of boards in a group',
            min=1, default=4,
            update=update
            )

    # tile specific
    tile_width = FloatProperty(
            name='Width',
            description='Width of the tiles',
            unit='LENGTH', subtype='DISTANCE',
            min=0.002,
            default=0.2,
            precision=2,
            update=update
            )
    tile_length = FloatProperty(
            name='Length',
            description='Length of the tiles',
            unit='LENGTH', subtype='DISTANCE',
            precision=2,
            min=0.02,
            default=0.3,
            update=update
            )

    # grout
    add_grout = BoolProperty(
            name='Add Grout',
            description='Add grout',
            default=False,
            update=update
            )
    mortar_depth = FloatProperty(
            name='Depth',
            description='The depth of the mortar from the surface of the tile',
            unit='LENGTH', subtype='DISTANCE',
            precision=2,
            step=0.005,
            min=0,
            default=0.001,
            update=update
            )

    # regular tile
    random_offset = BoolProperty(
            name='Random Offset',
            description='Random amount of offset for each row of tiles',
            update=update, default=False
            )
    offset = FloatProperty(
            name='Offset',
            description='How much to offset each row of tiles',
            min=0, max=100, default=0,
            precision=2,
            update=update
            )
    offset_variance = FloatProperty(
            name='Variance',
            description='How much to vary the offset each row of tiles',
            min=0.001, max=100, default=50,
            precision=2,
            update=update
            )

    # bevel
    bevel = BoolProperty(
            name='Bevel',
            update=update,
            default=False,
            description='Bevel upper faces'
            )
    bevel_amount = FloatProperty(
            name='Bevel',
            description='Bevel amount',
            unit='LENGTH', subtype='DISTANCE',
            min=0.0001, default=0.001,
            precision=2, step=0.05,
            update=update
            )
    solidify = BoolProperty(
            name="Solidify",
            default=True,
            update=update
            )
    vary_materials = BoolProperty(
            name="Random Material",
            default=True,
            description="Vary Material indexes",
            update=update)
    matid = IntProperty(
            name="#variations",
            min=1,
            max=10,
            default=7,
            description="Material index maxi",
            update=update)
    z = FloatProperty(
            name="dumb z",
            description="Dumb z for manipulator placeholder",
            default=0.01,
            options={'SKIP_SAVE'}
            )
    x_offset = FloatProperty(
            name='Offset',
            description='How much to offset boundary',
            default=0,
            precision=2,
            update=update
            )

    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )

    closed = BoolProperty(
            default=True,
            options={'SKIP_SAVE'}
            )
    always_closed = BoolProperty(
            default=True,
            options={'SKIP_SAVE'}
            )

    def get_generator(self, o=None):
        g = FloorGenerator(self, o)
        for part in self.parts:
            g.add_part(part)
        g.set_offset(self.x_offset)
        g.close(self.x_offset)
        return g

    def from_spline(self, context, wM, resolution, spline):

        o = self.find_in_selection(context)

        if o is None:
            return

        pts = self.coords_from_spline(
            spline,
            wM,
            resolution,
            ccw=True,
            close=True
            )

        if len(pts) < 3:
            return

        # pretranslate
        o.matrix_world = Matrix.Translation(pts[0].copy())
        auto_update = self.auto_update
        self.auto_update = False
        self.from_points(pts)
        self.auto_update = auto_update

    def add_manipulator(self, name, pt1, pt2, pt3):
        m = self.manipulators.add()
        m.prop1_name = name
        m.set_pts([pt1, pt2, pt3])

    def update_manipulators(self):
        self.manipulators.clear()  # clear every time, add new ones
        self.add_manipulator("length", (0, 0, 0), (0, self.length, 0), (-0.4, 0, 0))
        self.add_manipulator("width", (0, 0, 0), (self.width, 0, 0), (0.4, 0, 0))

        z = self.thickness

        if self.pattern == "boards":
            self.add_manipulator("board_length", (0, 0, z), (0, self.board_length, z), (0.1, 0, z))
            self.add_manipulator("board_width", (0, 0, z), (self.board_width, 0, z), (-0.2, 0, z))
        elif self.pattern == "square_parquet":
            self.add_manipulator("short_board_length", (0, 0, z), (0, self.short_board_length, z), (-0.2, 0, z))
        elif self.pattern in ("herringbone", "herringbone_parquet"):
            dia = self.short_board_length * cos(radians(45))
            dia2 = self.board_width * cos(radians(45))
            self.add_manipulator("short_board_length", (0, 0, z), (dia, dia, z), (0, 0, z))
            self.add_manipulator("board_width", (dia, 0, z), (dia - dia2, dia2, z), (0, 0, z))
        else:
            tl = self.tile_length
            tw = self.tile_width

            if self.pattern in ("regular_tile", "hopscotch", "stepping_stone"):
                self.add_manipulator("tile_width", (0, tl, z), (tw, tl, z), (0, 0, z))
                self.add_manipulator("tile_length", (0, 0, z), (0, tl, z), (0, 0, z))
            elif self.pattern == "hexagon":
                self.add_manipulator("tile_width", (tw / 2 + self.spacing, 0, z), (tw * 1.5 + self.spacing, 0, z),
                                     (0, 0, 0))
            elif self.pattern == "windmill":
                self.add_manipulator("tile_width", (0, 0, z), (tw, 0, 0), (0, 0, z))
                self.add_manipulator("tile_length", (0, tl / 2 + self.spacing, z), (0, tl * 1.5 + self.spacing, z),
                                     (0, 0, z))

    def setup_manipulators(self):

        if len(self.manipulators) < 1:
            s = self.manipulators.add()
            s.type_key = "SIZE"
            s.prop1_name = "z"
            s.normal = Vector((0, 1, 0))

        self.setup_parts_manipulators('z')

    def update(self, context, manipulable_refresh=False):

        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        tim = time.time()

        # clean up manipulators before any data model change
        if manipulable_refresh:
            self.manipulable_disable(context)

        self.update_parts()

        if len(self.parts) < 3:
            self.restore_context(context)
            return

        throttle.add(context, o, self, 0.1)

        g = self.get_generator()
        g.locate_manipulators()

        for i, seg in enumerate(g.segs):
            g.segs[i] = seg.line

        g.cut(context, o, self)
        g.floor(context, o, self, throttle.is_active(o.name))

        self.update_dimensions(context, o)

        # enable manipulators rebuild
        if manipulable_refresh:
            self.manipulable_refresh = True

        logger.debug("Floor.update() %s :%.4f seconds", o.name, time.time() - tim)
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

        # locate manipulators
        self.get_generator()
        self.manipulable_setup(context)
        self.manipulate_mode = True

        self._manipulable_invoke(context)

        return True


class archipack_floor_cutter_segment(ArchipackCutterPart, PropertyGroup):
    manipulators = CollectionProperty(type=archipack_manipulator)

    def get_datablock(self, o):
        return archipack_floor_cutter.datablock(o)


class archipack_floor_cutter(ArchipackCutter, ArchipackObject, Manipulable, DimensionProvider, PropertyGroup):
    parts = CollectionProperty(type=archipack_floor_cutter_segment)

    def update_points(self, context, o, pts, update_parent=False):
        """
         Create boundary from points
        """
        self.auto_update = False
        self.from_points(pts)
        self.auto_update = True
        if update_parent:
            self.update_parent(context, o)

    def update_parent(self, context, o):
        if o is not None:
            d = archipack_floor.datablock(o.parent)
            if d is not None:
                self.select_object(context, o.parent, True)
                d.update(context)
                self.unselect_object(o.parent)
            self.select_object(context, o, True)


# ------------------------------------------------------------------
# Define panel class to show object parameters in ui panel (N)
# ------------------------------------------------------------------


class ARCHIPACK_PT_floor(Panel):
    bl_idname = "ARCHIPACK_PT_floor"
    bl_label = "Flooring"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Archipack"

    @classmethod
    def poll(cls, context):
        # ensure your object panel only show when active object is the right one
        return archipack_floor.filter(context.active_object)

    def draw(self, context):
        o = context.active_object
        if not archipack_floor.filter(o):
            return
        layout = self.layout
        # retrieve datablock of your object
        props = archipack_floor.datablock(o)
        # manipulate
        row = layout.row(align=True)
        row.operator("archipack.manipulate", icon="HAND")
        row.operator("archipack.floor", text="Delete", icon='ERROR').mode = 'DELETE'
        box = layout.box()
        row = box.row(align=True)
        not_realtime = props.pattern != "realtime"
        # Presets operators
        row.operator("archipack.floor_preset_menu",
                     text=bpy.types.ARCHIPACK_OT_floor_preset_menu.bl_label)
        row.operator("archipack.floor_preset",
                      text="",
                      icon='ZOOMIN')
        row.operator("archipack.floor_preset",
                      text="",
                      icon='ZOOMOUT').remove_active = True

        box = layout.box()
        box.operator('archipack.floor_cutter').parent = o.name

        box = layout.box()
        expand = props.template_user_path(context, box)
        if expand:
            box.prop(props, 'x_offset')

        props.template_parts(context, layout)

        layout.separator()
        box = layout.box()
        box.prop(props, 'pattern', text="")

        if not_realtime:
            box.prop(props, 'rotation')

        # thickness
        box.separator()
        box.prop(props, 'thickness')
        if not_realtime:
            box.prop(props, 'vary_thickness', icon='RNDCURVE')
            if props.vary_thickness:
                box.prop(props, 'thickness_variance')
            box.separator()
            box.prop(props, 'solidify', icon='MOD_SOLIDIFY')
            box.separator()
            if props.pattern == 'boards':
                box.prop(props, 'board_length')
                box.prop(props, 'vary_length', icon='RNDCURVE')
                if props.vary_length:
                    box.prop(props, 'length_variance')
                    box.prop(props, 'max_boards')
                box.separator()

                # width
                box.prop(props, 'board_width')
                # vary width
                box.prop(props, 'vary_width', icon='RNDCURVE')
                if props.vary_width:
                    box.prop(props, 'width_variance')
                box.separator()
                box.prop(props, 'length_spacing')
                box.prop(props, 'width_spacing')

            elif props.pattern in {'square_parquet', 'herringbone_parquet', 'herringbone'}:
                box.prop(props, 'short_board_length')

                if props.pattern != "square_parquet":
                    box.prop(props, "board_width")
                box.prop(props, "spacing")

                if props.pattern == 'square_parquet':
                    box.prop(props, 'boards_in_group')
            elif props.pattern in {'regular_tile', 'hopscotch', 'stepping_stone', 'hexagon', 'windmill'}:
                # width and length and mortar
                if props.pattern != "hexagon":
                    box.prop(props, "tile_length")
                box.prop(props, "tile_width")
                box.prop(props, "spacing")

            if props.pattern in {"regular_tile", "boards"}:
                box.separator()
                box.prop(props, "random_offset", icon="RNDCURVE")
                if props.random_offset:
                    box.prop(props, "offset_variance")
                else:
                    box.prop(props, "offset")

            # grout
            box.separator()
            box.prop(props, 'add_grout', icon='MESH_GRID')
            if props.add_grout:
                box.prop(props, 'mortar_depth')

            # bevel
            box.separator()
            box.prop(props, 'bevel', icon='MOD_BEVEL')
            if props.bevel:
                box.prop(props, 'bevel_amount')

            box.separator()
            box.prop(props, "vary_materials", icon="MATERIAL")
            if props.vary_materials:
                box.prop(props, "matid")


class ARCHIPACK_PT_floor_cutter(Panel):
    bl_idname = "ARCHIPACK_PT_floor_cutter"
    bl_label = "Floor Cutter"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_floor_cutter.poll(context.active_object)

    def draw(self, context):
        d = archipack_floor_cutter.datablock(context.active_object)
        if d is None:
            return
        layout = self.layout
        d.draw(context, layout, draw_offset=True)


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_floor(ArchipackCreateTool, Operator):
    bl_idname = "archipack.floor"
    bl_label = "Floor"
    bl_description = "Floor"
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
        """
            expose only basic params in operator
            use object property for other params
        """
        m = bpy.data.meshes.new("Floor")
        o = bpy.data.objects.new("Floor", m)
        d = m.archipack_floor.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        angle_90 = pi / 2
        x, y, = 4, 4
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
        # Link object into scene
        self.link_object_to_scene(context, o)
        # select and make active
        self.select_object(context, o, True)
        self.load_preset(d)
        self.add_material(o)
        return o

    def delete(self, context):
        o = context.active_object
        if archipack_floor.filter(o):
            bpy.ops.archipack.disable_manipulate()
            self.delete_object(context, o)

    def execute(self, context):
        if context.mode == "OBJECT":
            if self.mode == 'CREATE':
                bpy.ops.object.select_all(action="DESELECT")
                o = self.create(context)
                o.location = context.scene.cursor_location
                # activate manipulators at creation time
                # select and make active
                self.select_object(context, o, True)
                self.manipulate()
            else:
                self.delete(context)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_floor_from_curve(ArchipackCreateTool, Operator):
    bl_idname = "archipack.floor_from_curve"
    bl_label = "Floor curve"
    bl_description = "Create a floor from a curve"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

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
        bpy.ops.archipack.floor(auto_manipulate=self.auto_manipulate, filepath=self.filepath)
        o = context.active_object
        d = archipack_floor.datablock(o)
        d.user_defined_path = curve.name
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


class ARCHIPACK_OT_floor_from_wall(ArchipackCreateTool, Operator):
    bl_idname = "archipack.floor_from_wall"
    bl_label = "->Floor"
    bl_description = "Create floor(s) from a wall"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o is not None and o.data is not None and 'archipack_wall2' in o.data

    def create(self, context):
        m = bpy.data.meshes.new("Floor")
        o = bpy.data.objects.new("Floor", m)
        d = m.archipack_floor.add()
        d.manipulable_selectable = True
        d.auto_update = False
        self.link_object_to_scene(context, o)
        self.select_object(context, o, True)
        bpy.ops.script.python_file_run(filepath=self.filepath)
        # Link object into scene
        # select and make active
        return o

    def floor_from_wall(self, context, w, wd):
        """
         Create flooring from surrounding wall
         Use slab cutters, windows and doors, T childs walls
        """
        tim = time.time()
        # wall is either a single or collection of polygons
        io, wall, childs = wd.as_geom(context, w, 'FLOORS', [], [], [])

        # find slab holes if any
        cutters = []
        if w.parent:
            ref = w.parent
            for c in ref.children:
                if c.data and "archipack_slab" in c.data:
                    for cutter in c.children:
                        cd = cutter.data
                        if cd and "archipack_slab_cutter" in cd:
                            coords = cd.archipack_slab_cutter[0].get_coords()
                            poly = io.coords_to_polygon(cutter.matrix_world, coords)
                            cutters.append(poly)
        o = None
        # Multipolygon
        if wall.type_id == 6:
            polys = wall.geoms
        else:
            polys = [wall]
        sel = []
        logger.debug("floor_from_wall() curves :%.4f seconds", time.time() - tim)

        for poly in polys:

            boundary = io._to_curve(poly.exterior, "{}-boundary".format(w.name), '2D')
            boundary.location.z = w.matrix_world.translation.z - wd.z_offset
            logger.debug("floor_from_wall() boundary :%.4f seconds", time.time() - tim)
            o = self.create(context)
            sel.append(o)
            o.matrix_world = w.matrix_world.copy()
            d = archipack_floor.datablock(o)
            d.auto_update = False
            d.thickness = wd.z_offset
            logger.debug("floor_from_wall() create :%.4f seconds", time.time() - tim)
            d.user_defined_path = boundary.name
            logger.debug("floor_from_wall() user_defined_path :%.4f seconds", time.time() - tim)
            self.delete_object(context, boundary)
            logger.debug("floor_from_wall() delete_object :%.4f seconds", time.time() - tim)
            d.user_defined_path = ""
            logger.debug("floor_from_wall() floor :%.4f seconds", time.time() - tim)

            for hole in poly.interiors:
                curve = io._to_curve(hole, "{}-cut".format(o.name), '3D')
                bpy.ops.archipack.floor_cutter(auto_manipulate=False, parent=o.name, curve=curve.name)
                c = context.active_object
                cd = archipack_floor_cutter.datablock(c)
                cd.user_defined_path = ""
                self.delete_object(context, curve)
                self.unselect_object(c)
            for cutter in cutters:
                if poly.intersects(cutter):
                    curve = io._to_curve(cutter, "{}-cut".format(o.name), '3D')
                    bpy.ops.archipack.floor_cutter(auto_manipulate=False, parent=o.name, curve=curve.name)
                    c = context.active_object
                    cd = archipack_floor_cutter.datablock(c)
                    cd.user_defined_path = ""
                    self.delete_object(context, curve)
                    self.unselect_object(c)
            logger.debug("floor_from_wall() cutters :%.4f seconds", time.time() - tim)

            # select and make active
            self.select_object(context, o, True)
            d.auto_update = True
            self.unselect_object(o)
            logger.debug("floor_from_wall() %s :%.4f seconds", o.name, time.time() - tim)

        for obj in sel:
            self.select_object(context, obj)

        self.select_object(context, w, True)
        bpy.ops.archipack.add_reference_point()
        self.unselect_object(w)
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            # same behiavour as molding from wall
            # wich require moldings to be selected to be editable
            bpy.ops.object.select_all(action="DESELECT")
            wall = context.active_object
            wd = wall.data.archipack_wall2[0]
            o = self.floor_from_wall(context, wall, wd)
            # select and make active
            self.select_object(context, o, True)
            # if self.auto_manipulate:
            #    bpy.ops.archipack.manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_floor_cutter(ArchipackCreateTool, Operator):
    bl_idname = "archipack.floor_cutter"
    bl_label = "Floor Cutter"
    bl_description = "Floor Cutter"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    parent = StringProperty("")
    curve = StringProperty("")

    def create(self, context):
        m = bpy.data.meshes.new("Floor Cutter")
        o = bpy.data.objects.new("Floor Cutter", m)
        d = m.archipack_floor_cutter.add()
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
        # Link object into scene
        self.link_object_to_scene(context, o)

        # select and make active
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
            # select and make active
            self.select_object(context, o, True)
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to load presets (with polls)
# ------------------------------------------------------------------


class ARCHIPACK_OT_floor_preset_create(PresetMenuOperator, Operator):
    bl_description = "Show Floor presets and create object at cursor location"
    bl_idname = "archipack.floor_preset_create"
    bl_label = "Floor"
    preset_subdir = "archipack_floor"


class ARCHIPACK_OT_floor_preset_menu(PresetMenuOperator, Operator):
    bl_description = "Create Floor from presets"
    bl_idname = "archipack.floor_preset_menu"
    bl_label = "Floor preset"
    preset_subdir = "archipack_floor"


class ARCHIPACK_OT_floor_preset_from_wall(PresetMenuOperator, Operator):
    bl_description = "Create floor(s) from a wall"
    bl_idname = "archipack.floor_preset_from_wall"
    bl_label = "-> Floor"
    preset_subdir = "archipack_floor"
    preset_operator = StringProperty(
        options={'SKIP_SAVE'},
        default="archipack.floor_from_wall"
    )

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o and o.data and "archipack_wall2" in o.data


class ARCHIPACK_OT_floor_preset_from_curve(PresetMenuOperator, Operator):
    bl_description = "Create a Floor from curve"
    bl_idname = "archipack.floor_preset_from_curve"
    bl_label = "-> Floor"
    preset_subdir = "archipack_floor"
    preset_operator = StringProperty(
        options={'SKIP_SAVE'},
        default="archipack.floor_from_curve"
    )

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o and o.type == 'CURVE'


class ARCHIPACK_OT_floor_preset(ArchipackPreset, Operator):
    """Add a Floor Preset"""
    bl_idname = "archipack.floor_preset"
    bl_label = "Add Floor preset"
    preset_menu = "ARCHIPACK_OT_floor_preset_menu"

    @property
    def blacklist(self):
        return ['manipulators', 'parts', 'thickness', 'n_parts', 'user_defined_path', 'user_defined_resolution']


def register():
    bpy.utils.register_class(archipack_floor_cutter_segment)
    bpy.utils.register_class(archipack_floor_cutter)
    Mesh.archipack_floor_cutter = CollectionProperty(type=archipack_floor_cutter)
    bpy.utils.register_class(ARCHIPACK_OT_floor_cutter)
    bpy.utils.register_class(ARCHIPACK_PT_floor_cutter)

    bpy.utils.register_class(archipack_floor_part)
    bpy.utils.register_class(archipack_floor)
    Mesh.archipack_floor = CollectionProperty(type=archipack_floor)
    bpy.utils.register_class(ARCHIPACK_PT_floor)
    bpy.utils.register_class(ARCHIPACK_OT_floor)
    bpy.utils.register_class(ARCHIPACK_OT_floor_preset_from_wall)
    bpy.utils.register_class(ARCHIPACK_OT_floor_preset_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_floor_preset_menu)
    bpy.utils.register_class(ARCHIPACK_OT_floor_preset_create)
    bpy.utils.register_class(ARCHIPACK_OT_floor_preset)
    bpy.utils.register_class(ARCHIPACK_OT_floor_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_floor_from_wall)


def unregister():
    bpy.utils.unregister_class(archipack_floor_cutter_segment)
    bpy.utils.unregister_class(archipack_floor_cutter)
    del Mesh.archipack_floor_cutter
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_cutter)
    bpy.utils.unregister_class(ARCHIPACK_PT_floor_cutter)

    bpy.utils.unregister_class(archipack_floor_part)
    bpy.utils.unregister_class(archipack_floor)
    del Mesh.archipack_floor
    bpy.utils.unregister_class(ARCHIPACK_PT_floor)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_preset_from_wall)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_preset_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_preset_create)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_from_wall)
