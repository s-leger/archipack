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
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import (
    FloatProperty,
    BoolProperty,
    CollectionProperty,
    IntProperty,
    EnumProperty
    )
from math import sin, cos, tan, ceil, floor, pi, asin, acos, radians, atan2, sqrt
from mathutils import Vector, Matrix
from .bmesh_utils import BmeshEdit as bmed
from .archipack_manipulator import Manipulable
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_object import ArchipackCreateTool, ArchipackObject
import logging
logger = logging.getLogger("archipack")

def update(self, context):
    self.update(context)


    
    

class archipack_blind(ArchipackObject, Manipulable, PropertyGroup):
    width = FloatProperty(
            name='Width',
            min=0.10, default=1, precision=3,
            description='Total width', update=update,
            )
    frame_height = FloatProperty(
            name='Height',
            min=0.01, default=0.2, precision=3,
            description='Frame height',
            update=update,
            )
    altitude = FloatProperty(
            name='Altitude',
            min=0, default=1.0, precision=3,
            description='Altitude',
            update=update,
            )
    frame_depth = FloatProperty(
            name='Frame depth', min=0.02, default=0.04,
            precision=3,
            description='Frame depth', update=update,
            )
    height = FloatProperty(
            name='Height',
            min=0.20, max=10, default=1.7, precision=3,
            description='Total height',
            update=update,
            )
    depth = FloatProperty(
            name='Depth', min=0.02, default=0.04,
            precision=3,
            description='Slat depth/width', update=update,
            )
    panels = IntProperty(
            name="Panels",
            min=1,
            default=2,
            description='#panels', update=update
            )
    angle = FloatProperty(
            name='Angle',
            min=-pi,
            max=pi,
            default=0,
            precision=1,
            subtype='ANGLE', unit='ROTATION',
            description='Angle of the slats', update=update,
            )
    blind_angle = FloatProperty(
            name='Angle',
            min=0,
            max=0.5*pi,
            default=0.5*pi,
            precision=1,
            subtype='ANGLE', unit='ROTATION',
            description='Angle of the slats', update=update,
            )
    ratio = FloatProperty(
            name='Extend', min=0, max=100, default=100,
            description='% of extension (100 full extend)', update=update,
            )
    style = EnumProperty(
        items=(
            ('VENITIAN', 'Venitian', 'Venitian', 0),  # /  /  ///
            ('SLAT', 'Slat', 'Slat', 1), # / / /
            ('ROLLER', 'Roller', 'Roller', 2),  # __________
            ('BLADES', 'Blades', 'Blades', 3), # _ _ ___
            ('PLATED', 'Plated', 'Plated', 4),  # Vertical \/\/        -> regular
            # ('CELLULAR', 'Cellular blind', 'Cellular blind', 5),    # OOOO         -> venitian like variant
            ('JAPANESE', 'Japanese', 'Japanese', 6), # _-  -> horizontal
            ('VERTICAL_SLOTTED', 'Vertical slotted', 'Vertical slotted', 7)  # / / -> horizontal
        ),
        default='VENITIAN', update=update
    )
    auto_update = BoolProperty(
            # Wont save auto_update state in any case
            options={'SKIP_SAVE'},
            default=True,
            update=update
            )

    def frame(self, verts, faces, matids):
        nv = len(verts)
        x = 0.5 * self.width
        y = 0.5 * self.frame_depth
        z0 = self.altitude + self.height + self.frame_height
        z1 = self.altitude + self.height
        verts.extend([(-x, -y, z0),
                (-x, y, z0),
                (x, y, z0),
                (x, -y, z0),
                (-x, -y, z1),
                (-x, y, z1),
                (x, y, z1),
                (x, -y, z1)])
        faces.extend([tuple([idx + nv for idx in f]) for f in [
            (0, 3, 2, 1),
            (0, 1, 5, 4),
            (1, 2, 6, 5),
            (2, 3, 7, 6),
            (5, 6, 7, 4),
            (0, 4, 7, 3)]])
        matids.extend([
            0,0,0,0,0,0
            ])

    def cylinder(self, radius, size, x, y, z, axis, verts, faces, matids):
        seg = 12
        deg = 2 * pi / seg
        nv = len(verts)

        if axis == 'X':
            tM = Matrix([
                [0, 0, size,  x],
                [0, radius, 0, y],
                [radius, 0, 0, z],
                [0, 0, 0, 1]
            ])

        elif axis == 'Y':
            tM = Matrix([
                [radius, 0, 0, x],
                [0, 0, size, y],
                [0, radius, 0, z],
                [0, 0, 0, 1]
            ])

        else:
            tM = Matrix([
                [radius, 0, 0, x],
                [0, radius, 0, y],
                [0, 0, size, z],
                [0, 0, 0, 1]
            ])

        verts.extend([tM * Vector((sin(deg * a), cos(deg * a), 0)) for a in range(seg)])
        verts.extend([tM * Vector((sin(deg * a), cos(deg * a), 1)) for a in range(seg)])

        faces.extend([tuple([nv + i + f for f in (0, 1, seg + 1, seg)]) for i in range(seg - 1)])
        faces.append((nv + seg - 1, nv, nv + seg, nv + 2 * seg - 1))
        matids.extend([2 for i in range(seg)])

    def make_blade(self, verts, faces, matids, posz, offset_z, angle):
        # ------------------------------------
        # Mesh data
        # ------------------------------------
        c = cos(angle)
        s = sin(angle)
        x = 0.5 * self.width
        y = 0.5 * self.depth

        tM0 = Matrix([
            [1, 0, 0, 0],
            [0, c, -s, c * posz],
            [0, s, c, s * posz + offset_z],
            [0, 0, 0, 1]
        ])
        tM1 = Matrix([
            [-1, 0, 0, 0],
            [0, c, -s, c * posz],
            [0, s, c, s * posz + offset_z],
            [0, 0, 0, 1]
        ])

        nv = len(verts)

        # blade shape angle
        a = tan(radians(3))

        # y from center
        y0 = -0.0025
        y1 = y - 0.00195

        # z from center to out
        z0 = 0
        z1 = -a * y1
        z2 = -a * y

        # x from out to center
        x1 = x - 0.0045
        x2 = x - 0.0017

        shape =[(x2, -y1, z1),
                (x2, y1, z1),
                (x1, -y, z2),
                (x1, y, z2),
                (x, -y0, z0),
                (x, y0, z0),
                (x1, -y0, z0),
                (x1, y0, z0)]

        # Vertex
        verts.extend([tM0 * Vector(v) for v in shape])
        verts.extend([tM1 * Vector(v) for v in shape])

        # Faces
        faces.extend([tuple([i + nv for i in f]) for f in [
            (6, 4, 1, 3), (7, 5, 4, 6), (2, 0, 5, 7),
            (14, 6,  3, 11), (15, 7, 6, 14), (10, 2, 7, 15),
            (12, 14, 11, 9), (13, 15, 14, 12), (8, 10, 15, 13)]])

        matids.extend([1 for i in range(9)])
        """"""

    def roller_blades(self, verts, faces, matids):
        gap = 0.01
        # total available space
        height = self.height
        bottom = self.altitude
        # numslats
        numslats = ceil(height / self.depth)
        separation = gap + self.depth

        # upper blade may overflow inside frame
        # top when all blades are collapsed
        absolute_top = bottom + numslats * self.depth

        # when up, there are gaps between blades
        unfold_height = numslats * separation

        # fully open upper blade
        unfold_top = absolute_top + unfold_height

        # altitude of upper blade
        upper_altitude = unfold_top - unfold_height * self.ratio / 100

        # available space, to compute collapsed one
        available_space = min(unfold_height, upper_altitude - bottom)
        # space left for gaps
        available_gaps = available_space - numslats * self.depth
        # regular including gaps
        regular_slats = floor(available_gaps / gap)
        regular_space = separation * regular_slats
        collapsed_slats = numslats - regular_slats
        collapsed_space = self.depth * collapsed_slats
        
        top = upper_altitude - (0.5 * self.depth + gap) - (self.altitude + self.height)
        
        for i in range(regular_slats):
            if top <= 0.5 * self.depth:
                self.make_blade(verts, faces, matids, top, self.altitude + self.height, self.blind_angle)
            top -= separation

        if collapsed_slats > 0:
            top -= available_space - (collapsed_space + regular_space) - gap

        for i in range(collapsed_slats):
            self.make_blade (verts, faces, matids, top, self.altitude + self.height, self.blind_angle)
            top -= self.depth
    
    def make_slat(self, verts, faces, matids, posz, angle):
        # ------------------------------------
        # Mesh data
        # ------------------------------------
        x = 0.5 * self.width
        y = 0.5 * self.depth
        c = cos(angle)
        s = sin(angle)

        tM0 = Matrix([
            [1, 0, 0, 0],
            [0, c, s, 0],
            [0, -s, c, posz],
            [0, 0, 0, 1]
        ])
        tM1 = Matrix([
            [-1, 0, 0, 0],
            [0, c, s, 0],
            [0, -s, c, posz],
            [0, 0, 0, 1]
        ])

        # half gap width
        gap = 0.0025

        # gap distance from border
        if self.width < 0.60:
            sep = 0.06
        else:
            sep = 0.15

        nv = len(verts)

        # slat angle
        a = tan(radians(3))

        # y from center
        y0 = gap
        y1 = y - 0.00195

        # z from center to out
        z0 = 0
        z1 = -a * y1
        z2 = -a * y

        # x from out to center
        x1 = x - 0.0045
        x2 = x - 0.0017
        # border gap
        x3 = x - sep - gap
        x4 = x - sep + gap
        # center gap
        x5 = gap

        shape =[(x2, -y1, z1),
                (x2, y1, z1),
                (x1, -y, z2),
                (x1, y, z2),
                (x, -y0, z0),
                (x, y0, z0),
                (x1, -y0, z0),
                (x1, y0, z0),
                (x5, -y, z2),
                (x5, y, z2),
                (x5, -y0, z0),
                (x5, y0, z0),
                (x3, -y, z2),
                (x3, y, z2),
                (x3, -y0, z0),
                (x3, y0, z0),
                (x4, -y, z2),
                (x4, y, z2),
                (x4, -y0, z0),
                (x4, y0, z0)]

        # Vertex
        verts.extend([tM0 * Vector(v) for v in shape])
        verts.extend([tM1 * Vector(v) for v in shape])

        # Faces
        faces.extend([tuple([idx + nv for idx in f]) for f in [
            (7, 5, 1, 3), (6, 4, 5, 7), (2, 0, 4, 6),
            (19, 7, 3, 17), (16, 2, 6, 18), (18, 6, 7, 19),
            (11, 15, 13, 9), (8, 12, 14, 10),
            (10, 14, 15, 11), (15, 19, 17, 13), (12, 16, 18, 14),
            (39, 35, 33, 37), (34, 38, 36, 32),
            (27, 23, 21, 25), (27, 25, 24, 26), (24, 20, 22, 26),
            (39, 37, 23, 27), (39, 27, 26, 38), (38, 26, 22, 36),
            (30, 34, 32, 28), (34, 30, 31, 35), (35, 31, 29, 33),
            (11, 9, 29, 31),    (8, 10, 30, 28)]])

        matids.extend([1 for i in range(24)])

    def roller_slats(self, verts, faces, matids):
        # total available space
        height = self.height
        # numslats (slats +1)
        numslats = ceil(height / (self.depth * 0.85))
        # separation
        separation = height / numslats
        half_separation = 0.5 * separation
        available = height * self.ratio / 100
        regular_slats = available / separation
        top = self.altitude + self.height - (regular_slats % 1 - 2) * separation
        absolute_top = self.altitude + self.height + 0.5 * self.depth

        angle = min(0.47222 * pi, max(-0.47222 * pi, self.angle))

        for i in range(int(regular_slats) + 2):
            top -= half_separation
            if top < absolute_top:
                self.make_slat(verts, faces, matids, top, angle)
            top -= half_separation

        radius = 0.00025
        bottom = top + half_separation
        z = self.altitude + self.height
        y = 0
        if self.width < 0.60:
            sep = 0.06
        else:
            sep = 0.15

        s = bottom - z

        self.cylinder(radius, s, 0, y, z, 'Z', verts, faces, matids)
        self.cylinder(radius, s, 0.5 * self.width - sep, y, z, 'Z', verts, faces, matids)
        self.cylinder(radius, s, sep - 0.5 * self.width, y, z, 'Z', verts, faces, matids)

    def venitian_slats(self, verts, faces, matids):
        gap = 0.0015
        # total available space
        height = self.height
        top = self.altitude + self.height
        angle = min(0.47222 * pi, max(-0.47222 * pi, self.angle))
        # numslats
        numslats = ceil(height / (self.depth * 0.85))
        # separation
        separation = height / numslats
        half_separation = 0.5 * separation
        # available space
        available = (height - numslats * gap) * self.ratio / 100
        collapsed_parts = (available - numslats * (separation - gap)) / (gap - separation)
        collapsed_slats = floor(collapsed_parts)
        collapsed_space = collapsed_slats * gap
        regular_slats = floor(numslats - collapsed_parts)
        regular_space = regular_slats * separation
        rotated_space = available + numslats * gap - (regular_space + collapsed_space)
        rotated_height = abs(sin(angle) * 0.5 * self.depth)
        """
        logger.debug("n_c:%s, n_r:%s, n_s:%s, c:%s, r:%s, rot:%s, a:%s",
            collapsed_slats,
            regular_slats,
            numslats,
            collapsed_space,
            regular_space,
            rotated_space,
            available
        )
        """
        for i in range(regular_slats):
            top -= half_separation
            self.make_slat(verts, faces, matids, top, angle)
            top -= half_separation

        # rotated slat
        if collapsed_slats < numslats and regular_slats < numslats:
            # space available between collapsed and regular
            #
            #  \|
            #   |\ _
            #  _|_ _
            before = min(half_separation, rotated_space)

            # space in use for collapsed
            if before < half_separation:
                rotated_angle = 0
            else:
                rotated_angle = asin(2 * (rotated_space - half_separation) / self.depth)
                if self.angle < 0:
                    rotated_angle = max(self.angle, -rotated_angle)
                else:
                    rotated_angle = min(self.angle, rotated_angle)
                    
            top -= before
            self.make_slat(verts, faces, matids, top, rotated_angle)
            top -= 0.5 * gap + max(0, rotated_space - before)


        for i in range(collapsed_slats):
            top -= 0.5 * gap
            self.make_slat(verts, faces, matids, top, 0)
            top -= 0.5 * gap

        radius = 0.0003
        bottom = top

        if collapsed_slats > 0:
            bottom += 0.5 * gap
        else:
            bottom = max(
                top + 0.5 * gap, 
                self.altitude + self.height - separation * (numslats - 0.5)
                )

        z = self.altitude + self.height

        y = 0
        if self.width < 0.60:
            sep = 0.06
        else:
            sep = 0.15

        s = bottom - z
        self.cylinder(radius, s, 0, y, z, 'Z', verts, faces, matids)
        self.cylinder(radius, s, 0.5 * self.width - sep, y, z, 'Z', verts, faces, matids)
        self.cylinder(radius, s, sep - 0.5 * self.width, y, z, 'Z', verts, faces, matids)

        # curvature of slats (z)
        ca = cos(angle)
        sa = sin(angle)
        slat_depth = tan(radians(3)) * 0.5 * self.depth
        curvature = slat_depth * sa

        rot = 0

        if self.ratio == 100:
            rot = sa * 0.5 * self.depth
            s -= slat_depth

        y = 0.5 * self.depth * ca - curvature + radius

        self.cylinder(radius, s - rot, 0, y, z, 'Z', verts, faces, matids)
        self.cylinder(radius, s - rot, 0.5 * self.width - sep, y, z, 'Z', verts, faces, matids)
        self.cylinder(radius, s - rot, sep - 0.5 * self.width, y, z, 'Z', verts, faces, matids)

        y = -0.5 * self.depth * ca - curvature - radius

        self.cylinder(radius, s + rot, 0, y, z, 'Z', verts, faces, matids)
        self.cylinder(radius, s + rot, 0.5 * self.width - sep, y, z, 'Z', verts, faces, matids)
        self.cylinder(radius, s + rot, sep - 0.5 * self.width, y, z, 'Z', verts, faces, matids)

    def roller_curtain(self, verts, faces, matids):
        x = 0.5 * self.width
        y = 0.5 * self.depth
        c = cos(self.blind_angle)
        s = sin(self.blind_angle)
        z0 = self.altitude + self.height
        ext = -self.height * self.ratio / 100
        z1 = z0 + s * ext
        y = c * ext
        nv = len(verts)
        verts.extend([Vector(v) for v in [
            [x, 0, z0],
            [-x, 0, z0],
            [-x, y, z1],
            [x, y, z1]
            ]])
        faces.append(tuple([nv + f for f in range(4)]))
        matids.append(1)

    def plated_blind(self, verts, faces, matids):
        x = 0.5 * self.width
        a = atan2(0.002, 0.5 * self.depth)
        zmax = self.depth * cos(a)
        num = ceil(self.height / zmax)
        z = self.height / num * self.ratio / 100
        a = asin(z / self.depth)
        y = 0.5 * cos(a) * self.depth
        z0 = self.altitude + self.height

        nv = len(verts)
        for yi, zi in [[(2 * (i % 2) - 1) * y, z0 - z * i] for i in range(num + 1)]:
            verts.extend([Vector((x, yi, zi)), Vector((-x, yi, zi))])
        if num % 2 == 1:
            y = -y

        verts.extend([Vector((x, y, z0 - num * z)), Vector((-x, y, z0 - num * z))])
        faces.extend([tuple([nv + f, nv + 1 + f, nv + 3 + f, nv + 2 + f]) for f in range(0, 2 * (num + 1), 2)])
        matids.extend([1 for i in range(num + 1)])

        radius = 0.0003

        s = -num * z
        z = self.altitude + self.height
        y = 0
        if self.width < 0.60:
            sep = 0.06
        else:
            sep = 0.15

        self.cylinder(radius, s, 0.5 * self.width - sep, y, z, 'Z', verts, faces, matids)
        self.cylinder(radius, s, sep - 0.5 * self.width, y, z, 'Z', verts, faces, matids)
    
    def cathenary(self, num, x, y, z, r, h, w, verts, faces, matids):
        # pts = [round(v.co.z, 4) for v in C.object.data.splines[0].points]
        zs = [-0.0, -0.3056, -0.5556, -0.75, -0.8889, -0.9722]
        seg = 6
        deg = 2 * pi / seg
        da = pi / 24
        nv = len(verts)
        for j in range(num):
            # rotate on y axis
            ay = -pi / 4
            for i, zi in enumerate(zs):
                xi = 1 / 12 * i
                tM = Matrix([
                    [r * sin(ay), 0, r * cos(ay), x + xi * w],
                    [0, r, 0, y],
                    [-r * cos(ay), 0, r * sin(ay), z + zi * h],
                    [0, 0, 0, 1]
                ])
                ay += da
                verts.extend([tM * Vector((sin(deg * a), cos(deg * a), 0)) for a in range(seg)])
            
            # lower vert (center)
            tM = Matrix([
                    [r * sin(ay), 0, r * cos(ay), x + 0.5 * w],
                    [0, r, 0, y],
                    [-r * cos(ay), 0, r * sin(ay), z - h],
                    [0, 0, 0, 1]
                ])
            ay += da
            verts.extend([tM * Vector((sin(deg * a), cos(deg * a), 0)) for a in range(seg)])    
            
            for i, zi in enumerate(reversed(zs)):
                xi = 0.5 + 1 / 12 * (i + 1)    
                tM = Matrix([
                    [r * sin(ay), 0, r * cos(ay), x + xi * w],
                    [0, r, 0, y],
                    [-r * cos(ay), 0, r * sin(ay), z + zi * h],
                    [0, 0, 0, 1]
                ])
                ay += da
                verts.extend([tM * Vector((sin(deg * a), cos(deg * a), 0)) for a in range(seg)])
            
            x += w
            
            for s in range(12):
                faces.extend([tuple([nv + i + f for f in (0, 1, seg + 1, seg)]) for i in range(seg - 1)])
                faces.append((nv + seg - 1, nv, nv + seg, nv + 2 * seg - 1))
                nv += seg
                matids.extend([2 for i in range(seg)])
            nv += seg
            
    def vertical_slotted(self, verts, faces, matids):
        c = cos(self.angle)
        s = sin(self.angle)
        x = 0.5 * self.depth
        z = self.altitude
        h = self.height
        nv = len(verts)
        num = ceil(self.width / (0.85 * self.depth))
        spacing = max(0.002, (self.width / num) * self.ratio / 100)
        x0 = 0.5 * spacing - 0.5 * self.width
        dx = c * x
        dy = s * x
        for i in range(num):
            tM = Matrix([
                [dx, -dy, 0, x0 + i * spacing],
                [dy, dx, 0, 0],
                [0, 0, h, z],
                [0, 0, 0, 1]
            ])
            verts.extend([tM * Vector(v) for v in [
                [1, 0, 1],
                [-1, 0, 1],
                [-1, 0, 0],
                [1, 0, 0]
            ]])
            faces.append(tuple([nv + 4 * i + j for j in range(4)]))
            matids.append(1)
        
        # cathenary scaling estimation
        h = (self.width / num) * sin(acos(self.ratio / 101.003) * 1.01003)
        w = self.ratio / 100
        # pts = [round(v.co.z, 4) for v in C.object.data.splines[0].points]
        
        r = 0.00025
        self.cathenary(num - 1, x0 + dx, dy, z, r, h, spacing, verts, faces, matids)
        self.cathenary(num - 1, x0 - dx, -dy, z, r, h, spacing, verts, faces, matids)
                 
    def japanese(self, verts, faces, matids):
        z = self.altitude
        nv = len(verts)
        num = self.panels
        w  = self.width / num
        h = self.height
        dx = self.ratio / 100 * w
        for i in range(num):
            x = -0.5 * self.width + i * dx
            y = 0.01 * (2 * (i % 2) - 1)
            verts.extend([Vector(v) for v in [
                [x + w, y, h + z],
                [x, y, h + z],
                [x, y, z],
                [x + w, y, z]
            ]])
            faces.append(tuple([nv + 4 * i + j for j in range(4)]))
            matids.append(1)
    
    def setup_manipulators(self):
        if len(self.manipulators) < 1:
            # add manipulator for x property
            s = self.manipulators.add()
            s.prop1_name = "width"
            s.prop2_name = "x"
            s.type_key = 'SNAP_SIZE_LOC'

            # add manipulator for y property
            s = self.manipulators.add()
            s.prop1_name = "height"
            s.type_key = 'SIZE'
            s.normal = Vector((0, 1, 0))

            # add manipulator for z property
            s = self.manipulators.add()
            s.prop1_name = "altitude"
            s.type_key = 'SIZE'
            # draw this one on xz plane
            s.normal = Vector((0, 1, 0))

    def update(self, context):

        # provide support for "copy to selected"
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # dynamically create manipulators when needed
        self.setup_manipulators()
        verts = []
        faces = []
        matids = []
        self.frame(verts, faces, matids)

        if self.style == 'VENITIAN':
            self.venitian_slats(verts, faces, matids)
        elif self.style == 'SLAT':
            self.roller_slats(verts, faces, matids)
        elif self.style == 'BLADES':
            self.roller_blades(verts, faces, matids)
        elif self.style == 'ROLLER':
            self.roller_curtain(verts, faces, matids)
        elif self.style == 'PLATED':
            self.plated_blind(verts, faces, matids)
        elif self.style == 'VERTICAL_SLOTTED':
            self.vertical_slotted(verts, faces, matids)
        elif self.style == 'JAPANESE':
            self.japanese(verts, faces, matids)
            
        # update your mesh from parameters
        bmed.buildmesh(context,
                       o,
                       verts,
                       faces,
                       matids=matids,
                       weld=False)
        
        # update manipulators location (3d location in object coordsystem)
        x, y = 0.5 * self.width, 0
        self.manipulators[0].set_pts([(-x, -y, 0), (x, -y, 0), (1, 0, 0)])
        self.manipulators[1].set_pts([(x, -y, self.altitude), (x, -y, self.altitude + self.height), (-1, 0, 0)])
        self.manipulators[2].set_pts([(x, -y, 0), (x, -y, self.altitude), (-1, 0, 0)])

        # always restore context
        self.restore_context(context)
        

class ARCHIPACK_PT_blind(Panel):
    bl_idname = "ARCHIPACK_PT_blind"
    bl_label = "Blind"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        # ensure your object panel only show when active object is the right one
        return archipack_blind.filter(context.active_object)

    def draw(self, context):
        o = context.active_object
        if not archipack_blind.filter(o):
            return
        layout = self.layout

        # retrieve datablock of your object
        props = archipack_blind.datablock(o)

        # Manipulate mode operator
        layout.operator('archipack.blind_manipulate', icon='HAND')

        box = layout.box()
        row = box.row(align=True)

        # Presets operators
        row.operator("archipack.blind_preset_menu",
                     text=bpy.types.ARCHIPACK_OT_blind_preset_menu.bl_label)
        row.operator("archipack.blind_preset",
                      text="",
                      icon='ZOOMIN')
        row.operator("archipack.blind_preset",
                      text="",
                      icon='ZOOMOUT').remove_active = True

        box = layout.box()
        box.prop(props, 'style')
        box.prop(props, 'width')
        box.prop(props, 'height')
        box.prop(props, 'altitude')
        box.prop(props, 'ratio')
        
        box = layout.box()
        box.label(text="Slat")
        
        if props.style == 'JAPANESE':
            box.prop(props, 'panels')
            
        else:
            if props.style != 'ROLLER':
                box.prop(props, 'depth')
            
            if props.style in {'ROLLER', 'BLADES'}:
                box.prop(props, 'blind_angle')
            
            elif props.style != 'PLATED':
                box.prop(props, 'angle')
            
        box = layout.box()
        box.label(text="Frame")
        box.prop(props, 'frame_height')
        box.prop(props, 'frame_depth')


class ARCHIPACK_OT_blind(ArchipackCreateTool, Operator):
    bl_idname = "archipack.blind"
    bl_label = "Blind"
    bl_description = "Create Blind"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    
    width = FloatProperty(
            name='Width',
            min=0.10, default=1, precision=3,
            description='Total width'
            ) 
    height = FloatProperty(
            name='Height',
            min=0.1, default=1.2, precision=3,
            description='Height'
            ) 
    altitude = FloatProperty(
            name='Altitude',
            min=0, default=0, precision=3,
            description='Altitude'
            )
    style = EnumProperty(
        items=(
            ('VENITIAN', 'Venitian', 'Venitian', 0),  # /  /  ///
            ('SLAT', 'Slat', 'Slat', 1), # / / /
            ('ROLLER', 'Roller', 'Roller', 2),  # __________
            ('BLADES', 'Blades', 'Blades', 3), # _ _ ___
            ('PLATED', 'Plated', 'Plated', 4),  # Vertical \/\/        -> regular
            # ('CELLULAR', 'Cellular blind', 'Cellular blind', 5),    # OOOO         -> venitian like variant
            ('JAPANESE', 'Japanese', 'Japanese', 6), # _-  -> horizontal
            ('VERTICAL_SLOTTED', 'Vertical slotted', 'Vertical slotted', 7)  # / / -> horizontal
        ),
        default='SLAT'
        )
    
    def create(self, context):

        # Create an empty mesh datablock
        m = bpy.data.meshes.new("Blind")

        # Create an object using the mesh datablock
        o = bpy.data.objects.new("Blind", m)

        # Add your properties on mesh datablock
        d = m.archipack_blind.add()
        
        d.width = self.width
        d.height = self.height
        d.altitude = self.altitude
        d.style = self.style
        
        # Link object into scene
        context.scene.objects.link(o)

        # select and make active
        o.select = True
        context.scene.objects.active = o

        # Load preset into datablock
        self.load_preset(d)

        # add a material
        self.add_material(o)
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = bpy.context.scene.cursor_location
            o.select = True
            context.scene.objects.active = o

            # Start manipulate mode
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_blind_preset_menu(PresetMenuOperator, Operator):
    bl_idname = "archipack.blind_preset_menu"
    bl_label = "Blind preset"
    preset_subdir = "archipack_blind"


class ARCHIPACK_OT_blind_preset(ArchipackPreset, Operator):
    """Add a Blind Preset"""
    bl_idname = "archipack.blind_preset"
    bl_label = "Add Blind preset"
    preset_menu = "ARCHIPACK_OT_blind_preset_menu"

    @property
    def blacklist(self):
        return ['manipulators']


class ARCHIPACK_OT_blind_manipulate(Operator):
    bl_idname = "archipack.blind_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        return archipack_blind.filter(context.active_object)

    def invoke(self, context, event):
        d = archipack_blind.datablock(context.active_object)
        d.manipulable_invoke(context)
        return {'FINISHED'}


def register():
    bpy.utils.register_class(archipack_blind)
    Mesh.archipack_blind = CollectionProperty(type=archipack_blind)
    bpy.utils.register_class(ARCHIPACK_PT_blind)
    bpy.utils.register_class(ARCHIPACK_OT_blind)
    bpy.utils.register_class(ARCHIPACK_OT_blind_preset_menu)
    bpy.utils.register_class(ARCHIPACK_OT_blind_preset)
    bpy.utils.register_class(ARCHIPACK_OT_blind_manipulate)


def unregister():
    bpy.utils.unregister_class(archipack_blind)
    del Mesh.archipack_blind
    bpy.utils.unregister_class(ARCHIPACK_PT_blind)
    bpy.utils.unregister_class(ARCHIPACK_OT_blind)
    bpy.utils.unregister_class(ARCHIPACK_OT_blind_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_OT_blind_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_blind_manipulate)
