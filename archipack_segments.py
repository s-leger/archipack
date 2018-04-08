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
from math import pi, atan2, cos, sin
from mathutils import Vector, Matrix
import bpy
from bpy.types import Operator
from bpy.props import (
    FloatProperty, IntProperty, BoolProperty,
    EnumProperty
    )
from .archipack_2d import Line, Arc
from .archipack_object import ArchipackGenericOperator


class Segment():

    def __init__(self):
        pass
    
    def set_offset(self, offset, last=None):
        """
            Offset line and compute intersection point
            between segments
        """
        if offset == 0:
            self.line = self
        else:
            self.line = self.make_offset(offset, last)


class StraightSegment(Segment, Line):

    def __init__(self, p, v):
        Line.__init__(self, p, v)
        Segment.__init__(self)
    
    

class CurvedSegment(Segment, Arc):

    def __init__(self, c, radius, a0, da):
        Arc.__init__(self, c, radius, a0, da)
        Segment.__init__(self)
    

        
class OpeningGenerator():
    """
     Provide relocate compatibility 
     to handle openings / customs
     the same way as other objects
    """
    def __init__(self, d, o=None):
        
        if o is None:
            # delta angle between local and world absolute
            tM = Matrix.Rotation(0, 3, Vector((0, 0, 1))).transposed()
            p = Vector((0, 0))
        else:
            if type(o).__name__ == "Matrix":
                tM = o.transposed()
                p = o.translation.to_2d()
            else:
                tM = o.matrix_world.transposed()
                p = o.matrix_world.translation.to_2d()
        y = -tM[1].to_2d() * d.y
        x = tM[0].to_2d() * 0.5 * d.x
        s0 = StraightSegment(p - 0.5 * y + x, y)
        s1 = StraightSegment(p - 0.5 * y - x, y)
        self.segs = [s0, s1]
        
    
class Generator():
    """
     Main class to inherit for generators
     provide:
      - io from datamodel
      - support for offset
      - update parts manipulators location
      - basic coordsys transformations
    """
    def __init__(self, d, o=None):
        """
         d: archipack datablock
         o: base object, or matrix
            when defined the generator is in world coordsys
            when None the generator is in local coordsys
        """
        self.d = d
        if o is None:
            # delta angle between local and world absolute
            self.rot = Matrix.Rotation(0, 3, Vector((0, 0, 1)))
            self.location = Vector((0, 0))
        else:
            if type(o).__name__ == "Matrix":
                tM = o
            else:
                tM = o.matrix_world

            loc, rot, scale = tM.decompose()
            self.rot = rot.to_matrix()
            self.location = loc.to_2d()

        self.segs = []

    @property
    def a0(self):
        """
         Return a0 for first segment in local coordsys
        """
        s = self.segs[0]
        if "C_" in self.parts[0].type:
            dp = s.p0 - s.c
            v = self.rot.inverted() * dp.to_3d()
        else:
            v = self.rot.inverted() * s.v.to_3d()
        return atan2(v.y, v.x)

    def next(self, idx, change_side):
        """
         Return next
         - index, part, segment in generator
        """
        next, p, s = None, None, None
        
        if idx is not None:
            if change_side == 'RIGHT':
                if idx + 1 < self.n_parts:
                    next = idx + 1
                elif self.closed:
                    next = 0
            else:    
                if idx - 1 > -1:
                    next = idx - 1
                elif self.closed:
                    next = -1
        
            if next is not None:
                p, s = self.parts[next], self.segs[next]
                
        return next, p, s

    @property
    def n_parts(self):
        n_parts = len(self.parts)
        # when open, parts use a hidden last segment
        # to handle manipulator at end point
        # when closed, the hidden part hold
        # closing segment properties
        if not self.closed:
            n_parts -= 1
        return n_parts

    @property
    def origin(self):
        return (self.rot * self.d.origin).to_2d()

    @property
    def closed(self):
        return self.d.closed

    @property
    def parts(self):
        return self.d.parts

    def _add_part(self, last, typ, length, a, radius, da):

        if 'S_' in typ:
            if last is None:
                p = self.location
                v = (self.rot * Vector((length, 0, 0))).to_2d()
                s = StraightSegment(p, v).rotate(a)
            else:
                seg = last.straight(length).rotate(a)
                s = StraightSegment(seg.p, seg.v)
        else:
            if last is None:
                c = self.location - (self.rot * (radius * Vector((cos(a), sin(a), 0)))).to_2d()
                s = CurvedSegment(c, radius, a, da)
            else:
                n = last.normal(1).rotate(a).scale(radius)
                if da < 0:
                    n.v = -n.v
                a0 = n.angle
                c = n.p - n.v
                s = CurvedSegment(c, radius, a0, da)

        self.segs.append(s)

    def add_part(self, part):

        last = None

        if len(self.segs) > 0:
            last = self.segs[-1]

        self._add_part(last, part.type, part.length, part.a0, part.radius, part.da)

        self.last_type = part.type

    def set_offset(self, offset):
        n_segs = self.n_parts
        last = None
        for i, seg in enumerate(self.segs):
            if i < n_segs:
                seg.set_offset(offset + self.parts[i].offset, last)
                last = seg.line
            else:
                seg.line = Line(last.p1, last.straight(0, 1).p1)

    def close(self, offset):
        # Make last segment implicit closing one
        if self.closed:

            _segs = self.segs
            _parts = self.parts

            part = _parts[-1]
            s0 = _segs[-1]
            s1 = _segs[0]
            dp = s1.p0 - s0.p0

            if "C_" in part.type:
                ds = (s0.p1 - s0.p0)
                s0.r = part.radius / ds.length * dp.length
                # angle pt - p0        - angle p0 p1
                da = atan2(dp.y, dp.x) - atan2(ds.y, ds.x)
                a0 = s0.a0 + da
                if a0 > pi:
                    a0 -= 2 * pi
                if a0 < -pi:
                    a0 += 2 * pi
                s0.a0 = a0
                # update angle and radius of closing segment to real values
                part.radius = s0.r
                part.a0 = s0.a0
            else:
                s0.v = dp
                # update angle and size of closing segment to real values
                part.length = s0.length
                if len(_segs) > 1:
                    part.a0 = s0.delta_angle(_segs[-2])

            if len(_segs) > 1:
                s0.line = s0.make_offset(offset + part.offset, _segs[-2].line)

            p1 = s1.line.p1
            s1.line = s1.make_offset(offset + _parts[0].offset, s0.line)
            s1.line.p1 = p1
        
    def locate_manipulators(self, side=1):
        """
            setup manipulators
        """
        _segs = self.segs
        _parts = self.parts
        n_parts = len(_parts)

        for i, f in enumerate(_segs):
            part = _parts[i]
            manipulators = part.manipulators
            p0 = f.p0.to_3d()
            p1 = f.p1.to_3d()
            # angle from last to current segment
            if i > 0:
                v0 = _segs[i - 1].straight(-side, 1).v.to_3d()
                v1 = f.straight(side, 0).v.to_3d()
                manipulators[0].type_key = 'DUAL_ANGLE'
                manipulators[0].prop1_name = "a_ui"
                manipulators[0].prop2_name = "change_side"
                manipulators[0].set_pts([p0, v0, v1])

            if 'S_' in part.type:
                # segment length
                manipulators[1].type_key = 'DUAL_SIZE'
                manipulators[1].prop1_name = "l_ui"
                manipulators[1].prop2_name = "change_side"
                manipulators[1].set_pts([p0, p1, (side, 0, 0)])
            else:
                # segment radius + angle
                v0 = (f.p0 - f.c).to_3d()
                v1 = (f.p1 - f.c).to_3d()
                manipulators[1].type_key = 'ARC_ANGLE_RADIUS'
                manipulators[1].prop1_name = "da_ui"
                manipulators[1].prop2_name = "r_ui"
                manipulators[1].set_pts([f.c.to_3d(), v0, v1])

            # snap manipulator, dont change index !
            manipulators[2].set_pts([p0, p1, (side, 0, 0)])
            # dumb segment id
            manipulators[3].set_pts([p0, p1, (side, 0, 0)])

            """
            # offset (currently not in use)
            manipulators[4].set_pts([
                p0,
                p0 + f.sized_normal(0, max(0.0001, self.parts[i].offset)).v.to_3d(),
                (0.5, 0, 0)
            ])
            """
            # Dimensions points
            self.d.add_dimension_point(part.uid, p0)
            if i == n_parts:
                self.d.add_dimension_point(part.uid + 1, p1)

    def draw(self, context):
        """
            draw generator using gl
        """
        for seg in self.segs:
            seg.draw(context, render=False)

    def get_verts(self, verts):
        verts.extend([s.p0.to_3d() for s in self.segs])


def update_type(self, context):

    i, o, d = self.find_datablock_in_selection(context)

    if d is not None and d.auto_update:

        d.auto_update = False
        g = Generator(d)
        _parts = d.parts
        a0 = 0
        w0 = None
        if i > 0:
            part0 = _parts[i - 1]
            g._add_part(None, part0.type, part0.length, part0.a0, part0.radius, part0.da)
            w0 = g.segs[0]
            a0 = w0.straight(1).angle

        if "C_" in self.type:
            typ = 'S_SEG'
        else:
            typ = 'C_SEG'

        g._add_part(w0, typ, self.length, self.a0, self.radius, self.da)
        w = g.segs[-1]

        # w0 - w - w1
        if d.closed and i == d.n_parts:
            dp = -w.p0
        else:
            dp = w.p1 - w.p0

        if "C_" in self.type:
            self.radius = 0.5 * dp.length
            self.dangle = pi
            a0 = atan2(dp.y, dp.x) - pi / 2 - a0
        else:
            self.length = dp.length
            a0 = atan2(dp.y, dp.x) - a0

        if a0 > pi:
            a0 -= 2 * pi
        if a0 < -pi:
            a0 += 2 * pi
        self.a0 = a0

        if i + 1 < d.n_parts:
            # adjust rotation of next part
            part1 = _parts[i + 1]
            if "C_" in self.type:
                a0 = part1.a0 - pi / 2
            else:
                a0 = part1.a0 + w.straight(1).angle - atan2(dp.y, dp.x)

            if a0 > pi:
                a0 -= 2 * pi
            if a0 < -pi:
                a0 += 2 * pi
            part1.a0 = a0

        d.auto_update = True


def update(self, context):
    self.update(context)


def update_angle(self, context):
    self.update_angle(context)


def set_a(self, value):
    self.last_a = self.a0
    self.a0 = value
    return None


def get_a(self):
    return self.a0


def update_length(self, context):
    self.update_length(context)


def get_l(self):
    return self.length


def set_l(self, value):
    self.last_l = self.length
    self.length = value
    return None


def update_radius(self, context):
    self.update_radius(context)


def get_r(self):
    return self.radius


def set_r(self, value):
    self.last_r = self.radius
    self.radius = value
    return None


def update_da(self, context):
    self.update_dangle(context)


def get_da(self):
    return self.da


def set_da(self, value):
    self.last_da = self.da
    self.da = value
    return None


class ArchipackSegment():
    """
        Base class for all archipack's line based objects parts
        wall, fence, slab, floor, cutter, roof
        Childs MUST implements
        - get_datablock(o)

        Use 3 vars for each param (2 unsaved one)
        where:
        - main store current value (saved)
        - last_ store previous value on change (not saved)
        - _ui is getter/setter and trigger auto update

        _ui is ment to be in use in ui and simple manipulators

        main dosent auto update so use when
        change does affect neighboors segments
        and call update by hand
    """
    type = EnumProperty(
            items=(
                ('S_SEG', 'Straight', '', 0),
                ('C_SEG', 'Curved', '', 1),
                ),
            default='S_SEG',
            update=update_type
            )
    
    change_side = EnumProperty(
        description="Side the changes apply, eg: size may change in left or right side of segment",
        items=(
            ('RIGHT', 'Right', 'Right'),
            ('LEFT', 'Left', 'Left')
            ),
        default='RIGHT',
        options={'SKIP_SAVE'}
        )

    # Segment length
    last_l = FloatProperty(
            description="Store last length between segments on update",
            default=2.0,
            unit='LENGTH', subtype='DISTANCE',
            options={'SKIP_SAVE'}
            )
    length = FloatProperty(
            description="Store segment length",
            min=0.001,
            default=2.0,
            unit='LENGTH', subtype='DISTANCE'
            )
    l_ui = FloatProperty(
            name="Length",
            description="Length of segment",
            min=0.001,
            default=2.0,
            set=set_l,
            get=get_l,
            update=update_length,
            unit='LENGTH', subtype='DISTANCE',
            options={'SKIP_SAVE'}
            )

    # Curved parts
    last_r = FloatProperty(
            description="Store last radius between segments on update",
            default=0.7,
            unit='LENGTH', subtype='DISTANCE',
            options={'SKIP_SAVE'}
            )
    radius = FloatProperty(
            description="Store segment radius",
            default=0.7,
            unit='LENGTH', subtype='DISTANCE'
            )
    r_ui = FloatProperty(
            name="Radius",
            description="Radius of curved segments",
            min=0.001,
            default=0.7,
            update=update_radius,
            get=get_r,
            set=set_r,
            unit='LENGTH', subtype='DISTANCE',
            options={'SKIP_SAVE'}
            )

    last_da = FloatProperty(
            description="Store last angle of curved  segment on update",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            options={'SKIP_SAVE'}
            )
    da = FloatProperty(
            description="Store angle of curved segments",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION'
            )
    da_ui = FloatProperty(
            name="Angle",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            set=set_da,
            get=get_da,
            update=update_da,
            options={'SKIP_SAVE'}
            )

    # Angle between segments
    last_a = FloatProperty(
            description="Store last angle between segments on update",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            options={'SKIP_SAVE'}
            )
    a0 = FloatProperty(
            description="Store angle between segments",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION'
            )
    a_ui = FloatProperty(
            name="Start angle",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            set=set_a,
            get=get_a,
            update=update_angle,
            options={'SKIP_SAVE'}
            )

    # Lateral offset
    offset = FloatProperty(
            name="Offset",
            description="Side offset of segment",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    auto_update = BoolProperty(
        default=True,
        options={'SKIP_SAVE'}
        )

    # DimensionProvider
    uid = IntProperty(default=0)

    def get_datablock(self, o):
        """
         Must return base object's datablock
        """
        raise NotImplementedError

    def find_datablock_in_selection(self, context):
        """
         Return selected object this instance belongs to
         provide support for "copy to selected"
        """
        selected = context.selected_objects[:]
        for o in selected:
            d = self.get_datablock(o)
            if d:
                for idx, part in enumerate(d.parts):
                    if part == self:
                        return idx, o, d
        return -1, None, None

    def generator(self, d, idx, changed):
        """
         Init a generator in last state and apply changes
         idx: changed segment
         changed: changed property to retrieve
        """
        # flag to prevent update when limit are reached
        res = True

        # Init generator in last state (state without current changes)
        g = Generator(d)
        last = None
        for i, part in enumerate(d.parts):
            if i == idx:
                if changed == "L":
                    g._add_part(last, self.type, self.last_l, self.a0, self.radius, self.da)
                elif changed == "A":
                    g._add_part(last, self.type, self.length, self.last_a, self.radius, self.da)
                elif changed == "R":
                    g._add_part(last, self.type, self.length, self.a0, self.last_r, self.da)
                elif changed == "DA":
                    g._add_part(last, self.type, self.length, self.a0, self.radius, self.last_da)
            else:
                g.add_part(part)
            last = g.segs[-1]

        # neighboors segments
        s0 = g.segs[idx]

        idx1, part1, s1 = g.next(idx, self.change_side)
        idx2, part2, s2 = g.next(idx1, self.change_side)

        # apply changes to segs
        if changed == "L":
            if self.change_side == 'RIGHT':
                s0.scale(self.length)
                if s1:
                    s1.p0 = s0.p1
            else:
                p1 = s0.p1.copy()
                s0.scale(self.length)
                s0.p = p1 - s0.v
                if s1:
                    s1.p1 = s0.p0
                
        elif changed == "A":
            if self.change_side == 'RIGHT':
                s0.rotate(self.a0 - self.last_a)
                if s1:
                    it, p, u, v = s0.intersect_ext(s1)
                    res = v < 1
                    s0.p1 = p
                    s1.p0 = p
            else:
                # rotate -1 part       
                if s1:
                    s1.rotate(self.last_a - self.a0)
                    s1.p = s0.p0 - s1.v
                if s2:
                    it, p, u, v = s1.intersect_ext(s2)
                    res = v > 0
                    s1.p0 = p
                    s2.p1 = p
                
        elif changed == "R":
            s0.scale(self.radius)
            if s1:
                s1.p0 = s0.p1

        elif changed == "DA":
            s0.da = self.da
            if s1:
                s1.p0 = s0.p1
  
        return g, res, s0, s1, s2, part1, part2, idx1, idx2

    def update_angle(self, context):
        """
         Update neighboor segments on angle change
        """
        i, o, d = self.find_datablock_in_selection(context)

        if d is not None:

            g, res, s0, s1, s2, part1, part2, idx1, idx2 = self.generator(d, i, "A")
            
            if self.change_side == 'RIGHT':
                if s1:
                    a = s1.delta_angle(s0)

                    if abs(a) < pi and res:
                        # when s1 is first segment: move object
                        if idx1 == 0:
                            # angle wont change here
                            d.move_object(o, o.matrix_world * s1.p0.to_3d())
                        else:
                            # compute start angle of next segment
                            # first part angle wont change
                            part1.a0 = a

                        if 'C_' in part1.type:
                            part1.radius = s1.r
                            part1.da = s1.da
                        else:
                            part1.length = s1.length

                        self.length = s0.length
                        # compute length of both segments
                        self.update(context)
                    else:
                        self.a0 = self.last_a
                else:
                    self.update(context)
            else:
                if s1:
                    # s2 - rotated s1 - s0
                    a = s0.delta_angle(s1)
                    print(a)
                    if abs(a) < pi and res:
                        self.a0 = a
                        
                        # when s1 is first segment: move object
                        if idx1 == 0:
                            part1.a0 = g.a0
                            d.move_object(o, o.matrix_world * s1.p0.to_3d())
                        
                        elif s2:
                            part1.a0 = s1.delta_angle(s2)

                        if 'C_' in part1.type:
                            part1.radius = s1.r
                            part1.da = s1.da
                        else:
                            part1.length = s1.length
                        
                        if s2:
                            part2.length = s2.length
                        
                        # compute length of both segments
                        self.update(context)
                    else:
                        self.a0 = self.last_a
                else:
                    self.update(context)

    def update_dangle(self, context):
        """
         Update neighboor segments on angle change
        """
        i, o, d = self.find_datablock_in_selection(context)

        if d is not None:
            g, res, s0, s1, s2, part1, part2, idx1, idx2 = self.generator(d, i, "DA")
            if s1:
                a = s1.delta_angle(s0)

                if abs(a) < pi:

                    # compute start angle of next segment
                    if idx1 == 0:
                        part1.a0 = s1.a0
                    else:
                        part1.a0 = a

                    if 'C_' in part1.type:
                        part1.radius = s1.r
                        part1.da = s1.da
                    else:
                        part1.length = s1.length

                    if s2:
                        if idx2 != 0:
                            part2.a = s2.delta_angle(s1)

                    # compute length of both segments
                    self.update(context)
                else:
                    self.da = self.last_da
            else:
                self.update(context)

    def update_length(self, context):
        """
         Update neighboor segments on length change
        """
        i, o, d = self.find_datablock_in_selection(context)

        if d is not None:

            g, res, s0, s1, s2, part1, part2, idx1, idx2 = self.generator(d, i, "L")
            if self.change_side == 'RIGHT':
                if s1:
                    if 'C_' in part1.type:
                        part1.radius = s1.r
                        part1.da = s1.da
                    else:
                        part1.length = s1.length

                    if idx1 == 0:
                        d.move_object(o, o.matrix_world * s1.p0.to_3d())
                        # seg angle in local object space
                        part1.a0 = s1.a0
                    else:
                        part1.a0 = s1.delta_angle(s0)

                if s2:
                    if idx2 != 0:
                        # angle for first seg dosen't change
                        part2.a0 = s2.delta_angle(s1)
            
            else:
                
                if i == 0:
                    d.move_object(o, o.matrix_world * s0.p0.to_3d())
                    self.a0 = g.a0
                    
                elif s1:
                    if 'C_' in part1.type:
                        part1.radius = s1.r
                        part1.da = s1.da
                    else:
                        part1.length = s1.length

                    if idx1 == 0:
                        # seg angle in local object space
                        part1.a0 = g.a0
                    
                    elif s2:
                        part1.a0 = s1.delta_angle(s2)
                        
                    self.a0 = s0.delta_angle(s1)

            self.update(context)

    def update_radius(self, context):
        """
         Update neighboor segments on length change
        """
        i, o, d = self.find_datablock_in_selection(context)

        if d is not None:

            g, res, s0, s1, s2, part1, part2, idx1, idx2 = self.generator(d, i, "R")

            if s1:
                if 'C_' in part1.type:
                    part1.radius = s1.r
                    part1.da = s1.da
                else:
                    part1.length = s1.length

                if idx1 == 0:
                    d.move_object(o, o.matrix_world * s1.p0.to_3d())
                    # seg angle in local object space
                    part1.a0 = s1.a0
                else:
                    part1.a0 = s1.delta_angle(s0)

            if s2:
                if idx2 != 0:
                    # angle for first seg dosen't change
                    part2.a0 = s2.delta_angle(s1)

            self.update(context)

    def update(self, context, manipulable_refresh=False):
        
        # Reset change side
        self.change_side = 'RIGHT'
        
        if self.auto_update:
            idx, o, d = self.find_datablock_in_selection(context)
            if d is not None:
                d.update(context, manipulable_refresh)

    def draw_insert(self, context, layout, index):
        """
            May implement draw for insert / remove segment operators
        """
        row = layout.row(align=True)
        row.operator("archipack.segment_insert", text="Split").index = index
        if index > 0:
            row.operator("archipack.segment_remove", text="Remove").index = index

    def draw(self, context, layout, index, draw_type=False):

        if draw_type:
            layout.prop(self, "type", text="Seg {}".format(str(index + 1)))
        else:
            layout.label(text="Seg {}".format(str(index + 1)))

        self.draw_insert(context, layout, index)

        if "C_" in self.type:
            layout.prop(self, "r_ui")
            layout.prop(self, "da_ui")
        else:
            layout.prop(self, "l_ui")
        layout.prop(self, "a_ui")


# ------------------------------------------------------------------
# Define operator class to manage parts
# ------------------------------------------------------------------


class ARCHIPACK_OT_segment_insert(ArchipackGenericOperator, Operator):
    bl_idname = "archipack.segment_insert"
    bl_label = "Insert"
    bl_description = "Insert segment splitting current in 2 equal parts"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    index = IntProperty(default=0)

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            d = self.datablock(o)
            if d is None:
                return {'CANCELLED'}
            d.insert_part(context, o, self.index)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_segment_remove(ArchipackGenericOperator, Operator):
    bl_idname = "archipack.segment_remove"
    bl_label = "Remove"
    bl_description = "Remove segment merging with previous one"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    index = IntProperty(default=0)

    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            d = self.datablock(o)
            if d is None:
                return {'CANCELLED'}
            d.remove_part(context, o, self.index)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


def register():
    bpy.utils.register_class(ARCHIPACK_OT_segment_insert)
    bpy.utils.register_class(ARCHIPACK_OT_segment_remove)


def unregister():
    bpy.utils.unregister_class(ARCHIPACK_OT_segment_insert)
    bpy.utils.unregister_class(ARCHIPACK_OT_segment_remove)
