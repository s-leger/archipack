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

import bpy
from bpy.types import Operator, PropertyGroup, Mesh, Curve, Panel
from bpy.props import (
    FloatProperty, CollectionProperty, StringProperty,
    BoolProperty, IntProperty, EnumProperty
    )
from mathutils import Vector, Matrix
from mathutils.geometry import interpolate_bezier
from math import cos, sin, pi, tan, atan2
from .archipack_2d import Line, Arc
from .archipack_manipulator import Manipulable, archipack_manipulator
from .archipack_preset import ArchipackPreset, PresetMenuOperator
from .archipack_object import ArchipackCreateTool, ArchipackObject
from .archipack_cutter import (
    CutAblePolygon, CutAbleGenerator,
    ArchipackCutter,
    ArchipackCutterPart
    )
from .archipack_gl import GlText
from .archipack_polylines import Io, Qtree, Envelope
import logging
logger = logging.getLogger("archipack")


class SeekPoint():
    """
     A location to find segments on tree
    """
    def __init__(self, p):
        self.p = p
        self.envelope = Envelope(p, p)

    def move(self, p):
        self.p = p
        self.envelope.initByPoints(p, p)


class SeekBox():
    def __init__(self):
        self.envelope = Envelope()

    def init_right(self, turtle):
        p0, v0 = turtle.front.p, turtle.right.v
        p1 = p0 + v0
        p2 = p0 + turtle.front.v
        p3 = p2 + v0
        self.envelope.initByPoints(p0, p1)
        self.envelope.expandToInclude(p2)
        self.envelope.expandToInclude(p3)

    def init_left(self, turtle, length):
        p0, v0 = turtle.front.p, turtle.left.v
        p1 = p0 + v0
        p2 = p0 + turtle.front.v * length
        p3 = p2 + v0
        self.envelope.initByPoints(p0, p1)
        self.envelope.expandToInclude(p2)
        self.envelope.expandToInclude(p3)

    def does_intersect(self, Q_segs, seg, skip):
        nb, found = Q_segs.intersects_ext(self, 0.001)
        u_min = 2
        v_min = 0
        i_idx = -1
        intersect = False
        for idx in found:
            it, pt, u, v = seg.intersect_ext(Q_segs._geoms[idx])
            # u is limited by tree
            if it and 1.25 >= u >= -0.25 and 1.0001 >= v >= -0.0001:
                intersect = True
                if u < u_min:
                    u_min = u
                    v_min = v
                    i_idx = idx
        return intersect, u_min, v_min, i_idx


class Seg():
    def __init__(self, p, v):
        # c0 c1 are Points
        self.p = p
        self.v = v
        self.envelope = Envelope(p, p + v)
        self.idx = -1

    def copy(self):
        return Seg(self.p.copy(), self.v.copy())

    def update_envelope(self):
        self.envelope.initByPoints(self.p, self.p + self.v)

    def t(self, other):
        """
            return param t of intersection on this seg
        """
        c = Vector((other.v.y, -other.v.x, 0))
        d = self.v * c
        if d == 0:
            # parallel
            return 0
        dp = other.p - self.p
        t = (c * dp) / d
        return t

    def intersect_ext(self, other):
        """
            intersect, return param t on both lines
        """
        c = Vector((other.v.y, -other.v.x, 0))
        d = self.v * c
        if d == 0:
            return False, self.p, 0, 0
        dp = other.p - self.p
        c2 = Vector((self.v.y, -self.v.x, 0))
        u = (c * dp) / d
        v = (c2 * dp) / d
        return True, self.p + self.v * u, u, v

    def point_on_seg(self, pt):
        dp = pt - self.p
        dl = self.v.length
        if dl == 0:
            return 0
        t = (self.v * dp) / (dl ** 2)
        return t

    def farest_point(self, other):
        """
         find param t on current seg
         and point on other seg
         of worst point of other seg
         return point on this segment
         and farest point on other one
        """
        p0 = other.p.copy()
        p1 = p0 + other.v
        t0 = self.point_on_seg(p0)
        t1 = self.point_on_seg(p1)
        if t1 > t0:
            p0 = p1
            t0 = t1
        return t0, p0

    def nearest_point(self, other):
        """
         find param t on current seg
         of closest point of other seg
        """
        p0 = other.p.copy()
        p1 = p0 + other.v
        t0 = self.point_on_seg(p0)
        t1 = self.point_on_seg(p1)
        if t1 < t0:
            p0 = p1
            t0 = t1
        return t0, p0

    def init(self, p, v):
        self.p = p
        self.v = v
        self.update_envelope()

    def output(self, context, gf, coordsys, name="Line"):
        p0 = self.p
        p1 = p0 + self.v
        if self.v.length > 0:
            line = gf.createLineString((p0, p1))
            Io.to_curve(context.scene, coordsys, line, name=name)

    def distance_pt(self, pt):
        dp = pt - self.p
        dl = self.v.length
        if dl == 0:
            return 0, 0
        d = (self.v.x * dp.y - self.v.y * dp.x) / dl
        t = (self.v * dp) / (dl ** 2)
        return d, t

    def minimal_dist(self, other):
        d0, t = self.distance_pt(other.p)
        d1, t = self.distance_pt(other.p + other.v)
        if d1 < d0:
            d0 = d1
        return d0


class Turtle():
    def __init__(self, p, v, cw=True):
        self.init_side(cw)
        side = v.cross(self.zAxis)
        self.p = p
        self.front = Seg(p, v)
        self.front2 = Seg(p, 2 * v)
        self.front3 = Seg(p, 3 * v)
        self.left = Seg(p, -side)
        self.left2 = Seg(p, 2 * -side)
        self.right = Seg(p, side)

    def init_side(self, cw):
        self.cw = cw
        if cw:
            self.up = 1
        else:
            self.up = -1
        self.zAxis = Vector((0, 0, self.up))

    def update_envelopes(self):
        self.front.update_envelope()
        self.front2.update_envelope()
        self.front3.update_envelope()
        self.left.update_envelope()
        self.left2.update_envelope()
        self.right.update_envelope()

    def rotate(self, v):
        """
          set absolute rotation
        """
        side = v.cross(self.zAxis)
        self.front.v = v
        self.front2.v = 2 * v
        self.front3.v = 3 * v
        self.left.v = -side
        self.left2.v = 2 * -side
        self.right.v = side
        self.update_envelopes()

    def move(self, v):
        self.p += v
        self.update_envelopes()

    def turn_right(self):
        self.rotate(self.right.v)

    def turn_left(self):
        self.rotate(self.left.v)

    def step_forward(self):
        self.move(self.front.v)

    def step_backward(self):
        self.move(-self.front.v)

    def reverse(self):
        self.init_side(not self.cw)
        self.rotate(self.front.v)

    def scale(self, factor):
        self.rotate(factor * self.front.v)

    def relocate(self, p, v):
        self.p += p - self.p
        self.rotate(v)


class Tree(Qtree):
    def __init__(self, coordsys):
        Qtree.__init__(self, coordsys)

    def newSegment(self, c0, c1):
        p0 = Vector((c0.coord.x, c0.coord.y, 0))
        p1 = Vector((c1.coord.x, c1.coord.y, 0))
        new_seg = Seg(p0, p1 - p0)
        self.insert(self.ngeoms, new_seg)
        return new_seg


        
class AbstractAnalyser():
    def __init__(self):
        self._result = None
        # segment to check for intersection
        self.seg = None
        # tree
        self.tree = None
    def reset(self):
        self._result = None
    def parse_kwargs(self, kwargs):
        keys = self.__dict__.keys()
        self.__dict__.update((key, value) for key, value in kwargs.items() if key in keys)
    def compute(self, **kwargs):
        if self._result is None:
            self.parse_kwargs(kwargs)
            self.analyse()
        
        
class AnalyseIntersection(AbstractAnalyser):
    """
      Find closest intersection in the tree
      for given seg
    """
    def __init__(self):
        AbstractAnalyser.__init__(self)        
    def _intersect(self):
        _seg = self.seg
        nb, found = self.tree.intersects(_seg)
        t = 1e32
        idx = -1
        intersect = False
        p0 = _seg.p0
        for i in found:
            f_seg = self.tree._geoms[i]
            # prevent intesection with touching segment
            if p0 != f_seg.p1:
                it, pt, u, v = _seg.intersect_ext(f_seg)                
                # u and v are limited by tree
                if it and u > 0 and 1.0001 >= v >= -0.0001:
                    if u < t:
                        intersect = True
                        t = u
                        idx = i
        return intersect, t, idx
    def analyse(self):
        # max param t for intersection on seek seg
        it, t, idx = self._intersect()
        self._result = it, t, idx
    def intersect(self, d):
        """
          return intersecting segment
          or None
        """
        seg = None
        it, t, idx = self.result
        if it:
            tmax = d / self.seg.length
            if t <= tmax:
                seg = self.tree._geoms[idx]
        return seg, t
        

class AnalyseObstacle(AbstractAnalyser):
    """
      Find closest cutter in the tree
      for given seg
    """
    def __init__(self):
        AbstractAnalyser.__init__(self)
    def _intersect(self):
        _seg = self.seg
        nb, found = self.tree.intersects(_seg)
        t = 1e32
        idx = -1
        intersect = False
        p0 = _seg.p0
        for i in found:
            f_seg = self.tree._geoms[i]
            # prevent intesection with touching segment
            if p0 != f_seg.p1:
                it, pt, u, v = _seg.intersect_ext(f_seg)                
                # u and v are limited by tree
                if it and u > 0 and 1.0001 >= v >= -0.0001:
                    if u < t:
                        intersect = True
                        t = u
                        idx = i
        return intersect, t, idx
    def analyse(self):
        # init SeekBox around given segment
        # find any cutter in that box
        # compute closest distance to cutter
        it, t, idx = self._intersect()
        self._result = it, t, idx
    def intersect(self, d):
        """
          return intersecting segment
          or None
        """
        seg = None
        it, t, idx = self.result
        if it:
            tmax = d / self.seg.length
            if t <= tmax:
                seg = self.tree._geoms[idx]
        return seg, t
        
"""    
a = AnalyseIntersection()
a.reset()
a.compute(seg=3, tree=1, d=2)      
seg, t = a.intersect(dist) 
"""       
        
class AbstractRule():    
    
    def __init__(self, apply_immediately):
        self.apply_immediately = apply_immediately
    
    def add_test(self, test):
        self.tests.append(test)
        
    def apply(self, tree, turtle):
        return 
    
    def run(self, tree, turtle):
        res = False
        for test in self.tests:
            res = test(tree, turtle)
            if res:
                self.apply(tree, turtle)
                break
        return res    
  
"""
Fast and reliabile way detect obstacles ?
Init a tree with cutters only !

Analysis:
- intersections in 3 directions
- obstacles
- stuck state 


Rule:
determine direction and size for next segment
-exception: apply immediately
-generic: let other rules try to apply


"""
     
class PathFinder():
    
    def __init__(self, tree, spacing, allow_backward, cw):
        self.tree = tree
        self.spacing = spacing
        self.segs = tree._geoms
        self.n_boundarys = tree.ngeoms
        # wall segment idx
        self.last = -1
        # helpers to seek in tree
        p = Vector((0, 0, 0))
        v = Vector((1, 0, 0))
        self.pt = SeekPoint(p)
        self.seg = Seg(p, v)
        self.left = Seg(p, v)
        self.right = Seg(p, v)
        self.turtle = Turtle(p, v * spacing, cw)
        # output coords
        self.coords = []
        # iteration current and max
        self.run = True
        self.allow_backward = allow_backward
        self.is_forward = True
        self.iter = 0
        self.max_iter = (tree.width * tree.height) / (spacing ** 2)

    def insert(self, p0, p1):
        seg = Seg(p0, p1 - p0)
        idx = self.tree.ngeoms
        seg.idx = idx
        self.tree.insert(idx, seg)
        return idx

    def intersect(self, seg):
        nb, found = self.tree.intersects_ext(seg, 0.005 * self.spacing)
        t = 1e32
        idx = -1
        intersect = False
        skip = self.tree.ngeoms - 1
        for i in found:
            if i != skip:
                it, pt, u, v = seg.intersect_ext(self.segs[i])
                # u is limited by tree
                if it and 1.0001 >= u > 0 and 1.0001 >= v >= -0.0001:
                    intersect = True
                    if u < t:
                        t = u
                        idx = i
        return intersect, t, idx

    def normal(self, seg):
        """
         return
         n : normal with same orientation or turtle front
         d : length of normal + turtle front normalized
        """
        n = seg.v.normalized().cross(self.turtle.zAxis)
        # direction must match turtle one
        # when both are in same direction length > 1
        v = self.turtle.front.v.normalized()
        d = (n + v).length
        if d < 1:
            n = -n
            d = (n + v).length
        return self.spacing * n, d

    def seg_at_pos(self, seg, p):
        """
         segment either start or end at location
        """
        v0 = seg.p - p
        return v0.length < 0.0001 or (v0 + seg.v).length < 0.0001

    def obstacle_seg(self, p):
        """
         find segment with largest normal at location
        """
        self.pt.move(p)
        d_max = -1
        idx = -1
        n = self.turtle.front.v
        nb, found = self.tree.intersects_ext(self.pt, 0.001)
        for i in found:
            seg = self.segs[i]
            if self.seg_at_pos(seg, p):
                ni, d = self.normal(seg)
                if d > d_max:
                    d_max = d
                    idx = i
                    n = ni
        return idx, n

    def wall_seg(self, p):
        """
         find segment != skip at location
        """
        self.pt.move(p)
        nb, found = self.tree.intersects_ext(self.pt, 0.001)
        for i in found:
            if i != self.last:
                seg = self.segs[i]
                if self.seg_at_pos(seg, p):
                    return i
        return -1

    def obstacle(self):
        """
         find obstacles for self.left
         return
         t param on turtle front
         n normal of hit segment
        """
        it, t, idx = self.intersect(self.left)
        if it:
            # find u of closest point on intersected vector
            seg = self.segs[idx]
            # self.left.output(context, gf, coordsys, name="Left_{}".format(self.iter))
            # seg.output(context, gf, coordsys, name="Hit_{}".format(self.iter))
            # self.turtle.front.output(context, gf, coordsys, name="Turtle_{}_front".format(self.iter))
            # self.turtle.left.output(context, gf, coordsys, name="Turtle_{}_left".format(self.iter))
            t, p = self.turtle.front.nearest_point(seg)
            next, n = self.obstacle_seg(p)
            next_seg = self.segs[next]
            # next_seg.output(context, gf, coordsys, name="next_seg_{}".format(self.iter))
            if next < self.n_boundarys:
                p -= 0.5 * n
            else:
                p -= n
            self.seg.init(p, next_seg.v)
            # self.seg.output(context, gf, coordsys, name="seek_seg_{}".format(self.iter))
            # print(next, p, next_seg.v)
            t = self.turtle.front.t(self.seg)
            return True, t, n, idx
        return False, 0, self.turtle.front.v, -1

    def start(self, p, v):
        self.turtle.relocate(p, v)
        self.coords.append(p.copy())

    def next(self):
        if self.is_forward:
            res = self.forward()
        else:
            if self.allow_backward:
                res = self.backward()
            else:
                return False
        if not res:
            self.is_forward = not self.is_forward
            logger.debug("%s forward:%s", self.iter, self.is_forward)
        if self.iter > self.max_iter:
            logger.debug("%s exit max iter", self.iter)
            return False
        return self.run

    def forward(self):
        self.iter += 1
        if self.iter > self.max_iter:
            logger.debug("%s exit max iter", self.iter)
            return False
        v = self.turtle.front.v
        # follow wall
        if self.last < 0:
            it, t, idx = self.intersect(self.turtle.right)
            logger.debug("%s hit right %s", self.iter, idx)
            self.last = idx
        
        #################
        # regular segment
        #################
        # estimate front t_max on wall segment
        # this is the size of seek area for
        # front and left intersections
        seg = self.segs[self.last]
        t, p0 = self.turtle.front.farest_point(seg)
        # store this one as last if nothing else hit
        next = self.wall_seg(p0)
        seg = self.segs[next]
        # identify the side: when d > 0 we are left side
        # multiply by turtle.up for reverse cases
        self.seg.init(p0, self.turtle.front.v)
        p1 = seg.p
        if (p0 - p1).length < 0.0001:
            p1 = seg.p + seg.v
        d, t1 = self.seg.distance_pt(p1)
        dir = 'RIGHT'
        side = 1
        if d * self.turtle.up > 0:
            dir = 'LEFT'
            side = -side
        # logger.debug("%s side:%s", self.iter, side)
        n, d = self.normal(seg)
        # logger.debug("%s last:%s next:%s", self.iter, self.last, next)
        # find intersection of turtle front vector
        # and segment parallel to seg
        # we know if the seg is
        # right (wall) or left (empty) side
        if next < self.n_boundarys:
            p = seg.p + side * 0.5 * n
        else:
            p = seg.p + side * n
        self.seg.init(p, seg.v)
        # t param for wall segment (if nothing else is hit)
        t = self.turtle.front.t(self.seg)

        if t == 0:
            print(p0)
            return False

        hit = False
        #################
        # obstacle on segment
        #################
        # segment parallel to seg
        # check for obstacle
        if not hit:
            # p1 = self.turtle.p + 0.5 * self.turtle.left.v
            p1 = self.turtle.p + self.turtle.left.v
            self.left.init(p1, v)
            # resize seek segment to max
            # parallel to next -+ 0.5 spacing
            self.seg.init(p + side * 0.5 * n, seg.v)
            t_max = self.left.t(self.seg)
            # logger.debug("%s t:%s t_max:%s", self.iter, t, t_max)
            if t_max < 0:
                t_max = 1 - t_max
            self.left.init(p1, t_max * v)
            it, o_t, o_n, idx = self.obstacle()
            # left hit seg
            if it and idx < self.n_boundarys and o_t < t:
                dir = 'LEFT'
                t, n = o_t, o_n
                next = idx
                logger.debug("%s hit left o %s t:%s", self.iter, idx, t)
                hit = True

        # check for pipe
        # dosent always work with such (1.5) space as it does hit |_|
        if not hit:
            p1 = self.turtle.p + 1.5 * self.turtle.left.v
            self.left.init(p1, v)
            # resize seek segment to max
            # parallel to next -+ spacing
            self.seg.init(p + side * 1.5 * n, seg.v)
            t_max = self.left.t(self.seg)
            # logger.debug("%s t:%s t_max:%s", self.iter, t, t_max)
            if t_max < 0:
                t_max = 1 - t_max
            self.left.init(p1, t_max * v)
            it, o_t, o_n, idx = self.obstacle()
            # left hit seg
            if it and idx >= self.n_boundarys and o_t < t:
                dir = 'LEFT'
                t, n = o_t, o_n
                next = idx
                logger.debug("%s hit left p %s t:%s", self.iter, idx, t)
                hit = True

        if next < self.n_boundarys:
            p = p0 - 0.5 * n
        else:
            p = p0 - n

        ###########################
        # forward past segment
        # as we might not go further
        # when we do left an obstacle
        ############################
        if not hit:
            # find intersection with wall seg line
            # resize seek segment to max + 2 to ensure there
            # is enougth space to come back
            self.seg.init(self.turtle.p, (t + 2.5) * v)
            # front hit seg
            it, t_f, idx = self.intersect(self.seg)
            if it:
                if idx != next:
                    # hit seg is not next one
                    next = idx
                    seg = self.segs[next]
                    n, d = self.normal(seg)
                    side = -1
                if next < self.n_boundarys:
                    p = seg.p + side * 0.5 * n
                else:
                    p = seg.p + side * n
                self.seg.init(p, seg.v)
                dir = 'LEFT'
                t = self.turtle.front.t(self.seg)
                logger.debug("%s hit front %s t:%s", self.iter, idx, t)
                hit = True

        ###########################
        # obstacle past segment
        # (is there enougth space to realy go right?)
        # When dir is right (wall side)
        # check for obstacle on right side
        # from t to t + 1  (1.5 for pipe ?)
        # if hit check for t < 0.5 bound t < 1 pipe
        #          __
        #    | _h_|o |
        #  __|    |__|
        #  ____t   n
        #     |0.5_|
        #     s
        # if hit, t is -1 / -0.5 from parallel to hit seg
        # 2 cases on hit
        # 1 obstacle seg is on same line as right one
        #   and obstacle is on right side of this line
        #   -> t becomes end obstacle
        # 2 dir becomes left next become obstacle one
        if not hit and dir == 'RIGHT':
            # build check seg on right side
            # start at intersection of parallel seg on turtle.right side
            # and next seg + 0.5
            # check for obstacle at 0.5 * right
            # to know if there is enougth space left
            self.seg.init(seg.p + 0.5 * n, seg.v)
            self.right.init(self.turtle.p + 0.5 * self.turtle.right.v, v)
            t0 = self.right.t(self.seg)
            self.right.init(self.right.p + t0 * v, 2 * v)
            it, o_t, idx = self.intersect(self.right)
            if it and idx < self.n_boundarys:
                # hit seg
                o_t, p = self.turtle.front.nearest_point(self.segs[idx])
                idx, o_n = self.obstacle_seg(p)
                # seg with greatest normal
                d = self.segs[idx].minimal_dist(seg)
                if d < 1.5 * self.spacing:
                    dir = 'LEFT'
                    seg = self.segs[idx]
                    next = idx
                    n = o_n
                    self.seg.init(seg.p - 0.5 * n, seg.v)
                    t = self.turtle.front.t(self.seg)
                    logger.debug("%s right obstacle %s t:%s", self.iter, idx, t)
                    hit = True

        if not hit and dir == 'RIGHT':
            # check for pipe at 1 * right
            # to know if there is enougth space left
            self.seg.init(seg.p + 0.5 * n, seg.v)
            self.right.init(self.turtle.p + self.turtle.right.v, v)
            t0 = self.right.t(self.seg)
            self.right.init(self.right.p + t0 * v, 2 * v)
            it, o_t, idx = self.intersect(self.right)
            if it and idx >= self.n_boundarys:
                # hit seg
                o_t, p = self.turtle.front.nearest_point(self.segs[idx])
                idx, o_n = self.obstacle_seg(p)
                # seg with greatest normal
                d = self.segs[idx].minimal_dist(seg)
                if d < 2 * self.spacing:
                    dir = 'LEFT'
                    seg = self.segs[idx]
                    next = idx
                    n = o_n
                    self.seg.init(seg.p - n, seg.v)
                    t = self.turtle.front.t(self.seg)
                    logger.debug("%s right pipe %s t:%s", self.iter, idx, t)
                    hit = True

        self.last = next
        logger.debug("%s t:%s", dir, t)
        if t < 0.5:
            # self.turtle.front.output(context, gf, coordsys, name="Turtle_{}".format(self.iter))
            self.turtle.rotate(self.segs[-1].v.normalized() * self.spacing)
            self.turtle.turn_right()
            self.turtle.move(0.5 * self.turtle.front.v)
            p0 = self.coords[-1]
            p1 = self.turtle.p.copy()
            self.coords.append(p1)
            self.insert(p0, p1)
            self.turtle.turn_right()
            it, t, idx = self.intersect(self.turtle.left)
            logger.debug("%s hit: last=%s", self.iter, idx)
            self.last = idx
            return False
        else:
            self.turtle.move(t * v)
            # rotate turtle to normal
            self.turtle.rotate(n)
            p0 = self.coords[-1]
            p1 = self.turtle.p.copy()
            self.coords.append(p1)
            self.insert(p0, p1)
            if dir == 'LEFT':
                self.turtle.turn_left()
            elif dir == 'RIGHT':
                self.turtle.turn_right()
        # self.turtle.front.output(context, gf, coordsys, name="Turtle_{}".format(self.iter))
        return True

    def backward(self):
        self.iter += 1
        backward = True
        v = self.turtle.front.v
        #################
        # wall found
        # follow til obstacle
        # or hit front
        #################
        # estimate front t_max on wall segment
        # this is the size of seek area for
        # front and left intersections
        seg = self.segs[self.last]
        t, p0 = self.turtle.front.farest_point(seg)
        # store this one as last if nothing else hit
        next = self.wall_seg(p0)
        dir = 'RIGHT'

        if next == -1:
            # end condition: found first segment
            self.run = False
            n, d = self.normal(seg)
        else:
            seg = self.segs[next]
            # identify the side: when d > 0 we are left side
            # multiply by turtle.up for reverse cases
            self.seg.init(p0, v)
            p1 = seg.p.copy()
            if (p0 - p1).length < 0.0001:
                p1 += seg.v
            d, t1 = self.seg.distance_pt(p1)
            side = -1
            if d * self.turtle.up > 0:
                dir = 'LEFT'
                side = -side
            n, d = self.normal(seg)
            # logger.debug("%s last:%s next:%s", self.iter, self.last, next)
            # find intersection of turtle front vector
            # and segment parallel to seg
            # so we have to know if the seg is
            # right (wall) or left (empty) side
            if next < self.n_boundarys:
                p = seg.p + side * 0.25 * n
            else:
                p = seg.p + side * 0.5 * n
            self.seg.init(p, seg.v)
            # t param for wall segment (if nothing else is hit)
            t = self.turtle.front.t(self.seg)
            logger.debug("%s side:%s, t:%s next:%s", self.iter, side, t, next)
            # found first segment

            # find any intersection along that segment
            # this occurs when angle < 90
            self.seg.init(self.turtle.p, v * t)
            it, o_t, idx = self.intersect(self.seg)
            if it:
                next = idx
                seg = self.segs[next]
                n, d = self.normal(seg)
                if next < self.n_boundarys:
                    p = seg.p - 0.25 * n
                else:
                    p = seg.p - 0.5 * n
                self.seg.init(p, seg.v)
                t = self.turtle.front.t(self.seg)
                logger.debug("%s hit along t:%s", self.iter, t)
                dir = 'RIGHT'
            else:
                #################
                # available
                # directions
                #################
                # front must hit unless we do have space
                p = self.turtle.p + v * t
                self.seg.init(p, v * 1.5)
                it, o_t, idx = self.intersect(self.seg)
                if it or True:
                    logger.debug("%s hit front", self.iter)
                    # other dir must hit
                    # check in next wall direction
                    # if there is space
                    # if hit something, use closest point on that segment - spacing as t
                    if dir == 'LEFT':
                        seg = self.turtle.left
                    else:
                        seg = self.turtle.right
                    self.seg.init(p, 0.75 * seg.v)
                    it, o_t, idx = self.intersect(self.seg)
                    if it:
                        seg = self.segs[idx]
                        # self.left.output(context, gf, coordsys, name="Left_{}".format(self.iter))
                        # seg.output(context, gf, coordsys, name="Hit_{}".format(self.iter))
                        # self.turtle.front.output(context, gf, coordsys, name="Turtle_{}_front".format(self.iter))
                        # self.turtle.left.output(context, gf, coordsys, name="Turtle_{}_left".format(self.iter))
                        t, p = self.turtle.front.nearest_point(seg)
                        next, n = self.obstacle_seg(p)
                        next_seg = self.segs[next]
                        # next_seg.output(context, gf, coordsys, name="next_seg_{}".format(self.iter))
                        # seg dosent cross -> -n
                        # seg cross +n
                        self.seg.init(self.turtle.p, 2 * seg.v.length * -n)
                        o_t = next_seg.t(self.seg)
                        side = 1
                        if 0 < o_t < 1:
                            side = -side
                        if next < self.n_boundarys:
                            p -= side * 0.25 * n
                        else:
                            p -= side * 0.5 * n
                        self.seg.init(p, next_seg.v)
                        n = -n
                        # self.seg.output(context, gf, coordsys, name="seek_seg_{}".format(self.iter))
                        # print(next, p, next_seg.v)
                        t = self.turtle.front.t(self.seg)
                        logger.debug("%s hit side t:%s", self.iter, t)

                    if abs(t) < 0.5:
                        # self.turtle.front.output(context, gf, coordsys, name="Turtle_{}".format(self.iter))
                        self.run = False
                        return False

                else:
                    # space in front, could we run forward ?
                    backward = False
        self.last = next
        self.turtle.move(t * v)
        # rotate turtle to normal
        self.turtle.rotate(n)
        p0 = self.coords[-1]
        p1 = self.turtle.p.copy()
        self.coords.append(p1)
        self.insert(p0, p1)
        if dir == 'LEFT':
            self.turtle.turn_left()
        elif dir == 'RIGHT':
            self.turtle.turn_right()
        # self.turtle.front.output(context, gf, coordsys, name="Turtle_{}".format(self.iter))
        return backward


class Floor():

    def __init__(self):
        # self.colour_inactive = (1, 1, 1, 1)
        pass

    def set_offset(self, offset, last=None):
        """
            Offset line and compute intersection point
            between segments
        """
        self.line = self.make_offset(offset, last)

    def straight_floor_heating(self, a0, length):
        s = self.straight(length).rotate(a0)
        return StraightFloor(s.p, s.v)

    def curved_floor_heating(self, a0, da, radius):
        n = self.normal(1).rotate(a0).scale(radius)
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedFloor(c, radius, a0, da)


class StraightFloor(Floor, Line):

    def __init__(self, p, v):
        Line.__init__(self, p, v)
        Floor.__init__(self)


class CurvedFloor(Floor, Arc):

    def __init__(self, c, radius, a0, da):
        Arc.__init__(self, c, radius, a0, da)
        Floor.__init__(self)


class FloorGenerator(CutAblePolygon, CutAbleGenerator):

    def __init__(self, parts):
        self.parts = parts
        self.segs = []
        self.holes = []
        self.convex = True
        self.xsize = 0
        # return values
        self.iter = 0
        self.pipe_len = 0
        self.area = 0

    def add_part(self, part):

        if len(self.segs) < 1:
            s = None
        else:
            s = self.segs[-1]
        # start a new floor
        if s is None:
            if part.type == 'S_SEG':
                p = Vector((0, 0))
                v = part.length * Vector((cos(part.a0), sin(part.a0)))
                s = StraightFloor(p, v)
            elif part.type == 'C_SEG':
                c = -part.radius * Vector((cos(part.a0), sin(part.a0)))
                s = CurvedFloor(c, part.radius, part.a0, part.da)
        else:
            if part.type == 'S_SEG':
                s = s.straight_floor_heating(part.a0, part.length)
            elif part.type == 'C_SEG':
                s = s.curved_floor_heating(part.a0, part.da, part.radius)

        self.segs.append(s)
        self.last_type = part.type

    def set_offset(self, offset):
        last = None
        for i, seg in enumerate(self.segs):
            seg.set_offset(offset + self.parts[i].offset, last)
            last = seg.line

    def close(self, offset):
        # Make last segment implicit closing one
        part = self.parts[-1]
        w = self.segs[-1]
        dp = self.segs[0].p0 - self.segs[-1].p0
        if "C_" in part.type:
            dw = (w.p1 - w.p0)
            w.r = part.radius / dw.length * dp.length
            # angle pt - p0        - angle p0 p1
            da = atan2(dp.y, dp.x) - atan2(dw.y, dw.x)
            a0 = w.a0 + da
            if a0 > pi:
                a0 -= 2 * pi
            if a0 < -pi:
                a0 += 2 * pi
            w.a0 = a0
        else:
            w.v = dp

        if len(self.segs) > 1:
            w.line = w.make_offset(offset + self.parts[-1].offset, self.segs[-2].line)

        p1 = self.segs[0].line.p1
        self.segs[0].line = self.segs[0].make_offset(offset + self.parts[0].offset, w.line)
        self.segs[0].line.p1 = p1

    def locate_manipulators(self):
        """
            setup manipulators
        """
        for i, f in enumerate(self.segs):

            manipulators = self.parts[i].manipulators
            p0 = f.p0.to_3d()
            p1 = f.p1.to_3d()
            # angle from last to current segment
            if i > 0:
                v0 = self.segs[i - 1].straight(-1, 1).v.to_3d()
                v1 = f.straight(1, 0).v.to_3d()
                manipulators[0].set_pts([p0, v0, v1])

            if type(f).__name__ == "StraightFloor":
                # segment length
                manipulators[1].type_key = 'SIZE'
                manipulators[1].prop1_name = "length"
                manipulators[1].set_pts([p0, p1, (1, 0, 0)])
            else:
                # segment radius + angle
                v0 = (f.p0 - f.c).to_3d()
                v1 = (f.p1 - f.c).to_3d()
                manipulators[1].type_key = 'ARC_ANGLE_RADIUS'
                manipulators[1].prop1_name = "da"
                manipulators[1].prop2_name = "radius"
                manipulators[1].set_pts([f.c.to_3d(), v0, v1])

            # snap manipulator, dont change index !
            manipulators[2].set_pts([p0, p1, (1, 0, 0)])
            # dumb segment id
            manipulators[3].set_pts([p0, p1, (1, 0, 0)])

    def get_verts(self, verts):
        for s in self.segs:
            if "Curved" in type(s).__name__:
                for i in range(16):
                    # x, y = floor.line.lerp(i / 16)
                    verts.append(s.lerp(i / 16).to_3d())
            else:
                # x, y = s.line.p0
                verts.append(s.p0.to_3d())
            """
            for i in range(33):
                x, y = floor.line.lerp(i / 32)
                verts.append((x, y, 0))
            """

    def rotate(self, idx_from, a):
        """
            apply rotation to all following segs
        """
        self.segs[idx_from].rotate(a)
        ca = cos(a)
        sa = sin(a)
        rM = Matrix([
            [ca, -sa],
            [sa, ca]
            ])
        # rotation center
        p0 = self.segs[idx_from].p0
        for i in range(idx_from + 1, len(self.segs)):
            seg = self.segs[i]
            # rotate seg
            seg.rotate(a)
            # rotate delta from rotation center to segment start
            dp = rM * (seg.p0 - p0)
            seg.translate(dp)

    def translate(self, idx_from, dp):
        """
            apply translation to all following segs
        """
        self.segs[idx_from].p1 += dp
        for i in range(idx_from + 1, len(self.segs)):
            self.segs[i].translate(dp)

    def draw(self, context):
        """
            draw generator using gl
        """
        for seg in self.segs:
            seg.draw(context, render=False)

    def limits(self):
        x_size = [s.p0.x for s in self.segs]
        y_size = [s.p0.y for s in self.segs]
        for s in self.segs:
            if "Curved" in type(s).__name__:
                x_size.append(s.c.x + s.r)
                x_size.append(s.c.x - s.r)
                y_size.append(s.c.y + s.r)
                y_size.append(s.c.y - s.r)

        self.xmin = min(x_size)
        self.xmax = max(x_size)
        self.xsize = self.xmax - self.xmin
        self.ymin = min(y_size)
        self.ymax = max(y_size)
        self.ysize = self.ymax - self.ymin

    def cut(self, context, o):
        """
            either external or holes cuts
        """
        self.limits()
        self.as_lines()
        self.is_convex()
        for b in o.children:
            d = archipack_floor_heating_cutter.datablock(b)
            if d is not None:
                g = d.ensure_direction()
                g.change_coordsys(b.matrix_world, o.matrix_world)
                self.slice(g)

    def _add_spline(self, curve, closed, coords):
        spline = curve.splines.new('POLY')
        spline.use_endpoint_u = False
        spline.use_cyclic_u = closed
        spline.points.add(len(coords) - 1)
        for i, coord in enumerate(coords):
            x, y, z = coord
            spline.points[i].co = (x, y, z, 1)

    def bevel(self, coords, radius, a):
        n_coords = len(coords) - 1
        verts = [coords[0]]
        for i, co in enumerate(coords):
            if i > 0 and i < n_coords:
                self.roundedCorner(co, coords[i - 1], coords[i + 1], radius, a, verts)
        verts.append(coords[-1])
        return verts

    def roundedCorner(self, p, p1, p2, radius, a, verts):
        # Vector 1
        u = p1 - p

        # Vector 2
        v = p2 - p

        # Angle between vector 1 and vector 2 divided by 2
        angle = atan2(u.x * v.y - u.y * v.x, u.x * v.x + u.y * v.y) / 2
        # The length of segment between angular point and the
        # points of intersection with the circle of a given radius
        tang = tan(abs(angle))
        segment = radius / tang

        # Check the segment
        length1 = u.length
        length2 = v.length
        length = 0.5 * min(length1, length2)

        if segment > length:
            segment = length
            radius = length * tang

        # Points of intersection are calculated by the proportion between
        # the coordinates of the vector, length of vector and the length of the segment.
        uc = (segment / length1) * u
        vc = (segment / length2) * v
        # Calculation of the coordinates of the circle
        # center by the addition of angular vectors.
        d = (segment ** 2 + radius ** 2) ** 0.5
        c = (uc + vc).normalized() * d
        n1 = uc - c
        n2 = vc - c
        center = p + c
        # StartAngle and EndAngle of arc
        startAngle = atan2(n1.y, n1.x)
        endAngle = atan2(n2.y, n2.x)

        # Sweep angle
        sweepAngle = endAngle - startAngle
        if sweepAngle > pi:
            sweepAngle -= 2 * pi

        if sweepAngle < -pi:
            sweepAngle += 2 * pi

        steps = int(abs(sweepAngle) / a)

        a0 = startAngle
        if steps == 0:
            da = 0
        else:
            da = sweepAngle / steps

        verts.extend([center + radius * Vector((cos(a0 + da * i), sin(a0 + da * i), 0)) for i in range(steps + 1)])

    def floor_heating(self, context, o, d):
        curve = o.data
        curve.splines.clear()
        verts = []
        self.get_verts(verts)
        self._add_spline(curve, True, verts)
        for hole in self.holes:
            verts = [s.p0.to_3d() for s in hole.segs if s.length > 0]
            if len(verts) > 0:
                self._add_spline(curve, True, verts)

        context.scene.update()

        child = None
        for c in o.children:
            if c.type == 'CURVE':
                child = c
                pipes = child.data

        if child is None:
            pipes = bpy.data.curves.new("Pipes", type='CURVE')
            pipes.dimensions = '3D'
            child = bpy.data.objects.new("Pipes", pipes)
            context.scene.objects.link(child)
            child.parent = o
            child.matrix_world = o.matrix_world.copy()

        pipes.splines.clear()
        curves = [o]
        resolution = 1
        # gf = GeometryFactory()
        # init boundarys
        coordsys = Io.getCoordsys(curves)
        tree = Tree(coordsys)
        Q_points = Qtree(coordsys)
        Io.add_curves(Q_points, tree, coordsys, curves, resolution)
        logger.debug("n geoms:%s", tree.ngeoms)
        cw = d.pattern == 'CW'
        start = 0
        center = Vector((0.5 * tree.width, 0.5 * tree.height, 0))
        for s in self.segs:
            length = s.length
            if start + length < d.start_location:
                start += length
            else:
                t = (d.start_location - start) / length
                t_min = 0.5 * d.spacing / length
                n = s.sized_normal(min(max(t_min, t), 1 - t_min), -0.5 * d.spacing)
                p = (n.p + n.v).to_3d() - center
                if cw:
                    v = d.spacing * s.v.normalized().to_3d()
                else:
                    v = -d.spacing * s.v.normalized().to_3d()
        
        print(p, v)      
        
        #Vector((-0.5 * tree.width, -0.5 * tree.height, 0))
        # p += d.spacing * Vector((0.5, 0.5, 0))
        # if cw:
        #    v = Vector((d.spacing, 0, 0))
        # else:
        #    v = Vector((0, d.spacing, 0))

        pf = PathFinder(tree, d.spacing, d.backward, cw)
        pf.start(p, v)
        pf.max_iter = d.max_iter
        run = True
        while run:
            run = pf.next()
        if d.enable_radius:
            a = 12 / 180 * pi
            coords = self.bevel(pf.coords, d.radius, a)
        else:
            coords = pf.coords
        self._add_spline(pipes, False, coords)
        child.location = center
        self.iter = pf.iter
        self.pipe_len = sum([(coords[i - 1] - co).length for i, co in enumerate(coords) if i > 0])
        self.area = 0

    def add_manipulator(self, name, pt1, pt2, pt3):
        m = self.manipulators.add()
        m.prop1_name = name
        m.set_pts([pt1, pt2, pt3])


def update(self, context):
    self.update(context)


def update_type(self, context):

    d = self.find_in_selection(context)

    if d is not None and d.auto_update:

        d.auto_update = False
        # find part index
        idx = 0
        for i, part in enumerate(d.parts):
            if part == self:
                idx = i
                break

        part = d.parts[idx]
        a0 = 0
        if idx > 0:
            g = d.get_generator()
            w0 = g.segs[idx - 1]
            a0 = w0.straight(1).angle
            if "C_" in self.type:
                w = w0.straight_floor_heating(part.a0, part.length)
            else:
                w = w0.curved_floor_heating(part.a0, part.da, part.radius)
        else:
            if "C_" in self.type:
                p = Vector((0, 0))
                v = self.length * Vector((cos(self.a0), sin(self.a0)))
                w = StraightFloor(p, v)
                a0 = pi / 2
            else:
                c = -self.radius * Vector((cos(self.a0), sin(self.a0)))
                w = CurvedFloor(c, self.radius, self.a0, pi)

        # w0 - w - w1
        if idx + 1 == d.n_parts:
            dp = - w.p0
        else:
            dp = w.p1 - w.p0

        if "C_" in self.type:
            part.radius = 0.5 * dp.length
            part.da = pi
            a0 = atan2(dp.y, dp.x) - pi / 2 - a0
        else:
            part.length = dp.length
            a0 = atan2(dp.y, dp.x) - a0

        if a0 > pi:
            a0 -= 2 * pi
        if a0 < -pi:
            a0 += 2 * pi
        part.a0 = a0

        if idx + 1 < d.n_parts:
            # adjust rotation of next part
            part1 = d.parts[idx + 1]
            if "C_" in part.type:
                a0 = part1.a0 - pi / 2
            else:
                a0 = part1.a0 + w.straight(1).angle - atan2(dp.y, dp.x)

            if a0 > pi:
                a0 -= 2 * pi
            if a0 < -pi:
                a0 += 2 * pi
            part1.a0 = a0

        d.auto_update = True


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


def update_path(self, context):
    self.update_path(context)


class archipack_floor_heating_part(PropertyGroup):

    """
        A single manipulable polyline like segment
        polyline like segment line or arc based
    """
    type = EnumProperty(
            items=(
                ('S_SEG', 'Straight', '', 0),
                ('C_SEG', 'Curved', '', 1),
                ),
            default='S_SEG',
            update=update_type
            )
    length = FloatProperty(
            name="Length",
            min=0.01,
            default=2.0,
            update=update
            )
    radius = FloatProperty(
            name="Radius",
            min=0.5,
            default=0.7,
            update=update
            )
    da = FloatProperty(
            name="Angle",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    a0 = FloatProperty(
            name="Start angle",
            min=-2 * pi,
            max=2 * pi,
            default=0,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    offset = FloatProperty(
            name="Offset",
            description="Side offset of segment",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    manipulators = CollectionProperty(type=archipack_manipulator)

    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = archipack_floor_heating.datablock(o)
            if props:
                for part in props.parts:
                    if part == self:
                        return props
        return None

    def update(self, context, manipulable_refresh=False):
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context, manipulable_refresh)

    def draw(self, context, layout, index):
        box = layout.box()
        # box.prop(self, "type", text=str(index + 1))
        box.label(text="#" + str(index + 1))
        if self.type in ['C_SEG']:
            box.prop(self, "radius")
            box.prop(self, "da")
        else:
            box.prop(self, "length")
        box.prop(self, "a0")


class archipack_floor_heating(ArchipackObject, Manipulable, PropertyGroup):
    n_parts = IntProperty(
            name="Parts",
            min=1,
            default=1, update=update_manipulators
            )
    parts = CollectionProperty(type=archipack_floor_heating_part)
    user_defined_path = StringProperty(
            name="User defined",
            update=update_path
            )
    user_defined_resolution = IntProperty(
            name="Resolution",
            min=1,
            max=128,
            default=12, update=update_path
            )

    # UI layout related
    parts_expand = BoolProperty(
            options={'SKIP_SAVE'},
            default=False
            )
    max_iter = IntProperty(
            name='Max iter',
            description='Maximum iteration',
            min=10,
            default=200,
            update=update
            )
    pattern = EnumProperty(
            name='Floor Pattern',
            items=(("CW", "CW", "regular up = 1"),
                    ("CCW", "CCW", "reverse up = -1 left right")),
            default="CW",
            update=update
            )
    spacing = FloatProperty(
            name='Spacing',
            description='The amount of space between pipes',
            unit='LENGTH', subtype='DISTANCE',
            min=0.01,
            default=0.3,
            precision=2,
            update=update
            )
    enable_radius = BoolProperty(
            name="Enable radius",
            default=True,
            update=update
            )
    radius = FloatProperty(
            name='radius',
            description='Radius',
            unit='LENGTH', subtype='DISTANCE',
            min=0.01,
            default=0.1,
            precision=2,
            update=update
            )
    backward = BoolProperty(
            name="backward",
            default=True,
            update=update
            )

    pipe_len = StringProperty()
    iter = StringProperty()
    area = StringProperty()
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )
    closed = BoolProperty(
            options={'SKIP_SAVE'},
            default=True
            )
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
    start_location = FloatProperty(
            name="Start location",
            default=0.1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
            
    def get_generator(self):
        g = FloorGenerator(self.parts)
        for part in self.parts:
            # type, radius, da, length
            g.add_part(part)

        g.set_offset(self.x_offset)

        g.close(self.x_offset)
        g.locate_manipulators()
        for i, seg in enumerate(g.segs):
            g.segs[i] = seg.line

        return g

    def update_parts(self, o):

        for i in range(len(self.parts), self.n_parts, -1):
            self.parts.remove(i - 1)

        # add rows
        for i in range(len(self.parts), self.n_parts):
            self.parts.add()

        self.setup_manipulators()

        g = self.get_generator()

        return g

    def is_cw(self, pts):
        p0 = pts[0]
        d = 0
        for p in pts[1:]:
            d += (p.x * p0.y - p.y * p0.x)
            p0 = p
        return d > 0

    def interpolate_bezier(self, pts, wM, p0, p1, resolution):
        # straight segment, worth testing here
        # since this can lower points count by a resolution factor
        # use normalized to handle non linear t
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

    def from_spline(self, context, wM, resolution, spline):
        pts = []
        if spline.type == 'POLY':
            pts = [wM * p.co.to_3d() for p in spline.points]
            if spline.use_cyclic_u:
                pts.append(pts[0])
        elif spline.type == 'BEZIER':
            points = spline.bezier_points
            for i in range(1, len(points)):
                p0 = points[i - 1]
                p1 = points[i]
                self.interpolate_bezier(pts, wM, p0, p1, resolution)
            if spline.use_cyclic_u:
                p0 = points[-1]
                p1 = points[0]
                self.interpolate_bezier(pts, wM, p0, p1, resolution)
                pts.append(pts[0])
            else:
                pts.append(wM * points[-1].co)

        # pretranslate
        o = self.find_in_selection(context, self.auto_update)
        o.matrix_world = self.from_points(pts)

    def from_points(self, pts):

        if self.is_cw(pts):
            pts = list(reversed(pts))

        self.auto_update = False

        self.n_parts = len(pts) - 1

        self.update_parts(None)

        p0 = pts.pop(0)
        tM = Matrix.Translation(p0.copy())
        a0 = 0
        for i, p1 in enumerate(pts):
            dp = p1 - p0
            da = atan2(dp.y, dp.x) - a0
            if da > pi:
                da -= 2 * pi
            if da < -pi:
                da += 2 * pi
            if i >= len(self.parts):
                break
            p = self.parts[i]
            p.length = dp.to_2d().length
            p.dz = dp.z
            p.a0 = da
            a0 += da
            p0 = p1

        self.closed = True
        self.auto_update = True
        return tM
        
    def update_path(self, context):
        user_def_path = context.scene.objects.get(self.user_defined_path)
        if user_def_path is not None and user_def_path.type == 'CURVE':
            self.from_spline(
                context,
                user_def_path.matrix_world,
                self.user_defined_resolution,
                user_def_path.data.splines[0])

    def add_manipulator(self, name, pt1, pt2, pt3):
        m = self.manipulators.add()
        m.prop1_name = name
        m.set_pts([pt1, pt2, pt3])

    def update_manipulators(self):
        self.manipulators.clear()  # clear every time, add new ones
        self.add_manipulator("length", (0, 0, 0), (0, self.length, 0), (-0.4, 0, 0))
        self.add_manipulator("width", (0, 0, 0), (self.width, 0, 0), (0.4, 0, 0))

    def setup_manipulators(self):

        if len(self.manipulators) < 1:
            s = self.manipulators.add()
            s.type_key = "SIZE"
            s.prop1_name = "z"
            s.normal = Vector((0, 1, 0))

        for i in range(self.n_parts):
            p = self.parts[i]
            n_manips = len(p.manipulators)
            if n_manips < 1:
                s = p.manipulators.add()
                s.type_key = "ANGLE"
                s.prop1_name = "a0"
            p.manipulators[0].type_key = 'ANGLE'
            if n_manips < 2:
                s = p.manipulators.add()
                s.type_key = "SIZE"
                s.prop1_name = "length"
            if n_manips < 3:
                s = p.manipulators.add()
                s.type_key = 'WALL_SNAP'
                s.prop1_name = str(i)
                s.prop2_name = 'z'
            if n_manips < 4:
                s = p.manipulators.add()
                s.type_key = 'DUMB_STRING'
                s.prop1_name = str(i + 1)
            p.manipulators[2].prop1_name = str(i)
            p.manipulators[3].prop1_name = str(i + 1)

        self.parts[-1].manipulators[0].type_key = 'DUMB_ANGLE'

    def text(self, context, value, type, precision=2):

        dimension = 1

        if type == 'AREA':
            dimension = 2
        unit_type = 'SIZE'
        unit_mode = 'AUTO'

        label = GlText(
            label="",
            value=value,
            precision=precision,
            unit_mode=unit_mode,
            unit_type=unit_type,
            dimension=dimension
            )
        return label.add_units(context)

    def update(self, context, manipulable_refresh=False):

        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # clean up manipulators before any data model change
        if manipulable_refresh:
            self.manipulable_disable(context)

        g = self.update_parts(o)

        g.cut(context, o)
        g.floor_heating(context, o, self)

        # enable manipulators rebuild
        if manipulable_refresh:
            self.manipulable_refresh = True

        self.iter = "Iter {}".format(g.iter)
        self.pipe_len = "Len {}".format(self.text(context, g.pipe_len, 'SIZE', 2))
        self.area = "Area {}".format(self.text(context, g.area, 'AREA', 2))

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

        self.manipulable_setup(context)
        self.manipulate_mode = True

        self._manipulable_invoke(context)

        return True


def update_hole(self, context):
    self.update(context, update_parent=True)


def update_operation(self, context):
    self.reverse(context, make_ccw=(self.operation == 'INTERSECTION'))


class archipack_floor_heating_cutter_segment(ArchipackCutterPart, PropertyGroup):
    manipulators = CollectionProperty(type=archipack_manipulator)
    type = EnumProperty(
        name="Type",
        items=(
            ('DEFAULT', 'Side', 'Side with rake', 0),
            ),
        default='DEFAULT',
        update=update_hole
        )

    def find_in_selection(self, context):
        selected = [o for o in context.selected_objects]
        for o in selected:
            d = archipack_floor_heating_cutter.datablock(o)
            if d:
                for part in d.parts:
                    if part == self:
                        return d
        return None

    def draw(self, layout, context, index):
        box = layout.box()
        box.label(text="Part:" + str(index + 1))
        # box.prop(self, "type", text=str(index + 1))
        box.prop(self, "length")
        box.prop(self, "a0")


class archipack_floor_heating_cutter(ArchipackCutter, ArchipackObject, Manipulable, PropertyGroup):
    parts = CollectionProperty(type=archipack_floor_heating_cutter_segment)

    def update_points(self, context, o, pts, update_parent=False):
        """
            Create boundary from roof
        """
        self.auto_update = False
        self.from_points(pts)
        self.auto_update = True
        if update_parent:
            self.update_parent(context, o)

    def update_parent(self, context, o):

        d = archipack_floor_heating.datablock(o.parent)
        if d is not None:
            o.parent.select = True
            context.scene.objects.active = o.parent
            d.update(context)
        o.parent.select = False
        context.scene.objects.active = o


# ------------------------------------------------------------------
# Define panel class to show object parameters in ui panel (N)
# ------------------------------------------------------------------


class ARCHIPACK_PT_floor_heating(Panel):
    bl_idname = "ARCHIPACK_PT_floor_heating"
    bl_label = "Floor heating"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Archipack"

    @classmethod
    def poll(cls, context):
        # ensure your object panel only show when active object is the right one
        return archipack_floor_heating.filter(context.active_object)

    def draw(self, context):
        o = context.active_object
        if not archipack_floor_heating.filter(o):
            return
        layout = self.layout
        scene = context.scene
        # retrieve datablock of your object
        props = archipack_floor_heating.datablock(o)
        # manipulate
        layout.operator("archipack.manipulate", icon="HAND")
        layout.separator()
        box = layout.box()
        row = box.row(align=True)

        # Presets operators
        row.operator("archipack.floor_heating_preset_menu",
                     text=bpy.types.ARCHIPACK_OT_floor_heating_preset_menu.bl_label)
        row.operator("archipack.floor_heating_preset",
                      text="",
                      icon='ZOOMIN')
        row.operator("archipack.floor_heating_preset",
                      text="",
                      icon='ZOOMOUT').remove_active = True

        box = layout.box()
        box.operator('archipack.floor_heating_cutter').parent = o.name

        box = layout.box()
        box.label(text="From curve")
        box.prop_search(props, "user_defined_path", scene, "objects", text="", icon='OUTLINER_OB_CURVE')
        if props.user_defined_path != "":
            box.prop(props, 'user_defined_resolution')
        box.prop(props, 'x_offset')
        box.prop(props, 'start_location')
        box = layout.box()
        row = box.row()
        if props.parts_expand:
            row.prop(props, 'parts_expand', icon="TRIA_DOWN", icon_only=True, text="Parts", emboss=False)
            box.prop(props, 'n_parts')
            # box.prop(prop, 'closed')
            for i, part in enumerate(props.parts):
                part.draw(context, layout, i)
        else:
            row.prop(props, 'parts_expand', icon="TRIA_RIGHT", icon_only=True, text="Parts", emboss=False)
        layout.separator()
        box = layout.box()
        box.prop(props, 'pattern', text="")
        box.prop(props, 'spacing')
        box.prop(props, 'enable_radius')
        if props.enable_radius:
            box.prop(props, 'radius')
        box.prop(props, 'max_iter')
        box.prop(props, 'backward')
        box.label(text=props.iter)
        box.label(text=props.pipe_len)
        box.label(text=props.area)
        # thickness


class ARCHIPACK_PT_floor_heating_cutter(Panel):
    bl_idname = "ARCHIPACK_PT_floor_heating_cutter"
    bl_label = "Floor heating obstacle"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'

    @classmethod
    def poll(cls, context):
        return archipack_floor_heating_cutter.filter(context.active_object)

    def draw(self, context):
        prop = archipack_floor_heating_cutter.datablock(context.active_object)
        if prop is None:
            return
        layout = self.layout
        scene = context.scene
        box = layout.box()
        box.operator('archipack.manipulate', icon='HAND')
        box.prop(prop, 'operation', text="")
        prop.draw(layout, context)


# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------


class ARCHIPACK_OT_floor_heating(ArchipackCreateTool, Operator):
    bl_idname = "archipack.floor_heating"
    bl_label = "Floor heating"
    bl_description = "Floor heating pipes"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    def create(self, context):
        """
            expose only basic params in operator
            use object property for other params
        """
        c = bpy.data.curves.new("Floor heating", type='CURVE')
        c.dimensions = '3D'
        o = bpy.data.objects.new("Floor heating", c)
        d = c.archipack_floor_heating.add()
        # make manipulators selectable
        d.manipulable_selectable = True
        angle_90 = pi / 2
        x, y, = 4, 4
        p = d.parts.add()
        p.a0 = 0
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
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        self.load_preset(d)
        self.add_material(o)
        return o

    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.location = context.scene.cursor_location
            # activate manipulators at creation time
            o.select = True
            context.scene.objects.active = o
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_floor_heating_from_curve(ArchipackCreateTool, Operator):
    bl_idname = "archipack.floor_heating_from_curve"
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
        bpy.ops.archipack.floor_heating(auto_manipulate=self.auto_manipulate, filepath=self.filepath)
        o = context.active_object
        d = archipack_floor_heating.datablock(o)
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


class ARCHIPACK_OT_floor_heating_from_wall(ArchipackCreateTool, Operator):
    bl_idname = "archipack.floor_heating_from_wall"
    bl_label = "->Floor"
    bl_description = "Create a floor from a wall"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(self, context):
        o = context.active_object
        return o is not None and o.data is not None and 'archipack_wall2' in o.data

    def create(self, context):
        wall = context.active_object
        wd = wall.data.archipack_wall2[0]
        bpy.ops.archipack.floor_heating(auto_manipulate=False, filepath=self.filepath)
        o = context.scene.objects.active
        d = archipack_floor_heating.datablock(o)
        d.auto_update = False
        d.closed = True
        d.parts.clear()
        d.n_parts = wd.n_parts + 1
        for part in wd.parts:
            p = d.parts.add()
            if "S_" in part.type:
                p.type = "S_SEG"
            else:
                p.type = "C_SEG"
            p.length = part.length
            p.radius = part.radius
            p.da = part.da
            p.a0 = part.a0

        side = 1
        if wd.flip:
            side = -1
        d.x_offset = -side * 0.5 * (1 - wd.x_offset) * wd.width

        d.auto_update = True
        # pretranslate
        o.matrix_world = wall.matrix_world.copy()
        return o

    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            o.select = True
            context.scene.objects.active = o
            if self.auto_manipulate:
                bpy.ops.archipack.floor_manipulate('INVOKE_DEFAULT')
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


class ARCHIPACK_OT_floor_heating_cutter(ArchipackCreateTool, Operator):
    bl_idname = "archipack.floor_heating_cutter"
    bl_label = "Floor Heating Cutter"
    bl_description = "Floor Heating Cutter"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}

    parent = StringProperty("")
    curve = StringProperty("")

    def create(self, context):
        m = bpy.data.meshes.new("Floor Cutter")
        o = bpy.data.objects.new("Floor Cutter", m)
        d = m.archipack_floor_heating_cutter.add()
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
            o.matrix_world = parent.matrix_world * Matrix([
                [1, 0, 0, -3 * x],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
                ])
            p = d.parts.add()
            p.a0 = - angle_90
            p.length = y
            p = d.parts.add()
            p.a0 = angle_90
            p.length = x
            p = d.parts.add()
            p.a0 = angle_90
            p.length = y
            d.n_parts = 3
            # d.close = True
            # pd = archipack_floor_heating.datablock(parent)
            # pd.boundary = o.name
        else:
            o.location = context.scene.cursor_location
        # make manipulators selectable
        d.manipulable_selectable = True
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o
        # self.add_material(o)
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
            o.select = True
            context.scene.objects.active = o
            self.manipulate()
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}


# ------------------------------------------------------------------
# Define operator class to manipulate object
# ------------------------------------------------------------------


class ARCHIPACK_OT_floor_heating_preset_menu(PresetMenuOperator, Operator):
    bl_description = "Show Floor heating presets"
    bl_idname = "archipack.floor_heating_preset_menu"
    bl_label = "Presets"
    preset_subdir = "archipack_floor_heating"


class ARCHIPACK_OT_floor_heating_preset(ArchipackPreset, Operator):
    """Add a Floor Heating Preset"""
    bl_idname = "archipack.floor_heating_preset"
    bl_label = "Add Floor heating preset"
    preset_menu = "ARCHIPACK_OT_floor_heating_preset_menu"

    @property
    def blacklist(self):
        return ['manipulators', 'parts', 'n_parts', 'user_defined_path', 'user_defined_resolution']


def register():
    bpy.utils.register_class(archipack_floor_heating_cutter_segment)
    bpy.utils.register_class(archipack_floor_heating_cutter)
    Mesh.archipack_floor_heating_cutter = CollectionProperty(type=archipack_floor_heating_cutter)
    bpy.utils.register_class(ARCHIPACK_OT_floor_heating_cutter)
    bpy.utils.register_class(ARCHIPACK_PT_floor_heating_cutter)
    bpy.utils.register_class(archipack_floor_heating_part)
    bpy.utils.register_class(archipack_floor_heating)
    Curve.archipack_floor_heating = CollectionProperty(type=archipack_floor_heating)
    bpy.utils.register_class(ARCHIPACK_PT_floor_heating)
    bpy.utils.register_class(ARCHIPACK_OT_floor_heating)
    bpy.utils.register_class(ARCHIPACK_OT_floor_heating_preset_menu)
    bpy.utils.register_class(ARCHIPACK_OT_floor_heating_preset)
    bpy.utils.register_class(ARCHIPACK_OT_floor_heating_from_curve)
    bpy.utils.register_class(ARCHIPACK_OT_floor_heating_from_wall)


def unregister():
    bpy.utils.unregister_class(archipack_floor_heating_cutter_segment)
    bpy.utils.unregister_class(archipack_floor_heating_cutter)
    del Mesh.archipack_floor_heating_cutter
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_heating_cutter)
    bpy.utils.unregister_class(ARCHIPACK_PT_floor_heating_cutter)
    bpy.utils.unregister_class(archipack_floor_heating_part)
    bpy.utils.unregister_class(archipack_floor_heating)
    del Curve.archipack_floor_heating
    bpy.utils.unregister_class(ARCHIPACK_PT_floor_heating)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_heating)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_heating_preset_menu)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_heating_preset)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_heating_from_curve)
    bpy.utils.unregister_class(ARCHIPACK_OT_floor_heating_from_wall)
