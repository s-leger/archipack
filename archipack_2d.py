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
from mathutils import Vector
from math import sin, cos, pi, atan2, sqrt, acos


class Projection():
    def proj_xy(self, t, next=None):
        """
            length of projection of sections at crossing line / circle intersections
            deformation unit vector for profil in xy axis
            so f(x_profile) = position of point in xy plane
        """
        if next is None:
            return self.normal(t).v.normalized(), 1
        v0 = self.normal(1).v.normalized()
        v1 = next.normal(0).v.normalized()
        direction = v0 + v1
        adj = (v0 * self.length) * (v1 * next.length)
        hyp = (self.length * next.length)
        c = min(1, max(-1, adj / hyp))
        size = 1 / cos(0.5 * acos(c))
        return direction.normalized(), size

    def proj_z(self, t, dz0, next=None, dz1=0):
        """
            length of projection along crossing line / circle
            deformation unit vector for profil in z axis at line / line intersection
            so f(y) = position of point in yz plane
        """
        return Vector((0, 1)), 1
        """
            NOTE (to myself):
              In theory this is how it has to be done so sections follow path,
              but in real world results are better when sections are z-up.
              So return a dumb 1 so f(y) = y
        """
        if next is None:
            dz = dz0 / self.length
        else:
            dz = (dz1 + dz0) / (self.length + next.length)
        return Vector((0, 1)), sqrt(1 + dz * dz)
        # 1 / sqrt(1 + (dz0 / self.length) * (dz0 / self.length))
        if next is None:
            return Vector((-dz0, self.length)).normalized(), 1
        v0 = Vector((self.length, dz0))
        v1 = Vector((next.length, dz1))
        direction = Vector((-dz0, self.length)).normalized() + Vector((-dz1, next.length)).normalized()
        adj = v0 * v1
        hyp = (v0.length * v1.length)
        c = min(1, max(-1, adj / hyp))
        size = -cos(pi - 0.5 * acos(c))
        return direction.normalized(), size


class Line(Projection):
    def __init__(self, p, v):
        self.p = p
        self.v = v
        self.v2 = v * v

    @property
    def length(self):
        return self.v.length

    @property
    def angle(self):
        return atan2(self.v.y, self.v.x)

    @property
    def angle_normal(self):
        return atan2(-self.v.x, self.v.y)

    @property
    def reversed(self):
        return Line(self.p, -self.v)

    @property
    def oposite(self):
        return Line(self.p + self.v, -self.v)

    @property
    def cross_z(self):
        return Vector((self.v.y, -self.v.x))

    def normal(self, t=0):
        # perpendicular on the right of segment
        return Line(self.lerp(t), self.cross_z)

    def sized_normal(self, t, size):
        return Line(self.lerp(t), size * self.cross_z.normalized())

    def lerp(self, t):
        return self.p + self.v * t

    def intersect(self, line):
        """ point_sur_segment return
            p: point d'intersection
            t: param t de l'intersection sur le segment courant
        """
        c = line.cross_z
        d = self.v * c
        if d == 0:
            return False, 0, 0
        t = (c * (line.p - self.p)) / d
        return True, self.lerp(t), t

    def steps(self, len):
        steps = max(1, round(self.length / len, 0))
        return 1 / steps, int(steps)

    def offset(self, offset):
        """
            offset > 0 on the right part
        """
        return Line(self.p + offset * self.cross_z.normalized(), self.v)

    def tangeant(self, t, da, radius):
        p = self.lerp(t)
        if da < 0:
            c = p + radius * self.cross_z.normalized()
        else:
            c = p - radius * self.cross_z.normalized()
        return Arc(c, radius, self.angle_normal, da)

    def straight(self, length, t=1):
        return Line(self.lerp(t), self.v.normalized() * length)

    def rotate(self, da):
        cs = cos(da)
        sn = sin(da)
        x, y = self.v
        self.v.x = x * cs - y * sn
        self.v.y = x * sn + y * cs
        return self

    def scale(self, length):
        self.v = length * self.v.normalized()
        return self

    def tangeant_unit_vector(self, t):
        return self.v.normalized()

    def point_sur_segment(self, pt):
        """ _point_sur_segment
            point: Vector 2d
            t: param t de l'intersection sur le segment courant
            d: distance laterale perpendiculaire positif a droite
        """
        dp = pt - self.p
        dl = self.length
        d = (self.v.x * dp.y - self.v.y * dp.x) / dl
        t = (self.v * dp) / (dl * dl)
        return t > 0 and t < 1, d, t


class Circle(Projection):
    def __init__(self, c, radius):
        self.r = radius
        self.r2 = radius * radius
        self.c = c

    def intersect(self, line):
        v = line.p - self.c
        A = line.v2
        B = 2 * v * line.v
        C = v * v - self.r2
        d = B * B - 4 * A * C
        if A <= 0.0000001 or d < 0:
            res, d, t = line.point_sur_segment(self.c)
            return False, line.lerp(t), t
        elif d == 0:
            t = -B / 2 * A
            return True, line.lerp(t), t
        else:
            AA = 2 * A
            dsq = sqrt(d)
            t0 = (-B + dsq) / AA
            t1 = (-B - dsq) / AA
            if abs(t0) < abs(t1):
                return True, line.lerp(t0), t0
            else:
                return True, line.lerp(t1), t1


class Arc(Circle):

    def __init__(self, c, radius, a0, da):
        """
            a0 and da arguments are in radians
            a0 = 0   on the right side
            a0 = pi on the left side
            da > 0 CCW contrary-clockwise
            da < 0 CW  clockwise
            stored internally as radians
        """
        Circle.__init__(self, c, radius)
        self.a0 = a0
        self.da = da

    @property
    def length(self):
        return self.r * abs(self.da)

    def normal(self, t=0):
        """
            always on the right side
        """
        p = self.lerp(t)
        if self.da < 0:
            return Line(p, self.c - p)
        else:
            return Line(p, p - self.c)

    def sized_normal(self, t, size):
        p = self.lerp(t)
        if self.da < 0:
            v = self.c - p
        else:
            v = p - self.c
        return Line(p, size * v.normalized())

    def lerp(self, t):
        a = self.a0 + t * self.da
        return self.c + Vector((self.r * cos(a), self.r * sin(a)))

    # not the same, wall does use angle instead of length..
    # maybe use a 2nd method here ?
    def steps(self, len):
        steps = max(1, round(self.length / len, 0))
        return 1.0 / steps, int(steps)

    # this is for wall
    def steps_by_angle(self, step_angle):
        steps = max(1, round(abs(self.da) / step_angle, 0))
        return 1.0 / steps, int(steps)

    def offset(self, offset):
        """
            offset > 0 on the right part
        """
        if self.da > 0:
            radius = self.r + offset
        else:
            radius = self.r - offset
        return Arc(self.c, radius, self.a0, self.da)

    def tangeant(self, t, length):
        a = self.a0 + t * self.da
        ca = cos(a)
        sa = sin(a)
        p = self.c + Vector((self.r * ca, self.r * sa))
        v = Vector((length * sa, -length * ca))
        if self.da > 0:
            v = -v
        return Line(p, v)

    def tangeant_unit_vector(self, t):
        a = self.a0 + t * self.da
        ca = cos(a)
        sa = sin(a)
        v = Vector((sa, -ca))
        if self.da > 0:
            v = -v
        return v

    def straight(self, length, t=1):
        s = self.tangeant(t, length)
        return s

    def point_sur_segment(self, pt):
        dp = pt - self.c
        d = dp.length - self.r
        a = atan2(dp.y, dp.x)
        t = (a - self.a0) / self.da
        return t > 0 and t < 1, d, t

    def rotate(self, da):
        cs = cos(da)
        sn = sin(da)
        x, y = self.v
        self.v.x = x * cs - y * sn
        self.v.y = x * sn + y * cs
        return self
