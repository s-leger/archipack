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

# ----------------------------------------------------------
# Author: Stephen Leger (s-leger)
#
# ----------------------------------------------------------
# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import FloatProperty, BoolProperty, IntProperty, CollectionProperty , EnumProperty
from .bmesh_utils import BmeshEdit as bmed
from .materialutils import MaterialUtils
from mathutils import Vector
from math import sin, cos, pi, atan2, sqrt

class Line():
    def __init__(self, p, v):
        self.p = p
        self.v = v
        self.v2 = v*v
    @property    
    def length(self):
        return self.v.length
    @property
    def angle(self):
        return atan2(self.v.y, self.v.x)
    @property
    def angle_normal(self):
        return atan2( - self.v.x, self.v.y )
    @property
    def reversed(self):
        return Line(self.p, -self.v)
    @property
    def cross_z(self):
        return Vector((self.v.y, - self.v.x))
    def normal(self, t=0):
        # perpendiculaire a droite du segment
        return Line(self.lerp(t), self.cross_z)
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
        t = (c * (line.p-self.p))/d
        return True, self.lerp(t), t
    def steps(self, len):
        steps = round(self.length / len, 0)
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
    def straight(self, length):
        return self.p+self.v, self.v.normalized()*length
      
class Circle():
    def __init__(self, c, radius):
        self.r = radius
        self.r2 = radius*radius
        self.c = c
    def intersect(self, line):
        v = line.p - self.c
        A = line.v2
        B = 2 * v * line.v
        C = v * v - self.r2
        d = B*B-4*A*C
        if A <= 0.0000001 or d < 0:
            return False, 0, 0
        elif d == 0:
            t = -B / 2*A
            return True, line.lerp(t), t
        else:
            AA = 2*A
            dsq = sqrt(d)
            t0 = (-B + dsq)/AA
            t1 = (-B - dsq)/AA
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
        super().__init__(c, radius)
        self.a0 = a0
        self.da = da
    @property
    def length(self):
        return self.r*abs(self.da)
    def normal(self, t=0):
        """
            always on the right side
        """
        p = self.lerp(t)
        if self.da < 0:
            return Line(p, self.c-p)
        else:
            return Line(p, p-self.c)
    def lerp(self, t):
        a = self.a0+t*self.da
        return self.c+Vector((self.r*cos(a),self.r*sin(a)))
    def steps(self, len):
        steps = round(self.length / len, 0)
        return 1.0 / steps, int(steps)
    def offset(self, offset):
        """
            offset > 0 on the right part
        """
        if self.da > 0:
            radius = self.r+offset
        else:
            radius = self.r-offset
        return Arc(self.c, radius, self.a0, self.da)
    def tangeant(self, t, length):
        a = self.a0+t*self.da
        ca = cos(a)
        sa = sin(a)
        p = self.c+Vector((self.r*ca,self.r*sa))
        v = Vector((length*sa, -length*ca)) 
        if self.da > 0:
            v = -v
        return Line(p, v)
    def straight(self, length):
        t = self.tangeant(1,length)
        return t.p, t.v
  
class Stair():
    def __init__(self, left_offset, right_offset, faces_type, z_mode, step_z, bottom_z):
        self.faces_type = faces_type
        self.next_type = 'NONE'
        self.last_type = 'NONE'
        self.z_mode = z_mode
        # depth of open step
        self.step_z = step_z
        # size under the step on bottom
        self.bottom_z = bottom_z
        self.left_offset = left_offset
        self.right_offset = right_offset
        self.last_height = 0
    def set_height(self, step_height, z0):
        self.step_height = step_height
        self.z0 = z0
    @property 
    def height(self):
        return self.n_step * self.step_height
    @property
    def top_offset(self):
        return self.t_step / self.step_depth
    @property
    def top(self):
        return self.z0 + self.height
    def step_size(self, step_depth):
        t_step, n_step = self.steps(step_depth) 
        self.n_step = n_step
        self.t_step = t_step
        self.step_depth = step_depth
        return n_step
    def p3d_bottom(self, verts, p2d, i, t):
        x, y = p2d
        if 'LINEAR' in self.z_mode:
            z = self.z0 + t * self.height
        else:    
            z = self.z0 + i * self.step_height
        if self.faces_type == 'FULL':
            zb = 0
        else:
            zb = max(0, z - self.bottom_z)
        verts.append((x, y, zb))
        verts.append((x, y, z))
    def p3d_top(self, verts, p2d, i, t):
        x, y = p2d
        if self.z_mode == 'LINEAR_TOP':
            z = self.z0 + t * self.height
        else:    
            z = self.z0 + (i+1) * self.step_height
        zb = z - self.step_z
        verts.append((x, y, zb))
        verts.append((x, y, z))
    def p3d_end(self, verts, p2d, i, t):
        x, y = p2d
        z3 = self.z0 + self.step_height
        z2 = z3 - self.bottom_z
        z1 = self.z0 - self.bottom_z + 0.5 * self.step_height
        if self.faces_type == 'CLOSED':
            z0 = max(0, self.z0 - self.bottom_z)
        else:
            z0 = 0
        verts.append((x, y, z0))
        verts.append((x, y, z1))
        verts.append((x, y, z2))
        verts.append((x, y, z3))
    def straight_stair(self, length):
        self.next_type = 'STAIR'
        p, v = self.straight(length)
        return StraightStair(p, v, self.left_offset, self.right_offset, self.faces_type, self.z_mode, self.step_z, self.bottom_z)
    def straight_palier(self, length, last_type='STAIR'):
        self.next_type = 'PALIER'
        p, v = self.straight(length)
        return StraightPalier(p, v, self.left_offset, self.right_offset, self.faces_type, self.z_mode, self.step_z, self.bottom_z, last_type=last_type)
    def curved_stair(self, da, radius, left_shape, right_shape, double_limit=pi):
        self.next_type = 'STAIR'
        n = self.normal(1)
        n.v = radius * n.v.normalized()
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedStair(c, radius, a0, da, self.left_offset, self.right_offset, self.faces_type, self.z_mode, self.step_z, self.bottom_z, left_shape, right_shape, double_limit=double_limit)
    def curved_palier(self, da, radius, left_shape, right_shape, double_limit=pi, last_type='STAIR'):
        self.next_type = 'PALIER'
        n = self.normal(1)
        n.v = radius * n.v.normalized()
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedPalier(c, radius, a0, da, self.left_offset, self.right_offset, self.faces_type, self.z_mode, self.step_z, self.bottom_z, left_shape, right_shape, double_limit=double_limit, last_type=last_type)
    
class StraightStair(Stair, Line):
    def __init__(self, p, v, left_offset, right_offset, faces_type, z_mode, step_z, bottom_z):
        Stair.__init__(self, left_offset, right_offset, faces_type, z_mode, step_z, bottom_z)
        Line.__init__(self, p, v)
        self.l_line = self.offset(left_offset)
        self.r_line = self.offset(-right_offset)
    def make_step(self, i, verts, faces, t_offset=0):
        tb = self.t_step * i
        tt = tb+t_offset
        nb = self.normal(tb)
        nt = self.normal(tt)
        f = len(verts)
        res, p, u = self.l_line.intersect(nb)
        self.p3d_bottom(verts, p, i, u)
        res, p, u = self.l_line.intersect(nt)
        self.p3d_top(verts, p, i, u)
        res, p, u = self.r_line.intersect(nb)
        self.p3d_bottom(verts, p, i, u)
        res, p, u = self.r_line.intersect(nt)
        self.p3d_top(verts, p, i, u)
        if i < self.n_step:
            # contre marche haut
            faces.append((f+2, f+3, f+7, f+6))
            # dessus marche rectangle
            faces.append((f+3, f+9, f+13, f+7))
            if self.faces_type in ['CLOSED', 'FULL']:
                # contre-marche bas
                faces.append((f+1, f+2, f+6, f+5))
                # dessous marche rectangle
                faces.append((f, f+4, f+12, f+8 ))  
                # side right top
                faces.append((f+1, f+9, f+3, f+2))
                # side left top
                faces.append((f+5, f+6, f+7, f+13))
                # side left bottom
                faces.append((f+4, f+5, f+13, f+12))
                # side right bottom
                faces.append((f, f+8, f+9, f+1))
            else:
                # side left top
                faces.append((f+6, f+7, f+13, f+12))
                # side right top
                faces.append((f+8, f+9, f+3, f+2))
                # dessous marche rectangle
                faces.append((f+6, f+12, f+8, f+2))
                # arriere
                faces.append((f+8, f+12, f+13, f+9))
        elif self.faces_type in ['CLOSED', 'FULL']:
            faces.append((f, f+4, f+5, f+1))       

class CurvedStair(Stair, Arc):
    def __init__(self, c, radius, a0, da, left_offset, right_offset, faces_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=pi):
        Stair.__init__(self, left_offset, right_offset, faces_type, z_mode, step_z, bottom_z)
        Arc.__init__(self, c, radius, a0, da)
        self.l_shape = left_shape 
        self.r_shape = right_shape
        self.edges_multiples = round(abs(da),6) > double_limit
        #left arc, tangeant at start and end
        self.l_arc = self.offset(left_offset)
        if left_shape != 'CIRCLE':
            self.l_t0 = self.l_arc.tangeant(0,1)
            self.l_t1 = self.l_arc.tangeant(1,1)
            if self.edges_multiples: 
                self.l_tc = self.l_arc.tangeant(0.5,1)
                i, p, t = self.l_t0.intersect(self.l_tc)
                self.l_tc.v *= 2*t
                self.l_tc.p = p
                i, p, t2 = self.l_tc.intersect(self.l_t1)
            else:
                i, p, t = self.l_t0.intersect(self.l_t1)
            self.l_t0.v *= t
            self.l_t1.p = p
            self.l_t1.v *= t
        #right arc, tangeant at start and end
        self.r_arc = self.offset(-right_offset)
        if right_shape != 'CIRCLE':
            self.r_t0 = self.r_arc.tangeant(0,1)
            self.r_t1 = self.r_arc.tangeant(1,1)
            if self.edges_multiples: 
                self.r_tc = self.r_arc.tangeant(0.5,1)
                i, p, t = self.r_t0.intersect(self.r_tc)
                self.r_tc.v *= 2*t
                self.r_tc.p = p
                i, p, t2 = self.r_tc.intersect(self.r_t1)
            else:
                i, p, t = self.r_t0.intersect(self.r_t1)
            self.r_t0.v *= t
            self.r_t1.p = p
            self.r_t1.v *= t
    def p3d_corner(self, verts, p2d, i, t):
        x, y = p2d
        if self.faces_type == 'FULL':
            z = 0
        else:
            if 'LINEAR'in self.z_mode:
                z = self.z0 + t * self.height - self.bottom_z
            else:    
                z = self.z0 + (i+0.5) * self.step_height - self.bottom_z
        verts.append((x, y, z))
    def make_faces(self, i, f, faces, e_right, e_left):
        if i < self.n_step:
            # contre marche haut
            faces.append((f+2, f+3, f+7, f+6))
            if self.faces_type in ['CLOSED', 'FULL']:
                # contre-marche bas
                faces.append((f+1, f+2, f+6, f+5))
                if e_right and e_left:
                    # dessus faces with two more edges
                    faces.append((f+3, f+10, f+15, f+19, f+13, f+7))
                    # side left bottom haut
                    faces.append((f+5, f+6, f+7, f+13, f+12))
                    # side left bottom 
                    faces.append((f+4, f+5, f+12, f+11))
                    # side left top haut
                    faces.append((f+12, f+13, f+19))
                    # side left top bas
                    faces.append((f+11, f+12, f+19, f+18))
                    # side right top bas
                    faces.append((f+8, f+14, f+15, f+9))
                    # side right top haut
                    faces.append((f+9, f+15, f+10))
                    # side right bottom bas
                    faces.append((f, f+8, f+9, f+1))
                    # side right bottom haut
                    faces.append((f+1, f+9, f+10, f+3, f+2))
                    # bottom
                    faces.append((f, f+4, f+11, f+8))
                    faces.append((f+11, f+18, f+14, f+8))
                elif e_right:
                    #dessus right with one more edge
                    faces.append((f+3, f+12, f+16, f+10, f+7))
                    # side right haut
                    faces.append((f+1, f+12, f+3, f+2))
                    # side right bas
                    faces.append((f, f+11, f+12, f+1))
                    # side left bottom bas
                    faces.append((f+4, f+5, f+9, f+8))
                    # side left bottom haut
                    faces.append((f+5, f+6, f+7, f+10, f+9))
                    # side left top haut
                    faces.append((f+9, f+10, f+16))
                    # side left top bas
                    faces.append((f+15, f+8, f+9, f+16))
                    # bottom
                    faces.append((f, f+4, f+8, f+15, f+11))
                elif e_left:
                    #dessus  left with one more edge 
                    faces.append((f+3, f+10, f+12, f+16, f+7))
                    # side left haut
                    faces.append((f+5, f+6, f+7, f+16))
                    # side left bas
                    faces.append((f+4, f+5, f+16, f+15))
                    # side right top haut
                    faces.append((f+9, f+12, f+10))
                    # side right top bas
                    faces.append((f+8, f+11, f+12, f+9))
                    # side right bottom bas
                    faces.append((f, f+8, f+9, f+1))
                    # side right bottom haut
                    faces.append((f+1, f+9, f+10, f+3, f+2))
                    # bottom
                    faces.append((f, f+4, f+15, f+11, f+8))
                else:
                    # dessus marche rectangle
                    faces.append((f+3, f+9, f+13, f+7))
                    # dessous marche rectangle
                    faces.append((f, f+4, f+12, f+8 )) 
                    # side right top
                    faces.append((f+1, f+9, f+3, f+2))
                    # side left top
                    faces.append((f+5, f+6, f+7, f+13))
                    # side right bottom
                    faces.append((f, f+8, f+9, f+1))
                    # side left bottom
                    faces.append((f+4, f+5, f+13, f+12))
            else:
                if e_right and e_left:
                    # dessus faces with two more edges
                    faces.append((f+3, f+10, f+15, f+19, f+13, f+7))
                    # dessous faces with two more edges
                    faces.append((f+6, f+12, f+18, f+14, f+9, f+2))
                    # side left bottom haut
                    faces.append((f+6, f+7, f+13, f+12))
                    # side left top haut
                    faces.append((f+12, f+13, f+19, f+18))
                    # side right bottom haut
                    faces.append((f+2, f+9, f+10, f+3))
                    # side right top haut
                    faces.append((f+9, f+14, f+15, f+10))                    
                elif e_right:
                    #dessus right with one more edge
                    faces.append((f+3, f+12, f+16, f+10, f+7))
                    #dessous right with one more edge
                    faces.append((f+6, f+9, f+15, f+11, f+2))
                    # side right haut
                    faces.append((f+11, f+12, f+3, f+2))
                    # side left bottom haut
                    faces.append((f+6, f+7, f+10, f+9))
                    # side left top haut
                    faces.append((f+9, f+10, f+16, f+15))
                elif e_left:
                    #dessus  left with one more edge 
                    faces.append((f+3, f+10, f+12, f+16, f+7))
                    #dessous  left with one more edge 
                    faces.append((f+6, f+15, f+11, f+9, f+2))
                    # side right top haut
                    faces.append((f+11, f+12, f+10, f+9))
                    # side right bottom haut
                    faces.append((f+9, f+10, f+3, f+2))
                    # side left haut
                    faces.append((f+6, f+7, f+16, f+15))                
                else:
                    # dessus marche rectangle
                    faces.append((f+3, f+9, f+13, f+7))
                    # dessous marche rectangle
                    faces.append((f+6, f+12, f+8, f+2))
                    # side right top
                    faces.append((f+8, f+9, f+3, f+2))
                    # side left top
                    faces.append((f+6, f+7, f+13, f+12))
                    # arriere
                    faces.append((f+8, f+12, f+13, f+9))
        elif self.faces_type in ['CLOSED', 'FULL']:
            # close last face
            faces.append((f, f+4, f+5, f+1))
    def make_step(self, i, verts, faces, t_offset=0):
        tb = self.t_step * i
        tt = tb+t_offset
        nbl = self.r_arc.normal(tb)
        ntl = self.r_arc.normal(tt)
        nbr = self.l_arc.normal(tb).reversed
        ntr = self.l_arc.normal(tt).reversed
        e_left  = False
        e_right = False
        f = len(verts)
        if self.l_shape == 'CIRCLE':
            res, p, u = self.l_arc.intersect(nbl)
            self.p3d_bottom(verts, p, i, tb)
            res, p, u = self.l_arc.intersect(ntl)
            self.p3d_top(verts, p, i, tt)
        else:
            if self.edges_multiples:    
                # two edges
                if tb < 0.25:
                    res, p, u = self.l_t0.intersect(nbl)
                    self.p3d_bottom(verts, p, i, 0.25 * u)
                elif tb < 0.75:
                    res, p, u = self.l_tc.intersect(nbl)
                    self.p3d_bottom(verts, p, i, 0.25 + 0.5 * u)
                else:
                    res, p, u = self.l_t1.intersect(nbl)
                    self.p3d_bottom(verts, p, i, 0.75 + 0.25 * u)
                if tb < 0.25:
                    res, p, u = self.l_t0.intersect(ntl)
                    self.p3d_top(verts, p, i, 0.25 * u)
                elif tb < 0.75:
                    res, p, u = self.l_tc.intersect(ntl)
                    self.p3d_top(verts, p, i, 0.25 + 0.5 * u)
                else:
                    res, p, u = self.l_t1.intersect(ntl)
                    self.p3d_top(verts, p, i, 0.75 + 0.25 * u)
            else:    
                if tb < 0.5:
                    res, p, u = self.l_t0.intersect(nbl)
                    self.p3d_bottom(verts, p, i, 0.5 * u)
                else:
                    res, p, u = self.l_t1.intersect(nbl)
                    self.p3d_bottom(verts, p, i, 0.5 + 0.5 * u)
                if tb < 0.5:
                    res, p, u = self.l_t0.intersect(ntl)
                    self.p3d_top(verts, p, i, 0.5 * u)
                else:    
                    res, p, u = self.l_t1.intersect(ntl)
                    self.p3d_top(verts, p, i, 0.5 + 0.5 * u)
        if self.r_shape == 'CIRCLE':
            res, p, u = self.r_arc.intersect(nbr)
            self.p3d_bottom(verts, p, i, tb)
            res, p, u = self.r_arc.intersect(ntr)
            self.p3d_top(verts, p, i, tt)
        else:
            if self.edges_multiples:
                # two edges
                if tb < 0.25:
                    res, p, u = self.r_t0.intersect(nbr)
                    self.p3d_bottom(verts, p, i, 0.25 * u)
                elif tb < 0.75:
                    res, p, u = self.r_tc.intersect(nbr)
                    self.p3d_bottom(verts, p, i, 0.25 + 0.5 * u)
                else:
                    res, p, u = self.r_t1.intersect(nbr)
                    self.p3d_bottom(verts, p, i, 0.75 + 0.25 * u)
                if tb < 0.25:
                    res, p, u = self.r_t0.intersect(ntr)
                    self.p3d_top(verts, p, i, 0.25 * u)
                elif tb < 0.75:
                    res, p, u = self.r_tc.intersect(ntr)
                    self.p3d_top(verts, p, i, 0.25 + 0.5 * u)
                else:
                    res, p, u = self.r_t1.intersect(ntr)
                    self.p3d_top(verts, p, i, 0.75 + 0.25 * u)
            else:
                # one edge
                if tb < 0.5:
                    res, p, u = self.r_t0.intersect(nbr)
                    self.p3d_bottom(verts, p, i, 0.5 * u)
                else:
                    res, p, u = self.r_t1.intersect(nbr)
                    self.p3d_bottom(verts, p, i, 0.5 + 0.5 * u)
                if tb < 0.5:
                    res, p, u = self.r_t0.intersect(ntr)
                    self.p3d_top(verts, p, i, 0.5 * u)
                else:
                    res, p, u = self.r_t1.intersect(ntr)
                    self.p3d_top(verts, p, i, 0.5 + 0.5 * u)
        # make edges verts after regular ones
        if self.edges_multiples:
            # edge 1
            if tb < 0.25 and tb + self.t_step > 0.25:     
                if self.l_shape != 'CIRCLE': 
                    e_left = True
                    # the step goes through the edge
                    self.p3d_corner(verts, self.l_tc.p, i, 0.25)
                    self.p3d_top(verts, self.l_tc.p, i, 0.25)                
                if self.r_shape != 'CIRCLE':
                    e_right = True
                    # the step goes through the edge
                    self.p3d_corner(verts,  self.r_tc.p, i, 0.25)
                    self.p3d_top(verts,  self.r_tc.p, i, 0.25)
            # edge 2
            if tb < 0.75 and tb + self.t_step > 0.75:     
                if self.l_shape != 'CIRCLE': 
                    e_left = True
                    # the step goes through the edge
                    self.p3d_corner(verts, self.l_t1.p, i, 0.75)
                    self.p3d_top(verts, self.l_t1.p, i, 0.75)                
                if self.r_shape != 'CIRCLE':
                    e_right = True
                    # the step goes through the edge
                    self.p3d_corner(verts,  self.r_t1.p, i, 0.75)
                    self.p3d_top(verts,  self.r_t1.p, i, 0.75)
        else:        
            if tb < 0.5 and tb + self.t_step > 0.5:     
                if self.l_shape != 'CIRCLE': 
                    e_left = True
                    # the step goes through the edge
                    self.p3d_corner(verts, self.l_t1.p, i, 0.5)
                    self.p3d_top(verts, self.l_t1.p, i, 0.5)                
                if self.r_shape != 'CIRCLE':
                    e_right = True
                    # the step goes through the edge
                    self.p3d_corner(verts,  self.r_t1.p, i, 0.5)
                    self.p3d_top(verts,  self.r_t1.p, i, 0.5)
        self.make_faces(i, f, faces, e_right, e_left)
               
class StraightPalier(StraightStair):
    def __init__(self, p, v, left_offset, right_offset, faces_type, z_mode, step_z, bottom_z, last_type='STAIR'):
        StraightStair.__init__(self, p, v, left_offset, right_offset, faces_type, z_mode, step_z, bottom_z)
        self.last_type=last_type
    @property
    def height(self):
        return 0
    @property
    def top_offset(self):
        return self.t_step / self.v.length
    @property
    def top(self):
        if self.next_type == 'PALIER':
            return self.z0
        else:
            return self.z0 + self.step_height
    def step_size(self, step_depth):
        self.n_step = 1
        self.t_step = 1
        self.step_depth = step_depth
        if self.last_type == 'PALIER':
            return 0
        else:
            return 1
    def make_step(self, i, verts, faces, t_offset=0):
        f = len(verts)
        j = 0
        if i == self.n_step and self.next_type != 'PALIER' and self.faces_type in ['CLOSED','FULL']:
            p = self.l_line.lerp(1)
            self.p3d_end(verts, p, j, 1)
            p = self.r_line.lerp(1)
            self.p3d_end(verts, p, j, 1)
        else:
            tb = self.t_step * i
            tt = tb+t_offset
            if i > 0 or self.last_type == 'PALIER':
                # disabled between flat parts
                tt = tb
            nb = self.normal(tb)
            nt = self.normal(tt)
            res, p, u = self.l_line.intersect(nb)
            self.p3d_bottom(verts, p, j, u)
            res, p, u = self.l_line.intersect(nt)
            self.p3d_top(verts, p, j, u)
            res, p, u = self.r_line.intersect(nb)
            self.p3d_bottom(verts, p, j, u)
            res, p, u = self.r_line.intersect(nt)
            self.p3d_top(verts, p, j, u)
                 
        if i < self.n_step:  
            # contre marche haut
            faces.append((f+2, f+3, f+7, f+6))
            # dessus marche rectangle
            faces.append((f+3, f+11, f+15, f+7))
            if self.faces_type in ['CLOSED','FULL']:
                # contre-marche bas
                faces.append((f+1, f+2, f+6, f+5))
                # side right top
                faces.append((f+1, f+9, f+10, f+11, f+3, f+2))
                # side right bottom
                faces.append((f, f+8, f+9, f+1))
                # side left top
                faces.append((f+5, f+6, f+7, f+15, f+14, f+13))
                # side left bottom
                faces.append((f+4, f+5, f+13, f+12))
                # dessous marche rectangle
                faces.append((f, f+4, f+12, f+8 ))
            else:
                # dessous marche rectangle
                faces.append((f+6, f+14, f+10, f+2))
                # side right top
                faces.append((f+10, f+11, f+3, f+2))
                # side left top
                faces.append((f+6, f+7, f+15, f+14))
        else:
            if self.faces_type in ['CLOSED','FULL']:
                if self.next_type == 'STAIR':
                    # close last part bottom
                    faces.append((f, f+4, f+5, f+6, f+2,f+1))
                elif self.next_type == 'NONE':
                    faces.append((f, f+4, f+5, f+6, f+7, f+3, f+2, f+1))
    def straight_palier(self, length):
        return Stair.straight_palier(self, length, last_type='PALIER')
    def curved_palier(self, da, radius, left_shape, right_shape, double_limit=pi):
        return Stair.curved_palier(self, da, radius, left_shape, right_shape, double_limit=double_limit, last_type='PALIER')
        
class CurvedPalier(CurvedStair):
    # @TODO: Implements RAMP mode
    def __init__(self, c, radius, a0, da, left_offset, right_offset, faces_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=pi, last_type='STAIR', ramp_mode='FLAT'):
        CurvedStair.__init__(self, c, radius, a0, da, left_offset, right_offset, faces_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=double_limit)
        self.last_type=last_type
        self.ramp_mode = ramp_mode
    @property
    def top_offset(self):
        if self.l_shape == 'CIRCLE' or self.r_shape == 'CIRCLE':
            return self.t_step / self.step_depth
        else:
            if self.edges_multiples:
                return 0.5 / self.length
            else:
                return 1 / self.length
    @property
    def height(self):
        if self.ramp_mode == 'RAMP':
            return self.n_step * self.step_height
        return 0
    @property
    def top(self):
        if self.next_type == 'PALIER':
            return self.z0 + self.height
        else:
            return self.z0 + self.step_height
    def step_size(self, step_depth):
        if self.l_shape == 'CIRCLE' or self.r_shape == 'CIRCLE':
            t_step, n_step = self.steps(step_depth) 
        else:
            if self.edges_multiples:
                t_step, n_step = 0.5, 2
            else:
                t_step, n_step = 1, 1
        self.n_step = n_step
        self.t_step = t_step
        self.step_depth = step_depth
        if self.last_type == 'PALIER':
            return 0
        else:
            return 1
    def p3d_corner(self, verts, p2d, i, t):
        x, y = p2d
        if self.faces_type == 'FULL':
            z = 0
        else:
            z = self.z0 - self.bottom_z
        verts.append((x, y, z))
    def make_faces(self, i, f, faces, e_right, e_left):
        if i < self.n_step:
            if i == 0 and self.last_type != 'PALIER':
                # only on first part
                # disabled between flat parts
                # contre marche haut
                faces.append((f+2, f+3, f+7, f+6))
                if self.faces_type in ['CLOSED', 'FULL']:
                    # contre-marche bas
                    faces.append((f+1, f+2, f+6, f+5))
            if self.faces_type in ['CLOSED', 'FULL']:
                if e_right and e_left:
                    # dessus faces with two more edges
                    faces.append((f+3, f+10, f+17, f+21, f+13, f+7))
                    # side left bottom bas
                    faces.append((f+4, f+5, f+6, f+12, f+11))
                    # side left bottom haut
                    faces.append((f+6, f+7, f+13, f+12))
                    # side left top bas
                    faces.append((f+11, f+12, f+20, f+19, f+18))
                    # side left top haut
                    faces.append((f+12, f+13, f+21, f+20))
                    # side right bottom bas
                    faces.append((f, f+8, f+9, f+2, f+1))
                    # side right bottom haut
                    faces.append((f+2, f+9, f+10, f+3))
                    # side right top bas
                    faces.append((f+8, f+14, f+15, f+16, f+9))
                    # side right top haut
                    faces.append((f+10, f+9, f+16, f+17 ))
                    # bottom
                    faces.append((f, f+4, f+11, f+8))
                    faces.append((f+11, f+18, f+14, f+8))
                elif e_right:
                    #dessus right with one more edge
                    faces.append((f+3, f+14, f+18, f+10, f+7))
                    # side right bas
                    faces.append((f, f+11, f+12, f+1))
                    # side right haut
                    faces.append((f+12, f+13, f+14, f+3, f+2, f+1))
                    # side left bottom haut
                    faces.append((f+6, f+7, f+10, f+9))
                    # side left bottom bas
                    faces.append((f+4, f+5, f+6, f+9, f+8))
                    # side left top bas
                    faces.append((f+8, f+9, f+17, f+16, f+15))
                    # side left top haut
                    faces.append((f+9, f+10, f+18, f+17))
                    # bottom
                    faces.append((f, f+4, f+8, f+15, f+11))
                elif e_left:
                    #dessus  left with one more edge 
                    faces.append((f+3, f+10, f+14, f+18, f+7))
                    # side right bottom bas
                    faces.append((f, f+8, f+9, f+2, f+1))
                    # side right bottom haut
                    faces.append((f+9, f+10, f+3, f+2))
                    # side right top bas
                    faces.append((f+8, f+11, f+12, f+13, f+9))
                    # side right top haut
                    faces.append((f+9, f+13, f+14, f+10))
                    # side left haut
                    faces.append((f+5, f+6, f+7, f+18, f+17, f+16))
                    # side left bas
                    faces.append((f+4, f+5, f+16, f+15))
                    # bottom
                    faces.append((f, f+4, f+15, f+11, f+8))
                else:
                    # dessus marche rectangle
                    faces.append((f+3, f+11, f+15, f+7))
                    # side right top
                    faces.append((f+1, f+9, f+10, f+11, f+3, f+2))
                    # side right bottom
                    faces.append((f, f+8, f+9, f+1))
                    # side left top
                    faces.append((f+5, f+6, f+7, f+15, f+14, f+13))
                    # side left bottom
                    faces.append((f+4, f+5, f+13, f+12))
                    # dessous marche rectangle
                    faces.append((f, f+4, f+12, f+8 ))
            else:
                if e_right and e_left:
                    # dessus faces with two more edges
                    faces.append((f+3, f+10, f+17, f+21, f+13, f+7))
                    # dessous faces with two more edges
                    faces.append((f+6, f+12, f+20, f+16, f+9, f+2))
                    # side left bottom haut
                    faces.append((f+6, f+7, f+13, f+12))
                    # side left top haut
                    faces.append((f+12, f+13, f+21, f+20))
                    # side right bottom haut
                    faces.append((f+2, f+9, f+10, f+3))
                    # side right top haut
                    faces.append((f+10, f+9, f+16, f+17 ))
                elif e_right:
                    #dessus right with one more edge
                    faces.append((f+3, f+12, f+16, f+10, f+7))
                    #dessous right with one more edge
                    faces.append((f+6, f+9, f+15, f+11, f+2))
                    # side right haut
                    faces.append((f+11, f+12, f+3, f+2))
                    # side left bottom haut
                    faces.append((f+6, f+7, f+12, f+11))
                    # side left top haut
                    faces.append((f+9, f+10, f+16, f+15))
                elif e_left:
                    #dessus  left with one more edge 
                    faces.append((f+3, f+10, f+14, f+18, f+7))
                    #dessous  left with one more edge 
                    faces.append((f+6, f+17, f+13, f+9, f+2))
                    # side right bottom haut
                    faces.append((f+9, f+10, f+3, f+2))
                    # side right top haut
                    faces.append((f+9, f+13, f+14, f+10))
                    # side left haut
                    faces.append((f+6, f+7, f+18, f+17))
                else:
                    # dessus marche rectangle
                    faces.append((f+3, f+11, f+15, f+7))
                    # dessous marche rectangle
                    faces.append((f+6, f+14, f+10, f+2))
                    # side right top
                    faces.append((f+10, f+11, f+3, f+2))
                    # side left top
                    faces.append((f+6, f+7, f+15, f+14))
    def make_step(self, i, verts, faces, t_offset=0): 
        f = len(verts)
        if self.ramp_mode == 'RAMP':
            j = i
        else:
            j = 0
        if i == self.n_step and self.faces_type in ['CLOSED','FULL']:
            if self.next_type == 'STAIR':
                p = self.l_arc.lerp(1)
                self.p3d_end(verts, p, j, 1)
                p = self.r_arc.lerp(1)
                self.p3d_end(verts, p, j, 1)
                # close last part bottom
                faces.append((f, f+4, f+5, f+6, f+2,f+1))
                return
            elif self.next_type == 'NONE':
                faces.append((f, f+4, f+5, f+6, f+7, f+3, f+2, f+1))
            
        tb = self.t_step * i
        tt = tb+t_offset
        nb = self.normal(tb)
        if i > 0 or self.last_type == 'PALIER':
            # disabled between flat parts
            tt = tb
        nbl = self.r_arc.normal(tb)
        ntl = self.r_arc.normal(tt)
        nbr = self.l_arc.normal(tb).reversed
        ntr = self.l_arc.normal(tt).reversed
        e_left  = False
        e_right = False
        f = len(verts)

        if self.l_shape == 'CIRCLE':
            res, p, u = self.l_arc.intersect(nbl)
            self.p3d_bottom(verts, p, j, tb)
            res, p, u = self.l_arc.intersect(ntl)
            self.p3d_top(verts, p, j, tt)
        else:
            if self.edges_multiples:
                # two edges
                if tb < 0.25:
                    res, p, u = self.l_t0.intersect(nbl)
                    self.p3d_bottom(verts, p, j, 0.25 * u)
                elif tb < 0.75:
                    res, p, u = self.l_tc.intersect(nbl)
                    self.p3d_bottom(verts, p, j, 0.25 + 0.5 * u)
                else:
                    res, p, u = self.l_t1.intersect(nbl)
                    self.p3d_bottom(verts, p, j, 0.75 + 0.25 * u)
                if tt < 0.25:
                    res, p, u = self.l_t0.intersect(ntl)
                    self.p3d_top(verts, p, j, 0.25 * u)
                elif tt < 0.75:
                    res, p, u = self.l_tc.intersect(ntl)
                    self.p3d_top(verts, p, j, 0.25 + 0.5 * u)
                else:
                    res, p, u = self.l_t1.intersect(ntl)
                    self.p3d_top(verts, p, j, 0.75 + 0.25 * u)
            else:  
                if tb < 0.5:
                    res, p, u = self.l_t0.intersect(nbl)
                    self.p3d_bottom(verts, p, j, 0.5 * u)
                else:
                    res, p, u = self.l_t1.intersect(nbl)
                    self.p3d_bottom(verts, p, j, 0.5 + 0.5 * u)
                if tt < 0.5:
                    res, p, u = self.l_t0.intersect(ntl)
                    self.p3d_top(verts, p, j, 0.5 * u)
                else:    
                    res, p, u = self.l_t1.intersect(ntl)
                    self.p3d_top(verts, p, j, 0.5 + 0.5 * u)  
        if self.r_shape == 'CIRCLE':
            res, p, u = self.r_arc.intersect(nbr)
            self.p3d_bottom(verts, p, j, tb)
            res, p, u = self.r_arc.intersect(ntr)
            self.p3d_top(verts, p, j, tt)
        else:
            if self.edges_multiples:
                # two edges
                if tb < 0.25:
                    res, p, u = self.r_t0.intersect(nbr)
                    self.p3d_bottom(verts, p, j, 0.25 * u)
                elif tb < 0.75:
                    res, p, u = self.r_tc.intersect(nbr)
                    self.p3d_bottom(verts, p, j, 0.25 + 0.5 * u)
                else:
                    res, p, u = self.r_t1.intersect(nbr)
                    self.p3d_bottom(verts, p, j, 0.75 + 0.25 * u)
                if tt < 0.25:
                    res, p, u = self.r_t0.intersect(ntr)
                    self.p3d_top(verts, p, j, 0.25 * u)
                elif tt < 0.75:
                    res, p, u = self.r_tc.intersect(ntr)
                    self.p3d_top(verts, p, j, 0.25 + 0.5 * u)
                else:
                    res, p, u = self.r_t1.intersect(ntr)
                    self.p3d_top(verts, p, j, 0.75 + 0.25 * u)                
            else:
                # one edge
                if tb < 0.5:
                    res, p, u = self.r_t0.intersect(nbr)
                    self.p3d_bottom(verts, p, j, 0.5 * u)
                else:
                    res, p, u = self.r_t1.intersect(nbr)
                    self.p3d_bottom(verts, p, j, 0.5 + 0.5 * u)
                if tt < 0.5:
                    res, p, u = self.r_t0.intersect(ntr)
                    self.p3d_top(verts, p, j, 0.5 * u)
                else:
                    res, p, u = self.r_t1.intersect(ntr)
                    self.p3d_top(verts, p, j, 0.5 + 0.5 * u)
        # make edges verts after regular ones
        if self.edges_multiples:
            # edge 1
            if tb < 0.25 and tb + self.t_step > 0.25:     
                if self.l_shape != 'CIRCLE': 
                    e_left = True
                    # the step goes through the edge
                    self.p3d_corner(verts, self.l_tc.p, j, 0.25)
                    self.p3d_top(verts, self.l_tc.p, j, 0.25)                
                if self.r_shape != 'CIRCLE':
                    e_right = True
                    # the step goes through the edge
                    self.p3d_corner(verts,  self.r_tc.p, j, 0.25)
                    self.p3d_top(verts,  self.r_tc.p, j, 0.25)
            # edge 2
            if tb < 0.75 and tb + self.t_step > 0.75:     
                if self.l_shape != 'CIRCLE': 
                    e_left = True
                    # the step goes through the edge
                    self.p3d_corner(verts, self.l_t1.p, j, 0.75)
                    self.p3d_top(verts, self.l_t1.p, j, 0.75)                
                if self.r_shape != 'CIRCLE':
                    e_right = True
                    # the step goes through the edge
                    self.p3d_corner(verts,  self.r_t1.p, j, 0.75)
                    self.p3d_top(verts,  self.r_t1.p, j, 0.75)
        else:        
            if tb < 0.5 and tb + self.t_step > 0.5:     
                if self.l_shape != 'CIRCLE': 
                    e_left = True
                    # the step goes through the edge
                    self.p3d_corner(verts, self.l_t1.p, j, 0.5)
                    self.p3d_top(verts, self.l_t1.p, j, 0.5)                
                if self.r_shape != 'CIRCLE':
                    e_right = True
                    # the step goes through the edge
                    self.p3d_corner(verts,  self.r_t1.p, j, 0.5)
                    self.p3d_top(verts,  self.r_t1.p, j, 0.5)
        
        self.make_faces(i, f, faces, e_right, e_left)
    def straight_palier(self, length):
        return Stair.straight_palier(self, length, last_type='PALIER')
    def curved_palier(self, da, radius, left_shape, right_shape, double_limit=pi):
        return Stair.curved_palier(self, da, radius, left_shape, right_shape, double_limit=double_limit,  last_type='PALIER')
        
class StairGenerator():
    def __init__(self):
        self.last_type = 'NONE'
        self.stairs = []
        self.faces_type = 'NONE'
    def add_part(self, type, faces_type, z_mode, step_z, bottom_z, center, radius, da, width_left, width_right, length, left_shape, right_shape):
        self.faces_type = faces_type
        if len(self.stairs) < 1:
            s = None
        else:
            s = self.stairs[-1]
        # start a new stair    
        if s is None:
            if type == 'S_STAIR':
                p = Vector((0,0))
                v = Vector((0,length))
                s = StraightStair(p, v, width_left, width_right, faces_type, z_mode, step_z, bottom_z)
            elif type == 'C_STAIR':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedStair(c, radius, 0, da, width_left, width_right, faces_type, z_mode, step_z, bottom_z, left_shape, right_shape)
            elif type == 'D_STAIR':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedStair(c, radius, 0, da, width_left, width_right, faces_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=0)
            elif type == 'S_PALIER':
                p = Vector((0,0))
                v = Vector((0,length))
                s = StraightPalier(p, v, width_left, width_right, faces_type, z_mode, step_z, bottom_z)
            elif type == 'C_PALIER':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedPalier(c, radius, 0, da, width_left, width_right, faces_type, z_mode, step_z, bottom_z, left_shape, right_shape)
            elif type == 'D_PALIER':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedPalier(c, radius, 0, da, width_left, width_right, faces_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=0)
        else:
            if type == 'S_STAIR':
                s = s.straight_stair(length)
            elif type == 'C_STAIR':
                s = s.curved_stair(da, radius, left_shape, right_shape)
            elif type == 'D_STAIR':
                s = s.curved_stair(da, radius, left_shape, right_shape, double_limit=0)
            elif type == 'S_PALIER':
                s = s.straight_palier(length)
            elif type == 'C_PALIER':
                s = s.curved_palier(da, radius, left_shape, right_shape)
            elif type == 'D_PALIER':
                s = s.curved_palier(da, radius, left_shape, right_shape, double_limit=0)
        self.stairs.append(s)
        self.last_type = type
    def n_steps(self, step_depth):
        n_steps = 0
        for stair in self.stairs:
            n_steps += stair.step_size(step_depth)
        return n_steps
    def set_height(self, step_height):
        z = 0
        for stair in self.stairs:
            stair.set_height(step_height, z)
            z = stair.top
    def make_stair(self, height, step_depth, verts, faces, t_offset=0):
        n_steps = self.n_steps(step_depth)
        self.set_height(height / n_steps)
        if self.faces_type in ['CLOSED', 'FULL'] and len(self.stairs) > 0:
            faces.append((0,1,5,4))
        for stair in self.stairs:
            for i in range(stair.n_step+1):
                t_o = stair.top_offset * t_offset
                stair.make_step(i, verts, faces, t_offset=-t_o)
  
def update(self, context):
    self.update(context)

def update_preset(self, context):
    self.auto_update = False
    if self.presets == 'STAIR_I':
        self.n_parts = 1
        self.update_parts()
        self.parts[0].type = 'S_STAIR'
    elif self.presets == 'STAIR_L':
        self.n_parts = 3
        self.update_parts()
        self.parts[0].type = 'S_STAIR'
        self.parts[1].type = 'C_STAIR'
        self.parts[2].type = 'S_STAIR'
        self.da = pi / 2
    elif self.presets == 'STAIR_U':
        self.n_parts = 3
        self.update_parts()
        self.parts[0].type = 'S_STAIR'
        self.parts[1].type = 'D_STAIR'
        self.parts[2].type = 'S_STAIR'
        self.da = pi
    elif self.presets == 'STAIR_O':
        self.n_parts = 2
        self.update_parts()
        self.parts[0].type = 'D_STAIR'
        self.parts[1].type = 'D_STAIR'
        self.da = pi
    self.auto_update = True
    self.update(context)
    
class StairPartProperty(PropertyGroup):
    type = EnumProperty(
            items=(
                ('S_STAIR', 'Straight stair','',0),
                ('C_STAIR', 'Curved stair','',1),
                ('D_STAIR', 'Dual Curved stair','',2),
                ('S_PALIER', 'Straight palier','',3),
                ('C_PALIER', 'Curved palier','',4),
                ('D_PALIER', 'Dual Curved palier','',5)
                ),
            default='S_STAIR',
            update=update
            )
    length = FloatProperty(
            name="length",
            min=0.5,
            max=100.0,
            default=2.0,
            update=update
            )
    radius = FloatProperty(
            name="radius",
            min=0.5,
            max=100.0,
            default=0.7,
            update=update
            )
    da = FloatProperty(
            name="angle",
            min=-pi,
            max=pi,
            default=pi/2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    left_shape = EnumProperty(
            items=(
                ('RECTANGLE', 'Straight','',0),
                ('CIRCLE', 'Curved ','',1)
                ),
            default='RECTANGLE',
            update=update
            )
    right_shape = EnumProperty(
            items=(
                ('RECTANGLE', 'Straight','',0),
                ('CIRCLE', 'Curved ','',1)
                ),
            default='RECTANGLE',
            update=update
            )
    
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = ARCHIPACK_PT_stair.params(o)
            if props:
                for part in props.parts:
                    if part == self:
                        return props
        return None
        
    def update(self, context):
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context)
        
    def draw(self, layout, context, index, user_mode):
        row = layout.row(align=True)
        row.label(text="Part "+str(index+1))
        if user_mode:
            row.prop(self, "type", text="")
        if self.type in ['C_STAIR', 'C_PALIER', 'D_STAIR', 'D_PALIER']: 
            if user_mode:
                row = layout.row()
                row.prop(self, "radius")
                row = layout.row()
                row.prop(self, "da")
        else:
            row = layout.row()
            row.prop(self, "length")
        if user_mode:
            row = layout.row(align=True)
            row.prop(self, "left_shape", text="")
            row.prop(self, "right_shape", text="")
        
class StairProperty(PropertyGroup):
    parts = CollectionProperty(type=StairPartProperty)
    n_parts = IntProperty(
            name="parts",
            min=1,
            max=32,
            default=1, update=update
            )
    step_depth = FloatProperty(
            name="step depth",
            min=0.2,
            max=2.0,
            default=0.25,
            update=update
            )
    width = FloatProperty(
            name="width",
            min=0.01,
            max=100.0,
            default=1.2,
            update=update
            )
    height= FloatProperty(
            name="Height",
            min=0.1, max=1000,
            default=2.4, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )    
    t_offset = FloatProperty(
        name="Top offset",
            min=0.0, max=1000,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )    
    x_offset = FloatProperty(
        name="x offset",
            min=-1000, max=1000,
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )    
    step_z = FloatProperty(
            name="Depth",
            min=0.001, max=1000,
            default=0.03, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )    
    bottom_z = FloatProperty(
            name="Stair bottom",
            min=0.001, max=1000,
            default=0.03, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    radius = FloatProperty(
            name="radius",
            min=0.5,
            max=100.0,
            default=0.7,
            update=update
            )
    da = FloatProperty(
            name="angle",
            min=-pi,
            max=pi,
            default=pi/2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    faces_type = EnumProperty(
            items=(
                ('CLOSED', 'Closed','',0),
                ('FULL', 'Full height','',1),
                ('OPEN', 'Open ','',2)
                ),
            default='CLOSED',
            update=update
            )
    left_shape = EnumProperty(
            items=(
                ('RECTANGLE', 'Straight','',0),
                ('CIRCLE', 'Curved ','',1)
                ),
            default='RECTANGLE',
            update=update
            )
    right_shape = EnumProperty(
            items=(
                ('RECTANGLE', 'Straight','',0),
                ('CIRCLE', 'Curved ','',1)
                ),
            default='RECTANGLE',
            update=update
            )
    z_mode = EnumProperty(
            items=(
                ('STANDARD', 'Standard','',0),
                ('LINEAR', 'Bottom Linear','',1),
                ('LINEAR_TOP', 'All Linear','',2)
                ),
            default='STANDARD',
            update=update
            )
    presets = EnumProperty(
        items=(
            ('STAIR_I','I stair', '', 0),
            ('STAIR_L','L stair', '', 1),
            ('STAIR_U','U stair', '', 2),
            ('STAIR_O','O stair', '', 3),
            ('STAIR_USER','User defined stair', '', 4),
            ),
        default='STAIR_I', update=update_preset
        )
    auto_update=BoolProperty(default=True)
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        active = context.active_object
        selected = [o for o in context.selected_objects]
        #selected = [o for o in context.scene.objects]
        for o in selected:
            if ARCHIPACK_PT_stair.params(o) == self:
                return active, selected, o
        return active, selected, None
    
    def update_parts(self):
        # remove rows
        for i in range(len(self.parts), self.n_parts, -1):
            self.parts.remove(i-1)
        # add rows
        for i in range(len(self.parts), self.n_parts):
            self.parts.add()
        
    def update(self, context):   
        
        if not self.auto_update:
            return
        active, selected, o = self.find_in_selection(context)
        
        if o is None:
            return
            
        self.update_parts()
        
        center = Vector((0,0))
        verts = []
        faces = []
        
        # depth at bottom
        bottom_z = self.bottom_z
        if self.faces_type == 'OPEN':
            # depth at front
            bottom_z = self.step_z
        
        width_left  = 0.5*self.width + self.x_offset
        width_right = 0.5*self.width - self.x_offset   
        
        g = StairGenerator()
        if self.presets == 'STAIR_USER':
            for part in self.parts:
                g.add_part(part.type, self.faces_type, self.z_mode, self.step_z, bottom_z, center, max(width_left+0.01, width_right+0.01, part.radius), part.da, width_left, width_right, part.length, part.right_shape, part.left_shape)
        else:
            for part in self.parts:
                g.add_part(part.type, self.faces_type, self.z_mode, self.step_z, bottom_z, center, max(width_left+0.01, width_right+0.01, self.radius), self.da,  width_left, width_right, part.length, self.right_shape, self.left_shape)
        g.make_stair(self.height, self.step_depth, verts, faces, t_offset=self.t_offset)
        
        
        
        bmed.buildmesh(context, o, verts, faces, matids=None, uvs=None, weld=True)
        
        # restore context
        try:
            for o in selected:
                o.select = True
        except:
            pass    
        
        active.select = True
        context.scene.objects.active = active  
 
class ARCHIPACK_PT_stair(Panel):
    bl_idname = "ARCHIPACK_PT_stair"
    bl_label = "Stair"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'
    
    def draw(self, context):
        o = context.object
        prop = ARCHIPACK_PT_stair.params(o)
        if prop is None:
            return 
        layout = self.layout
        row = layout.row(align=True)
        row.prop(prop, 'presets')
        box = layout.box()
        box.prop(prop, 'faces_type')
        box.prop(prop, 'n_parts')
        box.prop(prop, 'step_depth')
        box.prop(prop, 'width')
        box.prop(prop, 'x_offset')
        box.prop(prop, 'height')
        box.prop(prop, 'step_z')
        box.prop(prop, 'z_mode')
        box.prop(prop, 'bottom_z')
        box.prop(prop, 't_offset')
        if prop.presets != 'STAIR_USER':
            row = layout.row()
            row.prop(prop, "radius")
            row = layout.row()
            row.prop(prop, "da")                
            row = layout.row(align=True)
            row.prop(prop, "left_shape", text="")
            row.prop(prop, "right_shape", text="")
        for i, part in enumerate(prop.parts):
            box = layout.box()
            part.draw(box, context, i, prop.presets == 'STAIR_USER')
        
    @classmethod
    def params(cls, o):
        try:
            if 'StairProperty' not in o.data:
                return False
            else:
                return o.data.StairProperty[0]
        except:
            return False
    @classmethod
    def filter(cls, o):
        try:
            if 'StairProperty' not in o.data:
                return False
            else:
                return True
        except:
            return False
    @classmethod
    def poll(cls, context):
        o = context.object
        if o is None:
            return False
        return cls.filter(o)
 
# ------------------------------------------------------------------
# Define operator class to create object
# ------------------------------------------------------------------
class ARCHIPACK_OT_stair(Operator):
    bl_idname = "archipack.stair"
    bl_label = "Stair (alpha)"
    bl_description = "Stair (alpha)"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    
    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')
        
    def create(self, context):
        m = bpy.data.meshes.new("Stair")
        o = bpy.data.objects.new("Stair", m)
        d = m.StairProperty.add()
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o 
        d.update(context)
        #MaterialUtils.add_stair_materials(o)
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
            o = self.create(context)
            o.location = bpy.context.scene.cursor_location
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}
            
bpy.utils.register_class(StairPartProperty)
bpy.utils.register_class(StairProperty)
Mesh.StairProperty = CollectionProperty(type=StairProperty)
bpy.utils.register_class(ARCHIPACK_PT_stair)
bpy.utils.register_class(ARCHIPACK_OT_stair)



