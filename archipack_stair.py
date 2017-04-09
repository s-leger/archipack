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
from bpy.props import FloatProperty, BoolProperty, IntProperty, CollectionProperty , EnumProperty, FloatVectorProperty, IntVectorProperty
from .bmesh_utils import BmeshEdit as bmed
from .materialutils import MaterialUtils
from mathutils import Vector
from math import sin, cos, pi, atan2, sqrt, floor

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
        steps = max(1, round(self.length / len, 0))
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
    def __init__(self, left_offset, right_offset, steps_type, z_mode, step_z, bottom_z):
        self.steps_type = steps_type
        self.l_shape = None 
        self.r_shape = None
        self.next_type = 'NONE'
        self.last_type = 'NONE'
        # verts mode in STAIR LADDER POST
        self.verts_mode = 'STAIR' 
        self.z_mode = z_mode
        # depth of open step
        self.step_z = step_z
        # size under the step on bottom
        self.bottom_z = bottom_z
        self.left_offset = left_offset
        self.right_offset = right_offset
        self.last_height = 0
    def set_matids(self, matids):
        self.idmat_top, self.idmat_face_top, self.idmat_face_bottom, self.idmat_side, self.idmat_bottom = matids 
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
        if self.z_mode == 'LINEAR_TOP':
            z = self.z0 + t * self.height
        else:
            z = self.z0 + i * self.step_height
        if 'LINEAR' in self.z_mode:
            zb = self.z0 + t * self.height
        else:    
            zb = z
        if self.steps_type == 'FULL':
            zb = 0
        else:
            zb = max(0, zb - self.bottom_z)
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
        if self.z_mode == 'LINEAR_TOP':
            z3 = self.z0 + t * self.height
        else: 
            z3 = self.z0 + self.step_height
        z2 = z3 - self.bottom_z
        
        z1 = self.z0 - self.bottom_z + 0.5 * self.step_height
        if self.z_mode == 'LINEAR_TOP':
            z0 = self.z0 + t * self.height
        else:
            z0 = self.z0 
        if self.steps_type == 'CLOSED':
            z0 = max(0, z0 - self.bottom_z)
        else:
            z0 = 0
        verts.append((x, y, z0))
        verts.append((x, y, z1))
        verts.append((x, y, z2))
        verts.append((x, y, z3))
    def straight_stair(self, length):
        self.next_type = 'STAIR'
        p, v = self.straight(length)
        return StraightStair(p, v, self.left_offset, self.right_offset, self.steps_type, self.z_mode, self.step_z, self.bottom_z)
    def straight_palier(self, length, last_type='STAIR'):
        self.next_type = 'PALIER'
        p, v = self.straight(length)
        return StraightPalier(p, v, self.left_offset, self.right_offset, self.steps_type, self.z_mode, self.step_z, self.bottom_z, last_type=last_type)
    def curved_stair(self, da, radius, left_shape, right_shape, double_limit=pi):
        self.next_type = 'STAIR'
        n = self.normal(1)
        n.v = radius * n.v.normalized()
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedStair(c, radius, a0, da, self.left_offset, self.right_offset, self.steps_type, self.z_mode, self.step_z, self.bottom_z, left_shape, right_shape, double_limit=double_limit)
    def curved_palier(self, da, radius, left_shape, right_shape, double_limit=pi, last_type='STAIR'):
        self.next_type = 'PALIER'
        n = self.normal(1)
        n.v = radius * n.v.normalized()
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedPalier(c, radius, a0, da, self.left_offset, self.right_offset, self.steps_type, self.z_mode, self.step_z, self.bottom_z, left_shape, right_shape, double_limit=double_limit, last_type=last_type)
    def get_z(self, t, mode):
        if mode == 'LINEAR':
            return self.z0 + t * self.height
        else:
            step = 1+floor(t/self.t_step)
            return self.z0 + step * self.step_height
         
class StraightStair(Stair, Line):
    def __init__(self, p, v, left_offset, right_offset, steps_type, z_mode, step_z, bottom_z):
        Stair.__init__(self, left_offset, right_offset, steps_type, z_mode, step_z, bottom_z)
        Line.__init__(self, p, v)
        self.l_line = self.offset(left_offset)
        self.r_line = self.offset(-right_offset)
    def make_step(self, i, verts, faces, matids, t_offset=0):
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
            if self.z_mode != 'LINEAR_TOP':
                # contre marche haut
                faces.append((f+2, f+3, f+7, f+6))
                matids.append(self.idmat_face_top)
            # dessus marche rectangle
            faces.append((f+3, f+9, f+13, f+7))
            matids.append(self.idmat_top)
            if self.steps_type in ['CLOSED', 'FULL']:
                if self.z_mode != 'LINEAR_TOP':
                    # contre-marche bas
                    faces.append((f+1, f+2, f+6, f+5))
                    matids.append(self.idmat_face_bottom)
                # dessous marche rectangle
                faces.append((f, f+4, f+12, f+8 ))  
                matids.append(self.idmat_bottom)
                # side right top
                faces.append((f+1, f+9, f+3, f+2))
                matids.append(self.idmat_side)
                # side left top
                faces.append((f+5, f+6, f+7, f+13))
                matids.append(self.idmat_side)
                # side left bottom
                faces.append((f+4, f+5, f+13, f+12))
                matids.append(self.idmat_side)
                # side right bottom
                faces.append((f, f+8, f+9, f+1))
                matids.append(self.idmat_side)
            else:
                # side left top
                faces.append((f+6, f+7, f+13, f+12))
                matids.append(self.idmat_side)
                # side right top
                faces.append((f+8, f+9, f+3, f+2))
                matids.append(self.idmat_side)
                # dessous marche rectangle
                faces.append((f+6, f+12, f+8, f+2))
                matids.append(self.idmat_bottom)
                if self.z_mode != 'LINEAR_TOP':
                    # arriere
                    faces.append((f+8, f+12, f+13, f+9))
                    matids.append(self.idmat_bottom)
        elif self.steps_type in ['CLOSED', 'FULL'] and self.z_mode != 'LINEAR_TOP':
            faces.append((f, f+4, f+5, f+1))
            matids.append(self.idmat_bottom)       
    def get_post(self, posts, side, i, t_step, z_offset, respect_edges):
        t = i * t_step
        if side == 'LEFT':
            n = self.l_line.normal(t)
        else:
            n = self.r_line.normal(t)
        z0 = self.get_z(t, 'STEP')
        z1 = z_offset+self.get_z(t, 'LINEAR')
        posts.append((n, z0, z1))
    def n_posts(self, post_spacing, side, respect_edges):
        return self.steps(post_spacing)

class CurvedStair(Stair, Arc):
    def __init__(self, c, radius, a0, da, left_offset, right_offset, steps_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=pi):
        Stair.__init__(self, left_offset, right_offset, steps_type, z_mode, step_z, bottom_z)
        Arc.__init__(self, c, radius, a0, da)
        self.l_shape = left_shape 
        self.r_shape = right_shape
        self.edges_multiples = round(abs(da),6) > double_limit
        #left arc, tangeant at start and end
        self.l_arc, self.l_t0, self.l_t1, self.l_tc = self.set_offset(left_offset, left_shape)
        self.r_arc, self.r_t0, self.r_t1, self.r_tc = self.set_offset(-right_offset, right_shape)
    def set_offset(self, offset, shape):
        arc = self.offset(offset)
        t0, t1, tc = None, None, None
        if shape != 'CIRCLE':
            t0 = arc.tangeant(0,1)
            t1 = arc.tangeant(1,1)
            if self.edges_multiples: 
                tc = arc.tangeant(0.5,1)
                i, p, t = t0.intersect(tc)
                tc.v *= 2*t
                tc.p = p
                i, p, t2 = tc.intersect(t1)
            else:
                i, p, t = t0.intersect(t1)
            t0.v *= t
            t1.p = p
            t1.v *= t
        return arc, t0, t1, tc
    def p3d_corner(self, verts, p2d, i, t):
        x, y = p2d
        if self.steps_type == 'FULL':
            z = 0
        else:
            if 'LINEAR'in self.z_mode:
                z = self.z0 + t * self.height - self.bottom_z
            else:    
                z = self.z0 + (i+0.5) * self.step_height - self.bottom_z
        verts.append((x, y, z))
    def make_faces(self, i, f, faces, matids, e_right, e_left):
        if i < self.n_step:
            if self.z_mode != 'LINEAR_TOP':
                # self.idmat_top, self.idmat_face_top, self.idmat_face_bottom, self.idmat_side, self.idmat_bottom
                # contre marche haut
                faces.append((f+2, f+3, f+7, f+6))
                matids.append(self.idmat_face_top)
            if self.steps_type in ['CLOSED', 'FULL']:
                if self.z_mode != 'LINEAR_TOP':
                    # contre-marche bas
                    faces.append((f+1, f+2, f+6, f+5))
                    matids.append(self.idmat_face_bottom)
                if e_right and e_left:
                    # dessus faces with two more edges
                    faces.append((f+3, f+10, f+15, f+19, f+13, f+7))
                    matids.append(self.idmat_top)
                    # side left bottom haut
                    faces.append((f+5, f+6, f+7, f+13, f+12))
                    matids.append(self.idmat_side)
                    # side left bottom 
                    faces.append((f+4, f+5, f+12, f+11))
                    matids.append(self.idmat_side)
                    # side left top haut
                    faces.append((f+12, f+13, f+19))
                    matids.append(self.idmat_side)
                    # side left top bas
                    faces.append((f+11, f+12, f+19, f+18))
                    matids.append(self.idmat_side)
                    # side right top bas
                    faces.append((f+8, f+14, f+15, f+9))
                    matids.append(self.idmat_side)
                    # side right top haut
                    faces.append((f+9, f+15, f+10))
                    matids.append(self.idmat_side)
                    # side right bottom bas
                    faces.append((f, f+8, f+9, f+1))
                    matids.append(self.idmat_side)
                    # side right bottom haut
                    faces.append((f+1, f+9, f+10, f+3, f+2))
                    matids.append(self.idmat_side)
                    # bottom
                    faces.append((f, f+4, f+11, f+8)) 
                    matids.append(self.idmat_bottom)
                    faces.append((f+11, f+18, f+14, f+8)) 
                    matids.append(self.idmat_bottom)
                elif e_right:
                    #dessus right with one more edge
                    faces.append((f+3, f+12, f+16, f+10, f+7))
                    matids.append(self.idmat_top)
                    # side right haut
                    faces.append((f+1, f+12, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side right bas
                    faces.append((f, f+11, f+12, f+1))
                    matids.append(self.idmat_side)
                    # side left bottom bas
                    faces.append((f+4, f+5, f+9, f+8))
                    matids.append(self.idmat_side)
                    # side left bottom haut
                    faces.append((f+5, f+6, f+7, f+10, f+9))
                    matids.append(self.idmat_side)
                    # side left top haut
                    faces.append((f+9, f+10, f+16))
                    matids.append(self.idmat_side)
                    # side left top bas
                    faces.append((f+15, f+8, f+9, f+16))
                    matids.append(self.idmat_side)
                    # bottom
                    faces.append((f, f+4, f+8, f+15, f+11)) 
                    matids.append(self.idmat_bottom)
                elif e_left:
                    #dessus  left with one more edge 
                    faces.append((f+3, f+10, f+12, f+16, f+7))
                    matids.append(self.idmat_top)
                    # side left haut
                    faces.append((f+5, f+6, f+7, f+16))
                    matids.append(self.idmat_side)
                    # side left bas
                    faces.append((f+4, f+5, f+16, f+15))
                    matids.append(self.idmat_side)
                    # side right top haut
                    faces.append((f+9, f+12, f+10))
                    matids.append(self.idmat_side)
                    # side right top bas
                    faces.append((f+8, f+11, f+12, f+9))
                    matids.append(self.idmat_side)
                    # side right bottom bas
                    faces.append((f, f+8, f+9, f+1))
                    matids.append(self.idmat_side)
                    # side right bottom haut
                    faces.append((f+1, f+9, f+10, f+3, f+2))
                    matids.append(self.idmat_side)
                    # bottom
                    faces.append((f, f+4, f+15, f+11, f+8)) 
                    matids.append(self.idmat_bottom)
                else:
                    # dessus marche rectangle
                    faces.append((f+3, f+9, f+13, f+7))
                    matids.append(self.idmat_top)
                    # dessous marche rectangle
                    faces.append((f, f+4, f+12, f+8 )) 
                    matids.append(self.idmat_bottom)
                    # side right top
                    faces.append((f+1, f+9, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side left top
                    faces.append((f+5, f+6, f+7, f+13))
                    matids.append(self.idmat_side)
                    # side right bottom
                    faces.append((f, f+8, f+9, f+1))
                    matids.append(self.idmat_side)
                    # side left bottom
                    faces.append((f+4, f+5, f+13, f+12))
                    matids.append(self.idmat_side)
            else:
                if e_right and e_left:
                    # dessus faces with two more edges
                    faces.append((f+3, f+10, f+15, f+19, f+13, f+7))
                    matids.append(self.idmat_top)
                    # dessous faces with two more edges
                    faces.append((f+6, f+12, f+18, f+14, f+9, f+2))
                    matids.append(self.idmat_bottom)
                    # side left bottom haut
                    faces.append((f+6, f+7, f+13, f+12))
                    matids.append(self.idmat_side)
                    # side left top haut
                    faces.append((f+12, f+13, f+19, f+18))
                    matids.append(self.idmat_side)
                    # side right bottom haut
                    faces.append((f+2, f+9, f+10, f+3))
                    matids.append(self.idmat_side)
                    # side right top haut
                    faces.append((f+9, f+14, f+15, f+10))  
                    matids.append(self.idmat_side)                  
                elif e_right:
                    #dessus right with one more edge
                    faces.append((f+3, f+12, f+16, f+10, f+7))
                    matids.append(self.idmat_top)
                    #dessous right with one more edge
                    faces.append((f+6, f+9, f+15, f+11, f+2))
                    matids.append(self.idmat_bottom)
                    # side right haut
                    faces.append((f+11, f+12, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side left bottom haut
                    faces.append((f+6, f+7, f+10, f+9))
                    matids.append(self.idmat_side)
                    # side left top haut
                    faces.append((f+9, f+10, f+16, f+15))
                    matids.append(self.idmat_side)
                elif e_left:
                    #dessus  left with one more edge 
                    faces.append((f+3, f+10, f+12, f+16, f+7))
                    matids.append(self.idmat_top)
                    #dessous  left with one more edge 
                    faces.append((f+6, f+15, f+11, f+9, f+2))
                    matids.append(self.idmat_bottom)
                    # side right top haut
                    faces.append((f+11, f+12, f+10, f+9))
                    matids.append(self.idmat_side)
                    # side right bottom haut
                    faces.append((f+9, f+10, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side left haut
                    faces.append((f+6, f+7, f+16, f+15))
                    matids.append(self.idmat_side)                
                else:
                    # dessus marche rectangle
                    faces.append((f+3, f+9, f+13, f+7))
                    matids.append(self.idmat_top)
                    # dessous marche rectangle
                    faces.append((f+6, f+12, f+8, f+2))
                    matids.append(self.idmat_bottom)
                    # side right top
                    faces.append((f+8, f+9, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side left top
                    faces.append((f+6, f+7, f+13, f+12))
                    matids.append(self.idmat_side)
                    # arriere
                    faces.append((f+8, f+12, f+13, f+9))
                    matids.append(self.idmat_side)
        elif self.steps_type in ['CLOSED', 'FULL'] and self.z_mode != 'LINEAR_TOP':
            # close last face
            faces.append((f, f+4, f+5, f+1))
            matids.append(self.idmat_bottom)
    def make_step(self, i, verts, faces, matids, t_offset=0):
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
        self.make_faces(i, f, faces, matids, e_right, e_left)
    def get_post(self, posts, side, i, t_step, z_offset, respect_edges):
        t  = i * t_step
        # vect y
        if side == 'RIGHT':
            arc = self.r_arc
            shape = self.r_shape
            t0, t1, tc = self.r_t0, self.r_t1, self.r_tc
        else:
            arc = self.l_arc
            shape = self.l_shape
            t0, t1, tc = self.l_t0, self.l_t1, self.l_tc
        
        n = arc.normal(t)
        z0 = self.get_z(t, 'STEP')
        z1 = z_offset+self.get_z(t, 'LINEAR')
        
        if shape != 'CIRCLE':
            if self.edges_multiples:    
                # two edges
                if t < 0.25:
                    n = t0.normal(t*4)
                elif t < 0.75:
                    n = tc.normal((t-0.25)*2)
                else:
                    n = t1.normal((t-0.75)*4)
            else:    
                if t < 0.5:
                    n = t0.normal(t*2)
                else:
                    n = t1.normal((t-0.5)*2)
        posts.append((n, z0, z1))
        if shape != 'CIRCLE' and respect_edges:
            if self.edges_multiples:
                if t < 0.25 and t+t_step > 0.25:
                    z0 = self.get_z(0.25, 'STEP')
                    z1 = z_offset+self.get_z(0.25, 'LINEAR')
                    n = t0.normal(1)
                    posts.append((n, z0, z1))
                elif t < 0.75 and t+t_step > 0.75:
                    z0 = self.get_z(0.75, 'STEP')
                    z1 = z_offset+self.get_z(0.75, 'LINEAR')
                    n = tc.normal(1)
                    posts.append((n, z0, z1))
            elif t < 0.5 and t+t_step > 0.5:
                    z0 = self.get_z(0.5, 'STEP')
                    z1 = z_offset+self.get_z(0.5, 'LINEAR')
                    n = t0.normal(1)
                    posts.append((n, z0, z1))
    def n_posts(self, post_spacing, side, respect_edges):
        if side == 'LEFT':
            arc, t0, shape = self.l_arc, self.l_t0, self.l_shape
        else:
            arc, t0, shape = self.r_arc, self.r_t0, self.r_shape
        step_factor = 1
        if shape == 'CIRCLE':
            length = arc.length
        else:
            if self.edges_multiples: 
                if respect_edges:
                    step_factor = 2
                length = 4*t0.length
            else:
                length = 2*t0.length
        steps = step_factor * max(1, round(length / post_spacing, 0))
        print("respect_edges:%s t_step:%s n_step:%s" % (respect_edges, 1.0 / steps, int(steps)))
        return 1.0 / steps, int(steps)
        
class StraightPalier(StraightStair):
    def __init__(self, p, v, left_offset, right_offset, steps_type, z_mode, step_z, bottom_z, last_type='STAIR'):
        StraightStair.__init__(self, p, v, left_offset, right_offset, steps_type, z_mode, step_z, bottom_z)
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
    def make_step(self, i, verts, faces, matids, t_offset=0):
        f = len(verts)
        j = 0
        if i == self.n_step and self.next_type != 'PALIER' and self.steps_type in ['CLOSED','FULL']:
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
        # self.idmat_top, self.idmat_face_top, self.idmat_face_bottom, self.idmat_side, self.idmat_bottom
                        
        if i < self.n_step:  
            # contre marche haut
            if self.z_mode != 'LINEAR_TOP':
                faces.append((f+2, f+3, f+7, f+6))
                matids.append(self.idmat_face_top)
            # dessus marche rectangle
            faces.append((f+3, f+11, f+15, f+7))
            matids.append(self.idmat_top)
            if self.steps_type in ['CLOSED','FULL']:
                # contre-marche bas
                faces.append((f+1, f+2, f+6, f+5))
                matids.append(self.idmat_face_bottom)
                # side right top
                faces.append((f+1, f+9, f+10, f+11, f+3, f+2))
                matids.append(self.idmat_side)
                # side right bottom
                faces.append((f, f+8, f+9, f+1))
                matids.append(self.idmat_side)
                # side left top
                faces.append((f+5, f+6, f+7, f+15, f+14, f+13))
                matids.append(self.idmat_side)
                # side left bottom
                faces.append((f+4, f+5, f+13, f+12))
                matids.append(self.idmat_side)
                # dessous marche rectangle
                faces.append((f, f+4, f+12, f+8 ))
                matids.append(self.idmat_bottom)
            else:
                # dessous marche rectangle
                faces.append((f+6, f+14, f+10, f+2))
                matids.append(self.idmat_bottom)
                # side right top
                faces.append((f+10, f+11, f+3, f+2))
                matids.append(self.idmat_side)
                # side left top
                faces.append((f+6, f+7, f+15, f+14))
                matids.append(self.idmat_side)
        else:
            if self.steps_type in ['CLOSED','FULL']:
                if self.next_type == 'STAIR':
                    # close last part bottom
                    faces.append((f, f+4, f+5, f+6, f+2,f+1))
                    matids.append(self.idmat_bottom)
                elif self.next_type == 'NONE':
                    faces.append((f, f+4, f+5, f+6, f+7, f+3, f+2, f+1))
                    matids.append(self.idmat_bottom)
    def straight_palier(self, length):
        return Stair.straight_palier(self, length, last_type='PALIER')
    def curved_palier(self, da, radius, left_shape, right_shape, double_limit=pi):
        return Stair.curved_palier(self, da, radius, left_shape, right_shape, double_limit=double_limit, last_type='PALIER')
    def get_z(self, t, mode):
        if mode == 'STEP':
            return self.z0 + self.step_height
        else:
            return self.z0

class CurvedPalier(CurvedStair):
    def __init__(self, c, radius, a0, da, left_offset, right_offset, steps_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=pi, last_type='STAIR'):
        CurvedStair.__init__(self, c, radius, a0, da, left_offset, right_offset, steps_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=double_limit)
        self.last_type=last_type
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
        return 0
    @property
    def top(self):
        if self.next_type == 'PALIER':
            return self.z0
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
        if self.steps_type == 'FULL':
            z = 0
        else:
            z = self.z0 - self.bottom_z
        verts.append((x, y, z))
    def make_faces(self, i, f, faces, matids, e_right, e_left):
        if i < self.n_step:
            if i == 0 and self.last_type != 'PALIER':
                # only on first part
                # disabled between flat parts
                # contre marche haut
                if self.z_mode != 'LINEAR_TOP':
                    faces.append((f+2, f+3, f+7, f+6))
                    matids.append(self.idmat_face_top)
                if self.steps_type in ['CLOSED', 'FULL'] and self.z_mode != 'LINEAR_TOP':
                    # contre-marche bas
                    faces.append((f+1, f+2, f+6, f+5))
                    matids.append(self.idmat_face_bottom)
            if self.steps_type in ['CLOSED', 'FULL']:
                if e_right and e_left:
                    # dessus faces with two more edges
                    faces.append((f+3, f+10, f+17, f+21, f+13, f+7))
                    matids.append(self.idmat_top)
                    # side left bottom bas
                    faces.append((f+4, f+5, f+6, f+12, f+11))
                    matids.append(self.idmat_side)
                    # side left bottom haut
                    faces.append((f+6, f+7, f+13, f+12))
                    matids.append(self.idmat_side)
                    # side left top bas
                    faces.append((f+11, f+12, f+20, f+19, f+18))
                    matids.append(self.idmat_side)
                    # side left top haut
                    faces.append((f+12, f+13, f+21, f+20))
                    matids.append(self.idmat_side)
                    # side right bottom bas
                    faces.append((f, f+8, f+9, f+2, f+1))
                    matids.append(self.idmat_side)
                    # side right bottom haut
                    faces.append((f+2, f+9, f+10, f+3))
                    matids.append(self.idmat_side)
                    # side right top bas
                    faces.append((f+8, f+14, f+15, f+16, f+9))
                    matids.append(self.idmat_side)
                    # side right top haut
                    faces.append((f+10, f+9, f+16, f+17 ))
                    matids.append(self.idmat_side)
                    # bottom
                    faces.append((f, f+4, f+11, f+8))
                    matids.append(self.idmat_bottom)
                    faces.append((f+11, f+18, f+14, f+8))
                    matids.append(self.idmat_bottom)
                elif e_right:
                    #dessus right with one more edge
                    faces.append((f+3, f+14, f+18, f+10, f+7))
                    matids.append(self.idmat_top)
                    # side right bas
                    faces.append((f, f+11, f+12, f+1))
                    matids.append(self.idmat_side)
                    # side right haut
                    faces.append((f+12, f+13, f+14, f+3, f+2, f+1))
                    matids.append(self.idmat_side)
                    # side left bottom haut
                    faces.append((f+6, f+7, f+10, f+9))
                    matids.append(self.idmat_side)
                    # side left bottom bas
                    faces.append((f+4, f+5, f+6, f+9, f+8))
                    matids.append(self.idmat_side)
                    # side left top bas
                    faces.append((f+8, f+9, f+17, f+16, f+15))
                    matids.append(self.idmat_side)
                    # side left top haut
                    faces.append((f+9, f+10, f+18, f+17))
                    matids.append(self.idmat_side)
                    # bottom
                    faces.append((f, f+4, f+8, f+15, f+11))
                    matids.append(self.idmat_bottom)
                elif e_left:
                    #dessus  left with one more edge 
                    faces.append((f+3, f+10, f+14, f+18, f+7))
                    matids.append(self.idmat_top)
                    # side right bottom bas
                    faces.append((f, f+8, f+9, f+2, f+1))
                    matids.append(self.idmat_side)
                    # side right bottom haut
                    faces.append((f+9, f+10, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side right top bas
                    faces.append((f+8, f+11, f+12, f+13, f+9))
                    matids.append(self.idmat_side)
                    # side right top haut
                    faces.append((f+9, f+13, f+14, f+10))
                    matids.append(self.idmat_side)
                    # side left haut
                    faces.append((f+5, f+6, f+7, f+18, f+17, f+16))
                    matids.append(self.idmat_side)
                    # side left bas
                    faces.append((f+4, f+5, f+16, f+15))
                    matids.append(self.idmat_side)
                    # bottom
                    faces.append((f, f+4, f+15, f+11, f+8))
                    matids.append(self.idmat_bottom)
                else:
                    # dessus marche rectangle
                    faces.append((f+3, f+11, f+15, f+7))
                    matids.append(self.idmat_top)
                    # side right top
                    faces.append((f+1, f+9, f+10, f+11, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side right bottom
                    faces.append((f, f+8, f+9, f+1))
                    matids.append(self.idmat_side)
                    # side left top
                    faces.append((f+5, f+6, f+7, f+15, f+14, f+13))
                    matids.append(self.idmat_side)
                    # side left bottom
                    faces.append((f+4, f+5, f+13, f+12))
                    matids.append(self.idmat_side)
                    # dessous marche rectangle
                    faces.append((f, f+4, f+12, f+8 ))
                    matids.append(self.idmat_bottom)
            else:
                if e_right and e_left:
                    # dessus faces with two more edges
                    faces.append((f+3, f+10, f+17, f+21, f+13, f+7))
                    matids.append(self.idmat_top)
                    # dessous faces with two more edges
                    faces.append((f+6, f+12, f+20, f+16, f+9, f+2))
                    matids.append(self.idmat_bottom)
                    # side left bottom haut
                    faces.append((f+6, f+7, f+13, f+12))
                    matids.append(self.idmat_side)
                    # side left top haut
                    faces.append((f+12, f+13, f+21, f+20))
                    matids.append(self.idmat_side)
                    # side right bottom haut
                    faces.append((f+2, f+9, f+10, f+3))
                    matids.append(self.idmat_side)
                    # side right top haut
                    faces.append((f+10, f+9, f+16, f+17 ))
                    matids.append(self.idmat_side)
                elif e_right:
                    #dessus right with one more edge
                    faces.append((f+3, f+12, f+16, f+10, f+7))
                    matids.append(self.idmat_top)
                    #dessous right with one more edge
                    faces.append((f+6, f+9, f+15, f+11, f+2))
                    matids.append(self.idmat_bottom)
                    # side right haut
                    faces.append((f+11, f+12, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side left bottom haut
                    faces.append((f+6, f+7, f+12, f+11))
                    matids.append(self.idmat_side)
                    # side left top haut
                    faces.append((f+9, f+10, f+16, f+15))
                    matids.append(self.idmat_side)
                elif e_left:
                    #dessus  left with one more edge 
                    faces.append((f+3, f+10, f+14, f+18, f+7))
                    matids.append(self.idmat_top)
                    #dessous  left with one more edge 
                    faces.append((f+6, f+17, f+13, f+9, f+2))
                    matids.append(self.idmat_bottom)
                    # side right bottom haut
                    faces.append((f+9, f+10, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side right top haut
                    faces.append((f+9, f+13, f+14, f+10))
                    matids.append(self.idmat_side)
                    # side left haut
                    faces.append((f+6, f+7, f+18, f+17))
                    matids.append(self.idmat_side)
                else:
                    # dessus marche rectangle
                    faces.append((f+3, f+11, f+15, f+7))
                    matids.append(self.idmat_top)
                    # dessous marche rectangle
                    faces.append((f+6, f+14, f+10, f+2))
                    matids.append(self.idmat_bottom)
                    # side right top
                    faces.append((f+10, f+11, f+3, f+2))
                    matids.append(self.idmat_side)
                    # side left top
                    faces.append((f+6, f+7, f+15, f+14))
                    matids.append(self.idmat_side)
    def make_step(self, i, verts, faces, matids, t_offset=0): 
        f = len(verts)
        j = 0
        if i == self.n_step and self.steps_type in ['CLOSED','FULL']:
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
        
        self.make_faces(i, f, faces, matids, e_right, e_left)
    def straight_palier(self, length):
        return Stair.straight_palier(self, length, last_type='PALIER')
    def curved_palier(self, da, radius, left_shape, right_shape, double_limit=pi):
        return Stair.curved_palier(self, da, radius, left_shape, right_shape, double_limit=double_limit,  last_type='PALIER')
    def get_z(self, t, mode):
        if mode == 'STEP':
            return self.z0 + self.step_height
        else:
            return self.z0

class StairGenerator():
    def __init__(self):
        self.last_type = 'NONE'
        self.stairs = []
        self.steps_type = 'NONE'
    def add_part(self, type, steps_type, z_mode, step_z, bottom_z, center, radius, da, width_left, width_right, length, left_shape, right_shape):
        self.steps_type = steps_type
        if len(self.stairs) < 1:
            s = None
        else:
            s = self.stairs[-1]
        # start a new stair    
        if s is None:
            if type == 'S_STAIR':
                p = Vector((0,0))
                v = Vector((0,length))
                s = StraightStair(p, v, width_left, width_right, steps_type, z_mode, step_z, bottom_z)
            elif type == 'C_STAIR':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedStair(c, radius, 0, da, width_left, width_right, steps_type, z_mode, step_z, bottom_z, left_shape, right_shape)
            elif type == 'D_STAIR':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedStair(c, radius, 0, da, width_left, width_right, steps_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=0)
            elif type == 'S_PALIER':
                p = Vector((0,0))
                v = Vector((0,length))
                s = StraightPalier(p, v, width_left, width_right, steps_type, z_mode, step_z, bottom_z)
            elif type == 'C_PALIER':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedPalier(c, radius, 0, da, width_left, width_right, steps_type, z_mode, step_z, bottom_z, left_shape, right_shape)
            elif type == 'D_PALIER':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedPalier(c, radius, 0, da, width_left, width_right, steps_type, z_mode, step_z, bottom_z, left_shape, right_shape, double_limit=0)
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
    def make_stair(self, height, step_depth, verts, faces, matids, t_offset=0):
        n_steps = self.n_steps(step_depth)
        self.set_height(height / n_steps)
        if self.steps_type in ['CLOSED', 'FULL'] and len(self.stairs) > 0:
            faces.append((0,1,5,4))
            matids.append(self.stairs[-1].idmat_face_bottom)
        for stair in self.stairs:
            for i in range(stair.n_step+1):
                t_o = stair.top_offset * t_offset
                stair.make_step(i, verts, faces, matids, t_offset=-t_o)
    def get_post(self, post, post_x, post_y, id_mat, verts, faces, matids):
        n, z0, z1 = post
        vn = n.v.normalized()
        dx = post_x * Vector((vn.y, -vn.x))
        dy = post_y * vn
        x0, y0 = n.p + dx + dy
        x1, y1 = n.p + dx - dy 
        x2, y2 = n.p - dx - dy
        x3, y3 = n.p - dx + dy 
        f = len(verts)
        verts.extend([(x0, y0, z0),(x0, y0, z1),
                    (x1, y1, z0),(x1, y1, z1),
                    (x2, y2, z0),(x2, y2, z1),
                    (x3, y3, z0),(x3, y3, z1)])
        faces.extend([(f, f+1, f+3, f+2),
                    (f+2, f+3, f+5, f+4),
                    (f+4, f+5, f+7, f+6),
                    (f+6, f+7, f+1, f),
                    (f, f+2, f+4, f+6),
                    (f+7, f+5, f+3, f+1)])
        matids.extend([id_mat, id_mat, id_mat, id_mat, id_mat, id_mat])
    def get_panel(self, n0, n1, dist, z0, z1, panel_x, panel_z, id_mat, verts, faces, matids):
        dp = n1.p-n0.p
        c0 = dp.normalized()
        v0 = Vector((c0.y, -c0.x))
        dx = dist * c0
        dy = panel_x * v0
        dz = (z1-z0)/dp.length*dist
        x0, y0 = n0.p + dx + dy
        x1, y1 = n1.p - dx + dy
        x2, y2 = n1.p - dx - dy 
        x3, y3 = n0.p + dx - dy
        z0, z1 = z0+dz, z1-dz 
        z2, z3 = z0+panel_z, z1+panel_z
        f = len(verts)
        verts.extend([(x0, y0, z0),(x0, y0, z2),
                    (x3, y3, z0),(x3, y3, z2),
                    (x2, y2, z1),(x2, y2, z3),
                    (x1, y1, z1),(x1, y1, z3)])
        faces.extend([(f, f+1, f+3, f+2),
                    (f+2, f+3, f+5, f+4),
                    (f+4, f+5, f+7, f+6),
                    (f+6, f+7, f+1, f),
                    (f, f+2, f+4, f+6),
                    (f+7, f+5, f+3, f+1)])
        matids.extend([id_mat, id_mat, id_mat, id_mat, id_mat, id_mat])
    def make_post(self, height, step_depth, post_x, post_y, post_spacing, respect_edges, x_offset, z_offset, enable_r, enable_l, id_mat, enable_post, 
        enable_panel, panel_alt, panel_x, panel_z, panel_dist, panel_mat, verts, faces, matids):
        n_steps = self.n_steps(step_depth)
        self.set_height(height / n_steps)
        l_posts = []
        r_posts = []
        n_stairs = len(self.stairs)-1
        for s, stair in enumerate(self.stairs):
            if type(stair).__name__ in ['CurvedStair', 'CurvedPalier']:
                stair.l_arc, stair.l_t0, stair.l_t1, stair.l_tc = stair.set_offset(x_offset, stair.l_shape)
                stair.r_arc, stair.r_t0, stair.r_t1, stair.r_tc = stair.set_offset(-x_offset, stair.r_shape)
            else:
                stair.l_line = stair.offset(x_offset)
                stair.r_line = stair.offset(-x_offset)
            if enable_l:
                t_step, n_step = stair.n_posts(post_spacing, 'LEFT', respect_edges)
                if s == n_stairs:
                    n_step += 1
                for i in range(n_step):
                    stair.get_post(l_posts, 'LEFT', i, t_step, z_offset, respect_edges)
                    if s == n_stairs and i == n_step-1:
                        n, z0, z1 = l_posts[-1]
                        l_posts[-1] = (n, z0-stair.step_height, z1)
            if enable_r:  
                t_step, n_step = stair.n_posts(post_spacing, 'RIGHT', respect_edges)
                if s == n_stairs:
                    n_step += 1
                for i in range(n_step):
                    stair.get_post(r_posts, 'RIGHT', i, t_step, z_offset, respect_edges)
                    if s == n_stairs and i == n_step-1:
                        n, z0, z1 = r_posts[-1]
                        r_posts[-1] = (n, z0-stair.step_height, z1)
        if enable_post: 
            if enable_l:
                for i, post in enumerate(l_posts):
                    self.get_post(post, post_x, post_y, id_mat, verts, faces, matids)
            if enable_r:    
                for i, post in enumerate(r_posts):
                    self.get_post(post, post_x, post_y, id_mat, verts, faces, matids)
        """
        # panels between posts
        # this is near production ready, still need to be exposed
        """ 
        if enable_panel:
            if enable_l: 
                for i, p0 in enumerate(l_posts):
                    if i < len(l_posts)-1:
                        n0, z, z0 = p0
                        n1, z, z1 = l_posts[i+1]
                        self.get_panel( n0, n1, panel_dist, z0+panel_alt, z1+panel_alt, panel_x, panel_z, panel_mat, verts, faces, matids)
            if enable_r: 
                for i, p0 in enumerate(r_posts):
                    if i < len(r_posts)-1:
                        n0, z, z0 = p0
                        n1, z, z1 = r_posts[i+1]
                        self.get_panel( n0, n1, panel_dist, z0+panel_alt, z1+panel_alt, panel_x, panel_z, panel_mat, verts, faces, matids)
    def make_part(self, height, step_depth, part_x, part_z, x_offset, z_offset, z_mode, steps_type, verts, faces, matids):
        # NOTE: may allow FULL_HEIGHT
        params = [(stair.z_mode, stair.l_shape, stair.r_shape, stair.bottom_z, stair.steps_type) for stair in self.stairs]
        for stair in self.stairs:
            if x_offset < 0:
                stair.l_shape = stair.r_shape
            else:
                stair.r_shape = stair.l_shape
            
            stair.steps_type = steps_type
            stair.z_mode = z_mode 
            stair.bottom_z = part_z
            if type(stair).__name__ in ['CurvedStair', 'CurvedPalier']:
                stair.l_arc, stair.l_t0, stair.l_t1, stair.l_tc = stair.set_offset(x_offset+0.5*part_x, stair.l_shape)
                stair.r_arc, stair.r_t0, stair.r_t1, stair.r_tc = stair.set_offset(x_offset-0.5*part_x, stair.r_shape)
            else:
                stair.l_line = stair.offset(x_offset+0.5*part_x)
                stair.r_line = stair.offset(x_offset-0.5*part_x)
        n_steps = self.n_steps(step_depth)
        self.set_height(height / n_steps)
        for j, stair in enumerate(self.stairs):
            stair.z0 += z_offset+part_z
            stair.n_step *= 2
            stair.t_step /= 2
            stair.step_height /= 2
            for i in range(stair.n_step+1):
                stair.make_step(i, verts, faces, matids, t_offset=0)
            stair.n_step /= 2
            stair.t_step *= 2
            stair.step_height *= 2
            stair.z_mode = params[j][0]
            stair.l_shape = params[j][1]
            stair.r_shape = params[j][2]
            stair.bottom_z = params[j][3]
            stair.steps_type = params[j][4]
            stair.z0 -= z_offset
    def set_matids(self, id_materials):
        for stair in self.stairs:
            stair.set_matids(id_materials)

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
 
materials_enum = (
            ('0','Ceiling','',0),
            ('1','White','',1),
            ('2','Concrete','',2),
            ('3','Wood','',3),
            ('4','Metal','',4),
            ('5','Glass','',5)
            ) 
 
class StairMaterialProperty(PropertyGroup):
    index = EnumProperty(
        items=materials_enum,
        default='4',
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
                for part in props.ladder_mat:
                    if part == self:
                        return props
        return None
    def update(self, context):
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context)
     
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
        if user_mode:
            box = layout.box()
            row = box.row()
            row.prop(self, "type", text="")
            if self.type in ['C_STAIR', 'C_PALIER', 'D_STAIR', 'D_PALIER']:
                row = box.row()
                row.prop(self, "radius")
                row = box.row()
                row.prop(self, "da")
            else:
                row.prop(self, "length")
            row = box.row(align=True)
            row.prop(self, "left_shape", text="")
            row.prop(self, "right_shape", text="")
        else:
            if self.type in ['S_STAIR', 'S_PALIER']:
                box = layout.box()
                row = box.row()
                row.prop(self, "length")
        
class StairProperty(PropertyGroup):
    parts = CollectionProperty(type=StairPartProperty)
    n_parts = IntProperty(
            name="parts",
            min=1,
            max=32,
            default=1, update=update
            )
    step_depth = FloatProperty(
            name="Depth",
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
            name="Height",
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
    total_angle = FloatProperty(
            name="angle",
            min=-50*pi,
            max=50*pi,
            default=2*pi,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    steps_type = EnumProperty(
            name="Steps",
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
            name="Interp z",
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
    left_ladder=BoolProperty(
        name="left",
        update=update,
        default=False
        )
    right_ladder=BoolProperty(
        name="right",
        update=update,
        default=False
        )
    enable_post = BoolProperty(
            name='Posts',
            default=True,
            update=update
            )
    post_spacing = FloatProperty(
            name="spacing",
            min=0.001, max=1000,
            default=0.50, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    post_x = FloatProperty(
            name="width",
            min=0.001, max=1000,
            default=0.03, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    post_y = FloatProperty(
            name="length",
            min=0.001, max=1000,
            default=0.03, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    post_z = FloatProperty(
            name="height",
            min=0.001, max=1000,
            default=1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    post_offset = FloatProperty(
            name="offset",
            min=-100.0, max=100,
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    post_corners = BoolProperty(
        name="only on edges",
        update=update,
        default=False
        )
    post_mat = EnumProperty(
        name="Post",
        items=materials_enum,
        default='4',
        update=update
        )
    enable_panel = BoolProperty(
            name='Panels',
            default=False,
            update=update
            )
    panel_alt= FloatProperty(
            name="altitude",
            min=-100, max=1000,
            default=0.25, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    panel_x= FloatProperty(
            name="width",
            min=0.001, max=1000,
            default=0.01, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    panel_z= FloatProperty(
            name="height",
            min=0.001, max=1000,
            default=0.6, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    panel_dist= FloatProperty(
            name="space",
            min=0.001, max=1000,
            default=0.05, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    panel_mat = EnumProperty(
        name="Panels",
        items=materials_enum,
        default='5',
        update=update
        )
    ladder_n = IntProperty(
            name="number",
            default=1,
            min=0,
            max=31,
            update=update
            )
    ladder_x = FloatVectorProperty(
            name="width",
            default=[
                0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
                0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
                0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
                0.05,0.05,0.05,0.05,0.05,0.05,0.05
            ],
            size=31,
            min=0.001, max=1000,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            ) 
    ladder_z = FloatVectorProperty(
            name="height",
            default=[
                0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
                0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
                0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
                0.05,0.05,0.05,0.05,0.05,0.05,0.05
            ],
            size=31,
            min=0.001, max=1000,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            ) 
    ladder_offset = FloatVectorProperty(
            name="offset",
            default=[
                0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0,0,
                0,0,0,0,0,0,0
            ],
            size=31,
            min=-100, max=100,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            ) 
    ladder_alt = FloatVectorProperty(
            name="altitude",
            default=[
                1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                1.0,1.0,1.0,1.0,1.0,1.0,1.0
            ],
            size=31,
            min=-100, max=100,
            precision=2, step=1,
            unit='LENGTH',
            update=update
            ) 
    ladder_mat = CollectionProperty(type=StairMaterialProperty)
    auto_update=BoolProperty(default=True)
    idmat_bottom = EnumProperty(
        name="Bottom",
        items=materials_enum,
        default='1',
        update=update
        )
    idmat_face_bottom = EnumProperty(
        name="Front bot",
        items=materials_enum,
        default='0',
        update=update
        )
    idmat_face_top = EnumProperty(
        name="Front up",
        items=materials_enum,
        default='0',
        update=update
        )
    idmat_top = EnumProperty(
        name="Top",
        items=materials_enum,
        default='0',
        update=update
        )
    idmat_side = EnumProperty(
        name="Side",
        items=materials_enum,
        default='1',
        update=update
        )
    
    # UI layout related
    parts_expand = BoolProperty(
            default = False
            )
    ladder_expand = BoolProperty(
            default = False
            )
    idmats_expand = BoolProperty(
            default = False
            )
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
    
        # remove ladders materials
        for i in range(len(self.ladder_mat), self.ladder_n, -1):
            self.ladder_mat.remove(i-1)
        # add ladders
        for i in range(len(self.ladder_mat), self.ladder_n):
            self.ladder_mat.add()
        
        
        # remove parts
        for i in range(len(self.parts), self.n_parts, -1):
            self.parts.remove(i-1)
        # add parts
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
        matids = []
        
        #top, int(self.idmat_face_top), int(self.idmat_face_bottom), int(self.idmat_side), int(self.idmat_bottom) 
        id_materials = [int(self.idmat_top), int(self.idmat_face_top), int(self.idmat_face_bottom), int(self.idmat_side), int(self.idmat_bottom)]
        
        # depth at bottom
        bottom_z = self.bottom_z
        if self.steps_type == 'OPEN':
            # depth at front
            bottom_z = self.step_z
        
        width_left  = 0.5*self.width + self.x_offset
        width_right = 0.5*self.width - self.x_offset   
        
        g = StairGenerator()
        if self.presets == 'STAIR_USER':
            for part in self.parts:
                g.add_part(part.type, self.steps_type, self.z_mode, self.step_z, bottom_z, center, max(width_left+0.01, 
                        width_right+0.01, part.radius), part.da, width_left, width_right, part.length, part.right_shape, part.left_shape)
        elif self.presets == 'STAIR_O':
            n_parts = max(1, int(round(abs(self.total_angle)/pi,0)))
            if self.total_angle > 0:
                dir = 1
            else:
                dir = -1
            last_da = self.total_angle - dir * (n_parts-1) * pi
            if dir * last_da > pi:
                n_parts += 1
                last_da -= dir * pi
            abs_last = dir * last_da
            
            for part in range(n_parts-1):
                g.add_part('D_STAIR', self.steps_type, self.z_mode, self.step_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), dir*pi,  width_left, width_right, 1.0, self.right_shape, self.left_shape)
            print("dir%s, last_da:%s abs_last:%s n_parts:%s" % (dir, last_da, abs_last, n_parts) )
            if round(abs_last, 2) > 0:
                if abs_last > pi/2:
                    g.add_part('C_STAIR', self.steps_type, self.z_mode, self.step_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), dir*pi/2,  width_left, width_right, 1.0, self.right_shape, self.left_shape)
                    g.add_part('C_STAIR', self.steps_type, self.z_mode, self.step_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), last_da - dir*pi/2,  width_left, width_right, 1.0, self.right_shape, self.left_shape)
                else:
                    g.add_part('C_STAIR', self.steps_type, self.z_mode, self.step_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), last_da,  width_left, width_right, 1.0, self.right_shape, self.left_shape)
        else:
            for part in self.parts:
                g.add_part(part.type, self.steps_type, self.z_mode, self.step_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), self.da,  width_left, width_right, part.length, self.right_shape, self.left_shape)
        
        # Stair basis
        g.set_matids(id_materials)
        g.make_stair(self.height, self.step_depth, verts, faces, matids, t_offset=self.t_offset)
        
        # Ladder
        if self.left_ladder or self.right_ladder:
        
            offset_x = 0.5 * self.width - self.post_offset
            post_spacing = self.post_spacing
            if self.post_corners:
                post_spacing = 10000
            if self.enable_post or self.enable_panel:
                g.make_post(self.height, self.step_depth, 0.5*self.post_x, 0.5*self.post_y, post_spacing, self.post_corners, 
                            offset_x, self.post_z, self.left_ladder, self.right_ladder, int(self.post_mat), self.enable_post,
                            self.enable_panel, self.panel_alt-self.post_z, 0.5*self.panel_x, self.panel_z, self.panel_dist, int(self.panel_mat), verts, faces, matids)
            if self.right_ladder:
                for i in range(self.ladder_n):
                    id_materials = [int(self.ladder_mat[i].index) for j in range(5)]
                    g.set_matids(id_materials)
                    g.make_part(self.height, self.step_depth, self.ladder_x[i], self.ladder_z[i], 
                            offset_x + self.ladder_offset[i], self.ladder_alt[i], 'LINEAR_TOP', 'CLOSED', verts, faces, matids)
            if self.left_ladder:
                for i in range(self.ladder_n):
                    id_materials = [int(self.ladder_mat[i].index) for j in range(5)]
                    g.set_matids(id_materials)
                    g.make_part(self.height, self.step_depth, self.ladder_x[i], self.ladder_z[i], 
                            -offset_x - self.ladder_offset[i], self.ladder_alt[i], 'LINEAR_TOP', 'CLOSED', verts, faces, matids)
        
        bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=None, weld=True, clean=True)
        
        
        #bpy.ops.mesh.select_linked()
        #bpy.ops.mesh.faces_shade_smooth()
        
        
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
        box.prop(prop, 'width')
        box.prop(prop, 'height')
        box.prop(prop, 'bottom_z')
        box.prop(prop, 'x_offset')
        box.prop(prop, 'z_mode')
        box = layout.box()
        box.prop(prop, 'steps_type')
        box.prop(prop, 'step_depth')
        box.prop(prop, 'step_z')
        box.prop(prop, 't_offset')
        box = layout.box()
        row = box.row()
        if prop.parts_expand:
            row.prop(prop, 'parts_expand',icon="TRIA_DOWN", icon_only=True, text="Parts", emboss=False)
            if prop.presets == 'STAIR_USER': 
                box.prop(prop, 'n_parts')
            if prop.presets != 'STAIR_USER':
                row = box.row(align=True)
                row.prop(prop, "left_shape", text="")
                row.prop(prop, "right_shape", text="")
                row = box.row()
                row.prop(prop, "radius")
                row = box.row()
                if prop.presets == 'STAIR_O': 
                    row.prop(prop, 'total_angle')
                else:
                    row.prop(prop, 'da')                
            if prop.presets != 'STAIR_O': 
                for i, part in enumerate(prop.parts):
                    part.draw(layout, context, i, prop.presets == 'STAIR_USER')
        else:
            row.prop(prop, 'parts_expand',icon="TRIA_RIGHT", icon_only=True, text="Parts", emboss=False)
        box = layout.box()
        row = box.row()
        if prop.ladder_expand:
            row.prop(prop, 'ladder_expand',icon="TRIA_DOWN", icon_only=True, text="Ladder", emboss=False)
            row = box.row(align=True)
            row.prop(prop, 'left_ladder')
            row.prop(prop, 'right_ladder')            
            if prop.right_ladder or prop.left_ladder:
                box.prop(prop, 'ladder_n')
                for i in range(prop.ladder_n):
                    box = layout.box()
                    box.label(text="Ladder" + str(i+1))
                    box.prop(prop, 'ladder_x', index=i)
                    box.prop(prop, 'ladder_z', index=i)
                    box.prop(prop, 'ladder_alt', index=i)
                    box.prop(prop, 'ladder_offset', index=i)
                    box.prop(prop.ladder_mat[i], 'index', text="")

                box = layout.box()
                box.prop(prop, 'enable_post')
                if prop.enable_post:
                    box.prop(prop, 'post_spacing')
                    box.prop(prop, 'post_x')
                    box.prop(prop, 'post_y')
                    box.prop(prop, 'post_z')
                    box.prop(prop, 'post_offset')
                    box.prop(prop, 'post_corners')
                box = layout.box()
                box.prop(prop, 'enable_panel')
                if prop.enable_panel:
                    box.prop(prop, 'panel_dist')
                    box.prop(prop, 'panel_x')
                    box.prop(prop, 'panel_z')
                    box.prop(prop, 'panel_alt')
        else:
            row.prop(prop, 'ladder_expand',icon="TRIA_RIGHT", icon_only=True, text="Ladder", emboss=False)
        box = layout.box()
        row = box.row()
        if prop.idmats_expand:
            row.prop(prop, 'idmats_expand',icon="TRIA_DOWN", icon_only=True, text="Materials", emboss=False)
            box.prop(prop, 'idmat_top')
            box.prop(prop, 'idmat_face_top')
            box.prop(prop, 'idmat_face_bottom')
            box.prop(prop, 'idmat_bottom')
            box.prop(prop, 'idmat_side')
            box.prop(prop, 'panel_mat')
            box.prop(prop, 'post_mat')
        else:
            row.prop(prop, 'idmats_expand',icon="TRIA_RIGHT", icon_only=True, text="Materials", emboss=False)
        
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
        MaterialUtils.add_stair_materials(o)
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

bpy.utils.register_class(StairMaterialProperty)
bpy.utils.register_class(StairPartProperty)
bpy.utils.register_class(StairProperty)
Mesh.StairProperty = CollectionProperty(type=StairProperty)
bpy.utils.register_class(ARCHIPACK_PT_stair)
bpy.utils.register_class(ARCHIPACK_OT_stair)



