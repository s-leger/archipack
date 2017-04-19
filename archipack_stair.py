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
from .panel import Panel as Lofter
from mathutils import Vector, Matrix
from math import sin, cos, tan, pi, atan2, sqrt, floor, acos
from .archipack_manipulator import ManipulatorStack

class Project():
    def proj_xy(self, t, next=None):
        """
            length of projection along crossing line/circle
            deformation unit vector for profil in xy axis at line/line intersection
            so f(x) = position of point in xy plane
        """
        if next is None:
            return self.normal(t).v.normalized(), 1
        v0 = self.normal(1).v.normalized()  #self.tangeant_unit_vector(1) * self.length
        v1 = next.normal(0).v.normalized()  #next.tangeant_unit_vector(0) * next.length
        direction = v0+v1
        adj = (v0*self.length)*(v1*next.length)
        hyp = (self.length * next.length)
        c = min(1, max(-1, adj/hyp))
        size = 1/cos(0.5*acos(c)) 
        return direction.normalized(), size
    
    def proj_z(self, t, dz0, next=None, dz1=0):
        """
            length of projection along crossing line/circle
            deformation unit vector for profil in z axis at line/line intersection
            so f(y) = position of point in yz plane
        """
        return Vector((0 ,1)), 1
        if next is None:
            dz = dz0/self.length
        else:
            dz = (dz1+dz0)/(self.length+next.length)
        return Vector((0 ,1)), sqrt(1+dz*dz)
        #1/sqrt(1+(dz0/self.length)*(dz0/self.length))
        if next is None:
            return Vector((-dz0, self.length)).normalized(), 1
        v0 = Vector((self.length, dz0))
        v1 = Vector((next.length, dz1))
        direction = Vector((-dz0, self.length)).normalized()+Vector((-dz1, next.length)).normalized()
        adj = v0*v1
        hyp = (v0.length * v1.length)
        c = min(1, max(-1, adj/hyp))
        size = -cos(pi-0.5*acos(c))
        return direction.normalized(), size
    

class Line(Project):
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
    def tangeant_unit_vector(self, t):
        return self.v.normalized()
    def point_sur_segment(self, pt):
        """ _point_sur_segment 
            point: Vector 2d
            t: param t de l'intersection sur le segment courant
            d: distance laterale perpendiculaire positif a droite
        """
        dp = pt-self.p
        dl = self.length
        d = (self.v.x*dp.y-self.v.y*dp.x)/dl
        t = (self.v*dp)/(dl*dl)
        return t > 0 and t < 1, d, t 

class Circle(Project):
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
            res, d, t = line.point_sur_segment(self.c)
            return False, line.lerp(t), t
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
    def tangeant_unit_vector(self, t):
        a = self.a0+t*self.da
        ca = cos(a)
        sa = sin(a)
        v = Vector((sa, -ca)) 
        if self.da > 0:
            v = -v
        return v
    def straight(self, length):
        t = self.tangeant(1,length)
        return t.p, t.v
    def point_sur_segment(self, pt):
        dp = pt-self.c
        d = dp.length-self.r
        a = atan2(dp.y, dp.x)
        t = (a - self.a0)/self.da
        return t > 0 and t < 1, d, t
        
class Stair():
    def __init__(self, left_offset, right_offset, steps_type, nose_type, z_mode, nose_z, bottom_z):
        self.steps_type = steps_type
        self.nose_type = nose_type
        self.l_shape = None 
        self.r_shape = None
        self.next_type = 'NONE'
        self.last_type = 'NONE'
        self.z_mode = z_mode
        # depth of open step
        self.nose_z = nose_z
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
    @property
    def left_length(self):
        return self.get_length("LEFT")
    @property
    def right_length(self):
        return self.get_length("RIGHT")
     
    def step_size(self, step_depth):
        t_step, n_step = self.steps(step_depth) 
        self.n_step = n_step
        self.t_step = t_step
        self.step_depth = step_depth
        return n_step
    
    def p3d_left(self, verts, p2d, i, t, landing=False):
        x, y = p2d
        nose_z = min(self.step_height, self.nose_z)
        zl = self.z0 + t * self.height
        zs = self.z0 + i * self.step_height
        if self.z_mode == 'LINEAR':
            z0 = max(0, zl)
            z1 = z0 - self.nose_z
            verts.extend([(x, y, z0),(x, y, z1)])
        else:
            if "FULL" in self.steps_type:
                z0 = 0
            else:
                z0 = max(0, zl - nose_z - self.bottom_z)
            z3 = zs + max(0, self.step_height - nose_z)
            z4 = zs + self.step_height
            if landing:
                if "FULL" in self.steps_type:
                    z2 = 0
                    z1 = 0
                else:
                    z2 = max(0, min(z3, z3 - self.bottom_z))
                    z1 = z2
            else:    
                z1 = min(z3, max(z0, zl - nose_z))
                z2 = min(z3, max(z1, zl))
            verts.extend([(x, y, z0),
                        (x, y, z1),
                        (x, y, z2),
                        (x, y, z3),
                        (x, y, z4)])
            
    def p3d_right(self, verts, p2d, i, t, landing=False):
        x, y = p2d
        nose_z = min(self.step_height, self.nose_z)
        zl = self.z0 + t * self.height
        zs = self.z0 + i * self.step_height
        if self.z_mode == 'LINEAR':
            z0 = max(0, zl)
            z1 = z0 - self.nose_z
            verts.extend([(x, y, z1),(x, y, z0)])
        else:
            if "FULL" in self.steps_type:
                z0 = 0
            else:
                z0 = max(0, zl - nose_z - self.bottom_z)
            z3 = zs + max(0, self.step_height - nose_z)
            z4 = zs + self.step_height
            if landing:
                if "FULL" in self.steps_type:
                    z2 = 0
                    z1 = 0
                else:
                    z2 = max(0, min(z3, z3 - self.bottom_z))
                    z1 = z2
            else:    
                z1 = min(z3, max(z0, zl - nose_z))
                z2 = min(z3, max(z1, zl))
            verts.extend([(x, y, z4),
                          (x, y, z3),
                          (x, y, z2),
                          (x, y, z1),
                          (x, y, z0)])
            
    def p3d_cstep_left(self, verts, p2d, i, t):
        x, y = p2d
        nose_z = min(self.step_height, self.nose_z)
        zs = self.z0 + i * self.step_height
        z3 = zs + max(0, self.step_height - nose_z)
        z1 = min(z3, zs - nose_z)
        verts.append((x, y, z1))
        verts.append((x, y, z3))
            
    def p3d_cstep_right(self, verts, p2d, i, t):
        x, y = p2d
        nose_z = min(self.step_height, self.nose_z)
        zs = self.z0 + i * self.step_height
        z3 = zs + max(0, self.step_height - nose_z)
        z1 = min(z3, zs - nose_z)
        verts.append((x, y, z3))
        verts.append((x, y, z1))
     
    def straight_stair(self, length):
        self.next_type = 'STAIR'
        p, v = self.straight(length)
        return StraightStair(p, v, self.left_offset, self.right_offset, self.steps_type, self.nose_type, self.z_mode, self.nose_z, self.bottom_z)
    def straight_landing(self, length, last_type='STAIR'):
        self.next_type = 'LANDING'
        p, v = self.straight(length)
        return StraightLanding(p, v, self.left_offset, self.right_offset, self.steps_type, self.nose_type, self.z_mode, self.nose_z, self.bottom_z, last_type=last_type)
    def curved_stair(self, da, radius, left_shape, right_shape, double_limit=pi):
        self.next_type = 'STAIR'
        n = self.normal(1)
        n.v = radius * n.v.normalized()
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedStair(c, radius, a0, da, self.left_offset, self.right_offset, self.steps_type, self.nose_type, self.z_mode, self.nose_z, self.bottom_z, left_shape, right_shape, double_limit=double_limit)
    def curved_landing(self, da, radius, left_shape, right_shape, double_limit=pi, last_type='STAIR'):
        self.next_type = 'LANDING'
        n = self.normal(1)
        n.v = radius * n.v.normalized()
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedLanding(c, radius, a0, da, self.left_offset, self.right_offset, self.steps_type, self.nose_type, self.z_mode, self.nose_z, self.bottom_z, left_shape, right_shape, double_limit=double_limit, last_type=last_type)
    def get_z(self, t, mode):
        if mode == 'LINEAR':
            return self.z0 + t * self.height
        else:
            step = 1+floor(t/self.t_step)
            return self.z0 + step * self.step_height
    def make_profile(self, t, side, profile, verts, faces, matids, next=None, tnext=0):
        z0 = self.get_z(t, 'LINEAR')
        dz1 = 0
        t, part, dz0, shape = self.get_part(t, side)
        if next is not None:
            tnext, next, dz1, shape1 = next.get_part(tnext, side)
        xy, s = part.proj_xy(t, next)
        v_xy = s * xy.to_3d()
        z, s  = part.proj_z(t, dz0, next, dz1)
        v_z = s * Vector((-xy.y * z.x, xy.x * z.x, z.y)) 
        x, y = part.lerp(t)
        verts += [Vector((x, y, z0)) + v.x * v_xy + v.y * v_z for v in profile]
    def project_uv(self, rM, uvs, verts, indexes, up_axis='Z'):
        if up_axis == 'Z':
            uvs.append([(rM * Vector(verts[i])).to_2d() for i in indexes])
        elif up_axis == 'Y':
            uvs.append([(x, z) for x, y, z in [(rM * Vector(verts[i])) for i in indexes]])
        else:
            uvs.append([(y, z) for x, y, z in [(rM * Vector(verts[i])) for i in indexes]])
    def get_proj_matrix(self, part, t, nose_y):
        # a matrix to project verts 
        # into uv space for horizontal parts of this step
        # so uv = (rM * vertex).to_2d()
        no = part.normal(t).offset(nose_y)
        v = no.v.normalized()
        return Matrix([
            [-v.y,v.x, 0, no.p.x],
            [v.x, v.y, 0, no.p.y],
            [0,    0,  1, 0],
            [0,    0,  0, 1]
        ]).inverted()
    def _make_nose(self, i, s, verts, faces, matids, uvs, nose_y):
        
        t = self.t_step * i
        
        if 'Curved' in type(self).__name__: 
            nr = self.r_arc
            nl = self.l_arc
        else:
            nr = self.r_line
            nl = self.l_line
        
        # a matrix to project verts 
        # into uv space for horizontal parts of this step
        # so uv = (rM * vertex).to_2d()
        rM = self.get_proj_matrix(nl, t, nose_y)
        
        if self.z_mode == 'LINEAR':
            return rM
        
        
        f = len(verts)
            
        nr = nr.normal(t).offset(nose_y)
        nl = nl.normal(t).reversed.offset(-nose_y)
        
        t2, part, dz, shape = self.get_part(t, "LEFT")
        res, p, u = part.intersect(nl)
        self.p3d_left(verts, p, s, u)
        
        t2, part, dz, shape = self.get_part(t, "RIGHT")        
        res, p, u = part.intersect(nr)
        self.p3d_right(verts, p, s, u)
        
        start = 3
        end   = 6
        offset= 10
        
        # left, top, right 
        matids.extend([self.idmat_side,
             self.idmat_top,
             self.idmat_side])
             
        faces += [(f+j, f+j+1, f+j+offset+1, f+j+offset) for j in range(start, end)]
        
        u = nose_y
        v = (nr.p-nl.p).length
        w = verts[f+2][2]-verts[f+3][2]
        s = int((end-start)/2)
        uvs += [[(u, verts[f+j][2]), (u, verts[f+j+1][2]), (0,verts[f+j+1][2]), (0,verts[f+j][2])] for j in range(start, start+s)]
        uvs.append([(0,0), (0, v), (u, v), (u, 0)])
        uvs += [[(u, verts[f+j][2]), (u, verts[f+j+1][2]), (0,verts[f+j+1][2]), (0,verts[f+j][2])] for j in range(start+s+1, end)]
        
        if 'STRAIGHT' in self.nose_type or 'OPEN' in self.steps_type:
            # face bottom
            matids.append(self.idmat_bottom)
            faces.append((f+end,f+start,f+offset+start,f+offset+end))
            uvs.append([(u, v), (u, 0), (0,0), (0, v)])
        
        if self.steps_type != 'OPEN':
            if 'STRAIGHT' in self.nose_type:
                # front face bottom straight
                matids.append(self.idmat_face_bottom)
                faces.append((f+12, f+17, f+16, f+13))
                uvs.append([(0,w), (v,w), (v,0), (0,0)])
                
            elif 'OBLIQUE' in self.nose_type:        
                # front face bottom oblique
                matids.append(self.idmat_face_bottom)
                faces.append((f+12, f+17, f+6, f+3))
                
                uvs.append([(0,w), (v,w), (v,0), (0,0)])
                
                matids.append(self.idmat_side)
                faces.append((f+3, f+13, f+12))
                uvs.append([(0,0), (u,0), (u,w)])
                
                matids.append(self.idmat_side)
                faces.append((f+6, f+17, f+16))
                uvs.append([(0,0), (u,w), (u,0)])
                
        # front face top
        w = verts[f+3][2]-verts[f+4][2]
        matids.append(self.idmat_face_top)
        faces.append((f+4, f+3, f+6, f+5))
        uvs.append([(0,0), (0,w), (v,w), (v,0)])
        return rM
    def make_faces(self, f, rM, verts, faces, matids, uvs):
        
        if self.z_mode == 'LINEAR':
            start = 0
            end   = 3
            offset= 4
            matids.extend([self.idmat_side,
                 self.idmat_top,
                 self.idmat_side,
                 self.idmat_bottom])
        elif "OPEN" in self.steps_type:
            # faces dessus-dessous-lateral marches fermees
            start = 3
            end   = 6
            offset= 10
            matids.extend([self.idmat_side,
                 self.idmat_top,
                 self.idmat_side,
                 self.idmat_bottom])
        else:
            # faces dessus-dessous-lateral marches fermees
            start = 0
            end   = 9            
            offset= 10
            matids.extend([self.idmat_side,
                 self.idmat_side,
                 self.idmat_side,
                 self.idmat_side,
                 self.idmat_top,
                 self.idmat_side,
                 self.idmat_side,
                 self.idmat_side,
                 self.idmat_side,
                 self.idmat_bottom])
        
        u_l0 = 0
        u_l1 = self.t_step * self.left_length
        u_r0 = 0
        u_r1 = self.t_step * self.right_length
        
        w = 1
        s = int((end-start)/2)
        uvs += [[(u_l0, verts[f+j][2]), (u_l0, verts[f+j+1][2]), (u_l1, verts[f+j+offset+1][2]), (u_l1, verts[f+j+offset][2])] for j in range(start, start+s)]
        self.project_uv(rM, uvs, verts, [f+start+s, f+start+s+1, f+start+s+offset+1, f+start+s+offset])
        #uvs.append([(u_l0, 0), (u_r0, w), (u_r1, w), (u_l1,0)])
        uvs += [[(u_r0, verts[f+j][2]), (u_r0, verts[f+j+1][2]), (u_r1, verts[f+j+offset+1][2]), (u_r1, verts[f+j+offset][2])] for j in range(start+s+1, end)]
        #uvs.append([(u_r1, w), (u_l1,0), (u_l0, 0), (u_r0, w)])
        self.project_uv(rM, uvs, verts, [f+end, f+start, f+offset+start,f+offset+end])
        faces += [(f+j, f+j+1, f+j+offset+1, f+j+offset) for j in range(start, end)]
        faces.append((f+end,f+start,f+offset+start,f+offset+end))
    
        
    
class StraightStair(Stair, Line):
    def __init__(self, p, v, left_offset, right_offset, steps_type, nose_type, z_mode, nose_z, bottom_z):
        Stair.__init__(self, left_offset, right_offset, steps_type, nose_type, z_mode, nose_z, bottom_z)
        Line.__init__(self, p, v)
        self.l_line = self.offset(-left_offset)
        self.r_line = self.offset(right_offset)
    def make_step(self, i, verts, faces, matids, uvs, nose_y=0):
        
        rM = self._make_nose(i, i, verts, faces, matids, uvs, nose_y)
        
        t0 = self.t_step * i
        
        f = len(verts)
        
        p = self.l_line.lerp(t0)
        self.p3d_left(verts, p, i, t0)
        p = self.r_line.lerp(t0)
        self.p3d_right(verts, p, i, t0)
        
        t1 = t0 + self.t_step
        
        p = self.l_line.lerp(t1)
        self.p3d_left(verts, p, i, t1)
        p = self.r_line.lerp(t1)
        self.p3d_right(verts, p, i, t1)
            
        self.make_faces(f, rM, verts, faces, matids, uvs)
        
        if "OPEN" in self.steps_type:
            faces.append((f+13, f+14, f+15, f+16))
            matids.append(self.idmat_face_top)
            uvs.append([(0, 0), (0, 1), (1, 1), (1,0)])
    def get_length(self, side):
        return self.length
    def get_lerp_vect(self, posts, side, i, t_step, respect_edges):
        t0 = i * t_step
        t, part, dz, shape = self.get_part(t0, side)
        n = part.normal(t)
        z0 = self.get_z(t0, 'STEP')
        z1 = self.get_z(t0, 'LINEAR')
        posts.append((n, dz, z0, z1))
        return t0
        
    def n_posts(self, post_spacing, side, respect_edges):
        return self.steps(post_spacing)
    def get_part(self, t, side):
        if side == 'LEFT':
            part = self.l_line
        else:
            part = self.r_line
        return t, part, self.height, 'LINE'
    
class CurvedStair(Stair, Arc):
    def __init__(self, c, radius, a0, da, left_offset, right_offset, steps_type, nose_type, z_mode, nose_z, bottom_z, left_shape, right_shape, double_limit=pi):
        Stair.__init__(self, left_offset, right_offset, steps_type, nose_type, z_mode, nose_z, bottom_z)
        Arc.__init__(self, c, radius, a0, da)
        self.l_shape = left_shape 
        self.r_shape = right_shape
        self.edges_multiples = round(abs(da),6) > double_limit
        #left arc, tangeant at start and end
        self.l_arc, self.l_t0, self.l_t1, self.l_tc = self.set_offset(-left_offset, left_shape)
        self.r_arc, self.r_t0, self.r_t1, self.r_tc = self.set_offset(right_offset, right_shape)
    def set_offset(self, offset, shape):
        arc = self.offset(offset)
        t0 = arc.tangeant(0,1)
        t1 = arc.tangeant(1,1)
        tc = arc.tangeant(0.5,1)
        if self.edges_multiples: 
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
    
    def get_length(self, side):
        if side == 'RIGHT':
            arc = self.r_arc
            shape = self.r_shape
            t0 = self.r_t0
        else:
            arc = self.l_arc
            shape = self.l_shape
            t0 = self.l_t0
        if shape == 'CIRCLE':
            return arc.length
        else:
            if self.edges_multiples:    
                # two edges
                return t0.length * 4
            else:    
                return t0.length * 2
       
    def _make_step(self, t_step, i, s, verts, landing=False):
        
        tb = t_step * i
        
        f = len(verts)
               
        t, part, dz, shape = self.get_part(tb, "LEFT")
        p = part.lerp(t)
        self.p3d_left(verts, p, s, tb, landing)
        
        t, part, dz, shape = self.get_part(tb, "RIGHT")        
        p = part.lerp(t)
        self.p3d_right(verts, p, s, tb, landing)
        return f
    
    def _make_edge(self, t_step, i, j, f, rM, verts, faces, matids, uvs):
        tb = t_step * i
        # make edges verts after regular ones
        if self.l_shape != 'CIRCLE' or self.r_shape != 'CIRCLE':
            if self.edges_multiples:
                # edge 1
                if tb < 0.25 and tb + t_step > 0.25:
                    f0 = f
                    f = len(verts)
                    if self.l_shape == 'CIRCLE':
                        self.p3d_left(verts, self.l_arc.lerp(0.25), j, 0.25)
                    else:
                        self.p3d_left(verts, self.l_tc.p, j, 0.25)
                    if self.r_shape == 'CIRCLE':
                        self.p3d_right(verts, self.r_arc.lerp(0.25), j, 0.25)
                    else:
                        self.p3d_right(verts,  self.r_tc.p, j, 0.25)
                    self.make_faces(f0, rM, verts, faces, matids, uvs)
                # edge 2
                if tb < 0.75 and tb + t_step > 0.75:     
                    f0 = f
                    f = len(verts)
                    if self.l_shape == 'CIRCLE':
                        self.p3d_left(verts, self.l_arc.lerp(0.75), j, 0.75)
                    else:
                        self.p3d_left(verts, self.l_t1.p, j, 0.75)
                    if self.r_shape == 'CIRCLE':
                        self.p3d_right(verts, self.r_arc.lerp(0.75), j, 0.75)
                    else:
                        self.p3d_right(verts,  self.r_t1.p, j, 0.75)
                    self.make_faces(f0, rM, verts, faces, matids, uvs)
            else:        
                if tb < 0.5 and tb + t_step > 0.5:     
                    f0 = f
                    f = len(verts)
                    # the step goes through the edge
                    if self.l_shape == 'CIRCLE':
                        self.p3d_left(verts, self.l_arc.lerp(0.5), j, 0.5)
                    else:
                        self.p3d_left(verts, self.l_t1.p, j, 0.5)
                    if self.r_shape == 'CIRCLE':
                        self.p3d_right(verts, self.r_arc.lerp(0.5), j, 0.5)
                    else:
                        self.p3d_right(verts,  self.r_t1.p, j, 0.5)
                    self.make_faces(f0,  rM, verts, faces, matids, uvs)
        return f
        
    def make_step(self, i, verts, faces, matids, uvs, nose_y=0):
        
        # open stair with closed face
        
        # step nose
        rM = self._make_nose(i, i, verts, faces, matids, uvs, nose_y)
        f = 0
        if self.l_shape == 'CIRCLE' or self.r_shape == 'CIRCLE':
            # every 6 degree
            n_subs =  max(1, int(abs(self.da)/pi*30/self.n_step))
            t_step = self.t_step/n_subs
            for j in range(n_subs):
                f0 = f
                f = self._make_step(t_step, n_subs*i+j, i, verts)
                if j > 0:
                    self.make_faces(f0, rM, verts, faces, matids, uvs)
                f = self._make_edge(t_step, n_subs*i+j, i, f, rM, verts, faces, matids, uvs)
        else:
            f = self._make_step(self.t_step, i, i, verts)
            f = self._make_edge(self.t_step, i, i, f, rM, verts, faces, matids, uvs)
            
        self._make_step(self.t_step, i+1, i, verts)
        self.make_faces(f, rM, verts, faces, matids, uvs)
        
        if "OPEN" in self.steps_type and self.z_mode != 'LINEAR':
            # back face top
            faces.append((f+13, f+14, f+15, f+16))
            matids.append(self.idmat_face_top)
            uvs.append([(0, 0), (0, 1), (1, 1), (1,0)])
            
    def get_part(self, t, side):
        if side == 'RIGHT':
            arc = self.r_arc
            shape = self.r_shape
            t0, t1, tc = self.r_t0, self.r_t1, self.r_tc
        else:
            arc = self.l_arc
            shape = self.l_shape
            t0, t1, tc = self.l_t0, self.l_t1, self.l_tc
        if shape == 'CIRCLE':
            return t, arc, self.height, shape
        else:
            if self.edges_multiples:    
                # two edges
                if t <= 0.25:
                    return 4 * t, t0, 0.25*self.height, shape
                elif t <= 0.75:
                    return 2 * ( t - 0.25 ), tc, 0.5*self.height, shape
                else:
                    return 4 * (t - 0.75), t1, 0.25*self.height, shape
            else:    
                if t <= 0.5:
                    return 2 * t, t0, 0.5*self.height, shape
                else:
                    return 2 * (t - 0.5), t1, 0.5*self.height, shape 
                    
    def get_lerp_vect(self, posts, side, i, t_step, respect_edges):
        
        t0 = i * t_step
        t2 = t0
        t1 = t0 + t_step
        zs = self.get_z(t0, 'STEP')
        zl = self.get_z(t0, 'LINEAR')
        
        # vect normal
        t, part, dz, shape = self.get_part(t0, side)
        n = part.normal(t)
        posts.append((n, dz, zs, zl))
        
        if shape != 'CIRCLE' and respect_edges:
            if self.edges_multiples:
                if t0 < 0.25 and t1 > 0.25:
                    zs = self.get_z(0.25, 'STEP')
                    zl = self.get_z(0.25, 'LINEAR')
                    t, part, dz, shape = self.get_part(0.25, side)
                    n = part.normal(1)
                    posts.append((n, dz, zs, zl))
                    t2 = 0.25
                if t0 < 0.75 and t1 > 0.75:
                    zs = self.get_z(0.75, 'STEP')
                    zl = self.get_z(0.75, 'LINEAR')
                    t, part, dz, shape = self.get_part(0.75, side)
                    n = part.normal(1)
                    posts.append((n, dz, zs, zl))
                    t2 = 0.75
            elif t0 < 0.5 and t1 > 0.5:
                    zs = self.get_z(0.5, 'STEP')
                    zl = self.get_z(0.5, 'LINEAR')
                    t, part, dz, shape = self.get_part(0.5, side)
                    n = part.normal(1)
                    posts.append((n, dz, zs, zl))
                    t2 = 0.5
        return t2
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
        
class StraightLanding(StraightStair):
    def __init__(self, p, v, left_offset, right_offset, steps_type, nose_type, z_mode, nose_z, bottom_z, last_type='STAIR'):
        StraightStair.__init__(self, p, v, left_offset, right_offset, steps_type, nose_type, z_mode, nose_z, bottom_z)
        self.last_type=last_type
    @property
    def height(self):
        return 0
    @property
    def top_offset(self):
        return self.t_step / self.v.length
    @property
    def top(self):
        if self.next_type == 'LANDING':
            return self.z0
        else:
            return self.z0 + self.step_height
    def step_size(self, step_depth):
        self.n_step = 1
        self.t_step = 1
        self.step_depth = step_depth
        if self.last_type == 'LANDING':
            return 0
        else:
            return 1
    def make_step(self, i, verts, faces, matids, uvs, nose_y=0):
    
        if i == 0 and self.last_type != 'LANDING':
            rM = self._make_nose(i, 0, verts, faces, matids, uvs, nose_y)
        else:
            rM = self.get_proj_matrix(self.l_line, self.t_step * i, nose_y)
        
            
        f = len(verts)
        j = 0
        t0 = self.t_step * i
       
        p = self.l_line.lerp(t0)
        self.p3d_left(verts, p, j, t0)
        
        p = self.r_line.lerp(t0)
        self.p3d_right(verts, p, j, t0)
        
        t1 = t0 + self.t_step
        p = self.l_line.lerp(t1)
        self.p3d_left(verts, p, j, t1, self.next_type != 'LANDING')
        
        p = self.r_line.lerp(t1)
        self.p3d_right(verts, p, j, t1, self.next_type != 'LANDING')
        
        self.make_faces(f, rM, verts, faces, matids, uvs)
        
        if "OPEN" in self.steps_type and self.next_type != 'LANDING':
            faces.append((f+13, f+14, f+15, f+16))
            matids.append(self.idmat_face_top)
            uvs.append([(0, 0), (0, 1), (1, 1), (1,0)])
        
    def straight_landing(self, length):
        return Stair.straight_landing(self, length, last_type='LANDING')
    def curved_landing(self, da, radius, left_shape, right_shape, double_limit=pi):
        return Stair.curved_landing(self, da, radius, left_shape, right_shape, double_limit=double_limit, last_type='LANDING')
    def get_z(self, t, mode):
        if mode == 'STEP':
            return self.z0 + self.step_height
        else:
            return self.z0

class CurvedLanding(CurvedStair):
    def __init__(self, c, radius, a0, da, left_offset, right_offset, steps_type, nose_type, z_mode, nose_z, bottom_z, left_shape, right_shape, double_limit=pi, last_type='STAIR'):
        CurvedStair.__init__(self, c, radius, a0, da, left_offset, right_offset, steps_type, nose_type, z_mode, nose_z, bottom_z, left_shape, right_shape, double_limit=double_limit)
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
        if self.next_type == 'LANDING':
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
        if self.last_type == 'LANDING':
            return 0
        else:
            return 1
        
    def make_step(self, i, verts, faces, matids, uvs, nose_y=0): 
        """
        if i == self.n_step and self.steps_type in ['CLOSED','FULL']:
            if self.next_type == 'STAIR':
                p = self.l_arc.lerp(1)
                self.p3d_end(verts, p, j, 1)
                p = self.r_arc.lerp(1)
                self.p3d_end(verts, p, j, 1)
                # close last part bottom
                #faces.append((f, f+4, f+5, f+6, f+2,f+1))
                return
            elif self.next_type == 'NONE':
                #faces.append((f, f+4, f+5, f+6, f+7, f+3, f+2, f+1))
        """
        if i == 0 and self.last_type != 'LANDING':
            rM = self._make_nose(i, 0, verts, faces, matids, uvs, nose_y)
        else:
            rM = self.get_proj_matrix(self.l_arc, self.t_step * i, nose_y)
            
        f = len(verts)
        
        if self.l_shape == 'CIRCLE' or self.r_shape == 'CIRCLE':
            n_subs = max(1, int(abs(self.da/pi*30/self.n_step)))
            t_step = self.t_step/n_subs
            for j in range(n_subs):
                f0 = f
                f = self._make_step(t_step, n_subs*i+j, 0, verts)
                if j > 0:
                    self.make_faces(f0, rM,  verts, faces, matids, uvs)
                f = self._make_edge(t_step, n_subs*i+j, 0, f, rM, verts, faces, matids, uvs)
        else:
            f = self._make_step(self.t_step, i, 0, verts)
            f = self._make_edge(self.t_step, i, 0, f, rM, verts, faces, matids, uvs)
          
        self._make_step(self.t_step, i+1, 0, verts, i == self.n_step-1 and self.next_type != 'LANDING')
        self.make_faces(f, rM,  verts, faces, matids, uvs)
        
        if "OPEN" in self.steps_type and self.next_type != 'LANDING':
            faces.append((f+13, f+14, f+15, f+16))
            matids.append(self.idmat_face_top)
            uvs.append([(0, 0), (0, 1), (1, 1), (1,0)])
        
    def straight_landing(self, length):
        return Stair.straight_landing(self, length, last_type='LANDING')
    def curved_landing(self, da, radius, left_shape, right_shape, double_limit=pi):
        return Stair.curved_landing(self, da, radius, left_shape, right_shape, double_limit=double_limit,  last_type='LANDING')
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
    def add_part(self, type, steps_type, nose_type, z_mode, nose_z, bottom_z, center, radius, da, width_left, width_right, length, left_shape, right_shape):
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
                s = StraightStair(p, v, width_left, width_right, steps_type, nose_type, z_mode, nose_z, bottom_z)
            elif type == 'C_STAIR':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedStair(c, radius, 0, da, width_left, width_right, steps_type, nose_type, z_mode, nose_z, bottom_z, left_shape, right_shape)
            elif type == 'D_STAIR':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedStair(c, radius, 0, da, width_left, width_right, steps_type, nose_type, z_mode, nose_z, bottom_z, left_shape, right_shape, double_limit=0)
            elif type == 'S_LANDING':
                p = Vector((0,0))
                v = Vector((0,length))
                s = StraightLanding(p, v, width_left, width_right, steps_type, nose_type, z_mode, nose_z, bottom_z)
            elif type == 'C_LANDING':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedLanding(c, radius, 0, da, width_left, width_right, steps_type, nose_type, z_mode, nose_z, bottom_z, left_shape, right_shape)
            elif type == 'D_LANDING':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedLanding(c, radius, 0, da, width_left, width_right, steps_type, nose_type, z_mode, nose_z, bottom_z, left_shape, right_shape, double_limit=0)
        else:
            if type == 'S_STAIR':
                s = s.straight_stair(length)
            elif type == 'C_STAIR':
                s = s.curved_stair(da, radius, left_shape, right_shape)
            elif type == 'D_STAIR':
                s = s.curved_stair(da, radius, left_shape, right_shape, double_limit=0)
            elif type == 'S_LANDING':
                s = s.straight_landing(length)
            elif type == 'C_LANDING':
                s = s.curved_landing(da, radius, left_shape, right_shape)
            elif type == 'D_LANDING':
                s = s.curved_landing(da, radius, left_shape, right_shape, double_limit=0)
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
    def make_stair(self, height, step_depth, verts, faces, matids, uvs, nose_y=0):
        n_steps = self.n_steps(step_depth)
        self.set_height(height / n_steps)
        #if self.steps_type in ['CLOSED', 'FULL'] and len(self.stairs) > 0:
            #faces.append((0,1,5,4))
            #matids.append(self.stairs[-1].idmat_face_bottom)
        for s, stair in enumerate(self.stairs):
            for i in range(stair.n_step):
                stair.make_step(i, verts, faces, matids, uvs, nose_y=nose_y)
                if s < len(self.stairs)-1 and self.steps_type != 'OPEN' and 'Landing' in type(stair).__name__ and stair.next_type != "LANDING":
                    f = len(verts)-10
                    faces.append((f,f+1,f+8,f+9))
                    matids.append(self.stairs[-1].idmat_bottom)
                    u = verts[f+1][2]-verts[f][2]
                    v = (Vector(verts[f])-Vector(verts[f+9])).length
                    uvs.append([(0,0),(0,u),(v,u),(v,0)])
                    
        if self.steps_type != 'OPEN' and len(self.stairs) > 0:
            f = len(verts)-10
            faces.append((f,f+1,f+2,f+3,f+4,f+5,f+6,f+7,f+8,f+9))
            matids.append(self.stairs[-1].idmat_bottom)
            uvs.append([(0,0),(.1,0),(.2,0),(.3,0),(.4,0),(.4,1),(.3,1),(.2,1),(.1,1),(0,1)])
    def get_post(self, post, post_x, post_y, post_z, post_alt, id_mat, verts, faces, matids, uvs):
        n, dz, z0, z1 = post
        z0 += post_alt
        z1 += post_z+post_alt
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
        x = [(0,0), (0,post_z), (post_x,post_z), (post_x,0)]
        y = [(0,0), (0,post_z), (post_y,post_z), (post_y, 0)]
        z = [(0,0), (post_x,0),(post_x,post_y), (0,post_y)]
        uvs.extend([x, y, x, y, z, z]) 
       
    def get_panel(self, n0, n1, dist, z0, z1, panel_x, panel_z, id_mat, verts, faces, matids, uvs):
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
        uv = [(0,0), (0,1), (1,1), (1, 0)]
        uvs.extend([uv, uv, uv, uv, uv, uv]) 
       
    def make_post(self, height, step_depth, x, y, z, altitude, post_spacing, panel_dist, respect_edges, move_x, x_offset, enable_r, enable_l, mat, enable_post, 
        verts, faces, matids, uvs):
        n_steps = self.n_steps(step_depth)
        self.set_height(height / n_steps)
        l_posts = []
        r_posts = []
        n_stairs = len(self.stairs)-1
        for s, stair in enumerate(self.stairs):
            if type(stair).__name__ in ['CurvedStair', 'CurvedLanding']:
                stair.l_arc, stair.l_t0, stair.l_t1, stair.l_tc = stair.set_offset(move_x - x_offset, stair.l_shape)
                stair.r_arc, stair.r_t0, stair.r_t1, stair.r_tc = stair.set_offset(move_x + x_offset, stair.r_shape)
            else:
                stair.l_line = stair.offset(move_x - x_offset)
                stair.r_line = stair.offset(move_x + x_offset)
            if enable_l:
                t_step, n_step = stair.n_posts(post_spacing, 'LEFT', respect_edges)
                if s == n_stairs:
                    n_step += 1
                for i in range(n_step):
                    stair.get_lerp_vect(l_posts, 'LEFT', i, t_step, respect_edges)
                    if s == n_stairs and i == n_step-1:
                        n, dz, z0, z1 = l_posts[-1]
                        l_posts[-1] = (n, dz, z0-stair.step_height, z1)
            if enable_r:  
                t_step, n_step = stair.n_posts(post_spacing, 'RIGHT', respect_edges)
                if s == n_stairs:
                    n_step += 1
                for i in range(n_step):
                    stair.get_lerp_vect(r_posts, 'RIGHT', i, t_step, respect_edges)
                    if s == n_stairs and i == n_step-1:
                        n, dz, z0, z1 = r_posts[-1]
                        r_posts[-1] = (n, dz, z0-stair.step_height, z1)
        if enable_post: 
            if enable_l:
                for i, post in enumerate(l_posts):
                    self.get_post(post, x, y, z, altitude, mat, verts, faces, matids, uvs)
            if enable_r:    
                for i, post in enumerate(r_posts):
                    self.get_post(post, x, y, z, altitude, mat, verts, faces, matids, uvs)
        else:
            """
            # panels between posts
            # this is near production ready, still need to be exposed
            """ 
            if enable_l: 
                for i, p0 in enumerate(l_posts):
                    if i < len(l_posts)-1:
                        n0, dz, zx, z0 = p0
                        n1, dz, zx, z1 = l_posts[i+1]
                        self.get_panel( n0, n1, panel_dist, z0+altitude, z1+altitude, x, z, mat, verts, faces, matids, uvs)
            if enable_r: 
                for i, p0 in enumerate(r_posts):
                    if i < len(r_posts)-1:
                        n0, dz, zx, z0 = p0
                        n1, dz, zx, z1 = r_posts[i+1]
                        self.get_panel( n0, n1, panel_dist, z0+altitude, z1+altitude, x, z, mat, verts, faces, matids, uvs)
    def make_part(self, height, step_depth, part_x, part_z, x_offset, z_offset, z_mode, steps_type, verts, faces, matids, uvs):
        # NOTE: may allow FULL_HEIGHT
        params = [(stair.z_mode, stair.l_shape, stair.r_shape, stair.bottom_z, stair.steps_type) for stair in self.stairs]
        for stair in self.stairs:
            if x_offset > 0:
                stair.l_shape = stair.r_shape
            else:
                stair.r_shape = stair.l_shape
            stair.steps_type = steps_type
            stair.z_mode = "LINEAR" 
            stair.bottom_z = part_z
            if 'Curved' in type(stair).__name__:
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
            for i in range(stair.n_step):
                stair.make_step(i, verts, faces, matids, uvs, nose_y=0)
            stair.n_step /= 2
            stair.t_step *= 2
            stair.step_height *= 2
            stair.z_mode = params[j][0]
            stair.l_shape = params[j][1]
            stair.r_shape = params[j][2]
            stair.bottom_z = params[j][3]
            stair.steps_type = params[j][4]
            stair.z0 -= z_offset+part_z

    def make_profile(self, profile, idmat, side, slice, height, step_depth, x_offset, z_offset, extend, verts, faces, matids, uvs):
        for stair in self.stairs:
            if 'Curved' in type(stair).__name__:
                stair.l_arc, stair.l_t0, stair.l_t1, stair.l_tc = stair.set_offset(-x_offset, stair.l_shape)
                stair.r_arc, stair.r_t0, stair.r_t1, stair.r_tc = stair.set_offset(x_offset, stair.r_shape)
            else:
                stair.l_line = stair.offset(-x_offset)
                stair.r_line = stair.offset(x_offset)
        
        
        n_steps = self.n_steps(step_depth)  
        self.set_height(height/n_steps)
        
        n_stairs = len(self.stairs)-1
        
        if n_stairs < 0:
            return
            
        sections = []
        sections.append([])
        
        # first step
        if extend != 0:
            t = -extend/self.stairs[0].length
            self.stairs[0].get_lerp_vect(sections[-1], side, 1, t, True)
            
        for s, stair in enumerate(self.stairs):
            n_step = 1
            is_circle = False
            
            if 'Curved' in type(stair).__name__:
                if side == "LEFT":
                    part = stair.l_arc
                    is_circle = stair.l_shape == "CIRCLE"
                else:
                    part = stair.r_arc
                    is_circle = stair.r_shape == "CIRCLE"
            else:
                if side == "LEFT":
                    part = stair.l_line
                else:
                    part = stair.r_line
                    
            if is_circle:
                n_step = 3*stair.n_step
            
            t_step = 1/n_step
               
            last_t = 1.0
            
            # last section 1 step before stair
            if 'Landing' in type(stair).__name__ and stair.next_type == 'STAIR':
                if not slice:
                    line = stair.normal(1).offset(self.stairs[s+1].step_depth)
                    res, p, t_part = part.intersect(line)
                    res, d, last_t = part.point_sur_segment(p)
                        
            if s == n_stairs:
                n_step += 1
            cur_t = 0 
            for i in range(n_step):
                cur_t = stair.get_lerp_vect(sections[-1], side, i, t_step, True)
                if cur_t > last_t:
                    del sections[-1][-1]
                    break
        
            # last section 1 step before next stair start
            if 'Landing' in type(stair).__name__ and stair.next_type == 'STAIR':
                stair.get_lerp_vect(sections[-1], side, 1, last_t, True)
                if slice:
                    sections.append([])
                    if extend > 0:
                        t = -extend/self.stairs[s+1].length
                        self.stairs[s+1].get_lerp_vect(sections[-1], side, 1, t, True)
        
        t = 1+extend/self.stairs[-1].length
        self.stairs[-1].get_lerp_vect(sections[-1], side, 1, t, True)
            
        for cur_sect in sections:
            user_path_verts = len(cur_sect)   
            f = len(verts)
            if user_path_verts > 0:
                uv_v = 0
                user_path_uv_v = [] 
                n, dz, z0, z1 = cur_sect[-1]
                cur_sect[-1] = (n, dz, z0-stair.step_height, z1)
                n_sections = user_path_verts-1
                n, dz, zs, zl = cur_sect[0]
                p0 = n.p
                v0 = n.v.normalized()
                for s, section in enumerate(cur_sect):
                    n, dz, zs, zl = section
                    p1 = n.p
                    if s < n_sections:
                        v1 = cur_sect[s+1][0].v.normalized()
                    dir = (v0+v1).normalized()
                    scale = 1/cos(0.5*acos(min(1, max(-1, v0*v1)))) 
                    for p in profile:
                        x, y = n.p + scale * p.x * dir
                        z = zl + p.y + z_offset
                        verts.append((x,y,z))
                    if s > 0:
                        user_path_uv_v.append((p1 - p0).length)
                    p0 = p1
                    v0 = v1
                    
                # build faces using Panel    
                lofter = Lofter(
                    # closed_shape, index, x, y, idmat
                    True,
                    [i for i in range(len(profile))],
                    [p.x for p in profile],
                    [p.y for p in profile],
                    [idmat for i in range(len(profile))],
                    closed_path=False,
                    user_path_uv_v=user_path_uv_v,
                    user_path_verts=user_path_verts
                    )
                faces += lofter.faces(16, offset=f, path_type='USER_DEFINED')
                matids += lofter.mat(16, idmat, idmat, path_type='USER_DEFINED')
                v = Vector((0,0))
                uvs += lofter.uv(16, v, v, v, v, 0, v, 0, 0, path_type='USER_DEFINED')
                
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
                ('S_LANDING', 'Straight landing','',3),
                ('C_LANDING', 'Curved landing','',4),
                ('D_LANDING', 'Dual Curved landing','',5)
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
            if self.type in ['C_STAIR', 'C_LANDING', 'D_STAIR', 'D_LANDING']:
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
            if self.type in ['S_STAIR', 'S_LANDING']:
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
            name="Going",
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
    nose_y = FloatProperty(
        name="Depth",
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
    nose_z = FloatProperty(
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
    nose_type = EnumProperty(
            name="Nosing",
            items=(
                ('STRAIGHT', 'Straight','',0),
                ('OBLIQUE', 'Oblique','',1),
                ),
            default='OBLIQUE',
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
    left_post = BoolProperty(
            name='left',
            default=True,
            update=update
            )
    right_post = BoolProperty(
            name='right',
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
    post_alt= FloatProperty(
            name="altitude",
            min=-100, max=1000,
            default=0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    post_offset_x = FloatProperty(
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
    idmat_post = EnumProperty(
        name="Post",
        items=materials_enum,
        default='4',
        update=update
        )
    left_panel = BoolProperty(
            name='left',
            default=True,
            update=update
            )
    right_panel = BoolProperty(
            name='right',
            default=True,
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
    idmat_panel = EnumProperty(
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
    
    left_handrail=BoolProperty(
            name="left",
            update=update,
            default=False
            )
    right_handrail=BoolProperty(
            name="right",
            update=update,
            default=False
            )
    handrail_offset = FloatProperty(
            name="offset",
            min=-100.0, max=100,
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    handrail_alt = FloatProperty(
            name="altitude",
            min=-100, max=1000,
            default=1.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    handrail_extend = FloatProperty(
            name="extend",
            min=-100, max=1000,
            default=0.1, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    handrail_slice_left=BoolProperty(
            name='slice',
            default=True,
            update=update
            )
    handrail_slice_right=BoolProperty(
            name='slice',
            default=True,
            update=update
            )
    
    left_string=BoolProperty(
            name="left",
            update=update,
            default=False
            )
    right_string=BoolProperty(
            name="right",
            update=update,
            default=False
            )
    string_x = FloatProperty(
            name="width",
            min=-100.0, max=100,
            default=0.02, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    string_z = FloatProperty(
            name="height",
            min=-100.0, max=100,
            default=0.3, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    string_offset = FloatProperty(
            name="offset",
            min=-100.0, max=100,
            default=0.0, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    string_alt = FloatProperty(
            name="altitude",
            min=-100, max=1000,
            default=-0.04, precision=2, step=1,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            ) 
    
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
    idmat_handrail = EnumProperty(
        name="Handrail",
        items=materials_enum,
        default='4',
        update=update
        )
    idmat_string = EnumProperty(
        name="String",
        items=materials_enum,
        default='3',
        update=update
        )
   
    # UI layout related
    parts_expand = BoolProperty(
            default = False
            )
    balluster_expand = BoolProperty(
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
        uvs = []
        #top, int(self.idmat_face_top), int(self.idmat_face_bottom), int(self.idmat_side), int(self.idmat_bottom) 
        id_materials = [int(self.idmat_top), int(self.idmat_face_top), int(self.idmat_face_bottom), int(self.idmat_side), int(self.idmat_bottom)]
        
        # depth at bottom
        bottom_z = self.bottom_z
        if self.steps_type == 'OPEN':
            # depth at front
            bottom_z = self.nose_z
        
        width_left  = 0.5*self.width - self.x_offset
        width_right = 0.5*self.width + self.x_offset   
        
        g = StairGenerator()
        if self.presets == 'STAIR_USER':
            for part in self.parts:
                g.add_part(part.type, self.steps_type, self.nose_type, self.z_mode, self.nose_z, bottom_z, center, max(width_left+0.01,
                        width_right+0.01, part.radius), part.da, width_left, width_right, part.length, part.left_shape, part.right_shape)
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
                g.add_part('D_STAIR', self.steps_type, self.nose_type, self.z_mode, self.nose_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), dir*pi,  width_left, width_right, 1.0, self.left_shape, self.right_shape)
            print("dir%s, last_da:%s abs_last:%s n_parts:%s" % (dir, last_da, abs_last, n_parts) )
            if round(abs_last, 2) > 0:
                if abs_last > pi/2:
                    g.add_part('C_STAIR', self.steps_type, self.nose_type, self.z_mode, self.nose_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), dir*pi/2,  width_left, width_right, 1.0, self.left_shape, self.right_shape)
                    g.add_part('C_STAIR', self.steps_type, self.nose_type, self.z_mode, self.nose_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), last_da - dir*pi/2,  width_left, width_right, 1.0, self.left_shape, self.right_shape)
                else:
                    g.add_part('C_STAIR', self.steps_type, self.nose_type, self.z_mode, self.nose_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), last_da,  width_left, width_right, 1.0, self.left_shape, self.right_shape)
        else:
            for part in self.parts:
                g.add_part(part.type, self.steps_type, self.nose_type, self.z_mode, self.nose_z, bottom_z, center, max(width_left+0.01, 
                            width_right+0.01, self.radius), self.da,  width_left, width_right, part.length, self.left_shape, self.right_shape)
        
        # Stair basis
        g.set_matids(id_materials)
        g.make_stair(self.height, self.step_depth, verts, faces, matids, uvs, nose_y=self.nose_y)
        
        # Ladder
        offset_x =  0.5 * self.width - self.post_offset_x
        if self.left_post or self.right_post:
            post_spacing = self.post_spacing
            if self.post_corners:
                post_spacing = 10000
            g.make_post(self.height, self.step_depth, 0.5*self.post_x, 0.5*self.post_y, self.post_z, self.post_alt, post_spacing, 0, self.post_corners, 
                  self.x_offset, offset_x, self.right_post, self.left_post, int(self.idmat_post), True, verts, faces, matids, uvs)

        if self.left_panel or self.right_panel:
            post_spacing = self.post_spacing
            if self.post_corners:
                post_spacing = 10000
            g.make_post(self.height, self.step_depth, 0.5*self.panel_x, 0, self.panel_z, self.panel_alt, post_spacing, self.panel_dist, self.post_corners, 
                    self.x_offset, offset_x, self.right_panel, self.left_panel, int(self.idmat_panel), False, verts, faces, matids, uvs)

        if self.right_ladder:
            for i in range(self.ladder_n):
                id_materials = [int(self.ladder_mat[i].index) for j in range(5)]
                g.set_matids(id_materials)
                g.make_part(self.height, self.step_depth, self.ladder_x[i], self.ladder_z[i], 
                        offset_x + self.ladder_offset[i], self.ladder_alt[i], 'LINEAR', 'CLOSED', verts, faces, matids, uvs)
        
        if self.left_ladder:
            for i in range(self.ladder_n):
                id_materials = [int(self.ladder_mat[i].index) for j in range(5)]
                g.set_matids(id_materials)
                g.make_part(self.height, self.step_depth, self.ladder_x[i], self.ladder_z[i], 
                        -offset_x - self.ladder_offset[i], self.ladder_alt[i], 'LINEAR', 'CLOSED', verts, faces, matids, uvs)
        
        handrail = [3*Vector(x) for x in reversed([(-0.0038, 0.0375), (-0.0019, 0.0381), (0.0, 0.0383), 
            (0.002, 0.0381), (0.0038, 0.0375), (0.0056, 0.0366), (0.0071, 0.0354), (0.0083, 0.0339), 
            (0.0092, 0.0321), (0.0098, 0.0302), (0.01, 0.0283), (0.0098, 0.0263), (0.0092, 0.0245), 
            (0.0083, 0.0227), (0.0071, 0.0212), (0.0056, 0.02), (0.0051, 0.0185), (0.0066, 0.0171), 
            (0.01, 0.0171), (0.01, 0.0), (-0.01, 0.0), (-0.01, 0.0171), (-0.0066, 0.0171), (-0.0051, 0.0185), 
            (-0.0056, 0.02), (-0.0071, 0.0212), (-0.0083, 0.0227), (-0.0092, 0.0245), (-0.0098, 0.0263), 
            (-0.01, 0.0283), (-0.0098, 0.0302), (-0.0092, 0.0321), (-0.0083, 0.0339), (-0.0071, 0.0354), (-0.0056, 0.0366)])]
        
        
        if self.right_handrail:
            g.make_profile(handrail, int(self.idmat_handrail), "RIGHT", self.handrail_slice_right, self.height, self.step_depth, 
                self.x_offset + offset_x + self.handrail_offset, self.handrail_alt, self.handrail_extend, verts, faces, matids, uvs)
        
        if self.left_handrail:
            g.make_profile(handrail, int(self.idmat_handrail), "LEFT", self.handrail_slice_left, self.height, self.step_depth, 
                -self.x_offset + offset_x + self.handrail_offset, self.handrail_alt, self.handrail_extend, verts, faces, matids, uvs)
         
        w = 0.5*self.string_x
        h = self.string_z
        string = [Vector((-w, 0)),Vector((w, 0)),Vector((w, h)),Vector((-w, h))] 
         
        if self.right_string:
            g.make_profile(string, int(self.idmat_string), "RIGHT", False, self.height, self.step_depth, 
                self.x_offset + 0.5 * self.width + self.string_offset, self.string_alt, 0, verts, faces, matids, uvs)
        
        if self.left_string:
            g.make_profile(string, int(self.idmat_string), "LEFT", False, self.height, self.step_depth, 
                -self.x_offset + 0.5 * self.width + self.string_offset, self.string_alt, 0, verts, faces, matids, uvs)
        
         
         
        bmed.buildmesh(context, o, verts, faces, matids=matids, uvs=uvs, weld=False, clean=False)
        
        
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
        row.operator('archipack.stair_manipulate')
        row = layout.row(align=True)
        row.prop(prop, 'presets')
        box = layout.box()
        box.prop(prop, 'width')
        box.prop(prop, 'height')
        box.prop(prop, 'bottom_z')
        box.prop(prop, 'x_offset')
        #box.prop(prop, 'z_mode')
        box = layout.box()
        box.prop(prop, 'steps_type')
        box.prop(prop, 'step_depth')
        box.prop(prop, 'nose_type')
        box.prop(prop, 'nose_z')
        box.prop(prop, 'nose_y')
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
        if prop.balluster_expand:
            row.prop(prop, 'balluster_expand',icon="TRIA_DOWN", icon_only=True, text="Components", emboss=False)
            box = layout.box()
            row = box.row()
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
            box.label(text="Handrail")
            row = box.row(align=True)
            row.prop(prop, 'left_handrail')
            row.prop(prop, 'right_handrail')
            if prop.left_handrail or prop.right_handrail:
                box.prop(prop, 'handrail_alt')
                box.prop(prop, 'handrail_offset')
                box.prop(prop, 'handrail_extend')
                row = box.row(align=True)
                row.prop(prop, 'handrail_slice_left')
                row.prop(prop, 'handrail_slice_right')
            box = layout.box()
            box.label(text="String")
            row = box.row(align=True)
            row.prop(prop, 'left_string')
            row.prop(prop, 'right_string')
            if prop.left_string or prop.right_string:
                box.prop(prop, 'string_x')
                box.prop(prop, 'string_z')
                box.prop(prop, 'string_alt')
                box.prop(prop, 'string_offset')
            box = layout.box()
            box.label(text="Post")
            row = box.row(align=True)
            row.prop(prop, 'left_post')
            row.prop(prop, 'right_post')
            if prop.left_post or prop.right_post:
                box.prop(prop, 'post_spacing')
                box.prop(prop, 'post_x')
                box.prop(prop, 'post_y')
                box.prop(prop, 'post_z')
                box.prop(prop, 'post_alt')
                box.prop(prop, 'post_offset_x')
                box.prop(prop, 'post_corners')
            box = layout.box()
            box.label(text="Panels")
            row = box.row(align=True)
            row.prop(prop, 'left_panel')
            row.prop(prop, 'right_panel')
            if prop.left_panel or prop.right_panel:
                box.prop(prop, 'panel_dist')
                box.prop(prop, 'panel_x')
                box.prop(prop, 'panel_z')
                box.prop(prop, 'panel_alt')
        else:
            row.prop(prop, 'balluster_expand',icon="TRIA_RIGHT", icon_only=True, text="Components", emboss=False)
        box = layout.box()
        row = box.row()
        
        if prop.idmats_expand:
            row.prop(prop, 'idmats_expand',icon="TRIA_DOWN", icon_only=True, text="Materials", emboss=False)
            box.prop(prop, 'idmat_top')
            box.prop(prop, 'idmat_face_top')
            box.prop(prop, 'idmat_face_bottom')
            box.prop(prop, 'idmat_bottom')
            box.prop(prop, 'idmat_side')
            box.prop(prop, 'idmat_handrail')
            box.prop(prop, 'idmat_panel')
            box.prop(prop, 'idmat_post')
            box.prop(prop, 'idmat_string')
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
        o.location = bpy.context.scene.cursor_location
        o.select = True
        context.scene.objects.active = o 
        bpy.ops.archipack.stair_manipulate('INVOKE_DEFAULT')
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

# ------------------------------------------------------------------
# Define operator class to manipulate object
# ------------------------------------------------------------------       
class ARCHIPACK_OT_stair_manipulate(Operator):
    bl_idname = "archipack.stair_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}
    
    @classmethod
    def poll(self, context):
        return ARCHIPACK_PT_stair.filter(context.active_object)
        
    def modal(self, context, event):
        return self.manips.modal(context, event)
        
    def update(self, context):
        self.manips.exit()
        o = context.active_object
        # context, axis, location, size, size_x, size_y, action
        obj = o.data.StairProperty[0]
        tM = o.matrix_world
        self.manips.add(context, "X", "TOP", 15, 15, obj, "width", min=0.01, max=50, step_round=2, sensitive=2, label="width", pos = tM * Vector((0.5*obj.width, 0, 0)))
        if obj.presets == 'STAIR_O':
            pos = tM * Vector((obj.radius-0.5*obj.width, obj.radius, 0))
            self.manips.add(context, "Y", "TOP", 15, 15, obj, "radius", min=0.01, max=50, step_round=2, sensitive=1, label="radius", pos=pos )
        elif obj.presets == 'STAIR_I':
            pos = tM * Vector((0.5*obj.width, obj.parts[0].length, obj.height))
            self.manips.add(context, "Y", "TOP", 15, 15, obj.parts[0], "length", min=0.01, max=50, step_round=2, sensitive=1, label="length", pos=pos)
        elif obj.presets == 'STAIR_L':
            pos = tM * Vector((0.5*obj.width, obj.parts[0].length, 0.3*obj.height))
            self.manips.add(context, "Y", "TOP", 15, 15, obj.parts[0], "length", min=0.01, max=50, step_round=2, sensitive=1, label="length", pos=pos)
            pos = tM * Vector((0.5*obj.width, obj.parts[0].length+obj.radius+0.5*obj.width, 0.5*obj.height))
            self.manips.add(context, "Y", "TOP", 15, 15, obj, "radius", min=0.01, max=50, step_round=2, sensitive=2, label="radius", pos=pos)
            #pos = tM * Vector((0.5*obj.width, obj.parts[0].length, 0))
            #self.manips.add(context, "Y", "BOTTOM", 15, 15, obj.parts[2], "length", min=0.01, max=50, step_round=2, sensitive=2, label="length", pos=pos)
        elif obj.presets == 'STAIR_U':
            pos = tM * Vector((0.5*obj.width, obj.parts[0].length, 0.3*obj.height))
            self.manips.add(context, "Y", "TOP", 15, 15, obj.parts[0], "length", min=0.01, max=50, step_round=2, sensitive=1, label="length", pos=pos)
            pos = tM * Vector((0.5*obj.width, obj.parts[0].length+obj.radius+0.5*obj.width, 0.5*obj.height))
            self.manips.add(context, "Y", "TOP", 15, 15, obj, "radius", min=0.01, max=50, step_round=2, sensitive=1, label="radius", pos=pos)
            #pos = tM * Vector((0.5*obj.width, obj.parts[0].length, 0))
            #self.manips.add(context, "Y", "BOTTOM", 15, 15, obj.parts[2], "length", min=0.01, max=50, step_round=2, sensitive=2, label="length", pos=pos)
        self.manips.add(context, "Z", "TOP", 15, 15, obj, "height", min=0.01, max=50, step_round=2, sensitive=1, label="height")
        
    def invoke(self, context, event):
        if context.space_data.type == 'VIEW_3D':
            self.manips = ManipulatorStack(context, self.update)
            self.update(context)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'} 
             
             
            
bpy.utils.register_class(StairMaterialProperty)
bpy.utils.register_class(StairPartProperty)
bpy.utils.register_class(StairProperty)
Mesh.StairProperty = CollectionProperty(type=StairProperty)
bpy.utils.register_class(ARCHIPACK_PT_stair)
bpy.utils.register_class(ARCHIPACK_OT_stair)



