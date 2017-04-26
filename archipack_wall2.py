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

# ----------------------------------------------------------
# Author: Stephen Leger (s-leger)
#
# ----------------------------------------------------------
# noinspection PyUnresolvedReferences
import bpy
# noinspection PyUnresolvedReferences
from bpy.types import Operator, PropertyGroup, Mesh, Panel
from bpy.props import FloatProperty, BoolProperty, IntProperty, StringProperty, FloatVectorProperty, CollectionProperty , EnumProperty
from .bmesh_utils import BmeshEdit as bmed
from .materialutils import MaterialUtils
from mathutils import Vector, Matrix
from math import sin, cos, pi, atan2, sqrt
from .archipack_manipulator import Manipulable, ManipulatorProperty

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
    def oposite(self):
        return Line(self.p+self.v, -self.v)
    @property
    def cross_z(self):
        return Vector((self.v.y, - self.v.x))
    def normal(self, t=0):
        # perpendiculaire a droite du segment
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
        t = (c * (line.p-self.p))/d
        return True, self.lerp(t), t
    def steps(self, len):
        return 1, 1
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
        return Line(self.lerp(t), self.v.normalized()*length)
    def rotate(self, da):
        cs = cos(da)
        sn = sin(da)
        x,y = self.v
        self.v.x = x*cs - y*sn
        self.v.y = x*sn + y*cs
        return self
    def scale(self, length):
        self.v = length*self.v.normalized()
        return self
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
    def gl_pts(self):
        p0 = self.p.to_3d()
        p1 = (self.p+self.v).to_3d()
        return [p0, p1]
        
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
            return False, line.lerp(1), 1
        elif d == 0:
            t = -B / 2*A
            return True, line.lerp(t), t
        else:
            AA = 2*A
            dsq = sqrt(d)
            t0 = (-B + dsq)/AA
            t1 = (-B - dsq)/AA
            if abs(abs(t0)-1) < abs(abs(t1)-1):
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
    def sized_normal(self, t, size):
        p = self.lerp(t)
        if self.da < 0:
            v = self.c-p
        else:
            v = p-self.c
        return Line(p, length*v.normalized())
    def lerp(self, t):
        a = self.a0+t*self.da
        return self.c+Vector((self.r*cos(a),self.r*sin(a)))
    def steps(self, step_angle):
        steps = max(1, round(abs(self.da) / step_angle, 0))
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
    def straight(self, length, t=1):
        return self.tangeant(t,length)
    def rotate(self, da):
        cs = cos(da)
        sn = sin(da)
        x,y = self.v
        self.v.x = x*cs - y*sn
        self.v.y = x*sn + y*cs
        return self
    def point_sur_segment(self, pt):
        dp = pt-self.c
        d = dp.length-self.r
        a = atan2(dp.y, dp.x)
        t = (a - self.a0)/self.da
        return t > 0 and t < 1, d, t
    def gl_pts(self):
        n_pts = max(1, int(round(abs(self.da)/pi*30,0)))
        t_step = 1/n_pts
        return [self.lerp(i*t_step).to_3d() for i in range(n_pts+1)]
           
class Wall():
    def __init__(self, last, z, t, flip):
        self.z = z
        self.t = t
        self.flip = flip
        self.z_step = len(z)
    def get_z(self, t):
        t0 = self.t[0]
        z0 = self.z[0]
        for i in range(1, self.z_step):
            t1 = self.t[i]
            z1 = self.z[i]
            if t <= t1:
                return z0 + (t - t0) / (t1-t0) * (z1-z0)     
            t0, z0 = t1, z1
        return self.z[-1]
    def make_faces(self, i, f, faces):
        if i < self.n_step:
            # 1 3   5 7
            # 0 2   4 6
            if self.flip:
                faces.append((f+2, f, f+1, f+3))
            else:    
                faces.append((f, f+2, f+3, f+1))
    def p3d(self, verts, t):
        x, y = self.lerp(t)
        z = self.get_z(t)
        verts.append((x, y, 0))
        verts.append((x, y, z))
    def make_wall(self, i, verts, faces):
        t = self.t_step[i]
        f = len(verts)
        self.p3d(verts, t)
        self.make_faces(i, f, faces)
    def straight_wall(self, last, length, da, z, t):
        r = last.straight(length).rotate(da)
        return StraightWall(last, r.p, r.v, z, t, self.flip)
    def curved_wall(self, last, a0, da, radius, z, t):
        n = last.normal(1).rotate(a0).scale(radius)
        if da < 0:
            n.v = -n.v
        a0 = n.angle
        c = n.p - n.v
        return CurvedWall(last, c, radius, a0, da, z, t, self.flip)
        
class StraightWall(Wall, Line):
    def __init__(self, last, p, v, z, t, flip):
        Line.__init__(self, p, v)
        Wall.__init__(self, last, z, t, flip)
    
    def param_t(self, step_angle):
        self.t_step = self.t
        self.n_step = len(self.t)-1
        
class CurvedWall(Wall, Arc):
    def __init__(self, last, c, radius, a0, da, z, t, flip):
        Arc.__init__(self, c, radius, a0, da)
        Wall.__init__(self, last, z, t, flip)
    
    def param_t(self, step_angle):
        t_step, n_step = self.steps(step_angle) 
        self.t_step = list(sorted([i*t_step for i in range(1, n_step)] + self.t))
        self.n_step = len(self.t_step)-1
                 
class WallGenerator():
    def __init__(self, parts):
        self.last_type = 'NONE'
        self.walls = []
        self.parts = parts
        self.faces_type = 'NONE'
    def add_part(self, type, center, radius, a0, da, length, part_z, part_t, n_splits, flip):
        
        # TODO:
        # refactor this part (height manipulators)
        manip_index = []
        if len(self.walls) < 1:
            s = None
            z = [part_z[0]]
            manip_index.append(0)
        else:
            s = self.walls[-1]
            z = [s.z[-1]]
        
        t_cur = 0
        z_last = n_splits-1
        t = [0]
        
        for i in range(n_splits):
            t_try = t[-1] + part_t[i] 
            if t_try == t_cur:
                continue
            if t_try <= 1:
                t_cur = t_try
                t.append(t_cur)
                z.append(part_z[i])
                manip_index.append(i)
            else:
                z_last = i
                break
                
        if t_cur < 1:
            t.append(1)
            manip_index.append(z_last)
            z.append(part_z[z_last])
            
        # start a new stair    
        if s is None:
            if type == 'S_WALL':
                p = Vector((0,0))
                v = length*Vector((cos(a0),sin(a0)))
                s = StraightWall(s, p, v, z, t, flip)
            elif type == 'C_WALL':
                if da < 0:
                    c = Vector((radius,0))
                else:
                    c = Vector((-radius,0))
                s = CurvedWall(s, c, radius, a0, da, z, t, flip)
        else:
            if type == 'S_WALL':
                s = s.straight_wall(s, length, a0, z, t)
            elif type == 'C_WALL':
                s = s.curved_wall(s, a0, da, radius, z, t)
        self.walls.append(s)
        self.last_type = type
        return manip_index
        
    def make_wall(self, step_angle, verts, faces):
        for i, wall in enumerate(self.walls):
            manipulators = self.parts[i].manipulators
            # angle from last to current segment
            if i > 0:
                v0 = self.walls[i-1].straight(-1, 1).v.to_3d()
                v1 = wall.straight(1, 0).v.to_3d()
                manipulators[0].set_pts([wall.lerp(0).to_3d(), v0, v1])
            if type(wall).__name__ == "StraightWall":
                # segment length
                manipulators[1].type = 'SIZE'
                manipulators[1].prop1_name = "length"
                manipulators[1].set_pts([wall.lerp(0).to_3d(), wall.lerp(1).to_3d(), (1,0,0)])
            else:
                # segment radius + angle
                v0 = (wall.lerp(0)-wall.c).to_3d()
                v1 = (wall.lerp(1)-wall.c).to_3d()
                manipulators[1].type = 'ARC_ANGLE_RADIUS'
                manipulators[1].prop1_name = "da"
                manipulators[1].prop2_name = "radius"
                manipulators[1].set_pts([wall.c.to_3d(), v0, v1])
                
            wall.param_t(step_angle)
            for j in range(wall.n_step+1):
                wall.make_wall(j, verts, faces)
                
    def debug(self, verts):
        for wall in self.walls:
            for i in range(33):
                x, y = wall.lerp(i/32)
                verts.append((x,y,0))
            
def update(self, context):
    self.update(context)

def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)
 
def set_splits(self, value):
    if self.n_splits != value:
        self.auto_update = False
        self._set_t(value)  
        self.auto_update = True
        self.n_splits = value
    return None
    
def get_splits(self):
    return self.n_splits
  
class Wall2PartProperty(PropertyGroup):
    type = EnumProperty(
            items=(
                ('S_WALL', 'Straight','',0),
                ('C_WALL', 'Curved','',1)
                ),
            default='S_WALL',
            update=update_manipulators
            )
    length = FloatProperty(
            name="length",
            min=0.01,
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
    a0 = FloatProperty(
            name="angle depart",
            min=-pi,
            max=pi,
            default=pi/2,
            subtype='ANGLE', unit='ROTATION',
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
    z = FloatVectorProperty(
            name="height",
            min=0,
            max=1000,
            default=[
                    2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,
                    2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,
                    2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,
                    2.4,2.4,2.4,2.4,2.4,2.4,1
                ],
            size=31,
            update=update
            )
    t = FloatVectorProperty(
            name="position",
            min=0,
            max=1,
            default=[
                1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1
            ],
            size=31,
            update=update
            )
    splits = IntProperty(
        name="splits",
        default=1,
        min=1,
        max=31, 
        get=get_splits, set=set_splits
        )
    n_splits= IntProperty(
        name="splits",
        default=1,
        min=1,
        max=31, 
        update=update
        )
    auto_update=BoolProperty(default=True)
    manipulators = CollectionProperty(type=ManipulatorProperty)
    # ui related
    expand = BoolProperty(default=False)
    def _set_t(self, splits):
        t = 1/splits
        for i in range(splits):
            self.t[i] = t
    
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        selected = [o for o in context.selected_objects]
        for o in selected:
            props = ARCHIPACK_PT_wall2.params(o)
            if props:
                for part in props.parts:
                    if part == self:
                        return props
        return None
        
    def update(self, context, manipulable_refresh=False):
        if not self.auto_update:
            return 
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context, manipulable_refresh)
        
    def draw(self, layout, context, index):
        
        row = layout.row(align=True)
        if self.expand:
            row.prop(self, 'expand',icon="TRIA_DOWN", icon_only=True, text="Part "+str(index+1), emboss=False)
        else:
            row.prop(self, 'expand',icon="TRIA_RIGHT", icon_only=True, text="Part "+str(index+1), emboss=False)
        #row.label(text="Part "+str(index+1))
        row.prop(self, "type", text="")
            
        if self.expand:        
            row = layout.row(align=True)
            row.operator("archipack.wall2_insert", text="Split").index = index
            row.operator("archipack.wall2_remove", text="Remove").index = index
            if self.type == 'C_WALL':
                row = layout.row()
                row.prop(self, "radius")
                row = layout.row()
                row.prop(self, "da")
            else:
                row = layout.row()
                row.prop(self, "length")
            row = layout.row()
            row.prop(self, "a0")
            row = layout.row()
            row.prop(self, "splits")
            for split in range(self.n_splits):
                row = layout.row()
                row.prop(self, "z", text="alt", index=split)
                row.prop(self, "t", text="pos", index=split)

class Wall2ChildProperty(PropertyGroup):
    # Size  Loc
    # Delta Loc
    manipulators=CollectionProperty(type=ManipulatorProperty)
    name = StringProperty()
    wall_idx = IntProperty()
    pos = FloatVectorProperty(subtype='XYZ')
    flip = BoolProperty(default=False)
    
    def get_child(self, context):
        d = None
        child = context.scene.objects.get(self.name)
        if child.data is not None:
            if 'archipack_window' in child.data:
                d = child.data.archipack_window[0]
            elif 'archipack_door' in child.data:
                d = child.data.archipack_door[0]
        return child, d
       
class Wall2Property(Manipulable, PropertyGroup):
    parts = CollectionProperty(type=Wall2PartProperty)
    n_parts = IntProperty(
            name="parts",
            min=1,
            max=32,
            default=1, update=update_manipulators
            )
    step_angle = FloatProperty(
            name="step angle",
            min=1/180*pi,
            max=pi,
            default=6/180*pi,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
    width = FloatProperty(
            name="width",
            min=0.01,
            max=100.0,
            default=0.2,
            update=update
            )
    x_offset = FloatProperty(
        name="x offset",
            min=-1, max=1,
            default=-1, precision=2, step=1,
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
    flip=BoolProperty(default=False, update=update)
    closed=BoolProperty(
            default=False,
            name="Close",
            update=update
            )
    auto_update=BoolProperty(default=True)
    realtime=BoolProperty(default=True, name="RealTime", description="Relocate childs in realtime")
    # dumb manipulators to show sizes between childs
    childs_manipulators=CollectionProperty(type=ManipulatorProperty)
    childs=CollectionProperty(type=Wall2ChildProperty)
    
    def insert_part(self, context, where):
        self.manipulable_disable()
        self.auto_update = False
        part_0 = self.parts[where]
        part_0.length /=2
        part_0.da /=2
        p = self.parts.add()
        s = p.manipulators.add()
        s.type = "ANGLE"
        s.prop1_name = "a0"
        s = p.manipulators.add()
        s.prop1_name = "length"
        part_1 = self.parts[len(self.parts)-1] 
        part_1.type = part_0.type
        #part_1.z_left = part_0.z_left
        #part_1.z_right = part_0.z_right
        part_1.length = part_0.length
        part_1.da = part_0.da
        part_1.a0 = part_0.a0
        part_0.a0 = 0
        self.parts.move(len(self.parts)-1, where)
        self.n_parts += 1
        self.auto_update = True
        self.update(context, manipulable_refresh = True)
        
    def remove_part(self, context, where):
        self.manipulable_disable()
        self.auto_update = False
        self.parts.remove(where)
        self.n_parts -= 1
        self.auto_update = True
        self.update(context, manipulable_refresh = True)
       
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        active = context.active_object
        selected = [o for o in context.selected_objects]
        #selected = [o for o in context.scene.objects]
        for o in selected:
            if ARCHIPACK_PT_wall2.params(o) == self:
                return active, selected, o
        return active, selected, None
    
    def get_generator(self):
        # print("get_generator")
        center = Vector((0,0))    
        g = WallGenerator(self.parts)
        for part in self.parts:
            g.add_part(part.type, center, part.radius, part.a0, part.da, part.length, part.z, part.t, part.n_splits, self.flip)
        return g
        
    def update_parts(self, o):
        # print("update_parts")
        # remove rows
        row_change = False
        for i in range(len(self.parts), self.n_parts, -1):
            row_change = True
            self.parts.remove(i-1)
        
        # add rows
        for i in range(len(self.parts), self.n_parts):
            row_change = True
            p = self.parts.add()
            s = p.manipulators.add()
            s.type = "ANGLE"
            s.prop1_name = "a0"
            s = p.manipulators.add()
            s.type = "SIZE"
            s.prop1_name = "length"   
        
        g = self.get_generator()
        
        if row_change:
            self.setup_childs(o, g)
        
        return g
        
    def update(self, context, manipulable_refresh=False):   
        # print("update manipulable_refresh:%s" % (manipulable_refresh))
            
        if not self.auto_update:
            return
            
        active, selected, o = self.find_in_selection(context)
        
        if o is None:
            return
        
        if manipulable_refresh:
            # prevent crash by removing all manipulators refs to datablock before changes 
            self.manipulable_disable(context)
            
        verts = []
        faces = []
        
        g = self.update_parts(o)
        # print("make_wall")
        g.make_wall(self.step_angle, verts, faces)
        
        if self.closed:
            f = len(verts)
            faces.append((f-2, 0, 1, f-1))
            
        # print("buildmesh")    
        bmed.buildmesh(context, o, verts, faces, matids=None, uvs=None, weld=True)
        
        # Width
        self.manipulators[0].set_pts([(0,0,0),(-self.width,0,0),(-1,0,0)])
        # Parts COUNTER
        self.manipulators[1].set_pts([g.walls[-1].lerp(1.1).to_3d(),g.walls[-1].lerp(1.1+0.5/g.walls[-1].length).to_3d(),(-1,0,0)])

        if self.realtime:
            # update child location and size
            self.relocate_childs(context, o, g)
            # store gl points
            self.update_childs(context, o, g)
        else:
            bpy.ops.archipack.wall2_delay_update(name=o.name)
            
        modif = o.modifiers.get('Wall')
        if modif is None:
            modif = o.modifiers.new('Wall', 'SOLIDIFY')
            modif.use_quality_normals = True
            modif.use_even_offset = True
            modif.material_offset_rim = 2
            modif.material_offset = 1
        
        modif.thickness = self.width
        modif.offset = self.x_offset
        
        if manipulable_refresh:
            # print("manipulable_refresh=True")
            self.manipulable_refresh=True
        
        # restore context
        try:
            for o in selected:
                o.select = True
        except:
            pass    
        
        active.select = True
        context.scene.objects.active = active  
    
    # manipulable 
    def child_partition(self, array, begin, end):
        pivot = begin
        for i in range(begin+1, end+1):
            # wall idx
            if array[i][0] < array[begin][0]:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
            # param t on the wall
            elif array[i][0] == array[begin][0] and array[i][3] <= array[begin][3]:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
        array[pivot], array[begin] = array[begin], array[pivot]
        return pivot

    def sort_child(self, array, begin=0, end=None):
        # print("sort_child")
        if end is None:
            end = len(array) - 1
        def _quicksort(array, begin, end):
            if begin >= end:
                return
            pivot = self.child_partition(array, begin, end)
            _quicksort(array, begin, pivot-1)
            _quicksort(array, pivot+1, end)
        return _quicksort(array, begin, end)
    
    def add_child(self, name, wall_idx, pos, flip):
        # print("add_child %s %s" % (name, wall_idx))
        c = self.childs.add()
        c.name = name
        c.wall_idx = wall_idx
        c.pos = pos
        c.flip = flip
        m = c.manipulators.add()
        m.type = 'DELTA_LOC'
        m.prop1_name = "x"
        m = c.manipulators.add()
        m.type = 'SIZE_LOC'
        m.prop1_name = "x"
        m.prop2_name = "x"
    
    def setup_childs(self, o, g):
        """
            Store childs
            create manipulators
            call after a boolean oop
        """
        # print("setup_childs")
            
        self.childs.clear()
        self.childs_manipulators.clear()
        if o.parent is None:
            return
        wall_with_childs = [0 for i in range(self.n_parts)]    
        relocate = []
        dmax = 2*self.width   
        itM = o.matrix_world.inverted() * o.parent.matrix_world
        rM = itM.to_3x3()
        for child in o.parent.children:
            if (child != o and 
                'archipack_robusthole' not in child and 
                'archipack_hole' not in child):
                tM = child.matrix_world.to_3x3()
                pt = (itM * child.location).to_2d()
                dir_y = (rM * tM * Vector((0,1,0))).to_2d()
                for wall_idx, wall in enumerate(g.walls):
                    res, d, t = wall.point_sur_segment(pt)
                    dir = -wall.normal(t).v.normalized()
                    if res and t > 0 and t < 1 and abs(d) < dmax:
                        wall_with_childs[wall_idx] = 1
                        m = self.childs_manipulators.add()
                        m.type = 'DUMB_SIZE'
                        relocate.append((child.name, wall_idx, (t*wall.length, d, child.location.z), (dir-dir_y).length > 0.5))
        self.sort_child(relocate)
        for child in relocate:
            name, wall_idx, pos, flip = child
            self.add_child(name, wall_idx, pos, flip)
            
        for i in range(sum(wall_with_childs)):
            m = self.childs_manipulators.add()
            m.type = 'DUMB_SIZE'
        
    def relocate_childs(self, context, o, g):
        """
            Move and resize childs after wall edition
        """
        ## print("relocate_childs")
        
        w = self.width
        if self.flip:
            w = -w
        tM = o.matrix_world
        for child in self.childs:
            c, d = child.get_child(context)
            if c is None:
                continue
            t = child.pos.x/g.walls[child.wall_idx].length        
            n = g.walls[child.wall_idx].normal(t)
            rx, ry = -n.v.normalized()
            rx, ry = ry, -rx
            if child.flip:
                rx, ry = -rx, -ry
            #context.scene.objects.active = c
            if d is not None:
                c.select = True
                d.y = w
                d.update(context)
                c.select = False
                x, y = n.p - (0.5*w * n.v.normalized())
            else:    
                x, y = n.p - (child.pos.y * n.v.normalized())
            #c.select = False
            context.scene.objects.active = o
            # preTranslate
            c.matrix_world = tM * Matrix([
                [rx,-ry,0,x],
                [ry,rx,0,y],
                [0,0,1,child.pos.z],
                [0,0,0,1]
            ])
     
    def update_childs(self, context, o, g):
        """
            setup gl points for childs
        """
        # print("update_childs")
        
        if o.parent is None:
            return
            
        itM = o.matrix_world.inverted() * o.parent.matrix_world
        m_idx = 0
        for wall_idx, wall in enumerate(g.walls):
            p0 = wall.lerp(0)
            wall_has_childs = False
            for child in self.childs:
                if child.wall_idx == wall_idx:
                    c, d = child.get_child(context)
                    if d is not None:
                        wall_has_childs = True
                        dt = 0.5*d.x/wall.length
                        pt = (itM * c.location).to_2d()
                        res, y, t = wall.point_sur_segment(pt)
                        child.pos = (wall.length*t, y, child.pos.z)
                        p1 = wall.lerp(t-dt)
                        # dumb
                        self.childs_manipulators[m_idx].set_pts([(p0.x, p0.y, 0), (p1.x, p1.y, 0), (0.5,0,0)])
                        m_idx += 1
                        # delta loc
                        x, y = 0.5*d.x, 0.5*d.y
                        if child.flip:
                            side = -1
                        else:
                            side = 1      
                        child.manipulators[0].set_pts([(-x,side*-y,0),(x,side*-y,0),(side,0,0)])
                        # loc size
                        child.manipulators[1].set_pts([(-x,side*-y,0),(x,side*-y,0),(0.5*side,0,0)])
                        p0 = wall.lerp(t+dt)
            p1 = wall.lerp(1)  
            if wall_has_childs:
                self.childs_manipulators[m_idx].set_pts([(p0.x, p0.y, 0), (p1.x, p1.y, 0), (0.5,0,0)])
                m_idx += 1
                  
    def manipulate_childs(self, context):
        """
            setup child manipulators
        """
        # print("manipulate_childs")
            
        for wall_idx in range(self.n_parts):
            for child in self.childs:
                if child.wall_idx == wall_idx:
                    c, d = child.get_child(context)
                    if d is not None:
                        self.manip_stack.append(child.manipulators[0].setup(context, c, d))
                        self.manip_stack.append(child.manipulators[1].setup(context, c, d))
    
    def manipulable_manipulate(self, context, type):
        # print("manipulable_manipulate %s" % (type))
        if type in ['DeltaLocationManipulator', 'SizeLocationManipulator']:
            # update manipulators pos of childs 
            o = context.active_object
            if o.parent is None:
                return
            g = self.get_generator()
            itM = o.matrix_world.inverted() * o.parent.matrix_world
            for child in self.childs:
                c, d = child.get_child(context)
                if d is not None:
                    wall = g.walls[child.wall_idx]
                    pt = (itM * c.location).to_2d()
                    res, d, t = wall.point_sur_segment(pt)
                    child.pos = (t*wall.length, d, child.pos.z)
            if not self.realtime:
                self.update_childs(context, o, g)
            # trigger update by hand since those manipulations 
            # dosen't affect wall data
            self.update(context)
      
    def manipulable_setup(self, context):
        # print("manipulable_setup")
        self.manipulable_disable(context)
        o = context.active_object
        d = self
        
        # setup childs manipulators
        self.manipulate_childs(context)
        
        for i, part in enumerate(d.parts):
            if i >= d.n_parts:
                break
            # length / radius+angle
            self.manip_stack.append(part.manipulators[1].setup(context, o, part))
            if i > 0:
                # start angle 
                self.manip_stack.append(part.manipulators[0].setup(context, o, part))
        # width + counter
        for m in self.manipulators:
            self.manip_stack.append(m.setup(context, o, self))  
        # dumb between childs
        for m in self.childs_manipulators:
            self.manip_stack.append(m.setup(context, o, self))
            
    def manipulable_exit(self, context):
        return
    
    def manipulable_invoke(self, context):
        """
            call this in operator invoke() 
        """
        # print("manipulable_invoke")
        
        self.manip_stack = []
        o = context.active_object
        g = self.get_generator()
        # setup childs manipulators
        self.setup_childs(o, g)
        # store gl points
        self.update_childs(context, o, g)
        self.manipulable_release(context)
        self.manipulable_setup(context)
 
# use 2 globals to store a timer and state of update_action 
update_timer = None
update_timer_updating = False

class ARCHIPACK_OT_wall2_delay_update(Operator):
    bl_idname = "archipack.wall2_delay_update"
    bl_label = "Update childs with a delay"
    
    name = StringProperty()
    
    def modal(self, context, event):
        global update_timer_updating
        if event.type == 'TIMER' and not update_timer_updating:
            update_timer_updating = True
            o = context.scene.objects.get(self.name)
            #print("delay update of %s" % (self.name))
            if o is not None:
                o.select = True
                context.scene.objects.active = o
                d = o.data.archipack_wall2[0]
                g = d.get_generator()
                # update child location and size
                d.relocate_childs(context, o, g)
                # store gl points
                d.update_childs(context, o, g)
                self.cancel(context)
        return {'PASS_THROUGH'}

    def execute(self, context):
        global update_timer
        global update_timer_updating
        if update_timer is not None:
            if update_timer_updating:
                return {'CANCELLED'}
            # reset update_timer so it only occurs once 0.1s after last action
            context.window_manager.event_timer_remove(update_timer)
            update_timer = context.window_manager.event_timer_add(0.1, context.window)
            return {'CANCELLED'}
        update_timer_updating = False
        context.window_manager.modal_handler_add(self)
        update_timer = context.window_manager.event_timer_add(0.1, context.window)
        return {'RUNNING_MODAL'}

    def cancel(self, context):
        global update_timer
        context.window_manager.event_timer_remove(update_timer)
        update_timer = None
        return {'CANCELLED'}
   
class ARCHIPACK_PT_wall2(Panel):
    bl_idname = "ARCHIPACK_PT_wall2"
    bl_label = "Wall"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'
    
    def draw(self, context):
        o = context.object
        prop = ARCHIPACK_PT_wall2.params(o)
        if prop is None:
            return
        layout = self.layout
        row = layout.row(align=True)
        row.operator("archipack.wall2_manipulate")
        row.prop(prop, 'realtime')
        box = layout.box()
        box.prop(prop, 'n_parts')
        box.prop(prop, 'step_angle')
        box.prop(prop, 'width')
        box.prop(prop, 'flip')
        box.prop(prop, 'x_offset')
        row = layout.row()
        row.prop(prop, "closed")
        for i, part in enumerate(prop.parts):
            box = layout.box()
            part.draw(box, context, i)
        
    @classmethod
    def params(cls, o):
        try:
            if 'archipack_wall2' not in o.data:
                return False
            else:
                return o.data.archipack_wall2[0]
        except:
            return False
    @classmethod
    def filter(cls, o):
        try:
            if 'archipack_wall2' not in o.data:
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
class ARCHIPACK_OT_wall2(Operator):
    bl_idname = "archipack.wall2"
    bl_label = "Wall (alpha)"
    bl_description = "Wall (alpha)"
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
        m = bpy.data.meshes.new("Wall")
        o = bpy.data.objects.new("Wall", m)
        d = m.archipack_wall2.add()
        s = d.manipulators.add()
        s.prop1_name = "width"
        s = d.manipulators.add()
        s.prop1_name = "n_parts"
        s.type = 'COUNTER'
        context.scene.objects.link(o)
        o.select = True
        context.scene.objects.active = o 
        d.update(context)
        MaterialUtils.add_wall_materials(o)
        # select frame
        o.select = True
        context.scene.objects.active = o 
        bpy.ops.archipack.wall2_manipulate('INVOKE_DEFAULT')
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
# Define operator class to create object
# ------------------------------------------------------------------
class ARCHIPACK_OT_wall2_insert(Operator):
    bl_idname = "archipack.wall2_insert"
    bl_label = "Insert"
    bl_description = "Insert part"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    index = IntProperty(default=0)
    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')
    
    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            prop = ARCHIPACK_PT_wall2.params(o)
            if prop is None:
                return {'CANCELLED'}
            prop.insert_part(context, self.index)   
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}
   
class ARCHIPACK_OT_wall2_remove(Operator):
    bl_idname = "archipack.wall2_remove"
    bl_label = "Remove"
    bl_description = "Remove part"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    index = IntProperty(default=0)
    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')
    
    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            o = context.active_object
            prop = ARCHIPACK_PT_wall2.params(o)
            if prop is None:
                return {'CANCELLED'}
            prop.remove_part(context, self.index)   
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}
 
# ------------------------------------------------------------------
# Define operator class to manipulate object
# ------------------------------------------------------------------          
class ARCHIPACK_OT_wall2_manipulate(Operator):
    bl_idname = "archipack.wall2_manipulate"
    bl_label = "Manipulate"
    bl_description = "Manipulate"
    bl_options = {'REGISTER', 'UNDO'}
    
    @classmethod
    def poll(self, context):
        return ARCHIPACK_PT_wall2.filter(context.active_object)
    
    def modal(self, context, event):
        return self.d.manipulable_modal(context, event)
        
    def invoke(self, context, event):
        if context.space_data.type == 'VIEW_3D':
            o = context.active_object
            self.d = o.data.archipack_wall2[0]
            self.d.manipulable_invoke(context)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'} 
    
    def execute(self, context):
        """
            For use in boolean ops
        """
        if ARCHIPACK_PT_wall2.filter(context.active_object):
            o = context.active_object
            d = o.data.archipack_wall2[0]
            g = d.get_generator()
            d.setup_childs(o, g)
            d.update_childs(context, o, g)
            d.manipulable_release(context)
            d.update(context)
            o.select = True
            context.scene.objects.active = o
        return {'FINISHED'} 
        
bpy.utils.register_class(Wall2PartProperty)
bpy.utils.register_class(Wall2ChildProperty)
bpy.utils.register_class(Wall2Property)
Mesh.archipack_wall2 = CollectionProperty(type=Wall2Property)
bpy.utils.register_class(ARCHIPACK_PT_wall2)
bpy.utils.register_class(ARCHIPACK_OT_wall2)
bpy.utils.register_class(ARCHIPACK_OT_wall2_insert)
bpy.utils.register_class(ARCHIPACK_OT_wall2_remove)
bpy.utils.register_class(ARCHIPACK_OT_wall2_manipulate)
bpy.utils.register_class(ARCHIPACK_OT_wall2_delay_update)
