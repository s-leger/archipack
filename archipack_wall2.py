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
from .archipack_manipulator import Manipulator, position_2d_from_coord

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
    def straight(self, length):
        return Line(self.p+self.v, self.v.normalized()*length)
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
    def gl_pts(self, context, tM):
        p0 = self.p.to_3d()
        p1 = (self.p+self.v).to_3d()
        gl_0 = position_2d_from_coord(context, tM * p0)
        gl_1 = position_2d_from_coord(context, tM * p1)
        return [gl_0, gl_1]
        
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
    def straight(self, length):
        return self.tangeant(1,length)
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
    def gl_pts(self, context, tM):
        n_pts = int(round(self.da/pi*6,0))
        t_step = 1/n_pts
        return [position_2d_from_coord(context, tM * self.lerp(i*t_step).to_3d()) for i in range(n_pts+1)]
           
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
    def make_step(self, i, verts, faces):
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
    def __init__(self):
        self.last_type = 'NONE'
        self.walls = []
        self.faces_type = 'NONE'
    def add_part(self, type, center, radius, a0, da, length, part_z, part_t, n_splits, flip):
        
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
        for wall in self.walls:
            wall.param_t(step_angle)
            for i in range(wall.n_step+1):
                wall.make_step(i, verts, faces)
                
    def debug(self, verts):
        for wall in self.walls:
            for i in range(33):
                x, y = wall.lerp(i/32)
                verts.append((x,y,0))
            
def update(self, context):
    self.update(context)
 
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
            update=update
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
        
    def update(self, context):
        if not self.auto_update:
            return 
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context)
        
    def draw(self, layout, context, index):
        row = layout.row(align=True)
        row.label(text="Part "+str(index+1))
        row.prop(self, "type", text="")
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
        
class Wall2Property(PropertyGroup):
    parts = CollectionProperty(type=Wall2PartProperty)
    n_parts = IntProperty(
            name="parts",
            min=1,
            max=32,
            default=1, update=update
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
            default=1.2,
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
    
    def insert_part(self, context, where):
        self.auto_update = False
        part_0 = self.parts[where]
        part_0.length /=2
        part_0.da /=2
        self.parts.add()
        part_1 = self.parts[len(self.parts)-1] 
        part_1.type = part_0.type
        part_1.z_left = part_0.z_left
        part_1.z_right = part_0.z_right
        part_1.length = part_0.length
        part_1.da = part_0.da
        part_1.a0 = part_0.a0
        part_0.a0 = 0
        self.parts.move(len(self.parts)-1, where)
        self.n_parts += 1
        self.auto_update = True
        self.update(context)
        
    def remove_part(self, context, where):
        self.auto_update = False
        self.parts.remove(where)
        self.n_parts -= 1
        self.auto_update = True
        self.update(context)
       
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
        
        g = WallGenerator()
        for part in self.parts:
            g.add_part(part.type, center, part.radius, part.a0, part.da, part.length, part.z, part.t, part.n_splits, self.flip)
        g.make_wall(self.step_angle, verts, faces)
       
        if self.closed:
            f = len(verts)
            faces.append((f-2, 0, 1, f-1))
            
        bmed.buildmesh(context, o, verts, faces, matids=None, uvs=None, weld=True)
        
        modif = o.modifiers.get('Wall')
        if modif is None:
            modif = o.modifiers.new('Wall', 'SOLIDIFY')
            modif.use_quality_normals = True
            modif.use_even_offset = True
            modif.material_offset_rim = 2
            modif.material_offset = 1
        
        modif.thickness = self.width
        modif.offset = self.x_offset
        
        # restore context
        try:
            for o in selected:
                o.select = True
        except:
            pass    
        
        active.select = True
        context.scene.objects.active = active  
 
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
        layout.operator("archipack.wall2_manipulate")
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
            if 'Wall2Property' not in o.data:
                return False
            else:
                return o.data.Wall2Property[0]
        except:
            return False
    @classmethod
    def filter(cls, o):
        try:
            if 'Wall2Property' not in o.data:
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
        d = m.Wall2Property.add()
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
            else:
                pivot += 1
                array[i], array[pivot] = array[pivot], array[i]
        array[pivot], array[begin] = array[begin], array[pivot]
        return pivot

    def sort_child(self, array, begin=0, end=None):
        if end is None:
            end = len(array) - 1
        def _quicksort(array, begin, end):
            if begin >= end:
                return
            pivot = self.child_partition(array, begin, end)
            _quicksort(array, begin, pivot-1)
            _quicksort(array, pivot+1, end)
        return _quicksort(array, begin, end)
    
    def store_childs(self, o, g):
        dmax = 2*o.data.Wall2Property[0].width
        itM = o.matrix_world.inverted() * o.parent.matrix_world
        rM = itM.to_3x3()
        self.relocate = []
        for child in o.parent.children:
            if child != o:
                tM = child.matrix_world.to_3x3()
                pt = (itM * child.location).to_2d()
                dir_y = (rM * tM * Vector((0,1,0))).to_2d()
                for wall_idx, wall in enumerate(g.walls):
                    res, d, t = wall.point_sur_segment(pt)
                    dir = -wall.normal(t).v.normalized()
                    if res and abs(d) < dmax:
                        #print("%s dir:%s dir_y:%s d:%s len:%s" % (child.name, dir, dir_y, d, (dir-dir_y).length))
                        self.relocate.append((wall_idx, child.name, d, t, child.location.z, (dir-dir_y).length > 0.5))
        self.sort_child(self.relocate)
        
    def relocate_childs(self, context, o, g):
        if not hasattr(self, "relocate"):
            return
        tM = o.matrix_world
        for i, name, d, t, z, flip in self.relocate:
            child = context.scene.objects.get(name)
            if child is None:
                print("child is none")
                continue
            n = g.walls[i].normal(t)
            rx, ry = -n.v.normalized()
            rx, ry = ry, -rx
            if flip:
                rx, ry = -rx, -ry
                #print("flip %s rx:%s, ry:%s cxx:%s, cxy:%s cyx:%s, cyy:%s" %(child.name, rx, ry, child.matrix_world[0][0], child.matrix_world[0][1], child.matrix_world[1][0], child.matrix_world[1][1]) )
            #else:
            #    print("%s rx:%s, ry:%s cxx:%s, cxy:%s cyx:%s, cyy:%s" %(child.name, rx, ry, child.matrix_world[0][0], child.matrix_world[0][1], child.matrix_world[1][0], child.matrix_world[1][1]) )
            x, y = n.p - (d * n.v.normalized())
            # preTranslate
            child.matrix_world = tM * Matrix([
                [rx,-ry,0,x],
                [ry,rx,0,y],
                [0,0,1,z],
                [0,0,0,1]
            ])
            
    def update(self, context):
        size = 10
        tM = self.o.matrix_world
        props = self.o.data.Wall2Property[0]
        props.update_parts()
        center = Vector((0,0))
        for manip in self.manipulators:
            manip.exit()
        self.manipulators = []
        g = WallGenerator()
        manip_t = [g.add_part(part.type, center, part.radius, part.a0, part.da, part.length, part.z, part.t, part.n_splits, props.flip) for part in props.parts]
        
        self.relocate_childs(context, self.o, g)
        context.scene.update()
        self.store_childs(self.o, g)
        
        x, y = g.walls[-1].lerp(1.1)
        pos = tM * Vector((x, y, 1.5))
        manip = Manipulator(context, "X", "TOP", size, size, props, "n_parts", min=props.n_parts, max=31, step_round=0, sensitive=1, label="add parts", pos=pos)
        self.manipulators.append(manip)
                   
        for wall_idx, params_t in enumerate(manip_t):
            
            childs = [child for child in self.relocate if child[0] == wall_idx]
            
            part = props.parts[wall_idx]
            wall = g.walls[wall_idx]
            normal = wall.straight(1).v.to_3d()
            t0 = 0
            t_abs = []
            for t in part.t:
                t0 += t
                if t0 > 1:
                    t0 = 1
                t_abs.append(t0)
            # Heights
            for j, t in enumerate(params_t):
                if ((j==0 and i==0) or t_abs[j] != 0) and j < part.n_splits:
                    x, y = wall.lerp(t_abs[j])
                    z = part.z[j]
                    pos = tM * Vector((x, y, z))
                    manip = Manipulator(context, "Z", "TOP", size, size, part, "z", min=0.01, max=50, step_round=2, sensitive=1, label="height", index=j, pos=pos)
                    self.manipulators.append(manip)
                    pos = tM * Vector((x, y, 0.5+z))
                    manip = Manipulator(context, "X", "TOP", size, size, part, "t", min=0.0, max=1, step_round=4, sensitive=0.1, label="position", index=j, pos=pos, normal=normal, label_factor=g.walls[i].length)
                    self.manipulators.append(manip)
            x, y = wall.lerp(0.5)
            z = wall.get_z(0.5)
            pos = tM * Vector((x, y,1.0+z))
            manip = Manipulator(context, "X", "TOP", size, size, part, "splits", min=1, max=31, step_round=0, sensitive=1, label="split", pos=pos, normal=normal)
            self.manipulators.append(manip)
            
            # Length and radius
            x, y = wall.lerp(1)
            pos = tM * Vector((x, y, 1.0))
            """
            child = context.scene.objects.get(c[1])
            x = child.data.WindowProperty[0].x
            dt = 0.5*x/wall.length
            c[3]-dt
            c[3]+dt
            """
            """
            o = wall.offset(1).gl_pts(context, tM)
            n0 = wall.sized_normal(0, 1.1).gl_pts(context, tM)
            n1 = wall.sized_normal(0, 1.1).gl_pts(context, tM)
            """
            if part.type == 'S_WALL':
                manip = Manipulator(context, "X", "TOP", size, size, part, "length", min=0.0, max=100, step_round=2, sensitive=1, label="length", pos=pos, normal=normal)
            else:
                manip = Manipulator(context, "X", "TOP", size, size, part, "radius", min=0.0, max=100, step_round=2, sensitive=0.5, label="radius", pos=pos)
                self.manipulators.append(manip)
                x, y = wall.lerp(0.5)
                pos = tM * Vector((x, y, 0.0))
                manip = Manipulator(context, "X", "TOP", size, size, part, "da", min=-180, max=180, step_round=1, sensitive=10, label="angle", pos=pos, value_factor=180/pi)
            # start angle
            self.manipulators.append(manip)
            x, y = wall.lerp(0)
            pos = tM * Vector((x, y, 0.0))
            manip = Manipulator(context, "X", "TOP", size, size, part, "a0", min=-180, max=180, step_round=1, sensitive=10, label="start angle", pos=pos, value_factor=180/pi)
            self.manipulators.append(manip)
            
    def modal(self, context, event):
        context.area.tag_redraw()
        if self.o != context.active_object:
            return {'FINISHED'}
        if event.type == 'RIGHTMOUSE' and event.value == 'PRESS':
            for manip in self.manipulators:
                manip.exit()
            return {'FINISHED'}
        for manip in self.manipulators:
            if manip.modal(context, event):
                return {'RUNNING_MODAL'}
        active = False
        for manip in self.manipulators:
            if manip.active or manip.hover:
                active = True
        if active == False:
            self.update(context)
        return {'PASS_THROUGH'} 
        
    def invoke(self, context, event):
        if context.space_data.type == 'VIEW_3D':
            self.manipulators = []
            self.relocate = []
            self.o = context.active_object
            self.update(context)
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "Active space must be a View3d")
            return {'CANCELLED'} 
 
bpy.utils.register_class(Wall2PartProperty)
bpy.utils.register_class(Wall2Property)
Mesh.Wall2Property = CollectionProperty(type=Wall2Property)
bpy.utils.register_class(ARCHIPACK_PT_wall2)
bpy.utils.register_class(ARCHIPACK_OT_wall2)
bpy.utils.register_class(ARCHIPACK_OT_wall2_insert)
bpy.utils.register_class(ARCHIPACK_OT_wall2_remove)
bpy.utils.register_class(ARCHIPACK_OT_wall2_manipulate)

