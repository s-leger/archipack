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
from bpy.props import FloatProperty, IntProperty, CollectionProperty, EnumProperty 
from mathutils import Vector, Matrix
from .bmesh_utils import BmeshEdit as bmed
from .panel import Panel as DoorPanel
from .materialutils import MaterialUtils
from .archipack_handle import create_handle, door_handle_horizontal_01

class DoorPanelProperty(PropertyGroup):
    x = FloatProperty(
            name='width',
            min=0.25, max=10000,
            default=100.0, precision=2,
            description='Width'
            )
    y = FloatProperty(
            name='Depth',
            min=0.001, max=10000,
            default=0.02, precision=2,
            description='depth'
            )    
    z = FloatProperty(
            name='height',
            min=0.1, max=10000,
            default=2.0, precision=2,
            description='height'
            )
    direction = IntProperty(
            name="Direction", 
            min=0,
            max=1,
            description="open direction"
            )
    model = IntProperty(
            name="model",
            min=0,
            max=3,
            default=0,
            description="Model"
            )
    chanfer =FloatProperty(
            name='chanfer',
            min=0.001, max=10000,
            default=0.005, precision=3,
            description='chanfer'
            )
    panel_spacing =FloatProperty(
            name='spacing',
            min=0.001, max=10000,
            default=0.1, precision=2,
            description='distance between panels'
            )
    panel_bottom =FloatProperty(
            name='bottom',
            min=0.0, max=10000,
            default=0.0, precision=2,
            description='distance from bottom'
            )
    panel_border =FloatProperty(
            name='border',
            min=0.001, max=10000,
            default=0.2, precision=2,
            description='distance from border'
            )
    panels_x= IntProperty(
            name="panels h",
            min=1,
            max=50,
            default=1,
            description="panels h"
            )
    panels_y= IntProperty(
            name="panels v",
            min=1,
            max=50,
            default=1,
            description="panels v"
            )
    panels_distrib=EnumProperty(
            name='distribution',
            items=(
                ('REGULAR', 'Regular','',0),
                ('ONE_THIRD', '1/3 2/3', '', 1)
                ),
            default = 'REGULAR'
            )
    handle=EnumProperty(
            name='Shape',
            items=(
                ('NONE','No handle','',0),
                ('BOTH','Inside and outside','',1)
                ),
            default='BOTH'
            )
    
    @property
    def panels(self):        
        
        # subdivide side to weld panels
        subdiv_x = self.panels_x-1
        
        if self.panels_distrib == 'REGULAR':
            subdiv_y = self.panels_y-1
        else:
            subdiv_y = 2
            
        #  __ y0
        # |__ y1
        # x0 x1
        y0 = -self.y
        y1 = 0
        x0 = 0
        x1 = max(0.001, self.panel_border-0.5*self.panel_spacing)

        side = DoorPanel(
            False,               # profil closed
            [1,0,0,1],           # x index
            [x0, x1], 
            [y0, y0, y1, y1],
            [0,1,1,1],           # material index
            closed_path=True,    #
            subdiv_x=subdiv_x,
            subdiv_y=subdiv_y
            )
        
        face = None
        back = None
           
        if self.model == 1:
            #     /   y2-y3
            #  __/    y1-y0
            #   x2 x3
            x2 = 0.5*self.panel_spacing
            x3 = x2+self.chanfer
            y2 = y1+self.chanfer
            y3 = y0-self.chanfer

            face = DoorPanel(
                False,              # profil closed
                [0,1,2],            # x index
                [0, x2, x3], 
                [y1, y1, y2],
                [1,1,1],             # material index
                side_cap_front=2,    # cap index
                closed_path=True
                )
                
            back = DoorPanel(
                False,              # profil closed
                [0,1,2],            # x index
                [x3, x2, 0], 
                [y3, y0, y0],
                [0,0,0],             # material index
                side_cap_back=0,     # cap index
                closed_path=True
                )
                
        elif self.model == 2:     
            #               /   y2-y3
            #  ___    _____/    y1-y0
            #     \  /
            #      \/           y4-y5
            # 0 x2 x4 x5 x6 x3
            x2 = 0.5*self.panel_spacing
            x4 = x2+self.chanfer
            x5 = x4+self.chanfer
            x6 = x5+4*self.chanfer
            x3 = x6+self.chanfer
            y2 = y1-self.chanfer
            y4 = y1+self.chanfer
            y3 = y0+self.chanfer
            y5 = y0-self.chanfer
            face = DoorPanel(
                False,                    # profil closed
                [0,1,2,3,4,5],            # x index
                [0,  x2, x4, x5, x6, x3], 
                [y1, y1, y4, y1, y1, y2],
                [1,1,1,1,1,1],            # material index
                side_cap_front=5,          # cap index
                closed_path=True
                )
                
            back = DoorPanel(
                False,                    # profil closed
                [0,1,2,3,4,5],            # x index
                [x3, x6, x5, x4, x2, 0], 
                [y3, y0, y0, y5, y0, y0],
                [0,0,0,0,0,0],             # material index
                side_cap_back=0,          # cap index
                closed_path=True
                ) 
        
        elif self.model == 3:         
            #      _____      y2-y3
            #     /     \     y4-y5
            #  __/            y1-y0
            # 0 x2 x3 x4 x5 
            x2 = 0.5*self.panel_spacing
            x3 = x2+self.chanfer
            x4 = x3+4*self.chanfer
            x5 = x4+2*self.chanfer
            y2 = y1-self.chanfer
            y3 = y0+self.chanfer
            y4 = y2+self.chanfer
            y5 = y3-self.chanfer
            face = DoorPanel(
                False,              # profil closed
                [0,1,2,3,4],            # x index
                [0, x2, x3, x4, x5], 
                [y1, y1, y2, y2, y4],
                [1,1,1,1,1],             # material index
                side_cap_front=4,    # cap index
                closed_path=True
                )
                
            back = DoorPanel(
                False,              # profil closed
                [0,1,2,3,4],            # x index
                [x5, x4, x3, x2, 0], 
                [y5, y3, y3, y0, y0],
                [0,0,0,0,0],             # material index
                side_cap_back=0,     # cap index
                closed_path=True
                ) 
        
        else:
            side.side_cap_front = 3
            side.side_cap_back = 0
           
        return side, face, back
            
    @property
    def verts(self):
        
        subdiv_x = self.panels_x-1
        if self.panels_distrib == 'REGULAR':
            subdiv_y = self.panels_y-1
        else:
            subdiv_y = 2
            
        radius = Vector((0.8,0.5,0))
        center = Vector((0,self.z-radius.x,0))
        
        if self.direction == 0:
            pivot = 1
        else:
            pivot = -1
        
        path_type = 'RECTANGLE'
        curve_steps = 16
        side, face, back = self.panels
        
        x1 = max(0.001, self.panel_border-0.5*self.panel_spacing)
        bottom_z = self.panel_bottom
        shape_z = [0,bottom_z,bottom_z,0]
        origin = Vector((-pivot*0.5*self.x, 0 ,0))
        offset = Vector((0, 0 ,0))
        size = Vector((self.x, self.z, 0))
        verts = side.vertices(curve_steps, offset, center, origin, size, radius, 0, pivot, shape_z=shape_z, path_type=path_type)
        if face is not None:
            p_radius = radius.copy()
            p_radius.x -= x1
            p_radius.y -= x1
            if self.panels_distrib == 'REGULAR':
                p_size = Vector(((self.x-2*x1)/self.panels_x, (self.z-2*x1-bottom_z)/self.panels_y, 0))
                for i in range(self.panels_x):
                    for j in range(self.panels_y):
                        if j < subdiv_y:
                            shape = 'RECTANGLE'
                        else:
                            shape = path_type
                        offset = Vector(((pivot*0.5*self.x)+p_size.x*(i+0.5)-0.5*size.x+x1, bottom_z+p_size.y*j+x1,0))
                        origin = Vector((p_size.x*(i+0.5)-0.5*size.x+x1, bottom_z+p_size.y*j+x1,0))
                        verts += face.vertices(curve_steps, offset, center, origin, p_size, p_radius, 0, 0, shape_z=None, path_type=shape)
                        if back is not None:
                            verts += back.vertices(curve_steps, offset, center, origin, p_size, p_radius, 0, 0, shape_z=None, path_type=shape)
            else:
                ####################################
                # Ratio vertical panels 1/3 - 2/3
                ####################################
                p_size = Vector(((self.x-2*x1)/self.panels_x, (self.z-2*x1-bottom_z)/3, 0))
                p_size_2x = Vector((p_size.x, p_size.y*2,0))
                for i in range(self.panels_x):
                    j = 0
                    offset = Vector(((pivot*0.5*self.x)+p_size.x*(i+0.5)-0.5*size.x+x1, bottom_z+p_size.y*j+x1,0))
                    origin = Vector((p_size.x*(i+0.5)-0.5*size.x+x1, bottom_z+p_size.y*j+x1,0))
                    shape = 'RECTANGLE'
                    face.subdiv_y = 0
                    verts += face.vertices(curve_steps, offset, center, origin, p_size, p_radius, 0, 0, shape_z=None, path_type=shape)
                    if back is not None:
                        back.subdiv_y = 0
                        verts += back.vertices(curve_steps, offset, center, origin, p_size, p_radius, 0, 0, shape_z=None, path_type=shape)
                    j = 1
                    offset = Vector(((pivot*0.5*self.x)+p_size.x*(i+0.5)-0.5*size.x+x1, bottom_z+p_size.y*j+x1,0))
                    origin = Vector((p_size.x*(i+0.5)-0.5*size.x+x1, bottom_z+p_size.y*j+x1,0))
                    shape = path_type
                    face.subdiv_y = 1
                    verts += face.vertices(curve_steps, offset, center, origin, p_size_2x, p_radius, 0, 0, shape_z=None, path_type=path_type)
                    if back is not None:
                        back.subdiv_y = 1
                        verts += back.vertices(curve_steps, offset, center, origin, p_size_2x, p_radius, 0, 0, shape_z=None, path_type=path_type)
                    
        return verts    
        
    @property
    def faces(self): 
        
        subdiv_x = self.panels_x-1
        if self.panels_distrib == 'REGULAR':
            subdiv_y = self.panels_y-1
        else:
            subdiv_y = 2
            
        path_type = 'RECTANGLE'
        curve_steps = 16
        side, face, back = self.panels
        
        faces = side.faces(curve_steps, path_type=path_type)
        faces_offset = side.n_verts(curve_steps, path_type=path_type)
        
        if face is not None:
            if self.panels_distrib == 'REGULAR':
                for i in range(self.panels_x):
                    for j in range(self.panels_y):
                        if j < subdiv_y:
                            shape = 'RECTANGLE'
                        else:
                            shape = path_type
                        faces += face.faces(curve_steps, path_type=shape, offset=faces_offset)   
                        faces_offset += face.n_verts(curve_steps, path_type=shape)
                        if back is not None:
                            faces += back.faces(curve_steps, path_type=shape, offset=faces_offset)   
                            faces_offset += back.n_verts(curve_steps, path_type=shape)
            else:
                ####################################
                # Ratio vertical panels 1/3 - 2/3
                ####################################
                for i in range(self.panels_x):
                    j = 0
                    shape = 'RECTANGLE'
                    face.subdiv_y = 0
                    faces += face.faces(curve_steps, path_type=shape, offset=faces_offset)
                    faces_offset += face.n_verts(curve_steps, path_type=shape)
                    if back is not None:
                        back.subdiv_y = 0
                        faces += back.faces(curve_steps, path_type=shape, offset=faces_offset)
                        faces_offset += back.n_verts(curve_steps, path_type=shape)
                    j = 1
                    shape = path_type
                    face.subdiv_y = 1
                    faces += face.faces(curve_steps, path_type=path_type, offset=faces_offset)
                    faces_offset += face.n_verts(curve_steps, path_type=path_type)
                    if back is not None:
                        back.subdiv_y = 1
                        faces += back.faces(curve_steps, path_type=path_type, offset=faces_offset)
                        faces_offset += back.n_verts(curve_steps, path_type=path_type)
        
        print("faces %s" % (faces) )
        return faces
        
    @property
    def uvs(self):
        
        subdiv_x = self.panels_x-1
        if self.panels_distrib == 'REGULAR':
            subdiv_y = self.panels_y-1
        else:
            subdiv_y = 2
        
        radius = Vector((0.8,0.5,0))
        center = Vector((0,self.z-radius.x,0))
        
        if self.direction == 0:
            pivot = 1
        else:
            pivot = -1
            
        path_type = 'RECTANGLE'
        curve_steps = 16
        side, face, back = self.panels 
        
        x1 = max(0.001, self.panel_border-0.5*self.panel_spacing)
        bottom_z = self.panel_bottom
        origin = Vector((-pivot*0.5*self.x, 0 ,0))
        size = Vector((self.x, self.z, 0))
        uvs = side.uv(curve_steps, center, origin, size, radius, 0, pivot, 0, self.panel_border, path_type=path_type)
        if face is not None:
            p_radius = radius.copy()
            p_radius.x -= x1
            p_radius.y -= x1
            if self.panels_distrib == 'REGULAR':
                p_size = Vector(((self.x-2*x1)/self.panels_x, (self.z-2*x1-bottom_z)/self.panels_y, 0))
                for i in range(self.panels_x):
                    for j in range(self.panels_y):
                        if j < subdiv_y:
                            shape = 'RECTANGLE'
                        else:
                            shape = path_type
                        origin = Vector((p_size.x*(i+0.5)-0.5*size.x+x1, bottom_z+p_size.y*j+x1,0))
                        uvs += face.uv(curve_steps, center, origin, p_size, p_radius, 0, 0, 0, 0, path_type=shape)
                        if back is not None:
                            uvs += back.uv(curve_steps, center, origin, p_size, p_radius, 0, 0, 0, 0, path_type=shape)
            else:
                ####################################
                # Ratio vertical panels 1/3 - 2/3
                ####################################
                p_size = Vector(((self.x-2*x1)/self.panels_x, (self.z-2*x1-bottom_z)/3, 0))
                p_size_2x = Vector((p_size.x, p_size.y*2,0))
                for i in range(self.panels_x):
                    j = 0
                    origin = Vector((p_size.x*(i+0.5)-0.5*size.x+x1, bottom_z+p_size.y*j+x1,0))
                    shape = 'RECTANGLE'
                    face.subdiv_y = 0
                    uvs += face.uv(curve_steps, center, origin, p_size, p_radius, 0, 0, 0, 0, path_type=shape)
                    if back is not None:
                        back.subdiv_y = 0
                        uvs += back.uv(curve_steps, center, origin, p_size, p_radius, 0, 0, 0, 0, path_type=shape)
                    j = 1
                    origin = Vector((p_size.x*(i+0.5)-0.5*size.x+x1, bottom_z+p_size.y*j+x1,0))
                    shape = path_type
                    face.subdiv_y = 1
                    uvs += face.uv(curve_steps, center, origin, p_size_2x, p_radius, 0, 0, 0, 0, path_type=path_type)
                    if back is not None:
                        back.subdiv_y = 1
                        uvs += back.uv(curve_steps, center, origin, p_size_2x, p_radius, 0, 0, 0, 0, path_type=path_type)
        
        
        return uvs

    @property
    def matids(self): 
        if self.panels_distrib == 'REGULAR':
            subdiv_y = self.panels_y-1
        else:
            subdiv_y = 2
        
        path_type = 'RECTANGLE'
        curve_steps = 16
        side, face, back = self.panels
        
        mat = side.mat(curve_steps, 1, 0, path_type=path_type)
        
        if face is not None:
            if self.panels_distrib == 'REGULAR':
                for i in range(self.panels_x):
                    for j in range(self.panels_y):
                        if j < subdiv_y:
                            shape = 'RECTANGLE'
                        else:
                            shape = path_type
                        mat += face.mat(curve_steps, 1, 1, path_type=shape)
                        if back is not None:
                            mat += back.mat(curve_steps, 0, 0, path_type=shape)
            else:
                ####################################
                # Ratio vertical panels 1/3 - 2/3
                ####################################
                for i in range(self.panels_x):
                    j = 0
                    shape = 'RECTANGLE'
                    face.subdiv_y = 0
                    mat += face.mat(curve_steps, 1, 1, path_type=shape)
                    if back is not None:
                        back.subdiv_y = 0
                        mat += back.mat(curve_steps, 0, 0, path_type=shape)
                    j = 1
                    shape = path_type
                    face.subdiv_y = 1
                    mat += face.mat(curve_steps, 1, 1, path_type=shape)
                    if back is not None:
                        back.subdiv_y = 1
                        mat += back.mat(curve_steps, 0, 0, path_type=shape)
        
        print("mat %s" % (mat) )
        return mat
	
    def find_in_selection(self, context):
        """
            find witch selected object this instance belongs to
            provide support for "copy to selected"
        """
        active = context.active_object
        selected = [o for o in context.selected_objects]
        for o in selected:
            c, props = ARCHIPACK_PT_door_panel.params(o)
            if props == self:
                return active, selected, o
        return active, selected, None    
    
    def find_handle(self, o):
        for child in o.children:
            if 'archipack_handle' in child:
                return child
        return None
     
    def update_handle(self, context, o):
        handle = self.find_handle(o)
        if handle is None:
            m = bpy.data.meshes.new("Handle")
            handle = create_handle(context, o, m)
            MaterialUtils.add_handle_materials(handle)
        verts, faces = door_handle_horizontal_01(self.direction, 1)
        b_verts, b_faces = door_handle_horizontal_01(self.direction, 0, offset=len(verts)) 
        b_verts = [(v[0], v[1]-self.y, v[2]) for v in b_verts]
        handle_y = 0.07
        handle.location = ((1-self.direction*2)*(self.x-handle_y), 0, 0.5*self.z)
        bmed.buildmesh(context, handle, verts+b_verts, faces+b_faces)
            
    def remove_handle(self, context, o):
        handle = self.find_handle(o)
        if handle is not None:
            context.scene.objects.unlink(handle)
            bpy.data.objects.remove(handle, do_unlink=True)
     
    def update(self, context):
        active, selected, o = self.find_in_selection(context)
        if o is None:
            return
            
        bmed.buildmesh(context, o, self.verts, self.faces, matids=self.matids, uvs=self.uvs, weld=True)
        
        if self.handle == 'NONE':
            self.remove_handle(context, o)
        else:
            self.update_handle(context, o)
         
        active.select = True
        context.scene.objects.active = active
        
class ARCHIPACK_PT_door_panel(Panel):
    bl_idname = "ARCHIPACK_PT_door_panel"
    bl_label = "Door"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'ArchiPack'
    def draw(self, context):
        """
        o = context.object
        o, prop = ARCHIPACK_PT_door_panel.params(o)
        if prop is None:
            return 
        """
        layout = self.layout
        """
        layout.prop(prop, 'direction')
        layout.prop(prop, 'model')
        layout.prop(prop, 'chanfer')
        layout.prop(prop, 'panel_bottom')
        layout.prop(prop, 'panel_spacing')
        layout.prop(prop, 'panel_border')
        layout.prop(prop, 'panels_x')
        layout.prop(prop, 'panels_y')
        layout.prop(prop, 'panels_distrib')
        layout.label(text="Show properties")
        """
        layout.operator("archipack.select_parent")
        return
        
    @classmethod
    def params(cls, o):
        if cls.filter(o):
            if 'archipack_doorpanel' in o.data:
                return o, o.data.archipack_doorpanel[0]
            else:
                for child in o.children:
                    o, props = cls.params(child)
                    if props is not None:
                        return o, props
        return o, None
    @classmethod
    def filter(cls, o):
        try:
            return bool('archipack_doorpanel' in o.data)
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
class ARCHIPACK_OT_door_panel(Operator):
    bl_idname = "archipack.door_panel"
    bl_label = "Door model 1"
    bl_description = "Door model 1"
    bl_category = 'Archipack'
    bl_options = {'REGISTER', 'UNDO'}
    x = FloatProperty(
            name='width',
            min=0.1, max=10000,
            default=0.80, precision=2,
            description='Width'
            )
    z = FloatProperty(
            name='height',
            min=0.1, max=10000,
            default=2.0, precision=2,
            description='height'
            )
    y = FloatProperty(
            name='depth',
            min=0.001, max=10000,
            default=0.02, precision=2,
            description='Depth'
            )
    direction = IntProperty(
            name="direction", 
            min=0,
            max=1,
            description="open direction"
            )
    model = IntProperty(
            name="model", 
            min=0,
            max=3,
            description="panel type"
            )
    chanfer =FloatProperty(
            name='chanfer',
            min=0.001, max=10000,
            default=0.005, precision=3,
            description='chanfer'
            )
    panel_spacing =FloatProperty(
            name='spacing',
            min=0.001, max=10000,
            default=0.1, precision=2,
            description='distance between panels'
            )
    panel_bottom =FloatProperty(
            name='bottom',
            min=0.0, max=10000,
            default=0.0, precision=2,
            description='distance from bottom'
            )
    panel_border =FloatProperty(
            name='border',
            min=0.001, max=10000,
            default=0.2, precision=2,
            description='distance from border'
            )
    panels_x= IntProperty(
            name="panels h",
            min=1,
            max=50,
            default=1,
            description="panels h"
            )
    panels_y= IntProperty(
            name="panels v",
            min=1,
            max=50,
            default=1,
            description="panels v"
            )
    panels_distrib=EnumProperty(
            name='distribution',
            items=(
                ('REGULAR', 'Regular','',0),
                ('ONE_THIRD', '1/3 2/3', '', 1)
                ),
            default = 'REGULAR'
            )
    handle=EnumProperty(
            name='Shape',
            items=(
                ('NONE','No handle','',0),
                ('BOTH','Inside and outside','',1)
                ),
            default='BOTH'
            )
    # -----------------------------------------------------
    # Draw (create UI interface)
    # -----------------------------------------------------
    # noinspection PyUnusedLocal
    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label("Use Properties panel (N) to define parms", icon='INFO')
    
    def create(self, context):
        """
            expose only basic params in operator
            use object property for other params
        """
        m = bpy.data.meshes.new("Door")
        o = bpy.data.objects.new("Door", m)
        d = m.archipack_doorpanel.add()
        d.x = self.x
        d.y = self.y
        d.z = self.z
        d.model = self.model
        d.direction = self.direction
        d.chanfer = self.chanfer
        d.panel_border = self.panel_border
        d.panel_bottom = self.panel_bottom
        d.panel_spacing = self.panel_spacing
        d.panels_distrib = self.panels_distrib
        d.panels_x = self.panels_x
        d.panels_y = self.panels_y 
        d.handle = self.handle 
        context.scene.objects.link(o)
        o.lock_location[0] = True
        o.lock_location[1] = True
        o.lock_location[2] = True
        o.lock_rotation[0] = True
        o.lock_rotation[1] = True
        o.lock_scale[0] = True
        o.lock_scale[1] = True
        o.lock_scale[2] = True
        o.select = True
        context.scene.objects.active = o 
        d.update(context)
        MaterialUtils.add_door_materials(o)
        o.select = True
        context.scene.objects.active = o 
        o.lock_rotation[0] = True
        o.lock_rotation[1] = True
        return o
    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            bpy.ops.object.select_all(action="DESELECT")
            o = self.create(context)
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}

class ARCHIPACK_OT_select_parent(Operator):
    bl_idname = "archipack.select_parent"
    bl_label = "Edit parameters"
    bl_description = "Edit parameters located on parent"
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
    
    # -----------------------------------------------------
    # Execute
    # -----------------------------------------------------
    def execute(self, context):
        if context.mode == "OBJECT":
            if context.active_object is not None and context.active_object.parent is not None:
                bpy.ops.object.select_all(action="DESELECT")
                context.active_object.parent.select = True
                context.scene.objects.active = context.active_object.parent
            return {'FINISHED'}
        else:
            self.report({'WARNING'}, "Archipack: Option only valid in Object mode")
            return {'CANCELLED'}
            
            
bpy.utils.register_class(DoorPanelProperty)
Mesh.archipack_doorpanel = CollectionProperty(type=DoorPanelProperty)
bpy.utils.register_class(ARCHIPACK_PT_door_panel)
bpy.utils.register_class(ARCHIPACK_OT_door_panel)            
bpy.utils.register_class(ARCHIPACK_OT_select_parent) 
