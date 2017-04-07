
import bpy, bgl, blf
from mathutils import Vector
from mathutils.geometry import intersect_line_plane
from bpy_extras import view3d_utils

def position_2d_from_coord(context, coord):
    """ coord given in local input coordsys
    """
    scene = context.scene
    region = context.region
    rv3d = context.region_data
    loc = view3d_utils.location_3d_to_region_2d(region, rv3d, coord)
    return loc

class Manipulator():
    """
        pick_tools actions
    """
    def __init__(self, context, axis, location, size_x, size_y, obj, attr, min=0, max=1, step_round=2, sensitive=1, value_factor=1, label="value", label_factor=1, index=None, pos=None, normal=None):
        """
            size_x, size_y on screen size in pixels
            location in 3d
        """
        self.active_color = (1.0, 0.0, 0.0, 1.0)
        self.hover_color  = (1.0, 1.0, 0.0, 1.0)
        self.basis_color  = (1.0, 1.0, 1.0, 1.0)
        self.manipulator_size = (0.5 * size_x, 0.5 * size_y)
        self.startPoint = (0,0)
        self.endPoint = (0,0)
        self.active = False
        self.o = context.active_object
        self.tM =  self.o.matrix_world
        self.itM = self.tM.inverted()
        self.axis = axis
        self.location = location
        self.obj = obj
        self.pos = pos
        self.attr = attr
        self.index = index
        self.label = label
        self.delta = 0
        self.max = max
        self.min = min
        self.sensitive = sensitive
        self.step_round = step_round
        self.label_factor = label_factor
        self.value_factor = value_factor
        if normal is None:
            if axis == 'X':
                self.normal = Vector((0,1,0))
            elif axis == 'Y':
                self.normal = Vector((1,0,0))
            else:
                self.normal = Vector((1,0,0))
        else:
            self.normal = normal
        self.get_value()
        self.draw_location(context)
        # Post selection actions
        args = (self, context)
        self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_callback, args, 'WINDOW', 'POST_PIXEL')
    
    def get_value(self):
        if self.index is not None:
            self.value = getattr(self.obj, self.attr)[self.index]*self.value_factor
        else:    
            self.value = getattr(self.obj, self.attr)*self.value_factor
            
    def set_value(self, value):
        if self.index is not None:
            getattr(self.obj, self.attr)[self.index] = value / self.value_factor
        else:
            setattr(self.obj, self.attr, value / self.value_factor)
            
    def exit(self):
        bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')
    
    def _position_3d_from_coord(self, context, coord):
        """return point in local input coordsys
        """
        scene = context.scene
        region = context.region
        rv3d = context.region_data
        view_vector_mouse = view3d_utils.region_2d_to_vector_3d(region, rv3d, coord)
        ray_origin_mouse = view3d_utils.region_2d_to_origin_3d(region, rv3d, coord)
        loc = intersect_line_plane(ray_origin_mouse, ray_origin_mouse+view_vector_mouse, self.origin, self.normal, False)
        return self.itM * loc
    
    def _contains(self, coord):
        x, y = coord
        lx, ly = self.manipulator_coords
        sx, sy = self.manipulator_size
        return (
            x > lx - sx and x < lx + sx and 
            y > ly - sy and y < ly + sy)
    
    @property
    def hover(self):
        return self._contains(self.endPoint)
    
    def object_center(self):
        self.tM = self.o.matrix_world
        self.itM = self.tM.inverted()
        x, y, z = self.o.bound_box[0]
        min = self.tM * Vector((x, y, z))
        x, y, z = self.o.bound_box[6]
        max = self.tM * Vector((x, y, z))
        return min + 0.5 * (max - min)
    
    def draw_location(self, context):
        if self.pos is not None:
            self.origin = self.pos
        else:    
            self.origin = self.object_center()
            if self.location == 'TOP':
                x, y, z = 0.5 * self.o.dimensions
            elif self.location == 'BOTTOM':
                x, y, z = -0.5 * self.o.dimensions
            else:
                size = 0
            if self.axis == 'X':
                self.origin += self.tM.to_3x3() * Vector((x, 0, 0))
            elif self.axis == 'Y':
                self.origin += self.tM.to_3x3() * Vector((0, y, 0))
            else:
                self.origin += Vector((0, 0, z))
        """
            origin is 3d location of handle
        """
        self.manipulator_coords = position_2d_from_coord(context, self.origin)  
        
    def complete(self, context):
        print("Manipulator.complete()")
        x, y, z = self._position_3d_from_coord(context, self.startPoint) - self._position_3d_from_coord(context, self.endPoint)
        if self.axis == 'X':
            v = round(x*self.sensitive, self.step_round)
        elif self.axis == 'Y':
            v = round(y*self.sensitive, self.step_round)
        else:
            v = round(z*self.sensitive, self.step_round)
        self.delta = min(self.max, max(self.min, self.value - v)) - self.value
        self.set_value(self.value + self.delta)
            
    def modal(self, context, event):
        if event.type == 'LEFTMOUSE' and event.value == 'PRESS':
            if self.hover:
                self.get_value()
                self.startPoint = (event.mouse_region_x, event.mouse_region_y)
                self.endPoint = (event.mouse_region_x, event.mouse_region_y)
                self.active = True
                return True
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            if self.active:
                self.endPoint = (event.mouse_region_x, event.mouse_region_y)
                self.draw_location(context)
                self.complete(context)
                self.active = False
        if event.type == 'MOUSEMOVE':
            self.endPoint = (event.mouse_region_x, event.mouse_region_y)
            self.draw_location(context)
            if self.active:
                self.complete(context)
                return True
        return False
    
    def _end(self):
        bgl.glEnd()
        bgl.glPopAttrib()
        bgl.glLineWidth(1)
        bgl.glDisable(bgl.GL_BLEND)
        bgl.glColor4f(0.0, 0.0, 0.0, 1.0)
    
    def string(self, text, cx, cy, size, colour):
        ''' my_string : the text we want to print
            pos_x, pos_y : coordinates in integer values
            size : font height.
            colour : used for definining the colour'''
        dpi, font_id = 72, 0 # dirty fast assignment
        bgl.glColor4f(*colour)
        blf.position(font_id, cx, cy, 0)
        blf.size(font_id, size, dpi)
        blf.draw(font_id, text)
        
    def triangle(self, cx, cy, pts, colour, width=1, style=bgl.GL_LINE):
        bgl.glPushAttrib(bgl.GL_ENABLE_BIT)
        bgl.glEnable(bgl.GL_BLEND)
        bgl.glColor4f(*colour)    
        bgl.glBegin(bgl.GL_POLYGON)
        for pt in pts:
            x, y = pt
            bgl.glVertex2f(cx+x, cy+y) 
        self._end()
        return
        self._start_line(colour, width=width, style=style) 
        for pt in pts:
            x, y = pt
            bgl.glVertex2f(cx+x, cy+y)
        x, y = pts[0]    
        bgl.glVertex2f(cx+x, cy+y)
        self._end()
        
    def draw_callback(self, _self, context):
        """
            draw on screen feedback using gl.
            
        """
        if self.active:
            color = self.active_color
        elif self.hover:
            color = self.hover_color
        else:
            color = self.basis_color
        cx, cy = self.manipulator_coords
        x, y = self.manipulator_size
        if self.axis != 'Z':
            if self.location == 'TOP':
                pts = [(-x, y), (x, 0), (-x, -y)]
            else:    
                pts = [(x, y), (-x, 0), (x, -y)]
        else:
            if self.location == 'TOP':
                pts = [(-x, -y), (0, y), (x, -y)]
            else:    
                pts = [(-x, y), (0, -y), (x, y)]
        self.triangle(cx, cy, pts, color)
        self.string(self.label+":"+str(round(self.label_factor*(self.value+self.delta),2)), x+cx, cy, 10, color)
  