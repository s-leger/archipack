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
import bpy
import os
from bl_operators.presets import AddPresetBase
from mathutils import Vector
from bpy.props import StringProperty
from .archipack_gl import (
    ThumbHandle, Screen, GlRect
)

preset_paths = bpy.utils.script_paths("presets")
addons_paths = bpy.utils.script_paths("addons")


class PresetMenuItem():
    def __init__(self, thumbsize, image, preset):
        name = bpy.path.display_name_from_filepath(preset)
        self.preset = preset
        self.handle = ThumbHandle(thumbsize, image, name, draggable=True)
        self.enable = True

    def set_pos(self, context, pos):
        self.handle.set_pos(context, pos)

    def check_hover(self, mouse_pos):
        self.handle.check_hover(mouse_pos)

    def mouse_press(self):
        if self.handle.hover:
            self.handle.hover = False
            self.handle.active = True
            return True
        return False

    def draw(self, context):
        if self.enable:
            self.handle.draw(context)


class PresetMenu():
    def __init__(self, context, category, thumbsize=Vector((150, 100))):
        self.imageList = []
        self.menuItems = []
        self.thumbsize = thumbsize
        file_list = self.scan_files(category)
        self.default_image = None
        self.load_default_image()
        for filepath in file_list:
            self.make_menuitem(filepath)
        self.margin = 50
        self.y_scroll = 0
        self.scroll_max = 1000
        self.spacing = Vector((25, 25))
        self.screen = Screen(self.margin)
        self.mouse_pos = Vector((0, 0))
        self.bg = GlRect(colour=(0, 0, 0, 0.5))
        self.set_pos(context)

    def load_default_image(self):
        img_idx = bpy.data.images.find("missing.png")
        if img_idx > -1:
            self.default_image = bpy.data.images[img_idx]
            self.imageList.append(self.default_image.filepath_raw)
            return
        sub_path = "archipack" + os.path.sep + "presets" + os.path.sep + "missing.png"
        for path in addons_paths:
            filepath = os.path.join(path, sub_path)
            if os.path.exists(filepath) and os.path.isfile(filepath):
                self.default_image = bpy.data.images.load(filepath=filepath)
                self.imageList.append(self.default_image.filepath_raw)
                break
        if self.default_image is None:
            raise EnvironmentError("archipack/presets/missing.png not found")

    def scan_files(self, category):
        file_list = []
        # load default presets
        sub_path = "archipack" + os.path.sep + "presets" + os.path.sep + category
        for path in addons_paths:
            presets_path = os.path.join(path, sub_path)
            if os.path.exists(presets_path):
                file_list += [presets_path + os.path.sep + f[:-3]
                    for f in os.listdir(presets_path)
                    if f.endswith('.py') and
                    not f.startswith('.')]
        # load user def presets
        for path in preset_paths:
            presets_path = os.path.join(path, category)
            if os.path.exists(presets_path):
                file_list += [presets_path + os.path.sep + f[:-3]
                    for f in os.listdir(presets_path)
                    if f.endswith('.py') and
                    not f.startswith('.')]

        file_list.sort()
        return file_list

    def clearImages(self):
        for image in bpy.data.images:
            if image.filepath_raw in self.imageList:
                # image.user_clear()
                bpy.data.images.remove(image, do_unlink=True)
        self.imageList.clear()

    def make_menuitem(self, filepath):
        image = None
        img_idx = bpy.data.images.find(os.path.basename(filepath) + '.png')
        if img_idx > -1:
            image = bpy.data.images[img_idx]
            self.imageList.append(image.filepath_raw)
        elif os.path.exists(filepath + '.png') and os.path.isfile(filepath + '.png'):
            image = bpy.data.images.load(filepath=filepath + '.png')
            self.imageList.append(image)
        if image is None:
            image = self.default_image
        item = PresetMenuItem(self.thumbsize, image, filepath + '.py')
        self.menuItems.append(item)

    def set_pos(self, context):

        x_min, x_max, y_min, y_max = self.screen.size(context)
        self.bg.set_pos([Vector((x_min, y_min)), Vector((x_max, y_max))])
        x_min += 0.5 * self.thumbsize.x + 0.5 * self.margin
        x_max -= 0.5 * self.thumbsize.x + 0.5 * self.margin
        y_max -= 0.5 * self.thumbsize.y + 0.5 * self.margin
        y_min += 0.5 * self.margin
        x = x_min
        y = y_max + self.y_scroll
        n_rows = 0
        for item in self.menuItems:
            if y > y_max or y < y_min:
                item.enable = False
            else:
                item.enable = True
            item.set_pos(context, Vector((x, y)))
            x += self.thumbsize.x + self.spacing.x
            if x > x_max:
                n_rows += 1
                x = x_min
                y -= self.thumbsize.y + self.spacing.y
        
        self.scroll_max = (n_rows - 1) * (self.thumbsize.y + self.spacing.y)
                
    def draw(self, context):
        self.bg.draw(context)
        for item in self.menuItems:
            item.draw(context)

    def mouse_press(self, context, event):
        self.mouse_position(event)
        for item in self.menuItems:
            if item.mouse_press():
                # load item preset
                return item.preset
        return None

    def mouse_position(self, event):
        self.mouse_pos.x, self.mouse_pos.y = event.mouse_region_x, event.mouse_region_y

    def mouse_move(self, context, event):
        self.mouse_position(event)
        for item in self.menuItems:
            item.check_hover(self.mouse_pos)
    
    def scroll_up(self, context, event):
        self.y_scroll = max(0, self.y_scroll - (self.thumbsize.y + self.spacing.y))
        self.set_pos(context)
        print("scroll_up %s" % (self.y_scroll))
        
    def scroll_down(self, context, event):
        self.y_scroll = min( self.scroll_max, self.y_scroll + (self.thumbsize.y + self.spacing.y))
        self.set_pos(context)
        print("scroll_down %s" % (self.y_scroll))
        
class PresetMenuOperator():

    preset_operator = StringProperty(
        options={'SKIP_SAVE'},
        default="script.execute_preset"
    )

    def __init__(self):
        self.menu = None
        self._handle = None

    def exit(self, context):
        self.menu.clearImages()
        bpy.types.SpaceView3D.draw_handler_remove(self._handle, 'WINDOW')

    def draw_handler(self, _self, context):
        self.menu.draw(context)

    def modal(self, context, event):
        context.area.tag_redraw()
        if event.type == 'MOUSEMOVE':
            self.menu.mouse_move(context, event)
        elif event.type == 'WHEELUPMOUSE':
            self.menu.scroll_up(context, event)
        elif event.type == 'WHEELDOWNMOUSE':
            self.menu.scroll_down(context, event)
        elif event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
            preset = self.menu.mouse_press(context, event)
            if preset is not None:
                self.exit(context)
                po = self.preset_operator.split(".")
                op = getattr(getattr(bpy.ops, po[0]), po[1])
                if self.preset_operator == 'script.execute_preset':
                    # call from preset menu
                    # ensure right active_object class
                    d = getattr(bpy.types, self.preset_subdir).datablock(context.active_object)
                    if d is not None:
                        d.auto_update = False
                        op(filepath=preset, menu_idname=self.bl_idname)
                        d.auto_update = True
                else:
                    # call draw operator
                    if op.poll():
                        op('INVOKE_DEFAULT', filepath=preset)
                    else:
                        print("Poll failed")

                return {'FINISHED'}

        elif event.type in {'RIGHTMOUSE', 'ESC'}:
            self.exit(context)
            return {'CANCELLED'}

        return {'RUNNING_MODAL'}

    def invoke(self, context, event):
        if context.area.type == 'VIEW_3D':
            self.menu = PresetMenu(context, self.preset_subdir)

            # the arguments we pass the the callback
            args = (self, context)
            # Add the region OpenGL drawing callback
            # draw in view space with 'POST_VIEW' and 'PRE_VIEW'
            self._handle = bpy.types.SpaceView3D.draw_handler_add(self.draw_handler, args, 'WINDOW', 'POST_PIXEL')
            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "View3D not found, cannot show preset flinger")
            return {'CANCELLED'}


class ArchipackPreset(AddPresetBase):

    @classmethod
    def poll(cls, context):
        o = context.active_object
        return o is not None and \
            o.data is not None and \
            "archipack_" + cls.__name__[13:-7] in o.data

    @property
    def preset_subdir(self):
        return "archipack_" + self.__class__.__name__[13:-7]

    @property
    def blacklist(self):
        """
            properties black list for presets
            may override on addon basis
        """
        return []

    @property
    def preset_values(self):
        blacklist = self.blacklist
        blacklist.extend(bpy.types.Mesh.bl_rna.properties.keys())
        d = getattr(bpy.context.active_object.data, self.preset_subdir)[0]
        props = d.rna_type.bl_rna.properties.items()
        ret = []
        for prop_id, prop in props:
            if prop_id not in blacklist:
                if not (prop.is_hidden or prop.is_skip_save):
                    ret.append("d.%s" % prop_id)
        return ret

    @property
    def preset_defines(self):
        return ["d = bpy.context.active_object.data." + self.preset_subdir + "[0]"]

    def pre_cb(self, context):
        return

    def remove(self, context, filepath):
        # remove preset
        os.remove(filepath)
        # remove thumb
        os.remove(filepath[:-3] + ".png")

    def post_cb(self, context):

        if not self.remove_active:

            name = self.name.strip()
            if not name:
                return

            filename = self.as_filename(name)
            target_path = os.path.join("presets", self.preset_subdir)
            target_path = bpy.utils.user_resource('SCRIPTS',
                                                  target_path,
                                                  create=True)

            filepath = os.path.join(target_path, filename) + ".png"

            # render thumb
            scene = context.scene
            render = scene.render
            render.resolution_x = 150
            render.resolution_y = 100
            render.resolution_percentage = 100
            render.filepath = filepath
            render.use_file_extension = True
            render.use_overwrite = True
            render.image_settings.file_format = 'PNG'
            render.image_settings.color_mode = 'RGBA'
            render.image_settings.color_depth = '8'
            render.use_compositing = False
            render.use_sequencer = False
            bpy.ops.render.render(animation=False, write_still=True, use_viewport=False)
            return
