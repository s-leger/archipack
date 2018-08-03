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
# Author: Stephen Leger (s-leger)
# ----------------------------------------------------------
from bpy.utils import previews
import bpy
import os


"""
from archipack.archipack_iconmanager import icons

im = icons['archipack_curves']['Molding_1']
p = bpy.data.images.new("Temp", im.image_size[0], im.image_size[1])
p.pixels = im.image_pixels_float[:]
p.scale(128, 128)

im.icon_size = [128, 128]
im.icon_pixels_float = p.pixels[:]
res = bpy.data.images.new("Temp", im.icon_size[0], im.icon_size[1])
res.pixels = im.icon_pixels_float[:]
"""


class IconsCollectionManager(dict):
    """
     Global icons management for archipack ui and presets
     Handle icons for ui
     Provide 
      - loading icon from folder 
      - reload single icon after change
      - remove icon / category
      - new icons loading after render
      - watch for icon changes on folder
     Icons for presets:
      - Icons for gl menuitems
      - Icons for dynamic enum for use in template_icon_view
    """
    def __init__(self):
        dict.__init__(self)
        # category files remove
        self._refresh = {}
        # category files add/change watcher
        self._watcher = {}
        # hold reference to icon enums
        self._enums = {}
        self._folders = {}
        self.allowed_ext = {
            '.png',
            '.jpg'
            }

    def add(self, category):
        """
         Add a category
        """
        if category not in self:
            self._enums[category] = []
            self._folders[category] = {}
            self[category] = previews.new()

    def remove(self, category):
        """
         Remove a category
        """
        if category in self:
            self._enums[category].clear()
            self._folders[category].clear()
            self[category].close()
            previews.remove(self[category])
            del self[category]
        if category in self._watcher:
            del self._watcher[category]
        if category in self._refresh:
            del self._refresh[category]

    def add_files(self, files, category):
        """
         Add files in list
         @folder: folder with images
         @files: array of (full file path, filename.ext)
         @category: category name
        """
        folders = self._folders[category]
        icons = self[category]
        for full_path, file in files:
            name, ext = os.path.splitext(file)
            if name not in icons and ext in self.allowed_ext:
                if name not in folders:
                    folders[name] = {
                        'label': bpy.path.display_name(name),
                        'file_name': name,
                        'full_path': full_path
                        }
                icons.load(
                    name,
                    full_path,
                    'IMAGE')

    def load(self, folder, category):
        """
         Load files from folder
         @folder: folder with images
         @category: category name
        """
        self.add(category)
        if os.path.exists(folder):
            files = [(os.path.join(folder, file), file) for file in os.listdir(folder)]
            print("load_presets", category, len(files), folder)
            self.add_files(files, category)
        else:
            print("load_presets folder not found:", category, folder)

    def load_presets(self, category):
        """
         Load thumbs from archipack's preset dir,
         and all presets/category sub folders
         User presets override factory one
        """
        folders = bpy.utils.script_paths("presets")
        # append after user script so user override factory
        if category not in self:
            # only load factory on startup
            # so user is able to remove and override
            folders.append(os.path.join(
                os.path.dirname(__file__),
                "presets"))

        for folder in folders:
            self.load(os.path.join(folder, category), category)

    def add_preset(self, category, name):
        """
          Add / Replace a preset thumb (from user path only)
          call this when a preset is saved
          rescan folder until the count of images grow
        """
        if category not in self._watcher:
            self._watcher[category] = {}

        # refresh enum until file is loaded
        self._watcher[category][name] = True

        if category not in self:
            self.load_presets(category)

        if name in self[category]:
            del self._folders[category][name]
            self[category].pop(name)

    def remove_preset(self, category, name):
        """
          Remove a preset thumb (from user path only)
          call this when a preset is saved so enum will refresh from files
        """
        if name in self[category]:
            del self._folders[category][name]
            del self._enums[category][name]
            self[category].pop(name)

        if category in self._watcher and name in self._watcher[category]:
            del self._watcher[category][name]

        # ensure enum rescan folder once on next run
        self._refresh[category] = True

    def cleanup(self):
        """
         Clear all icons in all category
        """
        self._enums.clear()
        self._folders.clear()
        self._watcher.clear()
        self._refresh.clear()
        for icons in self.values():
            previews.remove(icons)
        self.clear()

    def enum(self, context, category):
        """
         Build and return icon enumerator
        """

        # load missing category
        if category not in self:
            self.load_presets(category)

        # refresh enum on file delete
        if category in self._refresh:
            del self._refresh[category]

        # reload until watched images are found
        if category in self._watcher:
            kill_watcher = True
            for name in self._watcher[category]:
                if name not in self[category]:
                    kill_watcher = False
            if kill_watcher:
                for name in self._watcher[category]:
                    print("reload", category, name)
                    self[category][name].reload()
                del self._watcher[category]
            else:
                self.load_presets(category)

        self._enums[category].clear()
        folders = self._folders[category]
        items = list(self[category].items())
        items.sort(key=lambda i: i[0].lower())
        for i, item in enumerate(items):
            name, preview = item

            # use preview image
            """
            if not preview.is_icon_custom:
                preview.icon_size = preview.image_size[:]
                preview.icon_pixels_float = preview.image_pixels_float[:]
                preview.icon_pixels = preview.image_pixels[:]
            """

            self._enums[category].append((
                folders[name]['file_name'],
                folders[name]['label'],
                folders[name]['label'],
                preview.icon_id,
                i
                ))

        return self._enums[category]

    def as_blender_image(self, preview):
        """
         Convert ImagePreview into regular Image
         allow fastest bgl display
        """
        w, h = preview.image_size
        image = bpy.data.images.new("Icon_{}".format(preview.icon_id), width=w, height=h)
        # copy pixels for full image trick: convert to a list make it fast!! [:]
        image.pixels = preview.image_pixels_float[:]
        return image

    def menuitems(self, category):
        """
         Build preset menu thumb compatible enum
         Using Image
        """
        # refresh thumbs
        self.load_presets(category)

        folders = self._folders[category]
        items = list(self[category].items())
        items.sort(key=lambda i: i[0].lower())
        return [(
                folders[name]['full_path'],
                folders[name]['label'],
                self.as_blender_image(preview)
                )
                for name, preview in items]


icons = IconsCollectionManager()


def unregister():
    global icons
    icons.cleanup()
