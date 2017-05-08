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
from bpy.ops import op_as_string


def operator_exists(idname):
    try:
        op_as_string(idname)
        return True
    except:
        return False


class Shortcuts:
    """
        define shortcuts mapping available
        check for specific event in modal
        build gl label including name, type and modifiers
    """
    def __init__(self, context):
        self._test = {}
        self._keys = {}

    def add_key(self, context, identifier, keymap_name, idname, key_type=None):
        """
            Add keymap by identifier
            identifier : string your identifier for this keymap
            keymap_name : string keymaps name as category in user preferences dialog
            idname : keymap_items idname as in user preferences dialog
            key_type : allow override of key.type eg: for MOUSEWHEEL events
        """
        key = self.get_keymap_item(context, keymap_name, idname)
        if key is not None:
            if key_type is None:
                key_type = key.type
            if identifier not in self._test.keys():
                self._test[identifier] = set()
                self._keys[identifier] = set()
            self._test[identifier].append((key.alt, key.ctrl, key.shift, key_type, key.value))
            self._keys[identifier].append(key)

    def check(self, event, identifier):
        """
            check event for specific identifier
            return Boolean
                True when event match
                False when not match or not found
        """
        evkey = (event.alt, event.ctrl, event.shift, event.type, event.value)
        if identifier in self._test.keys():
            return evkey in self._test[identifier]
        return False

    def get_key(self, identifier):
        """
            get stored keymap by identifier
            return set of keymap / empty set
        """
        if identifier in self._keys.keys():
            return self._keys[identifier]
        return set()

    def gl_info(self, identifier, label=None):
        """
            Return gl info text based on key
            label: string info label, default use key.name
        """
        keys = self.get_key(identifier)
        for key in keys:
            key.name, key.alt, key.shift, key.ctrl

    def get_keymap_item(self, context, keymap_name, idname):
        """
            retrieve a keymap item
            keymap_name : string keymaps name as category in user preferences dialog
            idname : keymap_items idname as in user preferences dialog
            return key or None if not found
        """
        km = context.window_manager.keyconfigs.user.keymaps
        if keymap_name in km.keys():
            idx = km[keymap_name].keymap_items.find(idname)
            if idx > -1:
                return km[keymap_name].keymap_items[idx]
        return None

    def dump_keys(self, context, filename="c:\\tmp\\keymap.txt"):
        """
            Dump all keymaps to a file
            filename : string a file path to dump keymaps
        """
        str = ""
        km = context.window_manager.keyconfigs.user.keymaps
        for key in km.keys():
            str += "\n\n#--------------------------------\n{}:\n#--------------------------------\n\n".format(key)
            for sub in km[key].keymap_items.keys():
                k = km[key].keymap_items[sub]
                str += "alt:{} ctrl:{} shift:{} type:{} value:{}  idname:{} name:{}\n".format(
                    k.alt, k.ctrl, k.shift, k.type, k.value, sub, k.name)
        file = open(filename, "w")
        file.write(str)
        file.close()
