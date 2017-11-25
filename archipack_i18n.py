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
#
# ----------------------------------------------------------
import os
import bpy
import codecs
import csv


def GetTranslationDict():
    dict = {}
    path = os.path.join(os.path.dirname(__file__), "i18n.csv")
    try:
        with codecs.open(path, 'r', 'utf-8') as f:
            reader = csv.reader(f)
            dict['fr_FR'] = {}
            for row in reader:
                for context in bpy.app.translations.contexts:
                    dict['fr_FR'][(context, row[0])] = row[1]
    except:
        print("Archipack: Error reading i18n translations")
        pass
    return dict


def register():
    translation_dict = GetTranslationDict()
    bpy.app.translations.register(__name__.split('.')[0], translation_dict)


def unregister():
    bpy.app.translations.unregister(__name__.split('.')[0])
