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
import numpy as np

class BitArray():

    def __init__(self, bitSize, fill=False):
        self.bitArray = np.array((fill,) * bitSize , dtype='bool_')
        
    def __str__(self):
        return str(self.list)

    def test(self, bit_num):
        return self.bitArray[bit_num]

    def set(self, bit_num):
        self.bitArray[bit_num] = True

    def clear(self, bit_num):
        self.bitArray[bit_num] = False

    def toggle(self, bit_num):
        self.bitArray[bit_num] = not self.bitArray[bit_num]

    @property
    def copy(self):
        copy = BitArray(0)
        copy.bitArray = self.bitArray.copy()
        return copy

    @property
    def list(self):
        return self.bitArray.nonzero()[0].tolist()

    def none(self):
        self.bitArray.fill(False)

    def reverse(self):
        self.bitArray = np.array([not b for b in self.bitArray], dtype='bool_')

    def all(self):
        self.bitArray.fill(True)

    def equals(self, other):
        return np.array_equal(self.bitArray, other)
               