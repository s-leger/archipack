# -*- coding:utf-8 -*-

# ##### BEGIN LGPL LICENSE BLOCK #####
# GEOS - Geometry Engine Open Source
# http://geos.osgeo.org
#
# Copyright (C) 2011 Sandro Santilli <strk@kbt.io>
# Copyright (C) 2005 2006 Refractions Research Inc.
# Copyright (C) 2001-2002 Vivid Solutions Inc.
# Copyright (C) 1995 Olivier Devillers <Olivier.Devillers@sophia.inria.fr>
#
# This is free software you can redistribute and/or modify it under
# the terms of the GNU Lesser General Public Licence as published
# by the Free Software Foundation.
# See the COPYING file for more information.
#
# ##### END LGPL LICENSE BLOCK #####

# <pep8 compliant>

# ----------------------------------------------------------
# Partial port (version 3.7.0) by: Stephen Leger (s-leger)
#
# ----------------------------------------------------------

from .constants import (
    Coordinate,
    CoordinateFilter
    )


class CommonBits():
    """
     * Determines the maximum number of common most-significant
     * bits in the mantissa of one or numbers.
     *
     * Can be used to compute the double-precision number which
     * is represented by the common bits.
     * If there are no common bits, the number computed is 0.0.
     *
    """
    def __init__(self):
        self.isFirst = True
        self.commonMantissaBitsCount = 53
        self.common = 0
        self.commonSignExp = 0

    @staticmethod
    def signExpBits(num: int) -> int:
        """
         * Computes the bit pattern for the sign and exponent of a
         * double-precision number.
         *
         * @param num
         * @return the bit pattern for the sign and exponent
        """
        return num >> 52

    @staticmethod
    def numCommonMostSigMantissaBits(num1: int, num2: int) -> int:
        """
         * This computes the number of common most-significant
         * bits in the mantissas of two double-precision numbers.
         *
         * It does not count the hidden bit, which is always 1.
         * It does not determine whether the numbers have the same
         * exponent - if they do not, the value computed by this
         * function is meaningless.
         * @param db
         * @return the number of common most-significant mantissa bits
        """
        count = 0
        for i in range(52, -1, -1):
            if CommonBits.getBit(num1, i) != CommonBits.getBit(num2, i):
                return count
            count += 1
        return 52

    @staticmethod
    def zeroLowerBits(bits: int, nBits: int) -> int:
        """
         * Zeroes the lower n bits of a bitstring.
         *
         * @param bits the bitstring to alter
         * @param i the number of bits to zero
         * @return the zeroed bitstring
        """
        invMask = (1 << nBits) - 1
        mask = ~ invMask
        zeroed = bits & mask
        return zeroed

    @staticmethod
    def getBit(bits: int, i: int) -> int:
        """
         * Extracts the i'th bit of a bitstring.
         *
         * @param bits the bitstring to extract from
         * @param i the bit to extract
         * @return the value of the extracted bit
        """
        mask = (1 << i)
        if (bits & mask) != 0:
            return 1
        else:
            return 0

    def add(self, num: float):
        numBits = int(num)
        if self.isFirst:
            self.common = numBits
            self.commonSignExp = CommonBits.signExpBits(numBits)
            self.isFirst = False
            return
        numSignExp = CommonBits.signExpBits(numBits)
        if numSignExp != self.commonSignExp:
            self.common = 0
            return
        self.commonMantissaBitsCount = CommonBits.numCommonMostSigMantissaBits(self.commonBits, numBits)
        self.common = CommonBits.zeroLowerBits(self.common, 64 - (12 + self.commonMantissaBitsCount))


class Translater(CoordinateFilter):

    def __init__(self, trans):
        CoordinateFilter.__init__(self)
        self.trans = trans

    def filter_ro(self, coord) -> None:
        pass

    def filter_rw(self, coord) -> None:
        coord.x += self.trans.x
        coord.y += self.trans.y


class CommonCoordinateFilter(CoordinateFilter):
    def __init__(self):
        CoordinateFilter.__init__(self)
        self.commonBitsX = CommonBits()
        self.commonBitsY = CommonBits()

    def filter_ro(self, coord) -> None:
        self.commonBitsX.add(coord.x)
        self.commonBitsY.add(coord.x)

    def filter_rw(self, coord) -> None:
        pass

    @property
    def common(self):
        return Coordinate(self.commonBitsX.common, self.commonBitsY.common)


class CommonBitsRemover():
    """
     * Allow computing and removing common mantissa bits from one or
     * more Geometries.
    """
    def __init__(self):
        # Coordinate
        self.common = None
        self._disabled = False
        self.ccFilter = CommonCoordinateFilter()

    def add(self, geom) -> None:
        """
         * Add a geometry to the set of geometries whose common bits are
         * being computed.  After this method has executed the
         * common coordinate reflects the common bits of all added
         * geometries.
         *
         * @param geom a Geometry to test for common bits
        """
        if self._disabled:
            return
        geom.apply_ro(self.ccFilter)
        self.common = self.ccFilter.common
        self._disabled = self._disabled or (self.common.x == 0.0 and self.common.y == 0.0)

    def removeCommonBits(self, geom):
        """
         * Removes the common coordinate bits from a Geometry.
         * The coordinates of the Geometry are changed.
         *
         * @param geom the Geometry from which to remove the common
         *             coordinate bits
         * @return the shifted Geometry
        """
        if self._disabled:
            return geom
        invCoord = Coordinate(-self.common.x, -self.common.y)
        trans = Translater(invCoord)
        geom.apply_rw(trans)
        geom.geometryChanged()
        return geom

    def addCommonBits(self, geom):
        """
         * Adds the common coordinate bits back into a Geometry.
         * The coordinates of the Geometry are changed.
         *
         * @param geom the Geometry to which to add the common coordinate bits
         * @return the shifted Geometry
        """
        if self._disabled:
            return geom
        trans = Translater(self.common)
        geom.apply_rw(trans)
        geom.geometryChanged()
        return geom
