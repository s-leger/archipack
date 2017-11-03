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

from .algorithms import BoundaryNodeRule
from .constants import (
    logger,
    TopologyException
)
from .precision import CommonBitsRemover
from .op_overlay import (
    GeometrySnapper
    )
from .op_simple import IsSimpleOp
from .op_valid import IsValidOp


CBR_BEFORE_SNAPPING = False


def check_valid(geom, label: str, doThrow: bool=False, validOnly: bool=False) -> bool:

    if geom.geometryTypeId in [1, 5]:
        # Lineal geoms
        sop = IsSimpleOp(geom, BoundaryNodeRule.getBoundaryEndPoint())
        if not sop.is_simple:
            if doThrow:
                raise TopologyException("%s is not simple".format(label))
            return False
        else:
            ivo = IsValidOp(geom)
            if not ivo.is_valid:
                if doThrow:
                    raise TopologyException("%s is invalid".format(label))
                return False

    return True


def fix_self_intersections(geom, label: str):
    # Only multi-components can be fixed by UnaryUnion
    if geom.geometryTypeId < 4:
        return geom

    ivo = IsValidOp(geom)
    # poygon is valid, nothing to do
    if ivo.is_valid:
        return geom

    # Not all invalidities can be fixed by this code
    return geom


def SnapOp(geom0, geom1, _Op):
    snapTolerance = GeometrySnapper.computeOverlaySnapTolerance(geom0, geom1)
    if CBR_BEFORE_SNAPPING:
        cbr = CommonBitsRemover()
        cbr.add(geom0)
        cbr.add(geom1)
        rG0 = cbr.removeCommonBits(geom0.clone())
        rG1 = cbr.removeCommonBits(geom1.clone())
    else:
        rG0, rG1 = geom0, geom1

    snapper0 = GeometrySnapper(rG0)
    snapG0 = snapper0.snapTo(rG1, snapTolerance)

    snapper1 = GeometrySnapper(rG1)
    snapG1 = snapper1.snapTo(snapG0, snapTolerance)

    result = _Op(snapG0, snapG1)
    check_valid(result, "SNAP: result (before common-bits addition")

    if CBR_BEFORE_SNAPPING:
        cbr.addCommonBits(result)

    return result


def BinaryOp(geom0, geom1, _Op):

    origException = None

    try:
        res = _Op.execute(geom0, geom1)
        check_valid(res, "Overlay result between original inputs", True, True)
        logger.debug("Attempt with original input succeeded")
        return res
    except TopologyException as ex:
        origException = ex
        pass

    check_valid(geom0, "Input geom 0", True, True)
    check_valid(geom1, "Input geom 1", True, True)

    # USE_COMMONBITS_POLICY
    # USE_SNAPPING_POLICY
    # USE_PRECISION_REDUCTION_POLICY
    # USE_TP_SIMPLIFY_POLICY
    raise origException
