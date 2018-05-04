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
# Cutter / CutAble shared by roof, slab, and floor
# ----------------------------------------------------------
from mathutils import Vector, Matrix
import bmesh
from random import uniform
from bpy.props import (
    FloatProperty, BoolProperty,
    EnumProperty
    )
from .archipack_2d import Line
from .archipack_curveman import ArchipackUserDefinedPath
from .archipack_segments import StraightSegment, Generator, ArchipackSegment


class CutterSegment(StraightSegment):

    def __init__(self, p, v, side_type='DEFAULT'):
        StraightSegment.__init__(self, p, v)
        self.side_type = side_type
        self.is_hole = True

    @property
    def copy(self):
        return CutterSegment(self.p.copy(), self.v.copy(), self.side_type)

    def straight(self, length, t=1):
        s = self.copy
        s.p = self.lerp(t)
        s.v = self.v.normalized() * length
        return s

    def offset(self, offset):
        s = self.copy
        s.p += self.sized_normal(0, offset).v
        return s

    @property
    def oposite(self):
        s = self.copy
        s.p += s.v
        s.v = -s.v
        return s


class CutterGenerator(Generator):
    """
      Generator for cutter objects
    """
    def __init__(self, d, o=None):
        Generator.__init__(self, d, o)
        self.operation = d.operation

    def add_part(self, part):

        if len(self.segs) < 1:
            last = None
        else:
            last = self.segs[-1]

        if hasattr(part, "side_type"):
            side_type = part.side_type
        else:
            side_type = 'DEFAULT'

        # start a new Cutter
        if last is None:
            p = self.location
            v = (self.rot * Vector((part.length, 0, 0))).to_2d()
        else:
            p = last.p1
            v = last.v.normalized() * part.length

        s = CutterSegment(p, v, side_type).rotate(part.a0)
        self.segs.append(s)

    def locate_manipulators(self):
        if self.operation == 'DIFFERENCE':
            side = -1
        else:
            side = 1
        Generator.locate_manipulators(self, side=side)

    def get_index(self, index):
        n_segs = len(self.segs)
        if index >= n_segs:
            index -= n_segs
        return index

    def next_seg(self, index):
        idx = self.get_index(index + 1)
        return self.segs[idx]

    def last_seg(self, index):
        return self.segs[index - 1]

    def get_coords(self, verts):
        for s in self.segs:
            verts.append(s.line.p0.to_3d())

    def get_verts(self, verts, edges):

        n_segs = len(self.segs) - 1

        self.get_coords(verts)

        for i in range(n_segs):
            edges.append([i, i + 1])


class CutAblePolygon():
    """
        Simple boolean operations (mod 2 rule based)
        Cutable generator / polygon
        Object MUST have properties
        - segs
        - holes
        - convex
    """
    def as_lines(self, step_angle=0.104):
        """
            Convert curved segments to straight lines
            Use offset ones when set
        """
        segs = []
        for s in self.segs:
            if "Curved" in type(s).__name__:
                dt, steps = s.steps_by_angle(step_angle)
                segs.extend(s.as_lines(steps))
            else:
                if s.line is None:
                    segs.append(s)
                else:
                    segs.append(s.line)
        self.segs = segs

    def inside(self, pt, segs=None):
        """
            Point inside poly (raycast method)
            support concave polygons
            TODO:
            make s1 angle different than all othr segs
        """
        s1 = Line(pt, Vector((min(10000, 100 * self.xsize), uniform(-0.5, 0.5))))
        counter = 0
        if segs is None:
            segs = self.segs
        for s in segs:
            res, p, t, u = s.intersect_ext(s1)
            if res:
                counter += 1
        return counter % 2 == 1

    def get_index(self, index):
        n_segs = len(self.segs)
        if index >= n_segs:
            index -= n_segs
        return index

    def is_convex(self):
        n_segs = len(self.segs)
        self.convex = True
        sign = False
        s0 = self.segs[-1]
        for i in range(n_segs):
            s1 = self.segs[i]
            if "Curved" in type(s1).__name__:
                self.convex = False
                return
            c = s0.v.cross(s1.v)
            if i == 0:
                sign = (c > 0)
            elif sign != (c > 0):
                self.convex = False
                return
            s0 = s1

    def get_intersections(self, border, cutter, s_start, segs, start_by_hole):
        """
            Detect all intersections
            for boundary: store intersection point, t, idx of segment, idx of cutter
            sort by t
        """
        s_segs = border.segs
        b_segs = cutter.segs
        s_nsegs = len(s_segs)
        b_nsegs = len(b_segs)
        inter = []

        # find all intersections
        for idx in range(s_nsegs):
            s_idx = border.get_index(s_start + idx)
            s = s_segs[s_idx]
            for b_idx, b in enumerate(b_segs):
                res, p, u, v = s.intersect_ext(b)
                if res:
                    inter.append((s_idx, u, b_idx, v, p))

        # print("%s" % (self.side))
        # print("%s" % (inter))

        if len(inter) < 1:
            return True

        # sort by seg and param t of seg
        inter.sort()

        # reorder so we realy start from s_start
        for i, it in enumerate(inter):
            if it[0] >= s_start:
                order = i
                break

        inter = inter[order:] + inter[:order]

        # print("%s" % (inter))
        p0 = border.segs[s_start].p0

        n_inter = len(inter) - 1

        for i in range(n_inter):
            s_end, u, b_start, v, p = inter[i]
            s_idx = border.get_index(s_start)
            s = s_segs[s_idx].copy
            s.is_hole = not start_by_hole
            segs.append(s)
            idx = s_idx
            max_iter = s_nsegs
            # walk through s_segs until intersection
            while s_idx != s_end and max_iter > 0:
                idx += 1
                s_idx = border.get_index(idx)
                s = s_segs[s_idx].copy
                s.is_hole = not start_by_hole
                segs.append(s)
                max_iter -= 1
            segs[-1].p1 = p

            s_start, u, b_end, v, p = inter[i + 1]
            b_idx = cutter.get_index(b_start)
            s = b_segs[b_idx].copy
            s.is_hole = start_by_hole
            segs.append(s)
            idx = b_idx
            max_iter = b_nsegs
            # walk through b_segs until intersection
            while b_idx != b_end and max_iter > 0:
                idx += 1
                b_idx = cutter.get_index(idx)
                s = b_segs[b_idx].copy
                s.is_hole = start_by_hole
                segs.append(s)
                max_iter -= 1
            segs[-1].p1 = p

        # add part between last intersection and start point
        idx = s_start
        s_idx = border.get_index(s_start)
        s = s_segs[s_idx].copy
        s.is_hole = not start_by_hole
        segs.append(s)
        max_iter = s_nsegs
        # go until end of segment is near start of first one
        while (s_segs[s_idx].p1 - p0).length > 0.0001 and max_iter > 0:
            idx += 1
            s_idx = border.get_index(idx)
            s = s_segs[s_idx].copy
            s.is_hole = not start_by_hole
            segs.append(s)
            max_iter -= 1

        if len(segs) > s_nsegs + b_nsegs + 1:
            # print("slice failed found:%s of:%s" % (len(segs), s_nsegs + b_nsegs))
            return False

        for i, s in enumerate(segs):
            s.p0 = segs[i - 1].p1

        return True

    def slice(self, cutter):
        """
            Simple 2d Boolean between boundary and roof part
            Dosen't handle slicing roof into multiple parts

            4 cases:
            1 pitch has point in boundary -> start from this point
            2 boundary has point in pitch -> start from this point
            3 no points inside -> find first crossing segment
            4 not points inside and no crossing segments
        """
        # print("************")

        # keep inside or cut inside
        # keep inside must be CCW
        # cut inside must be CW
        keep_inside = (cutter.operation == 'INTERSECTION')

        start = -1

        f_segs = self.segs
        c_segs = cutter.segs
        store = []

        slice_res = True
        is_inside = False

        # find if either a cutter or
        # cutter intersects
        # (at least one point of any must be inside other one)

        # find if all points of the object are inside cutter
        # so there is no need to cut when mode is intersection
        if keep_inside:
            one_outside = False
            for i, s in enumerate(f_segs):
                res = self.inside(s.p0, c_segs)
                if not res:
                    one_outside = True
                    break
            if not one_outside:
                return True

        # find a point of this pitch inside cutter
        for i, s in enumerate(f_segs):
            res = self.inside(s.p0, c_segs)
            if res:
                is_inside = True
            if res == keep_inside:
                start = i
                # print("pitch pt %sside f_start:%s %s" % (in_out, start, self.side))
                slice_res = self.get_intersections(self, cutter, start, store, True)
                break

        # seek for point of cutter inside pitch
        for i, s in enumerate(c_segs):
            res = self.inside(s.p0)
            if res:
                is_inside = True
            # no pitch point found inside cutter
            if start < 0 and res == keep_inside:
                start = i
                # print("cutter pt %sside c_start:%s %s" % (in_out, start, self.side))
                # swap cutter / pitch so we start from cutter
                slice_res = self.get_intersections(cutter, self, start, store, False)
                break

        # no points found at all
        if start < 0:
            # print("no pt inside")
            return not keep_inside

        if not slice_res:
            # print("slice fails")
            # found more segments than input
            # cutter made more than one loop
            return True

        if len(store) < 1:
            if is_inside:
                # print("not touching, add as hole")
                if keep_inside:
                    self.segs = cutter.segs
                else:
                    self.holes.append(cutter)

            return True

        self.segs = store
        self.is_convex()

        return True


class CutAbleGenerator(Generator):
    """
     Generator for cutable objects
    """
    def __init__(self, d, o=None):
        Generator.__init__(self, d, o)
        self.holes = []
        self.convex = True
        self.xsize = 0

    def limits(self):
        x_size = [s.p0.x for s in self.segs]
        self.xsize = max(x_size) - min(x_size)

    def bissect(self, bm,
            plane_co,
            plane_no,
            dist=0.001,
            use_snap_center=False,
            clear_outer=True,
            clear_inner=False
            ):
        geom = bm.verts[:]
        geom.extend(bm.edges[:])
        geom.extend(bm.faces[:])

        bmesh.ops.bisect_plane(bm,
            geom=geom,
            dist=dist,
            plane_co=plane_co,
            plane_no=plane_no,
            use_snap_center=False,
            clear_outer=clear_outer,
            clear_inner=clear_inner
            )

    def cut_holes(self, bm, cutable, offset={'DEFAULT': 0}):
        o_keys = offset.keys()
        has_offset = len(o_keys) > 1 or offset['DEFAULT'] != 0
        # cut holes
        for hole in cutable.holes:

            if has_offset:

                for s in hole.segs:
                    if s.length > 0:
                        if s.side_type in o_keys:
                            of = offset[s.side_type]
                        else:
                            of = offset['DEFAULT']
                        n = s.sized_normal(0, 1).v
                        p0 = s.p0 + n * of
                        self.bissect(bm, p0.to_3d(), n.to_3d(), clear_outer=False)

                # compute boundary with offset
                new_s = None
                segs = []
                for s in hole.segs:
                    if s.length > 0:
                        if s.side_type in o_keys:
                            of = offset[s.side_type]
                        else:
                            of = offset['DEFAULT']
                        new_s = s.make_offset(of, new_s)
                        segs.append(new_s)
                # last / first intersection
                if len(segs) > 0:
                    res, p0, t = segs[0].intersect(segs[-1])
                    if res:
                        segs[0].p0 = p0
                        segs[-1].p1 = p0

            else:
                for s in hole.segs:
                    if s.length > 0:
                        n = s.sized_normal(0, 1).v
                        self.bissect(bm, s.p0.to_3d(), n.to_3d(), clear_outer=False)
                # use hole boundary
                segs = hole.segs

            if len(segs) > 0:
                # when hole segs are found clear parts inside hole
                f_geom = [f for f in bm.faces
                    if cutable.inside(
                        f.calc_center_median().to_2d(),
                        segs=segs)]
                if len(f_geom) > 0:
                    bmesh.ops.delete(bm, geom=f_geom, context=5)

    def cut_boundary(self, bm, cutable, offset={'DEFAULT': 0}):
        o_keys = offset.keys()
        has_offset = len(o_keys) > 1 or offset['DEFAULT'] != 0
        # cut outside parts
        if has_offset:
            for s in cutable.segs:
                if s.length > 0:
                    if s.side_type in o_keys:
                        of = offset[s.side_type]
                    else:
                        of = offset['DEFAULT']
                    n = s.sized_normal(0, 1).v
                    p0 = s.p0 + n * of
                    self.bissect(bm, p0.to_3d(), n.to_3d(), clear_outer=cutable.convex)
        else:
            for s in cutable.segs:
                if s.length > 0:
                    n = s.sized_normal(0, 1).v
                    self.bissect(bm, s.p0.to_3d(), n.to_3d(), clear_outer=cutable.convex)

        if not cutable.convex:
            f_geom = [f for f in bm.faces
                if not cutable.inside(f.calc_center_median().to_2d())]
            if len(f_geom) > 0:
                bmesh.ops.delete(bm, geom=f_geom, context=5)


class ArchipackCutterPart(ArchipackSegment):
    """
        Cutter segment PropertyGroup
        Childs MUST implements
        -get_datablock
    """

    def update(self, context, update_parent=True):
        idx, o, d = self.find_datablock_in_selection(context)
        if d is not None:
            d.update(context, update_parent=update_parent)


def update_operation(self, context):
    g = self.get_generator()
    pts = [seg.p0.to_3d() for seg in g.segs]
    if self.is_cw(pts) != (self.operation == 'INTERSECTION'):
        return
    o = self.find_in_selection(context)
    if o is None:
        return
    self.reverse(context, o)


def update_path(self, context):
    self.update_path(context)


def update(self, context):
    self.update(context, update_parent=True)


def update_manipulators(self, context):
    self.update(context, manipulable_refresh=True)


class ArchipackCutter(ArchipackUserDefinedPath):

    z = FloatProperty(
            name="dumb z",
            description="Dumb z for manipulator placeholder",
            default=0.01,
            options={'SKIP_SAVE'}
            )
    operation = EnumProperty(
            items=(
                ('DIFFERENCE', 'Difference', 'Cut inside part', 0),
                ('INTERSECTION', 'Intersection', 'Keep inside part', 1)
                ),
            default='DIFFERENCE',
            update=update_operation
            )
    auto_update = BoolProperty(
            options={'SKIP_SAVE'},
            default=True,
            update=update_manipulators
            )
    closed = BoolProperty(
            options={'SKIP_SAVE'},
            description="keep closed to be wall snap manipulator compatible",
            default=True
            )
    always_closed = BoolProperty(
            options={'SKIP_SAVE'},
            name="Always closed geometry",
            description="Flag indicate whenever geometry path is always closed",
            default=True
            )
    offset = FloatProperty(
            default=0,
            name="Offset",
            description="Lateral offset of cutter",
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )

    def draw(self, context, layout, draw_offset=False, draw_type=False):
        box = layout.box()
        box.operator('archipack.manipulate', icon='HAND')
        box.prop(self, 'operation', text="")
        box = layout.box()
        self.template_user_path(context, box)
        if draw_offset:
            box.prop(self, 'offset')
        self.template_parts(context, layout, draw_type=draw_type)

    def update_parent(self, context):
        raise NotImplementedError

    def setup_manipulators(self):
        self.setup_parts_manipulators('z')

    def get_generator(self, o=None):
        g = CutterGenerator(self, o)
        for part in self.parts:
            g.add_part(part)
        g.set_offset(self.offset)
        g.close(self.offset)
        return g

    def ensure_direction(self, o=None):
        # get segs ensure they are cw or ccw depending on operation
        # whatever the user do with points
        g = self.get_generator(o)
        pts = [seg.p0.to_3d() for seg in g.segs]
        if self.is_cw(pts) != (self.operation == 'INTERSECTION'):
            return g

        g.segs = [s.oposite for s in reversed(g.segs)]

        return g

    def from_spline(self, context, wM, resolution, spline):
        o = self.find_in_selection(context, self.auto_update)
        if o is None:
            return

        pts = self.coords_from_spline(
            spline,
            wM,
            resolution,
            ccw=(self.operation == 'INTERSECTION'),
            cw=(self.operation != 'INTERSECTION'),
            close=True
            )

        if len(pts) < 3:
            return

        # pretranslate
        o.matrix_world = Matrix.Translation(pts[0].copy())
        auto_update = self.auto_update
        self.auto_update = False
        self.from_points(pts)
        self.auto_update = auto_update
        self.update_parent(context, o)

    def after_reverse(self, context, o):
        self.auto_update = True
        self.update_parent(context, o)

    def make_surface(self, o, verts, edges):
        bm = bmesh.new()
        for v in verts:
            bm.verts.new(v)
        bm.verts.ensure_lookup_table()
        for ed in edges:
            bm.edges.new((bm.verts[ed[0]], bm.verts[ed[1]]))
        bm.edges.new((bm.verts[-1], bm.verts[0]))
        bm.edges.ensure_lookup_table()
        bm.to_mesh(o.data)
        bm.free()

    def get_coords(self):
        """
         return coordinates in object coordsys
        """
        self.update_parts()
        verts = []
        g = self.get_generator()
        g.get_coords(verts)
        return verts

    def update(self, context, manipulable_refresh=False, update_parent=True):
        """
         Does update parent make sense at all ?
         as cutter changes must always update parent
        """
        o = self.find_in_selection(context, self.auto_update)

        if o is None:
            return

        # clean up manipulators before any data model change
        if manipulable_refresh:
            self.manipulable_disable(context)

        self.update_parts()

        verts = []
        edges = []

        g = self.get_generator()
        g.locate_manipulators()

        # vertex index in order to build axis
        g.get_verts(verts, edges)

        if len(verts) > 2:
            self.make_surface(o, verts, edges)

        # enable manipulators rebuild
        if manipulable_refresh:
            self.manipulable_refresh = True

        # update parent on direct edit
        if manipulable_refresh or update_parent:
            self.update_parent(context, o)

        self.update_dimensions(context, o)

        # restore context
        self.restore_context(context)

    def manipulable_setup(self, context):

        self.manipulable_disable(context)
        o = context.active_object

        n_parts = self.n_parts + 1

        self.setup_manipulators()

        for i, part in enumerate(self.parts):
            if i < n_parts:

                if i > 0:
                    # start angle
                    self.manip_stack.append(part.manipulators[0].setup(context, o, part))

                # length
                self.manip_stack.append(part.manipulators[1].setup(context, o, part))
                # index
                self.manip_stack.append(part.manipulators[3].setup(context, o, self))
                # offset
                # self.manip_stack.append(part.manipulators[4].setup(context, o, part))

            # snap point
            self.manip_stack.append(part.manipulators[2].setup(context, o, self))
