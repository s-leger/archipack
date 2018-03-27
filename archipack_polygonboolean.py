
# See archipack_cutter

from .archipack_2d import Line


class CutterSegment(Line):

    def __init__(self, p, v, type='DEFAULT'):
        Line.__init__(self, p, v)
        self.type = type
    
    @property
    def copy(self):
        s = Cutter(self.p.copy(), self.v.copy(), self.type)
        return s

    def straight(self, length, t=1):
        s = self.copy
        s.p = self.lerp(t)
        s.v = self.v.normalized() * length
        return s

    def set_offset(self, offset, last=None):
        """
            Offset line and compute intersection point
            between segments
        """
        self.line = self.make_offset(offset, last)

    def offset(self, offset):
        s = self.copy
        s.p += offset * self.cross_z.normalized()
        return s

    @property
    def oposite(self):
        s = self.copy
        s.p += s.v
        s.v = -s.v
        return s

   
class CutAblePolygon():
    """
        Simple boolean operations
        Cutable generator
        Object MUST have properties
        - segs
        - holes
        - convex
    """
    
    def inside(self, pt, segs=None):
        """
            Point inside poly (raycast method)
            support concave polygons
            TODO:
            make s1 angle different than all othr segs
        """
        s1 = Line(pt, Vector((100 * self.xsize, 0.1)))
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
            print("slice failed found:%s of:%s" % (len(segs), s_nsegs + b_nsegs))
            return False

        for i, s in enumerate(segs):
            s.p0 = segs[i - 1].p1

        return True

    def slice(self, boundary):
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
        keep_inside = (boundary.operation == 'INTERSECTION')

        start = -1
        cutter = boundary
        roof = self

        f_segs = self.segs
        c_segs = boundary.segs
        store = []

        slice_res = True
        is_inside = False

        # find if either a boundary or
        # cutter intersects
        # (at least one point of any must be inside other one)

        # find a point of this pitch inside boundary
        for i, s in enumerate(f_segs):
            res = self.inside(s.p0, c_segs)
            if res:
                is_inside = True
            if res == keep_inside:
                start = i
                # print("pitch pt %sside f_start:%s %s" % (in_out, start, self.side))
                slice_res = self.get_intersections(roof, cutter, start, store, True)
                break

        # seek for point of boundary inside pitch
        for i, s in enumerate(c_segs):
            res = self.inside(s.p0)
            if res:
                is_inside = True
            # no pitch point found inside boundary
            if start < 0 and res == keep_inside:
                start = i
                # print("boundary pt %sside c_start:%s %s" % (in_out, start, self.side))
                # swap boundary / pitch so we start from boundary
                slice_res = self.get_intersections(cutter, roof, start, store, False)
                break

        # no points found at all
        if start < 0:
            # print("no pt inside")
            return

        if not slice_res:
            # print("slice fails")
            # found more segments than input
            # cutter made more than one loop
            return

        if len(store) < 1:
            if is_inside:
                # print("not touching, add as hole")
                self.holes.append(boundary)
            return

        self.segs = store
        self.is_convex()
        self.limits()

        
class CutterGenerator():

    def __init__(self, d):
        self.parts = d.parts
        self.operation = d.operation
        self.segs = []
    
    def add_part(self, part):
        
        if len(self.segs) < 1:
            s = None
        else:
            s = self.segs[-1]

        # start a new Cutter
        if s is None:
            v = part.length * Vector((cos(part.a0), sin(part.a0)))
            s = CutterSegment(Vector((0, 0)), v, part.type)
        else:
            s = s.straight(part.length).rotate(part.a0)
            s.type = part.type
            
        self.segs.append(s)
    
    def add_part(self, part):
        raise NotImplementedError

    def set_offset(self):
        last = None
        for i, seg in enumerate(self.segs):
            seg.set_offset(self.parts[i].offset, last)
            last = seg.line

    def close(self):
        # Make last segment implicit closing one
        s0 = self.segs[-1]
        s1 = self.segs[0]
        dp = s1.p0 - s0.p0
        s0.v = dp

        if len(self.segs) > 1:
            s0.line = s0.make_offset(self.parts[-1].offset, self.segs[-2].line)

        p1 = s1.line.p1
        s1.line = s1.make_offset(self.parts[0].offset, s0.line)
        s1.line.p1 = p1

    def locate_manipulators(self):
        if self.operation == 'DIFFERENCE':
            side = -1
        else:
            side = 1
        for i, f in enumerate(self.segs):

            manipulators = self.parts[i].manipulators
            p0 = f.p0.to_3d()
            p1 = f.p1.to_3d()
            # angle from last to current segment
            if i > 0:

                if i < len(self.segs) - 1:
                    manipulators[0].type_key = 'ANGLE'
                else:
                    manipulators[0].type_key = 'DUMB_ANGLE'

                v0 = self.segs[i - 1].straight(-side, 1).v.to_3d()
                v1 = f.straight(side, 0).v.to_3d()
                manipulators[0].set_pts([p0, v0, v1])

            # segment length
            manipulators[1].type_key = 'SIZE'
            manipulators[1].prop1_name = "length"
            manipulators[1].set_pts([p0, p1, (side, 0, 0)])

            # snap manipulator, dont change index !
            manipulators[2].set_pts([p0, p1, (side, 0, 0)])
            # dumb segment id
            manipulators[3].set_pts([p0, p1, (side, 0, 0)])

            # offset
            manipulators[4].set_pts([
                p0,
                p0 + f.sized_normal(0, max(0.0001, self.parts[i].offset)).v.to_3d(),
                (0.5, 0, 0)
            ])

    def change_coordsys(self, fromTM, toTM):
        """
            move shape fromTM into toTM coordsys
        """
        dp = (toTM.inverted() * fromTM.translation).to_2d()
        da = toTM.row[1].to_2d().angle_signed(fromTM.row[1].to_2d())
        ca = cos(da)
        sa = sin(da)
        rM = Matrix([
            [ca, -sa],
            [sa, ca]
            ])
        for s in self.segs:
            tp = (rM * s.p0) - s.p0 + dp
            s.rotate(da)
            s.translate(tp)

    def get_index(self, index):
        n_segs = len(self.segs)
        if index >= n_segs:
            index -= n_segs
        return index


