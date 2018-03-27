class Segment():

    def __init__(self):
        pass

    def set_offset(self, offset, last=None):
        """
            Offset line and compute intersection point
            between segments
        """
        self.line = self.make_offset(offset, last)
    

class StraightSegment(Segment, Line):

    def __init__(self, p, v):
        Line.__init__(self, p, v)
        Segment.__init__(self)


class CurvedSegment(Segment, Arc):

    def __init__(self, c, radius, a0, da):
        Arc.__init__(self, c, radius, a0, da)
        Segment.__init__(self)

# CutAblePolygon, CutAbleGenerator
class Generator():

    def __init__(self, d):
        # Shared
        self.d = d
        self.segs = []
        
        # Cutable
        self.holes = []
        self.convex = True
        # Private
        self.xsize = 0
    
    @property
    def closed(self):
        return self.d.closed
        
    @property
    def parts(self):
        return self.d.parts
        
    def add_part(self, part):

        last = None
        
        if len(self.segs) > 0:
            last = self.segs[-1]
        
        if 'S_' in part.type:
            if last is None:
                p = Vector((0, 0))
                v = Vector((part.length, 0))
                s = StraightSegment(p, v).rotate(part.a)
            else:
                seg = last.straight(part.length).rotate(part.a)
                s = StraightSegment(seg.p, seg.v)
        else:
            if last is None:
                c = -part.radius * Vector((cos(part.a), sin(part.a)))
                s = CurvedSegment(c, part.radius, part.a, part.da)
            else:    
                n = last.normal(1).rotate(part.a).scale(part.radius)
                if part.da < 0:
                    n.v = -n.v
                a = n.angle
                c = n.p - n.v
                s = CurvedSegment(c, part.radius, a, part.da)
        
        self.segs.append(s)
        
        self.last_type = part.type

    def set_offset(self, offset):
        last = None
        for i, seg in enumerate(self.segs):
            seg.set_offset(offset + self.parts[i].offset, last)
            last = seg.line

    def close(self, offset):
        # Make last segment implicit closing one
        if self.closed:
            _parts = self.parts
            _segs = self.segs
            
            part = _parts[-1]
            s0 = _segs[-1]
            s1 = _segs[0]
            dp = _s1.p0 - s0.p0
            
            if "C_" in part.type:
                ds = (s0.p1 - s0.p0)
                s0.r = part.radius / ds.length * dp.length
                # angle pt - p0        - angle p0 p1
                da = atan2(dp.y, dp.x) - atan2(ds.y, ds.x)
                a0 = s0.a0 + da
                if a0 > pi:
                    a0 -= 2 * pi
                if a0 < -pi:
                    a0 += 2 * pi
                s0.a0 = a0
            else:
                s0.v = dp

            if len(_segs) > 1:
                s0.line = s0.make_offset(offset + part.offset, _segs[-2].line)
                
            p1 = s1.line.p1
            s1.line = s1.make_offset(offset + _parts[0].offset, s0.line)
            s1.line.p1 = p1

    def locate_manipulators(self, side=1):
        """
            setup manipulators
        """
        _segs = self.segs
        _parts = self.parts
        
        for i, f in enumerate(_segs):
            part = _parts[i]
            manipulators = part.manipulators
            p0 = f.p0.to_3d()
            p1 = f.p1.to_3d()
            # angle from last to current segment
            if i > 0:
                v0 = _segs[i - 1].straight(-side, 1).v.to_3d()
                v1 = f.straight(side, 0).v.to_3d()
                manipulators[0].set_pts([p0, v0, v1])

            if 'S_' in part.type:
                # segment length
                manipulators[1].type_key = 'SIZE'
                manipulators[1].prop1_name = "length"
                manipulators[1].set_pts([p0, p1, (side, 0, 0)])
            else:
                # segment radius + angle
                v0 = (f.p0 - f.c).to_3d()
                v1 = (f.p1 - f.c).to_3d()
                manipulators[1].type_key = 'ARC_ANGLE_RADIUS'
                manipulators[1].prop1_name = "da"
                manipulators[1].prop2_name = "radius"
                manipulators[1].set_pts([f.c.to_3d(), v0, v1])

            # snap manipulator, dont change index !
            manipulators[2].set_pts([p0, p1, (side, 0, 0)])
            # dumb segment id
            manipulators[3].set_pts([p0, p1, (side, 0, 0)])
            
            """
            # offset
            manipulators[4].set_pts([
                p0,
                p0 + f.sized_normal(0, max(0.0001, self.parts[i].offset)).v.to_3d(),
                (0.5, 0, 0)
            ])
            """
            
            # Dimensions points
            self.d.add_dimension_point(part.uid, p0)

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
    
    def rotate(self, idx_from, a):
        """
            apply rotation to all following segs
        """
        self.segs[idx_from].rotate(a)
        ca = cos(a)
        sa = sin(a)
        rM = Matrix([
            [ca, -sa],
            [sa, ca]
            ])
        # rotation center
        p0 = self.segs[idx_from].p0
        for i in range(idx_from + 1, len(self.segs)):
            seg = self.segs[i]
            # rotate seg
            seg.rotate(a)
            # rotate delta from rotation center to segment start
            dp = rM * (seg.p0 - p0)
            seg.translate(dp)

    def translate(self, idx_from, dp):
        """
            apply translation to all following segs
        """
        self.segs[idx_from].p1 += dp
        for i in range(idx_from + 1, len(self.segs)):
            self.segs[i].translate(dp)

    def draw(self, context):
        """
            draw generator using gl
        """
        for seg in self.segs:
            seg.draw(context, render=False)

    def get_verts(self, verts):
        verts.extend([s.p0.to_3d() for s in self.segs])

    def limits(self):
        x_size = [s.p0.x for s in self.segs]
        self.xsize = max(x_size) - min(x_size)

    def cut(self, context, o):
        """
            either external or holes cuts
        """
        # use offset segs as base
        self.as_lines(step_angle=0.0502)
        self.limits()

        for b in o.children:
            d = archipack_slab_cutter.datablock(b)
            if d is not None:
                g = d.ensure_direction()
                g.change_coordsys(b.matrix_world, o.matrix_world)
                self.slice(g)

    def slab(self, context, o, d):

        verts = []
        self.get_verts(verts)
        if len(verts) > 2:

            # ensure verts are CCW
            if d.is_cw(verts):
                verts = list(reversed(verts))

            bm = bmesh.new()

            for v in verts:
                bm.verts.new(v)
            bm.verts.ensure_lookup_table()
            for i in range(1, len(verts)):
                bm.edges.new((bm.verts[i - 1], bm.verts[i]))
            bm.edges.new((bm.verts[-1], bm.verts[0]))
            bm.edges.ensure_lookup_table()
            bmesh.ops.contextual_create(bm, geom=bm.edges)

            self.cut_holes(bm, self)

            bmesh.ops.dissolve_limit(bm,
                        angle_limit=0.01,
                        use_dissolve_boundaries=False,
                        verts=bm.verts,
                        edges=bm.edges,
                        delimit=1)

            bm.to_mesh(o.data)
            bm.free()
            # geom = bm.faces[:]
            # verts = bm.verts[:]
            # bmesh.ops.solidify(bm, geom=geom, thickness=d.z)

            # merge with object
            # bmed.bmesh_join(context, o, [bm], normal_update=True)

            bpy.ops.object.mode_set(mode='OBJECT')


def update_type(self, context):
    
    # TODO:
    # Use Arc and Line
    # Instead of generator's one
    
    d = self.find_datablock_in_selection(context)

    if d is not None and d.auto_update:

        d.auto_update = False
        idx = 0
        for i, part in enumerate(d.parts):
            if part == self:
                idx = i
                break
        a0 = 0
        if idx > 0:
            g = d.get_generator()
            w0 = g.segs[idx - 1]
            a0 = w0.straight(1).angle
            if "C_" in self.type:
                w = w0.straight_wall(self.a0, self.length, d.z, self.z, self.t)
            else:
                w = w0.curved_wall(self.a0, self.da, self.radius, d.z, self.z, self.t)
        else:
            if "C_" in self.type:
                p = Vector((0, 0))
                v = self.length * Vector((cos(self.a0), sin(self.a0)))
                w = StraightWall(p, v, d.z, self.z, self.t, d.flip)
                a0 = pi / 2
            else:
                c = -self.radius * Vector((cos(self.a0), sin(self.a0)))
                w = CurvedWall(c, self.radius, self.a0, pi, d.z, self.z, self.t, d.flip)

        # w0 - w - w1
        if d.closed and idx == d.n_parts:
            dp = - w.p0
        else:
            dp = w.p1 - w.p0

        if "C_" in self.type:
            self.radius = 0.5 * dp.length
            self.da = pi
            a0 = atan2(dp.y, dp.x) - pi / 2 - a0
        else:
            self.length = dp.length
            a0 = atan2(dp.y, dp.x) - a0

        if a0 > pi:
            a0 -= 2 * pi
        if a0 < -pi:
            a0 += 2 * pi
        self.a0 = a0

        if idx + 1 < d.n_parts:
            # adjust rotation of next part
            part1 = d.parts[idx + 1]
            if "C_" in self.type:
                a0 = part1.a0 - pi / 2
            else:
                a0 = part1.a0 + w.straight(1).angle - atan2(dp.y, dp.x)

            if a0 > pi:
                a0 -= 2 * pi
            if a0 < -pi:
                a0 += 2 * pi
            part1.a0 = a0

        d.auto_update = True

 

def update_angle(self, context):
    self.update_angle(context)
        
        
def set_a(self, value):
    self.last_a = self.a
    self.a = value 
    return None

    
def get_a(self):
    return self.a


class ArchipackSegment():
    """
        A single manipulable polyline like segment
        polyline like segment line or arc based
        @TODO: share this base class with
        stair, wall, fence, slab
    """
    type = EnumProperty(
            items=(
                ('S_SEG', 'Straight', '', 0),
                ('C_SEG', 'Curved', '', 1),
                ),
            default='S_SEG',
            update=update_type
            )
            
    length = FloatProperty(
            name="Length",
            min=0.001,
            default=2.0,
            update=update
            )
            
    # Curved parts
    radius = FloatProperty(
            name="Radius",
            min=0.5,
            default=0.7,
            update=update
            )
    da = FloatProperty(
            name="Angle",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update
            )
            
    # Angle between segments       
    last_a = FloatProperty(
            description="Store last angle between segments on update",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            options={'SKIP_SAVE'}
            )
    a = FloatProperty(
            description="Store angle between segments",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            update=update_angle
            )
    a0 = FloatProperty(
            name="Start angle",
            min=-pi,
            max=pi,
            default=pi / 2,
            subtype='ANGLE', unit='ROTATION',
            set=set_a,
            get=get_a,
            options={'SKIP_SAVE'}
            )
            
    offset = FloatProperty(
            name="Offset",
            description="Add to current segment offset",
            default=0,
            unit='LENGTH', subtype='DISTANCE',
            update=update
            )
    linked_idx = IntProperty(default=-1)

    # @TODO:
    # flag to handle wall's  x_offset
    # when set add wall offset value to segment offset
    # pay attention at allowing per wall segment offset

    manipulators = CollectionProperty(type=archipack_manipulator)

    def find_in_selection(self, context):
        raise NotImplementedError

    def update(self, context, manipulable_refresh=False):
        props = self.find_in_selection(context)
        if props is not None:
            props.update(context, manipulable_refresh)

    def draw_insert(self, context, layout, index):
        """
            May implement draw for insert / remove segment operators
        """
        pass

    def draw(self, context, layout, index):
        box = layout.box()
        box.prop(self, "type", text=str(index + 1))
        self.draw_insert(context, box, index)
        if self.type in ['C_SEG']:
            box.prop(self, "radius")
            box.prop(self, "da")
        else:
            box.prop(self, "length")
        box.prop(self, "a0")
        
