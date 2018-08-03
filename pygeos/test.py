"""
from archipack.pygeos import test as test
tb = test.TestBuffer()
tb.test_buffer0()
tb.test_poly7()


ta = test.TestAlgorithms()
ta.test_isPointInRing()
ta.test_lineIntersector()
ta.test_pointLocator()

tp = ta = test.TestPrecision()
tp.test_reduce()
tp.test_collapse()


"""

import unittest
from .geom import GeometryFactory
from .algorithms import (
    CGAlgorithms,
    PointLocator,
    LineIntersector
    )
from .shared import (
    GeomTypeId,
    Location,
    PrecisionModel,
    logger
    )
from .precision import (
    GeometryPrecisionReducer
    )
from .op_buffer import BufferParameters, BufferOp
        

gf = GeometryFactory()


def createCoord(co):
    return gf.createCoordinate(co)


def createCoords(coords):
    return gf.coordinateSequenceFactory.create([createCoord(co) for co in coords])


def createLinarRing(coords):
    return gf.createLinearRing(createCoords(coords))


def createLineString(coords):
    return gf.createLineString(createCoords(coords))


def createCollection(geoms):
    return gf.buildGeometry(geoms)


def createPolygon(coords):
    return gf.createPolygon(createLinarRing(coords))


def createPoint(coord):
    return gf.createPoint(createCoord((10, 10)))


def point_equals(p0, p1, tolerance):
    dist = p0.distance(p1)
    logger.debug("point_equals %s" , dist)
    return dist <= tolerance


def test_lineIntersector(tester, p0, p1, q0, q1, expectedIntersectionNum, intPt=[], tolerance=0):

    li = LineIntersector()
    li.computeLinesIntersection(p0, p1, q0, q1)

    tester.assertEqual(li.intersections, expectedIntersectionNum)

    if len(intPt) == 0:
        return

    tester.assertEqual(li.intersections, len(intPt))

    if li.intersections == 1:
        logger.debug("%s == %s", intPt[0], li.intersectionPts[0])
        tester.assertTrue(point_equals(intPt[0], li.intersectionPts[0], tolerance))

    if li.intersections == 2:
        logger.debug("%s == %s", intPt[0], li.intersectionPts[0])
        logger.debug("%s == %s", intPt[1], li.intersectionPts[1])
        eq0 = point_equals(intPt[0], li.intersectionPts[0], tolerance) or \
            point_equals(intPt[1], li.intersectionPts[0], tolerance)
        eq1 = point_equals(intPt[0], li.intersectionPts[1], tolerance) or \
            point_equals(intPt[1], li.intersectionPts[1], tolerance)
        tester.assertTrue(eq0 and eq1)


class TestAlgorithms(unittest.TestCase):

    """Check algorithms"""
    def test_isPointInRing(self):
        """IsPointInRing"""
        coords = [(0, 0), (0, 20), (20, 20), (20, 0), (0, 0)]
        poly = createCoords(coords)
        pt = createCoord((10, 10))
        res = CGAlgorithms.isPointInRing(pt, poly)
        self.assertTrue(res)

        coords = [(-40, 80), (-40, -80), (20, 0), (20, -100), (40, 40), (80, -80),
            (100, 80), (140, -20), (120, 140), (40, 180), (60, 40), (0, 120), (-20, -20), (-40, 80)]
        poly = createCoords(coords)
        pt = createCoord((0, 0))
        res = CGAlgorithms.isPointInRing(pt, poly)
        self.assertTrue(res)

    def test_pointLocator(self):
        """pointLocator"""
        locator = PointLocator()
        coords = [(-40, 80), (-40, -80), (20, 0), (20, -100), (40, 40), (80, -80),
            (100, 80), (140, -20), (120, 140), (40, 180), (60, 40), (0, 120), (-20, -20), (-40, 80)]
        poly = createPolygon(coords)
        c0 = createCoord((0, 0))
        loc = locator.locate(c0, poly)
        self.assertEqual(loc, Location.INTERIOR)

        ls = createLineString([(0, 0), (10, 10)])
        lr = createLinarRing([(10, 10), (10, 20), (20, 10), (10, 10)])
        coll = createCollection([ls, lr])
        loc = locator.locate(c0, coll)
        self.assertEqual(loc, Location.BOUNDARY)

        c1 = createCoord((11, 11))
        loc = locator.locate(c1, lr)
        self.assertEqual(loc, Location.INTERIOR)

        # fails probably due to object type
        c1 = createCoord((1, 1))
        pt0 = createPoint(c0)
        pt1 = createPoint(c1)
        coll = createCollection([pt0, pt1])
        loc = locator.locate(c1, coll)
        self.assertEqual(loc, Location.INTERIOR)

    def test_lineIntersector(self):
        tolerance = 1e-12
        
        p0 = createCoord((0, 0))
        p1 = createCoord((1, 0))
        q0 = createCoord((0, 0))
        q1 = createCoord((0, 1))
        pt = createCoord((0, 0))
        test_lineIntersector(self, p0, p1, q0, q1, 1, [pt], tolerance)
        
        p0 = createCoord((0, 0))
        p1 = createCoord((1, 0))
        q0 = createCoord((0, 0))
        q1 = createCoord((-1, 0))
        pt = createCoord((0, 0))
        test_lineIntersector(self, p0, p1, q0, q1, 1, [pt], tolerance)

        # fails in geos
        p0 = createCoord((588750.7429703881, 4518950.493668233))
        p1 = createCoord((588748.2060409798, 4518933.9452804085))
        q0 = createCoord((588745.824857241, 4518940.742239175))
        q1 = createCoord((588748.2060437313, 4518933.9452791475))
        pt = createCoord((588748.2060416829, 4518933.945284994))
        test_lineIntersector(self, p0, p1, q0, q1, 1, [pt], tolerance)

        # fails in geos
        p0 = createCoord((588743.626135934, 4518924.610969561))
        p1 = createCoord((588732.2822865889, 4518925.4314047815))
        q0 = createCoord((588739.1191384895, 4518927.235700594))
        q1 = createCoord((588731.7854614238, 4518924.578370095))
        pt = createCoord((588733.8306132929, 4518925.319423238))
        test_lineIntersector(self, p0, p1, q0, q1, 1, [pt], tolerance)


        # fails to find 2nd point
        # almost parallel (d=0.0003147125)
        p0 = createCoord((2089426.5233462777, 1180182.3877339689))
        p1 = createCoord((2085646.6891757075, 1195618.7333999649))
        q0 = createCoord((1889281.8148903656, 1997547.0560044837))
        q1 = createCoord((2259977.3672235999, 483675.17050843034))
        pt0 = createCoord((2089426.5233462777, 1180182.3877339689))
        pt1 = createCoord((2085646.6891757075, 1195618.7333999649))
        test_lineIntersector(self, p0, p1, q0, q1, 2, [pt0, pt1], tolerance)

        p0 = createCoord((-42.0, 163.2))
        p1 = createCoord((21.2, 265.2))
        q0 = createCoord((-26.2, 188.7))
        q1 = createCoord((37.0, 290.7))
        test_lineIntersector(self, p0, p1, q0, q1, 2, [p1, q0], tolerance)

        # non almost paralel lines
        p0 = createCoord((305690.0434123494, 254176.46578338774))
        p1 = createCoord((305601.9999843455, 254243.19999846347))
        q0 = createCoord((305689.6153764265, 254177.33102743194))
        q1 = createCoord((305692.4999844298, 254171.4999983967))
        pt = createCoord((305690.0434123494, 254176.46578338774))
        test_lineIntersector(self, p0, p1, q0, q1, 1, [pt], tolerance)

        p0 = createCoord((163.81867067, -211.31840378))
        p1 = createCoord((165.9174252, -214.1665075))
        q0 = createCoord((2.84139601, -57.95412726))
        q1 = createCoord((469.59990601, -502.63851732))
        pt = createCoord((163.81867067, -211.31840378))
        test_lineIntersector(self, p0, p1, q0, q1, 1, [pt], tolerance)

        p0 = createCoord((-58.00593335955, -1.43739086465))
        p1 = createCoord((-513.86101637525, -457.29247388035))
        q0 = createCoord((-215.22279674875, -158.65425425385))
        q1 = createCoord((-218.1208801283, -160.68343590235))
        pt = createCoord((-215.22279674875, -158.65425425385))
        test_lineIntersector(self, p0, p1, q0, q1, 1, [pt], tolerance)

        
class TestPrecision(unittest.TestCase):

    def test_reduce(self):
        pm = PrecisionModel(1)
        reducer = GeometryPrecisionReducer(precisionModel=pm)
        g1 = createPolygon([(0, 0), (0, 1.4), (1.4, 1.4), (1.4, 0), (0, 0)])
        g2 = createPolygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
        res = reducer._reduce(g1)
        self.assertTrue(g2.almost_equals(res))

    def test_collapse(self):
        pm = PrecisionModel(1)
        reducer = GeometryPrecisionReducer(precisionModel=pm)
        g1 = createPolygon([(0, 0), (0, 0.4), (0.4, 0.4), (0.4, 0), (0, 0)])
        res = reducer._reduce(g1)
        self.assertTrue(res.is_empty)

        
class TestBuffer(unittest.TestCase):   

    def test_poly7(self):
        coords = ((150213.068985, 524020.273985), (150226.206985, 524020.218985), (150245.019985, 524022.421985),
            (150248.570985, 524017.633985), (150243.628985, 523991.402985), (150233.942985, 523969.968985), 
            (150203.852985, 523929.546985), (150189.509985, 523905.746985), (150179.578985, 523905.795985),
            (150163.996985, 523891.573985), (150150.529985, 523872.591985), (150158.960985, 523848.710985), 
            (150150.727985, 523827.265985), (150129.075985, 523796.394985), (150110.126985, 523782.119985), 
            (150064.853985, 523787.149985), (150051.774985, 523791.993985), (150035.273985, 523796.784985), 
            (150034.124985, 523803.948985), (150047.317985, 523842.088985), (150048.538985, 523846.850985), 
            (150048.758985, 523856.362985), (150044.002985, 523858.774985), (150033.285985, 523861.216985), 
            (150022.584985, 523866.044985), (150013.239985, 523875.626985), (150010.897985, 523882.801985), 
            (150007.322985, 523904.295985), (150015.725985, 523913.798985), (150028.883985, 523920.894985), 
            (150036.292985, 523937.604985), (150033.171985, 523964.012985), (150028.264985, 524013.973985), 
            (150020.417985, 524042.804985), (150014.532985, 524064.428985), (150004.476985, 524083.491985), 
            (149987.717985, 524115.262985), (149981.881985, 524139.242985), (149991.382985, 524146.196985), 
            (150012.547985, 524165.288985), (150017.553385, 524169.126585), (150024.575985, 524166.982985), 
            (150037.645985, 524157.385985), (150054.301985, 524147.777985), (150067.231985, 524142.754985), 
            (150080.313985, 524135.548985), (150096.808985, 524132.911985), (150108.662985, 524120.938985), 
            (150113.586985, 524111.551985), (150113.285985, 524097.054985), (150114.403985, 524085.116985), 
            (150121.501985, 524075.543985), (150134.308985, 524061.036985), (150143.802985, 524053.844985), 
            (150159.042985, 524051.270985), (150177.151985, 524046.558985), (150188.764985, 524039.234985), 
            (150195.842985, 524027.285985), (150213.068985, 524020.273985))
        g0 = createPolygon(coords)
        segments = 32
        params = BufferParameters(segments)
        op = BufferOp(g0, params)
        res = op.getResultGeometry(-75.0)
        self.assertTrue(not res.is_empty)
        self.assertTrue(res.is_valid)
        self.assertTrue(res.type_id == GeomTypeId.GEOS_POLYGON)
        self.assertTrue(len(res.coords) == 24)
    
    def test_buffer0(self):
        coords = [(0, 0), (0, 2), (1, 1), (2, 2), (2, 0), (1, 1), (0, 0)]
        bowtie = createPolygon(coords)
        self.assertTrue(not bowtie.is_valid)
        clean = bowtie.buffer(0)
        # clean cant be valid as polygons toches
        self.assertTrue(not clean.is_valid)
        
        
if __name__ == 'main':
    unittest.main()
