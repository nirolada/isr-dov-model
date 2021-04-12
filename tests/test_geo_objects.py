import unittest
from geo_objects import (
    ECEFPoint, Point2D, Point3D, GeodeticPoint
)


class TestPoint3D(unittest.TestCase):
    def test_valid(self):
        vals = 1, 2.5, 3
        p = Point3D(*vals)
        self.assertTupleEqual(p.tup, vals)
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            Point3D((1,), '2', 3)


class TestGeodeticPoint(unittest.TestCase):
    def test_valid(self):
        lon, lat, h = 0.602138592, 0.593411946, 156.755
        p = GeodeticPoint(lon, lat, h)
        self.assertEqual(p.tup, (lon, lat, h))
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            GeodeticPoint('a', (3,), 360.5)
            
    def test_invalid_lon(self):
        with self.assertRaises(ValueError):
            GeodeticPoint(-4, 0.593411946, 156.755)
    
    def test_invalid_lat(self):
        with self.assertRaises(ValueError):
            GeodeticPoint(0.602138592, 2, 156.755)


class TestECEFPoint(unittest.TestCase):
    def test_valid(self):
        x, y, z = 4449038, 3094664, 3351722
        p = ECEFPoint(x, y, z)
        self.assertEqual(p.tup, (x, y, z))
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            ECEFPoint('a', [2], -345)


class TestPoint2D(unittest.TestCase):
    def test_valid(self):
        vals = 1, 2
        p = Point2D(*vals)
        self.assertEqual(p.tup, vals)
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            Point2D('a', 3)


class TestUTMPoint(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()