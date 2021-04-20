import unittest
from collections import namedtuple

from geo_objects import (
    ECEFPoint, Point2D, Point3D, GeodeticPoint, UTMPoint
)


class TestPoint3D(unittest.TestCase):
    def test_valid(self):
        valid = 1, 2.5, 3
        p = Point3D(*valid)
        self.assertTupleEqual(p.tup, valid)
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            Point3D((1,), '2', 3)


class TestGeodeticPoint(unittest.TestCase):
    Point = namedtuple('Point', ['lon', 'lat', 'h'])
    valid = Point(0.602138592, 0.593411946, 156.755)
    
    def test_valid(self):
        p = GeodeticPoint(*self.valid)
        self.assertEqual(p.tup, self.valid)
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            GeodeticPoint('a', (3,), self.valid.h)
            
    def test_invalid_lon(self):
        with self.assertRaises(ValueError):
            GeodeticPoint(-4, self.valid.lat, self.valid.h)
    
    def test_invalid_lat(self):
        with self.assertRaises(ValueError):
            GeodeticPoint(self.valid.lon, 2, self.valid.h)


class TestECEFPoint(unittest.TestCase):
    Point = namedtuple('Point', ['x', 'y', 'z'])
    valid = Point(4449038.5, 3094664, 3351722)
    
    def test_valid(self):
        p = ECEFPoint(*self.valid)
        self.assertEqual(p.tup, self.valid)
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            ECEFPoint('a', [2], self.valid.z)


class TestPoint2D(unittest.TestCase):
    def test_valid(self):
        vals = 1, 2
        p = Point2D(*vals)
        self.assertEqual(p.tup, vals)
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            Point2D('a', 3)


class TestUTMPoint(unittest.TestCase):
    Point = namedtuple('Point', ['zone', 'east', 'north'])
    valid = Point(36, 683383.341, 3630033.034)
    
    def test_valid(self):
        p = UTMPoint(*self.valid)
        self.assertEqual(p.tup, self.valid)
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            UTMPoint('a', (3,), self.valid.north)
    
    def test_invalid_zone(self):
        with self.assertRaises(ValueError):
            UTMPoint(67, self.valid.east, self.valid.north)
    
    def test_invalid_east(self):
        with self.assertRaises(ValueError):
            UTMPoint(self.valid.zone, 100000, self.valid.north)
    
    def test_invalid_north(self):
        with self.assertRaises(ValueError):
            UTMPoint(self.valid.zone, self.valid.east, 12000000)


if __name__ == '__main__':
    unittest.main()