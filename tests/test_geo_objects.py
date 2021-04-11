import unittest
from geo_objects import (
    Point3D, GeodeticPoint
)


class TestPoint3D(unittest.TestCase):
    def test_valid(self):
        vals = (1, 2.5, 3)
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
        lon, lat, h = 'a', (3,), 360.5
        with self.assertRaises(TypeError):
            GeodeticPoint(lon, lat, h)
    def test_invalid_lon(self):
        lon, lat, h = -4, 0.593411946, 156.755
        with self.assertRaises(ValueError):
            GeodeticPoint(lon, lat, h)
    
    def test_invalid_lat(self):
        lon, lat, h = 0.602138592, 2, 156.755
        with self.assertRaises(ValueError):
            GeodeticPoint(lon, lat, h)


class TestECEFPoint(unittest.TestCase):
    def test_valid(self):
        pass


if __name__ == '__main__':
    unittest.main()