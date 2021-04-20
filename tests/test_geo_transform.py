import unittest
# from collections import namedtuple

from geo_transform import geo_to_ecef, geo_to_utm
from geo_objects import (
    ECEFPoint, GeodeticPoint, UTMPoint, WGS84Ellipsoid
)


class TestGeoToECEF(unittest.TestCase):
    ellipsoid = WGS84Ellipsoid()
    geo_point = GeodeticPoint(0.602138592, 0.593411946, 156.755)
    ecef_point = None
    
    def test_valid(self):
        self.assertTupleEqual(p.tup, valid)
    
    def test_invalid(self):
        with self.assertRaises(TypeError):
            Point3D((1,), '2', 3)