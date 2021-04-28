import unittest

from geo_transform import geo_to_ecef, geo_to_utm
from geo_objects import (
    ECEFPoint, GeodeticPoint, UTMPoint, wgs84
)


class TestGeoToECEF(unittest.TestCase):
    geo_point = GeodeticPoint(0.602138592, 0.593411946, 156.755)
    ecef_point = ECEFPoint(4362419.917, 2998208.142, 3546534.220)
    
    def test_valid(self):
        converted_point = geo_to_ecef(self.geo_point, wgs84)
        self.assertAlmostEqual(converted_point.x, self.ecef_point.x, places=2)
        self.assertAlmostEqual(converted_point.y, self.ecef_point.y, places=2)
        self.assertAlmostEqual(converted_point.z, self.ecef_point.z, places=2)


class TestGeoToUTM(unittest.TestCase):
    geo_point = GeodeticPoint(0.602138592, 0.593411946, 156.755)
    utm_point = UTMPoint(36, 638527.725, 3763170.114)
    
    def test_valid(self):
        converted_point = geo_to_utm(self.geo_point, wgs84)
        self.assertEqual(converted_point.zone, self.utm_point.zone)
        self.assertAlmostEqual(converted_point.east, self.utm_point.east, places=2)
        self.assertAlmostEqual(converted_point.north, self.utm_point.north, places=2)