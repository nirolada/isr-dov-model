from collections import namedtuple
import math
from typing import TypeVar, NamedTuple

Number = TypeVar('Number', int, float)


class GeodeticPoint:
    def __init__(self, lon_rad: Number, lat_rad: Number, h_m: Number) -> None:
        if lon_rad < -math.pi or lon_rad > math.pi:
            raise ValueError('Longitude in radians must be between [-π, π]!')
        if lat_rad < -math.pi / 2 or lat_rad > math.pi / 2:
            raise ValueError(
                'Latitude in radians must be between [-π/2, π/2]!')
        self.lon, self.lat, self.h = lon_rad, lat_rad, h_m

class ECEFPoint(NamedTuple):
    # _m - in meters
    x_m: Number
    y_m: Number
    z_m: Number


class UTMPoint:
    def __init__(self):
        pass


class Ellipsoid:
    def __init__(self, name: str, a_m: Number, f: Number) -> None:
        self.name = name
        self.a = a_m  # semi major axis in meters
        self.f = f  # flattening
        self.e = math.sqrt(2*f - f**2)  # eccentricity
        self.et = math.sqrt(self.e**2 / (1 - self.e**2))  # second eccentricity

    def calc_help_param(self, lat_rad: Number) -> None:
        self.t = math.tan(lat_rad)
        
        # radius of curvature in the prime vertical
        self.rn = self.a / math.sqrt(1 - self.e**2 * math.sin(lat_rad)**2)
