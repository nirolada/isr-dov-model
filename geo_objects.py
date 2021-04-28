import math
from typing import Tuple, TypeVar

Number = TypeVar('Number', int, float)


class Point3D:
    def __init__(self, val1: Number, val2: Number, val3: Number) -> None:
        if not all(map(lambda v: isinstance(v, (int, float)), (val1, val2, val3))):
            raise TypeError('All values must be integers or floats!')
        self._val1, self._val2, self._val3 = val1, val2, val3
    
    @property
    def tup(self) -> Tuple[Number, Number, Number]:
        return self._val1, self._val2, self._val3


class GeodeticPoint(Point3D):
    def __init__(self, lon_rad: Number, lat_rad: Number, h_m: Number) -> None:
        # _rad - radians  # _m - meters
        if lon_rad < -math.pi or lon_rad > math.pi:
            raise ValueError('Longitude in radians must be between [-π, π]!')
        if lat_rad < -math.pi / 2 or lat_rad > math.pi / 2:
            raise ValueError(
                'Latitude in radians must be between [-π/2, π/2]!')
        super().__init__(lon_rad, lat_rad, h_m)

    @property
    def lon(self) -> Number:
        return self._val1

    @property
    def lat(self) -> Number:
        return self._val2

    @property
    def h(self) -> Number:
        return self._val3


class ECEFPoint(Point3D):
    def __init__(self, x_m: Number, y_m: Number, z_m: Number) -> None:
        # _m - meter
        super().__init__(x_m, y_m, z_m)
    
    @property
    def x(self) -> Number:
        return self._val1

    @property
    def y(self)-> Number:
        return self._val2

    @property
    def z(self)-> Number:
        return self._val3


class Point2D:
    def __init__(self, val1: Number, val2: Number) -> None:
        if not all(map(lambda v: isinstance(v, (int, float)), (val1, val2))):
            raise TypeError('All values must be integers or floats!')
        self._val1, self._val2 = val1, val2
    
    @property
    def tup(self)-> Tuple[Number, Number]:
        return self._val1, self._val2


class UTMPoint(Point2D):
    def __init__(self, zone: int, east_m: Number, north_m: Number):
        # _m - meter
        if not isinstance(zone, int):
            raise TypeError
        
        super().__init__(east_m, north_m)
        
        if not 1 <= zone <= 60:
            raise ValueError('Zone number must be between [1, 60]!')
        if not 200000 <= east_m <= 800000:
            raise ValueError('East value must be between [200K, 800K] meters!')
        if not 0 <= north_m <= 10000000:
            raise ValueError('North value must be between [0, 10M] meters!')
        
        self._zone = zone
        

    @property
    def east(self) -> Number:
        return self._val1

    @property
    def north(self) -> Number:
        return self._val2

    @property
    def zone(self) -> int:
        return self._zone

    @property
    def tup(self) -> Tuple[int, Number, Number]:
        return (self._zone,) + super().tup


class WGS84Ellipsoid:
    def __init__(self) -> None:
        self._name = 'WGS84'
        self._a = 6378137  # semi major axis in meters
        self._f = 1/298.257223563  # flattening
        self._e = math.sqrt(2*self._f - self._f**2)  # eccentricity
        self._et = math.sqrt(self._e**2 / (1 - self._e**2))  # 2nd eccentricity

    # Calculate help parameters for a given latitude in radians.
    def calc_help_param(self, lat_rad: Number):
        t: Number = math.tan(lat_rad)
        
        # radius of curvature in the prime vertical
        rn: Number =  self._a / math.sqrt(1 - self._e**2 * math.sin(lat_rad)**2)
        return t, rn                      

    @property
    def name(self) -> str:
        return self._name

    @property
    def a(self) -> Number:
        return self._a

    @property
    def f(self) -> Number:
        return self._f

    @property
    def e(self) -> Number:
        return self._e

    @property
    def et(self) -> Number:
        return self._et
    
    
wgs84 = WGS84Ellipsoid()
