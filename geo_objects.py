import math
from validators import validate_input


class Ellipsoid:
    def __init__(self, name: str, a: float, f: float) -> None:
        validate_input(self.__init__, name=name, a=a, f=f)
        self.name = name
        self.a = a  # semi major axis
        self.f = f  # flattening
        self.b = self.a * (1 - self.f)  # semi minor axis
        self.e = math.sqrt(2*self.f - self.f**2)  # eccentricity
        self.et = math.sqrt(self.e**2 / (1 - self.e**2))  # second eccentricity
        
    def calc_help_param(self, lat_rad: float) -> None:
        validate_input(self.calc_help_param, lat_rad=lat_rad)
        self.eta = self.et * math.cos(lat_rad)
        self.t = math.tan(lat_rad)
        self.v = math.sqrt(1 + self.eta**2)
        self.w = math.sqrt(1 - self.e**2 * math.sin(lat_rad)**2)
        self.rn = self.a / self.w
        self.rm = self.a * (1 - self.e**2) / self.w**3
        self.rg = math.sqrt(self.rn * self.rm)


class CommonEllipsoids:
    WGS84 = Ellipsoid(name='WGS84', a=6378137.0, f=1/298.257223563)
    GRS80 = Ellipsoid(name='GRS80', a=6378137.0, f=1/298.257222101)
    WGS72 = Ellipsoid(name='WGS72', a=6378135.0, f=1/298.26)
    INT1924 = Ellipsoid(name='International1924', a=6378388.0, f=1/297.0)
    CLARKE1880 = Ellipsoid(name='Clarke1880', a=6378300.79, f=1/293.466307656)


class GeoPoint:
    def __init__(self, lon_rad: float, lat_rad: float, h: float) -> None:
        validate_input(self.__init__, lon_rad=lon_rad, lat_rad=lat_rad, h=h)
        if lon_rad < -math.pi or lon_rad > math.pi:
            raise ValueError('Longitude in radians must be between [-π, π]!')
        if lat_rad < -math.pi / 2 or lat_rad > math.pi / 2:
            raise ValueError('Latitude in radians must be between [-π/2, π/2]!')
        self.lon, self.lat, self.h = lon_rad, lat_rad, h


class ECEFPoint:
    def __init__(self, x: float, y: float, z: float) -> None:
        validate_input(self.__init__, x=x, y=y, z=z)
        self.x, self.y, self.z = x, y, z
    
    def calc_distance(self, other) -> float:
        if not isinstance(other, type(self)):
            raise TypeError('Parameter must be an instance of ECEFPoint as well!')
        return math.sqrt()
        