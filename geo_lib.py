import math
from validators import validate_input


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


class Ellipsoid:
    def __init__(self, name: str, a: float, f: float) -> None:
        validate_input(self.__init__, name=name, a=a, f=f)
        self.name = name
        self.a = a  # semi major axis
        self.f = f  # flattening
        self.b = self.a * (1 - self.f)  # semi minor axis
        self.e = math.sqrt(2*self.f - self.f**2)  # eccentricity
        self.et = math.sqrt(self.e**2 / (1 - self.e**2))  # second eccentricity
        
    def calc_help_param(self, pnt: GeoPoint) -> None:
        validate_input(self.calc_help_param, pnt=pnt)
        self.eta = self.et * math.cos(pnt.lat)
        self.t = math.tan(pnt.lat)
        self.v = math.sqrt(1 + self.eta**2)
        self.w = math.sqrt(1 - self.e**2 * math.sin(pnt.lat)**2)
        self.rn = self.a / self.w
        self.rm = self.a * (1 - self.e**2) / self.w**3
        self.rg = math.sqrt(self.rn * self.rm)


class CommonEllipsoids:
    WGS84 = Ellipsoid(name='WGS84', a=6378137.0, f=1/298.257223563)
    GRS80 = Ellipsoid(name='GRS80', a=6378137.0, f=1/298.257222101)
    WGS72 = Ellipsoid(name='WGS72', a=6378135.0, f=1/298.26)
    INT1924 = Ellipsoid(name='International1924', a=6378388.0, f=1/297.0)
    CLARKE1880 = Ellipsoid(name='Clarke1880', a=6378300.79, f=1/293.466307656)


def ecef_distance(pnt1: ECEFPoint, pnt2: ECEFPoint):
    validate_input(ecef_distance, pnt1=pnt1, pnt2=pnt2)
    return math.sqrt((pnt2.x - pnt1.x)**2 + (pnt2.y - pnt1.y)**2 + (pnt2.z - pnt1.z)**2)


def geo_to_ecef(pnt: GeoPoint, ell: Ellipsoid) -> tuple:
    validate_input(geo_to_ecef, pnt=pnt, ell=ell)
    ell.calc_help_param(pnt.lat)
    x = (ell.rn + pnt.h) * math.cos(pnt.lat) * math.cos(pnt.lon)
    y = (ell.rn + pnt.h) * math.cos(pnt.lat) * math.sin(pnt.lon)
    z = (ell.rn * (1 - ell.e**2) + pnt.h) * math.sin(pnt.lat)
    return ECEFPoint(x, y, z)


def utm_zone_find(pnt: GeoPoint):
    if pnt.lon < math.pi:
        zone = int((math.degrees(pnt.lon) + 180) / 6) + 1
    elif pnt.lon == math.pi:
        zone = 60
    center_lon = -180 + 6*zone - 3 # center meridian longitude
    return (zone, math.radians(center_lon))


def geo2utm(pnt: GeoPoint, ell: Ellipsoid):
    # need: k0, rn, A, T, C, et, M, M0
    ell.calc_help_param(pnt)
    zone, center_lon = utm_zone_find(pnt)
    
    k0 = 0.9996
    N = ell.rn
    T = ell.t ** 2
    C = ell.et**2 * math.cos(pnt.lat)**2
    A = (pnt.lon - center_lon) * math.cos(pnt.lat)
    M0 = 0
    M = ell.a * (+(1 - ell.e**2/4 - 3*ell.e**4/64 - 5*ell.e**6/256) * pnt.lat 
                 -(3*ell.e**2/8 + 3*ell.e**4/32 + 45*ell.e**6/1024) * math.sin(2 * pnt.lat)
            +(15*e**4/256 +45*e**6/1024)*math.sin(4*pnt.lat)
           -(35*e**6/3072)*math.sin(6*pnt.lat))
    del b, f, v, w, rm, rn, rg, a, eta
    
    x = k0*N*(A 
               +(1 -T +C)*A**3/6
               +(5 -18*T +T**2 +72*C -58*et**2)*A**5/120)
    y = k0*(M - M0 +N*t*( A**2/2 
                         +(5 -T +9*C +4*C**2)*A**4/24
                         +(61 -58*T +T**2 +600*C -330*et**2)*A**6/720))
    east = 500000 + x
    north = y
    return (east, north, zone)
  
    
def distUTM(p1, p2):
    if len(p1) == len(p2) == 3:
        return math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)
    else:
        print("Error! These are not 2D points.")
        return None
