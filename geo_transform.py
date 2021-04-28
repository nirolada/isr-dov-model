import math

from geo_objects import GeodeticPoint, ECEFPoint, UTMPoint, WGS84Ellipsoid
from geo_functions import find_utm_zone


def geo_to_ecef(p: GeodeticPoint, ell: WGS84Ellipsoid) -> ECEFPoint:
    _, rn = ell.calc_help_param(p.lat)
    x = (rn + p.h) * math.cos(p.lat) * math.cos(p.lon)
    y = (rn + p.h) * math.cos(p.lat) * math.sin(p.lon)
    z = (rn * (1 - ell.e**2) + p.h) * math.sin(p.lat)
    return ECEFPoint(x, y, z)


def geo_to_utm(p: GeodeticPoint, ell: WGS84Ellipsoid):
    t, rn = ell.calc_help_param(p.lat)
    zone, center_lon = find_utm_zone(p.lon)

    k0 = 0.9996
    N = rn
    T = t ** 2
    C = ell.et**2 * math.cos(p.lat)**2
    A = (p.lon - center_lon) * math.cos(p.lat)
    M0 = 0
    M = ell.a * (
        + (1 - ell.e**2/4 - 3*ell.e**4/64 - 5*ell.e**6/256) * p.lat
        - (3*ell.e**2/8 + 3*ell.e**4/32 + 45 *
           ell.e**6/1024) * math.sin(2 * p.lat)
        + (15*ell.e**4/256 + 45*ell.e**6/1024) * math.sin(4 * p.lat)
        - (35 * ell.e**6/3072) * math.sin(6 * p.lat))

    x = k0 * N * (
        A + (1 - T + C) * A**3/6 + (5 - 18*T + T**2 + 72*C - 58*ell.et**2) * A**5/120)
    y = k0 * (
        M - M0 + N*t*(
            + A**2/2
            + (5 - T + 9*C + 4*C**2) * A**4/24
            + (61 - 58*T + T**2 + 600*C - 330*ell.et**2) * A**6/720))
    east = 500000 + x
    north = y
    return UTMPoint(zone, east, north)