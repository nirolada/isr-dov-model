import math

from geo_objects import GeodeticPoint, ECEFPoint, UTMPoint, Ellipsoid
from geo_functions import find_utm_zone


def geo_to_ecef(pnt: GeodeticPoint, ell: Ellipsoid) -> ECEFPoint:
    ell.calc_help_param(pnt.lat)
    x = (ell.rn + pnt.h) * math.cos(pnt.lat) * math.cos(pnt.lon)
    y = (ell.rn + pnt.h) * math.cos(pnt.lat) * math.sin(pnt.lon)
    z = (ell.rn * (1 - ell.e**2) + pnt.h) * math.sin(pnt.lat)
    return ECEFPoint(x, y, z)


def geo_to_utm(pnt: GeodeticPoint, ell: Ellipsoid):
    ell.calc_help_param(pnt)
    zone, center_lon = find_utm_zone(pnt)

    k0 = 0.9996
    N = ell.rn
    T = ell.t ** 2
    C = ell.et**2 * math.cos(pnt.lat)**2
    A = (pnt.lon - center_lon) * math.cos(pnt.lat)
    M0 = 0
    M = ell.a * (
        + (1 - ell.e**2/4 - 3*ell.e**4/64 - 5*ell.e**6/256) * pnt.lat
        - (3*ell.e**2/8 + 3*ell.e**4/32 + 45 *
           ell.e**6/1024) * math.sin(2 * pnt.lat)
        + (15*ell.e**4/256 + 45*ell.e**6/1024) * math.sin(4 * pnt.lat)
        - (35 * ell.e**6/3072) * math.sin(6 * pnt.lat))

    x = k0 * N * (
        A + (1 - T + C) * A**3/6 + (5 - 18*T + T**2 + 72*C - 58*ell.et**2) * A**5/120)
    y = k0 * (
        M - M0 + N*ell.t*(
            + A**2/2
            + (5 - T + 9*C + 4*C**2) * A**4/24
            + (61 - 58*T + T**2 + 600*C - 330*ell.et**2) * A**6/720))
    east = 500000 + x
    north = y
    return UTMPoint(east, north, zone)