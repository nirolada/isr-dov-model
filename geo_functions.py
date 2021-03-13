import math

from geo_objects import GeodeticPoint, ECEFPoint


def find_utm_zone(pnt: GeodeticPoint):
    if pnt.lon < math.pi:
        zone = int((math.degrees(pnt.lon) + 180) / 6) + 1
    elif pnt.lon == math.pi:
        zone = 60
    center_lon = -180 + 6*zone - 3  # center meridian longitude
    return (zone, math.radians(center_lon))


def distance_utm(p1, p2):
    if len(p1) == len(p2) == 3:
        return math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)
    else:
        print("Error! These are not 2D points.")
        return None


def distance_ecef(pnt1: ECEFPoint, pnt2: ECEFPoint):
    return math.sqrt((pnt2.x - pnt1.x)**2 + (pnt2.y - pnt1.y)**2 + (pnt2.z - pnt1.z)**2)