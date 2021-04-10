import math
from typing import Optional, Tuple, NewType

from geo_objects import Number, GeodeticPoint, ECEFPoint, UTMPoint


def find_utm_zone(lon_rad: Number) -> Tuple[int, Number]:
    if -math.pi <= lon_rad < math.pi:
        zone = int((math.degrees(lon_rad) + 180) / 6) + 1
    elif lon_rad == math.pi:
        zone = 60
    else:
        raise ValueError('Longitude in radians must be between [-π, π]!')
    
    center_lon = -180 + 6*zone - 3  # center meridian longitude
    return zone, math.radians(center_lon)


def distance_utm(p1: UTMPoint, p2: UTMPoint) -> Optional[Number]:
    if len(p1) == len(p2) == 3:
        return math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)
    else:
        print("Error! These are not 2D points.")
        return None
