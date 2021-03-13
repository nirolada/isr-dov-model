from geo_objects import Ellipsoid


WGS84 = Ellipsoid(name='WGS84', a=6378137.0, f=1/298.257223563)
GRS80 = Ellipsoid(name='GRS80', a=6378137.0, f=1/298.257222101)
WGS72 = Ellipsoid(name='WGS72', a=6378135.0, f=1/298.26)
INT1924 = Ellipsoid(name='International1924', a=6378388.0, f=1/297.0)
CLARKE1880 = Ellipsoid(name='Clarke1880', a=6378300.79, f=1/293.466307656)