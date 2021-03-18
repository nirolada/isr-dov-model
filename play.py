import geo_objects
import math

ep1 = geo_objects.ECEFPoint(1000, 1000, 1000)
ep2 = geo_objects.ECEFPoint(2000, 2000, 2000)
print(math.dist(ep1, ep2))
print(math.sqrt((ep2.x_m - ep1.x_m)**2 + (ep2.y_m - ep1.y_m)**2 + (ep2.z_m - ep1.z_m)**2))
