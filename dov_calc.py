#import tkinter as tk
#from tkinter import filedialog
import math
import numpy as np
# import time as t
# import imageio


class Ellipsoid:
    def __init__(self, name: str, a: float, f: float):
        self.name = name
        self.a = a  # semi major axis
        self.f = f  # flattening
        self.b = self.a * (1 - self.f)  # semi minor axis
        self.e = math.sqrt(2*self.f - self.f**2)  # eccentricity
        self.et = math.sqrt(self.e**2 / (1 - self.e**2))  # second eccentricity
        
    def calc_help_param(self, lat_rad):
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


def geo2ecef(lon_rad: float, lat_rad: float, h: float, ellipsoid: Ellipsoid) -> tuple:
    ellipsoid.calc_help_param(lat_rad)
    x = (ellipsoid.rn + h) * math.cos(lat_rad) * math.cos(lon_rad)
    y = (ellipsoid.rn + h) * math.cos(lat_rad) * math.sin(lon_rad)
    z = (ellipsoid.rn * (1 - ellipsoid.e**2) + h) * math.sin(lat_rad)
    return (x, y, z)


def dist3D(p1, p2):
    if len(p1) == len(p2) == 3:
        return math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2 + (p2[2] - p1[2]) ** 2)
    else:
        print("Error! These are not 3D points.")
        return None


def utm_zone_find(lon_rad):
    if lon_rad < math.pi:
        n = int((math.degrees(lon_rad) +180)/6) +1 # zone number
    elif lon_rad == math.pi:
        n = 60
    clon = -180 + n*6 -3 # center meridian longitude
    return (n, math.radians(clon))


def geo2utm(lon_rad, lat_rad, datum_name):
    """
    :param lon_rad:
    :param lat_rad:
    :param datum_name:
    """
    # need: k0, rn, A, T, C, et, M, M0
    (a, b, f, e, et) = ellipsoid_param(datum_name)
    (eta, t, v, w, rn, rm, rg) = help_param(lat_rad, datum_name)
    (n, clon) = utm_zone_find(lon_rad)
    
    k0 = 0.9996
    N = rn
    T = t**2
    C = et**2*math.cos(lat_rad)**2
    A = (lon_rad - clon)*math.cos(lat_rad)
    M0 = 0
    M = a*(+(1 -e**2/4 -3*e**4/64 -5*e**6/256)*lat_rad 
           -(3*e**2/8 +3*e**4/32 +45*e**6/1024)*math.sin(2*lat_rad)
           +(15*e**4/256 +45*e**6/1024)*math.sin(4*lat_rad)
           -(35*e**6/3072)*math.sin(6*lat_rad))
    del b, f, v, w, rm, rn, rg, a, eta
    
    x = k0*N*(A 
               +(1 -T +C)*A**3/6
               +(5 -18*T +T**2 +72*C -58*et**2)*A**5/120)
    y = k0*(M - M0 +N*t*( A**2/2 
                         +(5 -T +9*C +4*C**2)*A**4/24
                         +(61 -58*T +T**2 +600*C -330*et**2)*A**6/720))
    east = 500000 + x
    north = y
    return (east, north, n)
  
    
def distUTM(p1, p2):
    if len(p1) == len(p2) == 3:
        return math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)
    else:
        print("Error! These are not 2D points.")
        return None
    
    
def xyz2grids(xyz_file_path):
    """
    Converts an xyz text file to 3 proper numpy grids (numpy matrix).
    
    Args:
        xyz_file_path (string): file path to the xyz file. 
                                the lines in the file should be in this format - "longitude,latitude,undulation"
        
    Returns:
        lon_grid (numpy matrix): a grid containing the longitude values
        lat_grid (numpy matrix): a grid containing the latitude values
        und_grid (numpy matrix): a grid containing the undulation values
        
    """    
    # open the file the user selected
    xyz_file = open(xyz_file_path) 
    # initiate lists
    lon = [] # longitude list
    lat = [] # latitude list
    und = [] # undulation list
    # fill the lists with values
    while True:
        line = xyz_file.readline()
        if not line:
            break
        line_vals = line.split(',')
        lon.append(float(line_vals[0]))
        lat.append(float(line_vals[1]))
        und.append(float(line_vals[2]))
    del line, line_vals
    xyz_file.close()
    # convert lists to numpy arrays
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    und = np.asarray(und)
    # calculate the output grids dimensions
    n = np.unique(lon).size
    m = np.unique(lat).size
    # reshape the numpy arrays to grids
    lon_grid = np.reshape(lon, (n,m))
    lon_grid = np.transpose(lon_grid)
    np.save("lon_grid.npy", lon_grid)
    
    lat_grid = np.reshape(lat, (n,m))
    lat_grid = np.transpose(lat_grid)
    np.save("lat_grid.npy", lat_grid)
    
    und_grid = np.reshape(und, (n,m))
    und_grid = np.transpose(und_grid)
    np.save("und_grid.npy", und_grid)
    
    return None


def dov_grids(lon_grid, lat_grid, und_grid):
    """
    
    """
    (m, n) = und_grid.shape
    
    xi_grid = np.empty((m - 2, n - 2), dtype=float)
    eta_grid = np.empty((m - 2, n - 2), dtype=float)
    
    print('DOV grids processing status:')
    for i in range(1, m - 1):    
        for j in range(1, n - 1):        
            # calc of the DOV in the North direction, Xi
            p_here = geo2utm(math.radians(lon_grid[i, j]), 
                              math.radians(lat_grid[i, j]), 
                              "WGS84")
            p_north = geo2utm(math.radians(lon_grid[i - 1, j]), 
                               math.radians(lat_grid[i - 1, j]), 
                               "WGS84")

            d_north = distUTM(p_here, p_north)
            d_und_north = und_grid[i - 1, j] - und_grid[i, j]           
            xi_grid[i - 1, j - 1] = math.degrees(-d_und_north / d_north) * 3600
    
            # calc of the DOV in the East direction, Eta
            p_east = geo2utm(math.radians(lon_grid[i, j + 1]), 
                              math.radians(lat_grid[i, j + 1]), 
                              "WGS84")
            
            d_east = distUTM(p_here, p_east)
            d_und_east = und_grid[i, j + 1] - und_grid[i, j]       
            eta_grid[i - 1, j - 1] = math.degrees(-d_und_east / d_east) * 3600

        if i % 409 == 0:
            print(int(i/(m-2)*100), end='%.. ')
    print('100%!\nDOV grids ready! saved them as xi_grid.npy and eta_grid.npy')
    np.save("xi_grid.npy", xi_grid)
    np.save("eta_grid.npy", eta_grid)


def grid_bilinear_interp(x_mesh_grid, y_mesh_grid, z_grid, points):
    """
    Calculate the z value of input points inside the grid using bilinear interpolation
    
    Args:
        x_mesh_grid (numpy 2d array): a mesh grid of the x-longitude values (top = high, bot = low)
        y_mesh_grid (numpy 2d array): a mesh grid of the y-latitude values (right = high, left = low)
        z_grid (numpy 2d array): a grid of the z values
        points (list) of (tuple): a list of points inside the grid for which we want
                                  to calculate z value.
                                  required format: [(x1,y1), (x2, y2), ...(xn, yn)]
    
    Returns:
        interp_z (list): a list of the interpolated z values for the input points
    """
    # initiate returned list
    interp_z = []
    # get input grid size
    (m, n) = z_grid.shape
    # interpolate z val for input points
    for k in range(len(points)):
        (x, y) = (points[k][0], points[k][1])
        # find the 4 surrounding grid points
        for i in range(m):
            if(y < y_mesh_grid[i, 0]):
                continue
            else:
                y1 = y_mesh_grid[i, 0]
                y2 = y_mesh_grid[i-1, 0]
                break
        for j in range(n):        
            if(x > x_mesh_grid[0, j]):
                continue
            else:
                x1 = x_mesh_grid[0, j-1]
                x2 = x_mesh_grid[0, j]
                break
        
        q11 = z_grid[i, j-1] # x1, y1
        q12 = z_grid[i-1, j-1] # x1, y2
        q21 = z_grid[i, j] # x2, y1
        q22 = z_grid[i-1, j] # x2, y2
        
        fxy1 = (x2-x)/(x2-x1)*q11 + (x-x1)/(x2-x1)*q21
        fxy2 = (x2-x)/(x2-x1)*q12 + (x-x1)/(x2-x1)*q22
        fxy = (y2-y)/(y2-y1)*fxy1 + (y-y1)/(y2-y1)*fxy2
        # print("x1:", x1, "x:", x, "x2:", x2)
        # print("y1:", y1, "y:", y, "y2:", y2)
        # print("Q11:", q11, "Q12:", q12, "Q21:", q21, "Q22:", q22)
        # print("interpolated z:", fxy,'\n')
        interp_z.append(fxy)
    
    return interp_z

# ask user to select the xyz file
# root = tk.Tk()
# xyz_file_path = filedialog.askopenfilename(filetypes = [("XYZ File", ".xyz")])
# root.withdraw()
# create grids and save them as npy files
# xyz2grids(xyz_file_path)


# display grids
# und_grid_norm = (und_grid - und_grid.min())/(und_grid.max() - und_grid.min())
# und_grid_norm = 255*und_grid_norm
# und_grid_uint8 = und_grid_norm.astype(np.uint8)
# imageio.imwrite('und.tiff',und_grid_uint8)

# load grids
lon_grid = np.load("lon_grid.npy")
lat_grid = np.load("lat_grid.npy")
und_grid = np.load("und_grid.npy")

# create DOV grids
# dov_grids(lon_grid, lat_grid, und_grid)


# load DOV grids
xi_grid = np.load("xi_grid.npy")
eta_grid = np.load("eta_grid.npy")

# load Laplas stations data - laplace points geographic and astronomic position
points_geo_pos = []
points_astro_pos = []
points_id = []
with open("laplace_stations.csv", 'r') as f:
    while True:
        line = f.readline()
        if not line:
            break
        vals = line.split(',')
        points_id.append(vals[0])
        
        lon_geo = float(vals[1]) + float(vals[2]) / 60 + float(vals[3]) / 3600
        lat_geo = float(vals[4]) + float(vals[5]) / 60 + float(vals[6]) / 3600
        points_geo_pos.append([lon_geo, lat_geo])
        
        lon_astro = float(vals[7]) + float(vals[8]) / 60 + float(vals[9]) / 3600
        lat_astro = float(vals[10]) + float(vals[11]) / 60 + float(vals[12]) / 3600
        points_astro_pos.append([lon_astro, lat_astro])
    del line, vals, lon_geo, lat_geo, lon_astro, lat_astro

# calculate xi and eta by astronomic and geodetic coordinates
# xi = astro_lat - geo_lat
# eta = (astro_lon - geo_lon) * cos(astro_lat)    
points_laplace_xi = []
points_laplace_eta = []
for k in range(len(points_geo_pos)):
    astro_lon = points_astro_pos[k][0]
    geo_lon = points_geo_pos[k][0]
    astro_lat = points_astro_pos[k][1]
    geo_lat = points_geo_pos[k][1]
    xi = (astro_lat - geo_lat) * 3600
    eta = (astro_lon - geo_lon) * math.cos(math.radians(astro_lat)) * 3600
    points_laplace_xi.append(xi)
    points_laplace_eta.append(eta)
del astro_lon, astro_lat, geo_lon, geo_lat, xi, eta

# calculate xi and eta by interpolation on our generated xi and eta grids
points_model_xi = grid_bilinear_interp(lon_grid[1 : -1, 1 : -1], 
                                       lat_grid[1 : -1, 1 : -1], 
                                       xi_grid, 
                                       points_geo_pos)

points_model_eta = grid_bilinear_interp(lon_grid[1 : -1, 1 : -1], 
                                        lat_grid[1 : -1, 1 : -1], 
                                        eta_grid, 
                                        points_geo_pos)

# make a comparison csv file to show the differences between our 
# model xi and eta and the calculated xi and eta
with open("model_vs_laplace.csv", 'w') as f:
    header = ','.join(['id', 'lon', 'lat', 'model xi', 'model eta', 'laplace xi', 'laplace eta', 'xi diff', 'eta diff']) + '\n'
    f.write(header)
    for k in range(len(points_id)):
        pid = points_id[k]
        lon = str(round(points_geo_pos[k][0], 8))
        lat = str(round(points_geo_pos[k][1], 8))
        laplace_xi = str(round(points_laplace_xi[k], 1))
        laplace_eta = str(round(points_laplace_eta[k], 1))
        model_xi = str(round(points_model_xi[k], 1))
        model_eta = str(round(points_model_eta[k], 1))
        diff_xi = str(round(points_model_xi[k] - points_laplace_xi[k], 1))
        diff_eta = str(round(points_model_eta[k] - points_laplace_eta[k], 1))
        line = ','.join([pid, lon, lat, model_xi, model_eta, laplace_xi, laplace_eta, diff_xi, diff_eta]) + '\n'
        f.write(line)
