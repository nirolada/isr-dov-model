
#import tkinter as tk
#from tkinter import filedialog
import math
import numpy as np
# import time as t
# import imageio
from geo_objects import Ellipsoid, GeoPoint, ECEFPoint
    
    
def xyz2grids(xyz_file_path):
    """
    Converts an xyz file to 3 proper numpy grids.
    Saves the corresponding grids as '.npy' files.
    
    @type xyz_file_path: str
    @param xyz_file_path: A path to a xyz file. 
                          Expected line format: <longitude>,<latitude>,<undulation>
    """    
    xyz_file = open(xyz_file_path) 
    lon, lat, und = list(), list(), list()
    # TO BE CONTINUED FROM HERE!!!
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
