import geopandas as gpd
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.ticker as mticker
import numpy as np

def read_average(filename,bbox,lat_grid,lon_grid):
    [w, e, s, n] = bbox
    lonsize = int((e-w)/lon_grid)
    latsize = int((n-s)/lat_grid)
    average_z = np.zeros((latsize, lonsize))
    lines = []
    with open(filename,'r') as f:
        for line in f:
            lines.append(line)
        counter = 0
        for i in range(latsize):
            for j in range(lonsize):
                average_z[i,j] = float(lines[counter])
                counter += 1
        f.close()
    return average_z

def plot_figure(inputdata,shapefile,bbox,lat_grid,lon_grid,figtitle):
    
    #create the grid
    [w, e, s, n] = bbox
    xx, yy = np.meshgrid(np.linspace(w, e, int((e-w)/lon_grid)), np.linspace(s, n, int((n-s)/lat_grid)))

    # load shape file
    shape_gpd = gpd.read_file(shapefile)

    # initilize the plot, plot coastline
    fig, ax = plt.subplots(figsize = (15,15),subplot_kw={'projection': ccrs.PlateCarree()})
    ax.coastlines(resolution='10m')
    ax.set_extent(bbox, crs=ccrs.PlateCarree())
    
    # plot average data
    im = ax.pcolormesh(xx, yy, np.squeeze(inputdata), cmap = 'RdYlBu_r', vmin=0, vmax = 80,transform=ccrs.PlateCarree())

    #plot highway on top of it
    shape_gpd.plot(ax = ax, color = 'black')

    #refine plot properties.
    
    cb = plt.colorbar(im, pad = 0.1, orientation = 'horizontal', shrink = 0.6)
    
    z = 'NO2'

    if z == 'CH4':
       units = 'ppb'
    if z == 'NO2':
       units = 'μmol m-2'
    if z == 'SO2':
       units = 'mol m-2'
    if z == 'CO':
       units = 'mol m-2'
    cb.set_label(units)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels = True, color = 'gray', linestyle='--')
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    # gl.xlocator = mticker.FixedLocator([xmin, xmin + lon_dif, xmin + (lon_dif * 2), xmin + (lon_dif * 3), xmin + (lon_dif * 4), xmin + (lon_dif * 5)])
    # gl.ylocator = mticker.FixedLocator([ymax, ymax - lat_dif, ymax - (lat_dif * 2), ymax - (lat_dif * 3), ymax - (lat_dif * 4), ymax - (lat_dif * 5)])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    plt.savefig(figtitle + '.png')
