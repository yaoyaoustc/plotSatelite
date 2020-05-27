import plotfigure

bbox = [144, 149, -39, -36]
lat_grid, lon_grid = 0.02, 0.02

filename = 'MelbourneNO220200408-20200408.txt'
shapefilename = './shape/vic_lg_roadmap_shape.shp'

# read in txt data
data = plotfigure.read_average(filename,bbox,lat_grid,lon_grid)

# plot and save the fig as png
plotfigure.plot_figure(data,shapefilename,bbox,lat_grid,lon_grid,filename)
