#
#this routine reads in Tropomi data and avrages the data for defined time periods the orginal data is approx 4x7km but this routine defines a new grid that is lower resolution and remaps the data. A number of variables need to be defined at the top of the file for example the region to be plotted and input and output directories the lat/lon grid onto which the data will be averaged
#

# Import All Necessary Programs
import matplotlib
matplotlib.use('Agg')
import netCDF4 as nc4
import math
from netCDF4 import Dataset
import numpy as np
import os
import os.path
from os import path
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs
import cartopy
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import sys

# If wishing to plot, leave as True, otherwise switch to False
plot_on = True

# Change these dates to those you wish to focus on this is done outside the script

start_year = int(sys.argv[1])
start_month = int(sys.argv[2])
start_day = int(sys.argv[3])
end_year = int(sys.argv[4])
end_month = int(sys.argv[5])
end_day = int(sys.argv[6])

# Enter here the path to the netCDF files
path_ = '/scratch/k10/cp2786/PRODUCT/'
print(start_year,end_year,start_month,end_month,start_day,end_day)

# Simply uncomment below whichever gas you would like to look at 
#z = 'CH4'
#z = 'NO2'
#z = 'CO'
#z = 'SO2'
z = sys.argv[7]

# Set the lonmin lonmas latmin latmax boundary for the plot 
xmin, xmax, ymin, ymax = 0, 0, 0, 0


#Melbourne
xmin, xmax, ymin, ymax = 144, 149, -39, -36
city='Melbourne '

#Victoria and NSW region
#xmin, xmax, ymin, ymax = 144, 155, -40, -30
#city='Victoria and NSW fires '


print('city',city)
# Set your image title and what you want to name the output file, as well as the destination of the output file 
# E.g. 'Methane in Australia June 2018 - June 2019', 'aus_methane_june_18_june_19' and 'scratch/k10/lp4439/'
title = city + z + ' '+str(start_year)+str(start_month).zfill(2)+ str(start_day).zfill(2)+'-'+str(end_year)+str(end_month).zfill(2)+ str(end_day).zfill(2)

output_file=title.replace(' ','')
output_destination = '/scratch/k10/cp2786/'
# Set your min and max for the colour bar
v_min = 0
v_max = 0
if z == 'CH4':
    v_min = 1770
    v_max = 1850
elif z == 'NO2':
    v_min = 0
    v_max = 80
#   v_max = 150
elif z == 'CO':
    v_min = 0
    v_max = 0.25 #changed from .35 better for monthly
    
#This creates a fine grid over the box created by the points
# Lat and lon length is how wide each pixel will be in the final grid, i.e. 0.1 = 0.1 of a degree
#approx 10km grid
lat_length = 0.1
lon_length = 0.1
#approx 2km grid
lat_length = 0.02
lon_length = 0.02

#create the grid
xx, yy = np.meshgrid(np.linspace(xmin, xmax, (xmax-xmin)/lon_length), np.linspace(ymin, ymax, (ymax-ymin)/lat_length))

num_rows = len(xx)
num_cols = len(xx[0])

#Create a zero array where the sum of z will be stored
average_z = np.zeros((num_rows, num_cols))
#Create a zero array which will store the number of times a measurement is added to a cell in average_z
counter = np.zeros((num_rows, num_cols))


# Set up some quality control if looking at Methane set to true
# Set alebdo to False, as CH4 is the only object that looks at albedo 
albedo = False


if z == 'CH4':
    filename_check = 'CH4'
elif z == 'NO2':
    filename_check = 'NO2'
elif z == 'CO':
    filename_check = 'CO_'
elif z == 'SO2':
    filename_check = 'SO2'

#this loop seraches through the files
    
if not path.exists(output_destination + output_file + '.txt'):
    #Loop through the files
    for filename in os.listdir(path_):
        #Read files one at a time
        # Split the filename from the extension so can check if netcdf or not
        file_name, file_extension = os.path.splitext(filename)

        if file_extension == '.nc' and filename[13:16] == filename_check:
        
            year = int(filename[20:24])
            month = int(filename[24:26])
            day = int(filename[26:28])

            # Check if in date range 
            if dt.date(start_year, start_month, start_day) <= dt.date(year, month, day) <= dt.date(end_year, end_month, end_day):

                # Read the file 
                sat_file = path_  + filename
                print('sat_file',sat_file)
                fh = Dataset(sat_file, mode='r')
            
                # Read the data depending on what z is
                if z == 'CH4':
            
                    z_main = fh.groups['PRODUCT'].variables['methane_mixing_ratio_bias_corrected'][0,:,:]
                    albedo_SWIR = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].variables['surface_albedo_SWIR'][0,:,:]
                    albedo_NIR = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].variables['surface_albedo_NIR'][0,:,:]
                    albedo = True 
            
                elif z == 'NO2':
                    z_main = fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric_column'][0,:,:] * 1000000
                
                elif z == 'CO':
                    z_main = fh.groups['PRODUCT'].variables['carbonmonoxide_total_column'][0,:,:]
                
                elif z == 'SO2':
                    z_main = fh.groups['PRODUCT'].variables['sulfurdioxide_total_vertical_column'][0,:,:]

                # Read in the longitude, latitude and qa values 
                lons = fh.groups['PRODUCT'].variables['longitude'][:][0,:,:]
                lats = fh.groups['PRODUCT'].variables['latitude'][:][0,:,:]
                qa = fh.groups['PRODUCT'].variables['qa_value'][0,:,:]
            
                #Work out the satellite resolution by finding the difference between two adjacent points  this is a bit experimental and could be a better way of doing this
                #Produces the best plot when latitude and longitude resolution are the same 
                #lat_sat_res = (abs(lats[0,0])-abs(lats[0,1])) 
                #lon_sat_res = lat_sat_res 
                #lat_half_res = abs(lat_sat_res)/2
                #lon_half_res = abs(lon_sat_res)/2
                #lat_sat_res = (abs(lats[0,0])-abs(lats[0,1])) #orig
                lat_sat_res = (abs(lats[0,0])-abs(lats[1,0]))
                
                lon_sat_res = (abs(lons[0,0])-abs(lons[0,1]))
                lat_half_res = abs(lat_sat_res)/2
                lon_half_res = abs(lon_sat_res)/2
                
                # Find the number of rows and columns in the original CH4 data            
                sat_rows = len(z_main)
                sat_cols = len(z_main[0])
                
                #Create a mask, that will mask out the points that are outside the area we are looking at 
                satmask = np.full((sat_rows, sat_cols), True, dtype=bool)
                satmask = (lons >  (xmin-lon_sat_res-2))  & (lons < (xmax+lon_sat_res+1)) & (lats > (ymin-lat_sat_res-1)) & (lats <(ymax+lat_sat_res+1)) 
            
                # Apply the mask to the data 
                short_z_main = z_main[satmask]
                new_qa = qa[satmask]
                sat_len = len(short_z_main)
                new_lat = lats[satmask]
                new_lon = lons[satmask]
                if albedo is True:
                    new_SWIR = albedo_SWIR[satmask]
                    new_NIR = albedo_NIR[satmask]

                for b in range(sat_len):
                    if new_qa[b] < 0.5:
                        short_z_main[b] = float('nan')
                    
                #Count how many nonzero elements are in ch4
                nn = np.count_nonzero(short_z_main)
                nn1 = np.count_nonzero(~np.isnan(short_z_main))
                
                #Check if there are any ch4 measurements in the specified area, if so then enter the loop 
                if nn > 0 and nn1 > 0:

                    # Loop through the finer grid 
                    for i in range(num_rows):
                        for j in range(num_cols):
                        
                            # Mask out data that isn't close the point in each loop 
                            resv=.1
                            maskshort= (new_lat-resv < yy[i,j]) & (new_lat+resv > yy[i,j]) & (new_lon-resv < xx[i,j]) & (new_lon+resv > xx[i,j])
                            
                            new_lat_2 = new_lat[maskshort]
                            new_lon_2 = new_lon[maskshort]
                            short_z_main_2 = short_z_main[maskshort]
                            sat_len_2 = len(short_z_main_2)
                            if albedo is True:
                                new_SWIR_2 = new_SWIR[maskshort]
                                new_NIR_2 = new_NIR[maskshort]

                            # Create two zero variables to store the average and counter for each point 
                            average = 0 
                            short_counter = 0

                            # Loop through the coarse data 
                            for k in range(sat_len_2):
                    
                                # If the latitude and longitude is within range, the value is not nan, the qa valye and/or the albedo is within range then add to the average and increase counter 
                                #if (yy[i,j] > new_lat_2[k]-(lat_half_res * 2)) & (yy[i,j] < new_lat_2[k]+(lat_half_res * 2)) & (xx[i,j] > new_lon_2[k]-(lon_half_res * 2)) & (xx[i,j] < new_lon_2[k]+(lon_half_res * 2)):
                                if (yy[i,j] > new_lat_2[k]-(lat_half_res * 2) - 0.014) & (yy[i,j] < new_lat_2[k]+(lat_half_res * 2) + 0.007) & (xx[i,j] > new_lon_2[k]-(lon_half_res * 2) - 0.02) & (xx[i,j] < new_lon_2[k]+(lon_half_res * 2) +0.008):   
                                    if albedo is True:
                                        if not math.isnan(short_z_main_2[k]) and new_SWIR_2[k] > 0.05 and new_NIR_2[k] > 0.05:
                                            average += short_z_main_2[k]
                                            short_counter += 1
                                        else:
                                            average += 0
                                    
                                    else:
                                        if not math.isnan(short_z_main_2[k]):
                                            average += short_z_main_2[k]
                                            short_counter += 1
                                        else:
                                            average += 0
                            # If the average isn't zero at that point then divide by the counter and add that value to the finer grid and increase the counter for that point by 1 
                            if average != 0:
                                average_z[i,j] += average / short_counter  
                                counter[i,j] += 1

                # close the file 
                fh.close()


    # terate through the final average array 
    for i in range(num_rows):
        for j in range(num_cols):

            #If average ch4 is still 0, nothing has been added to that cell and so set it to NaN for the plot
            if average_z[i,j] == 0:
                average_z[i,j] = float('nan')
                
                #Otherwise divide the sum by the counter to get the average for each cell
            else:
                average_z[i,j] = (average_z[i,j] / counter[i,j])

    data_file = open(output_destination + output_file + '.txt', 'w')
    for l in range(num_rows):
        for m in range(num_cols):
            data_file.write(str(average_z[l,m]))
            data_file.write('\n')
    data_file.close()

else:
    average_z = np.zeros((num_rows, num_cols))
    lines = []
    read_file = open(output_destination + output_file + '.txt')
    for line in read_file:
        lines.append(line)
    counter = 0 
    for i in range(num_rows):
        for j in range(num_cols):
            average_z[i,j] = float(lines[counter])
            counter += 1

    read_file.close()
    
# If plot_on is True then map the data
if plot_on is True:
    
    # Extent = the area we want to look at 
    extent = [xmin, xmax, ymin, ymax]
    lon_dif = (xmax - xmin) / 4
    lat_dif = (ymax - ymin) / 4
    # Set the projection type and the extent 
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    # Add coastlines and title 
    ax.coastlines()

    # Plot data as well as a colour bar and its label 
    im = plt.pcolormesh(xx, yy, np.squeeze(average_z), cmap = 'RdYlBu_r', vmin= v_min, vmax = v_max)
    cb = plt.colorbar(im, pad = 0.1, orientation = 'horizontal', shrink = 0.6)

    if z == 'CH4':
        units = 'ppb'
    if z == 'NO2':
        units = 'μmol m-2'
    if z == 'SO2':
        units = 'mol m-2'
    if z == 'CO':
        units = 'mol m-2'
        
    cb.set_label(units)
    
    # Plot the gridlines 
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels = True, color = 'gray', linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlocator = mticker.FixedLocator([xmin, xmin + lon_dif, xmin + (lon_dif * 2), xmin + (lon_dif * 3), xmin + (lon_dif * 4), xmin + (lon_dif * 5)])
    gl.ylocator = mticker.FixedLocator([ymax, ymax - lat_dif, ymax - (lat_dif * 2), ymax - (lat_dif * 3), ymax - (lat_dif * 4), ymax - (lat_dif * 5)])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Save the figure
    plt.title(title)
    plt.savefig(output_destination + output_file + '.png')
    #plt.show()
