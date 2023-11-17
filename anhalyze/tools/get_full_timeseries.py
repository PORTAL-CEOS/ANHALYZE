#!/usr/bin/env python
# coding: utf-8

# Data-related libraries
import numpy as np

# Project custom-made libraries
import anhalyze_utils as au

# global variable selection.

Hudson_bay = False           # Boolean if using Hudson Bay vs James Bay locations.
depth = 0                    # z axis location from 50 unit "depth" 

if Hudson_bay: 
    hudson_east = -75
    hudson_west = -95
    hudson_north = 65
    hudson_south = 50
    loc = 'hudson_bay'

else:
    hudson_east = -78.5
    hudson_west = -82.5
    hudson_north = 54.7
    hudson_south = 51
    loc = 'james_bay'

lat_range = (hudson_south, hudson_north)
lon_range = (hudson_west, hudson_east)
    

# Making year list
year_list = [str(y) for y in np.arange(1958, 2010, 1)]
    
# Make file list from year list    
file_list = au.get_file_list(years=year_list)

# Get timeseries
timeseries_var = au.get_timeseries(file_list,  lat_range, lon_range, no_min_max=False)

# Save data
timeseries_var.to_csv(loc+'_timeseries_data.csv', index=False)
