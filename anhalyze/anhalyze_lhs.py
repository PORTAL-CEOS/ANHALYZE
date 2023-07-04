#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 2023

@author: Luiz Henrique da Silva
based on Dr. Emilio Enriquez codes

"""

# Data-related libraries
import numpy as np
import netCDF4 as nc
import datetime
import pandas as pd

# Plotting-related libraries
import matplotlib
import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# OS-specific libraries
import sys
import os, glob




def get_date(filename, how=''):
    """  Get date from filename.
         Assuming filename format: */*/ANHA4-EPM111_y1998m04d05_gridB.nc
         Multiple output formats possible.
    """

    # Get full date from filename
    date = filename.split('_')[-2]

    # Return format specific date info
    if how == 'ymd':
        return int(date[1:5]), int(date[6:8]), int(date[9:11])
    elif how == 'y':
        return int(date[1:5])
    elif how == 'm':
        return int(date[6:8])
    else:
        return date



def get_row_col_range(data, grid, lat_range, lon_range):
    """ Get the row and col range given lat and lon range.  """
    
    if grid == 'gridT':
        # Get all lat-lon data
        lat = data['nav_lat_grid_T'][:]
        lon = data['nav_lon_grid_T'][:]

    else:
        # Get all lat-lon data
        lat = data['nav_lat'][:]
        lon = data['nav_lon'][:]
        
    # Create mask given lat lon values.
    lat_mask = np.ma.filled((lat.data >= lat_range[0]) & (lat.data <= lat_range[1]))
    lon_mask = np.ma.filled((lon.data >= lon_range[0]) & (lon.data <= lon_range[1]))

    # Apply masks to data
    mask = lat
    mask[~(lat_mask & lon_mask)] = 0

    # Find the row,col range by collapsing each axis.
    row_ranges = np.where(mask.data.sum(axis=1) > 0)[0]
    col_ranges = np.where(mask.data.sum(axis=0) > 0)[0]

    # Select range
    row_range = (row_ranges[0], row_ranges[-1]+1)
    col_range = (col_ranges[0], col_ranges[-1]+1)

    return row_range, col_range



def get_paths(run):
    """ Get paths (in portal) to data and mask from a given simulation.
    
    Input
    run = Simulation name (e.g. 'ANHA4-WJM004'); type: string
    
    """

    run = run.upper()
    
    if sys.platform == "linux" or sys.platform == "linux2":
        # linux
        # setup paths
        # Set directory and sections mask file
        
        if any([x in run.split('-')[1][1:3].upper() for x in ["JM","PM"]]):
    
            mask_path = "/mnt/storage0/jmarson/ANALYSES/MASKS/"
            data_path = '/mnt/storage0/jmarson/NEMO/ANHA4/{}-S/'.format(run)
        
        
        elif run.split('-')[1][1:3].upper() == 'MC':
            user = os.environ['USER']
            
            mask_path = "/mnt/storage0/{}/ANALYSES/MASKS/".format(user)
            data_path = '/mnt/storage0/madhurima/NEMO/ANHA4/{}-S/'.format(run)
        
        else:
            data_path = ''
            mask_path = ''
        
    else:
        data_path = ''
        mask_path = ''
        # raise ValueError("Platform not recognized.")

    return data_path, mask_path



def get_file_list(run,grid, years_list, month_list=None, one_per_month=False, monthly_mean=False):
        """  Returns file list given a list of years, a grid type,
             and either all the days in a month, or the first one.
             
             Inputs:
             run = Simulation name (e.g. 'ANHA4-WJM004'); type: string
             grid = Simulation output grid. (e.g. gridT, gridU); type: string
             years_list = List of years to be analysed; type: string
             month_list = List of months. If empty, return all 12 months; type: string
             one_per_month = True or False; If True return only the first file of each month
             monthly_mean = True or False; If not using the default 5 days averaged outputs, 
                             the function ask the user to insert the monthly mean path.
             
        """

        # Setup
        selected_file_list = []
        monthly_file_list = []
        if month_list is None:
            month_list = []

            
        # Get paths
        data_path, mask_path = get_paths(run=run)
        
        if monthly_mean:
            print('---------------------------')
            data_path = str(input('Insert monthly averaged data path: '))
        
        # Get complete file list from path
        file_list = os.listdir(data_path)

        # Selecting list of files given params
        for year in years_list:
            selected_file_list += (sorted([f for f in file_list if '_y'+year in f and '_'+grid in f]))

        # Selecting first day on given month
        if one_per_month:
            if not month_list:
                month_list = [get_date(filename, how='m') for filename in selected_file_list]

            for year in years_list:
                for month in month_list:
                    file_name_stump = 'y{}m{}'.format(year, month)
                    file_month_name = [f for f in selected_file_list if file_name_stump in f][0]
                    monthly_file_list.append(file_month_name)
        else:
            # Make month selection
            if month_list:
                for year in years_list:
                    for month in month_list:
                        file_name_stump = 'y{}m{}'.format(year, month)
                        file_month_names = [f for f in selected_file_list if file_name_stump in f]
                        monthly_file_list += file_month_names

        if monthly_file_list:
            selected_file_list = monthly_file_list

        # Adding full path to filenames
        selected_file_list = [data_path+filename for filename in selected_file_list]

        return selected_file_list
    
    
    
def getIndex_sec(sectName):
    """Return the number associated to each section.
    One can check the respective numbers and sections on the Lab Guide.
   
    Inputs
    sectName = Section name; type: String
   
    Options:
    Bering Strait
    Lancaster Sound
    Jones Sound
    Nares Strait
    Davis Strait
    Fram Strait

"""

# Todo, add the tilted sections. So far, the section functions only work with
# zonal or meridional oriented sections.

    if sectName == 'Bering Strait':
        sect = 3180
    elif sectName=='Lancaster Sound':
        sect = 4180
    elif sectName=='Jones Sound':
        sect = 4180
    elif sectName=='Nares Strait':
        sect = 4270
    elif sectName=='Davis Strait':
        sect = 9360
    elif sectName=='Fram Strait':
        sect = 5360
    else:
        print("Section name not found.\nCall help/doc to check the sections available in this function.")
        
    return sect



def getMask_region(run,depth,row_range,col_range):
    """ Return the land/ocean ANHA4 mask given lat and lon sliced area
    
    Inputs
    run = Simulation name (e.g.. 'ANHA4-WJM004'); type: string
    depth = Number of depth layers including surface (e.g. Whole water column = 50)
    row_range = Row max and min indices
    col_ranfe = Col max and min indices
    
    """

        # Get paths
    data_path, mask_path = get_paths(run=run)
    mask_file = mask_path + '/ANHA4_mask.nc'
                
    # Using the rows and cols indices to extract the mask
    # over the region of interest.
    mask = nc.Dataset(mask_file)
    mask_region = mask['tmask'][:,:depth,row_range[0]:row_range[1],col_range[0]:col_range[1]]
        
    return mask_region[:].data



def getVar_region(run, grid, depth, lon_range,lat_range,var,years_list, month_list=None, one_per_month=False,monthly_mean=False,masked=True,cardinal=True):
    
    """
    
    Extract a variable in an pre-defined area from ANHA4 domain.
    
    Inputs:
    run = Simulation (i.e. ANHA4-WJM004); type: string
    grid = Variable model grid (i.e. griT, gridB, icemod); type: string
    depth = Number of depth levels, including surface; type: integer
    lon_range = A list with west and east region limits; type: integer/float
    lat_range = A list with south and north region limits; type: integer/float
    var = Variable name (i.e. votemper, vosaline); type: string
    years_list = List made of all years of interest; type: string
    month_list = List made of all years of interest; type: string
    one_per_month = If True, extract just the first day of each month
    monthly_mean = If True, ask the user to input the monthly mean directory, since this is not
                   a standard PORTAL output format.
    
    """
    
    # Get mask path
    data_path, mask_path = get_paths(run=run)
    
    # Get file names list
    file_list = get_file_list(run=run,grid=grid, years_list=years_list, month_list=month_list, one_per_month=one_per_month, monthly_mean=monthly_mean)
    
    # Get depth levels
    depth_levels = nc.Dataset(file_list[0])['deptht'][:depth]
    
    # Get latitude and longitude limits and indices
    if cardinal:
        row_range, col_range = get_row_col_range(nc.Dataset(file_list[0]),grid, lat_range, lon_range)
    else:
        row_range, col_range = lat_range, lon_range

    # Get the region land/ocean mask
    mask_region = getMask_region(run,row_range,col_range)

    # For tracers
    if grid == 'gridT':
        
        # Get latitude and longitude
        lon = nc.Dataset(file_list[0])['nav_lon_grid_T'][row_range[0]:row_range[1],col_range[0]:col_range[1]]
        lon[lon == 0] = np.nan
        lat = nc.Dataset(file_list[0])['nav_lat_grid_T'][row_range[0]:row_range[1],col_range[0]:col_range[1]]
        lat[lat == 0] = np.nan
        
        # Prealocating variable
        var_region = np.empty([1, depth, np.size(lat,0),np.size(lon,1)])
        
        # Extracting data
        for filename in file_list:
            ds = nc.Dataset(filename)
            ds_var = ds[var][:,:depth,row_range[0]:row_range[1],col_range[0]:col_range[1]]
            ds_var[mask_region == 0] = np.nan
            var_region = np.append(var_region,ds_var.data,axis=0)
                                              
        
        # Remove empty first time dimension
        var_region = np.squeeze(var_region[1:,:,:,:])
        
    # For velocities
    elif any([g in grid for g in ['gridU','gridV','gridW']]):
    
        # Pre alocating variable
        var_region = np.empty([1, depth, np.size(lat,0),np.size(lon,1)])
        
        # Get latitude and longitude
        lon = nc.Dataset(file_list[0])['nav_lon'][row_range[0]:row_range[1],col_range[0]:col_range[1]]
        lon[lon == 0] = np.nan
        lat = nc.Dataset(file_list[0])['nav_lat'][row_range[0]:row_range[1],col_range[0]:col_range[1]]
        lat[lat == 0] = np.nan
        # Extracting data
        for filename in file_list:
            ds = nc.Dataset(filename)
            ds_var = ds[var][:,:depth,row_range[0]:row_range[1],col_range[0]:col_range[1]]
            ds_var[mask_region == 0] = np.nan
            var_region = np.append(var_region,ds_var.data,axis=0)
        
        # Remove empty first time dimension
        var_region = var_region[1:,:,:,:]
        
    else:
        # Get latitude and longitude
        lon = nc.Dataset(file_list[0])['nav_lon'][row_range[0]:row_range[1],col_range[0]:col_range[1]]
        lon[lon == 0] = np.nan
        lat = nc.Dataset(file_list[0])['nav_lat'][row_range[0]:row_range[1],col_range[0]:col_range[1]]
        lat[lat == 0] = np.nan
        
        # Prealocating variable
        var_region = np.empty([1, np.size(lat,0),np.size(lon,1)])

        # Extracting data
        for filename in file_list:
            ds = nc.Dataset(filename)
            ds_var = ds[var][:,row_range[0]:row_range[1],col_range[0]:col_range[1]]
            ds_var[mask_region[0,:,:] == 0] = np.nan
            var_region = np.append(var_region,ds_var.data,axis=0)
                                              
        
        # Remove empty first time dimension
        var_region = var_region[1:,:,:]
        
    return var_region, lon, lat, depth_levels



def getClim_region(run, grid, depth, lon_range, lat_range, var, years_list,monthly_mean=False):
    
    """ Get the climatologic year of a choosen variable from an especific region """

    # Month of the year list
    month_list = ['{:02d}'.format(mm) for mm in np.arange(1,13,1)]
    
    
    # Climatologic averaged
    
    for mm in month_list:
        
        if mm == '01':
            
            var_region, lon, lat, depth_levels = getVar_region(run,grid,depth, lon_range, lat_range, var, years_list=years_list,month_list=[mm],monthly_mean=monthly_mean)
            # Average the data
            var_region_avg = np.nanmean(var_region,axis=0,keepdims=True)
            # Create climatologic np array
            var_clim = var_region_avg
        else:
            
            # Get var data
            var_region, lon, lat, depth_levels = getVar_region(run,grid,depth, lon_range, lat_range, var, years_list=years_list,month_list=[mm],monthly_mean=monthly_mean)
        
            # Average the data
            var_region_avg = np.nanmean(var_region,axis=0,keepdims=True)
        
            # Concatenate to the climatologic array
            var_clim = np.append(var_clim,var_region_avg,axis=0)
    
    return var_clim, lon, lat, depth_levels



def getMask_sec(run,depth,sectName,i,j):

    # Get paths
    data_path, mask_path = get_paths(run=run)
    mask_file = mask_path + '/ANHA4_mask.nc'
        
    if any([n in sectName for n in ["Davis","Lancaster","Jones","Bering"]]):
    
        mask = nc.Dataset(mask_file)
        mask_sect = np.squeeze(mask['tmask'][:,:depth,i[0]:i[-1]+1,j[0]:j[-1]+1],axis=2)
                                              
    else:           
        mask = nc.Dataset(mask_file)
        mask_sect = np.squeeze(mask['tmask'][:,:depth,i[0]:i[-1]+1,j[0]:j[-1]+1],axis=3)
        
    return mask_sect[:].data
    
    
    
def getVar_sec(run, sectName, grid, depth, var,years_list, month_list=None, one_per_month=False,monthly_mean=False):
    
    """ var_sec, lon, lat, depth_levels = getVar_sec(run, sectName, grid, depth, years_list, month_list=None, var, one_per_month=False,monthly_mean=False):


    Get variable in one of the ANHA4 sections.

    run = Name of the run (e.g. ANHA4-WJM004), type: string
    var = Variable to be extract (e.g. votemper), type: string
    grid = What kind of grid your are looking for?
    (e.g. griT, gridU(V),icemod),type: string
    secName = The name of the section, type: string

    Bering Strait = 3180
    Lancaster Sound = 4180
    Jones Sound = 4180
    Nares Strait = 4270
    Davis Strait = 9360
    Fram Strait = 5360

    depth = Number of vertical levels (include surface), type: integer
    years_list: , List of years; type: str list
    month_list: List of month, if empty will get all 12 months of the year; type: str list

    Obs:
    For the outputs of gridU and gridV, you are getting 2 sections.
    Its up to you to keep that way or average to fit into
    the gridT coordinates.

    """

    # Get mask path
    data_path, mask_path = get_paths(run=run)
    
    # Get file names list
    file_list = get_file_list(run=run,grid=grid, years_list=years_list, month_list=month_list, one_per_month=one_per_month, monthly_mean=monthly_mean)
    
    # Get depth levels
    depth_levels = nc.Dataset(file_list[0])['deptht'][:depth]
    
    # Get the Index from the transect tracer mask
    fileMask = mask_path+'ANHA4_trc_sec_mask_Nov2022.nc';

    if os.path.isfile(fileMask):
        # Get tmask data
        mask = nc.Dataset(fileMask);
        tmask = mask['tmask'][:].data
    
        # Get the i's and j's indices
        sect = getIndex_sec(sectName=sectName)
    
        i = np.where(tmask == sect)[0]
        j = np.where(tmask == sect)[1]
    
    
    # Get the vertical section mask
    mask_sect = getMask_sec(run=run,sectName=sectName,i=i,j=j)

    # For tracers
    if any([g in grid for g in ['gridT','gridB']]):
        
        # Prealocating variable
        var_sec = np.empty([1, depth, len(i)])
        
        # Get latitude and longitude
        lon = np.squeeze(nc.Dataset(file_list[0])['nav_lon_grid_T'][i[0]:i[-1]+1,j[0]:j[-1]+1])
        lat = np.squeeze(nc.Dataset(file_list[0])['nav_lat_grid_T'][i[0]:i[-1]+1,j[0]:j[-1]+1])
        
        # Extracting data
        if any([n in sectName for n in ["Davis","Lancaster","Jones","Bering"]]):
    
            for filename in file_list:
                ds = nc.Dataset(filename)
                ds_var = np.squeeze(ds[var][:,:depth,i[0]:i[-1]+1,j[0]:j[-1]+1],axis=2)
                ds_var[mask_sect == 0] = np.nan
                var_sec = np.append(var_sec,ds_var.data,axis=0)
                                              
        else:           
            for filename in file_list:
                ds = nc.Dataset(filename)
                ds_var = np.squeeze(ds[var][:,:depth,i[0]:i[-1]+1,j[0]:j[-1]+1],axis=3)
                ds_var[mask_sect == 0] = np.nan
                var_sec = np.append(var_sec,ds_var,axis=0)
        
        # Remove empty first time dimension
        var_sec = var_sec[1:,:,:]
        
    # For velocities
    elif grid == 'gridU':
    
        # Pre alocating variable
        var_sec = np.empty([1, depth, len(i), 2])
        
        # Get latitude and longitude
        lon = nc.Dataset(file_list[0])['nav_lon'][i[0]:i[-1]+1,j[0]:j[-1]+1]
        lat = nc.Dataset(file_list[0])['nav_lat'][i[0]:i[-1]+1,j[0]:j[-1]+1]
        
        # Extracting data
        for filename in file_list:
            ds = nc.Dataset(filename)
            ds_var = ds[var][:,:depth,i[0]:i[-1]+1,j[0]-1:j[-1]+1]
            ds_var[mask_sect == 0] = np.nan
            var_sec = np.append(var_sec,ds_var.data,axis=0)
        
        # Remove empty first time dimension
        var_sec = var_sec[1:,:,:,:]
    
    elif grid == 'gridV':
        
        # Prealocating variable
        var_sec = np.empty([1, depth, 2, len(i)])
        
        # Get latitude and longitude
        lon = nc.Dataset(file_list[0])['nav_lon'][i[0]:i[-1]+1,j[0]:j[-1]+1]
        lat = nc.Dataset(file_list[0])['nav_lat'][i[0]:i[-1]+1,j[0]:j[-1]+1]
        
        # Extracting data
        for filename in file_list:
            ds = nc.Dataset(filename)
            ds_var = ds[var][:,:depth,i[0]-1:i[-1]+1,j[0]:j[-1]+1]
            ds_var[mask_sect == 0] = np.nan
            var_sec = np.append(var_sec,ds_var.data,axis=0)
        
        # Remove empty first time dimension   
        var_sec = var_sec[1:,:,:,:]
        
    return var_sec, lon, lat, depth_levels



def getClim_sec(run, sectName, grid, depth, var, years_list, monthly_mean=False):
    
    """ Get the climatologic year of a choosen variable from a predefined ANHA section 

    run = Name of the run (e.g. ANHA4-WJM004), type: string
    var = Variable to be extract (e.g. votemper), type: string
    grid = What kind of grid your are looking for 
    (e.g. griT, gridU(V),icemod),type: string
    secName = The name of the section, type: string

    Bering Strait
    Lancaster Sound
    Jones Sound
    Nares Strait
    Davis Strait
    Fram Strait

    depth = Number of vertical levels (include surface), type: integer
    years_list: , List of yearstype: str list
    month_list: List of month, if empty will get all 12 months of the year, type: str list

    Obs:
    For the outputs of gridU and gridV, you are getting 2 sections.
    Its up to you to keep that way or average to fit into
    the gridT coordinates.

    

    """

    # Month of the year list
    month_list = ['{:02d}'.format(mm) for mm in np.arange(1,13,1)]
    
    
    # Loop through months to create the climatologic year
    
    for mm in month_list:
        
            if mm == '01':
            
                var_sec, lon, lat, depth_levels = getVar_sec(run,sectName,grid,depth, var, years_list=years_list,month_list=[mm],monthly_mean=monthly_mean)
                # Average the data
                var_sec_avg = np.nanmean(var_sec,axis=0,keepdims=True)
                # Create climatologic np array
                var_clim_sec = var_sec_avg
            
            else:

                # Get var data
                var_sec, lon, lat, depth_levels = getVar_sec(run,sectName,grid,depth, var, years_list=years_list,month_list=[mm],monthly_mean=monthly_mean)
        
                # Average the data
                var_sec_avg = np.nanmean(var_sec,axis=0,keepdims=True)
        
                # Concatenate to the climatologic matrix
                var_clim_sec = np.append(var_clim_sec,var_sec_avg,axis=0)
    
    return var_clim_sec, lon, lat, depth_levels


