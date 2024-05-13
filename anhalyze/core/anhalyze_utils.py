#!/usr/bin/env python3
# coding: utf-8

# Data-related libraries
import numpy as np
import netCDF4 as nc

# Plotting-related libraries

# OS-specific libraries
import os
import socket
import deprecation

# Project custom made libraries
import anhalyze.core.anhalyze_plot_utils as apu
from anhalyze.core.anhalyze_geo import getIndex_sec


@deprecation.deprecated()
def get_paths(run_name=None, environ_paths=False):
    """ Get paths to data and mask standard locations.

    Parameters
    ----------
    run_name :  str, optional
        Simulation name (e.g. 'ANHA4-WJM004')
    environ_paths :  bool, optional
        True if using environment variables
    Returns
    -------
    data_path, mask_path: str
        directory paths for data and mask
    """

    if 'portal' not in socket.gethostname() or not run_name or environ_paths:
        # Try getting paths from environment variables
        try:
            mask_path = os.environ['MASK_PATH']
            data_path = os.environ['DATA_PATH']

        except KeyError:
            message = "Relevant paths not defined, "
            message += "please add in your .bash_profile (or .bashrc, etc..) something like this: \n\n"
            message += "#-------------------------------------------------------------\n"
            message += "# ANHALIZE setup\n"
            message += "#-------------------------------------------------------------\n"
            message += "export MASK_PATH='/root_path/user/ANALYSES/MASKS/'\n"
            message += "export DATA_PATH='/root_path/user/NEMO/ANHA4/ANHA4-EPM111-S/'\n"
            message += "#-------------------------------------------------------------\n"

            print(message)

    else:

        # Asserting input format
        assert '-' in run_name
        assert 'ANHA' in run_name
        assert len(run_name) == len('ANHA4-WJM000')

        # Get simulation info
        model_path = run_name.split('-')[0]
        user_initials = run_name.split('-')[1][1:3].upper()
        if any([x in user_initials for x in ["JM", "PM"]]):
            user_path = 'jmarson'
        elif user_initials == 'MC':
            user_path = 'madhurima'
        elif user_initials == 'EE':
            user_path = 'emilio'
        else:
            raise ValueError('Incorrect run_name.')

        data_path = f'/mnt/storage0/{user_path}/NEMO/{model_path}/{run_name}-S/'
        mask_path = f'/mnt/storage0/{user_path}/ANALYSES/MASKS/'

    return data_path, mask_path


# TODO need to move bunch of functions out of here
#   also need to separate io, timeseries, and whatever else needed


@deprecation.deprecated()
def get_lat_lon(data, lat_range, lon_range, cartesian=True):
    """  Getting Latitude and Longitude """

    # Given data selection range in lat-lon or row-col
    if cartesian:
        row_range, col_range = get_row_col_range(data, lat_range, lon_range)
    else:
        row_range, col_range = lat_range, lon_range

    lat = data['nav_lat_grid_T'][row_range[0]:row_range[1], col_range[0]:col_range[1]]
    lon = data['nav_lon_grid_T'][row_range[0]:row_range[1], col_range[0]:col_range[1]]

    return lat, lon


@deprecation.deprecated()
def get_mask(data, lat_range, lon_range, depth=0, cartesian=True):
    """  Getting Mask given Latitude and Longitude """

    # Get paths
    data_path, mask_path = get_paths()

    # Reading ANHA4 mask data
    mask = nc.Dataset(mask_path + '/'+ "ANHA4_mask.nc")

    # Extracting mask data for given depth
    tmask = mask['tmask'][0][depth]

    # Given data selection range in lat-lon or row-col
    if cartesian:
        row_range, col_range = get_row_col_range(data, lat_range, lon_range)
    else:
        row_range, col_range = lat_range, lon_range

    # Extracting mask data for given range
    surf_mask = tmask[row_range[0]:row_range[1], col_range[0]:col_range[1]]

    return surf_mask


@deprecation.deprecated()
def get_var_data(data, lat_range, lon_range, depth=0, var='votemper', masked=True, cartesian=True):
    """  Getting Data Latitude and Longitude

    depth : z axis location from 50 unit "depth"

    """

    # Get var data
    var_data = data[var][:]

    # Given data selection range in lat-lon or row-col
    if cartesian:
        row_range, col_range = get_row_col_range(data, lat_range, lon_range)
    else:
        row_range, col_range = lat_range, lon_range

    # Extracting data given lat-lon selection and depth
    var_data = var_data[0, depth, row_range[0]:row_range[1], col_range[0]:col_range[1]]

    if masked:
        # Get mask
        surf_mask = get_mask(data, lat_range, lon_range, depth=depth)

        # Mask data
        var_data.data[~np.ma.filled((1 == surf_mask.data))] = np.nan

    return var_data


@deprecation.deprecated()
def get_row_col_range(data, lat_range, lon_range, grid='gridT'):
    """ Get the row and col range given lat and lon range.  """

    # Get all lat-lon data
    if grid == 'gridT':
        lat = data['nav_lat_grid_T'][:]
        lon = data['nav_lon_grid_T'][:]
    else:
        # Get all lat-lon data
        lat = data['nav_lat'][:]
        lon = data['nav_lon'][:]

    # Create mask given lat lon values.
    lat_mask = np.ma.filled((lat.data > lat_range[0]) & (lat.data < lat_range[1]))
    lon_mask = np.ma.filled((lon.data > lon_range[0]) & (lon.data < lon_range[1]))

    # Apply masks to data
    mask = lat
    mask[~(lat_mask & lon_mask)] = 0

    # Find the row,col range by collapsing each axis.
    row_ranges = np.where(mask.data.sum(axis=1) > 0)[0]
    col_ranges = np.where(mask.data.sum(axis=0) > 0)[0]

    # Select range
    row_range = (row_ranges[0], row_ranges[-1])
    col_range = (col_ranges[0], col_ranges[-1])

    return row_range, col_range


@deprecation.deprecated()
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


@deprecation.deprecated()
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
    data_path, mask_path = get_paths(run_name=run)
    
    # Get file names list
    file_list = get_file_list(run=run, grid=grid, years_list=years_list, month_list=month_list, one_per_month=one_per_month, monthly_mean=monthly_mean)
    
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


@deprecation.deprecated()
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


@deprecation.deprecated()
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


@deprecation.deprecated()
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


@deprecation.deprecated()
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


@deprecation.deprecated()
def show_var_data_map(data, lat_range, lon_range, depth=0, var='votemper'):
    """ Displays map of given var in lat-lon range and depth.
        Note: depth has not been tested.
    """

    apu.show_var_data_map(data, lat_range, lon_range, depth=depth, var=var)


@deprecation.deprecated()
def calc_stats_var_data(data, lat_range, lon_range, depth=0, no_min_max=True, var='votemper'):
    """  """

    # Get var data
    var_data = get_var_data(data, lat_range, lon_range, depth=depth, var=var)

    # Calculating stats
    var_mean = np.nanmean(var_data)
    var_std = np.nanstd(var_data)

    if no_min_max:
        return var_mean, var_std

    else:
        var_min = np.nanmin(var_data)
        var_max = np.nanmax(var_data)

        return var_mean, var_std, var_min, var_max


