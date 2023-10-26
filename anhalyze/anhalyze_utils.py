#!/usr/bin/env python3
# coding: utf-8

# Data-related libraries
import numpy as np
import netCDF4 as nc
import datetime
import pandas as pd

# Plotting-related libraries
import matplotlib
from cartopy import feature as cfeature

# OS-specific libraries
import os
import socket

# Project custom made libraries
import anhalyze_plot_utils as apu


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


class Anhalyze:
    """ This class will do analysis of ANHA4 data, for now it initializes the location of files.
   ...

    Attributes
    ----------
        run_name : str
            Simulation name (e.g. 'ANHA4-WJM004')  (default is ??)
        grid : str
            Simulation output grid. (e.g. gridT, gridU)   (default is T)
        years : str
            List of years to be analysed  (default is [1998])
        month_list : str
            List of months to be analysed. If None(default), return all 12 months
        one_per_month : bool
            If True return only the first file of each month (default is False)
        verbose : bool
            Print outputs (default is True)

    Methods
    -------
    get_file_list(month_list=None, one_per_month=False, monthly_mean=False)
        Returns file list given a list of years, a grid type,
        and either all the days in a month, or the first one.

    """

    def __init__(self, run_name=None, grid=None, years=None, month_list=None, one_per_month=False, verbose=True):
        """ Initializing class.

        Parameters
        ----------
        run_name : str, optional
            Simulation name (e.g. 'ANHA4-WJM004')  (default is ??)
        grid : str, optional
            Simulation output grid. (e.g. gridT, gridU)   (default is T)
        years : str, optional
            List of years to be analysed  (default is [1998])
        month_list : str, optional
            List of months to be analysed. If None(default), return all 12 months
        one_per_month : bool, optional
            If True return only the first file of each month (default is False)
        verbose : bool, optional
            Print outputs (default is True)

         TODO: where is the nc.Dataset to get the data?
        """

        # -------
        # The global wild west
        self.run_name = run_name
        if grid is None:
            self.grid = 'T'
        else:
            self.grid = grid
        # TODO: add assert to grid values
        if years is None:
            self.years = ['1998']
        else:
            self.years = years
        # TODO: add assert to years
        self.month_list = month_list
        self.one_per_month = one_per_month
        self.verbose = verbose

        # Init file list given conditions
        self.file_list = self.get_file_list(one_per_month=one_per_month, month_list=month_list)

        if verbose:
            print(self.file_list)

    def get_file_list(self, month_list=None, one_per_month=False, monthly_mean=False):
        """  Returns file list given a list of years, a grid type,
             and either all the days in a month, or the first one.

        Parameters
        ----------
        month_list :  str, optional
            List of months to be analysed. If None(default), return all 12 months
        one_per_month :  bool, optional
            If True return only the first file of each month (default is False)
        monthly_mean :  bool, optional
            If not using the default 5 days averaged outputs,
            The function ask the user to insert the monthly mean path.

        Returns
        -------
        selected_file_list: list
            List of files for analysis

        """

        # Setup
        selected_file_list = []
        monthly_file_list = []
        if month_list is None:
            month_list = []

        # Get paths
        data_path, mask_path = get_paths(self.run_name)

        # TODO: update code below (from luiz)
        if monthly_mean:
            print('---------------------------')
            data_path = str(input('Insert monthly averaged data path: '))

        # Get complete file list from path
        file_list = os.listdir(data_path)

        # Selecting list of files given params
        for year in self.years:
            selected_file_list += (sorted([f for f in file_list if '_y' + year in f and '_grid' + self.grid in f]))

        # Selecting first day on given month
        if one_per_month:
            if not month_list:
                month_list = [get_date(filename, how='m') for filename in selected_file_list]

            for year in self.years:
                for month in month_list:
                    file_name_stump = 'y{}m{}'.format(year, month)
                    file_month_name = [f for f in selected_file_list if file_name_stump in f][0]
                    monthly_file_list.append(file_month_name)
        else:
            # Make month selection
            if month_list:
                for year in self.years:
                    for month in month_list:
                        file_name_stump = 'y{}m{}'.format(year, month)
                        file_month_names = [f for f in selected_file_list if file_name_stump in f]
                        monthly_file_list += file_month_names

        if monthly_file_list:
            selected_file_list = monthly_file_list

        # Adding full path to filenames
        selected_file_list = [data_path + filename for filename in selected_file_list]

        return selected_file_list


def get_lat_lon(data, lat_range, lon_range, cardinal=True):
    """  Getting Latitude and Longitude """

    # Given data selection range in lat-lon or row-col
    if cardinal:
        row_range, col_range = get_row_col_range(data, lat_range, lon_range)
    else:
        row_range, col_range = lat_range, lon_range

    lat = data['nav_lat_grid_T'][row_range[0]:row_range[1], col_range[0]:col_range[1]]
    lon = data['nav_lon_grid_T'][row_range[0]:row_range[1], col_range[0]:col_range[1]]

    return lat, lon


def get_mask(data, lat_range, lon_range, depth=0, cardinal=True):
    """  Getting Mask given Latitude and Longitude """

    # Get paths
    data_path, mask_path = get_paths()

    # Reading ANHA4 mask data
    mask = nc.Dataset(mask_path + "ANHA4_mask.nc")

    # Extracting mask data for given depth
    tmask = mask['tmask'][0][depth]

    # Given data selection range in lat-lon or row-col
    if cardinal:
        row_range, col_range = get_row_col_range(data, lat_range, lon_range)
    else:
        row_range, col_range = lat_range, lon_range

    # Extracting mask data for given range
    surf_mask = tmask[row_range[0]:row_range[1], col_range[0]:col_range[1]]

    return surf_mask


def get_var_data(data, lat_range, lon_range, depth=0, var='votemper', masked=True, cardinal=True):
    """  Getting Data Latitude and Longitude

    depth : z axis location from 50 unit "depth"

    """

    # Get var data
    var_data = data[var][:]

    # Given data selection range in lat-lon or row-col
    if cardinal:
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
    data_path, mask_path = get_paths(run_name=run)
    
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




def get_feature_mask(feature='land', resolution='50m'):
    """   """

    # Select facecolor
    if 'land' in feature:
        facecolor = matplotlib.colors.to_hex('wheat')
    elif 'ocean' in feature:
        facecolor = '#000066'
    else:
        facecolor = matplotlib.colors.to_hex('gray')

        # Construct feature mask
    feature_mask = cfeature.NaturalEarthFeature('physical', feature,
                                                scale=resolution,
                                                edgecolor='face',
                                                facecolor=facecolor)

    return feature_mask


def show_var_data_map(data, lat_range, lon_range, depth=0, var='votemper'):
    """ Displays map of given var in lat-lon range and depth.
        Note: depth has not been tested.
    """

    apu.show_var_data_map(data, lat_range, lon_range, depth=depth, var=var)


def init_location(hudson_bay=True):
    """
    Setting up model parameters based on Hudson Bay and James Bay locations.

    Hudson_bay : Boolean if using Hudson Bay vs James Bay locations.
    """

    if hudson_bay:
        # proj_size = 2
        region = 'HudsonBay'

        # Setting up James Bay location
        east = -75
        west = -93
        north = 65
        south = 50
        standard_parallels = (52.5, 62.5)
        central_longitude = -80

        lat_range = (south, north)
        lon_range = (west, east)

    else:

        # proj_size = 4
        region = 'JamesBay'

        # Setting up James Bay location
        east = -78.5
        west = -82.5
        north = 54.7
        south = 51
        standard_parallels = (52, 53)
        central_longitude = -80

        lat_range = (south, north)
        lon_range = (west, east)

    location_info = {'lat_range': lat_range,
                     'lon_range': lon_range,
                     'region': region,
                     'standard_parallels': standard_parallels,
                     'central_longitude': central_longitude,
                     # 'proj_size': proj_size,
                     }

    return location_info


def show_var_data_maps(file_list, lat_range, lon_range, depth=0, var='votemper'):
    """ Displays map of given var in lat-lon range and depth.
        Note: depth has not been tested.
    """

    apu.show_var_data_maps(file_list, lat_range, lon_range, depth=depth, var=var)


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


def get_timeseries(file_list, lat_range, lon_range, depth=0, no_min_max=True, var='votemper'):
    """    """

    # Setting up
    dates = []
    var_means = var_stds = []
    var_mins = var_maxs = []
    var_min = var_max = None

    # Get datapoints by looping over files
    for filename in file_list:
        # Get ANHA4 data
        data = nc.Dataset(filename)

        # Calculate var_data mean
        if no_min_max:
            var_mean, var_std = calc_stats_var_data(data, lat_range, lon_range, depth=depth,
                                                    no_min_max=no_min_max, var=var)
        else:
            var_mean, var_std, var_min, var_max = calc_stats_var_data(data, lat_range, lon_range, depth=depth,
                                                                      no_min_max=no_min_max, var=var)

        # Save date from filename/time step
        y, m, d = get_date(filename, how='ymd')
        date = datetime.date(y, m, d)

        # Append data to timeseries
        var_means.append(var_mean)
        var_stds.append(var_std)
        dates.append(date)
        if not no_min_max:
            var_mins.append(var_min)
            var_maxs.append(var_max)

    # Create timeseries df
    if no_min_max:
        timeseries_var = pd.DataFrame({'date': dates,
                                       'var_mean': var_means,
                                       'var_std': var_stds})
    else:
        timeseries_var = pd.DataFrame({'date': dates,
                                       'var_mean': var_means,
                                       'var_std': var_stds,
                                       'var_min': var_mins,
                                       'var_max': var_maxs})

    return timeseries_var


def calc_timeseries(timeseries, action='g_mean', n_year=None):
    """  Generalized code that calculates one specific operation, depending on how it is called.
    """

    # Making timeseries copy
    timeseries = timeseries.copy()

    # Calculating full period climatology stats
    if 'g_' in action:

        # Grouping timeseries
        grouped_timeseries = timeseries.groupby('wrap_day')[['var_mean', 'var_std']]

        if 'mean' in action:
            action_timeseries = grouped_timeseries.mean().reset_index()
        elif 'quantile90' in action:
            action_timeseries = grouped_timeseries.quantile(.9).reset_index()
        elif 'quantile10' in action:
            action_timeseries = grouped_timeseries.quantile(.1).reset_index()
        elif 'median' in action:
            action_timeseries = grouped_timeseries.median().reset_index()
        elif 'max' in action:
            action_timeseries = grouped_timeseries.max().reset_index()
        else:
            action_timeseries = None

    # Adding climatology and MHW references to timeseries dataframe for easy comparison.
    elif 'add_' in action:

        # Adding climatology stats to timeseries dataframe for easy comparison.
        if 'long' in action:
            action_timeseries = np.array(timeseries['var_mean'].to_list() * n_year)

        # Calculating MHW categories to timeseries dataframe.
        else:
            delta_t = timeseries['var_mean_quantile'] - timeseries['var_mean_mean']

            if '2T' in action:
                action_timeseries = timeseries['var_mean_mean'] + 2*delta_t
            elif '3T' in action:
                action_timeseries = timeseries['var_mean_mean'] + 3*delta_t
            elif '4T' in action:
                action_timeseries = timeseries['var_mean_mean'] + 4*delta_t
            else:
                action_timeseries = None
    elif 'remove_' in action:

        # Calculating Marine Cold Spells categories to timeseries dataframe.
        delta_t = timeseries['var_mean_mean'] - timeseries['var_mean_quantile']

        if '2T' in action:
            action_timeseries = timeseries['var_mean_mean'] - 2*delta_t
        elif '3T' in action:
            action_timeseries = timeseries['var_mean_mean'] - 3*delta_t
        elif '4T' in action:
            action_timeseries = timeseries['var_mean_mean'] - 4*delta_t
        else:
            action_timeseries = None

    else:
        action_timeseries = None

    return action_timeseries


def anhalize_timeseries(raw_timeseries, mhw=True):
    """
    """

    # Set year vars
    year_standard = 2000
    year_min = 1958
    year_max = 2009
    n_year = year_max - year_min + 1

    if mhw:
        actions = ['g_quantile90', 'add_2T', 'add_3T', 'add_4T']
    else:
        actions = ['g_quantile10', 'remove_2T', 'remove_3T', 'remove_4T']

    # Make copy of raw data
    anhalyzed_timeseries = raw_timeseries.copy()

    # Change date formatting
    anhalyzed_timeseries['date'] = pd.to_datetime(anhalyzed_timeseries['date'], format='%Y-%m-%d')

    # Add Year column
    anhalyzed_timeseries['year'] = anhalyzed_timeseries.date.dt.year

    # Add Month column
    anhalyzed_timeseries['month'] = anhalyzed_timeseries.date.dt.month

    # Add day column
    anhalyzed_timeseries['day'] = anhalyzed_timeseries.date.dt.day

    # Change date formatting
    anhalyzed_timeseries['date2'] = anhalyzed_timeseries['date'].dt.date
    anhalyzed_timeseries.drop('date', axis=1, inplace=True)
    anhalyzed_timeseries.rename(columns={'date2': 'date'}, inplace=True)

    # Add wrap-day column
    anhalyzed_timeseries['wrap_day'] = anhalyzed_timeseries.apply(
        lambda row: datetime.date(year_standard, row.month, row.day), axis=1)

    # Folding data yearly to get day stats.
    timeseries_year_mean = calc_timeseries(anhalyzed_timeseries, action='g_mean')
    timeseries_year_quantile = calc_timeseries(anhalyzed_timeseries, action=actions[0])
    # timeseries_year_max = calc_timeseries(anhalyzed_timeseries, action='g_max')
    # timeseries_year_median = calc_timeseries(anhalyzed_timeseries, action='g_median')

    # Calculating timeseries mean/quantile values and adding them to timeseries
    anhalyzed_timeseries['var_mean_mean'] = calc_timeseries(timeseries_year_mean,
                                                            action='add_long',
                                                            n_year=n_year)
    anhalyzed_timeseries['var_mean_quantile'] = calc_timeseries(timeseries_year_quantile,
                                                                action='add_long',
                                                                n_year=n_year)

    # Calculating timeseries MHW categories and adding them to timeseries
    anhalyzed_timeseries['var_mean_2T'] = calc_timeseries(anhalyzed_timeseries, action=actions[1])
    anhalyzed_timeseries['var_mean_3T'] = calc_timeseries(anhalyzed_timeseries, action=actions[2])
    anhalyzed_timeseries['var_mean_4T'] = calc_timeseries(anhalyzed_timeseries, action=actions[3])

    return anhalyzed_timeseries


def plot_timeseries(timeseries_var, data_variables, lat_range, lon_range, var='votemper'):
    """  Wrapper to plot timeseries function  """

    apu.plot_timeseries(timeseries_var, data_variables, lat_range, lon_range, var=var)


def plot_mhw(anhalyzed_timeseries, year=1998, remove_mean=True, show_cat4=False, region="James Bay", mhw=True):
    """  Wrapper to plot mhw function  """

    apu.plot_mhw(anhalyzed_timeseries, year=year, remove_mean=remove_mean, show_cat4=show_cat4, region=region, mhw=mhw)


def find_mhw_info(anhalyzed_timeseries, mhw=True):
    """  placeholder, to find dates and deltaTs/categories on mhws or mcs """

#    if mhw:

    return None
