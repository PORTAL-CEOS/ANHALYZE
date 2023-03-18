#!/usr/bin/env python
# coding: utf-8

# Data-related libraries
import numpy as np
import netCDF4 as nc
import datetime
import pandas as pd

# Plotting-related libraries
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cartopy import crs as ccrs, feature as cfeature

# OS-specific libraries
from sys import platform
import os


def get_paths():
    """ Get paths to data and mask standard locations."""

    if platform == "linux" or platform == "linux2":
        # linux
        # setup paths
        mask_path = '/mnt/storage0/jmarson/ANALYSES/MASKS/'
        data_path = '/mnt/storage0/jmarson/NEMO/ANHA4/ANHA4-EPM111-S/'
    elif platform == "darwin":
        # OS X
        # setup paths
        data_path = '/Users/jeenriquez/Documents/CEOS/ANHA4/Test_Data/'
        mask_path = '/Users/jeenriquez/Documents/CEOS/ANHA4/Test_Data/'
    else:
        data_path = ''
        mask_path = ''
        # raise ValueError("Platform not recognized.")

    return data_path, mask_path


def get_file_list(years=['1998'], grid='T', one_per_month=False, month_list=[]):
    """  Returns file list given a list of years, a grid type,
         and ether all the days in a month, or the first one.
    """

    # Get paths
    data_path, mask_path = get_paths()

    # Get complete file list from path
    file_list = os.listdir(data_path)

    selected_file_list = []

    # Selecting list of files given params
    for year in years:
        selected_file_list += (sorted([f for f in file_list if 'y'+year in f and 'grid'+grid in f]))

    # Selecting first day on given month.
    if one_per_month:
        if not month_list:
            month_list = [get_date(filename, how='m') for filename in selected_file_list]

        monthly_file_list = []

        for month in month_list:
            file_name_stump = 'y{}m{}'.format(years[0], month)
            file_month_name = [f for f in selected_file_list if file_name_stump in f][0]
            monthly_file_list.append(file_month_name)

        selected_file_list = monthly_file_list

    # Adding full path to filenames
    selected_file_list = [data_path+filename for filename in selected_file_list]

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


def get_row_col_range(data, lat_range, lon_range):
    """ Get the row col range given lat lon range.  """

    # Get all lat-lon data
    lat = data['nav_lat_grid_T'][:]
    lon = data['nav_lon_grid_T'][:]

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

    # Get var data
    var_data = get_var_data(data, lat_range, lon_range, depth=depth, var=var)

    # getting lat and lon
    lat, lon = get_lat_lon(data, lat_range, lon_range)

    # EE: should move this elsewhere
    # Setting plotting vars.
    levels = 42
    vmax = 20.
    vmin = -20.
    standard_parallels = (55, 60)
    central_longitude = -80
    cmap = 'coolwarm'

    # Set up plot
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1,
                         projection=ccrs.LambertConformal(central_longitude=central_longitude,
                                                          standard_parallels=standard_parallels))
    ax.set_extent([lon_range[0], lon_range[1], lat_range[0], lat_range[1]])

    # Adding ocean and land features
    ax.add_feature(get_feature_mask())
    ax.add_feature(get_feature_mask(feature='ocean'))

    # Plotting var data as filled contour regions
    im = ax.contourf(lon, lat, var_data, levels=levels, cmap=cmap,
                     vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), zorder=2)

    # Plotting data contour lines
    ax.contour(lon, lat, var_data, levels=levels, cmap='Greys', linewidths=.2, transform=ccrs.PlateCarree())

    # Create grid-line labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False,
                      y_inline=False, color='k', alpha=.3, linewidths=.01)
    gl.right_labels = gl.top_labels = False

    # Set Color-bar
    axins = inset_axes(ax, width="5%", height="100%", loc='right', borderpad=-1)
    label = '%s [%s]' % (data.variables[var].long_name.title(), data.variables[var].units.title())
    fig.colorbar(im, cax=axins, orientation="vertical", label=label, format='%.1f')


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


def get_date(filename, how=''):
    """  Get date from filename. Multiple formats possible.  """

    date = filename.split('_')[-2]

    if how == 'ymd':
        return int(date[1:5]), int(date[6:8]), int(date[9:11])
    elif how == 'y':
        return int(date[1:5])
    elif how == 'm':
        return int(date[6:8])
    else:
        return date


def get_timeseries(file_list, lat_range, lon_range, depth=0, var='votemper'):
    """    """

    var_means = []
    var_stds = []
    dates = []

    # Get datapoints by looping over files
    for filename in file_list:
        # Get Anha4 data
        data = nc.Dataset(filename)

        # Calculate var_data mean
        var_mean, var_std = calc_stats_var_data(data, lat_range, lon_range, depth=depth, var=var)

        # Save date from filename/time step
        date = datetime.date(get_date(filename, how='ymd'))

        data = nc.Dataset(filename)  # Append data to timeseries
        var_means.append(var_mean)
        var_stds.append(var_std)
        dates.append(date)

    # Create timeseries df
    timeseries_var = pd.DataFrame({'date': dates, 'var_mean': var_means, 'var_std': var_stds})

    return timeseries_var
