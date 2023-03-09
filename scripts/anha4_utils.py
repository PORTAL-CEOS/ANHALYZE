#!/usr/bin/env python
# coding: utf-8

# Data-related libraries
import numpy as np
import netCDF4 as nc

# Plotting-related libraries
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy
import cartopy.crs as ccrs

# OS-specific libraries
from sys import platform


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


def get_var_data(data, lat_range, lon_range, depth=0, var='votemper', cardinal=True):
    """  Getting Data Latitude and Longitude """

    # Get var data
    var_data = data[var][:]

    # Given data selection range in lat-lon or row-col
    if cardinal:
        row_range, col_range = get_row_col_range(data, lat_range, lon_range)
    else:
        row_range, col_range = lat_range, lon_range

    # Extracting data given lat-lon selection and depth
    var_data = var_data[0, depth, row_range[0]:row_range[1], col_range[0]:col_range[1]]

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


def plot_var_data(data, lat_range, lon_range, depth=0, var='votemper'):
    """  """

    # Get var data
    var_data = get_var_data(data, lat_range, lon_range, depth=depth)

    # Get mask
    surf_mask = get_mask(data, lat_range, lon_range, depth=depth)

    # Mask data
    var_data.data[~np.ma.filled((1 == surf_mask.data))] = np.nan
    surf_mask.data[np.ma.filled((1 == surf_mask.data))] = np.nan

    # getting lat and lon
    lat, lon = get_lat_lon(data, lat_range, lon_range)

    # Setting plotting vars.
    levels = 42
    vmax = 20.
    vmin = -20.
    standard_parallels = (55, 60)
    central_longitude = -80
    cmap = 'coolwarm'
    ocean_color = '#000066'
    land_color = matplotlib.colors.to_hex('wheat')
    
    # Set up plot
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1,
                         projection=ccrs.LambertConformal(central_longitude=central_longitude,
                                                          standard_parallels=standard_parallels))
    ax.set_extent([lon_range[0], lon_range[1], lat_range[0], lat_range[1]])

    # Adding ocean and land features
    ax.add_feature(cartopy.feature.OCEAN, facecolor=ocean_color)
    ax.add_feature(cartopy.feature.LAND, facecolor=land_color)

    # Plotting var data as filled countour regions
    im = ax.contourf(lon, lat, var_data, levels=levels, cmap=cmap,
                     vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

    # Create grid-line labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False,
                      y_inline=False, color='k', alpha=.3, linewidths=.01)
    gl.right_labels = gl.top_labels = False

    # Set Color-bar
    axins = inset_axes(ax, width="5%", height="100%", loc='right', borderpad=-1)
    label = data.variables[var].long_name + ' [' + data.variables[var].units + ']'
    fig.colorbar(im, cax=axins, orientation="vertical", label=label, format='%.1f')

