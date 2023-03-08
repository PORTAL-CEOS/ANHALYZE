#!/usr/bin/env python
# coding: utf-8

# Data related libraries
import numpy as np
import netCDF4 as nc
import glob
import os

# Plotting related libraries
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Machine specific libraries
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

    return data_path,mask_path

def get_lat_lon(data, x_range, y_range):
    """  Getting Latitude and Longitude """

    lat = data['nav_lat_grid_T'][x_range[0]:x_range[1], y_range[0]:y_range[1]]
    lon = data['nav_lon_grid_T'][x_range[0]:x_range[1], y_range[0]:y_range[1]]

    return lat, lon

def get_mask(x_range, y_range, depth):
    """  Getting Mask given Latitude and Longitude """

    # Get paths
    data_path, mask_path = get_paths()

    # Reading ANHA4 mask data
    mask = nc.Dataset(mask_path + "ANHA4_mask.nc")

    # Extracting mask data for given depth
    tmask = mask['tmask'][0][depth]

    # Extracting mask data for given range
    surf_mask = tmask[x_range[0]:x_range[1], y_range[0]:y_range[1]]

    return surf_mask

def get_data(filename, row_range, col_range, depth=0, var='votemper'):
    """  Getting Data Latitude and Longitude """

    # Get data
    data = nc.Dataset(filename)
    temp = data[var][:]

    # Extracting data given lat-lon selection and depth
    temp = temp[0, depth, row_range[0]:row_range[1], col_range[0]:col_range[1]]

    return temp

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