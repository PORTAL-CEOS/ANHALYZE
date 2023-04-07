#!/usr/bin/env python
# coding: utf-8

# Data-related libraries
import numpy as np
import netCDF4 as nc

# OS-specific libraries
from sys import platform

# Plotting-related libraries
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cartopy import crs as ccrs
import seaborn as sns

# Project custom made libaries
import anha4_utils as au


def show_var_data_map(data, lat_range, lon_range, depth=0, var='votemper'):
    """ Displays map of given var in lat-lon range and depth.
        Note: depth has not been tested.
    """

    # Get var data
    var_data = au.get_var_data(data, lat_range, lon_range, depth=depth, var=var)

    # getting lat and lon
    lat, lon = au.get_lat_lon(data, lat_range, lon_range)

    # EE: should move this elsewhere
    # Setting plotting vars.
    levels = 42
    line_levels = 11
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
    ax.add_feature(au.get_feature_mask())
    ax.add_feature(au.get_feature_mask(feature='ocean'))

    # Plotting var data as filled contour regions
    im = ax.contourf(lon, lat, var_data, levels=levels, cmap=cmap,
                     vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), zorder=2)

    # Plotting data contour lines
    ax.contour(lon, lat, var_data, levels=line_levels, cmap='Greys', linewidths=.2, transform=ccrs.PlateCarree())

    # Create grid-line labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False,
                      y_inline=False, color='k', alpha=.3, linewidths=.01)
    gl.right_labels = gl.top_labels = False

    # Set Color-bar
    axins = inset_axes(ax, width="5%", height="100%", loc='right', borderpad=-1)
    label = '%s [%s]' % (data.variables[var].long_name.title(), data.variables[var].units.title())
    fig.colorbar(im, cax=axins, orientation="vertical", label=label, format='%.1f')


def show_var_data_maps(file_list, lat_range, lon_range, depth=0, var='votemper'):
    """ Displays map of given var in lat-lon range and depth.
        Note: depth has not been tested.
    """

    # Setting plotting vars.
    levels = 42
    line_levels = 11
    vmax = 20.
    vmin = -20.
    standard_parallels = (52.5, 62.5)
    central_longitude = -80
    cmap = 'coolwarm'

    region, proj_size = au.init_location()

    # Get
    date_start = au.get_date(file_list[0])
    date_end = au.get_date(file_list[-1])

    # Setting up figure
    f_size = 4
    ncols = 3
    nrows = int(np.ceil(len(file_list) / ncols))
    fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=[ncols * f_size, nrows * f_size + 1],
                           subplot_kw={'projection': ccrs.LambertConformal(central_longitude=central_longitude,
                                                                           standard_parallels=standard_parallels)})

    # Looping over files
    for i, xx in enumerate(fig.axes):

        # Stop loop when number of subplots is larger than the number of files.
        if i == len(file_list):
            break

        # Get data
        data = nc.Dataset(file_list[i])

        # Get var data
        var_data = au.get_var_data(data, lat_range, lon_range, depth=depth, var=var)

        # getting lat and lon
        if i == 0:
            lat, lon = au.get_lat_lon(data, lat_range, lon_range)

        # Plotting data
        # Set map extent
        xx.set_extent([lon_range[0], lon_range[1], lat_range[0], lat_range[1]])

        # Set map features
        xx.add_feature(au.get_feature_mask())
        xx.add_feature(au.get_feature_mask(feature='ocean'))

        # Plotting var data as filled contour regions
        im = xx.contourf(lon, lat, var_data, levels=levels, cmap=cmap,
                         vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), zorder=2)

        # Plotting data contour lines
        xx.contour(lon, lat, var_data, levels=line_levels, cmap='Greys', linewidths=.2, transform=ccrs.PlateCarree())

        # Create grid-line labels
        gl = xx.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False,
                          y_inline=False, color='k', alpha=.3, linewidths=.01)
        gl.right_labels = gl.top_labels = False

        # Set axis labels
        # To show x label only at the bottom most row
        if i > ncols*(nrows-1)-1:
            xx.text(0.5, -0.22, 'Longitude', va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=xx.transAxes)
        # To show y label only at left most column
        if i % ncols == 0:
            xx.text(-0.25, 0.55, 'Latitude', va='bottom', ha='center',
                    rotation='vertical', rotation_mode='anchor',
                    transform=xx.transAxes)

        # Set title
        xx.set_title(file_list[i].split('/')[-1].split('_')[1])

        # Set Color-bar
        axins = inset_axes(xx, width="5%", height="100%", loc='right', borderpad=-1)
        label = '%s [%s]' % (data.variables[var].long_name.title(), data.variables[var].units.title())
        fig.colorbar(im, cax=axins, orientation="vertical", label=label, format='%.1f')

    # Temp fix to plot position given number of plots.
    if len(fig.axes) == 18:
        hspace = .27
    else:
        hspace = .1

    fig.suptitle(data.variables[var].standard_name.replace('_', ' '), fontsize=16)
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.07, right=.92, wspace=.65, hspace=hspace)

    # Save figure to file.
    if platform == "linux" or platform == "linux2":
        output_fig_name = '../figs/%s_%s_%s-%s.png' % (region, data.variables[var].long_name, date_start, date_end)
        plt.savefig(output_fig_name, bbox_inches='tight', dpi=500)


def plot_timeseries(timeseries_var, data_variables, lat_range, lon_range, var='votemper'):
    """    """

    # Setting up seaborn defaults
    sns.set()

    # Setting up timeseries plot
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(9, 5))

    # Plotting data with error bars
    ax.errorbar(timeseries_var['date'], timeseries_var['var_mean'], yerr=timeseries_var['var_std'], fmt='o')

    # Labels and axis
    ax.set_ylabel('%s [%s]'%(data_variables[var].long_name.title(), data_variables[var].units))
    ax.set_xlabel('Time')
    ax.set_title('Mean values taken from region: Lat %s , Lon %s' % (str(lat_range), str(lon_range)))
    ax.xaxis.set_tick_params(rotation=30, labelsize=10)

    # Returning to matplotlib defaults
    sns.reset_orig()