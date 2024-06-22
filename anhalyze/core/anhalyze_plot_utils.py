#!/usr/bin/env python
# coding: utf-8

# Data-related libraries
import matplotlib
import numpy as np
import os

# Plotting-related libraries
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cartopy import crs as ccrs, feature as cfeature

# Project custom made libraries
import anhalyze.core.anhalyze
import anhalyze.core.anhalyze_geo

# Setting plotting variables as global constants for now
LEVELS = 42
LINE_LEVELS = 11


def get_plot_config(var, var_data, color_range='physical'):
    """ Return var-dependent plotting information

        Parameters
        ----------
        var : str
            Variable name.
        var_data : np.array
            numpy array with var data.
        color_range : str
            Color range either `physical` limits, or `relative` values.
    """

    if var == 'votemper':  # Temperature
        # cmap = 'coolwarm'
        # vrange = [-20, 20]   # color map based values
        cmap = 'magma'  # Other possible colors: 'plasma', 'magma'
        vrange = [-2, 35]    # Physical based values
    elif var == 'vosaline':  # Salinity
        cmap = 'viridis'  # Other possible colors: 'winter'
        vrange = [25, 39]    # Physical based values
    elif var == 'chl':  # Chlorophyll
        cmap = 'BuGn'
        vrange = None  # Placeholder for physical based values
    else:
        # cmap = 'cividis'
        cmap = 'spring'
        vrange = None

    if not vrange or color_range == 'relative':
        vrange = [np.nanmin(var_data), np.nanmax(var_data)]
        print(f'  vrange: {vrange}')
    else:
        vrange = vrange

    return cmap, vrange


def get_feature_mask(feature='land', resolution='50m'):
    """   """

    # Select face color
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



def get_projection(proj='LambertConformal', proj_info=None):
    """
    Select Cartopy projections option and configurations based on users choice
    of projection and coordinates info from `AnhaDataset`.
    
    Parameters
    ----------
    proj : str
        Projection name from Cartopy list.
    proj_info : dict
        Information for projection calculated by get_projection_info.
    """

    if proj_info is None:
        raise "Argument proj_info is None. Use get_projection_info do obtain projection information from `AnhaDataset`."

    proj_list = {'PlateCarree': ccrs.PlateCarree(central_longitude=proj_info['central_longitude']),
                 'LambertAzimuthalEqualArea': ccrs.LambertAzimuthalEqualArea(
                     central_longitude=proj_info['central_longitude'],
                     central_latitude=proj_info['central_latitude']),
                 'AlbersEqualArea': ccrs.AlbersEqualArea(
                     central_longitude=proj_info['central_longitude'],
                     central_latitude=proj_info['central_latitude'],
                     standard_parallels=proj_info['standard_parallels']),
                 'NorthPolarStereo': ccrs.NorthPolarStereo(central_longitude=proj_info['central_longitude']),
                 'Orthographic': ccrs.Orthographic(central_longitude=proj_info['central_longitude'],
                                                   central_latitude=proj_info['central_latitude']),
                 'Robinson': ccrs.Robinson(central_longitude=proj_info['central_longitude']),
                 'LambertConformal': ccrs.LambertConformal(central_longitude=proj_info['central_longitude'],
                                                           standard_parallels=proj_info['standard_parallels']),
                 'Mercartor': ccrs.Mercator(central_longitude=proj_info['central_longitude'],
                                            min_latitude=proj_info['lat_range'][0],
                                            max_latitude=proj_info['lat_range'][1]),
                 'AzimuthalEquidistant': ccrs.AzimuthalEquidistant(central_longitude=proj_info['central_longitude'],
                                                                   central_latitude=proj_info['central_latitude']),
                 }

    assert proj in list(proj_list.keys()), \
        f'[anhalyze_plot_utils] Projection {proj} not found in list of projections available: {list(proj_list.keys())}'

    proj_config = proj_list[proj]

    return proj_config


def get_projection_info(data):
    """
    Calculate information used to set figure projection.
    
    Parameters
    ----------
    data: `AnhaDataset`
        Dataset to be plotted
    """

    # Setting up user's region
    east = data.attrs['coord_lon_range'][1]
    west = data.attrs['coord_lon_range'][0]
    north = data.attrs['coord_lat_range'][1]
    south = data.attrs['coord_lat_range'][0]

    # 1/6th law to calculate the standard parallels. We calculate 1/6 of the
    # distance in degrees from south to north, then add it to the southern
    # figure limit and subtract from the northern limit.
    #
    law16 = (north - south) / 6
    standard_parallels = (south + law16, north - law16)

    # Central longitude and latitude are halfway from east to west
    # and from south to north, respectively
    # Central longitude
    midway_lon = (east - west) / 2
    central_longitude = east - midway_lon

    # Central latitude
    midway_lat = (north - south) / 2
    central_latitude = north - midway_lat

    lat_range = (south, north)
    lon_range = (west, east)

    proj_info = {'lat_range': lat_range,
                 'lon_range': lon_range,
                 'standard_parallels': standard_parallels,
                 'central_longitude': central_longitude,
                 'central_latitude': central_latitude,
                 }

    return proj_info


def show_var_data_map(data, var='', idepth=0, proj='', color_range='physical', savefig=None):
    """ Displays map of given parameter (var) in lat-lon range and depth.
        Note: depth has not been tested.

        Parameters
        ----------
        data: AnhaDataset
            In `AnhaDataset` or `xarray.Dataset` formats.
        var : str
            Variable name.
        idepth: int
            dep level; default first level (0)
        color_range : str
            Color range either `physical` limits, or `relative` values.
        savefig : str
            Filename to save figure including path.
    """

    assert var in list(data.data_vars), \
        f'[anhalyze_plot_utils] Variable {var} not found in data_vars: {list(data.data_vars)}'

    # Calculate projection information (e.g. Standard parallels) based on the dataset lat and lon limits	
    proj_info = get_projection_info(data)

    # Select figure projection	 
    proj_config = get_projection(proj, proj_info)

    # Get var data
    if len(data.data_vars[var].dims) == 4:
        if 'depth' in data.data_vars[var].dims[3]:
            var_data = data.data_vars[var].data[0, :, :, 0]
        else:
            var_data = data.data_vars[var].data[0, idepth, :]
    elif len(data.data_vars[var].dims) == 3:
        var_data = data.data_vars[var].data[0, :]
    else:
        raise ValueError(f"Variable {var} should be a 2D or 3D variable.")

    # getting lat and lon
    lat, lon = data.coords[data.attrs['coord_lat']].data, data.coords[data.attrs['coord_lon']].data

    # Set up figure and projection
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1,
                         projection=proj_config)

    ax.set_extent([data.attrs['coord_lon_range'][0],
                   data.attrs['coord_lon_range'][1],
                   data.attrs['coord_lat_range'][0],
                   data.attrs['coord_lat_range'][1]],
                  crs=ccrs.PlateCarree(),
                  )

    # Adding ocean and land features
    ax.add_feature(get_feature_mask())
    ax.add_feature(get_feature_mask(feature='ocean'))

    # Get var-dependent plotting information
    cmap, vrange = get_plot_config(var, var_data, color_range=color_range)

    # Plotting var data as filled contour regions
    im = ax.contourf(lon, lat, var_data, levels=LEVELS, cmap=cmap, extend='both',
                     vmin=vrange[0], vmax=vrange[1], transform=ccrs.PlateCarree(), zorder=2)

    # Plotting var data contour lines
    ax.contour(lon, lat, var_data, levels=LINE_LEVELS, cmap='Greys', linewidths=.2, transform=ccrs.PlateCarree())

    # Create grid-line labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False,
                      y_inline=False, color='k', alpha=.3, linewidths=.01)
    gl.right_labels = gl.top_labels = False

    # Set Color-bar
    axins = inset_axes(ax, width="2.5%", height="100%", loc='right', borderpad=-1)
    label = '%s [%s]' % (data.data_vars[var].attrs['long_name'].title(), data.data_vars[var].attrs['units'])
    fig.colorbar(im, cax=axins, orientation="vertical", label=label, extend='both')

    # Save plot if filename is provided
    if savefig:
        # Including path from original file if not given
        if not os.path.dirname(savefig):
            savefig = os.path.join(data.attrs["filepath"], savefig)

        print(f'[Anhalyze] Saving figure: {savefig}')

        fig.savefig(savefig)
