#!/usr/bin/env python
# coding: utf-8

# Data-related libraries
import matplotlib
import numpy as np

# OS-specific libraries

# Plotting-related libraries
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cartopy import crs as ccrs, feature as cfeature

# Project custom made libraries
import anhalyze.core.anhalyze
import anhalyze.core.anhalyze_geo

# TODO: check best practices for these global variables
# Setting plotting variables
levels = 42
line_levels = 11


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

    if var == 'votemper':
        # cmap = 'coolwarm'
        # vrange = [-20, 20]   # color map based values
        cmap = 'plasma'  # 'magma'
        vrange = [-2, 35]    # Physical based values
    elif var == 'vosaline':
        cmap = 'YlGn'
        vrange = [25, 39]    # Physical based values
    else:
        cmap = 'cividis'
        vrange = None

    if not vrange or color_range == 'relative':
        vrange = [np.nanmin(var_data), np.nanmax(var_data)]
        print(f'vrange: {vrange}')
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


def show_var_data_map(data, idepth=0, var='', color_range='physical', savefig=None):
    """ Displays map of given parameter (var) in lat-lon range and depth.
        Note: depth has not been tested.

        Parameters
        ----------
        var : str
            Variable name.
        data: AnhaDataset
            In AnhaDataset or xarray.Dataset formats.
        idepth: int
            dep level; default first level (0)
        color_range : str
            Color range either `physical` limits, or `relative` values.
        savefig : str
            Filename to save figure including path.
    """
    # TODO add projection options

    assert var in list(data.data_vars), \
        f'[anhalyze_plot_utils] Variable {var} not found in data_vars: {list(data.data_vars)}'

    # TODO better use of init_location, if location set then load it, otherwise calculate it.
    # Init location info
    location_info = anhalyze.core.anhalyze_geo.init_location()

    # Get var data
    var_data = data.data_vars[var].data[0, idepth, :]

    # getting lat and lon
    lat, lon = data.coords[data.attrs['coord_lat']].data, data.coords[data.attrs['coord_lon']].data

    # Set up figure and projection
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1,
                         projection=ccrs.LambertConformal(central_longitude=location_info['central_longitude'],
                                                          standard_parallels=location_info['standard_parallels']))
    ax.set_extent([data.attrs['coord_lon_range'][0],
                   data.attrs['coord_lon_range'][1],
                   data.attrs['coord_lat_range'][0],
                   data.attrs['coord_lat_range'][1]])

    # Adding ocean and land features
    ax.add_feature(get_feature_mask())
    ax.add_feature(get_feature_mask(feature='ocean'))

    # Get var-dependent plotting information
    cmap, vrange = get_plot_config(var, var_data, color_range=color_range)

    # Plotting var data as filled contour regions
    im = ax.contourf(lon, lat, var_data, levels=levels, cmap=cmap,
                     vmin=vrange[0], vmax=vrange[1], transform=ccrs.PlateCarree(), zorder=2)

    # Plotting var data contour lines
    ax.contour(lon, lat, var_data, levels=line_levels, cmap='Greys', linewidths=.2, transform=ccrs.PlateCarree())

    # Create grid-line labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False,
                      y_inline=False, color='k', alpha=.3, linewidths=.01)
    gl.right_labels = gl.top_labels = False

    # Set Color-bar
    axins = inset_axes(ax, width="5%", height="100%", loc='right', borderpad=-1)
    label = '%s [%s]' % (data.data_vars[var].attrs['long_name'].title(), data.data_vars[var].attrs['units'])
    fig.colorbar(im, cax=axins, orientation="vertical", label=label)

    # Save plot if filename is provided
    if savefig:
        fig.savefig(savefig)
