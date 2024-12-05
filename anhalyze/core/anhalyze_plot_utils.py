#!/usr/bin/env python
# coding: utf-8

# Data-related libraries
import matplotlib
import numpy as np
import os

# Plotting-related libraries
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean.cm as cmo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cartopy import crs as ccrs, feature as cfeature

# Project custom made libraries

# Setting plotting variables as global constants for now
LEVELS = 21
LINE_LEVELS = 11


def get_plot_config(var, var_data, attrs, color_range='default'):
    """ Return var-dependent plotting information

        Parameters
        ----------
        var : str
            Variable name.
        var_data : np.array
            numpy array with var data.
        attrs : dict
            Attributes from `AnhaDataset`
        color_range : str
            Color range either `default` limits, `local` data values or a two items list [vmin, vmax].
            Color range options:
             default: Limits decide by Anhalyze developers. It is based on the more
                      likely limits the user can find in a ANHA4 outputs.
             local: Color range based on the values within the area selected by the user.
             [vmin, vmax] = List of color range limits chosen by the user.
             
    """

    color_range_options = ['default', 'local']
    
    assert_message = f"[anhalyze_plot_utils] Color range option '{color_range}' "
    assert_message += f"not found. Color range should be either '{color_range_options[0]}'," \
                      f" '{color_range_options[1]}', or a two-item list ([vmin, vmax])."
    assert color_range in color_range_options or isinstance(color_range, list), assert_message

    # Selection of cmap and vrange given var.
    if var == 'votemper':  # Temperature
        cmap = cmo.thermal  # Other possible colors: 'plasma', 'magma'
        vrange = [-2, 28]    # Default values
    elif var == 'vosaline':  # Salinity
        cmap = cmo.haline  # Other possible colors: 'winter'
        vrange = [20, 40]    # Default based values
    elif var == 'ileadfra':  # Sea ice concentration
        cmap = cmo.ice
        vrange = [0, 1]  # Default based values
    elif var == 'chl':  # Chlorophyll
        cmap = cmo.algae
        vrange = [10, 1000]  # Placeholder for default based values
    elif (attrs['grid'] in ['gridU', 'gridV']) or (var in ['iicevelu', 'iicevelv']):
        vrange = [-1.5, 1.5]  # Default based values
        cmap = cmo.balance
    elif attrs['grid'] == 'icebergs':
        cmap = cmo.thermal
        vrange = [0.00000001, 0.001]
    else:
        cmap = 'spring'
        vrange = None
        
    #
    if not vrange or color_range == 'local':
        if cmap == cmo.balance:
            vdistmax = np.nanmax(np.abs(var_data))
            vrange = [-vdistmax, vdistmax]
        else:    
            vrange = [np.nanmin(var_data), np.nanmax(var_data)]
            print(f'  vrange: {vrange}')
    else:
        vrange = vrange
        
    # Color range manually input by user
    if isinstance(color_range, list):
        vrange = sorted(color_range)
        
    # Colorbar boundaries based on the vranges and LEVELS
    if attrs['grid'] == 'icebergs' or var == 'chl':
        if 0 in vrange:
            print('[Anhalyze] Local variable range are either min or max equal to 0.\n '
                  'The value was replaced by the data value closest to 0.')
            newv = np.nanmin(np.abs(var_data[np.nonzero(var_data)]))[0]
            print(f'[Anhalyze] New vmin: {newv}')
            cnorm = mcolors.LogNorm(vmin=newv, vmax=vrange[1])
        else:
            cnorm = mcolors.LogNorm(vmin=vrange[0], vmax=vrange[1])
    else:
        bounds = np.linspace(vrange[0], vrange[1], LEVELS)
        cnorm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)

    return cmap, vrange, cnorm


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


def get_projection(proj_name='LambertConformal', proj_info=None):
    """
    Select Cartopy projections option and configurations based on users choice
    of projection and coordinates info from `AnhaDataset`.
    
    Parameters
    ----------
    proj_name : str
        Projection name from Cartopy list. The projections available are: 'PlateCarree', 'LambertAzimuthalEqualArea',
        'AlbersEqualArea', 'NorthPolarStereo', 'Orthographic', 'Robinson', 'LambertConformal',
        'Mercator', and 'AzimuthalEquidistant'.
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
                 'Robinson': ccrs.Robinson(central_longitude=0),
                 'LambertConformal': ccrs.LambertConformal(central_longitude=proj_info['central_longitude'],
                                                           standard_parallels=proj_info['standard_parallels']),
                 'Mercator': ccrs.Mercator(central_longitude=0,
                                           min_latitude=proj_info['lat_range'][0],
                                           max_latitude=proj_info['lat_range'][1]),
                 'AzimuthalEquidistant': ccrs.AzimuthalEquidistant(central_longitude=proj_info['central_longitude'],
                                                                   central_latitude=proj_info['central_latitude']),
                 }

    assert_message = f'[anhalyze_plot_utils] Projection {proj_name} '
    assert_message += f'not found in list of projections available: {list(proj_list.keys())}'
    assert proj_name in list(proj_list.keys()), assert_message

    proj_config = proj_list[proj_name]

    return proj_config


def get_projection_info(attrs):
    """ Calculate information used to set figure projection.
    
        Parameters
        ----------
        attrs : dict
            Attributes from `AnhaDataset`
    """

    # Setting up user's region
    east = attrs['coord_lon_range'][1]
    west = attrs['coord_lon_range'][0]
    north = attrs['coord_lat_range'][1]
    south = attrs['coord_lat_range'][0]

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


def set_colorbar_ticks(cbar, var_data, color_range, vrange):
    """ Set colorbar ticks adjusted to each type of color_range in Anhalyze.
    
        Parameters
        ----------
        cbar : object
            Colorbar object.
        var_data : AnhaDataset
            AnhaDataset to be plotted.
        color_range : str
            Color range either `default` limits, or `local` values.
        vrange : list
            Range of values returned from get_plot_config
    """
    
    if color_range == 'default':
            
        if np.nanmin(var_data) >= vrange[0] and np.nanmax(var_data) >= vrange[1]:
            step = (vrange[1] - vrange[0])/(LEVELS-1)
            cbar.ax.set_ylim(np.nanmin(var_data), vrange[1])
            cbar.ax.set_yticks(np.arange(vrange[0], np.nanmax(var_data)+step/2, step*2))
            
        elif np.nanmin(var_data) <= vrange[0] and np.nanmax(var_data) <= vrange[1]:
            step = (vrange[1] - vrange[0])/(LEVELS-1)
            cbar.ax.set_ylim(vrange[0], np.nanmax(var_data))
            cbar.ax.set_yticks(np.arange(vrange[0], np.nanmax(var_data)+step/2, step*2))
            
        elif np.nanmin(var_data) >= vrange[0] and np.nanmax(var_data) <= vrange[1]:
            
            cbar.ax.set_ylim(np.nanmin(var_data), np.nanmax(var_data))
                
    else:
        step = (vrange[1] - vrange[0])/(LEVELS-1)
        cbar.ax.set_yticks(np.arange(vrange[0], vrange[1]+step/2, step*2))
            
    return cbar


def show_var_data_map(var_da, attrs, color_range='default', savefig=None, proj_name=''):
    """ Displays map of given parameter (var) in lat-lon range and depth.

        Parameters
        ----------
        var_da: xarray.DataArray
            xarray.DataArray for given var(var_da.name).
        attrs : dict
            Attributes from `AnhaDataset`
        color_range : str
            Color range either `default` limits, `local` data values or a two
             items list [vmin, vmax].
             
            Color range options:
                 default: Color range limits predefined based on the more
                 likely values the user can find in a ANHA4 outputs.
             
                 local: Color range based on the values within the area selected by the user.
                 
                 [vmin, vmax] = List of color range limits chosen by the user.
        savefig : str
            Filename to save figure including path.
        proj_name : str
            Projection name from Cartopy list. The projections available are: 'PlateCarree',
            'LambertAzimuthalEqualArea', 'AlbersEqualArea', 'NorthPolarStereo', 'Orthographic', 'Robinson',
            'LambertConformal', 'Mercator', and 'AzimuthalEquidistant'.
        """

    # Calculate projection information (e.g. Standard parallels) based on the dataset lat and lon limits
    proj_info = get_projection_info(attrs)

    # Select figure projection
    proj_config = get_projection(proj_name, proj_info)

    # getting lat and lon
    lat, lon = np.squeeze(var_da.coords[attrs['coord_lat']].data), np.squeeze(var_da.coords[attrs['coord_lon']].data)

    # Get var data to 2D for plotting.
    var_data = np.squeeze(var_da.data)

    # Set up figure and projection
    fig = plt.figure(num=var_da.name)
    ax = fig.add_subplot(1, 1, 1,
                         projection=proj_config)

    #
    ax.set_extent([attrs['coord_lon_range'][0],
                   attrs['coord_lon_range'][1],
                   attrs['coord_lat_range'][0],
                   attrs['coord_lat_range'][1]],
                  crs=ccrs.PlateCarree(),
                  )

    # Adding ocean and land features
    ax.add_feature(get_feature_mask(), zorder=1)
    ax.add_feature(get_feature_mask(feature='ocean'), zorder=0)

    # Get var-dependent plotting information
    cmap, vrange, cnorm = get_plot_config(var_da.name, var_data, attrs, color_range=color_range)

    # Plotting var data as filled contour regions
    im = plt.pcolormesh(lon, lat, var_data, cmap=cmap,
                        norm=cnorm, transform=ccrs.PlateCarree(), zorder=2)
    
    # Plotting var data contour lines
    ax.contour(lon, lat, var_data, levels=LINE_LEVELS, cmap='Greys', linewidths=.2, transform=ccrs.PlateCarree())

    # Create grid-line labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False,
                      y_inline=False, color='k', alpha=.3)
    gl.right_labels = gl.top_labels = False

    # Set Color-bar
    # Setting color bar feature
    if color_range != 'default':
        bar_extend = 'both'
    else:
        bar_extend = None

    # Color-bar axis and properties
    axins = inset_axes(ax, width="3%", height="100%", loc='right', borderpad=-1)
    label = '%s [%s]' % (var_da.attrs['long_name'].title(), var_da.attrs['units'])
    cbar = fig.colorbar(im, cax=axins, orientation="vertical", label=label, extend=bar_extend)
    
    # Set colorbar ticks
    if attrs['grid'] != 'icebergs' and var_da.name != 'chl':
        new_cbar = set_colorbar_ticks(cbar, var_data, color_range, vrange)
        cbar = new_cbar
    
    # Display map
    if get_ipython().__class__.__name__ != 'ZMQInteractiveShell':
        plt.ion()
        fig.show()

    # Save plot if filename is provided
    if savefig:
        # Including path from original file if not given
        if not os.path.dirname(savefig):
            savefig = os.path.join(attrs["filepath"], savefig)

        print(f'[Anhalyze] Saving figure: {savefig}')

        fig.savefig(savefig)
