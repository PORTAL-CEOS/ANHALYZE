import netCDF4 as nc
import numpy as np
import seaborn as sns
from cartopy import crs as ccrs
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import anhalyze.core
from anhalyze.core import anhalyze_utils as au
from anhalyze.core.anhalyze_plot_utils import get_feature_mask, levels, cmap, vmin, vmax, line_levels


def show_var_data_map_old(data, lat_range, lon_range, depth=0, var='votemper'):
    """ Displays map of given parameter (var) in lat-lon range and depth.
        Note: depth has not been tested.

       data: ANHA4 data in NCDF format. [nc]
       lat_range: latitute range [tupple]
       lon_range: longitude range [tupple]
       depth: dep level [int]
       var: variable name [str]
    """

    # Init location info
    location_info = anhalyze.core.anhalyze_geo.init_location()

    # Get var data
    var_data = au.get_var_data(data, lat_range, lon_range, depth=depth, var=var)

    # getting lat and lon
    lat, lon = au.get_lat_lon(data, lat_range, lon_range)

    # Set up figure and projection
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1,
                         projection=ccrs.LambertConformal(central_longitude=location_info['central_longitude'],
                                                          standard_parallels=location_info['standard_parallels']))
    ax.set_extent([lon_range[0], lon_range[1], lat_range[0], lat_range[1]])

    # Adding ocean and land features
    ax.add_feature(get_feature_mask())
    ax.add_feature(get_feature_mask(feature='ocean'))

    # Plotting var data as filled contour regions
    im = ax.contourf(lon, lat, var_data, levels=levels, cmap=cmap,
                     vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), zorder=2)

    # Plotting var data contour lines
    ax.contour(lon, lat, var_data, levels=line_levels, cmap='Greys', linewidths=.2, transform=ccrs.PlateCarree())

    # Create grid-line labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False,
                      y_inline=False, color='k', alpha=.3, linewidths=.01)
    gl.right_labels = gl.top_labels = False

    # Set Color-bar
    axins = inset_axes(ax, width="5%", height="100%", loc='right', borderpad=-1)
    label = '%s [%s]' % (data.variables[var].long_name.title(), data.variables[var].units.title())
    fig.colorbar(im, cax=axins, orientation="vertical", label=label, format='%.1f')


def show_var_data_maps(file_list, lat_range, lon_range, depth=0, var='votemper', save_fig=False):
    """ Displays maps of given parameter (var) in lat-lon range and depth, and date selection.
        Note: depth has not been tested.

       file_list: list of ANHA4 data files in NCDF format. [list]
       lat_range: latitute range [tupple]
       lon_range: longitude range [tupple]
       depth: dep level [int]
       var: variable name [str]
    """

    # Setting loop variables.
    lat = None
    lon = None
    data = None

    # Init location info
    # region, proj_size = au.init_location()
    location_info = anhalyze.core.anhalyze_geo.init_location()

    # Get date
    date_start = anhalyze.core.anhalyze.get_date(file_list[0])
    date_end = anhalyze.core.anhalyze.get_date(file_list[-1])

    # Setting up figure
    f_size = 4
    ncols = 3
    nrows = int(np.ceil(len(file_list) / ncols))
    fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=[ncols * f_size, nrows * f_size + 1],
                           subplot_kw={'projection': ccrs.LambertConformal(
                                                          central_longitude=location_info['central_longitude'],
                                                          standard_parallels=location_info['standard_parallels'])})

    # Looping over ANHA4 model files
    for i, xx in enumerate(fig.axes):

        # Stop loop when current subplot number is larger than the number of files in list.
        if i == len(file_list):
            break

        # Get ncdf data
        data = nc.Dataset(file_list[i])

        # Get var data from ncdf data
        var_data = au.get_var_data(data, lat_range, lon_range, depth=depth, var=var)

        # Getting lat and lon
        if i == 0:
            lat, lon = au.get_lat_lon(data, lat_range, lon_range)

        # Plotting data
        # Set map extent
        xx.set_extent([lon_range[0], lon_range[1], lat_range[0], lat_range[1]])

        # Set map features
        xx.add_feature(get_feature_mask())
        xx.add_feature(get_feature_mask(feature='ocean'))

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
        # TODO: change to anha.file_date. If/when moving into class.
        xx.set_title(file_list[i].split('/')[-1].split('_')[1])

        # Set Color-bar
        axins = inset_axes(xx, width="5%", height="100%", loc='right', borderpad=-1)
        label = '%s [%s]' % (data.variables[var].long_name.title(), data.variables[var].units.title())
        fig.colorbar(im, cax=axins, orientation="vertical", label=label, format='%.1f')

    # TODO: generalize for any number of plots. If/when moving into class.
    # Temp fix to plot position given number of plots.
    if len(fig.axes) == 18:
        hspace = .27
    else:
        hspace = .1

    # Add titles to eah figure
    fig.suptitle(data.variables[var].standard_name.replace('_', ' '), fontsize=16)
    # Set up plots positions
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.07, right=.92, wspace=.65, hspace=hspace)

    # Save figure
    if save_fig:
        output_fig_name = '../figs/%s_%s_%s-%s.png' % (location_info['region'],
                                                       data.variables[var].long_name,
                                                       date_start,
                                                       date_end)
        plt.savefig(output_fig_name, bbox_inches='tight', dpi=500)


def plot_timeseries(timeseries_var, data_variables, lat_range, lon_range, var='votemper'):
    """  Displays timeseries of given parameter (var) in lat-lon range and depth, and date selection.
        Note: no depth option available, currently .

       file_list: list of ANHA4 data files in NCDF format. [list]
       lat_range: latitute range [tupple]
       lon_range: longitude range [tupple]
       depth: dep level [int]
       var: variable name [str]  """

    # TODO: update funtion with newer versions found in ipynb reports.
    # Setting up seaborn defaults
    sns.set(rc={'axes.facecolor': 'whitesmoke'})

    # Setting up timeseries plot
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(9, 5))

    # Plotting data with error bars
    ax.errorbar(timeseries_var['date'], timeseries_var['var_mean'], yerr=timeseries_var['var_std'], fmt='o')

    # Labels and axis
    ax.set_ylabel('%s [%s]' % (data_variables[var].long_name.title(), data_variables[var].units))
    ax.set_xlabel('Time')
    ax.set_title('Mean values taken from region: Lat %s , Lon %s' % (str(lat_range), str(lon_range)))
    ax.xaxis.set_tick_params(rotation=30, labelsize=10)

    # Returning to matplotlib defaults
    sns.reset_orig()


def plot_mhw(anhalyzed_timeseries, year=1998, remove_mean=True, show_cat4=False, region="James Bay", mhw=True):
    """
        Plot time series data for given year, along with Marine Heat Wave or Marine Cold Spells curve categories.
        Categories are defined in Hobday et. al. (2008).

        anhalyzed_timeseries: calculated time series for a given region
        year: year [int]
        remove_mean: remove overall mean from data  [bool]
        show_cat4: show Category 4 curve [bool]
        region: region name [str]
        mhw: MHW or MCS [bool]
    """

    # Setting up seaborn defaults
    sns.set(rc={'axes.facecolor': 'whitesmoke'})

    # Selecting year from timeseries
    anhalyzed_timeseries_year = anhalyzed_timeseries[anhalyzed_timeseries['year'] == year].copy()

    # Prep data/plot if removing mean values or not
    if remove_mean:
        labels = ['SST - T$_{c}$', r'$\Delta$T', r'2$\Delta$T', r'3$\Delta$T', r'4$\Delta$T']

        # Get freezing line
        freezing_line = anhalyzed_timeseries_year['var_mean_mean'] * -1

        # Removing the yearly average T$_{c}$
        anhalyzed_timeseries_year['var_mean'] = anhalyzed_timeseries_year.apply(
            lambda row: row.var_mean - row.var_mean_mean, axis=1)
        anhalyzed_timeseries_year['var_mean_quantile'] = anhalyzed_timeseries_year.apply(
            lambda row: row.var_mean_quantile - row.var_mean_mean, axis=1)
        anhalyzed_timeseries_year['var_mean_2T'] = anhalyzed_timeseries_year.apply(
            lambda row: row.var_mean_2T - row.var_mean_mean, axis=1)
        anhalyzed_timeseries_year['var_mean_3T'] = anhalyzed_timeseries_year.apply(
            lambda row: row.var_mean_3T - row.var_mean_mean, axis=1)
        anhalyzed_timeseries_year['var_mean_4T'] = anhalyzed_timeseries_year.apply(
            lambda row: row.var_mean_4T - row.var_mean_mean, axis=1)

    else:
        if mhw:
            labels = ['SST', 'T$_{90}$', r'T$_{c}$+2$\Delta$T', r'T$_{c}$+3$\Delta$T', r'T$_{c}$+4$\Delta$T']
        else:
            labels = ['SST', 'T$_{90}$', r'T$_{c}$-2$\Delta$T', r'T$_{c}$-3$\Delta$T', r'T$_{c}$-4$\Delta$T']
        freezing_line = anhalyzed_timeseries_year['var_mean_mean'] * 0. - 2.

    # Setting up MHW/MCS related plotting variables.
    if mhw:
        colors = ['gold', 'orange', 'red', 'maroon']
        alphas = [0.2, 0.3, 0.1]
        mhw_title = 'MHW'
    else:
        colors = ['DeepSkyBlue', 'dodgerblue', 'blue', 'navy']
        alphas = [0.1, 0.5, 0.1]
        mhw_title = 'MCS'

    # Setting up figure
    plt.figure(figsize=(9, 5))

    # Plotting data and categories
    if not remove_mean:
        plt.plot(anhalyzed_timeseries['date'],
                 anhalyzed_timeseries.var_mean_mean,
                 c='gray', alpha=.5, label='T$_{c}$')
    else:
        plt.plot(anhalyzed_timeseries['date'],
                 anhalyzed_timeseries.var_mean_mean*0,
                 c='gray', alpha=.5, label='T$_{c}$')

    plt.scatter(anhalyzed_timeseries_year['date'],
                anhalyzed_timeseries_year.var_mean,
                c='k', alpha=.3, label=labels[0])

    plt.plot(anhalyzed_timeseries_year['date'],
             anhalyzed_timeseries_year.var_mean_quantile,
             c=colors[0], alpha=alphas[1], label=labels[1])
    plt.plot(anhalyzed_timeseries_year['date'],
             anhalyzed_timeseries_year.var_mean_2T,
             c=colors[1], alpha=alphas[1], label=labels[2])
    plt.plot(anhalyzed_timeseries_year['date'],
             anhalyzed_timeseries_year.var_mean_3T,
             c=colors[2], alpha=alphas[1], label=labels[3])
    if show_cat4:
        plt.plot(anhalyzed_timeseries_year['date'],
                 anhalyzed_timeseries_year.var_mean_4T,
                 c=colors[3], alpha=alphas[1], label=labels[4])

    # Fill in categories colors
    plt.fill_between(anhalyzed_timeseries_year['date'],
                     anhalyzed_timeseries_year.var_mean_quantile,
                     y2=anhalyzed_timeseries_year.var_mean_2T,
                     alpha=alphas[2], facecolor=colors[0])

    plt.fill_between(anhalyzed_timeseries_year['date'],
                     anhalyzed_timeseries_year.var_mean_2T,
                     y2=anhalyzed_timeseries_year.var_mean_3T,
                     alpha=alphas[2], facecolor=colors[1])

    if show_cat4:
        plt.fill_between(anhalyzed_timeseries_year['date'],
                         anhalyzed_timeseries_year.var_mean_3T,
                         y2=anhalyzed_timeseries_year.var_mean_4T,
                         alpha=.1, facecolor=colors[2])

    if remove_mean:
        plt.fill_between(anhalyzed_timeseries_year['date'], freezing_line, alpha=alphas[0],
                         y2=freezing_line-20,
                         label='Freezing Zone', hatch='/', edgecolor='navy')
    else:
        plt.fill_between(anhalyzed_timeseries_year['date'], freezing_line, alpha=alphas[0],
                         label='Freezing Zone', hatch='/', edgecolor='navy')

    # Set axis limits and labels
    if remove_mean:
        if mhw:
            plt.ylim(bottom=0)
        else:
            plt.ylim(top=2, bottom=-12)
    else:
        plt.ylim(bottom=-2)

    plt.xlim([anhalyzed_timeseries_year['date'].iloc[0], anhalyzed_timeseries_year['date'].iloc[-1]])
    plt.ylabel(' Residual\n Temperature [\N{DEGREE SIGN}C]', fontsize=15)
    plt.title('%s %i %s' % (region, year, mhw_title), fontsize=18)
    plt.xticks(fontsize=14, rotation=30)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)

    plt.tight_layout()
    plt.show()

    # Returning to matplotlib defaults
    sns.reset_orig()
