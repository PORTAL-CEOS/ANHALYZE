#!/usr/bin/env python
# coding: utf-8

# General libraries
import datetime
import os

import netCDF4 as nc
import numpy as np
import pandas as pd
import deprecation

# project specific
import anhalyze.core.anhalyze_plot_utils_depr
from anhalyze.core import anhalyze_plot_utils as apu
from anhalyze.core.anhalyze_utils import calc_stats_var_data, get_paths
from anhalyze.core.anhalyze import get_date


# Possible names AnhalyzeProject, AnhaProject, AnhalyzeTimeseries, AnhaTimeseries
@deprecation.deprecated()
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


@deprecation.deprecated()
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


@deprecation.deprecated()
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


@deprecation.deprecated()
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


@deprecation.deprecated()
def plot_timeseries(timeseries_var, data_variables, lat_range, lon_range, var='votemper'):
    """  Wrapper to plot timeseries function  """

    anhalyze.core.anhalyze_plot_utils_depr.plot_timeseries(timeseries_var, data_variables, lat_range, lon_range, var=var)


@deprecation.deprecated()
def show_var_data_maps(file_list, lat_range, lon_range, depth=0, var='votemper'):
    """ Displays map of given var in lat-lon range and depth.
        Note: depth has not been tested.
    """

    anhalyze.core.anhalyze_plot_utils_depr.show_var_data_maps(file_list, lat_range, lon_range, depth=depth, var=var)


@deprecation.deprecated()
def plot_mhw(anhalyzed_timeseries, year=1998, remove_mean=True, show_cat4=False, region="James Bay", mhw=True):
    """  Wrapper to plot mhw function  """

    anhalyze.core.anhalyze_plot_utils_depr.plot_mhw(anhalyzed_timeseries, year=year, remove_mean=remove_mean, show_cat4=show_cat4, region=region, mhw=mhw)


@deprecation.deprecated()
def find_mhw_info(anhalyzed_timeseries, mhw=True):
    """  placeholder, to find dates and deltaTs/categories on mhws or mcs """

#    if mhw:

    return None
