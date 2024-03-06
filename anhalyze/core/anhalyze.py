#!/usr/bin/env python3
# coding: utf-8

# System-related libraries
import os
#import pdb

# Data-related libraries
import netCDF4 as nc

# Project-related libraries


# Possible other names AnhaModelData, Dataset, AnhaData, AnhaDataframe, AnhaReader
class AnhaDataset:
    """ This wrapper class opens netCDF4 files with ANHA specific properties.

    Attributes
    ----------
        filename : str
            Filename given in format ANHA?-??????_y????m??d??_grid?.nc

    Methods
    -------


    """

    def __init__(self, filename, load_data=True, cartesian=True):
        """ Initializing class.

        Parameters
        ----------
        filename : str
            Filename given in format ANHA?-??????_y????m??d??_grid?.nc
            or *_mask*.nc
        load_data : bool, optional
            Bool for loading data (Default is False)
        cartesian : bool, optional
            Bool for using cartesian coordinates (Default is True)
            Otherwise use geographical coordinates and SI units

        """

        # Setting up filename
        assert os.path.isfile(filename)
        if os.path.dirname(filename):
            self.filename = os.path.basename(filename)
            self.filepath = os.path.dirname(filename)
        else:
            self.filename = filename
            self.filepath = ''

        # Initialize data properties
        # TODO update to use xarray to be able to load file info without loading data
        if load_data:
            self.anha_data = nc.Dataset(os.path.join(self.filepath, self.filename))
            # self.depth = self.anha_data.dimensions['deptht'].size
            # self.xdim = self.anha_data.dimensions['x'].size
            # self.ydim = self.anha_data.dimensions['x'].size

            #TODO need to update, placeholder for now
            self.nc_variables = self.anha_data.variables
            self.nc_dimensions = self.anha_data.dimensions
            self.var = 'votemper'

        else:
            self.anha_data = None

        # Initialize grid type from filename
        self._init_grid()

        # Initialize model properties from filename
        if '_grid' in self.filename:
            self.model_run = self.filename.split('_')[0]
            self.year = self.filename.split('y')[-1][:4]
            self.month = self.filename.split('m')[-1][:2]
            self.day = self.filename.split('d')[1][:2]

            assert 'ANHA' in self.model_run, f'Filename format not recognized: {self.model_run}'

            self.configuration = self.filename.split('-')[0]

        else:
            pass

        # Initialize unit properties
        self.cartesian = cartesian

        # Initialize selection
        self._setup_selection_range(init=True)

    def _init_grid(self):
        """Setup grid related values
        """

        if '_grid' in self.filename:
            # Init grid type
            self.grid = self.filename.split('_grid')[-1][0]
            grid_values = 'TBUVW'
            assert self.grid in grid_values, 'File type not recognized'
            self.grid = 'grid'+self.grid
            self.is_mask = False

            # Init grid dimensions var names
            self.depth_var_name = [var for var in list(self.nc_dimensions.keys()) if 'depth' in var][0]
            # Need index 1 below, since dimensions has also key 'axis_nbounds'
            self.x_var_name = [var for var in list(self.nc_dimensions.keys()) if 'x' in var][1]
            self.y_var_name = [var for var in list(self.nc_dimensions.keys()) if 'y' in var][0]

            # Init grid geocoordinates var names
            self.lat_var_name = [var for var in list(self.nc_variables.keys()) if 'nav_lat' in var][0]
            self.lon_var_name = [var for var in list(self.nc_variables.keys()) if 'nav_lon' in var][0]

            # TODO question: is the grid_W in the gridT files the same as the one in the gridW files?

        elif '_mask' in self.filename:
            # Init grid type
            self.grid = 'mask'
            self.is_mask = True

            # Init grid dimensions var names
            self.x_var_name = [var for var in list(self.nc_dimensions.keys()) if 'x' in var][0]
            self.y_var_name = [var for var in list(self.nc_dimensions.keys()) if 'y' in var][0]
            self.depth_var_name = [var for var in list(self.nc_dimensions.keys()) if 'z' in var][0]

            # Init grid geocoordinates var names
            self.lat_var_name = [var for var in list(self.nc_variables.keys()) if 'nav_lat' in var][0]
            self.lon_var_name = [var for var in list(self.nc_variables.keys()) if 'nav_lon' in var][0]

        else:
            self.grid = ''
            self.is_mask = None

    def _setup_selection_range(self, lat_range=None, lon_range=None, i_range=None, j_range=None, init=False):
        """Setup data selection range,


        """
        # TODO could, set an absolut minimum for lat at: -20.07611 , or min in file.

        if init:

            # For cartesian coordinates
            self.i_begin = 0
            self.i_end = self.nc_dimensions[self.x_var_name].size
            self.j_begin = 0
            self.j_end = self.nc_dimensions[self.y_var_name].size
            self.k_begin = 0
            self.k_end = self.nc_dimensions[self.depth_var_name].size

            self.i_range = [self.i_begin, self.i_end]
            self.j_range = [self.j_begin, self.j_end]
            self.k_range = [self.k_begin, self.k_end]

            # For geographical coordinates
            self.lat_range = [self.anha_data[self.lat_var_name][:].min(),
                              self.anha_data[self.lat_var_name][:].max()]
            self.lon_range = [self.anha_data[self.lon_var_name][:].min(),
                              self.anha_data[self.lon_var_name][:].max()]
            self.depth_range = [self.anha_data[self.depth_var_name][:].min(),
                                self.anha_data[self.depth_var_name][:].max()]

        else:

            #TODO add assest if ranges are correct

            if lat_range:
                self.lat_range = lat_range
            if lon_range:
                self.lon_range = lon_range

            if i_range:
                self.i_range = i_range
            if j_range:
                self.j_range = j_range

    def show_var_data_map(self, var=''):
        """ Displays map of given var.
        """
        import anhalyze.core.anhalyze_plot_utils as apu

        if var:
            self.var = var

        apu.show_var_data_map(self.anha_data,
                              self.lat_range,
                              self.lon_range,
                              depth=self.depth,
                              var=self.var)


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
