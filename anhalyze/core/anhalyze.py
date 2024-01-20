#!/usr/bin/env python3
# coding: utf-8

# System-related libraries
import os

# Data-related libraries
import netCDF4 as nc

# Project-related libraries
import anhalyze.core.anhalyze_plot_utils as apu


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
            Otherwise use geographic coordinates and SI units

        """

        assert os.path.isfile(filename)
        self.filename = filename

        # Initialize data properties
        # TODO update to use xarray to be able to load file info without loading data
        if load_data:
            self.anha_data = nc.Dataset(filename)
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
        if '_grid' in filename:
            self.model_run = filename.split('_')[0]
            self.year = filename.split('y')[-1][:4]
            self.month = filename.split('m')[-1][:2]
            self.day = filename.split('d')[1][:2]

            assert 'ANHA' in self.model_run, 'Filename format not recognized.'

            self.configuration = self.filename.split('-')[0]

        else:
            pass


        # Initialize unit properties
        self.cartesian = cartesian

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
        # TODO set an absolut minimum for lat at: -20.07611 , or min in file.

        if init:
            self.i_begin = 0
            self.i_end = self.anha_data.dimensions['x'].size
            self.j_begin = 0
            self.j_end = self.anha_data.dimensions['y'].size
            self.k_begin = 0
            self.k_end = self.anha_data.dimensions['deptht'].size

            self.lat_range = [self.anha_data[self.lat_var_name][:].min(),
                              self.anha_data[self.lat_var_name][:].max()]
            self.lon_range = [self.anha_data[self.lon_var_name][:].min(),
                              self.anha_data[self.lon_var_name][:].max()]
            self.i_range = [0, self.anha_data.dimensions['x'].size]
            self.j_range = [0, self.anha_data.dimensions['y'].size]
        else:
            if lat_range:
                self.lat_range = lat_range
            if lon_range:
                self.lon_range = lon_range
            #TODO add assest if ranges are correct
            self.i_range = [0, self.anha_data.dimensions['x'].size]
            self.j_range = [0, self.anha_data.dimensions['y'].size]


    def show_var_data_map(self, var=''):
        """ Displays map of given var.
        """

        if var:
            self.var = var

        apu.show_var_data_map(self.anha_data,
                              self.lat_range,
                              self.lon_range,
                              depth=self.depth,
                              var=self.var)
