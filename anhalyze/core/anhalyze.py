#!/usr/bin/env python3
# coding: utf-8

# System-related libraries
import os
import numpy as np

# Data-related libraries
import xarray as xr

# Project-related libraries


# Possible other names AnhaModelData, Dataset, AnhaData, AnhaDataframe, AnhaReader, AnhaGrid
class AnhaDataset:
    """ This class opens a single netCDF file with ANHA/NEMO specific properties.

    Attributes
    ----------
        filename : str
            Filename given in format ANHA?-??????_y????m??d??_grid?.nc

    Methods
    -------


    """

    def __repr__(self):
        """ Return string representation of object
        """
        # TODO return own dims/coords/data_vars, instead of the ones from _xr_dataset
        # return xr.core.formatting.dataset_repr(self._anha_dataset)  # placeholder, may help create own version.
        return f'Filename: {self.attrs["filename"]} \n'+str(self._xr_dataset)

    def __init__(self, filename, load_data=True, xr_dataset=None):
        """ Initializing object.

        Parameters
        ----------
        filename : str
            Filename given in format ANHA?-??????_y????m??d??_grid?.nc
            or *_mask*.nc
        load_data : bool, optional
            Bool for loading data (Default is False)
        # cartesian : bool, optional
        #     Bool for using cartesian coordinates (Default is True)
        #     Otherwise use geographical coordinates and SI units

        """

        # Initialize info from filename
        self._init_filename_attrs(filename)

        # Loading data.
        if xr_dataset:
            #TODO assert xr_data is of xarray.Dataset type.
            self._xr_dataset = xr_dataset
        else:
            # Open dataset
            if load_data:
                self._xr_dataset = xr.open_dataset(os.path.join(self.attrs['filepath'], self.attrs['filename']))
            else:
                self._xr_dataset = xr.open_dataset(os.path.join(self.attrs['filepath'], self.attrs['filename']), decode_cf=False)

        # Initialize file metadata
        self._load_data = load_data
        self._init_metadata()

    def _init_coords(self):
        """ Initialize coordinates
        """

        # Load all keys from data into lists
        dims_list = list(self._xr_dataset.dims)
        vars_list = list(self._xr_dataset.data_vars)
        coords_list = list(self._xr_dataset.coords)

        # Select appropriate list
        if self._load_data:
            geocoords_list = coords_list
        else:
            geocoords_list = vars_list

        # Init grid geocoordinates var names
        self.attrs['coord_lat'] = [var for var in geocoords_list if 'nav_lat' in var][0]
        self.attrs['coord_lon'] = [var for var in geocoords_list if 'nav_lon' in var][0]
        self.attrs['coord_depth'] = [var for var in dims_list if 'depth' in var][0]

        # TODO this may still be an issue if we want coords to reflect the  "real" ones (when load_data=T)
        return self._xr_dataset.coords

    def _init_data_vars(self):
        """ Initialize data variables
        """
        return self._xr_dataset.data_vars

    def _init_xr_attrs(self):
        """ Initialize data attributes
        """
        self.attrs['xr_attrs'] = self._xr_dataset.attrs

    def _init_dims(self):
        """ Initialize dimensions
        """

        dims_list = list(self._xr_dataset.dims)

        # Init grid dimensions var names
        # Need to exclude 'axis_nbounds'
        self.attrs['dim_x'] = [var for var in dims_list if 'x' in var and 'axis' not in var][0]
        self.attrs['dim_y'] = [var for var in dims_list if 'y' in var][0]

        return self._xr_dataset.dims

    def _init_filename_attrs(self, filename):
        """ Initialize properties from filename
        """

        # Setting up filename
        assert os.path.isfile(filename)

        if os.path.dirname(filename):
            filepath = os.path.dirname(filename)
            filename = os.path.basename(filename)
        else:
            filepath = ''

        # Init attributes
        self.attrs = {'filename': filename,
                      'filepath': filepath}

        if '_grid' in self.attrs['filename']:

            #TODO add a few asserts here in filename format

            # Initialize model config
            self.attrs['model_run'] = self.attrs['filename'].split('_')[0]
            assert 'ANHA' in self.attrs['model_run'], f'Model run format not recognized: {self.attrs["model_run"]}'
            self.attrs['model_config'] = self.attrs['filename'].split('-')[0]
            self.attrs['model_case'] = self.attrs['filename'].split('-')[1].split('_')[0]

            # Init grid type
            self.attrs['grid'] = self.attrs['filename'].split('_grid')[-1][0]
            grid_value_options = 'TBUVW'
            assert self.attrs['grid'] in grid_value_options, f'Grid type not recognized: {self.attrs["grid"]}'
            self.attrs['grid'] = 'grid'+self.attrs['grid']
            self.attrs['is_mask'] = False

            # Initialize time
            self.attrs['year'] = self.attrs['filename'].split('y')[-1][:4]
            self.attrs['month'] = self.attrs['filename'].split('m')[-1][:2]
            self.attrs['day'] = self.attrs['filename'].split('d')[1][:2]
            self.attrs['date'] = get_date(self.attrs['filename'])

        elif '_mask' in self.attrs['filename']:
            raise NotImplementedError("Loading mask data not implemented yet.")

            # # Init grid type.
            # self.attrs['grid'] = 'mask'
            # self.attrs['is_mask'] = True

        else:
            raise NotImplementedError("Loading non_grid data not implemented yet.")
            # TODO ask Luiz again about data types.
            # self.attrs['grid'] = ''
            # self.attrs['is_mask'] = False

    def _init_metadata(self):
        """ Initialize model properties from filename
        """

        # Initialize xarray main attributes
        self.data_vars = self._init_data_vars()
        self.coords = self._init_coords()
        self.dims = self._init_dims()
        self._init_xr_attrs()

        if '_grid' in self.attrs['filename']:
            pass
            # # Init grid dimensions var names
            # # Need to exclude 'axis_nbounds'
            # self.x_var_name = [var for var in dims_list if 'x' in var and 'axis' not in var][0]
            # self.y_var_name = [var for var in dims_list if 'y' in var][0]

        elif '_mask' in self.attrs['filename']:
            raise NotImplementedError("Loading mask data not implemented yet.")

            # # Init grid dimensions var names
            # self.x_var_name = [var for var in dims_list if 'x' in var][0]
            # self.y_var_name = [var for var in dims_list if 'y' in var][0]
            #
            # # Init grid geocoordinates var names
            # self.lat_var_name = [var for var in coords_list if 'nav_lat' in var][0]
            # self.lon_var_name = [var for var in coords_list if 'nav_lon' in var][0]
            # self.depth_var_name = [var for var in coords_list if 'gdep' in var][0]

        else:
            pass

    def _setup_selection_range(self, lat_range=None, lon_range=None, i_range=None, j_range=None, init=False):
        """Setup data selection range,


        """
        # TODO could, set an absolut minimum for lat at: -20.07611 , or min in file.

        if init:

            # For cartesian coordinates
            self.i_begin = 0
            self.i_end = self.dims[self.x_var_name].size
            self.j_begin = 0
            self.j_end = self.dims[self.y_var_name].size
            self.k_begin = 0
#            self.k_end = self.dims[self.depth_var_name].size

            self.i_range = [self.i_begin, self.i_end]
            self.j_range = [self.j_begin, self.j_end]
#            self.k_range = [self.k_begin, self.k_end]

            # For geographical coordinates
            self.lat_range = [self._xr_dataset[self.lat_var_name][:].min(),
                              self._xr_dataset[self.lat_var_name][:].max()]
            self.lon_range = [self._xr_dataset[self.lon_var_name][:].min(),
                              self._xr_dataset[self.lon_var_name][:].max()]
            self.depth_range = [self._xr_dataset[self.depth_var_name][:].min(),
                                self._xr_dataset[self.depth_var_name][:].max()]

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

    def _get_row_col_range(self, lat_range, lon_range):
        """ Get the row and col range given lat and lon range.  """

        # TODO add asserts to handle lat and lon ranges

        if not self._load_data:
            raise NotImplementedError("Need load_data=True, otherwise not implemented yet.")

        # Get all lat-lon data in file
        lat = self.coords[self.attrs['coord_lat']].data[:]
        lon = self.coords[self.attrs['coord_lon']].data[:]

        # Create mask given lat lon values.
        lat_mask = np.ma.filled((lat > lat_range[0]) & (lat < lat_range[1]))
        lon_mask = np.ma.filled((lon > lon_range[0]) & (lon < lon_range[1]))

        # Apply masks to data
        mask = lat
        mask[~(lat_mask & lon_mask)] = 0

        # Find the row,col range by collapsing each axis.
        row_ranges = np.where(mask.sum(axis=1) > 0)[0]
        col_ranges = np.where(mask.sum(axis=0) > 0)[0]

        # Select range
        row_range = (row_ranges[0], row_ranges[-1])
        col_range = (col_ranges[0], col_ranges[-1])

        return row_range, col_range

    def sel(self, lat_range=None, lon_range=None):
        """
        """

        # Find row,col ranges from lat,lon values
        row_range, col_range = self._get_row_col_range(lat_range, lon_range)

        # Selection of xarray instance
        _xr_dataset = self._xr_dataset.isel({self.attrs['dim_y']: slice(row_range[0], row_range[1]),
                                            self.attrs['dim_x']: slice(col_range[0], col_range[1])})

        return AnhaDataset(self.attrs['filename'], load_data=self._load_data, xr_dataset=_xr_dataset)

    def isel(self, lat_range=None, lon_range=None):
        """
        """

        # Set row,col ranges from lat,lon values
        row_range, col_range = lat_range, lon_range

        return self._xr_dataset.isel({self.attrs['dim_y']: slice(row_range[0], row_range[1]),
                                      self.attrs['dim_x']: slice(col_range[0], col_range[1])})

    def show_var_data_map(self, var=''):
        """ Displays map of given var.
        """
        import anhalyze.core.anhalyze_plot_utils as apu

        if var:
            self.var = var

        apu.show_var_data_map(self._xr_dataset,
                              self.lat_range,
                              self.lon_range,
                              depth=self.depth,
                              var=self.var)

#     def set_range(self, lat_range=None, lon_range=None, depth_range=None,
#                   i_range=None, j_range=None, k_range=None):
#         """ Set data range, and update both geographical and cartesian coordinates.
#             #TODO clarify the difference between this and _setup_selection_range
#         """
#         # TODO  set location with default values
#         # TODO set lat, long, depth ranges
#         # TODO set i, j, k ranges
#
#         if lat_range:
#                 self.lat_range = lat_range
# #                self.i_range =   # TODO move get_row_col_range to this class, and use it here.
#
#         if lon_range:
#                 self.lon_range = lon_range
#
#         if depth_range:
#                 self.depth_range = depth_range
#
#         if i_range:
#             self.i_range = i_range
#
#         if j_range:
#             self.j_range = j_range
#
#         if k_range:
#             self.k_range = k_range
#
#     def set_location(self):
#         """
#
#         """


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
