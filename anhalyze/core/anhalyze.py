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
    """ Wrapper :py:class: for `xarray.Dataset` with specific implementation
    for netCDF files created with ANHA/NEMO ocean models. Similarly to a
    `xarray.Dataset`, it organizes netCDF data into data_vars, coords,
    attrs and dims.

    Parameters
    ----------
    filename : str
        Filename given with format */*/ANHA?-??????_y????m??d??_grid?.nc
    load_data : bool, optional
        Bool for loading data (Default: True)

    Returns
    -------
    dataset : AnhaDataset
        `AnhaDataset` instance from given filename.

    """

    def __repr__(self):
        """ Return string representation of object
        """
        # TODO return own dims/coords/data_vars, instead of the ones from _xr_dataset
        # return xr.core.formatting.dataset_repr(self._anha_dataset)  # placeholder, may help create own version.
        return f'[Anhalyze] Filename: {self.attrs["filename"]} \n'+str(self._xr_dataset)

    def __init__(self, filename, load_data=True, _xr_dataset=None, _attrs=None):
        """ Initializing object.

        Parameters
        ----------
        _xr_dataset : xarray.Dataset, optional
            Instance of xarray.Dataset object.
        _attrs : dict, optional
            Dict of attributes, use internally.
        """

        # Initialize info from filename
        if _attrs:
            self.attrs = _attrs
        else:
            self._init_filename_attrs(filename)

        # Loading data.
        if _xr_dataset:
            assert type(_xr_dataset) == xr.core.dataset.Dataset, \
                TypeError('[Anhalyze] Parameter xr_dataset incorrect type.')
            self._xr_dataset = _xr_dataset
        else:
            # Open dataset
            if load_data:
                self._xr_dataset = xr.open_dataset(os.path.join(self.attrs['filepath'], self.attrs['filename']))
            else:
                self._xr_dataset = xr.open_dataset(os.path.join(self.attrs['filepath'], self.attrs['filename']), decode_cf=False)

        # Initialize file metadata
        self._load_data = load_data
        self._init_metadata()

        # TODO update verbose with logging level
        # Initialize other attrs
        self.attrs['verbose'] = True

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

        # For some reason there are W dims in a gridT file, removing them since not used.
        if 'T' in self.attrs['grid']:
            del_list = [dim for dim in dims_list if 'w' in dim or 'W' in dim]
            self._xr_dataset = self._xr_dataset.drop_dims(del_list)

        # Init grid dimensions var names
        # Need to exclude 'axis_nbounds'
        self.attrs['dim_x'] = [var for var in dims_list if 'x' in var and 'axis' not in var][0]
        self.attrs['dim_y'] = [var for var in dims_list if 'y' in var][0]
        self.attrs['dim_z'] = [var for var in dims_list if 'depth' in var][0]

        return self._xr_dataset.dims

    def _init_filename_attrs(self, filename):
        """ Initialize properties from filename.
            Filename given in format ANHA?-??????_y????m??d??_grid?.nc
        """

        # Setting up filename
        assert os.path.isfile(filename), f'[Anhalyze] File {filename} not found.'

        if os.path.dirname(filename):
            filepath = os.path.dirname(filename)
            filename = os.path.basename(filename)
        else:
            filepath = ''

        # Init attributes
        self.attrs = {'filename': filename,
                      'filepath': filepath}

        # Making sure filename format is correct
        assert '.nc' in self.attrs['filename'], IOError('[Anhalyze] Incorrect file format.')
        assert '_grid' in self.attrs['filename'], f'[Anhalyze] Filename {filename} does not contain "_grid".'
        assert 'ANHA' in self.attrs['filename'], f'[Anhalyze] Filename {filename} does not contain "ANHA".'

        # Initialize model config
        self.attrs['model_run'] = self.attrs['filename'].split('_')[0]
        assert 'ANHA' in self.attrs['model_run'],\
            f'[Anhalyze] model_run format not recognized: {self.attrs["model_run"]} should include ANHA'
        self.attrs['model_config'] = self.attrs['filename'].split('-')[0]
        self.attrs['model_case'] = self.attrs['filename'].split('-')[1].split('_')[0]

        # Init grid type
        self.attrs['grid'] = self.attrs['filename'].split('_grid')[-1][0]
        grid_value_options = 'TBUVW'
        assert self.attrs['grid'] in grid_value_options,\
            f'[Anhalyze] Grid type not recognized: {self.attrs["grid"]}'
        self.attrs['grid'] = 'grid'+self.attrs['grid']

        # Initialize time
        self.attrs['year'] = self.attrs['filename'].split('y')[-1][:4]
        self.attrs['month'] = self.attrs['filename'].split('m')[-1][:2]
        self.attrs['day'] = self.attrs['filename'].split('d')[1][:2]
        self.attrs['date'] = get_date(self.attrs['filename'])

        # Initialize other attrs
        if 'mask_applied' not in self.attrs.keys():
            self.attrs['mask_applied'] = False

    def _init_metadata(self):
        """ Initialize model properties from filename
        """

        # Initialize attributes
        self.dims = self._init_dims()  # Note this should init before coords and data_vars
        self.data_vars = self._init_data_vars()
        self.coords = self._init_coords()
        if 'xr_attrs' not in self.attrs.keys():
            self._init_xr_attrs()
        self._init_range()

    def _init_range(self):
        """ Initialize boundary values.
        """

        # Get all lat-lon data in file
        lat = self.coords[self.attrs['coord_lat']].data.copy()
        lon = self.coords[self.attrs['coord_lon']].data.copy()

        lat[(lat == 0)] = np.nan
        lon[(lon == 0)] = np.nan

        # Init grid geocoordinates range
        self.attrs['coord_lat_range'] = [np.nanmin(lat), np.nanmax(lat)]
        self.attrs['coord_lon_range'] = [np.nanmin(lon), np.nanmax(lon)]
        self.attrs['coord_depth_range'] = [self.coords[self.attrs['coord_depth']].data.min(),
                                           self.coords[self.attrs['coord_depth']].data.max()]

        # Init grid dims range
        self.attrs['dim_x_range'] = [0, self._xr_dataset.sizes[self.attrs['dim_x']]]
        self.attrs['dim_y_range'] = [0, self._xr_dataset.sizes[self.attrs['dim_y']]]
        self.attrs['dim_z_range'] = [0, self._xr_dataset.sizes[self.attrs['dim_z']]]

    def _update_range(self, coord_name, coord_range, mode='loose'):
        """ Setup data selection range. Assert values are in order, within range and valid.
        """

        # Make sure range has correct format ( a list with two values).
        assert isinstance(coord_range, list), f"[Anhalyze] The variable {coord_name} is not a list."
        assert len(coord_range) == 2, '[Anhalyze] Coordinate range size need to be equal to two.'

        # Make sure the range is ordered from lower to higher.
        if coord_range[0] > coord_range[1]:
            coord_range = [coord_range[1], coord_range[0]]

        # Get full range
        full_range = self.attrs[coord_name+'_range']

        # Make sure values are within full range.
        if 'rigid' in mode:
            # Only accepting values within range
            for value in coord_range:
                assert value >= full_range[0], f"[Anhalyze] {coord_name} value {value} is out of range {full_range}"
                assert value <= full_range[1], f"[Anhalyze] {coord_name} value {value} is out of range {full_range}"

        else:
            # Use edges if values outside range
            if coord_range[0] < full_range[0]:
                mode_message = f'[Anhalyze] Warning: {coord_name} '
                mode_message += f'using edge value {full_range[0]} since given value {coord_range[0]} is out of bounds.'
                coord_range[0] = full_range[0]

                if self.attrs['verbose']:
                    print(mode_message)

            if coord_range[1] > full_range[1]:
                mode_message = f'[Anhalyze] Warning: {coord_name} '
                mode_message += f'using edge value {full_range[1]} since given value {coord_range[1]} is out of bounds.'
                coord_range[1] = full_range[1]

                if self.attrs['verbose']:
                    print(mode_message)

        return coord_range

    def _get_row_col_range(self, lat_range, lon_range):
        """ Get the row AND col range given lat AND lon range.  """

        if not self._load_data:
            raise NotImplementedError("[Anhalyze] Need load_data=True, otherwise not implemented yet.")

        # Get all lat-lon data in file
        lat = self.coords[self.attrs['coord_lat']].data.copy()
        lon = self.coords[self.attrs['coord_lon']].data.copy()

        # Create mask given lat lon values.
        lat_mask = np.ma.filled((lat > lat_range[0]) & (lat < lat_range[1]))
        lon_mask = np.ma.filled((lon > lon_range[0]) & (lon < lon_range[1]))

        # Apply masks to data
        mask = lat
        mask[~(lat_mask & lon_mask)] = np.nan

        # Find the row,col range by collapsing each axis.
        row_ranges = np.where(np.nansum(mask, axis=1) > 0)[0]
        col_ranges = np.where(np.nansum(mask, axis=0) > 0)[0]

        # Select range
        row_range = (row_ranges[0], row_ranges[-1])
        col_range = (col_ranges[0], col_ranges[-1])

        return row_range, col_range

    def _get_row_or_col_range(self, coord_range, coord_name):
        """ Get the row OR col range given lat OR lon range.  """

        if not self._load_data:
            raise NotImplementedError("[Anhalyze] Need load_data=True, otherwise not implemented yet.")

        # Get all lat-lon data in file
        coord = self.coords[coord_name].data.copy()

        # Create mask given lat lon values.
        coord_mask = np.ma.filled((coord > coord_range[0]) & (coord < coord_range[1]))

        # Apply masks to data
        mask = coord
        mask[~coord_mask] = np.nan
        # TODO there maybe a bug at the equator or the meridian.
        #  could try using nans instead of zeros, probably not here but at the end?
        #  I guess not since i'm not using copy() Need to plot to see...

        if 'lat' in coord_name:
            axis = 1
        elif 'lon' in coord_name:
            axis = 0
        else:
            raise ValueError('[Anhalyze] coord_name should be "lat" or "lon".')

        # Find the row,col range by collapsing each axis.
        coord_ranges = np.where(np.nansum(mask, axis=axis) > 0)[0]

        # Select range
        coord_range = (coord_ranges[0], coord_ranges[-1])

        return coord_range

    def sel(self, lat_range=None, lon_range=None, depth_range=None):
        """
        Returns a new `AnhaDataset` with each data array indexed
        along the specified coordinate(s) in `AnhaDataset.coords`.

        In contrast to `AnhaDataset.isel`, indexers for this method should use
        'geographical' values, instead of 'cartesian' integers.

        Parameters
        ----------
        lat_range : list
            Two element list containing min and max Latitude values for selection. [in degrees]
        lon_range : list
            Two element list containing min and max Longitude values for selection. [in degrees]
        depth_range : list
            Two element list containing min and max depth values for selection. [in meters]

        Returns
        -------
        out : AnhaDataset
            An `AnhaDataset` given lat, lon and/or depth range.

        """

        # TODO figure out selecting by location

        # Setting up dict
        dict_range = {}

        # Make copy of xarray
        _xr_dataset = self._xr_dataset.copy()

        # Populating dict for lat, lon selection
        if lat_range and lon_range:
            lat_range = self._update_range('coord_lat', lat_range)
            lon_range = self._update_range('coord_lon', lon_range)

            if self.attrs['verbose']:
                print(f'[Anhalyze] Selecting Latitude range: {lat_range}')
                print(f'[Anhalyze] Selecting Longitude range: {lon_range}')

            # Find row and col ranges from lat or lon values
            row_range, col_range = self._get_row_col_range(lat_range, lon_range)
            dict_range.update({self.attrs['dim_x']: slice(col_range[0], col_range[1]),
                               self.attrs['dim_y']: slice(row_range[0], row_range[1])})

        else:
            if lat_range:
                lat_range = self._update_range('coord_lat', lat_range)

                if self.attrs['verbose']:
                    print(f'[Anhalyze] Selecting Latitude range: {lat_range}')

                # Find row ranges from lat values
                row_range = self._get_row_or_col_range(lat_range, self.attrs['coord_lat'])
                dict_range.update({self.attrs['dim_y']: slice(row_range[0], row_range[1])})
            if lon_range:
                lon_range = self._update_range('coord_lon', lon_range)

                if self.attrs['verbose']:
                    print(f'[Anhalyze] Selecting Longitude range: {lon_range}')

                # Find col ranges from lon values
                col_range = self._get_row_or_col_range(lon_range, self.attrs['coord_lon'])
                dict_range.update({self.attrs['dim_x']: slice(col_range[0], col_range[1])})

        # Lat/lon selection
        if dict_range:
            _xr_dataset = _xr_dataset.isel(dict_range)

        # Populating dict for lat, lon selection
        if depth_range:
            depth_range = self._update_range('coord_depth', depth_range)

            if self.attrs['verbose']:
                print(f'[Anhalyze] Selecting Depth range: {depth_range}')

            dict_range = {self.attrs['dim_z']: slice(depth_range[0], depth_range[1])}

            # Depth selection
            _xr_dataset = _xr_dataset.sel(dict_range)

        # TODO Placeholder for:
        #      I could mask the data outside the lat,lon range region asked here.
        #      Will wait for now, until I have a way to mask data.

        # TODO should add something in the filename just to mark is not the original file.
        return AnhaDataset(os.path.join(self.attrs['filepath'], self.attrs['filename']),
                           load_data=self._load_data, _xr_dataset=_xr_dataset, _attrs=self.attrs)

    def isel(self, x_range=None, y_range=None, z_range=None):
        """
        Returns a new `AnhaDataset` with each data array indexed
        along the specified dimension(s) in `AnhaDataset.dims`.
        Using 'cartesian' integer values.

        Parameters
        ----------
        x_range : list
            Two element list containing min and max x values for selection.
        y_range : list
            Two element list containing min and max y values for selection.
        z_range : list
            Two element list containing min and max z values for selection.

        Returns
        -------
        out : AnhaDataset
            An AnhaDataset given x, y and/or z range.

        """

        # Setting up dict
        dict_range = {}

        # Populating dict for selection
        if x_range:
            x_range = self._update_range('dim_x', x_range)
            if self.attrs['verbose']:
                print(f'[Anhalyze] Selecting x range: {x_range}')
            dict_range.update({self.attrs['dim_x']: slice(x_range[0], x_range[1])})
        if y_range:
            y_range = self._update_range('dim_y', y_range)
            if self.attrs['verbose']:
                print(f'[Anhalyze] Selecting y range: {y_range}')
            dict_range.update({self.attrs['dim_y']: slice(y_range[0], y_range[1])})
        if z_range:
            z_range = self._update_range('dim_z', z_range)
            if self.attrs['verbose']:
                print(f'[Anhalyze] Selecting z range: {z_range}')
            dict_range.update({self.attrs['dim_z']: slice(z_range[0], z_range[1])})

        # Selection of xarray instance
        _xr_dataset = self._xr_dataset.isel(dict_range)

        return AnhaDataset(os.path.join(self.attrs['filepath'], self.attrs['filename']),
                           load_data=self._load_data, _xr_dataset=_xr_dataset, _attrs=self.attrs)

    def show_var_data_map(self, var, idepth=0, proj='AzimuthalConformal', color_range='physical'):
        """ Displays a map for given var in `AnhaDataset.data_vars`.

        Parameters
        ----------
        var : str
            Variable name.
        idepth : int
            Depth value from z_range.

        """
        import anhalyze.core.anhalyze_plot_utils as apu

        assert var in list(self.data_vars), f'[anhalyze] Variable {var} not found in data_vars: {list(self.data_vars)}'

        apu.show_var_data_map(self, var=var, idepth=idepth, proj=proj, color_range=color_range)

    def apply_mask(self, mask_filename=None):
        """ Applies mask to `AnhaDataset.data_vars`.

        Parameters
        ----------
        mask_filename : str
            Mask filename

        """

        if not mask_filename:
            mask_filename = 'ANHA4_mask.nc'

        # Getting mask data
        mask = xr.open_dataset(mask_filename).tmask.data

        # Adding mask data to coords
        self.coords['mask'] = ((self.attrs['dim_z'], self.attrs['dim_y'], self.attrs['dim_x']),
                               mask[0, :, :, :])

        # TODO not exactly what I want. for now this works if used before using sel(), since not returning new instance
        # Applying mask data
        self._xr_dataset = self._xr_dataset.where(self.coords['mask'] == 1)
        self.attrs['mask_applied'] = True


#     def set_location(self):
#         """
#
#         """
#         # TODO  set location with default values


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
