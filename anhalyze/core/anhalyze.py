#!/usr/bin/env python3
# coding: utf-8

# System-related libraries
import os
import numpy as np

# Data-related libraries
import xarray as xr

# Project-related libraries
import anhalyze
import anhalyze.config as config


#
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
    mask_filename : str, optional
        Mask filename(Default: )

    Returns
    -------
    dataset : AnhaDataset
        `AnhaDataset` instance from given filename.

    """

    def __repr__(self):
        """ Return string representation of object
        """

        anhalyze_repr = f'[Anhalyze] AnhaDataset. \n'
        anhalyze_repr += f'[Anhalyze] Filename: {self.attrs["filename"]}\n'
        anhalyze_repr += f'[Anhalyze] Description: {self.attrs["description"]}\n'
        anhalyze_repr += f'[Anhalyze] File category: {self.attrs["file_category"]}\n'

        xarray_repr = str(self._xr_dataset)

        anhalyze_warning = '\n[Anhalyze] Note: Above we show the xarray repr of this file, ' \
                           'it should be mostly complete, ' \
                           'but use `self.attrs` for the full set of Attributes.'

        return "{0}{1}{2}".format(anhalyze_repr, xarray_repr, anhalyze_warning)

    def _repr_html_(self):
        """ Returns html representation of object
        """

        anhalyze_repr = f'[Anhalyze] AnhaDataset. <br/>'
        anhalyze_repr += f'[Anhalyze] Filename: {self.attrs["filename"]}<br/>'
        anhalyze_repr += f'[Anhalyze] Description: {self.attrs["description"]}<br/>'
        anhalyze_repr += f'[Anhalyze] File category: {self.attrs["file_category"]}<br/>'

        xarray_repr = self._xr_dataset._repr_html_()

        anhalyze_warning = '[Anhalyze] Note: Above we show the xarray repr of this file, ' \
                           'it should be mostly complete, ' \
                           'but use `self.attrs` for the full set of Attributes.'

        return "{0}{1}{2}".format(anhalyze_repr, xarray_repr, anhalyze_warning)

    def __init__(self, filename, load_data=True, mask_filename=None, _xr_dataset=None, _attrs=None):
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
            assert _xr_dataset, '[Anhalyze] Parameter _xr_dataset needs to be provided with _attrs'
            self.attrs = _attrs
        else:
            self._init_filename_attrs(filename)
            self.attrs['file_category'] = 'original'

        # Loading data
        if _xr_dataset:
            # Check correct xarray.Dataset type.
            assert type(_xr_dataset) == xr.core.dataset.Dataset, \
                TypeError('[Anhalyze] Parameter _xr_dataset incorrect type.')
            # Updating attrs
            _xr_dataset.attrs = _attrs
            # Setting internal xarray.Dataset
            self._xr_dataset = _xr_dataset
        else:
            # Open dataset
            if load_data:
                self._xr_dataset = xr.open_dataset(os.path.join(self.attrs['filepath'], self.attrs['filename']))
            else:
                raise FutureWarning("[Anhalyze] load_data=false option hasn't been fully developed.")
                # self._xr_dataset = xr.open_dataset(os.path.join(self.attrs['filepath'], self.attrs['filename']),
                #                                    decode_cf=False)

        # Loading mask data
        if not _xr_dataset:
            # Get mask from filename
            if 'mask' not in list(self._xr_dataset.data_vars):
                self._mask_filename = mask_filename

        # Initialize file metadata
        self._load_data = load_data
        self._init_metadata()

        # Initialize other attrs
        # TODO could replace verbose with logging levels
        self._verbose = True

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
        if any('depth' in dim for dim in dims_list):
            self.attrs['coord_depth'] = [var for var in dims_list if 'depth' in var][0]

        # TODO this may still be an issue if we want coords to reflect the "real" ones (when load_data=T)
        return self._xr_dataset.coords

    def _init_data_vars(self):
        """ Initialize data variables
        """
        # TODO need to add case of load_data = False

        # TODO May need a way to force to add mask even if already in file. In case one wants to update mask.
        #      Could create a child class that 'fixes' files.

        # Getting mask data if not already in data_vars
        if 'mask' not in list(self._xr_dataset.data_vars):
            self._get_mask(mask_filename=self._mask_filename)

        return self._xr_dataset.data_vars

    def _init_xr_attrs(self):
        """ Initialize data attributes
        """

        self.attrs |= self._xr_dataset.attrs

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
        if any('depth' in dim for dim in dims_list):
            self.attrs['dim_z'] = [var for var in dims_list if 'depth' in var][0]

        return self._xr_dataset.dims

    def _init_filename_attrs(self, filename):
        """ Initialize properties from filename.
            Filename given in format ANHA?-??????_y????m??d??_grid?.nc
        """

        # Setting up filename
        assert os.path.isfile(filename), f'[Anhalyze] File {filename} not found.'

        if os.path.dirname(filename):
            filepath = os.path.dirname(os.path.realpath(filename))
            filename = os.path.basename(filename)
        else:
            filepath = ''

        # Init attributes
        self.attrs = {'filename': filename,
                      'filepath': filepath}

        # Making sure filename format is correct
        assert '.nc' in self.attrs['filename'], IOError('[Anhalyze] Incorrect file format.')
        assert '_grid' or '_ice' in self.attrs['filename'], \
            f'[Anhalyze] Filename {filename} does not contain neither "_grid" nor "_ice".'
        assert 'ANHA' in self.attrs['filename'], f'[Anhalyze] Filename {filename} does not contain "ANHA".'

        # Initialize model config
        self.attrs['model_run'] = self.attrs['filename'].split('_')[0]
        assert 'ANHA' in self.attrs['model_run'], \
            f'[Anhalyze] model_run format not recognized: {self.attrs["model_run"]} should include ANHA'
        self.attrs['model_config'] = self.attrs['filename'].split('-')[0]
        self.attrs['model_case'] = self.attrs['filename'].split('-')[1].split('_')[0]

        # Init grid type
        grid_value_options = ['gridT', 'gridB',
                              'gridU', 'gridV', 'gridW',
                              'icemod', 'icebergs']
        self.attrs['grid'] = [grid for grid in grid_value_options if grid in self.attrs['filename']][0]
        assert self.attrs['grid'] in grid_value_options,\
            f'[Anhalyze] Grid type not recognized: {self.attrs["grid"]}'

        # Initialize time
        self.attrs['date'] = get_date(self.attrs['filename'])
        self.attrs['year'] = self.attrs['date'][1:5]
        self.attrs['month'] = self.attrs['date'][6:8]
        self.attrs['day'] = self.attrs['date'][9:11]

        # Initialize other attrs

    def _init_metadata(self):
        """ Initialize model properties from filename
        """

        # Initialize attributes
        self.dims = self._init_dims()  # Note this should init before coords and data_vars
        self.data_vars = self._init_data_vars()
        self.coords = self._init_coords()
        if 'description' not in self.attrs.keys():
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
        
        # Init grid dims range
        self.attrs['dim_x_range'] = [0, self._xr_dataset.sizes[self.attrs['dim_x']]]
        self.attrs['dim_y_range'] = [0, self._xr_dataset.sizes[self.attrs['dim_y']]]
        
        # Only for 3 dimension variable
        if 'coord_depth' in self.attrs.keys():
            self.attrs['coord_depth_range'] = [self.coords[self.attrs['coord_depth']].data.min(),
                                               self.coords[self.attrs['coord_depth']].data.max()]
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
        full_range = self.attrs[coord_name + '_range']

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
                mode_message += f'using edge value {full_range[0]:.2f} ' \
                                f'since given value {coord_range[0]} is out of bounds.'
                coord_range[0] = full_range[0]

                if self._verbose:
                    print(mode_message)

            if coord_range[1] > full_range[1]:
                mode_message = f'[Anhalyze] Warning: {coord_name} '
                mode_message += f'using edge value {full_range[1]:.2f} ' \
                                f'since given value {coord_range[1]} is out of bounds.'
                coord_range[1] = full_range[1]

                if self._verbose:
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

    def _get_var_data_array(self, var='', mask_data=True):
        """ Returns DataArray for given var
        """

        # Get DataArray for given var
        var_da = self.data_vars[var]

        # Mask data
        if mask_data:
            var_da = self._apply_mask(var_da)

        return var_da

    def _get_mask(self, mask_filename=None):
        """ Get mask from given mask_filename or default location.

        Parameters
        ----------
        mask_filename : str
            Mask filename

        """

        if not mask_filename:
            # Check if a mask is declared as environment variable
            if os.environ.get('MASK_PATH_FILENAME'):
                mask_filename = os.environ.get('MASK_PATH_FILENAME')
            else:
                # Get mask from standard location.
                mask_filename = os.path.join(anhalyze.PACKAGE_DATA_DIR, 'ANHA4_mask.nc')

                # Download mask if not present and if config allows.
                if not os.path.isfile(mask_filename) and config.package_data['mask']['autodownload_file']:
                    print('[Anhalyze] No mask file found, downloading ...')

                    from anhalyze.core.downloader import download_mask
                    download_mask()

        # Check if there is a mask file
        assert_message = '[Anhalyze] No mask file found, '
        assert_message += 'please provide correct path, or see README for options.'
        assert os.path.isfile(mask_filename), assert_message

        # Get mask
        if mask_filename:

            # Getting mask data
            if 'gridT' in self.attrs['grid']:
                mask = xr.open_dataset(mask_filename).tmask.data
            elif 'gridU' in self.attrs['grid']:
                mask = xr.open_dataset(mask_filename).umask.data
            elif 'gridV' in self.attrs['grid']:
                mask = xr.open_dataset(mask_filename).vmask.data
            elif 'gridW' in self.attrs['grid']:
                # TODO use tmask and add a warning
                raise NotImplementedError('TODO: Need to figure this out')
                # mask = xr.open_dataset(mask_filename).tmask.data
            else:
                mask = xr.open_dataset(mask_filename).tmask.data

            # TODO: for icemod,  there are u and v data variables that need to have their exceptions
            #       (with in the same file)
            #       rest gridT.  For  icebergs is all gridT.
            # TODO: need to get an icemod, and an iceberg file for testing.

            # TODO: add assertion of dimensions. mask dimensions need to match

            # TODO may want to update this to save mask dataArray instead of numpy array (.data)
            # Adding mask data to data_vars
            
            if 'dim_z' not in self.attrs.keys():
                self._xr_dataset = self._xr_dataset.assign({'mask': ((self.attrs['dim_y'],
                                                                      self.attrs['dim_x']),
                                                                     mask[0, 0, :, :])})
                
            else:
                self._xr_dataset = self._xr_dataset.assign({'mask': ((self.attrs['dim_z'],
                                                                      self.attrs['dim_y'],
                                                                      self.attrs['dim_x']),
                                                                     mask[0, :, :, :])})

            # Add mask filename to attrs
            self._xr_dataset.attrs['mask_filename'] = mask_filename

        else:
            raise OSError('[Anhalyze] No mask/mesh file found.')

    def _apply_mask(self, var_data, at_top_layer=False):
        """ Applies mask to single var in `AnhaDataset.data_vars`.
        """

        # Applying mask data
        if at_top_layer:
            var_data[~np.ma.filled((1 == self.data_vars['mask'][0, :]))] = np.nan
        else:
            # var_data[~np.ma.filled((1 == self.data_vars['mask']))] = np.nan
            var_data = var_data.where(self.data_vars['mask'] == 1)

        # previous versions.
        # self._xr_dataset = self._xr_dataset.where(self.coords['mask'] == 1)
        # self.data_vars = self._xr_dataset.where(self.coords['mask'] == 1)

        return var_data

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

            if self._verbose:
                print(f'[Anhalyze] Selecting Latitude range: {lat_range}')
                print(f'[Anhalyze] Selecting Longitude range: {lon_range}')

            # Find row and col ranges from lat or lon values
            row_range, col_range = self._get_row_col_range(lat_range, lon_range)
            dict_range.update({self.attrs['dim_x']: slice(col_range[0], col_range[1]),
                               self.attrs['dim_y']: slice(row_range[0], row_range[1])})

        else:
            if lat_range:
                lat_range = self._update_range('coord_lat', lat_range)

                if self._verbose:
                    print(f'[Anhalyze] Selecting Latitude range: {lat_range}')

                # Find row ranges from lat values
                row_range = self._get_row_or_col_range(lat_range, self.attrs['coord_lat'])
                dict_range.update({self.attrs['dim_y']: slice(row_range[0], row_range[1])})
            if lon_range:
                lon_range = self._update_range('coord_lon', lon_range)

                if self._verbose:
                    print(f'[Anhalyze] Selecting Longitude range: {lon_range}')

                # Find col ranges from lon values
                col_range = self._get_row_or_col_range(lon_range, self.attrs['coord_lon'])
                dict_range.update({self.attrs['dim_x']: slice(col_range[0], col_range[1])})

        # Lat/lon selection
        if dict_range:
            _xr_dataset = _xr_dataset.isel(dict_range)

        # Populating dict for depth selection
        if depth_range:
            depth_range = self._update_range('coord_depth', depth_range)

            if self._verbose:
                print(f'[Anhalyze] Selecting Depth range: {depth_range}')

            dict_range = {self.attrs['dim_z']: slice(depth_range[0], depth_range[1])}

            # Depth selection
            _xr_dataset = _xr_dataset.sel(dict_range)

        # Set attrs
        _attrs = self.attrs.copy()
        _attrs['file_category'] = 'regional'
        # TODO could add section/transect or something specific like this.

        return AnhaDataset('', load_data=self._load_data, _xr_dataset=_xr_dataset, _attrs=_attrs)

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
            if self._verbose:
                print(f'[Anhalyze] Selecting x range: {x_range}')
            dict_range.update({self.attrs['dim_x']: slice(x_range[0], x_range[1])})
        if y_range:
            y_range = self._update_range('dim_y', y_range)
            if self._verbose:
                print(f'[Anhalyze] Selecting y range: {y_range}')
            dict_range.update({self.attrs['dim_y']: slice(y_range[0], y_range[1])})
        if z_range:
            z_range = self._update_range('dim_z', z_range)
            if self._verbose:
                print(f'[Anhalyze] Selecting z range: {z_range}')
            dict_range.update({self.attrs['dim_z']: slice(z_range[0], z_range[1])})

        # Selection of xarray instance
        _xr_dataset = self._xr_dataset.isel(dict_range)

        # Set attrs
        _attrs = self.attrs.copy()
        _attrs['file_category'] = 'regional'
        # TODO could add section/transect or something specific like this.

        return AnhaDataset('', load_data=self._load_data, _xr_dataset=_xr_dataset, _attrs=_attrs)

    def show_var_data_map(self, var, color_range='default', savefig=None, projection_name='LambertConformal'):
        """ Displays a map for given var in `AnhaDataset.data_vars`.

        Parameters
        ----------
        var : str
            Variable name.
        color_range : str
            Color range either `default` limits, `local` data values or a two items list [vmin, vmax].
            Color range options:
             default: Limits decide by Anhalyze developers. It is based on the more
                      likely limits the user can find in a ANHA4 outputs.
             local: Color range based on the values within the area selected by the user.
             [vmin, vmax] = List of color range limits chosen by the user.
        savefig : str
            Filename to save figure including path.
        projection_name : str
            Projection name from Cartopy list. The projections available are: 'PlateCarree',
            'LambertAzimuthalEqualArea', 'AlbersEqualArea', 'NorthPolarStereo', 'Orthographic', 'Robinson',
            'LambertConformal', 'Mercator', and 'AzimuthalEquidistant'.
        """

        import anhalyze.core.anhalyze_plot_utils as apu

        assert var in list(self.data_vars), f'[anhalyze] Variable {var} not found in data_vars: {list(self.data_vars)}'

        # Get DataArray for given var
        var_da = self._get_var_data_array(var=var)

        # Select top layer
        if 'dim_z' in self.attrs.keys():
            var_da = var_da.isel(indexers={self.attrs['dim_z']: [0]})

        # Show var data map
        apu.show_var_data_map(var_da,
                              attrs=self.attrs,
                              color_range=color_range,
                              savefig=savefig,
                              proj_name=projection_name)

    def to_netcdf(self, path=None, filename=None, suffix='_CutRegion', **kwargs):
        """ Writes `AnhaDataset` contents to netCDF file.
            For additional options see: `self._xr_dataset.to_netcdf`

            Note: Behaviour differs from `xarray.Dataset.to_netcdf` since here we avoid overwriting files.

        Parameters
        ----------
        path : str, optional
            Path to which to save this  `AnhaDatabase`.
        filename : str, optional
            Filename to which to save this `AnhaDatabase`.
        suffix : str, default: _CutRegion.nc
            Suffix added to filename to avoid overwriting.

        """

        # set path
        if not path:
            path = self.attrs['filepath']
            # TODO could change default to current location

        # set filename
        if filename:
            assert '.nc' in filename, ValueError('[Anhalyze] Filename should be .nc type.')
        else:
            filename = self.attrs['filename']

        # Making sure the filename ends in .nc
        if '.nc' not in suffix[-3:]:
            suffix += '.nc'

        # Setting up new full filename
        new_full_filename = os.path.join(path, filename.replace('.nc', suffix))

        # Avoiding overwriting files by adding extra suffix continuously until available.
        while os.path.isfile(new_full_filename):
            print(f'[Anhalyze] Warning, file exists: {new_full_filename}')
            new_full_filename = new_full_filename.replace('.nc', '_copy.nc')

        if self._verbose:
            print(f'[Anhalyze] Saving: {new_full_filename}')

        # Updating filename
        self._xr_dataset.attrs['filename'] = os.path.basename(new_full_filename)

        # Saving new file
        self._xr_dataset.to_netcdf(new_full_filename, **kwargs)


def get_date(filename, how=None):
    # TODO find/apply naming convention documentation.
    """  Get date information from filename.
         Assuming filename format: */*/ANHA?-??????_y????m??d??_{grid}_*.nc

        Parameters
        ----------
        filename: str
            Filename
        how : str
            What to output. Multiple options and output formats possible:
                None/other: str: 'y????m??d??'
                ymd: tuple (y,m,d)
                y: int: ????
                m: int: ??
                d: int: ??

    """

    # Get full date from filename
    date = filename.split('_')[1]
    assert date[0] == 'y' and date[5] == 'm' and date[8] == 'd', '[Anhalyze] Filename format not supported.'

    # Return format specific date info
    if how == 'ymd':
        return int(date[1:5]), int(date[6:8]), int(date[9:11])
    elif how == 'y':
        return int(date[1:5])
    elif how == 'm':
        return int(date[6:8])
    elif how == 'd':
        return int(date[9:11])
    else:
        return date
