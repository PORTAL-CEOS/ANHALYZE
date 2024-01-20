#!/usr/bin/env python3
# coding: utf-8

import os

# Data-related libraries
import netCDF4 as nc

# Possible other names AnhaModelData, Dataset, AnhaData, AnhaDataframe
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

        # Initialize file type from filename
        if '_grid' in filename:
            self.grid = filename.split('_grid')[-1][0]
            self.is_mask = False
        elif '_mask' in filename:
            self.grid = 'K'  # used for masks
            self.is_mask = True
        else:
            self.grid = ''  # used for masks
            self.is_mask = None

        grid_values = 'TBUVWK'
        assert self.grid in grid_values, 'File type not recognized'

        # Initialize data properties
        if load_data:
            self.anha_data = nc.Dataset(filename)
            self.depth = self.anha_data.dimensions['deptht'].size
            self.i_range = [0, self.anha_data.dimensions['x'].size]
            self.j_range = [0, self.anha_data.dimensions['y'].size]
            self.xdim = self.anha_data.dimensions['x'].size
            self.ydim = self.anha_data.dimensions['x'].size
        else:
            self.anha_data = None

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