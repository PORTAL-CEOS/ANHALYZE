#!/usr/bin/env python3
# coding: utf-8

import os

# Data-related libraries
import numpy as np
import netCDF4 as nc


# Possible other names AnhaModelData, Dataset, AnhaData, AnhaDataframe
class AnhaDataset:
    """ This class will do analysis of ANHA4 data, for now it initializes the location of files.
   ...

    Attributes
    ----------
        filename : str
            Filename given in format ANHA?-??????_y????m??d??_grid?.nc

    Methods
    -------

    """

    def __init__(self, filename, load_data=True):
        """ Initializing class.

        Parameters
        ----------
        filename : str
            Filename given in format ANHA?-??????_y????m??d??_grid?.nc
            or *_mask*.nc
        load_data : bool, optional
            Bool for loading data (Default is False)

        """

        assert os.path.isfile(filename)
        self.filename = filename

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

        if load_data:
            self.anha_data = nc.Dataset(filename)
        else:
            self.anha_data = None


