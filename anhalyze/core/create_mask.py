#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Masking region from ANHA4 domain using Polygon Selector

Interactively allows the user to create a mask file drawing a polygon over the ANHA4 domain.
The user can zoom in the figure to get a more accurate polygon.

"""
# %%
# Importing libraries
import os
import matplotlib
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Polygon related libraries
from PIL import Image, ImageDraw
from matplotlib.widgets import PolygonSelector

# Project-related libraries
import anhalyze.config as config
import anhalyze

matplotlib.use('TkAgg')

# %%


def create_mask(mask_source=None, grid='tmask', path=None, suffix='_CutMask.nc', verbose=True):
    """ Get mask from given mask_source or default location.

    Parameters
    ----------
    mask_source : str
        Original ANHA4 mask file.

    grid : str
        Grid in which the user wants to create the new mask from.
        Options: 'tmask', 'umask', and 'vmask',

    path : str
        Full path where anhalyze will save the new mask file as NetCDF.
        The default path is anhalyze.config directory.

    suffix : str
        Indicate a suffix to be added to the new mask file name

    verbose : bool
        Activate the verbose mode.
    """
    # Get mask file source. If the user doesn't insert mask source file name,
    # the code will check if there are any mask file already in the config folder.
    # If it doesn't find anything, will use anhalyze 'downloader' to download a default mask file.
    if not mask_source:
        # Check if a mask is declared as environment variable
        if os.environ.get('MASK_PATH_FILENAME'):
            mask_source = os.environ.get('MASK_PATH_FILENAME')
        else:
            # Get mask from standard location.
            mask_source = os.path.join(anhalyze.PACKAGE_DATA_DIR, 'ANHA4_mask.nc')

            # Download mask if not present and if config allows.
            if not os.path.isfile(mask_source) and config.package_data['mask']['autodownload_file']:
                print('[Anhalyze] No mask file found, downloading ...')

                from anhalyze.core.downloader import download_mask
                download_mask()

    # Check wheter is a netcdf file
    assert '.nc' in mask_source[-3:], (f'[Anhalyze] Mask source file {mask_source.split("/")[-1]}'
                                       f' must be a NetCDF file type.')

    # Check if there is a mask file
    assert_message = '[Anhalyze] No mask file found, '
    assert_message += 'please provide correct path, or see README for options.'
    assert os.path.isfile(mask_source), assert_message

    # Get mask
    if mask_source:

        # Getting mask data depending on the grid
        if grid == 'tmask':
            mask = xr.open_dataset(mask_source).tmask.data[0, 0, :, :]
            mask_var = 'tmask'
        elif grid == 'umask':
            mask = xr.open_dataset(mask_source).umask.data[0, 0, :, :]
            mask_var = 'umask'
        elif grid == 'vmask':
            mask = xr.open_dataset(mask_source).vmask.data[0, 0, :, :]
            mask_var = 'vmask'
        elif grid == 'wmask':
            # TODO use tmask and add a warning
            raise NotImplementedError('TODO: Need to figure this out')
            # mask = xr.open_dataset(mask_filename).tmask.data
    else:
        raise OSError('[Anhalyze] No mask/mesh file found.')

    # Extract latitude and longitude information
    lat = xr.open_dataset(mask_source).nav_lat.data
    lon = xr.open_dataset(mask_source).nav_lon.data

    # %%
    # Create polygon
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.pcolormesh(mask)

    def onselect(verts):
        pass

    # PolygonSelector open a window using the previous figure.
    # The user can draw the vertices around the region of interes.
    selector = PolygonSelector(ax, onselect)

    print("Click on the figure to create a polygon. Select as many vertices as you need until close the polyogon.")
    print("Press 'esc' ley to start a new polygon.")
    print("Close the window to save the new mask as NetCDF.")
    plt.show(block=True)

    # %%
    # Set everything out of the polygon area as zero.

    vert = np.array([np.intc(np.round(vertices, 0)) for vertices in selector.verts], 'int32')
    vert = [tuple(x) for x in vert]

    # Create a new numpy array with the ANHA4 mask dimensions
    mask_sel = np.zeros(mask.shape)

    # Mask dimensions
    x = mask_sel.shape[1]
    y = mask_sel.shape[0]

    # Image and ImageDraw.Draw().Polygon to fill up the area within the Polygon with ones
    img = Image.new('L', (x, y), 0)
    ImageDraw.Draw(img).polygon(vert, fill=1)

    # Convert the image object to a numpy object
    mask_sel = np.array(img)

    # Multiply the selected area by ANHA4 mask so we get back the land within the polygon.
    mask_new = mask_sel*mask

    # %%
    # Create a xarray dataset from the new mask.
    mask_ds = xr.Dataset(
                data_vars={
                    mask_var: (['lat', 'lon'], mask_new,
                               {'long_name': f'ANHA4_{mask_var}'})
                },
                coords={
                    "latitude": (["lat", "lon"], lat, {"long_name": "Latitude", "units": "degrees_north"}),
                    "longitude": (["lat", "lon"], lon, {"long_name": "Longitude", "units": "degrees_east"})
                },
                attrs={
                    "ANHALYZE": f"ANHA4 subregion {mask_var} created by using anhalyze create_mask tool. "
                                f"Refers to https://github.com/PORTAL-CEOS/ANHALYZE for more information."
                }
    )

    # Fix bug that creates gigantic values for some unknown reason.
    mask_ds.where(mask_ds == 1, 0)

    # %%
    # Save the new mask array as a new NetCDF file
    # Set destiny path to save the new mask file.
    if not path:
        path = anhalyze.PACKAGE_DATA_DIR+'/'

    # New mask file name
    old_mask_name = mask_source.split('/')[-1]
    new_mask_name = old_mask_name.replace('.nc', suffix)

    # Making sure the filename ends in .nc
    if '.nc' not in new_mask_name[-3:]:
        new_mask_name += '.nc'

    # Avoiding overwriting files by adding an extra suffix continuously until available.
    while os.path.isfile(new_mask_name):
        print(f'[Anhalyze Warning, file exists: {new_mask_name}]')
        new_mask_name = new_mask_name.replace('.nc', '_copy.nc')

    if verbose:
        print(f'[Anhalyze] Saving new mask: {path+new_mask_name}')

    # Saving new mask file
    mask_ds.to_netcdf(path+new_mask_name)
