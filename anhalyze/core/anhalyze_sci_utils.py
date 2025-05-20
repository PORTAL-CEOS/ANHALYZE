#!/usr/bin/env python
# coding: utf-8

import numpy as np
import xarray as xr
import os
import anhalyze
import warnings


def get_ohc(temp_da, e3t_da, tref=0, rho_ref=1030, cp=4218):
    """
     Calculate the heat content relative to a temperature reference.
     """

    ohc_da = rho_ref * cp * np.sum((temp_da - tref) * e3t_da, axis=1)

    return ohc_da


def get_thermheig(temp_da, e3t_da, tref=0, sref=34.8, rho_ref=None, cp=4218):
    """
    Calculates thermosteric height relative to a reference temperature
    """
    import gsw as gsw

    if not rho_ref:
        rho_ref = gsw.density.rho(sref, tref, 1)

    # Alpha thermal expansion coefficient under a reference salinity
    alpha = gsw.density.alpha(sref, temp_da, 1)

    # Calculate the depth integrated ocean heat content
    ohc_da = get_ohc(temp_da=temp_da, e3t_da=e3t_da, tref=tref, rho_ref=rho_ref, cp=cp)

    # Estimate the thermosteric height from the total heat content
    thermheig = np.nanmean(alpha / (rho_ref * rho_ref * cp), axis=1) * ohc_da

    return thermheig


def get_fwc(sal_da, e3t_da, sref=34.8):
    """
    Calculate the freshwater content thickness relative to a salinity reference.
    """

    fwc_da = np.sum((((sref - sal_da) / sref) * e3t_da), axis=1)

    return fwc_da


def get_halosheig(sal_da, e3t_da, sref=34.8, tref=0, rho_ref=None):
    import gsw as gsw

    if not rho_ref:
        rho_ref = gsw.density.rho(sref, tref, 1)

    # Beta salinity contraction coefficient under a reference temperature
    beta = gsw.density.beta(sal_da, tref, 1)

    # Calculate the depth integrated  freshwater content
    fwc_da = get_fwc(sal_da=sal_da, e3t_da=e3t_da, sref=sref)

    # Estimate the thermosteric height from the total fwc content
    halosheig = np.nanmean(beta * sref, axis=1) * fwc_da

    return halosheig


def get_mke(gridu_file, varu, gridv_file, varv):
    """
    Calculate Mean Kinect Energy from gridU and gridV velocities variables.

    """

    # Getting the AnhaDataset object
    adsu = anhalyze.AnhaDataset(gridu_file)
    adsv = anhalyze.AnhaDataset(gridv_file)

    # Velocities components are oriented in the i and j directions within the model grid.
    # To calculate the actual Mean Kinect Energy, two steps must be done:
    # 1 - Average the components into the grid T
    # 2 - Rotate the vectors to change the orientation reference to the North Pole

    # Rotating the vector
    adsu, adsv = rot_vec(adsu, varu, adsv, varv)

    # Calculating Mean Kinect Energy
    varu_t_r_da = adsu._get_var_data_array(var=f"{varu}_t_r")
    varv_t_r_da = adsv._get_var_data_array(var=f"{varv}_t_r")

    # Calculating the Mean Kinect Energy
    mke_da = 0.5 * ((varu_t_r_da ** 2) + (varv_t_r_da ** 2))

    # Update variable name in the AnhaDataset _xr_dataset
    mke = adsu._xr_dataset.rename({varu: 'mke'})

    # Replacing the velocities values with mke values
    mke.data_vars['mke'].data = mke_da

    # Update variable unit
    mke.attrs['standard_name'] = 'sea_water_mke'
    mke.attrs['long_name'] = 'ocean mean kinect energy'
    mke.attrs['units'] = 'm².s⁻²'

    return mke


def rot_vec(adsu, varu_t, adsv, varv_t):
    """
    Rotating vectors from model grid to Earth surface reference

    varu_t : string
        Name of the ector variable component i-axis oriented
    varv_t : string
        Name of Vector variable component j-axis oriented
    """

    # Verify if variables are in grid T
    if 'gridu2gridt' not in adsu._xr_dataset[varu_t].attrs:
        warnings.warn(f"[anhalyze_sci_utils] Vector component {varu_t} not in grid T."
                      f"\nAnhalyze will average the variable into the grid T.")
        adsu = adsu.gridu2gridt(var=varu_t)
        varu_t = f'{varu_t}_t'

    if 'gridv2gridt' not in adsu._xr_dataset[varu_t].attrs:
        warnings.warn(f"[anhalyze_sci_utils] Vector component {varv_t} not in grid T."
                      f"\nAnhalyze will average the variable into the grid T.")
        adsv = adsv.gridv2gridt(var=varv_t)
        varv_t = f'{varv_t}_t'

    # Verify whether the grid configurations are the same for both components
    assert adsu.attrs['model_config'] == adsv.attrs['model_config'], \
        (f"[Anhalyze] Model configuration doesn't match for U ({adsu.attrs['model_config']})"
         f" and V ({adsv.attrs['model_config']}) components.")

    # Verify whether the inputs have the same dimension: in X
    assert adsu.dims['x'] == adsv.dims['x'], \
        (f"[Anhalyze] Input 'x' dims doesn't match for U ({adsu.dims['x']})"
         f" and V ({adsv.dims['x']}).")

    # Verify whether the inputs have the same dimension: in Y
    assert adsu.dims['y'] == adsv.dims['y'], \
        (f"[Anhalyze] Input 'y' dims doesn't match for U ({adsu.dims['y']})"
         f" and V ({adsv.dims['y']}).")

    # Verify whether the inputs have the same dimension: in Z
    assert adsu.dims['depthu'] == adsv.dims['depthv'], \
        (f"[Anhalyze] Input dims doesn't match for U ({adsu.dims['depthu']})"
         f" and V ({adsv.dims['depthv']}).")

    # Get copies of AnhaDataset information for both components
    _load_data_u, _load_data_v = adsu._load_data, adsv._load_data
    _xr_dataset_u, _xr_dataset_v = adsu._xr_dataset.copy(), adsv._xr_dataset.copy()

    # Check if any other method had already created '_attrs' attribute
    if not hasattr([adsu, adsv], '_attrs'):
        _attrs_u, _attrs_v = adsu.attrs.copy(), adsv.attrs.copy()
    elif not hasattr(adsu, '_attrs') & hasattr(adsv, '_attrs'):
        _attrs_u, _attrs_v = adsu.attrs.copy(), adsv._attrs.copy()
    elif not hasattr(adsv, '_attrs') & hasattr(adsu, '_attrs'):
        _attrs_u, _attrs_v = adsu._attrs.copy(), adsv.attrs.copy()
    else:
        _attrs_u, _attrs_v = adsu._attrs.copy(), adsv._attrs.copy()

    # Extract configuration information from either files since is the same in both
    config = adsu.attrs['model_config']

    # Get sin(θ) and cons(θ) variables for a given grid configuration
    rotangle = get_rotatedangle(config=config)
    gsint = rotangle['gsint']
    gcost = rotangle['gcost']

    # Slice the angle if necessary
    if hasattr(adsu, '_attrs') and 'regional' in adsu._attrs['file_catogory']:
        # Get all lat-lon data in file
        lat = rotangle['nav_lat'].data.copy()
        lon = rotangle['nav_lon'].data.copy()

        # Extract lat_range and lon_range from dataset
        lat_range = adsu.attrs['coord_lat_range']
        lon_range = adsu.attrs['coord_lon_range']

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

        # Slicing the sin and cos to match the input region selection
        gsint = gsint.isel({'x': slice(col_range[0], col_range[1]),
                            'y': slice(row_range[0], row_range[1])})

        gcost = gcost.isel({'x': slice(col_range[0], col_range[1]),
                            'y': slice(row_range[0], row_range[1])})

    # Rearrange angle array dimensions to match the AnhaDataset's
    gsint = gsint.expand_dims(dim={'time_counter': adsu.dims['time_counter'], 'deptht': adsu.dims['depthu']},
                              axis=[0, 1])
    gcost = gcost.expand_dims(dim={'time_counter': adsu.dims['time_counter'], 'deptht': adsu.dims['depthu']},
                              axis=[0, 1])

    # Rotate the components
    varu_t_rot = ((_xr_dataset_u[varu_t].data * gcost.data) -
                  (_xr_dataset_v[varv_t].data * gsint.data))
    varv_t_rot = ((_xr_dataset_v[varv_t].data * gcost.data) +
                  (_xr_dataset_u[varu_t].data * gsint.data))

    # Add to the data set as a new variable
    _xr_dataset_u[f'{varu_t}_r'] = (_xr_dataset_u[varu_t].dims, varu_t_rot)
    _xr_dataset_v[f'{varv_t}_r'] = (_xr_dataset_v[varv_t].dims, varv_t_rot)

    # Add a parameter to inform that these AnhaDataset components are rotated
    _xr_dataset_u[f'{varu_t}_r'].attrs = _xr_dataset_u[varu_t].attrs
    _xr_dataset_u[f'{varu_t}_r'].attrs['Rotated'] = True
    _xr_dataset_v[f'{varv_t}_r'].attrs = _xr_dataset_v[varv_t].attrs
    _xr_dataset_v[f'{varv_t}_r'].attrs['Rotated'] = True

    return (anhalyze.AnhaDataset('', load_data=_load_data_u, _xr_dataset=_xr_dataset_u, _attrs=_attrs_u),
            anhalyze.AnhaDataset('', load_data=_load_data_v, _xr_dataset=_xr_dataset_v, _attrs=_attrs_v))


def get_rotatedangle(config='ANHA4'):
    """ Get mask from given mask_filename or default location.

    Parameters
    ----------
    config : str
        Model configuration. Default: ANHA4

    """
    # Rotated angle file name. Each grid configuration have they own file.
    rotatedangle_filename = f'RotatedAngle_{config}.nc'

    # Get rotating file for a given configuration (ANHA4 is the default)
    rot_file = os.path.join(anhalyze.PACKAGE_DATA_DIR, rotatedangle_filename)

    # Extracting sin(θ) and cons(θ) variables
    rotangle = xr.open_dataset(rot_file)

    return rotangle

#
#def get_e1t(lat_range=None, lon_range=None, config='ANHA4'):

#def get_e2t():

#def get_e3t(lat_range=None, lon_range=None, depth_range, config='ANHA4'):
