#!/usr/bin/env python
# coding: utf-8

import numpy as np
import xarray as xr
import anhalyze
import os


def get_ohc(temp_da, e3t_da, tref=0, rho_ref=1030, cp=4218):
    """
     Calculate the heat content relative to a temperature reference.
     """

    ohc_da = rho_ref * cp * np.sum((temp_da - tref) * e3t_da, axis=1)

    return ohc_da


def get_thermheig(temp_da, e3t_da, tref=0, sref=34.8, rho_ref=None, cp=4218):
    import gsw as gsw

    if not rho_ref:
        rho_ref = gsw.density.rho(sref, tref, 1)

    # Alpha thermal expansion coefficient under a reference salinity
    alpha = gsw.density.alpha(sref, temp_da, 1)

    # Calculate the depth integrated ocean heat content
    ohc_da = get_ohc(temp_da=temp_da, e3t_da=e3t_da, tref=tref, rho_ref=rho_ref, cp=cp)

    # Estimate the thermosteric height from the total heat content
    thermheig = np.nanmean(alpha / (rho_ref*rho_ref*cp), axis=1) * ohc_da

    return thermheig
#
#def get_e1t(lat_range=None, lon_range=None, config='ANHA4'):

#def get_e2t():

#def get_e3t(lat_range=None, lon_range=None, depth_range, config='ANHA4'):
