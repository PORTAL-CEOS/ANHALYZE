#!/usr/bin/env python
# coding: utf-8


# TODO can make these sections into its own class
def init_location(hudson_bay=True):
    """
    Setting up model parameters based on Hudson Bay and James Bay locations.

    Hudson_bay : Boolean if using Hudson Bay vs James Bay locations.
    """

    if hudson_bay:
        # proj_size = 2
        region = 'HudsonBay'

        # Setting up James Bay location
        east = -75
        west = -93
        north = 65
        south = 50
        standard_parallels = (52.5, 62.5)  # TODO could calculate this instead
        central_longitude = -80  # TODO Could calculate from range (like mid point), or as keyword in functions

        lat_range = (south, north)
        lon_range = (west, east)

    else:

        # proj_size = 4
        region = 'JamesBay'

        # Setting up James Bay location
        east = -78.5
        west = -82.5
        north = 54.7
        south = 51
        standard_parallels = (52, 53)
        central_longitude = -80

        lat_range = (south, north)
        lon_range = (west, east)

    location_info = {'lat_range': lat_range,
                     'lon_range': lon_range,
                     'region': region,
                     'standard_parallels': standard_parallels,
                     'central_longitude': central_longitude,
                     # 'proj_size': proj_size,
                     }

    return location_info
