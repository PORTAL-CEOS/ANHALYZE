#!/usr/bin/env python
# coding: utf-8

import deprecation


# TODO can make these sections into its own class
@deprecation.deprecated()
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


# TODO can make these sections into its own class
@deprecation.deprecated()
def getIndex_sec(sectName):
    """Return the number associated to each section.
    One can check the respective numbers and sections on the Lab Guide.

    Inputs
    sectName = Section name; type: String

    Options:
    Bering Strait
    Lancaster Sound
    Jones Sound
    Nares Strait
    Davis Strait
    Fram Strait

"""

# Todo, add the tilted sections. So far, the section functions only work with
# zonal or meridional oriented sections.

    if sectName == 'Bering Strait':
        sect = 3180
    elif sectName=='Lancaster Sound':
        sect = 4180
    elif sectName=='Jones Sound':
        sect = 4180
    elif sectName=='Nares Strait':
        sect = 4270
    elif sectName=='Davis Strait':
        sect = 9360
    elif sectName=='Fram Strait':
        sect = 5360
    else:
        sect = None
        print("Section name not found.\nCall help/doc to check the sections available in this function.")

    return sect
