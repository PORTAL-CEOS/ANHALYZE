#!/usr/bin/env python
# coding: utf-8

from anhalyze import AnhaDataset
import anhalyze

# Class to work using temperature
class AnhalyzeTemp(AnhaDataset):
    def __init__(self, filename, load_data=True, mask_filename=None, _xr_dataset=None, _attrs=None):
        super().__init__(filename, load_data=True, mask_filename=None, _xr_dataset=None, _attrs=None)

    def get_ohc(self, tref=0, rho_ref=1030, cp=4218, load_data=None, _xr_dataset=None, attrs=None):
        """
        Calculate the heat content relative to a temperature reference.
        """

        import anhalyze.core.anhalyze_sci_utils as asu

        # Get DataArray for temperature variable and grid cell thickness
        temp_da = self._get_var_data_array(var='votemper')
        e3t_da = self._get_var_data_array(var='e3t')

        # Calculating OHC
        ohc = asu.get_ohc(
            temp_da=temp_da,
            e3t_da=e3t_da,
            tref=tref,
            rho_ref=rho_ref,
            cp=cp,
        )

        # Adding ohc data to data_vars
        self._xr_dataset = self._xr_dataset.assign({'ohc': ((self.attrs['dim_t'],
                                                             self.attrs['dim_y'],
                                                             self.attrs['dim_x']),
                                                            ohc.data)})

        self._xr_dataset.data_vars['ohc'].attrs = {
            'standard_name': 'heat_content_depth_integrated',
            'temperature_reference': f'{tref} ºC',
            'ocean_density_reference': f'{rho_ref} Kg.m⁻³',
            'heat_capacity': f'{cp} J.ºC⁻¹.Kg⁻¹',
            'long_name': 'depth_integrated_ocean_heat_content',
            'units': 'J.m⁻²',
            'online_operation': 'average',
            'interval_operation': '1080 s',
            'interval_write': '5 d',
            'cell_methods': 'time: mean (interval: 1080 s)'
        }

    def get_thermheig(self, tref=0, sref=34.8, cp=4218):
        """
        Calculates thermosteric height
        """
        import anhalyze.core.anhalyze_sci_utils as asu

        # Get DataArray for temperature variable and grid cell thickness
        temp_da = self._get_var_data_array(var='votemper')
        e3t_da = self._get_var_data_array(var='e3t')

        # Calculating Thermosteric Height
        thermheig = asu.get_thermheig(
            temp_da=temp_da,
            e3t_da=e3t_da,
            tref=tref,
            sref=sref,
            rho_ref=None,
            cp=cp
        )

        # Adding thermosteric height data to data_vars
        self._xr_dataset = self._xr_dataset.assign({'thermheig': ((self.attrs['dim_t'],
                                                                   self.attrs['dim_y'],
                                                                   self.attrs['dim_x']),
                                                                  thermheig.data)})

        self._xr_dataset.data_vars['thermheig'].attrs = {
            'standard_name': 'heat_content_depth_integrated',
            'temperature_reference': f'{tref} ºC',
            'salinity_reference': f'{sref}',
            'ocean_density_reference': f'{rho_ref} Kg.m⁻³',
            'heat_capacity': f'{cp} J.ºC⁻¹.Kg⁻¹',
            'long_name': 'depth_integrated_thermosteric_height',
            'units': 'm',
            'online_operation': 'average',
            'interval_operation': '1080 s',
            'interval_write': '5 d',
            'cell_methods': 'time: mean (interval: 1080 s)'
        }

#    TODO: def get_ohf_sec(self):

# Class to work using salinity
class AnhalyzeSal(AnhaDataset):
    def __init__(self, filename, load_data=True, mask_filename=None, _xr_dataset=None, _attrs=None):
        super().__init__(filename, load_data=True, mask_filename=None, _xr_dataset=None, _attrs=None)

    def get_fwc(self, sref=34.8, load_data=None, _xr_dataset=None, attrs=None):
        """
        Calculate the freshwater content relative to a salinity reference.
        """

        import anhalyze.core.anhalyze_sci_utils as asu

        # Get DataArray for salinity variable and grid cell thickness
        sal_da = self._get_var_data_array(var='vosaline')
        e3t_da = self._get_var_data_array(var='e3t')

        # Calculating FWC
        fwc = asu.get_fwc(
            sal_da=sal_da,
            e3t_da=e3t_da,
            sref=sref,
        )

        # Adding fwc data to data_vars
        self._xr_dataset = self._xr_dataset.assign({'fwc': ((self.attrs['dim_t'],
                                                             self.attrs['dim_y'],
                                                             self.attrs['dim_x']),
                                                            fwc.data)})

        self._xr_dataset.data_vars['fwc'].attrs = {
            'standard_name': 'freshwater_content_depth_integrated',
            'temperature_reference': f'{sref}',
            'long_name': 'depth_integrated_freshwater_content',
            'units': 'meters',
            'online_operation': 'average',
            'interval_operation': '1080 s',
            'interval_write': '5 d',
            'cell_methods': 'time: mean (interval: 1080 s)'
        }

    def get_halosheig(self, sref=34.8, tref=0, rho_ref=None):
        """
        Calculates halosteric height
        """
        import anhalyze.core.anhalyze_sci_utils as asu

        # Get DataArray for temperature variable and grid cell thickness
        sal_da = self._get_var_data_array(var='vosaline')
        e3t_da = self._get_var_data_array(var='e3t')

        # Calculating Thermosteric Height
        halosheig = asu.get_halosheig(
            sal_da=sal_da,
            e3t_da=e3t_da,
            sref=sref,
            tref=tref,
            rho_ref=None,
        )

        # Adding halossteric height data to data_vars
        self._xr_dataset = self._xr_dataset.assign({'halosheig': ((self.attrs['dim_t'],
                                                                   self.attrs['dim_y'],
                                                                   self.attrs['dim_x']),
                                                                  halosheig.data)})

        self._xr_dataset.data_vars['halosheig'].attrs = {
            'standard_name': 'freshwater_content_depth_integrated',
            'temperature_reference': f'{tref} ºC',
            'salinity_reference': f'{sref}',
            'ocean_density_reference': f'{rho_ref} Kg.m⁻³',
            'long_name': 'depth_integrated_halosteric_height',
            'units': 'm',
            'online_operation': 'average',
            'interval_operation': '1080 s',
            'interval_write': '5 d',
            'cell_methods': 'time: mean (interval: 1080 s)'
        }

    #TODO: def get_ohf_sec(self):
