#!/usr/bin/env python3
# coding: utf-8

# System-related libraries
import requests
import os

# Project-related libraries
import anhalyze.config as config
import anhalyze as ah
from anhalyze.core.anhalyze import get_date


def download_sharepoint_file(download_url, download_destination):
    """ Download a large file given url into destination

        Parameters
        ----------
        download_url : str
            URL location
        download_destination : str
            Path and filename where to save file.

        Return
        ----------
        Saves file in destination.

    """

    # Code from:
    # https://www.geeksforgeeks.org/how-to-download-large-file-in-python-with-requests/
    try:
        with requests.get(download_url, stream=True) as response:
            response.raise_for_status()
            with open(download_destination, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
        print("[Anhalyze.Downloader] File downloaded successfully!")
    except requests.exceptions.RequestException as e:
        print("[Anhalyze.Downloader] Error downloading file:", e)


def download_mask():
    """ Downloads standard mask.
    """

    # URL location for mask.
    mask_url = config.package_data['mask']['url']

    # Download mask to standard location.
    mask_destination = os.path.join(ah.PACKAGE_DATA_DIR, 'ANHA4_mask.nc')

    # Downloading mask.
    print("[Anhalyze.Downloader] Downloading mask file.")
    download_sharepoint_file(mask_url, mask_destination)


def download_example(file_type='gridT'):
    """ Downloads Anha file example.

        Parameters
        ----------
        file_type : str
            File type either 'gridT'(default) or 'icebergs'

    """

    assert file_type in ['gridT', 'icebergs'], '[Anhalyze.Downloader] Incorrect file_type.'

    # URL location.
    file_url = config.package_data[file_type]['url']

    # Download mask to standard location.
    file_destination = os.path.join(ah.PACKAGE_DATA_DIR, config.package_data[file_type]['filename'])

    # Downloading mask.
    print(f"[Anhalyze.Downloader] Downloading {file_type} file.")
    download_sharepoint_file(file_url, file_destination)

    # Check filename reflects correct ymd
    test_filename(file_destination)


def test_filename(file_destination):
    """ Checks filename contains correct ymd values.
    """

    # Extract filename
    filename = os.path.basename(file_destination)
    # Open file
    file_example = ah.AnhaDataset(file_destination)

    # Get ymd values from file
    year = file_example._xr_dataset.coords['time_counter'].data[0].year
    month = file_example._xr_dataset.coords['time_counter'].data[0].month
    # Filename shows day at end of timestep, file shows at middle of timestep
    day = file_example._xr_dataset.coords['time_counter'].data[0].day + 2

    # Assert values in filename correspond to values in file
    assert get_date(file_example.attrs['filename'], how='y') == year, \
        f'[Anhalyze.Downloader] Filename {filename} does not contain year value found in file: {year}.'
    assert get_date(file_example.attrs['filename'], how='m') == month, \
        f'[Anhalyze.Downloader] Filename {filename} does not contain month value found in file: {month}.'
    assert get_date(file_example.attrs['filename'], how='d') == day, \
        f'[Anhalyze.Downloader] Filename {filename} does not contain day value found in file: {day}.'


def download_tutorial():
    """ Downloads tutorial in html.
    """

    # Download tutorial to standard location
    tutorial_destination = os.path.join(ah.PACKAGE_DATA_DIR.replace('package_data', 'tutorials'),
                                        'anhalyze_tutorial.html')

    # Downloading tutorial
    print(f"[Anhalyze.Downloader] Downloading tutorial version: {config.package_data['tutorial']['version']}, "
          f"here:{tutorial_destination}.")
    download_sharepoint_file(config.package_data['tutorial']['url'], tutorial_destination)


if __name__ == '__main__':

    download_mask()
