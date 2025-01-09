#!/usr/bin/env python3
# coding: utf-8

# System-related libraries
import requests
import os

# Project-related libraries
import anhalyze.config as config
import anhalyze as ah


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
