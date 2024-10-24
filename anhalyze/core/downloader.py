#!/usr/bin/env python3
# coding: utf-8

# System-related libraries
import requests
import os

# Project-related libraries
import anhalyze.config as config
import anhalyze


def download_nc_file(download_url, download_destination):
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
    mask_destination = os.path.join(anhalyze.PACKAGE_DATA_DIR, 'ANHA4_mask.nc')

    # Downloading mask.
    print("[Anhalyze.Downloader] Downloading mask file.")
    download_nc_file(mask_url, mask_destination)


if __name__ == '__main__':

    download_mask()
