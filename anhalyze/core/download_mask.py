#!/usr/bin/env python3
# coding: utf-8

import requests

MASK_URL = 'https://bit.ly/3TkYlGL'


def download_large_file(url, destination):
    """ Download a large file given url into destination


    :param url: str
        URL location
    :param destination: str
        Path and filename where to save file.
    :return:
        Saves file in destination.
    Original code from:
    https://www.geeksforgeeks.org/how-to-download-large-file-in-python-with-requests/

    """

    try:
        with requests.get(url, stream=True) as response:
            response.raise_for_status()
            with open(destination, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
        print("File downloaded successfully!")
    except requests.exceptions.RequestException as e:
        print("Error downloading the file:", e)


def test_download():
    """ Testing download_large_file
    """

    # TODO update to smaller file and move to tests
    # Tested with this 3Gb file, and it works well.
    url = 'https://releases.ubuntu.com/20.04.4/ubuntu-20.04.4-desktop-amd64.iso'

    destination = 'ubuntu-20.04.4-desktop-amd64.iso'

    download_large_file(url, destination)


if __name__ == '__main__':

    # URL location for mask.
    url = MASK_URL

    # Filename
    destination = 'ANHA4_mask_test.nc'

    download_large_file(url, destination)
