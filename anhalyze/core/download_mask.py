#!/usr/bin/env python3
# coding: utf-8

import requests


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
    url = 'https://drive.google.com/file/d/17f6Lb0qy_mB5zSw2ZVEZt8PiqcnFu-S7/view?usp=sharing'
    # url2 = 'https://umanitoba-my.sharepoint.com/personal/dasilvlh_myumanitoba_ca/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fdasilvlh%5Fmyumanitoba%5Fca%2FDocuments%2FANHALIZE%2FRESOURCES%2FANHA4%5Fmask%2Enc&parent=%2Fpersonal%2Fdasilvlh%5Fmyumanitoba%5Fca%2FDocuments%2FANHALIZE%2FRESOURCES&ga=1'
    # url2 = 'https://umanitoba-my.sharepoint.com/personal/dasilvlh_myumanitoba_ca/_layouts/15/download.aspx?UniqueId=a07a5b33-7508-41a7-ad89-d63e84169096'

    # Filename
    destination = 'ANHA4_mask_test.nc'

    download_large_file(url, destination)
