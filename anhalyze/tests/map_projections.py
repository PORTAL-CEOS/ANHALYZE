

# Library imports
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import unittest


class MyTestCase(unittest.TestCase):
    def test_cartopy(self):
        """ Testing if current installation can make a plot using cartopy projections. """

        # self.assertEqual(True, False)  # add assertion here

        # Arbitrary values to plot hudson bay
        hudson_east = -75
        hudson_west = -95
        hudson_north = 65
        hudson_south = 50
        standard_parallels = (55, 60)
        central_longitude = -80

        # Set up plot
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1,
                             projection=ccrs.LambertConformal(central_longitude=central_longitude,
                                                              standard_parallels=standard_parallels))
        ax.set_extent([hudson_west, hudson_east, hudson_south, hudson_north])
        ax.stock_img()
        ax.coastlines()
        ax.gridlines()

        ax.scatter([-80], [55], np.array([100]), color='Orange', transform=ccrs.PlateCarree())

        # Adding grid-line labels
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False)
        gl.right_labels = gl.top_labels = False


if __name__ == '__main__':
    unittest.main()
