"""
Masking region from ANHA4 domain using Polygon Selector
================


"""

import numpy as np
import netCDF4 as nc
import cv2

import matplotlib.pyplot as plt
from matplotlib.widgets import PolygonSelector


# %%
# to do: Use the ANHALYZE get_paths
mask_file = "/mnt/storage0/luiz/ANALYSES/MASKS/ANHA4_mesh_zgr.nc"

# Extract ANHA4 mask
mask_anha4 = nc.Dataset(mask_file)
tmask_anha4 = mask_anha4['tmask'][0][0]



# %%
# Create polygon

fig, ax = plt.subplots()
ax.contourf(tmask)
fig.show()

def onselect(vertices):
    pass

selector = PolygonSelector(ax, onselect)

print("Click on the figure to create a polygon.")
print("Press the 'esc' key to start a new polygon.")

# %%
# Set everything out of the polygon as zeros

# New variable with the vertices indices. Rounded because the
# PolygonSelector returns float numbers 
vert = [tuple(np.round(vertices,0)) for vertices in selector.verts]

# Create a new numpy array with the ANHA4 mask dimensions
mask_region = np.zeros(mask_anha4.shape) # to do

# Set every value within the polygon as 1
cv2.fillPoly(mask_region,polys,1)
mask_region = mask_region.astype(bool)

# to do: Multiply the masked region by the ANHA4 mask

plt.imshow(mask_region)
