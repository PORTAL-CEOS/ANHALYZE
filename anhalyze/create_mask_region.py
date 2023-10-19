def create_mask_region(run_name=None, environ_paths=True)

    """
    Masking region from ANHA4 domain using Polygon Selector
    
    Interactively creates a Mask file drawing a polygon over the ANHA4 domain.
    The user can zoom in in the figure to get a more accurate polygon.
    
    """

    import numpy as np
    import netCDF4 as nc
    import cv2
    import anhalyze_utils as au
    
    import matplotlib.pyplot as plt
    from matplotlib.widgets import PolygonSelector


    # %%
    # Use the ANHALYZE get_paths
    #au.get_paths(run_name=run_name, environ_paths=True)

    # Mas path to test
    mask_file = "/mnt/storage0/luiz/ANALYSES/MASKS/ANHA4_mesh_zgr.nc"

    # Extract ANHA4 mask
    mask_anha4 = nc.Dataset(mask_file)
    tmask_anha4 = mask_anha4['tmask'][0][0]

    # %%
    # Create polygon

    fig, ax = plt.subplots()
    ax.contourf(tmask_anha4)
    fig.show()

    def onselect(vertices):
        pass

    selector = PolygonSelector(ax, onselect)

    print("Click on the figure to create a polygon.")
    print("Press the 'esc' key to start a new polygon.")

    # %%
    # Set everything out of the polygon as zeros
    
    # New variable with the vertices indices. Rounded because the
    # PolygonSelector returns float numbers.  
    vert = np.array([np.intc(np.round(vertices,0)) for vertices in selector.verts],'int32')
    
    # Create a new numpy array with the ANHA4 mask dimensions
    mask_region = np.zeros(tmask_anha4.shape,dtype=np.uint8)
    
    # Set every value within the polygon as 1
    cv2.fillConvexPoly(mask_region,vert,1) 
    
    # Multiply the masked region by the ANHA4 mask so we get back the land within 
    # the new mask
    mask = mask_region*tmask_anha4
    
    # Todo SAVE MASK AS NETCDF
    