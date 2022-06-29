.. _Reading and saving raster:

Reading and saving raster data
==============================

For reading raster data (DEMs) from files (GeoTIFFs) to a numpy array we suggest using the ``rvt.default`` module (which uses ``gdal``).
You can also use ``rasterio``, ``gdal`` or any other module for reading and saving geo rasters.

Reading raster data
-------------------

**Example**

To read a raster with ``rvt.default``:

.. code-block:: python

    # import the module
    import rvt.default
    
    # change path to your GeoTIFF
    dem_path = r"C:/data/dem.tif"  
    
    # call the function rvt.default.get_raster_arr() to return a dictionary with keys:
    # array (contains numpy array of raster),
    # resolution (contains the tuple(x resolution, y resolution)),
    # no_data (contains the value of no_data)
    dem_dict = rvt.default.get_raster_arr(dem_path) 
    
    # create numpy array
    dem_arr = dem_dict["array"]
    
    # the resolution tuple (x-direction resolution, y-direction resolution)
    dem_resolution_tuple = dem_dict["resolution"] 
    
    # the first element of the resolution tuple (the x-direction resolution)
    dem_x_resolution = dem_resolution_tuple[0]  
    
    # the second element of resolution tuple (the y-direction resolution)
    dem_y_resolution = dem_resolution_tuple[1] 
    
    # the value of no_data stored in DEM
    dem_no_data = dem_dict["no_data"] 

Saving raster data
------------------

**Example**

Let's say we wanted to use a DEM stored in ``dem_path`` to compute (using ``rvt.vis``) a hillshade stored in ``hillshade_arr``, and then save this hillshade visualization to ``hillshade_path``:

.. code-block:: python

    # import the required modules
    import rvt.default
    import numpy as np
    
    # call the function rvt.default.save_raster() and define the function parameters:
    # src_raster_path: source raster path (dem_path) to copy metadata, 
    # out_raster_path: path to new file (visualization tif), 
    # out_raster_arr: vizualization numpy array, 
    # no_data: value of no_data (visualizations return no data as np.nan)
    rvt.default.save_raster(src_raster_path=dem_path, out_raster_path=hillshade_path, out_raster_arr=hillshade_arr, no_data=np.nan)
    
.. seealso:: Find out more about defining default values in :ref:`rvt.default`.
