.. _Reading and saving raster:

Reading and saving raster data
==============================

For reading raster data (DEMs) from files (GeoTIFFs) to a numpy array we suggest using the ``rvt.default`` module (which uses ``gdal``).
You can also use ``rasterio``, ``gdal`` or any other module for reading and saving geo rasters.

To read rasters with ``rvt.default``, first import the module. Then call the function ``rvt.default.get_raster_arr()`` to get a dictionary with keys array (contains numpy array of raster), resolution (contains tuple(x resolution, y resolution)) and no_data (contains value of no_data). See example below:

.. code-block:: python

    import rvt.default

    dem_path = r"C:/data/dem.tif"  # change path to your GeoTIFF
    dem_dict = rvt.default.get_raster_arr(dem_path)  # returns dictionary: {"array": array, "resolution": (x_res, y_res), "no_data": no_data}
    dem_arr = dem_dict["array"]  # numpy array
    dem_resolution_tuple = dem_dict["resolution"]  # resolution tuple (x direction resolution, y direction resolution)
    dem_x_resolution = dem_resolution_tuple[0]  # first element of resolution tuple is x direction resolution
    dem_y_resolution = dem_resolution_tuple[1]  # second element of resolution tuple is y direction resolution
    dem_no_data = dem_dict["no_data"]  # returns value of no_data stored in DEM


To save a raster with ``rvt.default``, import the module then call the function ``rvt.default.save_raster()``. You have to define the function parameters: ``src_raster_path``: source raster path (dem_path) to copy metadata, ``out_raster_path``: path to new file (visualization tif), ``out_raster_arr``: vizualization numpy array, ``no_data``: value of no_data (visualizations return no data as np.nan).

For example, the example below would compute a hillshade (with rvt.vis) in hillshade_arr from a DEM stored in dem_path location, storing this hillshade visualization to hillshade_path:

.. code-block:: python

    import rvt.default
    import numpy as np

    rvt.default.save_raster(src_raster_path=dem_path, out_raster_path=hillshade_path, out_raster_arr=hillshade_arr, no_data=np.nan)
