.. _arcgis:

ArcGIS Pro
==========

Relief Visualisation Toolbox can be used in ArcGIS Pro as raster functions.

As described in `Raster functions—ArcGIS Pro | Documentation <https://pro.arcgis.com/en/pro-app/help/data/imagery/raster-functions.htm>`_, raster functions are operations that apply processing directly to the pixels of imagery and raster datasets, as opposed to geoprocessing tools, which write out a new raster to disk. Calculations are applied to the pixels of the original data as the raster is displayed, so only pixels that are visible on your screen are processed. As you zoom and pan around, the calculations are performed on the fly. Since no intermediate datasets are created, processes can be applied quickly, as opposed to the time it would take to create a processed file on disk.

If you would like to use use Relief Visualization Toolbox in ArcGIS Pro you have to activate each function individually.

#. Select ``Imagery → Raster Functions``

   .. image:: ./figures/rvt_esri_ribbon.png

#. Select ``≡ → Open Python Raster Function``.

   .. image:: ./figures/rvt_esri_toolbox.png

#. Then select Python Module ``rvt_esri_*.py`` and Class Name (there is only  one class in every function).

   .. image:: ./figures/rvt_esri_open.png

#. In toolbox set the function parameters.

   .. image:: ./figures/rvt_esri_menu.png

#. Click ``Create new layer``

The result is computed ad displayed in viewer.

.. image:: ./figures/rvt_esri_svf.png

The modules are the same as described in :ref:`usage`.

For ``rvt_esri_blender.py`` you will need to install ``rasterio`` into Python ArcGIS Pro conda environment. To do that open ArcGIS Pro, click ``Project → Python`` then click ``Manage Environments`` and clone the default environment. After that activate the environment and click ``OK``. Then click ``Add Packages`` search for ``rasterio`` and install it. In case you are having problems try older version of rasterio.
