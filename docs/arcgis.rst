.. _arcgis:

ArcGIS Pro
==========

Relief Visualisation Toolbox can be used in ArcGIS Pro as Raster Functions.

As described in `Raster functions—ArcGIS Pro | Documentation <https://pro.arcgis.com/en/pro-app/help/data/imagery/raster-functions.htm>`_, raster functions are operations that apply processing directly to the pixels of imagery and raster datasets, as opposed to geoprocessing tools, which write out a new raster to disk. Calculations are applied to the pixels of the original data as displayed, so only pixels that are visible on your screen are processed. As you zoom and/or pan around, the calculations are performed on the fly. Because no intermediate datasets are created, processes can be applied quickly, as opposed to the time it would take to create a processed file on disk.

#. Start ArcGIS Pro and open a project. Add a DEM layer.

   .. image:: ./figures/ArcGISPro_dem.png

#. Select ``Imagery → Raster Functions``.

   .. image:: ./figures/ArcGISPro_ribbon_raster_functions.png

#. Select ``Custom`` and open the group ``rvt-arcgis-pro``.

   .. image:: ./figures/ArcGISPro_raster_functions.png

#. Select the appropriate function, e.g. ``svf``. Specify processing parameters. Click ``Create new layer``.

   .. image:: ./figures/ArcGISPro_svf_parameters.png

#. The visualization is computed and displayd as layer in main window.

   .. image:: ./figures/ArcGISPro_result_svf.jpg

The visualizations are described in the `Relief Visualization Toolbox — Relief Visualization Toolbox in Python documentation <https://rvt-py.readthedocs.io>`_.
