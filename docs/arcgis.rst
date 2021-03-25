.. _arcgis:

ArcGIS Pro
==========

Relief Visualisation Toolbox can be used in ArcGIS Pro as Raster Functions.

As described in `Raster functions—ArcGIS Pro | Documentation <https://pro.arcgis.com/en/pro-app/help/data/imagery/raster-functions.htm>`_, raster functions are operations that apply processing directly to the pixels of imagery and raster datasets, as opposed to geoprocessing tools, which write out a new raster to disk. Calculations are applied to the pixels of the original data as displayed, so only pixels that are visible on your screen are processed. As you zoom and/or pan around, the calculations are performed on the fly. Because no intermediate datasets are created, processes can be applied quickly, as opposed to the time it would take to create a processed file on disk.

The development version of the plugin is available in `Github EarthObservation/rvt-arcgis-pro <https://github.com/EarthObservation/rvt-arcgis-pro>`_.

.. toctree::

    arcgis_install
    arcgis_usage
