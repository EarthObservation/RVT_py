.. _qgis_releases:

QGIS Plugin Release history
===========================

v 0.8.0
-------

*   Luminosity blending bug fix.
*   Added soft light blending mode.
*   Added enhanced Multi-Scale Topographic Position (e3MSTP).
*   Fixed summed area table (used in MSTP, SLRM, MSRM) bug when DEM contains nodata.
*   Removed fill no-data option from visualizations (is still available under Other tab).

v 0.7.1
-------

*   Buttons "Select all" and "Select none" to select/deselect all DEMs.


v 0.7.0
-------

*   Enabled qgis_process command line utility.


v 0.6.4
-------

*   Blending with Multiple directions hillshade bug fixed. Now MHS 8-bit is used.
*   Fixed SVF, ASVF, OPNS file name, added noise remove parameter in the output name.
*   Added Multi-scale topographic position (MSTP) processing algorithm (function).


v 0.6.3
-------

*   Negative-Openness 8bit reversed colors bug fix.
*   Added Multi-scale topographic position (MSTP) visualization.


v 0.6.2.1
---------

*   Changed RVT QGIS plugin documentation link in About.


v 0.6.2
-------

*   Added fill no-data methods (Inverse Distance Weighting, K-D Tree, Nearest Neighbour).


v 0.6.1
-------

*   Multi-scale relief model fixed and added back. Sky illumination is still not working as it should (will be fixed soon). Read the Docs sites (RVT python library rvt_py, RVT QGIS plugin, RVT ArcGIS raster fn) were merged into one.


v 0.5.3
-------

*   Multi-scale relief model and Sky illumination visualizations temporarily removed, because they don't work as they should. They will be fixed soon.

v 0.5.2
-------

*   8-bit no data values changed from 0 to 255 (white).
*   Added from osgeo import gdal to default module.

v 0.5.1
-------

*   Blender tab Blend images button position changed.

v 0.5.0
-------

*   Added Other tab where you can cut-off raster values, normalize raster and change raster to 8 bit.
*   Plugin saves all output raster files as LZW compressed GeoTIFFs (previously it was saving without compression).
