.. _qgis_releases:

QGIS Plugin Release history
===========================

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
