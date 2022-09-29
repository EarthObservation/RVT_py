.. _releases_qgis:

QGIS plugin release history
===========================

v0.9.4
------
*   Fixed paths to processing functions.

v0.9.3
------
*   Added MSTP float and 8bit option.

v0.9.2
------
*   Try to install scipy if it doesn't exist.
*   Added processing functions for filling no-data.


v0.9.1
------
*   Added 1 pixel edge padding before the calculation of hillshade and slope, to avoid no data edge in the final results.


v0.9.0
------
*   Added tiling module, visualizations on huge rasters are now calculated tile by tile.
*   Changed bytescale to 8bit parameters of all the visualizations to value mode (value ranges are different on each tile, this is why percent mode is not suitable).


v0.8.1
------
*   Hillshade negative values set to 0.
*   Changed vertical exaggeration factor limit from [-1000, 1000] to [-10000, 10000].


v0.8.0
------

*   Luminosity blending bug fix.
*   Added soft light blending mode.
*   Added enhanced Multi-Scale Topographic Position (e3MSTP).
*   Fixed summed area table (used in MSTP, SLRM, MSRM) bug when DEM contains nodata.
*   Removed fill no-data option from visualizations (is still available under Other tab).

v0.7.1
------

*   Buttons "Select all" and "Select none" to select/deselect all DEMs.


v0.7.0
------

*   Enabled qgis_process command line utility.


v0.6.4
------

*   Blending with Multiple directions hillshade bug fixed. Now MHS 8-bit is used.
*   Fixed SVF, ASVF, OPNS file name, added noise remove parameter in the output name.
*   Added Multi-scale topographic position (MSTP) processing algorithm (function).


v0.6.3
------

*   Negative-Openness 8bit reversed colors bug fix.
*   Added Multi-scale topographic position (MSTP) visualization.


v0.6.2.1
--------

*   Changed RVT QGIS plugin documentation link in About.


v0.6.2
------

*   Added fill no-data methods (Inverse Distance Weighting, K-D Tree, Nearest Neighbour).


v0.6.1
------

*   Multi-scale relief model fixed and added back. Sky illumination is still not working as it should (will be fixed soon). Read the Docs sites (RVT python library rvt_py, RVT QGIS plugin, RVT ArcGIS raster fn) were merged into one.


v0.5.3
------

*   Multi-scale relief model and Sky illumination visualizations temporarily removed, because they don't work as they should. They will be fixed soon.

v0.5.2
------

*   8-bit no data values changed from 0 to 255 (white).
*   Added from osgeo import gdal to default module.

v0.5.1
------

*   Blender tab Blend images button position changed.

v0.5.0
------

*   Added Other tab where you can cut-off raster values, normalize raster and change raster to 8 bit.
*   Plugin saves all output raster files as LZW compressed GeoTIFFs (previously it was saving without compression).
