.. _releases_rvtpy:

Python package release history
==============================

UNRELEASED
----------

2.2.3
-----
*   Added enhance Multi-Scale Topographic Position (e4MSTP) blending combination visualization.

RELEASE

Jul 4, 2025


2.2.2
-----
*   Replace np.cast with np.asarray .
*   Added BIGTIFF=IF_SAFER to tiff creation properties in default.py.
*   Added dependencies: pandas, geopandas, rasterio, jupyter.

RELEASE

Feb 26, 2025


2.2.1
-----
*   Fixed bug in blending.

RELEASE

May 23, 2023


2.2.0
-----
*   Added 1 pixel edge padding before the calculation of hillshade and slope, to avoid no data edge in the final results.
*   Added float option for MSTP visualization.
*   Changed default parameters of MSTP visualization.
*   Changed nodata handling when calculating slope visualization. When calculating slope in the specific pixel if any of
    neighbour pixels are nodata use middle value in the calculation instead.

RELEASE

May 22, 2023


2.1.0
-----
*   Changed 8bit (bytescale) parameters of some visualizations (all changed to value mode, to avoid tiling effect when using tile module).
*   Float to 8bit bug fix.

RELEASE

March 6, 2022


2.0.0
-----
*   Module multiproc.py replaced with tile.py.

RELEASE

February 5, 2022


1.0.0
-----
*   Added soft light blending mode.
*   Added Color Relief Image Map (CRIM) blending combination visualization.
*   Added enhance Multi-Scale Topographic Position (e3MSTP) blending combination visualization.
*   Fixed luminosity blending bug.
*   Fixed summed area table algorithm (used in SLRM, MSRM, MSTP) no data bug.

RELEASE

September 14, 2021


1.0.0a11
--------

*   Blending with Multiple directions hillshade bug fixed. Now MHS 8-bit is used.
*   Fixed SVF, ASVF, OPNS file name, added noise remove parameter in the output name.

PRE-RELEASE

May 20, 2021


1.0.0a10
--------

*   Added fill no data methods (IDW, Nearest neighbor, K-D Tree)
*   Fixed Negative Openness 8bit image (reverted colors).
*   Added Multi-scale topographic position (MSTP)

PRE-RELEASE

Apr 16, 2021


1.0.0a9
-------

PRE-RELEASE

Mar 9, 2021


1.0.0a8
-------

PRE-RELEASE

Jan 29, 2021


1.0.0a7
-------

PRE-RELEASE

Jan 19, 2021


1.0.0a6
-------

PRE-RELEASE

Jan 11, 2021


1.0.0a5
-------

PRE-RELEASE

Jan 10, 2021


1.0.0a4
-------

PRE-RELEASE

Jan 8, 2021


1.0.0a3
-------

PRE-RELEASE

Jan 8, 2021


1.0.0a2
-------

PRE-RELEASE

Jan 8, 2021


1.0.0a1
-------

PRE-RELEASE

Jan 8, 2021
