ArcGIS Pro
==========

If you would like to use visualization functions in ArcGIS Pro you have to click Analysis->Raster functions->"three parallel lines"->Open Python Raster Function. Then select Python Module (rvt_esri_*.py) and Class Name (it is only one).

For rvt_esri_blender.py you will need to install rasterio into Python ArcGIS Pro conda environment. To do that open ArcGIS Pro click Python then click Manage Environments and Clone the default environment. After you clone default environment set new env to active and click ok. Then click Add Packages search for rasterio and install it. In case you are having problems try older version of rasterio.
