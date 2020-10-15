.. _arcgis:

ArcGIS Pro
==========

If you would like to use visualization functions in ArcGIS Pro you have to activate it.

#. Select ``Imagery → Raster Functions``

   .. image:: ./figures/rvt_esri_ribbon.png

#. Select→ ≡ → Open Python Raster Function``.

   .. image:: ./figures/rvt_esri_toolbox.png

#. Then select Python Module ``rvt_esri_*.py`` and Class Name (it is only one).

   .. image:: ./figures/rvt_esri_open.png

#. In toolbox set the function parameters.

   .. image:: ./figures/rvt_esri_menu.png

#. Click ``Create new layer``

The result is computed ad displayed in viewer.

.. image:: ./figures/rvt_esri_svf.png

The modules are the same as described in :ref:`usage`.

For ``rvt_esri_blender.py`` you will need to install ``rasterio`` into Python ArcGIS Pro conda environment. To do that open ArcGIS Pro, click ``Project → Python`` then click ``Manage Environments`` and clone the default environment. After that activate the environment and click ``OK``. Then click ``Add Packages`` search for ``rasterio`` and install it. In case you are having problems try older version of rasterio.
