.. _qgis:

QGIS Plugin
===========

.. # TODO Describe, web portal

QGIS plugin will be available at the QGIS Plugins web portal soon. Please stay tuned.

The development version of the plugin is available in `Github <https://github.com/EarthObservation/rvt-qgis>`_. There are also installation instructions.

If you are in hurry, just download the file `<https://github.com/EarthObservation/RVT_py/blob/master/rvt_qgis_plugin/Q-RVTpy.zip?raw=true>`_ and install it from QGIS Plugins menu. Select ``Plugins → Manage and Install Plugins`` and follow the instructions in ``Install from ZIP``.

Once the plugin is installed, it can be accessed from the ``Raster menu``. First, check that the plugin is available in ``Plugins → Manage and Install Plugins``.

   .. image:: ./figures/rvt_qgis_plugins.png

#. Open a DEM file to be visualized.

   .. image:: ./figures/rvt_qgis_dem.png

#. Select ``Raster → Relief Visualization Toolbox``.

   .. image:: ./figures/rvt_qgis_menu.png

#. Choose the options and visualizations.

   .. image:: ./figures/rvt_qgis_toolbox.png

#. Click ``Start``.

The visualizations are computed and displayed in the main QGIS window. Files are stored in the same folder as the input file.

   .. image:: ./figures/rvt_qgis_svf.png

The visualizations are described in :ref:`usage`.