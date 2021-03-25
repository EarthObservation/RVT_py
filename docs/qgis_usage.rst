.. _qgis_usage:

Using RVT
=========

#. Open a DEM file to be visualized.

   .. image:: ./figures/rvt_qgis_dem.png

#. Select ``Raster → Relief Visualization Toolbox`` or click the RVT icon.

   .. image:: ./figures/rvt_qgis_menu.png

Computing visualizations
------------------------

#. Chose the DEM in ``List of currently selected files:``, then choose the Visualizations tab. In the Visualization tab select preferred visualizations and set their parameters (options).

   .. image:: ./figures/rvt_qgis_toolbox.png

#. Click ``Start`` to calculate the visualizations.

The visualizations are stored as GeoTIFFs in the same folder as the input file or to a custom location (if ``Save to raster location`` check box is unchecked and directory is set in the line edit next to it).
Visualizations are also added to the main window of QGIS  if ``Add to QGIS`` check box is checked. If the ``Overwrite`` check box is checked, the program overwrites the already existing visualization files.

   .. image:: ./figures/rvt_qgis_svf.png

The visualizations are described in the `RVT Python core library documentation <https://rvt-py.readthedocs.io>`_.

Using the blender
-----------------

#. Chose a DEM in ``List of currently selected files:``, then choose the Blender tab. In the Blender tab select your ``Blend combination:`` or build your own in layers.

You can add your own custom combination to the list. Write a name in the ``Combination name`` line edit and click ``Add``. To remove it, just select it in the (``Blend combination`` list) and click ``Remove``.
You can also save a specific combination to a JSON file (if you want to share it, for example). To do this, input its name (in the ``Combination name`` line edit) and click ``Save ...`` (then select the location and name of the file).
Saved JSON combinations can be added by clicking the ``Load ...`` button (select file). You can change the parameters for each visualization method in the blend combination in the Visualizations tab.
If you check the ``Use preset values for terrain type`` it applies the selected terrain type settings (this changes the normalization min and max, and visualizations parameters). If you check ``Save visualizations``, all the visualization parameters used in the blender combination will be saved.

   .. image:: ./figures/rvt_qgis_blender.png

#. Click ``Blend images`` to calculate the blended (fused, combined) image.

   .. image:: ./figures/rvt_qgis_vat.png

The blended image is stored as a GeoTIFF in the same folder as the input file or to a custom location (if ``Save to raster location`` check box is unchecked and a directory is set in the line edit next to it).

Using the processing functions
------------------------------

#. In QGIS go to the ``Processing Toolbox → Relief visualization toolbox`` where all the Relief Visualization Toolbox visualization functions are accessible.

   .. image:: ./figures/rvt_qgis_processing_toolbox.png
