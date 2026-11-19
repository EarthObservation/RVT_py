.. _install_arcgis:

ArcGIS installation
===================

To use RVT in ArcGIS Pro, `download the ArcGIS Raster Functions repository <https://github.com/EarthObservation/rvt-arcgis-pro>`_ by selecting ``Code → Download ZIP``.

Unzip the downloaded repository folder and *rename* it to ``rvt-arcgis-pro``, then copy the whole repository folder to: ``<ArcGIS Pro install path>/Resources/Raster/Functions/Custom``

Usually the path is: ``c:/Program Files/ArcGIS/Pro/Resources/Raster/Functions/Custom``

For ArcGIS Server use, copy the whole repository folder (``rvt-arcgis-pro``) to every federated server machine of your enterprise setup: ``<ArcGIS Server install path>/framework/runtime/ArcGIS/Resources/Raster/Functions/Custom``

Open or restart ArcGIS Pro. Select ``Imagery → Raster Functions`` to open the Raster Functions pane. 

In the ``Raster Functions`` pane, select the ``Custom`` tab to access the ``rvt-arcgis-pro`` group containing the raster functions.

.. image:: ./figures/ArcGISPro_raster_functions.png
   :width: 450px

Importing from a Local Folder (No Admin Rights Required)
--------------------------------------------------------

1. Download the `RVT ArcGIS Pro repository <https://github.com/EarthObservation/rvt-arcgis-pro>`_ from GitHub by selecting **Code → Download ZIP**.
2. Extract the ZIP file to a folder on your computer.
3. Open ArcGIS Pro and select **Imagery → Raster Functions**.
4. In the **Raster Functions** pane, open the menu (three horizontal bars) and select **Add New Sub-Category**.
5. In the newly created subcategory, click the orange **Import Functions** icon (downward arrow).
6. Browse to the extracted ``rvt-arcgis-pro`` folder.
7. Select the desired RVT raster functions (``.rft.xml`` files) and import them.

The imported functions will appear in your custom subcategory.
