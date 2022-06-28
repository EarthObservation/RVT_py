.. _install:

Installation
============

Requirements
------------

Required libraries (specified versions have been tested, other versions may also work):

*   numpy 1.19.2
*   scipy 1.5.2
*   gdal 3.0.2

We recommend using Python 3.6 or higher and a conda environment (this works best with gdal).

Conda
-----

The rvt_py package is available from the `Anaconda Cloud repository <https://anaconda.org/zmigyyy/rvt_py>`_. Using conda to install the rvt_py package will include all required libraries.

To use this method, first `install Anaconda and conda <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_.

Then open Anaconda Prompt (Windows) or Terminal (MacOS) and run:

``conda install -c zmigyyy rvt_py``

Pypi
----

Another option is to install the rvt-py package and required libraries using the `Python Package Index <https://pypi.org/project/rvt-py>`_ (pypi).

Pypi usually has problems installing gdal, so `install gdal first <https://pypi.org/project/GDAL/>`_ to use this method.

Then open Command Prompt (Windows) or Terminal (MacOS) and run:

``pip install rvt-py``

ArcGIS
------

To use RVT in ArcGIS Pro, download the `ArcGIS Raster Functions repository <https://github.com/EarthObservation/rvt-arcgis-pro>`_ by selecting ``Code → Download ZIP``.

Unzip the downloaded repository folder and rename it to ``rvt-arcgis-pro``, then copy the whole repository folder to:

``<ArcGIS Pro install path>/Resources/Raster/Functions/Custom``

Usually the path is:

``c:/Program Files/ArcGIS/Pro/Resources/Raster/Functions/Custom``

For ArcGIS Server use, copy the whole repository folder (``rvt-arcgis-pro``) to every federated server machine of your enterprise setup:

``<ArcGIS Server install path>/framework/runtime/ArcGIS/Resources/Raster/Functions/Custom``

Open or restart ArcGIS Pro. Select ``Imagery → Raster Functions`` to open the Raster Functions pane.

In the ``Raster Functions`` pane, select the ``Custom`` tab to access the ``rvt-arcgis-pro`` group containing the raster functions.

   .. image:: ./figures/ArcGISPro_raster_functions.png

QGIS
----

The RVT QGIS plugin uses the RVT python (core) library. To use RVT in QGIS, first open QGIS and select ``Plugins → Manage and Install Plugins → All``.

Search for ``RVT`` or ``Relief Visualization Toolbox`` and select ``Install Plugin``.

Once the plugin is installed, it can be accessed from the ``Raster menu`` or by selecting the icon that should appear on the toolbar. All the visualization functions are also available as processing functions in the ``Processing toolbox``.

The RVT QGIS plugin can also be downloaded from the `QGIS Python Plugins Repository <https://plugins.qgis.org/plugins/rvt-qgis/>`_.

   .. image:: ./figures/rvt_qgis_plugins.png
 
.. toctree::

    install
