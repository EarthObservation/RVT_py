.. # Relief Visualization Toolbox in Python documentation master file, created by
   sphinx-quickstart on Mon Oct 12 08:11:32 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Relief Visualization Toolbox
============================

.. image:: /figures/RVT_head.png

Relief Visualization Toolbox was produced to help scientist visualize raster elevation model datasets. We have narrowed down the selection to include techniques that have proven to be effective for identification of small scale features. Default settings therefore assume working with high resolution digital elevation models, derived from airborne laser scanning missions (lidar).

Despite this, techniques are also used for different other purposes. Sky-view factor, for example, can be efficiently used in numerous studies where digital elevation model visualizations and automatic feature extraction techniques are indispensable, e.g. in geography, archeology,  geomorphology, cartography, hydrology, glaciology, forestry and disaster management. It can be used even in engineering applications, such as, predicting the availability of the GPS signal in urban areas.

Methods currently implemented are:

*   hillshading,
*   hillshading from multiple directions,
*   slope gradient,
*   simple local relief model,
*   sky illumination,
*   sky-view factor (as developed by our team),
*   anisotropic sky-view factor,
*   positive and negative openness,
*   local dominance.

.. image:: /figures/RVT_visualisations.jpg

For a more detailed description see references given at each method in the manual and a comparative paper describing them (e.g. Kokalj and Hesse 2017, see below).

RVT python library called rvt contains 3 modules: vis (rvt.vis), blend (rvt.blend) and default (rvt.default). Modules contains:
* vis       -   visualization functions (mentioned above), for computing visualizations;
* blend     -   blender (mixer), for blending visualizations;
* default   -   default values, class for defining default parameters with methods to compute and save visualization functions using set parameters.

References
----------

When using the tools, please cite:

*   Kokalj, Ž., Somrak, M. 2019. Why Not a Single Image? Combining Visualizations to Facilitate Fieldwork and On-Screen Mapping. Remote Sensing 11(7): 747.
*   Zakšek, K., Oštir, K., Kokalj, Ž. 2011. Sky-View Factor as a Relief Visualization Technique. Remote Sensing 3: 398-415.
*   Kokalj, Ž., Zakšek, K., Oštir, K. 2011. Application of Sky-View Factor for the Visualization of Historic Landscape Features in Lidar-Derived Relief Models. Antiquity 85, 327: 263-273.

Contributing
------------

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please report any bugs and suggestions for improvements.

Acknowledgment
--------------

Development of RVT was part financed by the European Commission's Culture Programme through the ArchaeoLandscapes Europe project and by the Slovenian Research Agency core funding No. P2-0406, and by research projects No. J6-7085 and No. J6-9395.

License
-------

This project is licensed under the terms of the `Apache License <http://www.apache.org/licenses/LICENSE-2.0>`_.


Copyright 2020 ZRC SAZU and University of Ljubljana

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   install
   usage
   examples
   arcgis
