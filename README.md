# Relief Visualization Toolbox in Python

![](./docs/figures/RVT_head.png)

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

For a more detailed description see references given at each method in the manual and a comparative paper describing them (e.g. Kokalj and Hesse 2017, see below).

RVT python library called rvt contains 3 modules: vis (rvt.vis), blend (rvt.blend) and default (rvt.default). Modules contains:
* vis       -   visualization functions (mentioned above), for computing visualizations;
* blend     -   blender (mixer), for blending visualizations;
* default   -   default values, class for defining default parameters with methods to compute and save visualization functions using set parameters.

For every visualization function directory also contains Python Esri raster functions for ArcGIS Pro (rvt_esri_*.py).

## References

When using the tools, please cite:

*   Kokalj, Ž., Somrak, M. 2019. Why Not a Single Image? Combining Visualizations to Facilitate Fieldwork and On-Screen Mapping. Remote Sensing 11(7): 747.
*   Zakšek, K., Oštir, K., Kokalj, Ž. 2011. Sky-View Factor as a Relief Visualization Technique. Remote Sensing 3: 398-415.

## Installation

### conda

You can install all required libraries and rvt-py with Anaconda environment from Anaconda cloud ([conda rvt_py](https://anaconda.org/zmigyyy/rvt_py)). To do that open Anaconda Prompt (activate conda environment) and run:

`conda install -c zmigyyy rvt_py`

### pypi

Another option is to install required libraries and rvt-py with Python Package Index, pypi ([pypi rvt-py](https://pypi.org/project/rvt-py)). To do that open command prompt (terminal) and run:

`pip install rvt-py`

This might not work, because pypi usually has problems installing gdal. To solve that first try to install gdal then run above command.

### requirements

We suggest using an Anaconda environment (easiest gdal installation) and Python 3.6 or higher. Required libraries with tested versions (could also work with other versions):

*   numpy 1.19.2
*   scipy 1.5.2
*   gdal 3.0.2


You can also clone the repository ([github rvt_py](https://github.com/EarthObservation/RVT_py)).

Library `rvt-py` can be used in Python scripts, Jupyter Notebooks and in ArcGIS Pro.

## Documentation

Documentation of the package and its usage is available at [Read the Docs](https://rvt-py.readthedocs.io/).

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please report any bugs and suggestions for improvements.

## Acknowledgment

Development of RVT Python scripts was part financed by the Slovenian Research Agency core funding No. P2-0406, and by research project No. J6-9395.

## License
This project is licensed under the terms of the [Apache License](LICENSE).

