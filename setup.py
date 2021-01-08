import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rvt_py",
    version="1.0.0a3",
    author="ZRC SAZU and University of Ljubljana",
    author_email="ziga.kokalj@zrc-sazu.si",
    description="Relief Visualization Toolbox python library. "
                "It helps scientist visualize raster elevation model datasets. ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/EarthObservation/RVT_py",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent"
    ],
    keywords="relief_visualization_toolbox relief_visualization relief rvt raster raster_visualization visualization",
    python_requires='>=3.6',
    project_urls={
        'Documentation': 'https://rvt-py.readthedocs.io/en/latest/',
        'Source': 'https://github.com/EarthObservation/RVT_py',
        'Old RVT': 'https://iaps.zrc-sazu.si/en/rvt#v',
        'ArcGIS Pro': 'https://github.com/EarthObservation/rvt-arcgis-pro',
        "QGIS plugin": 'https://github.com/EarthObservation/rvt-qgis'''
    },
    data_files=[
        ('settings', ['rvt/settings/blender_custom_layers.json', 'rvt/settings/blender_file_example.json',
                      'rvt/settings/default_blender_combinations.json', 'rvt/settings/default_settings.json',
                      'rvt/settings/default_terrains_settings.json'])
    ],
    install_requires=['numpy', 'scipy', 'gdal']
)
