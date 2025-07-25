import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rvt_py",
    version="2.2.3",
    author="ZRC SAZU and University of Ljubljana (UL FGG)",
    author_email="ziga.kokalj@zrc-sazu.si",
    description="Relief Visualization Toolbox Python library. "
                "It helps scientist visualize raster elevation model datasets. ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/EarthObservation/RVT_py",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent"
    ],
    keywords="relief_visualization_toolbox relief_visualization relief rvt raster raster_visualization visualization",
    python_requires='>=3.6, <3.12',
    project_urls={
        'Documentation': 'https://rvt-py.readthedocs.io/en/latest/',
        'Source': 'https://github.com/EarthObservation/RVT_py',
        'Old RVT': 'https://iaps.zrc-sazu.si/en/rvt#v',
        'ArcGIS Pro': 'https://github.com/EarthObservation/rvt-arcgis-pro',
        "QGIS plugin": 'https://github.com/EarthObservation/rvt-qgis'''
    },
    install_requires=['numpy', 'scipy', 'gdal', 'matplotlib', 'pandas', 'geopandas', 'jupyter', 'rasterio']
)
