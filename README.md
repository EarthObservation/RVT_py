[![PyPI](https://img.shields.io/pypi/v/RVT_py?style=flat-square)](https://pypi.org/project/rvt-py/)

# Relief Visualization Toolbox Python library

Relief Visualization Toolbox (RVT) helps scientists visualize raster elevation model datasets. The selection of techniques is narrowed down to those that have proven effective for identification of small-scale features. The default settings assume high-resolution digital elevation models derived from airborne laser scanning (lidar), but RVT methods can be used for other purposes as well.

Version 3 is a ground-up rewrite: a typed, rasterio-based API verified against the previous implementation with a golden-master test suite. The legacy API (`rvt.vis`, `rvt.blend`, `rvt.default`) is gone; version 2 documentation remains valid for version 2 releases.

Implemented visualizations:

* hillshading, and hillshading from multiple directions,
* slope gradient and aspect,
* simple local relief model,
* multi-scale relief model (MSRM),
* multi-scale topographic position (MSTP),
* sky-view factor and anisotropic sky-view factor,
* positive and negative openness,
* sky illumination and shadow horizon,
* local dominance,

and blends of them, including the shipped presets: visualization for archaeological topography (VAT), color relief image map (CRIM), e3MSTP and e4MSTP.

## Quick start

Every visualization is a frozen config dataclass with a `compute` method over a float32 DEM array (NaN as nodata). Raster I/O is rasterio-backed; `rvt.compute` bridges the two:

```python
import rvt
import rvt.io
from rvt.visualization import Hillshade, SkyViewFactor

dem = rvt.io.read("dem.tif")
hillshade = rvt.compute(Hillshade(sun_azimuth=315, sun_elevation=35), dem)
rvt.io.write(hillshade, "dem_hillshade.tif")

result = SkyViewFactor(include_positive_openness=True).compute(dem.array, dem.resolution)
```

Batch processing with parameter-encoded output names replaces the old `DefaultValues.save_*` workflow:

```python
import rvt
from rvt.visualization import Hillshade, SimpleLocalReliefModel, SkyViewFactor

session = rvt.Session(
    dem_path="dem.tif",
    tasks=[
        rvt.SessionTask(Hillshade(), save_8bit=True),
        rvt.SessionTask(SimpleLocalReliefModel()),
        rvt.SessionTask(SkyViewFactor(include_positive_openness=True)),
    ],
)
written_paths = session.run()  # dem_HS_A315_H35.tif, dem_SLRM_R20.tif, ...
```

Blends are typed layer stacks; the shipped presets are factory functions:

```python
import rvt.io
from rvt.blending import presets

dem = rvt.io.read("dem.tif")
vat = presets.vat()
image = vat.render(dem.array, dem.resolution)
```

## Development

The project uses [uv](https://docs.astral.sh/uv/). Common tasks:

```bash
uv sync
uv run pytest -n auto
uv run ruff check .
uv run ruff format .
uv run pyrefly check
```

The test suite verifies the implementation against golden-master fixtures in `tests/golden/`, generated with the legacy (version 2) code; see `tests/golden/manifest.json` for exact parameters and provenance. The legacy implementation is temporarily kept under `legacy/` until the rewrite fully replaces it.

## References

When using the tools, please cite:

* Kokalj, Ž., Somrak, M. 2019. Why Not a Single Image? Combining Visualizations to Facilitate Fieldwork and On-Screen Mapping. Remote Sensing 11(7): 747.
* Zakšek, K., Oštir, K., Kokalj, Ž. 2011. Sky-View Factor as a Relief Visualization Technique. Remote Sensing 3: 398-415.

## License

See [LICENSE](LICENSE).

## Acknowledgment

Development of RVT was part financed by the European Commission's Culture Programme through the ArchaeoLandscapes Europe project and by the Slovenian Research Agency core funding No. P2-0406, and projects No. J6-7085 and No. J6-9395.

© Copyright 2010-2026 ZRC SAZU and University of Ljubljana
