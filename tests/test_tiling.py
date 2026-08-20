from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
import rasterio
from affine import Affine
from rasterio.crs import CRS

import rvt
from rvt.core import Resolution
from rvt.io import read
from rvt.session import Session, SessionTask
from rvt.tiling import compute_tiled
from rvt.visualization import (
    Hillshade,
    LocalDominance,
    MultiHillshade,
    MultiScaleReliefModel,
    MultiScaleTopographicPosition,
    NegativeOpenness,
    ShadowHorizon,
    SimpleLocalReliefModel,
    SkyIllumination,
    SkyViewFactor,
    SlopeAspect,
)
from tests.conftest import GOLDEN_DIRECTORY

DEM_PATH = GOLDEN_DIRECTORY / "inputs" / "TM1_564_146_crop.tif"
NODATA_DEM_PATH = GOLDEN_DIRECTORY / "inputs" / "TM1_564_146_nodata_crop.tif"

MSRM_CONFIG = MultiScaleReliefModel(
    feature_minimum=1.0, feature_maximum=20.0, scaling_factor=2
)


class _Output:
    def __init__(self, path: Path) -> None:
        with rasterio.open(path) as source:
            array = source.read()
            self.array: np.ndarray = array[0] if array.shape[0] == 1 else array
            self.transform: Affine = source.transform
            self.crs: CRS | None = source.crs
            self.nodata: float | None = source.nodata


@pytest.mark.parametrize("tile_size", [(200, 200), (190, 190)])
@pytest.mark.parametrize(
    "config",
    [
        Hillshade(),
        SimpleLocalReliefModel(),
        NegativeOpenness(),
        LocalDominance(),
        MSRM_CONFIG,
    ],
    ids=lambda config: type(config).__name__,
)
def test_tiled_equals_untiled(
    config: rvt.SingleImageVisualization,
    tile_size: tuple[int, int],
    tmp_path: Path,
) -> None:
    untiled = rvt.compute(config, read(DEM_PATH))
    output_path = tmp_path / "tiled.tif"

    compute_tiled(config, DEM_PATH, output_path, tile_size=tile_size)

    output = _Output(output_path)
    np.testing.assert_array_equal(output.array, untiled.array)
    assert output.transform == untiled.transform
    assert output.crs == untiled.crs
    assert output.nodata is not None and np.isnan(output.nodata)


def test_mstp_overlap_exceeding_raster_equals_untiled(tmp_path: Path) -> None:
    # The default broad-scale overlap (501 px) exceeds the 500 px test DEM, so
    # every expanded tile read clamps to the whole raster; the multi-band
    # output must still equal the untiled compute exactly.
    config = MultiScaleTopographicPosition()
    untiled = rvt.compute(config, read(DEM_PATH))
    output_path = tmp_path / "tiled.tif"

    compute_tiled(config, DEM_PATH, output_path, tile_size=(190, 190))

    output = _Output(output_path)
    assert output.array.shape == (3, 500, 500)
    np.testing.assert_array_equal(output.array, untiled.array)


@pytest.mark.parametrize(
    "config",
    [Hillshade(), SimpleLocalReliefModel()],
    ids=lambda config: type(config).__name__,
)
def test_tiled_equals_untiled_with_nodata(
    config: rvt.SingleImageVisualization, tmp_path: Path
) -> None:
    untiled = rvt.compute(config, read(NODATA_DEM_PATH))
    output_path = tmp_path / "tiled.tif"

    compute_tiled(config, NODATA_DEM_PATH, output_path, tile_size=(190, 190))

    np.testing.assert_array_equal(_Output(output_path).array, untiled.array)


def test_workers_produce_identical_output(tmp_path: Path) -> None:
    serial_path = tmp_path / "serial.tif"
    parallel_path = tmp_path / "parallel.tif"

    compute_tiled(Hillshade(), DEM_PATH, serial_path, tile_size=(190, 190))
    compute_tiled(
        Hillshade(), DEM_PATH, parallel_path, tile_size=(190, 190), workers=2
    )

    np.testing.assert_array_equal(
        _Output(parallel_path).array, _Output(serial_path).array
    )


@pytest.mark.parametrize(
    "config",
    [SkyIllumination(), MultiHillshade()],
    ids=lambda config: type(config).__name__,
)
def test_non_tileable_config_raises(config: object, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="cannot be computed tiled"):
        compute_tiled(cast(Any, config), DEM_PATH, tmp_path / "out.tif")


@pytest.mark.parametrize(
    "config",
    [SlopeAspect(), SkyViewFactor(), ShadowHorizon()],
    ids=lambda config: type(config).__name__,
)
def test_multi_output_config_raises(config: object, tmp_path: Path) -> None:
    with pytest.raises(TypeError, match="multiple outputs"):
        compute_tiled(cast(Any, config), DEM_PATH, tmp_path / "out.tif")


def test_invalid_tile_size_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="tile_size"):
        compute_tiled(Hillshade(), DEM_PATH, tmp_path / "out.tif", tile_size=(0, 200))


def test_invalid_workers_raises(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="workers"):
        compute_tiled(Hillshade(), DEM_PATH, tmp_path / "out.tif", workers=0)


ONE_METRE = Resolution(1.0, 1.0)


@pytest.mark.parametrize(
    "config, resolution, expected",
    [
        (Hillshade(), ONE_METRE, 1),
        (SlopeAspect(), ONE_METRE, 1),
        (MultiHillshade(), ONE_METRE, 1),
        (SimpleLocalReliefModel(radius_cells=30), ONE_METRE, 31),
        (SkyViewFactor(max_radius=15), ONE_METRE, 15),
        (NegativeOpenness(), ONE_METRE, 10),
        (LocalDominance(maximum_radius=25), ONE_METRE, 25),
        (ShadowHorizon(), ONE_METRE, 100),
        (SkyIllumination(maximum_fine_radius=250), ONE_METRE, 250),
        (MultiScaleTopographicPosition(), ONE_METRE, 501),
        # ceil(((20 - 1) / 2) ** (1 / 2)) ** 2 + 1 = 17: the kernel radius in
        # pixels, not the legacy meters-as-pixels value (20).
        (MSRM_CONFIG, ONE_METRE, 17),
        # At half-metre cells the honest pixel radius exceeds the legacy value:
        # ceil(((20 - 0.5) / 1) ** (1 / 2)) ** 2 + 1 = 26.
        (MSRM_CONFIG, Resolution(0.5, 0.5), 26),
    ],
    ids=lambda value: type(value).__name__ if hasattr(value, "compute") else None,
)
def test_required_overlap(
    config: Any, resolution: Resolution, expected: int
) -> None:
    assert config.required_overlap(resolution) == expected


def _run_session(
    dem_path: Path, output_directory: Path, **session_kwargs: Any
) -> list[Path]:
    session = Session(
        dem_path,
        tasks=[
            SessionTask(Hillshade(), save_8bit=True),
            SessionTask(SimpleLocalReliefModel()),
            SessionTask(SlopeAspect()),
        ],
        output_directory=output_directory,
        **session_kwargs,
    )
    return session.run()


def test_session_switches_to_tiled_processing(tmp_path: Path) -> None:
    untiled_paths = _run_session(DEM_PATH, tmp_path / "untiled")
    tiled_paths = _run_session(
        DEM_PATH,
        tmp_path / "tiled",
        automatic_tiling_threshold=10_000,
        tile_size=(190, 190),
    )

    assert [path.name for path in tiled_paths] == [
        path.name for path in untiled_paths
    ]
    for tiled_path, untiled_path in zip(tiled_paths, untiled_paths, strict=True):
        tiled, untiled = _Output(tiled_path), _Output(untiled_path)
        np.testing.assert_array_equal(tiled.array, untiled.array)
        assert tiled.transform == untiled.transform
        assert tiled.crs == untiled.crs
