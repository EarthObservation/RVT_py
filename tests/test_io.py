from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio

from rvt.io import Raster, read, write
from tests.conftest import GOLDEN_DIRECTORY

REAL_CROP = GOLDEN_DIRECTORY / "inputs" / "TM1_564_146_crop.tif"
NODATA_CROP = GOLDEN_DIRECTORY / "inputs" / "TM1_564_146_nodata_crop.tif"


def test_read_returns_float32_with_derived_resolution() -> None:
    raster = read(REAL_CROP)

    assert raster.array.dtype == np.float32
    assert raster.array.shape == (500, 500)
    assert raster.resolution.x == 1.0
    assert raster.resolution.y == 1.0
    assert raster.crs is not None


def test_read_converts_nodata_to_nan() -> None:
    raster = read(NODATA_CROP)

    with rasterio.open(NODATA_CROP) as source:
        raw = source.read(1)
        nodata = source.nodata
        assert nodata is not None
        expected_nan_count = (
            int(np.isnan(raw).sum()) if np.isnan(nodata) else int((raw == nodata).sum())
        )

    assert int(np.isnan(raster.array).sum()) == expected_nan_count
    assert expected_nan_count > 0


def test_write_read_round_trip(tmp_path: Path) -> None:
    raster = read(NODATA_CROP)
    output_path = tmp_path / "round_trip.tif"

    write(raster, output_path)
    round_tripped = read(output_path)

    np.testing.assert_array_equal(round_tripped.array, raster.array)
    assert round_tripped.transform == raster.transform
    assert round_tripped.crs == raster.crs
    with rasterio.open(output_path) as source:
        assert source.nodata is not None and np.isnan(source.nodata)
        assert source.profile["compress"] == "lzw"


def test_write_uint8_has_no_nodata_tag(tmp_path: Path) -> None:
    template = read(REAL_CROP)
    eight_bit = Raster(
        array=np.zeros((500, 500), dtype=np.uint8),
        transform=template.transform,
        crs=template.crs,
    )
    output_path = tmp_path / "eight_bit.tif"

    write(eight_bit, output_path)

    with rasterio.open(output_path) as source:
        assert source.nodata is None
        assert source.dtypes[0] == "uint8"


def test_write_multiband(tmp_path: Path) -> None:
    template = read(REAL_CROP)
    three_band = Raster(
        array=np.random.default_rng(0).random((3, 500, 500)).astype(np.float32),
        transform=template.transform,
        crs=template.crs,
    )
    output_path = tmp_path / "three_band.tif"

    write(three_band, output_path)

    with rasterio.open(output_path) as source:
        assert source.count == 3
        np.testing.assert_array_equal(source.read(), three_band.array)
