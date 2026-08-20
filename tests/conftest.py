from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import rasterio
from numpy.typing import NDArray

GOLDEN_DIRECTORY = Path(__file__).parent / "golden"


@pytest.fixture(scope="session")
def golden_manifest() -> dict[str, Any]:
    with (GOLDEN_DIRECTORY / "manifest.json").open() as manifest_file:
        return json.load(manifest_file)["fixtures"]


def read_golden_array(relative_path: str) -> NDArray[Any]:
    """Read a golden fixture raster as a bare array (2D, or (bands, rows, cols))."""
    with rasterio.open(GOLDEN_DIRECTORY / relative_path) as source:
        array = source.read()
    if array.shape[0] == 1:
        return array[0]
    return array


def read_input_dem(name: str) -> NDArray[np.float32]:
    """Read a golden input DEM as float32 with nodata converted to NaN."""
    with rasterio.open(GOLDEN_DIRECTORY / "inputs" / name) as source:
        array = source.read(1).astype(np.float32)
        if source.nodata is not None and not np.isnan(source.nodata):
            array[array == source.nodata] = np.nan
    return array
