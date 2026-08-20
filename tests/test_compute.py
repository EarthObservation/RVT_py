from __future__ import annotations

import numpy as np

import rvt
from rvt.io import read
from rvt.visualization import Hillshade, MultiScaleTopographicPosition
from tests.conftest import GOLDEN_DIRECTORY, read_golden_array


def test_compute_wraps_result_with_georeferencing() -> None:
    dem = read(GOLDEN_DIRECTORY / "inputs" / "TM1_564_146_crop.tif")

    hillshade = rvt.compute(Hillshade(), dem)

    np.testing.assert_allclose(
        hillshade.array.astype(np.float32),
        read_golden_array("vis/hillshade/defaults__float.tif"),
        rtol=1e-7,
    )
    assert hillshade.transform == dem.transform
    assert hillshade.crs == dem.crs
    assert hillshade.resolution == dem.resolution


def test_compute_supports_multiband_results() -> None:
    dem = read(GOLDEN_DIRECTORY / "inputs" / "TM1_564_146_crop.tif")

    mstp = rvt.compute(MultiScaleTopographicPosition(), dem)

    assert mstp.band_count == 3
    assert mstp.transform == dem.transform
