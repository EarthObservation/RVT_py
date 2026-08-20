from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from rvt.io import read
from rvt.session import Session, SessionTask, VisualizationConfig
from rvt.visualization import (
    Hillshade,
    LocalDominance,
    MultiHillshade,
    MultiScaleReliefModel,
    MultiScaleTopographicPosition,
    NegativeOpenness,
    NoiseRemoval,
    ShadowHorizon,
    SimpleLocalReliefModel,
    SkyIllumination,
    SkyViewFactor,
    SlopeAspect,
    SlopeUnits,
)
from tests.conftest import GOLDEN_DIRECTORY, read_golden_array

DEM_PATH = GOLDEN_DIRECTORY / "inputs" / "TM1_564_146_crop.tif"
DEM_NAME = "TM1_564_146_crop"


def test_session_writes_param_encoded_names(tmp_path: Path) -> None:
    session = Session(
        DEM_PATH,
        tasks=[
            SessionTask(Hillshade()),
            SessionTask(SlopeAspect(units=SlopeUnits.DEGREE)),
            SessionTask(
                SkyViewFactor(include_anisotropic=True, include_positive_openness=True)
            ),
            SessionTask(SimpleLocalReliefModel()),
        ],
        output_directory=tmp_path,
    )

    written = session.run()

    expected_names = {
        f"{DEM_NAME}_HS_A315_H35.tif",
        f"{DEM_NAME}_SLOPE.tif",
        f"{DEM_NAME}_ASPECT.tif",
        f"{DEM_NAME}_SVF_R10_D16.tif",
        f"{DEM_NAME}_SVF-A_R10_D16_A315_ALlow.tif",
        f"{DEM_NAME}_OPEN-POS_R10_D16.tif",
        f"{DEM_NAME}_SLRM_R20.tif",
    }
    assert {path.name for path in written} == expected_names
    assert all(path.exists() for path in written)
    assert all(path.parent == tmp_path for path in written)


@pytest.mark.parametrize(
    "config, file_name",
    [
        (
            Hillshade(sun_azimuth=90.0, sun_elevation=10.0),
            f"{DEM_NAME}_HS_A90_H10.tif",
        ),
        (
            SkyViewFactor(max_radius=20, directions=8, noise_removal=NoiseRemoval.HIGH),
            f"{DEM_NAME}_SVF_R20_D8_NRhigh.tif",
        ),
        (
            NegativeOpenness(noise_removal=NoiseRemoval.HIGH),
            f"{DEM_NAME}_OPEN-NEG_R10_D16_NRhigh.tif",
        ),
        (MultiHillshade(), f"{DEM_NAME}_MULTI-HS_D16_H35.tif"),
        (
            LocalDominance(),
            f"{DEM_NAME}_LD_R_M10-20_DI1_A15_OH1.7.tif",
        ),
        (
            MultiScaleReliefModel(
                feature_minimum=0.0, feature_maximum=20.0, scaling_factor=2
            ),
            f"{DEM_NAME}_MSRM_F_M0-20_S2.tif",
        ),
        (
            MultiScaleTopographicPosition(),
            f"{DEM_NAME}_MSTP_5_50_500_L1.2.tif",
        ),
        (
            SkyIllumination(),
            f"{DEM_NAME}_SIM_overcast_D32_100px.tif",
        ),
    ],
)
def test_output_file_names_are_legacy_shaped(
    config: VisualizationConfig, file_name: str
) -> None:
    dem_path = Path("/data") / f"{DEM_NAME}.tif"
    session = Session(dem_path, tasks=[SessionTask(config)])

    assert {path.name for path in session.planned_paths()} >= {file_name}


def test_shadow_horizon_member_names() -> None:
    session = Session(DEM_PATH, tasks=[SessionTask(ShadowHorizon())])

    assert {path.name for path in session.planned_paths()} == {
        f"{DEM_NAME}_shadow_A315_H35.tif",
        f"{DEM_NAME}_horizon_A315_H35.tif",
    }


def test_float_outputs_round_trip_against_golden(tmp_path: Path) -> None:
    session = Session(
        DEM_PATH,
        tasks=[
            SessionTask(Hillshade()),
            SessionTask(SimpleLocalReliefModel()),
            SessionTask(SkyViewFactor()),
        ],
        output_directory=tmp_path,
    )

    session.run()

    for file_name, golden_relative_path in [
        (f"{DEM_NAME}_HS_A315_H35.tif", "vis/hillshade/defaults__float.tif"),
        (f"{DEM_NAME}_SLRM_R20.tif", "vis/slrm/defaults__float.tif"),
        (f"{DEM_NAME}_SVF_R10_D16.tif", "vis/sky_view_factor/defaults_svf__float.tif"),
    ]:
        raster = read(tmp_path / file_name)
        assert raster.array.dtype == np.float32
        np.testing.assert_allclose(
            raster.array, read_golden_array(golden_relative_path), rtol=1e-7
        )

    dem_raster = read(DEM_PATH)
    hillshade_raster = read(tmp_path / f"{DEM_NAME}_HS_A315_H35.tif")
    assert hillshade_raster.transform == dem_raster.transform
    assert hillshade_raster.crs == dem_raster.crs


def test_8bit_outputs_match_golden(tmp_path: Path) -> None:
    session = Session(
        DEM_PATH,
        tasks=[
            SessionTask(Hillshade(), save_float=False, save_8bit=True),
            SessionTask(
                SlopeAspect(units=SlopeUnits.DEGREE), save_float=False, save_8bit=True
            ),
            SessionTask(SimpleLocalReliefModel(), save_float=False, save_8bit=True),
        ],
        output_directory=tmp_path,
    )

    written = session.run()

    # Aspect has no legacy 8-bit rendering, so only three files are written.
    assert {path.name for path in written} == {
        f"{DEM_NAME}_HS_A315_H35_8bit.tif",
        f"{DEM_NAME}_SLOPE_8bit.tif",
        f"{DEM_NAME}_SLRM_R20_8bit.tif",
    }
    for file_name, golden_relative_path in [
        (f"{DEM_NAME}_HS_A315_H35_8bit.tif", "vis/hillshade/defaults__8bit.tif"),
        (f"{DEM_NAME}_SLOPE_8bit.tif", "vis/slope_aspect/defaults_slope__8bit.tif"),
        (f"{DEM_NAME}_SLRM_R20_8bit.tif", "vis/slrm/defaults__8bit.tif"),
    ]:
        raster = read(tmp_path / file_name)
        np.testing.assert_array_equal(
            raster.array.astype(np.uint8), read_golden_array(golden_relative_path)
        )


def test_skip_existing_honours_overwrite(tmp_path: Path) -> None:
    tasks = [SessionTask(Hillshade()), SessionTask(SimpleLocalReliefModel())]
    session = Session(DEM_PATH, tasks=tasks, output_directory=tmp_path)

    first_run = session.run()
    second_run = session.run()

    assert len(first_run) == 2
    assert second_run == []

    # A missing file is recreated without touching the existing one.
    (tmp_path / f"{DEM_NAME}_SLRM_R20.tif").unlink()
    third_run = session.run()
    assert [path.name for path in third_run] == [f"{DEM_NAME}_SLRM_R20.tif"]

    overwriting = Session(
        DEM_PATH, tasks=tasks, output_directory=tmp_path, overwrite=True
    )
    assert len(overwriting.run()) == 2


def test_per_task_name_override(tmp_path: Path) -> None:
    session = Session(
        DEM_PATH,
        tasks=[
            SessionTask(Hillshade(), name="custom_shading", save_8bit=True),
            SessionTask(SlopeAspect(), name="terrain"),
        ],
        output_directory=tmp_path,
    )

    written = session.run()

    assert {path.name for path in written} == {
        "custom_shading.tif",
        "custom_shading_8bit.tif",
        "terrain_slope.tif",
        "terrain_aspect.tif",
    }


def test_session_wide_vertical_exaggeration(tmp_path: Path) -> None:
    session = Session(
        DEM_PATH,
        tasks=[SessionTask(Hillshade(sun_azimuth=90.0, sun_elevation=10.0))],
        output_directory=tmp_path,
        vertical_exaggeration=2.0,
    )

    written = session.run()

    assert [path.name for path in written] == [f"{DEM_NAME}_HS_A90_H10.tif"]
    np.testing.assert_allclose(
        read(written[0]).array,
        read_golden_array("vis/hillshade/az90_el10_ve2__float.tif"),
        rtol=1e-7,
    )


def test_default_output_directory_is_dem_directory(tmp_path: Path) -> None:
    dem_copy = tmp_path / DEM_PATH.name
    shutil.copy(DEM_PATH, dem_copy)

    written = Session(dem_copy, tasks=[SessionTask(Hillshade())]).run()

    assert written == [tmp_path / f"{DEM_NAME}_HS_A315_H35.tif"]


def test_missing_output_directory_is_created(tmp_path: Path) -> None:
    output_directory = tmp_path / "outputs" / "nested"

    written = Session(
        DEM_PATH, tasks=[SessionTask(Hillshade())], output_directory=output_directory
    ).run()

    assert written == [output_directory / f"{DEM_NAME}_HS_A315_H35.tif"]


def test_multi_hillshade_8bit_is_not_implemented() -> None:
    with pytest.raises(NotImplementedError, match="[Mm]ulti-hillshade"):
        Session(DEM_PATH, tasks=[SessionTask(MultiHillshade(), save_8bit=True)])


def test_task_without_any_output_raises() -> None:
    with pytest.raises(ValueError, match="save_float"):
        SessionTask(Hillshade(), save_float=False, save_8bit=False)
