"""
Every visualization config satisfies the shared :class:`rvt.core.Visualization`
contract — the base class enforces the methods, these tests enforce the parts
it cannot: that the outputs are well-formed and match the compute result.
"""

from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution, Visualization, VisualizationOutput
from rvt.stretch import Stretch
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

DEM_NAME = "dem"

# Every config the session can run, with the optional sky-view members on so
# their conditional outputs are covered too.
ALL_CONFIGS: list[Visualization] = [
    Hillshade(),
    LocalDominance(),
    MultiHillshade(),
    MultiScaleReliefModel(feature_minimum=0.0, feature_maximum=20.0, scaling_factor=2),
    MultiScaleTopographicPosition(),
    NegativeOpenness(),
    ShadowHorizon(),
    SimpleLocalReliefModel(),
    SkyIllumination(),
    SkyViewFactor(include_anisotropic=True, include_positive_openness=True),
    SlopeAspect(),
]

# Members with no legacy 8-bit form, plus the not-yet-implemented 3-band RGB
# multi-hillshade rendering.
FLOAT_ONLY = {
    ("SlopeAspect", "aspect"),
    ("ShadowHorizon", "horizon"),
    ("MultiHillshade", None),
}


def _identify(config: Visualization) -> str:
    return type(config).__name__


@pytest.mark.parametrize("config", ALL_CONFIGS, ids=_identify)
def test_config_exposes_the_contract(config: Visualization) -> None:
    assert isinstance(config, Visualization)
    assert isinstance(config.TILEABLE, bool)
    assert isinstance(config.vertical_exaggeration, float)
    assert isinstance(config.required_overlap(Resolution(1.0, 1.0)), int)


def test_subclass_must_declare_tileable() -> None:
    with pytest.raises(TypeError, match="TILEABLE"):

        class _NoTileable(Visualization[np.ndarray]):
            vertical_exaggeration: float = 1.0


def test_incomplete_subclass_cannot_be_instantiated() -> None:
    class _Incomplete(Visualization[np.ndarray]):
        TILEABLE = True
        vertical_exaggeration: float = 1.0

    with pytest.raises(TypeError, match="abstract"):
        _Incomplete()  # pyrefly: ignore  # deliberately abstract


@pytest.mark.parametrize("config", ALL_CONFIGS, ids=_identify)
def test_outputs_are_well_formed(config: Visualization) -> None:
    outputs = config.outputs(DEM_NAME)

    assert isinstance(outputs, tuple)
    assert outputs, "a config must produce at least one output"
    for output in outputs:
        assert isinstance(output, VisualizationOutput)
        assert output.stem.startswith(f"{DEM_NAME}_") or output.stem == DEM_NAME
        assert output.member is None or isinstance(output.member, str)
        if (_identify(config), output.member) in FLOAT_ONLY:
            assert output.eight_bit_stretch is None
        else:
            assert isinstance(output.eight_bit_stretch, Stretch)


@pytest.mark.parametrize("config", ALL_CONFIGS, ids=_identify)
def test_output_stems_are_distinct(config: Visualization) -> None:
    outputs = config.outputs(DEM_NAME)

    assert len({output.stem for output in outputs}) == len(outputs)


@pytest.mark.parametrize("config", ALL_CONFIGS, ids=_identify)
def test_output_members_match_the_compute_result(config: Visualization) -> None:
    """
    A single-output config returns the array itself; a multi-output one
    returns a result object carrying exactly the named members."""
    dem = np.zeros((24, 24), dtype=np.float32)

    result = config.compute(dem, Resolution(1.0, 1.0))

    members = [output.member for output in config.outputs(DEM_NAME)]
    if members == [None]:
        assert isinstance(result, np.ndarray)
        return
    for member in members:
        assert member is not None, "a multi-output config must name every member"
        assert getattr(result, member) is not None


def test_optional_sky_view_members_are_omitted_when_not_requested() -> None:
    members = [output.member for output in SkyViewFactor().outputs(DEM_NAME)]

    assert members == ["sky_view_factor"]
