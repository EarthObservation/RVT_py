from __future__ import annotations

import json
from collections.abc import Callable

import numpy as np
import pytest

from rvt.blending import (
    BlendCombination,
    BlendedSource,
    BlendLayer,
    BlendMode,
    ColormappedSource,
    ColorRamp,
    ImageSource,
    MultiDirectionalHillshadeRgbSource,
    NegativeOpennessSource,
    Normalization,
    NormalizationMode,
    OpennessDifferenceSource,
    PositiveOpennessSource,
    SlopeGradientSource,
    VisualizationSource,
    combination_from_json,
    combination_to_json,
    presets,
)
from rvt.visualization import (
    LocalDominance,
    MultiScaleTopographicPosition,
    NegativeOpenness,
    NoiseRemoval,
    SkyViewFactor,
    SlopeAspect,
    SlopeUnits,
)

VAT_EXPECTED_JSON = """\
{
  "layers": [
    {
      "blend_mode": "multiply",
      "normalization": {
        "maximum": 1.0,
        "minimum": 0.7,
        "mode": "value"
      },
      "opacity": 25.0,
      "source": {
        "config": {
          "anisotropy_direction": 315.0,
          "anisotropy_level": 1,
          "directions": 16,
          "include_anisotropic": false,
          "include_positive_openness": false,
          "max_radius": 10,
          "noise_removal": "none",
          "type": "SkyViewFactor",
          "vertical_exaggeration": 1.0
        },
        "inverts_normalization": false,
        "kind": "sky_view_factor"
      }
    },
    {
      "blend_mode": "overlay",
      "normalization": {
        "maximum": 93.0,
        "minimum": 68.0,
        "mode": "value"
      },
      "opacity": 50.0,
      "source": {
        "config": {
          "anisotropy_direction": 315.0,
          "anisotropy_level": 1,
          "directions": 16,
          "include_anisotropic": false,
          "include_positive_openness": true,
          "max_radius": 10,
          "noise_removal": "none",
          "type": "SkyViewFactor",
          "vertical_exaggeration": 1.0
        },
        "inverts_normalization": false,
        "kind": "positive_openness"
      }
    },
    {
      "blend_mode": "luminosity",
      "normalization": {
        "maximum": 50.0,
        "minimum": 0.0,
        "mode": "value"
      },
      "opacity": 50.0,
      "source": {
        "config": {
          "type": "SlopeAspect",
          "units": "degree",
          "vertical_exaggeration": 1.0
        },
        "inverts_normalization": true,
        "kind": "slope_gradient"
      }
    },
    {
      "blend_mode": "normal",
      "normalization": {
        "maximum": 1.0,
        "minimum": 0.0,
        "mode": "value"
      },
      "opacity": 100.0,
      "source": {
        "config": {
          "sun_azimuth": 315.0,
          "sun_elevation": 35.0,
          "type": "Hillshade",
          "vertical_exaggeration": 1.0
        },
        "inverts_normalization": false,
        "kind": "visualization"
      }
    }
  ],
  "name": "VAT - Archaeological",
  "schema": "rvt-blend/1"
}"""


def _every_variant_combination() -> BlendCombination:
    """A hand-built combination exercising every source kind, blend mode, and normalization mode."""
    nested = BlendCombination(
        name="nested",
        layers=(
            BlendLayer(
                source=MultiDirectionalHillshadeRgbSource(
                    sun_elevation=40.0, vertical_exaggeration=2.0
                ),
                normalization=Normalization(NormalizationMode.PERCENT, 1.0, 2.0),
                blend_mode=BlendMode.SCREEN,
                opacity=30.0,
            ),
            BlendLayer(
                source=NegativeOpennessSource(
                    config=NegativeOpenness(
                        directions=8, noise_removal=NoiseRemoval.MEDIUM
                    ),
                    inverts_normalization=False,
                ),
                normalization=Normalization(NormalizationMode.VALUE, 60.0, 95.0),
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )
    return BlendCombination(
        name="every variant",
        layers=(
            BlendLayer(
                source=VisualizationSource(
                    MultiScaleTopographicPosition(local_scale=(2, 8, 2), lightness=0.9)
                ),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
                blend_mode=BlendMode.OVERLAY,
                opacity=90.0,
            ),
            BlendLayer(
                source=BlendedSource(nested, inverts_normalization=True),
                normalization=Normalization(NormalizationMode.VALUE, -0.5, 0.5),
                blend_mode=BlendMode.SOFT_LIGHT,
                opacity=70.0,
            ),
            BlendLayer(
                source=OpennessDifferenceSource(
                    positive_config=SkyViewFactor(
                        max_radius=20, include_positive_openness=True
                    ),
                    negative_config=NegativeOpenness(max_radius=20),
                ),
                normalization=Normalization(NormalizationMode.PERCENT, 2.0, 3.0),
                blend_mode=BlendMode.LUMINOSITY,
                opacity=50.0,
            ),
            BlendLayer(
                source=PositiveOpennessSource(
                    config=SkyViewFactor(include_positive_openness=True, directions=32)
                ),
                normalization=Normalization(NormalizationMode.VALUE, 60.0, 95.0),
                blend_mode=BlendMode.MULTIPLY,
                opacity=25.0,
            ),
            BlendLayer(
                source=ColormappedSource(
                    source=SlopeGradientSource(
                        config=SlopeAspect(units=SlopeUnits.PERCENT),
                        inverts_normalization=False,
                    ),
                    normalization=Normalization(NormalizationMode.VALUE, 0.0, 0.8),
                    color_ramp=ColorRamp("viridis"),
                ),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
            BlendLayer(
                source=VisualizationSource(LocalDominance(minimum_radius=5)),
                normalization=Normalization(NormalizationMode.VALUE, 0.5, 1.8),
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )


@pytest.mark.parametrize(
    "factory", [presets.vat, presets.crim, presets.e3mstp, presets.e4mstp]
)
def test_preset_round_trips(factory: Callable[[], BlendCombination]) -> None:
    combination = factory()

    assert combination_from_json(combination_to_json(combination)) == combination


def test_every_variant_combination_round_trips() -> None:
    combination = _every_variant_combination()

    assert combination_from_json(combination_to_json(combination)) == combination


def test_blend_combination_methods_delegate_to_module_functions() -> None:
    combination = presets.vat()

    text = combination.to_json()

    assert text == combination_to_json(combination)
    assert BlendCombination.from_json(text) == combination


def test_vat_json_is_stable() -> None:
    assert combination_to_json(presets.vat()) == VAT_EXPECTED_JSON


def test_image_source_is_rejected() -> None:
    combination = BlendCombination(
        name="with image",
        layers=(
            BlendLayer(
                source=ImageSource(np.zeros((2, 2), dtype=np.float32)),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
            ),
        ),
    )

    with pytest.raises(ValueError, match=r"layer 0 of combination 'with image'"):
        combination_to_json(combination)


def test_colormapped_image_source_is_rejected() -> None:
    combination = BlendCombination(
        name="with colormapped image",
        layers=(
            BlendLayer(
                source=ColormappedSource(
                    source=ImageSource(np.zeros((2, 2), dtype=np.float32)),
                    normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
                    color_ramp=ColorRamp("viridis"),
                ),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
            ),
        ),
    )

    with pytest.raises(ValueError, match="ImageSource"):
        combination_to_json(combination)


def test_nested_image_source_is_rejected() -> None:
    inner = BlendCombination(
        name="inner",
        layers=(
            BlendLayer(
                source=ImageSource(np.zeros((2, 2), dtype=np.float32)),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
            ),
        ),
    )
    combination = BlendCombination(
        name="outer",
        layers=(
            BlendLayer(
                source=BlendedSource(inner),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
            ),
        ),
    )

    with pytest.raises(ValueError, match=r"layer 0 of combination 'inner'"):
        combination_to_json(combination)


def _vat_document() -> dict:
    return json.loads(combination_to_json(presets.vat()))


def test_unknown_schema_version_is_rejected() -> None:
    document = _vat_document()
    document["schema"] = "rvt-blend/99"

    with pytest.raises(ValueError, match=r"'rvt-blend/99'"):
        combination_from_json(json.dumps(document))


def test_missing_schema_is_rejected() -> None:
    document = _vat_document()
    del document["schema"]

    with pytest.raises(ValueError, match="schema"):
        combination_from_json(json.dumps(document))


def test_unknown_source_kind_is_rejected() -> None:
    document = _vat_document()
    document["layers"][0]["source"]["kind"] = "hologram"

    with pytest.raises(ValueError, match=r"'hologram'"):
        combination_from_json(json.dumps(document))


def test_unknown_top_level_field_is_rejected() -> None:
    document = _vat_document()
    document["legacy_terrain"] = "flat"

    with pytest.raises(ValueError, match=r"'legacy_terrain'"):
        combination_from_json(json.dumps(document))


def test_unknown_layer_field_is_rejected() -> None:
    document = _vat_document()
    document["layers"][0]["visibility"] = True

    with pytest.raises(ValueError, match=r"'visibility'"):
        combination_from_json(json.dumps(document))


def test_unknown_config_field_is_rejected() -> None:
    document = _vat_document()
    document["layers"][0]["source"]["config"]["sun_azimuth"] = 315.0

    with pytest.raises(ValueError, match=r"'sun_azimuth'"):
        combination_from_json(json.dumps(document))


def test_unknown_config_type_is_rejected() -> None:
    document = _vat_document()
    document["layers"][3]["source"]["config"] = {"type": "Hologram"}

    with pytest.raises(ValueError, match=r"'Hologram'"):
        combination_from_json(json.dumps(document))


def test_invalid_blend_mode_is_rejected() -> None:
    document = _vat_document()
    document["layers"][0]["blend_mode"] = "sparkle"

    with pytest.raises(ValueError, match=r"'sparkle'"):
        combination_from_json(json.dumps(document))


def test_invalid_config_enum_value_is_rejected() -> None:
    document = _vat_document()
    document["layers"][0]["source"]["config"]["noise_removal"] = "extreme"

    with pytest.raises(ValueError, match=r"'extreme'"):
        combination_from_json(json.dumps(document))


def test_mismatched_typed_config_is_rejected() -> None:
    document = _vat_document()
    document["layers"][2]["source"]["config"]["type"] = "Hillshade"

    with pytest.raises(ValueError, match=r"'Hillshade'"):
        combination_from_json(json.dumps(document))
