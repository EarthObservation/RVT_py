"""
JSON serialization of blend combinations (schema ``rvt-blend/1``).

The schema is new and intentionally simple: enums are serialized as their
string values, and layers are flat objects referencing visualization configs
by class name plus their dataclass fields. There is no reader for the legacy
``blender_*.json`` format.
"""

from __future__ import annotations

import json
from dataclasses import fields, is_dataclass
from enum import Enum, StrEnum, auto
from typing import Any, get_origin, get_type_hints

from rvt.blending.blend_modes import BlendMode
from rvt.blending.colormap import ColorRamp
from rvt.blending.combination import (
    BlendCombination,
    BlendedSource,
    BlendLayer,
    ColormappedSource,
    ImageSource,
    LayerSource,
    MultiDirectionalHillshadeRgbSource,
    NegativeOpennessSource,
    OpennessDifferenceSource,
    PositiveOpennessSource,
    SkyViewFactorSource,
    SlopeGradientSource,
    VisualizationSource,
)
from rvt.blending.normalization import Normalization, NormalizationMode
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

SCHEMA = "rvt-blend/1"

_CONFIG_TYPES: dict[str, type] = {
    config_type.__name__: config_type
    for config_type in (
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
}


class _SourceKind(StrEnum):
    """Stable schema identifier for each serializable :data:`LayerSource` type."""

    VISUALIZATION = auto()
    SLOPE_GRADIENT = auto()
    SKY_VIEW_FACTOR = auto()
    POSITIVE_OPENNESS = auto()
    NEGATIVE_OPENNESS = auto()
    OPENNESS_DIFFERENCE = auto()
    MULTI_DIRECTIONAL_HILLSHADE_RGB = auto()
    COLORMAPPED = auto()
    BLENDED = auto()


def combination_to_json(combination: BlendCombination) -> str:
    """
    Serialize a blend combination to a human-readable JSON document
    (schema ``rvt-blend/1``, sorted keys, two-space indentation).

    :param combination: The combination to serialize.
    :returns: The JSON text.
    :raises ValueError: If any layer (also in nested combinations) uses
        :class:`ImageSource`, whose raw array cannot be serialized.
    """
    payload = {"schema": SCHEMA, **_combination_to_dict(combination)}
    return json.dumps(payload, indent=2, sort_keys=True)


def combination_from_json(text: str) -> BlendCombination:
    """
    Deserialize a blend combination from its ``rvt-blend/1`` JSON document.

    :param text: JSON text produced by :func:`combination_to_json`.
    :returns: The combination, equal to the one serialized.
    :raises ValueError: If the schema version, a layer source kind, an enum
        value, or any field is unknown, or the document is malformed.
    """
    data = json.loads(text)
    if not isinstance(data, dict):
        raise ValueError("The top-level JSON value must be an object!")
    schema = data.get("schema")
    if schema != SCHEMA:
        raise ValueError(f"Unknown blend schema '{schema}', expected '{SCHEMA}'!")
    _check_keys(data, {"schema", "name", "layers"}, "the top-level object")
    return _combination_from_dict(data, "the top-level object")


def _combination_to_dict(combination: BlendCombination) -> dict[str, Any]:
    return {
        "name": combination.name,
        "layers": [
            _layer_to_dict(layer, f"layer {index} of combination '{combination.name}'")
            for index, layer in enumerate(combination.layers)
        ],
    }


def _combination_from_dict(data: dict[str, Any], context: str) -> BlendCombination:
    layers = _required(data, "layers", context)
    if not isinstance(layers, list):
        raise ValueError(f"'layers' must be a list in {context}!")
    return BlendCombination(
        name=_required(data, "name", context),
        layers=tuple(
            _layer_from_dict(
                _object(layer, f"layer {index} in {context}"),
                f"layer {index} in {context}",
            )
            for index, layer in enumerate(layers)
        ),
    )


def _layer_to_dict(layer: BlendLayer, layer_label: str) -> dict[str, Any]:
    return {
        "source": _source_to_dict(layer.source, layer_label),
        "normalization": _normalization_to_dict(layer.normalization),
        "blend_mode": layer.blend_mode.value,
        "opacity": layer.opacity,
    }


def _layer_from_dict(data: dict[str, Any], context: str) -> BlendLayer:
    _check_keys(data, {"source", "normalization", "blend_mode", "opacity"}, context)
    return BlendLayer(
        source=_source_from_dict(
            _object(_required(data, "source", context), f"'source' of {context}"),
            context,
        ),
        normalization=_normalization_from_dict(
            _object(
                _required(data, "normalization", context),
                f"'normalization' of {context}",
            ),
            context,
        ),
        blend_mode=_parse_enum(
            BlendMode, _required(data, "blend_mode", context), context
        ),
        opacity=_required(data, "opacity", context),
    )


def _source_to_dict(source: LayerSource, layer_label: str) -> dict[str, Any]:
    inverts = {"inverts_normalization": source.inverts_normalization}
    match source:
        case ImageSource():
            raise ValueError(
                f"Cannot serialize {layer_label}: ImageSource holds a raw image"
                " array; use a visualization-config source instead!"
            )
        case VisualizationSource():
            return {
                "kind": _SourceKind.VISUALIZATION.value,
                "config": _config_to_dict(source.config, layer_label),
                **inverts,
            }
        case SlopeGradientSource():
            return {
                "kind": _SourceKind.SLOPE_GRADIENT.value,
                "config": _config_to_dict(source.config, layer_label),
                **inverts,
            }
        case SkyViewFactorSource():
            return {
                "kind": _SourceKind.SKY_VIEW_FACTOR.value,
                "config": _config_to_dict(source.config, layer_label),
                **inverts,
            }
        case PositiveOpennessSource():
            return {
                "kind": _SourceKind.POSITIVE_OPENNESS.value,
                "config": _config_to_dict(source.config, layer_label),
                **inverts,
            }
        case NegativeOpennessSource():
            return {
                "kind": _SourceKind.NEGATIVE_OPENNESS.value,
                "config": _config_to_dict(source.config, layer_label),
                **inverts,
            }
        case OpennessDifferenceSource():
            return {
                "kind": _SourceKind.OPENNESS_DIFFERENCE.value,
                "positive_config": _config_to_dict(source.positive_config, layer_label),
                "negative_config": _config_to_dict(source.negative_config, layer_label),
                **inverts,
            }
        case MultiDirectionalHillshadeRgbSource():
            return {
                "kind": _SourceKind.MULTI_DIRECTIONAL_HILLSHADE_RGB.value,
                "sun_elevation": source.sun_elevation,
                "vertical_exaggeration": source.vertical_exaggeration,
                **inverts,
            }
        case ColormappedSource():
            return {
                "kind": _SourceKind.COLORMAPPED.value,
                "source": _source_to_dict(source.source, layer_label),
                "normalization": _normalization_to_dict(source.normalization),
                "colormap": _color_ramp_to_dict(source.color_ramp),
                **inverts,
            }
        case BlendedSource():
            return {
                "kind": _SourceKind.BLENDED.value,
                "combination": _combination_to_dict(source.combination),
                **inverts,
            }


def _source_from_dict(data: dict[str, Any], context: str) -> LayerSource:
    kind_value = _required(data, "kind", context)
    try:
        kind = _SourceKind(kind_value)
    except ValueError as error:
        raise ValueError(
            f"Unknown layer source kind '{kind_value}' in {context}!"
        ) from error

    inverts = data.get("inverts_normalization", False)
    config_keys = {"kind", "config", "inverts_normalization"}

    match kind:
        case _SourceKind.VISUALIZATION:
            _check_keys(data, config_keys, context)
            return VisualizationSource(
                config=_config_from_dict(_required(data, "config", context), context),
                inverts_normalization=inverts,
            )
        case _SourceKind.SLOPE_GRADIENT:
            _check_keys(data, config_keys, context)
            return SlopeGradientSource(
                config=_typed_config_from_dict(data, "config", SlopeAspect, context),
                inverts_normalization=inverts,
            )
        case _SourceKind.SKY_VIEW_FACTOR:
            _check_keys(data, config_keys, context)
            return SkyViewFactorSource(
                config=_typed_config_from_dict(data, "config", SkyViewFactor, context),
                inverts_normalization=inverts,
            )
        case _SourceKind.POSITIVE_OPENNESS:
            _check_keys(data, config_keys, context)
            return PositiveOpennessSource(
                config=_typed_config_from_dict(data, "config", SkyViewFactor, context),
                inverts_normalization=inverts,
            )
        case _SourceKind.NEGATIVE_OPENNESS:
            _check_keys(data, config_keys, context)
            return NegativeOpennessSource(
                config=_typed_config_from_dict(
                    data, "config", NegativeOpenness, context
                ),
                inverts_normalization=inverts,
            )
        case _SourceKind.OPENNESS_DIFFERENCE:
            _check_keys(
                data,
                {"kind", "positive_config", "negative_config", "inverts_normalization"},
                context,
            )
            return OpennessDifferenceSource(
                positive_config=_typed_config_from_dict(
                    data, "positive_config", SkyViewFactor, context
                ),
                negative_config=_typed_config_from_dict(
                    data, "negative_config", NegativeOpenness, context
                ),
                inverts_normalization=inverts,
            )
        case _SourceKind.MULTI_DIRECTIONAL_HILLSHADE_RGB:
            _check_keys(
                data,
                {
                    "kind",
                    "sun_elevation",
                    "vertical_exaggeration",
                    "inverts_normalization",
                },
                context,
            )
            return MultiDirectionalHillshadeRgbSource(
                sun_elevation=_required(data, "sun_elevation", context),
                vertical_exaggeration=_required(data, "vertical_exaggeration", context),
                inverts_normalization=inverts,
            )
        case _SourceKind.COLORMAPPED:
            _check_keys(
                data,
                {
                    "kind",
                    "source",
                    "normalization",
                    "colormap",
                    "inverts_normalization",
                },
                context,
            )
            return ColormappedSource(
                source=_source_from_dict(
                    _object(
                        _required(data, "source", context), f"'source' of {context}"
                    ),
                    context,
                ),
                normalization=_normalization_from_dict(
                    _object(
                        _required(data, "normalization", context),
                        f"'normalization' of {context}",
                    ),
                    context,
                ),
                color_ramp=_color_ramp_from_dict(
                    _object(
                        _required(data, "colormap", context), f"'colormap' of {context}"
                    ),
                    context,
                ),
                inverts_normalization=inverts,
            )
        case _SourceKind.BLENDED:
            _check_keys(data, {"kind", "combination", "inverts_normalization"}, context)
            return BlendedSource(
                combination=_combination_from_dict(
                    _object(
                        _required(data, "combination", context),
                        f"'combination' of {context}",
                    ),
                    context,
                ),
                inverts_normalization=inverts,
            )


def _config_to_dict(config: object, layer_label: str) -> dict[str, Any]:
    type_name = type(config).__name__
    if (
        type_name not in _CONFIG_TYPES
        or not is_dataclass(config)
        or isinstance(config, type)
    ):
        raise ValueError(
            f"Cannot serialize {layer_label}: '{type_name}' is not a known"
            " visualization config!"
        )
    result: dict[str, Any] = {"type": type_name}
    for config_field in fields(config):
        result[config_field.name] = _plain_value(getattr(config, config_field.name))
    return result


def _config_from_dict(data: Any, context: str) -> Any:
    config_data = _object(data, f"'config' of {context}")
    type_name = _required(config_data, "type", context)
    config_type = _CONFIG_TYPES.get(type_name)
    if config_type is None:
        raise ValueError(
            f"Unknown visualization config type '{type_name}' in {context}!"
        )
    return _build_config(config_type, config_data, context)


def _typed_config_from_dict(
    data: dict[str, Any], key: str, config_type: type, context: str
) -> Any:
    if key not in data:
        return config_type()
    config_data = _object(data[key], f"'{key}' of {context}")
    type_name = config_data.get("type", config_type.__name__)
    if type_name != config_type.__name__:
        raise ValueError(
            f"'{key}' of {context} must be a '{config_type.__name__}' config,"
            f" got '{type_name}'!"
        )
    return _build_config(config_type, config_data, context)


def _build_config(config_type: type, config_data: dict[str, Any], context: str) -> Any:
    field_types = get_type_hints(config_type)
    field_names = {config_field.name for config_field in fields(config_type)}
    keyword_arguments: dict[str, Any] = {}
    for key, value in config_data.items():
        if key == "type":
            continue
        if key not in field_names:
            raise ValueError(
                f"Unknown field '{key}' for config '{config_type.__name__}'"
                f" in {context}!"
            )
        keyword_arguments[key] = _typed_value(field_types[key], value, context)
    return config_type(**keyword_arguments)


def _normalization_to_dict(normalization: Normalization) -> dict[str, Any]:
    return {
        "mode": normalization.mode.value,
        "minimum": normalization.minimum,
        "maximum": normalization.maximum,
    }


def _normalization_from_dict(data: dict[str, Any], context: str) -> Normalization:
    _check_keys(data, {"mode", "minimum", "maximum"}, f"'normalization' of {context}")
    return Normalization(
        mode=_parse_enum(NormalizationMode, _required(data, "mode", context), context),
        minimum=_required(data, "minimum", context),
        maximum=_required(data, "maximum", context),
    )


def _color_ramp_to_dict(color_ramp: ColorRamp) -> dict[str, Any]:
    return {
        "name": color_ramp.name,
        "minimum_cut": color_ramp.minimum_cut,
        "maximum_cut": color_ramp.maximum_cut,
    }


def _color_ramp_from_dict(data: dict[str, Any], context: str) -> ColorRamp:
    _check_keys(
        data, {"name", "minimum_cut", "maximum_cut"}, f"'colormap' of {context}"
    )
    return ColorRamp(
        name=_required(data, "name", context),
        minimum_cut=data.get("minimum_cut"),
        maximum_cut=data.get("maximum_cut"),
    )


def _plain_value(value: object) -> object:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return list(value)
    return value


def _typed_value(field_type: Any, value: Any, context: str) -> Any:
    if isinstance(field_type, type) and issubclass(field_type, Enum):
        return _parse_enum(field_type, value, context)
    if get_origin(field_type) is tuple and isinstance(value, list):
        return tuple(value)
    return value


def _parse_enum(enum_type: type[Any], value: Any, context: str) -> Any:
    try:
        return enum_type(value)
    except ValueError as error:
        raise ValueError(
            f"'{value}' is not a valid {enum_type.__name__} in {context}!"
        ) from error


def _object(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be a JSON object, got '{value}'!")
    return value


def _required(data: dict[str, Any], key: str, context: str) -> Any:
    if key not in data:
        raise ValueError(f"Missing field '{key}' in {context}!")
    return data[key]


def _check_keys(data: dict[str, Any], allowed: set[str], context: str) -> None:
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ValueError(f"Unknown field '{unknown[0]}' in {context}!")
