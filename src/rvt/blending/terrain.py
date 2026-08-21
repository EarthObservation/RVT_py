"""
Terrain settings: per-terrain-type parameter and normalization overrides that
can be applied to any blend combination.

The three shipped terrains ("general", "flat", "steep") carry the values of the
legacy ``settings/default_terrains_settings.json`` as code-level constants; the
settings JSON files are not read any more.

The legacy :class:`rvt.blend.TerrainSettings` assigned the "Sky illumination"
and the "Local dominance" stretch to the same field as the "Openness -
Negative" stretch, so those two stretches were parsed but never applied. Here
every visualization's normalization bounds go to its own field, which is the
intended behaviour and therefore a deliberate deviation from version 2.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import StrEnum, auto

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
    ShadowSource,
    SkyViewFactorSource,
    SlopeGradientSource,
    VisualizationSource,
)
from rvt.blending.normalization import Normalization, NormalizationMode
from rvt.visualization import (
    Hillshade,
    LocalDominance,
    MultiHillshade,
    NegativeOpenness,
    NoiseRemoval,
    SimpleLocalReliefModel,
    SkyIllumination,
    SkyViewFactor,
)


def _value_normalization(minimum: float, maximum: float) -> Normalization:
    """Return an absolute-value normalization, as every shipped terrain uses."""
    return Normalization(NormalizationMode.VALUE, minimum, maximum)


@dataclass(frozen=True)
class TerrainSettings:
    """
    Visualization parameters and normalization bounds recommended for one
    terrain type. Every override is optional; ``None`` leaves the value the
    combination already carries untouched.

    The sky-view factor parameters are also applied to the openness sources,
    because openness is computed by the same horizon search (as in version 2,
    where all of them read ``svf_r_max`` and ``svf_noise``).

    :param name: Name of the terrain type.
    :param hillshade_sun_elevation: Solar vertical angle in degrees for
        :class:`rvt.visualization.Hillshade` layers.
    :param multi_hillshade_sun_elevation: Solar vertical angle in degrees for
        :class:`rvt.visualization.MultiHillshade` and
        :class:`rvt.blending.MultiDirectionalHillshadeRgbSource` layers.
    :param simple_local_relief_model_radius_cells: Trend removal radius in
        pixels for :class:`rvt.visualization.SimpleLocalReliefModel` layers.
    :param sky_view_factor_max_radius: Horizon search radius in pixels for
        sky-view factor and openness layers.
    :param sky_view_factor_noise_removal: Noise removal level for sky-view
        factor and openness layers.
    :param local_dominance_minimum_radius: Smallest observation distance in
        pixels for :class:`rvt.visualization.LocalDominance` layers.
    :param local_dominance_maximum_radius: Largest observation distance in
        pixels for :class:`rvt.visualization.LocalDominance` layers.
    :param local_dominance_observer_height: Observer height in meters for
        :class:`rvt.visualization.LocalDominance` layers.
    :param slope_gradient_normalization: Bounds for slope gradient layers.
    :param hillshade_normalization: Bounds for hillshade layers.
    :param multi_hillshade_normalization: Bounds for multiple directions
        hillshade layers.
    :param simple_local_relief_model_normalization: Bounds for simple local
        relief model layers.
    :param sky_view_factor_normalization: Bounds for sky-view factor layers.
    :param anisotropic_sky_view_factor_normalization: Bounds for sky-view
        factor layers whose config requests anisotropy.
    :param positive_openness_normalization: Bounds for positive openness
        layers.
    :param negative_openness_normalization: Bounds for negative openness
        layers.
    :param sky_illumination_normalization: Bounds for sky illumination layers.
    :param local_dominance_normalization: Bounds for local dominance layers.
    """

    name: str

    hillshade_sun_elevation: float | None = None
    multi_hillshade_sun_elevation: float | None = None
    simple_local_relief_model_radius_cells: int | None = None
    sky_view_factor_max_radius: int | None = None
    sky_view_factor_noise_removal: NoiseRemoval | None = None
    local_dominance_minimum_radius: int | None = None
    local_dominance_maximum_radius: int | None = None
    local_dominance_observer_height: float | None = None

    slope_gradient_normalization: Normalization | None = None
    hillshade_normalization: Normalization | None = None
    multi_hillshade_normalization: Normalization | None = None
    simple_local_relief_model_normalization: Normalization | None = None
    sky_view_factor_normalization: Normalization | None = None
    anisotropic_sky_view_factor_normalization: Normalization | None = None
    positive_openness_normalization: Normalization | None = None
    negative_openness_normalization: Normalization | None = None
    sky_illumination_normalization: Normalization | None = None
    local_dominance_normalization: Normalization | None = None

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("A terrain settings name must not be empty!")
        if (
            self.local_dominance_minimum_radius is not None
            and self.local_dominance_maximum_radius is not None
            and self.local_dominance_minimum_radius
            > self.local_dominance_maximum_radius
        ):
            raise ValueError(
                "local_dominance_minimum_radius must not exceed"
                f" local_dominance_maximum_radius, got"
                f" '{self.local_dominance_minimum_radius}' and"
                f" '{self.local_dominance_maximum_radius}'!"
            )

    def apply(self, combination: BlendCombination) -> BlendCombination:
        """
        Apply this terrain to a blend combination.

        Every layer whose visualization kind this terrain has settings for gets
        its visualization parameters and its normalization bounds overridden;
        all other layers, and the input combination itself, are left untouched.
        Nested combinations (:class:`rvt.blending.BlendedSource`) are adjusted
        recursively. Inside a :class:`rvt.blending.ColormappedSource` only the
        wrapped visualization parameters are overridden: its normalization is
        part of the color mapping recipe, not of the terrain.

        :param combination: The combination to adjust. Not modified.
        :returns: A new combination with the terrain's overrides applied.
        """
        return replace(
            combination,
            layers=tuple(self._apply_to_layer(layer) for layer in combination.layers),
        )

    def _apply_to_layer(self, layer: BlendLayer) -> BlendLayer:
        normalization = self._normalization_for(layer.source)
        return replace(
            layer,
            source=self._apply_to_source(layer.source),
            normalization=normalization
            if normalization is not None
            else layer.normalization,
        )

    def _apply_to_source(self, source: LayerSource) -> LayerSource:
        match source:
            case VisualizationSource():
                return replace(source, config=self._apply_to_config(source.config))
            case MultiDirectionalHillshadeRgbSource():
                if self.multi_hillshade_sun_elevation is None:
                    return source
                return replace(source, sun_elevation=self.multi_hillshade_sun_elevation)
            case SkyViewFactorSource():
                return replace(source, config=self._apply_to_config(source.config))
            case PositiveOpennessSource():
                return replace(source, config=self._apply_to_config(source.config))
            case NegativeOpennessSource():
                return replace(source, config=self._apply_to_config(source.config))
            case OpennessDifferenceSource():
                return replace(
                    source,
                    positive_config=self._apply_to_config(source.positive_config),
                    negative_config=self._apply_to_config(source.negative_config),
                )
            case ColormappedSource():
                return replace(source, source=self._apply_to_source(source.source))
            case BlendedSource():
                return replace(source, combination=self.apply(source.combination))
            case ImageSource() | SlopeGradientSource() | ShadowSource():
                return source

    def _apply_to_config[ConfigT](self, config: ConfigT) -> ConfigT:
        """Override the parameters of one visualization config, by its kind."""
        match config:
            case Hillshade():
                if self.hillshade_sun_elevation is None:
                    return config
                return replace(config, sun_elevation=self.hillshade_sun_elevation)
            case MultiHillshade():
                if self.multi_hillshade_sun_elevation is None:
                    return config
                return replace(config, sun_elevation=self.multi_hillshade_sun_elevation)
            case SimpleLocalReliefModel():
                if self.simple_local_relief_model_radius_cells is None:
                    return config
                return replace(
                    config,
                    radius_cells=self.simple_local_relief_model_radius_cells,
                )
            case SkyViewFactor() | NegativeOpenness():
                return replace(
                    config,
                    max_radius=(
                        self.sky_view_factor_max_radius
                        if self.sky_view_factor_max_radius is not None
                        else config.max_radius
                    ),
                    noise_removal=(
                        self.sky_view_factor_noise_removal
                        if self.sky_view_factor_noise_removal is not None
                        else config.noise_removal
                    ),
                )
            case LocalDominance():
                return replace(
                    config,
                    minimum_radius=(
                        self.local_dominance_minimum_radius
                        if self.local_dominance_minimum_radius is not None
                        else config.minimum_radius
                    ),
                    maximum_radius=(
                        self.local_dominance_maximum_radius
                        if self.local_dominance_maximum_radius is not None
                        else config.maximum_radius
                    ),
                    observer_height=(
                        self.local_dominance_observer_height
                        if self.local_dominance_observer_height is not None
                        else config.observer_height
                    ),
                )
            case _:
                return config

    def _normalization_for(self, source: LayerSource) -> Normalization | None:
        """Return the bounds this terrain prescribes for a layer source, if any."""
        match source:
            case SlopeGradientSource():
                return self.slope_gradient_normalization
            case MultiDirectionalHillshadeRgbSource():
                return self.multi_hillshade_normalization
            case SkyViewFactorSource():
                if source.config.include_anisotropic:
                    return self.anisotropic_sky_view_factor_normalization
                return self.sky_view_factor_normalization
            case PositiveOpennessSource():
                return self.positive_openness_normalization
            case NegativeOpennessSource():
                return self.negative_openness_normalization
            case VisualizationSource():
                return self._normalization_for_config(source.config)
            case (
                ImageSource()
                | ShadowSource()
                | OpennessDifferenceSource()
                | ColormappedSource()
                | BlendedSource()
            ):
                return None

    def _normalization_for_config(self, config: object) -> Normalization | None:
        """Return the bounds this terrain prescribes for a visualization config, if any."""
        match config:
            case Hillshade():
                return self.hillshade_normalization
            case MultiHillshade():
                return self.multi_hillshade_normalization
            case SimpleLocalReliefModel():
                return self.simple_local_relief_model_normalization
            case SkyIllumination():
                return self.sky_illumination_normalization
            case LocalDominance():
                return self.local_dominance_normalization
            case _:
                return None


_GENERAL_TERRAIN_SETTINGS = TerrainSettings(
    name="general",
    hillshade_sun_elevation=35.0,
    multi_hillshade_sun_elevation=35.0,
    simple_local_relief_model_radius_cells=20,
    sky_view_factor_max_radius=10,
    sky_view_factor_noise_removal=NoiseRemoval.NONE,
    local_dominance_minimum_radius=10,
    local_dominance_maximum_radius=20,
    local_dominance_observer_height=1.7,
    slope_gradient_normalization=_value_normalization(0.0, 50.0),
    hillshade_normalization=_value_normalization(0.0, 1.0),
    multi_hillshade_normalization=_value_normalization(0.0, 1.0),
    simple_local_relief_model_normalization=_value_normalization(-1.0, 1.0),
    sky_view_factor_normalization=_value_normalization(0.7, 1.0),
    anisotropic_sky_view_factor_normalization=_value_normalization(0.65, 1.0),
    positive_openness_normalization=_value_normalization(68.0, 93.0),
    negative_openness_normalization=_value_normalization(60.0, 95.0),
    sky_illumination_normalization=_value_normalization(0.25, 0.3),
    local_dominance_normalization=_value_normalization(0.5, 1.8),
)

_FLAT_TERRAIN_SETTINGS = TerrainSettings(
    name="flat",
    hillshade_sun_elevation=15.0,
    multi_hillshade_sun_elevation=15.0,
    simple_local_relief_model_radius_cells=10,
    sky_view_factor_max_radius=20,
    sky_view_factor_noise_removal=NoiseRemoval.HIGH,
    local_dominance_minimum_radius=10,
    local_dominance_maximum_radius=20,
    local_dominance_observer_height=1.7,
    slope_gradient_normalization=_value_normalization(0.0, 15.0),
    sky_view_factor_normalization=_value_normalization(0.9, 1.0),
    positive_openness_normalization=_value_normalization(85.0, 93.0),
    negative_openness_normalization=_value_normalization(75.0, 95.0),
    local_dominance_normalization=_value_normalization(0.5, 3.0),
)

_STEEP_TERRAIN_SETTINGS = TerrainSettings(
    name="steep",
    hillshade_sun_elevation=55.0,
    multi_hillshade_sun_elevation=55.0,
    simple_local_relief_model_radius_cells=10,
    sky_view_factor_max_radius=10,
    sky_view_factor_noise_removal=NoiseRemoval.NONE,
    local_dominance_minimum_radius=10,
    local_dominance_maximum_radius=10,
    local_dominance_observer_height=16.0,
    slope_gradient_normalization=_value_normalization(0.0, 60.0),
    sky_view_factor_normalization=_value_normalization(0.55, 1.0),
    positive_openness_normalization=_value_normalization(55.0, 95.0),
    negative_openness_normalization=_value_normalization(45.0, 95.0),
    local_dominance_normalization=_value_normalization(0.55, 0.95),
)


class Terrain(StrEnum):
    """
    The terrain types RVT ships settings for: a diverse landscape of mountains
    and plains (``GENERAL``), very flat terrain (``FLAT``), and steep terrain
    (``STEEP``).
    """

    GENERAL = auto()
    FLAT = auto()
    STEEP = auto()

    def settings(self) -> TerrainSettings:
        """
        Return the shipped settings of this terrain type (the values of the
        version 2 ``default_terrains_settings.json``).

        :returns: The terrain settings, ready to apply to a combination.
        """
        match self:
            case Terrain.GENERAL:
                return _GENERAL_TERRAIN_SETTINGS
            case Terrain.FLAT:
                return _FLAT_TERRAIN_SETTINGS
            case Terrain.STEEP:
                return _STEEP_TERRAIN_SETTINGS
