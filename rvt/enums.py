from enum import Enum


class RVTVisualization(Enum):
    """Visualizations supported in RVT."""

    SLOPE = "slope"
    HILLSHADE = "hillshade"
    SHADOW = "shadow"
    MULTIPLE_DIRECTIONS_HILLSHADE = "multiple_directions_hillshade"
    SIMPLE_LOCAL_RELIEF_MODEL = "simple_local_relief_model"
    SKY_VIEW_FACTOR = "sky_view_factor"
    ANISOTROPIC_SKY_VIEW_FACTOR = "anisotropic_sky_view_factor"
    POSITIVE_OPENNESS = "positive_openness"
    NEGATIVE_OPENNESS = "negative_openness"
    SKY_ILLUMINATION = "sky_illumination"
    LOCAL_DOMINANCE = "local_dominance"
    MULTI_SCALE_RELIEF_MODEL = "multi_scale_relief_model"
    MULTI_SCALE_TOPOGRAPHIC_POSITION = "multi_scale_topographic_position"


class SlopeOutputUnit(Enum):
    PERCENT = "percent"
    DEGREE = "degree"
    RADIAN = "radian"


class SvfNoiseRemove(Enum):
    NO_REMOVE = "no_remove"
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"

    @property
    def maximum_radius_ignore_percentage(self):
        """
        The portion (percent) of the maximal search radius to ignore in horizon estimation, for selected noise level.
        """
        if self == SvfNoiseRemove.NO_REMOVE:
            return 0.0
        elif self == SvfNoiseRemove.LOW:
            return 10.0
        elif self == SvfNoiseRemove.MEDIUM:
            return 20.0
        elif self == SvfNoiseRemove.HIGH:
            return 40.0
        else:
            raise ValueError(
                f"Percentage for selected noise level ({self.name}) is unspecified!"
            )


class AnisotropyLevel(Enum):
    LOW = "low"
    HIGH = "high"

    @property
    def polynomial_level(self):
        """Level of polynomial that determines the anisotropy, for selected anisotropy level."""
        if self == AnisotropyLevel.LOW:
            return 4
        elif self == AnisotropyLevel.HIGH:
            return 8

    @property
    def minimal_weight(self):
        """Minimal weight for selected anisotropy level."""
        if self == AnisotropyLevel.LOW:
            return 0.4
        elif self == AnisotropyLevel.HIGH:
            return 0.1


class OpennessType(Enum):
    POSITIVE = "positive"
    NEGATIVE = "negative"


class NormalizationMode(Enum):
    PERCENT = "percent"
    VALUE = "value"
