from enum import Enum


class RVTVisualization(Enum):
    """Visualizations supported in RVT and their abbreviations."""

    SLOPE = "slp"
    HILLSHADE = "hs"
    SHADOW = "shd"
    MULTIPLE_DIRECTIONS_HILLSHADE = "mhs"
    SIMPLE_LOCAL_RELIEF_MODEL = "slrm"
    SKY_VIEW_FACTOR = "svf"
    ANISOTROPIC_SKY_VIEW_FACTOR = "asvf"
    POSITIVE_OPENNESS = "pos_opns"
    NEGATIVE_OPENNESS = "neg_opns"
    SKY_ILLUMINATION = "sim"
    LOCAL_DOMINANCE = "ld"
    MULTI_SCALE_RELIEF_MODEL = "msrm"
    MULTI_SCALE_TOPOGRAPHIC_POSITION = "mstp"


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
