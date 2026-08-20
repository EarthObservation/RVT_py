"""
Visualization configs: one frozen dataclass per visualization, each with a
``compute`` method over a float32 NaN-nodata DEM array.
"""

from rvt.visualization.hillshade import Hillshade
from rvt.visualization.local_dominance import LocalDominance
from rvt.visualization.msrm import MultiScaleReliefModel
from rvt.visualization.mstp import MultiScaleTopographicPosition
from rvt.visualization.multi_hillshade import MultiHillshade
from rvt.visualization.shadow import ShadowHorizon, ShadowHorizonResult
from rvt.visualization.sky_illumination import SkyIllumination, SkyModel
from rvt.visualization.sky_view_factor import (
    NegativeOpenness,
    NoiseRemoval,
    SkyViewFactor,
    SkyViewFactorResult,
)
from rvt.visualization.slope import SlopeAspect, SlopeAspectResult, SlopeUnits
from rvt.visualization.slrm import SimpleLocalReliefModel

__all__ = [
    "Hillshade",
    "LocalDominance",
    "MultiHillshade",
    "MultiScaleReliefModel",
    "MultiScaleTopographicPosition",
    "NegativeOpenness",
    "NoiseRemoval",
    "ShadowHorizon",
    "ShadowHorizonResult",
    "SimpleLocalReliefModel",
    "SkyIllumination",
    "SkyModel",
    "SkyViewFactor",
    "SkyViewFactorResult",
    "SlopeAspect",
    "SlopeAspectResult",
    "SlopeUnits",
]
