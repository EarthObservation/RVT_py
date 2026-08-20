"""
Visualization configs: one frozen dataclass per visualization, each with a
``compute`` method over a float32 NaN-nodata DEM array.
"""

from rvt.visualization.hillshade import Hillshade
from rvt.visualization.local_dominance import LocalDominance
from rvt.visualization.msrm import MultiScaleReliefModel
from rvt.visualization.multi_hillshade import MultiHillshade
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
    "NegativeOpenness",
    "NoiseRemoval",
    "SimpleLocalReliefModel",
    "SkyViewFactor",
    "SkyViewFactorResult",
    "SlopeAspect",
    "SlopeAspectResult",
    "SlopeUnits",
]
