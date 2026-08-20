"""
Visualization configs: one frozen dataclass per visualization, each with a
``compute`` method over a float32 NaN-nodata DEM array.
"""

from rvt.visualization.hillshade import Hillshade
from rvt.visualization.local_dominance import LocalDominance
from rvt.visualization.msrm import MultiScaleReliefModel
from rvt.visualization.multi_hillshade import MultiHillshade
from rvt.visualization.slope import SlopeAspect, SlopeAspectResult, SlopeUnits
from rvt.visualization.slrm import SimpleLocalReliefModel

__all__ = [
    "Hillshade",
    "LocalDominance",
    "MultiHillshade",
    "MultiScaleReliefModel",
    "SimpleLocalReliefModel",
    "SlopeAspect",
    "SlopeAspectResult",
    "SlopeUnits",
]
