"""
Visualization configs: one frozen dataclass per visualization, each with a
``compute`` method over a float32 NaN-nodata DEM array.
"""

from rvt.visualization.hillshade import Hillshade
from rvt.visualization.multi_hillshade import MultiHillshade
from rvt.visualization.slope import SlopeAspect, SlopeAspectResult, SlopeUnits
from rvt.visualization.slrm import SimpleLocalReliefModel

__all__ = [
    "Hillshade",
    "MultiHillshade",
    "SimpleLocalReliefModel",
    "SlopeAspect",
    "SlopeAspectResult",
    "SlopeUnits",
]
