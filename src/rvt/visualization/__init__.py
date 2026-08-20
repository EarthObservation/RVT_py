"""
Visualization configs: one frozen dataclass per visualization, each with a
``compute`` method over a float32 NaN-nodata DEM array.
"""

from rvt.visualization.hillshade import Hillshade
from rvt.visualization.slope import SlopeAspect, SlopeAspectResult, SlopeUnits

__all__ = [
    "Hillshade",
    "SlopeAspect",
    "SlopeAspectResult",
    "SlopeUnits",
]
