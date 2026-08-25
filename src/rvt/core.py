"""
Shared core types for the RVT library.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar, NamedTuple

import numpy as np
from numpy.typing import NDArray

from rvt.stretch import Stretch

MAXIMUM_VERTICAL_EXAGGERATION = 10000.0


class Resolution(NamedTuple):
    """
    Raster cell size in map units.

    :param x: Cell size in the X direction (always positive).
    :param y: Cell size in the Y direction (always positive).
    """

    x: float
    y: float


def format_parameter(value: float) -> str:
    """
    Format a numeric parameter for a file name: integral values render without
    a decimal part (``315.0`` becomes ``315``), others as-is (``22.5``).

    :param value: The parameter value to encode.
    :returns: The compact string form of the value.
    """
    if value == int(value):
        return str(int(value))
    return str(value)


def validate_vertical_exaggeration(vertical_exaggeration: float) -> None:
    """
    Raise ``ValueError`` when the exaggeration factor is outside ±10000.

    :param vertical_exaggeration: Vertical exaggeration factor applied to the DEM.
    :raises ValueError: If the value is outside the supported range.
    """
    if not (
        -MAXIMUM_VERTICAL_EXAGGERATION
        <= vertical_exaggeration
        <= MAXIMUM_VERTICAL_EXAGGERATION
    ):
        raise ValueError(
            f"vertical_exaggeration must be between -{MAXIMUM_VERTICAL_EXAGGERATION:g}"
            f" and {MAXIMUM_VERTICAL_EXAGGERATION:g}, got '{vertical_exaggeration}'!"
        )


def validate_dem_and_resolution(dem_ndim: int, resolution: Resolution) -> None:
    """
    Raise ``ValueError`` for a non-2D DEM or a non-positive resolution.

    :param dem_ndim: Number of dimensions of the DEM array.
    :param resolution: Cell size of the DEM.
    :raises ValueError: If the DEM is not 2D or the resolution is not positive.
    """
    if dem_ndim != 2:
        raise ValueError(f"dem has to be a 2D array, got {dem_ndim} dimensions!")
    if resolution.x <= 0 or resolution.y <= 0:
        raise ValueError(f"resolution must be positive, got '{resolution}'!")


class VisualizationOutput(NamedTuple):
    """
    One file-producing member of a visualization: where it goes and how it is
    rendered to 8-bit.

    :param member: Attribute of the ``compute`` result this output reads, or
        ``None`` when ``compute`` returns the array itself.
    :param stem: Parameter-encoded output file stem, without extension.
    :param eight_bit_stretch: Stretch mapping the float member onto 0-255, or
        ``None`` for a member with no legacy 8-bit form (aspect, horizon).
    """

    member: str | None
    stem: str
    eight_bit_stretch: Stretch | None


class Visualization[ResultT](ABC):
    """
    The contract every visualization config satisfies. To add a visualization,
    subclass this as a frozen dataclass and provide:

    - the parameters as fields, validated in ``__post_init__``, with
      ``vertical_exaggeration: float = 1.0`` last (reuse
      :func:`validate_vertical_exaggeration`);
    - ``TILEABLE``, a ``ClassVar[bool]``: ``True`` only when a pixel's result
      depends on nothing outside :meth:`required_overlap` — no data-dependent
      normalisation (e.g. dividing by the array's own maximum) and no output
      shape change, so tiled output equals untiled output;
    - :meth:`outputs`, one :class:`VisualizationOutput` per file produced;
    - :meth:`required_overlap`, the neighbourhood radius in pixels;
    - :meth:`compute`, taking a 2D float32 DEM with NaN as nodata (not
      modified; call :func:`validate_dem_and_resolution` first) and
      propagating NaN into the float32 result.

    ``ResultT`` is what ``compute`` returns: an array for a single-image
    config, a typed :class:`~typing.NamedTuple` for a multi-output one.

    Then register the new config: export it from ``rvt.visualization``, add it
    to ``VisualizationConfig`` in :mod:`rvt.session` (and to
    ``_MULTI_OUTPUT_CONFIGS`` in :mod:`rvt.tiling` if it is multi-output), and
    add it to ``ALL_CONFIGS`` in ``tests/test_visualization_contract.py``.
    """

    TILEABLE: ClassVar[bool]

    @property
    @abstractmethod
    def vertical_exaggeration(self) -> float:
        """Return the factor the DEM heights are multiplied by."""

    def __init_subclass__(cls, **kwargs: object) -> None:
        """Reject a concrete subclass that forgot to declare ``TILEABLE``."""
        super().__init_subclass__(**kwargs)
        if not isinstance(getattr(cls, "TILEABLE", None), bool):
            raise TypeError(
                f"{cls.__name__} must declare TILEABLE: ClassVar[bool] —"
                " whether its result is independent of the tiling."
            )

    @abstractmethod
    def outputs(self, dem_name: str) -> tuple[VisualizationOutput, ...]:
        """Return the outputs this config produces, in file-name order."""

    @abstractmethod
    def required_overlap(self, resolution: Resolution) -> int:
        """Return the neighbourhood radius in pixels a tile must be expanded by."""

    @abstractmethod
    def compute(self, dem: NDArray[np.float32], resolution: Resolution) -> ResultT:
        """Compute the visualization of a DEM."""
