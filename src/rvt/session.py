"""
Batch file-in/file-out sessions: read a DEM once, compute a list of
visualizations, and write each output as a GeoTIFF with a parameter-encoded
file name (e.g. ``dem_HS_A315_H35.tif``), replacing the legacy
``DefaultValues.save_*`` workflow.

Tiled processing of DEMs too large for memory is forthcoming; sessions
currently read the whole DEM at once.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import NamedTuple, Sequence

import numpy as np

from rvt.io import Raster, read, write
from rvt.stretch import Stretch, to_8bit
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

VisualizationConfig = (
    Hillshade
    | LocalDominance
    | MultiHillshade
    | MultiScaleReliefModel
    | MultiScaleTopographicPosition
    | NegativeOpenness
    | ShadowHorizon
    | SimpleLocalReliefModel
    | SkyIllumination
    | SkyViewFactor
    | SlopeAspect
)


@dataclass(frozen=True)
class SessionTask:
    """
    One visualization to compute and save within a :class:`Session`.

    :param config: The visualization config to compute.
    :param save_float: Write the float32 GeoTIFF (NaN nodata).
    :param save_8bit: Also write the 8-bit rendering (``_8bit`` suffix), using
        the config's ``EIGHT_BIT_STRETCH``. Multi-output members without a
        legacy 8-bit form (aspect, horizon) are written as float only.
    :param name: Override for the parameter-encoded output file stem. For
        multi-output configs the member name is appended, e.g. ``name_slope``.
    :raises ValueError: If neither ``save_float`` nor ``save_8bit`` is set.
    """

    config: VisualizationConfig
    save_float: bool = True
    save_8bit: bool = False
    name: str | None = None

    def __post_init__(self) -> None:
        if not self.save_float and not self.save_8bit:
            raise ValueError(
                "at least one of save_float and save_8bit must be enabled!"
            )


class _PlannedOutput(NamedTuple):
    """One file-producing member of a task: where it goes and how to render it."""

    member: str | None
    stem: str
    stretch: Stretch | None


def _planned_outputs(task: SessionTask, dem_name: str) -> list[_PlannedOutput]:
    """Return the outputs a task produces, honouring its name override."""
    config = task.config
    match config:
        case SlopeAspect():
            outputs = [
                _PlannedOutput(
                    "slope",
                    config.slope_output_stem(dem_name),
                    SlopeAspect.EIGHT_BIT_STRETCH,
                ),
                _PlannedOutput("aspect", config.aspect_output_stem(dem_name), None),
            ]
        case SkyViewFactor():
            outputs = [
                _PlannedOutput(
                    "sky_view_factor",
                    config.sky_view_factor_output_stem(dem_name),
                    SkyViewFactor.EIGHT_BIT_STRETCH,
                )
            ]
            if config.include_anisotropic:
                outputs.append(
                    _PlannedOutput(
                        "anisotropic",
                        config.anisotropic_output_stem(dem_name),
                        SkyViewFactor.ANISOTROPIC_EIGHT_BIT_STRETCH,
                    )
                )
            if config.include_positive_openness:
                outputs.append(
                    _PlannedOutput(
                        "positive_openness",
                        config.positive_openness_output_stem(dem_name),
                        SkyViewFactor.POSITIVE_OPENNESS_EIGHT_BIT_STRETCH,
                    )
                )
        case ShadowHorizon():
            outputs = [
                _PlannedOutput(
                    "shadow",
                    config.shadow_output_stem(dem_name),
                    ShadowHorizon.EIGHT_BIT_STRETCH,
                ),
                _PlannedOutput("horizon", config.horizon_output_stem(dem_name), None),
            ]
        case MultiHillshade():
            # The legacy 8-bit multi-hillshade is a special 3-band RGB rendering
            # from the DEM, not a stretch of the float stack; not implemented yet.
            outputs = [_PlannedOutput(None, config.output_stem(dem_name), None)]
        case _:
            outputs = [
                _PlannedOutput(
                    None, config.output_stem(dem_name), config.EIGHT_BIT_STRETCH
                )
            ]

    if task.name is None:
        return outputs
    return [
        output._replace(
            stem=task.name if output.member is None else f"{task.name}_{output.member}"
        )
        for output in outputs
    ]


class _Target(NamedTuple):
    """One file to write: its path, source output, and rendering kind."""

    path: Path
    output: _PlannedOutput
    is_8bit: bool


class Session:
    """
    A batch run: read a DEM, compute visualizations, write GeoTIFFs.

    Output file names encode the visualization parameters in the legacy scheme
    (e.g. ``dem_HS_A315_H35.tif``, ``dem_SVF_R10_D16.tif``, ``_8bit`` suffix
    for 8-bit renderings), overridable per task via :class:`SessionTask`'s
    ``name``. Existing files are skipped unless ``overwrite`` is set.

    Tiled processing of DEMs too large for memory is forthcoming; the whole
    DEM is currently read and computed in memory.

    :param dem_path: Path of the DEM raster to visualize.
    :param tasks: The visualizations to compute and save.
    :param output_directory: Directory the outputs are written to, created if
        missing. Defaults to the DEM's directory.
    :param overwrite: Rewrite outputs whose file already exists.
    :param vertical_exaggeration: When given, replaces every config's
        ``vertical_exaggeration`` at construction; ``None`` leaves the configs
        as they are.
    :raises ValueError: If a task requests no output at all.
    :raises NotImplementedError: If a task requests the 8-bit multi-hillshade
        (the legacy 3-band RGB rendering), which is not implemented yet.
    """

    def __init__(
        self,
        dem_path: str | Path,
        tasks: Sequence[SessionTask],
        *,
        output_directory: str | Path | None = None,
        overwrite: bool = False,
        vertical_exaggeration: float | None = None,
    ) -> None:
        self.dem_path = Path(dem_path)
        if vertical_exaggeration is not None:
            tasks = [
                replace(
                    task,
                    config=replace(
                        task.config, vertical_exaggeration=vertical_exaggeration
                    ),
                )
                for task in tasks
            ]
        self.tasks = tuple(tasks)
        self.output_directory = (
            Path(output_directory)
            if output_directory is not None
            else self.dem_path.parent
        )
        self.overwrite = overwrite
        for task in self.tasks:
            if task.save_8bit and isinstance(task.config, MultiHillshade):
                raise NotImplementedError(
                    "The 8-bit multi-hillshade is a special legacy rendering"
                    " (3-band RGB computed from the DEM, not a stretch of the"
                    " float stack) and is not implemented yet; save the float"
                    " stack instead."
                )

    def _task_targets(self, task: SessionTask) -> list[_Target]:
        """Return the files a task would write, ignoring skip-if-exists."""
        targets = []
        for output in _planned_outputs(task, self.dem_path.stem):
            if task.save_float:
                targets.append(
                    _Target(self.output_directory / f"{output.stem}.tif", output, False)
                )
            if task.save_8bit and output.stretch is not None:
                targets.append(
                    _Target(
                        self.output_directory / f"{output.stem}_8bit.tif", output, True
                    )
                )
        return targets

    def planned_paths(self) -> list[Path]:
        """
        Return the paths of every file the session would produce, in task
        order, regardless of whether the files already exist.

        :returns: The planned output file paths.
        """
        return [
            target.path for task in self.tasks for target in self._task_targets(task)
        ]

    def run(self) -> list[Path]:
        """
        Read the DEM, compute each task, and write its outputs.

        Outputs whose file already exists are skipped (not an error) unless
        the session was created with ``overwrite``; a task none of whose files
        need writing is not computed at all.

        :returns: The paths of the files actually written, in task order.
        """
        raster = read(self.dem_path)
        self.output_directory.mkdir(parents=True, exist_ok=True)

        written = []
        for task in self.tasks:
            targets = [
                target
                for target in self._task_targets(task)
                if self.overwrite or not target.path.exists()
            ]
            if not targets:
                continue
            dem = raster.array.astype(np.float32, copy=False)
            result = task.config.compute(dem, raster.resolution)
            for target in targets:
                array = (
                    result
                    if target.output.member is None
                    else getattr(result, target.output.member)
                )
                array = np.asarray(array).astype(np.float32)
                if target.is_8bit:
                    assert target.output.stretch is not None
                    array = to_8bit(array, target.output.stretch)
                write(
                    Raster(array=array, transform=raster.transform, crs=raster.crs),
                    target.path,
                )
                written.append(target.path)
        return written
