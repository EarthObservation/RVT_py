"""
Batch file-in/file-out sessions: read a DEM once, compute a list of
visualizations, and write each output as a GeoTIFF with a parameter-encoded
file name (e.g. ``dem_HS_A315_H35.tif``), replacing the legacy
``DefaultValues.save_*`` workflow.

DEMs above :data:`DEFAULT_TILING_THRESHOLD` pixels are processed tile by tile
(see :mod:`rvt.tiling`) rather than read whole.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import NamedTuple, Sequence, cast

import numpy as np
import rasterio

from rvt.io import Raster, read, write
from rvt.stretch import Stretch, to_8bit
from rvt.tiling import DEFAULT_TILE_SIZE, TileableVisualization, compute_tiled
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

# The legacy switch-over point: 10000 x 10000 pixels.
DEFAULT_TILING_THRESHOLD = 100_000_000


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

    DEMs larger than ``automatic_tiling_threshold`` pixels are processed tile
    by tile through :func:`rvt.tiling.compute_tiled` instead of being read
    whole. Tasks whose visualization cannot be tiled (a multi-output config, or
    one declaring ``TILEABLE = False``) always read the whole DEM, as does the
    8-bit rendering of a tiled result, whose stretch needs the finished image.

    :param dem_path: Path of the DEM raster to visualize.
    :param tasks: The visualizations to compute and save.
    :param output_directory: Directory the outputs are written to, created if
        missing. Defaults to the DEM's directory.
    :param overwrite: Rewrite outputs whose file already exists.
    :param vertical_exaggeration: When given, replaces every config's
        ``vertical_exaggeration`` at construction; ``None`` leaves the configs
        as they are.
    :param automatic_tiling_threshold: DEM size in pixels above which tiled
        processing is used; ``None`` always computes in memory.
    :param tile_size: Tile size as ``(rows, columns)`` for tiled processing.
    :param workers: Number of worker processes for tiled processing.
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
        automatic_tiling_threshold: int | None = DEFAULT_TILING_THRESHOLD,
        tile_size: tuple[int, int] = DEFAULT_TILE_SIZE,
        workers: int = 1,
    ) -> None:
        self.dem_path = Path(dem_path)
        self.automatic_tiling_threshold = automatic_tiling_threshold
        self.tile_size = tile_size
        self.workers = workers
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

    def _dem_pixel_count(self) -> int:
        """Return the DEM's pixel count without reading its data."""
        with rasterio.open(self.dem_path) as source:
            return source.width * source.height

    def _tiles_task(self, task: SessionTask) -> bool:
        """
        Return whether a task's float output is computed tile by tile: only a
        tileable, single-output config qualifies.
        """
        if not getattr(task.config, "TILEABLE", False):
            return False
        return all(
            output.member is None
            for output in _planned_outputs(task, self.dem_path.stem)
        )

    def _write_8bit(
        self, float_raster: Raster, targets: Sequence[_Target]
    ) -> list[Path]:
        """Stretch a computed float result to 8-bit and write it."""
        written = []
        array = np.asarray(float_raster.array).astype(np.float32)
        for target in targets:
            assert target.output.stretch is not None
            write(
                Raster(
                    array=to_8bit(array, target.output.stretch),
                    transform=float_raster.transform,
                    crs=float_raster.crs,
                ),
                target.path,
            )
            written.append(target.path)
        return written

    def _run_task_tiled(
        self, task: SessionTask, targets: Sequence[_Target]
    ) -> list[Path]:
        """
        Compute one task tile by tile. The 8-bit rendering is stretched from
        the finished float raster, so it is read back into memory.
        """
        written = []
        float_targets = [target for target in targets if not target.is_8bit]
        eight_bit_targets = [target for target in targets if target.is_8bit]
        # :meth:`_tiles_task` has established that the config is tileable and
        # single-output, which the union type cannot express.
        config = cast(TileableVisualization, task.config)

        float_path = None
        for target in float_targets:
            compute_tiled(
                config,
                self.dem_path,
                target.path,
                tile_size=self.tile_size,
                workers=self.workers,
            )
            float_path = target.path
            written.append(target.path)

        if eight_bit_targets:
            with TemporaryDirectory() as directory:
                if float_path is None:
                    float_path = Path(directory) / f"{self.dem_path.stem}_float.tif"
                    compute_tiled(
                        config,
                        self.dem_path,
                        float_path,
                        tile_size=self.tile_size,
                        workers=self.workers,
                    )
                written.extend(self._write_8bit(read(float_path), eight_bit_targets))
        return written

    def _run_task_in_memory(
        self, task: SessionTask, targets: Sequence[_Target], raster: Raster
    ) -> list[Path]:
        """Compute one task from the whole DEM held in memory."""
        written = []
        dem = raster.array.astype(np.float32, copy=False)
        result = task.config.compute(dem, raster.resolution)
        for target in targets:
            array = (
                result
                if target.output.member is None
                else getattr(result, target.output.member)
            )
            member = Raster(
                array=np.asarray(array).astype(np.float32),
                transform=raster.transform,
                crs=raster.crs,
            )
            if target.is_8bit:
                written.extend(self._write_8bit(member, [target]))
            else:
                write(member, target.path)
                written.append(target.path)
        return written

    def run(self) -> list[Path]:
        """
        Compute each task and write its outputs, tile by tile where the DEM is
        large enough and the visualization allows it.

        Outputs whose file already exists are skipped (not an error) unless
        the session was created with ``overwrite``; a task none of whose files
        need writing is not computed at all.

        :returns: The paths of the files actually written, in task order.
        """
        self.output_directory.mkdir(parents=True, exist_ok=True)
        tiling_wanted = (
            self.automatic_tiling_threshold is not None
            and self._dem_pixel_count() > self.automatic_tiling_threshold
        )
        raster = None

        written = []
        for task in self.tasks:
            targets = [
                target
                for target in self._task_targets(task)
                if self.overwrite or not target.path.exists()
            ]
            if not targets:
                continue
            if tiling_wanted and self._tiles_task(task):
                written.extend(self._run_task_tiled(task, targets))
                continue
            if raster is None:
                raster = read(self.dem_path)
            written.extend(self._run_task_in_memory(task, targets, raster))
        return written
