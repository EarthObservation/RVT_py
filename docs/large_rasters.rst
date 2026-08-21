Large rasters and batch runs
============================

Two facilities cover work that does not fit a single in-memory
``compute`` call: :class:`rvt.Session`, which computes and saves a list of
visualizations for one DEM, and :func:`rvt.tiling.compute_tiled`, which
processes a DEM window by window. A session uses tiled processing on its own
when the DEM is large enough.

Sessions
--------

A session reads the DEM once, computes each task, and writes every output as a
GeoTIFF next to the DEM (or into ``output_directory``):

.. code-block:: python

   import rvt
   from rvt.visualization import Hillshade, SimpleLocalReliefModel, SkyViewFactor

   session = rvt.Session(
       dem_path="dem.tif",
       tasks=[
           rvt.SessionTask(Hillshade(), save_8bit=True),
           rvt.SessionTask(SimpleLocalReliefModel()),
           rvt.SessionTask(SkyViewFactor(include_positive_openness=True)),
       ],
       output_directory="visualizations",
   )
   written_paths = session.run()

Each :class:`rvt.SessionTask` chooses whether the float raster, the 8-bit
rendering, or both are written, and can override the output name. Outputs whose
file already exists are skipped rather than recomputed, unless the session was
created with ``overwrite=True``; :meth:`rvt.Session.planned_paths` reports what
a run would produce without computing anything.

``vertical_exaggeration`` given to the session replaces the field on every
config, so a whole batch can be re-run with one exaggeration factor without
rebuilding the configs.

Parameter-encoded file names
````````````````````````````

Output names follow the version 2 scheme: the DEM's stem plus a short
visualization code and its parameters, with ``_8bit`` appended for the 8-bit
rendering. For example:

=================================================  =========================================================
Config                                             File name
=================================================  =========================================================
``Hillshade()``                                    ``dem_HS_A315_H35.tif``
``MultiHillshade()``                               ``dem_MULTI-HS_D16_H35.tif``
``SlopeAspect()``                                  ``dem_SLOPE.tif``, ``dem_ASPECT.tif``
``SimpleLocalReliefModel()``                       ``dem_SLRM_R20.tif``
``SkyViewFactor()``                                ``dem_SVF_R10_D16.tif``
``SkyViewFactor(include_anisotropic=True)``        ``dem_SVF-A_R10_D16_A315_ALlow.tif``
``SkyViewFactor(include_positive_openness=True)``  ``dem_OPEN-POS_R10_D16.tif``
``NegativeOpenness()``                             ``dem_OPEN-NEG_R10_D16.tif``
``LocalDominance()``                               ``dem_LD_R_M10-20_DI1_A15_OH1.7.tif``
``SkyIllumination()``                              ``dem_SIM_overcast_D32_100px.tif``
``MultiScaleReliefModel(0, 20, 2)``                ``dem_MSRM_F_M0-20_S2.tif``
``MultiScaleTopographicPosition()``                ``dem_MSTP_5_50_500_L1.2.tif``
``ShadowHorizon()``                                ``dem_shadow_A315_H35.tif``, ``dem_horizon_A315_H35.tif``
=================================================  =========================================================

The naming lives on the configs themselves (``output_stem``, or one method per
member for the multi-output configs), so a name can be predicted without a
session.

.. autoclass:: rvt.Session
   :members:
   :no-index:

.. autoclass:: rvt.SessionTask
   :members:
   :no-index:

Tiled processing
----------------

:func:`rvt.tiling.compute_tiled` reads a DEM window by window, computes a
visualization per tile, and writes each result into a fixed window of the output
raster, so a DEM far larger than memory can be processed:

.. code-block:: python

   from rvt.tiling import compute_tiled
   from rvt.visualization import SimpleLocalReliefModel

   compute_tiled(
       SimpleLocalReliefModel(),
       "big_dem.tif",
       "big_dem_slrm.tif",
       tile_size=(4000, 4000),
       workers=4,
   )

Tiled output equals untiled output. Every window is read expanded by the
config's ``required_overlap(resolution)`` pixels — its neighbourhood radius —
clamped to the raster bounds, and the overlap is sliced off before the core
window is written. Tile borders therefore see real data, and raster borders fall
back on each visualization's own edge padding, exactly as an untiled compute
does. With ``workers`` above 1 the tiles are computed by a process pool; the
result is identical either way, because each tile goes to a fixed window of the
pre-created output raster.

``required_overlap`` replaces version 2's hardcoded overlap table. Note that the
multi-scale relief model derives its overlap from the cell size instead of
treating ``feature_maximum`` (map units) as a pixel count, which is what version
2 did — that buffer was too small below one-metre cells and produced visible
seams.

What can be tiled
`````````````````

Every config declares a ``TILEABLE`` class attribute. Tiling is refused, rather
than silently producing seams, when it cannot be exact:

*   :class:`rvt.visualization.SkyIllumination` — the overcast model divides the
    result by the array's own maximum, and the horizon pyramid pads with the
    array's own minimum elevation; both differ between a tile and the whole
    raster.
*   :class:`rvt.visualization.ShadowHorizon` — the horizon pyramid pads with the
    tile's own minimum elevation.
*   :class:`rvt.visualization.MultiHillshade` — the output shape differs from
    the input (a one-pixel border is trimmed), so it cannot be written back into
    fixed windows.

These raise ``ValueError``. The multi-output configs
(:class:`rvt.visualization.SlopeAspect`,
:class:`rvt.visualization.SkyViewFactor`,
:class:`rvt.visualization.ShadowHorizon`) raise ``TypeError``: ``compute_tiled``
writes one image, so compute those untiled through their own ``compute``.

Percentile-based conversions are data-dependent in the same way. A tiled float
result is stretched to 8-bit from the finished raster, not per tile.

.. autofunction:: rvt.tiling.compute_tiled
   :no-index:

.. autoclass:: rvt.tiling.TileableVisualization
   :members:
   :no-index:

Automatic switchover in a session
`````````````````````````````````

A session checks the DEM's pixel count before reading it. Above
``automatic_tiling_threshold`` (default :data:`rvt.session.DEFAULT_TILING_THRESHOLD`,
100 000 000 pixels — the version 2 limit of 10 000 × 10 000) every task whose
config is tileable and single-output is computed through
:func:`rvt.tiling.compute_tiled` instead of reading the whole DEM. The remaining
tasks still read the DEM whole, as does the 8-bit rendering of a tiled result,
whose stretch needs the finished image.

Pass ``automatic_tiling_threshold=None`` to always compute in memory, or a lower
number to force tiling earlier; ``tile_size`` and ``workers`` are forwarded to
the tiled runner.
