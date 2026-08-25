Migrating from version 2
========================

Version 3 replaces the three version 2 modules — ``rvt.vis`` (array functions),
``rvt.blend`` (blending) and ``rvt.default`` (the ``DefaultValues`` bag of
parameters and its ``save_*`` workflow) — with a typed API: one frozen config
class per visualization, a rasterio-backed I/O layer, an explicit blend layer
model, and :class:`rvt.Session` for file-in/file-out batches. The numeric
results are verified against version 2 with a golden-master test suite.

Nothing of the version 2 import surface survives, so version 2 scripts have to
be rewritten rather than adjusted. The tables below map every entry point.

Visualizations
--------------

Version 2 functions took a DEM array plus loose keyword parameters. Version 3
puts the parameters on a frozen config object and computes with
``config.compute(dem, resolution)``, where ``resolution`` is a
:class:`rvt.core.Resolution` (both cell sizes, always positive).

==========================================  ===============================================================
Version 2                                   Version 3
==========================================  ===============================================================
``rvt.vis.slope_aspect(...)``               :class:`rvt.visualization.SlopeAspect`
``rvt.vis.hillshade(...)``                  :class:`rvt.visualization.Hillshade`
``rvt.vis.multi_hillshade(...)``            :class:`rvt.visualization.MultiHillshade`
``rvt.vis.slrm(...)``                       :class:`rvt.visualization.SimpleLocalReliefModel`
``rvt.vis.msrm(...)``                       :class:`rvt.visualization.MultiScaleReliefModel`
``rvt.vis.mstp(...)``                       :class:`rvt.visualization.MultiScaleTopographicPosition`
``rvt.vis.sky_view_factor(...)``            :class:`rvt.visualization.SkyViewFactor`
``rvt.vis.sky_view_factor(compute_asvf)``   ``SkyViewFactor(include_anisotropic=True)``
``rvt.vis.sky_view_factor(compute_opns)``   ``SkyViewFactor(include_positive_openness=True)``
``rvt.default.DefaultValues.get_neg_opns``  :class:`rvt.visualization.NegativeOpenness`
``rvt.vis.sky_illumination(...)``           :class:`rvt.visualization.SkyIllumination`
``rvt.vis.shadow_horizon(...)``             :class:`rvt.visualization.ShadowHorizon`
``rvt.vis.byte_scale(...)``                 :func:`rvt.stretch.to_8bit` with a :class:`rvt.stretch.Stretch`
``rvt.vis.fill_where_nan(...)``             *not ported*
==========================================  ===============================================================

Parameter renames follow the config fields: ``ve_factor`` is
``vertical_exaggeration``, ``sun_azimuth``/``sun_elevation`` keep their names,
``nr_directions``/``svf_n_dir`` become ``directions``, ``svf_r_max`` becomes
``max_radius``, ``radius_cell`` becomes ``radius_cells``,
``feature_min``/``feature_max`` become ``feature_minimum``/``feature_maximum``,
and the loose integer/string choices become enums
(:class:`rvt.visualization.SlopeUnits`,
:class:`rvt.visualization.NoiseRemoval`,
:class:`rvt.visualization.SkyModel`).

Dict returns are gone. The multi-output visualizations return typed results
instead: :class:`rvt.visualization.SlopeAspectResult` (``.slope``, ``.aspect``),
:class:`rvt.visualization.SkyViewFactorResult` (``.sky_view_factor``,
``.anisotropic``, ``.positive_openness``) and
:class:`rvt.visualization.ShadowHorizonResult` (``.shadow``, ``.horizon``).

.. code-block:: python

   # version 2
   import rvt.vis

   result = rvt.vis.sky_view_factor(
       dem=dem, resolution=0.5, compute_svf=True, compute_opns=True, svf_r_max=10
   )
   svf = result["svf"]
   opns = result["opns"]

   # version 3
   from rvt.core import Resolution
   from rvt.visualization import SkyViewFactor

   result = SkyViewFactor(max_radius=10, include_positive_openness=True).compute(
       dem, Resolution(0.5, 0.5)
   )
   svf = result.sky_view_factor
   opns = result.positive_openness

Reading and writing rasters
---------------------------

The GDAL helpers in ``rvt.default`` are replaced by :mod:`rvt.io`, which is
rasterio-based. :func:`rvt.io.read` returns a :class:`rvt.io.Raster`: a float32
array with nodata already converted to NaN, plus the transform and CRS, and the
cell size as the ``resolution`` property. :func:`rvt.io.write` writes an LZW
compressed, internally tiled GeoTIFF (float32 with NaN nodata, or ``uint8``
without a nodata tag).

=====================================  ================================
Version 2                              Version 3
=====================================  ================================
``rvt.default.get_raster_arr(path)``   :func:`rvt.io.read`
``rvt.default.get_raster_size(path)``  ``rasterio.open(path)`` directly
``rvt.default.save_raster(...)``       :func:`rvt.io.write`
=====================================  ================================

:func:`rvt.compute` combines the two: it applies a single-image config to a
:class:`rvt.io.Raster` and returns a raster with the same georeferencing.

Saving visualizations
---------------------

``DefaultValues`` — one mutable object carrying ~90 attributes, the per-output
``*_compute``/``*_save_float``/``*_save_8bit`` flags, the ``*_bytscl`` stretches
and the ``save_*`` methods — is gone. Its job is split between the configs
(which carry their own parameters, their ``EIGHT_BIT_STRETCH`` and their
parameter-encoded output names) and :class:`rvt.Session` (which runs a list of
:class:`rvt.SessionTask` and writes the files).

=================================================  ===============================================================================
Version 2                                          Version 3
=================================================  ===============================================================================
``DefaultValues()`` attributes                     fields of the individual configs
``DefaultValues.ve_factor``                        ``vertical_exaggeration`` per config, or ``Session(vertical_exaggeration=...)``
``DefaultValues.get_hillshade(dem_arr, ...)``      ``Hillshade(...).compute(dem, resolution)``
``DefaultValues.save_hillshade(dem_path)``         ``Session(dem_path, [SessionTask(Hillshade())]).run()``
``DefaultValues.save_visualizations(dem_path)``    one :class:`rvt.SessionTask` per visualization
``DefaultValues.*_save_float`` / ``*_save_8bit``   ``SessionTask(save_float=..., save_8bit=...)``
``DefaultValues.*_bytscl``                         ``Config.EIGHT_BIT_STRETCH``
``DefaultValues.float_to_8bit(...)``               :func:`rvt.stretch.to_8bit`
``DefaultValues.get_<vis>_file_name(dem_path)``    ``config.outputs(dem_name)``, or ``config.output_stem(dem_name)``
``DefaultValues.overwrite``                        ``Session(overwrite=...)``
``DefaultValues.tile_size_limit`` / ``tile_size``  ``Session(automatic_tiling_threshold=..., tile_size=...)``
``DefaultValues.create_log_file(...)``             *not ported*
=================================================  ===============================================================================

Blending
--------

Blends are still layer stacks rendered bottom-up, but the layers are typed
objects and the presets are functions returning inspectable combinations. See
:doc:`blending`.

================================================  ============================================
Version 2                                         Version 3
================================================  ============================================
``rvt.blend.BlenderLayer``                        :class:`rvt.blending.BlendLayer`
``rvt.blend.BlenderCombination``                  :class:`rvt.blending.BlendCombination`
``BlenderCombination.render_all_images(...)``     ``BlendCombination.render(dem, resolution)``
``rvt.blend.BlenderCombinations``                 *not ported* (keep your own list)
``rvt.blend.TerrainSettings``                     *not ported* (build the configs you need)
``rvt.blend_func.normalize_image(...)``           :class:`rvt.blending.Normalization`
``rvt.blend_func.blend_images(...)``              :class:`rvt.blending.BlendMode`
``rvt.blend_func.gray_scale_to_color_ramp(...)``  :class:`rvt.blending.ColorRamp`
``settings/blender_VAT.json``                     :func:`rvt.blending.presets.vat`
``rvt.blend.color_relief_image_map(...)``         :func:`rvt.blending.presets.crim`
``rvt.blend.e3mstp(...)``                         :func:`rvt.blending.presets.e3mstp`
``rvt.blend.e4mstp(...)``                         :func:`rvt.blending.presets.e4mstp`
================================================  ============================================

The visualization name strings that identified a layer's content
(``"sky-view factor"``, ``"openness - negative"``, ...) are replaced by typed
layer sources; the blender's internal three-azimuth
``"multiple directions hillshade"`` RGB image is now the explicitly named
:class:`rvt.blending.MultiDirectionalHillshadeRgbSource`.

Configuration files
-------------------

Version 2 read parameters and blends from JSON (and legacy ``.txt``) files:
``default_settings.json``, ``default_blender_combinations.json``,
``default_terrains_settings.json``, ``blender_*.json``. All of them are gone —
sessions and combinations are configured in code.

The one file format that survives is a new one:
:func:`rvt.blending.combination_to_json` and
:func:`rvt.blending.combination_from_json` round-trip a combination through the
``rvt-blend/1`` schema. **There is no reader for the version 2
``blender_*.json`` files**; old blend definitions have to be rebuilt as code (or
built once and re-serialized to the new schema).

Deliberate behaviour changes
----------------------------

*   **Nodata is NaN in the compute layer.** The ``no_data`` parameter is gone
    from every computation: :func:`rvt.io.read` converts a file's nodata value
    to NaN, and configs consume float32 arrays with NaN nodata and never mutate
    the caller's array (version 2 wrote NaN into the array it was passed).
*   **Validation actually happens.** Configs validate in ``__post_init__`` and
    raise ``ValueError``; version 2's opacity check was a no-op, and several
    tiling guards never raised.
*   **Tiled output equals untiled output.** Tiling is refused for
    visualizations whose result depends on the whole array — sky illumination
    (data-dependent overcast scaling and pyramid padding) and shadow/horizon
    (pyramid padding) cannot be tiled — instead of producing seams. The
    multi-scale relief model's tiling overlap is derived from the cell size,
    fixing version 2's treatment of ``feature_max`` in map units as a pixel
    count. See :doc:`large_rasters`.
*   **The negative-openness noise suffix is ``_NRhigh``**, consistently with
    the other sky-view visualizations; version 2 wrote ``_NRstrong`` for that
    one output.
*   **GPU support is dropped.** The experimental ``cupy`` code paths are not
    part of version 3.
*   **The 8-bit multi-hillshade is not implemented in a session yet.** Version 2
    rendered a special three-band RGB image from the DEM (azimuths 315, 22.5 and
    90) rather than stretching the float stack; a
    :class:`rvt.SessionTask` requesting it raises ``NotImplementedError``. The
    float stack saves normally, and the RGB image is available for blends as
    :class:`rvt.blending.MultiDirectionalHillshadeRgbSource`.
*   **Blend rendering is pure.** The in-place mutation of inputs is gone. The
    version 2 formulas are kept as they are, including the quirk that overlay
    and soft light effectively ignore the layer opacity, so results stay
    comparable.
*   **Sky-view-factor-bearing blends are not bit-identical.** Version 2
    rendered from float64 sky-view factor arrays, version 3's compute layer is
    float32; parity for those blends is verified at a relative tolerance of
    1e-7 rather than exactly.
*   **Python 3.12 and rasterio.** GDAL's Python bindings are no longer a
    dependency; see :doc:`installation`.
