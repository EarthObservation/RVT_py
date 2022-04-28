# TODO

This is the markdown todo file for feature/dask branch of Github Repo [RVT_py](https://github.com/EarthObservation/RVT_py).

### Content

* Modify rvt to support parallel processing and streaming computation on rasters that don't fit into memory. 
* Heavy use of dask functions `map_overlap` and `map_blocks`.
* Dask analysis workflow. Generate task graph: 
  * Load raster. 
    * Apply visualization. 
    * Apply normalization. 
    * Apply blending. 
    * Apply rendering. 
    * ... repeat ...
  * Save raster. 
  * Execute tasks by calling .compute() at the end.

* Mapping and chaining of these functions across all dask blocks is done in a same fashion as described in the third section of  [napari tutorials](https://napari.org/tutorials/processing/dask.html). Multiple cycles of visualisation -> normalization -> blending -> rendering restults in long tasks, taking one chunk "from start to finish".\*  

\*With [large inputs](https://github.com/dask/dask/issues/3514) takes a long time to start calculation (Is there a more efficient way of saying "this input chunk should map to this output chunk"?)

### Todo

- [ ] Additional code refactoring. 
- [ ] Dask memory issues.
  - [ ] Get available memory and compute / set memory (GB) per Dask worker. 
  - [ ] Get / compute optimal chunk size. 

##### _Notes on Todo_
- If chunks are too small: huge amount of tasks and a lot of time spent developing the task graph - task scheduler hangs
- If chunks are too big:  _"Unable to allocate XX MiB for an array with shape (dimx, dimy) and data type float32."_ [Error encountered.](https://stackoverflow.com/questions/62839068/memoryerror-unable-to-allocate-mib-for-an-array-with-shape-and-data-type-when) 
- Raster datasets and are usually stored in blocks (tiles). It is better to  make dask chunks **N * original tile dimensions** to avoid bringing up more data than is needed each time the data is accesed. 
- If overlap depth is greater than any chunk along a particular axis, then the array is rechunked -> Strange behavior, error at blending step.
- Runtime Warnings encountered:
  - `RuntimeWarning: All-NaN slice encountered if np.nanmin(image_chunk) < 0 or np.nanmax(image_chunk) > 1`
  - `RuntimeWarning: overflow encountered in long_scalars maxsize = math.ceil(nbytes / (other_numel * itemsize))`

### In Progress

- [ ] Work on restoring some of the previously existing non-dask functionality - keep both options.
- [ ] Fix numpy runtime errors / suppress warnings.
- [ ] Additional testing of visualization functions.

### Done ✓

- [x] Wrap exising vis functions and map over dask array with some overlap. 
- [x] Wrap exising blend_func functions and map over dask array. 
- [x] Read (lazy load) raster in chunks.
- [x] Save raster in chunk by chunk (.tif and .zarr). Parallel writes. 
- [x] Read multiband data. 