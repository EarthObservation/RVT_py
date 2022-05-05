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

- [ ] Dask memory issues.
  - [ ] Get available memory and compute / set memory (GB) per Dask worker. 
  - [ ] Get / compute optimal chunk size. 
- [ ] Fix most outter padding in some visualization functions in original `vis.py`.

### In Progress

- [ ] Work on restoring some of the previously existing non-dask functionality - keep both options.
- [ ] Additional testing of visualization functions (for i/o functionality, writing and reading raster data.)
- [ ] Fix failing test with multiband (`rvt.blend_func.blend_images`) blend mode =  "Luminosity".

### Done ✓

- [x] Wrap exising vis functions and map over dask array with some overlap. 
- [x] Wrap exising blend_func functions and map over dask array. 
- [x] Read (lazy load) raster in chunks.
- [x] Save raster in chunk by chunk (.tif and .zarr). Parallel writes. 
- [x] Read multiband data. 
- [x] Fix numpy runtime errors / suppress warnings (keep, not really relevant).

##### _Notes on Todo and Done_
- If chunks are too small: huge amount of tasks and a lot of time spent developing the task graph - task scheduler hangs.
- If chunks are too big: `Exception has occurred: _ArrayMemoryError`, `Unable to allocate xx MiB for an array with shape (dimx, dimy) and data type float32.` [Error encountered.](https://stackoverflow.com/questions/62839068/memoryerror-unable-to-allocate-mib-for-an-array-with-shape-and-data-type-when). This happens only with very large rasters (e.g. 30 GB). For example, computation is successful on 2 GB raster with chunk size ranging from 8 MiB to 512 MiB. The same computation on 30 GB raster is successful only for chunk sizes 16 MiB to 64 MiB (without pausing the workers, which increases computation time significantly). Example of smooth and paused task stream below. Same computation, on the same 30 GB raster, 16 MiB chunks. The only difference is the number of workers and threads per worker. 1 worker & 12 threads in the case of "ok looking" task stream and 2 workers & 6 threads per worker in the case of workers pausing.

![Smooth task stream, computation takes less than 40 mins.](./docs/bmarks/not_pausing_w1_t12.png)

![Smooth task stream, computation takes more than 50 mins.](./docs/bmarks/pausing_w2_t6.png)
Possible issues: Is dask scheduler "too eager" and overcommiting tasks, is computation faster than writing to disk and chunks keep accumulating in memory (issue persists if writing to tif or to zarr), waiting to be written but eventually fail, memory leak, some expensive memory operation (e.g. rechunking in axis = 0 direction when visualizations produce 3D array) in the middle of workflow, is data load faster than it can be computed downstream,..? [Memory issues](https://alimanfoo.github.io/dask/2021/03/22/dask-memory-thought.html).
Try tweaking some settings in `distributed.yaml` file for smoother process:
1. `default-task-durations`,
2. set target worker  memory fractions to stay below: `target, spil, pause, terminate`,
3. set `MALLOC_TRIM_THRESHOLD_` to trim unmanaged memory (only works on linux os),
4. change `optimization.fuse.ave-width` (to fuse neighbouring tasks together).
- Raster datasets are usually stored in blocks (tiles). It is better to  make dask chunks **N * original tile dimensions** to avoid bringing up more data than is needed each time the data is accesed. If `chunks = True` when loading data with `rioxarray.open_rasterio ` automatic chunkig is done in a way that takes into account default tiling and **N** is determined in a way that each chunk size is 128MiB (_Dask default, better performance and higher success/failed to finish ratio was recorded with smaller chunk size, but it depends on the hardware and data size...scalability?_). Set chunk size `dask.config.set({"array.chunk-size": limit})`. Graphs below were generated as results of computations with dask destributed, running Local Cluster on a machine with total of 32 GB RAM (`psutil.virtual_memory().total = 32 GB`) and 16 cores (`sum(client.ncores().values()) = 16`). Aprox. 64 MiB equals chunk of size (4096, 4096) = 32 * (128, 128). For reference, "whole raster in memory, pure numpy  calcualtion" takes around 380 s.

![Comparision for chunk_size 8MiB, 16Mib, 32Mib, 64MiB, 128MiB, 256MiB, 512MiB. Calculation VAT_Combined, float32. Save tif to disk.](./docs/bmarks/csize_wt_ratio.png)
- Can't "replicate" entire graph above running `VAT_combined_(dask).py` on very large rasters (e.g. 30 GB). [Problems](https://github.com/dask/distributed/issues/2602) if chunks are bigger and/or number of workers is higher. See possible memory issues above. Apply some sort of queue to limit chunks or tasks being processed at the same time, for smoother processing/without pausing workers so much? Why is there such a difference in processing 2 GB vs 30 GB file with same chunksizes and worker/threads ratios?

![Comparision for chunk_size 16MiB, 32Mib, 64MiB. Calculation VAT_Combined, float32. Save tif to disk.](./docs/bmarks/csize_wt_ratio_large.png)
- If overlap depth is greater than any chunk along a particular axis, then the array is rechunked -> Strange behavior, may result in error at blending step.
- _Multi_hillshade_ (if saving directly from `vis.py`)results is very large (chunk depth and final) output of file size: _nr_directions * original raster file size_. Calculation in each direction is additional raster band.
- Metadata and the file naming convention when reading _.tif_ raster with `rioxarray.open_rasterio`. Understanding how the data is structured (e.g. assumes dims[1:] are named 'x' and 'y'. Name(s) of the bands (not always the same)? Indexing position always true (e.g. `dims[0]= 'band'`, `dims[1] = 'y'`, `dims[2] = 'x'`)?
- How are [no_data values](https://corteva.github.io/rioxarray/stable/getting_started/nodata_management.html) represented (no_data = data_set.rio.nodata, no_data = data_set.rio.encoded_nodata)? When saving to 8-bit `Error: Python int too large to convert to C long`. Temporary exclude copying attributes from orignal to saved dataset. Fix later. [Issue](https://github.com/corteva/rioxarray/issues/113).
- Runtime Warnings encountered:
  - `RuntimeWarning: All-NaN slice encountered if np.nanmin(image_chunk) < 0 or np.nanmax(image_chunk) > 1`
  - `RuntimeWarning: overflow encountered in long_scalars maxsize = math.ceil(nbytes / (other_numel * itemsize))`
  - `RuntimeWarning: overflow encountered in square tan_slope = np.sqrt(dzdx ** 2 + dzdy ** 2)`
- [Scaling](https://www.jennakwon.page/2021/03/benchmarks-dask-distributed-vs-ray-for.html).
- [Benchmarks](https://matthewrocklin.com/blog/work/2017/07/03/scaling).
- When reading/ writing very large arrays to .tif  sometimes error (random, can't reproduce): `RasterioIOError: file used by other process` followed by `distributed.worker - WARNING - Compute Failed`. Current workaround to avoid computation failure is to catch the exception in `rasterio.__init__.py", line 230, in open`, wait for 5 s and try to acquire the lock again. Change the lines 230 - 232 to:
```diff
- s = get_writer_for_path(path, driver=driver)(
-  path, mode, driver=driver, sharing=sharing, **kwargs
-  )
+ try:
+     s = get_writer_for_path(path, driver=driver)(
+     path, mode, driver=driver, sharing=sharing, **kwargs
+     )
+ except Exception as ex:
+     #print("!!EXCEPTION OCCURRED!!")
+     print(str(ex))
+     import time
+     time.sleep(5)
+     s = get_writer_for_path(path, driver=driver)(
+     path, mode, driver=driver, sharing=sharing, **kwargs
+     )
```