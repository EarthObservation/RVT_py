# TODO

This is the markdown todo file for feature/dask branch of Github Repo [RVT_py](https://github.com/EarthObservation/RVT_py).

### Content

- [ ] Modify rvt to support parallel processing and streaming computation on rasters that don't fit into memory. 
- [ ] Dask analysis workflow. Generate task graph: 
  - [ ] Load raster. 
    - [ ] Apply visualization. 
    - [ ] Apply normalization. 
    - [ ] Apply blending. 
    - [ ] Apply rendering. 
    - [ ] ... repeat ...
  - [ ] Save raster. 
  - [ ] Execute tasks by calling .compute() at the end.
- [ ] Mapping and chaining of these functions across all dask blocks is done in a same fashion as described in the third section of  [napari tutorials](https://napari.org/tutorials/processing/dask.html). Multiple cycles of visualisation -> normalization -> blending -> rendering restults in long tasks, taking one chunk "from start to finish".

### Todo

- [ ] Additional code refactoring. 
- [ ] Dask memory issues.
  - [ ] Get available memory and compute / set memory (GB) per Dask worker. 
  - [ ] Get / compute optimal chunk size. 
- [ ] Read multiband data. 

### In Progress

- [ ] Work on restoring some of the previously existing non-dask functionality - keep both options.
- [ ] Fix numpy runtime errors / suppress warnings.
- [ ] Additional testing of visualization functions.

### Done ✓

- [x] Wrap exising vis functions and map over dask array with some overlap. 
- [x] Wrap exising blend_func functions and map over dask array. 
- [x] Read (lazy load) raster in chunks.
- [x] Save raster in chunk by chunk (.tif and .zarr). Parallel writes. 