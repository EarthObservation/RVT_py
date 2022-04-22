# TODO

This is the markdown todo file for feature/dask.

### Todo

- [ ] Fix some blend_func functions not passing tests (blend_overlay, blend_soft_light). 
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