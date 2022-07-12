from pathlib import Path
import rioxarray
import rvt.vis
import rvt.blend_func
import rvt.blend_func_dask
import rvt.default
import numpy as np
import dask.array as da
import xarray as xr
from nptyping import NDArray
import dask
import copy
import pytest

# Test dask and numpy array equality. Do not run tests on very large rasters. Data must fit in memory.

input_dem_path = Path(r"test_data/TM1_564_146.tif")
# input_dem_path = Path(r"test_data/synthetic_dem15_0.50.tif")
CHUNKSIZE = {'x': 100, 'y': 100}

def get_source_img(src_path):
    """Load test dem and reshape to 3D"""
    input_arr: xr.DataArray = rioxarray.open_rasterio(src_path, chunks = CHUNKSIZE, cache = False, lock = False).data.squeeze() 
    stckd_arr = da.stack((input_arr, input_arr, input_arr))
    np_stck = np.array(stckd_arr)
    da_stck = da.from_array(np_stck, chunks = {0: -1, 1: CHUNKSIZE['y'], 2: CHUNKSIZE['x']})
    return da_stck, np_stck

def get_artifical_img(src_path):
    """Generated 2D data same of the same x and y shape as loaded test dem."""
    input_arr: xr.DataArray = rioxarray.open_rasterio(src_path, chunks = CHUNKSIZE, cache = False, lock = False).data.squeeze() 
    dimx, dimy = input_arr.shape
    np_2d = np.array([[i for i in range(dimy)] for j in range(dimx)]).astype('float32')
    da_2d = da.from_array(np_2d, chunks =(CHUNKSIZE['y'], CHUNKSIZE['x']))
    return da_2d, np_2d


@pytest.mark.parametrize("norm", ["Value", "Percent", None])
@pytest.mark.parametrize("minn, maxn", [ (0.2, 0.7), (0, 1), (63, 98)])
def test_normalize_eq(norm, minn, maxn):
    da_stacked, np_stacked = get_source_img(input_dem_path)
    da_arr_3d = rvt.blend_func_dask.dask_normalize_image(image = da_stacked, visualization= "Sky-View Factor", min_norm= minn, 
                                            max_norm = maxn, normalization = norm).compute()
    np_arr_3d = rvt.blend_func.normalize_image(image = np_stacked, visualization= "Sky-View Factor", min_norm= minn, 
                                            max_norm = maxn, normalization = norm)

    np.testing.assert_array_equal(da_arr_3d, np_arr_3d)


@pytest.mark.parametrize("blend", ["Luminosity", "Multiply", "Normal", "Overlay", "Screen", "Soft_light"])
@pytest.mark.parametrize("minc, maxc", [(None, None), (0, 15), (0.7, 1), (63, 98)])
def test_blend_eq(blend, minc, maxc):
    da_stacked, np_stacked = get_source_img(input_dem_path)
    da_2d, np_2d = get_artifical_img(input_dem_path)
    da_arr_3d = rvt.blend_func_dask.dask_blend_images(active = da_2d, background = da_stacked, blend_mode = blend,
                                                    min_c = minc, max_c = maxc).compute()
    np_arr_3d = rvt.blend_func.blend_images(active = np_2d, background = np_stacked, blend_mode = blend,
                                                    min_c = minc, max_c = maxc) 

    np.testing.assert_array_equal(da_arr_3d, np_arr_3d)

@pytest.mark.parametrize("opac", [25, 75, 0, 100, -5])
def test_render_eq(opac):
    ## if active and background are both 2d, unexpected results
    da_stacked, np_stacked = get_source_img(input_dem_path)
    da_2d, np_2d = get_artifical_img(input_dem_path)
    da_arr_3d = rvt.blend_func_dask.dask_render_images(active = da_2d, background = da_stacked, opacity = opac).compute()
    np_arr_3d = rvt.blend_func.render_images(active = np_2d, background = np_stacked, opacity = opac)  

    if len(np_arr_3d.shape) ==2:
        #np.testing.assert_array_equal(da_arr_3d[:,0:99], np_arr_3d[:,0:99])
        #np.testing.assert_array_equal(da_arr_3d[:,200:], np_arr_3d[:,200:])
        # np.testing.assert_array_equal(da_arr_3d[:,99:199], np_arr_3d[:,99:199])
        np.testing.assert_array_equal(da_arr_3d, np_arr_3d)
    else: 
        np.testing.assert_array_equal(da_arr_3d, np_arr_3d)

      
