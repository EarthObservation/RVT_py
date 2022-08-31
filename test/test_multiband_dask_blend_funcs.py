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
    numpy_stck = np.array(stckd_arr)
    dask_stck = da.from_array(numpy_stck, chunks = {0: -1, 1: CHUNKSIZE['y'], 2: CHUNKSIZE['x']})
    return dask_stck, numpy_stck

def get_artifical_img(src_path):
    """Generated 2D data same of the same x and y shape as loaded test dem."""
    input_arr: xr.DataArray = rioxarray.open_rasterio(src_path, chunks = CHUNKSIZE, cache = False, lock = False).data.squeeze() 
    dimx, dimy = input_arr.shape
    numpy_2d = np.array([[i for i in range(dimy)] for j in range(dimx)]).astype('float32')
    dask_2d = da.from_array(numpy_2d, chunks =(CHUNKSIZE['y'], CHUNKSIZE['x']))
    return dask_2d, numpy_2d


@pytest.mark.parametrize("vis", ["Sky-View Factor", "Slope gradient", "mhs",None])
@pytest.mark.parametrize("norm", ["Value", "Perc", None])
@pytest.mark.parametrize("minn, maxn", [ (0.2, 0.7), (0, 1), (63, 98)])
def test_normalize_eq(vis, norm, minn, maxn):
    da_stck, np_stck = get_source_img(input_dem_path)
    da_stacked = copy.deepcopy(da_stck)
    np_stacked = copy.deepcopy(np_stck)
   
    da_arr_3d = rvt.blend_func_dask.dask_normalize_image(image = da_stacked, visualization= vis, min_norm= minn, 
                                            max_norm = maxn, normalization = norm).compute()
    np_arr_3d = rvt.blend_func.normalize_image(image = np_stacked, visualization= vis, min_norm= minn, 
                                            max_norm = maxn, normalization = norm)
    np.testing.assert_array_equal(da_arr_3d, np_arr_3d)


@pytest.mark.parametrize("blend", ["Luminosity", "Multiply", "Normal", "Overlay", "Screen", "Soft_light"])
@pytest.mark.parametrize("minc, maxc", [(None, None), (0, 15), (0.7, 1), (63, 98)])
def test_blend_eq(blend, minc, maxc):
    da_stck, np_stck = get_source_img(input_dem_path)
    da_2dim, np_2dim = get_artifical_img(input_dem_path)
    da_stacked = copy.deepcopy(da_stck)
    np_stacked = copy.deepcopy(np_stck)
    da_2d = copy.deepcopy(da_2dim)
    np_2d = copy.deepcopy(np_2dim)    

    da_arr_3d = rvt.blend_func_dask.dask_blend_images(active = da_2d, background = da_stacked, blend_mode = blend,
                                                    min_c = minc, max_c = maxc).compute()
    np_arr_3d = rvt.blend_func.blend_images(active = np_2d, background = np_stacked, blend_mode = blend,
                                                    min_c = minc, max_c = maxc) 

    np.testing.assert_array_equal(da_arr_3d, np_arr_3d)

@pytest.mark.parametrize("opac", [25, 75, 0, 100, -5, 0.5])
def test_render_eq(opac):
    da_stck, np_stck = get_source_img(input_dem_path)
    da_2dim, np_2dim = get_artifical_img(input_dem_path)
    da_stacked = copy.deepcopy(da_stck)
    np_stacked = copy.deepcopy(np_stck)
    da_2d = copy.deepcopy(da_2dim)
    np_2d = copy.deepcopy(np_2dim)    

    da_arr_3d = rvt.blend_func_dask.dask_render_images(active = da_2d, background = da_stacked, opacity = opac).compute()
    np_arr_3d = rvt.blend_func.render_images(active = np_2d, background = np_stacked, opacity = opac, bkg_abs_max = np.nanmax(np_stacked), bkg_abs_min = np.nanmin(np_stacked), act_abs_max = np.nanmax(np_2d), act_abs_min = np.nanmin(np_2d))  

    np.testing.assert_array_equal(da_arr_3d, np_arr_3d)
