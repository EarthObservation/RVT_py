from pathlib import Path
import rioxarray
import rvt.vis
import rvt.vis_dask
import rvt.default
import numpy as np
import dask.array as da
import xarray as xr
from nptyping import NDArray
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

    
@pytest.mark.parametrize("c_min, c_max", [(0, 1), (30, 60)])
@pytest.mark.parametrize("high, low", [(200, 0), (100, 50)])
def test_vis_byte_scale(c_min, c_max, high, low):
    da_stck, np_stck = get_source_img(input_dem_path)
    da_2d = copy.deepcopy(da_stck[0])
    np_2d = copy.deepcopy(np_stck[0])    

    da_arr_2d = rvt.vis_dask.dask_byte_scale(data = da_2d, c_min = c_min, c_max = c_max, high=high, low = low).compute()
    np_arr_2d = rvt.vis.byte_scale(data = np_2d, c_min = c_min, c_max = c_max, high=high, low = low)  

    np.testing.assert_array_equal(da_arr_2d, np_arr_2d)

@pytest.mark.parametrize("high, low", [(200, 0), (100, 50), (255, 1)])
def test_vis_byte_scale_None_case(high, low):
    da_stck, np_stck = get_source_img(input_dem_path)
    da_2d = copy.deepcopy(da_stck[0])
    np_2d = copy.deepcopy(np_stck[0])    

    da_arr_2d_NoneNone = rvt.vis_dask.dask_byte_scale(data = da_2d, c_min = None, c_max = None, high=high, low = low).compute()
    np_arr_2d_NoneNone = rvt.vis.byte_scale(data = np_2d, c_min = np.nanmin(np_2d), c_max = np.nanmax(np_2d), high=high, low = low)  

    # c_min = np.nanmin(np_2d)/c_max = np.nanmax(np_2d) are given because "if c_min/max == None" condition is commented out in rvt.vis.byte_scale
    da_arr_2d_c_maxNone = rvt.vis_dask.dask_byte_scale(data = da_2d, c_min = None, c_max = 500, high=high, low = low).compute()
    np_arr_2d_c_maxNone= rvt.vis.byte_scale(data = np_2d, c_min = np.nanmin(np_2d), c_max = 500, high=high, low = low)

    da_arr_2d_c_minNone = rvt.vis_dask.dask_byte_scale(data = da_2d, c_min = 1, c_max = None, high=high, low = low).compute()
    np_arr_2d_c_minNone= rvt.vis.byte_scale(data = np_2d, c_min = 1, c_max = np.nanmax(np_2d), high=high, low = low)  

    np.testing.assert_array_equal(da_arr_2d_NoneNone, np_arr_2d_NoneNone)
    np.testing.assert_array_equal(da_arr_2d_c_maxNone, np_arr_2d_c_maxNone)
    np.testing.assert_array_equal(da_arr_2d_c_minNone, np_arr_2d_c_minNone)