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
import pytest

# Test dask and numpy array equality. Do not run tests on very large rasters. Data must fit in memory.

input_dem_path = Path(r"test_data/TM1_564_146.tif")
CHUNKSIZE = {'x': 300, 'y': 200}
## Load test dem and reshape to 3D. 
input_arr: xr.DataArray = rioxarray.open_rasterio(input_dem_path, chunks = CHUNKSIZE, cache = False, lock = False).data.squeeze() 
stacked_arr = da.stack((input_arr, input_arr, input_arr))
np_stacked = np.array(stacked_arr)
da_stacked = da.from_array(np_stacked, chunks = {0: -1, 1: CHUNKSIZE['y'], 2: CHUNKSIZE['x']})
##Random 2D data
dimx, dimy = input_arr.shape
np_2d = np.array([[i for i in range(dimy)] for j in range(dimx)]).astype('float32')
da_2d = da.from_array(np_2d, chunks =(CHUNKSIZE['y'], CHUNKSIZE['x']))

@pytest.mark.parametrize("norm", ["Value", "Percent", None])
@pytest.mark.parametrize("minn, maxn", [ (0.2, 0.7), (0, 1), (63, 98)])
def test_normalize_eq(norm, minn, maxn):
    da_arr_3d = rvt.blend_func_dask.dask_normalize_image(image = da_stacked, visualization= "Sky-View Factor", min_norm= minn, 
                                            max_norm = maxn, normalization = norm).compute()
    np_arr_3d = rvt.blend_func.normalize_image(image = np_stacked, visualization= "Sky-View Factor", min_norm= minn, 
                                            max_norm = maxn, normalization = norm)

    assert len(np_arr_3d.shape) == 3
    assert len(da_arr_3d.shape) == 3

    np_3d_inner = np_arr_3d[:, 1:-1, 1:-1]
    da_3d_inner = da_arr_3d[:, 1:-1, 1:-1]
    np_3d_edges = [np_arr_3d[:, 0, :],
                     np_arr_3d[:, -1, :],
                     np_arr_3d[:, :, 0],
                     np_arr_3d[:, :, -1]]
    da_3d_edges = [da_arr_3d[:, 0, :],
                     da_arr_3d[:, -1, :],
                     da_arr_3d[:, :, 0],
                     da_arr_3d[:, :, -1]]

    np.testing.assert_array_equal(da_3d_inner, np_3d_inner)
    for i in range(len(da_3d_edges)):
        np.testing.assert_array_equal(da_3d_edges[i][0], np_3d_edges[i][0])
        np.testing.assert_array_equal(da_3d_edges[i][1], np_3d_edges[i][1])


@pytest.mark.parametrize("blend", ["Luminosity", "Multiply", "Normal", "Overlay", "Screen", "Soft_light"])
@pytest.mark.parametrize("minc, maxc", [(None, None), (0, 15), (0.7, 1), (63, 98)])
def test_blend_eq(blend, minc, maxc):
    da_arr_3d = rvt.blend_func_dask.dask_blend_images(active = da_stacked, background = da_2d, blend_mode = blend,
                                                    min_c = minc, max_c = maxc).compute()
    np_arr_3d = rvt.blend_func.blend_images(active = np_stacked, background = np_2d, blend_mode = blend,
                                                    min_c = minc, max_c = maxc)  
    # further parametrize for different shapes of active and background                                              
    np_3d_inner = np_arr_3d[:, 1:-1, 1:-1]
    da_3d_inner = da_arr_3d[:, 1:-1, 1:-1]
    np_3d_edges = [np_arr_3d[:, 0, :],
                     np_arr_3d[:, -1, :],
                     np_arr_3d[:, :, 0],
                     np_arr_3d[:, :, -1]]
    da_3d_edges = [da_arr_3d[:, 0, :],
                     da_arr_3d[:, -1, :],
                     da_arr_3d[:, :, 0],
                     da_arr_3d[:, :, -1]]

    np.testing.assert_array_equal(da_3d_inner, np_3d_inner)
    for i in range(len(da_3d_edges)):
        np.testing.assert_array_equal(da_3d_edges[i][0], np_3d_edges[i][0])
        np.testing.assert_array_equal(da_3d_edges[i][1], np_3d_edges[i][1])


@pytest.mark.parametrize("opac", [25, 75, 0, 100, -5])
def test_render_eq(opac):
    da_arr_3d = rvt.blend_func_dask.dask_render_images(active = da_stacked, background = da_2d, opacity = opac).compute()
    np_arr_3d = rvt.blend_func.render_images(active = np_stacked, background = np_2d, opacity = opac)  
    # further parametrize for different shapes of active and background  
    
    np_3d_inner = np_arr_3d[:, 1:-1, 1:-1]
    da_3d_inner = da_arr_3d[:, 1:-1, 1:-1]
    np_3d_edges = [np_arr_3d[:, 0, :],
                     np_arr_3d[:, -1, :],
                     np_arr_3d[:, :, 0],
                     np_arr_3d[:, :, -1]]
    da_3d_edges = [da_arr_3d[:, 0, :],
                     da_arr_3d[:, -1, :],
                     da_arr_3d[:, :, 0],
                     da_arr_3d[:, :, -1]]

    np.testing.assert_array_equal(da_3d_inner, np_3d_inner)
    for i in range(len(da_3d_edges)):
        np.testing.assert_array_equal(da_3d_edges[i][0], np_3d_edges[i][0])
        np.testing.assert_array_equal(da_3d_edges[i][1], np_3d_edges[i][1])
      
