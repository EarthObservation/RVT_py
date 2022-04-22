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
import pytest

# Test dask and numpy array equality. Do not run tests on very large rasters. Data must fit in memory.
# TEST DATA : input 2 dems, test "normalize_image", "blend_images",  "render_images" and gray_scale_to_color_ramp functions with different parameters

input_dem_path = Path(r"test_data/TM1_564_146.tif")
CHUNKSIZE = {'x': 150, 'y':150}
## first input dem
input_arr_1: xr.DataArray = rioxarray.open_rasterio(input_dem_path, chunks = CHUNKSIZE, cache = False, lock = False) 

def get_dask_result() -> da.Array: 
    input_da_arr = input_arr_1.data[0] #dask array 2D
    x_res = abs(input_arr_1.rio.resolution()[0])
    y_res = abs(input_arr_1.rio.resolution()[1])
    no_data = input_arr_1.rio.nodata 
    ve_factor = 1
    arr_svf = rvt.vis_dask.dask_sky_view_factor(input_dem=input_da_arr, resolution=x_res,compute_svf=True,
                                                compute_asvf=False, compute_opns=False,
                                                svf_n_dir=2, svf_r_max=10,
                                                svf_noise=0, asvf_dir=315,
                                                asvf_level=1, ve_factor=ve_factor,
                                                no_data = no_data)
    return arr_svf

## second input dem
da_input_arr = get_dask_result()
np_input_arr = get_dask_result().compute()


@pytest.mark.parametrize("norm", ["Value", "Percent", None])
@pytest.mark.parametrize("minn, maxn", [ (0.2, 0.7), (0, 1)])
def test_normalize_eq(norm, minn, maxn):
    np_arr = rvt.blend_func.normalize_image(image = np_input_arr, visualization= "Sky-View Factor", min_norm=minn, 
                                            max_norm =  maxn, normalization = norm)
    da_arr = rvt.blend_func_dask.dask_normalize_image(image = da_input_arr, visualization= "Sky-View Factor", min_norm= minn, 
                                            max_norm = maxn, normalization = norm).compute()
    da_inner = da_arr[ 1:-1, 1:-1]
    np_inner = np_arr[ 1:-1, 1:-1]
    da_edges = [da_arr[0, :],
                da_arr[-1, :],
                da_arr[:, 0],
                da_arr[:, -1]]
    np_edges = [np_arr[0, :],
                np_arr[-1, :],
                np_arr[:, 0],
                np_arr[:, -1]]
    np.testing.assert_array_equal(np_inner, da_inner)
    np.testing.assert_array_equal(np_edges[0], da_edges[0])
    np.testing.assert_array_equal(np_edges[1], da_edges[1])
    np.testing.assert_array_equal(np_edges[2], da_edges[2])
    np.testing.assert_array_equal(np_edges[3], da_edges[3])


@pytest.mark.parametrize("blend", ["Luminosity", "Multiply", "Normal", "Overlay", "Screen", "Soft_light"])
@pytest.mark.parametrize("minc, maxc", [(None, None), (0, 15), (0.7, 1)])
def test_blend_eq(blend, minc, maxc):
    np_arr = rvt.blend_func.blend_images(active = np_input_arr , background = np.array(input_arr_1.data[0]), blend_mode = blend,
                                         min_c =  minc, max_c = maxc )
    da_arr = rvt.blend_func_dask.dask_blend_images(active = da_input_arr, background = input_arr_1.data[0], blend_mode = blend,
                                                    min_c = minc, max_c = maxc).compute()
    da_inner = da_arr[ 1:-1, 1:-1]
    np_inner = np_arr[ 1:-1, 1:-1]
    da_edges = [da_arr[0, :],
                da_arr[-1, :],
                da_arr[:, 0],
                da_arr[:, -1]]
    np_edges = [np_arr[0, :],
                np_arr[-1, :],
                np_arr[:, 0],
                np_arr[:, -1]]
    np.testing.assert_array_equal(np_inner, da_inner)
    np.testing.assert_array_equal(np_edges[0], da_edges[0])
    np.testing.assert_array_equal(np_edges[1], da_edges[1])
    np.testing.assert_array_equal(np_edges[2], da_edges[2])
    np.testing.assert_array_equal(np_edges[3], da_edges[3])


@pytest.mark.parametrize("opac", [25, 75, 0, 100, -5])
def test_render_eq(opac):
    np_arr = rvt.blend_func.render_images(active = np_input_arr , background = np.array(input_arr_1.data[0]), opacity = opac )
    da_arr = rvt.blend_func_dask.dask_render_images(active = da_input_arr, background = input_arr_1.data[0], opacity = opac).compute()
    da_inner = da_arr[ 1:-1, 1:-1]
    np_inner = np_arr[ 1:-1, 1:-1]
    da_edges = [da_arr[0, :],
                da_arr[-1, :],
                da_arr[:, 0],
                da_arr[:, -1]]
    np_edges = [np_arr[0, :],
                np_arr[-1, :],
                np_arr[:, 0],
                np_arr[:, -1]]
    np.testing.assert_array_equal(np_inner, da_inner)
    np.testing.assert_array_equal(np_edges[0], da_edges[0])
    np.testing.assert_array_equal(np_edges[1], da_edges[1])
    np.testing.assert_array_equal(np_edges[2], da_edges[2])
    np.testing.assert_array_equal(np_edges[3], da_edges[3])
  

@pytest.mark.parametrize("cmap, min_cmap_cut, max_cmap_cut", [("OrRd", 0.2, 1), ("Blues", 0, 0.7)])
@pytest.mark.parametrize("alph", [False, True])
@pytest.mark.parametrize("output_8", [False, True])
def test_gray_to_color(cmap, min_cmap_cut, max_cmap_cut, alph, output_8):
    np_arr = rvt.blend_func.gray_scale_to_color_ramp(gray_scale = np_input_arr, colormap = cmap, min_colormap_cut = min_cmap_cut, 
                                                    max_colormap_cut = max_cmap_cut, alpha = alph, output_8bit = output_8 )
    da_arr = rvt.blend_func_dask.dask_gray_scale_to_color_ramp(gray_scale = da_input_arr, colormap = cmap, min_colormap_cut = min_cmap_cut, 
                                                             max_colormap_cut = max_cmap_cut, alpha = alph, 
                                                             output_8bit = output_8).compute()
    if alph == True:
        assert da_arr.shape[0] == 4
        assert np_arr.shape[0] == 4
    else:
        assert da_arr.shape[0] == 3
        assert np_arr.shape[0] == 3
    
    da_inner = da_arr[:, 1:-1, 1:-1]
    np_inner = np_arr[:, 1:-1, 1:-1]
    da_edges = [da_arr[:,0, :],
                da_arr[:,-1, :],
                da_arr[:, :, 0],
                da_arr[:, :, -1]]
    np_edges = [np_arr[:, 0, :],
                np_arr[:, -1, :],
                np_arr[:, :, 0],
                np_arr[:, :, -1]]
    np.testing.assert_array_equal(np_inner, da_inner)
    np.testing.assert_array_equal(np_edges[0], da_edges[0])
    np.testing.assert_array_equal(np_edges[1], da_edges[1])
    np.testing.assert_array_equal(np_edges[2], da_edges[2])
    np.testing.assert_array_equal(np_edges[3], da_edges[3])
