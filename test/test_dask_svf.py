from pathlib import Path
import rioxarray
import rvt.vis
import rvt.vis_dask
import rvt.default
import numpy as np
import dask.array as da
import xarray as xr
from nptyping import NDArray

# pytest test dask and numpy array equality

input_dem_path = Path(r"test_data/TM1_564_146.tif")
# default_values = rvt.default.DefaultValues()
CHUNKSIZE = {'x': 100, 'y':100}
input_arr: xr.DataArray = rioxarray.open_rasterio(input_dem_path, chunks = CHUNKSIZE, cache = False, lock = False)
comp_svf, comp_asvf, comp_opns = [False, True, False]

def get_np_result():
    input_np_arr = np.array(input_arr.data[0]) #numpy array 2D
    x_res = abs(input_arr.rio.resolution()[0])
    y_res = abs(input_arr.rio.resolution()[1])
    no_data = input_arr.rio.nodata
    ve_factor = 1
    dict_svf_asvf_opns = rvt.vis.sky_view_factor(dem=input_np_arr, resolution=x_res, compute_svf=comp_svf,
                                                compute_asvf=comp_asvf, compute_opns=comp_opns,
                                                svf_n_dir=16, svf_r_max=10,
                                                svf_noise=0, asvf_dir=315,
                                                asvf_level=1, ve_factor=ve_factor,
                                                no_data = no_data)
    arr_svf, = dict_svf_asvf_opns.values() 
    return arr_svf

def get_dask_result(): 
    input_da_arr = input_arr.data[0] #dask array 2D
    x_res = abs(input_arr.rio.resolution()[0])
    y_res = abs(input_arr.rio.resolution()[1])
    no_data = input_arr.rio.nodata 
    ve_factor = 1
    arr_svf = rvt.vis_dask.dask_sky_view_factor(input_dem=input_da_arr, resolution=x_res,compute_svf=comp_svf,
                                                compute_asvf=comp_asvf, compute_opns=comp_opns,
                                                svf_n_dir=16, svf_r_max=10,
                                                svf_noise=0, asvf_dir=315,
                                                asvf_level=1, ve_factor=ve_factor,
                                                no_data = no_data)
    return arr_svf

numpy_result = get_np_result()
dask_to_be_computed = get_dask_result()
dask_result = dask_to_be_computed.compute()


def test_general_io_checks():
    assert isinstance(dask_result, type(input_arr.data.compute()))
    assert dask_result.dtype == np.float32
    # check shape and other attributes 2D -> 3D
    assert input_arr.shape[1:]== dask_to_be_computed.shape
    assert input_arr.data.chunksize[1:] == dask_to_be_computed.chunksize

def test_arr_edges_svf():
    np_svf_edges = [numpy_result[0, :],
                    numpy_result[-1, :],
                    numpy_result[:, 0],
                    numpy_result[:, -1]]
    da_svf_edges = [dask_result[ 0,:] ,
                    dask_result[-1,:],
                    dask_result[:,0] ,
                    dask_result[:,-1]]
    np.testing.assert_array_equal(np_svf_edges[:], da_svf_edges[:])

def test_arr_inner():
    np_svf_inner = numpy_result[1:-1, 1:-1]
    da_svf_inner = dask_result[ 1:-1, 1:-1]
    np.testing.assert_array_equal(np_svf_inner, da_svf_inner)
    assert np_svf_inner.shape  == da_svf_inner.shape

