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

NR_DIR = 16

def get_np_result(nr_directions=NR_DIR, sun_elevation=35):
    input_np_arr = np.array(input_arr.data[0]) #numpy array 2D
    x_res = abs(input_arr.rio.resolution()[0])
    y_res = abs(input_arr.rio.resolution()[1])
    no_data = input_arr.rio.nodata
    ve_factor = 1
    arr_mhs = rvt.vis.multi_hillshade(dem=input_np_arr, resolution_x=x_res,
                                        resolution_y=y_res, nr_directions=nr_directions,
                                        sun_elevation=sun_elevation)
    return arr_mhs

def get_dask_result(nr_directions=NR_DIR, sun_elevation=35): 
    input_da_arr = input_arr.data[0] #dask array 2D
    x_res = abs(input_arr.rio.resolution()[0])
    y_res = abs(input_arr.rio.resolution()[1])
    no_data = input_arr.rio.nodata 
    ve_factor = 1
    arr_mhs = rvt.vis_dask.dask_multi_hillshade(input_dem=input_da_arr, resolution_x=x_res,
                                        resolution_y=y_res, nr_directions=nr_directions,
                                        sun_elevation=sun_elevation, ve_factor = ve_factor)
    return arr_mhs

# numpy_result = get_np_result(output_units = 'radian')
# dask_to_be_computed = get_dask_result(output_units = 'radian')
numpy_result = get_np_result()
dask_to_be_computed = get_dask_result()
dask_result = dask_to_be_computed.compute()


def test_general_io_checks():
    assert isinstance(dask_result, type(input_arr.data.compute()))
    assert dask_result.dtype == np.float32
    # check shape and other attributes 2D -> 3D
    assert input_arr.shape[1:] == dask_to_be_computed.shape[1:]
    assert input_arr.data.chunksize[0] == dask_to_be_computed.chunksize[0] / NR_DIR
    assert input_arr.data.chunksize[1:] == dask_to_be_computed.chunksize[1:]

def test_arr_edges_mhs():
    np_mhs_edges = [numpy_result[:,0, :],
                    numpy_result[:,-1, :],
                    numpy_result[:,:, 0],
                    numpy_result[:,:, -1]]
    da_mhs_edges = [dask_result[:, 0,:] ,
                    dask_result[:,-1,:],
                    dask_result[:,:,0] ,
                    dask_result[:,:,-1]]
    ## not passing
    np.testing.assert_array_equal(np_mhs_edges[0], da_mhs_edges[0])
    np.testing.assert_array_equal(np_mhs_edges[1], da_mhs_edges[1])
    np.testing.assert_array_equal(np_mhs_edges[2], da_mhs_edges[2])
    np.testing.assert_array_equal(np_mhs_edges[3], da_mhs_edges[3])


def test_arr_inner():
    np_mhs_inner = numpy_result[:,1:-1, 1:-1]
    da_mhs_inner = dask_result[:, 1:-1, 1:-1]
    np.testing.assert_array_equal(np_mhs_inner, da_mhs_inner)
    assert np_mhs_inner.shape[0] == NR_DIR
    assert da_mhs_inner.shape[0] == NR_DIR
