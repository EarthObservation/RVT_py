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
# CHUNKSIZE = True
CHUNKSIZE = {'x': 100, 'y':100}
input_arr: xr.DataArray = rioxarray.open_rasterio(input_dem_path, chunks = CHUNKSIZE, cache = False, lock = False)
 

def get_np_result(output_units = 'degree'):
    input_np_arr = np.array(input_arr.data[0]) #numpy array 2D
    x_res = abs(input_arr.rio.resolution()[0])
    y_res = abs(input_arr.rio.resolution()[1])
    no_data = input_arr.rio.nodata
    ve_factor = 1
    np_result = rvt.vis.slope_aspect(dem=input_np_arr, resolution_x=x_res,
                                    resolution_y=y_res, ve_factor=ve_factor,
                                    output_units=output_units, no_data= no_data)
    return np_result

def get_dask_result(output_units = 'degree'): 
    input_da_arr = input_arr.data[0] #dask array 2D
    x_res = abs(input_arr.rio.resolution()[0])
    y_res = abs(input_arr.rio.resolution()[1])
    no_data = input_arr.rio.nodata 
    ve_factor = 1
    arr_slp_asp = rvt.vis_dask.dask_slope_aspect(input_dem=input_da_arr, resolution_x=x_res,
                                        resolution_y=y_res, ve_factor=ve_factor,
                                        output_units=output_units, no_data = no_data)
    return arr_slp_asp


# numpy_result = get_np_result(output_units = 'radian')
# dask_to_be_computed = get_dask_result(output_units = 'radian')
numpy_result = get_np_result()
dask_to_be_computed = get_dask_result()
dask_result = dask_to_be_computed.compute()

def test_general_io_checks():
    assert isinstance(dask_result, type(input_arr.data.compute()))
    # check shape and other attributes 2D -> 3D
    assert input_arr.shape[1:] == dask_to_be_computed.shape[1:]
    assert input_arr.data.chunksize[0] == dask_to_be_computed.chunksize[0] / 2
    assert input_arr.data.chunksize[1:] == dask_to_be_computed.chunksize[1:]


def test_arr_edges_slope():
    np_slope_edges = [numpy_result['slope'][0, :],
                      numpy_result['slope'][-1, :],
                      numpy_result['slope'][:, 0],
                      numpy_result['slope'][:, -1]]
    da_slope_edges = [dask_result[0, 0,:] ,
                      dask_result[0,-1,:],
                      dask_result[0,:,0] ,
                      dask_result[0,:,-1]]
    ##test nan - np passing test
    for np_slope_edge in np_slope_edges:
        np.testing.assert_array_equal(np_slope_edge, np.nan)
    # for da_slope_edge in da_slope_edges:
        # np.testing.assert_array_equal(da_slope_edge, np.nan)

    ## not passing
    np.testing.assert_array_equal(np_slope_edges[0], da_slope_edges[0])
    np.testing.assert_array_equal(np_slope_edges[1], da_slope_edges[1])
    np.testing.assert_array_equal(np_slope_edges[2], da_slope_edges[2])
    np.testing.assert_array_equal(np_slope_edges[3], da_slope_edges[3])

def test_arr_edges_aspect():
    np_aspect_edges = [numpy_result['aspect'][0, :],
                      numpy_result['aspect'][-1, :],
                      numpy_result['aspect'][:, 0],
                      numpy_result['aspect'][:, -1]]
    da_aspect_edges = [dask_result[1,0,:],
                      dask_result[1,-1,:],
                      dask_result[1,:,0] ,
                      dask_result[1,:,-1]]

    np.testing.assert_array_equal(np_aspect_edges[0], da_aspect_edges[0])
    np.testing.assert_array_equal(np_aspect_edges[1], da_aspect_edges[1])
    np.testing.assert_array_equal(np_aspect_edges[2], da_aspect_edges[2])
    np.testing.assert_array_equal(np_aspect_edges[3], da_aspect_edges[3])

def test_arr_inner():
    np_slope_inner = numpy_result['slope'][1:-1, 1:-1]
    np_aspect_inner = numpy_result['aspect'][1:-1, 1:-1]
    da_slope_inner = dask_result[0, 1:-1, 1:-1]
    da_aspect_inner = dask_result[1, 1:-1, 1:-1]

    np.testing.assert_array_equal(np_slope_inner, da_slope_inner)
    np.testing.assert_array_equal(np_aspect_inner, da_aspect_inner)
