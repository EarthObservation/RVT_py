from pathlib import Path
import rvt.tile
import rvt.visualizations
import rvt.default
import numpy as np

# pytest rvt.tile.save_visualization_tile_by_tile

dem_path = Path(r"test_data\TM1_564_146.tif")
tile_size_x = 100
tile_size_y = 100


def test_slope_tile_by_tile() -> None:
    out_slope_path = Path(r"test_data\TM1_564_146_test_tile_slope.tif")
    # if resolution, resolution_x, resolution_y, no_data are None in function_parameters dict,
    # they are taken from dem (dem_path)
    output_units = "degree"
    function_parameters = {"resolution_x": None, "resolution_y": None, "no_data": None, "output_units": output_units}
    rvt.tile.save_visualization_tile_by_tile(
        visualization_function=rvt.vis.slope_aspect,
        function_parameters=function_parameters,
        dem_path=dem_path,
        overlap=1,
        tile_size_x=tile_size_x,
        tile_size_y=tile_size_y,
        out_raster_path=out_slope_path,
        out_raster_e_type=6,
        out_raster_nr_of_bands=1,
        out_visualization_dict_key="slope"
    )
    slope_tile_by_tile_arr = rvt.default.get_raster_arr(out_slope_path.as_posix())["array"]
    dem_arr_dict = rvt.default.get_raster_arr(dem_path.as_posix())
    slope_arr = rvt.vis.slope_aspect(
        dem=dem_arr_dict["array"],
        resolution_x=dem_arr_dict["resolution"][0],
        resolution_y=dem_arr_dict["resolution"][1],
        no_data=dem_arr_dict["no_data"],
        output_units=output_units
    )["slope"]
    assert np.array_equal(slope_tile_by_tile_arr, slope_arr, equal_nan=True)


def test_hillshade_tile_by_tile() -> None:
    out_hillshade_path = Path(r"test_data\TM1_564_146_test_tile_hillshade.tif")
    # if resolution, resolution_x, resolution_y, no_data are None in function_parameters dict,
    # they are taken from dem (dem_path)
    function_parameters = {"resolution_x": None, "resolution_y": None, "no_data": None}
    rvt.tile.save_visualization_tile_by_tile(
        visualization_function=rvt.vis.hillshade,
        function_parameters=function_parameters,
        dem_path=dem_path,
        overlap=1,
        tile_size_x=tile_size_x,
        tile_size_y=tile_size_y,
        out_raster_path=out_hillshade_path,
        out_raster_e_type=6,
        out_raster_nr_of_bands=1
    )
    hillshade_tile_by_tile_arr = rvt.default.get_raster_arr(out_hillshade_path.as_posix())["array"]
    dem_arr_dict = rvt.default.get_raster_arr(dem_path.as_posix())
    hillshade_arr = rvt.vis.hillshade(
        dem=dem_arr_dict["array"],
        resolution_x=dem_arr_dict["resolution"][0],
        resolution_y=dem_arr_dict["resolution"][1],
        no_data=dem_arr_dict["no_data"]
    )
    assert np.array_equal(hillshade_tile_by_tile_arr, hillshade_arr, equal_nan=True)


def test_svf_tile_by_tile() -> None:
    out_svf_path = Path(r"test_data\TM1_564_146_test_tile_svf.tif")
    # if resolution, resolution_x, resolution_y, no_data are None in function_parameters dict,
    # they are taken from dem (dem_path)
    svf_r_max = 10
    function_parameters = {"resolution": None, "no_data": None, "svf_r_max": svf_r_max}
    rvt.tile.save_visualization_tile_by_tile(
        visualization_function=rvt.vis.sky_view_factor,
        function_parameters=function_parameters,
        dem_path=dem_path,
        overlap=svf_r_max,
        tile_size_x=tile_size_x,
        tile_size_y=tile_size_y,
        out_raster_path=out_svf_path,
        out_raster_e_type=6,
        out_raster_nr_of_bands=1,
        out_visualization_dict_key="svf"
    )
    svf_tile_by_tile_arr = rvt.default.get_raster_arr(out_svf_path.as_posix())["array"]
    dem_arr_dict = rvt.default.get_raster_arr(dem_path.as_posix())
    svf_arr = rvt.vis.sky_view_factor(
        dem=dem_arr_dict["array"],
        resolution=dem_arr_dict["resolution"][0],
        no_data=dem_arr_dict["no_data"],
        svf_r_max=10
    )["svf"]
    assert np.array_equal(svf_tile_by_tile_arr, svf_arr, equal_nan=True)


def test_multi_hillshade_tile_by_tile() -> None:
    out_multi_hillshade_path = Path(r"test_data\TM1_564_146_test_tile_multi_hillshade.tif")
    # if resolution, resolution_x, resolution_y, no_data are None in function_parameters dict,
    # they are taken from dem (dem_path)
    nr_directions = 8
    function_parameters = {"resolution_x": None, "resolution_y": None, "no_data": None, "nr_directions": nr_directions}
    rvt.tile.save_visualization_tile_by_tile(
        visualization_function=rvt.vis.multi_hillshade,
        function_parameters=function_parameters,
        dem_path=dem_path,
        overlap=1,
        tile_size_x=tile_size_x,
        tile_size_y=tile_size_y,
        out_raster_path=out_multi_hillshade_path,
        out_raster_e_type=6,
        out_raster_nr_of_bands=nr_directions
    )
    multi_hillshade_tile_by_tile_arr = rvt.default.get_raster_arr(out_multi_hillshade_path.as_posix())["array"]
    dem_arr_dict = rvt.default.get_raster_arr(dem_path.as_posix())
    multi_hillshade_arr = rvt.vis.multi_hillshade(
        dem=dem_arr_dict["array"],
        resolution_x=dem_arr_dict["resolution"][0],
        resolution_y=dem_arr_dict["resolution"][1],
        no_data=dem_arr_dict["no_data"],
        nr_directions=nr_directions
    )
    assert np.array_equal(multi_hillshade_tile_by_tile_arr, multi_hillshade_arr, equal_nan=True)
