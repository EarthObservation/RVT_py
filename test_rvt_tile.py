from pathlib import Path

import rvt.enums
import rvt.tile
import rvt.visualizations
import rvt.default
import numpy as np

# pytest rvt.tile.save_rvt_visualization_tile_by_tile

dem_path = Path(r"test_data\TM1_564_146.tif")
default_values = rvt.default.DefaultValues()
default_values.tile_size = (100, 100)
default_values.tile_size_limit = 100 * 100


def test_rvt_slope_tile_by_tile() -> None:
    out_slope_float_path = Path(default_values.get_slope_path(dem_path))
    out_slope_8bit_path = Path(default_values.get_slope_path(dem_path, True))
    rvt.tile.save_rvt_visualization_tile_by_tile(
        rvt_visualization=rvt.enums.RVTVisualization.SLOPE,
        rvt_default=default_values,
        dem_path=dem_path,
        output_dir_path=dem_path.parent,
        save_float=True,
        save_8bit=True,
    )
    slope_float_tile_by_tile_arr = rvt.default.get_raster_arr(
        out_slope_float_path.as_posix()
    )["array"]
    slope_8bit_tile_by_tile_arr = rvt.default.get_raster_arr(
        out_slope_8bit_path.as_posix()
    )["array"]
    dem_arr_dict = rvt.default.get_raster_arr(dem_path.as_posix())
    slope_arr = default_values.get_slope(
        dem_arr=dem_arr_dict["array"],
        resolution_x=dem_arr_dict["resolution"][0],
        resolution_y=dem_arr_dict["resolution"][1],
        no_data=dem_arr_dict["no_data"],
    )
    slope_8bit_arr = default_values.float_to_8bit(
        float_arr=slope_arr, visualization=rvt.enums.RVTVisualization.SLOPE
    )
    assert np.array_equal(slope_float_tile_by_tile_arr, slope_arr, equal_nan=True)
    assert np.array_equal(slope_8bit_tile_by_tile_arr, slope_8bit_arr, equal_nan=True)


def test_rvt_hillshade_tile_by_tile() -> None:
    out_hillshade_float_path = Path(default_values.get_hillshade_path(dem_path))
    out_hillshade_8bit_path = Path(default_values.get_hillshade_path(dem_path, True))
    rvt.tile.save_rvt_visualization_tile_by_tile(
        rvt_visualization=rvt.enums.RVTVisualization.HILLSHADE,
        rvt_default=default_values,
        dem_path=dem_path,
        output_dir_path=dem_path.parent,
        save_float=True,
        save_8bit=True,
    )
    hillshade_float_tile_by_tile_arr = rvt.default.get_raster_arr(
        out_hillshade_float_path.as_posix()
    )["array"]
    hillshade_8bit_tile_by_tile_arr = rvt.default.get_raster_arr(
        out_hillshade_8bit_path.as_posix()
    )["array"]
    dem_arr_dict = rvt.default.get_raster_arr(dem_path.as_posix())
    hillshade_arr = default_values.get_hillshade(
        dem_arr=dem_arr_dict["array"],
        resolution_x=dem_arr_dict["resolution"][0],
        resolution_y=dem_arr_dict["resolution"][1],
        no_data=dem_arr_dict["no_data"],
    )
    hillshade_8bit_arr = default_values.float_to_8bit(
        float_arr=hillshade_arr, visualization=rvt.enums.RVTVisualization.HILLSHADE
    )
    assert np.array_equal(
        hillshade_float_tile_by_tile_arr, hillshade_arr, equal_nan=True
    )
    assert np.array_equal(
        hillshade_8bit_tile_by_tile_arr, hillshade_8bit_arr, equal_nan=True
    )


def test_rvt_svf_tile_by_tile() -> None:
    out_svf_float_path = Path(default_values.get_svf_path(dem_path))
    out_svf_8bit_path = Path(default_values.get_svf_path(dem_path, True))
    rvt.tile.save_rvt_visualization_tile_by_tile(
        rvt_visualization=rvt.enums.RVTVisualization.SKY_VIEW_FACTOR,
        rvt_default=default_values,
        dem_path=dem_path,
        output_dir_path=dem_path.parent,
        save_float=True,
        save_8bit=True,
    )
    svf_float_tile_by_tile_arr = rvt.default.get_raster_arr(
        out_svf_float_path.as_posix()
    )["array"]
    svf_8bit_tile_by_tile_arr = rvt.default.get_raster_arr(
        out_svf_8bit_path.as_posix()
    )["array"]
    dem_arr_dict = rvt.default.get_raster_arr(dem_path.as_posix())
    svf_arr = default_values.get_sky_view_factor(
        dem_arr=dem_arr_dict["array"],
        resolution=dem_arr_dict["resolution"][0],
        no_data=dem_arr_dict["no_data"],
    )["svf"]
    svf_8bit_arr = default_values.float_to_8bit(
        float_arr=svf_arr, visualization=rvt.enums.RVTVisualization.SKY_VIEW_FACTOR
    )
    assert np.array_equal(svf_float_tile_by_tile_arr, svf_arr, equal_nan=True)
    assert np.array_equal(svf_8bit_tile_by_tile_arr, svf_8bit_arr, equal_nan=True)


def test_rvt_multi_hillshade_tile_by_tile() -> None:
    out_mhs_float_path = Path(default_values.get_multi_hillshade_path(dem_path))
    out_mhs_8bit_path = Path(default_values.get_multi_hillshade_path(dem_path, True))
    rvt.tile.save_rvt_visualization_tile_by_tile(
        rvt_visualization=rvt.enums.RVTVisualization.MULTIPLE_DIRECTIONS_HILLSHADE,
        rvt_default=default_values,
        dem_path=dem_path,
        output_dir_path=dem_path.parent,
        save_float=True,
        save_8bit=True,
    )
    mhs_float_tile_by_tile_arr = rvt.default.get_raster_arr(
        out_mhs_float_path.as_posix()
    )["array"]
    mhs_8bit_tile_by_tile_arr = rvt.default.get_raster_arr(
        out_mhs_8bit_path.as_posix()
    )["array"]
    dem_arr_dict = rvt.default.get_raster_arr(dem_path.as_posix())
    mhs_arr = default_values.get_multi_hillshade(
        dem_arr=dem_arr_dict["array"],
        resolution_x=dem_arr_dict["resolution"][0],
        resolution_y=dem_arr_dict["resolution"][1],
        no_data=dem_arr_dict["no_data"],
    )
    mhs_8bit_arr = default_values.float_to_8bit(
        float_arr=dem_arr_dict["array"],
        visualization=rvt.enums.RVTVisualization.MULTIPLE_DIRECTIONS_HILLSHADE,
        x_res=dem_arr_dict["resolution"][0],
        y_res=dem_arr_dict["resolution"][1],
        no_data=dem_arr_dict["no_data"],
    )
    assert np.array_equal(mhs_float_tile_by_tile_arr, mhs_arr, equal_nan=True)
    assert np.array_equal(mhs_8bit_tile_by_tile_arr, mhs_8bit_arr, equal_nan=True)
