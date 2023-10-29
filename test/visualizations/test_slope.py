"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from pathlib import Path

from rvt.enums import SlopeOutputUnit, RVTVisualizationName
from rvt.factory import RVTVisualizationFactory, Slope
import pytest


@pytest.mark.parametrize(
    "rvt_vis_factory, directory_path, dem_file_name, is_8bit, expected_path",
    [
        (
            RVTVisualizationFactory(
                slope=Slope(output_unit=SlopeOutputUnit.RADIAN),
            ),
            Path("/some/directory"),
            "some_dem",
            False,
            Path("/some/directory/some_dem_SLOPE.tif"),
        ),
        (
            RVTVisualizationFactory(
                slope=Slope(output_unit=SlopeOutputUnit.DEGREE),
            ),
            Path("/results/"),
            "DEM_123",
            True,
            Path("/results/DEM_123_SLOPE_8bit.tif"),
        ),
    ],
)
def test_get_slope_path(
    rvt_vis_factory: RVTVisualizationFactory,
    directory_path: Path,
    dem_file_name: str,
    is_8bit: bool,
    expected_path: Path,
) -> None:
    path_from_factory_specific_method = rvt_vis_factory.get_slope_path(
        directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
    )
    path_from_factory_general_method = rvt_vis_factory.get_visualization_path(
        visualization_name=RVTVisualizationName.SLOPE,
        directory_path=directory_path,
        dem_file_name=dem_file_name,
        is_8bit=is_8bit,
    )
    path_from_rvt_visualization_method = rvt_vis_factory.slope.get_visualization_path(
        directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
    )
    assert path_from_factory_specific_method == expected_path
    assert path_from_factory_general_method == expected_path
    assert path_from_rvt_visualization_method == expected_path
