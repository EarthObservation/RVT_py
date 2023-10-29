"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from pathlib import Path

import pytest

from rvt.enums import RVTVisualizationName
from rvt.factory import RVTVisualizationFactory, MultipleDirectionsHillshade


@pytest.mark.parametrize(
    "rvt_vis_factory, directory_path, dem_file_name, is_8bit, expected_path",
    [
        (
            RVTVisualizationFactory(
                multiple_directions_hillshade=MultipleDirectionsHillshade(
                    number_of_directions=12, sun_elevation=40.5
                )
            ),
            Path("/some/directory"),
            "some_dem",
            False,
            Path("/some/directory/some_dem_MULTI-HS_D12_H40.5.tif"),
        ),
        (
            RVTVisualizationFactory(
                multiple_directions_hillshade=MultipleDirectionsHillshade(
                    number_of_directions=16, sun_elevation=50
                )
            ),
            Path("/results/"),
            "DEM_123",
            True,
            Path("/results/DEM_123_MULTI-HS_D16_H50_8bit.tif"),
        ),
    ],
)
def test_get_multiple_directions_hillshade_path(
    rvt_vis_factory: RVTVisualizationFactory,
    directory_path: Path,
    dem_file_name: str,
    is_8bit: bool,
    expected_path: Path,
) -> None:
    path_from_factory_specific_method = (
        rvt_vis_factory.get_multiple_directions_hillshade_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    path_from_factory_general_method = rvt_vis_factory.get_visualization_path(
        visualization_name=RVTVisualizationName.MULTIPLE_DIRECTIONS_HILLSHADE,
        directory_path=directory_path,
        dem_file_name=dem_file_name,
        is_8bit=is_8bit,
    )
    path_from_rvt_visualization_method = (
        rvt_vis_factory.multiple_directions_hillshade.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    assert path_from_factory_specific_method == expected_path
    assert path_from_factory_general_method == expected_path
    assert path_from_rvt_visualization_method == expected_path
