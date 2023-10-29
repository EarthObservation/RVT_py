"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from pathlib import Path

import pytest

from rvt.enums import RVTVisualizationName
from rvt.factory import RVTVisualizationFactory, LocalDominance


@pytest.mark.parametrize(
    "rvt_vis_factory, directory_path, dem_file_name, is_8bit, expected_path",
    [
        (
            RVTVisualizationFactory(
                local_dominance=LocalDominance(
                    minimum_radius=10,
                    maximum_radius=20,
                    radial_distance_step=1,
                    angular_step=15,
                    observer_height=1.7,
                ),
            ),
            Path("/some/directory"),
            "some_dem",
            False,
            Path("/some/directory/some_dem_LD_R_M10-20_DI1_A15_OH1.7.tif"),
        ),
        (
            RVTVisualizationFactory(
                local_dominance=LocalDominance(
                    minimum_radius=10,
                    maximum_radius=40,
                    radial_distance_step=2,
                    angular_step=25,
                    observer_height=1.9,
                ),
            ),
            Path("/results/"),
            "DEM_123",
            True,
            Path("/results/DEM_123_LD_R_M10-40_DI2_A25_OH1.9_8bit.tif"),
        ),
    ],
)
def test_get_local_dominance_path(
    rvt_vis_factory: RVTVisualizationFactory,
    directory_path: Path,
    dem_file_name: str,
    is_8bit: bool,
    expected_path: Path,
) -> None:
    path_from_factory_specific_method = rvt_vis_factory.get_local_dominance_path(
        directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
    )
    path_from_factory_general_method = rvt_vis_factory.get_visualization_path(
        visualization_name=RVTVisualizationName.LOCAL_DOMINANCE,
        directory_path=directory_path,
        dem_file_name=dem_file_name,
        is_8bit=is_8bit,
    )
    path_from_rvt_visualization_method = (
        rvt_vis_factory.local_dominance.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    assert path_from_factory_specific_method == expected_path
    assert path_from_factory_general_method == expected_path
    assert path_from_rvt_visualization_method == expected_path
