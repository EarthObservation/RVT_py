"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from pathlib import Path

import pytest

from rvt.enums import RVTVisualizationName
from rvt.factory import RVTVisualizationFactory, MultiScaleTopographicPosition


@pytest.mark.parametrize(
    "rvt_vis_factory, directory_path, dem_file_name, is_8bit, expected_path",
    [
        (
            RVTVisualizationFactory(
                multi_scale_topographic_position=MultiScaleTopographicPosition(
                    local_scale=(3, 21, 2),
                    meso_scale=(23, 203, 18),
                    broad_scale=(223, 2023, 180),
                    lightness=1.2,
                )
            ),
            Path("/some/directory"),
            "some_dem",
            False,
            Path("/some/directory/some_dem_MSTP_21_203_2023_L1.2.tif"),
        ),
        (
            RVTVisualizationFactory(
                multi_scale_topographic_position=MultiScaleTopographicPosition(
                    local_scale=(10, 20, 1),
                    meso_scale=(50, 250, 50),
                    broad_scale=(250, 750, 100),
                    lightness=1.4,
                )
            ),
            Path("/results/"),
            "DEM_123",
            True,
            Path("/results/DEM_123_MSTP_20_250_750_L1.4_8bit.tif"),
        ),
    ],
)
def test_get_multi_scale_topographic_position_path(
    rvt_vis_factory: RVTVisualizationFactory,
    directory_path: Path,
    dem_file_name: str,
    is_8bit: bool,
    expected_path: Path,
) -> None:
    path_from_factory_specific_method = (
        rvt_vis_factory.get_multi_scale_topographic_position_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    path_from_factory_general_method = rvt_vis_factory.get_visualization_path(
        visualization_name=RVTVisualizationName.MULTI_SCALE_TOPOGRAPHIC_POSITION,
        directory_path=directory_path,
        dem_file_name=dem_file_name,
        is_8bit=is_8bit,
    )
    path_from_rvt_visualization_method = (
        rvt_vis_factory.multi_scale_topographic_position.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    assert path_from_factory_specific_method == expected_path
    assert path_from_factory_general_method == expected_path
    assert path_from_rvt_visualization_method == expected_path
