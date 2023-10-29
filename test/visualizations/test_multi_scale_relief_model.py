"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from pathlib import Path

import pytest

from rvt.enums import RVTVisualizationName
from rvt.factory import RVTVisualizationFactory, MultiScaleReliefModel


@pytest.mark.parametrize(
    "rvt_vis_factory, directory_path, dem_file_name, is_8bit, expected_path",
    [
        (
            RVTVisualizationFactory(
                multi_scale_relief_model=MultiScaleReliefModel(
                    minimum_feature_size=0.0,
                    maximum_feature_size=20.0,
                    scaling_factor=2,
                )
            ),
            Path("/some/directory"),
            "some_dem",
            False,
            Path("/some/directory/some_dem_MSRM_FS0.0-20.0_S2.tif"),
        ),
        (
            RVTVisualizationFactory(
                multi_scale_relief_model=MultiScaleReliefModel(
                    minimum_feature_size=10,
                    maximum_feature_size=30,
                    scaling_factor=3,
                )
            ),
            Path("/results/"),
            "DEM_123",
            True,
            Path("/results/DEM_123_MSRM_FS10-30_S3_8bit.tif"),
        ),
    ],
)
def test_get_multi_scale_relief_model_path(
    rvt_vis_factory: RVTVisualizationFactory,
    directory_path: Path,
    dem_file_name: str,
    is_8bit: bool,
    expected_path: Path,
) -> None:
    path_from_factory_specific_method = (
        rvt_vis_factory.get_multi_scale_relief_model_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    path_from_factory_general_method = rvt_vis_factory.get_visualization_path(
        visualization_name=RVTVisualizationName.MULTI_SCALE_RELIEF_MODEL,
        directory_path=directory_path,
        dem_file_name=dem_file_name,
        is_8bit=is_8bit,
    )
    path_from_rvt_visualization_method = (
        rvt_vis_factory.multi_scale_relief_model.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    assert path_from_factory_specific_method == expected_path
    assert path_from_factory_general_method == expected_path
    assert path_from_rvt_visualization_method == expected_path
