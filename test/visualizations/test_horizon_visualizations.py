"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from pathlib import Path
from typing import Dict

import pytest

from rvt.enums import SvfNoiseRemove, AnisotropyLevel, RVTVisualizationName
from rvt.factory import RVTVisualizationFactory, HorizonVisualizations


@pytest.mark.parametrize(
    "rvt_vis_factory, directory_path, dem_file_name, is_8bit, expected_paths",
    [
        (
            RVTVisualizationFactory(
                horizon_visualizations=HorizonVisualizations(
                    number_of_directions=16,
                    maximum_search_radius=10,
                    noise_remove=SvfNoiseRemove.NO_REMOVE,
                    direction_of_anisotropy=300.5,
                    anisotropy_level=AnisotropyLevel.LOW,
                )
            ),
            Path("/some/directory"),
            "some_dem",
            False,
            {
                "sky_view_factor": Path(
                    "/some/directory/some_dem_SVF_R10_D16_NRno_remove.tif"
                ),
                "anisotropic_sky_view_factor": Path(
                    "/some/directory/some_dem_SVF-A_R10_D16_NRno_remove_A300.5_ALlow.tif"
                ),
                "positive_openness": Path(
                    "/some/directory/some_dem_OPEN-POS_R10_D16_NRno_remove.tif"
                ),
                "negative_openness": Path(
                    "/some/directory/some_dem_OPEN-NEG_R10_D16_NRno_remove.tif"
                ),
            },
        ),
        (
            RVTVisualizationFactory(
                horizon_visualizations=HorizonVisualizations(
                    number_of_directions=12,
                    maximum_search_radius=5,
                    noise_remove=SvfNoiseRemove.HIGH,
                    direction_of_anisotropy=255,
                    anisotropy_level=AnisotropyLevel.HIGH,
                )
            ),
            Path("/results/"),
            "DEM_123",
            True,
            {
                "sky_view_factor": Path("/results/DEM_123_SVF_R5_D12_NRhigh_8bit.tif"),
                "anisotropic_sky_view_factor": Path(
                    "/results/DEM_123_SVF-A_R5_D12_NRhigh_A255_ALhigh_8bit.tif"
                ),
                "positive_openness": Path(
                    "/results/DEM_123_OPEN-POS_R5_D12_NRhigh_8bit.tif"
                ),
                "negative_openness": Path(
                    "/results/DEM_123_OPEN-NEG_R5_D12_NRhigh_8bit.tif"
                ),
            },
        ),
    ],
)
def test_get_horizon_visualizations_paths(
    rvt_vis_factory: RVTVisualizationFactory,
    directory_path: Path,
    dem_file_name: str,
    is_8bit: bool,
    expected_paths: Dict[str, Path],
) -> None:
    # Sky View Factor
    svf_path_from_factory_specific_method = rvt_vis_factory.get_sky_view_factor_path(
        directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
    )
    svf_path_from_factory_general_method = rvt_vis_factory.get_visualization_path(
        visualization_name=RVTVisualizationName.SKY_VIEW_FACTOR,
        directory_path=directory_path,
        dem_file_name=dem_file_name,
        is_8bit=is_8bit,
    )
    svf_path_from_rvt_visualization_method = (
        rvt_vis_factory.horizon_visualizations.sky_view_factor.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    assert svf_path_from_factory_specific_method == expected_paths["sky_view_factor"]
    assert svf_path_from_factory_general_method == expected_paths["sky_view_factor"]
    assert svf_path_from_rvt_visualization_method == expected_paths["sky_view_factor"]

    # Anisotropic Sky View Factor
    asvf_path_from_factory_specific_method = (
        rvt_vis_factory.get_anisotropic_sky_view_factor_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    asvf_path_from_factory_general_method = rvt_vis_factory.get_visualization_path(
        visualization_name=RVTVisualizationName.ANISOTROPIC_SKY_VIEW_FACTOR,
        directory_path=directory_path,
        dem_file_name=dem_file_name,
        is_8bit=is_8bit,
    )
    asvf_path_from_rvt_visualization_method = rvt_vis_factory.horizon_visualizations.anisotropic_sky_view_factor.get_visualization_path(
        directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
    )
    assert (
        asvf_path_from_factory_specific_method
        == expected_paths["anisotropic_sky_view_factor"]
    )
    assert (
        asvf_path_from_factory_general_method
        == expected_paths["anisotropic_sky_view_factor"]
    )
    assert (
        asvf_path_from_rvt_visualization_method
        == expected_paths["anisotropic_sky_view_factor"]
    )

    # Positive Openness
    pos_opns_path_from_factory_specific_method = (
        rvt_vis_factory.get_positive_openness_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    pos_opns_path_from_factory_general_method = rvt_vis_factory.get_visualization_path(
        visualization_name=RVTVisualizationName.POSITIVE_OPENNESS,
        directory_path=directory_path,
        dem_file_name=dem_file_name,
        is_8bit=is_8bit,
    )
    pos_opns_path_from_rvt_visualization_method = (
        rvt_vis_factory.horizon_visualizations.positive_openness.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    assert (
        pos_opns_path_from_factory_specific_method
        == expected_paths["positive_openness"]
    )
    assert (
        pos_opns_path_from_factory_general_method == expected_paths["positive_openness"]
    )
    assert (
        pos_opns_path_from_rvt_visualization_method
        == expected_paths["positive_openness"]
    )

    # Negative Openness
    neg_opns_path_from_factory_specific_method = (
        rvt_vis_factory.get_negative_openness_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    neg_opns_path_from_factory_general_method = rvt_vis_factory.get_visualization_path(
        visualization_name=RVTVisualizationName.NEGATIVE_OPENNESS,
        directory_path=directory_path,
        dem_file_name=dem_file_name,
        is_8bit=is_8bit,
    )
    neg_opns_path_from_rvt_visualization_method = (
        rvt_vis_factory.horizon_visualizations.negative_openness.get_visualization_path(
            directory_path=directory_path, dem_file_name=dem_file_name, is_8bit=is_8bit
        )
    )
    assert (
        neg_opns_path_from_factory_specific_method
        == expected_paths["negative_openness"]
    )
    assert (
        neg_opns_path_from_factory_general_method == expected_paths["negative_openness"]
    )
    assert (
        neg_opns_path_from_rvt_visualization_method
        == expected_paths["negative_openness"]
    )
