"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from dataclasses import fields, dataclass, is_dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

import pytest

from rvt.factory import (
    RVTVisualizationFactory,
    Visualization,
    HorizonVisualizations,
    Hillshade,
    MultipleDirectionsHillshade,
)


def _assert_rvt_vis_factories(
    expected_rvt_vis_factory: RVTVisualizationFactory,
    actual_rvt_vis_factory: RVTVisualizationFactory,
) -> bool:
    """Function to compare 2 dataclasses."""
    assert type(actual_rvt_vis_factory) == type(expected_rvt_vis_factory)

    for field in fields(expected_rvt_vis_factory):
        field_name = field.name
        expected_field_value = getattr(expected_rvt_vis_factory, field_name)
        actual_field_value = getattr(actual_rvt_vis_factory, field_name)

        if issubclass(type(expected_field_value), Visualization):
            assert actual_field_value.__dict__ == expected_field_value.__dict__

        else:
            assert actual_field_value == expected_field_value


@pytest.mark.parametrize(
    "rvt_vis_factory",
    (
        RVTVisualizationFactory(),
        RVTVisualizationFactory(
            multiple_directions_hillshade=MultipleDirectionsHillshade(sun_elevation=10),
            hillshade=Hillshade(sun_azimuth=60.123),
            horizon_visualizations=HorizonVisualizations(direction_of_anisotropy=300.5),
        ),
    ),
)
def test_rvt_visualization_factory_parameters_json(
    rvt_vis_factory: RVTVisualizationFactory,
) -> None:
    with NamedTemporaryFile(suffix=".json", delete=False) as temp_file:
        temp_file_path = Path(temp_file.name)
        rvt_vis_factory.save_parameters_to_file(file_path=temp_file_path)
        from_file_rvt_vis_factory = RVTVisualizationFactory.read_parameters_from_file(
            file_path=temp_file_path
        )
        _assert_rvt_vis_factories(
            expected_rvt_vis_factory=rvt_vis_factory,
            actual_rvt_vis_factory=from_file_rvt_vis_factory,
        )
