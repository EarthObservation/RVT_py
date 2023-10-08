"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from pathlib import Path
from tempfile import NamedTemporaryFile

import pytest

from rvt.default import DefaultValues


# TODO: Rename from default


@pytest.mark.parametrize("default", ((DefaultValues(),)))
def test_parameters_json(default: DefaultValues) -> None:
    with NamedTemporaryFile(suffix=".json", delete=False) as temp_file:
        temp_file_path = Path(temp_file.name)
        default.save_to_file(file_path=temp_file_path)
        from_file_default = DefaultValues.read_from_file(file_path=temp_file_path)
        assert default == from_file_default
