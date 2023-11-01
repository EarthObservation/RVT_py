"""
Copyright:
    2010-2023 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2023 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from typing import Any

import numpy as np
import pytest
import numpy.typing as npt
from numpy.testing import assert_array_equal

from rvt.blender_functions import normalize_lin


@pytest.mark.parametrize(
    "input_image_array, minimum, maximum, expected_image_array",
    [
        (
            np.array([[1, 2, 3], [2, 2, 1]]),
            1,
            3,
            np.array([[0, 0.5, 1], [0.5, 0.5, 0]]),
        ),
        (
            np.array(
                [[[1, 2, 3], [2, 2, 1]], [[1, 2, 4], [2, 2, 1]], [[1, 2, 3], [1, 1, 1]]]
            ),
            1,
            3,
            np.array(
                [
                    [[0, 0.5, 1], [0.5, 0.5, 0]],
                    [[0, 0.5, 1], [0.5, 0.5, 0]],
                    [[0, 0.5, 1], [0, 0, 0]],
                ]
            ),
        ),
    ],
)
def test_normalize_lin(
    input_image_array: npt.NDArray[Any],
    minimum: float,
    maximum: float,
    expected_image_array: npt.NDArray[Any],
) -> None:
    assert_array_equal(
        normalize_lin(image=input_image_array, minimum=minimum, maximum=maximum),
        expected_image_array,
    )
