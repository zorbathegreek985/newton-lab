"""Tests for numerical input validation."""

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.utils.validation import validate_finite_array


def test_valid_array_is_float64() -> None:
    result = validate_finite_array([1, 2.5, 3], ndim=1)

    assert result.dtype == np.float64
    np.testing.assert_array_equal(result, np.array([1.0, 2.5, 3.0]))


@pytest.mark.parametrize("values", [[1, float("nan")], [float("inf")], [-float("inf")]])
def test_rejects_non_finite_values(values: list[float]) -> None:
    with pytest.raises(ScientificValidationError, match="finite"):
        validate_finite_array(values, name="measurements")


def test_rejects_non_numeric_values() -> None:
    with pytest.raises(ScientificValidationError, match="numeric"):
        validate_finite_array(["not a number"])


def test_rejects_wrong_dimension() -> None:
    with pytest.raises(ScientificValidationError, match="2 dimensions"):
        validate_finite_array([[1, 2]], ndim=1)
