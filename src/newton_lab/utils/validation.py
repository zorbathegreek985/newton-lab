"""Validation helpers for numerical inputs."""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from newton_lab.exceptions import ScientificValidationError


def validate_finite_scalar(value: float, *, name: str = "value") -> float:
    """Return a scalar as ``float`` after checking it is numeric and finite.

    Boolean values are rejected even though Python treats them as integers.

    Raises:
        ScientificValidationError: If the value is boolean, non-numeric, or
            non-finite.
    """
    if isinstance(value, bool):
        raise ScientificValidationError(f"{name} must be a finite number")
    try:
        converted = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ScientificValidationError(f"{name} must be a finite number") from exc
    if not np.isfinite(converted):
        raise ScientificValidationError(f"{name} must be finite")
    return converted


def validate_finite_array(
    values: ArrayLike,
    *,
    name: str = "values",
    ndim: int | None = None,
) -> NDArray[np.float64]:
    """Convert input to float64 and require all entries to be finite.

    Args:
        values: Numeric data accepted by ``numpy.asarray``.
        name: Name used in validation error messages.
        ndim: Optional required number of dimensions.

    Raises:
        ScientificValidationError: If conversion, dimensionality, or finiteness
            validation fails.
    """
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ScientificValidationError(f"{name} must contain numeric values") from exc

    if ndim is not None and array.ndim != ndim:
        expected_unit = "dimension" if ndim == 1 else "dimensions"
        actual_unit = "dimension" if array.ndim == 1 else "dimensions"
        raise ScientificValidationError(
            f"{name} must have {ndim} {expected_unit}; got {array.ndim} {actual_unit}"
        )
    if not np.isfinite(array).all():
        raise ScientificValidationError(f"{name} must contain only finite values")
    return array
