"""Tests for causality and the Phase 30 exponential-filter study."""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import cast

import numpy as np
import pytest

from newton_lab.causal_diffusion_filter import (
    exponential_filter,
    frequency_response,
    run_causal_filter_study,
)
from newton_lab.exceptions import ScientificValidationError


def test_filter_does_not_anticipate_changes_after_a_prefix() -> None:
    first = np.array([0.0, 1.0, -0.5, 0.25, 900.0, 2.0])
    second = first.copy()
    second[4:] = [-700.0, 1200.0]
    through = 3
    left = exponential_filter(first, 0.83)
    right = exponential_filter(second, 0.83)
    np.testing.assert_array_equal(left[: through + 1], right[: through + 1])


def test_beta_zero_is_identity_and_initial_sample_is_preserved() -> None:
    values = [3.0, -1.0, 4.0, 0.0]
    np.testing.assert_array_equal(exponential_filter(values, 0.0), values)
    assert exponential_filter(values, 0.7)[0] == values[0]
    assert values == [3.0, -1.0, 4.0, 0.0]


def test_constant_input_remains_constant_and_outputs_are_finite() -> None:
    output = exponential_filter([2.5] * 100, 0.97)
    np.testing.assert_array_equal(output, np.full(100, 2.5))
    assert np.isfinite(output).all()


@pytest.mark.parametrize("beta", [-0.01, 1.0, math.inf, math.nan, True, "bad"])
def test_invalid_beta_is_rejected(beta: object) -> None:
    with pytest.raises(ScientificValidationError):
        exponential_filter([1.0, 2.0], cast(float, beta))


@pytest.mark.parametrize("values", [[], [[1.0, 2.0]], [1.0, math.nan], [math.inf]])
def test_invalid_input_is_rejected(values: object) -> None:
    with pytest.raises(ScientificValidationError):
        exponential_filter(cast(Sequence[float], values), 0.5)


def test_analytic_response_and_group_delay_match_finite_difference() -> None:
    beta = 0.8
    omega = 0.4
    response = frequency_response(beta, omega)
    expected_magnitude = (1.0 - beta) / math.sqrt(
        1.0 + beta**2 - 2.0 * beta * math.cos(omega)
    )
    assert response.magnitude == pytest.approx(expected_magnitude)
    assert response.phase_rad < 0.0
    step = 1e-6
    phase_plus = frequency_response(beta, omega + step).phase_rad
    phase_minus = frequency_response(beta, omega - step).phase_rad
    numeric_group_delay = -(phase_plus - phase_minus) / (2.0 * step)
    assert response.group_delay_samples == pytest.approx(numeric_group_delay, abs=1e-9)
    assert frequency_response(beta, 0.0).group_delay_samples == pytest.approx(
        beta / (1.0 - beta)
    )


def test_frequency_checks_match_analytic_magnitude_and_phase() -> None:
    study = run_causal_filter_study(beta_sweep=(0.0, 0.7), sample_count=201)
    assert len(study.frequency_checks) == 3
    assert max(row.magnitude_absolute_error for row in study.frequency_checks) < 0.003
    assert max(row.phase_absolute_error_rad for row in study.frequency_checks) < 0.003


def test_sweep_preserves_order_repeats_invalids_and_records_causality() -> None:
    study = run_causal_filter_study(
        beta_sweep=(0.0, 0.8, 0.8, -0.1, 1.0), sample_count=201
    )
    assert study.beta_sweep == (0.0, 0.8, 0.8, -0.1, 1.0)
    assert [row.sweep_index for row in study.sweep_rows] == list(range(5))
    assert [row.status for row in study.sweep_rows] == [
        "succeeded",
        "succeeded",
        "succeeded",
        "failed",
        "failed",
    ]
    assert study.sweep_rows[0].filtered_rmse_to_truth is not None
    assert study.sweep_rows[0].filtered_rmse_to_truth == pytest.approx(
        study.sweep_rows[0].input_rmse_to_truth
    )
    assert study.sweep_rows[0].step_pretransition_max_abs == 0.0
    assert study.diffusion_mapping[0].matched_beta != pytest.approx(
        study.diffusion_mapping[-1].matched_beta
    )


def test_higher_beta_increases_step_width_and_low_frequency_delay() -> None:
    betas = (0.0, 0.5, 0.9, 0.98)
    study = run_causal_filter_study(beta_sweep=betas, sample_count=301)
    widths = [row.step_10_90_width_samples for row in study.sweep_rows]
    delays = [frequency_response(beta, 0.0).group_delay_samples for beta in betas]
    assert all(width is not None for width in widths)
    measured_widths = [width for width in widths if width is not None]
    assert len(measured_widths) == len(widths)
    assert measured_widths == sorted(measured_widths)
    assert delays == sorted(delays)
