"""Tests for the focused Phase 25 pendulum approximation study."""

from __future__ import annotations

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.pendulum_approximation import (
    AMPLITUDES_RAD,
    PendulumApproximationResearch,
    nonlinear_period_reference,
    run_pendulum_approximation_research,
)


@pytest.fixture(scope="module")
def research_results() -> PendulumApproximationResearch:
    return run_pendulum_approximation_research()


def test_elliptic_period_reference_has_small_angle_limit_and_amplitude_growth() -> None:
    tiny = nonlinear_period_reference(1e-4, length_m=1.0, gravity_m_per_s2=9.81)
    reference = 2.0 * np.pi * np.sqrt(1.0 / 9.81)
    assert tiny == pytest.approx(reference, rel=1e-9)
    periods = [
        nonlinear_period_reference(angle, length_m=1.0, gravity_m_per_s2=9.81)
        for angle in AMPLITUDES_RAD
    ]
    assert all(right > left for left, right in zip(periods, periods[1:], strict=False))


@pytest.mark.parametrize(
    ("angle", "length", "gravity"),
    [(0.0, 1.0, 9.81), (np.pi, 1.0, 9.81), (0.5, 0.0, 9.81), (0.5, 1.0, -1.0)],
)
def test_period_reference_rejects_invalid_physical_domain(
    angle: float, length: float, gravity: float
) -> None:
    with pytest.raises(ScientificValidationError):
        nonlinear_period_reference(angle, length_m=length, gravity_m_per_s2=gravity)


def test_amplitude_sweep_separates_analytical_and_sampled_periods(
    research_results: PendulumApproximationResearch,
) -> None:
    results = research_results
    assert results.damping_coefficient_kg_m2_per_s == 0.0
    assert results.initial_angular_velocity_rad_per_s == 0.0
    assert len(results.amplitude_comparisons) == len(AMPLITUDES_RAD)
    for item in results.amplitude_comparisons:
        assert item.nonlinear_numerical_period_s == pytest.approx(
            item.nonlinear_analytical_period_s, abs=2e-7
        )
        assert item.nonlinear_analytical_period_s >= item.small_angle_period_s
        assert item.function_evaluations > 0
        assert tuple(window.window_small_angle_periods for window in item.windows) == (
            1,
            3,
            5,
            10,
        )


def test_phase_and_displacement_errors_grow_with_window_and_amplitude(
    research_results: PendulumApproximationResearch,
) -> None:
    results = research_results
    low = results.amplitude_comparisons[2]
    high = results.amplitude_comparisons[-1]
    assert high.small_angle_period_difference_s > low.small_angle_period_difference_s
    assert high.windows[-1].accumulated_phase_difference_rad > (
        low.windows[-1].accumulated_phase_difference_rad
    )
    assert low.windows[-1].accumulated_phase_difference_rad > (
        low.windows[0].accumulated_phase_difference_rad
    )
    assert low.windows[-1].maximum_absolute_angle_error_rad > (
        low.windows[0].maximum_absolute_angle_error_rad
    )


def test_tolerance_results_are_grid_limited_and_solver_checks_are_quantified(
    research_results: PendulumApproximationResearch,
) -> None:
    results = research_results
    assert len(results.numerical_sensitivities) == 3
    assert all(
        item.maximum_shared_sample_angle_difference_rad >= 0.0
        and item.numerical_period_difference_s is not None
        for item in results.numerical_sensitivities
    )
    assert (
        results.numerical_sensitivities[0].maximum_shared_sample_angle_difference_rad
        < 1e-7
    )
    assert (
        results.numerical_sensitivities[1].maximum_shared_sample_angle_difference_rad
        == 0.0
    )
