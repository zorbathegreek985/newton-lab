"""Scientific tests for Phase 26's velocity-aware pendulum study."""

from __future__ import annotations

from math import pi, sqrt

import numpy as np
import pytest
from scipy.special import ellipk  # type: ignore[import-untyped]

import newton_lab.pendulum_initial_velocity as phase26
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.pendulum import simulate_damped_pendulum
from newton_lab.pendulum_approximation import (
    nonlinear_period_reference as phase25_period,
)
from newton_lab.pendulum_initial_velocity import (
    AccuracyCriteria,
    CaseStatus,
    InitialConditionSpec,
    InitialStateClass,
    classify_initial_condition,
    estimate_nonlinear_period,
    linear_period_reference,
    nonlinear_period_reference,
    render_initial_velocity_report,
    run_initial_velocity_research,
    turning_point_amplitude,
)


@pytest.fixture(scope="module")
def research_results() -> phase26.PendulumInitialVelocityResearch:
    return run_initial_velocity_research()


def _classify(angle: float, velocity: float) -> phase26.InitialStateClassification:
    return classify_initial_condition(
        angle,
        velocity,
        length_m=1.0,
        gravity_m_per_s2=9.81,
    )


def test_turning_point_amplitude_for_release_from_rest_matches_initial_angle() -> None:
    assert turning_point_amplitude(
        0.8, 0.0, length_m=1.0, gravity_m_per_s2=9.81
    ) == pytest.approx(0.8, abs=1e-14)


def test_nonzero_velocity_increases_turning_amplitude_by_energy() -> None:
    angle, speed = 0.5, 3.0
    actual = turning_point_amplitude(angle, speed, length_m=1.0, gravity_m_per_s2=9.81)
    expected = 2.0 * np.arcsin(
        np.sqrt((2.0 * np.sin(angle / 2.0) ** 2 + speed**2 / (2.0 * 9.81)) / 2.0)
    )
    assert actual == pytest.approx(expected)
    assert actual - angle > 0.5


def test_turning_amplitude_and_period_are_symmetric_under_velocity_reversal() -> None:
    positive = turning_point_amplitude(0.5, 3.0, length_m=1.0, gravity_m_per_s2=9.81)
    negative = turning_point_amplitude(0.5, -3.0, length_m=1.0, gravity_m_per_s2=9.81)
    assert positive == pytest.approx(negative, abs=1e-14)
    assert nonlinear_period_reference(
        0.5, 3.0, length_m=1.0, gravity_m_per_s2=9.81
    ) == pytest.approx(
        nonlinear_period_reference(0.5, -3.0, length_m=1.0, gravity_m_per_s2=9.81)
    )


@pytest.mark.parametrize(
    ("angle", "velocity", "length", "gravity"),
    [
        (0.2, 0.0, 0.0, 9.81),
        (0.2, 0.0, 1.0, -9.81),
        (pi + 0.01, 0.0, 1.0, 9.81),
        (float("nan"), 0.0, 1.0, 9.81),
        (0.2, float("inf"), 1.0, 9.81),
    ],
)
def test_classification_rejects_invalid_parameters_and_initial_conditions(
    angle: float, velocity: float, length: float, gravity: float
) -> None:
    with pytest.raises(ScientificValidationError):
        classify_initial_condition(
            angle, velocity, length_m=length, gravity_m_per_s2=gravity
        )


def test_near_boundary_libration_is_classified_and_not_clamped() -> None:
    state = _classify(3.12, 0.0)
    assert state.classification is InitialStateClass.NEAR_SEPARATRIX
    assert state.turning_point_amplitude_rad is not None
    assert state.turning_point_amplitude_rad == pytest.approx(3.12, abs=1e-12)
    assert state.turning_point_amplitude_rad < pi
    assert state.distance_below_separatrix > 0.0
    assert state.distance_below_separatrix < phase26.NEAR_SEPARATRIX_ENERGY_GAP


def test_roundoff_separatrix_and_rotation_never_receive_libration_period() -> None:
    separatrix_speed = sqrt(4.0 * 9.81)
    state = _classify(0.0, separatrix_speed)
    rotation = _classify(0.0, separatrix_speed * (1.0 + 1e-6))
    assert state.classification is InitialStateClass.SEPARATRIX
    assert rotation.classification is InitialStateClass.ROTATING
    assert _classify(pi, 1e-7).classification is InitialStateClass.SEPARATRIX
    assert _classify(pi, 1e-6).classification is InitialStateClass.ROTATING
    with pytest.raises(ScientificValidationError, match="period is unavailable"):
        nonlinear_period_reference(
            0.0, separatrix_speed, length_m=1.0, gravity_m_per_s2=9.81
        )
    with pytest.raises(ScientificValidationError, match="exceeds the separatrix"):
        nonlinear_period_reference(0.0, 7.0, length_m=1.0, gravity_m_per_s2=9.81)
    with pytest.raises(ScientificValidationError, match="not a libration amplitude"):
        turning_point_amplitude(0.0, 7.0, length_m=1.0, gravity_m_per_s2=9.81)


def test_elliptic_reference_uses_scipy_parameter_not_modulus() -> None:
    angle, speed = 0.5, 3.0
    amplitude = turning_point_amplitude(
        angle, speed, length_m=1.0, gravity_m_per_s2=9.81
    )
    parameter = np.sin(amplitude / 2.0) ** 2
    expected = 4.0 / sqrt(9.81) * ellipk(parameter)
    actual = nonlinear_period_reference(
        angle, speed, length_m=1.0, gravity_m_per_s2=9.81
    )
    assert actual == pytest.approx(expected, rel=1e-14)
    wrong_modulus_convention = 4.0 / sqrt(9.81) * ellipk(np.sqrt(parameter))
    assert abs(actual - wrong_modulus_convention) > 1e-3


def test_nonlinear_numerical_period_matches_analytical_for_libration(
    research_results: phase26.PendulumInitialVelocityResearch,
) -> None:
    analyzed = [
        case for case in research_results.cases if case.status is CaseStatus.ANALYZED
    ]
    assert analyzed
    for case in analyzed:
        assert case.nonlinear_numerical_period_s == pytest.approx(
            case.nonlinear_analytical_period_s, abs=2e-5
        )
        assert case.period_estimate is not None
        assert case.period_estimate.crossing_count >= 2


def test_zero_velocity_matches_phase25_and_small_amplitude_linear_limit() -> None:
    phase26_period = nonlinear_period_reference(
        0.5, 0.0, length_m=1.0, gravity_m_per_s2=9.81
    )
    assert phase26_period == pytest.approx(
        phase25_period(0.5, length_m=1.0, gravity_m_per_s2=9.81), abs=1e-14
    )
    tiny_period = nonlinear_period_reference(
        1e-5, 0.0, length_m=1.0, gravity_m_per_s2=9.81
    )
    assert tiny_period == pytest.approx(
        linear_period_reference(length_m=1.0, gravity_m_per_s2=9.81), rel=1e-10
    )


def test_phase_origin_records_initial_direction_and_phase_lag_sign(
    research_results: phase26.PendulumInitialVelocityResearch,
) -> None:
    positive = next(
        case
        for case in research_results.cases
        if case.case_id == "grid_theta_0.5_omega_+1"
    )
    negative = next(
        case
        for case in research_results.cases
        if case.case_id == "grid_theta_0.5_omega_-1"
    )
    assert positive.initial_linear_phase_rad is not None
    assert negative.initial_linear_phase_rad is not None
    assert positive.initial_linear_phase_rad < 0.0
    assert negative.initial_linear_phase_rad > 0.0
    assert positive.windows[0].signed_accumulated_phase_difference_rad > 0.0
    assert negative.windows[0].signed_accumulated_phase_difference_rad > 0.0


def test_period_estimate_reports_incomplete_trajectory_explicitly() -> None:
    short = simulate_damped_pendulum(
        1.0,
        1.0,
        0.0,
        9.81,
        0.5,
        0.0,
        0.1,
        num_points=21,
    )
    estimate = estimate_nonlinear_period(short)
    assert estimate.period_s is None
    assert estimate.status == "insufficient_complete_cycles"
    assert "before a complete period" in estimate.explanation


def test_zero_angle_initial_velocity_period_estimator_counts_initial_direction() -> (
    None
):
    trajectory = simulate_damped_pendulum(
        1.0,
        1.0,
        0.0,
        9.81,
        0.0,
        0.5,
        4.2,
        num_points=4201,
    )
    estimate = estimate_nonlinear_period(trajectory)
    assert estimate.crossing_direction == "positive"
    assert estimate.crossing_count >= 2
    assert estimate.period_s == pytest.approx(
        nonlinear_period_reference(0.0, 0.5, length_m=1.0, gravity_m_per_s2=9.81),
        abs=1e-5,
    )


def test_solver_failure_is_retained_as_failure_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_simulation(*args: object, **kwargs: object) -> None:
        raise IntegrationError("injected solver failure")

    monkeypatch.setattr(phase26, "simulate_damped_pendulum", fail_simulation)
    result = phase26.analyze_initial_condition(
        InitialConditionSpec(
            case_id="solver_failure_test",
            initial_angle_rad=0.5,
            initial_angular_velocity_rad_per_s=0.0,
        ),
        sweep_index=0,
        mass_kg=1.0,
        length_m=1.0,
        gravity_m_per_s2=9.81,
        criteria=AccuracyCriteria(),
    )
    assert result.status is CaseStatus.SOLVER_FAILURE
    assert result.failure_reason == "injected solver failure"
    assert result.nonlinear_numerical_period_s is None


def test_thresholds_are_observable_specific_and_conditional(
    research_results: phase26.PendulumInitialVelocityResearch,
) -> None:
    by_observable = {
        (item.observable, item.initial_angular_velocity_rad_per_s): item
        for item in research_results.threshold_slices
    }
    period_rest = by_observable[("relative_period_error", 0.0)]
    angle_rest = by_observable[("sampled_angle_error", 0.0)]
    phase_rest = by_observable[("absolute_accumulated_phase_difference", 0.0)]
    assert period_rest.largest_tested_passing_absolute_initial_angle_rad is not None
    assert angle_rest.largest_tested_passing_absolute_initial_angle_rad is not None
    assert phase_rest.largest_tested_passing_absolute_initial_angle_rad is not None
    assert period_rest.largest_tested_passing_absolute_initial_angle_rad != (
        phase_rest.largest_tested_passing_absolute_initial_angle_rad
    )
    assert period_rest.next_higher_point_failed is True
    assert "conditional" in period_rest.interpretation
    assert any(
        case.status is CaseStatus.ROTATION_EXCLUDED for case in research_results.cases
    )
    assert any(
        case.status is CaseStatus.SEPARATRIX_EXCLUDED for case in research_results.cases
    )
    assert any(
        case.status is CaseStatus.INVALID_INPUT for case in research_results.cases
    )


def test_tolerance_and_grid_sensitivity_include_small_moderate_and_near_boundary(
    research_results: phase26.PendulumInitialVelocityResearch,
) -> None:
    assert {item.comparison for item in research_results.numerical_sensitivities} == {
        "looser_tolerances",
        "coarser_output_grid",
    }
    assert {item.case_id for item in research_results.numerical_sensitivities} == {
        "grid_theta_0.1_omega_+1",
        "grid_theta_0.5_omega_+3",
        "near_boundary_3_12_rest",
    }
    near_tolerance = next(
        item
        for item in research_results.numerical_sensitivities
        if item.case_id == "near_boundary_3_12_rest"
        and item.comparison == "looser_tolerances"
    )
    small_tolerance = next(
        item
        for item in research_results.numerical_sensitivities
        if item.case_id == "grid_theta_0.1_omega_+1"
        and item.comparison == "looser_tolerances"
    )
    assert near_tolerance.period_estimate_difference_s is not None
    assert small_tolerance.period_estimate_difference_s is not None
    assert abs(near_tolerance.period_estimate_difference_s) > abs(
        small_tolerance.period_estimate_difference_s
    )


def test_sweep_order_and_result_provenance_are_preserved(
    research_results: phase26.PendulumInitialVelocityResearch,
) -> None:
    assert tuple(case.case_id for case in research_results.cases) == tuple(
        spec.case_id for spec in research_results.sweep_order
    )
    assert tuple(case.sweep_index for case in research_results.cases) == tuple(
        range(len(research_results.cases))
    )
    for case, spec in zip(
        research_results.cases, research_results.sweep_order, strict=True
    ):
        assert case.initial_angle_rad == spec.initial_angle_rad
        assert case.initial_angular_velocity_rad_per_s == (
            spec.initial_angular_velocity_rad_per_s
        )


def test_report_table_is_complete_and_has_consistent_columns(
    research_results: phase26.PendulumInitialVelocityResearch,
) -> None:
    report = render_initial_velocity_report(research_results)
    assert "Initial linear phase (rad)" in report
    assert "separatrix_period_not_finite" in report
    assert "rotation_excluded_from_libration_period" in report
    lines = report.splitlines()
    heading = lines.index("## Computed initial-condition sweep")
    header = lines[heading + 2]
    separator = lines[heading + 3]
    first_data_row = lines[heading + 4]
    assert header.count("|") == separator.count("|") == first_data_row.count("|")
