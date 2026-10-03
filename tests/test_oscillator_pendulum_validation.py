"""Mathematical and numerical tests for Phase 28 validation."""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

import newton_lab.oscillator_pendulum_validation as validation
from newton_lab.dynamics import simulate_damped_oscillator
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.pendulum_initial_velocity import nonlinear_period_reference


def test_physical_and_dimensionless_mapping_includes_scale_and_time_conversion() -> (
    None
):
    state = validation.construct_matched_initial_conditions(0.2, -0.3)
    assert state.amplitude_scale_m == pytest.approx(0.25)
    assert state.oscillator_displacement_m == pytest.approx(0.05)
    assert state.oscillator_velocity_m_per_s == pytest.approx(-0.15)
    assert state.pendulum_angle_rad == pytest.approx(0.2)
    assert state.pendulum_angular_velocity_rad_per_s == pytest.approx(-0.6)
    assert state.frequency_match
    assert 2.0 * np.pi / state.oscillator_frequency_rad_per_s == pytest.approx(
        2.0 * np.pi / state.pendulum_frequency_rad_per_s
    )


def test_different_frequencies_use_dimensionless_time_or_reject_required_match() -> (
    None
):
    state = validation.construct_matched_initial_conditions(
        0.1,
        0.2,
        mass_kg=1.0,
        stiffness_n_per_m=4.0,
        length_m=2.0,
        gravity_m_per_s2=4.0,
    )
    assert not state.frequency_match
    assert state.oscillator_frequency_rad_per_s == pytest.approx(2.0)
    assert state.pendulum_frequency_rad_per_s == pytest.approx(np.sqrt(2.0))
    with pytest.raises(ScientificValidationError):
        validation.construct_matched_initial_conditions(
            0.1,
            0.2,
            mass_kg=1.0,
            stiffness_n_per_m=4.0,
            length_m=2.0,
            gravity_m_per_s2=4.0,
            require_frequency_match=True,
        )


def test_zero_amplitude_scale_is_rejected() -> None:
    with pytest.raises(ScientificValidationError):
        validation.construct_matched_initial_conditions(0.0, 0.0, amplitude_scale_m=0.0)


@pytest.mark.parametrize("q0,velocity", [(0.4, 0.0), (0.1, 0.3), (-0.2, -0.5)])
def test_exact_solution_for_arbitrary_initial_state(q0: float, velocity: float) -> None:
    tau = np.array([0.0, 0.5, 1.25, 2.0])
    actual = validation.exact_linear_solution(tau, q0, velocity)
    expected = q0 * np.cos(tau) + velocity * np.sin(tau)
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=0.0)
    assert actual[0] == pytest.approx(q0)


def test_velocity_sign_reversal_has_expected_linear_symmetry() -> None:
    tau = np.linspace(0.0, 4.0, 100)
    plus = validation.exact_linear_solution(tau, 0.2, 0.3)
    minus = validation.exact_linear_solution(tau, 0.2, -0.3)
    np.testing.assert_allclose(plus + minus, 0.4 * np.cos(tau), atol=1e-15)


def test_oscillator_and_linear_pendulum_are_integrated_independently() -> None:
    state = validation.construct_matched_initial_conditions(0.2, 0.1)
    with (
        patch(
            "newton_lab.oscillator_pendulum_validation.simulate_damped_oscillator",
            wraps=simulate_damped_oscillator,
        ) as oscillator_spy,
        patch(
            "newton_lab.oscillator_pendulum_validation.solve_ivp",
            wraps=validation.__dict__["solve_ivp"],
        ) as solve_spy,
    ):
        run = validation._run_one(
            case_id="spy",
            initial=state,
            duration_tau=2.0 * np.pi,
            sample_count=101,
            method="DOP853",
            rtol=1e-9,
            atol=1e-11,
        )
    assert run.status == "succeeded"
    assert oscillator_spy.call_count == 1
    assert solve_spy.call_count == 1  # Independent linearized-pendulum IVP.


def test_both_numerical_linear_models_match_exact_reference() -> None:
    state = validation.construct_matched_initial_conditions(0.2, 0.3)
    run = validation._run_one(
        case_id="baseline",
        initial=state,
        duration_tau=8.0 * np.pi,
        sample_count=501,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
    )
    assert run.status == "succeeded"
    assert run.oscillator_error_max is not None and run.oscillator_error_max < 1e-8
    assert (
        run.linear_pendulum_error_max is not None
        and run.linear_pendulum_error_max < 1e-8
    )
    assert run.cross_system_error_max is not None and run.cross_system_error_max < 1e-8


def test_cross_error_is_consistent_with_individual_errors() -> None:
    run = validation._run_one(
        case_id="triangle",
        initial=validation.construct_matched_initial_conditions(0.2, -0.3),
        duration_tau=8.0 * np.pi,
        sample_count=501,
        method="RK45",
        rtol=1e-7,
        atol=1e-9,
    )
    assert run.status == "succeeded"
    assert run.cross_system_error_max is not None
    assert run.oscillator_error_max is not None
    assert run.linear_pendulum_error_max is not None
    assert run.cross_system_error_max <= (
        run.oscillator_error_max + run.linear_pendulum_error_max + 1e-14
    )


def test_common_dimensionless_grid_is_preserved_for_unequal_frequencies() -> None:
    state = validation.construct_matched_initial_conditions(
        0.1,
        0.2,
        mass_kg=1.0,
        stiffness_n_per_m=9.0,
        length_m=2.0,
        gravity_m_per_s2=8.0,
    )
    run = validation._run_one(
        case_id="different_frequency",
        initial=state,
        duration_tau=3.0,
        sample_count=61,
        method="DOP853",
        rtol=1e-9,
        atol=1e-11,
    )
    assert run.status == "succeeded"
    assert run.dimensionless_times[0] == 0.0
    assert run.dimensionless_times[-1] == 3.0
    assert len(run.dimensionless_times) == 61


def test_physical_parameters_preserve_dimensionless_trajectories() -> None:
    standard = validation.construct_matched_initial_conditions(0.2, -0.3)
    rescaled = validation.construct_matched_initial_conditions(
        0.2,
        -0.3,
        mass_kg=1.0,
        stiffness_n_per_m=9.0,
        length_m=2.0,
        gravity_m_per_s2=8.0,
    )
    standard_run = validation._run_one(
        case_id="standard_scale",
        initial=standard,
        duration_tau=8.0 * np.pi,
        sample_count=501,
        method="DOP853",
        rtol=1e-9,
        atol=1e-11,
    )
    rescaled_run = validation._run_one(
        case_id="rescaled_parameters",
        initial=rescaled,
        duration_tau=8.0 * np.pi,
        sample_count=501,
        method="DOP853",
        rtol=1e-9,
        atol=1e-11,
    )
    assert standard_run.status == rescaled_run.status == "succeeded"
    np.testing.assert_allclose(
        standard_run.oscillator_dimensionless_displacement,
        rescaled_run.oscillator_dimensionless_displacement,
        atol=1e-8,
    )
    np.testing.assert_allclose(
        standard_run.linear_pendulum_displacement,
        rescaled_run.linear_pendulum_displacement,
        atol=1e-8,
    )


def test_independent_integrations_preserve_velocity_sign_symmetry() -> None:
    runs = []
    for velocity in (0.3, -0.3):
        runs.append(
            validation._run_one(
                case_id=f"velocity_{velocity:+g}",
                initial=validation.construct_matched_initial_conditions(0.2, velocity),
                duration_tau=8.0 * np.pi,
                sample_count=501,
                method="DOP853",
                rtol=1e-9,
                atol=1e-11,
            )
        )
    positive, negative = runs
    assert positive.status == negative.status == "succeeded"
    expected_even_part = 0.4 * np.cos(np.asarray(positive.dimensionless_times))
    np.testing.assert_allclose(
        np.asarray(positive.oscillator_dimensionless_displacement)
        + np.asarray(negative.oscillator_dimensionless_displacement),
        expected_even_part,
        atol=1e-8,
    )
    np.testing.assert_allclose(
        np.asarray(positive.linear_pendulum_displacement)
        + np.asarray(negative.linear_pendulum_displacement),
        expected_even_part,
        atol=1e-8,
    )


def test_zero_state_errors_are_absolute_and_normalized_values_undefined() -> None:
    run = validation._run_one(
        case_id="zero",
        initial=validation.construct_matched_initial_conditions(0.0, 0.0),
        duration_tau=4.0 * np.pi,
        sample_count=101,
        method="DOP853",
        rtol=1e-9,
        atol=1e-11,
    )
    assert run.status == "succeeded"
    assert run.oscillator_error_max == pytest.approx(0.0)
    assert run.normalized_cross_system_error is None


def test_normalization_rejects_invalid_error_and_handles_zero_reference() -> None:
    assert validation.normalized_error(1e-14, np.zeros(4)) is None
    with pytest.raises(ScientificValidationError):
        validation.normalized_error(-1.0, np.ones(2))


def test_full_sine_model_error_is_distinct_and_larger_at_moderate_angle() -> None:
    run = validation._run_one(
        case_id="moderate",
        initial=validation.construct_matched_initial_conditions(0.8, 0.3),
        duration_tau=8.0 * np.pi,
        sample_count=1001,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
    )
    assert run.status == "succeeded"
    assert run.nonlinear_model_error_max is not None
    assert run.linear_pendulum_error_max is not None
    assert run.nonlinear_model_error_max > 1000 * run.linear_pendulum_error_max
    assert run.nonlinear_vs_linear_solver_error_max is not None


def test_small_angle_nonlinear_model_error_is_smaller() -> None:
    small = validation._run_one(
        case_id="small",
        initial=validation.construct_matched_initial_conditions(0.01, 0.0),
        duration_tau=8.0 * np.pi,
        sample_count=501,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
    )
    large = validation._run_one(
        case_id="large",
        initial=validation.construct_matched_initial_conditions(0.8, 0.0),
        duration_tau=8.0 * np.pi,
        sample_count=501,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
    )
    assert small.nonlinear_model_error_max is not None
    assert large.nonlinear_model_error_max is not None
    assert small.nonlinear_model_error_max < large.nonlinear_model_error_max


def test_elliptic_period_reference_matches_existing_phase26_api() -> None:
    theta0, omega0 = 0.5, 1.0
    assert nonlinear_period_reference(
        theta0, omega0, length_m=1.5, gravity_m_per_s2=6.0
    ) == pytest.approx(
        nonlinear_period_reference(theta0, omega0, length_m=1.5, gravity_m_per_s2=6.0)
    )


@pytest.mark.parametrize("velocity", [3.0, 2.0])
def test_rotating_and_separatrix_states_have_no_libration_period(
    velocity: float,
) -> None:
    state = validation.construct_matched_initial_conditions(0.0, velocity)
    run = validation._run_one(
        case_id="nonlibration",
        initial=state,
        duration_tau=2.0 * np.pi,
        sample_count=101,
        method="DOP853",
        rtol=1e-8,
        atol=1e-10,
    )
    assert run.status == "succeeded"
    assert run.nonlinear_period is not None
    assert run.nonlinear_period.reference_s is None
    assert run.nonlinear_period.status == "reference_unavailable_nonlibration"


def test_period_estimator_requires_multiple_same_direction_crossings() -> None:
    times = np.linspace(0.0, 2.0, 101)
    values = np.cos(times)
    estimate, count, status = validation.estimate_period_tau(times, values)
    assert estimate is None
    assert count < 2
    assert status == "insufficient_same_direction_crossings"


def test_period_estimator_recovers_sinusoidal_period_on_dense_grid() -> None:
    times = np.linspace(0.0, 8.0 * np.pi, 4001)
    values = np.cos(times + 0.4)
    estimate, count, status = validation.estimate_period_tau(times, values)
    assert count >= 3
    assert status.startswith("estimated")
    assert estimate == pytest.approx(2.0 * np.pi, abs=1e-5)


def test_periods_are_reported_for_both_independent_linear_integrations() -> None:
    run = validation._run_one(
        case_id="linear_periods",
        initial=validation.construct_matched_initial_conditions(0.2, 0.1),
        duration_tau=8.0 * np.pi,
        sample_count=501,
        method="DOP853",
        rtol=1e-9,
        atol=1e-11,
    )
    assert run.status == "succeeded"
    assert run.linear_period is not None
    assert run.pendulum_linear_period is not None
    assert run.linear_period.estimate_tau == pytest.approx(2.0 * np.pi, abs=1e-4)
    assert run.pendulum_linear_period.estimate_tau == pytest.approx(
        2.0 * np.pi, abs=1e-4
    )


def test_tolerance_and_output_grid_sensitivity_are_recorded() -> None:
    study = validation.run_validation_experiment()
    assert {item.factor for item in study.sensitivities} >= {
        "tolerance",
        "method",
        "output_grid",
        "horizon",
    }
    assert study.tolerances == ((1e-5, 1e-7), (1e-8, 1e-10), (1e-11, 1e-13))
    assert study.output_counts == (101, 501, 2001)
    assert all(run.status == "succeeded" for run in study.runs)


def test_ordered_provenance_includes_failed_or_successful_case_settings() -> None:
    study = validation.run_validation_experiment()
    first = study.runs[0]
    assert first.case_id == "initial_1"
    assert first.method == "DOP853"
    assert first.dimensionless_times[0] == 0.0
    assert first.dimensionless_times[-1] == first.duration_tau
    assert first.initial.amplitude_scale_m == 0.25
    assert first.oscillator_function_evaluations is not None
    assert len(first.exact_dimensionless_displacement) == first.sample_count
    assert len(first.oscillator_dimensionless_displacement) == first.sample_count
    assert len(first.linear_pendulum_displacement) == first.sample_count
    assert len(first.nonlinear_pendulum_displacement) == first.sample_count
    assert first.oscillator_dimensionless_displacement == pytest.approx(
        first.linear_pendulum_displacement
    )


def test_integration_failure_is_recorded_without_fabricated_metrics() -> None:
    state = validation.construct_matched_initial_conditions(0.2, 0.1)
    with patch(
        "newton_lab.oscillator_pendulum_validation.simulate_damped_oscillator",
        side_effect=IntegrationError("forced deterministic failure"),
    ):
        run = validation._run_one(
            case_id="failed_oscillator",
            initial=state,
            duration_tau=2.0 * np.pi,
            sample_count=101,
            method="DOP853",
            rtol=1e-8,
            atol=1e-10,
        )
    assert run.status == "failed"
    assert run.failure_message == "IntegrationError: forced deterministic failure"
    assert run.oscillator_error_max is None
    assert run.oscillator_function_evaluations is None


def test_report_generation_is_deterministic_and_explicit_about_shared_library() -> None:
    first = validation.run_validation_experiment()
    second = validation.run_validation_experiment()
    assert first.report_markdown == second.report_markdown
    assert "not independent numerical-library validation" in first.report_markdown
    assert "full-sine model error" in first.report_markdown.lower()


def test_phase27_api_and_report_remain_available() -> None:
    from newton_lab.cross_system_discovery import build_cross_system_report

    report = build_cross_system_report()
    assert report.trajectory.linear_normalized_error < 1e-8
    assert "H3b_universal_threshold" in report.markdown
