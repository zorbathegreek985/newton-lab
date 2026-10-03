"""Tests for the nonlinear pendulum's Phase 11 experiment adapter."""

import numpy as np
import pytest

from newton_lab.experiments import (
    AnalyticalReferenceRequest,
    ExperimentFailureKind,
    ExperimentRecordStatus,
    ExperimentRunStatus,
    SensitivityFailure,
    SensitivityResult,
    make_pendulum_experiment,
    run_experiment,
    run_sensitivity_analysis,
)
from newton_lab.pendulum import simulate_damped_pendulum
from newton_lab.simulation import ODESolverConfiguration, OutputSampling


def test_adapter_uses_existing_nonlinear_pendulum_and_extracts_metrics() -> None:
    specification = make_pendulum_experiment(
        mass_kg=0.4,
        length_m=0.8,
        damping_coefficient_kg_m2_per_s=0.02,
        gravity_m_per_s2=9.81,
        initial_angle_rad=0.7,
        initial_angular_velocity_rad_per_s=0.1,
        duration_s=2.0,
        sweep_values=None,
        metrics=(
            "final_angle_rad",
            "final_angular_velocity_rad_per_s",
            "maximum_absolute_angle_rad",
            "maximum_absolute_angular_velocity_rad_per_s",
        ),
        output_sampling=OutputSampling(point_count=121),
    )
    run = run_experiment(specification).runs[0]
    direct = simulate_damped_pendulum(
        mass_kg=0.4,
        length_m=0.8,
        damping_coefficient_kg_m2_per_s=0.02,
        gravity_m_per_s2=9.81,
        initial_angle_rad=0.7,
        initial_angular_velocity_rad_per_s=0.1,
        duration_s=2.0,
        num_points=121,
    )

    assert run.status == ExperimentRunStatus.SUCCEEDED
    metrics = {metric.metric_id: metric for metric in run.metrics}
    assert metrics["final_angle_rad"].unit == "rad"
    assert metrics["final_angle_rad"].value == pytest.approx(direct.angle_rad[-1])
    assert metrics["final_angular_velocity_rad_per_s"].unit == "rad/s"
    assert metrics["final_angular_velocity_rad_per_s"].value == pytest.approx(
        direct.angular_velocity_rad_per_s[-1]
    )
    assert metrics["maximum_absolute_angle_rad"].value == pytest.approx(
        np.max(np.abs(direct.angle_rad))
    )
    assert metrics[
        "maximum_absolute_angular_velocity_rad_per_s"
    ].value == pytest.approx(np.max(np.abs(direct.angular_velocity_rad_per_s)))
    assert run.reproducibility.parameters == specification.parameters
    assert run.reproducibility.output_point_count == 121
    assert run.function_evaluations == direct.function_evaluations


def test_angle_sweep_preserves_order_and_does_not_mutate_shared_inputs() -> None:
    specification = make_pendulum_experiment(
        initial_angle_rad=0.2,
        duration_s=1.5,
        sweep_parameter_name="initial_angle_rad",
        sweep_values=(0.5, 0.05, 0.2),
        metrics=("final_angle_rad", "maximum_absolute_angle_rad"),
        output_sampling=OutputSampling(point_count=81),
    )
    original_parameters = specification.parameters

    first = run_experiment(specification)
    second = run_experiment(specification)

    assert first.status == ExperimentRecordStatus.SUCCEEDED
    assert [run.swept_parameter_value for run in first.runs] == [0.5, 0.05, 0.2]
    assert all(
        run.parameters == other.parameters
        for run, other in zip(first.runs, second.runs, strict=True)
    )
    assert [run.metrics for run in first.runs] == [run.metrics for run in second.runs]
    assert specification.parameters == original_parameters
    assert all(run.solver_method == "RK45" for run in first.runs)


def test_expected_invalid_pendulum_run_is_recorded_and_later_runs_continue() -> None:
    specification = make_pendulum_experiment(
        sweep_parameter_name="length_m",
        sweep_values=(0.75, 0.0, 1.0),
        duration_s=0.8,
        metrics=("final_angle_rad",),
        output_sampling=OutputSampling(point_count=31),
    )

    record = run_experiment(specification)

    assert record.status == ExperimentRecordStatus.COMPLETED_WITH_FAILURES
    assert [run.swept_parameter_value for run in record.runs] == [0.75, 0.0, 1.0]
    assert [run.status for run in record.runs] == [
        ExperimentRunStatus.SUCCEEDED,
        ExperimentRunStatus.FAILED,
        ExperimentRunStatus.SUCCEEDED,
    ]
    assert record.runs[1].failure_kind == ExperimentFailureKind.MODEL_VALIDATION
    assert "length_m must be positive" in record.runs[1].message


def test_pendulum_sensitivity_uses_base_and_symmetric_parameter_runs() -> None:
    specification = make_pendulum_experiment(
        initial_angle_rad=0.3,
        damping_coefficient_kg_m2_per_s=0.01,
        sweep_values=None,
        duration_s=1.2,
        metrics=("final_angle_rad",),
        output_sampling=OutputSampling(point_count=81),
        solver_configuration=ODESolverConfiguration(
            relative_tolerance=1e-9, absolute_tolerance=1e-11
        ),
    )

    result = run_sensitivity_analysis(
        specification,
        parameter_name="initial_angle_rad",
        metric_id="final_angle_rad",
        perturbation=0.001,
    )

    assert isinstance(result, SensitivityResult)
    assert result.parameter_unit == "rad"
    assert result.metric_unit == "rad"
    assert result.derivative_unit == "rad/rad"
    assert result.lower_parameter_value == pytest.approx(0.299)
    assert result.upper_parameter_value == pytest.approx(0.301)
    assert tuple(run.run_index for run in result.runs) == (0, 1, 2)
    assert np.isfinite(result.derivative_estimate)


def test_invalid_perturbed_pendulum_parameter_is_reported() -> None:
    specification = make_pendulum_experiment(
        mass_kg=0.2,
        sweep_values=None,
        duration_s=0.8,
        metrics=("final_angle_rad",),
        output_sampling=OutputSampling(point_count=31),
    )

    result = run_sensitivity_analysis(
        specification,
        parameter_name="mass_kg",
        metric_id="final_angle_rad",
        perturbation=0.3,
    )

    assert isinstance(result, SensitivityFailure)
    assert result.failed_run_indices == (1,)
    assert result.runs[1].failure_kind == ExperimentFailureKind.MODEL_VALIDATION
    assert "mass_kg must be positive" in result.runs[1].message


def test_small_angle_reference_compares_same_final_angle_and_states_limit() -> None:
    specification = make_pendulum_experiment(
        damping_coefficient_kg_m2_per_s=0.0,
        initial_angle_rad=0.05,
        duration_s=1.0,
        sweep_values=None,
        metrics=("final_angle_rad",),
        output_sampling=OutputSampling(point_count=101),
        solver_configuration=ODESolverConfiguration(
            relative_tolerance=1e-10, absolute_tolerance=1e-12
        ),
        analytical_reference=AnalyticalReferenceRequest(
            metric_id="final_angle_rad", absolute_tolerance=1e-4
        ),
    )

    run = run_experiment(specification).runs[0]

    assert run.reference_comparison is not None
    assert run.reference_comparison.status == "passed"
    assert run.reference_comparison.unit == "rad"
    assert run.reference_comparison.absolute_error is not None
    assert run.reference_comparison.reference_description.startswith(
        "Undamped small-angle linearized pendulum"
    )
    assert "not an exact solution" in run.reference_comparison.reference_description


def test_larger_angle_can_exceed_small_angle_reference_tolerance() -> None:
    specification = make_pendulum_experiment(
        damping_coefficient_kg_m2_per_s=0.0,
        initial_angle_rad=0.8,
        duration_s=2.0,
        sweep_values=None,
        metrics=("final_angle_rad",),
        output_sampling=OutputSampling(point_count=101),
        analytical_reference=AnalyticalReferenceRequest(
            metric_id="final_angle_rad", absolute_tolerance=1e-5
        ),
    )

    comparison = run_experiment(specification).runs[0].reference_comparison

    assert comparison is not None
    assert comparison.status == "failed"
    assert comparison.absolute_error is not None
    assert comparison.absolute_error > comparison.absolute_tolerance


def test_damped_pendulum_has_no_undamped_small_angle_reference() -> None:
    specification = make_pendulum_experiment(
        damping_coefficient_kg_m2_per_s=0.01,
        sweep_values=None,
        metrics=("final_angle_rad",),
        analytical_reference=AnalyticalReferenceRequest(
            metric_id="final_angle_rad", absolute_tolerance=1e-5
        ),
    )

    comparison = run_experiment(specification).runs[0].reference_comparison

    assert comparison is not None
    assert comparison.status == "unavailable"
