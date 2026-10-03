"""Tests for ordered experiments, scalar metrics, and finite differences."""

import numpy as np
import pytest
from pydantic import ValidationError

from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.experiments import (
    AnalyticalReferenceRequest,
    ExperimentFailureKind,
    ExperimentRecordStatus,
    ExperimentRunStatus,
    ExperimentSpecification,
    FiniteDifferenceEstimate,
    MetricValue,
    SensitivityFailure,
    SensitivityResult,
    estimate_central_difference,
    make_oscillator_experiment,
    run_experiment,
    run_sensitivity_analysis,
)
from newton_lab.simulation import ODESolverConfiguration, OutputSampling


def _spec(
    *,
    initial_displacement_m: float = 0.1,
    initial_velocity_m_per_s: float = 0.0,
    duration_s: float = 2.0,
    metrics: tuple[str, ...] = ("final_displacement_m",),
) -> ExperimentSpecification:
    return make_oscillator_experiment(
        experiment_id="test_oscillator_run",
        sweep_values=None,
        initial_displacement_m=initial_displacement_m,
        initial_velocity_m_per_s=initial_velocity_m_per_s,
        duration_s=duration_s,
        metrics=metrics,
        output_sampling=OutputSampling(point_count=101),
    )


def test_sweep_preserves_value_order_and_records_each_point() -> None:
    specification = make_oscillator_experiment(
        sweep_values=(0.3, 0.0, 0.1),
        duration_s=1.0,
        output_sampling=OutputSampling(point_count=51),
        metrics=("final_displacement_m",),
    )

    record = run_experiment(specification)

    assert record.status == ExperimentRecordStatus.SUCCEEDED
    assert [run.swept_parameter_value for run in record.runs] == [0.3, 0.0, 0.1]
    assert [run.run_index for run in record.runs] == [0, 1, 2]
    expected_parameter_order = tuple(item.name for item in specification.parameters)
    assert all(
        tuple(item.name for item in run.parameters) == expected_parameter_order
        for run in record.runs
    )
    assert all(run.status == ExperimentRunStatus.SUCCEEDED for run in record.runs)
    assert all(run.function_evaluations is not None for run in record.runs)


def test_sweep_continues_after_expected_invalid_parameter_failure() -> None:
    specification = make_oscillator_experiment(
        sweep_values=(0.0, -0.1, 0.2),
        duration_s=1.0,
        output_sampling=OutputSampling(point_count=31),
        metrics=("final_displacement_m",),
    )

    record = run_experiment(specification)

    assert record.status == ExperimentRecordStatus.COMPLETED_WITH_FAILURES
    assert [run.swept_parameter_value for run in record.runs] == [0.0, -0.1, 0.2]
    assert [run.status for run in record.runs] == [
        ExperimentRunStatus.SUCCEEDED,
        ExperimentRunStatus.FAILED,
        ExperimentRunStatus.SUCCEEDED,
    ]
    failed = record.runs[1]
    assert failed.failure_kind == ExperimentFailureKind.MODEL_VALIDATION
    assert "non-negative" in failed.message


def test_oscillator_metrics_are_defined_over_requested_samples() -> None:
    specification = _spec(
        initial_displacement_m=0.15,
        metrics=(
            "maximum_absolute_displacement_m",
            "final_displacement_m",
            "final_velocity_m_per_s",
            "maximum_absolute_velocity_m_per_s",
        ),
    )
    record = run_experiment(specification)
    run = record.runs[0]

    assert run.status == ExperimentRunStatus.SUCCEEDED
    metrics = {metric.metric_id: metric for metric in run.metrics}
    assert metrics["maximum_absolute_displacement_m"].unit == "m"
    assert metrics["maximum_absolute_displacement_m"].value == pytest.approx(0.15)
    assert metrics["final_displacement_m"].unit == "m"
    assert metrics["final_velocity_m_per_s"].unit == "m/s"
    assert metrics["maximum_absolute_velocity_m_per_s"].unit == "m/s"


def test_repeated_experiments_are_deterministic_and_keep_inputs() -> None:
    specification = make_oscillator_experiment(
        sweep_values=(0.0, 0.2),
        duration_s=1.0,
        output_sampling=OutputSampling(point_count=41),
        metrics=("final_displacement_m",),
    )
    original = specification.parameters
    first = run_experiment(specification)
    second = run_experiment(specification)

    assert specification.parameters == original
    assert [run.metrics for run in first.runs] == [run.metrics for run in second.runs]
    assert [run.reproducibility for run in first.runs] == [
        run.reproducibility for run in second.runs
    ]
    assert first.runs[0].reproducibility.random_seed is None
    assert first.runs[0].reproducibility.output_point_count == 41


def test_undamped_analytical_reference_comparison() -> None:
    specification = make_oscillator_experiment(
        damping_coefficient_kg_per_s=0.0,
        sweep_values=None,
        initial_displacement_m=0.2,
        initial_velocity_m_per_s=0.1,
        duration_s=1.7,
        metrics=("final_displacement_m",),
        solver_configuration=ODESolverConfiguration(
            relative_tolerance=1e-10, absolute_tolerance=1e-12
        ),
        output_sampling=OutputSampling(point_count=71),
        analytical_reference=AnalyticalReferenceRequest(
            metric_id="final_displacement_m", absolute_tolerance=1e-8
        ),
    )

    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.SUCCEEDED
    assert run.reference_comparison is not None
    assert run.reference_comparison.status == "passed"
    assert run.reference_comparison.absolute_error is not None
    assert run.reference_comparison.absolute_error <= 1e-8


def test_reference_unavailable_does_not_claim_physical_failure() -> None:
    specification = make_oscillator_experiment(
        damping_coefficient_kg_per_s=0.1,
        sweep_values=None,
        duration_s=1.0,
        metrics=("final_displacement_m",),
        analytical_reference=AnalyticalReferenceRequest(
            metric_id="final_displacement_m", absolute_tolerance=1e-6
        ),
    )
    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.SUCCEEDED
    assert run.reference_comparison is not None
    assert run.reference_comparison.status == "unavailable"


def test_analytical_reference_reports_tolerance_failure() -> None:
    specification = make_oscillator_experiment(
        damping_coefficient_kg_per_s=0.0,
        sweep_values=None,
        duration_s=1.0,
        metrics=("final_displacement_m",),
        analytical_reference=AnalyticalReferenceRequest(
            metric_id="final_displacement_m", absolute_tolerance=1e-16
        ),
    )
    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.SUCCEEDED
    assert run.reference_comparison is not None
    assert run.reference_comparison.status == "failed"
    assert run.reference_comparison.absolute_error is not None
    assert run.reference_comparison.absolute_error > 1e-16


def test_solver_failures_are_recorded_and_sweep_continues(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import newton_lab.experiments as experiments

    calls = 0

    def fail_solver(*_args: object, **_kwargs: object) -> None:
        nonlocal calls
        calls += 1
        raise IntegrationError("forced integrator failure")

    monkeypatch.setattr(experiments, "simulate_damped_oscillator", fail_solver)
    specification = make_oscillator_experiment(
        sweep_values=(0.0, 0.1, 0.2),
        duration_s=1.0,
        metrics=("final_displacement_m",),
    )

    record = run_experiment(specification)

    assert calls == 3
    assert record.status == ExperimentRecordStatus.COMPLETED_WITH_FAILURES
    assert [run.failure_kind for run in record.runs] == [
        ExperimentFailureKind.SOLVER,
        ExperimentFailureKind.SOLVER,
        ExperimentFailureKind.SOLVER,
    ]
    assert all(run.solver_success is False for run in record.runs)


def test_metric_extraction_failure_is_distinct_from_solver_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import newton_lab.experiments as experiments
    from newton_lab.dynamics import DampedOscillatorResult

    def nonfinite_result(*_args: object, **_kwargs: object) -> DampedOscillatorResult:
        return DampedOscillatorResult(
            time_s=np.asarray([0.0, 1.0]),
            displacement_m=np.asarray([0.1, float("nan")]),
            velocity_m_per_s=np.asarray([0.0, 0.0]),
            mass_kg=1.0,
            damping_coefficient_kg_per_s=0.1,
            stiffness_n_per_m=4.0,
            initial_displacement_m=0.1,
            initial_velocity_m_per_s=0.0,
            method="RK45",
            relative_tolerance=1e-8,
            absolute_tolerance=1e-10,
            function_evaluations=5,
        )

    monkeypatch.setattr(experiments, "simulate_damped_oscillator", nonfinite_result)
    run = run_experiment(
        make_oscillator_experiment(sweep_values=None, metrics=("final_displacement_m",))
    ).runs[0]

    assert run.status == ExperimentRunStatus.FAILED
    assert run.failure_kind == ExperimentFailureKind.METRIC_EXTRACTION
    assert run.solver_success is True


def test_finite_difference_matches_known_quadratic_derivative() -> None:
    estimate = estimate_central_difference(lambda value: value**2, 3.0, 0.01)

    assert isinstance(estimate, FiniteDifferenceEstimate)
    assert estimate.method == "central_difference"
    assert estimate.derivative_estimate == pytest.approx(6.0)


def test_invalid_finite_difference_perturbations_are_rejected() -> None:
    for perturbation in (0.0, -0.1, float("nan"), float("inf")):
        with pytest.raises(ScientificValidationError):
            estimate_central_difference(lambda value: value, 1.0, perturbation)


def test_oscillator_sensitivity_records_base_and_symmetric_runs() -> None:
    specification = make_oscillator_experiment(
        sweep_values=None,
        duration_s=1.2,
        output_sampling=OutputSampling(point_count=81),
        metrics=("final_displacement_m",),
    )

    result = run_sensitivity_analysis(
        specification,
        parameter_name="stiffness_n_per_m",
        metric_id="final_displacement_m",
        perturbation=0.01,
    )

    assert isinstance(result, SensitivityResult)
    assert result.method == "central_difference"
    assert result.lower_parameter_value == pytest.approx(3.99)
    assert result.upper_parameter_value == pytest.approx(4.01)
    assert result.derivative_unit == "m/(N/m)"
    assert tuple(run.run_index for run in result.runs) == (0, 1, 2)
    assert np.isfinite(result.derivative_estimate)


def test_sensitivity_reports_invalid_perturbation_run_without_one_sided_fallback() -> (
    None
):
    specification = make_oscillator_experiment(
        damping_coefficient_kg_per_s=0.0,
        sweep_values=None,
        duration_s=1.0,
        metrics=("final_displacement_m",),
    )

    result = run_sensitivity_analysis(
        specification,
        parameter_name="damping_coefficient_kg_per_s",
        metric_id="final_displacement_m",
        perturbation=0.1,
    )

    assert isinstance(result, SensitivityFailure)
    assert result.failed_run_indices == (1,)
    assert result.runs[1].failure_kind == ExperimentFailureKind.MODEL_VALIDATION
    assert "No one-sided estimate" in result.message


def test_invalid_sensitivity_request_and_metric_values_are_rejected() -> None:
    with pytest.raises(ScientificValidationError):
        run_sensitivity_analysis(
            _spec(),
            parameter_name="damping_coefficient_kg_per_s",
            metric_id="final_displacement_m",
            perturbation=0.0,
        )
    with pytest.raises(ValidationError):
        MetricValue(metric_id="bad", value=float("nan"), unit="m")


def test_sweep_values_must_be_finite_and_reference_metric_requested() -> None:
    with pytest.raises(ValidationError):
        make_oscillator_experiment(sweep_values=(0.0, float("nan")))
    with pytest.raises(ValidationError):
        make_oscillator_experiment(
            metrics=("final_velocity_m_per_s",),
            analytical_reference=AnalyticalReferenceRequest(
                metric_id="final_displacement_m", absolute_tolerance=1e-6
            ),
        )


def test_experiment_rejects_unknown_sweep_parameter_and_duplicate_metrics() -> None:
    specification = make_oscillator_experiment(sweep_values=None)
    payload = specification.model_dump(mode="python")
    payload["sweep"] = {"parameter_name": "unknown_parameter", "values": (1.0,)}
    with pytest.raises(ValidationError, match="sweep parameter"):
        type(specification).model_validate(payload)
    with pytest.raises(ValidationError, match="metric IDs"):
        make_oscillator_experiment(
            sweep_values=None,
            metrics=("final_displacement_m", "final_displacement_m"),
        )


def test_solver_api_compatibility_is_preserved() -> None:
    from newton_lab.algebraic import (
        ScipyHybridRootSolver,
        make_nonlinear_spring_equilibrium_problem,
    )
    from newton_lab.boundary_value import (
        BoundaryValueFailure,
        ScipyBoundaryValueSolver,
        make_linear_heat_conduction_problem,
    )
    from newton_lab.dynamics import simulate_damped_oscillator
    from newton_lab.pendulum import simulate_damped_pendulum

    assert simulate_damped_oscillator(1.0, 0.0, 1.0, 0.1, 0.0, 0.1).time_s.size == 1001
    assert (
        simulate_damped_pendulum(1.0, 1.0, 0.0, 9.81, 0.1, 0.0, 0.1).time_s.size == 1001
    )
    assert ScipyHybridRootSolver().solve(make_nonlinear_spring_equilibrium_problem())
    bvp = ScipyBoundaryValueSolver().solve(make_linear_heat_conduction_problem())
    assert not isinstance(bvp, BoundaryValueFailure)


def test_ode_trajectory_retention_is_opt_in_and_copies_requested_samples() -> None:
    default_record = run_experiment(_spec())
    assert default_record.runs[0].ode_trajectory is None

    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="retained_trajectory",
            sweep_values=None,
            duration_s=8.0,
            initial_displacement_m=0.25,
            initial_velocity_m_per_s=-0.1,
            output_sampling=OutputSampling(point_count=401),
            retain_ode_trajectory=True,
        )
    )
    trajectory = record.runs[0].ode_trajectory
    assert trajectory is not None
    assert trajectory.model_id == "damped_harmonic_oscillator"
    assert trajectory.state_names == ("displacement", "velocity")
    assert trajectory.state_units == ("m", "m/s")
    assert len(trajectory.time_s) == 401
    assert len(trajectory.state_values[0]) == len(trajectory.time_s)
    assert trajectory.state_values[0][0] == pytest.approx(0.25)
    assert trajectory.state_values[1][0] == pytest.approx(-0.1)
