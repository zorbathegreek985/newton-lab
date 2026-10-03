"""Tests for spatial boundary-value experiments and recorded evidence."""

import numpy as np
import pytest

import newton_lab.boundary_value as boundary_value
from newton_lab.boundary_value import (
    BoundaryValueOutcome,
    BoundaryValueProblem,
    BoundaryValueResult,
    BoundaryValueSolverConfiguration,
    ScipyBoundaryValueSolver,
)
from newton_lab.experiments import (
    BVPReferenceRequest,
    ExperimentFailureKind,
    ExperimentRecordStatus,
    ExperimentRunStatus,
    SensitivityFailure,
    SensitivityResult,
    make_heat_conduction_experiment,
    run_experiment,
    run_sensitivity_analysis,
)


def test_successful_bvp_run_preserves_mesh_profile_and_spatial_metrics() -> None:
    specification = make_heat_conduction_experiment(
        length_m=2.0,
        left_temperature_k=290.0,
        right_temperature_k=350.0,
        sweep_values=None,
        metrics=(
            "sampled_minimum:temperature",
            "sampled_minimum_location:temperature",
            "sampled_maximum:temperature",
            "sampled_maximum_location:temperature",
            "left_endpoint:temperature",
            "right_endpoint:temperature",
            "boundary_residual_max_abs",
        ),
    )

    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.SUCCEEDED
    assert run.algebraic_diagnostics is None
    diagnostics = run.bvp_diagnostics
    assert diagnostics is not None
    assert diagnostics.mesh[0] == 0.0 and diagnostics.mesh[-1] == 2.0
    assert len(diagnostics.solution[0]) == len(diagnostics.mesh)
    assert diagnostics.state_names == ("temperature", "temperature_gradient")
    assert diagnostics.solver_success
    assert diagnostics.boundary_conditions_satisfied
    assert diagnostics.numerical_accepted
    assert diagnostics.physical_validity_assessed is False
    metrics = {metric.metric_id: metric.value for metric in run.metrics}
    assert metrics["sampled_minimum:temperature"] == pytest.approx(290.0)
    assert metrics["sampled_minimum_location:temperature"] == pytest.approx(0.0)
    assert metrics["sampled_maximum:temperature"] == pytest.approx(350.0)
    assert metrics["sampled_maximum_location:temperature"] == pytest.approx(2.0)
    assert metrics["left_endpoint:temperature"] == pytest.approx(290.0)
    assert metrics["right_endpoint:temperature"] == pytest.approx(350.0)
    assert metrics["boundary_residual_max_abs"] < 1e-6
    assert (
        next(
            metric
            for metric in run.metrics
            if metric.metric_id == "boundary_residual_max_abs"
        ).unit
        == "K"
    )
    assert run.reproducibility.bvp_domain == (0.0, 2.0)
    assert specification.bvp_problem is not None
    assert run.reproducibility.bvp_mesh == specification.bvp_problem.mesh


def test_heat_reference_is_compared_on_the_actual_returned_mesh() -> None:
    specification = make_heat_conduction_experiment(
        left_temperature_k=300.0,
        right_temperature_k=400.0,
        sweep_values=None,
        bvp_reference=BVPReferenceRequest(maximum_absolute_tolerance=1e-6),
    )

    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.SUCCEEDED
    comparison = run.bvp_reference_comparison
    assert comparison is not None
    assert comparison.status == "passed"
    assert comparison.method == "direct_evaluation_on_returned_mesh"
    assert comparison.maximum_absolute_error is not None
    assert comparison.maximum_absolute_error <= 1e-6
    assert comparison.rms_error is not None


def test_reference_is_unavailable_for_another_bvp_model() -> None:
    specification = make_heat_conduction_experiment(
        sweep_values=None,
        bvp_reference=BVPReferenceRequest(maximum_absolute_tolerance=1e-6),
    )
    assert specification.bvp_problem is not None
    incompatible_problem = specification.bvp_problem.model_copy(
        update={"problem_id": "other_spatial_problem"}
    )
    incompatible_specification = specification.model_copy(
        update={"bvp_problem": incompatible_problem}
    )

    run = run_experiment(incompatible_specification).runs[0]

    assert run.status == ExperimentRunStatus.SUCCEEDED
    assert run.bvp_reference_comparison is not None
    assert run.bvp_reference_comparison.status == "unavailable"


def test_boundary_temperature_sweep_preserves_order_and_validates_domain() -> None:
    specification = make_heat_conduction_experiment(
        sweep_parameter_name="right_temperature_k",
        sweep_values=(350.0, -1.0, 450.0),
        metrics=("right_endpoint:temperature",),
    )
    original_parameters = specification.parameters

    record = run_experiment(specification)

    assert record.status == ExperimentRecordStatus.COMPLETED_WITH_FAILURES
    assert [run.swept_parameter_value for run in record.runs] == [350.0, -1.0, 450.0]
    assert [run.parameters[-1].value for run in record.runs] == [350.0, -1.0, 450.0]
    assert record.runs[0].status == ExperimentRunStatus.SUCCEEDED
    assert record.runs[1].failure_kind == ExperimentFailureKind.MODEL_VALIDATION
    assert "non-negative" in (record.runs[1].message or "")
    assert record.runs[2].status == ExperimentRunStatus.SUCCEEDED
    assert record.runs[0].metrics[0].value == pytest.approx(350.0)
    assert record.runs[2].metrics[0].value == pytest.approx(450.0)
    assert specification.parameters == original_parameters
    assert all(
        run.reproducibility.parameters[-1].value == value
        for run, value in zip(record.runs, (350.0, -1.0, 450.0), strict=True)
    )


def test_solver_success_and_numerical_acceptance_are_separate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class Candidate:
        x = np.asarray([0.0, 0.5, 1.0])
        y = np.asarray([[300.0, 350.0, 400.0], [100.0, 100.0, 100.0]])
        rms_residuals = np.asarray([0.1, 0.1])
        success = True
        status = 0
        message = "solver reported success"
        niter = 1

    monkeypatch.setattr(
        boundary_value, "solve_bvp", lambda *args, **kwargs: Candidate()
    )
    specification = make_heat_conduction_experiment(
        sweep_values=None,
        metrics=("right_endpoint:temperature",),
        solver_configuration=BoundaryValueSolverConfiguration(
            acceptance_tolerance=1e-5
        ),
    )

    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.FAILED
    assert run.failure_kind == ExperimentFailureKind.RESIDUAL_ACCEPTANCE
    assert run.solver_success is True
    assert run.numerical_acceptance is False
    assert run.bvp_diagnostics is not None
    assert run.bvp_diagnostics.solver_success is True
    assert run.bvp_diagnostics.numerical_accepted is False
    assert run.metrics[0].value == pytest.approx(400.0)


def test_solver_failure_is_preserved_and_sweep_continues(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class Candidate:
        x = np.asarray([0.0, 1.0])
        y = np.asarray([[300.0, 400.0], [100.0, 100.0]])
        rms_residuals = np.asarray([0.0])
        success = False
        status = 1
        message = "maximum nodes exceeded"
        niter = 2

    monkeypatch.setattr(
        boundary_value, "solve_bvp", lambda *args, **kwargs: Candidate()
    )
    specification = make_heat_conduction_experiment(
        sweep_values=(350.0, 400.0), metrics=("right_endpoint:temperature",)
    )

    record = run_experiment(specification)

    assert record.status == ExperimentRecordStatus.COMPLETED_WITH_FAILURES
    assert all(run.failure_kind == ExperimentFailureKind.SOLVER for run in record.runs)
    assert all(run.solver_success is False for run in record.runs)
    assert all(run.metrics[0].value == pytest.approx(400.0) for run in record.runs)
    assert all(
        run.bvp_diagnostics is not None
        and run.bvp_diagnostics.message == "maximum nodes exceeded"
        for run in record.runs
    )


def test_bvp_sensitivity_is_central_and_domain_failures_are_retained() -> None:
    specification = make_heat_conduction_experiment(
        left_temperature_k=5.0,
        right_temperature_k=10.0,
        sweep_values=None,
        metrics=("right_endpoint:temperature",),
    )
    result = run_sensitivity_analysis(
        specification,
        parameter_name="left_temperature_k",
        metric_id="right_endpoint:temperature",
        perturbation=1.0,
    )
    assert isinstance(result, SensitivityResult)
    assert result.method == "central_difference"
    assert result.estimate_kind == "mesh_dependent_metric_difference"
    assert (result.lower_parameter_value, result.upper_parameter_value) == (4.0, 6.0)
    assert len(result.runs) == 3

    constrained_specification = make_heat_conduction_experiment(
        left_temperature_k=0.0,
        right_temperature_k=10.0,
        sweep_values=None,
        metrics=("right_endpoint:temperature",),
    )
    constrained = run_sensitivity_analysis(
        constrained_specification,
        parameter_name="left_temperature_k",
        metric_id="right_endpoint:temperature",
        perturbation=1.0,
    )
    assert isinstance(constrained, SensitivityFailure)
    assert constrained.failed_run_indices == (1,)
    assert constrained.runs[1].failure_kind == ExperimentFailureKind.MODEL_VALIDATION


def test_malformed_nonfinite_result_is_recorded_as_metric_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    solve = ScipyBoundaryValueSolver.solve

    def return_nonfinite_candidate(
        solver: ScipyBoundaryValueSolver,
        problem: BoundaryValueProblem,
        configuration: BoundaryValueSolverConfiguration,
    ) -> BoundaryValueOutcome:
        result = solve(solver, problem, configuration)
        if isinstance(result, BoundaryValueResult):
            result.solution[0, 0] = np.nan
        return result

    monkeypatch.setattr(ScipyBoundaryValueSolver, "solve", return_nonfinite_candidate)
    specification = make_heat_conduction_experiment(
        sweep_values=None, metrics=("right_endpoint:temperature",)
    )

    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.FAILED
    assert run.failure_kind == ExperimentFailureKind.METRIC_EXTRACTION
    assert run.bvp_diagnostics is None
