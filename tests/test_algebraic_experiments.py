"""Tests for local nonlinear algebraic experiments and diagnostics."""

from collections.abc import Mapping
from types import SimpleNamespace

import numpy as np
import pytest
from pydantic import ValidationError

import newton_lab.algebraic as algebraic
from newton_lab.algebraic import (
    AlgebraicParameter,
    AlgebraicSolverConfiguration,
    AlgebraicUnknown,
    NonlinearAlgebraicProblem,
)
from newton_lab.experiments import (
    AnalyticalReferenceRequest,
    ExperimentFailureKind,
    ExperimentRecordStatus,
    ExperimentRunStatus,
    ExperimentSpecification,
    ParameterSweep,
    SensitivityResult,
    make_nonlinear_spring_experiment,
    run_experiment,
    run_sensitivity_analysis,
)
from newton_lab.simulation import NamedValue, ODESolverConfiguration


def test_known_nonlinear_spring_root_and_reference_comparison() -> None:
    # Independently evaluate 100*x + 1000*x**3 - 1.001 at x=0.01.
    reference_displacement_m = 0.01
    force_balance_n = (
        100.0 * reference_displacement_m + 1000.0 * reference_displacement_m**3 - 1.001
    )
    assert force_balance_n == pytest.approx(0.0, abs=1e-14)

    specification = make_nonlinear_spring_experiment(
        sweep_values=None,
        initial_displacement_m=0.005,
        metrics=("solution:displacement", "residual:force_balance", "residual_norm"),
        analytical_reference=AnalyticalReferenceRequest(
            metric_id="solution:displacement", absolute_tolerance=1e-8
        ),
        reference_solution=(reference_displacement_m,),
    )
    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.SUCCEEDED
    metrics = {metric.metric_id: metric for metric in run.metrics}
    assert metrics["solution:displacement"].value == pytest.approx(
        reference_displacement_m, abs=1e-8
    )
    assert metrics["solution:displacement"].unit == "m"
    assert metrics["residual:force_balance"].value == pytest.approx(0.0, abs=1e-8)
    assert metrics["residual:force_balance"].unit == "N"
    assert metrics["residual_norm"].unit == "N"
    assert run.reference_comparison is not None
    assert run.reference_comparison.status == "passed"
    assert run.reference_comparison.numerical_value == pytest.approx(0.01, abs=1e-8)
    assert run.reference_comparison.reference_value == pytest.approx(0.01)
    assert run.reference_comparison.absolute_error is not None
    assert run.reference_comparison.absolute_error <= 1e-8
    assert run.reference_comparison.absolute_tolerance == 1e-8
    diagnostics = run.algebraic_diagnostics
    assert diagnostics is not None
    assert diagnostics.solver_converged is True
    assert diagnostics.residual_accepted is True
    assert diagnostics.physical_validity == "not_assessed"
    assert diagnostics.uniqueness_assessed is False
    assert diagnostics.solution_set_completeness_assessed is False


def test_force_sweep_preserves_order_and_fixed_initial_guess() -> None:
    specification = make_nonlinear_spring_experiment(
        applied_force_n=1.1,
        initial_displacement_m=0.005,
        sweep_values=(2.008, 0.8, 1.1),
        metrics=("solution:displacement", "residual_norm"),
    )
    assert specification.algebraic_problem is not None
    initial_guess = specification.algebraic_problem.initial_guess.copy()

    record = run_experiment(specification)

    assert record.status == ExperimentRecordStatus.SUCCEEDED
    assert [run.swept_parameter_value for run in record.runs] == [2.008, 0.8, 1.1]
    assert [run.parameters[-1].value for run in record.runs] == [2.008, 0.8, 1.1]
    assert all(
        run.algebraic_diagnostics is not None
        and run.algebraic_diagnostics.initial_guess == (0.005,)
        for run in record.runs
    )
    assert specification.algebraic_problem is not None
    assert np.array_equal(specification.algebraic_problem.initial_guess, initial_guess)
    assert all(run.reproducibility.output_point_count is None for run in record.runs)
    assert all(
        run.reproducibility.maximum_function_evaluations == 1000 for run in record.runs
    )


def test_solver_nonconvergence_preserves_candidate_and_diagnostics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    specification = make_nonlinear_spring_experiment(
        sweep_values=None, metrics=("solution:displacement", "residual_norm")
    )

    def nonconverged(*_args: object, **_kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(
            success=False,
            x=np.asarray([0.02]),
            message="maximum evaluations reached",
            nfev=4,
        )

    monkeypatch.setattr(algebraic, "root", nonconverged)
    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.FAILED
    assert run.failure_kind == ExperimentFailureKind.SOLVER
    assert run.solver_success is False
    assert run.algebraic_diagnostics is not None
    assert run.algebraic_diagnostics.solution_candidate == pytest.approx((0.02,))
    assert run.algebraic_diagnostics.termination_reason == "maximum evaluations reached"
    assert {metric.metric_id for metric in run.metrics} == {
        "solution:displacement",
        "residual_norm",
    }


def test_solver_convergence_and_residual_acceptance_remain_distinct(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    specification = make_nonlinear_spring_experiment(
        sweep_values=None,
        metrics=("solution:displacement",),
        solver_configuration=AlgebraicSolverConfiguration(residual_tolerance=1e-12),
    )

    def converged_but_inaccurate(*_args: object, **_kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(
            success=True,
            x=np.asarray([0.02]),
            message="converged by step tolerance",
            nfev=4,
        )

    monkeypatch.setattr(algebraic, "root", converged_but_inaccurate)
    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.FAILED
    assert run.failure_kind == ExperimentFailureKind.RESIDUAL_ACCEPTANCE
    assert run.solver_success is True
    assert run.numerical_acceptance is False
    assert run.algebraic_diagnostics is not None
    assert run.algebraic_diagnostics.solver_converged is True
    assert run.algebraic_diagnostics.residual_accepted is False


def test_sweep_continues_after_expected_solver_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fails_every_run(*_args: object, **_kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(
            success=False,
            x=np.asarray([0.005]),
            message="controlled nonconvergence",
            nfev=2,
        )

    monkeypatch.setattr(algebraic, "root", fails_every_run)
    specification = make_nonlinear_spring_experiment(
        sweep_values=(0.8, 1.1, 2.008), metrics=("solution:displacement",)
    )
    record = run_experiment(specification)

    assert record.status == ExperimentRecordStatus.COMPLETED_WITH_FAILURES
    assert [run.swept_parameter_value for run in record.runs] == [0.8, 1.1, 2.008]
    assert all(run.failure_kind == ExperimentFailureKind.SOLVER for run in record.runs)


def test_initial_guess_selects_candidates_without_claiming_uniqueness() -> None:
    def square_roots(
        unknowns: np.ndarray, parameters: Mapping[str, float]
    ) -> np.ndarray:
        return np.asarray([unknowns[0] ** 2 - parameters["target"]])

    def make_problem(initial_guess: float) -> ExperimentSpecification:
        problem = NonlinearAlgebraicProblem(
            problem_id="two_root_polynomial",
            name="Two-root polynomial",
            unknowns=(
                AlgebraicUnknown(name="x", unit="1", initial_guess=initial_guess),
            ),
            residual_function=square_roots,
            residual_names=("polynomial",),
            residual_units=("1",),
            parameters=(AlgebraicParameter(name="target", value=1.0, unit="1"),),
        )
        return ExperimentSpecification(
            experiment_id=(
                f"root_from_{'positive' if initial_guess > 0 else 'negative'}"
            ),
            name="Initial-guess dependence demonstration",
            model_id="nonlinear_algebraic_system",
            parameters=(NamedValue(name="target", value=1.0, unit="1"),),
            metrics=("solution:x",),
            solver_configuration=AlgebraicSolverConfiguration(),
            algebraic_problem=problem,
        )

    positive = run_experiment(make_problem(0.5)).runs[0]
    negative = run_experiment(make_problem(-0.5)).runs[0]

    assert positive.metrics[0].value == pytest.approx(1.0)
    assert negative.metrics[0].value == pytest.approx(-1.0)
    assert positive.algebraic_diagnostics is not None
    assert negative.algebraic_diagnostics is not None
    assert positive.algebraic_diagnostics.initial_guess == (0.5,)
    assert negative.algebraic_diagnostics.initial_guess == (-0.5,)
    assert positive.algebraic_diagnostics.uniqueness_assessed is False
    assert positive.algebraic_diagnostics.solution_set_completeness_assessed is False


def test_algebraic_sensitivity_is_labeled_candidate_difference() -> None:
    specification = make_nonlinear_spring_experiment(
        applied_force_n=1.1,
        initial_displacement_m=0.005,
        sweep_values=None,
        metrics=("solution:displacement",),
    )

    result = run_sensitivity_analysis(
        specification,
        parameter_name="applied_force_n",
        metric_id="solution:displacement",
        perturbation=0.01,
    )

    assert isinstance(result, SensitivityResult)
    assert result.estimate_kind == "candidate_difference"
    assert result.candidate_branch_consistency == "not_established"
    assert result.derivative_unit == "m/N"
    assert len(result.runs) == 3
    assert all(
        run.algebraic_diagnostics is not None
        and run.algebraic_diagnostics.initial_guess == (0.005,)
        for run in result.runs
    )


def test_algebraic_adapter_rejects_ode_configuration() -> None:
    specification = make_nonlinear_spring_experiment(sweep_values=None)
    with pytest.raises(ValidationError):
        ExperimentSpecification.model_validate(
            {
                **specification.model_dump(mode="python"),
                "solver_configuration": ODESolverConfiguration(),
            }
        )


def test_spring_parameter_domain_is_checked_for_each_sweep_point() -> None:
    specification = make_nonlinear_spring_experiment(
        sweep_values=None, metrics=("solution:displacement",)
    ).model_copy(
        update={
            "sweep": ParameterSweep(
                parameter_name="linear_stiffness_n_per_m", values=(-1.0, 100.0)
            )
        }
    )

    record = run_experiment(specification)

    assert record.status == ExperimentRecordStatus.COMPLETED_WITH_FAILURES
    assert record.runs[0].failure_kind == ExperimentFailureKind.MODEL_VALIDATION
    assert "must be positive" in (record.runs[0].message or "")
    assert record.runs[1].status == ExperimentRunStatus.SUCCEEDED


def test_reference_comparison_tolerance_is_separate_from_residual_acceptance(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    specification = make_nonlinear_spring_experiment(
        sweep_values=None,
        metrics=("solution:displacement",),
        solver_configuration=AlgebraicSolverConfiguration(residual_tolerance=0.01),
        analytical_reference=AnalyticalReferenceRequest(
            metric_id="solution:displacement", absolute_tolerance=1e-8
        ),
        reference_solution=(0.01,),
    )

    def converged_with_accepted_residual(
        *_args: object, **_kwargs: object
    ) -> SimpleNamespace:
        return SimpleNamespace(
            success=True,
            x=np.asarray([0.01001]),
            message="converged",
            nfev=4,
        )

    monkeypatch.setattr(algebraic, "root", converged_with_accepted_residual)
    run = run_experiment(specification).runs[0]

    assert run.status == ExperimentRunStatus.SUCCEEDED
    assert run.algebraic_diagnostics is not None
    assert run.algebraic_diagnostics.residual_accepted is True
    assert run.reference_comparison is not None
    assert run.reference_comparison.status == "failed"
