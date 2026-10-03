"""Tests for nonlinear algebraic problem contracts and the SciPy root family."""

from collections.abc import Mapping
from types import SimpleNamespace

import numpy as np
import pytest
from pydantic import ValidationError

import newton_lab.algebraic as algebraic
from newton_lab.algebraic import (
    AlgebraicFailure,
    AlgebraicPhysicalValidity,
    AlgebraicResult,
    AlgebraicSolverConfiguration,
    AlgebraicUnknown,
    NonlinearAlgebraicProblem,
    ResidualFunction,
    ScipyHybridRootSolver,
    make_nonlinear_spring_equilibrium_problem,
)


def test_cubic_spring_finds_known_equilibrium_and_reports_residual() -> None:
    problem = make_nonlinear_spring_equilibrium_problem()
    result = ScipyHybridRootSolver().solve(problem)

    assert isinstance(result, AlgebraicResult)
    assert result.success
    assert result.solver_converged
    assert result.residual_accepted
    assert result.variable_names == ("displacement",)
    assert result.variable_units == ("m",)
    assert result.residual_units == ("N",)
    assert result.solution_candidate[0] == pytest.approx(0.1, abs=1e-10)
    assert result.residual_vector[0] == pytest.approx(0.0, abs=1e-10)
    assert result.residual_norm == pytest.approx(
        np.linalg.norm(result.residual_vector, ord=2)
    )
    assert result.physical_validity == AlgebraicPhysicalValidity.NOT_ASSESSED
    assert not result.uniqueness_assessed
    assert not result.solution_set_completeness_assessed
    assert result.iteration_count is None
    assert result.function_evaluations > 0


def test_square_two_equation_system_with_unknown_units() -> None:
    problem = NonlinearAlgebraicProblem(
        problem_id="unit_circle_diagonal_intersection",
        name="Unit circle and diagonal intersection",
        unknowns=(
            AlgebraicUnknown(name="x", unit=None, initial_guess=0.7),
            AlgebraicUnknown(name="y", unit=None, initial_guess=0.6),
        ),
        residual_function=lambda values, _parameters: np.asarray(
            [values[0] ** 2 + values[1] ** 2 - 1.0, values[0] - values[1]]
        ),
        residual_names=("unit_circle", "diagonal"),
        residual_units=(None, None),
    )
    result = ScipyHybridRootSolver().solve(problem)

    assert isinstance(result, AlgebraicResult)
    assert result.success
    expected = np.sqrt(0.5)
    np.testing.assert_allclose(
        result.solution_candidate, (expected, expected), atol=1e-8
    )
    assert result.variable_units == (None, None)
    assert result.residual_units == (None, None)


@pytest.mark.parametrize(
    "residual",
    [
        lambda _values, _parameters: 1.0,
        lambda _values, _parameters: np.asarray([np.nan]),
    ],
)
def test_invalid_residual_shape_or_nonfinite_values_are_structured_failures(
    residual: ResidualFunction,
) -> None:
    problem = NonlinearAlgebraicProblem(
        problem_id="invalid_residual_example",
        name="Invalid residual test",
        unknowns=(AlgebraicUnknown(name="x", unit="m", initial_guess=1.0),),
        residual_function=residual,
    )

    result = ScipyHybridRootSolver().solve(problem)

    assert isinstance(result, AlgebraicFailure)
    assert result.last_candidate == (1.0,)
    assert result.function_evaluations == 1
    assert "residual function" in result.message


def test_problem_and_solver_reject_invalid_guesses_dimensions_and_tolerances() -> None:
    with pytest.raises(ValidationError, match="finite real number"):
        AlgebraicUnknown(name="x", initial_guess=float("inf"))
    with pytest.raises(ValidationError, match="one residual name per unknown"):
        NonlinearAlgebraicProblem(
            problem_id="mismatched_residual_names",
            name="Mismatched residual metadata",
            unknowns=(AlgebraicUnknown(name="x", initial_guess=0.0),),
            residual_function=lambda values, _parameters: values,
            residual_names=("r1", "r2"),
        )
    with pytest.raises(ValidationError, match="residual units must be identical"):
        NonlinearAlgebraicProblem(
            problem_id="mixed_residual_units",
            name="Mixed residual units",
            unknowns=(
                AlgebraicUnknown(name="x", initial_guess=0.0),
                AlgebraicUnknown(name="y", initial_guess=0.0),
            ),
            residual_function=lambda values, _parameters: values,
            residual_units=("N", "m"),
        )
    with pytest.raises(ValidationError):
        AlgebraicUnknown(name="boolean_guess", initial_guess=True)
    with pytest.raises(ValidationError, match="method='hybr' only"):
        AlgebraicSolverConfiguration(method="lm")
    with pytest.raises(ValidationError):
        AlgebraicSolverConfiguration(residual_tolerance=0.0)
    with pytest.raises(ValidationError):
        AlgebraicSolverConfiguration(maximum_function_evaluations=0)
    with pytest.raises(ValidationError):
        AlgebraicSolverConfiguration(step_tolerance=True)


def test_initial_guess_and_residual_evaluation_do_not_mutate_problem_inputs() -> None:
    observed_guesses: list[float] = []

    def residual(values: np.ndarray, _parameters: Mapping[str, float]) -> np.ndarray:
        observed_guesses.append(float(values[0]))
        values[0] = 123.0
        return np.asarray([values[0] - 0.1])

    problem = NonlinearAlgebraicProblem(
        problem_id="copy_isolation_test",
        name="Copy isolation",
        unknowns=(AlgebraicUnknown(name="x", initial_guess=0.0),),
        residual_function=residual,
    )
    original_guess = problem.initial_guess
    problem.initial_guess[0] = 55.0
    result = ScipyHybridRootSolver().solve(problem)

    assert original_guess[0] == 0.0
    assert problem.initial_guess[0] == 0.0
    assert isinstance(result, AlgebraicResult)
    assert observed_guesses


def test_nonconverged_solver_candidate_is_preserved_and_residual_checked(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    problem = make_nonlinear_spring_equilibrium_problem()

    def nonconverged(*_args: object, **_kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(
            success=False,
            x=np.asarray([0.1]),
            message="maximum evaluations reached",
            nfev=4,
        )

    monkeypatch.setattr(algebraic, "root", nonconverged)
    result = ScipyHybridRootSolver().solve(problem)

    assert isinstance(result, AlgebraicResult)
    assert not result.solver_converged
    assert result.residual_accepted
    assert not result.success
    assert result.solution_candidate[0] == pytest.approx(0.1)
    assert result.termination_reason == "maximum evaluations reached"


def test_iteration_limit_is_reported_without_claiming_a_root() -> None:
    problem = make_nonlinear_spring_equilibrium_problem(initial_displacement_m=10.0)
    result = ScipyHybridRootSolver().solve(
        problem,
        AlgebraicSolverConfiguration(
            maximum_function_evaluations=1,
            residual_tolerance=1e-12,
        ),
    )

    assert isinstance(result, AlgebraicResult)
    assert not result.success
    assert not result.solver_converged
    assert not result.residual_accepted
    assert "maxfev" in result.termination_reason


def test_capabilities_and_reproducibility_are_explicit() -> None:
    solver = ScipyHybridRootSolver()
    problem = make_nonlinear_spring_equilibrium_problem()
    result = solver.solve(problem)

    assert "nonlinear_algebraic_system" in solver.capabilities.problem_families
    assert solver.capabilities.supported_methods == ("hybr",)
    assert not solver.capabilities.supports_bounds_or_constraints
    assert isinstance(result, AlgebraicResult)
    assert result.reproducibility.problem_id == problem.problem_id
    assert result.reproducibility.solver_id == solver.capabilities.solver_id
    assert result.reproducibility.initial_guess == (0.08,)
    assert result.reproducibility.scipy_version
