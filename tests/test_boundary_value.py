"""Tests for the distinct boundary-value solver family."""

import numpy as np
import pytest

from newton_lab.boundary_value import (
    BoundaryValueFailure,
    BoundaryValueProblem,
    BoundaryValueSolverConfiguration,
    ScipyBoundaryValueSolver,
    make_linear_heat_conduction_problem,
)


def test_linear_heat_solution_matches_analytical_profile() -> None:
    problem = make_linear_heat_conduction_problem(
        length_m=2.0, left_temperature_k=290.0, right_temperature_k=350.0
    )
    result = ScipyBoundaryValueSolver().solve(problem)

    assert not isinstance(result, BoundaryValueFailure)
    expected = 290.0 + 30.0 * result.mesh
    assert result.solution.shape == (2, len(result.mesh))
    assert np.all(np.isfinite(result.solution))
    assert np.allclose(result.solution[0], expected, atol=1e-7, rtol=1e-9)
    assert result.solution[0, 0] == pytest.approx(290.0)
    assert result.solution[0, -1] == pytest.approx(350.0)
    assert result.solver_success and result.numerical_accepted
    assert result.boundary_conditions_satisfied
    assert not result.physical_validity_assessed
    assert not result.uniqueness_assessed
    assert not result.solution_set_completeness_assessed


def test_inputs_are_not_mutated_and_repeatable() -> None:
    problem = make_linear_heat_conduction_problem()
    mesh_before = problem.mesh
    guess_before = problem.initial_guess
    solver = ScipyBoundaryValueSolver()
    first, second = solver.solve(problem), solver.solve(problem)
    assert problem.mesh == mesh_before
    assert problem.initial_guess == guess_before
    assert not isinstance(first, BoundaryValueFailure)
    assert not isinstance(second, BoundaryValueFailure)
    assert np.array_equal(first.mesh, second.mesh)
    assert np.array_equal(first.solution, second.solution)


@pytest.mark.parametrize(
    "changes",
    [
        {"domain": (1.0, 0.0)},
        {"domain": (0.0, 2.0)},
        {"mesh": (0.0, 1.0, 0.5)},
        {"initial_guess": ((float("nan"),) * 5, (0.0,) * 5)},
        {"initial_guess": ((0.0,) * 4, (0.0,) * 4)},
    ],
)
def test_invalid_problem_grid_is_rejected(changes: dict[str, object]) -> None:
    problem = make_linear_heat_conduction_problem()
    values = problem.model_dump()
    values.update(changes)
    with pytest.raises((ValueError, TypeError)):
        BoundaryValueProblem.model_validate(values)


def test_invalid_solver_settings_are_rejected() -> None:
    with pytest.raises(ValueError):
        BoundaryValueSolverConfiguration(tolerance=0.0)
    with pytest.raises(ValueError):
        BoundaryValueSolverConfiguration(maximum_nodes=1)


def test_invalid_callback_output_returns_clear_failure() -> None:
    problem = make_linear_heat_conduction_problem()
    invalid = problem.model_copy(
        update={"ode_function": lambda x, y, p: np.full((2, len(x)), np.nan)}
    )
    result = ScipyBoundaryValueSolver().solve(invalid)
    assert isinstance(result, BoundaryValueFailure)
    assert result.error_type == "ScientificValidationError"
    assert "finite" in result.message


def test_solver_nonconvergence_is_reported(monkeypatch: pytest.MonkeyPatch) -> None:
    import newton_lab.boundary_value as module

    class Candidate:
        x = np.asarray([0.0, 1.0])
        y = np.asarray([[300.0, 400.0], [100.0, 100.0]])
        rms_residuals = np.asarray([0.0])
        success = False
        status = 1
        message = "maximum nodes exceeded"
        niter = 2

    monkeypatch.setattr(module, "solve_bvp", lambda *args, **kwargs: Candidate())
    result = ScipyBoundaryValueSolver().solve(make_linear_heat_conduction_problem())
    assert not isinstance(result, BoundaryValueFailure)
    assert not result.solver_success
    assert not result.numerical_accepted
    assert result.message == "maximum nodes exceeded"
