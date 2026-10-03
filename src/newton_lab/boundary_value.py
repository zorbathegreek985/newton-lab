"""Boundary-value problems solved by SciPy's collocation solver.

This family represents spatial boundary-value candidates, separately from
initial-value time-series simulations and finite-dimensional root problems.
Callbacks are trusted Python functions; descriptive equation text is never run.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from types import MappingProxyType
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from pydantic import Field, field_validator, model_validator
from scipy.integrate import solve_bvp  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.simulation import NamedValue, SimulationModel, StateVariable

FloatArray = NDArray[np.float64]
ODEFunction = Callable[[FloatArray, FloatArray, Mapping[str, float]], ArrayLike]
BoundaryFunction = Callable[[FloatArray, FloatArray, Mapping[str, float]], ArrayLike]


class BoundaryValueProblem(SimulationModel):
    """First-order system ``y'=f(x,y)`` with two-end boundary residuals.

    The current family accepts exactly two state variables and two boundary
    residuals, the standard first-order form of a scalar second-order problem.
    Callbacks receive arrays and an immutable named-parameter mapping.
    """

    problem_id: str = Field(pattern=r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
    name: str = Field(min_length=1)
    description: str = ""
    domain: tuple[float, float]
    domain_unit: str = Field(default="unknown", min_length=1)
    states: tuple[StateVariable, StateVariable]
    boundary_residual_units: tuple[str | None, str | None] = (None, None)
    parameters: tuple[NamedValue, ...] = ()
    mesh: tuple[float, ...]
    initial_guess: tuple[tuple[float, ...], tuple[float, ...]]
    ode_function: ODEFunction
    boundary_function: BoundaryFunction
    boundary_conditions: tuple[str, ...] = ()
    assumptions: tuple[str, ...] = ()
    validity_limits: tuple[str, ...] = ()
    equation_record_ids: tuple[str, ...] = ()
    source_reference_ids: tuple[str, ...] = ()

    @field_validator("domain")
    @classmethod
    def valid_domain(cls, value: tuple[float, float]) -> tuple[float, float]:
        if len(value) != 2 or not np.isfinite(value).all() or value[1] <= value[0]:
            raise ValueError("domain must have finite increasing endpoints")
        return value

    @field_validator("mesh")
    @classmethod
    def valid_mesh(cls, value: tuple[float, ...]) -> tuple[float, ...]:
        if (
            len(value) < 2
            or not np.isfinite(value).all()
            or not np.all(np.diff(value) > 0)
        ):
            raise ValueError("mesh must contain at least two finite increasing points")
        return value

    @model_validator(mode="after")
    def validate_grid_and_metadata(self) -> BoundaryValueProblem:
        if self.mesh[0] != self.domain[0] or self.mesh[-1] != self.domain[1]:
            raise ValueError("mesh endpoints must match the problem domain")
        guess = np.asarray(self.initial_guess, dtype=np.float64)
        if guess.shape != (2, len(self.mesh)) or not np.isfinite(guess).all():
            raise ValueError("initial_guess must be finite with shape (2, len(mesh))")
        if len({p.name for p in self.parameters}) != len(self.parameters):
            raise ValueError("parameter names must be unique")
        if len({state.name for state in self.states}) != 2:
            raise ValueError("state names must be unique")
        return self


class BoundaryValueSolverConfiguration(SimulationModel):
    """Explicit collocation and numerical acceptance settings."""

    method: Literal["fourth_order_collocation"] = "fourth_order_collocation"
    tolerance: float = Field(default=1e-6, gt=0.0)
    boundary_tolerance: float = Field(default=1e-6, gt=0.0)
    acceptance_tolerance: float = Field(default=1e-5, gt=0.0)
    maximum_nodes: int = Field(default=10000, ge=2, strict=True)

    @field_validator("tolerance", "boundary_tolerance", "acceptance_tolerance")
    @classmethod
    def finite_tolerance(cls, value: float) -> float:
        if not np.isfinite(value):
            raise ValueError("tolerances must be finite")
        return value


class BoundaryValueSolverCapabilities(SimulationModel):
    """Declared support and limitations of the SciPy BVP implementation."""

    solver_id: str = "scipy_solve_bvp"
    solver_version: str
    problem_family: str = "boundary_value_problem"
    method: str = "fourth_order_collocation"
    limitations: tuple[str, ...] = (
        "Supports two-state first-order systems with endpoint boundary residuals.",
        "Does not establish physical validity, uniqueness, or completeness.",
    )


@dataclass(frozen=True, slots=True)
class BoundaryValueResult:
    """Spatial solution candidate and separate numerical acceptance evidence."""

    problem: BoundaryValueProblem
    configuration: BoundaryValueSolverConfiguration
    solver_capabilities: BoundaryValueSolverCapabilities
    mesh: FloatArray
    solution: FloatArray
    differential_residuals: FloatArray
    boundary_residuals: FloatArray
    solver_success: bool
    numerical_accepted: bool
    boundary_conditions_satisfied: bool
    physical_validity_assessed: bool
    uniqueness_assessed: bool
    solution_set_completeness_assessed: bool
    status_code: int
    message: str
    iteration_count: int

    def __post_init__(self) -> None:
        if (
            self.mesh.ndim != 1
            or len(self.mesh) < 2
            or not np.isfinite(self.mesh).all()
        ):
            raise ScientificValidationError(
                "result mesh must be finite and one-dimensional"
            )
        if (
            self.solution.shape != (2, len(self.mesh))
            or not np.isfinite(self.solution).all()
        ):
            raise ScientificValidationError(
                "solution must be finite with shape (2, len(mesh))"
            )
        if (
            self.differential_residuals.ndim != 1
            or not np.isfinite(self.differential_residuals).all()
        ):
            raise ScientificValidationError("differential residuals must be finite")
        if (
            self.boundary_residuals.shape != (2,)
            or not np.isfinite(self.boundary_residuals).all()
        ):
            raise ScientificValidationError(
                "boundary residuals must contain two finite values"
            )
        if not np.all(np.diff(self.mesh) > 0):
            raise ScientificValidationError("result mesh must be strictly increasing")
        if self.numerical_accepted != (
            self.solver_success
            and self.boundary_conditions_satisfied
            and len(self.differential_residuals) > 0
            and float(np.max(self.differential_residuals))
            <= self.configuration.acceptance_tolerance
        ):
            raise ScientificValidationError(
                "numerical acceptance flags are inconsistent"
            )


@dataclass(frozen=True, slots=True)
class BoundaryValueFailure:
    """Failure before SciPy returned a valid numerical candidate."""

    problem_id: str
    solver_id: str
    error_type: str
    message: str


BoundaryValueOutcome = BoundaryValueResult | BoundaryValueFailure
_DEFAULT_CONFIGURATION = BoundaryValueSolverConfiguration()


def _scipy_version() -> str:
    try:
        return version("scipy")
    except PackageNotFoundError:
        return "unknown"


class ScipyBoundaryValueSolver:
    """Solve two-state endpoint BVPs with ``scipy.integrate.solve_bvp``."""

    capabilities = BoundaryValueSolverCapabilities(solver_version=_scipy_version())

    def solve(
        self,
        problem: BoundaryValueProblem,
        configuration: BoundaryValueSolverConfiguration = _DEFAULT_CONFIGURATION,
    ) -> BoundaryValueOutcome:
        """Return solver outcome, retaining nonconverged candidates when possible."""
        parameters = MappingProxyType(
            {item.name: item.value for item in problem.parameters}
        )
        x = np.asarray(problem.mesh, dtype=np.float64).copy()
        y = np.asarray(problem.initial_guess, dtype=np.float64).copy()

        def fun(mesh: FloatArray, values: FloatArray) -> FloatArray:
            candidate = np.asarray(
                problem.ode_function(mesh.copy(), values.copy(), parameters),
                dtype=np.float64,
            )
            if candidate.shape != values.shape or not np.isfinite(candidate).all():
                raise ScientificValidationError(
                    "ODE callback must return a finite array matching the state shape"
                )
            return candidate

        def bc(left: FloatArray, right: FloatArray) -> FloatArray:
            residual = np.asarray(
                problem.boundary_function(left.copy(), right.copy(), parameters),
                dtype=np.float64,
            )
            if residual.shape != (2,) or not np.isfinite(residual).all():
                raise ScientificValidationError(
                    "boundary callback must return two finite residuals"
                )
            return residual

        try:
            # Validate callback contracts before entering SciPy's iteration.
            fun(x, y)
            bc(y[:, 0], y[:, -1])
            solution = solve_bvp(
                fun,
                bc,
                x,
                y,
                tol=configuration.tolerance,
                bc_tol=configuration.boundary_tolerance,
                max_nodes=configuration.maximum_nodes,
                verbose=0,
            )
            result_mesh = np.asarray(solution.x, dtype=np.float64).copy()
            result_y = np.asarray(solution.y, dtype=np.float64).copy()
            residuals = np.asarray(solution.rms_residuals, dtype=np.float64).copy()
            boundary_residuals = bc(result_y[:, 0], result_y[:, -1]).copy()
            if (
                result_mesh.ndim != 1
                or len(result_mesh) < 2
                or result_y.shape != (2, len(result_mesh))
                or not np.isfinite(result_mesh).all()
                or not np.isfinite(result_y).all()
                or residuals.ndim != 1
                or not np.isfinite(residuals).all()
            ):
                raise ScientificValidationError(
                    "SciPy returned malformed or non-finite solution data"
                )
        except (
            ArithmeticError,
            RuntimeError,
            ValueError,
            TypeError,
            ScientificValidationError,
        ) as exc:
            return BoundaryValueFailure(
                problem.problem_id,
                self.capabilities.solver_id,
                type(exc).__name__,
                str(exc),
            )

        boundary_ok = bool(
            np.max(np.abs(boundary_residuals)) <= configuration.acceptance_tolerance
        )
        accepted = bool(
            solution.success
            and boundary_ok
            and len(residuals) > 0
            and np.max(residuals) <= configuration.acceptance_tolerance
        )
        return BoundaryValueResult(
            problem=problem,
            configuration=configuration,
            solver_capabilities=self.capabilities,
            mesh=result_mesh,
            solution=result_y,
            differential_residuals=residuals,
            boundary_residuals=boundary_residuals,
            solver_success=bool(solution.success),
            numerical_accepted=accepted,
            boundary_conditions_satisfied=boundary_ok,
            physical_validity_assessed=False,
            uniqueness_assessed=False,
            solution_set_completeness_assessed=False,
            status_code=int(solution.status),
            message=str(solution.message),
            iteration_count=int(solution.niter),
        )


def make_linear_heat_conduction_problem(
    *,
    length_m: float = 1.0,
    left_temperature_k: float = 300.0,
    right_temperature_k: float = 400.0,
    mesh_points: int = 5,
) -> BoundaryValueProblem:
    """Build steady 1D constant-conductivity heat conduction with fixed ends."""
    if not np.isfinite(length_m) or length_m <= 0:
        raise ScientificValidationError("length_m must be finite and positive")
    if not np.isfinite(left_temperature_k) or not np.isfinite(right_temperature_k):
        raise ScientificValidationError("boundary temperatures must be finite")
    if (
        isinstance(mesh_points, bool)
        or not isinstance(mesh_points, int)
        or mesh_points < 2
    ):
        raise ScientificValidationError(
            "mesh_points must be an integer of at least two"
        )
    mesh = np.linspace(0.0, length_m, mesh_points)

    def heat_equation(
        x: FloatArray, y: FloatArray, _p: Mapping[str, float]
    ) -> FloatArray:
        return np.vstack((y[1], np.zeros_like(x)))

    def fixed_temperatures(
        left: FloatArray, right: FloatArray, p: Mapping[str, float]
    ) -> FloatArray:
        return np.asarray(
            [left[0] - p["left_temperature_k"], right[0] - p["right_temperature_k"]]
        )

    return BoundaryValueProblem(
        problem_id="steady_linear_heat_conduction",
        name="Steady one-dimensional heat conduction",
        description=(
            "Steady conduction in a homogeneous rod with fixed endpoint temperatures."
        ),
        domain=(0.0, length_m),
        domain_unit="m",
        states=(
            StateVariable(name="temperature", unit="K"),
            StateVariable(name="temperature_gradient", unit="K/m"),
        ),
        boundary_residual_units=("K", "K"),
        parameters=(
            NamedValue(name="left_temperature_k", value=left_temperature_k, unit="K"),
            NamedValue(name="right_temperature_k", value=right_temperature_k, unit="K"),
        ),
        mesh=tuple(float(v) for v in mesh),
        initial_guess=(
            (left_temperature_k,) * mesh_points,
            ((right_temperature_k - left_temperature_k) / length_m,) * mesh_points,
        ),
        ode_function=heat_equation,
        boundary_function=fixed_temperatures,
        boundary_conditions=(
            f"T(0) = {left_temperature_k} K",
            f"T(L) = {right_temperature_k} K",
        ),
        assumptions=(
            "Steady state.",
            "One-dimensional homogeneous material with constant thermal conductivity.",
            "No internal heat generation.",
        ),
        validity_limits=(
            "A continuum model; contact resistance, material variation, and "
            "multidimensional effects are excluded.",
        ),
    )
