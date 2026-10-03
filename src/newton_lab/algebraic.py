"""Nonlinear real algebraic systems and a local SciPy root solver family.

Problems provide trusted Python residual functions. Registry expression text is
descriptive metadata and is never parsed or executed here.
"""

from __future__ import annotations

import platform
import sys
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from enum import StrEnum
from importlib.metadata import PackageNotFoundError, version
from types import MappingProxyType
from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray
from pydantic import Field, field_validator, model_validator
from scipy.optimize import root  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.simulation import SimulationModel

ResidualFunction = Callable[[NDArray[np.float64], Mapping[str, float]], ArrayLike]


class AlgebraicUnknown(SimulationModel):
    """One real-valued unknown, with optional unit and finite initial guess."""

    name: str = Field(min_length=1)
    unit: str | None = None
    initial_guess: float = Field(strict=True)

    @field_validator("initial_guess")
    @classmethod
    def guess_must_be_finite(cls, value: float) -> float:
        if isinstance(value, bool) or not np.isfinite(value):
            raise ValueError("initial_guess must be a finite real number")
        return value


class AlgebraicParameter(SimulationModel):
    """Scalar parameter with optional unit metadata."""

    name: str = Field(min_length=1)
    value: float = Field(strict=True)
    unit: str | None = None

    @field_validator("value")
    @classmethod
    def value_must_be_finite(cls, value: float) -> float:
        if not np.isfinite(value):
            raise ValueError("parameter value must be finite")
        return value


class NonlinearAlgebraicProblem(SimulationModel):
    """A square real system ``F(x)=0`` with an explicit trusted residual.

    Residual functions receive a copy of the unknown vector and a read-only
    parameter mapping. Each must return exactly one finite residual per unknown.
    The solver does not enforce descriptive physical constraints.
    """

    problem_id: str = Field(pattern=r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
    name: str = Field(min_length=1)
    description: str = ""
    unknowns: tuple[AlgebraicUnknown, ...]
    residual_function: ResidualFunction
    residual_names: tuple[str, ...] | None = None
    residual_units: tuple[str | None, ...] | None = None
    parameters: tuple[AlgebraicParameter, ...] = ()
    assumptions: tuple[str, ...] | None = None
    constraints: tuple[str, ...] | None = None
    validity_limits: tuple[str, ...] | None = None
    equation_record_ids: tuple[str, ...] = ()
    source_reference_ids: tuple[str, ...] = ()
    model_version: str = "1"

    @model_validator(mode="after")
    def system_dimensions_are_valid(self) -> NonlinearAlgebraicProblem:
        dimension = len(self.unknowns)
        if dimension == 0:
            raise ValueError("a nonlinear algebraic problem needs at least one unknown")
        if len({item.name for item in self.unknowns}) != dimension:
            raise ValueError("unknown variable names must be unique")
        if len({item.name for item in self.parameters}) != len(self.parameters):
            raise ValueError("parameter names must be unique")
        if self.residual_names is not None:
            if len(self.residual_names) != dimension:
                raise ValueError("a square system needs one residual name per unknown")
            if len(set(self.residual_names)) != dimension:
                raise ValueError("residual names must be unique")
        if self.residual_units is not None and len(self.residual_units) != dimension:
            raise ValueError("a square system needs one residual unit per unknown")
        if self.residual_units is not None:
            declared_units = {unit for unit in self.residual_units if unit is not None}
            if len(declared_units) > 1 or (
                declared_units and any(unit is None for unit in self.residual_units)
            ):
                raise ValueError(
                    "residual units must be identical or all unknown; "
                    "implicit scaling across units is unsupported"
                )
        if len(set(self.equation_record_ids)) != len(self.equation_record_ids):
            raise ValueError("equation record IDs must be unique")
        if len(set(self.source_reference_ids)) != len(self.source_reference_ids):
            raise ValueError("source reference IDs must be unique")
        return self

    @property
    def initial_guess(self) -> NDArray[np.float64]:
        """Return a new float64 array so callers cannot mutate problem data."""
        return np.asarray(
            [unknown.initial_guess for unknown in self.unknowns], dtype=np.float64
        )

    @property
    def parameter_values(self) -> Mapping[str, float]:
        """Return immutable scalar model parameters keyed by stable names."""
        return MappingProxyType({item.name: item.value for item in self.parameters})

    def evaluate_residual(self, values: ArrayLike) -> NDArray[np.float64]:
        """Evaluate and validate the user-supplied residual vector."""
        try:
            raw_vector = np.asarray(values)
            if raw_vector.dtype.kind == "b":
                raise ValueError("boolean unknowns are not real-valued inputs")
            vector = np.asarray(values, dtype=np.float64)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ScientificValidationError(
                "unknown vector must contain finite real values"
            ) from exc
        if vector.shape != (len(self.unknowns),):
            raise ScientificValidationError(
                f"unknown vector must have shape ({len(self.unknowns)},)"
            )
        if not np.isfinite(vector).all():
            raise ScientificValidationError(
                "unknown vector must contain only finite values"
            )
        try:
            raw_residual = self.residual_function(vector.copy(), self.parameter_values)
            if np.asarray(raw_residual).dtype.kind == "b":
                raise ValueError("boolean residuals are not real-valued equations")
            residual = np.asarray(raw_residual, dtype=np.float64)
        except Exception as exc:
            raise ScientificValidationError(
                f"residual function raised {type(exc).__name__}: {exc}"
            ) from exc
        if residual.shape != (len(self.unknowns),):
            raise ScientificValidationError(
                "residual function must return a one-dimensional vector with "
                f"shape ({len(self.unknowns)},); got {residual.shape}"
            )
        if not np.isfinite(residual).all():
            raise ScientificValidationError(
                "residual function must return only finite values"
            )
        return residual.copy()


class AlgebraicSolverConfiguration(SimulationModel):
    """Settings for the supported local Powell hybrid root method."""

    method: str = "hybr"
    step_tolerance: float = Field(default=1e-10, gt=0.0, strict=True)
    residual_tolerance: float = Field(default=1e-8, gt=0.0, strict=True)
    maximum_function_evaluations: int = Field(default=1000, ge=1, strict=True)

    @field_validator("step_tolerance", "residual_tolerance")
    @classmethod
    def tolerances_must_be_finite(cls, value: float) -> float:
        if not np.isfinite(value):
            raise ValueError("solver tolerances must be finite")
        return value

    @model_validator(mode="after")
    def supported_method(self) -> AlgebraicSolverConfiguration:
        if self.method != "hybr":
            raise ValueError("this solver family currently supports method='hybr' only")
        return self


class AlgebraicSolverCapabilities(SimulationModel):
    """Scope and known limitations for the SciPy hybrid root solver."""

    solver_id: str = "scipy_optimize_root_hybr"
    solver_version: str
    problem_families: tuple[str, ...] = ("nonlinear_algebraic_system",)
    supported_methods: tuple[str, ...] = ("hybr",)
    supports_square_systems: bool = True
    supports_bounds_or_constraints: bool = False
    supports_sparse_jacobians: bool = False
    limitations: tuple[str, ...] = (
        "Local method; results depend on the initial guess.",
        "Requires a square system of finite-dimensional real equations.",
        "Does not enforce physical constraints or establish physical validity.",
        "Does not establish uniqueness or enumerate all roots.",
        "No automatic variable or residual scaling is applied.",
    )


class AlgebraicPhysicalValidity(StrEnum):
    NOT_ASSESSED = "not_assessed"
    ASSESSED_VALID = "assessed_valid"
    ASSESSED_INVALID = "assessed_invalid"


@dataclass(frozen=True, slots=True)
class AlgebraicReproducibilityRecord:
    """Environment and inputs relevant to rerunning a nonlinear root solve."""

    python_version: str
    platform: str
    numpy_version: str
    scipy_version: str
    problem_id: str
    model_version: str
    initial_guess: tuple[float, ...]
    variable_names: tuple[str, ...]
    variable_units: tuple[str | None, ...]
    parameters: tuple[AlgebraicParameter, ...]
    solver_id: str
    solver_version: str
    method: str
    step_tolerance: float
    residual_tolerance: float
    maximum_function_evaluations: int
    equation_record_ids: tuple[str, ...]
    source_reference_ids: tuple[str, ...]
    assumptions: tuple[str, ...] | None
    constraints: tuple[str, ...] | None
    validity_limits: tuple[str, ...] | None


@dataclass(frozen=True, slots=True)
class AlgebraicResult:
    """Candidate root and distinct solver, residual, and physical assessments."""

    problem_id: str
    solver_id: str
    method: str
    solution_candidate: NDArray[np.float64]
    variable_names: tuple[str, ...]
    variable_units: tuple[str | None, ...]
    residual_vector: NDArray[np.float64]
    residual_names: tuple[str, ...] | None
    residual_units: tuple[str | None, ...] | None
    residual_norm: float
    residual_norm_definition: str
    solver_converged: bool
    residual_accepted: bool
    physical_validity: AlgebraicPhysicalValidity
    uniqueness_assessed: bool
    solution_set_completeness_assessed: bool
    termination_reason: str
    iteration_count: int | None
    function_evaluations: int
    configuration: AlgebraicSolverConfiguration
    parameters: tuple[AlgebraicParameter, ...]
    assumptions: tuple[str, ...] | None
    constraints: tuple[str, ...] | None
    validity_limits: tuple[str, ...] | None
    reproducibility: AlgebraicReproducibilityRecord

    @property
    def success(self) -> bool:
        """True only when the solver converged and the residual passed threshold."""
        return self.solver_converged and self.residual_accepted


@dataclass(frozen=True, slots=True)
class AlgebraicFailure:
    """Structured failure before a usable numerical candidate was produced."""

    problem_id: str
    solver_id: str
    message: str
    error_type: str
    last_candidate: tuple[float, ...] | None
    residual_vector: tuple[float, ...] | None
    function_evaluations: int
    configuration: AlgebraicSolverConfiguration
    reproducibility: AlgebraicReproducibilityRecord


AlgebraicOutcome = AlgebraicResult | AlgebraicFailure


class AlgebraicSolver(Protocol):
    """Protocol implemented by a solver for nonlinear algebraic systems."""

    capabilities: AlgebraicSolverCapabilities

    def solve(
        self,
        problem: NonlinearAlgebraicProblem,
        configuration: AlgebraicSolverConfiguration,
    ) -> AlgebraicOutcome:
        """Solve one explicitly compatible residual system."""


_DEFAULT_ALGEBRAIC_CONFIGURATION = AlgebraicSolverConfiguration()


def _scipy_version() -> str:
    try:
        return version("scipy")
    except PackageNotFoundError:
        return "unknown"


class ScipyHybridRootSolver:
    """Solve square systems using SciPy's local Powell hybrid method."""

    capabilities = AlgebraicSolverCapabilities(solver_version=_scipy_version())

    def solve(
        self,
        problem: NonlinearAlgebraicProblem,
        configuration: AlgebraicSolverConfiguration = _DEFAULT_ALGEBRAIC_CONFIGURATION,
    ) -> AlgebraicOutcome:
        """Return a candidate with separate solver and residual acceptance flags.

        Uses ``scipy.optimize.root(method='hybr')``. ``step_tolerance`` is
        passed as the method's relative iterate-change tolerance (``xtol``);
        ``residual_tolerance`` is independently applied to the final L2 norm.
        ``maximum_function_evaluations`` maps to SciPy's ``maxfev`` option.
        """
        reproducibility = _make_reproducibility(
            problem, configuration, self.capabilities
        )
        evaluation_count = 0
        initial_guess = problem.initial_guess
        evaluation_count += 1
        try:
            problem.evaluate_residual(initial_guess)
        except ScientificValidationError as exc:
            return AlgebraicFailure(
                problem_id=problem.problem_id,
                solver_id=self.capabilities.solver_id,
                message=str(exc),
                error_type=type(exc).__name__,
                last_candidate=tuple(float(value) for value in initial_guess),
                residual_vector=None,
                function_evaluations=evaluation_count,
                configuration=configuration,
                reproducibility=reproducibility,
            )

        last_valid_residual: NDArray[np.float64] | None = None
        last_evaluated_values = initial_guess.copy()

        def residual_function(values: NDArray[np.float64]) -> NDArray[np.float64]:
            nonlocal evaluation_count, last_evaluated_values, last_valid_residual
            evaluation_count += 1
            last_evaluated_values = np.asarray(values, dtype=np.float64).copy()
            last_valid_residual = problem.evaluate_residual(values)
            return last_valid_residual

        try:
            solution = root(
                residual_function,
                initial_guess,
                method=configuration.method,
                options={
                    "xtol": configuration.step_tolerance,
                    "maxfev": configuration.maximum_function_evaluations,
                },
            )
            candidate = np.asarray(solution.x, dtype=np.float64)
            final_residual = problem.evaluate_residual(candidate)
            evaluation_count += 1
        except ScientificValidationError as exc:
            return AlgebraicFailure(
                problem_id=problem.problem_id,
                solver_id=self.capabilities.solver_id,
                message=str(exc),
                error_type=type(exc).__name__,
                last_candidate=tuple(float(value) for value in last_evaluated_values),
                residual_vector=(
                    None
                    if last_valid_residual is None
                    else tuple(float(value) for value in last_valid_residual)
                ),
                function_evaluations=evaluation_count,
                configuration=configuration,
                reproducibility=reproducibility,
            )
        except (ArithmeticError, RuntimeError, ValueError, TypeError) as exc:
            return AlgebraicFailure(
                problem_id=problem.problem_id,
                solver_id=self.capabilities.solver_id,
                message=f"SciPy root solve failed: {exc}",
                error_type=type(exc).__name__,
                last_candidate=tuple(float(value) for value in last_evaluated_values),
                residual_vector=(
                    None
                    if last_valid_residual is None
                    else tuple(float(value) for value in last_valid_residual)
                ),
                function_evaluations=evaluation_count,
                configuration=configuration,
                reproducibility=reproducibility,
            )

        residual_norm = float(np.linalg.norm(final_residual, ord=2))
        return AlgebraicResult(
            problem_id=problem.problem_id,
            solver_id=self.capabilities.solver_id,
            method=configuration.method,
            solution_candidate=candidate.copy(),
            variable_names=tuple(unknown.name for unknown in problem.unknowns),
            variable_units=tuple(unknown.unit for unknown in problem.unknowns),
            residual_vector=final_residual.copy(),
            residual_names=problem.residual_names,
            residual_units=problem.residual_units,
            residual_norm=residual_norm,
            residual_norm_definition=(
                "Euclidean L2 norm of raw residual entries; no unit conversion "
                "or scaling is applied"
            ),
            solver_converged=bool(solution.success),
            residual_accepted=residual_norm <= configuration.residual_tolerance,
            physical_validity=AlgebraicPhysicalValidity.NOT_ASSESSED,
            uniqueness_assessed=False,
            solution_set_completeness_assessed=False,
            termination_reason=str(solution.message),
            # scipy.optimize.root(method='hybr') does not expose iteration count.
            iteration_count=None,
            function_evaluations=evaluation_count,
            configuration=configuration,
            parameters=problem.parameters,
            assumptions=problem.assumptions,
            constraints=problem.constraints,
            validity_limits=problem.validity_limits,
            reproducibility=reproducibility,
        )


def _make_reproducibility(
    problem: NonlinearAlgebraicProblem,
    configuration: AlgebraicSolverConfiguration,
    capabilities: AlgebraicSolverCapabilities,
) -> AlgebraicReproducibilityRecord:
    return AlgebraicReproducibilityRecord(
        python_version=sys.version,
        platform=platform.platform(),
        numpy_version=np.__version__,
        scipy_version=_scipy_version(),
        problem_id=problem.problem_id,
        model_version=problem.model_version,
        initial_guess=tuple(float(value) for value in problem.initial_guess),
        variable_names=tuple(unknown.name for unknown in problem.unknowns),
        variable_units=tuple(unknown.unit for unknown in problem.unknowns),
        parameters=problem.parameters,
        solver_id=capabilities.solver_id,
        solver_version=capabilities.solver_version,
        method=configuration.method,
        step_tolerance=configuration.step_tolerance,
        residual_tolerance=configuration.residual_tolerance,
        maximum_function_evaluations=configuration.maximum_function_evaluations,
        equation_record_ids=problem.equation_record_ids,
        source_reference_ids=problem.source_reference_ids,
        assumptions=problem.assumptions,
        constraints=problem.constraints,
        validity_limits=problem.validity_limits,
    )


def make_nonlinear_spring_equilibrium_problem(
    *,
    linear_stiffness_n_per_m: float = 100.0,
    cubic_stiffness_n_per_m3: float = 1000.0,
    applied_force_n: float = 11.0,
    initial_displacement_m: float = 0.08,
) -> NonlinearAlgebraicProblem:
    """Describe equilibrium of a one-dimensional cubic spring under load.

    The residual is ``k*x + a*x**3 - f``. With the defaults, ``x=0.1 m`` is
    an independently checkable root: ``100*0.1 + 1000*0.1**3 = 11 N``.
    Positive ``k`` and ``a`` make this particular scalar force law monotonic,
    but this property is specific to these stated parameter signs and does not
    imply general root-solver uniqueness guarantees.
    """
    for name, value in (
        ("linear_stiffness_n_per_m", linear_stiffness_n_per_m),
        ("cubic_stiffness_n_per_m3", cubic_stiffness_n_per_m3),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ScientificValidationError(f"{name} must be finite and positive")
    if not np.isfinite(applied_force_n):
        raise ScientificValidationError("applied_force_n must be finite")
    if not np.isfinite(initial_displacement_m):
        raise ScientificValidationError("initial_displacement_m must be finite")

    def force_balance(
        unknowns: NDArray[np.float64], parameters: Mapping[str, float]
    ) -> NDArray[np.float64]:
        displacement = unknowns[0]
        residual = (
            parameters["linear_stiffness_n_per_m"] * displacement
            + parameters["cubic_stiffness_n_per_m3"] * displacement**3
            - parameters["applied_force_n"]
        )
        return np.asarray([residual], dtype=np.float64)

    return NonlinearAlgebraicProblem(
        problem_id="cubic_spring_equilibrium",
        name="Loaded cubic-spring equilibrium",
        description="Static force-balance candidate for a one-dimensional spring.",
        unknowns=(
            AlgebraicUnknown(
                name="displacement",
                unit="m",
                initial_guess=initial_displacement_m,
            ),
        ),
        residual_function=force_balance,
        residual_names=("force_balance",),
        residual_units=("N",),
        parameters=(
            AlgebraicParameter(
                name="linear_stiffness_n_per_m",
                value=linear_stiffness_n_per_m,
                unit="N/m",
            ),
            AlgebraicParameter(
                name="cubic_stiffness_n_per_m3",
                value=cubic_stiffness_n_per_m3,
                unit="N/m^3",
            ),
            AlgebraicParameter(name="applied_force_n", value=applied_force_n, unit="N"),
        ),
        assumptions=(
            "One-dimensional quasi-static force balance.",
            "Restoring force is k*x + a*x^3 with constant coefficients.",
        ),
        constraints=("Displacement must lie within the calibrated spring range.",),
        validity_limits=(
            "The cubic constitutive law is illustrative and must be calibrated "
            "before use for a real spring.",
        ),
        model_version="1",
    )
