"""Reproducible experiment records, deterministic sweeps, and sensitivities.

The framework stores typed scalar metrics and adapter-level solver evidence.
Adapters retain model-specific parameter and metric semantics; current
end-to-end adapters cover the damped oscillator, nonlinear pendulum, and
nonlinear algebraic root problem.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from enum import StrEnum
from importlib.metadata import PackageNotFoundError, version
from typing import Literal, Protocol

import numpy as np
from pydantic import Field, field_validator, model_validator

from newton_lab.algebraic import (
    AlgebraicFailure,
    AlgebraicParameter,
    AlgebraicResult,
    AlgebraicSolver,
    AlgebraicSolverConfiguration,
    NonlinearAlgebraicProblem,
    ScipyHybridRootSolver,
    make_nonlinear_spring_equilibrium_problem,
)
from newton_lab.boundary_value import (
    BoundaryValueFailure,
    BoundaryValueProblem,
    BoundaryValueResult,
    BoundaryValueSolverConfiguration,
    ScipyBoundaryValueSolver,
    make_linear_heat_conduction_problem,
)
from newton_lab.dynamics import DampedOscillatorResult, simulate_damped_oscillator
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.pendulum import DampedPendulumResult, simulate_damped_pendulum
from newton_lab.simulation import (
    NamedValue,
    ODESolverConfiguration,
    OutputSampling,
    SimulationModel,
)

ExperimentSolverConfiguration = (
    ODESolverConfiguration
    | AlgebraicSolverConfiguration
    | BoundaryValueSolverConfiguration
)


class ParameterSweep(SimulationModel):
    """Ordered values for one named parameter; failures do not reorder points."""

    parameter_name: str = Field(min_length=1)
    values: tuple[float, ...] = Field(min_length=1)

    @field_validator("values")
    @classmethod
    def finite_values(cls, values: tuple[float, ...]) -> tuple[float, ...]:
        if not np.isfinite(values).all():
            raise ValueError("sweep values must be finite")
        return values


class AnalyticalReferenceRequest(SimulationModel):
    """Request an adapter-provided analytical scalar comparison."""

    metric_id: str = Field(min_length=1)
    absolute_tolerance: float = Field(gt=0.0)

    @field_validator("absolute_tolerance")
    @classmethod
    def finite_tolerance(cls, value: float) -> float:
        if not np.isfinite(value):
            raise ValueError("absolute_tolerance must be finite")
        return value


class AlgebraicParameterConstraint(SimulationModel):
    """A declared domain constraint for one algebraic experiment parameter."""

    parameter_name: str = Field(min_length=1)
    domain: Literal["positive", "non_negative"]


class BVPParameterConstraint(SimulationModel):
    """A declared domain constraint for one BVP experiment parameter."""

    parameter_name: str = Field(min_length=1)
    domain: Literal["positive", "non_negative"]


class BVPReferenceRequest(SimulationModel):
    """Request a qualified spatial comparison with an available BVP reference."""

    maximum_absolute_tolerance: float = Field(gt=0.0)

    @field_validator("maximum_absolute_tolerance")
    @classmethod
    def finite_tolerance(cls, value: float) -> float:
        if not np.isfinite(value):
            raise ValueError("maximum_absolute_tolerance must be finite")
        return value


class ExperimentSpecification(SimulationModel):
    """Validated experiment inputs independent of a solver's result type."""

    experiment_id: str = Field(pattern=r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
    name: str = Field(min_length=1)
    description: str = ""
    model_id: str = Field(min_length=1)
    configuration_version: Literal["1"] = "1"
    parameters: tuple[NamedValue, ...] = Field(min_length=1)
    sweep: ParameterSweep | None = None
    metrics: tuple[str, ...] = Field(min_length=1)
    solver_configuration: ExperimentSolverConfiguration = ODESolverConfiguration()
    output_sampling: OutputSampling = OutputSampling()
    algebraic_problem: NonlinearAlgebraicProblem | None = None
    algebraic_parameter_constraints: tuple[AlgebraicParameterConstraint, ...] = ()
    algebraic_reference_solution: tuple[float, ...] | None = None
    bvp_problem: BoundaryValueProblem | None = None
    bvp_parameter_constraints: tuple[BVPParameterConstraint, ...] = ()
    bvp_reference: BVPReferenceRequest | None = None
    analytical_reference: AnalyticalReferenceRequest | None = None
    assumptions: tuple[str, ...] = ()
    equation_record_ids: tuple[str, ...] = ()
    source_reference_ids: tuple[str, ...] = ()
    retain_ode_trajectory: bool = False

    @model_validator(mode="after")
    def validate_experiment_references(self) -> ExperimentSpecification:
        parameter_names = tuple(item.name for item in self.parameters)
        if len(parameter_names) != len(set(parameter_names)):
            raise ValueError("experiment parameter names must be unique")
        if len(self.metrics) != len(set(self.metrics)):
            raise ValueError("requested metric IDs must be unique")
        if self.sweep is not None and self.sweep.parameter_name not in parameter_names:
            raise ValueError("sweep parameter must be present in experiment parameters")
        if self.bvp_problem is not None:
            if self.algebraic_problem is not None:
                raise ValueError("an experiment can configure only one solver family")
            if self.algebraic_parameter_constraints:
                raise ValueError("algebraic constraints require an algebraic problem")
            if self.model_id != BVPExperimentAdapter.model_id:
                raise ValueError("BVP problems require the BVP model ID")
            if self.algebraic_reference_solution is not None:
                raise ValueError(
                    "algebraic reference solution requires an algebraic problem"
                )
            if self.analytical_reference is not None:
                raise ValueError(
                    "BVP experiments use the spatial bvp_reference comparison request"
                )
            if not isinstance(
                self.solver_configuration, BoundaryValueSolverConfiguration
            ):
                raise ValueError(
                    "BVP experiments require BoundaryValueSolverConfiguration"
                )
            bvp_parameters = self.bvp_problem.parameters
            if parameter_names != tuple(item.name for item in bvp_parameters):
                raise ValueError(
                    "experiment parameters must match BVP problem parameter order"
                )
            constrained_names = tuple(
                constraint.parameter_name
                for constraint in self.bvp_parameter_constraints
            )
            if len(constrained_names) != len(set(constrained_names)) or any(
                name not in parameter_names for name in constrained_names
            ):
                raise ValueError(
                    "BVP parameter constraints must uniquely reference "
                    "declared parameters"
                )
            if any(
                named.value != parameter.value
                or (parameter.unit is not None and named.unit != parameter.unit)
                for named, parameter in zip(
                    self.parameters, bvp_parameters, strict=True
                )
            ):
                raise ValueError(
                    "experiment parameter values and units must match BVP problem"
                )
        elif self.algebraic_problem is None:
            if (
                self.bvp_parameter_constraints
                or self.bvp_reference is not None
                or self.algebraic_parameter_constraints
            ):
                raise ValueError("family-specific settings require their problem type")
            if not isinstance(self.solver_configuration, ODESolverConfiguration):
                raise ValueError("ODE experiments require ODESolverConfiguration")
            if self.algebraic_reference_solution is not None:
                raise ValueError(
                    "algebraic reference solution requires an algebraic problem"
                )
        else:
            if self.bvp_parameter_constraints or self.bvp_reference is not None:
                raise ValueError("BVP settings require a BVP problem")
            if self.model_id != "nonlinear_algebraic_system":
                raise ValueError("algebraic problems require the algebraic model ID")
            if not isinstance(self.solver_configuration, AlgebraicSolverConfiguration):
                raise ValueError(
                    "algebraic experiments require AlgebraicSolverConfiguration"
                )
            algebraic_parameters = self.algebraic_problem.parameters
            if parameter_names != tuple(item.name for item in algebraic_parameters):
                raise ValueError(
                    "experiment parameters must match algebraic problem parameter order"
                )
            constrained_names = tuple(
                constraint.parameter_name
                for constraint in self.algebraic_parameter_constraints
            )
            if len(constrained_names) != len(set(constrained_names)) or any(
                name not in parameter_names for name in constrained_names
            ):
                raise ValueError(
                    "algebraic parameter constraints must uniquely reference "
                    "declared parameters"
                )
            if any(
                named.value != parameter.value
                or (parameter.unit is not None and named.unit != parameter.unit)
                for named, parameter in zip(
                    self.parameters, algebraic_parameters, strict=True
                )
            ):
                raise ValueError(
                    "experiment parameter values and units must match algebraic problem"
                )
            if self.algebraic_reference_solution is not None and (
                len(self.algebraic_reference_solution)
                != len(self.algebraic_problem.unknowns)
                or not np.isfinite(self.algebraic_reference_solution).all()
            ):
                raise ValueError(
                    "algebraic reference solution must be finite and match "
                    "unknown count"
                )
        if (
            self.analytical_reference is not None
            and self.analytical_reference.metric_id not in self.metrics
        ):
            raise ValueError("analytical reference metric must be requested")
        for label, values in (
            ("equation record IDs", self.equation_record_ids),
            ("source reference IDs", self.source_reference_ids),
        ):
            if len(values) != len(set(values)):
                raise ValueError(f"{label} must be unique")
        if self.retain_ode_trajectory and self.model_id not in {
            OscillatorExperimentAdapter.model_id,
            PendulumExperimentAdapter.model_id,
        }:
            raise ValueError(
                "ODE trajectory retention is available only for ODE models"
            )
        return self


class ExperimentFailureKind(StrEnum):
    MODEL_VALIDATION = "model_validation_failure"
    SOLVER = "solver_failure"
    METRIC_EXTRACTION = "metric_extraction_failure"
    RESIDUAL_ACCEPTANCE = "residual_acceptance_failure"


class ExperimentRunStatus(StrEnum):
    SUCCEEDED = "succeeded"
    FAILED = "failed"


class ExperimentRecordStatus(StrEnum):
    SUCCEEDED = "succeeded"
    COMPLETED_WITH_FAILURES = "completed_with_failures"


class MetricValue(SimulationModel):
    """A finite scalar metric with its declared unit."""

    metric_id: str = Field(min_length=1)
    value: float
    unit: str = Field(min_length=1)

    @field_validator("value")
    @classmethod
    def finite_metric(cls, value: float) -> float:
        if isinstance(value, bool) or not np.isfinite(value):
            raise ValueError("metric value must be finite")
        return value


class AnalyticalReferenceComparison(SimulationModel):
    """Scalar comparison; a pass is not a physical-validity assessment."""

    status: Literal["passed", "failed", "unavailable"]
    metric_id: str
    unit: str
    numerical_value: float | None = None
    reference_value: float | None = None
    absolute_error: float | None = None
    absolute_tolerance: float
    reference_description: str = ""
    message: str = ""

    @field_validator(
        "numerical_value", "reference_value", "absolute_error", "absolute_tolerance"
    )
    @classmethod
    def finite_comparison_values(cls, value: float | None) -> float | None:
        if value is not None and not np.isfinite(value):
            raise ValueError("reference comparison values must be finite")
        return value


class ExperimentReproducibility(SimulationModel):
    """Configuration/environment record for interpreting and reproducing a run."""

    experiment_id: str
    configuration_version: str
    model_id: str
    parameters: tuple[NamedValue, ...]
    assumptions: tuple[str, ...]
    equation_record_ids: tuple[str, ...]
    source_reference_ids: tuple[str, ...]
    requested_metrics: tuple[str, ...]
    solver_id: str
    solver_method: str
    problem_id: str | None = None
    initial_guess: tuple[float, ...] | None = None
    relative_tolerance: float | None = None
    absolute_tolerance: float | None = None
    maximum_step: float | None = None
    step_tolerance: float | None = None
    residual_tolerance: float | None = None
    maximum_function_evaluations: int | None = None
    output_point_count: int | None = None
    bvp_domain: tuple[float, float] | None = None
    bvp_mesh: tuple[float, ...] | None = None
    bvp_initial_guess: tuple[tuple[float, ...], tuple[float, ...]] | None = None
    bvp_solver_tolerance: float | None = None
    bvp_boundary_tolerance: float | None = None
    bvp_acceptance_tolerance: float | None = None
    bvp_maximum_nodes: int | None = None
    python_version: str
    numpy_version: str
    scipy_version: str
    random_seed: int | None = None


class AlgebraicRunDiagnostics(SimulationModel):
    """Candidate and solver evidence specific to an algebraic root solve."""

    problem_id: str
    initial_guess: tuple[float, ...]
    variable_names: tuple[str, ...]
    variable_units: tuple[str | None, ...]
    solution_candidate: tuple[float, ...] | None = None
    residual_vector: tuple[float, ...] | None = None
    residual_names: tuple[str, ...] | None = None
    residual_units: tuple[str | None, ...] | None = None
    residual_norm: float | None = None
    residual_norm_definition: str
    solver_converged: bool | None = None
    residual_accepted: bool | None = None
    residual_tolerance: float
    physical_validity: Literal["not_assessed"] = "not_assessed"
    uniqueness_assessed: Literal[False] = False
    solution_set_completeness_assessed: Literal[False] = False
    termination_reason: str
    iteration_count: int | None = None
    function_evaluations: int

    @model_validator(mode="after")
    def candidate_dimensions_match(self) -> AlgebraicRunDiagnostics:
        dimension = len(self.initial_guess)
        if dimension == 0 or len(self.variable_names) != dimension:
            raise ValueError("algebraic diagnostics need one name per initial guess")
        if len(self.variable_units) != dimension:
            raise ValueError("variable units must match the initial guess dimension")
        if (
            self.solution_candidate is not None
            and len(self.solution_candidate) != dimension
        ):
            raise ValueError("solution candidate must match unknown dimension")
        if self.residual_vector is not None and len(self.residual_vector) != dimension:
            raise ValueError("residual vector must match square system dimension")
        numeric = (
            *self.initial_guess,
            *(self.solution_candidate or ()),
            *(self.residual_vector or ()),
        )
        if not np.isfinite(numeric).all() or (
            self.residual_norm is not None and not np.isfinite(self.residual_norm)
        ):
            raise ValueError("algebraic diagnostics values must be finite")
        return self


class BVPRunDiagnostics(SimulationModel):
    """Serializable spatial candidate and solver evidence for one BVP run."""

    problem_id: str
    domain: tuple[float, float]
    domain_unit: str
    mesh: tuple[float, ...]
    state_names: tuple[str, str]
    state_units: tuple[str | None, str | None]
    solution: tuple[tuple[float, ...], tuple[float, ...]]
    boundary_residuals: tuple[float, float]
    boundary_residual_units: tuple[str | None, str | None]
    differential_residuals: tuple[float, ...]
    solver_success: bool
    boundary_conditions_satisfied: bool
    numerical_accepted: bool
    status_code: int
    message: str
    iteration_count: int
    solver_tolerance: float
    boundary_tolerance: float
    acceptance_tolerance: float
    physical_validity_assessed: bool = False
    uniqueness_assessed: bool = False
    solution_set_completeness_assessed: bool = False

    @model_validator(mode="after")
    def validate_spatial_data(self) -> BVPRunDiagnostics:
        if len(self.mesh) < 2 or not np.isfinite(self.mesh).all():
            raise ValueError("BVP diagnostic mesh must be finite and nonempty")
        if not np.all(np.diff(self.mesh) > 0):
            raise ValueError("BVP diagnostic mesh must be strictly increasing")
        if (
            not np.isfinite(self.domain).all()
            or self.domain[1] <= self.domain[0]
            or self.mesh[0] != self.domain[0]
            or self.mesh[-1] != self.domain[1]
        ):
            raise ValueError("BVP mesh endpoints must match its finite domain")
        if any(len(component) != len(self.mesh) for component in self.solution):
            raise ValueError("BVP solution components must match the spatial mesh")
        numeric = (
            *self.domain,
            *self.mesh,
            *self.solution[0],
            *self.solution[1],
            *self.boundary_residuals,
            *self.differential_residuals,
        )
        if not np.isfinite(numeric).all():
            raise ValueError("BVP diagnostics must contain finite numerical data")
        expected_boundary_status = (
            max(abs(value) for value in self.boundary_residuals)
            <= self.acceptance_tolerance
        )
        expected_acceptance = (
            self.solver_success
            and expected_boundary_status
            and self.boundary_conditions_satisfied
            and bool(self.differential_residuals)
            and max(self.differential_residuals) <= self.acceptance_tolerance
        )
        if (
            self.boundary_conditions_satisfied != expected_boundary_status
            or self.numerical_accepted != expected_acceptance
        ):
            raise ValueError("BVP numerical acceptance flags are inconsistent")
        return self


class BVPReferenceComparison(SimulationModel):
    """Comparison with the known linear profile for the heat example only."""

    status: Literal["passed", "failed", "unavailable"]
    method: Literal["direct_evaluation_on_returned_mesh"]
    state_name: str
    unit: str
    maximum_absolute_error: float | None = None
    rms_error: float | None = None
    maximum_absolute_tolerance: float
    reference_description: str
    message: str

    @field_validator(
        "maximum_absolute_error", "rms_error", "maximum_absolute_tolerance"
    )
    @classmethod
    def finite_reference_metrics(cls, value: float | None) -> float | None:
        if value is not None and not np.isfinite(value):
            raise ValueError("BVP reference comparison values must be finite")
        return value


class ODETrajectorySnapshot(SimulationModel):
    """Immutable copy of the requested ODE output samples, with SI units."""

    model_id: str
    time_s: tuple[float, ...] = Field(min_length=2)
    state_names: tuple[str, ...] = Field(min_length=1)
    state_units: tuple[str, ...] = Field(min_length=1)
    state_values: tuple[tuple[float, ...], ...] = Field(min_length=1)
    solver_method: str

    @model_validator(mode="after")
    def validate_aligned_samples(self) -> ODETrajectorySnapshot:
        if len(self.state_names) != len(self.state_units) or len(
            self.state_names
        ) != len(self.state_values):
            raise ValueError("trajectory state metadata and values must align")
        if not np.isfinite(self.time_s).all() or not np.all(np.diff(self.time_s) > 0.0):
            raise ValueError("trajectory times must be finite and strictly increasing")
        if any(len(values) != len(self.time_s) for values in self.state_values):
            raise ValueError("trajectory state samples must match the time grid")
        if any(not np.isfinite(values).all() for values in self.state_values):
            raise ValueError("trajectory state samples must be finite")
        return self


class ExperimentRun(SimulationModel):
    """One sweep point, retaining either metrics or a typed expected failure."""

    run_index: int = Field(ge=0)
    swept_parameter_name: str | None = None
    swept_parameter_value: float | None = None
    parameters: tuple[NamedValue, ...]
    status: ExperimentRunStatus
    metrics: tuple[MetricValue, ...] = ()
    failure_kind: ExperimentFailureKind | None = None
    error_type: str | None = None
    message: str = ""
    solver_id: str
    solver_method: str
    solver_success: bool | None = None
    numerical_acceptance: bool | None = None
    function_evaluations: int | None = Field(default=None, ge=0)
    reference_comparison: AnalyticalReferenceComparison | None = None
    algebraic_diagnostics: AlgebraicRunDiagnostics | None = None
    bvp_diagnostics: BVPRunDiagnostics | None = None
    bvp_reference_comparison: BVPReferenceComparison | None = None
    ode_trajectory: ODETrajectorySnapshot | None = None
    reproducibility: ExperimentReproducibility

    @model_validator(mode="after")
    def status_matches_outcome(self) -> ExperimentRun:
        if self.status == ExperimentRunStatus.SUCCEEDED:
            if self.failure_kind is not None or self.error_type is not None:
                raise ValueError("successful runs cannot carry failure fields")
        elif self.failure_kind is None or self.error_type is None:
            raise ValueError("failed runs require failure kind and error type")
        if len({metric.metric_id for metric in self.metrics}) != len(self.metrics):
            raise ValueError("run metric IDs must be unique")
        return self

    @field_validator("swept_parameter_value")
    @classmethod
    def finite_sweep_value(cls, value: float | None) -> float | None:
        if value is not None and not np.isfinite(value):
            raise ValueError("swept parameter value must be finite")
        return value


class ExperimentRecord(SimulationModel):
    """Ordered aggregate for one experiment execution."""

    specification: ExperimentSpecification
    runs: tuple[ExperimentRun, ...] = Field(min_length=1)
    status: ExperimentRecordStatus

    @model_validator(mode="after")
    def status_matches_runs(self) -> ExperimentRecord:
        expected = (
            ExperimentRecordStatus.SUCCEEDED
            if all(run.status == ExperimentRunStatus.SUCCEEDED for run in self.runs)
            else ExperimentRecordStatus.COMPLETED_WITH_FAILURES
        )
        if self.status != expected:
            raise ValueError("experiment status does not match its run outcomes")
        if tuple(run.run_index for run in self.runs) != tuple(range(len(self.runs))):
            raise ValueError("run indices must preserve contiguous execution order")
        for run in self.runs:
            trajectory = run.ode_trajectory
            if trajectory is not None and (
                trajectory.model_id != self.specification.model_id
                or not self.specification.retain_ode_trajectory
                or run.status != ExperimentRunStatus.SUCCEEDED
            ):
                raise ValueError(
                    "ODE trajectories must match retained successful ODE runs"
                )
            if (
                self.specification.retain_ode_trajectory
                and run.status == ExperimentRunStatus.SUCCEEDED
                and trajectory is None
            ):
                raise ValueError(
                    "successful retained ODE runs require their trajectory"
                )
        return self


class AdapterSuccess(SimulationModel):
    """Family adapter's extracted metrics and available solver diagnostics."""

    metrics: tuple[MetricValue, ...]
    solver_id: str
    solver_method: str
    solver_success: bool
    numerical_acceptance: bool | None = None
    function_evaluations: int | None = Field(default=None, ge=0)
    algebraic_diagnostics: AlgebraicRunDiagnostics | None = None
    bvp_diagnostics: BVPRunDiagnostics | None = None
    ode_trajectory: ODETrajectorySnapshot | None = None


class AdapterFailure(SimulationModel):
    """An expected failure classified by the model/solver adapter."""

    failure_kind: ExperimentFailureKind
    error_type: str
    message: str
    solver_id: str
    solver_method: str
    solver_success: bool | None = None
    numerical_acceptance: bool | None = None
    function_evaluations: int | None = Field(default=None, ge=0)
    metrics: tuple[MetricValue, ...] = ()
    algebraic_diagnostics: AlgebraicRunDiagnostics | None = None
    bvp_diagnostics: BVPRunDiagnostics | None = None


AdapterOutcome = AdapterSuccess | AdapterFailure
_DEFAULT_EXPERIMENT_SOLVER_CONFIGURATION = ODESolverConfiguration()
_DEFAULT_EXPERIMENT_OUTPUT_SAMPLING = OutputSampling()
_DEFAULT_EXPERIMENT_ALGEBRAIC_CONFIGURATION = AlgebraicSolverConfiguration()
_DEFAULT_EXPERIMENT_BVP_CONFIGURATION = BoundaryValueSolverConfiguration()


class AdapterAnalyticalReference(SimulationModel):
    """Adapter reference scalar and a note describing its interpretation."""

    value: float
    unit: str
    description: str

    @field_validator("value")
    @classmethod
    def finite_reference(cls, value: float) -> float:
        if not np.isfinite(value):
            raise ValueError("analytical reference must be finite")
        return value


class ExperimentAdapter(Protocol):
    """Narrow adapter protocol that does not force solver result interchangeability."""

    model_id: str
    solver_id: str

    def execute(
        self,
        parameters: tuple[NamedValue, ...],
        metrics: tuple[str, ...],
        solver_configuration: ExperimentSolverConfiguration,
        output_sampling: OutputSampling,
    ) -> AdapterOutcome:
        """Run one model point and return scalar metrics or a classified failure."""

    def analytical_reference(
        self, parameters: tuple[NamedValue, ...], metric_id: str
    ) -> AdapterAnalyticalReference | None:
        """Return a trusted scalar reference and unit when this adapter has one."""


class OscillatorExperimentAdapter:
    """Experiment adapter for the existing damped oscillator function."""

    model_id = "damped_harmonic_oscillator"
    solver_id = "scipy_solve_ivp"
    _parameter_names = (
        "mass_kg",
        "damping_coefficient_kg_per_s",
        "stiffness_n_per_m",
        "initial_displacement_m",
        "initial_velocity_m_per_s",
        "duration_s",
    )
    _metrics = {
        "maximum_absolute_displacement_m": "m",
        "final_displacement_m": "m",
        "final_velocity_m_per_s": "m/s",
        "maximum_absolute_velocity_m_per_s": "m/s",
    }

    def execute(
        self,
        parameters: tuple[NamedValue, ...],
        metrics: tuple[str, ...],
        solver_configuration: ExperimentSolverConfiguration,
        output_sampling: OutputSampling,
        *,
        retain_trajectory: bool = False,
    ) -> AdapterOutcome:
        """Call the existing model, classifying expected input and solver errors."""
        if not isinstance(solver_configuration, ODESolverConfiguration):
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type="ConfigurationTypeError",
                message="oscillator experiments require ODESolverConfiguration",
                solver_id=self.solver_id,
                solver_method=solver_configuration.method,
            )
        solver_method = solver_configuration.method
        values = {item.name: item.value for item in parameters}
        try:
            if set(values) != set(self._parameter_names):
                missing = set(self._parameter_names) - set(values)
                extra = set(values) - set(self._parameter_names)
                raise ScientificValidationError(
                    f"oscillator parameters mismatch; missing={sorted(missing)}, "
                    f"extra={sorted(extra)}"
                )
            unknown_metrics = set(metrics) - set(self._metrics)
            if unknown_metrics:
                raise ScientificValidationError(
                    f"unsupported oscillator metrics: {sorted(unknown_metrics)}"
                )
        except ScientificValidationError as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_method,
            )

        try:
            result = simulate_damped_oscillator(
                mass_kg=values["mass_kg"],
                damping_coefficient_kg_per_s=values["damping_coefficient_kg_per_s"],
                stiffness_n_per_m=values["stiffness_n_per_m"],
                initial_displacement_m=values["initial_displacement_m"],
                initial_velocity_m_per_s=values["initial_velocity_m_per_s"],
                duration_s=values["duration_s"],
                num_points=output_sampling.point_count,
                method=solver_method,
                relative_tolerance=solver_configuration.relative_tolerance,
                absolute_tolerance=solver_configuration.absolute_tolerance,
                maximum_step_s=solver_configuration.maximum_step,
            )
        except ScientificValidationError as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_method,
            )
        except IntegrationError as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.SOLVER,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_method,
                solver_success=False,
            )

        try:
            extracted = _extract_oscillator_metrics(result, metrics)
        except (ScientificValidationError, ValueError, FloatingPointError) as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.METRIC_EXTRACTION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_method,
                solver_success=True,
            )
        return AdapterSuccess(
            metrics=extracted,
            solver_id=self.solver_id,
            solver_method=solver_method,
            solver_success=True,
            # The legacy direct solver exposes local tolerances, not a separate
            # experiment-level numerical acceptance flag.
            numerical_acceptance=None,
            function_evaluations=result.function_evaluations,
            ode_trajectory=(
                _oscillator_trajectory(result) if retain_trajectory else None
            ),
        )

    def analytical_reference(
        self, parameters: tuple[NamedValue, ...], metric_id: str
    ) -> AdapterAnalyticalReference | None:
        """Exact undamped oscillator displacement at the final simulation time."""
        if metric_id != "final_displacement_m":
            return None
        values = {item.name: item.value for item in parameters}
        if values["damping_coefficient_kg_per_s"] != 0.0:
            return None
        angular_frequency = np.sqrt(values["stiffness_n_per_m"] / values["mass_kg"])
        final_time = values["duration_s"]
        displacement = values["initial_displacement_m"] * np.cos(
            angular_frequency * final_time
        ) + values["initial_velocity_m_per_s"] / angular_frequency * np.sin(
            angular_frequency * final_time
        )
        return AdapterAnalyticalReference(
            value=float(displacement),
            unit="m",
            description="Exact analytical solution for the undamped linear oscillator.",
        )


class PendulumExperimentAdapter:
    """Adapter for the existing full-sine damped nonlinear pendulum model."""

    model_id = "damped_nonlinear_pendulum"
    solver_id = "scipy_solve_ivp"
    _parameter_names = (
        "mass_kg",
        "length_m",
        "damping_coefficient_kg_m2_per_s",
        "gravity_m_per_s2",
        "initial_angle_rad",
        "initial_angular_velocity_rad_per_s",
        "duration_s",
    )
    _metrics = {
        "final_angle_rad": "rad",
        "final_angular_velocity_rad_per_s": "rad/s",
        "maximum_absolute_angle_rad": "rad",
        "maximum_absolute_angular_velocity_rad_per_s": "rad/s",
    }

    def execute(
        self,
        parameters: tuple[NamedValue, ...],
        metrics: tuple[str, ...],
        solver_configuration: ExperimentSolverConfiguration,
        output_sampling: OutputSampling,
        *,
        retain_trajectory: bool = False,
    ) -> AdapterOutcome:
        """Delegate to the existing nonlinear pendulum simulation function."""
        if not isinstance(solver_configuration, ODESolverConfiguration):
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type="ConfigurationTypeError",
                message="pendulum experiments require ODESolverConfiguration",
                solver_id=self.solver_id,
                solver_method=solver_configuration.method,
            )
        solver_method = solver_configuration.method
        values = {item.name: item.value for item in parameters}
        try:
            if set(values) != set(self._parameter_names):
                missing = set(self._parameter_names) - set(values)
                extra = set(values) - set(self._parameter_names)
                raise ScientificValidationError(
                    f"pendulum parameters mismatch; missing={sorted(missing)}, "
                    f"extra={sorted(extra)}"
                )
            unknown_metrics = set(metrics) - set(self._metrics)
            if unknown_metrics:
                raise ScientificValidationError(
                    f"unsupported pendulum metrics: {sorted(unknown_metrics)}"
                )
        except ScientificValidationError as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_method,
            )

        try:
            result = simulate_damped_pendulum(
                mass_kg=values["mass_kg"],
                length_m=values["length_m"],
                damping_coefficient_kg_m2_per_s=values[
                    "damping_coefficient_kg_m2_per_s"
                ],
                gravity_m_per_s2=values["gravity_m_per_s2"],
                initial_angle_rad=values["initial_angle_rad"],
                initial_angular_velocity_rad_per_s=values[
                    "initial_angular_velocity_rad_per_s"
                ],
                duration_s=values["duration_s"],
                num_points=output_sampling.point_count,
                method=solver_method,
                relative_tolerance=solver_configuration.relative_tolerance,
                absolute_tolerance=solver_configuration.absolute_tolerance,
                maximum_step_s=solver_configuration.maximum_step,
            )
        except ScientificValidationError as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_method,
            )
        except IntegrationError as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.SOLVER,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_method,
                solver_success=False,
            )

        try:
            extracted = _extract_pendulum_metrics(result, metrics)
        except (ScientificValidationError, ValueError, FloatingPointError) as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.METRIC_EXTRACTION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_method,
                solver_success=True,
            )
        return AdapterSuccess(
            metrics=extracted,
            solver_id=self.solver_id,
            solver_method=solver_method,
            solver_success=True,
            numerical_acceptance=None,
            function_evaluations=result.function_evaluations,
            ode_trajectory=(
                _pendulum_trajectory(result) if retain_trajectory else None
            ),
        )

    def analytical_reference(
        self, parameters: tuple[NamedValue, ...], metric_id: str
    ) -> AdapterAnalyticalReference | None:
        """Return the undamped small-angle linearized final angle, if applicable."""
        if metric_id != "final_angle_rad":
            return None
        values = {item.name: item.value for item in parameters}
        if values["damping_coefficient_kg_m2_per_s"] != 0.0:
            return None
        angular_frequency = np.sqrt(values["gravity_m_per_s2"] / values["length_m"])
        final_time = values["duration_s"]
        angle = values["initial_angle_rad"] * np.cos(
            angular_frequency * final_time
        ) + values["initial_angular_velocity_rad_per_s"] / angular_frequency * np.sin(
            angular_frequency * final_time
        )
        return AdapterAnalyticalReference(
            value=float(angle),
            unit="rad",
            description=(
                "Undamped small-angle linearized pendulum reference; it is an "
                "approximation, not an exact solution of the nonlinear model."
            ),
        )


class BVPExperimentAdapter:
    """Experiment adapter that preserves a spatial BVP candidate and evidence."""

    model_id = "boundary_value_problem"

    def __init__(
        self,
        problem: BoundaryValueProblem,
        configuration: BoundaryValueSolverConfiguration,
        *,
        parameter_constraints: tuple[BVPParameterConstraint, ...] = (),
        solver: ScipyBoundaryValueSolver | None = None,
    ) -> None:
        self.problem = problem
        self.configuration = configuration
        self.parameter_constraints = parameter_constraints
        self.solver = solver or ScipyBoundaryValueSolver()
        self.solver_id = self.solver.capabilities.solver_id

    def execute(
        self,
        parameters: tuple[NamedValue, ...],
        metrics: tuple[str, ...],
        solver_configuration: ExperimentSolverConfiguration,
        output_sampling: OutputSampling,
    ) -> AdapterOutcome:
        del output_sampling  # BVP mesh settings belong to the problem itself.
        if not isinstance(solver_configuration, BoundaryValueSolverConfiguration):
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type="ConfigurationTypeError",
                message="BVP experiments require BoundaryValueSolverConfiguration",
                solver_id=self.solver_id,
                solver_method="unknown",
            )
        parameter_names = tuple(item.name for item in parameters)
        expected_names = tuple(item.name for item in self.problem.parameters)
        if parameter_names != expected_names:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type="ScientificValidationError",
                message="BVP parameter names and order must match the problem",
                solver_id=self.solver_id,
                solver_method=solver_configuration.method,
            )
        values = {item.name: item.value for item in parameters}
        try:
            for constraint in self.parameter_constraints:
                value = values[constraint.parameter_name]
                valid = value > 0.0 if constraint.domain == "positive" else value >= 0.0
                if not valid:
                    relation = (
                        "positive"
                        if constraint.domain == "positive"
                        else "non-negative"
                    )
                    raise ScientificValidationError(
                        f"parameter {constraint.parameter_name!r} must be {relation}"
                    )
            run_problem = self.problem.model_copy(update={"parameters": parameters})
            self._validate_metric_ids(run_problem, metrics)
        except (ScientificValidationError, KeyError) as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_configuration.method,
            )

        outcome = self.solver.solve(run_problem, solver_configuration)
        if isinstance(outcome, BoundaryValueFailure):
            failure_kind = (
                ExperimentFailureKind.MODEL_VALIDATION
                if outcome.error_type == "ScientificValidationError"
                else ExperimentFailureKind.SOLVER
            )
            return AdapterFailure(
                failure_kind=failure_kind,
                error_type=outcome.error_type,
                message=outcome.message,
                solver_id=outcome.solver_id,
                solver_method=solver_configuration.method,
                solver_success=(
                    False if failure_kind == ExperimentFailureKind.SOLVER else None
                ),
            )

        diagnostics: BVPRunDiagnostics | None = None
        try:
            diagnostics = _bvp_diagnostics(outcome)
            extracted = _extract_bvp_metrics(outcome, metrics)
        except (ScientificValidationError, ValueError, FloatingPointError) as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.METRIC_EXTRACTION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_configuration.method,
                solver_success=outcome.solver_success,
                numerical_acceptance=outcome.numerical_accepted,
                metrics=(),
                bvp_diagnostics=diagnostics,
            )
        if not outcome.solver_success:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.SOLVER,
                error_type="BVPNonconvergence",
                message=outcome.message,
                solver_id=self.solver_id,
                solver_method=solver_configuration.method,
                solver_success=False,
                numerical_acceptance=outcome.numerical_accepted,
                metrics=extracted,
                bvp_diagnostics=diagnostics,
            )
        if not outcome.numerical_accepted:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.RESIDUAL_ACCEPTANCE,
                error_type="BVPNumericalAcceptanceFailure",
                message=(
                    "BVP candidate did not satisfy the solver's configured "
                    "boundary and differential residual acceptance criteria"
                ),
                solver_id=self.solver_id,
                solver_method=solver_configuration.method,
                solver_success=True,
                numerical_acceptance=False,
                metrics=extracted,
                bvp_diagnostics=diagnostics,
            )
        return AdapterSuccess(
            metrics=extracted,
            solver_id=self.solver_id,
            solver_method=solver_configuration.method,
            solver_success=True,
            numerical_acceptance=outcome.numerical_accepted,
            bvp_diagnostics=diagnostics,
        )

    @staticmethod
    def _validate_metric_ids(
        problem: BoundaryValueProblem, metrics: tuple[str, ...]
    ) -> None:
        supported = {"domain_start", "domain_end", "mesh_point_count"}
        supported.update(
            {
                f"{prefix}:{state.name}"
                for prefix in (
                    "sampled_minimum",
                    "sampled_minimum_location",
                    "sampled_maximum",
                    "sampled_maximum_location",
                    "left_endpoint",
                    "right_endpoint",
                )
                for state in problem.states
            }
        )
        supported.update(
            {
                "boundary_residual_max_abs",
                "boundary_residual_l2_norm",
                "boundary_residual_component_0",
                "boundary_residual_component_1",
            }
        )
        unsupported = set(metrics) - supported
        if unsupported:
            raise ScientificValidationError(
                f"unsupported BVP metrics: {sorted(unsupported)}"
            )

    def analytical_reference(
        self, parameters: tuple[NamedValue, ...], metric_id: str
    ) -> AdapterAnalyticalReference | None:
        del parameters, metric_id
        return None


class AlgebraicExperimentAdapter:
    """Adapter around the existing local nonlinear algebraic root solver."""

    model_id = "nonlinear_algebraic_system"

    def __init__(
        self,
        problem: NonlinearAlgebraicProblem,
        configuration: AlgebraicSolverConfiguration,
        *,
        reference_solution: tuple[float, ...] | None = None,
        parameter_constraints: tuple[AlgebraicParameterConstraint, ...] = (),
        solver: AlgebraicSolver | None = None,
    ) -> None:
        self.problem = problem
        self.configuration = configuration
        self.reference_solution = reference_solution
        self.parameter_constraints = parameter_constraints
        self.solver = solver or ScipyHybridRootSolver()
        self.solver_id = self.solver.capabilities.solver_id

    def execute(
        self,
        parameters: tuple[NamedValue, ...],
        metrics: tuple[str, ...],
        solver_configuration: ExperimentSolverConfiguration,
        output_sampling: OutputSampling,
    ) -> AdapterOutcome:
        del output_sampling  # Algebraic solves have no sampled independent variable.
        if not isinstance(solver_configuration, AlgebraicSolverConfiguration):
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type="ConfigurationTypeError",
                message="algebraic experiments require AlgebraicSolverConfiguration",
                solver_id=self.solver_id,
                solver_method="unknown",
            )
        parameter_names = tuple(item.name for item in parameters)
        expected_names = tuple(item.name for item in self.problem.parameters)
        solver_method = solver_configuration.method
        if parameter_names != expected_names:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type="ScientificValidationError",
                message=(
                    "algebraic parameter names and order must match the configured "
                    "problem"
                ),
                solver_id=self.solver_id,
                solver_method=solver_method,
            )
        try:
            run_problem = self._problem_for_parameters(parameters)
            self._validate_metric_ids(run_problem, metrics)
        except ScientificValidationError as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.MODEL_VALIDATION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=self.solver_id,
                solver_method=solver_method,
            )

        outcome = self.solver.solve(run_problem, solver_configuration)
        if isinstance(outcome, AlgebraicFailure):
            diagnostics = AlgebraicRunDiagnostics(
                problem_id=run_problem.problem_id,
                initial_guess=tuple(float(item) for item in run_problem.initial_guess),
                variable_names=tuple(item.name for item in run_problem.unknowns),
                variable_units=tuple(item.unit for item in run_problem.unknowns),
                solution_candidate=outcome.last_candidate,
                residual_vector=outcome.residual_vector,
                residual_names=run_problem.residual_names,
                residual_units=run_problem.residual_units,
                residual_norm=(
                    None
                    if outcome.residual_vector is None
                    else float(np.linalg.norm(outcome.residual_vector, ord=2))
                ),
                residual_norm_definition=(
                    "Euclidean L2 norm of raw residual entries; no scaling is applied"
                ),
                residual_tolerance=solver_configuration.residual_tolerance,
                termination_reason=outcome.message,
                function_evaluations=outcome.function_evaluations,
            )
            kind = (
                ExperimentFailureKind.MODEL_VALIDATION
                if outcome.error_type == "ScientificValidationError"
                else ExperimentFailureKind.SOLVER
            )
            return AdapterFailure(
                failure_kind=kind,
                error_type=outcome.error_type,
                message=outcome.message,
                solver_id=outcome.solver_id,
                solver_method=solver_method,
                solver_success=False if kind == ExperimentFailureKind.SOLVER else None,
                function_evaluations=outcome.function_evaluations,
                algebraic_diagnostics=diagnostics,
            )

        diagnostics = self._diagnostics(run_problem, outcome)
        try:
            extracted = _extract_algebraic_metrics(
                run_problem,
                tuple(float(item) for item in outcome.solution_candidate),
                tuple(float(item) for item in outcome.residual_vector),
                outcome.residual_norm,
                metrics,
            )
        except (ScientificValidationError, ValueError, FloatingPointError) as exc:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.METRIC_EXTRACTION,
                error_type=type(exc).__name__,
                message=str(exc),
                solver_id=outcome.solver_id,
                solver_method=solver_method,
                solver_success=outcome.solver_converged,
                numerical_acceptance=outcome.residual_accepted,
                function_evaluations=outcome.function_evaluations,
                metrics=(),
                algebraic_diagnostics=diagnostics,
            )

        if not outcome.solver_converged:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.SOLVER,
                error_type="AlgebraicNonconvergence",
                message=outcome.termination_reason,
                solver_id=outcome.solver_id,
                solver_method=solver_method,
                solver_success=False,
                numerical_acceptance=outcome.residual_accepted,
                function_evaluations=outcome.function_evaluations,
                metrics=extracted,
                algebraic_diagnostics=diagnostics,
            )
        if not outcome.residual_accepted:
            return AdapterFailure(
                failure_kind=ExperimentFailureKind.RESIDUAL_ACCEPTANCE,
                error_type="ResidualAcceptanceFailure",
                message=(
                    f"residual norm {outcome.residual_norm:g} exceeds configured "
                    f"tolerance {solver_configuration.residual_tolerance:g}"
                ),
                solver_id=outcome.solver_id,
                solver_method=solver_method,
                solver_success=True,
                numerical_acceptance=False,
                function_evaluations=outcome.function_evaluations,
                metrics=extracted,
                algebraic_diagnostics=diagnostics,
            )
        return AdapterSuccess(
            metrics=extracted,
            solver_id=outcome.solver_id,
            solver_method=solver_method,
            solver_success=outcome.solver_converged,
            numerical_acceptance=outcome.residual_accepted,
            function_evaluations=outcome.function_evaluations,
            algebraic_diagnostics=diagnostics,
        )

    def _problem_for_parameters(
        self, parameters: tuple[NamedValue, ...]
    ) -> NonlinearAlgebraicProblem:
        values = {item.name: item for item in parameters}
        for constraint in self.parameter_constraints:
            value = values[constraint.parameter_name].value
            valid = value > 0.0 if constraint.domain == "positive" else value >= 0.0
            if not valid:
                relation = (
                    "positive" if constraint.domain == "positive" else "non-negative"
                )
                raise ScientificValidationError(
                    f"parameter {constraint.parameter_name!r} must be {relation}"
                )
        updated = tuple(
            AlgebraicParameter(
                name=parameter.name,
                value=values[parameter.name].value,
                unit=parameter.unit,
            )
            for parameter in self.problem.parameters
        )
        return self.problem.model_copy(update={"parameters": updated})

    @staticmethod
    def _validate_metric_ids(
        problem: NonlinearAlgebraicProblem, metrics: tuple[str, ...]
    ) -> None:
        supported = {"residual_norm"}
        supported.update(f"solution:{item.name}" for item in problem.unknowns)
        residual_names = problem.residual_names or tuple(
            f"equation_{index}" for index in range(len(problem.unknowns))
        )
        supported.update(f"residual:{name}" for name in residual_names)
        unsupported = set(metrics) - supported
        if unsupported:
            raise ScientificValidationError(
                f"unsupported algebraic metrics: {sorted(unsupported)}"
            )

    def _diagnostics(
        self, problem: NonlinearAlgebraicProblem, result: AlgebraicResult
    ) -> AlgebraicRunDiagnostics:
        return AlgebraicRunDiagnostics(
            problem_id=problem.problem_id,
            initial_guess=tuple(float(item) for item in problem.initial_guess),
            variable_names=result.variable_names,
            variable_units=result.variable_units,
            solution_candidate=tuple(float(item) for item in result.solution_candidate),
            residual_vector=tuple(float(item) for item in result.residual_vector),
            residual_names=result.residual_names,
            residual_units=result.residual_units,
            residual_norm=result.residual_norm,
            residual_norm_definition=result.residual_norm_definition,
            solver_converged=result.solver_converged,
            residual_accepted=result.residual_accepted,
            residual_tolerance=result.configuration.residual_tolerance,
            termination_reason=result.termination_reason,
            iteration_count=result.iteration_count,
            function_evaluations=result.function_evaluations,
        )

    def analytical_reference(
        self, parameters: tuple[NamedValue, ...], metric_id: str
    ) -> AdapterAnalyticalReference | None:
        """Return a configured candidate only when it satisfies this run's residual."""
        if self.reference_solution is None:
            return None
        problem = self._problem_for_parameters(parameters)
        self._validate_metric_ids(problem, (metric_id,))
        if not metric_id.startswith("solution:"):
            return None
        variable_name = metric_id.removeprefix("solution:")
        index = next(
            (
                i
                for i, unknown in enumerate(problem.unknowns)
                if unknown.name == variable_name
            ),
            None,
        )
        if index is None:
            return None
        residual = problem.evaluate_residual(np.asarray(self.reference_solution))
        if (
            float(np.linalg.norm(residual, ord=2))
            > self.configuration.residual_tolerance
        ):
            return None
        unit = problem.unknowns[index].unit or "unknown"
        return AdapterAnalyticalReference(
            value=self.reference_solution[index],
            unit=unit,
            description=(
                "Caller-supplied independently checkable algebraic root; this "
                "comparison does not establish uniqueness or completeness."
            ),
        )


def _extract_oscillator_metrics(
    result: DampedOscillatorResult, metric_ids: tuple[str, ...]
) -> tuple[MetricValue, ...]:
    """Extract explicitly defined metrics from requested output samples."""
    metric_values: dict[str, float] = {
        "maximum_absolute_displacement_m": float(np.max(np.abs(result.displacement_m))),
        "final_displacement_m": float(result.displacement_m[-1]),
        "final_velocity_m_per_s": float(result.velocity_m_per_s[-1]),
        "maximum_absolute_velocity_m_per_s": float(
            np.max(np.abs(result.velocity_m_per_s))
        ),
    }
    units = OscillatorExperimentAdapter._metrics
    try:
        return tuple(
            MetricValue(
                metric_id=metric_id,
                value=metric_values[metric_id],
                unit=units[metric_id],
            )
            for metric_id in metric_ids
        )
    except (KeyError, ValueError) as exc:
        raise ScientificValidationError(
            f"could not extract requested metric: {exc}"
        ) from exc


def _oscillator_trajectory(result: DampedOscillatorResult) -> ODETrajectorySnapshot:
    """Copy the requested oscillator output grid into an immutable record."""
    return ODETrajectorySnapshot(
        model_id=OscillatorExperimentAdapter.model_id,
        time_s=tuple(float(value) for value in result.time_s),
        state_names=("displacement", "velocity"),
        state_units=("m", "m/s"),
        state_values=(
            tuple(float(value) for value in result.displacement_m),
            tuple(float(value) for value in result.velocity_m_per_s),
        ),
        solver_method=result.method,
    )


def _pendulum_trajectory(result: DampedPendulumResult) -> ODETrajectorySnapshot:
    """Copy the requested pendulum output grid into an immutable record."""
    return ODETrajectorySnapshot(
        model_id=PendulumExperimentAdapter.model_id,
        time_s=tuple(float(value) for value in result.time_s),
        state_names=("angle", "angular_velocity"),
        state_units=("rad", "rad/s"),
        state_values=(
            tuple(float(value) for value in result.angle_rad),
            tuple(float(value) for value in result.angular_velocity_rad_per_s),
        ),
        solver_method=result.method,
    )


def _extract_pendulum_metrics(
    result: DampedPendulumResult, metric_ids: tuple[str, ...]
) -> tuple[MetricValue, ...]:
    """Extract pendulum metrics from the inclusive requested output grid."""
    metric_values: dict[str, float] = {
        "final_angle_rad": float(result.angle_rad[-1]),
        "final_angular_velocity_rad_per_s": float(
            result.angular_velocity_rad_per_s[-1]
        ),
        "maximum_absolute_angle_rad": float(np.max(np.abs(result.angle_rad))),
        "maximum_absolute_angular_velocity_rad_per_s": float(
            np.max(np.abs(result.angular_velocity_rad_per_s))
        ),
    }
    try:
        return tuple(
            MetricValue(
                metric_id=metric_id,
                value=metric_values[metric_id],
                unit=PendulumExperimentAdapter._metrics[metric_id],
            )
            for metric_id in metric_ids
        )
    except (KeyError, ValueError) as exc:
        raise ScientificValidationError(
            f"could not extract requested pendulum metric: {exc}"
        ) from exc


def _bvp_diagnostics(result: BoundaryValueResult) -> BVPRunDiagnostics:
    """Copy BVP grids and diagnostics into a serializable experiment record."""
    first_state, second_state = result.problem.states
    return BVPRunDiagnostics(
        problem_id=result.problem.problem_id,
        domain=result.problem.domain,
        domain_unit=result.problem.domain_unit,
        mesh=tuple(float(value) for value in result.mesh),
        state_names=(first_state.name, second_state.name),
        state_units=(first_state.unit, second_state.unit),
        solution=(
            tuple(float(value) for value in result.solution[0]),
            tuple(float(value) for value in result.solution[1]),
        ),
        boundary_residuals=(
            float(result.boundary_residuals[0]),
            float(result.boundary_residuals[1]),
        ),
        boundary_residual_units=result.problem.boundary_residual_units,
        differential_residuals=tuple(
            float(value) for value in result.differential_residuals
        ),
        solver_success=result.solver_success,
        boundary_conditions_satisfied=result.boundary_conditions_satisfied,
        numerical_accepted=result.numerical_accepted,
        status_code=result.status_code,
        message=result.message,
        iteration_count=result.iteration_count,
        solver_tolerance=result.configuration.tolerance,
        boundary_tolerance=result.configuration.boundary_tolerance,
        acceptance_tolerance=result.configuration.acceptance_tolerance,
    )


def _extract_bvp_metrics(
    result: BoundaryValueResult, metric_ids: tuple[str, ...]
) -> tuple[MetricValue, ...]:
    """Extract sampled metrics from aligned spatial nodes, not continuous extrema."""
    mesh = np.asarray(result.mesh, dtype=np.float64)
    profile = np.asarray(result.solution, dtype=np.float64)
    if (
        mesh.ndim != 1
        or len(mesh) < 2
        or profile.shape != (2, len(mesh))
        or not np.isfinite(mesh).all()
        or not np.isfinite(profile).all()
        or not np.all(np.diff(mesh) > 0.0)
    ):
        raise ScientificValidationError(
            "BVP result must contain finite aligned profiles on an increasing mesh"
        )
    BVPExperimentAdapter._validate_metric_ids(result.problem, metric_ids)
    boundary_units = result.problem.boundary_residual_units
    distinct_boundary_units = {unit for unit in boundary_units}
    common_boundary_unit = "unknown"
    if len(distinct_boundary_units) == 1 and None not in distinct_boundary_units:
        common_boundary_unit = boundary_units[0] or "unknown"
    values: dict[str, tuple[float, str]] = {
        "domain_start": (float(mesh[0]), result.problem.domain_unit),
        "domain_end": (float(mesh[-1]), result.problem.domain_unit),
        "mesh_point_count": (float(len(mesh)), "count"),
        "boundary_residual_max_abs": (
            float(np.max(np.abs(result.boundary_residuals))),
            common_boundary_unit,
        ),
        "boundary_residual_l2_norm": (
            float(np.linalg.norm(result.boundary_residuals, ord=2)),
            common_boundary_unit,
        ),
        "boundary_residual_component_0": (
            float(result.boundary_residuals[0]),
            boundary_units[0] or "unknown",
        ),
        "boundary_residual_component_1": (
            float(result.boundary_residuals[1]),
            boundary_units[1] or "unknown",
        ),
    }
    state_index = {
        state.name: index for index, state in enumerate(result.problem.states)
    }
    for metric_id in metric_ids:
        prefix, separator, state_name = metric_id.partition(":")
        if not separator:
            continue
        index = state_index[state_name]
        component = profile[index]
        unit = result.problem.states[index].unit or "unknown"
        if prefix == "sampled_minimum":
            values[metric_id] = (float(np.min(component)), unit)
        elif prefix == "sampled_maximum":
            values[metric_id] = (float(np.max(component)), unit)
        elif prefix == "sampled_minimum_location":
            values[metric_id] = (
                float(mesh[int(np.argmin(component))]),
                result.problem.domain_unit,
            )
        elif prefix == "sampled_maximum_location":
            values[metric_id] = (
                float(mesh[int(np.argmax(component))]),
                result.problem.domain_unit,
            )
        elif prefix == "left_endpoint":
            values[metric_id] = (float(component[0]), unit)
        elif prefix == "right_endpoint":
            values[metric_id] = (float(component[-1]), unit)
    try:
        return tuple(
            MetricValue(
                metric_id=metric_id,
                value=values[metric_id][0],
                unit=values[metric_id][1],
            )
            for metric_id in metric_ids
        )
    except (KeyError, ValueError) as exc:
        raise ScientificValidationError(
            f"could not extract requested BVP metric: {exc}"
        ) from exc


def _extract_algebraic_metrics(
    problem: NonlinearAlgebraicProblem,
    solution: tuple[float, ...],
    residual: tuple[float, ...],
    residual_norm: float,
    metric_ids: tuple[str, ...],
) -> tuple[MetricValue, ...]:
    """Extract named candidate/residual scalars and the raw residual L2 norm."""
    solution_values = {
        f"solution:{unknown.name}": (
            value,
            unknown.unit or "unknown",
        )
        for unknown, value in zip(problem.unknowns, solution, strict=True)
    }
    residual_names = problem.residual_names or tuple(
        f"equation_{index}" for index in range(len(problem.unknowns))
    )
    residual_units = problem.residual_units or tuple(
        None for _ in range(len(problem.unknowns))
    )
    residual_values = {
        f"residual:{name}": (value, unit or "unknown")
        for name, unit, value in zip(
            residual_names, residual_units, residual, strict=True
        )
    }
    declared_residual_units = {unit for unit in residual_units if unit is not None}
    norm_unit = (
        next(iter(declared_residual_units))
        if len(declared_residual_units) == 1
        else "unknown"
    )
    values = {
        **solution_values,
        **residual_values,
        "residual_norm": (residual_norm, norm_unit),
    }
    try:
        return tuple(
            MetricValue(
                metric_id=metric_id,
                value=values[metric_id][0],
                unit=values[metric_id][1],
            )
            for metric_id in metric_ids
        )
    except (KeyError, ValueError) as exc:
        raise ScientificValidationError(
            f"could not extract requested algebraic metric: {exc}"
        ) from exc


def _replace_parameter(
    parameters: tuple[NamedValue, ...], name: str, value: float
) -> tuple[NamedValue, ...]:
    """Return a fresh ordered parameter tuple with one substituted value."""
    return tuple(
        item.model_copy(update={"value": value}) if item.name == name else item
        for item in parameters
    )


def _scipy_version() -> str:
    try:
        return version("scipy")
    except PackageNotFoundError:
        return "unknown"


def _reproducibility(
    specification: ExperimentSpecification,
    parameters: tuple[NamedValue, ...],
    solver_id: str,
) -> ExperimentReproducibility:
    """Create run-specific configuration and environment metadata."""
    configuration = specification.solver_configuration
    if isinstance(configuration, ODESolverConfiguration):
        problem_id = None
        initial_guess = None
        step_tolerance = None
        residual_tolerance = None
        maximum_function_evaluations = None
        output_point_count = specification.output_sampling.point_count
        relative_tolerance = configuration.relative_tolerance
        absolute_tolerance = configuration.absolute_tolerance
        maximum_step = configuration.maximum_step
        bvp_metadata: dict[str, object] = {}
    elif isinstance(configuration, AlgebraicSolverConfiguration):
        problem = specification.algebraic_problem
        problem_id = None if problem is None else problem.problem_id
        initial_guess = (
            None
            if problem is None
            else tuple(float(item) for item in problem.initial_guess)
        )
        step_tolerance = configuration.step_tolerance
        residual_tolerance = configuration.residual_tolerance
        maximum_function_evaluations = configuration.maximum_function_evaluations
        output_point_count = None
        relative_tolerance = None
        absolute_tolerance = None
        maximum_step = None
        bvp_metadata = {}
    else:
        bvp_problem = specification.bvp_problem
        problem_id = None if bvp_problem is None else bvp_problem.problem_id
        initial_guess = None
        step_tolerance = None
        residual_tolerance = None
        maximum_function_evaluations = None
        output_point_count = None
        relative_tolerance = None
        absolute_tolerance = None
        maximum_step = None
        bvp_metadata = (
            {}
            if bvp_problem is None
            else {
                "bvp_domain": bvp_problem.domain,
                "bvp_mesh": bvp_problem.mesh,
                "bvp_initial_guess": bvp_problem.initial_guess,
                "bvp_solver_tolerance": configuration.tolerance,
                "bvp_boundary_tolerance": configuration.boundary_tolerance,
                "bvp_acceptance_tolerance": configuration.acceptance_tolerance,
                "bvp_maximum_nodes": configuration.maximum_nodes,
            }
        )
    return ExperimentReproducibility(
        experiment_id=specification.experiment_id,
        configuration_version=specification.configuration_version,
        model_id=specification.model_id,
        parameters=parameters,
        assumptions=specification.assumptions,
        equation_record_ids=specification.equation_record_ids,
        source_reference_ids=specification.source_reference_ids,
        requested_metrics=specification.metrics,
        solver_id=solver_id,
        solver_method=configuration.method,
        problem_id=problem_id,
        initial_guess=initial_guess,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        maximum_step=maximum_step,
        step_tolerance=step_tolerance,
        residual_tolerance=residual_tolerance,
        maximum_function_evaluations=maximum_function_evaluations,
        output_point_count=output_point_count,
        **bvp_metadata,
        python_version=sys.version.split()[0],
        numpy_version=np.__version__,
        scipy_version=_scipy_version(),
        random_seed=None,
    )


def _compare_reference(
    specification: ExperimentSpecification,
    parameters: tuple[NamedValue, ...],
    metrics: tuple[MetricValue, ...],
    adapter: ExperimentAdapter,
) -> AnalyticalReferenceComparison | None:
    request = specification.analytical_reference
    if request is None:
        return None
    value = adapter.analytical_reference(parameters, request.metric_id)
    if value is None:
        return AnalyticalReferenceComparison(
            status="unavailable",
            metric_id=request.metric_id,
            unit="unknown",
            absolute_tolerance=request.absolute_tolerance,
            message=(
                "No analytical reference is defined for this model point and metric."
            ),
        )
    reference_value = value.value
    reference_unit = value.unit
    numerical_metric = next(
        item for item in metrics if item.metric_id == request.metric_id
    )
    if reference_unit != numerical_metric.unit:
        return AnalyticalReferenceComparison(
            status="failed",
            metric_id=request.metric_id,
            unit=numerical_metric.unit,
            numerical_value=numerical_metric.value,
            reference_value=reference_value,
            absolute_tolerance=request.absolute_tolerance,
            reference_description=value.description,
            message=(
                f"Reference unit {reference_unit!r} does not match metric unit "
                f"{numerical_metric.unit!r}."
            ),
        )
    error = abs(numerical_metric.value - reference_value)
    passed = error <= request.absolute_tolerance
    return AnalyticalReferenceComparison(
        status="passed" if passed else "failed",
        metric_id=request.metric_id,
        unit=numerical_metric.unit,
        numerical_value=numerical_metric.value,
        reference_value=reference_value,
        absolute_error=error,
        absolute_tolerance=request.absolute_tolerance,
        reference_description=value.description,
        message=(
            f"{value.description} "
            + (
                "Within absolute tolerance."
                if passed
                else "Absolute error exceeds tolerance."
            )
        ),
    )


def _compare_bvp_reference(
    specification: ExperimentSpecification,
    parameters: tuple[NamedValue, ...],
    diagnostics: BVPRunDiagnostics,
) -> BVPReferenceComparison | None:
    """Compare the heat example with its exact linear profile on returned nodes."""
    request = specification.bvp_reference
    if request is None:
        return None
    problem = specification.bvp_problem
    values = {item.name: item.value for item in parameters}
    is_heat_reference = (
        problem is not None
        and problem.problem_id == "steady_linear_heat_conduction"
        and diagnostics.state_names[0] == "temperature"
        and diagnostics.state_units[0] == "K"
        and "left_temperature_k" in values
        and "right_temperature_k" in values
        and diagnostics.domain_unit == "m"
    )
    if not is_heat_reference:
        return BVPReferenceComparison(
            status="unavailable",
            method="direct_evaluation_on_returned_mesh",
            state_name="temperature",
            unit="K",
            maximum_absolute_tolerance=request.maximum_absolute_tolerance,
            reference_description=(
                "Exact linear heat profile for steady one-dimensional constant-"
                "conductivity conduction with no internal heat generation."
            ),
            message="No independently justified reference is available for this BVP.",
        )
    x = np.asarray(diagnostics.mesh, dtype=np.float64)
    left, right = diagnostics.domain
    expected = values["left_temperature_k"] + (
        values["right_temperature_k"] - values["left_temperature_k"]
    ) * (x - left) / (right - left)
    actual = np.asarray(diagnostics.solution[0], dtype=np.float64)
    errors = np.abs(actual - expected)
    maximum_error = float(np.max(errors))
    rms_error = float(np.sqrt(np.mean(np.square(errors))))
    passed = maximum_error <= request.maximum_absolute_tolerance
    return BVPReferenceComparison(
        status="passed" if passed else "failed",
        method="direct_evaluation_on_returned_mesh",
        state_name="temperature",
        unit="K",
        maximum_absolute_error=maximum_error,
        rms_error=rms_error,
        maximum_absolute_tolerance=request.maximum_absolute_tolerance,
        reference_description=(
            "Exact linear profile for constant-conductivity steady 1D heat "
            "conduction with no internal generation and fixed endpoint temperatures."
        ),
        message=(
            "Profile error is within tolerance."
            if passed
            else "Maximum sampled profile error exceeds tolerance."
        ),
    )


def _run_point(
    specification: ExperimentSpecification,
    parameters: tuple[NamedValue, ...],
    adapter: ExperimentAdapter,
    *,
    run_index: int,
    swept_parameter_name: str | None = None,
    swept_parameter_value: float | None = None,
) -> ExperimentRun:
    if specification.retain_ode_trajectory:
        if not isinstance(
            adapter, (OscillatorExperimentAdapter, PendulumExperimentAdapter)
        ):
            raise ScientificValidationError(
                "trajectory retention requires a built-in ODE experiment adapter"
            )
        outcome = adapter.execute(
            parameters,
            specification.metrics,
            specification.solver_configuration,
            specification.output_sampling,
            retain_trajectory=True,
        )
    else:
        outcome = adapter.execute(
            parameters,
            specification.metrics,
            specification.solver_configuration,
            specification.output_sampling,
        )
    reproducibility = _reproducibility(specification, parameters, outcome.solver_id)
    if isinstance(outcome, AdapterFailure):
        bvp_reference = (
            None
            if outcome.bvp_diagnostics is None
            else _compare_bvp_reference(
                specification, parameters, outcome.bvp_diagnostics
            )
        )
        return ExperimentRun(
            run_index=run_index,
            swept_parameter_name=swept_parameter_name,
            swept_parameter_value=swept_parameter_value,
            parameters=parameters,
            status=ExperimentRunStatus.FAILED,
            failure_kind=outcome.failure_kind,
            error_type=outcome.error_type,
            message=outcome.message,
            solver_id=outcome.solver_id,
            solver_method=outcome.solver_method,
            solver_success=outcome.solver_success,
            numerical_acceptance=outcome.numerical_acceptance,
            function_evaluations=outcome.function_evaluations,
            metrics=outcome.metrics,
            algebraic_diagnostics=outcome.algebraic_diagnostics,
            bvp_diagnostics=outcome.bvp_diagnostics,
            bvp_reference_comparison=bvp_reference,
            ode_trajectory=None,
            reproducibility=reproducibility,
        )
    reference = _compare_reference(specification, parameters, outcome.metrics, adapter)
    bvp_reference = (
        None
        if outcome.bvp_diagnostics is None
        else _compare_bvp_reference(specification, parameters, outcome.bvp_diagnostics)
    )
    return ExperimentRun(
        run_index=run_index,
        swept_parameter_name=swept_parameter_name,
        swept_parameter_value=swept_parameter_value,
        parameters=parameters,
        status=ExperimentRunStatus.SUCCEEDED,
        metrics=outcome.metrics,
        solver_id=outcome.solver_id,
        solver_method=outcome.solver_method,
        solver_success=outcome.solver_success,
        numerical_acceptance=outcome.numerical_acceptance,
        function_evaluations=outcome.function_evaluations,
        reference_comparison=reference,
        algebraic_diagnostics=outcome.algebraic_diagnostics,
        bvp_diagnostics=outcome.bvp_diagnostics,
        bvp_reference_comparison=bvp_reference,
        ode_trajectory=outcome.ode_trajectory,
        reproducibility=reproducibility,
    )


def _select_adapter(
    specification: ExperimentSpecification, adapter: ExperimentAdapter | None
) -> ExperimentAdapter:
    """Choose one of the explicitly supported model adapters."""
    if adapter is not None:
        selected = adapter
    elif specification.model_id == OscillatorExperimentAdapter.model_id:
        selected = OscillatorExperimentAdapter()
    elif specification.model_id == PendulumExperimentAdapter.model_id:
        selected = PendulumExperimentAdapter()
    elif (
        specification.model_id == BVPExperimentAdapter.model_id
        and specification.bvp_problem is not None
        and isinstance(
            specification.solver_configuration, BoundaryValueSolverConfiguration
        )
    ):
        selected = BVPExperimentAdapter(
            specification.bvp_problem,
            specification.solver_configuration,
            parameter_constraints=specification.bvp_parameter_constraints,
        )
    elif (
        specification.model_id == AlgebraicExperimentAdapter.model_id
        and specification.algebraic_problem is not None
        and isinstance(specification.solver_configuration, AlgebraicSolverConfiguration)
    ):
        selected = AlgebraicExperimentAdapter(
            specification.algebraic_problem,
            specification.solver_configuration,
            reference_solution=specification.algebraic_reference_solution,
            parameter_constraints=specification.algebraic_parameter_constraints,
        )
    else:
        raise ScientificValidationError(
            f"no experiment adapter supports {specification.model_id!r}"
        )
    if specification.model_id != selected.model_id:
        raise ScientificValidationError(
            f"adapter {selected.solver_id!r} does not support model "
            f"{specification.model_id!r}"
        )
    return selected


def run_experiment(
    specification: ExperimentSpecification,
    adapter: ExperimentAdapter | None = None,
) -> ExperimentRecord:
    """Execute an ordered sweep, recording expected failures and continuing.

    The default adapters support the oscillator, pendulum, algebraic, and BVP
    models. An expected model-validation, solver, or metric failure is recorded
    for that point; subsequent points still run. Unexpected programming
    exceptions propagate.
    """
    selected_adapter = _select_adapter(specification, adapter)
    points: tuple[tuple[str | None, float | None, tuple[NamedValue, ...]], ...]
    if specification.sweep is None:
        points = ((None, None, specification.parameters),)
    else:
        sweep = specification.sweep
        points = tuple(
            (
                sweep.parameter_name,
                value,
                _replace_parameter(
                    specification.parameters, sweep.parameter_name, value
                ),
            )
            for value in sweep.values
        )
    runs = tuple(
        _run_point(
            specification,
            parameters,
            selected_adapter,
            run_index=index,
            swept_parameter_name=parameter_name,
            swept_parameter_value=value,
        )
        for index, (parameter_name, value, parameters) in enumerate(points)
    )
    status = (
        ExperimentRecordStatus.SUCCEEDED
        if all(run.status == ExperimentRunStatus.SUCCEEDED for run in runs)
        else ExperimentRecordStatus.COMPLETED_WITH_FAILURES
    )
    return ExperimentRecord(specification=specification, runs=runs, status=status)


class FiniteDifferenceEstimate(SimulationModel):
    """Central finite-difference estimate of an absolute scalar derivative."""

    method: Literal["central_difference"] = "central_difference"
    base_value: float
    perturbation: float
    lower_parameter_value: float
    upper_parameter_value: float
    lower_metric_value: float
    upper_metric_value: float
    derivative_estimate: float

    @field_validator(
        "base_value",
        "perturbation",
        "lower_parameter_value",
        "upper_parameter_value",
        "lower_metric_value",
        "upper_metric_value",
        "derivative_estimate",
    )
    @classmethod
    def finite_estimate_values(cls, value: float) -> float:
        if not np.isfinite(value):
            raise ValueError("finite-difference values must be finite")
        return value


def estimate_central_difference(
    function: Callable[[float], float], base_value: float, perturbation: float
) -> FiniteDifferenceEstimate:
    """Estimate ``df/dp`` using an explicit symmetric perturbation.

    The evaluator's errors propagate; no one-sided replacement is attempted.
    This is an absolute derivative, with units of output divided by input.
    """
    if not np.isfinite(base_value):
        raise ScientificValidationError("base_value must be finite")
    if not np.isfinite(perturbation) or perturbation <= 0.0:
        raise ScientificValidationError("perturbation must be finite and positive")
    lower = base_value - perturbation
    upper = base_value + perturbation
    lower_metric = function(lower)
    upper_metric = function(upper)
    if not np.isfinite(lower_metric) or not np.isfinite(upper_metric):
        raise ScientificValidationError("function evaluations must be finite")
    estimate = (upper_metric - lower_metric) / (2.0 * perturbation)
    return FiniteDifferenceEstimate(
        base_value=base_value,
        perturbation=perturbation,
        lower_parameter_value=lower,
        upper_parameter_value=upper,
        lower_metric_value=lower_metric,
        upper_metric_value=upper_metric,
        derivative_estimate=estimate,
    )


class SensitivityResult(SimulationModel):
    """Successful central-difference experiment with recorded evaluation runs."""

    experiment_id: str
    parameter_name: str
    parameter_unit: str
    metric_id: str
    metric_unit: str
    base_value: float
    perturbation: float
    lower_parameter_value: float
    upper_parameter_value: float
    base_metric_value: float
    lower_metric_value: float
    upper_metric_value: float
    derivative_estimate: float
    derivative_unit: str
    method: Literal["central_difference"] = "central_difference"
    estimate_kind: Literal[
        "metric_sensitivity",
        "candidate_difference",
        "mesh_dependent_metric_difference",
    ] = "metric_sensitivity"
    candidate_branch_consistency: Literal["not_applicable", "not_established"] = (
        "not_applicable"
    )
    runs: tuple[ExperimentRun, ExperimentRun, ExperimentRun]

    @field_validator(
        "base_value",
        "perturbation",
        "lower_parameter_value",
        "upper_parameter_value",
        "base_metric_value",
        "lower_metric_value",
        "upper_metric_value",
        "derivative_estimate",
    )
    @classmethod
    def finite_sensitivity_values(cls, value: float) -> float:
        if not np.isfinite(value):
            raise ValueError("sensitivity values must be finite")
        return value


class SensitivityFailure(SimulationModel):
    """Central difference that could not be completed, with all run outcomes."""

    experiment_id: str
    parameter_name: str
    metric_id: str
    perturbation: float
    runs: tuple[ExperimentRun, ...]
    failed_run_indices: tuple[int, ...]
    message: str


SensitivityOutcome = SensitivityResult | SensitivityFailure


def run_sensitivity_analysis(
    specification: ExperimentSpecification,
    *,
    parameter_name: str,
    metric_id: str,
    perturbation: float,
    adapter: ExperimentAdapter | None = None,
) -> SensitivityOutcome:
    """Run base, lower, and upper experiment points for central sensitivity.

    The specification must not also define a sweep. Invalid perturbed values
    become explicit model-validation run failures. All three evaluations are
    attempted in base/lower/upper order; no one-sided estimate is substituted.
    """
    if not np.isfinite(perturbation) or perturbation <= 0.0:
        raise ScientificValidationError("perturbation must be finite and positive")
    if specification.sweep is not None:
        raise ScientificValidationError(
            "sensitivity analysis requires a specification without a sweep"
        )
    if metric_id not in specification.metrics:
        raise ScientificValidationError("sensitivity metric must be requested")
    base_parameter = next(
        (item for item in specification.parameters if item.name == parameter_name), None
    )
    if base_parameter is None:
        raise ScientificValidationError(
            f"unknown sensitivity parameter {parameter_name!r}"
        )
    selected_adapter = _select_adapter(specification, adapter)
    values = (
        base_parameter.value,
        base_parameter.value - perturbation,
        base_parameter.value + perturbation,
    )
    if not np.isfinite(values).all():
        raise ScientificValidationError("perturbed parameter values must remain finite")
    runs = tuple(
        _run_point(
            specification,
            _replace_parameter(specification.parameters, parameter_name, value),
            selected_adapter,
            run_index=index,
        )
        for index, value in enumerate(values)
    )
    failed_indices = tuple(
        run.run_index for run in runs if run.status == ExperimentRunStatus.FAILED
    )
    if failed_indices:
        return SensitivityFailure(
            experiment_id=specification.experiment_id,
            parameter_name=parameter_name,
            metric_id=metric_id,
            perturbation=perturbation,
            runs=runs,
            failed_run_indices=failed_indices,
            message=(
                "Central difference requires successful base, lower, and upper "
                f"runs; failed indices: {failed_indices}. No one-sided estimate used."
            ),
        )
    base_run, lower_run, upper_run = runs
    metric_values = tuple(
        next(item for item in run.metrics if item.metric_id == metric_id)
        for run in runs
    )
    if len({item.unit for item in metric_values}) != 1:
        return SensitivityFailure(
            experiment_id=specification.experiment_id,
            parameter_name=parameter_name,
            metric_id=metric_id,
            perturbation=perturbation,
            runs=runs,
            failed_run_indices=(),
            message="Metric units differ across sensitivity evaluations.",
        )
    base_metric, lower_metric, upper_metric = metric_values
    derivative = (upper_metric.value - lower_metric.value) / (2.0 * perturbation)
    return SensitivityResult(
        experiment_id=specification.experiment_id,
        parameter_name=parameter_name,
        parameter_unit=base_parameter.unit,
        metric_id=metric_id,
        metric_unit=base_metric.unit,
        base_value=base_parameter.value,
        perturbation=perturbation,
        lower_parameter_value=values[1],
        upper_parameter_value=values[2],
        base_metric_value=base_metric.value,
        lower_metric_value=lower_metric.value,
        upper_metric_value=upper_metric.value,
        derivative_estimate=derivative,
        derivative_unit=(
            f"{base_metric.unit}/({base_parameter.unit})"
            if "/" in base_parameter.unit
            else f"{base_metric.unit}/{base_parameter.unit}"
        ),
        estimate_kind=(
            "candidate_difference"
            if isinstance(selected_adapter, AlgebraicExperimentAdapter)
            else (
                "mesh_dependent_metric_difference"
                if isinstance(selected_adapter, BVPExperimentAdapter)
                else "metric_sensitivity"
            )
        ),
        candidate_branch_consistency=(
            "not_established"
            if isinstance(selected_adapter, AlgebraicExperimentAdapter)
            else "not_applicable"
        ),
        runs=(base_run, lower_run, upper_run),
    )


def make_oscillator_experiment(
    *,
    experiment_id: str = "oscillator_damping_sweep",
    mass_kg: float = 1.0,
    damping_coefficient_kg_per_s: float = 0.2,
    stiffness_n_per_m: float = 4.0,
    initial_displacement_m: float = 0.1,
    initial_velocity_m_per_s: float = 0.0,
    duration_s: float = 10.0,
    sweep_values: tuple[float, ...] | None = (0.0, 0.1, 0.2, 0.4),
    metrics: tuple[str, ...] = (
        "maximum_absolute_displacement_m",
        "final_displacement_m",
        "final_velocity_m_per_s",
        "maximum_absolute_velocity_m_per_s",
    ),
    solver_configuration: ODESolverConfiguration = (
        _DEFAULT_EXPERIMENT_SOLVER_CONFIGURATION
    ),
    output_sampling: OutputSampling = _DEFAULT_EXPERIMENT_OUTPUT_SAMPLING,
    analytical_reference: AnalyticalReferenceRequest | None = None,
    retain_ode_trajectory: bool = False,
) -> ExperimentSpecification:
    """Build the documented oscillator demonstration without hiding its inputs."""
    parameters = (
        NamedValue(name="mass_kg", value=mass_kg, unit="kg"),
        NamedValue(
            name="damping_coefficient_kg_per_s",
            value=damping_coefficient_kg_per_s,
            unit="kg/s",
        ),
        NamedValue(name="stiffness_n_per_m", value=stiffness_n_per_m, unit="N/m"),
        NamedValue(
            name="initial_displacement_m", value=initial_displacement_m, unit="m"
        ),
        NamedValue(
            name="initial_velocity_m_per_s", value=initial_velocity_m_per_s, unit="m/s"
        ),
        NamedValue(name="duration_s", value=duration_s, unit="s"),
    )
    sweep = (
        None
        if sweep_values is None
        else ParameterSweep(
            parameter_name="damping_coefficient_kg_per_s", values=sweep_values
        )
    )
    return ExperimentSpecification(
        experiment_id=experiment_id,
        name="Damped oscillator parameter experiment",
        description=(
            "Deterministic parameter study using the existing oscillator model."
        ),
        model_id=OscillatorExperimentAdapter.model_id,
        parameters=parameters,
        sweep=sweep,
        metrics=metrics,
        solver_configuration=solver_configuration,
        output_sampling=output_sampling,
        analytical_reference=analytical_reference,
        retain_ode_trajectory=retain_ode_trajectory,
        assumptions=(
            "Linear spring with constant stiffness.",
            "Linear viscous damping with constant coefficient.",
            "No external driving force.",
        ),
    )


def make_pendulum_experiment(
    *,
    experiment_id: str = "pendulum_initial_angle_sweep",
    mass_kg: float = 0.2,
    length_m: float = 0.75,
    damping_coefficient_kg_m2_per_s: float = 0.01,
    gravity_m_per_s2: float = 9.81,
    initial_angle_rad: float = 0.2,
    initial_angular_velocity_rad_per_s: float = 0.0,
    duration_s: float = 8.0,
    sweep_parameter_name: str = "initial_angle_rad",
    sweep_values: tuple[float, ...] | None = (0.05, 0.2, 0.5),
    metrics: tuple[str, ...] = (
        "final_angle_rad",
        "final_angular_velocity_rad_per_s",
        "maximum_absolute_angle_rad",
        "maximum_absolute_angular_velocity_rad_per_s",
    ),
    solver_configuration: ODESolverConfiguration = (
        _DEFAULT_EXPERIMENT_SOLVER_CONFIGURATION
    ),
    output_sampling: OutputSampling = _DEFAULT_EXPERIMENT_OUTPUT_SAMPLING,
    analytical_reference: AnalyticalReferenceRequest | None = None,
    retain_ode_trajectory: bool = False,
) -> ExperimentSpecification:
    """Build an experiment for the existing full-nonlinear pendulum solver.

    Damping ``b`` follows the model equation and is measured in kg*m^2/s;
    angles are radians. Parameter domains are validated by the existing model
    on every independent run.
    """
    parameters = (
        NamedValue(name="mass_kg", value=mass_kg, unit="kg"),
        NamedValue(name="length_m", value=length_m, unit="m"),
        NamedValue(
            name="damping_coefficient_kg_m2_per_s",
            value=damping_coefficient_kg_m2_per_s,
            unit="kg*m^2/s",
        ),
        NamedValue(name="gravity_m_per_s2", value=gravity_m_per_s2, unit="m/s^2"),
        NamedValue(name="initial_angle_rad", value=initial_angle_rad, unit="rad"),
        NamedValue(
            name="initial_angular_velocity_rad_per_s",
            value=initial_angular_velocity_rad_per_s,
            unit="rad/s",
        ),
        NamedValue(name="duration_s", value=duration_s, unit="s"),
    )
    sweep = (
        None
        if sweep_values is None
        else ParameterSweep(parameter_name=sweep_parameter_name, values=sweep_values)
    )
    return ExperimentSpecification(
        experiment_id=experiment_id,
        name="Nonlinear pendulum parameter experiment",
        description=(
            "Parameter study using the existing full-sine damped pendulum model."
        ),
        model_id=PendulumExperimentAdapter.model_id,
        parameters=parameters,
        sweep=sweep,
        metrics=metrics,
        solver_configuration=solver_configuration,
        output_sampling=output_sampling,
        analytical_reference=analytical_reference,
        retain_ode_trajectory=retain_ode_trajectory,
        assumptions=(
            "Point mass on a massless rigid rod with a fixed pivot.",
            "Uniform gravitational field and constant model parameters.",
            "Linear viscous damping torque and no external driving torque.",
            "The nonlinear restoring term sin(theta) is retained.",
        ),
    )


def make_nonlinear_spring_experiment(
    *,
    experiment_id: str = "cubic_spring_force_sweep",
    linear_stiffness_n_per_m: float = 100.0,
    cubic_stiffness_n_per_m3: float = 1000.0,
    applied_force_n: float = 1.001,
    initial_displacement_m: float = 0.005,
    sweep_values: tuple[float, ...] | None = (0.8, 1.001, 2.008),
    metrics: tuple[str, ...] = (
        "solution:displacement",
        "residual:force_balance",
        "residual_norm",
    ),
    solver_configuration: AlgebraicSolverConfiguration = (
        _DEFAULT_EXPERIMENT_ALGEBRAIC_CONFIGURATION
    ),
    analytical_reference: AnalyticalReferenceRequest | None = None,
    reference_solution: tuple[float, ...] | None = None,
) -> ExperimentSpecification:
    """Build an experiment around the existing cubic-spring problem factory.

    The sweep values alter the applied force. Every solve uses the fixed
    initial guess supplied here; previous candidates are never reused.
    """
    problem = make_nonlinear_spring_equilibrium_problem(
        linear_stiffness_n_per_m=linear_stiffness_n_per_m,
        cubic_stiffness_n_per_m3=cubic_stiffness_n_per_m3,
        applied_force_n=applied_force_n,
        initial_displacement_m=initial_displacement_m,
    )
    parameters = tuple(
        NamedValue(
            name=parameter.name,
            value=parameter.value,
            unit=parameter.unit or "unknown",
        )
        for parameter in problem.parameters
    )
    sweep = (
        None
        if sweep_values is None
        else ParameterSweep(parameter_name="applied_force_n", values=sweep_values)
    )
    return ExperimentSpecification(
        experiment_id=experiment_id,
        name="Cubic-spring nonlinear equilibrium experiment",
        description=(
            "Local root-solving experiment for a linear-plus-cubic spring force law."
        ),
        model_id=AlgebraicExperimentAdapter.model_id,
        parameters=parameters,
        sweep=sweep,
        metrics=metrics,
        solver_configuration=solver_configuration,
        algebraic_problem=problem,
        algebraic_parameter_constraints=(
            AlgebraicParameterConstraint(
                parameter_name="linear_stiffness_n_per_m", domain="positive"
            ),
            AlgebraicParameterConstraint(
                parameter_name="cubic_stiffness_n_per_m3", domain="positive"
            ),
        ),
        algebraic_reference_solution=reference_solution,
        analytical_reference=analytical_reference,
        assumptions=problem.assumptions or (),
        equation_record_ids=problem.equation_record_ids,
        source_reference_ids=problem.source_reference_ids,
    )


def make_heat_conduction_experiment(
    *,
    experiment_id: str = "heat_conduction_boundary_sweep",
    length_m: float = 1.0,
    left_temperature_k: float = 300.0,
    right_temperature_k: float = 400.0,
    mesh_points: int = 5,
    sweep_parameter_name: str = "right_temperature_k",
    sweep_values: tuple[float, ...] | None = (350.0, 400.0, 450.0),
    metrics: tuple[str, ...] = (
        "domain_start",
        "domain_end",
        "mesh_point_count",
        "sampled_minimum:temperature",
        "sampled_minimum_location:temperature",
        "sampled_maximum:temperature",
        "sampled_maximum_location:temperature",
        "left_endpoint:temperature",
        "right_endpoint:temperature",
        "boundary_residual_max_abs",
    ),
    solver_configuration: BoundaryValueSolverConfiguration = (
        _DEFAULT_EXPERIMENT_BVP_CONFIGURATION
    ),
    bvp_reference: BVPReferenceRequest | None = None,
) -> ExperimentSpecification:
    """Build a sweepable experiment around the existing steady heat BVP."""
    problem = make_linear_heat_conduction_problem(
        length_m=length_m,
        left_temperature_k=left_temperature_k,
        right_temperature_k=right_temperature_k,
        mesh_points=mesh_points,
    )
    parameters = problem.parameters
    sweep = (
        None
        if sweep_values is None
        else ParameterSweep(parameter_name=sweep_parameter_name, values=sweep_values)
    )
    return ExperimentSpecification(
        experiment_id=experiment_id,
        name="Steady heat-conduction boundary-value experiment",
        description=(
            "Spatial parameter study using the existing constant-conductivity "
            "steady heat-conduction problem."
        ),
        model_id=BVPExperimentAdapter.model_id,
        parameters=parameters,
        sweep=sweep,
        metrics=metrics,
        solver_configuration=solver_configuration,
        bvp_problem=problem,
        bvp_parameter_constraints=(
            BVPParameterConstraint(
                parameter_name="left_temperature_k", domain="non_negative"
            ),
            BVPParameterConstraint(
                parameter_name="right_temperature_k", domain="non_negative"
            ),
        ),
        bvp_reference=bvp_reference,
        assumptions=problem.assumptions,
        equation_record_ids=problem.equation_record_ids,
        source_reference_ids=problem.source_reference_ids,
    )
