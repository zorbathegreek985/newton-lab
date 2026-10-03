"""Reusable, family-aware contracts for numerical simulation workflows.

The core representations are deliberately independent of SciPy and of any one
physical model. They describe problems and results; solver adapters live beside
the relevant model implementations.
"""

from __future__ import annotations

import platform
import sys
from dataclasses import dataclass
from enum import StrEnum
from importlib.metadata import PackageNotFoundError, version
from typing import TYPE_CHECKING, Literal, Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from newton_lab.exceptions import ScientificValidationError

if TYPE_CHECKING:
    from newton_lab.knowledge.registry import EquationRecord


class SimulationModel(BaseModel):
    """Immutable validated metadata shared by simulation records."""

    model_config = ConfigDict(frozen=True, extra="forbid", str_strip_whitespace=True)


class NamedValue(SimulationModel):
    """A named physical value and its declared unit (without conversion)."""

    name: str = Field(min_length=1)
    value: float
    unit: str = Field(min_length=1)

    @field_validator("value")
    @classmethod
    def finite_value(cls, value: float) -> float:
        if isinstance(value, bool) or not np.isfinite(value):
            raise ValueError("value must be finite")
        return value


class StateVariable(SimulationModel):
    """Name and unit for one simulated state variable."""

    name: str = Field(min_length=1)
    unit: str = Field(min_length=1)


class OutputSampling(SimulationModel):
    """Requested output grid; distinct from the solver's internal steps."""

    point_count: int = Field(default=1001, ge=2, strict=True)
    policy: Literal["uniform_including_endpoints"] = "uniform_including_endpoints"


class KnowledgeMappingKind(StrEnum):
    """How executable behavior relates to a descriptive equation record."""

    EXACT_MODEL_IMPLEMENTATION = "exact_model_implementation"
    APPROXIMATION = "approximation_under_conditions"
    NUMERICAL_REALIZATION = "numerical_realization_of_stated_model"


class KnowledgeMappingStatus(StrEnum):
    CANDIDATE = "candidate_requires_review"
    REVIEWED = "reviewed"


class VariableCorrespondence(SimulationModel):
    """Declared correspondence between an equation variable and model state."""

    equation_symbol: str = Field(min_length=1)
    model_state: str = Field(min_length=1)
    equation_unit: str | None = None
    model_unit: str = Field(min_length=1)


class ParameterCorrespondence(SimulationModel):
    """Declared correspondence between an equation parameter and model input."""

    equation_symbol: str = Field(min_length=1)
    model_parameter: str = Field(min_length=1)
    equation_unit: str | None = None
    model_unit: str = Field(min_length=1)


class KnowledgeModelMapping(SimulationModel):
    """Explicit model-to-equation correspondence, independent of verification."""

    equation_record_id: str = Field(pattern=r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
    kind: KnowledgeMappingKind
    variable_correspondences: tuple[VariableCorrespondence, ...] = ()
    parameter_correspondences: tuple[ParameterCorrespondence, ...] = ()
    assumptions: tuple[str, ...] = ()
    validity_limits: tuple[str, ...] = ()
    source_reference_ids: tuple[str, ...] = ()
    status: KnowledgeMappingStatus = KnowledgeMappingStatus.CANDIDATE
    reviewed_by: str | None = None
    review_notes: str | None = None

    @model_validator(mode="after")
    def review_is_auditable(self) -> KnowledgeModelMapping:
        if self.status == KnowledgeMappingStatus.REVIEWED and (
            not self.reviewed_by
            or not self.reviewed_by.strip()
            or not self.review_notes
            or not self.review_notes.strip()
        ):
            raise ValueError("reviewed mappings require reviewer and review notes")
        if len(self.source_reference_ids) != len(set(self.source_reference_ids)):
            raise ValueError("mapping source reference IDs must be unique")
        return self

    def record_review(self, *, reviewer: str, notes: str) -> KnowledgeModelMapping:
        """Return a reviewed mapping with a named reviewer and rationale."""
        if self.status != KnowledgeMappingStatus.CANDIDATE:
            raise ScientificValidationError("only candidate mappings can be reviewed")
        return KnowledgeModelMapping.model_validate(
            {
                **self.model_dump(mode="python"),
                "status": KnowledgeMappingStatus.REVIEWED,
                "reviewed_by": reviewer,
                "review_notes": notes,
            }
        )


class SimulationProblem(SimulationModel):
    """Model-independent description of a concrete simulation request.

    ``family`` is an open string so later ODE, PDE, algebraic, optimization,
    stochastic, and event-driven problem types can be added without changing
    this common envelope. Family-specific contracts belong in their adapters.
    """

    model_id: str = Field(pattern=r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
    name: str = Field(min_length=1)
    description: str = ""
    family: str = Field(min_length=1)
    model_version: str = Field(min_length=1)
    independent_variable: StateVariable | None = None
    interval: tuple[float, float] | None = None
    states: tuple[StateVariable, ...] = ()
    parameters: tuple[NamedValue, ...] = ()
    initial_conditions: tuple[NamedValue, ...] = ()
    sampling: OutputSampling | None = None
    boundary_conditions: tuple[str, ...] = ()
    assumptions: tuple[str, ...] = ()
    validity_limits: tuple[str, ...] = ()
    equation_record_ids: tuple[str, ...] = ()
    source_reference_ids: tuple[str, ...] = ()
    knowledge_mappings: tuple[KnowledgeModelMapping, ...] = ()

    @field_validator("interval")
    @classmethod
    def valid_interval(
        cls, value: tuple[float, float] | None
    ) -> tuple[float, float] | None:
        if value is None:
            return None
        if len(value) != 2 or not np.isfinite(value).all() or value[1] <= value[0]:
            raise ValueError("interval must contain finite increasing endpoints")
        return value

    @model_validator(mode="after")
    def unique_names_and_ids(self) -> SimulationProblem:
        for label, items in (
            ("state names", tuple(state.name for state in self.states)),
            ("parameter names", tuple(item.name for item in self.parameters)),
            (
                "initial-condition names",
                tuple(item.name for item in self.initial_conditions),
            ),
            ("equation record IDs", self.equation_record_ids),
            ("source reference IDs", self.source_reference_ids),
            (
                "knowledge mapping equation IDs",
                tuple(item.equation_record_id for item in self.knowledge_mappings),
            ),
        ):
            if len(items) != len(set(items)):
                raise ValueError(f"{label} must be unique")
        state_names = {state.name for state in self.states}
        unknown_initials = {
            value.name for value in self.initial_conditions
        } - state_names
        if state_names and unknown_initials:
            raise ValueError(
                "initial conditions refer to unknown states: "
                f"{sorted(unknown_initials)}"
            )
        unmapped_equations = {
            item.equation_record_id for item in self.knowledge_mappings
        } - set(self.equation_record_ids)
        if unmapped_equations:
            raise ValueError(
                "knowledge mappings must be listed in equation_record_ids: "
                f"{sorted(unmapped_equations)}"
            )
        state_names = {state.name for state in self.states}
        parameter_names = {parameter.name for parameter in self.parameters}
        for mapping in self.knowledge_mappings:
            unknown_states = {
                item.model_state for item in mapping.variable_correspondences
            } - state_names
            unknown_parameters = {
                item.model_parameter for item in mapping.parameter_correspondences
            } - parameter_names
            if unknown_states or unknown_parameters:
                raise ValueError(
                    "knowledge mapping refers to unknown model quantities: "
                    f"states={sorted(unknown_states)}, "
                    f"parameters={sorted(unknown_parameters)}"
                )
        return self


def validate_knowledge_mappings(
    problem: SimulationProblem,
    equation_records: tuple[EquationRecord, ...],
) -> None:
    """Validate mapping targets and provenance without promoting candidates."""
    records_by_id = {record.equation_id: record for record in equation_records}
    for mapping in problem.knowledge_mappings:
        record = records_by_id.get(mapping.equation_record_id)
        if record is None:
            raise ScientificValidationError(
                f"knowledge mapping refers to unknown equation "
                f"{mapping.equation_record_id!r}"
            )
        source_ids = {source.source_id for source in record.source_references}
        unknown_source_ids = set(mapping.source_reference_ids) - source_ids
        if unknown_source_ids:
            raise ScientificValidationError(
                "knowledge mapping refers to unknown source IDs: "
                f"{sorted(unknown_source_ids)}"
            )


class ODESolverConfiguration(SimulationModel):
    """Explicit settings for the supported SciPy initial-value ODE adapter."""

    method: Literal["RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"] = "RK45"
    relative_tolerance: float = Field(default=1e-8, gt=0.0)
    absolute_tolerance: float = Field(default=1e-10, gt=0.0)
    maximum_step: float | None = Field(default=None, gt=0.0)

    @field_validator("relative_tolerance", "absolute_tolerance", "maximum_step")
    @classmethod
    def finite_setting(cls, value: float | None) -> float | None:
        if value is not None and not np.isfinite(value):
            raise ValueError("solver settings must be finite")
        return value


class SolverCapabilities(SimulationModel):
    """Declared scope and limitations of a concrete solver implementation."""

    solver_id: str
    solver_version: str
    problem_families: tuple[str, ...]
    supported_methods: tuple[str, ...]
    adaptive_internal_steps: bool
    supports_requested_output_grid: bool
    supports_events: bool
    stiffness_support: str
    limitations: tuple[str, ...]


class SimulationStatus(StrEnum):
    SUCCESS = "success"
    FAILED = "failed"


class ReproducibilityRecord(SimulationModel):
    """Environment and configuration details useful for later reproduction."""

    python_version: str
    platform: str
    numpy_version: str
    scipy_version: str
    model_id: str
    model_version: str
    solver_id: str
    solver_version: str
    solver_method: str
    solver_relative_tolerance: float
    solver_absolute_tolerance: float
    solver_maximum_step: float | None
    interval: tuple[float, float] | None
    requested_output_sampling: OutputSampling | None
    model_parameters: tuple[NamedValue, ...]
    initial_conditions: tuple[NamedValue, ...]
    state_variables: tuple[StateVariable, ...]
    boundary_conditions: tuple[str, ...]
    assumptions: tuple[str, ...]
    validity_limits: tuple[str, ...]
    equation_record_ids: tuple[str, ...]
    source_reference_ids: tuple[str, ...]
    random_seed: int | None = None


class IntegrationDiagnostics(SimulationModel):
    """Solver outcome details shared by supported ODE adapters."""

    function_evaluations: int = Field(ge=0)
    accepted_internal_steps: int | None = Field(default=None, ge=0)
    message: str


@dataclass(frozen=True, slots=True)
class StateSeries:
    """Sampled numerical values for one named state in its declared unit."""

    name: str
    unit: str
    values: NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class SimulationResult:
    """Common result envelope; arrays use requested output times only."""

    independent_variable: NDArray[np.float64]
    independent_variable_name: str
    independent_variable_unit: str
    states: tuple[StateSeries, ...]
    parameters: tuple[NamedValue, ...]
    problem: SimulationProblem
    solver_configuration: ODESolverConfiguration
    solver_capabilities: SolverCapabilities
    status: SimulationStatus
    diagnostics: IntegrationDiagnostics
    reproducibility: ReproducibilityRecord
    assumptions: tuple[str, ...]
    validity_limits: tuple[str, ...]

    def __post_init__(self) -> None:
        samples = self.independent_variable
        if (
            self.problem.sampling is None
            or self.problem.interval is None
            or self.problem.independent_variable is None
        ):
            raise ScientificValidationError(
                "time-series results require interval, sampling, and "
                "independent-variable metadata"
            )
        if samples.ndim != 1 or not np.isfinite(samples).all():
            raise ScientificValidationError(
                "result independent_variable must be a finite one-dimensional array"
            )
        if len(samples) != self.problem.sampling.point_count:
            raise ScientificValidationError(
                "result sample count must match the problem output grid"
            )
        if len(self.states) != len(self.problem.states):
            raise ScientificValidationError(
                "result states must match problem state metadata"
            )
        if (
            samples[0] != self.problem.interval[0]
            or samples[-1] != self.problem.interval[1]
            or not np.all(np.diff(samples) > 0)
        ):
            raise ScientificValidationError(
                "result samples must cover the declared interval in increasing order"
            )
        if any(
            (series.name, series.unit) != (metadata.name, metadata.unit)
            for series, metadata in zip(self.states, self.problem.states, strict=True)
        ):
            raise ScientificValidationError(
                "result state names and units must match the problem"
            )
        if any(
            state.values.ndim != 1
            or len(state.values) != len(samples)
            or not np.isfinite(state.values).all()
            for state in self.states
        ):
            raise ScientificValidationError(
                "each successful state series must be finite and match the output grid"
            )
        if self.status != SimulationStatus.SUCCESS:
            raise ScientificValidationError(
                "failed simulations are represented by SimulationFailure"
            )


@dataclass(frozen=True, slots=True)
class SimulationFailure:
    """Failure outcome for adapter workflows; legacy functions still raise."""

    problem: SimulationProblem
    solver_configuration: ODESolverConfiguration
    error_type: str
    message: str
    reproducibility: ReproducibilityRecord
    status: Literal[SimulationStatus.FAILED] = SimulationStatus.FAILED


SimulationOutcome = SimulationResult | SimulationFailure


@runtime_checkable
class SimulationAdapter(Protocol):
    """Protocol for model-specific adapters over compatible solver families."""

    capabilities: SolverCapabilities

    def run(
        self,
        problem: SimulationProblem,
        configuration: ODESolverConfiguration,
    ) -> SimulationOutcome:
        """Run one validated problem without interpreting equation text."""


def _library_version(distribution: str) -> str:
    try:
        return version(distribution)
    except PackageNotFoundError:
        return "unknown"


def make_reproducibility_record(
    problem: SimulationProblem,
    configuration: ODESolverConfiguration,
    capabilities: SolverCapabilities,
    *,
    random_seed: int | None = None,
) -> ReproducibilityRecord:
    """Capture current environment and declared numerical settings."""
    return ReproducibilityRecord(
        python_version=sys.version,
        platform=platform.platform(),
        numpy_version=np.__version__,
        scipy_version=_library_version("scipy"),
        model_id=problem.model_id,
        model_version=problem.model_version,
        solver_id=capabilities.solver_id,
        solver_version=capabilities.solver_version,
        solver_method=configuration.method,
        solver_relative_tolerance=configuration.relative_tolerance,
        solver_absolute_tolerance=configuration.absolute_tolerance,
        solver_maximum_step=configuration.maximum_step,
        interval=problem.interval,
        requested_output_sampling=problem.sampling,
        model_parameters=problem.parameters,
        initial_conditions=problem.initial_conditions,
        state_variables=problem.states,
        boundary_conditions=problem.boundary_conditions,
        assumptions=problem.assumptions,
        validity_limits=problem.validity_limits,
        equation_record_ids=problem.equation_record_ids,
        source_reference_ids=problem.source_reference_ids,
        random_seed=random_seed,
    )
