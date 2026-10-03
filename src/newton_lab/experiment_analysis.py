"""Read-only summaries and comparisons for existing experiment records.

Analysis results are projections of experiment records. They preserve source
identity and family-specific diagnostics without changing solver results or
assigning a common success meaning across solver families.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import isfinite, sqrt
from typing import Literal

import numpy as np
from pydantic import Field, field_validator

from newton_lab.experiments import (
    AlgebraicRunDiagnostics,
    AnalyticalReferenceComparison,
    BVPReferenceComparison,
    BVPRunDiagnostics,
    ExperimentRecord,
    ExperimentRecordStatus,
    ExperimentReproducibility,
    ExperimentRun,
    ExperimentRunStatus,
    MetricValue,
)
from newton_lab.simulation import NamedValue, SimulationModel

ProblemFamily = Literal["ode", "algebraic", "boundary_value", "unknown"]
MetricStatus = Literal["valid", "missing", "non_finite", "invalid"]
ObservationStatus = Literal[
    "valid", "execution_failed", "missing", "non_finite", "invalid"
]
ComparisonStatus = Literal["direct", "conditional", "not_comparable"]


class AnalysisCount(SimulationModel):
    """A count associated with one explicit category."""

    category: str
    count: int = Field(ge=0)


class DuplicateExperimentIdentifier(SimulationModel):
    """Duplicate experiment identifier retained rather than silently merged."""

    experiment_id: str
    record_count: int = Field(ge=2)
    record_indices: tuple[int, ...] = Field(min_length=2)


class ParameterConfigurationCount(SimulationModel):
    """Number of run records with one exact ordered model configuration."""

    model_id: str
    parameters: tuple[NamedValue, ...]
    run_count: int = Field(ge=1)


class ExperimentCollectionSummary(SimulationModel):
    """Counts of source records and runs; duplicate records remain separate."""

    record_count: int = Field(ge=0)
    duplicate_experiment_ids: tuple[DuplicateExperimentIdentifier, ...]
    run_count: int = Field(ge=0)
    successful_run_count: int = Field(ge=0)
    failed_run_count: int = Field(ge=0)
    counts_by_family: tuple[AnalysisCount, ...]
    counts_by_model_id: tuple[AnalysisCount, ...]
    counts_by_execution_status: tuple[AnalysisCount, ...]
    solver_success_count: int = Field(ge=0)
    solver_failure_count: int = Field(ge=0)
    solver_success_unavailable_count: int = Field(ge=0)
    numerical_acceptance_count: int = Field(ge=0)
    numerical_rejection_count: int = Field(ge=0)
    numerical_acceptance_unavailable_count: int = Field(ge=0)
    parameter_configurations: tuple[ParameterConfigurationCount, ...]


class RunAnalysisProjection(SimulationModel):
    """Common read-only view that retains the run's family-specific evidence."""

    record_index: int = Field(ge=0)
    experiment_id: str
    record_status: ExperimentRecordStatus
    model_id: str
    problem_family: ProblemFamily
    run_index: int = Field(ge=0)
    parameters: tuple[NamedValue, ...]
    swept_parameter_name: str | None
    swept_parameter_value: float | None
    execution_status: ExperimentRunStatus
    metrics: tuple[MetricValue, ...]
    solver_id: str
    solver_method: str
    solver_success: bool | None
    numerical_acceptance: bool | None
    function_evaluations: int | None
    failure_kind: str | None
    error_type: str | None
    message: str
    reproducibility: ExperimentReproducibility | None
    algebraic_diagnostics: AlgebraicRunDiagnostics | None
    bvp_diagnostics: BVPRunDiagnostics | None
    reference_comparison: AnalyticalReferenceComparison | None
    bvp_reference_comparison: BVPReferenceComparison | None


class MetricObservation(SimulationModel):
    """One metric lookup tied to its source record and run."""

    record_index: int = Field(ge=0)
    experiment_id: str
    record_status: ExperimentRecordStatus
    model_id: str
    problem_id: str | None
    problem_family: ProblemFamily
    run_index: int = Field(ge=0)
    parameters: tuple[NamedValue, ...]
    swept_parameter_name: str | None
    swept_parameter_value: float | None
    metric_id: str
    value: float | None = None
    unit: str | None = None
    definition: str | None = None
    metric_status: MetricStatus
    observation_status: ObservationStatus
    execution_status: ExperimentRunStatus
    solver_id: str
    solver_method: str
    solver_success: bool | None
    numerical_acceptance: bool | None
    function_evaluations: int | None
    failure_kind: str | None
    error_type: str | None
    message: str
    reproducibility: ExperimentReproducibility | None
    algebraic_diagnostics: AlgebraicRunDiagnostics | None
    bvp_diagnostics: BVPRunDiagnostics | None
    reference_comparison: AnalyticalReferenceComparison | None
    bvp_reference_comparison: BVPReferenceComparison | None

    @field_validator("value")
    @classmethod
    def finite_value(cls, value: float | None) -> float | None:
        if value is not None and not isfinite(value):
            raise ValueError("analysis observations store only finite numeric values")
        return value


class MetricSourceDeclaration(SimulationModel):
    """One source metric's meaning in an explicit compatibility declaration."""

    model_id: str = Field(min_length=1)
    problem_id: str | None = None
    metric_id: str = Field(min_length=1)
    unit: str = Field(min_length=1)
    source_definition: str = Field(min_length=1)

    @field_validator("unit")
    @classmethod
    def known_unit(cls, value: str) -> str:
        if value == "unknown":
            raise ValueError("compatibility declarations require an explicit unit")
        return value


class MetricCompatibilityDeclaration(SimulationModel):
    """Caller-declared mapping needed to compare different model metrics."""

    shared_definition: str = Field(min_length=1)
    unit: str = Field(min_length=1)
    sources: tuple[MetricSourceDeclaration, ...] = Field(min_length=1)
    assumptions: tuple[str, ...] = ()

    @field_validator("unit")
    @classmethod
    def shared_unit_known(cls, value: str) -> str:
        if value == "unknown":
            raise ValueError("compatibility declarations require an explicit unit")
        return value

    @field_validator("sources")
    @classmethod
    def unique_sources(
        cls, values: tuple[MetricSourceDeclaration, ...]
    ) -> tuple[MetricSourceDeclaration, ...]:
        keys = tuple(
            (item.model_id, item.problem_id, item.metric_id) for item in values
        )
        if len(keys) != len(set(keys)):
            raise ValueError("compatibility declaration sources must be unique")
        if any(item.unit != values[0].unit for item in values):
            raise ValueError("compatibility declaration cannot convert units")
        return values


class MetricComparison(SimulationModel):
    """Metric observations plus an explicit decision about comparability."""

    status: ComparisonStatus
    observations: tuple[MetricObservation, ...]
    unit: str | None = None
    shared_definition: str | None = None
    assumptions: tuple[str, ...] = ()
    explanation: str


class MetricStatistics(SimulationModel):
    """Descriptive statistics using population standard deviation (ddof=0)."""

    comparison_status: ComparisonStatus
    observation_count: int = Field(ge=0)
    valid_count: int = Field(ge=0)
    excluded_count: int = Field(ge=0)
    execution_failed_count: int = Field(ge=0)
    missing_metric_count: int = Field(ge=0)
    non_finite_metric_count: int = Field(ge=0)
    invalid_metric_count: int = Field(ge=0)
    unit: str | None = None
    mean: float | None = None
    minimum: float | None = None
    maximum: float | None = None
    population_standard_deviation: float | None = None
    status: Literal["computed", "no_valid_observations", "not_comparable"]


class SweepPointAnalysis(SimulationModel):
    """One point in the exact execution order of a recorded sweep."""

    run_index: int = Field(ge=0)
    parameter_value: float | None = None
    repeated_parameter_value: bool
    observation: MetricObservation


class ParameterSweepAnalysis(SimulationModel):
    """Ordered sweep projection; unreliable sweep metadata remains explicit."""

    experiment_id: str
    parameter_name: str
    metric_id: str
    ordering_reliable: bool
    message: str
    points: tuple[SweepPointAnalysis, ...]


def _problem_family(model_id: str) -> ProblemFamily:
    if model_id in {"damped_harmonic_oscillator", "damped_nonlinear_pendulum"}:
        return "ode"
    if model_id == "nonlinear_algebraic_system":
        return "algebraic"
    if model_id == "boundary_value_problem":
        return "boundary_value"
    return "unknown"


def _metric_definition(
    model_id: str, problem_id: str | None, metric_id: str
) -> str | None:
    """Return only definitions established by the existing adapters."""
    oscillator_definitions = {
        "maximum_absolute_displacement_m": (
            "Maximum absolute displacement over requested ODE samples."
        ),
        "final_displacement_m": "Displacement at the final requested ODE sample.",
        "final_velocity_m_per_s": "Velocity at the final requested ODE sample.",
        "maximum_absolute_velocity_m_per_s": (
            "Maximum absolute velocity over requested ODE samples."
        ),
    }
    pendulum_definitions = {
        "final_angle_rad": "Angle at the final requested pendulum sample.",
        "final_angular_velocity_rad_per_s": (
            "Angular velocity at the final requested pendulum sample."
        ),
        "maximum_absolute_angle_rad": (
            "Maximum absolute angle over requested pendulum samples."
        ),
        "maximum_absolute_angular_velocity_rad_per_s": (
            "Maximum absolute angular velocity over requested pendulum samples."
        ),
    }
    if model_id == "damped_harmonic_oscillator":
        return oscillator_definitions.get(metric_id)
    if model_id == "damped_nonlinear_pendulum":
        return pendulum_definitions.get(metric_id)
    if model_id == "nonlinear_algebraic_system":
        if metric_id == "residual_norm":
            return "Euclidean L2 norm of raw, unscaled residual entries."
        if metric_id.startswith("solution:"):
            return (
                f"Solver candidate component {metric_id.removeprefix('solution:')!r}."
            )
        if metric_id.startswith("residual:"):
            return (
                "Signed raw residual component "
                f"{metric_id.removeprefix('residual:')!r}."
            )
        return None
    if model_id == "boundary_value_problem":
        if metric_id in {"domain_start", "domain_end", "mesh_point_count"}:
            return {
                "domain_start": "First returned spatial mesh coordinate.",
                "domain_end": "Last returned spatial mesh coordinate.",
                "mesh_point_count": "Number of returned spatial mesh coordinates.",
            }[metric_id]
        if metric_id.startswith("sampled_minimum:"):
            return "Minimum sampled nodal value for the named BVP state."
        if metric_id.startswith("sampled_maximum:"):
            return "Maximum sampled nodal value for the named BVP state."
        if metric_id.startswith("sampled_minimum_location:"):
            return "First returned coordinate attaining the sampled state minimum."
        if metric_id.startswith("sampled_maximum_location:"):
            return "First returned coordinate attaining the sampled state maximum."
        if metric_id.startswith(("left_endpoint:", "right_endpoint:")):
            return "State value at the corresponding returned domain endpoint."
        if metric_id.startswith("boundary_residual"):
            return "Raw boundary callback residual metric."
    del problem_id
    return None


def _metric_for_run(run: ExperimentRun, metric_id: str) -> tuple[MetricValue, ...]:
    return tuple(metric for metric in run.metrics if metric.metric_id == metric_id)


def _observation(
    record: ExperimentRecord,
    record_index: int,
    run: ExperimentRun,
    metric_id: str,
) -> MetricObservation:
    matches = _metric_for_run(run, metric_id)
    metric = matches[0] if len(matches) == 1 else None
    metric_status: MetricStatus
    value: float | None = None
    unit: str | None = None
    if len(matches) > 1:
        metric_status = "invalid"
    elif metric is None:
        metric_status = "missing"
    else:
        raw_value = metric.value
        raw_unit = metric.unit
        unit = raw_unit if isinstance(raw_unit, str) else None
        try:
            numeric_value = float(raw_value)
        except (TypeError, ValueError, OverflowError):
            metric_status = "invalid"
        else:
            if not isfinite(numeric_value):
                metric_status = "non_finite"
            elif unit is None or not unit:
                metric_status = "invalid"
            else:
                metric_status = "valid"
                value = numeric_value
    if run.status == ExperimentRunStatus.FAILED:
        observation_status: ObservationStatus = "execution_failed"
    elif metric_status == "valid":
        observation_status = "valid"
    elif metric_status == "missing":
        observation_status = "missing"
    elif metric_status == "non_finite":
        observation_status = "non_finite"
    else:
        observation_status = "invalid"
    problem_id = run.reproducibility.problem_id
    return MetricObservation(
        record_index=record_index,
        experiment_id=record.specification.experiment_id,
        record_status=record.status,
        model_id=record.specification.model_id,
        problem_id=problem_id,
        problem_family=_problem_family(record.specification.model_id),
        run_index=run.run_index,
        parameters=run.parameters,
        swept_parameter_name=run.swept_parameter_name,
        swept_parameter_value=run.swept_parameter_value,
        metric_id=metric_id,
        value=value,
        unit=unit,
        definition=(
            None
            if metric is None
            else _metric_definition(
                record.specification.model_id, problem_id, metric_id
            )
        ),
        metric_status=metric_status,
        observation_status=observation_status,
        execution_status=run.status,
        solver_id=run.solver_id,
        solver_method=run.solver_method,
        solver_success=run.solver_success,
        numerical_acceptance=run.numerical_acceptance,
        function_evaluations=run.function_evaluations,
        failure_kind=None if run.failure_kind is None else run.failure_kind.value,
        error_type=run.error_type,
        message=run.message,
        reproducibility=run.reproducibility,
        algebraic_diagnostics=run.algebraic_diagnostics,
        bvp_diagnostics=run.bvp_diagnostics,
        reference_comparison=run.reference_comparison,
        bvp_reference_comparison=run.bvp_reference_comparison,
    )


def summarize_experiments(
    records: Sequence[ExperimentRecord],
) -> ExperimentCollectionSummary:
    """Summarize records without merging duplicate experiment identifiers."""
    items = tuple(records)
    by_id: dict[str, list[int]] = {}
    family_counts: dict[str, int] = {}
    model_counts: dict[str, int] = {}
    status_counts: dict[str, int] = {}
    config_counts: dict[
        tuple[object, ...], tuple[str, tuple[NamedValue, ...], int]
    ] = {}
    succeeded = failed = solver_true = solver_false = solver_unavailable = 0
    acceptance_true = acceptance_false = acceptance_unavailable = 0
    run_count = 0
    for record_index, record in enumerate(items):
        experiment_id = record.specification.experiment_id
        by_id.setdefault(experiment_id, []).append(record_index)
        model_id = record.specification.model_id
        family = _problem_family(model_id)
        family_counts[family] = family_counts.get(family, 0) + len(record.runs)
        model_counts[model_id] = model_counts.get(model_id, 0) + len(record.runs)
        for run in record.runs:
            run_count += 1
            status = run.status.value
            status_counts[status] = status_counts.get(status, 0) + 1
            if run.status == ExperimentRunStatus.SUCCEEDED:
                succeeded += 1
            else:
                failed += 1
            if run.solver_success is True:
                solver_true += 1
            elif run.solver_success is False:
                solver_false += 1
            else:
                solver_unavailable += 1
            if run.numerical_acceptance is True:
                acceptance_true += 1
            elif run.numerical_acceptance is False:
                acceptance_false += 1
            else:
                acceptance_unavailable += 1
            parameter_key = tuple(
                (item.name, item.value, item.unit) for item in run.parameters
            )
            config_key = (model_id, parameter_key)
            if config_key in config_counts:
                old_model, old_parameters, count = config_counts[config_key]
                config_counts[config_key] = (old_model, old_parameters, count + 1)
            else:
                config_counts[config_key] = (model_id, run.parameters, 1)
    duplicates = tuple(
        DuplicateExperimentIdentifier(
            experiment_id=experiment_id,
            record_count=len(indices),
            record_indices=tuple(indices),
        )
        for experiment_id, indices in by_id.items()
        if len(indices) > 1
    )
    return ExperimentCollectionSummary(
        record_count=len(items),
        duplicate_experiment_ids=duplicates,
        run_count=run_count,
        successful_run_count=succeeded,
        failed_run_count=failed,
        counts_by_family=tuple(
            AnalysisCount(category=key, count=count)
            for key, count in family_counts.items()
        ),
        counts_by_model_id=tuple(
            AnalysisCount(category=key, count=count)
            for key, count in model_counts.items()
        ),
        counts_by_execution_status=tuple(
            AnalysisCount(category=key, count=count)
            for key, count in status_counts.items()
        ),
        solver_success_count=solver_true,
        solver_failure_count=solver_false,
        solver_success_unavailable_count=solver_unavailable,
        numerical_acceptance_count=acceptance_true,
        numerical_rejection_count=acceptance_false,
        numerical_acceptance_unavailable_count=acceptance_unavailable,
        parameter_configurations=tuple(
            ParameterConfigurationCount(
                model_id=model_id, parameters=parameters, run_count=count
            )
            for model_id, parameters, count in config_counts.values()
        ),
    )


def project_experiment_runs(
    records: Sequence[ExperimentRecord],
) -> tuple[RunAnalysisProjection, ...]:
    """Create a common view while retaining family-specific diagnostics."""
    return tuple(
        RunAnalysisProjection(
            record_index=record_index,
            experiment_id=record.specification.experiment_id,
            record_status=record.status,
            model_id=record.specification.model_id,
            problem_family=_problem_family(record.specification.model_id),
            run_index=run.run_index,
            parameters=run.parameters,
            swept_parameter_name=run.swept_parameter_name,
            swept_parameter_value=run.swept_parameter_value,
            execution_status=run.status,
            metrics=run.metrics,
            solver_id=run.solver_id,
            solver_method=run.solver_method,
            solver_success=run.solver_success,
            numerical_acceptance=run.numerical_acceptance,
            function_evaluations=run.function_evaluations,
            failure_kind=None if run.failure_kind is None else run.failure_kind.value,
            error_type=run.error_type,
            message=run.message,
            reproducibility=run.reproducibility,
            algebraic_diagnostics=run.algebraic_diagnostics,
            bvp_diagnostics=run.bvp_diagnostics,
            reference_comparison=run.reference_comparison,
            bvp_reference_comparison=run.bvp_reference_comparison,
        )
        for record_index, record in enumerate(records)
        for run in record.runs
    )


def extract_metric_observations(
    records: Sequence[ExperimentRecord], metric_id: str
) -> tuple[MetricObservation, ...]:
    """Extract one named metric in caller record and run order."""
    if not metric_id:
        raise ValueError("metric_id must not be empty")
    return tuple(
        _observation(record, record_index, run, metric_id)
        for record_index, record in enumerate(records)
        for run in record.runs
    )


def _source_key(observation: MetricObservation) -> tuple[str, str | None, str]:
    return (observation.model_id, observation.problem_id, observation.metric_id)


def compare_metric_observations(
    observations: Sequence[MetricObservation],
    *,
    compatibility: MetricCompatibilityDeclaration | None = None,
) -> MetricComparison:
    """Check metric comparability without unit conversion or implicit mapping."""
    items = tuple(observations)
    if len(items) < 2:
        return MetricComparison(
            status="not_comparable",
            observations=items,
            explanation="At least two metric observations are required for comparison.",
        )
    origins = tuple(dict.fromkeys(_source_key(item) for item in items))
    present = tuple(
        item for item in items if item.metric_status in {"valid", "non_finite"}
    )
    units = {item.unit for item in present}
    if compatibility is not None:
        declared = {
            (source.model_id, source.problem_id, source.metric_id): source
            for source in compatibility.sources
        }
        if not all(origin in declared for origin in origins):
            return MetricComparison(
                status="not_comparable",
                observations=items,
                explanation=(
                    "The explicit compatibility declaration does not cover every "
                    "observed model, problem, and metric."
                ),
            )
        declared_units = {declared[origin].unit for origin in origins}
        if len(declared_units) != 1 or compatibility.unit not in declared_units:
            return MetricComparison(
                status="not_comparable",
                observations=items,
                explanation=(
                    "Unit conversion is not supported; all declared source units "
                    "must exactly match the shared unit."
                ),
            )
        for item in present:
            source = declared[_source_key(item)]
            if item.unit != source.unit:
                return MetricComparison(
                    status="not_comparable",
                    observations=items,
                    explanation="A metric unit conflicts with its declaration.",
                )
            if (
                item.definition is not None
                and item.definition != source.source_definition
            ):
                return MetricComparison(
                    status="not_comparable",
                    observations=items,
                    explanation=(
                        "A declared source definition conflicts with the "
                        "adapter-defined metric meaning."
                    ),
                )
        return MetricComparison(
            status="conditional",
            observations=items,
            unit=compatibility.unit,
            shared_definition=compatibility.shared_definition,
            assumptions=compatibility.assumptions,
            explanation=(
                "Comparison uses the caller's explicit mapping and assumptions; "
                "the analysis layer does not independently validate that mapping."
            ),
        )
    if len(origins) != 1:
        return MetricComparison(
            status="not_comparable",
            observations=items,
            explanation=(
                "Different model, problem, or metric sources require an explicit "
                "compatibility declaration; matching units alone are insufficient."
            ),
        )
    if len(units) != 1 or None in units or "unknown" in units:
        return MetricComparison(
            status="not_comparable",
            observations=items,
            explanation="A single known, consistent unit is required for comparison.",
        )
    definitions = {item.definition for item in present}
    if len(definitions) != 1 or None in definitions:
        return MetricComparison(
            status="not_comparable",
            observations=items,
            explanation="The metric definition is unavailable or inconsistent.",
        )
    if any(item.metric_status == "invalid" for item in items):
        return MetricComparison(
            status="not_comparable",
            observations=items,
            explanation="Malformed metric observations prevent comparison.",
        )
    return MetricComparison(
        status="direct",
        observations=items,
        unit=next(iter(units)),
        shared_definition=next(iter(definitions)),
        explanation=(
            "Observations share the same adapter-defined model, problem, metric, "
            "definition, and unit. Failed and unavailable values remain marked."
        ),
    )


def summarize_metric(comparison: MetricComparison) -> MetricStatistics:
    """Summarize compatible finite successful observations with population SD."""
    items = comparison.observations
    valid = tuple(
        item
        for item in items
        if item.execution_status == ExperimentRunStatus.SUCCEEDED
        and item.metric_status == "valid"
        and item.value is not None
    )
    failed = sum(item.execution_status == ExperimentRunStatus.FAILED for item in items)
    missing = sum(item.metric_status == "missing" for item in items)
    non_finite = sum(item.metric_status == "non_finite" for item in items)
    invalid = sum(item.metric_status == "invalid" for item in items)
    base = {
        "comparison_status": comparison.status,
        "observation_count": len(items),
        "valid_count": len(valid),
        "excluded_count": len(items) - len(valid),
        "execution_failed_count": failed,
        "missing_metric_count": missing,
        "non_finite_metric_count": non_finite,
        "invalid_metric_count": invalid,
        "unit": comparison.unit,
    }
    if comparison.status == "not_comparable":
        return MetricStatistics(
            **base,
            status="not_comparable",
        )
    if not valid:
        return MetricStatistics(
            **base,
            status="no_valid_observations",
        )
    values = np.asarray([item.value for item in valid], dtype=np.float64)
    mean = float(np.mean(values))
    return MetricStatistics(
        **base,
        mean=mean,
        minimum=float(np.min(values)),
        maximum=float(np.max(values)),
        population_standard_deviation=float(sqrt(float(np.mean((values - mean) ** 2)))),
        status="computed",
    )


def analyze_parameter_sweep(
    record: ExperimentRecord,
    *,
    parameter_name: str,
    metric_id: str,
) -> ParameterSweepAnalysis:
    """Analyze recorded sweep points in stored order without sorting or filling."""
    declared_sweep = record.specification.sweep
    reliable = (
        declared_sweep is not None and declared_sweep.parameter_name == parameter_name
    )
    message = "Sweep order follows the stored run order."
    if declared_sweep is None:
        message = "The experiment record does not declare a parameter sweep."
    elif declared_sweep.parameter_name != parameter_name:
        message = (
            f"The declared sweep parameter is {declared_sweep.parameter_name!r}, "
            f"not {parameter_name!r}."
        )
    points: list[SweepPointAnalysis] = []
    values_seen: set[float] = set()
    if reliable and declared_sweep is not None:
        recorded_values = tuple(run.swept_parameter_value for run in record.runs)
        names_match = all(
            run.swept_parameter_name == parameter_name
            and run.swept_parameter_value is not None
            for run in record.runs
        )
        values_match = recorded_values == declared_sweep.values
        if (
            len(record.runs) != len(declared_sweep.values)
            or not names_match
            or not values_match
        ):
            reliable = False
            message = (
                "Stored run sweep metadata does not match the declared ordered "
                "sweep; values are shown where recorded but ordering is not certified."
            )
    for _record_index, run in enumerate(record.runs):
        parameter_value = run.swept_parameter_value
        if parameter_value is None or run.swept_parameter_name != parameter_name:
            parameter_value = next(
                (item.value for item in run.parameters if item.name == parameter_name),
                None,
            )
        repeated = parameter_value is not None and parameter_value in values_seen
        if parameter_value is not None:
            values_seen.add(parameter_value)
        observation = _observation(record, 0, run, metric_id)
        points.append(
            SweepPointAnalysis(
                run_index=run.run_index,
                parameter_value=parameter_value,
                repeated_parameter_value=repeated,
                observation=observation.model_copy(update={"record_index": 0}),
            )
        )
    return ParameterSweepAnalysis(
        experiment_id=record.specification.experiment_id,
        parameter_name=parameter_name,
        metric_id=metric_id,
        ordering_reliable=reliable,
        message=message,
        points=tuple(points),
    )
