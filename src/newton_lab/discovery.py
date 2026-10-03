"""Evidence-aware discovery from existing experiment records.

The first workflow characterizes recorded scalar metrics and parameter-sweep
endpoint changes, compares explicitly declared equation structures, and
produces a bounded application hypothesis for research. It performs no solver
execution, symbolic reasoning, or automatic registry mutation.
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import StrEnum
from typing import Literal

import numpy as np
from pydantic import Field

from newton_lab.experiment_analysis import (
    ExperimentCollectionSummary,
    MetricObservation,
    analyze_parameter_sweep,
    extract_metric_observations,
    summarize_experiments,
)
from newton_lab.experiments import (
    ExperimentRecord,
    ExperimentRunStatus,
    ODETrajectorySnapshot,
)
from newton_lab.knowledge import (
    EquationRecord,
    EquationRelationship,
    EvidenceCategory,
    KnowledgeBase,
    Linearity,
    RelationshipReviewStatus,
    RelationshipType,
    build_example_knowledge_base,
    compare_equations,
    propose_relationships,
)
from newton_lab.knowledge.structure import StructureComparison
from newton_lab.simulation import NamedValue, SimulationModel
from newton_lab.visualization import build_experiment_report


class ObservationBasis(StrEnum):
    """How a discovery descriptor relates to recorded numerical results."""

    COMPUTED = "computed_from_solver_output"
    DERIVED = "derived_from_recorded_observations"
    ASSUMED = "assumed_from_model_metadata"
    UNAVAILABLE = "unavailable_or_insufficient_data"


class ObservationAvailability(StrEnum):
    """Whether a descriptor has usable numerical evidence."""

    AVAILABLE = "available"
    FAILED_RUN = "failed_run"
    MISSING_METRIC = "missing_metric"
    NON_FINITE_METRIC = "non_finite_metric"
    INVALID_METRIC = "invalid_metric"
    UNRELIABLE_SWEEP = "unreliable_sweep_metadata"
    INSUFFICIENT_ENDPOINTS = "insufficient_valid_endpoints"
    INSUFFICIENT_SAMPLES = "insufficient_or_unretained_trajectory"


class DiscoveryEvidenceStatus(StrEnum):
    """Qualified evidence levels used by discovery outputs."""

    ESTABLISHED_RELATIONSHIP = "established_relationship"
    CANDIDATE_STRUCTURAL_CONNECTION = "candidate_structural_connection"
    PROPOSED_APPLICATION_HYPOTHESIS = "proposed_application_hypothesis"
    VALIDATED_APPLICATION = "validated_application"


class BehavioralObservation(SimulationModel):
    """One computed or derived behavioral descriptor with source provenance."""

    observation_id: str
    descriptor: str
    source_experiment_id: str
    source_record_index: int = Field(ge=0)
    source_run_ids: tuple[str, ...] = Field(min_length=1)
    model_id: str
    equation_record_ids: tuple[str, ...] = ()
    equation_identity_method: str
    parameters: tuple[NamedValue, ...] = ()
    metric_id: str
    swept_parameter_name: str | None = None
    swept_parameter_values: tuple[float, ...] = ()
    swept_parameter_unit: str | None = None
    metric_definition: str | None = None
    value: float | None = None
    unit: str | None = None
    basis: ObservationBasis
    availability: ObservationAvailability
    solver_success: bool | None = None
    numerical_acceptance: bool | None = None
    method: str
    limitations: tuple[str, ...] = ()
    state_name: str | None = None
    time_interval_s: tuple[float, float] | None = None
    sample_count: int | None = Field(default=None, ge=0)


class DiscoveryStructureMatch(SimulationModel):
    """Explainable, non-mutating match between recorded equation structures."""

    match_id: str
    source_experiment_id: str
    source_record_index: int = Field(ge=0)
    source_equation_id: str
    target_equation_id: str
    evidence_status: DiscoveryEvidenceStatus
    relationship_type: RelationshipType
    relationship_source_equation_id: str | None = None
    relationship_target_equation_id: str | None = None
    relationship_id: str | None = None
    relationship_review_status: RelationshipReviewStatus | None = None
    relationship_source_reference_ids: tuple[str, ...] = ()
    shared_features: tuple[str, ...] = ()
    supporting_observation_ids: tuple[str, ...] = Field(min_length=1)
    conditions: tuple[str, ...] = ()
    explanation: str
    limitations: tuple[str, ...]


class ApplicationHypothesis(SimulationModel):
    """A testable cross-domain proposal that is never auto-validated."""

    hypothesis_id: str
    statement: str
    source_experiment_id: str
    source_record_index: int = Field(ge=0)
    source_model_id: str
    source_equation_id: str
    target_domain: str
    mathematical_structure: str
    source_behavior_observation_ids: tuple[str, ...] = Field(min_length=1)
    supporting_match_ids: tuple[str, ...] = Field(min_length=1)
    proposed_mechanism: str
    assumptions: tuple[str, ...] = Field(min_length=1)
    similarities: tuple[str, ...] = Field(min_length=1)
    differences: tuple[str, ...] = Field(min_length=1)
    test_plan: tuple[str, ...] = Field(min_length=1)
    evidence_status: Literal[
        DiscoveryEvidenceStatus.PROPOSED_APPLICATION_HYPOTHESIS
    ] = DiscoveryEvidenceStatus.PROPOSED_APPLICATION_HYPOTHESIS
    validated_application: Literal[False] = False


class DiscoveryResult(SimulationModel):
    """Complete discovery output and a deterministic Phase 15 Markdown report."""

    observations: tuple[BehavioralObservation, ...]
    structure_matches: tuple[DiscoveryStructureMatch, ...]
    hypotheses: tuple[ApplicationHypothesis, ...]
    collection_summary: ExperimentCollectionSummary
    limitations: tuple[str, ...]
    report_markdown: str


def _equation_identity(
    record: ExperimentRecord,
    equations: dict[str, EquationRecord],
) -> tuple[tuple[str, ...], str, tuple[str, ...]]:
    declared = record.specification.equation_record_ids
    if declared:
        known = tuple(identifier for identifier in declared if identifier in equations)
        unknown = tuple(
            identifier for identifier in declared if identifier not in equations
        )
        notes = (
            (f"Unrecognized equation record IDs: {', '.join(unknown)}.",)
            if unknown
            else ()
        )
        return known, "explicit experiment equation_record_ids", notes
    model_id = record.specification.model_id
    if model_id in equations:
        return (model_id,), "exact model_id equals registry equation_id", ()
    return (
        (),
        "unresolved; no explicit equation ID and no exact model/equation ID match",
        ("No equation structure was associated; discovery did not guess from names.",),
    )


def _metric_availability(item: MetricObservation) -> ObservationAvailability:
    if item.execution_status == ExperimentRunStatus.FAILED:
        return ObservationAvailability.FAILED_RUN
    if item.metric_status == "valid" and item.value is not None:
        return ObservationAvailability.AVAILABLE
    if item.metric_status == "missing":
        return ObservationAvailability.MISSING_METRIC
    if item.metric_status == "non_finite":
        return ObservationAvailability.NON_FINITE_METRIC
    return ObservationAvailability.INVALID_METRIC


def _metric_observation(
    item: MetricObservation,
    *,
    record_index: int,
    equation_ids: tuple[str, ...],
    identity_method: str,
    identity_notes: tuple[str, ...],
) -> BehavioralObservation:
    availability = _metric_availability(item)
    available = availability == ObservationAvailability.AVAILABLE
    limitations = list(identity_notes)
    limitations.append(
        "Endpoint difference is descriptive; it does not establish a causal or "
        "monotonic parameter response."
    )
    if item.definition and "sample" in item.definition.lower():
        limitations.append(
            "Metric is defined over returned samples or mesh nodes; behavior "
            "between them is not established."
        )
    if item.solver_success is not True:
        limitations.append("Solver success is not recorded as true for this run.")
    if item.numerical_acceptance is not True:
        limitations.append(
            "Numerical acceptance is false or unavailable; interpret the "
            "metric accordingly."
        )
    return BehavioralObservation(
        observation_id=(
            f"{record_index}:{item.experiment_id}#{item.run_index}:{item.metric_id}"
        ),
        descriptor="recorded_metric",
        source_experiment_id=item.experiment_id,
        source_record_index=record_index,
        source_run_ids=(f"{record_index}:{item.experiment_id}#{item.run_index}",),
        model_id=item.model_id,
        equation_record_ids=equation_ids,
        equation_identity_method=identity_method,
        parameters=item.parameters,
        metric_id=item.metric_id,
        swept_parameter_name=item.swept_parameter_name,
        swept_parameter_values=(
            () if item.swept_parameter_value is None else (item.swept_parameter_value,)
        ),
        swept_parameter_unit=next(
            (
                parameter.unit
                for parameter in item.parameters
                if parameter.name == item.swept_parameter_name
            ),
            None,
        ),
        metric_definition=item.definition,
        value=item.value if available else None,
        unit=item.unit,
        basis=(
            ObservationBasis.COMPUTED if available else ObservationBasis.UNAVAILABLE
        ),
        availability=availability,
        solver_success=item.solver_success,
        numerical_acceptance=item.numerical_acceptance,
        method="Reuse Phase 14 metric extraction; no metric was recalculated.",
        limitations=tuple(limitations),
    )


def _temporal_observation(
    record: ExperimentRecord,
    run_index: int,
    record_index: int,
    equation_ids: tuple[str, ...],
    identity_method: str,
    identity_notes: tuple[str, ...],
    descriptor: str,
    trajectory: ODETrajectorySnapshot | None,
) -> BehavioralObservation:
    run = record.runs[run_index]
    run_id = f"{record_index}:{record.specification.experiment_id}#{run_index}"
    default_state = (
        "angle"
        if record.specification.model_id == "damped_nonlinear_pendulum"
        else "displacement"
    )
    default_unit = "rad" if default_state == "angle" else "m"
    state_name = default_state if trajectory is None else trajectory.state_names[0]
    unit = default_unit if trajectory is None else trajectory.state_units[0]
    limits = list(identity_notes)
    value: float | None = None
    availability = ObservationAvailability.INSUFFICIENT_SAMPLES
    method = "No retained sampled trajectory was available for this ODE run."
    interval: tuple[float, float] | None = None
    sample_count: int | None = None
    if run.status == ExperimentRunStatus.FAILED:
        availability = ObservationAvailability.FAILED_RUN
        method = "The source ODE run failed; no temporal descriptor was computed."
    elif trajectory is not None:
        sample_count = len(trajectory.time_s)
        interval = (trajectory.time_s[0], trajectory.time_s[-1])
        times = np.asarray(trajectory.time_s, dtype=np.float64)
        states = np.asarray(trajectory.state_values[0], dtype=np.float64)
        state_name = trajectory.state_names[0]
        unit = trajectory.state_units[0]
        if descriptor == "observed_period":
            crossings = [
                float(
                    times[index]
                    + (times[index + 1] - times[index])
                    * (-states[index])
                    / (states[index + 1] - states[index])
                )
                for index in range(len(states) - 1)
                if states[index] < 0.0
                and states[index + 1] >= 0.0
                and states[index + 1] != states[index]
            ]
            periods = np.diff(crossings)
            if len(periods) >= 2 and np.isfinite(periods).all():
                value = float(np.median(periods))
                availability = ObservationAvailability.AVAILABLE
                method = (
                    "Median interval between linearly interpolated upward crossings "
                    "of the model's zero-equilibrium state."
                )
            else:
                method = (
                    "Period unavailable: at least three upward zero crossings "
                    "(two observed cycles) are required."
                )
        elif descriptor == "amplitude_change":
            if len(states) >= 10:
                window_size = max(3, int(np.ceil(len(states) * 0.2)))
                first_amplitude = float(np.max(np.abs(states[:window_size])))
                last_amplitude = float(np.max(np.abs(states[-window_size:])))
                value = last_amplitude - first_amplitude
                availability = ObservationAvailability.AVAILABLE
                method = (
                    "Difference between maximum absolute sampled state in the first "
                    "and last equal-count 20% windows. This is not a fitted decay law."
                )
            else:
                method = (
                    "Amplitude comparison unavailable: at least ten samples are "
                    "required."
                )
        else:
            if len(states) >= 10:
                window_size = max(3, int(np.ceil(len(states) * 0.2)))
                first_amplitude = float(np.max(np.abs(states[:window_size])))
                last_amplitude = float(np.max(np.abs(states[-window_size:])))
                if first_amplitude > 0.0:
                    value = float(last_amplitude <= 0.05 * first_amplitude)
                    availability = ObservationAvailability.AVAILABLE
                    method = (
                        "Finite-window indicator: 1 when the last 20% window's "
                        "maximum absolute state is at most 5% of the first 20% "
                        "window maximum; otherwise 0. This threshold is descriptive."
                    )
                else:
                    method = (
                        "Return indicator unavailable: first-window amplitude is zero."
                    )
            else:
                method = (
                    "Return indicator unavailable: at least ten samples are required."
                )
        if (
            sample_count is not None
            and sample_count < 12
            and descriptor == "observed_period"
        ):
            availability = ObservationAvailability.INSUFFICIENT_SAMPLES
            value = None
            method = "Period unavailable: at least 12 returned samples are required."
    if availability != ObservationAvailability.AVAILABLE:
        limits.append(
            "Unavailable means returned samples did not support this descriptor; "
            "it is not evidence of non-oscillation or instability."
        )
    limits.append(
        "Descriptors summarize returned numerical samples, not the exact "
        "continuous solution or independent physical validation."
    )
    return BehavioralObservation(
        observation_id=f"{run_id}:temporal:{descriptor}",
        descriptor=descriptor,
        source_experiment_id=record.specification.experiment_id,
        source_record_index=record_index,
        source_run_ids=(run_id,),
        model_id=record.specification.model_id,
        equation_record_ids=equation_ids,
        equation_identity_method=identity_method,
        parameters=run.parameters,
        metric_id=descriptor,
        swept_parameter_name=run.swept_parameter_name,
        swept_parameter_values=(
            () if run.swept_parameter_value is None else (run.swept_parameter_value,)
        ),
        swept_parameter_unit=next(
            (
                item.unit
                for item in run.parameters
                if item.name == run.swept_parameter_name
            ),
            None,
        ),
        value=value,
        unit=(
            "s"
            if descriptor == "observed_period"
            else "1"
            if descriptor == "finite_window_return_indicator"
            else unit
        ),
        basis=ObservationBasis.COMPUTED
        if value is not None
        else ObservationBasis.UNAVAILABLE,
        availability=availability,
        solver_success=run.solver_success,
        numerical_acceptance=run.numerical_acceptance,
        method=method,
        limitations=tuple(limits),
        state_name=state_name,
        time_interval_s=interval,
        sample_count=sample_count,
    )


def _sweep_change_observation(
    record: ExperimentRecord,
    *,
    record_index: int,
    metric_id: str,
    equation_ids: tuple[str, ...],
    identity_method: str,
    identity_notes: tuple[str, ...],
) -> BehavioralObservation:
    assert record.specification.sweep is not None
    parameter_name = record.specification.sweep.parameter_name
    analysis = analyze_parameter_sweep(
        record, parameter_name=parameter_name, metric_id=metric_id
    )
    source_run_ids = tuple(
        f"{record_index}:{record.specification.experiment_id}#{point.run_index}"
        for point in analysis.points
    )
    first = analysis.points[0].observation if analysis.points else None
    last = analysis.points[-1].observation if analysis.points else None
    unit = None if first is None else first.unit
    parameters = () if first is None else first.parameters
    limitations = list(identity_notes)
    if not analysis.ordering_reliable:
        availability = ObservationAvailability.UNRELIABLE_SWEEP
        limitations.append(f"Sweep ordering unreliable: {analysis.message}")
    elif (
        first is None
        or last is None
        or first.metric_status != "valid"
        or last.metric_status != "valid"
        or first.execution_status != ExperimentRunStatus.SUCCEEDED
        or last.execution_status != ExperimentRunStatus.SUCCEEDED
        or first.value is None
        or last.value is None
    ):
        availability = ObservationAvailability.INSUFFICIENT_ENDPOINTS
        limitations.append(
            "The first and last recorded sweep points must both have successful, "
            "finite metrics."
        )
    else:
        availability = ObservationAvailability.AVAILABLE
    value = (
        last.value - first.value
        if availability == ObservationAvailability.AVAILABLE
        and first is not None
        and last is not None
        and first.value is not None
        and last.value is not None
        else None
    )
    if any(point.repeated_parameter_value for point in analysis.points):
        limitations.append(
            "Repeated parameter values are retained in their original run order."
        )
    if any(
        point.observation.execution_status == ExperimentRunStatus.FAILED
        for point in analysis.points
    ):
        limitations.append(
            "One or more sweep runs failed; this endpoint difference does not "
            "summarize the failed point or imply a continuous response."
        )
    return BehavioralObservation(
        observation_id=(
            f"{record_index}:{record.specification.experiment_id}:sweep:"
            f"{parameter_name}:{metric_id}"
        ),
        descriptor="sweep_endpoint_change",
        source_experiment_id=record.specification.experiment_id,
        source_record_index=record_index,
        source_run_ids=source_run_ids,
        model_id=record.specification.model_id,
        equation_record_ids=equation_ids,
        equation_identity_method=identity_method,
        parameters=parameters,
        metric_id=metric_id,
        swept_parameter_name=parameter_name,
        swept_parameter_values=tuple(
            point.parameter_value
            for point in analysis.points
            if point.parameter_value is not None
        ),
        swept_parameter_unit=(
            next(
                (
                    parameter.unit
                    for parameter in first.parameters
                    if parameter.name == parameter_name
                ),
                None,
            )
            if first is not None
            else None
        ),
        metric_definition=(None if first is None else first.definition),
        value=value,
        unit=unit,
        basis=(
            ObservationBasis.DERIVED
            if value is not None
            else ObservationBasis.UNAVAILABLE
        ),
        availability=availability,
        solver_success=(
            True
            if first is not None
            and last is not None
            and first.solver_success is True
            and last.solver_success is True
            else None
        ),
        numerical_acceptance=None,
        method=(
            f"Last minus first {metric_id} in recorded {parameter_name} sweep order; "
            "uses Phase 14 observations without sorting or interpolation."
        ),
        limitations=tuple(limitations),
    )


def _relationship_status(
    relationship: EquationRelationship,
) -> DiscoveryEvidenceStatus:
    if (
        relationship.review_status == RelationshipReviewStatus.ACCEPTED
        and relationship.source_reference_ids
    ):
        return DiscoveryEvidenceStatus.ESTABLISHED_RELATIONSHIP
    return DiscoveryEvidenceStatus.CANDIDATE_STRUCTURAL_CONNECTION


def _structure_matches(
    records: tuple[ExperimentRecord, ...],
    observations: tuple[BehavioralObservation, ...],
    knowledge_base: KnowledgeBase,
) -> tuple[DiscoveryStructureMatch, ...]:
    equations = {
        entry.record.equation_id: entry.record for entry in knowledge_base.entries
    }
    proposed = propose_relationships(knowledge_base)
    all_relationships = (*knowledge_base.relationships, *proposed)
    observations_by_record: dict[int, tuple[BehavioralObservation, ...]] = {}
    for record_index, _record in enumerate(records):
        observations_by_record[record_index] = tuple(
            item for item in observations if item.source_record_index == record_index
        )

    matches: list[DiscoveryStructureMatch] = []
    seen: set[tuple[str, str, str]] = set()
    structural_features = {"classifications", "ode_order", "derivative_operators"}
    for record_index, record in enumerate(records):
        source_ids = tuple(
            dict.fromkeys(
                identifier
                for observation in observations_by_record.get(record_index, ())
                for identifier in observation.equation_record_ids
            )
        )
        for source_id in source_ids:
            source = equations[source_id]
            relevant_observations = tuple(
                item.observation_id
                for item in observations_by_record.get(record_index, ())
                if item.availability == ObservationAvailability.AVAILABLE
                and item.equation_record_ids
            )
            if not relevant_observations:
                continue
            for target_id, target in equations.items():
                if target_id == source_id:
                    continue
                comparison: StructureComparison = compare_equations(
                    source, target, knowledge_base=knowledge_base
                )
                shared = tuple(
                    item
                    for item in comparison.feature_comparisons
                    if item.feature_name in structural_features
                    and item.status.value in {"shared", "partially_shared"}
                    and (item.shared_values or item.left_value == item.right_value)
                )
                equation_relation = next(
                    (
                        item
                        for item in source.related_equations
                        if item.target_equation_id == target_id
                    ),
                    None,
                )
                if equation_relation is None:
                    equation_relation = next(
                        (
                            item
                            for item in target.related_equations
                            if item.target_equation_id == source_id
                        ),
                        None,
                    )
                established_equation_relation = (
                    equation_relation is not None
                    and equation_relation.evidence.category
                    in {
                        EvidenceCategory.MATHEMATICALLY_DERIVED,
                        EvidenceCategory.ESTABLISHED_PHYSICAL_RESULT,
                    }
                )
                if len(shared) < 2 and not established_equation_relation:
                    continue
                identity = (str(record_index), source_id, target_id)
                if identity in seen:
                    continue
                seen.add(identity)
                relation = next(
                    (
                        item
                        for item in all_relationships
                        if (
                            item.source_equation_id == source_id
                            and item.target_equation_id == target_id
                            if item.directional
                            else {
                                item.source_equation_id,
                                item.target_equation_id,
                            }
                            == {source_id, target_id}
                        )
                    ),
                    None,
                )
                conditions: tuple[str, ...]
                if established_equation_relation and equation_relation is not None:
                    relation_type = equation_relation.relationship
                    if equation_relation in source.related_equations:
                        relationship_source_id = source_id
                        relationship_target_id = target_id
                    else:
                        relationship_source_id = target_id
                        relationship_target_id = source_id
                    status = DiscoveryEvidenceStatus.ESTABLISHED_RELATIONSHIP
                    relationship_id = None
                    review_status = None
                    relationship_source_ids = (
                        equation_relation.evidence.source_reference_ids
                    )
                    explanation = (
                        "Existing EquationRecord relationship with evidence "
                        f"category {equation_relation.evidence.category.value}: "
                        f"{equation_relation.evidence.statement}"
                    )
                    conditions = ()
                elif relation is not None:
                    relation_type = relation.relationship_type
                    relationship_source_id = relation.source_equation_id
                    relationship_target_id = relation.target_equation_id
                    status = _relationship_status(relation)
                    relationship_id = relation.relationship_id
                    review_status = relation.review_status
                    relationship_source_ids = relation.source_reference_ids
                    explanation = (
                        f"Registry relationship {relation.relationship_id!r}: "
                        f"{relation.justification}"
                    )
                    conditions = relation.conditions
                else:
                    shared_operators: tuple[str, ...] = next(
                        (
                            item.shared_values
                            for item in shared
                            if item.feature_name == "derivative_operators"
                        ),
                        (),
                    )
                    relation_type = (
                        RelationshipType.SHARES_MATHEMATICAL_OPERATOR
                        if shared_operators
                        else RelationshipType.SHARES_MATHEMATICAL_STRUCTURE
                    )
                    relationship_source_id = source_id
                    relationship_target_id = target_id
                    status = DiscoveryEvidenceStatus.CANDIDATE_STRUCTURAL_CONNECTION
                    relationship_id = None
                    review_status = RelationshipReviewStatus.CANDIDATE
                    relationship_source_ids = ()
                    explanation = (
                        "Candidate proposed from explicitly shared registry "
                        "structure metadata; no equivalence is inferred."
                    )
                    conditions = ()
                shared_features = tuple(
                    f"{item.feature_name}="
                    + (
                        ", ".join(item.shared_values)
                        if item.shared_values
                        else str(item.left_value)
                    )
                    for item in shared
                )
                limitations = [
                    "Shared structure is not proof of physical equivalence or an "
                    "application mapping.",
                    "Matching uses declared registry metadata only; equation text "
                    "is not parsed.",
                ]
                if not shared:
                    limitations.append(
                        "The registered relationship has no shared comparison "
                        "features recorded on both equation entries."
                    )
                if status != DiscoveryEvidenceStatus.ESTABLISHED_RELATIONSHIP:
                    limitations.append(
                        "This connection remains a review candidate and has not "
                        "been established for the target domain."
                    )
                matches.append(
                    DiscoveryStructureMatch(
                        match_id=(
                            f"{record_index}:{record.specification.experiment_id}:"
                            f"{source_id}:{target_id}"
                        ),
                        source_experiment_id=record.specification.experiment_id,
                        source_record_index=record_index,
                        source_equation_id=source_id,
                        target_equation_id=target_id,
                        evidence_status=status,
                        relationship_type=relation_type,
                        relationship_source_equation_id=relationship_source_id,
                        relationship_target_equation_id=relationship_target_id,
                        relationship_id=relationship_id,
                        relationship_review_status=review_status,
                        relationship_source_reference_ids=relationship_source_ids,
                        shared_features=shared_features,
                        supporting_observation_ids=relevant_observations,
                        conditions=conditions,
                        explanation=explanation,
                        limitations=tuple(limitations),
                    )
                )
    return tuple(matches)


def _application_hypotheses(
    records: tuple[ExperimentRecord, ...],
    observations: tuple[BehavioralObservation, ...],
    matches: tuple[DiscoveryStructureMatch, ...],
    knowledge_base: KnowledgeBase,
) -> tuple[ApplicationHypothesis, ...]:
    equation_map = {
        entry.record.equation_id: entry.record for entry in knowledge_base.entries
    }
    hypotheses: list[ApplicationHypothesis] = []
    for record_index, record in enumerate(records):
        equation_ids = tuple(
            dict.fromkeys(
                identifier
                for item in observations
                if item.source_record_index == record_index
                for identifier in item.equation_record_ids
            )
        )
        for equation_id in equation_ids:
            equation = equation_map[equation_id]
            if equation_id != "damped_harmonic_oscillator":
                continue
            behavior = tuple(
                item
                for item in observations
                if item.source_record_index == record_index
                and item.descriptor == "sweep_endpoint_change"
                and item.availability == ObservationAvailability.AVAILABLE
                and item.metric_id == "maximum_absolute_velocity_m_per_s"
                and item.swept_parameter_name == "damping_coefficient_kg_per_s"
            )
            complete_sweep = tuple(
                item
                for item in observations
                if item.source_record_index == record_index
                and item.descriptor == "recorded_metric"
                and item.metric_id == "maximum_absolute_velocity_m_per_s"
            )
            source_matches = tuple(
                item
                for item in matches
                if item.source_experiment_id == record.specification.experiment_id
                and item.source_equation_id == equation_id
            )
            if (
                not behavior
                or behavior[0].value is None
                or behavior[0].value >= 0.0
                or not source_matches
                or len(complete_sweep) != len(record.runs)
                or any(
                    item.availability != ObservationAvailability.AVAILABLE
                    for item in complete_sweep
                )
                or len(set(behavior[0].swept_parameter_values)) < 2
            ):
                continue
            structure = equation.mathematical_structure
            if (
                structure is None
                or structure.ode_order != 2
                or structure.linearity != Linearity.LINEAR
                or structure.autonomous is not True
                or not {
                    "first_time_derivative",
                    "second_time_derivative",
                }.issubset(structure.derivative_operators or ())
            ):
                continue
            linearity = structure.linearity.value
            structure_text = (
                f"{structure.ode_order}-order {linearity} autonomous ODE "
                "with explicit viscous damping parameter"
            )
            hypotheses.append(
                ApplicationHypothesis(
                    hypothesis_id=(
                        f"{record_index}:{record.specification.experiment_id}:"
                        "damping_risk_control"
                    ),
                    statement=(
                        "The recorded damping sweep's maximum sampled absolute "
                        f"velocity endpoint change was {behavior[0].value:.8g} "
                        f"{behavior[0].unit or 'unit unavailable'}. Test whether a "
                        "damping-like term in a second-order risk-control update "
                        "can reduce peak oscillatory adjustment rates under "
                        "matched information and risk constraints."
                    ),
                    source_experiment_id=record.specification.experiment_id,
                    source_record_index=record_index,
                    source_model_id=record.specification.model_id,
                    source_equation_id=equation_id,
                    target_domain="quantitative_finance_risk_control_research",
                    mathematical_structure=structure_text,
                    source_behavior_observation_ids=tuple(
                        item.observation_id for item in behavior
                    ),
                    supporting_match_ids=tuple(
                        item.match_id for item in source_matches
                    ),
                    proposed_mechanism=(
                        "Represent risk-control state and rate as a second-order "
                        "update with a tunable damping-like term; compare its "
                        "response against an otherwise matched baseline."
                    ),
                    assumptions=(
                        *equation.assumptions,
                        "A target risk-control process can be represented by "
                        "comparable state and rate variables.",
                        "The analogy is a research prompt; no finance-domain "
                        "equation or mapping is established in the registry.",
                    ),
                    similarities=(
                        "Both proposed systems would track a state and its rate "
                        "over ordered time steps.",
                        "A damping-like parameter could be varied and its "
                        "observed response compared.",
                    ),
                    differences=(
                        "The source is a deterministic, unforced mechanical "
                        "model; financial observations are stochastic and "
                        "nonstationary.",
                        "Mechanical viscous damping is not transaction cost, "
                        "market impact, or an economic mechanism.",
                        "No financial data, trading signal, profitability, or "
                        "risk reduction has been validated.",
                    ),
                    test_plan=(
                        "Specify the financial state equation and compare the "
                        "damping-like update with a simple prespecified baseline.",
                        "Use chronological training and untouched out-of-sample "
                        "periods, with transaction costs and market impact.",
                        "Report risk measures, drawdown, turnover, and uncertainty "
                        "across regimes; predefine rejection criteria and control "
                        "for overfitting.",
                    ),
                )
            )
    return tuple(hypotheses)


def _render_discovery_report(
    records: tuple[ExperimentRecord, ...],
    observations: tuple[BehavioralObservation, ...],
    matches: tuple[DiscoveryStructureMatch, ...],
    hypotheses: tuple[ApplicationHypothesis, ...],
    limitations: tuple[str, ...],
) -> str:
    base = build_experiment_report(records, title="Newton Lab discovery report")
    lines = [base.body.rstrip(), "", "## Behavioral observations", ""]
    lines.extend(
        [
            "| Descriptor | Source run(s) | Sweep value | Metric | Value | Unit | "
            "Availability | Basis |",
            "| --- | --- | --- | --- | ---: | --- | --- | --- |",
        ]
    )
    if not observations:
        lines.append("No metric observations were requested or available.")
    for observation in observations:
        value = (
            "unavailable" if observation.value is None else f"{observation.value:.8g}"
        )
        lines.append(
            f"| {observation.descriptor} | "
            f"{', '.join(observation.source_run_ids)} | "
            f"{observation.swept_parameter_values or '—'} "
            f"{observation.swept_parameter_unit or ''} | "
            f"{observation.metric_id} | {value} | "
            f"{observation.unit or 'unit unavailable'} | "
            f"{observation.availability.value} | {observation.basis.value} |"
        )
        lines.append(
            f"\nMethod: {observation.method} "
            f"Identity: {observation.equation_identity_method}."
        )
        if observation.metric_definition:
            lines.append(f" Metric definition: {observation.metric_definition}.")
        if observation.limitations:
            lines.extend(f"\n- Limitation: {note}" for note in observation.limitations)
    lines.extend(["", "## Candidate structure matches", ""])
    if not matches:
        lines.append(
            "No sufficiently supported structure matches were found from the "
            "available equation identities and explicit metadata."
        )
    for match in matches:
        review_status = (
            match.relationship_review_status.value
            if match.relationship_review_status
            else (
                "registry-qualified"
                if match.evidence_status
                == DiscoveryEvidenceStatus.ESTABLISHED_RELATIONSHIP
                else "not recorded"
            )
        )
        lines.extend(
            [
                f"### `{match.source_equation_id}` to `{match.target_equation_id}`",
                "",
                f"- Evidence status: **{match.evidence_status.value}**",
                f"- Relationship: `{match.relationship_type.value}` ({review_status})",
                "- Relationship direction: "
                f"`{match.relationship_source_equation_id}` to "
                f"`{match.relationship_target_equation_id}`",
                "- Shared features: "
                f"{', '.join(match.shared_features) or 'none recorded'}",
                "- Supporting observations: "
                f"{', '.join(match.supporting_observation_ids)}",
                f"- Explanation: {match.explanation}",
            ]
        )
        if match.relationship_source_reference_ids:
            lines.append(
                "- Relationship source references: "
                + ", ".join(
                    f"`{reference_id}`"
                    for reference_id in match.relationship_source_reference_ids
                )
            )
        lines.extend(f"- Condition: {condition}" for condition in match.conditions)
        lines.extend(f"- Limitation: {note}" for note in match.limitations)
        lines.append("")
    lines.extend(["## Cross-domain hypotheses", ""])
    if not hypotheses:
        lines.append(
            "No application hypothesis was emitted because the available "
            "observations and registry matches did not meet the bounded evidence "
            "requirements."
        )
    for hypothesis in hypotheses:
        lines.extend(
            [
                f"### {hypothesis.hypothesis_id}",
                "",
                f"- Status: **{hypothesis.evidence_status.value}; not validated**",
                f"- Hypothesis: {hypothesis.statement}",
                f"- Source: `{hypothesis.source_experiment_id}` / "
                f"`{hypothesis.source_equation_id}`",
                f"- Target: {hypothesis.target_domain}",
                f"- Structure: {hypothesis.mathematical_structure}",
                f"- Mechanism: {hypothesis.proposed_mechanism}",
                "- Supporting observations: "
                f"{', '.join(hypothesis.source_behavior_observation_ids)}",
                f"- Supporting matches: {', '.join(hypothesis.supporting_match_ids)}",
                "- Similarities:",
                *(f"  - {value}" for value in hypothesis.similarities),
                "- Differences:",
                *(f"  - {value}" for value in hypothesis.differences),
                "- Assumptions:",
                *(f"  - {value}" for value in hypothesis.assumptions),
                "- Proposed test:",
                *(f"  - {value}" for value in hypothesis.test_plan),
            ]
        )
    lines.extend(["", "## Discovery limitations", ""])
    lines.extend(f"- {note}" for note in limitations)
    lines.extend(
        [
            "- A successful simulation is not evidence of physical validity or "
            "target-domain usefulness.",
            "- A shared mathematical structure is not proof of equation "
            "equivalence or an application mapping.",
            "- No application is validated by this workflow; any financial "
            "hypothesis requires independent out-of-sample evaluation, baselines, "
            "costs, risk measures, and overfitting safeguards.",
            "",
        ]
    )
    return "\n".join(lines)


def run_discovery(
    records: Sequence[ExperimentRecord],
    *,
    knowledge_base: KnowledgeBase | None = None,
    selected_metrics: Sequence[str] | None = None,
) -> DiscoveryResult:
    """Characterize records, match declared structures, and render a report.

    Metric values are extracted by Phase 14. Retained ODE trajectories support
    qualified sampled period and finite-window amplitude descriptors. The
    workflow does not infer stability, exact periodicity, fitted decay laws,
    causality, or missing values. With no knowledge base, the built-in reviewed
    workspace is used.
    """
    record_tuple = tuple(records)
    kb = knowledge_base or build_example_knowledge_base()
    equations = {entry.record.equation_id: entry.record for entry in kb.entries}
    selected = (
        None if selected_metrics is None else tuple(dict.fromkeys(selected_metrics))
    )
    summary = summarize_experiments(record_tuple)
    observations: list[BehavioralObservation] = []
    limitations: list[str] = []
    selected_found: set[str] = set()
    for record_index, record in enumerate(record_tuple):
        equation_ids, identity_method, identity_notes = _equation_identity(
            record, equations
        )
        limitations.extend(identity_notes)
        metric_ids = tuple(
            metric
            for metric in (
                record.specification.metrics if selected is None else selected
            )
            if metric in record.specification.metrics
        )
        selected_found.update(metric_ids)
        if not metric_ids:
            limitations.append(
                f"Experiment {record.specification.experiment_id!r} has no selected "
                "requested metrics."
            )
        for metric_id in metric_ids:
            metric_observations = extract_metric_observations((record,), metric_id)
            observations.extend(
                _metric_observation(
                    item,
                    record_index=record_index,
                    equation_ids=equation_ids,
                    identity_method=identity_method,
                    identity_notes=identity_notes,
                )
                for item in metric_observations
            )
            if record.specification.sweep is not None:
                observations.append(
                    _sweep_change_observation(
                        record,
                        record_index=record_index,
                        metric_id=metric_id,
                        equation_ids=equation_ids,
                        identity_method=identity_method,
                        identity_notes=identity_notes,
                    )
                )
        if record.specification.model_id in {
            "damped_harmonic_oscillator",
            "damped_nonlinear_pendulum",
        } and (selected is None or bool(metric_ids)):
            for run_index, run in enumerate(record.runs):
                for descriptor in (
                    "observed_period",
                    "amplitude_change",
                    "finite_window_return_indicator",
                ):
                    observations.append(
                        _temporal_observation(
                            record,
                            run_index,
                            record_index,
                            equation_ids,
                            identity_method,
                            identity_notes,
                            descriptor,
                            run.ode_trajectory,
                        )
                    )
    if selected is not None:
        unmatched = set(selected) - selected_found
        if unmatched:
            limitations.append(
                "Selected metrics were not requested by any input record: "
                + ", ".join(sorted(unmatched))
                + "."
            )
    if not record_tuple:
        limitations.append("No experiment records were supplied.")
    observation_tuple = tuple(observations)
    matches = _structure_matches(record_tuple, observation_tuple, kb)
    hypotheses = _application_hypotheses(record_tuple, observation_tuple, matches, kb)
    if not matches:
        limitations.append(
            "The current registry supplied no sufficiently supported cross-record "
            "structure match for these observations."
        )
    if not hypotheses:
        limitations.append(
            "No hypothesis met the implemented evidence gate; no application "
            "was invented to fill the gap."
        )
    limitation_tuple = tuple(dict.fromkeys(limitations))
    report = _render_discovery_report(
        record_tuple, observation_tuple, matches, hypotheses, limitation_tuple
    )
    return DiscoveryResult(
        observations=observation_tuple,
        structure_matches=matches,
        hypotheses=hypotheses,
        collection_summary=summary,
        limitations=limitation_tuple,
        report_markdown=report,
    )
