"""Knowledge queries linking existing Newton Lab discoveries to equation records.

Retrospective experiment links are kept as immutable association records. The
original experiment specifications, discovery observations, and solver outputs
are not rewritten when a registry association is added.
"""

from __future__ import annotations

from collections import deque
from enum import StrEnum
from typing import Literal

from pydantic import Field, model_validator

from newton_lab.cross_domain import DampingTransferCaseStudy
from newton_lab.discovery import ApplicationHypothesis
from newton_lab.exceptions import RegistryError, ScientificValidationError
from newton_lab.experiments import ExperimentRecord, ExperimentRun
from newton_lab.heat_diffusion_case_study import HeatDiffusionSmoothingCaseStudy
from newton_lab.knowledge.examples import (
    build_example_knowledge_base,
    build_example_registry,
)
from newton_lab.knowledge.registry import (
    EquationRecord,
    EquationRegistry,
    EquationRelation,
    MathematicalStructure,
    RelationshipType,
    SourceReference,
)
from newton_lab.knowledge.workflow import (
    EquationRelationship,
    KnowledgeBase,
    RelationshipReviewStatus,
)
from newton_lab.simulation import SimulationModel


class ExperimentEquationAssociation(SimulationModel):
    """Explicit link to an experiment without changing its historical record."""

    association_id: str = Field(pattern=r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
    experiment_id: str
    experiment_model_id: str
    equation_id: str
    source_run_ids: tuple[str, ...] = Field(min_length=1)
    source_run_snapshots: tuple[ExperimentRun, ...] = Field(min_length=1)
    original_equation_record_ids: tuple[str, ...] = ()
    association_origin: Literal[
        "retrospective_compatibility_checked", "already_declared"
    ]
    association_basis: str = Field(min_length=1)
    case_study_id: str
    limitations: tuple[str, ...] = Field(min_length=1)

    @model_validator(mode="after")
    def run_references_match(self) -> ExperimentEquationAssociation:
        indexes = tuple(run.run_index for run in self.source_run_snapshots)
        if len(indexes) != len(set(indexes)):
            raise ValueError("associated run snapshots must have unique indices")
        expected_ids = tuple(
            f"0:{self.experiment_id}#{run.run_index}"
            for run in self.source_run_snapshots
        )
        if self.source_run_ids != expected_ids:
            raise ValueError("source run IDs must preserve discovery run identity")
        if any(
            run.reproducibility.experiment_id != self.experiment_id
            for run in self.source_run_snapshots
        ):
            raise ValueError("associated run snapshots must match the experiment ID")
        if self.association_origin == "retrospective_compatibility_checked" and (
            self.equation_id in self.original_equation_record_ids
        ):
            raise ValueError(
                "retrospective associations cannot claim prior equation IDs"
            )
        return self


class ObservationEquationProvenance(SimulationModel):
    """Original equation identity data, kept distinct from later associations."""

    observation_id: str
    original_equation_record_ids: tuple[str, ...]
    identity_method: str


class IntegratedDiscoveryLink(SimulationModel):
    """A case-study and its untouched source discovery identifiers."""

    case_study_id: str
    case_study_title: str
    experiment_id: str
    association_id: str
    equation_id: str
    original_experiment_equation_ids: tuple[str, ...]
    observation_equation_provenance: tuple[ObservationEquationProvenance, ...]
    application_hypotheses: tuple[ApplicationHypothesis, ...] = ()
    evidence_classification: str
    mathematical_relationship: str
    assumptions: tuple[str, ...] = Field(min_length=1)
    target_model_equation: str | None = None
    target_domain_status: str
    case_specific_details: tuple[str, ...] = ()
    limitations: tuple[str, ...] = Field(min_length=1)


class KnowledgeIntegrationIndex(SimulationModel):
    """Immutable Phase 18/19 integration built on the project registry and graph."""

    registry: EquationRegistry
    knowledge_base: KnowledgeBase
    experiment_associations: tuple[ExperimentEquationAssociation, ...] = Field(
        min_length=2
    )
    discoveries: tuple[IntegratedDiscoveryLink, ...] = Field(min_length=2)
    report_markdown: str

    @model_validator(mode="after")
    def links_are_consistent(self) -> KnowledgeIntegrationIndex:
        equation_ids = {record.equation_id for record in self.registry.equations}
        registry_by_id = {
            record.equation_id: record for record in self.registry.equations
        }
        knowledge_by_id = {
            entry.record.equation_id: entry.record
            for entry in self.knowledge_base.entries
        }
        if registry_by_id != knowledge_by_id:
            raise ValueError("knowledge workflow entries must match registry records")
        associations_by_id = {
            association.association_id: association
            for association in self.experiment_associations
        }
        if len(associations_by_id) != len(self.experiment_associations):
            raise ValueError("experiment association IDs must be unique")
        association_targets = tuple(
            (item.experiment_id, item.equation_id)
            for item in self.experiment_associations
        )
        if len(association_targets) != len(set(association_targets)):
            raise ValueError("experiment/equation associations must be unique")
        if any(
            association.equation_id not in equation_ids
            for association in self.experiment_associations
        ):
            raise ValueError("experiment association references an unknown equation")
        discovery_ids = {item.case_study_id for item in self.discoveries}
        if len(discovery_ids) != len(self.discoveries):
            raise ValueError("integrated case-study IDs must be unique")
        for discovery in self.discoveries:
            association = associations_by_id.get(discovery.association_id)
            if (
                association is None
                or association.case_study_id != discovery.case_study_id
                or association.experiment_id != discovery.experiment_id
                or association.equation_id != discovery.equation_id
            ):
                raise ValueError("discovery link must match its experiment association")
        return self


class QueryKind(StrEnum):
    EQUATION = "equation"
    DISCOVERY = "discovery"


class StructureReference(SimulationModel):
    equation_id: str
    structure: MathematicalStructure | None


class ResolvedRegistryRelationship(SimulationModel):
    """A registry-stored equation link and the target record it resolves to."""

    source_equation_id: str
    target_equation: EquationRecord
    relationship: EquationRelation
    evidence_classification: Literal["registry_mathematical_derivation"] = (
        "registry_mathematical_derivation"
    )


class DiscoveryKnowledgeQuery(SimulationModel):
    """Equation- or discovery-centered view of the integrated evidence graph."""

    query_kind: QueryKind
    query_id: str
    equation_records: tuple[EquationRecord, ...]
    mathematical_structures: tuple[StructureReference, ...]
    registry_relationships: tuple[ResolvedRegistryRelationship, ...]
    workflow_relationships: tuple[EquationRelationship, ...]
    experiment_associations: tuple[ExperimentEquationAssociation, ...]
    discoveries: tuple[IntegratedDiscoveryLink, ...]
    provenance: tuple[SourceReference, ...]
    evidence_classifications: tuple[str, ...]
    limitations: tuple[str, ...]


_OSCILLATOR_EQUATION_ID = "damped_harmonic_oscillator"
_HEAT_EQUATION_ID = "steady_heat_conduction_1d"
_TRANSIENT_HEAT_ID = "transient_heat_diffusion_1d"
_INTEGRATED_EQUATION_IDS = frozenset(
    {
        _OSCILLATOR_EQUATION_ID,
        _HEAT_EQUATION_ID,
        "steady_heat_conduction_divergence_1d",
        _TRANSIENT_HEAT_ID,
    }
)


def _compatibility_basis(record: ExperimentRecord, equation: EquationRecord) -> str:
    specification = record.specification
    if specification.equation_record_ids and (
        equation.equation_id not in specification.equation_record_ids
    ):
        raise ScientificValidationError(
            "experiment already declares an incompatible equation association"
        )
    declared = equation.equation_id in specification.equation_record_ids
    if declared and any(
        equation.equation_id not in run.reproducibility.equation_record_ids
        for run in record.runs
    ):
        raise ScientificValidationError(
            "experiment specification declares an equation absent from its "
            "original run provenance"
        )
    if specification.model_id == _OSCILLATOR_EQUATION_ID:
        if equation.equation_id != _OSCILLATOR_EQUATION_ID:
            raise ScientificValidationError(
                "oscillator experiments can only associate with the damped oscillator"
            )
        required = {
            "mass_kg": "kg",
            "damping_coefficient_kg_per_s": "kg/s",
            "stiffness_n_per_m": "N/m",
        }
        parameters = {value.name: value.unit for value in specification.parameters}
        if any(parameters.get(name) != unit for name, unit in required.items()):
            raise ScientificValidationError(
                "oscillator experiment parameter units do not match its equation"
            )
        if not any("Linear spring" in text for text in specification.assumptions):
            raise ScientificValidationError(
                "oscillator experiment is missing its linear-spring assumption"
            )
        basis = (
            "The stable experiment model ID matches the registry equation, its "
            "mass/damping/stiffness parameters use the declared SI units, and "
            "the experiment records a linear spring assumption."
        )
        return _association_basis_text(basis, declared)
    if specification.model_id == "boundary_value_problem":
        problem = specification.bvp_problem
        if equation.equation_id != _HEAT_EQUATION_ID:
            raise ScientificValidationError(
                "the existing heat BVP can only associate with the constant-"
                "conductivity steady heat equation"
            )
        if problem is None or problem.problem_id != "steady_linear_heat_conduction":
            raise ScientificValidationError(
                "BVP experiment is not the existing steady heat-conduction problem"
            )
        if problem.equation_record_ids != specification.equation_record_ids:
            raise ScientificValidationError(
                "BVP problem and experiment equation IDs are inconsistent"
            )
        if (
            problem.domain_unit != "m"
            or tuple((item.name, item.unit) for item in problem.states)
            != (("temperature", "K"), ("temperature_gradient", "K/m"))
            or "No internal heat generation." not in problem.assumptions
            or not any(
                "constant thermal conductivity" in text for text in problem.assumptions
            )
        ):
            raise ScientificValidationError(
                "heat BVP units or assumptions are incompatible with the equation"
            )
        basis = (
            "The existing heat BVP problem ID, spatial/temperature units, and "
            "recorded constant-conductivity, source-free assumptions match the "
            "reduced registry equation."
        )
        return _association_basis_text(basis, declared)
    raise ScientificValidationError(
        f"no compatibility rule exists for model {specification.model_id!r}"
    )


def _association_basis_text(basis: str, declared: bool) -> str:
    if declared:
        return (
            "Already declared association validated against model compatibility. "
            f"{basis}"
        )
    return f"Added retrospectively; source records remain unchanged. {basis}"


def create_experiment_equation_association(
    record: ExperimentRecord,
    equation: EquationRecord,
    *,
    case_study_id: str,
    limitations: tuple[str, ...],
) -> ExperimentEquationAssociation:
    """Create a retrospective or already-declared association without mutation."""
    basis = _compatibility_basis(record, equation)
    existing = equation.equation_id in record.specification.equation_record_ids
    return ExperimentEquationAssociation(
        association_id=f"{record.specification.experiment_id}_{equation.equation_id}_association",
        experiment_id=record.specification.experiment_id,
        experiment_model_id=record.specification.model_id,
        equation_id=equation.equation_id,
        source_run_ids=tuple(
            f"0:{record.specification.experiment_id}#{run.run_index}"
            for run in record.runs
        ),
        source_run_snapshots=record.runs,
        original_equation_record_ids=record.specification.equation_record_ids,
        association_origin=(
            "already_declared" if existing else "retrospective_compatibility_checked"
        ),
        association_basis=basis,
        case_study_id=case_study_id,
        limitations=limitations,
    )


def _oscillator_diffusion_candidate() -> EquationRelationship:
    return EquationRelationship(
        relationship_id="candidate_oscillator_heat_diffusion_shared_dynamics",
        source_equation_id=_OSCILLATOR_EQUATION_ID,
        target_equation_id=_TRANSIENT_HEAT_ID,
        relationship_type=RelationshipType.SHARES_MATHEMATICAL_STRUCTURE,
        description=(
            "Candidate shared structure between linear deterministic evolution "
            "equations with exponential solution components."
        ),
        justification=(
            "The registry explicitly classifies both equations as linear and "
            "deterministic. Oscillator characteristic roots can yield oscillatory "
            "or non-oscillatory responses depending on damping; diffusion modes "
            "decay without oscillation at rates set by spatial frequency."
        ),
        conditions=(
            "This is a mathematical-structure comparison only.",
            "No physical equivalence, common mechanism, or common units are asserted.",
            "The candidate remains subject to human review.",
        ),
        directional=False,
        review_status=RelationshipReviewStatus.CANDIDATE,
    )


def build_knowledge_integration(
    oscillator_case: DampingTransferCaseStudy,
    heat_case: HeatDiffusionSmoothingCaseStudy,
    *,
    registry: EquationRegistry | None = None,
) -> KnowledgeIntegrationIndex:
    """Link Phase 18/19 discoveries to existing and new registry identities.

    Source records and case-study discoveries are read only. Existing equation
    IDs are reused, and the oscillator/diffusion structural edge remains a
    review-required candidate.
    """
    equations = registry or build_example_registry()
    for equation_id in _INTEGRATED_EQUATION_IDS:
        try:
            equations.get(equation_id)
        except RegistryError as error:
            raise ScientificValidationError(
                f"integration registry is missing required equation {equation_id!r}"
            ) from error

    base_knowledge = (
        build_example_knowledge_base()
        if registry is None
        else KnowledgeBase.from_equation_registry(equations)
    )
    has_candidate = any(
        {edge.source_equation_id, edge.target_equation_id}
        == {_OSCILLATOR_EQUATION_ID, _TRANSIENT_HEAT_ID}
        and edge.relationship_type == RelationshipType.SHARES_MATHEMATICAL_STRUCTURE
        for edge in base_knowledge.relationships
    )
    if not has_candidate:
        base_knowledge = base_knowledge.add_relationship(
            _oscillator_diffusion_candidate()
        )

    oscillator_equation = equations.get(_OSCILLATOR_EQUATION_ID)
    heat_equation = equations.get(_HEAT_EQUATION_ID)
    oscillator_association = create_experiment_equation_association(
        oscillator_case.source_experiment,
        oscillator_equation,
        case_study_id=oscillator_case.case_study_id,
        limitations=oscillator_case.limitations,
    )
    heat_association = create_experiment_equation_association(
        heat_case.source_experiment,
        heat_equation,
        case_study_id=heat_case.case_study_id,
        limitations=heat_case.limitations,
    )
    oscillator_hypothesis = oscillator_case.source_hypothesis
    oscillator_observation_links = tuple(
        ObservationEquationProvenance(
            observation_id=observation_id,
            original_equation_record_ids=(oscillator_hypothesis.source_equation_id,),
            identity_method=(
                "preserved Phase 18 discovery hypothesis source equation identity"
            ),
        )
        for observation_id in oscillator_hypothesis.source_behavior_observation_ids
    )
    heat_observation_links = tuple(
        ObservationEquationProvenance(
            observation_id=observation.observation_id,
            original_equation_record_ids=observation.equation_record_ids,
            identity_method=observation.equation_identity_method,
        )
        for observation in heat_case.discovery.observations
    )
    discoveries = (
        IntegratedDiscoveryLink(
            case_study_id=oscillator_case.case_study_id,
            case_study_title="Oscillator-to-adjustment model case study",
            experiment_id=oscillator_case.source_experiment_id,
            association_id=oscillator_association.association_id,
            equation_id=oscillator_association.equation_id,
            original_experiment_equation_ids=(
                *oscillator_case.source_experiment.specification.equation_record_ids,
            ),
            observation_equation_provenance=oscillator_observation_links,
            application_hypotheses=(oscillator_hypothesis,),
            evidence_classification="proposed_application_hypothesis_unvalidated",
            mathematical_relationship=oscillator_case.mathematical_correspondence,
            assumptions=oscillator_case.assumptions,
            target_model_equation=oscillator_case.target_equation,
            target_domain_status=(
                f"{oscillator_case.target_domain_interpretation}; "
                f"real_world_application={oscillator_case.real_world_application_status}"
            ),
            case_specific_details=(
                f"State scale x=Aq with A={oscillator_case.amplitude_scale_m:g} m "
                "per dimensionless state unit.",
                f"Time constant tau={oscillator_case.time_constant_s:g} s; the "
                "parameter map is tau=sqrt(m/k), zeta=c/(2*sqrt(m*k)).",
            ),
            limitations=oscillator_case.limitations,
        ),
        IntegratedDiscoveryLink(
            case_study_id=heat_case.case_study_id,
            case_study_title="Heat conduction, diffusion, and smoothing case study",
            experiment_id=heat_case.source_experiment.specification.experiment_id,
            association_id=heat_association.association_id,
            equation_id=heat_association.equation_id,
            original_experiment_equation_ids=(
                *heat_case.source_experiment.specification.equation_record_ids,
            ),
            observation_equation_provenance=heat_observation_links,
            evidence_classification="analytical_model_result_and_synthetic_analogy",
            mathematical_relationship=heat_case.mathematical_relationship,
            assumptions=(
                "One-dimensional homogeneous continuum.",
                "Constant positive thermal diffusivity and fixed endpoint "
                "temperatures.",
                "No internal heat generation; the transient perturbation uses "
                "homogeneous Dirichlet boundaries.",
            ),
            target_model_equation=(
                "Synthetic signal: s(x,t)=sum(b_n sin(n*pi*x/L) "
                "exp(-alpha*(n*pi/L)^2*t)); arbitrary signal units."
            ),
            target_domain_status="synthetic_structural_analogy_only; unvalidated",
            case_specific_details=(
                f"L={heat_case.diffusion_model.length_m:g} m; "
                f"alpha={heat_case.diffusion_model.diffusivity_m2_per_s:g} m^2/s.",
                "The sample grid evaluates the continuous modal solution; it is "
                "not a discrete filter or market data.",
            ),
            limitations=heat_case.limitations,
        ),
    )
    shell = KnowledgeIntegrationIndex.model_construct(
        registry=equations,
        knowledge_base=base_knowledge,
        experiment_associations=(oscillator_association, heat_association),
        discoveries=discoveries,
        report_markdown="",
    )
    return KnowledgeIntegrationIndex(
        registry=equations,
        knowledge_base=base_knowledge,
        experiment_associations=(oscillator_association, heat_association),
        discoveries=discoveries,
        report_markdown=render_knowledge_integration_report(shell),
    )


def query_discovery_knowledge(
    index: KnowledgeIntegrationIndex, identifier: str
) -> DiscoveryKnowledgeQuery:
    """Query the connected equation graph by equation ID or case-study ID."""
    equation_ids = {record.equation_id for record in index.registry.equations}
    if identifier in equation_ids:
        query_kind = QueryKind.EQUATION
        roots = {identifier}
    else:
        linked = tuple(
            discovery
            for discovery in index.discoveries
            if discovery.case_study_id == identifier
        )
        if not linked:
            raise RegistryError(f"unknown equation or discovery ID: {identifier}")
        query_kind = QueryKind.DISCOVERY
        roots = {item.equation_id for item in linked}

    registry_edges = {
        (record.equation_id, relation.target_equation_id): relation
        for record in index.registry.equations
        for relation in record.related_equations
    }
    workflow_edges = index.knowledge_base.relationships
    # Traverse the declared registry graph, plus reviewed/candidate workflow
    # edges, without deriving any edge from names, units, or expression text.
    visited = set(roots)
    frontier = deque(roots)
    while frontier:
        current = frontier.popleft()
        for source_id, target_id in registry_edges:
            if current in {source_id, target_id}:
                other = target_id if current == source_id else source_id
                if other not in visited:
                    visited.add(other)
                    frontier.append(other)
        for edge in workflow_edges:
            if current in {edge.source_equation_id, edge.target_equation_id}:
                other = (
                    edge.target_equation_id
                    if current == edge.source_equation_id
                    else edge.source_equation_id
                )
                if other not in visited:
                    visited.add(other)
                    frontier.append(other)

    records = tuple(
        index.registry.get(equation_id)
        for equation_id in sorted(visited)
        if equation_id in equation_ids
    )
    resolved_relations: list[ResolvedRegistryRelationship] = []
    for record in records:
        for relation in record.related_equations:
            if relation.target_equation_id in visited:
                resolved_relations.append(
                    ResolvedRegistryRelationship(
                        source_equation_id=record.equation_id,
                        target_equation=index.registry.get(relation.target_equation_id),
                        relationship=relation,
                    )
                )
    related_workflow = tuple(
        edge
        for edge in workflow_edges
        if edge.source_equation_id in visited or edge.target_equation_id in visited
    )
    associations = tuple(
        association
        for association in index.experiment_associations
        if association.equation_id in visited
    )
    association_ids = {association.association_id for association in associations}
    discoveries = tuple(
        item for item in index.discoveries if item.association_id in association_ids
    )
    sources_by_id = {
        source.source_id: source
        for record in records
        for source in record.source_references
    }
    classifications = {"registry_equation_identity_and_structure"}
    if resolved_relations:
        classifications.add("registry_mathematical_relationship_with_conditions")
    if related_workflow:
        classifications.add("workflow_relationship_requires_review_or_is_reviewed")
    if associations:
        classifications.add(
            "experiment_equation_association_includes_retrospective_links"
        )
    if any(item.application_hypotheses for item in discoveries):
        classifications.add("cross_domain_application_hypothesis_unvalidated")
    limitations = tuple(
        dict.fromkeys(
            limitation
            for discovery in discoveries
            for limitation in discovery.limitations
        )
    )
    return DiscoveryKnowledgeQuery(
        query_kind=query_kind,
        query_id=identifier,
        equation_records=records,
        mathematical_structures=tuple(
            StructureReference(
                equation_id=record.equation_id,
                structure=record.mathematical_structure,
            )
            for record in records
        ),
        registry_relationships=tuple(resolved_relations),
        workflow_relationships=related_workflow,
        experiment_associations=associations,
        discoveries=discoveries,
        provenance=tuple(sources_by_id[key] for key in sorted(sources_by_id)),
        evidence_classifications=tuple(sorted(classifications)),
        limitations=limitations,
    )


def render_knowledge_integration_report(index: KnowledgeIntegrationIndex) -> str:
    """Render deterministic Markdown distinguishing registry and candidate links."""
    lines = [
        "# Newton Lab Phase 20 — Discovery knowledge integration",
        "",
        "This report links existing Phase 18 and Phase 19 records to the equation "
        "registry. Retrospective associations do not rewrite source experiment "
        "records, observations, or numerical results.",
        "",
        "## Equation records",
        "",
        "| Equation ID | Expression | Classifications | Source provenance |",
        "|---|---|---|---|",
    ]
    for equation in index.registry.equations:
        if equation.equation_id not in _INTEGRATED_EQUATION_IDS:
            continue
        lines.append(
            f"| `{equation.equation_id}` | {equation.mathematical_expression} | "
            f"{', '.join(item.value for item in equation.mathematical_classifications)}"
            " | "
            f"{', '.join(source.source_id for source in equation.source_references)} |"
        )
    sources = {
        source.source_id: source
        for equation in index.registry.equations
        if equation.equation_id in _INTEGRATED_EQUATION_IDS
        for source in equation.source_references
    }
    for source in sorted(sources.values(), key=lambda item: item.source_id):
        lines.append(
            f"- Provenance `{source.source_id}`: {source.title}, "
            f"{source.responsible_organization}, {source.year}; {source.url}."
        )
    lines.extend(["", "## Registry-backed conditional relationships", ""])
    for equation in index.registry.equations:
        for relation in equation.related_equations:
            if relation.target_equation_id not in {
                "steady_heat_conduction_divergence_1d",
                _HEAT_EQUATION_ID,
            }:
                continue
            lines.extend(
                [
                    f"### `{equation.equation_id}` → `{relation.target_equation_id}`",
                    "",
                    f"- Type: `{relation.relationship.value}`",
                    f"- Evidence: `{relation.evidence.category.value}`",
                    f"- Conditions and rationale: {relation.evidence.statement}",
                    "",
                ]
            )
    lines.extend(
        [
            "The transient-to-steady link sets the time derivative to zero and "
            "requires constant positive diffusivity, constant material properties, "
            "one-dimensional source-free conduction, and matching boundary "
            "conditions. Fixed endpoint temperatures select a linear profile. "
            "This is a conditional mathematical reduction, not identity of the "
            "transient and steady models.",
            "",
            "## Shared oscillator/diffusion structure candidate",
            "",
        ]
    )
    for edge in index.knowledge_base.relationships:
        if {edge.source_equation_id, edge.target_equation_id} != {
            _OSCILLATOR_EQUATION_ID,
            _TRANSIENT_HEAT_ID,
        }:
            continue
        lines.extend(
            [
                f"- `{edge.source_equation_id}` ↔ `{edge.target_equation_id}`: "
                f"{edge.description} Status: `{edge.review_status.value}`.",
                f"  Rationale: {edge.justification}",
                f"  Conditions: {'; '.join(edge.conditions)}",
            ]
        )
    lines.extend(
        [
            "",
            "The oscillator may be underdamped, critically damped, or overdamped; "
            "diffusion modes decay without oscillation, with rates proportional "
            "to spatial frequency squared. This candidate asserts no physical "
            "equivalence and requires review.",
            "",
            "## Retrospective experiment associations",
            "",
            "| Experiment | Equation ID | Association origin | "
            "Original equation IDs | Runs retained |",
            "|---|---|---|---|---:|",
        ]
    )
    for association in index.experiment_associations:
        lines.append(
            f"| `{association.experiment_id}` | `{association.equation_id}` | "
            f"`{association.association_origin}` | "
            f"{association.original_equation_record_ids or 'none'} | "
            f"{len(association.source_run_snapshots)} |"
        )
        lines.append(f"\n{association.association_basis}\n")
        for run in association.source_run_snapshots:
            solver_evidence = (
                f"solver_success={run.solver_success}; "
                f"numerical_acceptance={run.numerical_acceptance}"
            )
            diagnostics = run.bvp_diagnostics
            if diagnostics is not None:
                solver_evidence += (
                    f"; mesh_points={len(diagnostics.mesh)}; "
                    f"boundary_residuals={diagnostics.boundary_residuals}; "
                    f"max_differential_residual="
                    f"{max(diagnostics.differential_residuals, default=0.0):.6g}"
                )
            lines.append(
                f"- Preserved run `0:{association.experiment_id}#{run.run_index}`: "
                f"status=`{run.status.value}`, {solver_evidence}."
            )
    lines.extend(["", "## Linked discoveries and evidence status", ""])
    for discovery in index.discoveries:
        lines.extend(
            [
                f"### {discovery.case_study_title}",
                "",
                f"- Case: `{discovery.case_study_id}`; experiment: "
                f"`{discovery.experiment_id}`; equation: `{discovery.equation_id}`.",
                f"- Original experiment equation IDs: "
                f"{discovery.original_experiment_equation_ids or 'none'}.",
                f"- Original observation count: "
                f"{len(discovery.observation_equation_provenance)}.",
                f"- Evidence: `{discovery.evidence_classification}`.",
                "- Mathematical relationship: "
                f"`{discovery.mathematical_relationship}`.",
                f"- Target interpretation: `{discovery.target_domain_status}`.",
            ]
        )
        if discovery.target_model_equation is not None:
            lines.append(f"- Target model: `{discovery.target_model_equation}`.")
        lines.extend(f"- Assumption: {item}" for item in discovery.assumptions)
        lines.extend(
            f"- Case detail: {item}" for item in discovery.case_specific_details
        )
        if discovery.application_hypotheses:
            for hypothesis in discovery.application_hypotheses:
                lines.extend(
                    [
                        f"- Application hypothesis: {hypothesis.statement}",
                        f"  Status: `{hypothesis.evidence_status.value}`; "
                        f"validated_application={hypothesis.validated_application}.",
                        f"  Proposed mechanism: {hypothesis.proposed_mechanism}",
                        "  Similarities: " + " ".join(hypothesis.similarities),
                        "  Differences: " + " ".join(hypothesis.differences),
                    ]
                )
                lines.extend(f"  Assumption: {item}" for item in hypothesis.assumptions)
                lines.extend(f"  Test plan: {item}" for item in hypothesis.test_plan)
                if discovery.target_model_equation is not None:
                    lines.append(
                        f"- Target model: `{discovery.target_model_equation}`."
                    )
        else:
            lines.extend(f"- Assumption: {item}" for item in discovery.assumptions)
            lines.extend(
                f"- Case detail: {item}" for item in discovery.case_specific_details
            )
            if discovery.target_model_equation is not None:
                lines.append(f"- Target model: `{discovery.target_model_equation}`.")
        lines.extend(["- Limitations: " + " ".join(discovery.limitations), ""])
        for observation in discovery.observation_equation_provenance:
            lines.append(
                f"  - Observation `{observation.observation_id}` originally "
                f"recorded equation IDs: "
                f"{observation.original_equation_record_ids or 'none'} "
                f"({observation.identity_method})."
            )
        lines.append("")
    lines.extend(
        [
            "## Query example",
            "",
            "```python",
            "from newton_lab.knowledge.integration import (",
            "    build_knowledge_integration,",
            "    query_discovery_knowledge,",
            ")",
            "from newton_lab.cross_domain import run_damping_transfer_case_study",
            "from newton_lab.heat_diffusion_case_study import (",
            "    run_heat_diffusion_smoothing_case_study,",
            ")",
            "",
            "index = build_knowledge_integration(",
            "    run_damping_transfer_case_study(),",
            "    run_heat_diffusion_smoothing_case_study(),",
            ")",
            "query = query_discovery_knowledge(index, 'transient_heat_diffusion_1d')",
            "print(query.registry_relationships)",
            "print(query.workflow_relationships)  # candidate status is retained",
            "```",
            "",
            "A query by case-study ID returns the same connected equation and "
            "provenance view. Equation text is descriptive, not executable. "
            "Registry equations and conditional relations are recorded knowledge; "
            "the oscillator/diffusion edge is only a candidate; Phase 18's "
            "application hypothesis remains unvalidated. Numerical acceptance "
            "does not establish physical validation.",
            "",
        ]
    )
    return "\n".join(lines)
