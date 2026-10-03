"""Phase 20 equation, association, query, and evidence-boundary tests."""

from __future__ import annotations

import pytest

from newton_lab.cross_domain import (
    DampingTransferCaseStudy,
    run_damping_transfer_case_study,
)
from newton_lab.exceptions import RegistryError, ScientificValidationError
from newton_lab.heat_diffusion_case_study import (
    HeatDiffusionSmoothingCaseStudy,
    run_heat_diffusion_smoothing_case_study,
)
from newton_lab.knowledge.examples import build_example_registry
from newton_lab.knowledge.integration import (
    KnowledgeIntegrationIndex,
    build_knowledge_integration,
    create_experiment_equation_association,
    query_discovery_knowledge,
)
from newton_lab.knowledge.registry import RelationshipType
from newton_lab.knowledge.workflow import RelationshipReviewStatus


@pytest.fixture(scope="module")
def phase20() -> tuple[
    KnowledgeIntegrationIndex,
    DampingTransferCaseStudy,
    HeatDiffusionSmoothingCaseStudy,
]:
    oscillator = run_damping_transfer_case_study(output_point_count=201)
    heat = run_heat_diffusion_smoothing_case_study(sample_count=21)
    return build_knowledge_integration(oscillator, heat), oscillator, heat


def test_registry_distinguishes_general_steady_reduced_and_transient_equations() -> (
    None
):
    registry = build_example_registry()
    general = registry.get("steady_heat_conduction_divergence_1d")
    steady = registry.get("steady_heat_conduction_1d")
    transient = registry.get("transient_heat_diffusion_1d")

    assert general.mathematical_expression == "d/dx(k(x) dT/dx) = 0"
    assert steady.mathematical_expression == "d^2T/dx^2 = 0"
    assert transient.mathematical_expression == (
        "partial T/partial t = alpha partial^2 T/partial x^2"
    )
    assert "constant" not in general.parameters[0].meaning.casefold()
    assert steady.parameters[0].si_unit == "W/(m K)"
    alpha = next(item for item in transient.parameters if item.symbol == "alpha")
    assert alpha.si_unit == "m^2/s"
    assert "k/(rho cp)" in alpha.meaning
    assert transient.mathematical_structure is not None
    assert transient.mathematical_structure.pde_order == 2
    assert transient.mathematical_structure.linearity is not None
    assert transient.mathematical_structure.linearity.value == "linear"
    assert general.source_references[0].source_id == (
        "mit_ocw_intro_engineering_heat_transfer_2002"
    )
    assert "ocw.mit.edu" in (general.source_references[0].url or "")
    assert any(
        item.source_id == "mit_ocw_1d_thermal_diffusion_2003"
        for item in transient.source_references
    )


def test_registry_links_encode_conditional_reductions() -> None:
    registry = build_example_registry()
    steady = registry.get("steady_heat_conduction_1d")
    transient = registry.get("transient_heat_diffusion_1d")

    reduction = steady.related_equations[0]
    equilibrium = transient.related_equations[0]
    assert reduction.target_equation_id == "steady_heat_conduction_divergence_1d"
    assert reduction.relationship == RelationshipType.SPECIAL_CASE_OF
    assert "constant k" in reduction.evidence.statement
    assert equilibrium.target_equation_id == steady.equation_id
    assert equilibrium.relationship == RelationshipType.LIMIT_OF
    assert "partial T/partial t=0" in equilibrium.evidence.statement
    assert "Fixed endpoint temperatures" in equilibrium.evidence.statement
    assert all(
        relation.relationship != RelationshipType.EQUIVALENT_UNDER_ASSUMPTIONS
        for equation in registry.equations
        for relation in equation.related_equations
        if relation.target_equation_id
        in {
            "steady_heat_conduction_1d",
            "steady_heat_conduction_divergence_1d",
            "transient_heat_diffusion_1d",
        }
    )
    assert len({equation.equation_id for equation in registry.equations}) == len(
        registry.equations
    )


def test_integrated_associations_preserve_historical_records_and_run_evidence(
    phase20: tuple[
        KnowledgeIntegrationIndex,
        DampingTransferCaseStudy,
        HeatDiffusionSmoothingCaseStudy,
    ],
) -> None:
    index, oscillator, heat = phase20
    associations = {item.experiment_id: item for item in index.experiment_associations}
    heat_association = associations[heat.source_experiment.specification.experiment_id]
    oscillator_association = associations[oscillator.source_experiment_id]

    assert heat_association.equation_id == "steady_heat_conduction_1d"
    assert oscillator_association.equation_id == "damped_harmonic_oscillator"
    assert all(
        item.association_origin == "retrospective_compatibility_checked"
        for item in associations.values()
    )
    assert heat_association.original_equation_record_ids == ()
    assert heat.source_experiment.specification.equation_record_ids == ()
    assert heat.source_experiment.runs == heat_association.source_run_snapshots
    assert (
        oscillator.source_experiment.runs == oscillator_association.source_run_snapshots
    )
    diagnostics = heat_association.source_run_snapshots[0].bvp_diagnostics
    source_diagnostics = heat.source_experiment.runs[0].bvp_diagnostics
    assert diagnostics is not None
    assert source_diagnostics is not None
    assert diagnostics.mesh == source_diagnostics.mesh
    assert diagnostics.boundary_residuals == (source_diagnostics.boundary_residuals)
    assert diagnostics.numerical_accepted
    assert heat_association.source_run_ids == ("0:phase19_steady_heat_source#0",)
    assert heat_association.association_basis


def test_incompatible_retrospective_association_is_rejected(
    phase20: tuple[
        KnowledgeIntegrationIndex,
        DampingTransferCaseStudy,
        HeatDiffusionSmoothingCaseStudy,
    ],
) -> None:
    index, _, heat = phase20

    with pytest.raises(ScientificValidationError, match="heat BVP"):
        create_experiment_equation_association(
            heat.source_experiment,
            index.registry.get("damped_harmonic_oscillator"),
            case_study_id=heat.case_study_id,
            limitations=heat.limitations,
        )

    conflicting_specification = heat.source_experiment.specification.model_copy(
        update={"equation_record_ids": ("damped_harmonic_oscillator",)}
    )
    conflicting_record = heat.source_experiment.model_copy(
        update={"specification": conflicting_specification}
    )
    with pytest.raises(ScientificValidationError, match="incompatible equation"):
        create_experiment_equation_association(
            conflicting_record,
            index.registry.get("steady_heat_conduction_1d"),
            case_study_id=heat.case_study_id,
            limitations=heat.limitations,
        )


def test_query_returns_equations_structures_relationships_sources_and_discoveries(
    phase20: tuple[
        KnowledgeIntegrationIndex,
        DampingTransferCaseStudy,
        HeatDiffusionSmoothingCaseStudy,
    ],
) -> None:
    index, oscillator, heat = phase20
    query = query_discovery_knowledge(index, "transient_heat_diffusion_1d")
    equation_ids = {item.equation_id for item in query.equation_records}

    assert query.query_id == "transient_heat_diffusion_1d"
    assert "steady_heat_conduction_1d" in equation_ids
    assert "steady_heat_conduction_divergence_1d" in equation_ids
    assert any(
        item.relationship.relationship == RelationshipType.LIMIT_OF
        for item in query.registry_relationships
    )
    assert any(item.structure is not None for item in query.mathematical_structures)
    assert any(
        item.source_id == "mit_ocw_1d_thermal_diffusion_2003"
        for item in query.provenance
    )
    assert {item.case_study_id for item in query.discoveries} == {
        oscillator.case_study_id,
        heat.case_study_id,
    }
    assert any(
        item.review_status == RelationshipReviewStatus.CANDIDATE
        and item.relationship_type == RelationshipType.SHARES_MATHEMATICAL_STRUCTURE
        for item in query.workflow_relationships
    )

    by_discovery = query_discovery_knowledge(index, heat.case_study_id)
    assert by_discovery.query_id == heat.case_study_id
    heat_link = next(
        item
        for item in by_discovery.discoveries
        if item.case_study_id == heat.case_study_id
    )
    assert heat_link.original_experiment_equation_ids == ()
    assert all(
        not item.original_equation_record_ids
        for item in heat_link.observation_equation_provenance
    )


def test_oscillator_analogy_stays_candidate_and_financial_hypothesis_unvalidated(
    phase20: tuple[
        KnowledgeIntegrationIndex,
        DampingTransferCaseStudy,
        HeatDiffusionSmoothingCaseStudy,
    ],
) -> None:
    index, oscillator, _ = phase20
    query = query_discovery_knowledge(index, oscillator.case_study_id)
    candidate = next(
        item
        for item in query.workflow_relationships
        if item.relationship_id == "candidate_oscillator_heat_diffusion_shared_dynamics"
    )
    linked_hypothesis = next(
        hypothesis
        for discovery in query.discoveries
        for hypothesis in discovery.application_hypotheses
    )

    assert candidate.review_status == RelationshipReviewStatus.CANDIDATE
    assert "No physical equivalence" in " ".join(candidate.conditions)
    assert linked_hypothesis.source_equation_id == "damped_harmonic_oscillator"
    assert linked_hypothesis.validated_application is False
    assert "quantitative_finance" in linked_hypothesis.target_domain


def test_report_distinguishes_relationship_types_and_preserves_original_evidence(
    phase20: tuple[
        KnowledgeIntegrationIndex,
        DampingTransferCaseStudy,
        HeatDiffusionSmoothingCaseStudy,
    ],
) -> None:
    index, _, _ = phase20
    report = index.report_markdown

    assert "registry-backed conditional relationships" in report.casefold()
    assert "candidate_requires_review" in report
    assert "retrospective_compatibility_checked" in report
    assert "validated_application=False" in report
    assert "original observation" in report.casefold()
    assert "no physical equivalence" in report.casefold()
    assert report.count("### Oscillator-to-adjustment model case study") == 1
    assert report.count("### Heat conduction, diffusion, and smoothing case study") == 1
    assert "tau=sqrt(m/k)" in report
    assert "alpha=k/(rho cp)" in report
    restored = KnowledgeIntegrationIndex.model_validate_json(index.model_dump_json())
    assert restored == index


def test_unknown_query_and_incomplete_registry_fail_clearly(
    phase20: tuple[
        KnowledgeIntegrationIndex,
        DampingTransferCaseStudy,
        HeatDiffusionSmoothingCaseStudy,
    ],
) -> None:
    index, _, _ = phase20
    with pytest.raises(RegistryError, match="unknown equation or discovery"):
        query_discovery_knowledge(index, "not_a_known_equation")
    incomplete = build_example_registry().model_copy(
        update={
            "equations": tuple(
                record
                for record in build_example_registry().equations
                if record.equation_id != "transient_heat_diffusion_1d"
            )
        }
    )
    with pytest.raises(ScientificValidationError, match="missing required equation"):
        build_knowledge_integration(
            run_damping_transfer_case_study(output_point_count=201),
            run_heat_diffusion_smoothing_case_study(sample_count=21),
            registry=incomplete,
        )
