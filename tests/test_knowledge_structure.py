"""Tests for equation structure metadata, comparisons, and relationship proposals."""

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from newton_lab.knowledge import (
    EquationRelationship,
    Linearity,
    MathematicalStructure,
    RelationshipReviewStatus,
    RelationshipType,
    build_example_knowledge_base,
    compare_equations,
    propose_relationships,
)
from newton_lab.knowledge.registry import EquationRecord
from newton_lab.knowledge.storage import JsonKnowledgeStore
from newton_lab.knowledge.workflow import KnowledgeBase


def _example_records() -> tuple[KnowledgeBase, EquationRecord, EquationRecord]:
    workspace = build_example_knowledge_base()
    return (
        workspace,
        workspace.get_entry("damped_harmonic_oscillator").record,
        workspace.get_entry("damped_nonlinear_pendulum").record,
    )


def _edge(
    *,
    relationship_id: str = "edge_one",
    source: str = "damped_harmonic_oscillator",
    target: str = "damped_nonlinear_pendulum",
    relationship_type: RelationshipType = RelationshipType.APPROXIMATION_OF,
    directional: bool = True,
    source_reference_ids: tuple[str, ...] = (),
) -> EquationRelationship:
    return EquationRelationship(
        relationship_id=relationship_id,
        source_equation_id=source,
        target_equation_id=target,
        relationship_type=relationship_type,
        description="A relationship to review.",
        justification="This edge is a candidate only.",
        conditions=("Small angle approximation",),
        source_reference_ids=source_reference_ids,
        directional=directional,
    )


def test_structure_features_validate_values_and_unknowns_are_explicit() -> None:
    structure = MathematicalStructure(
        ode_order=2,
        linearity=Linearity.NONLINEAR,
        autonomous=None,
        derivative_operators=("time_derivative",),
        forcing_terms=(),
    )

    assert structure.autonomous is None
    assert structure.forcing_terms == ()
    with pytest.raises(ValidationError, match="must not contain duplicates"):
        MathematicalStructure(derivative_operators=("d_dt", "d_dt"))
    with pytest.raises(ValidationError, match="both ODE and PDE order"):
        MathematicalStructure(ode_order=2, pde_order=2)
    with pytest.raises(ValidationError):
        MathematicalStructure(ode_order=0)


def test_relationship_type_directionality_and_self_links_are_validated() -> None:
    edge = _edge()
    assert edge.relationship_type == RelationshipType.APPROXIMATION_OF
    with pytest.raises(ValidationError, match="directionality must be True"):
        _edge(directional=False)
    with pytest.raises(ValidationError, match="cannot link an equation to itself"):
        _edge(source="damped_harmonic_oscillator", target="damped_harmonic_oscillator")


def test_relationship_review_requires_reviewer_and_does_not_auto_promote() -> None:
    edge = _edge()
    with pytest.raises(ValidationError, match="reviewer and review notes"):
        edge.record_review(
            RelationshipReviewStatus.ACCEPTED,
            reviewer=" ",
            notes="Reviewed",
        )
    accepted = edge.record_review(
        RelationshipReviewStatus.ACCEPTED,
        reviewer="reviewer:physics",
        notes="Conditions checked.",
    )
    assert edge.review_status == RelationshipReviewStatus.CANDIDATE
    assert accepted.review_status == RelationshipReviewStatus.ACCEPTED


def test_workspace_validates_equation_references_source_ids_and_duplicate_edges() -> (
    None
):
    workspace, _, _ = _example_records()

    with pytest.raises(ValidationError, match="unknown equation IDs"):
        workspace.add_relationship(_edge(target="missing_equation"))
    with pytest.raises(ValidationError, match="unknown source IDs"):
        workspace.add_relationship(_edge(source_reference_ids=("missing_source",)))
    with pytest.raises(ValidationError, match="duplicate relationship identity"):
        workspace.add_relationship(_edge(relationship_id="another_approximation"))


def test_undirected_identity_deduplicates_endpoint_order() -> None:
    workspace, _, _ = _example_records()
    edge = _edge(
        relationship_type=RelationshipType.SHARES_MATHEMATICAL_OPERATOR,
        directional=False,
    )
    workspace = workspace.add_relationship(edge)
    with pytest.raises(ValidationError, match="duplicate relationship identity"):
        workspace.add_relationship(
            _edge(
                relationship_id="reversed_edge",
                source="damped_nonlinear_pendulum",
                target="damped_harmonic_oscillator",
                relationship_type=RelationshipType.SHARES_MATHEMATICAL_OPERATOR,
                directional=False,
            )
        )


def test_oscillator_pendulum_comparison_is_structural_not_equivalence() -> None:
    workspace, oscillator, pendulum = _example_records()
    comparison = compare_equations(oscillator, pendulum, knowledge_base=workspace)
    features = {item.feature_name: item for item in comparison.feature_comparisons}

    assert features["ode_order"].status.value == "shared"
    assert features["linearity"].status.value == "different"
    assert "ordinary_differential_equation" in features["classifications"].shared_values
    assert comparison.physical_equivalence_concluded is False
    assert comparison.documented_relationships[0].review_status == (
        RelationshipReviewStatus.CANDIDATE
    )
    assert "damped oscillator" not in comparison.shared_assumptions


def test_unknown_features_are_reported_and_not_guessed_from_expression_text() -> None:
    workspace, oscillator, _ = _example_records()
    unknown_record = oscillator.model_copy(
        update={
            "equation_id": "unknown_structure_example",
            "mathematical_structure": None,
            "mathematical_classifications": (),
        }
    )
    comparison = compare_equations(oscillator, unknown_record)

    assert "ode_order" in comparison.unknown_features
    assert "classifications" in comparison.unknown_features


def test_discovery_generates_review_candidates_with_reasoning_only() -> None:
    workspace, _, _ = _example_records()
    proposed = propose_relationships(workspace)

    assert proposed
    assert all(
        item.review_status == RelationshipReviewStatus.CANDIDATE for item in proposed
    )
    assert all(item.justification for item in proposed)
    assert all(not item.source_reference_ids for item in proposed)
    with_review = workspace.add_relationship(proposed[0])
    assert with_review.relationships[-1].review_status == (
        RelationshipReviewStatus.CANDIDATE
    )


def test_phase_five_json_document_loads_without_new_optional_fields(
    tmp_path: Path,
) -> None:
    workspace, _, _ = _example_records()
    legacy_payload = json.loads(workspace.model_dump_json())
    legacy_payload.pop("relationships")
    for entry in legacy_payload["entries"]:
        entry["record"].pop("mathematical_structure")
    path = tmp_path / "legacy-v1.json"
    path.write_text(
        json.dumps(
            {"schema_version": 1, "knowledge_base": legacy_payload},
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    loaded = JsonKnowledgeStore(path).load()

    assert isinstance(loaded, KnowledgeBase)
    assert loaded.relationships == ()
    assert loaded.entries[0].record.mathematical_structure is None


def test_relationship_review_persists_without_changing_schema_version(
    tmp_path: Path,
) -> None:
    workspace, _, _ = _example_records()
    reviewed = workspace.review_relationship(
        "small_angle_oscillator_approximation",
        RelationshipReviewStatus.REJECTED,
        reviewer="reviewer:physics",
        notes="Mapping conditions need a more precise parameter statement.",
    )
    store = JsonKnowledgeStore(tmp_path / "knowledge.json")
    store.save(reviewed)

    assert store.load() == reviewed
