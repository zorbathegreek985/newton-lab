"""Tests for candidate ingestion and the human review workflow."""

import json
from datetime import UTC, datetime
from pathlib import Path

import pytest
from pydantic import ValidationError

from newton_lab.exceptions import CandidateIngestionError, RegistryError
from newton_lab.knowledge.examples import build_example_registry
from newton_lab.knowledge.ingestion import (
    candidate_from_json_file,
    candidate_from_mapping,
    candidate_from_project_example,
    candidates_from_example_registry,
)
from newton_lab.knowledge.registry import (
    EquationRecord,
    EvidenceCategory,
    EvidenceClaim,
    VerificationStatus,
)
from newton_lab.knowledge.workflow import (
    CandidateEntry,
    IngestionOrigin,
    KnowledgeBase,
    ReviewChecklist,
    ReviewState,
)

FIXED_TIME = datetime(2025, 1, 2, 3, 4, 5, tzinfo=UTC)


def _newton_candidate() -> CandidateEntry:
    source_record = build_example_registry().get("newton_second_law")
    unassessed_record = EquationRecord.model_validate(
        {
            **source_record.model_dump(mode="python"),
            "verification_status": VerificationStatus.UNASSESSED,
        }
    )
    return candidate_from_project_example(unassessed_record, submitted_at=FIXED_TIME)


def _accept(candidate: CandidateEntry) -> CandidateEntry:
    awaiting_review = candidate.transition(
        ReviewState.NEEDS_REVIEW,
        actor_id="submitter:researcher",
        occurred_at=FIXED_TIME,
    )
    source_ids = tuple(
        source.source_id for source in awaiting_review.record.source_references
    )
    return awaiting_review.transition(
        ReviewState.ACCEPTED,
        actor_id="reviewer:physics",
        occurred_at=FIXED_TIME,
        notes="Reviewed the cited source, equation form, and stated assumptions.",
        reviewed_source_ids=source_ids,
        checklist=ReviewChecklist(
            source_content_inspected=True,
            mathematical_form_reviewed=True,
            physical_assumptions_reviewed=True,
        ),
    )


def test_source_metadata_supports_identifiers_dates_type_and_locator() -> None:
    record = build_example_registry().get("newton_second_law")
    source = type(record.source_references[0]).model_validate(
        {
            **record.source_references[0].model_dump(mode="python"),
            "publication_date": "1687-01-01",
            "edition": "First edition",
            "source_type": "archival_source",
            "accessed_date": "2025-01-02",
            "responsible_organization": "The Newton Project",
        }
    )

    assert source.source_id == "newton_principia_1687"
    assert source.publication_date == datetime(1687, 1, 1).date()
    assert source.source_type == "archival_source"
    assert source.accessed_date == datetime(2025, 1, 2).date()


def test_source_ids_must_be_unique_within_an_equation() -> None:
    source = build_example_registry().get("newton_second_law").source_references[0]

    with pytest.raises(ValidationError, match="source reference IDs must be unique"):
        EquationRecord(
            equation_id="duplicate_sources",
            canonical_name="Duplicate sources",
            branch_id="classical_mechanics",
            mathematical_expression="F = ma",
            source_references=(source, source),
        )


def test_multiple_equations_can_reference_the_same_source() -> None:
    examples = build_example_registry()
    shared_source = examples.get("newton_second_law").source_references[0]
    oscillator = EquationRecord.model_validate(
        {
            **examples.get("damped_harmonic_oscillator").model_dump(mode="python"),
            "source_references": (shared_source,),
        }
    )
    entries = candidates_from_example_registry(examples)
    records = tuple(
        candidate_from_project_example(oscillator, submitted_at=entry.submitted_at)
        if entry.record.equation_id == oscillator.equation_id
        else entry
        for entry in entries
    )

    knowledge_base = KnowledgeBase(
        taxonomy=examples.taxonomy,
        entries=records,
    )

    assert (
        len(knowledge_base.get_entry("newton_second_law").record.source_references) == 1
    )
    assert (
        knowledge_base.get_entry("damped_harmonic_oscillator").record.source_references[
            0
        ]
        == shared_source
    )


def test_manual_ingestion_preserves_original_structured_content() -> None:
    payload = {
        "equation_id": "manual_equation",
        "canonical_name": "  Manual equation  ",
        "branch_id": "classical_mechanics",
        "subfield_id": "newtonian_mechanics",
        "mathematical_expression": "x = y",
        "assumptions": ["input assumption"],
    }
    candidate = candidate_from_mapping(payload)

    assert candidate.origin == IngestionOrigin.MANUAL_STRUCTURED
    assert candidate.review_state == ReviewState.DRAFT
    assert candidate.record.canonical_name == "Manual equation"
    assert (
        json.loads(candidate.original_content)["canonical_name"]
        == "  Manual equation  "
    )


def test_invalid_manual_candidate_reports_actionable_validation_errors() -> None:
    with pytest.raises(CandidateIngestionError, match="equation_id") as caught:
        candidate_from_mapping(
            {
                "canonical_name": "Missing stable ID",
                "branch_id": "classical_mechanics",
                "mathematical_expression": "x = y",
            }
        )

    assert caught.value.raw_submission is not None


def test_local_json_ingestion_preserves_exact_submitted_text(tmp_path: Path) -> None:
    raw_json = (
        '{\n  "equation_id": "local_record",\n  "canonical_name": '
        '"Local record",\n  "branch_id": "classical_mechanics",\n  '
        '"mathematical_expression": "x = y"\n}\n'
    )
    path = tmp_path / "candidate.json"
    path.write_text(raw_json, encoding="utf-8")

    candidate = candidate_from_json_file(path, submitted_at=FIXED_TIME)

    assert candidate.origin == IngestionOrigin.LOCAL_JSON
    assert candidate.original_content == raw_json
    assert candidate.record.equation_id == "local_record"


def test_local_json_rejects_malformed_duplicate_or_nonobject_input(
    tmp_path: Path,
) -> None:
    path = tmp_path / "candidate.json"
    path.write_text(
        '{"equation_id": "first", "equation_id": "second"}', encoding="utf-8"
    )

    with pytest.raises(CandidateIngestionError, match="duplicate JSON property"):
        candidate_from_json_file(path)

    path.write_text("[1, 2]", encoding="utf-8")
    with pytest.raises(CandidateIngestionError, match="single object"):
        candidate_from_json_file(path)

    path.write_text("{", encoding="utf-8")
    with pytest.raises(CandidateIngestionError, match="Invalid JSON"):
        candidate_from_json_file(path)


def test_project_example_ingestion_is_draft_and_batch_import_preserves_links() -> None:
    registry = build_example_registry()
    candidate = candidate_from_project_example(registry.get("newton_second_law"))
    all_candidates = candidates_from_example_registry(registry)
    workspace = KnowledgeBase(taxonomy=registry.taxonomy).add_candidates(all_candidates)

    assert candidate.origin == IngestionOrigin.PROJECT_EXAMPLE
    assert candidate.review_state == ReviewState.DRAFT
    assert len(workspace.entries) == len(registry.equations)
    assert workspace.get_entry("damped_nonlinear_pendulum").record.related_equations


def test_review_transitions_record_metadata_and_require_acceptance_review() -> None:
    candidate = _newton_candidate()
    with pytest.raises(RegistryError, match="invalid review transition"):
        candidate.transition(ReviewState.ACCEPTED, actor_id="reviewer")

    awaiting_review = candidate.transition(
        ReviewState.NEEDS_REVIEW,
        actor_id="submitter",
        occurred_at=FIXED_TIME,
    )
    with pytest.raises(RegistryError, match="every attached source"):
        awaiting_review.transition(
            ReviewState.ACCEPTED,
            actor_id="reviewer",
            occurred_at=FIXED_TIME,
            notes="Reviewed.",
            checklist=ReviewChecklist(
                source_content_inspected=True,
                mathematical_form_reviewed=True,
                physical_assumptions_reviewed=True,
            ),
        )


def test_acceptance_requires_distinct_review_checklist_and_source_inspection() -> None:
    awaiting_review = _newton_candidate().transition(
        ReviewState.NEEDS_REVIEW,
        actor_id="submitter",
        occurred_at=FIXED_TIME,
    )
    source_ids = tuple(
        source.source_id for source in awaiting_review.record.source_references
    )
    with pytest.raises(RegistryError, match="source, mathematical-form"):
        awaiting_review.transition(
            ReviewState.ACCEPTED,
            actor_id="reviewer",
            occurred_at=FIXED_TIME,
            notes="Review without complete checks.",
            reviewed_source_ids=source_ids,
            checklist=ReviewChecklist(source_content_inspected=True),
        )


def test_record_edits_are_snapshotted_and_described_in_review_history() -> None:
    awaiting_review = _newton_candidate().transition(
        ReviewState.NEEDS_REVIEW,
        actor_id="submitter",
        occurred_at=FIXED_TIME,
    )
    revised_record = awaiting_review.record.model_copy(
        update={"description": "Description clarified during review."}
    )
    revised = awaiting_review.transition(
        ReviewState.DRAFT,
        actor_id="reviewer",
        occurred_at=FIXED_TIME,
        updated_record=revised_record,
        changes_made=("Clarified the record description.",),
    )

    assert revised.record.description == "Description clarified during review."
    assert revised.review_history[-1].record_snapshot == revised_record
    assert revised.review_history[-1].changes_made == (
        "Clarified the record description.",
    )


def test_acceptance_does_not_promote_scientific_verification() -> None:
    accepted = _accept(_newton_candidate())

    assert accepted.review_state == ReviewState.ACCEPTED
    assert accepted.record.verification_status == VerificationStatus.UNASSESSED
    assert accepted.review_history[-1].actor_id == "reviewer:physics"
    assert accepted.review_history[-1].checklist is not None


def test_empirical_claims_require_empirical_evidence_assessment() -> None:
    record = build_example_registry().get("newton_second_law")
    empirical_claim_record = EquationRecord.model_validate(
        {
            **record.model_dump(mode="python"),
            "physical_interpretations": (
                EvidenceClaim(
                    category=EvidenceCategory.EMPIRICALLY_INVESTIGATED,
                    statement="An empirical claim.",
                    source_reference_ids=("newton_principia_1687",),
                ),
            ),
        }
    )
    candidate = candidate_from_project_example(empirical_claim_record).transition(
        ReviewState.NEEDS_REVIEW,
        actor_id="submitter",
        occurred_at=FIXED_TIME,
    )
    source_ids = tuple(
        source.source_id for source in candidate.record.source_references
    )
    checklist = ReviewChecklist(
        source_content_inspected=True,
        mathematical_form_reviewed=True,
        physical_assumptions_reviewed=True,
    )

    with pytest.raises(RegistryError, match="empirical claims"):
        candidate.transition(
            ReviewState.ACCEPTED,
            actor_id="reviewer",
            occurred_at=FIXED_TIME,
            notes="Reviewed.",
            reviewed_source_ids=source_ids,
            checklist=checklist,
        )


def test_review_state_can_be_rejected_returned_to_draft_then_superseded() -> None:
    submitted = _newton_candidate().transition(
        ReviewState.NEEDS_REVIEW,
        actor_id="submitter",
        occurred_at=FIXED_TIME,
    )
    rejected = submitted.transition(
        ReviewState.REJECTED,
        actor_id="reviewer",
        occurred_at=FIXED_TIME,
        notes="Source passage needs clarification.",
    )
    returned = rejected.transition(
        ReviewState.DRAFT,
        actor_id="author",
        occurred_at=FIXED_TIME,
        notes="Candidate revised for resubmission.",
    )
    accepted = _accept(returned)
    superseded = accepted.transition(
        ReviewState.SUPERSEDED,
        actor_id="curator",
        occurred_at=FIXED_TIME,
        notes="Replaced by a corrected record.",
    )

    assert superseded.review_state == ReviewState.SUPERSEDED
    assert len(superseded.review_history) == 6
    with pytest.raises(RegistryError, match="invalid review transition"):
        superseded.transition(ReviewState.DRAFT, actor_id="author")


def test_workspace_prevents_drafts_from_accepted_registry() -> None:
    accepted = _accept(_newton_candidate())
    draft = candidate_from_mapping(
        {
            "equation_id": "another_draft",
            "canonical_name": "Another draft",
            "branch_id": "classical_mechanics",
            "mathematical_expression": "u = v",
        }
    )
    workspace = KnowledgeBase(
        entries=(accepted, draft),
    )

    assert workspace.get_accepted("newton_second_law") == accepted.record
    with pytest.raises(RegistryError, match="not accepted"):
        workspace.get_accepted("another_draft")
    assert tuple(
        record.equation_id for record in workspace.accepted_registry().equations
    ) == ("newton_second_law",)


def test_workspace_rejects_duplicate_candidate_ids() -> None:
    candidate = _newton_candidate()
    workspace = KnowledgeBase(entries=(candidate,))

    with pytest.raises(RegistryError, match="duplicate equation ID"):
        workspace.add_candidate(candidate)


def test_accepted_registry_requires_related_records_to_be_accepted() -> None:
    all_examples = candidates_from_example_registry(build_example_registry())
    pendulum_candidate = next(
        entry
        for entry in all_examples
        if entry.record.equation_id == "damped_nonlinear_pendulum"
    )
    newton_draft = next(
        entry
        for entry in all_examples
        if entry.record.equation_id == "newton_second_law"
    )
    pendulum_record = EquationRecord.model_validate(
        {
            **pendulum_candidate.record.model_dump(mode="python"),
            "source_references": newton_draft.record.source_references,
        }
    )
    accepted_pendulum = _accept(candidate_from_project_example(pendulum_record))
    workspace = KnowledgeBase(entries=(accepted_pendulum, newton_draft))

    with pytest.raises(RegistryError, match="non-accepted records"):
        workspace.accepted_registry()
