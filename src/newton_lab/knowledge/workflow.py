"""Candidate ingestion state, human review history, and accepted records."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from enum import StrEnum
from typing import Annotated, Self

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    JsonValue,
    StringConstraints,
    field_validator,
    model_validator,
)

from newton_lab.exceptions import RegistryError
from newton_lab.knowledge.registry import (
    EquationRecord,
    EquationRegistry,
    EvidenceCategory,
    RelationshipType,
)
from newton_lab.knowledge.taxonomy import PhysicsTaxonomy, default_physics_taxonomy


class WorkflowModel(BaseModel):
    """Common immutable settings for workflow data."""

    model_config = ConfigDict(frozen=True, extra="forbid", str_strip_whitespace=True)


class ReviewState(StrEnum):
    """Editorial workflow state, separate from scientific verification."""

    DRAFT = "draft"
    NEEDS_REVIEW = "needs_review"
    ACCEPTED = "accepted"
    REJECTED = "rejected"
    SUPERSEDED = "superseded"


class RelationshipReviewStatus(StrEnum):
    """Review status for a graph edge, separate from equation review."""

    CANDIDATE = "candidate_requires_review"
    ACCEPTED = "accepted_after_review"
    REJECTED = "rejected_after_review"


_DIRECTIONAL_RELATIONSHIPS = frozenset(
    {
        RelationshipType.DERIVED_FROM,
        RelationshipType.LIMIT_OF,
        RelationshipType.SPECIAL_CASE_OF,
        RelationshipType.APPROXIMATION_OF,
        RelationshipType.TRANSFORMABLE_INTO,
    }
)


class EquationRelationship(WorkflowModel):
    """A provenance-aware, typed, optionally reviewed equation-graph edge."""

    relationship_id: str = Field(pattern=r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
    source_equation_id: str = Field(pattern=r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
    target_equation_id: str = Field(pattern=r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
    relationship_type: RelationshipType
    description: str = Field(min_length=1)
    justification: str = Field(min_length=1)
    conditions: tuple[str, ...] = ()
    source_reference_ids: tuple[str, ...] = ()
    review_status: RelationshipReviewStatus = RelationshipReviewStatus.CANDIDATE
    directional: bool
    transformation_notes: str | None = None
    reviewed_by: str | None = None
    review_notes: str | None = None

    @model_validator(mode="after")
    def relationship_is_well_formed(self) -> Self:
        if self.source_equation_id == self.target_equation_id:
            raise ValueError("equation relationships cannot link an equation to itself")
        expected_directionality = self.relationship_type in _DIRECTIONAL_RELATIONSHIPS
        if self.directional != expected_directionality:
            raise ValueError(
                f"{self.relationship_type.value} directionality must be "
                f"{expected_directionality}"
            )
        if len(self.source_reference_ids) != len(set(self.source_reference_ids)):
            raise ValueError("relationship source reference IDs must be unique")
        if self.review_status != RelationshipReviewStatus.CANDIDATE and (
            not self.reviewed_by
            or not self.reviewed_by.strip()
            or not self.review_notes
            or not self.review_notes.strip()
        ):
            raise ValueError("reviewed relationships require reviewer and review notes")
        return self

    def record_review(
        self,
        status: RelationshipReviewStatus,
        *,
        reviewer: str,
        notes: str,
    ) -> EquationRelationship:
        """Return a reviewed copy; relationship claims are never auto-promoted."""
        if status == RelationshipReviewStatus.CANDIDATE:
            raise RegistryError(
                "a relationship review must accept or reject a candidate"
            )
        if self.review_status != RelationshipReviewStatus.CANDIDATE:
            raise RegistryError("only candidate relationships can be reviewed")
        return EquationRelationship.model_validate(
            {
                **self.model_dump(mode="python"),
                "review_status": status,
                "reviewed_by": reviewer,
                "review_notes": notes,
            }
        )


class IngestionOrigin(StrEnum):
    """Supported initial sources for equation candidates."""

    MANUAL_STRUCTURED = "manual_structured"
    PROJECT_EXAMPLE = "project_example"
    LOCAL_JSON = "local_json"


class ReviewChecklist(WorkflowModel):
    """Independent review dimensions; ``None`` means not assessed or not applicable."""

    source_content_inspected: bool | None = None
    mathematical_form_reviewed: bool | None = None
    physical_assumptions_reviewed: bool | None = None
    empirical_evidence_assessed: bool | None = None


class ReviewEvent(WorkflowModel):
    """One recorded state transition and the record version reviewed at that time."""

    from_state: ReviewState
    to_state: ReviewState
    actor_id: str = Field(min_length=1)
    occurred_at: datetime
    notes: str | None = None
    reviewed_source_ids: tuple[str, ...] = ()
    checklist: ReviewChecklist | None = None
    changes_made: tuple[str, ...] = ()
    record_snapshot: EquationRecord

    @field_validator("occurred_at")
    @classmethod
    def timestamp_must_include_timezone(cls, value: datetime) -> datetime:
        if value.tzinfo is None or value.utcoffset() is None:
            raise ValueError("review timestamps must include a timezone")
        return value

    @model_validator(mode="after")
    def reviewed_sources_are_present(self) -> Self:
        if len(self.reviewed_source_ids) != len(set(self.reviewed_source_ids)):
            raise ValueError("reviewed source IDs must be unique")
        source_ids = {
            source.source_id for source in self.record_snapshot.source_references
        }
        unknown_ids = set(self.reviewed_source_ids) - source_ids
        if unknown_ids:
            raise ValueError(
                f"review event refers to unknown sources: {sorted(unknown_ids)}"
            )
        if self.checklist and self.checklist.source_content_inspected is True:
            if not self.reviewed_source_ids:
                raise ValueError(
                    "source content cannot be marked inspected without "
                    "reviewed source IDs"
                )
        return self


_ALLOWED_TRANSITIONS: dict[ReviewState, frozenset[ReviewState]] = {
    ReviewState.DRAFT: frozenset({ReviewState.NEEDS_REVIEW}),
    ReviewState.NEEDS_REVIEW: frozenset(
        {ReviewState.DRAFT, ReviewState.ACCEPTED, ReviewState.REJECTED}
    ),
    ReviewState.ACCEPTED: frozenset({ReviewState.SUPERSEDED}),
    ReviewState.REJECTED: frozenset({ReviewState.DRAFT, ReviewState.SUPERSEDED}),
    ReviewState.SUPERSEDED: frozenset(),
}


class CandidateEntry(WorkflowModel):
    """An equation record with its original submission and review history."""

    record: EquationRecord
    origin: IngestionOrigin
    submitted_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    original_content: Annotated[str, StringConstraints(strip_whitespace=False)]
    review_state: ReviewState = ReviewState.DRAFT
    review_history: tuple[ReviewEvent, ...] = ()

    @field_validator("submitted_at")
    @classmethod
    def submission_timestamp_must_include_timezone(cls, value: datetime) -> datetime:
        if value.tzinfo is None or value.utcoffset() is None:
            raise ValueError("submission timestamps must include a timezone")
        return value

    @field_validator("original_content")
    @classmethod
    def original_content_must_be_a_json_object(cls, value: str) -> str:
        try:
            parsed: JsonValue = json.loads(value)
        except json.JSONDecodeError as exc:
            raise ValueError("original_content must contain valid JSON") from exc
        if not isinstance(parsed, dict):
            raise ValueError("original_content must contain a JSON object")
        return value

    @model_validator(mode="after")
    def history_must_match_current_state(self) -> Self:
        state = ReviewState.DRAFT
        for event in self.review_history:
            if event.from_state != state:
                raise ValueError(
                    "review history contains a discontinuous state transition"
                )
            if event.to_state not in _ALLOWED_TRANSITIONS[state]:
                raise ValueError(
                    f"invalid review transition: {event.from_state.value} -> "
                    f"{event.to_state.value}"
                )
            if event.record_snapshot.equation_id != self.record.equation_id:
                raise ValueError("review snapshots must preserve the equation ID")
            state = event.to_state
        if state != self.review_state:
            raise ValueError("review history does not end at review_state")
        return self

    def transition(
        self,
        to_state: ReviewState,
        *,
        actor_id: str,
        occurred_at: datetime | None = None,
        notes: str | None = None,
        reviewed_source_ids: tuple[str, ...] = (),
        checklist: ReviewChecklist | None = None,
        updated_record: EquationRecord | None = None,
        changes_made: tuple[str, ...] = (),
    ) -> CandidateEntry:
        """Return a new entry after a valid, auditable review transition.

        Acceptance requires the source content, mathematical form, and physical
        assumptions to have been reviewed, with every attached source inspected.
        This editorial acceptance does not change scientific verification status.
        """
        if to_state not in _ALLOWED_TRANSITIONS[self.review_state]:
            raise RegistryError(
                f"invalid review transition: {self.review_state.value} -> "
                f"{to_state.value}"
            )
        record = updated_record or self.record
        if record.equation_id != self.record.equation_id:
            raise RegistryError("review edits cannot change the stable equation ID")
        if (
            updated_record is not None
            and updated_record != self.record
            and not changes_made
        ):
            raise RegistryError(
                "record edits during review must be described in changes_made"
            )
        if len(reviewed_source_ids) != len(set(reviewed_source_ids)):
            raise RegistryError("reviewed source IDs must be unique")
        record_source_ids = {source.source_id for source in record.source_references}
        if set(reviewed_source_ids) - record_source_ids:
            raise RegistryError("reviewed source IDs must exist on the reviewed record")

        if to_state == ReviewState.ACCEPTED:
            self._validate_acceptance(
                record=record,
                actor_id=actor_id,
                notes=notes,
                reviewed_source_ids=reviewed_source_ids,
                checklist=checklist,
            )

        event = ReviewEvent(
            from_state=self.review_state,
            to_state=to_state,
            actor_id=actor_id,
            occurred_at=occurred_at or datetime.now(UTC),
            notes=notes,
            reviewed_source_ids=reviewed_source_ids,
            checklist=checklist,
            changes_made=changes_made,
            record_snapshot=record,
        )
        return CandidateEntry(
            record=record,
            origin=self.origin,
            submitted_at=self.submitted_at,
            original_content=self.original_content,
            review_state=to_state,
            review_history=(*self.review_history, event),
        )

    @staticmethod
    def _validate_acceptance(
        *,
        record: EquationRecord,
        actor_id: str,
        notes: str | None,
        reviewed_source_ids: tuple[str, ...],
        checklist: ReviewChecklist | None,
    ) -> None:
        if not actor_id.strip():
            raise RegistryError("acceptance requires a reviewer identifier")
        if not notes or not notes.strip():
            raise RegistryError("acceptance requires review notes")
        if not record.source_references:
            raise RegistryError("acceptance requires at least one recorded source")
        all_source_ids = {source.source_id for source in record.source_references}
        if set(reviewed_source_ids) != all_source_ids:
            raise RegistryError(
                "acceptance requires every attached source to be reviewed"
            )
        if checklist is None or any(
            value is not True
            for value in (
                checklist.source_content_inspected,
                checklist.mathematical_form_reviewed,
                checklist.physical_assumptions_reviewed,
            )
        ):
            raise RegistryError(
                "acceptance requires source, mathematical-form, and assumption review"
            )
        claims = (
            *record.physical_interpretations,
            *record.established_applications,
            *(relation.evidence for relation in record.related_equations),
            *(mapping.evidence for mapping in record.application_mappings),
        )
        empirical_categories = {
            EvidenceCategory.EMPIRICALLY_INVESTIGATED,
            EvidenceCategory.INDEPENDENTLY_VALIDATED,
        }
        needs_empirical_assessment = any(
            claim.category in empirical_categories for claim in claims
        )
        if (
            needs_empirical_assessment
            and checklist.empirical_evidence_assessed is not True
        ):
            raise RegistryError(
                "acceptance requires assessment of the record's empirical claims"
            )


class KnowledgeBase(WorkflowModel):
    """Candidate workspace integrating taxonomy, review workflow, and registry."""

    taxonomy: PhysicsTaxonomy = Field(default_factory=default_physics_taxonomy)
    entries: tuple[CandidateEntry, ...] = ()
    relationships: tuple[EquationRelationship, ...] = ()

    @model_validator(mode="after")
    def entries_form_a_consistent_registry(self) -> Self:
        identifiers = [entry.record.equation_id for entry in self.entries]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("candidate equation IDs must be unique")
        EquationRegistry(
            taxonomy=self.taxonomy,
            equations=tuple(entry.record for entry in self.entries),
        )
        entry_by_id = {entry.record.equation_id: entry for entry in self.entries}
        relationship_ids: set[str] = set()
        relationship_identities: set[tuple[str, str, RelationshipType]] = set()
        for relationship in self.relationships:
            if relationship.relationship_id in relationship_ids:
                raise ValueError("relationship IDs must be unique")
            relationship_ids.add(relationship.relationship_id)
            source_entry = entry_by_id.get(relationship.source_equation_id)
            target_entry = entry_by_id.get(relationship.target_equation_id)
            if source_entry is None or target_entry is None:
                raise ValueError(
                    "relationship references unknown equation IDs: "
                    f"{relationship.source_equation_id!r}, "
                    f"{relationship.target_equation_id!r}"
                )
            known_sources = {
                source.source_id
                for endpoint in (source_entry.record, target_entry.record)
                for source in endpoint.source_references
            }
            unknown_sources = set(relationship.source_reference_ids) - known_sources
            if unknown_sources:
                raise ValueError(
                    "relationship references unknown source IDs: "
                    f"{sorted(unknown_sources)}"
                )
            source_id, target_id = (
                (relationship.source_equation_id, relationship.target_equation_id)
                if relationship.directional
                else tuple(
                    sorted(
                        (
                            relationship.source_equation_id,
                            relationship.target_equation_id,
                        )
                    )
                )
            )
            identity = (source_id, target_id, relationship.relationship_type)
            if identity in relationship_identities:
                raise ValueError(
                    "duplicate relationship identity for endpoint pair and type"
                )
            relationship_identities.add(identity)
        return self

    @classmethod
    def from_equation_registry(cls, registry: EquationRegistry) -> KnowledgeBase:
        """Import project example records as drafts without inventing review events."""
        entries = tuple(
            CandidateEntry(
                record=record,
                origin=IngestionOrigin.PROJECT_EXAMPLE,
                original_content=json.dumps(
                    record.model_dump(mode="json"), ensure_ascii=False
                ),
            )
            for record in registry.equations
        )
        return cls(taxonomy=registry.taxonomy, entries=entries)

    def get_entry(self, equation_id: str) -> CandidateEntry:
        """Find a candidate regardless of workflow state."""
        for entry in self.entries:
            if entry.record.equation_id == equation_id:
                return entry
        raise RegistryError(f"unknown candidate equation ID: {equation_id}")

    def get_accepted(self, equation_id: str) -> EquationRecord:
        """Return an accepted record; drafts and rejected entries are inaccessible."""
        entry = self.get_entry(equation_id)
        if entry.review_state != ReviewState.ACCEPTED:
            raise RegistryError(
                f"equation {equation_id!r} is {entry.review_state.value}, not accepted"
            )
        return entry.record

    def accepted_registry(self) -> EquationRegistry:
        """Build an equation registry from accepted entries only."""
        accepted = tuple(
            entry.record
            for entry in self.entries
            if entry.review_state == ReviewState.ACCEPTED
        )
        accepted_ids = {record.equation_id for record in accepted}
        unresolved = [
            (record.equation_id, relation.target_equation_id)
            for record in accepted
            for relation in record.related_equations
            if relation.target_equation_id not in accepted_ids
        ]
        if unresolved:
            raise RegistryError(
                "accepted records have relationships to non-accepted records: "
                f"{unresolved}"
            )
        return EquationRegistry(taxonomy=self.taxonomy, equations=accepted)

    def add_candidate(self, entry: CandidateEntry) -> KnowledgeBase:
        """Return a workspace with a new candidate after checking IDs and references."""
        if any(
            item.record.equation_id == entry.record.equation_id for item in self.entries
        ):
            raise RegistryError(f"duplicate equation ID: {entry.record.equation_id}")
        return KnowledgeBase(
            taxonomy=self.taxonomy,
            entries=(*self.entries, entry),
            relationships=self.relationships,
        )

    def add_candidates(self, entries: tuple[CandidateEntry, ...]) -> KnowledgeBase:
        """Add a batch together, allowing references between new candidates."""
        existing_ids = {entry.record.equation_id for entry in self.entries}
        added_ids = [entry.record.equation_id for entry in entries]
        if len(added_ids) != len(set(added_ids)) or existing_ids.intersection(
            added_ids
        ):
            raise RegistryError("candidate batch contains duplicate equation IDs")
        return KnowledgeBase(
            taxonomy=self.taxonomy,
            entries=(*self.entries, *entries),
            relationships=self.relationships,
        )

    def add_relationship(self, relationship: EquationRelationship) -> KnowledgeBase:
        """Return a workspace with a validated relationship candidate or review."""
        return KnowledgeBase(
            taxonomy=self.taxonomy,
            entries=self.entries,
            relationships=(*self.relationships, relationship),
        )

    def review_relationship(
        self,
        relationship_id: str,
        status: RelationshipReviewStatus,
        *,
        reviewer: str,
        notes: str,
    ) -> KnowledgeBase:
        """Record a person's decision on a candidate graph edge."""
        found = False
        relationships: list[EquationRelationship] = []
        for relationship in self.relationships:
            if relationship.relationship_id == relationship_id:
                relationships.append(
                    relationship.record_review(status, reviewer=reviewer, notes=notes)
                )
                found = True
            else:
                relationships.append(relationship)
        if not found:
            raise RegistryError(f"unknown relationship ID: {relationship_id}")
        return KnowledgeBase(
            taxonomy=self.taxonomy,
            entries=self.entries,
            relationships=tuple(relationships),
        )

    def transition_entry(
        self,
        equation_id: str,
        to_state: ReviewState,
        *,
        actor_id: str,
        occurred_at: datetime | None = None,
        notes: str | None = None,
        reviewed_source_ids: tuple[str, ...] = (),
        checklist: ReviewChecklist | None = None,
        updated_record: EquationRecord | None = None,
        changes_made: tuple[str, ...] = (),
    ) -> KnowledgeBase:
        """Return a workspace with one entry transitioned and history recorded."""
        updated_entries = tuple(
            entry.transition(
                to_state,
                actor_id=actor_id,
                occurred_at=occurred_at,
                notes=notes,
                reviewed_source_ids=reviewed_source_ids,
                checklist=checklist,
                updated_record=updated_record,
                changes_made=changes_made,
            )
            if entry.record.equation_id == equation_id
            else entry
            for entry in self.entries
        )
        if updated_entries == self.entries:
            raise RegistryError(f"unknown candidate equation ID: {equation_id}")
        return KnowledgeBase(
            taxonomy=self.taxonomy,
            entries=updated_entries,
            relationships=self.relationships,
        )
