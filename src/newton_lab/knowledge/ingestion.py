"""Controlled candidate importers for structured data and local JSON files."""

from __future__ import annotations

import json
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from pydantic import ValidationError

from newton_lab.exceptions import CandidateIngestionError
from newton_lab.knowledge.registry import EquationRecord, EquationRegistry
from newton_lab.knowledge.workflow import CandidateEntry, IngestionOrigin


def _validation_summary(error: ValidationError) -> str:
    return "; ".join(
        f"{'.'.join(str(part) for part in issue['loc']) or '<record>'}: {issue['msg']}"
        for issue in error.errors()
    )


def candidate_from_mapping(
    payload: Mapping[str, Any], *, submitted_at: datetime | None = None
) -> CandidateEntry:
    """Validate a manually supplied structured equation and preserve its JSON form."""
    try:
        original_content = json.dumps(
            payload,
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
        )
    except (TypeError, ValueError) as exc:
        raise CandidateIngestionError(
            f"Candidate input must contain JSON-compatible values: {exc}"
        ) from exc
    try:
        record = EquationRecord.model_validate(payload)
    except ValidationError as exc:
        raise CandidateIngestionError(
            f"Invalid equation candidate: {_validation_summary(exc)}",
            raw_submission=original_content,
        ) from exc
    return CandidateEntry(
        record=record,
        origin=IngestionOrigin.MANUAL_STRUCTURED,
        submitted_at=submitted_at or datetime.now(UTC),
        original_content=original_content,
    )


def candidate_from_project_example(
    record: EquationRecord, *, submitted_at: datetime | None = None
) -> CandidateEntry:
    """Wrap an existing project example as an unreviewed draft candidate."""
    return CandidateEntry(
        record=record,
        origin=IngestionOrigin.PROJECT_EXAMPLE,
        submitted_at=submitted_at or datetime.now(UTC),
        original_content=json.dumps(
            record.model_dump(mode="json"), ensure_ascii=False, allow_nan=False
        ),
    )


def candidates_from_example_registry(
    registry: EquationRegistry,
    *,
    submitted_at: datetime | None = None,
) -> tuple[CandidateEntry, ...]:
    """Import all records from an example registry as drafts in one batch."""
    return tuple(
        candidate_from_project_example(record, submitted_at=submitted_at)
        for record in registry.equations
    )


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON property {key!r}")
        result[key] = value
    return result


def candidate_from_json_file(
    path: str | Path, *, submitted_at: datetime | None = None
) -> CandidateEntry:
    """Read one equation object from a local UTF-8 JSON file as a draft.

    The submitted text is retained verbatim in ``original_content`` (including
    whitespace) for comparison with the normalized Pydantic record.
    """
    candidate_path = Path(path)
    try:
        original_content = candidate_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise CandidateIngestionError(
            f"Could not read candidate file {candidate_path}: {exc}"
        ) from exc
    try:
        payload = json.loads(original_content, object_pairs_hook=_unique_object)
    except (json.JSONDecodeError, ValueError) as exc:
        raise CandidateIngestionError(
            f"Invalid JSON in {candidate_path}: {exc}",
            raw_submission=original_content,
        ) from exc
    if not isinstance(payload, dict):
        raise CandidateIngestionError(
            f"Candidate JSON in {candidate_path} must be a single object",
            raw_submission=original_content,
        )
    try:
        record = EquationRecord.model_validate(payload)
    except ValidationError as exc:
        raise CandidateIngestionError(
            f"Invalid equation candidate in {candidate_path}: "
            f"{_validation_summary(exc)}",
            raw_submission=original_content,
        ) from exc
    return CandidateEntry(
        record=record,
        origin=IngestionOrigin.LOCAL_JSON,
        submitted_at=submitted_at or datetime.now(UTC),
        original_content=original_content,
    )
