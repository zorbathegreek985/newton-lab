"""Tests for safe, versioned knowledge-workspace persistence."""

import json
from datetime import UTC, datetime
from pathlib import Path

import pytest

from newton_lab.exceptions import PersistenceError
from newton_lab.knowledge.examples import build_example_registry
from newton_lab.knowledge.storage import STORE_SCHEMA_VERSION, JsonKnowledgeStore
from newton_lab.knowledge.workflow import KnowledgeBase, ReviewState


def _workspace() -> KnowledgeBase:
    return KnowledgeBase.from_equation_registry(build_example_registry())


def test_store_round_trip_preserves_records_provenance_and_workflow(
    tmp_path: Path,
) -> None:
    workspace = _workspace()
    first_entry = workspace.entries[0]
    workspace = workspace.transition_entry(
        first_entry.record.equation_id,
        ReviewState.NEEDS_REVIEW,
        actor_id="submitter:researcher",
        occurred_at=datetime(2025, 1, 2, tzinfo=UTC),
    )
    path = tmp_path / "knowledge.json"
    store = JsonKnowledgeStore(path)

    store.save(workspace)
    restored = store.load()

    assert restored == workspace
    assert restored.entries[0].original_content
    assert restored.entries[0].record.source_references
    assert restored.entries[0].review_history[0].actor_id == "submitter:researcher"


def test_store_serialization_is_deterministic(tmp_path: Path) -> None:
    first = JsonKnowledgeStore(tmp_path / "first.json")
    second = JsonKnowledgeStore(tmp_path / "second.json")
    workspace = _workspace()

    first.save(workspace)
    second.save(workspace)

    assert first.path.read_bytes() == second.path.read_bytes()


def test_store_does_not_overwrite_existing_file_without_opt_in(
    tmp_path: Path,
) -> None:
    path = tmp_path / "knowledge.json"
    path.write_text("preserve this", encoding="utf-8")
    store = JsonKnowledgeStore(path)

    with pytest.raises(PersistenceError, match="already exists"):
        store.save(_workspace())

    assert path.read_text(encoding="utf-8") == "preserve this"

    workspace = _workspace()
    store.save(workspace, overwrite=True)
    assert store.load() == workspace


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("{", "malformed knowledge-store JSON"),
        ('{"schema_version": 1, "schema_version": 1}', "duplicate JSON property"),
        ("[]", "must be a JSON object"),
        ('{"schema_version": 0}', "unsupported knowledge-store schema version"),
        ('{"schema_version": 2}', "unsupported knowledge-store schema version"),
        (
            '{"schema_version": 1, "knowledge_base": {"entries": [null]}}',
            "invalid knowledge-store",
        ),
    ],
)
def test_store_rejects_malformed_or_unsupported_documents(
    tmp_path: Path, contents: str, message: str
) -> None:
    path = tmp_path / "knowledge.json"
    path.write_text(contents, encoding="utf-8")

    with pytest.raises(PersistenceError, match=message):
        JsonKnowledgeStore(path).load()


def test_store_document_declares_current_schema_version(tmp_path: Path) -> None:
    path = tmp_path / "knowledge.json"
    JsonKnowledgeStore(path).save(_workspace())

    assert json.loads(path.read_text(encoding="utf-8"))["schema_version"] == (
        STORE_SCHEMA_VERSION
    )


def test_store_reports_missing_parent_directory(tmp_path: Path) -> None:
    store = JsonKnowledgeStore(tmp_path / "missing" / "knowledge.json")

    with pytest.raises(PersistenceError, match="directory does not exist"):
        store.save(_workspace())
