"""Versioned, offline-first JSON persistence for a knowledge workspace."""

from __future__ import annotations

import errno
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Self

from pydantic import BaseModel, ConfigDict, ValidationError, model_validator

from newton_lab.exceptions import PersistenceError
from newton_lab.knowledge.workflow import KnowledgeBase

STORE_SCHEMA_VERSION = 1


class KnowledgeStoreDocument(BaseModel):
    """On-disk envelope; version changes require an explicit migration."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    schema_version: int
    knowledge_base: KnowledgeBase

    @model_validator(mode="after")
    def supported_version(self) -> Self:
        if self.schema_version != STORE_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported knowledge-store schema version {self.schema_version}; "
                f"this version supports {STORE_SCHEMA_VERSION}"
            )
        return self


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON property {key!r}")
        result[key] = value
    return result


class JsonKnowledgeStore:
    """Load and save a workspace using a versioned JSON document.

    Writes use a same-directory temporary file and atomic replacement. Existing
    files are never overwritten unless the caller explicitly opts in.
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)

    def save(self, knowledge_base: KnowledgeBase, *, overwrite: bool = False) -> None:
        """Persist a complete validated workspace without silent replacement."""
        if not self.path.parent.exists():
            raise PersistenceError(
                f"store directory does not exist: {self.path.parent}"
            )
        if self.path.exists() and not overwrite:
            raise PersistenceError(
                f"store already exists: {self.path}; pass overwrite=True to replace it"
            )

        document = KnowledgeStoreDocument(
            schema_version=STORE_SCHEMA_VERSION,
            knowledge_base=knowledge_base,
        )
        serialized = (
            json.dumps(
                document.model_dump(mode="json"),
                ensure_ascii=False,
                allow_nan=False,
                sort_keys=True,
                indent=2,
            )
            + "\n"
        )

        temporary_path: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="\n",
                dir=self.path.parent,
                prefix=f".{self.path.name}.",
                suffix=".tmp",
                delete=False,
            ) as temporary_file:
                temporary_path = Path(temporary_file.name)
                temporary_file.write(serialized)
                temporary_file.flush()
                os.fsync(temporary_file.fileno())

            if overwrite:
                os.replace(temporary_path, self.path)
            else:
                # Hard-link creation is atomic and fails if another writer created
                # the destination after the initial existence check.
                os.link(temporary_path, self.path)
                temporary_path.unlink()
        except FileExistsError as exc:
            raise PersistenceError(f"store already exists: {self.path}") from exc
        except OSError as exc:
            detail = (
                "atomic write was unavailable on this filesystem"
                if exc.errno in {errno.EPERM, errno.EOPNOTSUPP, errno.ENOTSUP}
                else str(exc)
            )
            raise PersistenceError(
                f"could not save knowledge store {self.path}: {detail}"
            ) from exc
        finally:
            if temporary_path is not None and temporary_path.exists():
                temporary_path.unlink()

    def load(self) -> KnowledgeBase:
        """Read and fully validate one supported store document."""
        try:
            contents = self.path.read_text(encoding="utf-8")
        except OSError as exc:
            raise PersistenceError(
                f"could not read knowledge store {self.path}: {exc}"
            ) from exc
        try:
            payload = json.loads(contents, object_pairs_hook=_unique_object)
        except (json.JSONDecodeError, ValueError) as exc:
            raise PersistenceError(f"malformed knowledge-store JSON: {exc}") from exc
        if not isinstance(payload, dict):
            raise PersistenceError("knowledge-store document must be a JSON object")
        version = payload.get("schema_version")
        if version != STORE_SCHEMA_VERSION:
            raise PersistenceError(
                f"unsupported knowledge-store schema version {version!r}; "
                f"this version supports {STORE_SCHEMA_VERSION}"
            )
        try:
            document = KnowledgeStoreDocument.model_validate(payload)
        except ValidationError as exc:
            raise PersistenceError(f"invalid knowledge-store document: {exc}") from exc
        return document.knowledge_base
