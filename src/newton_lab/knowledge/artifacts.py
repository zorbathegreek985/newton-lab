"""Versioned local persistence for completed discovery case studies.

Artifacts contain the complete validated Phase 18 or Phase 19 model as
JSON-compatible data. Loading validates data only; it never runs an experiment
or numerical solver.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from copy import deepcopy
from pathlib import Path
from typing import Any, cast

from pydantic import BaseModel, ConfigDict, ValidationError, model_validator

from newton_lab.boundary_value import (
    BoundaryValueProblem,
    make_linear_heat_conduction_problem,
)
from newton_lab.cross_domain import DampingTransferCaseStudy
from newton_lab.exceptions import PersistenceError
from newton_lab.heat_diffusion_case_study import HeatDiffusionSmoothingCaseStudy
from newton_lab.knowledge.integration import (
    KnowledgeIntegrationIndex,
    build_knowledge_integration,
)

ARTIFACT_SCHEMA_VERSION = 1
DEFAULT_ARTIFACT_DIRECTORY = Path(".newton_lab") / "artifacts"
_IDENTITY_PATTERN = re.compile(r"^[a-z0-9]+(?:_[a-z0-9]+)*$")
_ARTIFACT_TYPES = {
    "phase18_damping_transfer": DampingTransferCaseStudy,
    "phase19_heat_diffusion_smoothing": HeatDiffusionSmoothingCaseStudy,
}
CaseStudyArtifact = DampingTransferCaseStudy | HeatDiffusionSmoothingCaseStudy


def _known_heat_problem(bvp: dict[str, Any]) -> BoundaryValueProblem:
    """Construct the declarative built-in heat problem to bind its callbacks."""
    parameters = {item["name"]: item["value"] for item in bvp["parameters"]}
    if bvp.get("problem_id") != "steady_linear_heat_conduction":
        raise PersistenceError("unsupported BVP problem identity in artifact")
    return make_linear_heat_conduction_problem(
        length_m=bvp["domain"][1] - bvp["domain"][0],
        left_temperature_k=parameters["left_temperature_k"],
        right_temperature_k=parameters["right_temperature_k"],
        mesh_points=len(bvp["mesh"]),
    )


class DiscoveryArtifactDocument(BaseModel):
    """Strict JSON envelope around one typed case-study payload."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    schema_version: int
    artifact_id: str
    artifact_type: str
    payload: dict[str, Any]
    content_sha256: str

    @model_validator(mode="after")
    def supported_schema(self) -> DiscoveryArtifactDocument:
        if self.schema_version != ARTIFACT_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported discovery-artifact schema version "
                f"{self.schema_version}; this version supports "
                f"{ARTIFACT_SCHEMA_VERSION}"
            )
        if not _IDENTITY_PATTERN.fullmatch(self.artifact_id):
            raise ValueError("artifact_id must be a safe lowercase identifier")
        if self.artifact_type not in _ARTIFACT_TYPES:
            raise ValueError(f"unknown discovery artifact type {self.artifact_type!r}")
        return self


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON property {key!r}")
        result[key] = value
    return result


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _payload_hash(payload: dict[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _case_payload(artifact: CaseStudyArtifact) -> dict[str, Any]:
    """Dump the typed model, replacing only the known BVP callbacks by IDs."""
    raw = artifact.model_dump(mode="python")
    bvp = raw["source_experiment"]["specification"].get("bvp_problem")
    if bvp is not None:
        canonical = _known_heat_problem(bvp)
        if (
            bvp.get("problem_id") != "steady_linear_heat_conduction"
            or getattr(artifact, "case_study_id", None)
            != "phase19_heat_diffusion_smoothing"
            or getattr(bvp.get("ode_function"), "__code__", None)
            is not canonical.ode_function.__code__
            or getattr(bvp.get("boundary_function"), "__code__", None)
            is not canonical.boundary_function.__code__
        ):
            raise PersistenceError(
                "unsupported executable callback in discovery artifact"
            )
        bvp["ode_function"] = {"callback_id": "steady_heat_ode_v1"}
        bvp["boundary_function"] = {"callback_id": "fixed_endpoint_temperatures_v1"}
    try:
        return cast(
            dict[str, Any],
            json.loads(json.dumps(raw, allow_nan=False, ensure_ascii=False)),
        )
    except (TypeError, ValueError) as exc:
        raise PersistenceError(
            f"case study contains data that cannot be safely serialized: {exc}"
        ) from exc


def _restore_known_callbacks(
    artifact_type: str, payload: dict[str, Any]
) -> dict[str, Any]:
    """Rebind documented symbolic callbacks to the built-in steady heat model."""
    result = deepcopy(payload)
    bvp = (
        result.get("source_experiment", {}).get("specification", {}).get("bvp_problem")
    )
    if bvp is None:
        return result
    if artifact_type != "phase19_heat_diffusion_smoothing":
        raise PersistenceError("BVP callback fields are invalid for this artifact type")
    if bvp.get("ode_function") != {"callback_id": "steady_heat_ode_v1"} or bvp.get(
        "boundary_function"
    ) != {"callback_id": "fixed_endpoint_temperatures_v1"}:
        raise PersistenceError("unsupported or missing BVP callback identity")
    canonical = _known_heat_problem(bvp)
    bvp["ode_function"] = canonical.ode_function
    bvp["boundary_function"] = canonical.boundary_function
    return result


def serialize_case_study(artifact: CaseStudyArtifact) -> dict[str, Any]:
    """Return a complete JSON-compatible typed artifact document."""
    if isinstance(artifact, DampingTransferCaseStudy):
        artifact_type = "phase18_damping_transfer"
    elif isinstance(artifact, HeatDiffusionSmoothingCaseStudy):
        artifact_type = "phase19_heat_diffusion_smoothing"
    else:
        raise PersistenceError(f"unsupported case-study type {type(artifact).__name__}")
    payload = _case_payload(artifact)
    return DiscoveryArtifactDocument(
        schema_version=ARTIFACT_SCHEMA_VERSION,
        artifact_id=artifact.case_study_id,
        artifact_type=artifact_type,
        payload=payload,
        content_sha256=_payload_hash(payload),
    ).model_dump(mode="json")


def deserialize_case_study(document: dict[str, Any]) -> CaseStudyArtifact:
    """Validate envelope, content integrity, identity, and concrete case schema."""
    try:
        envelope = DiscoveryArtifactDocument.model_validate(document)
        try:
            actual_hash = _payload_hash(envelope.payload)
        except (TypeError, ValueError) as exc:
            raise PersistenceError(
                "discovery-artifact payload is not valid finite JSON data"
            ) from exc
        if actual_hash != envelope.content_sha256:
            raise PersistenceError("discovery-artifact content hash does not match")
        payload = _restore_known_callbacks(envelope.artifact_type, envelope.payload)
        if envelope.artifact_type == "phase18_damping_transfer":
            validated = DampingTransferCaseStudy.model_validate(payload)
            artifact: CaseStudyArtifact = validated.model_copy(
                update={"report_markdown": envelope.payload["report_markdown"]}
            )
        else:
            validated_heat = HeatDiffusionSmoothingCaseStudy.model_validate(payload)
            artifact = validated_heat.model_copy(
                update={"report_markdown": envelope.payload["report_markdown"]}
            )
    except PersistenceError:
        raise
    except ValidationError as exc:
        raise PersistenceError(f"invalid discovery artifact: {exc}") from exc
    except (KeyError, TypeError) as exc:
        raise PersistenceError(f"invalid discovery artifact structure: {exc}") from exc
    except ValueError as exc:
        raise PersistenceError(f"invalid discovery artifact: {exc}") from exc
    if artifact.case_study_id != envelope.artifact_id:
        raise PersistenceError(
            "artifact identity conflicts with the embedded case-study identity"
        )
    return artifact


class DiscoveryArtifactStore:
    """Store and load Phase 18/19 artifacts under one local directory.

    The default root is ``.newton_lab/artifacts`` relative to the current
    working directory. Writes use same-directory temporary files and an atomic
    create/replace operation. Existing artifacts require ``overwrite=True``.
    """

    def __init__(self, directory: str | Path = DEFAULT_ARTIFACT_DIRECTORY) -> None:
        self.directory = Path(directory)

    def _path(self, artifact_id: str) -> Path:
        if not _IDENTITY_PATTERN.fullmatch(artifact_id):
            raise PersistenceError("artifact identity must be a safe lowercase ID")
        root = self.directory.resolve()
        path = self.directory / f"{artifact_id}.json"
        if path.resolve().parent != root:
            raise PersistenceError("artifact path must remain inside its directory")
        return path

    def save(self, artifact: CaseStudyArtifact, *, overwrite: bool = False) -> Path:
        """Persist a validated case without replacing an existing identity."""
        path = self._path(artifact.case_study_id)
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise PersistenceError(
                f"could not create artifact directory {path.parent}: {exc}"
            ) from exc
        if path.exists() and not overwrite:
            raise PersistenceError(
                f"artifact already exists: {path}; pass overwrite=True to replace it"
            )
        document = serialize_case_study(artifact)
        serialized = (
            json.dumps(
                document,
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
                dir=path.parent,
                prefix=f".{path.name}.",
                suffix=".tmp",
                delete=False,
            ) as temporary_file:
                temporary_path = Path(temporary_file.name)
                temporary_file.write(serialized)
                temporary_file.flush()
                os.fsync(temporary_file.fileno())
            if overwrite:
                os.replace(temporary_path, path)
            else:
                os.link(temporary_path, path)
                temporary_path.unlink()
        except FileExistsError as exc:
            raise PersistenceError(f"artifact already exists: {path}") from exc
        except OSError as exc:
            raise PersistenceError(
                f"could not save discovery artifact {path}: {exc}"
            ) from exc
        finally:
            if temporary_path is not None and temporary_path.exists():
                temporary_path.unlink()
        return path

    def load(self, artifact_id: str) -> CaseStudyArtifact:
        """Load by safe artifact ID and validate without running solvers."""
        path = self._path(artifact_id)
        if path.resolve().parent != self.directory.resolve():
            raise PersistenceError("artifact path must remain inside its directory")
        try:
            contents = path.read_text(encoding="utf-8")
        except OSError as exc:
            raise PersistenceError(
                f"could not read discovery artifact {path}: {exc}"
            ) from exc
        try:
            document = json.loads(contents, object_pairs_hook=_unique_object)
        except (json.JSONDecodeError, ValueError) as exc:
            raise PersistenceError(f"malformed discovery-artifact JSON: {exc}") from exc
        if not isinstance(document, dict):
            raise PersistenceError("discovery-artifact document must be a JSON object")
        try:
            artifact = deserialize_case_study(document)
        except ValidationError as exc:
            raise PersistenceError(f"invalid discovery artifact: {exc}") from exc
        if artifact.case_study_id != artifact_id:
            raise PersistenceError(
                f"requested artifact ID {artifact_id!r} conflicts with embedded ID "
                f"{artifact.case_study_id!r}"
            )
        return artifact

    def list_artifacts(self) -> tuple[str, ...]:
        """Return valid-looking artifact IDs in stable filename order."""
        if not self.directory.is_dir():
            return ()
        return tuple(
            path.stem
            for path in sorted(
                self.directory.glob("*.json"), key=lambda item: item.name
            )
            if _IDENTITY_PATTERN.fullmatch(path.stem)
        )


def build_knowledge_integration_from_artifacts(
    store: DiscoveryArtifactStore,
    *,
    oscillator_artifact_id: str = "oscillator_to_second_order_adjustment",
    heat_artifact_id: str = "phase19_heat_diffusion_smoothing",
) -> KnowledgeIntegrationIndex:
    """Load both required artifacts and reuse the Phase 20 integration API."""
    oscillator = store.load(oscillator_artifact_id)
    heat = store.load(heat_artifact_id)
    if not isinstance(oscillator, DampingTransferCaseStudy):
        raise PersistenceError(
            f"artifact {oscillator_artifact_id!r} is not a Phase 18 case study"
        )
    if not isinstance(heat, HeatDiffusionSmoothingCaseStudy):
        raise PersistenceError(
            f"artifact {heat_artifact_id!r} is not a Phase 19 case study"
        )
    return build_knowledge_integration(oscillator, heat)
