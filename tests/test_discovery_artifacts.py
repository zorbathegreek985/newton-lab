"""Persistence tests for completed Phase 18 and Phase 19 case studies."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import newton_lab.cross_domain as cross_domain
import newton_lab.heat_diffusion_case_study as heat_study
from newton_lab.cross_domain import (
    AdjustmentSweepPoint,
    DampingTransferCaseStudy,
    run_damping_transfer_case_study,
)
from newton_lab.exceptions import PersistenceError
from newton_lab.heat_diffusion_case_study import (
    HeatDiffusionSmoothingCaseStudy,
    run_heat_diffusion_smoothing_case_study,
)
from newton_lab.knowledge.artifacts import (
    ARTIFACT_SCHEMA_VERSION,
    DiscoveryArtifactStore,
    build_knowledge_integration_from_artifacts,
    deserialize_case_study,
    serialize_case_study,
)
from newton_lab.knowledge.integration import query_discovery_knowledge


@pytest.fixture(scope="module")
def cases() -> tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy]:
    oscillator = run_damping_transfer_case_study(
        damping_ratios=(0.0, 0.025, 0.05),
        duration_s=6.0,
        output_point_count=101,
    )
    heat = run_heat_diffusion_smoothing_case_study(sample_count=11)
    return oscillator, heat


@pytest.mark.parametrize("case_index", [0, 1])
def test_case_study_json_round_trip_is_lossless_and_does_not_mutate(
    cases: tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy],
    case_index: int,
) -> None:
    original = cases[case_index]
    before = serialize_case_study(original)

    encoded = serialize_case_study(original)
    restored = deserialize_case_study(encoded)

    assert serialize_case_study(restored) == before
    assert serialize_case_study(original) == before
    assert encoded["schema_version"] == ARTIFACT_SCHEMA_VERSION
    assert encoded["content_sha256"]


def test_phase18_preserves_sweep_provenance_failures_and_evidence(
    cases: tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy],
) -> None:
    oscillator, _ = cases
    restored = deserialize_case_study(serialize_case_study(oscillator))
    assert isinstance(restored, DampingTransferCaseStudy)
    assert restored.source_experiment == oscillator.source_experiment
    assert restored.source_hypothesis == oscillator.source_hypothesis
    assert restored.real_world_application_status == "unvalidated"
    assert tuple(point.damping_ratio for point in restored.sweep_points) == (
        oscillator.target_damping_ratios
    )
    assert restored.sweep_points[0].status == oscillator.sweep_points[0].status
    assert restored.source_metric_unit == oscillator.source_metric_unit


def test_artifact_preserves_repeated_sweep_values_and_failed_point(
    cases: tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy],
) -> None:
    oscillator, _ = cases
    failed_point = AdjustmentSweepPoint(
        damping_ratio=0.025,
        status="failed",
        solver_method="test_failure_record",
        failure_type="IntegrationError",
        failure_message="recorded numerical failure",
    )
    repeated = oscillator.model_copy(
        update={
            "target_damping_ratios": (0.0, 0.025, 0.025, 0.05),
            "sweep_points": (
                oscillator.sweep_points[0],
                failed_point,
                oscillator.sweep_points[1],
                oscillator.sweep_points[2],
            ),
        }
    )
    restored = deserialize_case_study(serialize_case_study(repeated))
    assert isinstance(restored, DampingTransferCaseStudy)
    assert restored.target_damping_ratios == (0.0, 0.025, 0.025, 0.05)
    assert restored.sweep_points[1].status == "failed"
    assert restored.sweep_points[1].failure_message == "recorded numerical failure"


def test_phase19_preserves_analytical_samples_units_assumptions_and_limits(
    cases: tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy],
) -> None:
    _, heat = cases
    restored = deserialize_case_study(serialize_case_study(heat))
    assert isinstance(restored, HeatDiffusionSmoothingCaseStudy)
    assert (
        serialize_case_study(restored)["payload"]["source_experiment"]
        == serialize_case_study(heat)["payload"]["source_experiment"]
    )
    assert restored.discovery == heat.discovery
    assert restored.diffusion_model == heat.diffusion_model
    assert restored.evaluation == heat.evaluation
    assert restored.synthetic_signal == heat.synthetic_signal
    assert restored.real_world_or_financial_validation is False
    assert restored.limitations == heat.limitations
    assert restored.source_bvp_run_index == heat.source_bvp_run_index


def test_store_save_load_list_and_overwrite_policy(
    tmp_path: Path,
    cases: tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy],
) -> None:
    store = DiscoveryArtifactStore(tmp_path)
    oscillator, heat = cases
    path = store.save(oscillator)
    store.save(heat)
    assert serialize_case_study(
        store.load(oscillator.case_study_id)
    ) == serialize_case_study(oscillator)
    assert serialize_case_study(store.load(heat.case_study_id)) == serialize_case_study(
        heat
    )
    assert store.list_artifacts() == tuple(
        sorted((oscillator.case_study_id, heat.case_study_id))
    )
    original_bytes = path.read_bytes()
    with pytest.raises(PersistenceError, match="already exists"):
        store.save(oscillator)
    assert path.read_bytes() == original_bytes
    store.save(oscillator, overwrite=True)
    assert serialize_case_study(
        store.load(oscillator.case_study_id)
    ) == serialize_case_study(oscillator)


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("{", "malformed discovery-artifact JSON"),
        ("[]", "must be a JSON object"),
        ('{"schema_version":1,"schema_version":1}', "duplicate JSON property"),
    ],
)
def test_store_rejects_malformed_json(
    tmp_path: Path, contents: str, message: str
) -> None:
    store = DiscoveryArtifactStore(tmp_path)
    path = tmp_path / "bad_artifact.json"
    path.write_text(contents, encoding="utf-8")
    with pytest.raises(PersistenceError, match=message):
        store.load("bad_artifact")


def test_store_rejects_missing_unsupported_conflicting_and_corrupt_artifacts(
    tmp_path: Path,
    cases: tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy],
) -> None:
    store = DiscoveryArtifactStore(tmp_path)
    with pytest.raises(PersistenceError, match="could not read"):
        store.load("missing_case")

    document = serialize_case_study(cases[0])
    document["schema_version"] = 99
    (tmp_path / "unsupported.json").write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(PersistenceError, match="unsupported discovery-artifact"):
        store.load("unsupported")

    document = serialize_case_study(cases[0])
    document["artifact_type"] = "unknown_case"
    (tmp_path / "unknown.json").write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(PersistenceError, match="unknown discovery artifact type"):
        store.load("unknown")

    document = serialize_case_study(cases[0])
    document["payload"].pop("source_experiment")
    # A valid digest isolates typed payload validation from integrity checking.
    from newton_lab.knowledge.artifacts import _payload_hash

    document["content_sha256"] = _payload_hash(document["payload"])
    (tmp_path / f"{cases[0].case_study_id}.json").write_text(
        json.dumps(document), encoding="utf-8"
    )
    with pytest.raises(PersistenceError, match="invalid discovery artifact"):
        store.load(cases[0].case_study_id)

    document = serialize_case_study(cases[0])
    document["artifact_id"] = "conflicting_identity"
    with pytest.raises(PersistenceError, match="identity conflicts"):
        deserialize_case_study(document)

    document = serialize_case_study(cases[0])
    document["payload"]["initial_error"] = float("nan")
    document["content_sha256"] = ""
    with pytest.raises(PersistenceError, match="not valid finite JSON"):
        deserialize_case_study(document)


def test_store_rejects_path_traversal_and_does_not_leave_temp_on_failed_serialization(
    tmp_path: Path,
    cases: tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    store = DiscoveryArtifactStore(tmp_path)
    with pytest.raises(PersistenceError, match="safe lowercase"):
        store.load("../outside")

    def fail_serialization(_: object) -> dict[str, object]:
        raise PersistenceError("intentional serialization failure")

    monkeypatch.setattr(
        "newton_lab.knowledge.artifacts.serialize_case_study", fail_serialization
    )
    with pytest.raises(PersistenceError, match="intentional serialization"):
        store.save(cases[0])
    assert tuple(tmp_path.iterdir()) == ()


def test_integration_from_loaded_artifacts_does_not_rerun_experiments(
    tmp_path: Path,
    cases: tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    oscillator, heat = cases
    store = DiscoveryArtifactStore(tmp_path)
    store.save(oscillator)
    store.save(heat)

    def unexpected_run(*_: object, **__: object) -> object:
        raise AssertionError("artifact integration must not rerun experiments")

    monkeypatch.setattr(cross_domain, "run_damping_transfer_case_study", unexpected_run)
    monkeypatch.setattr(
        heat_study, "run_heat_diffusion_smoothing_case_study", unexpected_run
    )
    index = build_knowledge_integration_from_artifacts(store)

    assert index.experiment_associations[0].association_origin == (
        "retrospective_compatibility_checked"
    )
    query = query_discovery_knowledge(index, "transient_heat_diffusion_1d")
    assert "registry_mathematical_relationship_with_conditions" in (
        query.evidence_classifications
    )
    assert "candidate" in index.report_markdown.casefold()
    assert "unvalidated" in index.report_markdown.casefold()


def test_invalid_required_artifact_is_not_silently_dropped(
    tmp_path: Path,
    cases: tuple[DampingTransferCaseStudy, HeatDiffusionSmoothingCaseStudy],
) -> None:
    oscillator, heat = cases
    store = DiscoveryArtifactStore(tmp_path)
    store.save(oscillator)
    path = store.save(heat)
    path.write_text("{", encoding="utf-8")
    with pytest.raises(PersistenceError, match="malformed discovery-artifact"):
        build_knowledge_integration_from_artifacts(store)
