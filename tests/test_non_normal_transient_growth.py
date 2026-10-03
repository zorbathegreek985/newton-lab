"""Tests for the bounded Phase 64 non-normal transient-growth study."""

import csv
import hashlib
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.linalg import expm  # type: ignore[import-untyped]

import newton_lab.non_normal_transient_growth as study
from newton_lab.exceptions import ScientificValidationError

REPORT_DIRECTORY = Path("reports/phase_64_non_normal_transient_growth")


def _copy_frozen_inputs(destination: Path) -> None:
    destination.mkdir(parents=True)
    for filename in ("protocol.json", "design_exploration.json"):
        shutil.copy2(REPORT_DIRECTORY / filename, destination / filename)


def test_matrix_pair_has_matched_stable_eigenvalues_and_expected_normality() -> None:
    pair = study.build_matrix_pair(1.0, 2.0, 4.0)

    np.testing.assert_array_equal(np.linalg.eigvals(pair.normal), [-1.0, -2.0])
    np.testing.assert_array_equal(np.linalg.eigvals(pair.non_normal), [-1.0, -2.0])
    assert np.max(np.real(np.linalg.eigvals(pair.normal))) < 0.0
    assert np.max(np.real(np.linalg.eigvals(pair.non_normal))) < 0.0
    np.testing.assert_array_equal(
        pair.normal.T @ pair.normal, pair.normal @ pair.normal.T
    )
    assert not np.allclose(
        pair.non_normal.T @ pair.non_normal,
        pair.non_normal @ pair.non_normal.T,
    )


@pytest.mark.parametrize(
    ("a", "b", "kappa"),
    [(0.0, 2.0, 1.0), (1.0, -2.0, 1.0), (1.0, 1.0, 1.0), (1.0, 2.0, np.inf)],
)
def test_matrix_pair_rejects_invalid_parameters(
    a: float, b: float, kappa: float
) -> None:
    with pytest.raises(ScientificValidationError):
        study.build_matrix_pair(a, b, kappa)


def test_analytic_propagator_matches_independent_matrix_exponential() -> None:
    identity = study.analytic_propagator(1.0, 2.0, 4.0, 0.0)
    np.testing.assert_array_equal(identity, np.eye(2))
    pair = study.build_matrix_pair(1.0, 2.0, 4.0)
    for time in (0.1, 0.37, 1.2, 4.0):
        actual = study.analytic_propagator(1.0, 2.0, 4.0, time)
        np.testing.assert_allclose(actual, expm(pair.non_normal * time), atol=1e-14)


def test_propagator_satisfies_semigroup_and_matrix_differential_equation() -> None:
    pair = study.build_matrix_pair(1.0, 2.0, 4.0)
    first = study.analytic_propagator(1.0, 2.0, 4.0, 0.2)
    second = study.analytic_propagator(1.0, 2.0, 4.0, 0.7)
    combined = study.analytic_propagator(1.0, 2.0, 4.0, 0.9)
    np.testing.assert_allclose(combined, first @ second, rtol=1e-13, atol=1e-14)

    time = 0.8
    step = 1e-5
    derivative = (
        study.analytic_propagator(1.0, 2.0, 4.0, time + step)
        - study.analytic_propagator(1.0, 2.0, 4.0, time - step)
    ) / (2.0 * step)
    expected = pair.non_normal @ study.analytic_propagator(1.0, 2.0, 4.0, time)
    np.testing.assert_allclose(derivative, expected, rtol=1e-9, atol=1e-10)


def test_gain_curve_has_nonamplifying_normal_and_zero_coupling_baselines() -> None:
    times = np.linspace(0.0, 4.0, 401)
    normal = study.calculate_gain_curve(1.0, 2.0, 0.0, times)
    zero_coupling = study.calculate_gain_curve(1.0, 2.0, 0.0, times)
    strong_coupling = study.calculate_gain_curve(1.0, 2.0, 8.0, times)

    assert float(np.max(normal.gain)) == pytest.approx(1.0)
    np.testing.assert_array_equal(normal.gain, zero_coupling.gain)
    assert float(np.max(strong_coupling.gain[1:])) > 1.0


def test_gain_curve_rejects_invalid_time_grid() -> None:
    with pytest.raises(ScientificValidationError):
        study.calculate_gain_curve(1.0, 2.0, 1.0, np.asarray([0.0, 0.0]))


def test_failed_solver_result_is_retained_as_a_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        study,
        "solve_ivp",
        lambda *_args, **_kwargs: SimpleNamespace(
            success=False,
            message="controlled test failure",
            t=np.asarray([0.0]),
            y=np.asarray([[1.0], [0.0]]),
            nfev=1,
        ),
    )
    row, samples = study._integrate_condition(
        condition_id="failure_test",
        kappa=4.0,
        initial_state=np.asarray([1.0, 0.0]),
        times=np.asarray([0.0, 0.5, 1.0]),
        solver_data={
            "method": "DOP853",
            "relative_tolerance": 1e-10,
            "absolute_tolerance": 1e-12,
            "maximum_step": 0.025,
        },
        trajectory_acceptance_tolerance=1e-8,
        tight=False,
    )

    assert row["solver_status"] == "failed"
    assert row["solver_success"] is False
    assert row["numerical_accepted"] is False
    assert "controlled test failure" in row["message"]
    assert samples == []


def test_nonfinite_solver_state_is_retained_as_a_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        study,
        "solve_ivp",
        lambda *_args, **_kwargs: SimpleNamespace(
            success=True,
            message="returned a nonfinite state",
            t=np.asarray([0.0, 1.0]),
            y=np.asarray([[1.0, np.nan], [0.0, 0.0]]),
            nfev=2,
        ),
    )
    row, samples = study._integrate_condition(
        condition_id="nonfinite_test",
        kappa=4.0,
        initial_state=np.asarray([1.0, 0.0]),
        times=np.asarray([0.0, 1.0]),
        solver_data={
            "method": "DOP853",
            "relative_tolerance": 1e-10,
            "absolute_tolerance": 1e-12,
            "maximum_step": 0.025,
        },
        trajectory_acceptance_tolerance=1e-8,
        tight=False,
    )

    assert row["solver_status"] == "failed"
    assert row["solver_success"] is True
    assert row["numerical_accepted"] is False
    assert "nonfinite" in row["message"]
    assert samples == []


def test_primary_study_outputs_hypotheses_manifest_and_preserves_other_records(
    tmp_path: Path,
) -> None:
    unrelated_record = tmp_path / "source_experiment_record.json"
    unrelated_record.write_text('{"status":"unchanged"}\n', encoding="utf-8")
    original_bytes = unrelated_record.read_bytes()
    output = tmp_path / "phase64"
    _copy_frozen_inputs(output)

    metadata = study.run_non_normal_transient_growth_study(output)
    summary = {
        float(row["kappa"]): row
        for row in csv.DictReader((output / "gain_summary.csv").open(encoding="utf-8"))
    }
    hypotheses = json.loads((output / "hypothesis_results.json").read_text())
    manifest = json.loads((output / "sha256_manifest.json").read_text())

    assert metadata["trajectory_integrations_attempted"] == 6
    assert metadata["trajectory_integrations_successful"] == 6
    assert metadata["state_space"] == "abstract_dimensionless_R2"
    assert metadata["norm"] == "Euclidean_induced_matrix_2_norm"
    assert metadata["phase32_status"] == "BLOCKED_AUTHORIZATION"
    assert metadata["phase52_status"] == "NO-GO"
    assert float(summary[0.0]["sampled_max_gain_primary"]) == pytest.approx(1.0)
    assert float(summary[2.0]["sampled_max_gain_primary"]) == pytest.approx(1.0)
    assert float(summary[4.0]["sampled_max_gain_primary"]) == pytest.approx(
        1.17435318, abs=1e-7
    )
    assert float(summary[8.0]["sampled_max_gain_primary"]) == pytest.approx(
        2.08004555, abs=1e-7
    )
    assert hypotheses["H1_finite_time_amplification"]["classification"] == (
        "SUPPORTED_WITHIN_FROZEN_FAMILY"
    )
    assert hypotheses["H2_exceeds_matched_normal"]["classification"] == (
        "SUPPORTED_WITHIN_FROZEN_FAMILY"
    )
    assert hypotheses["H3_analytic_numerical_agreement"]["classification"] == (
        "SUPPORTED_WITHIN_FROZEN_TOLERANCES"
    )
    assert len(manifest) == 10
    for filename, expected_hash in manifest.items():
        actual_hash = (
            hashlib.sha256((output / filename).read_bytes()).hexdigest().upper()
        )
        assert actual_hash == expected_hash.upper()
    assert unrelated_record.read_bytes() == original_bytes
    assert set(path.name for path in output.iterdir()) == set(manifest) | {
        "sha256_manifest.json"
    }


def test_study_rejects_tampered_frozen_protocol(tmp_path: Path) -> None:
    output = tmp_path / "phase64"
    _copy_frozen_inputs(output)
    protocol_path = output / "protocol.json"
    protocol_path.write_text(protocol_path.read_text() + " ", encoding="utf-8")

    with pytest.raises(ScientificValidationError, match="protocol hash mismatch"):
        study.run_non_normal_transient_growth_study(output)


def test_study_refuses_to_overwrite_existing_results(tmp_path: Path) -> None:
    output = tmp_path / "phase64"
    _copy_frozen_inputs(output)
    sentinel = output / "existing_result.txt"
    sentinel.write_text("preserve me", encoding="utf-8")

    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        study.run_non_normal_transient_growth_study(output)
    assert sentinel.read_text(encoding="utf-8") == "preserve me"
