"""Tests for Phase 50's fold-controlled-system validation experiment."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path
from typing import cast

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.fold_control_validation import (
    estimate_fold_recovery_rate,
    fold_equilibria,
    fold_exact_perturbation,
    fold_local_recovery_rate,
    run_fold_control_validation,
)


def test_fold_equilibria_and_stability_eigenvalue() -> None:
    stable, unstable = fold_equilibria(0.25)
    assert stable == pytest.approx(0.5)
    assert unstable == pytest.approx(-0.5)
    assert fold_local_recovery_rate(0.25) == pytest.approx(1.0)


@pytest.mark.parametrize("margin", [0.0, -0.1, float("nan"), float("inf")])
def test_fold_boundary_and_invalid_margins_are_not_stable_conditions(
    margin: float,
) -> None:
    with pytest.raises(ScientificValidationError):
        fold_equilibria(margin)


def test_exact_fold_trajectory_has_expected_initial_state_and_decay() -> None:
    margin = 0.16
    times = np.linspace(0.0, 4.0, 201)
    small_perturbation = fold_exact_perturbation(margin, 0.01, times)
    larger_perturbation = fold_exact_perturbation(margin, 0.05, times)
    assert small_perturbation[0] == pytest.approx(0.01 * np.sqrt(margin))
    assert small_perturbation[-1] < small_perturbation[0]
    assert larger_perturbation[0] > small_perturbation[0]
    assert np.all(np.diff(small_perturbation) < 0.0)


def test_exact_fold_trajectory_rejects_boolean_perturbation() -> None:
    with pytest.raises(ScientificValidationError):
        fold_exact_perturbation(0.1, cast(float, True), np.array([0.0, 1.0]))


def test_recovery_rate_estimator_recovers_exponential_decay() -> None:
    times = np.linspace(0.0, 5.0, 101)
    measured_rate, r_squared = estimate_fold_recovery_rate(times, np.exp(-0.7 * times))
    assert measured_rate == pytest.approx(0.7)
    assert r_squared == pytest.approx(1.0)


def test_recovery_rate_estimator_rejects_invalid_observations() -> None:
    with pytest.raises(ScientificValidationError):
        estimate_fold_recovery_rate(np.array([0.0, 1.0]), np.array([1.0, 0.5]))
    with pytest.raises(ScientificValidationError):
        estimate_fold_recovery_rate(
            np.array([0.0, 1.0, 2.0]), np.array([1.0, 0.0, 0.2])
        )


def test_final_runner_respects_frozen_protocol_and_baseline(tmp_path: Path) -> None:
    project_root = Path(__file__).resolve().parents[1]
    protocol_source = (
        project_root / "reports" / "phase_50_cross_domain_validation" / "protocol.json"
    )
    shutil.copyfile(protocol_source, tmp_path / "protocol.json")
    metadata = run_fold_control_validation(tmp_path)
    assert metadata["protocol_status"] == "frozen_before_final_evaluation"
    assert metadata["record_count"] == 2200
    assert metadata["condition_count"] == 22
    assert metadata["ode_accuracy_criterion_passed"] is True
    assert metadata["primary_rate_measurement_criterion_passed"] is True
    assert metadata["large_perturbation_sensitivity_criterion_passed"] is True
    assert metadata["predictive_value_criterion_passed"] is True
    assert {path.name for path in tmp_path.iterdir()} == {
        "protocol.json",
        "per_condition_results.csv",
        "condition_summary.csv",
        "predictor_comparison.csv",
        "paired_comparison.csv",
        "held_out_recovery_rate_predictions.png",
        "held_out_prediction_error.png",
        "metadata.json",
    }
    with (tmp_path / "predictor_comparison.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        comparisons = list(csv.DictReader(stream))
    assert len(comparisons) == 2
    with (tmp_path / "paired_comparison.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        paired = next(csv.DictReader(stream))
    assert float(paired["proposed_to_baseline_error_ratio"]) < 1.0
    assert (
        json.loads((tmp_path / "metadata.json").read_text(encoding="utf-8"))[
            "protocol_sha256"
        ]
        == metadata["protocol_sha256"]
    )
