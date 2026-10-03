"""Tests for held-out Phase 44 robust abrupt-step detection."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from newton_lab.step_detection_robustness import (
    default_conditions,
    run_step_detection_study,
    select_component_threshold,
    summarize_predictions,
    write_step_detection_study,
)
from newton_lab.synthetic_timeseries_benchmark import generate_signal


def test_ground_truth_labels_and_composed_signal_components() -> None:
    conditions = {
        condition.condition_id: condition for condition in default_conditions()
    }

    assert conditions["stationary_null"].target_step_present is False
    assert conditions["trend_0.002"].target_step_present is False
    assert conditions["damped_0.002"].target_step_present is False
    assert conditions["mean_step_0.1"].target_step_present is True
    assert conditions["weak_step_plus_trend"].target_step_present is True
    assert conditions["weak_step_plus_damping"].target_step_present is True
    assert conditions["weak_step_plus_trend"].trend_component_present is True
    assert conditions["weak_step_plus_damping"].damped_component_present is True

    trend = generate_signal(
        "trend_noise",
        sample_count=512,
        seed=8,
        parameter_value=1.0,
        trend_slope=0.001,
        additional_mean_step=0.2,
    )
    damped = generate_signal(
        "damped_oscillation",
        sample_count=512,
        seed=8,
        parameter_value=1.0,
        damping_rate=0.002,
        additional_mean_step=0.2,
    )
    assert np.allclose(np.diff(trend.truth[:256]), 0.001)
    assert np.allclose(trend.truth[256:] - (0.001 * np.arange(256, 512)), 0.2)
    assert np.allclose(
        damped.truth[256:]
        - np.exp(-0.002 * np.arange(256, 512))
        * np.sin(2.0 * np.pi * 0.1 * np.arange(256, 512)),
        0.2,
    )


def test_confusion_calculation_uses_correct_class_denominators() -> None:
    summary = summarize_predictions(
        "step",
        "baseline",
        "overall",
        (True, True, False, False),
        (True, False, True, False),
    )

    assert (summary.true_positive, summary.false_positive) == (1, 1)
    assert (summary.false_negative, summary.true_negative) == (1, 1)
    assert summary.false_positive_rate == 0.5
    assert summary.false_negative_rate == 0.5


def test_calibration_and_evaluation_are_disjoint_and_thresholds_are_frozen() -> None:
    study = run_step_detection_study(calibration_replicates=5, evaluation_replicates=7)
    calibration = [row for row in study.records if row.split == "calibration"]
    evaluation = [row for row in study.records if row.split == "evaluation"]
    assert set(study.calibration_seeds).isdisjoint(study.evaluation_seeds)
    assert len(calibration) == len(study.conditions) * 5
    assert len(evaluation) == len(study.conditions) * 7

    expected_trend = select_component_threshold(
        tuple(row.estimated_global_slope for row in calibration),
        tuple(row.trend_component_present for row in calibration),
    )
    expected_envelope = select_component_threshold(
        tuple(row.envelope_decay_fraction for row in calibration),
        tuple(row.damped_component_present for row in calibration),
    )
    assert study.thresholds.trend_abs_slope == expected_trend[0]
    assert study.thresholds.trend_calibration_youden_j == expected_trend[1]
    assert study.thresholds.envelope_decay_fraction == expected_envelope[0]
    assert study.thresholds.envelope_calibration_youden_j == expected_envelope[1]


def test_writer_outputs_are_reproducible_and_include_weak_steps_and_nulls(
    tmp_path: Path,
) -> None:
    first = write_step_detection_study(
        tmp_path / "phase44", calibration_replicates=4, evaluation_replicates=6
    )
    output = tmp_path / "phase44"
    before = {
        name: (output / name).read_bytes()
        for name in ("per_record_results.csv", "summary.csv", "metadata.json")
    }
    second = write_step_detection_study(
        output, calibration_replicates=4, evaluation_replicates=6
    )
    after = {name: (output / name).read_bytes() for name in before}

    assert first == second
    assert before == after
    metadata = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
    assert set(metadata["calibration_seeds"]).isdisjoint(metadata["evaluation_seeds"])
    assert (output / "per_record_results.csv").is_file()
    assert (output / "summary.csv").is_file()
    assert (output / "rule_error_rates.png").is_file()
    assert (output / "weak_step_detection.png").is_file()
    assert any(condition.step_magnitude == 0.1 for condition in second.conditions)
    assert any(not condition.target_step_present for condition in second.conditions)
