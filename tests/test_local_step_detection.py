"""Tests for Phase 45 local midpoint step detection."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.local_step_detection import (
    LOCAL_METHODS,
    WINDOW_SIZES,
    local_mean_difference,
    run_local_step_study,
    trend_adjusted_local_difference,
    write_local_step_study,
)


def test_local_window_indices_and_difference_at_known_midpoint() -> None:
    values = np.zeros(20)
    values[7:10] = 2.0
    values[10:13] = 5.0

    assert local_mean_difference(values, 3) == 3.0


def test_linear_trend_adjustment_removes_known_drift_and_preserves_step() -> None:
    time = np.arange(32, dtype=np.float64)
    slope = 0.25
    step = 4.0
    values = slope * time
    values[16:] += step

    assert local_mean_difference(values, 4) == pytest.approx(step + slope * 4)
    assert trend_adjusted_local_difference(
        values, 4, slope_estimate=slope
    ) == pytest.approx(step)

    smooth = slope * time
    assert trend_adjusted_local_difference(
        smooth, 4, slope_estimate=slope
    ) == pytest.approx(0.0)


def test_local_window_validation_rejects_invalid_sizes_and_data() -> None:
    with pytest.raises(ScientificValidationError):
        local_mean_difference(np.zeros(16), 9)
    with pytest.raises(ScientificValidationError):
        local_mean_difference(np.asarray([[0.0, 1.0]]), 1)


def test_calibration_seeds_and_thresholds_are_separate_from_evaluation() -> None:
    study = run_local_step_study(calibration_replicates=5, evaluation_replicates=7)

    assert set(study.phase44_study.calibration_seeds).isdisjoint(
        study.phase44_study.evaluation_seeds
    )
    assert tuple(threshold.window_size for threshold in study.thresholds) == (
        16,
        16,
        32,
        32,
        64,
        64,
    )
    assert tuple(threshold.method for threshold in study.thresholds) == LOCAL_METHODS
    calibration = [row for row in study.records if row.split == "calibration"]
    for threshold in study.thresholds:
        selected = [row for row in calibration if row.method == threshold.method]
        assert len(selected) == len(study.phase44_study.conditions) * 5
        assert all(
            row.calibrated_threshold == threshold.absolute_difference_threshold
            for row in selected
        )


def test_phase44_baseline_is_preserved_and_summary_counts_match_records() -> None:
    # Use the Phase 44 documented calibration/evaluation design when asserting
    # its frozen default thresholds. Reduced replicate counts intentionally
    # produce different calibration thresholds.
    study = run_local_step_study()
    assert study.phase44_study.thresholds.trend_abs_slope == pytest.approx(
        0.0004527730452371508
    )
    assert study.phase44_study.thresholds.envelope_decay_fraction == pytest.approx(
        0.1563558789122506
    )
    baseline = next(
        row
        for row in study.summaries
        if row.method == "baseline" and row.stratum == "overall"
    )
    selected = [
        row
        for row in study.records
        if row.split == "evaluation" and row.method == "baseline"
    ]
    assert baseline.true_positive == sum(
        row.target_step_present and row.detected for row in selected
    )
    assert baseline.false_positive == sum(
        not row.target_step_present and row.detected for row in selected
    )
    assert baseline.false_negative == sum(
        row.target_step_present and not row.detected for row in selected
    )
    assert baseline.true_negative == sum(
        not row.target_step_present and not row.detected for row in selected
    )


def test_output_schema_and_repeated_run_are_reproducible(tmp_path: Path) -> None:
    output = tmp_path / "phase45"
    first = write_local_step_study(
        output, calibration_replicates=4, evaluation_replicates=6
    )
    names = (
        "per_record_results.csv",
        "summary.csv",
        "metadata.json",
    )
    before = {name: (output / name).read_bytes() for name in names}
    second = write_local_step_study(
        output, calibration_replicates=4, evaluation_replicates=6
    )

    assert first == second
    assert before == {name: (output / name).read_bytes() for name in names}
    metadata = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["local_window_sizes_per_side"] == list(WINDOW_SIZES)
    assert set(metadata["calibration_seeds"]).isdisjoint(metadata["evaluation_seeds"])
    assert (output / "method_error_rates.png").is_file()
    assert (output / "step_magnitude_sensitivity.png").is_file()
