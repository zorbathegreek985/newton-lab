"""Tests for the Phase 46 independently crossed noise study."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from newton_lab.noise_sensitivity import (
    LOCAL_RAW_W32_THRESHOLD,
    LOCAL_TREND_ADJUSTED_W32_THRESHOLD,
    LOCAL_TREND_ADJUSTED_W64_THRESHOLD,
    METHODS,
    NOISE_LEVELS,
    STEP_MAGNITUDES,
    NoiseCondition,
    factorial_conditions,
    generate_noise_condition_signal,
    run_noise_sensitivity_study,
    write_noise_sensitivity_study,
)


def test_noise_control_is_independent_and_reproducible() -> None:
    conditions = [
        condition
        for condition in factorial_conditions()
        if condition.family == "step_plus_trend" and condition.step_magnitude == 0.2
    ]
    low, middle, high = conditions
    first = generate_noise_condition_signal(low, seed=31)
    repeated = generate_noise_condition_signal(low, seed=31)
    medium = generate_noise_condition_signal(middle, seed=31)
    noisy = generate_noise_condition_signal(high, seed=31)

    assert np.array_equal(first.values, repeated.values)
    assert np.array_equal(first.truth, medium.truth)
    assert np.array_equal(medium.truth, noisy.truth)
    assert np.allclose(medium.values - medium.truth, 2.0 * (first.values - first.truth))
    assert np.allclose(noisy.values - noisy.truth, 3.0 * (first.values - first.truth))


def test_factorial_conditions_hold_non_noise_settings_fixed() -> None:
    conditions = factorial_conditions()
    assert len(conditions) == 30
    assert {condition.noise_sd for condition in conditions} == set(NOISE_LEVELS)
    assert {
        condition.step_magnitude
        for condition in conditions
        if condition.target_step_present
    } == set(STEP_MAGNITUDES)
    cells: dict[str, list[NoiseCondition]] = {}
    for condition in conditions:
        key = condition.cell_id.rsplit("_noise_", 1)[0]
        cells.setdefault(key, []).append(condition)
    assert len(cells) == 10
    for cell in cells.values():
        assert len(cell) == len(NOISE_LEVELS)
        assert len({condition.noise_sd for condition in cell}) == 3
        assert len({condition.step_magnitude for condition in cell}) == 1
        assert len({condition.trend_slope for condition in cell}) == 1
        assert len({condition.target_step_present for condition in cell}) == 1
    by_family = {
        family: [condition for condition in conditions if condition.family == family]
        for family in (
            "stationary_null",
            "gradual_trend",
            "pure_step",
            "step_plus_trend",
        )
    }
    assert all(not row.target_step_present for row in by_family["stationary_null"])
    assert all(not row.target_step_present for row in by_family["gradual_trend"])
    assert all(row.target_step_present for row in by_family["pure_step"])
    assert all(row.target_step_present for row in by_family["step_plus_trend"])
    assert {row.trend_slope for row in by_family["gradual_trend"]} == {0.001}
    assert {row.trend_slope for row in by_family["step_plus_trend"]} == {0.001}


def test_study_uses_disjoint_seeds_frozen_thresholds_and_correct_denominators() -> None:
    study = run_noise_sensitivity_study(
        calibration_replicates=3, evaluation_replicates=5
    )

    assert set(study.calibration_seeds).isdisjoint(study.evaluation_seeds)
    assert len(study.calibration_seeds) == 10 * 3
    assert len(study.evaluation_seeds) == 10 * 5
    assert {record.method for record in study.records} == set(METHODS)
    assert {record.threshold_source for record in study.records} == {
        "Phase 42/44 fixed threshold",
        "Phase 45 calibration threshold, frozen",
    }
    expected_thresholds = {
        "phase44_baseline": 3.0,
        "phase45_local_raw_w32": LOCAL_RAW_W32_THRESHOLD,
        "phase45_local_trend_adjusted_w32": LOCAL_TREND_ADJUSTED_W32_THRESHOLD,
        "phase45_local_trend_adjusted_w64": LOCAL_TREND_ADJUSTED_W64_THRESHOLD,
    }
    for method in METHODS:
        assert (
            len({row.threshold for row in study.records if row.method == method}) == 1
        )
        assert {row.threshold for row in study.records if row.method == method} == {
            expected_thresholds[method]
        }

    pure_step = next(
        row
        for row in study.summaries
        if row.method == "phase45_local_raw_w32"
        and row.stratum == "family_noise:pure_step|noise_sd=0.5"
    )
    assert pure_step.positive_count == 4 * 5
    assert pure_step.negative_count == 0
    assert pure_step.true_positive + pure_step.false_negative == 4 * 5
    assert pure_step.false_positive + pure_step.true_negative == 0
    pure_records = [
        row
        for row in study.records
        if row.split == "evaluation"
        and row.method == pure_step.method
        and row.family == "pure_step"
        and row.noise_sd == 0.5
    ]
    assert pure_step.true_positive == sum(
        row.target_step_present and row.detected for row in pure_records
    )
    assert pure_step.false_negative == sum(
        row.target_step_present and not row.detected for row in pure_records
    )

    trend = next(
        row
        for row in study.summaries
        if row.method == "phase44_baseline"
        and row.stratum == "family_noise:gradual_trend|noise_sd=1"
    )
    assert trend.positive_count == 0
    assert trend.negative_count == 5
    assert trend.false_positive + trend.true_negative == 5
    assert trend.false_positive_rate_low_95 is not None


def test_output_schema_and_artifacts_are_reproducible(tmp_path: Path) -> None:
    output = tmp_path / "phase46"
    write_noise_sensitivity_study(
        output, calibration_replicates=2, evaluation_replicates=3
    )
    names = (
        "per_record_results.csv",
        "summary.csv",
        "metadata.json",
        "error_rates_by_noise.png",
        "step_magnitude_sensitivity.png",
    )
    first = {name: (output / name).read_bytes() for name in names}
    write_noise_sensitivity_study(
        output, calibration_replicates=2, evaluation_replicates=3
    )
    assert first == {name: (output / name).read_bytes() for name in names}

    metadata = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["noise_sd_levels"] == list(NOISE_LEVELS)
    assert metadata["threshold_policy"].startswith("no Phase 46 thresholds fitted")
    assert len(metadata["conditions"]) == 30
    with (output / "per_record_results.csv").open(encoding="utf-8") as stream:
        header = stream.readline().strip().split(",")
    assert {"split", "method", "noise_sd", "target_step_present", "threshold"} <= set(
        header
    )
