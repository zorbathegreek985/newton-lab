"""Tests for Phase 43 controlled diagnostic robustness experiments."""

from __future__ import annotations

from pathlib import Path

from newton_lab.diagnostic_robustness import (
    RobustnessRow,
    run_robustness_study,
    summarize_confusion,
    write_robustness_study,
)


def test_study_is_reproducible_and_keeps_explicit_binary_classes() -> None:
    first = run_robustness_study()
    second = run_robustness_study()

    assert first == second
    assert len(first.rows) == 460
    assert {summary.diagnostic for summary in first.summaries} == {
        "mean_shift",
        "variance_shift",
    }
    assert all(summary.positive_count > 0 for summary in first.summaries)
    assert all(summary.negative_count > 0 for summary in first.summaries)
    assert all(0.0 <= summary.false_positive_rate <= 1.0 for summary in first.summaries)
    assert all(0.0 <= summary.false_negative_rate <= 1.0 for summary in first.summaries)


def test_parameter_levels_include_null_controls_and_multiple_seeds() -> None:
    study = run_robustness_study()
    trend = [
        row
        for row in study.rows
        if row.diagnostic == "mean_shift"
        and row.process == "trend_noise"
        and row.parameter_name == "trend_slope"
    ]
    damping = [
        row
        for row in study.rows
        if row.diagnostic == "variance_shift"
        and row.process == "damped_oscillation"
        and row.parameter_name == "damping_rate"
    ]

    assert {row.parameter_value for row in trend} == {
        0.0,
        0.0005,
        0.001,
        0.002,
        0.004,
    }
    assert {row.parameter_value for row in damping} == {0.0, 0.001, 0.002, 0.004}
    assert all(not row.target_present for row in trend + damping)
    assert len({row.seed for row in trend}) == 100


def test_confusion_metric_calculation_uses_class_denominators() -> None:
    rows = (
        RobustnessRow(
            "mean_shift", "mean_shift", True, "step", 1.0, "{}", 1, 512, 4.0, True
        ),
        RobustnessRow(
            "mean_shift", "mean_shift", True, "step", 0.1, "{}", 2, 512, 1.0, False
        ),
        RobustnessRow(
            "mean_shift", "trend_noise", False, "slope", 0.1, "{}", 3, 512, 4.0, True
        ),
        RobustnessRow(
            "mean_shift",
            "stationary_noise",
            False,
            "control",
            1.0,
            "{}",
            4,
            512,
            0.0,
            False,
        ),
    )

    summary = summarize_confusion(rows, "mean_shift", "midpoint step", 3.0)

    assert (summary.true_positive, summary.false_negative) == (1, 1)
    assert (summary.false_positive, summary.true_negative) == (1, 1)
    assert summary.false_positive_rate == 0.5
    assert summary.false_negative_rate == 0.5


def test_writer_creates_phase43_outputs_without_phase42_paths(tmp_path: Path) -> None:
    output = tmp_path / "phase43"
    study = write_robustness_study(output)

    assert len(study.rows) == 460
    assert (output / "parameter_results.csv").is_file()
    assert (output / "confusion_summary.csv").is_file()
    assert (output / "metadata.json").is_file()
    assert (output / "mean_diagnostic_sensitivity.png").is_file()
    assert (output / "variance_diagnostic_sensitivity.png").is_file()
