"""Controlled Phase 43 robustness and ambiguity studies for synthetic signals."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np

from newton_lab.exceptions import ScientificValidationError
from newton_lab.synthetic_timeseries_benchmark import (
    SignalKind,
    SyntheticSignal,
    generate_signal,
    measure_signal,
)

DiagnosticName = Literal["mean_shift", "variance_shift"]
REPLICATES = 20
SAMPLE_COUNT = 512
MEAN_Z_THRESHOLD = 3.0
VARIANCE_RATIO_THRESHOLD = 2.0


@dataclass(frozen=True, slots=True)
class RobustnessRow:
    """One generated record and the diagnostic classification it received."""

    diagnostic: DiagnosticName
    process: str
    target_present: bool
    parameter_name: str
    parameter_value: float
    fixed_settings: str
    seed: int
    sample_count: int
    measured_statistic: float
    detected: bool


@dataclass(frozen=True, slots=True)
class ConfusionSummary:
    """Binary confusion counts for one narrowly defined target property."""

    diagnostic: DiagnosticName
    positive_definition: str
    threshold: float
    positive_count: int
    negative_count: int
    true_positive: int
    false_positive: int
    false_negative: int
    true_negative: int
    false_positive_rate: float
    false_negative_rate: float


@dataclass(frozen=True, slots=True)
class RobustnessStudy:
    """Reproducible parameter-level rows and pooled binary summaries."""

    rows: tuple[RobustnessRow, ...]
    summaries: tuple[ConfusionSummary, ...]


def run_robustness_study() -> RobustnessStudy:
    """Run fixed-seed mean-step and variance-step robustness experiments.

    Twenty seeds are used for every configuration. The mean experiment has a
    positive class of midpoint step signals (magnitudes 0.1, 0.2, and 0.5) and
    negative controls comprising iid noise, linear trends, and AR(1). The
    variance experiment's positive class is a random variance step (post/pre
    SD multiplier greater than one); negative controls include equal-variance
    noise and damped sinusoids with varying damping/noise.
    """
    rows: list[RobustnessRow] = []
    mean_configurations: list[
        tuple[str, SignalKind, float, dict[str, float], bool]
    ] = []
    mean_configurations.extend(
        (
            "mean_shift",
            "mean_shift",
            magnitude,
            {"mean_displacement": magnitude, "noise_sd": 1.0},
            True,
        )
        for magnitude in (0.1, 0.2, 0.5)
    )
    mean_configurations.append(
        ("stationary_noise", "stationary_noise", 1.0, {"noise_sd": 1.0}, False)
    )
    mean_configurations.extend(
        (
            "trend_noise",
            "trend_noise",
            slope,
            {"noise_sd": 0.75, "trend_slope": slope},
            False,
        )
        for slope in (0.0, 0.0005, 0.001, 0.002, 0.004)
    )
    mean_configurations.extend(
        (
            "trend_noise",
            "trend_noise",
            noise_sd,
            {"trend_slope": 0.002, "noise_sd": noise_sd},
            False,
        )
        for noise_sd in (0.25, 1.5)
    )
    mean_configurations.append(
        (
            "autocorrelated",
            "autocorrelated",
            0.5,
            {"phi": 0.8, "observation_noise_sd": 0.5},
            False,
        )
    )
    for configuration_index, (process, kind, value, settings, positive) in enumerate(
        mean_configurations
    ):
        for replicate in range(REPLICATES):
            seed = 43000 + configuration_index * 100 + replicate
            signal = _generate(kind, value, settings, seed)
            measured = measure_signal(signal)
            rows.append(
                RobustnessRow(
                    "mean_shift",
                    process,
                    positive,
                    "mean_displacement"
                    if positive
                    else "trend_slope"
                    if process == "trend_noise" and settings["noise_sd"] == 0.75
                    else "noise_sd"
                    if process == "trend_noise"
                    else "control",
                    value,
                    json.dumps(settings, sort_keys=True),
                    seed,
                    SAMPLE_COUNT,
                    measured.mean_shift_z,
                    measured.mean_shift_flag,
                )
            )

    variance_configurations: list[
        tuple[str, SignalKind, float, dict[str, float], bool, str]
    ] = []
    variance_configurations.extend(
        (
            "variance_shift",
            "variance_shift",
            multiplier,
            {"post_change_sd_multiplier": multiplier},
            multiplier > 1.0,
            "post_change_sd_multiplier",
        )
        for multiplier in (1.0, 1.1, 1.4, 2.0)
    )
    variance_configurations.extend(
        (
            "damped_oscillation",
            "damped_oscillation",
            damping,
            {
                "noise_sd": 0.2,
                "damping_rate": damping,
                "amplitude": 1.0,
                "frequency_cycles_per_sample": 0.10,
            },
            False,
            "damping_rate",
        )
        for damping in (0.0, 0.001, 0.002, 0.004)
    )
    variance_configurations.extend(
        (
            "damped_oscillation",
            "damped_oscillation",
            noise_sd,
            {
                "noise_sd": noise_sd,
                "damping_rate": 0.002,
                "amplitude": 1.0,
                "frequency_cycles_per_sample": 0.10,
            },
            False,
            "noise_sd",
        )
        for noise_sd in (0.05, 0.5)
    )
    variance_configurations.append(
        (
            "stationary_noise",
            "stationary_noise",
            1.0,
            {"noise_sd": 1.0},
            False,
            "noise_sd",
        )
    )
    for configuration_index, (
        process,
        kind,
        value,
        settings,
        positive,
        parameter_name,
    ) in enumerate(variance_configurations):
        for replicate in range(REPLICATES):
            seed = 53000 + configuration_index * 100 + replicate
            signal = _generate(kind, value, settings, seed)
            measured = measure_signal(signal)
            rows.append(
                RobustnessRow(
                    "variance_shift",
                    process,
                    positive,
                    parameter_name,
                    value,
                    json.dumps(settings, sort_keys=True),
                    seed,
                    SAMPLE_COUNT,
                    measured.variance_ratio,
                    measured.variance_shift_flag,
                )
            )

    row_tuple = tuple(rows)
    return RobustnessStudy(
        row_tuple,
        (
            summarize_confusion(
                row_tuple,
                "mean_shift",
                "A nonzero abrupt midpoint mean step; trends and AR dependence "
                "are controls.",
                MEAN_Z_THRESHOLD,
            ),
            summarize_confusion(
                row_tuple,
                "variance_shift",
                "A random midpoint variance step with post/pre SD multiplier > 1.",
                VARIANCE_RATIO_THRESHOLD,
            ),
        ),
    )


def write_robustness_study(
    output_directory: str | Path = "reports/phase_43_diagnostic_robustness",
) -> RobustnessStudy:
    """Write parameter rows, confusion summaries, metadata, and two figures."""
    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    study = run_robustness_study()
    _write_csv(destination / "parameter_results.csv", study.rows)
    _write_csv(destination / "confusion_summary.csv", study.summaries)
    metadata = {
        "experiment": "Phase 43 diagnostic robustness and signal ambiguity",
        "data": "synthetic only",
        "sample_count": SAMPLE_COUNT,
        "replicates_per_configuration": REPLICATES,
        "mean_seed_rule": "43000 + configuration_index * 100 + replicate_index",
        "variance_seed_rule": "53000 + configuration_index * 100 + replicate_index",
        "thresholds": {
            "absolute_two_half_mean_z": MEAN_Z_THRESHOLD,
            "two_half_variance_ratio": VARIANCE_RATIO_THRESHOLD,
        },
        "fixed_oscillation_frequency_cycles_per_sample": 0.10,
        "fixed_oscillation_amplitude": 1.0,
        "fixed_ar_phi": 0.8,
        "result_rows": len(study.rows),
    }
    (destination / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    _write_plots(study.rows, destination)
    return study


def _generate(
    kind: SignalKind, value: float, settings: dict[str, float], seed: int
) -> SyntheticSignal:
    if kind == "trend_noise":
        return generate_signal(
            kind,
            seed=seed,
            parameter_value=settings["noise_sd"],
            trend_slope=settings["trend_slope"],
        )
    if kind == "damped_oscillation":
        return generate_signal(
            kind,
            seed=seed,
            parameter_value=settings["noise_sd"],
            damping_rate=settings["damping_rate"],
        )
    return generate_signal(kind, seed=seed, parameter_value=value)


def summarize_confusion(
    rows: tuple[RobustnessRow, ...],
    diagnostic: DiagnosticName,
    definition: str,
    threshold: float,
) -> ConfusionSummary:
    """Calculate confusion counts and class-conditional error rates."""
    selected = tuple(row for row in rows if row.diagnostic == diagnostic)
    tp = sum(row.target_present and row.detected for row in selected)
    fn = sum(row.target_present and not row.detected for row in selected)
    fp = sum(not row.target_present and row.detected for row in selected)
    tn = sum(not row.target_present and not row.detected for row in selected)
    positive_count, negative_count = tp + fn, fp + tn
    if positive_count == 0 or negative_count == 0:
        raise ScientificValidationError("confusion summaries require both classes")
    return ConfusionSummary(
        diagnostic,
        definition,
        threshold,
        positive_count,
        negative_count,
        tp,
        fp,
        fn,
        tn,
        fp / negative_count,
        fn / positive_count,
    )


def _write_csv(
    path: Path, rows: tuple[RobustnessRow, ...] | tuple[ConfusionSummary, ...]
) -> None:
    if not rows:
        raise ScientificValidationError("cannot write an empty study table")
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def _write_plots(rows: tuple[RobustnessRow, ...], destination: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mean_rows = tuple(row for row in rows if row.diagnostic == "mean_shift")
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    _plot_rate(
        axes[0],
        mean_rows,
        "trend_noise",
        "trend_slope",
        "Mean diagnostic vs trend slope",
    )
    _plot_rate(
        axes[1],
        mean_rows,
        "mean_shift",
        "mean_displacement",
        "Mean diagnostic vs step size",
    )
    for axis in axes:
        axis.set_ylabel("Mean-shift flag rate")
        axis.set_ylim(0.0, 1.0)
    figure.savefig(destination / "mean_diagnostic_sensitivity.png", dpi=140)
    plt.close(figure)

    variance_rows = tuple(row for row in rows if row.diagnostic == "variance_shift")
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    _plot_rate(
        axes[0],
        variance_rows,
        "damped_oscillation",
        "damping_rate",
        "Variance diagnostic vs damping",
    )
    _plot_rate(
        axes[1],
        variance_rows,
        "variance_shift",
        "post_change_sd_multiplier",
        "Variance diagnostic vs SD multiplier",
    )
    for axis in axes:
        axis.set_ylabel("Variance-shift flag rate")
        axis.set_ylim(0.0, 1.0)
    figure.savefig(destination / "variance_diagnostic_sensitivity.png", dpi=140)
    plt.close(figure)


def _plot_rate(
    axis: Any,
    rows: tuple[RobustnessRow, ...],
    process: str,
    parameter_name: str,
    title: str,
) -> None:
    selected = [
        row
        for row in rows
        if row.process == process and row.parameter_name == parameter_name
    ]
    settings = sorted({row.parameter_value for row in selected})
    rates = [
        np.mean([row.detected for row in selected if row.parameter_value == setting])
        for setting in settings
    ]
    axis.plot(settings, rates, marker="o")
    axis.set_title(title)
    axis.set_xlabel(parameter_name)


def main() -> None:
    """Run the fixed robustness study and write its dedicated outputs."""
    study = write_robustness_study()
    print(f"Wrote {len(study.rows)} synthetic robustness runs")


if __name__ == "__main__":
    main()
