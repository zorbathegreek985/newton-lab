"""Phase 45 local midpoint level-shift comparisons and held-out evaluation."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from newton_lab.exceptions import ScientificValidationError
from newton_lab.step_detection_robustness import (
    PRIMARY_RULES,
    StepCondition,
    StepDetectionStudy,
    StepRecordResult,
    generate_step_condition_signal,
    run_step_detection_study,
    select_component_threshold,
)

WINDOW_SIZES = (16, 32, 64)
LOCAL_METHODS = tuple(
    method
    for window in WINDOW_SIZES
    for method in (f"local_raw_w{window}", f"local_trend_adjusted_w{window}")
)
METHODS = PRIMARY_RULES + LOCAL_METHODS


@dataclass(frozen=True, slots=True)
class LocalThreshold:
    """Calibration-only threshold for one local score and window size."""

    method: str
    window_size: int
    absolute_difference_threshold: float
    calibration_youden_j: float


@dataclass(frozen=True, slots=True)
class LocalMethodRecord:
    """One signal/method result, retaining ground truth and Phase 44 baselines."""

    split: str
    method: str
    condition_id: str
    signal_group: str
    generator_kind: str
    target_step_present: bool
    step_magnitude: float
    trend_slope: float
    damping_rate: float
    noise_sd: float
    oscillation_amplitude: float
    seed: int
    sample_count: int
    window_size: int | None
    calibrated_threshold: float | None
    calibration_youden_j: float | None
    phase44_mean_step_z: float
    phase44_baseline_detected: bool
    phase44_trend_veto_detected: bool
    phase44_envelope_veto_detected: bool
    phase44_combined_veto_detected: bool
    local_difference: float | None
    local_score: float | None
    detected: bool


@dataclass(frozen=True, slots=True)
class LocalMethodSummary:
    """Confusion counts and paired error-rate changes by method/stratum."""

    target: str
    method: str
    stratum: str
    positive_count: int
    negative_count: int
    true_positive: int
    false_positive: int
    false_negative: int
    true_negative: int
    detection_rate: float | None
    false_positive_rate: float | None
    false_negative_rate: float | None
    false_positive_change_vs_baseline: float | None
    false_negative_change_vs_baseline: float | None


@dataclass(frozen=True, slots=True)
class LocalStepStudy:
    """Phase 44 matched records plus local methods and held-out summaries."""

    phase44_study: StepDetectionStudy
    thresholds: tuple[LocalThreshold, ...]
    records: tuple[LocalMethodRecord, ...]
    summaries: tuple[LocalMethodSummary, ...]


def local_mean_difference(values: ArrayLike, window_size: int) -> float:
    """Return mean(after midpoint window) minus mean(before midpoint window).

    For length ``n`` and window ``w``, the pre-window is ``[n//2-w:n//2]`` and
    the post-window is ``[n//2:n//2+w]``. The known midpoint is retrospective
    oracle information. This function does not search for a change location.
    """
    source = _validated_values(values, window_size)
    midpoint = len(source) // 2
    return float(
        np.mean(source[midpoint : midpoint + window_size])
        - np.mean(source[midpoint - window_size : midpoint])
    )


def global_linear_slope(values: ArrayLike) -> float:
    """Estimate a full-record OLS slope against integer sample index."""
    source = _validated_values(values, window_size=None)
    time = np.arange(len(source), dtype=np.float64)
    centered_time = time - np.mean(time)
    centered_values = source - np.mean(source)
    return float(
        np.dot(centered_time, centered_values) / np.dot(centered_time, centered_time)
    )


def trend_adjusted_local_difference(
    values: ArrayLike, window_size: int, *, slope_estimate: float | None = None
) -> float:
    """Subtract the linear-trend level change expected across local windows.

    The local window centers are ``window_size`` samples apart, so the
    adjustment is ``slope_estimate * window_size``. Unless supplied, slope is
    estimated by full-record OLS and can therefore be biased by the step.
    """
    difference = local_mean_difference(values, window_size)
    slope = global_linear_slope(values) if slope_estimate is None else slope_estimate
    if not np.isfinite(slope):
        raise ScientificValidationError("slope_estimate must be finite")
    return float(difference - slope * window_size)


def run_local_step_study(
    *, calibration_replicates: int = 40, evaluation_replicates: int = 100
) -> LocalStepStudy:
    """Reuse the Phase 44 seed split and baseline, calibrating local thresholds.

    All three predeclared windows (16, 32, 64 samples per side) are evaluated;
    none is selected based on evaluation results. Each local rule threshold is
    calibrated on Phase 44 calibration records only.
    """
    phase44 = run_step_detection_study(
        calibration_replicates=calibration_replicates,
        evaluation_replicates=evaluation_replicates,
    )
    condition_by_id = {
        condition.condition_id: condition for condition in phase44.conditions
    }
    source_rows = tuple(
        _local_scores(record, condition_by_id[record.condition_id])
        for record in phase44.records
    )
    calibration_rows = tuple(
        result for result in source_rows if result[0].split == "calibration"
    )
    thresholds: list[LocalThreshold] = []
    for window in WINDOW_SIZES:
        for method in (f"local_raw_w{window}", f"local_trend_adjusted_w{window}"):
            calibration_scores = (
                tuple(abs(row[2][window]) for row in calibration_rows)
                if method.startswith("local_raw")
                else tuple(abs(row[3][window]) for row in calibration_rows)
            )
            threshold, youden_j = select_component_threshold(
                calibration_scores,
                tuple(row[0].target_step_present for row in calibration_rows),
            )
            thresholds.append(LocalThreshold(method, window, threshold, youden_j))
    threshold_by_method = {item.method: item for item in thresholds}
    records: list[LocalMethodRecord] = []
    for phase44_record, condition, raw, adjusted in source_rows:
        base = phase44_record
        for method in PRIMARY_RULES:
            detected = bool(
                getattr(
                    phase44_record,
                    {
                        "baseline": "baseline_detected",
                        "trend_veto": "trend_veto_detected",
                        "envelope_veto": "envelope_veto_detected",
                        "combined_veto": "combined_veto_detected",
                    }[method],
                )
            )
            records.append(
                _record(
                    base,
                    method,
                    condition,
                    phase44_record,
                    None,
                    None,
                    None,
                    None,
                    detected,
                )
            )
        for method in LOCAL_METHODS:
            window = int(method.rsplit("w", maxsplit=1)[1])
            threshold_record = threshold_by_method[method]
            difference = (
                raw[window] if method.startswith("local_raw") else adjusted[window]
            )
            score = abs(difference)
            records.append(
                _record(
                    base,
                    method,
                    condition,
                    phase44_record,
                    window,
                    threshold_record.absolute_difference_threshold,
                    threshold_record.calibration_youden_j,
                    difference,
                    score >= threshold_record.absolute_difference_threshold,
                    score,
                )
            )
    record_tuple = tuple(records)
    summaries = _summarize(record_tuple)
    return LocalStepStudy(phase44, tuple(thresholds), record_tuple, summaries)


def write_local_step_study(
    output_directory: str | Path = "reports/phase_45_local_step_detection",
    *,
    calibration_replicates: int = 40,
    evaluation_replicates: int = 100,
) -> LocalStepStudy:
    """Write the Phase 45 results and plots to a dedicated directory."""
    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    study = run_local_step_study(
        calibration_replicates=calibration_replicates,
        evaluation_replicates=evaluation_replicates,
    )
    _write_csv(destination / "per_record_results.csv", study.records)
    _write_csv(destination / "summary.csv", study.summaries)
    metadata = {
        "experiment": "Phase 45 local level-shift detection with trend adjustment",
        "data": "synthetic only",
        "sample_count": 512,
        "calibration_replicates_per_condition": (
            study.phase44_study.calibration_replicates
        ),
        "evaluation_replicates_per_condition": (
            study.phase44_study.evaluation_replicates
        ),
        "calibration_seeds": list(study.phase44_study.calibration_seeds),
        "evaluation_seeds": list(study.phase44_study.evaluation_seeds),
        "seed_policy": (
            "reuses Phase 44 seed sets; same record is paired across methods"
        ),
        "phase44_baselines": {
            "baseline": "abs(two_half_mean_z) >= 3",
            "trend_veto": "Phase 44 baseline AND NOT calibrated global trend flag",
            "envelope_veto": "Phase 44 baseline AND NOT calibrated envelope flag",
            "combined_veto": "Phase 44 baseline AND neither Phase 44 confound flag",
            "trend_abs_slope_threshold": (
                study.phase44_study.thresholds.trend_abs_slope
            ),
            "envelope_decay_threshold": (
                study.phase44_study.thresholds.envelope_decay_fraction
            ),
        },
        "local_window_sizes_per_side": list(WINDOW_SIZES),
        "local_midpoint": "known oracle location at sample_count // 2",
        "local_rules": [asdict(threshold) for threshold in study.thresholds],
        "trend_adjustment": "local post-minus-pre mean difference minus "
        "Phase 44 full-record OLS slope times window size",
        "window_selection": (
            "all predeclared sizes reported; no held-out window selection"
        ),
        "conditions": [
            asdict(condition) for condition in study.phase44_study.conditions
        ],
    }
    (destination / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    _write_plots(study, destination)
    return study


def _local_scores(
    record: StepRecordResult, condition: StepCondition
) -> tuple[
    StepRecordResult,
    StepCondition,
    dict[int, float],
    dict[int, float],
]:
    signal = generate_step_condition_signal(condition, record.seed)
    raw = {
        window: local_mean_difference(signal.values, window) for window in WINDOW_SIZES
    }
    slope = record.estimated_global_slope
    adjusted = {window: raw[window] - slope * window for window in WINDOW_SIZES}
    return record, condition, raw, adjusted


def _record(
    base: StepRecordResult,
    method: str,
    condition: StepCondition,
    phase44_record: StepRecordResult,
    window: int | None,
    threshold: float | None,
    youden_j: float | None,
    difference: float | None,
    detected: bool,
    score: float | None = None,
) -> LocalMethodRecord:
    return LocalMethodRecord(
        base.split,
        method,
        base.condition_id,
        base.signal_group,
        base.generator_kind,
        base.target_step_present,
        base.step_magnitude,
        base.trend_slope,
        base.damping_rate,
        base.noise_sd,
        base.oscillation_amplitude,
        base.seed,
        base.sample_count,
        window,
        threshold,
        youden_j,
        phase44_record.mean_step_z,
        phase44_record.baseline_detected,
        phase44_record.trend_veto_detected,
        phase44_record.envelope_veto_detected,
        phase44_record.combined_veto_detected,
        difference,
        score,
        detected,
    )


def _summarize(
    records: tuple[LocalMethodRecord, ...],
) -> tuple[LocalMethodSummary, ...]:
    evaluation = tuple(record for record in records if record.split == "evaluation")
    strata: dict[str, tuple[LocalMethodRecord, ...]] = {
        "overall": evaluation,
        **{
            group: tuple(row for row in evaluation if row.signal_group == group)
            for group in dict.fromkeys(row.signal_group for row in evaluation)
        },
        **{
            f"condition:{condition_id}": tuple(
                row for row in evaluation if row.condition_id == condition_id
            )
            for condition_id in dict.fromkeys(row.condition_id for row in evaluation)
        },
    }
    step_sizes = sorted(
        {row.step_magnitude for row in evaluation if row.target_step_present}
    )
    strata.update(
        {
            f"step_magnitude={size:g}": tuple(
                row
                for row in evaluation
                if row.target_step_present and row.step_magnitude == size
            )
            for size in step_sizes
        }
    )
    results: list[LocalMethodSummary] = []
    for stratum, subset in strata.items():
        baseline = tuple(row for row in subset if row.method == "baseline")
        baseline_counts = _confusion(baseline)
        for method in METHODS:
            selected = tuple(row for row in subset if row.method == method)
            counts = _confusion(selected)
            tp, fp, fn, tn = counts
            pos, neg = tp + fn, fp + tn
            fpr = fp / neg if neg else None
            fnr = fn / pos if pos else None
            if method == "baseline" or not baseline:
                fpr_change = fnr_change = None
            else:
                btp, bfp, bfn, btn = baseline_counts
                baseline_fpr = bfp / (bfp + btn) if bfp + btn else None
                baseline_fnr = bfn / (btp + bfn) if btp + bfn else None
                fpr_change = (
                    None if fpr is None or baseline_fpr is None else fpr - baseline_fpr
                )
                fnr_change = (
                    None if fnr is None or baseline_fnr is None else fnr - baseline_fnr
                )
            results.append(
                LocalMethodSummary(
                    "abrupt_midpoint_mean_step",
                    method,
                    stratum,
                    pos,
                    neg,
                    tp,
                    fp,
                    fn,
                    tn,
                    tp / pos if pos else None,
                    fpr,
                    fnr,
                    fpr_change,
                    fnr_change,
                )
            )
    return tuple(results)


def _confusion(
    rows: tuple[LocalMethodRecord, ...],
) -> tuple[int, int, int, int]:
    tp = sum(row.target_step_present and row.detected for row in rows)
    fp = sum(not row.target_step_present and row.detected for row in rows)
    fn = sum(row.target_step_present and not row.detected for row in rows)
    tn = sum(not row.target_step_present and not row.detected for row in rows)
    return tp, fp, fn, tn


def _validated_values(
    values: ArrayLike, window_size: int | None
) -> NDArray[np.float64]:
    try:
        source = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ScientificValidationError("values must be finite numeric data") from error
    if source.ndim != 1 or source.size < 4 or not np.isfinite(source).all():
        raise ScientificValidationError(
            "values must be a finite 1-D array of length >= 4"
        )
    if window_size is not None and (
        isinstance(window_size, bool)
        or not isinstance(window_size, int)
        or window_size < 1
        or 2 * window_size > len(source)
    ):
        raise ScientificValidationError("window_size must be an integer in [1, n//2]")
    return source


def _write_csv(path: Path, rows: tuple[Any, ...]) -> None:
    if not rows:
        raise ScientificValidationError("cannot write an empty Phase 45 table")
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def _write_plots(study: LocalStepStudy, destination: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    overall = [
        row
        for row in study.summaries
        if row.stratum == "overall" and row.method in METHODS
    ]
    figure, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    axes[0].bar(
        [row.method for row in overall],
        [float(row.false_positive_rate or 0.0) for row in overall],
    )
    axes[0].set_title("Held-out false-positive rate")
    axes[0].set_ylim(0.0, 1.0)
    axes[1].bar(
        [row.method for row in overall],
        [float(row.false_negative_rate or 0.0) for row in overall],
    )
    axes[1].set_title("Held-out false-negative rate")
    axes[1].set_ylim(0.0, 1.0)
    for axis in axes:
        axis.tick_params(axis="x", labelrotation=70, labelsize=7)
    figure.savefig(destination / "method_error_rates.png", dpi=140)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(8, 4), constrained_layout=True)
    step_summaries = [
        row for row in study.summaries if row.stratum.startswith("step_magnitude=")
    ]
    magnitudes = sorted(
        {float(row.stratum.split("=", maxsplit=1)[1]) for row in step_summaries}
    )
    plotted_methods = (
        "baseline",
        "trend_veto",
        "local_raw_w16",
        "local_trend_adjusted_w16",
        "local_raw_w32",
        "local_trend_adjusted_w32",
        "local_raw_w64",
        "local_trend_adjusted_w64",
    )
    for method in plotted_methods:
        axis.plot(
            magnitudes,
            [
                _summary_rate(step_summaries, method, magnitude)
                for magnitude in magnitudes
            ],
            marker="o",
            label=method,
        )
    axis.set_xlabel("Abrupt step magnitude")
    axis.set_ylabel("Detection rate")
    axis.set_ylim(0.0, 1.0)
    axis.legend(fontsize=7)
    figure.savefig(destination / "step_magnitude_sensitivity.png", dpi=140)
    plt.close(figure)


def _summary_rate(
    summaries: list[LocalMethodSummary], method: str, magnitude: float
) -> float:
    summary = next(
        row
        for row in summaries
        if row.method == method and row.stratum == f"step_magnitude={magnitude:g}"
    )
    if summary.detection_rate is None:
        raise ScientificValidationError("step magnitude detection rate unavailable")
    return summary.detection_rate


def main() -> None:
    """Run the predeclared Phase 45 study and write its dedicated outputs."""
    study = write_local_step_study()
    print(
        f"Wrote {len(study.records)} method/record rows and "
        f"{len(study.summaries)} held-out summaries"
    )


if __name__ == "__main__":
    main()
