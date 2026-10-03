"""Phase 44 synthetic evaluation of veto checks for abrupt mean-step detection."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

from newton_lab.exceptions import ScientificValidationError
from newton_lab.synthetic_timeseries_benchmark import (
    SignalKind,
    SyntheticSignal,
    generate_signal,
    measure_signal,
)

SignalGroup = Literal[
    "stationary_noise",
    "gradual_trend",
    "damped_oscillation",
    "mean_step",
    "step_plus_trend",
    "step_plus_oscillation",
]
PRIMARY_RULES = ("baseline", "trend_veto", "envelope_veto", "combined_veto")
CALIBRATION_SEED_BASE = 44_000
EVALUATION_SEED_BASE = 100_000
SAMPLE_COUNT = 512
OSCILLATION_FREQUENCY_CYCLES_PER_SAMPLE = 0.10


@dataclass(frozen=True, slots=True)
class StepCondition:
    """One declared signal configuration and its structural ground truth."""

    condition_id: str
    group: SignalGroup
    generator_kind: SignalKind
    target_step_present: bool
    trend_component_present: bool
    damped_component_present: bool
    step_magnitude: float
    trend_slope: float
    damping_rate: float
    noise_sd: float
    oscillation_amplitude: float


@dataclass(frozen=True, slots=True)
class FrozenThresholds:
    """Candidate thresholds selected using calibration records only."""

    trend_abs_slope: float
    envelope_decay_fraction: float
    trend_calibration_youden_j: float
    envelope_calibration_youden_j: float
    trend_tie_break: str
    envelope_tie_break: str


@dataclass(frozen=True, slots=True)
class StepRecordResult:
    """One generated record, all diagnostics, and all step-rule decisions."""

    split: Literal["calibration", "evaluation"]
    condition_id: str
    signal_group: str
    generator_kind: str
    seed: int
    sample_count: int
    target_step_present: bool
    trend_component_present: bool
    damped_component_present: bool
    step_magnitude: float
    trend_slope: float
    damping_rate: float
    noise_sd: float
    oscillation_amplitude: float
    mean_step_z: float
    baseline_detected: bool
    estimated_global_slope: float
    trend_flag: bool
    envelope_decay_fraction: float
    envelope_flag: bool
    trend_veto_detected: bool
    envelope_veto_detected: bool
    combined_veto_detected: bool


@dataclass(frozen=True, slots=True)
class RuleSummary:
    """Confusion counts and error rates for one target, rule, and stratum."""

    target: str
    rule: str
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
class StepDetectionStudy:
    """Complete calibration and held-out evaluation results."""

    conditions: tuple[StepCondition, ...]
    thresholds: FrozenThresholds
    calibration_replicates: int
    evaluation_replicates: int
    calibration_seeds: tuple[int, ...]
    evaluation_seeds: tuple[int, ...]
    records: tuple[StepRecordResult, ...]
    summaries: tuple[RuleSummary, ...]


def default_conditions() -> tuple[StepCondition, ...]:
    """Return the fixed, compact Phase 44 condition grid."""
    conditions: list[StepCondition] = [
        StepCondition(
            "stationary_null",
            "stationary_noise",
            "stationary_noise",
            False,
            False,
            False,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
        )
    ]
    conditions.extend(
        StepCondition(
            f"trend_{slope:g}",
            "gradual_trend",
            "trend_noise",
            False,
            slope > 0.0,
            False,
            0.0,
            slope,
            0.0,
            0.75,
            0.0,
        )
        for slope in (0.0, 0.0005, 0.001, 0.002, 0.004)
    )
    conditions.extend(
        StepCondition(
            f"damped_{damping:g}",
            "damped_oscillation",
            "damped_oscillation",
            False,
            False,
            damping > 0.0,
            0.0,
            0.0,
            damping,
            0.2,
            1.0,
        )
        for damping in (0.0, 0.001, 0.002, 0.004)
    )
    conditions.extend(
        StepCondition(
            f"mean_step_{magnitude:g}",
            "mean_step",
            "mean_shift",
            True,
            False,
            False,
            magnitude,
            0.0,
            0.0,
            1.0,
            0.0,
        )
        for magnitude in (0.1, 0.2, 0.5, 1.0)
    )
    conditions.extend(
        (
            StepCondition(
                "weak_step_plus_trend",
                "step_plus_trend",
                "trend_noise",
                True,
                True,
                False,
                0.2,
                0.001,
                0.0,
                1.0,
                0.0,
            ),
            StepCondition(
                "weak_step_plus_damping",
                "step_plus_oscillation",
                "damped_oscillation",
                True,
                False,
                True,
                0.2,
                0.0,
                0.002,
                1.0,
                1.0,
            ),
        )
    )
    return tuple(conditions)


def run_step_detection_study(
    *, calibration_replicates: int = 40, evaluation_replicates: int = 100
) -> StepDetectionStudy:
    """Calibrate confounder checks on disjoint seeds, freeze, then evaluate.

    The trend check uses the absolute global ordinary-least-squares slope.
    The envelope check uses the fractional reduction between sinusoidal
    amplitudes fitted separately on the two known half-windows at the
    generator's known frequency (0.10 cycles/sample). Both thresholds maximize
    Youden's J for the corresponding component-presence label on calibration
    records only; exact ties select the largest threshold to favor fewer vetoes.
    The baseline remains the Phase 42 absolute two-half mean statistic >= 3.
    Veto rules preserve its detection only when the named checks do not flag.
    """
    _validate_replicates(calibration_replicates, "calibration_replicates")
    _validate_replicates(evaluation_replicates, "evaluation_replicates")
    conditions = default_conditions()
    calibration = tuple(
        _make_record(
            "calibration",
            condition,
            _seed(CALIBRATION_SEED_BASE, index, replicate),
            FrozenThresholds(0.0, 0.0, 0.0, 0.0, "", ""),
        )
        for index, condition in enumerate(conditions)
        for replicate in range(calibration_replicates)
    )
    trend_threshold, trend_j = select_component_threshold(
        tuple(record.estimated_global_slope for record in calibration),
        tuple(record.trend_component_present for record in calibration),
    )
    envelope_threshold, envelope_j = select_component_threshold(
        tuple(record.envelope_decay_fraction for record in calibration),
        tuple(record.damped_component_present for record in calibration),
    )
    thresholds = FrozenThresholds(
        trend_threshold,
        envelope_threshold,
        trend_j,
        envelope_j,
        "largest threshold among equal Youden-J maxima",
        "largest threshold among equal Youden-J maxima",
    )
    calibration_records = tuple(
        _make_record(
            "calibration",
            condition,
            _seed(CALIBRATION_SEED_BASE, index, replicate),
            thresholds,
        )
        for index, condition in enumerate(conditions)
        for replicate in range(calibration_replicates)
    )
    evaluation_records = tuple(
        _make_record(
            "evaluation",
            condition,
            _seed(EVALUATION_SEED_BASE, index, replicate),
            thresholds,
        )
        for index, condition in enumerate(conditions)
        for replicate in range(evaluation_replicates)
    )
    records = calibration_records + evaluation_records
    calibration_seeds = tuple(record.seed for record in calibration_records)
    evaluation_seeds = tuple(record.seed for record in evaluation_records)
    if set(calibration_seeds) & set(evaluation_seeds):
        raise ScientificValidationError("calibration and evaluation seeds overlap")
    summaries = _summarize_evaluation(evaluation_records)
    return StepDetectionStudy(
        conditions,
        thresholds,
        calibration_replicates,
        evaluation_replicates,
        calibration_seeds,
        evaluation_seeds,
        records,
        summaries,
    )


def write_step_detection_study(
    output_directory: str | Path = "reports/phase_44_step_detection",
    *,
    calibration_replicates: int = 40,
    evaluation_replicates: int = 100,
) -> StepDetectionStudy:
    """Write per-record results, evaluation summaries, metadata, and plots."""
    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    study = run_step_detection_study(
        calibration_replicates=calibration_replicates,
        evaluation_replicates=evaluation_replicates,
    )
    _write_rows(destination / "per_record_results.csv", study.records)
    _write_rows(destination / "summary.csv", study.summaries)
    metadata = {
        "experiment": "Phase 44 robust step detection under confounding signals",
        "data": "synthetic only",
        "sample_count": SAMPLE_COUNT,
        "calibration_replicates_per_configuration": calibration_replicates,
        "evaluation_replicates_per_configuration": evaluation_replicates,
        "calibration_seeds": list(study.calibration_seeds),
        "evaluation_seeds": list(study.evaluation_seeds),
        "threshold_selection": "maximize Youden J on calibration seeds only; "
        "ties select the largest threshold",
        "frozen_thresholds": asdict(study.thresholds),
        "baseline": "Phase 42 abs(two_half_mean_z) >= 3",
        "trend_check": "abs(OLS global slope) >= frozen trend threshold",
        "envelope_check": "fractional early-to-late sinusoid amplitude reduction "
        ">= frozen envelope threshold; known frequency 0.10 cycles/sample",
        "veto_rules": {
            "trend_veto": "baseline AND NOT trend_flag",
            "envelope_veto": "baseline AND NOT envelope_flag",
            "combined_veto": "baseline AND NOT trend_flag AND NOT envelope_flag",
        },
        "conditions": [asdict(condition) for condition in study.conditions],
    }
    (destination / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    _write_plots(study, destination)
    return study


def _validate_replicates(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 2:
        raise ScientificValidationError(f"{name} must be an integer >= 2")


def _seed(base: int, condition_index: int, replicate: int) -> int:
    return base + condition_index * 1000 + replicate


def generate_step_condition_signal(
    condition: StepCondition, seed: int
) -> SyntheticSignal:
    """Generate one declared Phase 44 condition for reuse by later studies."""
    if condition.generator_kind == "stationary_noise":
        return generate_signal("stationary_noise", seed=seed)
    if condition.generator_kind == "trend_noise":
        return generate_signal(
            "trend_noise",
            seed=seed,
            parameter_value=condition.noise_sd,
            trend_slope=condition.trend_slope,
            additional_mean_step=condition.step_magnitude,
        )
    if condition.generator_kind == "damped_oscillation":
        return generate_signal(
            "damped_oscillation",
            seed=seed,
            parameter_value=condition.noise_sd,
            oscillation_amplitude=condition.oscillation_amplitude,
            damping_rate=condition.damping_rate,
            additional_mean_step=condition.step_magnitude,
        )
    if condition.generator_kind == "mean_shift":
        return generate_signal(
            "mean_shift", seed=seed, parameter_value=condition.step_magnitude
        )
    raise ScientificValidationError("unsupported Phase 44 signal condition")


def _make_record(
    split: Literal["calibration", "evaluation"],
    condition: StepCondition,
    seed: int,
    thresholds: FrozenThresholds,
) -> StepRecordResult:
    signal = generate_step_condition_signal(condition, seed)
    measured = measure_signal(signal)
    slope = _global_slope(signal.values)
    envelope_decay = _envelope_decay_fraction(signal.values)
    trend_flag = abs(slope) >= thresholds.trend_abs_slope
    envelope_flag = envelope_decay >= thresholds.envelope_decay_fraction
    baseline = measured.mean_shift_flag
    return StepRecordResult(
        split,
        condition.condition_id,
        condition.group,
        condition.generator_kind,
        seed,
        SAMPLE_COUNT,
        condition.target_step_present,
        condition.trend_component_present,
        condition.damped_component_present,
        condition.step_magnitude,
        condition.trend_slope,
        condition.damping_rate,
        condition.noise_sd,
        condition.oscillation_amplitude,
        measured.mean_shift_z,
        baseline,
        slope,
        trend_flag,
        envelope_decay,
        envelope_flag,
        baseline and not trend_flag,
        baseline and not envelope_flag,
        baseline and not trend_flag and not envelope_flag,
    )


def _global_slope(values: NDArray[np.float64]) -> float:
    time = np.arange(len(values), dtype=np.float64)
    centered_time = time - np.mean(time)
    centered_values = values - np.mean(values)
    return float(
        np.dot(centered_time, centered_values) / np.dot(centered_time, centered_time)
    )


def _envelope_decay_fraction(values: NDArray[np.float64]) -> float:
    """Measure fitted sinusoid amplitude loss between midpoint half-windows."""
    midpoint = len(values) // 2
    indices = np.arange(len(values), dtype=np.float64)
    amplitudes = tuple(
        _sinusoid_amplitude(indices[start:end], values[start:end])
        for start, end in ((0, midpoint), (midpoint, len(values)))
    )
    early, late = amplitudes
    return float(1.0 - late / max(early, np.finfo(np.float64).eps))


def _sinusoid_amplitude(
    time: NDArray[np.float64], values: NDArray[np.float64]
) -> float:
    omega = 2.0 * np.pi * OSCILLATION_FREQUENCY_CYCLES_PER_SAMPLE
    design = np.column_stack(
        (np.ones(len(time)), np.sin(omega * time), np.cos(omega * time))
    )
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return float(np.hypot(coefficients[1], coefficients[2]))


def select_component_threshold(
    scores: tuple[float, ...], positives: tuple[bool, ...]
) -> tuple[float, float]:
    """Select score threshold and Youden J using labeled calibration values."""
    if len(scores) != len(positives) or not any(positives) or all(positives):
        raise ScientificValidationError("threshold calibration requires both classes")
    ordered = sorted(set(scores))
    candidates = [float(np.nextafter(ordered[0], -np.inf))]
    candidates.extend(
        (left + right) / 2.0 for left, right in zip(ordered, ordered[1:], strict=False)
    )
    candidates.append(float(np.nextafter(ordered[-1], np.inf)))
    best_threshold = candidates[0]
    best_j = -np.inf
    for threshold in candidates:
        predicted = tuple(score >= threshold for score in scores)
        tp = sum(
            actual and result
            for actual, result in zip(positives, predicted, strict=True)
        )
        fn = sum(
            actual and not result
            for actual, result in zip(positives, predicted, strict=True)
        )
        fp = sum(
            not actual and result
            for actual, result in zip(positives, predicted, strict=True)
        )
        tn = sum(
            not actual and not result
            for actual, result in zip(positives, predicted, strict=True)
        )
        sensitivity = tp / (tp + fn)
        specificity = tn / (tn + fp)
        youden_j = sensitivity + specificity - 1.0
        if youden_j > best_j or (
            np.isclose(youden_j, best_j) and threshold > best_threshold
        ):
            best_threshold = threshold
            best_j = youden_j
    return best_threshold, float(best_j)


def _summarize_evaluation(
    records: tuple[StepRecordResult, ...],
) -> tuple[RuleSummary, ...]:
    rows: list[RuleSummary] = []
    strata: dict[str, tuple[StepRecordResult, ...]] = {
        "overall": records,
        **{
            group: tuple(record for record in records if record.signal_group == group)
            for group in dict.fromkeys(record.signal_group for record in records)
        },
        **{
            f"condition:{condition_id}": tuple(
                record for record in records if record.condition_id == condition_id
            )
            for condition_id in dict.fromkeys(record.condition_id for record in records)
        },
    }
    step_sizes = sorted(
        {record.step_magnitude for record in records if record.target_step_present}
    )
    strata.update(
        {
            f"step_magnitude={step:g}": tuple(
                record
                for record in records
                if record.target_step_present and record.step_magnitude == step
            )
            for step in step_sizes
        }
    )
    for stratum, selected in strata.items():
        baseline_counts = _counts(
            tuple(record.target_step_present for record in selected),
            tuple(record.baseline_detected for record in selected),
        )
        for rule in PRIMARY_RULES:
            predictions = tuple(_rule_prediction(record, rule) for record in selected)
            rows.append(
                _summary_row(
                    "abrupt_mean_step",
                    rule,
                    stratum,
                    selected,
                    predictions,
                    baseline_counts if rule != "baseline" else None,
                )
            )
    for target, truth_attr, detection_attr in (
        ("trend_component", "trend_component_present", "trend_flag"),
        ("damped_envelope", "damped_component_present", "envelope_flag"),
    ):
        selected = records
        truth = tuple(bool(getattr(record, truth_attr)) for record in selected)
        prediction = tuple(bool(getattr(record, detection_attr)) for record in selected)
        rows.append(_summary_from_arrays(target, target, "overall", truth, prediction))
    return tuple(rows)


def _rule_prediction(record: StepRecordResult, rule: str) -> bool:
    return bool(getattr(record, f"{rule}_detected"))


def _counts(
    truth: tuple[bool, ...], predicted: tuple[bool, ...]
) -> tuple[int, int, int, int]:
    return (
        sum(actual and result for actual, result in zip(truth, predicted, strict=True)),
        sum(
            not actual and result
            for actual, result in zip(truth, predicted, strict=True)
        ),
        sum(
            actual and not result
            for actual, result in zip(truth, predicted, strict=True)
        ),
        sum(
            not actual and not result
            for actual, result in zip(truth, predicted, strict=True)
        ),
    )


def _summary_row(
    target: str,
    rule: str,
    stratum: str,
    selected: tuple[StepRecordResult, ...],
    predicted: tuple[bool, ...],
    baseline_counts: tuple[int, int, int, int] | None,
) -> RuleSummary:
    truth = tuple(record.target_step_present for record in selected)
    return _summary_from_arrays(
        target, rule, stratum, truth, predicted, baseline_counts
    )


def _summary_from_arrays(
    target: str,
    rule: str,
    stratum: str,
    truth: tuple[bool, ...],
    predicted: tuple[bool, ...],
    baseline_counts: tuple[int, int, int, int] | None = None,
) -> RuleSummary:
    tp, fp, fn, tn = _counts(truth, predicted)
    positive_count, negative_count = tp + fn, fp + tn
    detection_rate = tp / positive_count if positive_count else None
    fpr = fp / negative_count if negative_count else None
    fnr = fn / positive_count if positive_count else None
    if baseline_counts is None:
        fpr_delta = fnr_delta = None
    else:
        base_tp, base_fp, base_fn, base_tn = baseline_counts
        base_fpr = base_fp / (base_fp + base_tn) if base_fp + base_tn else None
        base_fnr = base_fn / (base_tp + base_fn) if base_tp + base_fn else None
        fpr_delta = None if fpr is None or base_fpr is None else fpr - base_fpr
        fnr_delta = None if fnr is None or base_fnr is None else fnr - base_fnr
    return RuleSummary(
        target,
        rule,
        stratum,
        positive_count,
        negative_count,
        tp,
        fp,
        fn,
        tn,
        detection_rate,
        fpr,
        fnr,
        fpr_delta,
        fnr_delta,
    )


def summarize_predictions(
    target: str,
    rule: str,
    stratum: str,
    truth: tuple[bool, ...],
    predicted: tuple[bool, ...],
) -> RuleSummary:
    """Return confusion counts and rates with class-specific denominators."""
    if len(truth) != len(predicted) or not truth:
        raise ScientificValidationError(
            "truth and predicted labels must be nonempty and equally sized"
        )
    return _summary_from_arrays(target, rule, stratum, truth, predicted)


def _write_rows(path: Path, rows: tuple[Any, ...]) -> None:
    if not rows:
        raise ScientificValidationError("cannot write empty Phase 44 results")
    dictionaries = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)


def _write_plots(study: StepDetectionStudy, destination: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    overall = [
        row
        for row in study.summaries
        if row.target == "abrupt_mean_step" and row.stratum == "overall"
    ]
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axes[0].bar(
        [row.rule for row in overall], [row.false_positive_rate for row in overall]
    )
    axes[0].set_title("Evaluation false-positive rate")
    axes[0].set_ylim(0.0, 1.0)
    axes[1].bar(
        [row.rule for row in overall], [row.false_negative_rate for row in overall]
    )
    axes[1].set_title("Evaluation false-negative rate")
    axes[1].set_ylim(0.0, 1.0)
    for axis in axes:
        axis.tick_params(axis="x", labelrotation=35)
    figure.savefig(destination / "rule_error_rates.png", dpi=140)
    plt.close(figure)

    step_summaries = [
        row
        for row in study.summaries
        if row.target == "abrupt_mean_step"
        and row.stratum.startswith("step_magnitude=")
    ]
    step_magnitudes = sorted(
        {float(row.stratum.split("=", maxsplit=1)[1]) for row in step_summaries}
    )
    figure, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
    for rule in PRIMARY_RULES:
        axis.plot(
            step_magnitudes,
            [
                _detection_rate_for_step(step_summaries, rule, step)
                for step in step_magnitudes
            ],
            marker="o",
            label=rule,
        )
    axis.set_xlabel("Abrupt step magnitude")
    axis.set_ylabel("Evaluation detection rate")
    axis.set_ylim(0.0, 1.0)
    axis.legend()
    figure.savefig(destination / "weak_step_detection.png", dpi=140)
    plt.close(figure)


def _detection_rate_for_step(
    summaries: list[RuleSummary], rule: str, step: float
) -> float:
    summary = next(
        row
        for row in summaries
        if row.rule == rule and row.stratum == f"step_magnitude={step:g}"
    )
    if summary.detection_rate is None:
        raise ScientificValidationError("step-specific detection rate is unavailable")
    return summary.detection_rate


def main() -> None:
    """Run Phase 44 with its frozen default design and write dedicated outputs."""
    study = write_step_detection_study()
    print(
        f"Wrote {len(study.records)} calibration/evaluation records and "
        f"{len(study.summaries)} evaluation summaries"
    )


if __name__ == "__main__":
    main()
