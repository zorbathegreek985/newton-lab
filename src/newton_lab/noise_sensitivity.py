"""Phase 46 factorial noise sensitivity study for fixed step detectors."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from newton_lab.exceptions import ScientificValidationError
from newton_lab.local_step_detection import (
    local_mean_difference,
    trend_adjusted_local_difference,
)
from newton_lab.synthetic_timeseries_benchmark import (
    SyntheticSignal,
    generate_signal,
    measure_signal,
)

SAMPLE_COUNT = 512
MIDPOINT = SAMPLE_COUNT // 2
NOISE_LEVELS = (0.5, 1.0, 1.5)
STEP_MAGNITUDES = (0.1, 0.2, 0.5, 1.0)
TREND_SLOPE = 0.001
CALIBRATION_SEED_BASE = 46_000
EVALUATION_SEED_BASE = 200_000
LOCAL_RAW_W32_THRESHOLD = 0.20193448253267382
LOCAL_TREND_ADJUSTED_W32_THRESHOLD = 0.15418153275303856
LOCAL_TREND_ADJUSTED_W64_THRESHOLD = 0.19161592877807834
METHODS = (
    "phase44_baseline",
    "phase45_local_raw_w32",
    "phase45_local_trend_adjusted_w32",
    "phase45_local_trend_adjusted_w64",
)
Family = Literal["stationary_null", "pure_step", "gradual_trend", "step_plus_trend"]


@dataclass(frozen=True, slots=True)
class NoiseCondition:
    """One cell in the independently crossed noise and signal grid."""

    cell_id: str
    family: Family
    noise_sd: float
    step_magnitude: float
    trend_slope: float
    target_step_present: bool


@dataclass(frozen=True, slots=True)
class NoiseRecord:
    """One prediction by one frozen method on a generated signal."""

    split: Literal["calibration", "evaluation"]
    method: str
    cell_id: str
    family: str
    noise_sd: float
    step_magnitude: float
    trend_slope: float
    seed: int
    sample_count: int
    midpoint: int
    target_step_present: bool
    statistic: float
    threshold: float
    threshold_source: str
    detected: bool


@dataclass(frozen=True, slots=True)
class NoiseSummary:
    """Confusion counts, class denominators, and Wilson rate intervals."""

    method: str
    stratum: str
    positive_count: int
    negative_count: int
    true_positive: int
    false_positive: int
    false_negative: int
    true_negative: int
    detection_rate: float | None
    detection_rate_low_95: float | None
    detection_rate_high_95: float | None
    false_positive_rate: float | None
    false_positive_rate_low_95: float | None
    false_positive_rate_high_95: float | None
    false_negative_rate: float | None
    false_negative_rate_low_95: float | None
    false_negative_rate_high_95: float | None


@dataclass(frozen=True, slots=True)
class NoiseSensitivityStudy:
    """Complete Phase 46 conditions, paired records, and evaluation summaries."""

    conditions: tuple[NoiseCondition, ...]
    calibration_replicates: int
    evaluation_replicates: int
    calibration_seeds: tuple[int, ...]
    evaluation_seeds: tuple[int, ...]
    records: tuple[NoiseRecord, ...]
    summaries: tuple[NoiseSummary, ...]


def factorial_conditions() -> tuple[NoiseCondition, ...]:
    """Cross three noise SDs with null, trend, step, and step-plus-trend cells."""
    cells: tuple[tuple[str, Family, float, float, bool], ...] = (
        ("stationary_null", "stationary_null", 0.0, 0.0, False),
        ("trend_0.001", "gradual_trend", 0.0, TREND_SLOPE, False),
        *tuple(
            (f"pure_step_{step:g}", "pure_step", step, 0.0, True)
            for step in STEP_MAGNITUDES
        ),
        *tuple(
            (f"step_plus_trend_{step:g}", "step_plus_trend", step, TREND_SLOPE, True)
            for step in STEP_MAGNITUDES
        ),
    )
    return tuple(
        NoiseCondition(
            f"{cell_id}_noise_{noise:g}",
            family,
            noise,
            step,
            slope,
            positive,
        )
        for cell_id, family, step, slope, positive in cells
        for noise in NOISE_LEVELS
    )


def generate_noise_condition_signal(
    condition: NoiseCondition, seed: int
) -> SyntheticSignal:
    """Generate one fixed-signal cell with independently scaled Gaussian noise.

    The existing Phase 42 generator supplies the latent trend/step and unit-SD
    Gaussian draw. Only the noise residual is scaled; the latent component and
    all generator defaults remain unchanged.
    """
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ScientificValidationError("seed must be a non-negative integer")
    if not np.isfinite(condition.noise_sd) or condition.noise_sd <= 0.0:
        raise ScientificValidationError("noise_sd must be finite and positive")
    if condition.family == "stationary_null":
        source = generate_signal("stationary_noise", seed=seed)
    elif condition.family == "pure_step":
        source = generate_signal(
            "mean_shift",
            seed=seed,
            parameter_value=condition.step_magnitude,
        )
    elif condition.family in ("gradual_trend", "step_plus_trend"):
        source = generate_signal(
            "trend_noise",
            seed=seed,
            parameter_value=1.0,
            trend_slope=condition.trend_slope,
            additional_mean_step=condition.step_magnitude,
        )
    else:
        raise ValueError(f"unsupported signal family: {condition.family}")
    residual = source.values - source.truth
    values = source.truth + condition.noise_sd * residual
    return SyntheticSignal(
        source.kind,
        seed,
        "noise_sd",
        condition.noise_sd,
        np.asarray(values, dtype=np.float64),
        np.array(source.truth, dtype=np.float64, copy=True),
    )


def run_noise_sensitivity_study(
    *, calibration_replicates: int = 40, evaluation_replicates: int = 100
) -> NoiseSensitivityStudy:
    """Evaluate the four Phase 44/45 rules with all thresholds frozen."""
    _validate_replicates(calibration_replicates, "calibration_replicates")
    _validate_replicates(evaluation_replicates, "evaluation_replicates")
    conditions = factorial_conditions()
    cells = tuple(
        dict.fromkeys(
            condition.cell_id.rsplit("_noise_", 1)[0] for condition in conditions
        )
    )
    cell_indices = {cell: index for index, cell in enumerate(cells)}
    records: list[NoiseRecord] = []
    calibration_seeds: set[int] = set()
    evaluation_seeds: set[int] = set()
    for condition in conditions:
        cell = condition.cell_id.rsplit("_noise_", 1)[0]
        cell_index = cell_indices[cell]
        runs: tuple[
            tuple[Literal["calibration", "evaluation"], int, int, set[int]], ...
        ] = (
            (
                "calibration",
                calibration_replicates,
                CALIBRATION_SEED_BASE,
                calibration_seeds,
            ),
            (
                "evaluation",
                evaluation_replicates,
                EVALUATION_SEED_BASE,
                evaluation_seeds,
            ),
        )
        for split, replicates, base, seed_set in runs:
            for replicate in range(replicates):
                seed = base + cell_index * 1_000 + replicate
                seed_set.add(seed)
                signal = generate_noise_condition_signal(condition, seed)
                records.extend(_evaluate_record(split, condition, seed, signal))
    if calibration_seeds & evaluation_seeds:
        raise RuntimeError("calibration and evaluation seeds overlap")
    record_tuple = tuple(records)
    evaluation = tuple(row for row in record_tuple if row.split == "evaluation")
    return NoiseSensitivityStudy(
        conditions,
        calibration_replicates,
        evaluation_replicates,
        tuple(sorted(calibration_seeds)),
        tuple(sorted(evaluation_seeds)),
        record_tuple,
        _summarize(evaluation),
    )


def write_noise_sensitivity_study(
    output_directory: str | Path = "reports/phase_46_noise_sensitivity",
    *,
    calibration_replicates: int = 40,
    evaluation_replicates: int = 100,
) -> NoiseSensitivityStudy:
    """Write matched per-record results, summaries, metadata, and figures."""
    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    study = run_noise_sensitivity_study(
        calibration_replicates=calibration_replicates,
        evaluation_replicates=evaluation_replicates,
    )
    _write_csv(destination / "per_record_results.csv", study.records)
    _write_csv(destination / "summary.csv", study.summaries)
    metadata = {
        "experiment": "Phase 46 independent noise-level sensitivity",
        "data": "synthetic only",
        "sample_count": SAMPLE_COUNT,
        "midpoint": MIDPOINT,
        "noise_sd_levels": list(NOISE_LEVELS),
        "step_magnitudes": list(STEP_MAGNITUDES),
        "trend_slope": TREND_SLOPE,
        "calibration_replicates_per_condition": study.calibration_replicates,
        "evaluation_replicates_per_condition": study.evaluation_replicates,
        "calibration_seed_base": CALIBRATION_SEED_BASE,
        "evaluation_seed_base": EVALUATION_SEED_BASE,
        "calibration_seeds": list(study.calibration_seeds),
        "evaluation_seeds": list(study.evaluation_seeds),
        "seed_policy": (
            "paired across noise levels within each family and step magnitude; "
            "calibration and evaluation ranges are disjoint"
        ),
        "threshold_policy": (
            "no Phase 46 thresholds fitted; all thresholds are frozen from "
            "Phase 44/45 before evaluation"
        ),
        "methods": {
            "phase44_baseline": {
                "statistic": "Phase 44 two-half mean z statistic",
                "decision": "abs(z) >= 3",
                "threshold_source": "fixed Phase 42/44 threshold",
            },
            "phase45_local_raw_w32": {
                "statistic": "abs(local post-mean minus pre-mean), w=32",
                "threshold": LOCAL_RAW_W32_THRESHOLD,
                "threshold_source": "Phase 45 calibration; frozen",
            },
            "phase45_local_trend_adjusted_w32": {
                "statistic": "abs(local difference minus global OLS slope * 32)",
                "threshold": LOCAL_TREND_ADJUSTED_W32_THRESHOLD,
                "threshold_source": "Phase 45 calibration; frozen",
            },
            "phase45_local_trend_adjusted_w64": {
                "statistic": "abs(local difference minus global OLS slope * 64)",
                "threshold": LOCAL_TREND_ADJUSTED_W64_THRESHOLD,
                "threshold_source": "Phase 45 calibration; frozen",
            },
        },
        "conditions": [asdict(condition) for condition in study.conditions],
        "uncertainty": "Wilson 95% intervals for held-out class rates",
    }
    (destination / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    _write_plots(study, destination)
    return study


def _evaluate_record(
    split: Literal["calibration", "evaluation"],
    condition: NoiseCondition,
    seed: int,
    signal: SyntheticSignal,
) -> tuple[NoiseRecord, ...]:
    baseline = measure_signal(signal)
    statistics = {
        "phase44_baseline": abs(baseline.mean_shift_z),
        "phase45_local_raw_w32": abs(local_mean_difference(signal.values, 32)),
        "phase45_local_trend_adjusted_w32": abs(
            trend_adjusted_local_difference(signal.values, 32)
        ),
        "phase45_local_trend_adjusted_w64": abs(
            trend_adjusted_local_difference(signal.values, 64)
        ),
    }
    thresholds = {
        "phase44_baseline": 3.0,
        "phase45_local_raw_w32": LOCAL_RAW_W32_THRESHOLD,
        "phase45_local_trend_adjusted_w32": (LOCAL_TREND_ADJUSTED_W32_THRESHOLD),
        "phase45_local_trend_adjusted_w64": (LOCAL_TREND_ADJUSTED_W64_THRESHOLD),
    }
    sources = {
        "phase44_baseline": "Phase 42/44 fixed threshold",
        "phase45_local_raw_w32": "Phase 45 calibration threshold, frozen",
        "phase45_local_trend_adjusted_w32": "Phase 45 calibration threshold, frozen",
        "phase45_local_trend_adjusted_w64": "Phase 45 calibration threshold, frozen",
    }
    return tuple(
        NoiseRecord(
            split,
            method,
            condition.cell_id,
            condition.family,
            condition.noise_sd,
            condition.step_magnitude,
            condition.trend_slope,
            seed,
            SAMPLE_COUNT,
            MIDPOINT,
            condition.target_step_present,
            statistic,
            thresholds[method],
            sources[method],
            statistic >= thresholds[method],
        )
        for method, statistic in statistics.items()
    )


def _summarize(records: tuple[NoiseRecord, ...]) -> tuple[NoiseSummary, ...]:
    strata: dict[str, tuple[NoiseRecord, ...]] = {"overall": records}
    families = tuple(dict.fromkeys(row.family for row in records))
    noise_levels = tuple(dict.fromkeys(row.noise_sd for row in records))
    for family in families:
        strata[f"family:{family}"] = tuple(
            row for row in records if row.family == family
        )
    for noise in noise_levels:
        strata[f"noise_sd:{noise:g}"] = tuple(
            row for row in records if row.noise_sd == noise
        )
    for family in families:
        for noise in noise_levels:
            strata[f"family_noise:{family}|noise_sd={noise:g}"] = tuple(
                row for row in records if row.family == family and row.noise_sd == noise
            )
    for family in ("pure_step", "step_plus_trend"):
        for magnitude in STEP_MAGNITUDES:
            for noise in noise_levels:
                strata[
                    f"step_magnitude:{family}|magnitude={magnitude:g}|noise_sd={noise:g}"
                ] = tuple(
                    row
                    for row in records
                    if row.family == family
                    and row.step_magnitude == magnitude
                    and row.noise_sd == noise
                )
    summaries: list[NoiseSummary] = []
    for stratum, subset in strata.items():
        for method in METHODS:
            selected = tuple(row for row in subset if row.method == method)
            tp = sum(row.target_step_present and row.detected for row in selected)
            fp = sum(not row.target_step_present and row.detected for row in selected)
            fn = sum(row.target_step_present and not row.detected for row in selected)
            tn = sum(
                not row.target_step_present and not row.detected for row in selected
            )
            positives = tp + fn
            negatives = fp + tn
            tpr = tp / positives if positives else None
            fpr = fp / negatives if negatives else None
            fnr = fn / positives if positives else None
            tpr_interval = _wilson_interval(tp, positives)
            fpr_interval = _wilson_interval(fp, negatives)
            fnr_interval = _wilson_interval(fn, positives)
            summaries.append(
                NoiseSummary(
                    method,
                    stratum,
                    positives,
                    negatives,
                    tp,
                    fp,
                    fn,
                    tn,
                    tpr,
                    _bound(tpr_interval, 0),
                    _bound(tpr_interval, 1),
                    fpr,
                    _bound(fpr_interval, 0),
                    _bound(fpr_interval, 1),
                    fnr,
                    _bound(fnr_interval, 0),
                    _bound(fnr_interval, 1),
                )
            )
    return tuple(summaries)


def _wilson_interval(successes: int, total: int) -> tuple[float, float] | None:
    if total == 0:
        return None
    if successes == 0:
        rate = 0.0
        z = 1.959963984540054
        upper = z**2 / (total + z**2)
        return rate, upper
    if successes == total:
        z = 1.959963984540054
        lower = total / (total + z**2)
        return lower, 1.0
    z = 1.959963984540054
    rate = successes / total
    denominator = 1.0 + z**2 / total
    center = (rate + z**2 / (2.0 * total)) / denominator
    margin = (
        z * np.sqrt(rate * (1.0 - rate) / total + z**2 / (4.0 * total**2)) / denominator
    )
    return max(0.0, center - margin), min(1.0, center + margin)


def _bound(interval: tuple[float, float] | None, index: int) -> float | None:
    return None if interval is None else interval[index]


def _write_csv(
    path: Path, rows: tuple[NoiseRecord, ...] | tuple[NoiseSummary, ...]
) -> None:
    if not rows:
        raise ValueError("cannot write an empty Phase 46 table")
    values = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def _write_plots(study: NoiseSensitivityStudy, destination: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metric_specs = (
        (
            "false_positive_rate",
            "False-positive rate",
            ("stationary_null", "gradual_trend"),
        ),
        (
            "false_negative_rate",
            "False-negative rate",
            ("pure_step", "step_plus_trend"),
        ),
    )
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for row_index, (metric, title, selected_families) in enumerate(metric_specs):
        for column_index, family in enumerate(selected_families):
            axis = axes[row_index, column_index]
            for method in METHODS:
                selected = [
                    item
                    for item in study.summaries
                    if item.method == method
                    and item.stratum.startswith(f"family_noise:{family}|")
                ]
                selected.sort(key=lambda item: float(item.stratum.rsplit("=", 1)[1]))
                values = [getattr(item, metric) for item in selected]
                lows = [getattr(item, f"{metric}_low_95") for item in selected]
                highs = [getattr(item, f"{metric}_high_95") for item in selected]
                if any(value is None for value in values):
                    continue
                noise = [float(item.stratum.rsplit("=", 1)[1]) for item in selected]
                errors = np.asarray(
                    [
                        [
                            float(value) - float(low)
                            for value, low in zip(values, lows, strict=True)
                        ],
                        [
                            float(high) - float(value)
                            for value, high in zip(values, highs, strict=True)
                        ],
                    ]
                )
                errors = np.maximum(errors, 0.0)
                axis.errorbar(noise, values, yerr=errors, marker="o", label=method)
            axis.set_title(f"{title}: {family}")
            axis.set_xlabel("Noise SD")
            axis.set_ylabel(title)
            axis.set_ylim(0.0, 1.0)
            axis.grid(alpha=0.25)
            axis.legend(fontsize=7)
    figure.savefig(destination / "error_rates_by_noise.png", dpi=150)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    for axis, family in zip(axes, ("pure_step", "step_plus_trend"), strict=True):
        for noise_level in NOISE_LEVELS:
            for method in METHODS:
                selected = [
                    item
                    for item in study.summaries
                    if item.method == method
                    and item.stratum.startswith(f"step_magnitude:{family}|magnitude=")
                    and item.stratum.endswith(f"noise_sd={noise_level:g}")
                ]
                selected.sort(
                    key=lambda item: float(
                        item.stratum.split("|magnitude=", 1)[1].split("|", 1)[0]
                    )
                )
                magnitudes = [
                    float(item.stratum.split("|magnitude=", 1)[1].split("|", 1)[0])
                    for item in selected
                ]
                rates = [item.detection_rate for item in selected]
                axis.plot(
                    magnitudes,
                    rates,
                    marker="o",
                    label=f"{method}, noise={noise_level:g}",
                )
        axis.set_title(f"Detection by magnitude: {family}")
        axis.set_xlabel("Step magnitude")
        axis.set_ylabel("Detection rate")
        axis.set_ylim(0.0, 1.0)
        axis.grid(alpha=0.25)
        axis.legend(fontsize=6, ncol=2)
    figure.savefig(destination / "step_magnitude_sensitivity.png", dpi=150)
    plt.close(figure)


def _validate_replicates(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ScientificValidationError(f"{name} must be a positive integer")


def main() -> None:
    """Run the predeclared Phase 46 factorial study."""
    study = write_noise_sensitivity_study()
    print(
        f"Wrote {len(study.records)} method/record rows and "
        f"{len(study.summaries)} held-out summaries"
    )


if __name__ == "__main__":
    main()
