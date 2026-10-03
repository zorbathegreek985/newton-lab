"""Deterministic synthetic time-series benchmark for Phase 42.

Signals and their known generating components are used only to assess simple
descriptive diagnostics and the existing causal exponential filter. They are
not market data and do not measure forecasting or economic performance.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from newton_lab.causal_diffusion_filter import exponential_filter
from newton_lab.exceptions import ScientificValidationError

FloatArray = NDArray[np.float64]
SignalKind = Literal[
    "stationary_noise",
    "trend_noise",
    "autocorrelated",
    "damped_oscillation",
    "variance_shift",
    "mean_shift",
]


@dataclass(frozen=True, slots=True)
class SyntheticSignal:
    """One observed sequence paired with known latent synthetic truth."""

    kind: SignalKind
    seed: int
    parameter_label: str
    parameter_value: float
    values: FloatArray
    truth: FloatArray


@dataclass(frozen=True, slots=True)
class BenchmarkRow:
    """Per-signal diagnostics; detector flags use fixed declared thresholds."""

    kind: str
    seed: int
    parameter_label: str
    parameter_value: float
    lag1_autocorrelation: float
    autocorrelation_flag: bool
    mean_shift_z: float
    mean_shift_flag: bool
    variance_ratio: float
    variance_shift_flag: bool
    raw_rmse_to_truth: float
    filtered_rmse_to_truth: float
    filter_improved_rmse: bool


def generate_signal(
    kind: SignalKind,
    *,
    sample_count: int = 512,
    seed: int = 4200,
    parameter_value: float = 1.0,
    trend_slope: float = 0.002,
    oscillation_amplitude: float = 1.0,
    damping_rate: float = 0.002,
    additional_mean_step: float = 0.0,
) -> SyntheticSignal:
    """Generate a reproducible process; ``parameter_value`` is case-specific.

    ``stationary_noise`` uses iid N(0, 1); ``trend_noise`` adds iid noise to a
    linear trend of slope 0.002/sample; ``autocorrelated`` uses a stationary
    AR(1), phi=0.8, after a 200-step burn-in. ``damped_oscillation`` is a
    deterministic decaying sinusoid plus iid noise. ``variance_shift`` changes
    Gaussian noise standard deviation from 1 to the given multiplier at the
    midpoint. ``mean_shift`` changes the latent mean from 0 to the given
    displacement at the midpoint. Case-specific parameter values set additive
    noise SD for trend/AR/oscillation, the post-change SD multiplier, or the
    mean displacement respectively. ``trend_slope``,
    ``oscillation_amplitude``, and ``damping_rate`` expose deterministic
    components for controlled studies. ``additional_mean_step`` composes a
    midpoint step with a trend or oscillation. The baseline has fixed unit
    noise.
    """
    if (
        isinstance(sample_count, bool)
        or not isinstance(sample_count, int)
        or sample_count < 64
    ):
        raise ScientificValidationError("sample_count must be an integer >= 64")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ScientificValidationError("seed must be a non-negative integer")
    if (
        isinstance(parameter_value, bool)
        or not np.isfinite(parameter_value)
        or parameter_value <= 0.0
    ):
        raise ScientificValidationError("parameter_value must be finite and positive")
    for value, name in (
        (trend_slope, "trend_slope"),
        (oscillation_amplitude, "oscillation_amplitude"),
        (damping_rate, "damping_rate"),
        (additional_mean_step, "additional_mean_step"),
    ):
        if isinstance(value, bool) or not np.isfinite(value) or value < 0.0:
            raise ScientificValidationError(f"{name} must be finite and non-negative")

    n = sample_count
    split = n // 2
    time = np.arange(n, dtype=np.float64)
    latent = np.zeros(n, dtype=np.float64)
    rng = np.random.default_rng(seed)
    if kind == "stationary_noise":
        label, value = "noise_sd", 1.0
        observed = rng.normal(0.0, value, n)
    elif kind == "trend_noise":
        label, value = "noise_sd", parameter_value
        latent = trend_slope * time
        observed = latent + rng.normal(0.0, value, n)
    elif kind == "autocorrelated":
        label, value = "observation_noise_sd", parameter_value
        innovations = rng.normal(0.0, 1.0, n + 200)
        ar = np.empty(n + 200, dtype=np.float64)
        ar[0] = innovations[0] / np.sqrt(1.0 - 0.8**2)
        for index in range(1, len(ar)):
            ar[index] = 0.8 * ar[index - 1] + innovations[index]
        latent = ar[200:]
        observed = latent + rng.normal(0.0, value, n)
    elif kind == "damped_oscillation":
        label, value = "noise_sd", parameter_value
        latent = (
            oscillation_amplitude
            * np.exp(-damping_rate * time)
            * np.sin(2.0 * np.pi * 0.10 * time)
        )
        observed = latent + rng.normal(0.0, value, n)
    elif kind == "variance_shift":
        label, value = "post_change_sd_multiplier", parameter_value
        scales = np.ones(n, dtype=np.float64)
        scales[split:] = value
        observed = rng.normal(0.0, scales, n)
    elif kind == "mean_shift":
        label, value = "mean_displacement", parameter_value
        latent[split:] = value
        observed = latent + rng.normal(0.0, 1.0, n)
    else:
        raise ScientificValidationError(f"unsupported synthetic signal kind: {kind}")
    if additional_mean_step > 0.0 and kind != "mean_shift":
        latent[split:] += additional_mean_step
        observed[split:] += additional_mean_step
    return SyntheticSignal(
        kind,
        seed,
        label,
        value,
        np.asarray(observed, dtype=np.float64),
        latent,
    )


def run_benchmark() -> tuple[BenchmarkRow, ...]:
    """Run the fixed 48-case seed/parameter grid and return ordered results."""
    configurations: tuple[tuple[SignalKind, tuple[float, ...]], ...] = (
        ("trend_noise", (0.25, 0.75, 1.5)),
        ("autocorrelated", (0.1, 0.5, 1.0)),
        ("damped_oscillation", (0.05, 0.2, 0.5)),
        ("variance_shift", (1.1, 1.4, 2.0)),
        ("mean_shift", (0.5, 1.0, 2.0)),
    )
    cases: list[tuple[SignalKind, float]] = [("stationary_noise", 1.0)]
    for kind, values in configurations:
        cases.extend((kind, value) for value in values)

    rows: list[BenchmarkRow] = []
    for case_index, (kind, parameter) in enumerate(cases):
        for replicate in range(3):
            seed = 42000 + case_index * 10 + replicate
            signal = generate_signal(kind, seed=seed, parameter_value=parameter)
            rows.append(_measure(signal))
    return tuple(rows)


def write_benchmark(
    output_directory: str | Path = "reports/phase_42_synthetic_benchmark",
) -> tuple[BenchmarkRow, ...]:
    """Write CSV, JSON metadata, and plots to a dedicated results directory."""
    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    rows = run_benchmark()
    csv_path = destination / "results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(rows[0])))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)
    metadata = {
        "experiment": "Phase 42 synthetic time-series benchmark",
        "sample_count": 512,
        "replicates_per_configuration": 3,
        "seeds": "42000 + configuration_index * 10 + replicate_index",
        "thresholds": {
            "absolute_lag1_autocorrelation": 0.25,
            "absolute_mean_shift_z": 3.0,
            "variance_ratio": 2.0,
        },
        "filter": "newton_lab.causal_diffusion_filter.exponential_filter(beta=0.8)",
        "results": "results.csv",
        "data": "synthetic only; no external or market data",
    }
    (destination / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    _write_plots(rows, destination)
    return rows


def main() -> None:
    """Run the fixed study and write results to the default report directory."""
    destination = Path("reports/phase_42_synthetic_benchmark")
    rows = write_benchmark(destination)
    print(f"Wrote {len(rows)} synthetic runs to {destination}")


def measure_signal(signal: SyntheticSignal) -> BenchmarkRow:
    """Compute Phase 42 statistics for one signal and its truth record."""
    values = signal.values
    midpoint = len(values) // 2
    centered = values - np.mean(values)
    denominator = float(np.dot(centered, centered))
    lag1 = (
        float(np.dot(centered[:-1], centered[1:]) / denominator)
        if denominator > 0.0
        else 0.0
    )
    first, second = values[:midpoint], values[midpoint:]
    se = float(
        np.sqrt(
            np.var(first, ddof=1) / len(first) + np.var(second, ddof=1) / len(second)
        )
    )
    mean_z = float((np.mean(second) - np.mean(first)) / se) if se > 0.0 else 0.0
    var_first = float(np.var(first, ddof=1))
    var_second = float(np.var(second, ddof=1))
    var_ratio = max(var_first, var_second) / min(var_first, var_second)
    filtered = exponential_filter(values, beta=0.8)
    raw_rmse = float(np.sqrt(np.mean(np.square(values - signal.truth))))
    filtered_rmse = float(np.sqrt(np.mean(np.square(filtered - signal.truth))))
    return BenchmarkRow(
        signal.kind,
        signal.seed,
        signal.parameter_label,
        signal.parameter_value,
        lag1,
        abs(lag1) >= 0.25,
        mean_z,
        abs(mean_z) >= 3.0,
        var_ratio,
        var_ratio >= 2.0,
        raw_rmse,
        filtered_rmse,
        filtered_rmse < raw_rmse,
    )


_measure = measure_signal


def _write_plots(rows: tuple[BenchmarkRow, ...], destination: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kinds = tuple(dict.fromkeys(row.kind for row in rows))
    measures = (
        ("autocorrelation_flag", "Lag-1 autocorrelation flag rate"),
        ("mean_shift_flag", "Mean-shift flag rate"),
        ("variance_shift_flag", "Variance-shift flag rate"),
    )
    figure, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for axis, (field, title) in zip(axes, measures, strict=True):
        rates = [
            np.mean([bool(getattr(row, field)) for row in rows if row.kind == kind])
            for kind in kinds
        ]
        axis.bar(range(len(kinds)), rates)
        axis.set_title(title)
        axis.set_ylim(0.0, 1.0)
        axis.set_ylabel("Fraction of deterministic runs flagged")
        axis.set_xticks(range(len(kinds)), kinds, rotation=70, ha="right")
    figure.savefig(destination / "detection_rates.png", dpi=140)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7, 5), constrained_layout=True)
    for kind in kinds:
        selected = [row for row in rows if row.kind == kind]
        axis.scatter(
            [row.lag1_autocorrelation for row in selected],
            [row.mean_shift_z for row in selected],
            label=kind,
            alpha=0.8,
        )
    axis.axvline(0.25, color="grey", linestyle="--", linewidth=1)
    axis.axhline(3.0, color="grey", linestyle="--", linewidth=1)
    axis.axhline(-3.0, color="grey", linestyle="--", linewidth=1)
    axis.set_xlabel("Sample lag-1 autocorrelation")
    axis.set_ylabel("Two-half mean difference / estimated standard error")
    axis.legend(fontsize=7)
    figure.savefig(destination / "diagnostic_overlap.png", dpi=140)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(8, 4), constrained_layout=True)
    positions = np.arange(len(kinds), dtype=np.float64)
    raw = [
        np.mean([r.raw_rmse_to_truth for r in rows if r.kind == kind]) for kind in kinds
    ]
    filtered = [
        np.mean([r.filtered_rmse_to_truth for r in rows if r.kind == kind])
        for kind in kinds
    ]
    axis.bar(positions - 0.18, raw, width=0.36, label="Observed input")
    axis.bar(positions + 0.18, filtered, width=0.36, label="Filtered (beta=0.8)")
    axis.set_xticks(positions, kinds, rotation=45, ha="right")
    axis.set_ylabel("RMSE to known synthetic latent component")
    axis.legend()
    figure.savefig(destination / "filter_error.png", dpi=140)
    plt.close(figure)


if __name__ == "__main__":
    main()
