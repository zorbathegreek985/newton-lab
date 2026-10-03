"""Phase 50 controlled transfer test for recovery rates near a fold."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import NDArray
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError


@dataclass(frozen=True)
class FoldConditionEstimate:
    """One noisy observed return trajectory and its estimated recovery rate."""

    margin: float
    partition: str
    perturbation_fraction: float
    replicate: int
    observation_seed: int
    exact_linearized_rate: float
    measured_rate: float
    relative_rate_error: float
    log_fit_r_squared: float
    observation_count: int
    fitting_window: float
    observation_noise_sd: float
    max_relative_ode_error: float


@dataclass(frozen=True)
class MarginSummary:
    """Replicate distribution for one margin and perturbation condition."""

    margin: float
    partition: str
    perturbation_fraction: float
    predicted_local_rate: float
    measured_mean_rate: float
    measured_sd_rate: float
    measured_rate_ci95_low: float
    measured_rate_ci95_high: float
    mean_relative_error_from_local_rate: float
    maximum_relative_ode_error: float
    replicate_count: int


@dataclass(frozen=True)
class PredictorSummary:
    """Held-out score for one exponent fixed before the test evaluation."""

    predictor: str
    exponent: float
    calibration_prefactor: float
    calibration_margin_count: int
    held_out_margin_count: int
    held_out_mean_absolute_log_error: float
    held_out_rmse_relative_rate_error: float
    bootstrap_mean_absolute_log_error_ci95_low: float
    bootstrap_mean_absolute_log_error_ci95_high: float


def fold_equilibria(margin: float) -> tuple[float, float]:
    """Return stable and unstable equilibria of ``dx/dt=mu-x**2``.

    The boundary ``margin == 0`` is rejected because it is not an ordinary
    stable operating point: the two equilibria coalesce and local recovery
    vanishes there.
    """
    if isinstance(margin, bool) or not math.isfinite(margin) or margin <= 0.0:
        raise ScientificValidationError("margin must be finite and positive")
    root = math.sqrt(margin)
    return root, -root


def fold_local_recovery_rate(margin: float) -> float:
    """Return ``-lambda`` for the stable eigenvalue ``lambda=-2*sqrt(mu)``."""
    stable, _unstable = fold_equilibria(margin)
    return 2.0 * stable


def fold_exact_perturbation(
    margin: float, perturbation_fraction: float, times: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Evaluate the exact deviation from the stable branch after a positive pulse."""
    stable, _unstable = fold_equilibria(margin)
    if (
        isinstance(perturbation_fraction, bool)
        or not math.isfinite(perturbation_fraction)
        or perturbation_fraction <= 0.0
    ):
        raise ScientificValidationError("perturbation_fraction must be positive")
    sample_times = np.asarray(times, dtype=np.float64)
    if (
        sample_times.ndim != 1
        or sample_times.size < 2
        or not np.isfinite(sample_times).all()
        or np.any(sample_times < 0.0)
        or np.any(np.diff(sample_times) <= 0.0)
    ):
        raise ScientificValidationError(
            "times must be finite, non-negative, ordered samples"
        )
    rate = 2.0 * stable
    initial_deviation = perturbation_fraction * stable
    decay = np.exp(-rate * sample_times)
    return initial_deviation * decay / (1.0 + initial_deviation / rate * (1.0 - decay))


def estimate_fold_recovery_rate(
    times: NDArray[np.float64], observed_deviation: NDArray[np.float64]
) -> tuple[float, float]:
    """Estimate positive recovery rate from a log-linear observed deviation."""
    sample_times = np.asarray(times, dtype=np.float64)
    values = np.asarray(observed_deviation, dtype=np.float64)
    if (
        sample_times.ndim != 1
        or values.ndim != 1
        or sample_times.size != values.size
        or sample_times.size < 3
        or not np.isfinite(sample_times).all()
        or not np.isfinite(values).all()
        or np.any(np.diff(sample_times) <= 0.0)
        or np.any(values <= 0.0)
    ):
        raise ScientificValidationError(
            "rate estimation requires positive finite observations at ordered times"
        )
    logged = np.log(values / values[0])
    slope, intercept = np.polyfit(sample_times, logged, 1)
    fitted = slope * sample_times + intercept
    residual = float(np.sum((logged - fitted) ** 2))
    total = float(np.sum((logged - np.mean(logged)) ** 2))
    if slope >= 0.0 or total == 0.0:
        raise ScientificValidationError("observed deviation did not show decay")
    return float(-slope), float(1.0 - residual / total)


def _read_protocol(protocol_path: Path) -> dict[str, Any]:
    try:
        protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificValidationError(
            f"cannot read frozen protocol: {error}"
        ) from error
    if not isinstance(protocol, dict) or protocol.get("protocol_status") != (
        "frozen_before_final_evaluation"
    ):
        raise ScientificValidationError("a frozen Phase 50 protocol is required")
    return protocol


def _integrate_condition(
    margin: float,
    perturbation_fraction: float,
    rtol: float,
    atol: float,
    samples: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64], float, float]:
    stable, _unstable = fold_equilibria(margin)
    rate = fold_local_recovery_rate(margin)
    duration = math.log(5.0) / rate
    times = np.linspace(0.0, duration, samples, dtype=np.float64)
    initial = stable * (1.0 + perturbation_fraction)
    solution = solve_ivp(
        lambda _time, state: np.array([margin - state[0] ** 2]),
        (0.0, duration),
        np.array([initial], dtype=np.float64),
        method="DOP853",
        t_eval=times,
        rtol=rtol,
        atol=atol,
    )
    if not solution.success or solution.y.shape != (1, samples):
        raise RuntimeError(f"fold-control integration failed: {solution.message}")
    trajectory = solution.y[0].copy()
    exact_deviation = fold_exact_perturbation(margin, perturbation_fraction, times)
    exact_state = stable + exact_deviation
    relative_error = float(
        np.max(np.abs(trajectory - exact_state)) / np.max(np.abs(exact_state))
    )
    return times, trajectory, duration, relative_error


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write an empty Phase 50 result table")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _predictor_prefactor(
    margins: NDArray[np.float64], rates: NDArray[np.float64], exponent: float
) -> float:
    return float(np.exp(np.mean(np.log(rates) - exponent * np.log(margins))))


def _bootstrap_comparison(
    calibration_rates: NDArray[np.float64],
    test_rates: NDArray[np.float64],
    calibration_margins: NDArray[np.float64],
    test_margins: NDArray[np.float64],
    exponents: tuple[float, float],
    bootstrap_count: int,
    seed: int,
) -> tuple[dict[float, tuple[float, float]], tuple[float, float], float]:
    rng = np.random.default_rng(seed)
    sample_count = calibration_rates.shape[1]
    scores: dict[float, list[float]] = {exponent: [] for exponent in exponents}
    differences: list[float] = []
    for _ in range(bootstrap_count):
        indices = rng.integers(0, sample_count, size=sample_count)
        boot_calibration = calibration_rates[:, indices].mean(axis=1)
        boot_test = test_rates[:, indices]
        replicate_scores: list[float] = []
        for exponent in exponents:
            prefactor = _predictor_prefactor(
                calibration_margins, boot_calibration, exponent
            )
            predictions = prefactor * test_margins**exponent
            error = np.abs(np.log(predictions[:, None] / boot_test))
            score = float(np.mean(np.mean(error, axis=0)))
            scores[exponent].append(score)
            replicate_scores.append(score)
        differences.append(replicate_scores[0] - replicate_scores[1])
    intervals = {
        exponent: (
            float(np.percentile(values, 2.5)),
            float(np.percentile(values, 97.5)),
        )
        for exponent, values in scores.items()
    }
    difference_interval = (
        float(np.percentile(differences, 2.5)),
        float(np.percentile(differences, 97.5)),
    )
    return intervals, difference_interval, float(np.mean(differences))


def _plot_results(
    output: Path,
    margins: list[float],
    partitions: list[str],
    summaries: list[MarginSummary],
    predictors: list[PredictorSummary],
    exponents: tuple[float, float],
) -> None:
    calibration_margins = np.array(
        [
            margin
            for margin, partition in zip(margins, partitions, strict=True)
            if partition == "calibration"
        ]
    )
    test_margins = np.array(
        [
            margin
            for margin, partition in zip(margins, partitions, strict=True)
            if partition == "held_out"
        ]
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    for partition, marker in (("calibration", "o"), ("held_out", "s")):
        points = [
            row
            for row in summaries
            if row.perturbation_fraction == 0.01 and row.partition == partition
        ]
        ax.errorbar(
            [row.margin for row in points],
            [row.measured_mean_rate for row in points],
            yerr=[
                1.96 * row.measured_sd_rate / math.sqrt(row.replicate_count)
                for row in points
            ],
            fmt=marker,
            color="#444444",
            capsize=3,
            label=f"measured {partition} mean ± 95% normal CI",
        )
    colors = ("#2463a6", "#c8582b")
    predictors_by_exp = {row.exponent: row for row in predictors}
    for exponent, color in zip(exponents, colors, strict=True):
        result = predictors_by_exp[exponent]
        all_margins = np.array(margins)
        ax.plot(
            all_margins,
            result.calibration_prefactor * all_margins**exponent,
            color=color,
            linestyle="--",
            label=f"{result.predictor} (exponent {exponent:g})",
        )
    ax.axvline(max(calibration_margins), color="black", alpha=0.25, linestyle=":")
    ax.axvline(min(test_margins), color="black", alpha=0.25, linestyle=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Normalized margin μ to fold")
    ax.set_ylabel("Measured positive recovery rate")
    ax.set_title("Held-out recovery-rate predictions in a fold control benchmark")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        output / "held_out_recovery_rate_predictions.png",
        dpi=160,
        metadata={"Software": "Newton Lab Phase 50"},
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    # Per-replicate scores are written to a CSV; this figure plots the held-out
    # mean scores and paired bootstrap intervals without treating time samples
    # as independent replicates.
    for index, predictor in enumerate(predictors):
        ax.bar(index, predictor.held_out_mean_absolute_log_error, color=colors[index])
        ax.errorbar(
            index,
            predictor.held_out_mean_absolute_log_error,
            yerr=[
                [
                    predictor.held_out_mean_absolute_log_error
                    - predictor.bootstrap_mean_absolute_log_error_ci95_low
                ],
                [
                    predictor.bootstrap_mean_absolute_log_error_ci95_high
                    - predictor.held_out_mean_absolute_log_error
                ],
            ],
            fmt="none",
            color="black",
            capsize=4,
        )
    ax.set_xticks(range(len(predictors)), [row.predictor for row in predictors])
    ax.set_ylabel("Held-out mean absolute log rate error")
    ax.set_title("Paired held-out recovery-rate prediction error")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        output / "held_out_prediction_error.png",
        dpi=160,
        metadata={"Software": "Newton Lab Phase 50"},
    )
    plt.close(fig)


def run_fold_control_validation(output_directory: Path) -> dict[str, Any]:
    """Run the protocol-frozen fold-control transfer test and write artifacts."""
    output = Path(output_directory)
    protocol_path = output / "protocol.json"
    protocol = _read_protocol(protocol_path)
    calibration_margins = np.asarray(protocol["calibration_margins"], dtype=np.float64)
    test_margins = np.asarray(protocol["held_out_margins"], dtype=np.float64)
    margins = np.concatenate((calibration_margins, test_margins))
    partitions = ["calibration"] * len(calibration_margins) + ["held_out"] * len(
        test_margins
    )
    perturbations = tuple(
        float(value)
        for value in protocol["perturbation_fractions_of_stable_equilibrium"]
    )
    primary_perturbation = float(protocol["primary_perturbation_fraction"])
    replication_count = int(
        protocol["independent_measurement_noise_replicates_per_margin"]
    )
    noise_sd = float(
        protocol["measurement_noise"]["standard_deviation_normalized_state"]
    )
    seed = int(protocol["random_seed"])
    solver = protocol["ode_solver"]
    samples = int(solver["output_samples"])
    rtol = float(solver["relative_tolerance"])
    atol = float(solver["absolute_tolerance"])
    bootstrap_count = int(protocol["bootstrap_replicates"])
    bootstrap_seed = int(protocol["bootstrap_seed"])
    exponent_proposed = 0.5
    exponent_baseline = 1.0
    exponents = (exponent_proposed, exponent_baseline)

    records: list[FoldConditionEstimate] = []
    trajectory_errors: dict[tuple[float, float], float] = {}
    for margin_index, (margin, partition) in enumerate(
        zip(margins, partitions, strict=True)
    ):
        for perturbation in perturbations:
            rate = fold_local_recovery_rate(float(margin))
            times, trajectory, duration, integration_error = _integrate_condition(
                float(margin), perturbation, rtol, atol, samples
            )
            stable, _unstable = fold_equilibria(float(margin))
            true_deviation = trajectory - stable
            trajectory_errors[(float(margin), perturbation)] = integration_error
            for replicate in range(replication_count):
                observation_seed = seed + margin_index * replication_count + replicate
                noise_rng = np.random.default_rng(observation_seed)
                noise = noise_rng.normal(0.0, noise_sd, size=samples)
                observed = true_deviation + noise
                measured_rate, r_squared = estimate_fold_recovery_rate(times, observed)
                records.append(
                    FoldConditionEstimate(
                        margin=float(margin),
                        partition=partition,
                        perturbation_fraction=perturbation,
                        replicate=replicate,
                        observation_seed=observation_seed,
                        exact_linearized_rate=rate,
                        measured_rate=measured_rate,
                        relative_rate_error=abs(measured_rate - rate) / rate,
                        log_fit_r_squared=r_squared,
                        observation_count=samples,
                        fitting_window=duration,
                        observation_noise_sd=noise_sd,
                        max_relative_ode_error=integration_error,
                    )
                )

    condition_summaries: list[MarginSummary] = []
    estimate_matrix: dict[float, dict[float, NDArray[np.float64]]] = {}
    for margin, partition in zip(margins, partitions, strict=True):
        estimate_matrix[float(margin)] = {}
        for perturbation in perturbations:
            selected = [
                row
                for row in records
                if row.margin == margin and row.perturbation_fraction == perturbation
            ]
            estimates = np.asarray([row.measured_rate for row in selected])
            mean = float(np.mean(estimates))
            standard_deviation = float(np.std(estimates, ddof=1))
            standard_error = standard_deviation / math.sqrt(len(estimates))
            predicted_rate = fold_local_recovery_rate(float(margin))
            condition_summaries.append(
                MarginSummary(
                    margin=float(margin),
                    partition=partition,
                    perturbation_fraction=perturbation,
                    predicted_local_rate=predicted_rate,
                    measured_mean_rate=mean,
                    measured_sd_rate=standard_deviation,
                    measured_rate_ci95_low=mean - 1.96 * standard_error,
                    measured_rate_ci95_high=mean + 1.96 * standard_error,
                    mean_relative_error_from_local_rate=(
                        abs(mean - predicted_rate) / predicted_rate
                    ),
                    maximum_relative_ode_error=trajectory_errors[
                        (float(margin), perturbation)
                    ],
                    replicate_count=len(estimates),
                )
            )
            estimate_matrix[float(margin)][perturbation] = estimates

    calibration_small = np.array(
        [
            estimate_matrix[float(margin)][primary_perturbation]
            for margin in calibration_margins
        ]
    )
    test_small = np.array(
        [
            estimate_matrix[float(margin)][primary_perturbation]
            for margin in test_margins
        ]
    )
    calibration_means = calibration_small.mean(axis=1)
    predictor_summaries: list[PredictorSummary] = []
    for exponent, predictor_name in (
        (exponent_proposed, "fold_theory_exponent_0_5"),
        (exponent_baseline, "phase49_linear_exponent_1_baseline"),
    ):
        prefactor = _predictor_prefactor(
            calibration_margins, calibration_means, exponent
        )
        predictions = prefactor * test_margins**exponent
        log_errors = np.abs(np.log(predictions[:, None] / test_small))
        per_replicate_log_error = np.mean(log_errors, axis=0)
        relative_errors = predictions[:, None] / test_small - 1.0
        mean_absolute_log_error = float(np.mean(per_replicate_log_error))
        rmse_relative_error = float(np.sqrt(np.mean(relative_errors**2)))
        predictor_summaries.append(
            PredictorSummary(
                predictor=predictor_name,
                exponent=exponent,
                calibration_prefactor=prefactor,
                calibration_margin_count=len(calibration_margins),
                held_out_margin_count=len(test_margins),
                held_out_mean_absolute_log_error=mean_absolute_log_error,
                held_out_rmse_relative_rate_error=rmse_relative_error,
                bootstrap_mean_absolute_log_error_ci95_low=0.0,
                bootstrap_mean_absolute_log_error_ci95_high=0.0,
            )
        )

    bootstrap_intervals, difference_interval, bootstrap_mean_difference = (
        _bootstrap_comparison(
            calibration_small,
            test_small,
            calibration_margins,
            test_margins,
            exponents,
            bootstrap_count,
            bootstrap_seed,
        )
    )
    predictor_summaries = [
        PredictorSummary(
            **{
                **asdict(summary),
                "bootstrap_mean_absolute_log_error_ci95_low": bootstrap_intervals[
                    summary.exponent
                ][0],
                "bootstrap_mean_absolute_log_error_ci95_high": bootstrap_intervals[
                    summary.exponent
                ][1],
            }
        )
        for summary in predictor_summaries
    ]
    proposed_score = predictor_summaries[0].held_out_mean_absolute_log_error
    baseline_score = predictor_summaries[1].held_out_mean_absolute_log_error
    criteria = protocol["acceptance_criteria"]
    ratio = proposed_score / baseline_score
    prediction_supported = bool(
        ratio
        <= float(criteria["fold_method_test_error_ratio_to_exponent_one_baseline_max"])
        and difference_interval[1]
        < float(
            criteria[
                "paired_95_percent_bootstrap_ci_upper_bound_for_fold_minus_baseline_error"
            ]
        )
    )
    primary_means = [
        row
        for row in condition_summaries
        if row.perturbation_fraction == primary_perturbation
    ]
    primary_rate_accepted = all(
        row.mean_relative_error_from_local_rate
        <= float(criteria["primary_mean_rate_relative_error_vs_linearized_rate_max"])
        for row in primary_means
    )
    large_perturbation = max(perturbations)
    large_means = [
        row
        for row in condition_summaries
        if row.perturbation_fraction == large_perturbation
    ]
    large_rate_accepted = all(
        row.mean_relative_error_from_local_rate
        <= float(
            criteria[
                "large_perturbation_mean_rate_relative_error_vs_linearized_rate_max"
            ]
        )
        for row in large_means
    )
    max_ode_error = max(row.max_relative_ode_error for row in records)
    ode_accepted = max_ode_error <= float(
        criteria["maximum_relative_ode_trajectory_error_against_exact_fold_solution"]
    )

    _write_csv(output / "per_condition_results.csv", [asdict(row) for row in records])
    _write_csv(
        output / "condition_summary.csv",
        [asdict(row) for row in condition_summaries],
    )
    _write_csv(
        output / "predictor_comparison.csv",
        [asdict(row) for row in predictor_summaries],
    )
    comparison = {
        "proposed_predictor": predictor_summaries[0].predictor,
        "baseline_predictor": predictor_summaries[1].predictor,
        "proposed_to_baseline_error_ratio": ratio,
        "relative_error_reduction": 1.0 - ratio,
        "paired_bootstrap_mean_error_difference": bootstrap_mean_difference,
        "paired_bootstrap_difference_ci95_low": difference_interval[0],
        "paired_bootstrap_difference_ci95_high": difference_interval[1],
        "prediction_value_criterion_passed": prediction_supported,
        "fold_law_prediction_added_value": prediction_supported,
    }
    _write_csv(output / "paired_comparison.csv", [comparison])
    _plot_results(
        output,
        [float(value) for value in margins],
        partitions,
        condition_summaries,
        predictor_summaries,
        exponents,
    )
    protocol_hash = hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    metadata = {
        "protocol_sha256": protocol_hash,
        "protocol_status": protocol["protocol_status"],
        "data_provenance": protocol["provenance"],
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "record_count": len(records),
        "condition_count": len(condition_summaries),
        "independent_noise_replicates": replication_count,
        "maximum_relative_ode_error_observed": max_ode_error,
        "ode_accuracy_criterion_passed": ode_accepted,
        "primary_rate_measurement_criterion_passed": primary_rate_accepted,
        "large_perturbation_sensitivity_criterion_passed": large_rate_accepted,
        "predictive_value_criterion_passed": prediction_supported,
        "proposed_to_baseline_error_ratio": ratio,
        "relative_error_reduction": 1.0 - ratio,
        "paired_bootstrap_difference_ci95": list(difference_interval),
        "bootstrap_mean_difference": bootstrap_mean_difference,
        "predictors": [asdict(row) for row in predictor_summaries],
        "limitation": (
            "Canonical saddle-node benchmark validates behavior only for its "
            "declared synthetic control model; no specific engineering plant "
            "is validated."
        ),
    }
    (output / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return metadata
