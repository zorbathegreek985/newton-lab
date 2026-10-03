"""Synthetic identifiability benchmark for canonical bifurcation recovery rates."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import NDArray
from scipy.optimize import (  # type: ignore[import-untyped]
    least_squares,
    minimize_scalar,
)

from newton_lab.exceptions import ScientificValidationError

MODEL_NAMES = ("fold", "transcritical", "supercritical_pitchfork")
INDETERMINATE = "indeterminate"
TRUE_MARGINS = (0.04, 0.09, 0.16, 0.25, 0.36)
CONTROL_BOUNDARY = 1.0
CONTROL_SETTINGS = tuple(CONTROL_BOUNDARY + margin for margin in TRUE_MARGINS)
NOISE_LEVELS = (0.0005, 0.002)
WINDOWS_S = (3.0, 8.0)
SAMPLE_INTERVAL_S = 0.1
PERTURBATION_FRACTION = 0.2
CALIBRATION_SEEDS = tuple(range(41000, 41020))
EVALUATION_SEEDS = tuple(range(51000, 51040))
MODEL_INFORMED_MIN_BIC_GAP = 6.0
AGNOSTIC_MIN_POSTERIOR = 0.60
AGNOSTIC_MIN_LOG_MARGIN = math.log(2.0)
BOOTSTRAP_REPLICATES = 500

_SIGNATURES: dict[str, tuple[float, float]] = {
    "fold": (0.5, 0.5),
    "transcritical": (1.0, 1.0),
    "supercritical_pitchfork": (0.5, 1.0),
}


@dataclass(frozen=True)
class ExponentialEstimate:
    """Free exponential fit; the offset estimates equilibrium from samples."""

    equilibrium: float
    equilibrium_standard_error: float
    recovery_rate: float
    recovery_rate_standard_error: float
    amplitude: float
    r_squared: float
    fit_rmse: float
    converged: bool


@dataclass(frozen=True)
class ModelTrajectoryFit:
    """One candidate normal-form fit to a trace using an estimated margin."""

    model: str
    estimated_margin: float
    estimated_equilibrium: float
    predicted_recovery_rate: float
    normalized_rss: float
    bic: float
    converged: bool


@dataclass(frozen=True)
class ScalingEstimate:
    """Joint model-agnostic power fit with an estimated control boundary."""

    equilibrium_exponent: float
    rate_exponent: float
    estimated_boundary: float
    equilibrium_log_rmse: float
    rate_log_rmse: float
    converged: bool


def stable_equilibrium(model: str, margin: float) -> float:
    """Return the selected stable equilibrium for positive ``margin``.

    The tested branches are fold ``+sqrt(mu)``, transcritical ``mu``, and
    positive supercritical-pitchfork ``+sqrt(mu)``. This function deliberately
    excludes the bifurcation point and all other branches/sides.
    """
    _validate_model_margin(model, margin)
    if model == "transcritical":
        return margin
    return math.sqrt(margin)


def local_recovery_rate(model: str, margin: float) -> float:
    """Return the positive negative-Jacobian eigenvalue on the tested branch."""
    equilibrium = stable_equilibrium(model, margin)
    if model == "fold":
        return 2.0 * equilibrium
    if model == "transcritical":
        return margin
    return 2.0 * margin


def exact_trajectory(
    model: str,
    margin: float,
    initial_state: float,
    times_s: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Evaluate the exact positive-branch trajectory for one canonical model."""
    _validate_model_margin(model, margin)
    times = np.asarray(times_s, dtype=np.float64)
    if (
        times.ndim != 1
        or times.size < 2
        or not np.isfinite(times).all()
        or times[0] < 0.0
        or np.any(np.diff(times) <= 0.0)
        or not math.isfinite(initial_state)
        or initial_state <= 0.0
    ):
        raise ScientificValidationError(
            "trajectory inputs require positive finite state and ordered times"
        )

    equilibrium = stable_equilibrium(model, margin)
    if model == "fold":
        root = equilibrium
        ratio = (initial_state - root) / (initial_state + root)
        decay = np.exp(-2.0 * root * times)
        denominator = 1.0 - ratio * decay
        if np.any(np.abs(denominator) < 1e-14):
            raise ScientificValidationError("fold trajectory is singular")
        return root * (1.0 + ratio * decay) / denominator
    if model == "transcritical":
        denominator = 1.0 + (margin / initial_state - 1.0) * np.exp(-margin * times)
        if np.any(denominator <= 0.0):
            raise ScientificValidationError("transcritical trajectory is invalid")
        return margin / denominator

    denominator = 1.0 + (margin / initial_state**2 - 1.0) * np.exp(
        -2.0 * margin * times
    )
    if np.any(denominator <= 0.0):
        raise ScientificValidationError("pitchfork trajectory is invalid")
    return np.sqrt(margin / denominator)


def generate_observation(
    model: str,
    margin: float,
    noise_sd: float,
    duration_s: float,
    seed: int,
    *,
    sample_interval_s: float = SAMPLE_INTERVAL_S,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Generate one seeded noisy observation; truth is not returned to estimators."""
    if not math.isfinite(noise_sd) or noise_sd <= 0.0:
        raise ScientificValidationError("noise_sd must be finite and positive")
    if not math.isfinite(duration_s) or duration_s <= 0.0:
        raise ScientificValidationError("duration_s must be finite and positive")
    if not math.isfinite(sample_interval_s) or sample_interval_s <= 0.0:
        raise ScientificValidationError("sample interval must be finite and positive")
    equilibrium = stable_equilibrium(model, margin)
    count = int(round(duration_s / sample_interval_s)) + 1
    times = np.linspace(0.0, duration_s, count, dtype=np.float64)
    initial_state = equilibrium * (1.0 + PERTURBATION_FRACTION)
    clean = exact_trajectory(model, margin, initial_state, times)
    rng = np.random.default_rng(seed)
    return times, clean + rng.normal(0.0, noise_sd, size=count)


def estimate_exponential(
    times_s: NDArray[np.float64],
    observed: NDArray[np.float64],
    noise_sd: float,
) -> ExponentialEstimate:
    """Fit ``b + a exp(-rho t)`` with free ``a``, ``b`` and positive ``rho``."""
    times = np.asarray(times_s, dtype=np.float64)
    values = np.asarray(observed, dtype=np.float64)
    if (
        times.ndim != 1
        or values.ndim != 1
        or times.size != values.size
        or times.size < 5
        or not np.isfinite(times).all()
        or not np.isfinite(values).all()
        or not math.isfinite(noise_sd)
        or noise_sd <= 0.0
        or np.any(np.diff(times) <= 0.0)
    ):
        raise ScientificValidationError("invalid data for exponential estimation")

    initial_amplitude = max(float(values[0] - values[-1]), noise_sd)
    initial = np.array(
        [math.log(initial_amplitude), math.log(0.2), float(values[-1])],
        dtype=np.float64,
    )

    def residuals(parameters: NDArray[np.float64]) -> NDArray[np.float64]:
        amplitude = math.exp(float(parameters[0]))
        rate = math.exp(float(parameters[1]))
        return np.asarray(
            amplitude * np.exp(-rate * times) + parameters[2] - values,
            dtype=np.float64,
        )

    def jacobian(parameters: NDArray[np.float64]) -> NDArray[np.float64]:
        amplitude = math.exp(float(parameters[0]))
        rate = math.exp(float(parameters[1]))
        decay = np.exp(-rate * times)
        return np.column_stack(
            (amplitude * decay, -amplitude * rate * times * decay, np.ones_like(times))
        )

    lower = np.array([math.log(1e-9), math.log(1e-5), -10.0])
    upper = np.array([math.log(10.0), math.log(10.0), 10.0])
    fit = least_squares(
        residuals,
        initial,
        jac=jacobian,
        bounds=(lower, upper),
        max_nfev=2000,
    )
    amplitude = math.exp(float(fit.x[0]))
    rate = math.exp(float(fit.x[1]))
    equilibrium = float(fit.x[2])
    residual = residuals(fit.x)
    rss = float(np.dot(residual, residual))
    dof = max(1, times.size - 3)
    try:
        covariance = np.linalg.inv(fit.jac.T @ fit.jac) * (rss / dof)
        rate_se = rate * math.sqrt(max(0.0, float(covariance[1, 1])))
        equilibrium_se = math.sqrt(max(0.0, float(covariance[2, 2])))
    except np.linalg.LinAlgError:
        rate_se = math.inf
        equilibrium_se = math.inf
    centered = values - np.mean(values)
    total = float(np.dot(centered, centered))
    r_squared = 1.0 - rss / total if total > 0.0 else 0.0
    return ExponentialEstimate(
        equilibrium=equilibrium,
        equilibrium_standard_error=equilibrium_se,
        recovery_rate=rate,
        recovery_rate_standard_error=rate_se,
        amplitude=amplitude,
        r_squared=r_squared,
        fit_rmse=math.sqrt(rss / times.size),
        converged=bool(fit.success and math.isfinite(rate) and rate > 0.0),
    )


def fit_canonical_trajectory(
    model: str,
    times_s: NDArray[np.float64],
    observed: NDArray[np.float64],
    noise_sd: float,
) -> ModelTrajectoryFit:
    """Fit unknown positive margin for a named model without truth inputs."""
    if model not in MODEL_NAMES:
        raise ScientificValidationError(f"unsupported model: {model}")
    times = np.asarray(times_s, dtype=np.float64)
    values = np.asarray(observed, dtype=np.float64)
    if (
        times.ndim != 1
        or values.ndim != 1
        or times.size != values.size
        or times.size < 5
        or not np.isfinite(times).all()
        or not np.isfinite(values).all()
        or not math.isfinite(noise_sd)
        or noise_sd <= 0.0
        or np.any(np.diff(times) <= 0.0)
        or values[0] <= 0.0
    ):
        raise ScientificValidationError("invalid data for canonical model fit")
    initial_state = float(values[0])
    sigma = float(noise_sd)

    def objective(log_margin: float) -> float:
        margin = math.exp(log_margin)
        try:
            prediction = exact_trajectory(model, margin, initial_state, times)
        except ScientificValidationError:
            return math.inf
        residual = (prediction - values) / sigma
        return float(np.dot(residual, residual))

    result = minimize_scalar(
        objective,
        method="bounded",
        bounds=(math.log(0.002), math.log(1.5)),
        options={"xatol": 1e-8, "maxiter": 300},
    )
    margin = math.exp(float(result.x))
    normalized_rss = objective(float(result.x))
    n = times.size
    bic = normalized_rss + math.log(n)
    equilibrium = stable_equilibrium(model, margin)
    rate = local_recovery_rate(model, margin)
    return ModelTrajectoryFit(
        model=model,
        estimated_margin=margin,
        estimated_equilibrium=equilibrium,
        predicted_recovery_rate=rate,
        normalized_rss=normalized_rss,
        bic=bic,
        converged=bool(result.success and math.isfinite(normalized_rss)),
    )


def estimate_scaling(
    control_settings: NDArray[np.float64],
    estimated_equilibria: NDArray[np.float64],
    estimated_rates: NDArray[np.float64],
) -> ScalingEstimate:
    """Estimate two power exponents and an unknown shared control boundary."""
    settings = np.asarray(control_settings, dtype=np.float64)
    equilibria = np.asarray(estimated_equilibria, dtype=np.float64)
    rates = np.asarray(estimated_rates, dtype=np.float64)
    if (
        settings.ndim != 1
        or equilibria.shape != settings.shape
        or rates.shape != settings.shape
        or settings.size < 4
        or not np.isfinite(settings).all()
        or not np.isfinite(equilibria).all()
        or not np.isfinite(rates).all()
        or np.any(equilibria <= 0.0)
        or np.any(rates <= 0.0)
        or np.any(np.diff(settings) <= 0.0)
    ):
        raise ScientificValidationError("scaling fit requires positive ordered data")

    log_x = np.log(equilibria)
    log_r = np.log(rates)
    minimum = float(np.min(settings))

    def residuals(parameters: NDArray[np.float64]) -> NDArray[np.float64]:
        log_cx, q, log_cr, p, log_gap = parameters
        boundary = minimum - math.exp(float(log_gap))
        log_distance = np.log(settings - boundary)
        return np.concatenate(
            (log_cx + q * log_distance - log_x, log_cr + p * log_distance - log_r)
        )

    gap0 = max(0.05, float(settings[-1] - settings[0]) / 2.0)
    initial = np.array(
        [float(np.mean(log_x)), 0.7, float(np.mean(log_r)), 0.7, math.log(gap0)]
    )
    fit = least_squares(
        residuals,
        initial,
        bounds=(
            np.array([-20.0, 0.05, -20.0, 0.05, math.log(1e-5)]),
            np.array([20.0, 2.0, 20.0, 2.0, math.log(5.0)]),
        ),
        max_nfev=3000,
    )
    residual = residuals(fit.x)
    n = settings.size
    eq_rmse = math.sqrt(float(np.mean(residual[:n] ** 2)))
    rate_rmse = math.sqrt(float(np.mean(residual[n:] ** 2)))
    return ScalingEstimate(
        equilibrium_exponent=float(fit.x[1]),
        rate_exponent=float(fit.x[3]),
        estimated_boundary=minimum - math.exp(float(fit.x[4])),
        equilibrium_log_rmse=eq_rmse,
        rate_log_rmse=rate_rmse,
        converged=bool(fit.success and np.isfinite(fit.x).all()),
    )


def run_identifiability_study(output_dir: Path) -> dict[str, Any]:
    """Run the frozen seeded benchmark and write reproducible CSV/JSON outputs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    protocol_path = output_dir / "protocol.json"
    try:
        protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificValidationError(
            f"cannot read frozen protocol: {error}"
        ) from error
    if protocol.get("protocol_status") != "frozen_before_implementation_and_evaluation":
        raise ScientificValidationError("a frozen Phase 53 protocol is required")
    protocol_hash = hashlib.sha256(protocol_path.read_bytes()).hexdigest()

    records: list[dict[str, Any]] = []
    noise_index = {value: index for index, value in enumerate(NOISE_LEVELS)}
    window_index = {value: index for index, value in enumerate(WINDOWS_S)}
    model_index = {value: index for index, value in enumerate(MODEL_NAMES)}
    for partition, seeds in (
        ("calibration", CALIBRATION_SEEDS),
        ("evaluation", EVALUATION_SEEDS),
    ):
        for seed in seeds:
            for model in MODEL_NAMES:
                for noise_sd in NOISE_LEVELS:
                    for duration_s in WINDOWS_S:
                        for margin, setting in zip(
                            TRUE_MARGINS, CONTROL_SETTINGS, strict=True
                        ):
                            sample_seed = np.random.SeedSequence(
                                [
                                    seed,
                                    model_index[model],
                                    noise_index[noise_sd],
                                    window_index[duration_s],
                                    TRUE_MARGINS.index(margin),
                                ]
                            )
                            # SeedSequence emits an isolated child seed without
                            # exposing the true margin to either estimator.
                            trace_seed = int(sample_seed.generate_state(1)[0])
                            times, observed = generate_observation(
                                model, margin, noise_sd, duration_s, trace_seed
                            )
                            agnostic = estimate_exponential(times, observed, noise_sd)
                            candidate_fits = {
                                candidate: fit_canonical_trajectory(
                                    candidate, times, observed, noise_sd
                                )
                                for candidate in MODEL_NAMES
                            }
                            ordered = sorted(
                                candidate_fits.values(),
                                key=lambda candidate: candidate.bic,
                            )
                            gap = ordered[1].bic - ordered[0].bic
                            informed_prediction = (
                                ordered[0].model
                                if gap >= MODEL_INFORMED_MIN_BIC_GAP
                                else INDETERMINATE
                            )
                            true_fit = candidate_fits[model]
                            records.append(
                                {
                                    "partition": partition,
                                    "seed": seed,
                                    "model": model,
                                    "noise_sd": noise_sd,
                                    "window_s": duration_s,
                                    "control_setting": setting,
                                    "true_margin": margin,
                                    "true_equilibrium": stable_equilibrium(
                                        model, margin
                                    ),
                                    "true_recovery_rate": local_recovery_rate(
                                        model, margin
                                    ),
                                    "estimated_equilibrium": agnostic.equilibrium,
                                    "equilibrium_se": (
                                        agnostic.equilibrium_standard_error
                                    ),
                                    "estimated_rate": agnostic.recovery_rate,
                                    "rate_se": agnostic.recovery_rate_standard_error,
                                    "exponential_r_squared": agnostic.r_squared,
                                    "exponential_fit_rmse": agnostic.fit_rmse,
                                    "exponential_converged": agnostic.converged,
                                    "true_model_margin_fit": true_fit.estimated_margin,
                                    "true_model_equilibrium_fit": (
                                        true_fit.estimated_equilibrium
                                    ),
                                    "true_model_bic": true_fit.bic,
                                    "model_informed_prediction": informed_prediction,
                                    "model_informed_bic_gap": gap,
                                    "model_agnostic_prediction": INDETERMINATE,
                                    **{
                                        f"{name}_bic": candidate_fits[name].bic
                                        for name in MODEL_NAMES
                                    },
                                    **{
                                        f"{name}_margin_fit": candidate_fits[
                                            name
                                        ].estimated_margin
                                        for name in MODEL_NAMES
                                    },
                                }
                            )

    _classify_agnostic(records)
    _write_csv(output_dir / "per_trajectory_estimates.csv", records)
    confusion = _confusion_table(records)
    _write_csv(output_dir / "classification_confusion.csv", confusion)
    metrics = _classification_metrics(records)
    _write_csv(output_dir / "classification_metrics.csv", metrics)
    scaling_rows, scaling_summary = _scaling_results(records)
    _write_csv(output_dir / "scaling_exponent_estimates.csv", scaling_rows)
    _write_csv(output_dir / "scaling_exponent_summary.csv", scaling_summary)
    parameter_summary = _parameter_recovery_summary(records)
    _write_csv(output_dir / "parameter_recovery_summary.csv", parameter_summary)

    evaluation = [record for record in records if record["partition"] == "evaluation"]
    primary = [
        row
        for row in metrics
        if row["noise_sd"] == min(NOISE_LEVELS)
        and row["window_s"] == max(WINDOWS_S)
        and row["distance"] == "all"
    ]
    criterion_results = _evaluate_frozen_criteria(
        evaluation, primary, scaling_summary, parameter_summary
    )
    metadata: dict[str, Any] = {
        "protocol_status": protocol["protocol_status"],
        "protocol_sha256": protocol_hash,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "record_count": len(records),
        "calibration_record_count": sum(
            r["partition"] == "calibration" for r in records
        ),
        "evaluation_record_count": len(evaluation),
        "calibration_seed_count": len(CALIBRATION_SEEDS),
        "evaluation_seed_count": len(EVALUATION_SEEDS),
        "noise_levels": list(NOISE_LEVELS),
        "windows_s": list(WINDOWS_S),
        "margins_for_generation_only": list(TRUE_MARGINS),
        "equations": {
            "fold": "dx/dt=mu-x^2; stable x=+sqrt(mu), rho=2sqrt(mu), mu>0",
            "transcritical": "dx/dt=mu*x-x^2; stable x=mu, rho=mu, mu>0",
            "supercritical_pitchfork": (
                "dx/dt=mu*x-x^3; stable x=+/-sqrt(mu), rho=2mu, mu>0"
            ),
        },
        "frozen_criteria": criterion_results,
        "interpretation_limit": (
            "All target trajectories are generated by the canonical equations under "
            "study. Results evaluate identifiability only for the frozen synthetic "
            "design."
        ),
    }
    metadata_path = output_dir / "metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = _hash_manifest(output_dir)
    (output_dir / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata


def _validate_model_margin(model: str, margin: float) -> None:
    if model not in MODEL_NAMES:
        raise ScientificValidationError(f"unsupported model: {model}")
    if isinstance(margin, bool) or not math.isfinite(margin) or margin <= 0.0:
        raise ScientificValidationError("margin must be finite and positive")


def _classify_agnostic(records: list[dict[str, Any]]) -> None:
    training: dict[
        tuple[str, float, float, float], tuple[NDArray[Any], NDArray[Any]]
    ] = {}
    for model in MODEL_NAMES:
        for noise_sd in NOISE_LEVELS:
            for duration_s in WINDOWS_S:
                for setting in CONTROL_SETTINGS:
                    samples = [
                        row
                        for row in records
                        if row["partition"] == "calibration"
                        and row["model"] == model
                        and row["noise_sd"] == noise_sd
                        and row["window_s"] == duration_s
                        and row["control_setting"] == setting
                        and row["estimated_equilibrium"] > 0.0
                        and row["estimated_rate"] > 0.0
                    ]
                    if not samples:
                        continue
                    features = np.array(
                        [
                            [
                                math.log(row["estimated_equilibrium"]),
                                math.log(row["estimated_rate"]),
                            ]
                            for row in samples
                        ],
                        dtype=np.float64,
                    )
                    mean = np.mean(features, axis=0)
                    variance = np.maximum(np.var(features, axis=0, ddof=1), 0.02**2)
                    training[(model, noise_sd, duration_s, setting)] = (mean, variance)

    for row in records:
        if row["partition"] != "evaluation":
            continue
        if row["estimated_equilibrium"] <= 0.0 or row["estimated_rate"] <= 0.0:
            continue
        feature = np.log(
            np.array([row["estimated_equilibrium"], row["estimated_rate"]])
        )
        log_likelihoods: list[tuple[str, float]] = []
        for model in MODEL_NAMES:
            stats = training.get(
                (model, row["noise_sd"], row["window_s"], row["control_setting"])
            )
            if stats is None:
                continue
            mean, variance = stats
            score = -0.5 * float(
                np.sum((feature - mean) ** 2 / variance + np.log(variance))
            )
            log_likelihoods.append((model, score))
        if len(log_likelihoods) != len(MODEL_NAMES):
            continue
        log_likelihoods.sort(key=lambda item: item[1], reverse=True)
        normalized = np.array([score for _, score in log_likelihoods], dtype=np.float64)
        posterior = np.exp(normalized - float(np.max(normalized)))
        posterior /= np.sum(posterior)
        log_margin = log_likelihoods[0][1] - log_likelihoods[1][1]
        prediction = (
            log_likelihoods[0][0]
            if posterior[0] >= AGNOSTIC_MIN_POSTERIOR
            and log_margin >= AGNOSTIC_MIN_LOG_MARGIN
            else INDETERMINATE
        )
        row["model_agnostic_prediction"] = prediction


def _confusion_table(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    evaluation = [r for r in records if r["partition"] == "evaluation"]
    for method in ("model_informed", "model_agnostic"):
        prediction_key = f"{method}_prediction"
        for noise_sd in NOISE_LEVELS:
            for duration_s in WINDOWS_S:
                for margin in TRUE_MARGINS:
                    setting = CONTROL_BOUNDARY + margin
                    subset = [
                        r
                        for r in evaluation
                        if r["noise_sd"] == noise_sd
                        and r["window_s"] == duration_s
                        and r["control_setting"] == setting
                    ]
                    for actual in MODEL_NAMES:
                        for predicted in (*MODEL_NAMES, INDETERMINATE):
                            count = sum(
                                r["model"] == actual and r[prediction_key] == predicted
                                for r in subset
                            )
                            rows.append(
                                {
                                    "method": method,
                                    "noise_sd": noise_sd,
                                    "window_s": duration_s,
                                    "margin_for_scoring_only": margin,
                                    "true_model": actual,
                                    "predicted_model": predicted,
                                    "count": count,
                                }
                            )
    return rows


def _classification_metrics(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    evaluation = [r for r in records if r["partition"] == "evaluation"]
    for method in ("model_informed", "model_agnostic"):
        prediction_key = f"{method}_prediction"
        for noise_sd in NOISE_LEVELS:
            for duration_s in WINDOWS_S:
                for distance, setting in zip(
                    TRUE_MARGINS, CONTROL_SETTINGS, strict=True
                ):
                    subset = [
                        r
                        for r in evaluation
                        if r["noise_sd"] == noise_sd
                        and r["window_s"] == duration_s
                        and r["control_setting"] == setting
                    ]
                    rows.append(
                        _metric_row(
                            method,
                            noise_sd,
                            duration_s,
                            distance,
                            subset,
                            prediction_key,
                        )
                    )
                rows.append(
                    _metric_row(
                        method,
                        noise_sd,
                        duration_s,
                        "all",
                        [
                            r
                            for r in evaluation
                            if r["noise_sd"] == noise_sd and r["window_s"] == duration_s
                        ],
                        prediction_key,
                    )
                )
    return rows


def _metric_row(
    method: str,
    noise_sd: float,
    duration_s: float,
    distance: float | str,
    subset: list[dict[str, Any]],
    prediction_key: str,
) -> dict[str, Any]:
    predicted = [row[prediction_key] for row in subset]
    classified = [
        (row["model"], row[prediction_key])
        for row in subset
        if row[prediction_key] != INDETERMINATE
    ]
    correct = sum(actual == guess for actual, guess in classified)
    false = len(classified) - correct
    by_model = {
        model: sum(actual == model and guess == model for actual, guess in classified)
        / max(1, sum(actual == model for actual, _ in classified))
        for model in MODEL_NAMES
    }
    fold_true = [row for row in subset if row["model"] == "fold"]
    linear_true = [row for row in subset if row["model"] != "fold"]
    fold_classified = [row for row in fold_true if row[prediction_key] != INDETERMINATE]
    linear_classified = [
        row for row in linear_true if row[prediction_key] != INDETERMINATE
    ]
    return {
        "method": method,
        "noise_sd": noise_sd,
        "window_s": duration_s,
        "distance": distance,
        "trajectory_count": len(subset),
        "correct_count": correct,
        "false_classification_count": false,
        "abstention_count": sum(value == INDETERMINATE for value in predicted),
        "accuracy_all": correct / max(1, len(subset)),
        "accuracy_non_abstained": correct / max(1, len(classified)),
        "false_classification_rate_non_abstained": false / max(1, len(classified)),
        "abstention_rate": sum(value == INDETERMINATE for value in predicted)
        / max(1, len(subset)),
        "balanced_accuracy_non_abstained": float(np.mean(list(by_model.values()))),
        "fold_sensitivity_non_abstained": sum(
            row[prediction_key] == "fold" for row in fold_classified
        )
        / max(1, len(fold_classified)),
        "linear_specificity_non_abstained": sum(
            row[prediction_key] != "fold" for row in linear_classified
        )
        / max(1, len(linear_classified)),
    }


def _scaling_results(
    records: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    estimates: list[dict[str, Any]] = []
    summary: list[dict[str, Any]] = []
    evaluation = [r for r in records if r["partition"] == "evaluation"]
    for model in MODEL_NAMES:
        for noise_sd in NOISE_LEVELS:
            for duration_s in WINDOWS_S:
                seed_rows: list[dict[str, Any]] = []
                for seed in EVALUATION_SEEDS:
                    group = sorted(
                        [
                            r
                            for r in evaluation
                            if r["model"] == model
                            and r["noise_sd"] == noise_sd
                            and r["window_s"] == duration_s
                            and r["seed"] == seed
                        ],
                        key=lambda row: row["control_setting"],
                    )
                    if len(group) != len(CONTROL_SETTINGS) or any(
                        row["estimated_equilibrium"] <= 0.0
                        or row["estimated_rate"] <= 0.0
                        for row in group
                    ):
                        estimate = ScalingEstimate(
                            math.nan, math.nan, math.nan, math.nan, math.nan, False
                        )
                    else:
                        estimate = estimate_scaling(
                            np.array([row["control_setting"] for row in group]),
                            np.array([row["estimated_equilibrium"] for row in group]),
                            np.array([row["estimated_rate"] for row in group]),
                        )
                    record = {
                        "model": model,
                        "noise_sd": noise_sd,
                        "window_s": duration_s,
                        "seed": seed,
                        "estimated_equilibrium_exponent": estimate.equilibrium_exponent,
                        "estimated_rate_exponent": estimate.rate_exponent,
                        "estimated_boundary": estimate.estimated_boundary,
                        "equilibrium_log_rmse": estimate.equilibrium_log_rmse,
                        "rate_log_rmse": estimate.rate_log_rmse,
                        "converged": estimate.converged,
                    }
                    estimates.append(record)
                    seed_rows.append(record)
                valid = [r for r in seed_rows if r["converged"]]
                if not valid:
                    summary.append(
                        {
                            "model": model,
                            "noise_sd": noise_sd,
                            "window_s": duration_s,
                            "valid_seed_count": 0,
                            "true_equilibrium_exponent": _SIGNATURES[model][0],
                            "mean_equilibrium_exponent": math.nan,
                            "equilibrium_exponent_ci95_low": math.nan,
                            "equilibrium_exponent_ci95_high": math.nan,
                            "true_rate_exponent": _SIGNATURES[model][1],
                            "mean_rate_exponent": math.nan,
                            "rate_exponent_ci95_low": math.nan,
                            "rate_exponent_ci95_high": math.nan,
                        }
                    )
                    continue
                q_values = np.array(
                    [r["estimated_equilibrium_exponent"] for r in valid]
                )
                p_values = np.array([r["estimated_rate_exponent"] for r in valid])
                rng = np.random.default_rng(
                    53053
                    + MODEL_NAMES.index(model) * 100
                    + NOISE_LEVELS.index(noise_sd) * 10
                    + WINDOWS_S.index(duration_s)
                )
                boot_q = np.empty(BOOTSTRAP_REPLICATES)
                boot_p = np.empty(BOOTSTRAP_REPLICATES)
                for index in range(BOOTSTRAP_REPLICATES):
                    draw = rng.integers(0, len(valid), size=len(valid))
                    boot_q[index] = np.mean(q_values[draw])
                    boot_p[index] = np.mean(p_values[draw])
                q_ci = np.percentile(boot_q, [2.5, 97.5])
                p_ci = np.percentile(boot_p, [2.5, 97.5])
                q_true, p_true = _SIGNATURES[model]
                summary.append(
                    {
                        "model": model,
                        "noise_sd": noise_sd,
                        "window_s": duration_s,
                        "valid_seed_count": len(valid),
                        "true_equilibrium_exponent": q_true,
                        "mean_equilibrium_exponent": float(np.mean(q_values)),
                        "median_absolute_equilibrium_exponent_error": float(
                            np.median(np.abs(q_values - q_true))
                        ),
                        "equilibrium_exponent_ci95_low": float(q_ci[0]),
                        "equilibrium_exponent_ci95_high": float(q_ci[1]),
                        "true_rate_exponent": p_true,
                        "mean_rate_exponent": float(np.mean(p_values)),
                        "median_absolute_rate_exponent_error": float(
                            np.median(np.abs(p_values - p_true))
                        ),
                        "rate_exponent_ci95_low": float(p_ci[0]),
                        "rate_exponent_ci95_high": float(p_ci[1]),
                        "mean_estimated_boundary": float(
                            np.mean([r["estimated_boundary"] for r in valid])
                        ),
                        "true_boundary_for_scoring_only": CONTROL_BOUNDARY,
                    }
                )
    return estimates, summary


def _parameter_recovery_summary(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    primary = [
        r
        for r in records
        if r["partition"] == "evaluation"
        and r["noise_sd"] == min(NOISE_LEVELS)
        and r["window_s"] == max(WINDOWS_S)
    ]
    for model in MODEL_NAMES:
        for margin in TRUE_MARGINS:
            subset = [
                r for r in primary if r["model"] == model and r["true_margin"] == margin
            ]
            margin_errors = [
                abs(r["true_model_margin_fit"] - margin) / margin for r in subset
            ]
            eq_errors = [
                abs(r["true_model_equilibrium_fit"] - r["true_equilibrium"])
                / r["true_equilibrium"]
                for r in subset
            ]
            ag_eq_errors = [
                abs(r["estimated_equilibrium"] - r["true_equilibrium"])
                / r["true_equilibrium"]
                for r in subset
            ]
            rate_errors = [
                abs(r["estimated_rate"] - r["true_recovery_rate"])
                / r["true_recovery_rate"]
                for r in subset
            ]
            rows.append(
                {
                    "model": model,
                    "margin_for_scoring_only": margin,
                    "evaluation_trajectory_count": len(subset),
                    "model_informed_margin_median_relative_error": float(
                        np.median(margin_errors)
                    ),
                    "model_informed_equilibrium_median_relative_error": float(
                        np.median(eq_errors)
                    ),
                    "model_agnostic_equilibrium_median_relative_error": float(
                        np.median(ag_eq_errors)
                    ),
                    "model_agnostic_rate_median_relative_error": float(
                        np.median(rate_errors)
                    ),
                    "model_informed_margin_p90_relative_error": float(
                        np.percentile(margin_errors, 90)
                    ),
                    "model_informed_equilibrium_p90_relative_error": float(
                        np.percentile(eq_errors, 90)
                    ),
                }
            )
    return rows


def _evaluate_frozen_criteria(
    evaluation: list[dict[str, Any]],
    primary_metrics: list[dict[str, Any]],
    exponent_summary: list[dict[str, Any]],
    parameter_summary: list[dict[str, Any]],
) -> dict[str, Any]:
    informed = next(row for row in primary_metrics if row["method"] == "model_informed")
    agnostic = next(row for row in primary_metrics if row["method"] == "model_agnostic")
    informed_pass = (
        informed["accuracy_non_abstained"] >= 0.70
        and informed["abstention_rate"] <= 0.25
    )
    agnostic_pass = (
        agnostic["balanced_accuracy_non_abstained"] >= 0.70
        and agnostic["abstention_rate"] <= 0.25
        and agnostic["fold_sensitivity_non_abstained"] >= 0.75
        and agnostic["linear_specificity_non_abstained"] >= 0.75
    )
    primary_exponents = [
        row
        for row in exponent_summary
        if row["noise_sd"] == min(NOISE_LEVELS) and row["window_s"] == max(WINDOWS_S)
    ]
    exponent_checks = [
        abs(row["mean_equilibrium_exponent"] - row["true_equilibrium_exponent"]) <= 0.20
        and abs(row["mean_rate_exponent"] - row["true_rate_exponent"]) <= 0.20
        and row["equilibrium_exponent_ci95_low"]
        <= row["true_equilibrium_exponent"]
        <= row["equilibrium_exponent_ci95_high"]
        and row["rate_exponent_ci95_low"]
        <= row["true_rate_exponent"]
        <= row["rate_exponent_ci95_high"]
        for row in primary_exponents
    ]
    primary_parameters = parameter_summary
    parameter_checks = [
        row["model_informed_margin_median_relative_error"] <= 0.25
        and row["model_informed_equilibrium_median_relative_error"] <= 0.25
        for row in primary_parameters
    ]
    return {
        "primary_condition": {
            "noise_sd": min(NOISE_LEVELS),
            "window_s": max(WINDOWS_S),
        },
        "model_informed_three_class": {
            "passed": bool(informed_pass),
            "accuracy_non_abstained": informed["accuracy_non_abstained"],
            "abstention_rate": informed["abstention_rate"],
            "thresholds": {"accuracy_min": 0.70, "abstention_max": 0.25},
        },
        "model_agnostic_three_class": {
            "passed": bool(agnostic_pass),
            "balanced_accuracy_non_abstained": agnostic[
                "balanced_accuracy_non_abstained"
            ],
            "abstention_rate": agnostic["abstention_rate"],
            "fold_sensitivity_non_abstained": agnostic[
                "fold_sensitivity_non_abstained"
            ],
            "linear_specificity_non_abstained": agnostic[
                "linear_specificity_non_abstained"
            ],
            "thresholds": {
                "balanced_accuracy_min": 0.70,
                "abstention_max": 0.25,
                "sensitivity_min": 0.75,
                "specificity_min": 0.75,
            },
        },
        "scaling_identifiability_all_classes": {
            "passed": bool(primary_exponents and all(exponent_checks)),
            "per_class_checks": exponent_checks,
            "absolute_mean_exponent_error_max": 0.20,
            "require_ci_cover_truth": True,
        },
        "parameter_recovery_all_classes_and_margins": {
            "passed": bool(primary_parameters and all(parameter_checks)),
            "per_model_margin_checks": parameter_checks,
            "median_relative_error_max": 0.25,
        },
        "overall_method_success": bool(
            informed_pass
            and agnostic_pass
            and primary_exponents
            and all(exponent_checks)
        ),
        "overall_success_does_not_imply_real_system_validity": True,
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ScientificValidationError(f"refusing to write empty result file: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _hash_manifest(directory: Path) -> dict[str, str]:
    return {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(directory.iterdir())
        if path.is_file() and path.name != "sha256_manifest.json"
    }
