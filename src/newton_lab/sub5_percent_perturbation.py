"""Phase 55 synthetic study of small perturbations and recovery identifiability."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import NDArray

from newton_lab.exceptions import ScientificValidationError
from newton_lab.recovery_scaling_identifiability import (
    CONTROL_BOUNDARY,
    MODEL_NAMES,
    TRUE_MARGINS,
    estimate_exponential,
    estimate_scaling,
    exact_trajectory,
    local_recovery_rate,
    stable_equilibrium,
)

SUB5_AMPLITUDES = (0.0025, 0.005, 0.01, 0.02, 0.04)
REFERENCE_AMPLITUDES = (0.05, 0.20)
AMPLITUDES = SUB5_AMPLITUDES + REFERENCE_AMPLITUDES
NOISE_ARMS = {"noisy": 0.01, "zero_noise": 0.0}
WINDOW_TAU = 3.0
SAMPLE_INTERVAL_TAU = 0.1
SAMPLE_COUNT = 31
CALIBRATION_SEEDS = tuple(range(56000, 56008))
EVALUATION_SEEDS = tuple(range(57000, 57050))
BOOTSTRAP_REPLICATES = 2000
PROTOCOL_VERSION = "1.2"
SIGNATURES = {
    "fold": (0.5, 0.5),
    "transcritical": (1.0, 1.0),
    "supercritical_pitchfork": (0.5, 1.0),
}


def validate_protocol(protocol: dict[str, Any]) -> None:
    """Reject a protocol that disagrees with the frozen study implementation."""
    if (
        protocol.get("protocol_version") != PROTOCOL_VERSION
        or protocol.get("status")
        != "amended_and_frozen_before_final_held_out_evaluation"
    ):
        raise ScientificValidationError("Phase 55 protocol is not frozen/version 1.2")
    amplitude = protocol.get("amplitude", {})
    if (
        tuple(amplitude.get("sub_5_percent_primary_levels", ())) != SUB5_AMPLITUDES
        or tuple(amplitude.get("phase_54_reference_levels", ())) != REFERENCE_AMPLITUDES
        or tuple(amplitude.get("trend_fit_levels", ())) != SUB5_AMPLITUDES
        or amplitude.get("definition") != "a=delta/x*, x(0)=x*(1+a), delta=a*x*"
    ):
        raise ScientificValidationError("Phase 55 amplitude protocol mismatch")
    models = protocol.get("models", {})
    if tuple(models) != MODEL_NAMES or any(
        (models[name].get("q"), models[name].get("p")) != SIGNATURES[name]
        for name in MODEL_NAMES
    ):
        raise ScientificValidationError("Phase 55 models/signatures mismatch Phase 53")
    if tuple(protocol.get("margins", ())) != TRUE_MARGINS:
        raise ScientificValidationError("Phase 55 margins mismatch Phase 53")
    config = protocol.get("primary_observation", {})
    if (
        config.get("noise_sd_fraction_of_equilibrium") != NOISE_ARMS["noisy"]
        or config.get("window_tau") != WINDOW_TAU
        or config.get("sampling_interval_tau") != SAMPLE_INTERVAL_TAU
        or config.get("sample_count") != SAMPLE_COUNT
        or config.get("information_condition") != "reference"
    ):
        raise ScientificValidationError("Phase 55 primary observation mismatch")
    diagnostic = protocol.get("diagnostic_arm", {})
    if (
        diagnostic.get("only_changed_setting") != "observation_noise_sd=0"
        or diagnostic.get("estimator_scale_input_fraction_of_equilibrium_both_arms")
        != NOISE_ARMS["noisy"]
        or diagnostic.get("noise_is_not_added_to_zero_noise_observations") is not True
    ):
        raise ScientificValidationError("Phase 55 diagnostic-arm protocol mismatch")
    design = protocol.get("design", {})
    calibration = design.get("calibration_seeds", [])
    evaluation = design.get("evaluation_seeds", [])
    if (
        calibration != [CALIBRATION_SEEDS[0], CALIBRATION_SEEDS[-1]]
        or evaluation != [EVALUATION_SEEDS[0], EVALUATION_SEEDS[-1]]
        or not set(CALIBRATION_SEEDS).isdisjoint(EVALUATION_SEEDS)
        or design.get("bootstrap_replicates") != BOOTSTRAP_REPLICATES
    ):
        raise ScientificValidationError("Phase 55 seed/bootstrap protocol mismatch")
    controls = protocol.get("controls", [])
    if len(controls) != len(TRUE_MARGINS) or any(
        not math.isclose(value, CONTROL_BOUNDARY + margin, abs_tol=1e-12)
        for value, margin in zip(controls, TRUE_MARGINS, strict=True)
    ):
        raise ScientificValidationError("Phase 55 controls mismatch Phase 53")


def amplitude_specification() -> list[dict[str, Any]]:
    """Return the exact displacement and displacement/margin for each setting."""
    rows: list[dict[str, Any]] = []
    for model, margin, amplitude in itertools.product(
        MODEL_NAMES, TRUE_MARGINS, AMPLITUDES
    ):
        equilibrium = stable_equilibrium(model, margin)
        displacement = amplitude * equilibrium
        rows.append(
            {
                "model": model,
                "margin": margin,
                "control_setting": CONTROL_BOUNDARY + margin,
                "amplitude_fraction_of_equilibrium": amplitude,
                "true_equilibrium": equilibrium,
                "initial_displacement": displacement,
                "displacement_over_margin": displacement / margin,
                "initial_state": equilibrium + displacement,
                "is_sub_5_percent_primary": amplitude in SUB5_AMPLITUDES,
                "is_phase_54_reference": amplitude in REFERENCE_AMPLITUDES,
            }
        )
    return rows


def signal_to_noise(
    clean: NDArray[np.float64], equilibrium: float, noise_sd: float
) -> float:
    """Compute RMS clean recovery excursion divided by observation-noise SD.

    The equilibrium and clean trajectory are generator truth used only for this
    diagnostic. Zero noise has no finite SNR and is represented by infinity.
    """
    values = np.asarray(clean, dtype=np.float64)
    if (
        values.ndim != 1
        or values.size == 0
        or not np.isfinite(values).all()
        or not math.isfinite(equilibrium)
        or not math.isfinite(noise_sd)
        or noise_sd < 0.0
    ):
        raise ScientificValidationError("invalid signal-to-noise inputs")
    if noise_sd == 0.0:
        return math.inf
    rms = float(np.sqrt(np.mean(np.square(values - equilibrium))))
    return rms / noise_sd


def paired_innovations(seed: int, model: str, margin: float) -> NDArray[np.float64]:
    """Generate a fixed innovation sequence reused across paired amplitudes."""
    if model not in MODEL_NAMES or margin not in TRUE_MARGINS:
        raise ScientificValidationError("unknown model or margin for paired noise")
    model_index = MODEL_NAMES.index(model)
    margin_index = TRUE_MARGINS.index(margin)
    child = int(
        np.random.SeedSequence([seed, model_index, margin_index]).generate_state(1)[0]
    )
    return np.random.default_rng(child).standard_normal(SAMPLE_COUNT, dtype=np.float64)


def paired_observation(
    model: str, margin: float, amplitude: float, arm: str, seed: int
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], float]:
    """Construct one clean/noisy trace under the fixed Phase 55 settings."""
    if model not in MODEL_NAMES or margin not in TRUE_MARGINS:
        raise ScientificValidationError("unknown Phase 53 model or margin")
    if amplitude not in AMPLITUDES or arm not in NOISE_ARMS:
        raise ScientificValidationError("unknown Phase 55 amplitude or noise arm")
    equilibrium = stable_equilibrium(model, margin)
    rate = local_recovery_rate(model, margin)
    tau = 1.0 / rate
    times = np.arange(SAMPLE_COUNT, dtype=np.float64) * SAMPLE_INTERVAL_TAU * tau
    initial = equilibrium * (1.0 + amplitude)
    clean = exact_trajectory(model, margin, initial, times)
    noise_sd = NOISE_ARMS[arm] * equilibrium
    observed = clean + noise_sd * paired_innovations(seed, model, margin)
    snr = signal_to_noise(clean, equilibrium, noise_sd)
    return times, clean, observed, snr


def estimate_rate(
    model: str, margin: float, amplitude: float, arm: str, seed: int
) -> dict[str, Any]:
    """Fit Phase 53's free-offset exponential without supplying true equilibrium."""
    times, clean, observed, snr = paired_observation(
        model, margin, amplitude, arm, seed
    )
    equilibrium = stable_equilibrium(model, margin)
    true_rate = local_recovery_rate(model, margin)
    noise_sd = NOISE_ARMS[arm] * equilibrium
    # Keep the Phase 54 estimator initialization scale identical in both
    # diagnostic arms. It is metadata for initialization, not added noise.
    estimator_noise_scale = NOISE_ARMS["noisy"] * equilibrium
    fit = estimate_exponential(times, observed, estimator_noise_scale)
    valid = bool(
        fit.converged and math.isfinite(fit.recovery_rate) and fit.recovery_rate > 0.0
    )
    initial_displacement = amplitude * equilibrium
    relative_signed_error = (
        (fit.recovery_rate - true_rate) / true_rate if valid else math.nan
    )
    return {
        "estimated_equilibrium": fit.equilibrium,
        "equilibrium_standard_error": fit.equilibrium_standard_error,
        "estimated_rate": fit.recovery_rate,
        "rate_standard_error": fit.recovery_rate_standard_error,
        "rate_fit_rmse": fit.fit_rmse,
        "rate_fit_rmse_over_initial_displacement": fit.fit_rmse / initial_displacement,
        "rate_fit_r_squared": fit.r_squared,
        "rate_valid": valid,
        "true_equilibrium": equilibrium,
        "true_rate": true_rate,
        "rate_signed_error": fit.recovery_rate - true_rate if valid else math.nan,
        "rate_absolute_error": abs(fit.recovery_rate - true_rate)
        if valid
        else math.nan,
        "rate_relative_signed_error": relative_signed_error,
        "rate_relative_absolute_error": abs(relative_signed_error)
        if valid
        else math.nan,
        "initial_state": equilibrium + initial_displacement,
        "initial_displacement": initial_displacement,
        "amplitude_fraction_of_equilibrium": amplitude,
        "displacement_over_margin": initial_displacement / margin,
        "noise_sd_state_units": noise_sd,
        "estimator_noise_scale_state_units": estimator_noise_scale,
        "snr_rms_clean_excursion_over_noise_sd": snr,
        "duration": WINDOW_TAU / true_rate,
        "sample_count": times.size,
    }


def fit_exponent_group(
    model: str, arm: str, amplitude: float, seed: int, records: list[dict[str, Any]]
) -> dict[str, Any]:
    """Estimate q and p from five reference-information rate/equilibrium fits."""
    ordered = sorted(records, key=lambda row: row["margin"])
    valid_rates = all(row["rate_valid"] for row in ordered)
    rate_valid_fraction = sum(row["rate_valid"] for row in ordered) / len(ordered)
    q = p = boundary = math.nan
    converged = False
    if valid_rates and len(ordered) == len(TRUE_MARGINS):
        settings = np.asarray(
            [row["control_setting"] for row in ordered], dtype=np.float64
        )
        equilibria = np.asarray(
            [row["estimated_equilibrium"] for row in ordered], dtype=np.float64
        )
        rates = np.asarray([row["estimated_rate"] for row in ordered], dtype=np.float64)
        try:
            fit = estimate_scaling(settings, equilibria, rates)
            q, p, boundary, converged = (
                fit.equilibrium_exponent,
                fit.rate_exponent,
                fit.estimated_boundary,
                fit.converged,
            )
        except (ScientificValidationError, ValueError, FloatingPointError):
            pass
    valid = bool(converged and math.isfinite(q) and math.isfinite(p))
    true_q, true_p = SIGNATURES[model]
    return {
        "partition": "",
        "replicate_seed": seed,
        "model": model,
        "noise_arm": arm,
        "amplitude_fraction_of_equilibrium": amplitude,
        "attempted_margin_count": len(ordered),
        "rate_valid_fraction": rate_valid_fraction,
        "valid": valid,
        "estimated_boundary": boundary,
        "mean_estimated_equilibrium": float(
            np.mean([r["estimated_equilibrium"] for r in ordered])
        ),
        "q": q if valid else math.nan,
        "true_q": true_q,
        "q_bias": q - true_q if valid else math.nan,
        "q_absolute_error": abs(q - true_q) if valid else math.nan,
        "p": p if valid else math.nan,
        "true_p": true_p,
        "p_bias": p - true_p if valid else math.nan,
        "p_absolute_error": abs(p - true_p) if valid else math.nan,
    }


def run_study(output_dir: Path) -> dict[str, Any]:
    """Run frozen calibration/evaluation traces and write deterministic artifacts."""
    output_dir.mkdir(parents=True, exist_ok=True)
    protocol_path = output_dir / "protocol.json"
    try:
        protocol = json.loads(protocol_path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificValidationError(
            f"cannot read frozen Phase 55 protocol: {error}"
        ) from error
    validate_protocol(protocol)
    protocol_hash = hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    trajectories: list[dict[str, Any]] = []
    groups: list[dict[str, Any]] = []
    partitions = (("calibration", CALIBRATION_SEEDS), ("evaluation", EVALUATION_SEEDS))

    for partition, seeds in partitions:
        for seed, model, amplitude, arm in itertools.product(
            seeds, MODEL_NAMES, AMPLITUDES, NOISE_ARMS
        ):
            block: list[dict[str, Any]] = []
            for margin in TRUE_MARGINS:
                result = estimate_rate(model, margin, amplitude, arm, seed)
                row = {
                    "partition": partition,
                    "replicate_seed": seed,
                    "model": model,
                    "margin": margin,
                    "control_setting": CONTROL_BOUNDARY + margin,
                    "noise_arm": arm,
                    **result,
                }
                trajectories.append(row)
                block.append(row)
            group = fit_exponent_group(model, arm, amplitude, seed, block)
            group["partition"] = partition
            groups.append(group)

    evaluation_trajectories = _evaluation_only(trajectories)
    evaluation_groups = _evaluation_only(groups)
    rate_summary = _rate_summary(evaluation_trajectories)
    exponent_summary = _exponent_summary(evaluation_groups)
    trend_summary, trend_draws = _trend_summary(
        evaluation_trajectories, evaluation_groups, exponent_summary
    )
    paired_summary = _paired_arm_summary(
        evaluation_trajectories, evaluation_groups, exponent_summary
    )
    model_contrasts = _model_trend_contrasts(trend_draws)

    _write_csv(output_dir / "amplitude_specification.csv", amplitude_specification())
    _write_csv(output_dir / "per_trajectory_results.csv", trajectories)
    _write_csv(output_dir / "per_group_results.csv", groups)
    _write_csv(output_dir / "rate_condition_summary.csv", rate_summary)
    _write_csv(output_dir / "exponent_condition_summary.csv", exponent_summary)
    _write_csv(output_dir / "amplitude_trend_summary.csv", trend_summary)
    _write_csv(output_dir / "paired_arm_contrasts.csv", paired_summary)
    _write_csv(output_dir / "model_trend_contrasts.csv", model_contrasts)
    metadata = {
        "phase": 55,
        "protocol_version": PROTOCOL_VERSION,
        "protocol_sha256": protocol_hash,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "models": list(MODEL_NAMES),
        "margins_per_model": len(TRUE_MARGINS),
        "amplitude_levels": len(AMPLITUDES),
        "sub_5_percent_levels": list(SUB5_AMPLITUDES),
        "reference_levels": list(REFERENCE_AMPLITUDES),
        "noise_arms": NOISE_ARMS,
        "calibration_seed_count": len(CALIBRATION_SEEDS),
        "evaluation_seed_count": len(EVALUATION_SEEDS),
        "calibration_trajectory_count": sum(
            row["partition"] == "calibration" for row in trajectories
        ),
        "evaluation_trajectory_count": len(evaluation_trajectories),
        "evaluation_group_count": len(evaluation_groups),
        "valid_evaluation_rate_fits": sum(
            row["rate_valid"] for row in evaluation_trajectories
        ),
        "valid_evaluation_exponent_groups": sum(
            row["valid"] for row in evaluation_groups
        ),
        "trajectory_count": len(trajectories),
        "interpretation_limit": (
            "Controlled synthetic method-capability results only; no real-system, "
            "engineering, financial, or predictive validation."
        ),
    }
    (output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = _manifest(output_dir)
    (output_dir / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata


def _evaluation_only(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Select held-out rows explicitly; calibration rows never enter summaries."""
    return [row for row in rows if row.get("partition") == "evaluation"]


def _rate_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for model, margin, amplitude, arm in itertools.product(
        MODEL_NAMES, TRUE_MARGINS, AMPLITUDES, NOISE_ARMS
    ):
        cell = [
            row
            for row in rows
            if row["model"] == model
            and row["margin"] == margin
            and row["amplitude_fraction_of_equilibrium"] == amplitude
            and row["noise_arm"] == arm
        ]
        valid = [row for row in cell if row["rate_valid"]]
        errors = [row["rate_relative_absolute_error"] for row in valid]
        signed = [row["rate_relative_signed_error"] for row in valid]
        estimates = [row["estimated_rate"] for row in valid]
        rmse = [row["rate_fit_rmse_over_initial_displacement"] for row in cell]
        result.append(
            {
                "model": model,
                "margin": margin,
                "control_setting": CONTROL_BOUNDARY + margin,
                "amplitude_fraction_of_equilibrium": amplitude,
                "noise_arm": arm,
                "attempted_trajectories": len(cell),
                "valid_rate_fits": len(valid),
                "rate_fit_valid_fraction": len(valid) / len(cell) if cell else 0.0,
                "true_rate": local_recovery_rate(model, margin),
                "mean_estimated_rate": _mean(estimates),
                "median_estimated_rate": _median(estimates),
                "rate_mean_ci95_low": _bootstrap_interval(
                    estimates, _stable_seed(model, margin, amplitude, arm, 1)
                )[0],
                "rate_mean_ci95_high": _bootstrap_interval(
                    estimates, _stable_seed(model, margin, amplitude, arm, 1)
                )[1],
                "mean_signed_relative_rate_bias": _mean(signed),
                "mean_absolute_relative_rate_error": _mean(errors),
                "rate_abs_error_ci95_low": _bootstrap_interval(
                    errors, _stable_seed(model, margin, amplitude, arm, 2)
                )[0],
                "rate_abs_error_ci95_high": _bootstrap_interval(
                    errors, _stable_seed(model, margin, amplitude, arm, 2)
                )[1],
                "mean_normalized_fit_rmse": _mean(rmse),
                "mean_rate_fit_standard_error": _mean(
                    [row["rate_standard_error"] for row in valid]
                ),
                "mean_snr": (
                    math.inf
                    if arm == "zero_noise"
                    else _mean(
                        [row["snr_rms_clean_excursion_over_noise_sd"] for row in cell]
                    )
                ),
            }
        )
    return result


def _exponent_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for model, amplitude, arm in itertools.product(MODEL_NAMES, AMPLITUDES, NOISE_ARMS):
        cell = [
            row
            for row in rows
            if row["model"] == model
            and row["amplitude_fraction_of_equilibrium"] == amplitude
            and row["noise_arm"] == arm
        ]
        valid = [row for row in cell if row["valid"]]
        q_values = [row["q"] for row in valid]
        p_values = [row["p"] for row in valid]
        q_ci = _bootstrap_interval(q_values, _stable_seed(model, amplitude, arm, 3))
        p_ci = _bootstrap_interval(p_values, _stable_seed(model, amplitude, arm, 4))
        true_q, true_p = SIGNATURES[model]
        q_mean, p_mean = _mean(q_values), _mean(p_values)
        valid_rate = len(valid) / len(cell) if cell else 0.0
        q_cover = _covers(q_ci, true_q)
        p_cover = _covers(p_ci, true_p)
        joint_pass = bool(
            valid_rate >= 0.8
            and math.isfinite(q_mean)
            and math.isfinite(p_mean)
            and abs(q_mean - true_q) <= 0.2
            and abs(p_mean - true_p) <= 0.2
            and q_cover
            and p_cover
        )
        p_pass = bool(
            valid_rate >= 0.8
            and math.isfinite(p_mean)
            and abs(p_mean - true_p) <= 0.2
            and p_cover
        )
        result.append(
            {
                "model": model,
                "amplitude_fraction_of_equilibrium": amplitude,
                "noise_arm": arm,
                "attempted_groups": len(cell),
                "valid_exponent_groups": len(valid),
                "valid_exponent_group_fraction": valid_rate,
                "true_q": true_q,
                "mean_q": q_mean,
                "q_bias": q_mean - true_q if math.isfinite(q_mean) else math.nan,
                "q_absolute_error": abs(q_mean - true_q)
                if math.isfinite(q_mean)
                else math.nan,
                "q_ci95_low": q_ci[0],
                "q_ci95_high": q_ci[1],
                "q_ci95_width": q_ci[1] - q_ci[0],
                "q_ci_covers_truth": q_cover,
                "true_p": true_p,
                "mean_p": p_mean,
                "p_bias": p_mean - true_p if math.isfinite(p_mean) else math.nan,
                "p_absolute_error": abs(p_mean - true_p)
                if math.isfinite(p_mean)
                else math.nan,
                "p_ci95_low": p_ci[0],
                "p_ci95_high": p_ci[1],
                "p_ci95_width": p_ci[1] - p_ci[0],
                "p_ci_covers_truth": p_cover,
                "ci_coverage_interpretable": arm == "noisy",
                "ci_interpretation": (
                    "seed-bootstrap coverage"
                    if arm == "noisy"
                    else (
                        "degenerate deterministic zero-noise seed interval; "
                        "coverage is not inferential"
                    )
                ),
                "mean_estimated_boundary": _mean(
                    [row["estimated_boundary"] for row in cell]
                ),
                "mean_equilibrium_rate_valid_fraction": _mean(
                    [row["rate_valid_fraction"] for row in cell]
                ),
                "identifiability_criterion_applicability": (
                    "primary" if arm == "noisy" else "diagnostic_only_not_scored"
                ),
                "phase_54_joint_identifiability_criterion": (
                    joint_pass if arm == "noisy" else None
                ),
                "p_specific_identifiability_criterion": p_pass
                if arm == "noisy"
                else None,
            }
        )
    return result


def _trend_summary(
    trajectories: list[dict[str, Any]],
    groups: list[dict[str, Any]],
    exponent_summary: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[tuple[str, str, str], NDArray[np.float64]]]:
    result: list[dict[str, Any]] = []
    slope_draws: dict[tuple[str, str, str], NDArray[np.float64]] = {}
    log_amplitudes = np.log10(np.asarray(SUB5_AMPLITUDES, dtype=np.float64))
    for model, arm in itertools.product(MODEL_NAMES, NOISE_ARMS):
        rates_by_key = {
            (
                r["replicate_seed"],
                r["amplitude_fraction_of_equilibrium"],
                r["margin"],
            ): r
            for r in trajectories
            if r["model"] == model and r["noise_arm"] == arm
        }
        groups_by_key = {
            (r["replicate_seed"], r["amplitude_fraction_of_equilibrium"]): r
            for r in groups
            if r["model"] == model and r["noise_arm"] == arm
        }
        metric_rows: dict[str, list[dict[float, float]]] = {
            metric: []
            for metric in (
                "rate_absolute_relative_error",
                "rate_signed_relative_bias",
                "rate_valid_fraction",
                "normalized_fit_rmse",
                "snr",
                "p_absolute_error",
                "p_valid_fraction",
                "p_ci95_width",
                "p_ci95_coverage",
            )
        }
        for seed in EVALUATION_SEEDS:
            metrics_for_seed: dict[str, dict[float, float]] = {
                name: {} for name in metric_rows
            }
            for amplitude in SUB5_AMPLITUDES:
                block = [
                    rates_by_key[(seed, amplitude, margin)] for margin in TRUE_MARGINS
                ]
                valid = [row for row in block if row["rate_valid"]]
                metrics_for_seed["rate_valid_fraction"][amplitude] = len(valid) / len(
                    block
                )
                for metric, field in (
                    ("rate_absolute_relative_error", "rate_relative_absolute_error"),
                    ("rate_signed_relative_bias", "rate_relative_signed_error"),
                    ("normalized_fit_rmse", "rate_fit_rmse_over_initial_displacement"),
                    ("snr", "snr_rms_clean_excursion_over_noise_sd"),
                ):
                    values = [
                        row[field] for row in (valid if metric != "snr" else block)
                    ]
                    metrics_for_seed[metric][amplitude] = _mean(values)
                group = groups_by_key[(seed, amplitude)]
                metrics_for_seed["p_valid_fraction"][amplitude] = float(group["valid"])
                metrics_for_seed["p_absolute_error"][amplitude] = group[
                    "p_absolute_error"
                ]
            for metric in metric_rows:
                metric_rows[metric].append(metrics_for_seed[metric])

        for metric, per_seed in metric_rows.items():
            if metric == "snr" and arm == "zero_noise":
                slope = math.nan
                draws = np.asarray([], dtype=np.float64)
                trend_seed_count = 0
            elif metric == "snr":
                # SNR is fixed by the generator condition, so replicate seeds
                # are not independent observations for an inferential interval.
                slope = _deterministic_snr_slope(trajectories, model)
                draws = np.asarray([], dtype=np.float64)
                trend_seed_count = 0
            elif metric in ("p_ci95_width", "p_ci95_coverage"):
                summary_by_amp = {
                    row["amplitude_fraction_of_equilibrium"]: row
                    for row in exponent_summary
                    if row["model"] == model and row["noise_arm"] == arm
                }
                trend_values = np.asarray(
                    [
                        summary_by_amp[a][
                            "p_ci95_width"
                            if metric == "p_ci95_width"
                            else "p_ci_covers_truth"
                        ]
                        for a in SUB5_AMPLITUDES
                    ],
                    dtype=np.float64,
                )
                slope = _linear_slope(log_amplitudes, trend_values)
                draws = _bootstrap_ci_metric_trend(
                    [
                        groups_by_key[(seed, a)]["p"]
                        if groups_by_key[(seed, a)]["valid"]
                        else math.nan
                        for seed in EVALUATION_SEEDS
                        for a in SUB5_AMPLITUDES
                    ],
                    len(EVALUATION_SEEDS),
                    len(SUB5_AMPLITUDES),
                    SIGNATURES[model][1],
                    metric,
                    _stable_seed(model, arm, metric, 90),
                )
                trend_seed_count = min(
                    sum(groups_by_key[(seed, a)]["valid"] for seed in EVALUATION_SEEDS)
                    for a in SUB5_AMPLITUDES
                )
            else:
                matrix = np.asarray(
                    [
                        [row.get(a, math.nan) for a in SUB5_AMPLITUDES]
                        for row in per_seed
                    ],
                    dtype=np.float64,
                )
                seed_slopes = np.asarray(
                    [_row_slope(log_amplitudes, row) for row in matrix],
                    dtype=np.float64,
                )
                trend_seed_count = int(np.isfinite(seed_slopes).sum())
                finite = seed_slopes[np.isfinite(seed_slopes)]
                slope = _mean(finite.tolist())
                draws = _bootstrap_mean_draws(
                    finite, _stable_seed(model, arm, metric, 91)
                )
            low, high = _draw_interval(draws)
            slope_draws[(model, arm, metric)] = draws
            result.append(
                {
                    "model": model,
                    "noise_arm": arm,
                    "metric": metric,
                    "trend_levels": ";".join(str(a) for a in SUB5_AMPLITUDES),
                    "slope_per_log10_amplitude": slope,
                    "slope_ci95_low": low,
                    "slope_ci95_high": high,
                    "contributing_seed_count": trend_seed_count,
                    "trend_status": (
                        "not_applicable_no_measurement_noise"
                        if metric == "snr" and arm == "zero_noise"
                        else (
                            "descriptive_deterministic_no_sampling_uncertainty"
                            if metric == "snr" or arm == "zero_noise"
                            else _effect_status(low, high)
                        )
                    ),
                    "interpretation": (
                        "SNR is undefined without measurement noise"
                        if metric == "snr" and arm == "zero_noise"
                        else (
                            "seed-level sub-5% trend; conditional error trends "
                            "retain validity denominators"
                        )
                    ),
                }
            )
    return result, slope_draws


def _deterministic_snr_slope(trajectories: list[dict[str, Any]], model: str) -> float:
    """Return the exact descriptive noisy-arm SNR trend across primary levels."""
    amplitude_means = np.asarray(
        [
            _mean(
                [
                    float(row["snr_rms_clean_excursion_over_noise_sd"])
                    for row in trajectories
                    if row["model"] == model
                    and row["noise_arm"] == "noisy"
                    and row["partition"] == "evaluation"
                    and row["amplitude_fraction_of_equilibrium"] == amplitude
                ]
            )
            for amplitude in SUB5_AMPLITUDES
        ],
        dtype=np.float64,
    )
    return _linear_slope(
        np.log10(np.asarray(SUB5_AMPLITUDES, dtype=np.float64)), amplitude_means
    )


def _paired_arm_summary(
    trajectories: list[dict[str, Any]],
    groups: list[dict[str, Any]],
    exponent_summary: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for model, amplitude in itertools.product(MODEL_NAMES, AMPLITUDES):
        noisy = {
            (r["replicate_seed"], r["margin"]): r
            for r in trajectories
            if r["model"] == model
            and r["noise_arm"] == "noisy"
            and r["amplitude_fraction_of_equilibrium"] == amplitude
        }
        quiet = {
            (r["replicate_seed"], r["margin"]): r
            for r in trajectories
            if r["model"] == model
            and r["noise_arm"] == "zero_noise"
            and r["amplitude_fraction_of_equilibrium"] == amplitude
        }
        metrics: dict[str, list[float]] = {
            "rate_absolute_relative_error_noisy_minus_zero": [],
            "rate_valid_fraction_noisy_minus_zero": [],
        }
        for seed in EVALUATION_SEEDS:
            pair = [(noisy[(seed, m)], quiet[(seed, m)]) for m in TRUE_MARGINS]
            paired_valid = [
                (a, b) for a, b in pair if a["rate_valid"] and b["rate_valid"]
            ]
            metrics["rate_absolute_relative_error_noisy_minus_zero"].append(
                _mean(
                    [
                        a["rate_relative_absolute_error"]
                        - b["rate_relative_absolute_error"]
                        for a, b in paired_valid
                    ]
                )
            )
            metrics["rate_valid_fraction_noisy_minus_zero"].append(
                sum(a["rate_valid"] - b["rate_valid"] for a, b in pair) / len(pair)
            )

        group_arm = {
            arm: {
                row["replicate_seed"]: row
                for row in groups
                if row["model"] == model
                and row["noise_arm"] == arm
                and row["amplitude_fraction_of_equilibrium"] == amplitude
            }
            for arm in NOISE_ARMS
        }
        exponent_diffs: list[float] = []
        valid_diffs: list[float] = []
        matched_exponents = 0
        for seed in EVALUATION_SEEDS:
            a = group_arm["noisy"][seed]
            b = group_arm["zero_noise"][seed]
            valid_diffs.append(float(a["valid"]) - float(b["valid"]))
            if a["valid"] and b["valid"]:
                exponent_diffs.append(a["p_absolute_error"] - b["p_absolute_error"])
                matched_exponents += 1
        metrics["p_absolute_error_noisy_minus_zero"] = exponent_diffs
        metrics["p_valid_fraction_noisy_minus_zero"] = valid_diffs
        for metric, values in metrics.items():
            low, high = _bootstrap_interval(
                values, _stable_seed(model, amplitude, metric, 77)
            )
            rows.append(
                {
                    "model": model,
                    "amplitude_fraction_of_equilibrium": amplitude,
                    "contrast": metric,
                    "estimate": _mean(values),
                    "ci95_low": low,
                    "ci95_high": high,
                    "matched_seed_count": sum(math.isfinite(value) for value in values),
                    "matched_exponent_seed_count": matched_exponents
                    if metric.startswith("p_")
                    else len(values),
                    "effect_status": _effect_status(low, high),
                }
            )

        # Width and coverage contrasts use paired seed resampling. Within each
        # resample the interval is the percentile interval of bootstrap means;
        # a normal approximation provides a deterministic, low-cost outer
        # bootstrap statistic for its width and truth coverage.
        for metric in (
            "p_ci95_width_noisy_minus_zero",
            "p_ci_coverage_noisy_minus_zero",
        ):
            draw = _bootstrap_ci_arm_contrast(
                group_arm["noisy"],
                group_arm["zero_noise"],
                SIGNATURES[model][1],
                metric,
                _stable_seed(model, amplitude, metric, 78),
            )
            low, high = _draw_interval(draw)
            noisy_row = next(
                r
                for r in exponent_summary
                if r["model"] == model
                and r["noise_arm"] == "noisy"
                and r["amplitude_fraction_of_equilibrium"] == amplitude
            )
            quiet_row = next(
                r
                for r in exponent_summary
                if r["model"] == model
                and r["noise_arm"] == "zero_noise"
                and r["amplitude_fraction_of_equilibrium"] == amplitude
            )
            field = "p_ci95_width" if "width" in metric else "p_ci_covers_truth"
            estimate = float(noisy_row[field]) - float(quiet_row[field])
            coverage_contrast = "coverage" in metric
            rows.append(
                {
                    "model": model,
                    "amplitude_fraction_of_equilibrium": amplitude,
                    "contrast": metric,
                    "estimate": math.nan if coverage_contrast else estimate,
                    "ci95_low": math.nan if coverage_contrast else low,
                    "ci95_high": math.nan if coverage_contrast else high,
                    "matched_seed_count": len(EVALUATION_SEEDS),
                    "matched_exponent_seed_count": len(EVALUATION_SEEDS),
                    "effect_status": (
                        "not_interpretable_degenerate_zero_noise_interval"
                        if coverage_contrast
                        else _effect_status(low, high)
                    ),
                }
            )
    return rows


def _model_trend_contrasts(
    draws: dict[tuple[str, str, str], NDArray[np.float64]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    metrics = (
        "rate_absolute_relative_error",
        "rate_valid_fraction",
        "p_absolute_error",
        "p_valid_fraction",
    )
    for arm, metric in itertools.product(NOISE_ARMS, metrics):
        for i, first in enumerate(MODEL_NAMES):
            for second in MODEL_NAMES[i + 1 :]:
                a = draws[(first, arm, metric)]
                b = draws[(second, arm, metric)]
                count = min(a.size, b.size)
                diff = a[:count] - b[:count]
                low, high = _draw_interval(diff)
                rows.append(
                    {
                        "noise_arm": arm,
                        "metric": metric,
                        "first_model": first,
                        "second_model": second,
                        "slope_difference_first_minus_second": _mean(diff.tolist()),
                        "ci95_low": low,
                        "ci95_high": high,
                        "resampling": "independent model-specific seed bootstrap",
                        "effect_status": (
                            "descriptive_deterministic_no_sampling_uncertainty"
                            if arm == "zero_noise"
                            else _effect_status(low, high)
                        ),
                    }
                )
    return rows


def _bootstrap_ci_arm_contrast(
    first: dict[int, dict[str, Any]],
    second: dict[int, dict[str, Any]],
    truth: float,
    metric: str,
    seed: int,
) -> NDArray[np.float64]:
    rng = np.random.default_rng(seed)
    seeds = list(EVALUATION_SEEDS)
    output = np.full(BOOTSTRAP_REPLICATES, np.nan, dtype=np.float64)
    for i in range(BOOTSTRAP_REPLICATES):
        selected = rng.integers(0, len(seeds), len(seeds))
        stats: list[float] = []
        for arm_rows in (first, second):
            values = np.asarray(
                [
                    arm_rows[seeds[j]]["p"]
                    for j in selected
                    if arm_rows[seeds[j]]["valid"]
                ],
                dtype=np.float64,
            )
            if values.size < 2:
                stats.append(math.nan)
                continue
            mean = float(np.mean(values))
            se = float(np.std(values, ddof=1) / math.sqrt(values.size))
            low, high = mean - 1.96 * se, mean + 1.96 * se
            stats.append(
                high - low if "width" in metric else float(low <= truth <= high)
            )
        if all(math.isfinite(value) for value in stats):
            output[i] = stats[0] - stats[1]
    return output[np.isfinite(output)]


def _bootstrap_ci_metric_trend(
    flattened: list[float],
    seed_count: int,
    amplitude_count: int,
    truth: float,
    metric: str,
    seed: int,
) -> NDArray[np.float64]:
    """Bootstrap seed-level trends in interval width/coverage diagnostics."""
    matrix = np.asarray(flattened, dtype=np.float64).reshape(
        seed_count, amplitude_count
    )
    rng = np.random.default_rng(seed)
    x = np.log10(np.asarray(SUB5_AMPLITUDES, dtype=np.float64))
    output = np.full(BOOTSTRAP_REPLICATES, np.nan, dtype=np.float64)
    for i in range(BOOTSTRAP_REPLICATES):
        sample = matrix[rng.integers(0, seed_count, seed_count)]
        vals: list[float] = []
        for col in range(amplitude_count):
            cell = sample[:, col]
            cell = cell[np.isfinite(cell)]
            if cell.size < 2:
                vals.append(math.nan)
                continue
            mean = float(np.mean(cell))
            se = float(np.std(cell, ddof=1) / math.sqrt(cell.size))
            lo, hi = mean - 1.96 * se, mean + 1.96 * se
            vals.append(hi - lo if "width" in metric else float(lo <= truth <= hi))
        output[i] = _linear_slope(x, np.asarray(vals, dtype=np.float64))
    return output[np.isfinite(output)]


def _bootstrap_mean_draws(
    values: NDArray[np.float64], seed: int
) -> NDArray[np.float64]:
    values = values[np.isfinite(values)]
    if values.size < 2:
        return np.asarray([], dtype=np.float64)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, values.size, (BOOTSTRAP_REPLICATES, values.size))
    return np.mean(values[indices], axis=1)


def _bootstrap_interval(values: list[float], seed: int) -> tuple[float, float]:
    array = np.asarray(values, dtype=np.float64)
    draws = _bootstrap_mean_draws(array, seed)
    return _draw_interval(draws)


def _draw_interval(draws: NDArray[np.float64]) -> tuple[float, float]:
    if draws.size < 2:
        return math.nan, math.nan
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(low), float(high)


def _stable_seed(*parts: Any) -> int:
    raw = json.dumps(parts, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return int.from_bytes(hashlib.sha256(raw).digest()[:4], "little")


def _mean(values: list[float]) -> float:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    return float(np.mean(finite)) if finite else math.nan


def _median(values: list[float]) -> float:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    return float(np.median(finite)) if finite else math.nan


def _covers(interval: tuple[float, float], truth: float) -> bool:
    low, high = interval
    return bool(math.isfinite(low) and math.isfinite(high) and low <= truth <= high)


def _row_slope(x: NDArray[np.float64], y: NDArray[np.float64]) -> float:
    valid = np.isfinite(y)
    return _linear_slope(x[valid], y[valid]) if int(valid.sum()) >= 3 else math.nan


def _linear_slope(x: NDArray[np.float64], y: NDArray[np.float64]) -> float:
    valid = np.isfinite(x) & np.isfinite(y)
    if int(valid.sum()) < 3:
        return math.nan
    return float(np.polyfit(x[valid], y[valid], 1)[0])


def _effect_status(low: float, high: float) -> str:
    if not math.isfinite(low) or not math.isfinite(high):
        return "indeterminate_insufficient_valid_data"
    return (
        "supported_interval_excludes_zero"
        if high < 0.0 or low > 0.0
        else "inconclusive_interval_includes_zero"
    )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ScientificValidationError(
            f"refusing to write empty Phase 55 CSV: {path.name}"
        )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _manifest(directory: Path) -> dict[str, str]:
    return {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(directory.iterdir())
        if path.is_file() and path.name != "sha256_manifest.json"
    }
