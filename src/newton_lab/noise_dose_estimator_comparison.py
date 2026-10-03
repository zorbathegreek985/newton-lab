"""Frozen synthetic comparison of observation-noise dose and rate estimators.

Only the model-informed diagnostic receives a declared model family. Primary
estimators use sampled times and values; generator truth is used only to score
their outputs after fitting.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import shutil
from collections import defaultdict
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
    fit_canonical_trajectory,
    local_recovery_rate,
    stable_equilibrium,
)

AMPLITUDES = (0.0025, 0.04)
NOISE_FRACTIONS = (0.0, 0.0025, 0.01, 0.05)
CALIBRATION_SEEDS = tuple(range(58000, 58008))
EVALUATION_SEEDS = tuple(range(58050, 58100))
TIME_STEP_S = 0.1
TIME_STOP_S = 75.0
TIME_COUNT = 751
BOOTSTRAP_REPLICATES = 2000
ESTIMATOR_NAMES = (
    "free_offset_exponential",
    "tail_offset_log_linear",
    "correct_normal_form",
)
PRIMARY_ESTIMATORS = ESTIMATOR_NAMES[:2]
SIGNATURES = {
    "fold": (0.5, 0.5),
    "transcritical": (1.0, 1.0),
    "supercritical_pitchfork": (0.5, 1.0),
}


def validate_protocol(protocol: dict[str, Any]) -> None:
    """Reject any executable protocol that differs from the frozen design."""
    if (
        protocol.get("phase") != 57
        or protocol.get("protocol_version") != "1.0"
        or protocol.get("protocol_status")
        != "frozen_before_implementation_and_held_out_evaluation"
        or tuple(protocol.get("margins", ())) != TRUE_MARGINS
        or tuple(protocol.get("amplitude_fraction_of_equilibrium", ())) != AMPLITUDES
        or tuple(protocol.get("noise_sd_fraction_of_equilibrium", ()))
        != NOISE_FRACTIONS
        or tuple(protocol.get("calibration_seeds", ())) != CALIBRATION_SEEDS
        or tuple(protocol.get("evaluation_seeds", ())) != EVALUATION_SEEDS
        or protocol.get("bootstrap_replicates") != BOOTSTRAP_REPLICATES
    ):
        raise ScientificValidationError("Phase 57 frozen protocol mismatch")
    grid = protocol.get("time_grid_seconds", {})
    if (
        grid.get("start") != 0.0
        or grid.get("stop") != TIME_STOP_S
        or grid.get("step") != TIME_STEP_S
        or grid.get("count") != TIME_COUNT
        or grid.get("same_for_all_conditions") is not True
    ):
        raise ScientificValidationError("Phase 57 fixed time grid mismatch")
    models = protocol.get("models", {})
    if tuple(models) != MODEL_NAMES or any(
        (models[name].get("q"), models[name].get("p")) != SIGNATURES[name]
        for name in MODEL_NAMES
    ):
        raise ScientificValidationError("Phase 57 normal forms/signatures mismatch")
    seeds = protocol["calibration_seeds"]
    eval_seeds = protocol["evaluation_seeds"]
    if not set(seeds).isdisjoint(eval_seeds):
        raise ScientificValidationError("calibration and evaluation seeds overlap")
    if protocol.get("phase32_status") != "BLOCKED_AUTHORIZATION":
        raise ScientificValidationError(
            "Phase 32 authorization state must be preserved"
        )
    if protocol.get("phase52_status") != "NO-GO":
        raise ScientificValidationError(
            "Phase 52 qualification state must be preserved"
        )


def theoretical_rate(model: str, margin: float) -> float:
    """Return the Phase 53 local recovery rate on its selected stable branch."""
    return local_recovery_rate(model, margin)


def theoretical_exponents(model: str) -> tuple[float, float]:
    """Return the equilibrium and recovery-rate exponents for a named branch."""
    try:
        return SIGNATURES[model]
    except KeyError as error:
        raise ScientificValidationError(f"unknown Phase 57 model: {model}") from error


def design() -> list[dict[str, Any]]:
    """Construct the stable, fully crossed generator conditions."""
    return [
        {
            "model": model,
            "margin": margin,
            "control_setting": CONTROL_BOUNDARY + margin,
            "amplitude_fraction": amplitude,
            "noise_fraction": noise,
        }
        for model in MODEL_NAMES
        for amplitude in AMPLITUDES
        for margin in TRUE_MARGINS
        for noise in NOISE_FRACTIONS
    ]


def _child_seed(seed: int, model: str, margin: float) -> int:
    return int(
        np.random.SeedSequence(
            [seed, MODEL_NAMES.index(model), TRUE_MARGINS.index(margin)]
        ).generate_state(1)[0]
    )


def standardized_innovations(
    seed: int, model: str, margin: float
) -> NDArray[np.float64]:
    """Return paired standard-normal innovations for all doses and amplitudes."""
    if model not in MODEL_NAMES or margin not in TRUE_MARGINS:
        raise ScientificValidationError("invalid key for paired Phase 57 innovations")
    rng = np.random.default_rng(_child_seed(seed, model, margin))
    return rng.standard_normal(TIME_COUNT, dtype=np.float64)


def observation(
    model: str,
    margin: float,
    amplitude: float,
    noise_fraction: float,
    seed: int,
) -> tuple[
    NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], dict[str, float]
]:
    """Generate one paired trace; returned truth metadata is not estimator input."""
    if model not in MODEL_NAMES or margin not in TRUE_MARGINS:
        raise ScientificValidationError("invalid model or margin")
    if amplitude not in AMPLITUDES or noise_fraction not in NOISE_FRACTIONS:
        raise ScientificValidationError("invalid Phase 57 amplitude or noise dose")
    times = np.arange(TIME_COUNT, dtype=np.float64) * TIME_STEP_S
    equilibrium = stable_equilibrium(model, margin)
    rate = theoretical_rate(model, margin)
    initial_state = equilibrium * (1.0 + amplitude)
    clean = exact_trajectory(model, margin, initial_state, times)
    sigma = noise_fraction * equilibrium
    observed = clean + sigma * standardized_innovations(seed, model, margin)
    metadata = {
        "true_equilibrium": equilibrium,
        "true_rate": rate,
        "true_initial_displacement": amplitude * equilibrium,
        "noise_sd": sigma,
        "snr_rms_excursion": (
            math.inf
            if sigma == 0.0
            else float(np.sqrt(np.mean(np.square(clean - equilibrium))) / sigma)
        ),
    }
    return times, clean, observed, metadata


def _observed_scale(observed: NDArray[np.float64]) -> float:
    """A numerical scale based only on samples, never on generator parameters."""
    return max(float(np.ptp(observed)) * 0.001, 1e-12)


def _estimate_free_exponential(
    times: NDArray[np.float64], observed: NDArray[np.float64]
) -> dict[str, Any]:
    fit = estimate_exponential(times, observed, _observed_scale(observed))
    return {
        "estimated_equilibrium": fit.equilibrium,
        "estimated_margin": math.nan,
        "estimated_rate": fit.recovery_rate,
        "rate_standard_error": fit.recovery_rate_standard_error,
        "fit_rmse": fit.fit_rmse,
        "converged": fit.converged,
        "included_points": int(observed.size),
        "failure_reason": "" if fit.converged else "optimizer_not_converged",
    }


def _estimate_tail_log_linear(
    times: NDArray[np.float64], observed: NDArray[np.float64]
) -> dict[str, Any]:
    tail_count = math.ceil(0.10 * observed.size)
    offset = float(np.mean(observed[-tail_count:]))
    residual = observed - offset
    positive = residual > 0.0
    count = int(np.count_nonzero(positive))
    base = {
        "estimated_equilibrium": offset,
        "estimated_margin": math.nan,
        "estimated_rate": math.nan,
        "rate_standard_error": math.nan,
        "fit_rmse": math.nan,
        "converged": False,
        "included_points": count,
        "failure_reason": "",
    }
    if count < 5:
        base["failure_reason"] = "fewer_than_five_positive_residuals"
        return base
    log_residual = np.log(residual[positive])
    design_matrix = np.column_stack((np.ones(count), times[positive]))
    coefficients, _, rank, _ = np.linalg.lstsq(design_matrix, log_residual, rcond=None)
    slope = float(coefficients[1])
    if rank < 2 or not math.isfinite(slope) or slope >= 0.0:
        base["failure_reason"] = "nonnegative_or_undefined_log_slope"
        return base
    residual_fit = log_residual - design_matrix @ coefficients
    rmse = float(np.sqrt(np.mean(np.square(residual_fit))))
    dof = max(1, count - 2)
    covariance = np.linalg.pinv(design_matrix.T @ design_matrix) * float(
        np.dot(residual_fit, residual_fit) / dof
    )
    slope_se = math.sqrt(max(0.0, float(covariance[1, 1])))
    rate = -slope
    if not math.isfinite(rate) or rate <= 0.0:
        base["failure_reason"] = "nonpositive_or_nonfinite_rate"
        return base
    base.update(
        estimated_rate=rate,
        rate_standard_error=slope_se,
        fit_rmse=rmse,
        converged=True,
    )
    return base


def _estimate_correct_form(
    model: str, times: NDArray[np.float64], observed: NDArray[np.float64]
) -> dict[str, Any]:
    fit = fit_canonical_trajectory(model, times, observed, _observed_scale(observed))
    return {
        "estimated_equilibrium": fit.estimated_equilibrium,
        "estimated_margin": fit.estimated_margin,
        "estimated_rate": fit.predicted_recovery_rate,
        "rate_standard_error": math.nan,
        "fit_rmse": math.sqrt(max(0.0, fit.normalized_rss))
        * _observed_scale(observed)
        / math.sqrt(observed.size),
        "converged": fit.converged,
        "included_points": int(observed.size),
        "failure_reason": "optimizer_not_converged" if not fit.converged else "",
    }


def _fit_reason(result: dict[str, Any]) -> str:
    if not result["converged"]:
        return str(result["failure_reason"] or "optimizer_not_converged")
    if not math.isfinite(float(result["estimated_rate"])):
        return "nonfinite_rate"
    if float(result["estimated_rate"]) <= 0.0:
        return "nonpositive_rate"
    if not math.isfinite(float(result["estimated_equilibrium"])):
        return "nonfinite_equilibrium"
    return ""


def _fit_record(
    partition: str,
    seed: int,
    condition: dict[str, Any],
    times: NDArray[np.float64],
    observed: NDArray[np.float64],
    truth: dict[str, float],
) -> list[dict[str, Any]]:
    model = str(condition["model"])
    methods = {
        PRIMARY_ESTIMATORS[0]: _estimate_free_exponential(times, observed),
        PRIMARY_ESTIMATORS[1]: _estimate_tail_log_linear(times, observed),
        "correct_normal_form": _estimate_correct_form(model, times, observed),
    }
    rows: list[dict[str, Any]] = []
    for estimator, result in methods.items():
        reason = _fit_reason(result)
        valid = not reason
        estimated_rate = float(result["estimated_rate"])
        relative_error = (
            (estimated_rate - truth["true_rate"]) / truth["true_rate"]
            if valid
            else math.nan
        )
        equilibrium = float(result["estimated_equilibrium"])
        rows.append(
            {
                "partition": partition,
                "seed": seed,
                **condition,
                "estimator": estimator,
                "true_equilibrium": truth["true_equilibrium"],
                "true_rate": truth["true_rate"],
                "true_margin": condition["margin"],
                "estimated_equilibrium": equilibrium,
                "estimated_margin": result["estimated_margin"],
                "equilibrium_relative_error": abs(
                    equilibrium - truth["true_equilibrium"]
                )
                / truth["true_equilibrium"],
                "estimated_rate": estimated_rate,
                "rate_standard_error": result["rate_standard_error"],
                "rate_signed_error": estimated_rate - truth["true_rate"]
                if valid
                else math.nan,
                "rate_relative_signed_error": relative_error,
                "rate_absolute_error": abs(estimated_rate - truth["true_rate"])
                if valid
                else math.nan,
                "rate_relative_absolute_error": abs(relative_error)
                if valid
                else math.nan,
                "fit_rmse": result["fit_rmse"],
                "included_points": result["included_points"],
                "valid": valid,
                "failure_reason": reason,
                "snr_rms_excursion": truth["snr_rms_excursion"],
                "noise_sd_state_units": truth["noise_sd"],
            }
        )
    return rows


def _seed_bootstrap_interval(
    values: list[float], label: str
) -> tuple[float | None, float | None]:
    finite = np.asarray([value for value in values if math.isfinite(value)])
    if finite.size < 2:
        return None, None
    seed = int.from_bytes(hashlib.sha256(label.encode()).digest()[:4], "big")
    rng = np.random.default_rng(seed)
    samples = rng.choice(finite, size=(BOOTSTRAP_REPLICATES, finite.size), replace=True)
    means = np.mean(samples, axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(low), float(high)


def _interval_excludes_zero(interval: tuple[float | None, float | None]) -> bool:
    """Return whether a complete interval is strictly on one side of zero."""
    low, high = interval
    return low is not None and high is not None and (low > 0.0 or high < 0.0)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ScientificValidationError(f"refusing to write empty table: {path.name}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _rate_summaries(
    records: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    eval_rows = [row for row in records if row["partition"] == "evaluation"]
    keys = ("model", "amplitude_fraction", "noise_fraction", "estimator")
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    per_margin: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in eval_rows:
        grouped[tuple(row[key] for key in keys)].append(row)
        per_margin[
            tuple(row[key] for key in (*keys[:3], "margin", "estimator"))
        ].append(row)

    summaries: list[dict[str, Any]] = []
    for key, rows in grouped.items():
        valid = [row for row in rows if row["valid"]]
        by_seed: dict[int, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            by_seed[int(row["seed"])].append(row)
        seed_bias = [
            float(
                np.mean([r["rate_relative_signed_error"] for r in values if r["valid"]])
            )
            if any(r["valid"] for r in values)
            else math.nan
            for values in by_seed.values()
        ]
        seed_mae = [
            float(
                np.mean(
                    [r["rate_relative_absolute_error"] for r in values if r["valid"]]
                )
            )
            if any(r["valid"] for r in values)
            else math.nan
            for values in by_seed.values()
        ]
        seed_validity = [
            float(np.mean([r["valid"] for r in values])) for values in by_seed.values()
        ]
        label = "|".join(map(str, key))
        bias_ci = (
            _seed_bootstrap_interval(seed_bias, "bias|" + label)
            if key[2]
            else (None, None)
        )
        mae_ci = (
            _seed_bootstrap_interval(seed_mae, "mae|" + label)
            if key[2]
            else (None, None)
        )
        valid_ci = (
            _seed_bootstrap_interval(seed_validity, "valid|" + label)
            if key[2]
            else (None, None)
        )
        abs_errors = [r["rate_relative_absolute_error"] for r in valid]
        summaries.append(
            {
                **dict(zip(keys, key, strict=True)),
                "attempted_fits": len(rows),
                "valid_fits": len(valid),
                "valid_proportion": len(valid) / len(rows),
                "valid_proportion_ci_low": valid_ci[0],
                "valid_proportion_ci_high": valid_ci[1],
                "mean_relative_signed_error": float(
                    np.mean([r["rate_relative_signed_error"] for r in valid])
                )
                if valid
                else math.nan,
                "mean_relative_signed_error_ci_low": bias_ci[0],
                "mean_relative_signed_error_ci_high": bias_ci[1],
                "mean_relative_absolute_error": float(np.mean(abs_errors))
                if valid
                else math.nan,
                "mean_relative_absolute_error_ci_low": mae_ci[0],
                "mean_relative_absolute_error_ci_high": mae_ci[1],
                "median_relative_absolute_error": float(np.median(abs_errors))
                if valid
                else math.nan,
                "p90_relative_absolute_error": float(np.quantile(abs_errors, 0.9))
                if valid
                else math.nan,
                "mean_equilibrium_relative_error": float(
                    np.mean([r["equilibrium_relative_error"] for r in rows])
                ),
                "failure_reasons": json.dumps(
                    {
                        reason: sum(r["failure_reason"] == reason for r in rows)
                        for reason in sorted(
                            {r["failure_reason"] for r in rows if r["failure_reason"]}
                        )
                    },
                    sort_keys=True,
                ),
                "zero_noise_interval_note": (
                    "no interval; seed-labelled observations repeat exactly"
                )
                if key[2] == 0.0
                else "seed-cluster percentile bootstrap",
            }
        )

    margin_summaries: list[dict[str, Any]] = []
    margin_keys = (
        "model",
        "amplitude_fraction",
        "noise_fraction",
        "margin",
        "estimator",
    )
    for key, rows in per_margin.items():
        valid = [r for r in rows if r["valid"]]
        errors = [r["rate_relative_absolute_error"] for r in valid]
        bias = [r["rate_relative_signed_error"] for r in valid]
        margin_summaries.append(
            {
                **dict(zip(margin_keys, key, strict=True)),
                "attempted_fits": len(rows),
                "valid_fits": len(valid),
                "valid_proportion": len(valid) / len(rows),
                "mean_relative_signed_error": float(np.mean(bias))
                if valid
                else math.nan,
                "mean_relative_absolute_error": float(np.mean(errors))
                if valid
                else math.nan,
                "median_relative_absolute_error": float(np.median(errors))
                if valid
                else math.nan,
                "p90_relative_absolute_error": float(np.quantile(errors, 0.9))
                if valid
                else math.nan,
                "mean_equilibrium_relative_error": float(
                    np.mean([r["equilibrium_relative_error"] for r in rows])
                ),
            }
        )
    return summaries, margin_summaries


def _exponent_groups(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_group: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in records:
        if (
            row["partition"] != "evaluation"
            or row["estimator"] not in PRIMARY_ESTIMATORS
        ):
            continue
        key = (
            row["seed"],
            row["model"],
            row["amplitude_fraction"],
            row["noise_fraction"],
            row["estimator"],
        )
        by_group[key].append(row)
    results: list[dict[str, Any]] = []
    for (seed, model, amplitude, noise, estimator), rows in by_group.items():
        rows.sort(key=lambda row: row["margin"])
        reasons: list[str] = []
        if len(rows) != len(TRUE_MARGINS):
            reasons.append("incomplete_margin_group")
        if any(not row["valid"] for row in rows):
            reasons.append("one_or_more_invalid_rates")
        if any(
            not math.isfinite(row["estimated_equilibrium"])
            or row["estimated_equilibrium"] <= 0.0
            for row in rows
        ):
            reasons.append("nonpositive_or_nonfinite_equilibrium")
        fit = None
        if not reasons:
            settings = np.asarray(
                [r["control_setting"] for r in rows], dtype=np.float64
            )
            equilibria = np.asarray(
                [r["estimated_equilibrium"] for r in rows], dtype=np.float64
            )
            rates = np.asarray([r["estimated_rate"] for r in rows], dtype=np.float64)
            try:
                fit = estimate_scaling(settings, equilibria, rates)
            except (ScientificValidationError, ValueError, FloatingPointError):
                reasons.append("scaling_fit_exception")
        if fit is not None and not fit.converged:
            reasons.append("scaling_optimizer_not_converged")
        valid = fit is not None and fit.converged and not reasons
        q_true, p_true = SIGNATURES[str(model)]
        results.append(
            {
                "seed": seed,
                "model": model,
                "amplitude_fraction": amplitude,
                "noise_fraction": noise,
                "estimator": estimator,
                "attempted_margin_count": len(rows),
                "rate_valid_fraction": sum(bool(r["valid"]) for r in rows) / len(rows)
                if rows
                else 0.0,
                "valid": bool(valid),
                "failure_reason": ";".join(reasons),
                "estimated_boundary": fit.estimated_boundary
                if valid and fit
                else math.nan,
                "estimated_q": fit.equilibrium_exponent if valid and fit else math.nan,
                "true_q": q_true,
                "q_bias": fit.equilibrium_exponent - q_true
                if valid and fit
                else math.nan,
                "q_absolute_error": abs(fit.equilibrium_exponent - q_true)
                if valid and fit
                else math.nan,
                "estimated_p": fit.rate_exponent if valid and fit else math.nan,
                "true_p": p_true,
                "p_bias": fit.rate_exponent - p_true if valid and fit else math.nan,
                "p_absolute_error": abs(fit.rate_exponent - p_true)
                if valid and fit
                else math.nan,
            }
        )
    return results


def _exponent_summaries(groups: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = ("model", "amplitude_fraction", "noise_fraction", "estimator")
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in groups:
        grouped[tuple(row[key] for key in keys)].append(row)
    summaries: list[dict[str, Any]] = []
    for key, rows in grouped.items():
        valid = [r for r in rows if r["valid"]]
        q_true, p_true = SIGNATURES[str(key[0])]
        q_values = [float(r["estimated_q"]) for r in valid]
        p_values = [float(r["estimated_p"]) for r in valid]
        label = "|".join(map(str, key))
        q_ci = (
            _seed_bootstrap_interval(q_values, "q|" + label) if key[2] else (None, None)
        )
        p_ci = (
            _seed_bootstrap_interval(p_values, "p|" + label) if key[2] else (None, None)
        )
        p_cover: bool | None = (
            None
            if key[2] == 0.0
            else bool(
                p_ci[0] is not None
                and p_ci[1] is not None
                and p_ci[0] <= p_true <= p_ci[1]
            )
        )
        mae_p = (
            float(np.mean([r["p_absolute_error"] for r in valid]))
            if valid
            else math.nan
        )
        valid_fraction = len(valid) / len(rows)
        identifiable: bool | None = (
            None
            if key[2] == 0.0
            else bool(valid_fraction >= 0.8 and mae_p <= 0.20 and p_cover)
        )
        summaries.append(
            {
                **dict(zip(keys, key, strict=True)),
                "attempted_groups": len(rows),
                "valid_groups": len(valid),
                "valid_group_fraction": valid_fraction,
                "mean_q": float(np.mean(q_values)) if valid else math.nan,
                "q_bias": float(np.mean([r["q_bias"] for r in valid]))
                if valid
                else math.nan,
                "q_mean_absolute_error": float(
                    np.mean([r["q_absolute_error"] for r in valid])
                )
                if valid
                else math.nan,
                "q_ci_low": q_ci[0],
                "q_ci_high": q_ci[1],
                "mean_p": float(np.mean(p_values)) if valid else math.nan,
                "p_bias": float(np.mean([r["p_bias"] for r in valid]))
                if valid
                else math.nan,
                "p_mean_absolute_error": mae_p,
                "p_ci_low": p_ci[0],
                "p_ci_high": p_ci[1],
                "p_ci_width": p_ci[1] - p_ci[0]
                if p_ci[0] is not None and p_ci[1] is not None
                else math.nan,
                "p_interval_covers_theory": p_cover,
                "p_identifiability_criterion_met": identifiable,
                "p_identifiability_criterion_status": (
                    "not_assessed_zero_noise_diagnostic"
                    if key[2] == 0.0
                    else "met"
                    if identifiable
                    else "not_met"
                ),
                "failure_reasons": json.dumps(
                    {
                        reason: sum(
                            reason in str(r["failure_reason"]).split(";") for r in rows
                        )
                        for reason in sorted(
                            {
                                part
                                for r in rows
                                for part in str(r["failure_reason"]).split(";")
                                if part
                            }
                        )
                    },
                    sort_keys=True,
                ),
                "zero_noise_interval_note": (
                    "no interval; seed-labelled trajectories repeat exactly"
                )
                if key[2] == 0.0
                else "seed-cluster percentile bootstrap",
            }
        )
    return summaries


def _paired_contrasts(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    evaluation = [r for r in records if r["partition"] == "evaluation"]
    index = {
        (
            r["seed"],
            r["model"],
            r["amplitude_fraction"],
            r["margin"],
            r["noise_fraction"],
            r["estimator"],
        ): r
        for r in evaluation
    }
    output: list[dict[str, Any]] = []
    for model in MODEL_NAMES:
        for amplitude in AMPLITUDES:
            for dose in NOISE_FRACTIONS[1:]:
                for estimand in (
                    "free_offset_exponential",
                    "tail_offset_log_linear",
                    "estimator_interaction",
                ):
                    seed_diffs: list[float] = []
                    matched_margin_count = 0
                    for seed in EVALUATION_SEEDS:
                        diffs: list[float] = []
                        for margin in TRUE_MARGINS:
                            if estimand == "estimator_interaction":
                                rows = [
                                    index.get((seed, model, amplitude, margin, d, est))
                                    for d in (dose, 0.0)
                                    for est in PRIMARY_ESTIMATORS
                                ]
                                present = [row for row in rows if row is not None]
                                if len(present) != 4 or any(
                                    not row["valid"] for row in present
                                ):
                                    continue
                                dose_gap = (
                                    present[1]["rate_relative_absolute_error"]
                                    - present[0]["rate_relative_absolute_error"]
                                )
                                zero_gap = (
                                    present[3]["rate_relative_absolute_error"]
                                    - present[2]["rate_relative_absolute_error"]
                                )
                                diffs.append(float(dose_gap - zero_gap))
                            else:
                                high = index.get(
                                    (seed, model, amplitude, margin, dose, estimand)
                                )
                                zero = index.get(
                                    (seed, model, amplitude, margin, 0.0, estimand)
                                )
                                if (
                                    high is None
                                    or zero is None
                                    or not high["valid"]
                                    or not zero["valid"]
                                ):
                                    continue
                                diffs.append(
                                    float(
                                        high["rate_relative_absolute_error"]
                                        - zero["rate_relative_absolute_error"]
                                    )
                                )
                        if diffs:
                            matched_margin_count += len(diffs)
                            seed_diffs.append(float(np.mean(diffs)))
                    label = f"{model}|{amplitude}|{dose}|{estimand}"
                    interval = _seed_bootstrap_interval(seed_diffs, "contrast|" + label)
                    output.append(
                        {
                            "model": model,
                            "amplitude_fraction": amplitude,
                            "noise_fraction": dose,
                            "contrast": "dose_by_estimator_interaction"
                            if estimand == "estimator_interaction"
                            else "dose_minus_zero_noise",
                            "estimator": estimand,
                            "paired_seed_count": len(seed_diffs),
                            "matched_margin_count": matched_margin_count,
                            "mean_paired_difference": float(np.mean(seed_diffs))
                            if seed_diffs
                            else math.nan,
                            "ci_low": interval[0],
                            "ci_high": interval[1],
                            "interval_excludes_zero": _interval_excludes_zero(interval),
                            "interpretation": "supported in observed direction"
                            if _interval_excludes_zero(interval)
                            else "inconclusive/not detected",
                        }
                    )
    # Primary estimator-form contrast paired within seed and margin.
    for model in MODEL_NAMES:
        for amplitude in AMPLITUDES:
            for dose in NOISE_FRACTIONS:
                seed_diffs = []
                matched_margin_count = 0
                for seed in EVALUATION_SEEDS:
                    diffs = []
                    for margin in TRUE_MARGINS:
                        ref = index.get(
                            (
                                seed,
                                model,
                                amplitude,
                                margin,
                                dose,
                                PRIMARY_ESTIMATORS[0],
                            )
                        )
                        alt = index.get(
                            (
                                seed,
                                model,
                                amplitude,
                                margin,
                                dose,
                                PRIMARY_ESTIMATORS[1],
                            )
                        )
                        if (
                            ref is not None
                            and alt is not None
                            and ref["valid"]
                            and alt["valid"]
                        ):
                            diffs.append(
                                float(
                                    alt["rate_relative_absolute_error"]
                                    - ref["rate_relative_absolute_error"]
                                )
                            )
                    if diffs:
                        matched_margin_count += len(diffs)
                        seed_diffs.append(float(np.mean(diffs)))
                interval = (
                    _seed_bootstrap_interval(
                        seed_diffs, f"estimator|{model}|{amplitude}|{dose}"
                    )
                    if dose
                    else (None, None)
                )
                output.append(
                    {
                        "model": model,
                        "amplitude_fraction": amplitude,
                        "noise_fraction": dose,
                        "contrast": "tail_minus_free_offset_absolute_relative_error",
                        "estimator": (
                            "tail_offset_log_linear-minus-free_offset_exponential"
                        ),
                        "paired_seed_count": len(seed_diffs),
                        "matched_margin_count": matched_margin_count,
                        "mean_paired_difference": float(np.mean(seed_diffs))
                        if seed_diffs
                        else math.nan,
                        "ci_low": interval[0],
                        "ci_high": interval[1],
                        "interval_excludes_zero": _interval_excludes_zero(interval),
                        "interpretation": "deterministic point contrast; no interval"
                        if dose == 0
                        else (
                            "supported in observed direction"
                            if _interval_excludes_zero(interval)
                            else "inconclusive/not detected"
                        ),
                    }
                )
    return output


def _manifest(output_dir: Path, names: list[str]) -> dict[str, str]:
    return {
        name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
        for name in names
    }


def _snr_summary(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Summarize the deterministic generator-only SNR over five margins."""
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in records:
        if row["partition"] == "evaluation":
            key = (row["model"], row["amplitude_fraction"], row["noise_fraction"])
            grouped[key].append(row)
    result: list[dict[str, Any]] = []
    for (model, amplitude, dose), rows in grouped.items():
        by_margin = {
            float(row["margin"]): float(row["snr_rms_excursion"]) for row in rows
        }
        values = list(by_margin.values())
        result.append(
            {
                "model": model,
                "amplitude_fraction": amplitude,
                "noise_fraction": dose,
                "margin_count": len(values),
                "snr_min_across_margins": min(values),
                "snr_max_across_margins": max(values),
                "snr_mean_across_margins": float(np.mean(values)),
                "interpretation": (
                    "generator-only descriptive diagnostic; no interval or causal claim"
                ),
            }
        )
    return result


def run_study(output_dir: Path) -> dict[str, Any]:
    """Run the frozen calibration/evaluation design and write auditable outputs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    protocol_path = output_dir / "protocol.json"
    try:
        protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificValidationError(
            f"cannot read frozen Phase 57 protocol: {error}"
        ) from error
    validate_protocol(protocol)
    protocol_hash = hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    times = np.arange(TIME_COUNT, dtype=np.float64) * TIME_STEP_S
    records: list[dict[str, Any]] = []
    for partition, seeds in (
        ("calibration", CALIBRATION_SEEDS),
        ("evaluation", EVALUATION_SEEDS),
    ):
        for seed in seeds:
            for condition in design():
                _, _, observed, truth = observation(
                    str(condition["model"]),
                    float(condition["margin"]),
                    float(condition["amplitude_fraction"]),
                    float(condition["noise_fraction"]),
                    seed,
                )
                records.extend(
                    _fit_record(partition, seed, condition, times, observed, truth)
                )

    rate_summary, margin_summary = _rate_summaries(records)
    exponent_groups = _exponent_groups(records)
    exponent_summary = _exponent_summaries(exponent_groups)
    contrasts = _paired_contrasts(records)
    snr_rows = _snr_summary(records)
    _write_csv(output_dir / "per_fit_results.csv", records)
    _write_csv(output_dir / "rate_condition_summary.csv", rate_summary)
    _write_csv(output_dir / "rate_margin_summary.csv", margin_summary)
    _write_csv(output_dir / "exponent_group_results.csv", exponent_groups)
    _write_csv(output_dir / "exponent_condition_summary.csv", exponent_summary)
    _write_csv(output_dir / "paired_contrasts.csv", contrasts)
    _write_csv(output_dir / "snr_summary.csv", snr_rows)
    metadata = {
        "phase": 57,
        "protocol_version": protocol["protocol_version"],
        "protocol_status": protocol["protocol_status"],
        "protocol_sha256": protocol_hash,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "time_count": TIME_COUNT,
        "trace_count_total": len(records) // len(ESTIMATOR_NAMES),
        "fit_attempt_count_total": len(records),
        "calibration_trace_count": len(CALIBRATION_SEEDS) * 120,
        "evaluation_trace_count": len(EVALUATION_SEEDS) * 120,
        "calibration_fit_attempt_count": sum(
            r["partition"] == "calibration" for r in records
        ),
        "evaluation_fit_attempt_count": sum(
            r["partition"] == "evaluation" for r in records
        ),
        "evaluation_exponent_group_count": len(exponent_groups),
        "seed_partition_disjoint": set(CALIBRATION_SEEDS).isdisjoint(EVALUATION_SEEDS),
        "estimators": list(ESTIMATOR_NAMES),
        "primary_estimators": list(PRIMARY_ESTIMATORS),
        "true_parameters_exposed_to_estimators": False,
        "model_info_exposed_only_to_diagnostic": True,
        "zero_noise_replicates_are_identical_by_seed": True,
        "phase32_status": "BLOCKED_AUTHORIZATION",
        "phase52_status": "NO-GO",
    }
    (output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest_names = [
        "protocol.json",
        "per_fit_results.csv",
        "rate_condition_summary.csv",
        "rate_margin_summary.csv",
        "exponent_group_results.csv",
        "exponent_condition_summary.csv",
        "paired_contrasts.csv",
        "snr_summary.csv",
        "metadata.json",
    ]
    manifest = _manifest(output_dir, manifest_names)
    (output_dir / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata


def verify_manifest(output_dir: Path) -> bool:
    """Return whether every listed output still matches its recorded SHA-256."""
    try:
        manifest = json.loads(
            (output_dir / "sha256_manifest.json").read_text(encoding="utf-8")
        )
        return bool(manifest) and all(
            hashlib.sha256((output_dir / name).read_bytes()).hexdigest() == digest
            for name, digest in manifest.items()
        )
    except (OSError, json.JSONDecodeError, TypeError):
        return False


def reproduce_to(output_dir: Path, repeat_dir: Path) -> bool:
    """Run in an isolated destination and compare every prescribed file hash."""
    repeat_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(output_dir / "protocol.json", repeat_dir / "protocol.json")
    run_study(repeat_dir)
    original = json.loads(
        (output_dir / "sha256_manifest.json").read_text(encoding="utf-8")
    )
    repeated = json.loads(
        (repeat_dir / "sha256_manifest.json").read_text(encoding="utf-8")
    )
    return bool(original == repeated)
