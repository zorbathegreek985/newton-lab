"""Controlled Phase 54 ablation of synthetic exponent identifiability."""

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
from scipy.optimize import least_squares  # type: ignore[import-untyped]

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

FACTOR_LEVELS: dict[str, tuple[Any, ...]] = {
    "amplitude_fraction": (0.05, 0.20, 0.50),
    "information_condition": (
        "reference",
        "oracle_equilibrium",
        "oracle_boundary",
        "oracle_both",
    ),
    "noise_fraction_of_equilibrium": (0.0, 0.01, 0.05),
    "window_tau": (1.0, 3.0),
}
SIGNATURES: dict[str, tuple[float, float]] = {
    "fold": (0.5, 0.5),
    "transcritical": (1.0, 1.0),
    "supercritical_pitchfork": (0.5, 1.0),
}
BOOTSTRAP_REPLICATES = 500
SEED_STARTS = {"calibration": 54000, "evaluation": 55000}
SEED_COUNTS = {"calibration": 8, "evaluation": 20}
INFORMATION_MODES = FACTOR_LEVELS["information_condition"]
_FACTOR_NAMES = tuple(FACTOR_LEVELS)


def factor_grid() -> list[dict[str, Any]]:
    """Return the frozen balanced 3x4x3x2 factorial in stable order."""
    return [
        dict(zip(_FACTOR_NAMES, values, strict=True))
        for values in itertools.product(
            *(FACTOR_LEVELS[name] for name in _FACTOR_NAMES)
        )
    ]


def validate_protocol(protocol: dict[str, Any]) -> None:
    """Reject protocol edits that would silently disagree with this implementation."""
    if protocol.get("protocol_status") != "frozen_before_implementation_and_evaluation":
        raise ScientificValidationError("Phase 54 protocol must be frozen")
    if protocol.get("protocol_version") != "1.0":
        raise ScientificValidationError("unsupported Phase 54 protocol version")
    factors = protocol.get("factors", {})
    expected = {
        "amplitude_fraction_of_equilibrium": list(FACTOR_LEVELS["amplitude_fraction"]),
        "information_condition": list(FACTOR_LEVELS["information_condition"]),
        "noise_sd_fraction_of_equilibrium": list(
            FACTOR_LEVELS["noise_fraction_of_equilibrium"]
        ),
        "window_multiples_of_local_recovery_time": list(FACTOR_LEVELS["window_tau"]),
    }
    if any(factors.get(name) != levels for name, levels in expected.items()):
        raise ScientificValidationError(
            "protocol factor levels do not match implementation"
        )
    design = protocol.get("design", {})
    calibration = design.get("calibration_seeds", [])
    evaluation = design.get("evaluation_seeds", [])
    if (
        len(calibration) != 2
        or len(evaluation) != 2
        or calibration[0] > calibration[1]
        or evaluation[0] > evaluation[1]
        or not set(range(calibration[0], calibration[1] + 1)).isdisjoint(
            range(evaluation[0], evaluation[1] + 1)
        )
    ):
        raise ScientificValidationError(
            "calibration/evaluation seeds must be disjoint ranges"
        )
    if calibration != [54000, 54007] or evaluation != [55000, 55019]:
        raise ScientificValidationError(
            "protocol seed ranges do not match implementation"
        )
    if design.get("factorial_cells_per_model") != 72:
        raise ScientificValidationError(
            "protocol factorial size does not match implementation"
        )
    controls = protocol.get("control_design", {})
    if (
        controls.get("true_boundary") != CONTROL_BOUNDARY
        or controls.get("margins_generator_only") != list(TRUE_MARGINS)
        or len(controls.get("controls", [])) != len(TRUE_MARGINS)
        or any(
            not math.isclose(actual, CONTROL_BOUNDARY + margin, abs_tol=1e-12)
            for actual, margin in zip(
                controls.get("controls", []), TRUE_MARGINS, strict=True
            )
        )
    ):
        raise ScientificValidationError("protocol controls do not match Phase 53")
    models = protocol.get("models", {})
    if tuple(models) != MODEL_NAMES or any(
        (models[name].get("q"), models[name].get("p")) != SIGNATURES[name]
        for name in MODEL_NAMES
    ):
        raise ScientificValidationError("protocol signatures do not match Phase 53")


def paired_noise_innovations(
    seed: int, model: str, margin: float
) -> NDArray[np.float64]:
    """Return the common standardized noise stream used across paired factor cells."""
    if model not in MODEL_NAMES or margin not in TRUE_MARGINS:
        raise ScientificValidationError(
            "noise pairing requires a declared model and margin"
        )
    model_index = MODEL_NAMES.index(model)
    margin_index = TRUE_MARGINS.index(margin)
    child_seed = int(
        np.random.SeedSequence([seed, model_index, margin_index]).generate_state(1)[0]
    )
    innovations = np.random.default_rng(child_seed).standard_normal(
        31, dtype=np.float64
    )
    return innovations


def _fixed_equilibrium_exponential(
    times: NDArray[np.float64],
    observed: NDArray[np.float64],
    equilibrium: float,
    noise_scale: float,
) -> tuple[float, float, bool, float]:
    """Estimate positive amplitude/rate with a supplied offset and no true rate."""
    amplitude0 = max(float(observed[0] - equilibrium), noise_scale, 1e-12)
    initial = np.array([math.log(amplitude0), math.log(1.0 / max(times[-1], 1e-12))])

    def residual(parameters: NDArray[np.float64]) -> NDArray[np.float64]:
        amplitude, rate = np.exp(parameters)
        return np.asarray(
            amplitude * np.exp(-rate * times) + equilibrium - observed,
            dtype=np.float64,
        )

    result = least_squares(
        residual,
        initial,
        bounds=(np.log([1e-12, 1e-8]), np.log([10.0, 1e8])),
        max_nfev=2000,
    )
    amplitude, rate = np.exp(result.x)
    errors = residual(result.x)
    rmse = float(np.sqrt(np.mean(np.square(errors))))
    return float(amplitude), float(rate), bool(result.success and rate > 0), rmse


def _ols_power(
    margins: NDArray[np.float64], values: NDArray[np.float64]
) -> tuple[float, float]:
    """Fit log-prefactor and exponent against known positive margins."""
    if (
        margins.ndim != 1
        or values.shape != margins.shape
        or margins.size < 4
        or not np.isfinite(margins).all()
        or not np.isfinite(values).all()
        or np.any(margins <= 0)
        or np.any(values <= 0)
    ):
        raise ScientificValidationError(
            "known-boundary power fit requires positive data"
        )
    slope, intercept = np.polyfit(np.log(margins), np.log(values), 1)
    return float(slope), float(intercept)


def _estimate_rate(
    model: str,
    margin: float,
    amplitude_fraction: float,
    noise_fraction: float,
    window_tau: float,
    information: str,
    seed: int,
) -> dict[str, Any]:
    """Generate one paired observation and estimate its local recovery rate."""
    times, clean, observed = _paired_observation(
        model,
        margin,
        amplitude_fraction,
        noise_fraction,
        window_tau,
        seed,
    )
    equilibrium = stable_equilibrium(model, margin)
    true_rate = local_recovery_rate(model, margin)
    tau = 1.0 / true_rate
    initial = equilibrium * (1.0 + amplitude_fraction)
    noise_sd = noise_fraction * equilibrium
    scale_for_estimator = max(noise_sd, 1e-12)

    if information in ("oracle_equilibrium", "oracle_both"):
        amplitude, rate, converged, rmse = _fixed_equilibrium_exponential(
            times, observed, equilibrium, scale_for_estimator
        )
        estimated_equilibrium = equilibrium
    else:
        estimate = estimate_exponential(times, observed, scale_for_estimator)
        amplitude = estimate.amplitude
        rate = estimate.recovery_rate
        converged = estimate.converged
        rmse = estimate.fit_rmse
        estimated_equilibrium = estimate.equilibrium
    valid = bool(converged and math.isfinite(rate) and rate > 0.0)
    relative_error = abs(rate - true_rate) / true_rate if valid else math.nan
    delta = initial - equilibrium
    return {
        "estimated_equilibrium": estimated_equilibrium,
        "equilibrium_relative_error": abs(estimated_equilibrium - equilibrium)
        / equilibrium,
        "estimated_rate": rate,
        "amplitude_estimate": amplitude,
        "rate_rmse": rmse,
        "rate_valid": valid,
        "true_rate": true_rate,
        "rate_relative_error": relative_error,
        "true_equilibrium": equilibrium,
        "initial_state": initial,
        "perturbation_displacement": delta,
        "perturbation_over_equilibrium": delta / equilibrium,
        "perturbation_over_margin": delta / margin,
        "tau_design_only": tau,
        "duration": window_tau * tau,
        "sample_count": times.size,
        "noise_sd_state": noise_sd,
    }


def _paired_observation(
    model: str,
    margin: float,
    amplitude_fraction: float,
    noise_fraction: float,
    window_tau: float,
    seed: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Construct clean and noisy paired traces for the declared design cell."""
    equilibrium = stable_equilibrium(model, margin)
    tau = 1.0 / local_recovery_rate(model, margin)
    count = int(round(window_tau / 0.1)) + 1
    times = np.arange(count, dtype=np.float64) * (0.1 * tau)
    clean = exact_trajectory(
        model, margin, equilibrium * (1.0 + amplitude_fraction), times
    )
    innovations = paired_noise_innovations(seed, model, margin)[:count]
    observed = clean + noise_fraction * equilibrium * innovations
    return times, clean, observed


def _fit_group(
    model: str,
    information: str,
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    """Estimate q and p for one seed/model/factor combination."""
    ordered = sorted(records, key=lambda row: row["control_setting"])
    margins = np.asarray([row["margin"] for row in ordered], dtype=np.float64)
    settings = np.asarray([row["control_setting"] for row in ordered], dtype=np.float64)
    true_q, true_p = SIGNATURES[model]
    rates_valid = all(row["rate_valid"] for row in ordered)
    rate_valid_fraction = sum(row["rate_valid"] for row in ordered) / len(ordered)
    rate_errors = [row["rate_relative_error"] for row in ordered if row["rate_valid"]]
    rate_error = float(np.mean(rate_errors)) if rate_errors else math.nan
    try:
        equilibria = np.asarray(
            [row["estimated_equilibrium"] for row in ordered], dtype=np.float64
        )
        if information in ("oracle_equilibrium", "oracle_both"):
            equilibria = np.asarray(
                [row["true_equilibrium"] for row in ordered], dtype=np.float64
            )
        rates = np.asarray([row["estimated_rate"] for row in ordered], dtype=np.float64)
        if not rates_valid:
            raise ScientificValidationError("one or more trajectory rate fits failed")
        if information in ("oracle_boundary", "oracle_both"):
            q, _ = _ols_power(margins, equilibria)
            p, _ = _ols_power(margins, rates)
            boundary = CONTROL_BOUNDARY
            converged = True
        else:
            fit = estimate_scaling(settings, equilibria, rates)
            q, p, boundary, converged = (
                fit.equilibrium_exponent,
                fit.rate_exponent,
                fit.estimated_boundary,
                fit.converged,
            )
        valid = bool(converged and math.isfinite(q) and math.isfinite(p))
    except (ScientificValidationError, ValueError, FloatingPointError):
        q = p = boundary = math.nan
        valid = False
    return {
        "valid": valid,
        "rate_valid_fraction": rate_valid_fraction,
        "rate_mean_relative_error": rate_error,
        "q": q,
        "p": p,
        "q_abs_error": abs(q - true_q) if valid else math.nan,
        "p_abs_error": abs(p - true_p) if valid else math.nan,
        "true_q": true_q,
        "true_p": true_p,
        "estimated_boundary": boundary,
        "mean_estimated_boundary": boundary,
        "attempted_margin_count": len(ordered),
    }


def run_ablation(output_dir: Path) -> dict[str, Any]:
    """Run the frozen Phase 54 design and write raw, summary, and hash outputs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    protocol_path = output_dir / "protocol.json"
    try:
        protocol = json.loads(protocol_path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScientificValidationError(
            f"cannot read Phase 54 protocol: {error}"
        ) from error
    validate_protocol(protocol)
    protocol_hash = hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    grid = factor_grid()
    raw: list[dict[str, Any]] = []
    groups: list[dict[str, Any]] = []
    model_ids = {name: index for index, name in enumerate(MODEL_NAMES)}
    margin_ids = {margin: index for index, margin in enumerate(TRUE_MARGINS)}

    for partition, count in SEED_COUNTS.items():
        seed_start = SEED_STARTS[partition]
        for seed in range(seed_start, seed_start + count):
            for model in MODEL_NAMES:
                for factors in grid:
                    block: list[dict[str, Any]] = []
                    for margin in TRUE_MARGINS:
                        noise_seed = int(
                            np.random.SeedSequence(
                                [seed, model_ids[model], margin_ids[margin]]
                            ).generate_state(1)[0]
                        )
                        result = _estimate_rate(
                            model,
                            margin,
                            factors["amplitude_fraction"],
                            factors["noise_fraction_of_equilibrium"],
                            factors["window_tau"],
                            factors["information_condition"],
                            noise_seed,
                        )
                        row = {
                            "partition": partition,
                            "replicate_seed": seed,
                            "model": model,
                            **factors,
                            "margin": margin,
                            "control_setting": CONTROL_BOUNDARY + margin,
                            **result,
                        }
                        raw.append(row)
                        block.append(row)
                    fitted = _fit_group(model, factors["information_condition"], block)
                    groups.append(
                        {
                            "partition": partition,
                            "replicate_seed": seed,
                            "model": model,
                            **factors,
                            **fitted,
                        }
                    )

    evaluation_groups = [row for row in groups if row["partition"] == "evaluation"]
    evaluation_raw = [row for row in raw if row["partition"] == "evaluation"]
    summaries = _cell_summaries(evaluation_groups)
    rate_summaries = _rate_summaries(evaluation_raw)
    contrasts = _factor_contrasts(evaluation_groups)
    _write_csv(output_dir / "per_trajectory_results.csv", raw)
    _write_csv(output_dir / "per_group_results.csv", groups)
    _write_csv(output_dir / "exponent_condition_summary.csv", summaries)
    _write_csv(output_dir / "rate_condition_summary.csv", rate_summaries)
    _write_csv(output_dir / "factor_contrasts.csv", contrasts)
    metadata = {
        "phase": 54,
        "protocol_version": protocol["protocol_version"],
        "protocol_sha256": protocol_hash,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "factor_cells_per_model": len(grid),
        "calibration_seed_count": SEED_COUNTS["calibration"],
        "evaluation_seed_count": SEED_COUNTS["evaluation"],
        "calibration_group_count": sum(
            row["partition"] == "calibration" for row in groups
        ),
        "evaluation_group_count": len(evaluation_groups),
        "trajectory_count": len(raw),
        "valid_evaluation_groups": sum(row["valid"] for row in evaluation_groups),
        "interpretation_limit": (
            "Synthetic method-ablation results only; no real-system or "
            "application validation."
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


def _cell_summaries(groups: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = _FACTOR_NAMES + ("model",)
    rows: list[dict[str, Any]] = []
    for values in itertools.product(
        *(FACTOR_LEVELS[name] for name in _FACTOR_NAMES), MODEL_NAMES
    ):
        condition = dict(zip(keys, values, strict=True))
        cell = [
            row
            for row in groups
            if all(row[key] == value for key, value in condition.items())
        ]
        valid = [row for row in cell if row["valid"]]
        q_values = np.asarray([row["q"] for row in valid], dtype=np.float64)
        p_values = np.asarray([row["p"] for row in valid], dtype=np.float64)
        q_ci = _bootstrap_interval(q_values, _cell_seed(condition, 0))
        p_ci = _bootstrap_interval(p_values, _cell_seed(condition, 1))
        true_q, true_p = SIGNATURES[condition["model"]]
        q_mean = float(np.mean(q_values)) if q_values.size else math.nan
        p_mean = float(np.mean(p_values)) if p_values.size else math.nan
        validity = len(valid) / len(cell) if cell else 0.0
        q_covers = q_ci[0] - 1e-12 <= true_q <= q_ci[1] + 1e-12
        p_covers = p_ci[0] - 1e-12 <= true_p <= p_ci[1] + 1e-12
        criterion = (
            bool(
                validity >= 0.8
                and abs(q_mean - true_q) <= 0.2
                and abs(p_mean - true_p) <= 0.2
                and q_covers
                and p_covers
            )
            if valid
            else False
        )
        rows.append(
            {
                **condition,
                "attempted_groups": len(cell),
                "valid_groups": len(valid),
                "valid_group_rate": validity,
                "mean_q": q_mean,
                "true_q": true_q,
                "mean_q_abs_error": abs(q_mean - true_q),
                "q_ci95_low": q_ci[0],
                "q_ci95_high": q_ci[1],
                "q_ci95_width": q_ci[1] - q_ci[0],
                "q_ci_covers_truth": q_covers,
                "mean_p": p_mean,
                "true_p": true_p,
                "mean_p_abs_error": abs(p_mean - true_p),
                "p_ci95_low": p_ci[0],
                "p_ci95_high": p_ci[1],
                "p_ci95_width": p_ci[1] - p_ci[0],
                "p_ci_covers_truth": p_covers,
                "mean_estimated_boundary": _mean_finite(
                    [row["mean_estimated_boundary"] for row in cell]
                ),
                "mean_rate_fit_valid_fraction": float(
                    np.mean([row["rate_valid_fraction"] for row in cell])
                )
                if cell
                else math.nan,
                "mean_rate_relative_error": _mean_finite(
                    [row["rate_mean_relative_error"] for row in cell]
                ),
                "identifiability_criterion_passed": criterion,
            }
        )
    return rows


def _rate_summaries(raw: list[dict[str, Any]]) -> list[dict[str, Any]]:
    factors = _FACTOR_NAMES
    rows: list[dict[str, Any]] = []
    combinations = itertools.product(
        *(FACTOR_LEVELS[name] for name in factors), MODEL_NAMES, TRUE_MARGINS
    )
    for values in combinations:
        condition = dict(zip((*factors, "model", "margin"), values, strict=True))
        cell = [
            row
            for row in raw
            if all(row[key] == value for key, value in condition.items())
        ]
        good = [row for row in cell if row["rate_valid"]]
        errors = [row["rate_relative_error"] for row in good]
        rows.append(
            {
                **condition,
                "attempted_trajectories": len(cell),
                "valid_rate_fits": len(good),
                "rate_fit_valid_fraction": len(good) / len(cell) if cell else 0.0,
                "mean_relative_rate_error": _mean_finite(errors),
                "median_relative_rate_error": (
                    float(np.median(errors)) if errors else math.nan
                ),
                "mean_relative_equilibrium_error": _mean_finite(
                    [row["equilibrium_relative_error"] for row in cell]
                ),
                "median_relative_equilibrium_error": float(
                    np.median([row["equilibrium_relative_error"] for row in cell])
                )
                if cell
                else math.nan,
            }
        )
    return rows


def _bootstrap_interval(values: NDArray[np.float64], seed: int) -> tuple[float, float]:
    if values.size < 2:
        return math.nan, math.nan
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, values.size, (BOOTSTRAP_REPLICATES, values.size))
    means = np.mean(values[indices], axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(low), float(high)


def _cell_seed(condition: dict[str, Any], suffix: int) -> int:
    payload = json.dumps(condition, sort_keys=True).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little") + suffix


def _mean_finite(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return float(np.mean(finite)) if finite else math.nan


def _factor_contrasts(groups: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Compute paired-seed factorial contrasts and interaction differences."""
    # Each seed-level value first averages its independent five-margin groups.
    group_metrics: dict[tuple[Any, ...], dict[str, float]] = {}
    for row in groups:
        key = (
            row["model"],
            row["replicate_seed"],
            *(row[name] for name in _FACTOR_NAMES),
        )
        group_metrics[key] = {
            "valid_group_probability": float(row["valid"]),
            "absolute_q_error": row["q_abs_error"],
            "absolute_p_error": row["p_abs_error"],
            "relative_rate_error": row["rate_mean_relative_error"],
            "all_rates_valid": float(row["rate_valid_fraction"] == 1.0),
        }
    results: list[dict[str, Any]] = []
    metrics = (
        "valid_group_probability",
        "absolute_q_error",
        "absolute_p_error",
        "relative_rate_error",
    )
    for model in MODEL_NAMES:
        for metric in metrics:
            for contrast_name, contrast_fn in _contrast_functions(metric):
                seed_effects: list[float] = []
                seed_ids = sorted({key[1] for key in group_metrics if key[0] == model})
                for seed in seed_ids:
                    value = contrast_fn(group_metrics, model, seed, metric)
                    if value is not None and math.isfinite(value):
                        seed_effects.append(value)
                ci = _bootstrap_interval(
                    np.asarray(seed_effects, dtype=np.float64),
                    _cell_seed(
                        {"model": model, "contrast": contrast_name, "metric": metric}, 8
                    ),
                )
                mean = float(np.mean(seed_effects)) if seed_effects else math.nan
                results.append(
                    {
                        "model": model,
                        "metric": metric,
                        "contrast": contrast_name,
                        "paired_seed_effect": mean,
                        "ci95_low": ci[0],
                        "ci95_high": ci[1],
                        "matched_seeds": len(seed_effects),
                        "effect_status": (
                            "supported"
                            if len(seed_effects)
                            >= math.ceil(SEED_COUNTS["evaluation"] / 2)
                            and (ci[1] < 0.0 or ci[0] > 0.0)
                            else "inconclusive"
                        ),
                    }
                )
    return results


def _condition_value(
    records: dict[tuple[Any, ...], dict[str, float]],
    model: str,
    seed: int,
    factors: dict[str, Any],
    metric: str,
) -> float | None:
    key = (model, seed, *(factors[name] for name in _FACTOR_NAMES))
    row = records.get(key)
    if row is None:
        return None
    if metric == "relative_rate_error" and not row["all_rates_valid"]:
        return None
    value = row[metric]
    return value if math.isfinite(value) else None


def _contrast_functions(metric: str) -> list[tuple[str, Any]]:
    functions: list[tuple[str, Any]] = []
    contexts = [name for name in _FACTOR_NAMES]

    def main(factor: str, left: Any, right: Any, label: str) -> Any:
        def effect(
            records: dict[tuple[Any, ...], dict[str, float]],
            model: str,
            seed: int,
            outcome: str,
        ) -> float | None:
            diffs: list[float] = []
            other = [name for name in contexts if name != factor]
            values = [FACTOR_LEVELS[name] for name in other]
            for combo in itertools.product(*values):
                base = dict(zip(other, combo, strict=True))
                a, b = dict(base), dict(base)
                a[factor], b[factor] = right, left
                av = _condition_value(records, model, seed, a, outcome)
                bv = _condition_value(records, model, seed, b, outcome)
                if av is not None and bv is not None:
                    diffs.append(av - bv)
            return float(np.mean(diffs)) if diffs else None

        functions.append((label, effect))
        return effect

    main("amplitude_fraction", 0.05, 0.50, "amplitude_0.50_minus_0.05")
    for info in INFORMATION_MODES[1:]:
        main(
            "information_condition",
            "reference",
            info,
            f"information_{info}_minus_reference",
        )
    main("noise_fraction_of_equilibrium", 0.0, 0.05, "noise_0.05_minus_0")
    main("window_tau", 1.0, 3.0, "window_3tau_minus_1tau")

    def interaction(
        first: str,
        levels_a: tuple[Any, Any],
        second: str,
        levels_b: tuple[Any, Any],
        label: str,
    ) -> None:
        def effect(
            records: dict[tuple[Any, ...], dict[str, float]],
            model: str,
            seed: int,
            outcome: str,
        ) -> float | None:
            other = [name for name in contexts if name not in (first, second)]
            diffs: list[float] = []
            for combo in itertools.product(*(FACTOR_LEVELS[name] for name in other)):
                base = dict(zip(other, combo, strict=True))
                cells: dict[tuple[int, int], float] = {}
                for i, a_level in enumerate(levels_a):
                    for j, b_level in enumerate(levels_b):
                        factors = dict(base)
                        factors[first], factors[second] = a_level, b_level
                        value = _condition_value(records, model, seed, factors, outcome)
                        if value is not None:
                            cells[i, j] = value
                if len(cells) == 4:
                    diffs.append(
                        (cells[1, 1] - cells[0, 1]) - (cells[1, 0] - cells[0, 0])
                    )
            return float(np.mean(diffs)) if diffs else None

        functions.append((label, effect))

    interaction(
        "noise_fraction_of_equilibrium",
        (0.0, 0.05),
        "window_tau",
        (1.0, 3.0),
        "interaction_noise_x_window",
    )
    interaction(
        "amplitude_fraction",
        (0.05, 0.50),
        "noise_fraction_of_equilibrium",
        (0.0, 0.05),
        "interaction_amplitude_x_noise",
    )
    interaction(
        "information_condition",
        ("reference", "oracle_both"),
        "window_tau",
        (1.0, 3.0),
        "interaction_information_x_window",
    )
    return functions


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ScientificValidationError(f"refusing to write empty CSV: {path.name}")
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
