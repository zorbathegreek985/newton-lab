"""Phase 61 finite-size and coupling-resolution study for Kuramoto onset."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import NDArray

from newton_lab.exceptions import ScientificValidationError
from newton_lab.kuramoto import simulate_kuramoto

FROZEN_PROTOCOL_SHA256 = (
    "ED307A5E8E7D6619F3D34CA1D2B5BC5B2EF04399DBB0F914EE0352C5FEC5DBF8"
)
PROTOCOL_NAME = "protocol.json"
POPULATION_SIZES = (64, 256, 512, 1024)
REPLICATE_SEEDS = tuple(range(61001, 61009))
COUPLING_FACTORS = (
    0.75,
    0.85,
    0.90,
    0.925,
    0.95,
    0.975,
    1.0,
    1.025,
    1.05,
    1.075,
    1.10,
    1.125,
    1.15,
    1.175,
    1.20,
    1.25,
    1.35,
)
BOOTSTRAP_REPLICATES = 2000
LORENTZIAN_GAMMA = 1.0
LORENTZIAN_CUTOFF = 5.0


def truncated_lorentzian_normalization_constant(
    gamma_rad_per_s: float, cutoff_rad_per_s: float
) -> float:
    """Return C making C*gamma/(pi*(omega^2+gamma^2)) integrate to one.

    The unnormalized density is restricted to the symmetric interval
    [-cutoff_rad_per_s, cutoff_rad_per_s] and is zero outside.
    """
    gamma = float(gamma_rad_per_s)
    cutoff = float(cutoff_rad_per_s)
    if not math.isfinite(gamma) or gamma <= 0.0:
        raise ScientificValidationError("gamma_rad_per_s must be finite and positive")
    if not math.isfinite(cutoff) or cutoff <= 0.0:
        raise ScientificValidationError("cutoff_rad_per_s must be finite and positive")
    return math.pi / (2.0 * math.atan(cutoff / gamma))


def truncated_lorentzian_density(
    omega_rad_per_s: NDArray[np.float64] | float,
    gamma_rad_per_s: float = LORENTZIAN_GAMMA,
    cutoff_rad_per_s: float = LORENTZIAN_CUTOFF,
) -> NDArray[np.float64] | float:
    """Evaluate the normalized symmetric truncated-Lorentzian density."""
    normalization = truncated_lorentzian_normalization_constant(
        gamma_rad_per_s, cutoff_rad_per_s
    )
    omega = np.asarray(omega_rad_per_s, dtype=np.float64)
    gamma = float(gamma_rad_per_s)
    cutoff = float(cutoff_rad_per_s)
    unnormalized = gamma / (math.pi * (omega**2 + gamma**2))
    density = np.where(np.abs(omega) <= cutoff, normalization * unnormalized, 0.0)
    if np.ndim(omega_rad_per_s) == 0:
        return float(density)
    return np.asarray(density, dtype=np.float64)


def truncated_lorentzian_critical_coupling(
    gamma_rad_per_s: float = LORENTZIAN_GAMMA,
    cutoff_rad_per_s: float = LORENTZIAN_CUTOFF,
) -> float:
    """Return continuum Kc=2/(pi*g(0)) for the normalized truncated density."""
    density_at_zero = float(
        truncated_lorentzian_density(0.0, gamma_rad_per_s, cutoff_rad_per_s)
    )
    return 2.0 / (math.pi * density_at_zero)


def _load_protocol(protocol_path: Path) -> tuple[dict[str, Any], str]:
    raw = protocol_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest().upper()
    if digest != FROZEN_PROTOCOL_SHA256:
        raise ScientificValidationError("Phase 61 frozen protocol hash mismatch")
    protocol = json.loads(raw)
    if (
        protocol.get("phase") != 61
        or protocol.get("protocol_version") != "1.0"
        or protocol.get("protocol_status") != "frozen_before_final_evaluation"
        or protocol.get("phase32_status") != "BLOCKED_AUTHORIZATION"
        or protocol.get("phase52_status") != "NO-GO"
    ):
        raise ScientificValidationError("Phase 61 protocol status is invalid")
    return protocol, digest


def _nested_population_sample(
    seed: int, maximum_size: int
) -> tuple[np.ndarray, np.ndarray]:
    frequency_rng = np.random.default_rng(
        np.random.SeedSequence([seed, maximum_size, 1])
    )
    phase_rng = np.random.default_rng(np.random.SeedSequence([seed, maximum_size, 2]))
    cutoff_angle = math.atan(LORENTZIAN_CUTOFF / LORENTZIAN_GAMMA)
    angles = frequency_rng.uniform(-cutoff_angle, cutoff_angle, maximum_size)
    frequencies = LORENTZIAN_GAMMA * np.tan(angles)
    phases = phase_rng.uniform(-math.pi, math.pi, maximum_size)
    return (
        np.asarray(frequencies, dtype=np.float64),
        np.asarray(phases, dtype=np.float64),
    )


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _estimate_transition(
    curves: np.ndarray,
    coupling_values: np.ndarray,
    *,
    bootstrap_seed: int,
) -> dict[str, Any]:
    """Estimate a max-slope grid bracket and seed-bootstrap uncertainty."""
    mean_curve = np.mean(curves, axis=0)
    slopes = np.diff(mean_curve) / np.diff(coupling_values)
    index = int(np.argmax(slopes))
    bracket_low = float(coupling_values[index])
    bracket_high = float(coupling_values[index + 1])
    boundary = index == 0 or index == len(coupling_values) - 2
    positive = bool(slopes[index] > 0.0)
    if not positive:
        return {
            "transition_estimate_rad_per_s": None,
            "censoring": "no_positive_slope",
            "grid_bracket_low_rad_per_s": bracket_low,
            "grid_bracket_high_rad_per_s": bracket_high,
            "bootstrap_low_rad_per_s": None,
            "bootstrap_high_rad_per_s": None,
            "bootstrap_boundary_fraction": None,
            "bootstrap_interval_censored": True,
            "seed_interval_width_rad_per_s": None,
            "selected_grid_interval_width_rad_per_s": bracket_high - bracket_low,
        }

    rng = np.random.default_rng(bootstrap_seed)
    samples = rng.integers(
        0, curves.shape[0], size=(BOOTSTRAP_REPLICATES, curves.shape[0])
    )
    sample_curves = np.mean(curves[samples, :], axis=1)
    sample_slopes = np.diff(sample_curves, axis=1) / np.diff(coupling_values)
    sample_indices = np.argmax(sample_slopes, axis=1)
    midpoints = (
        coupling_values[sample_indices] + coupling_values[sample_indices + 1]
    ) / 2.0
    boundary_fraction = float(
        np.mean((sample_indices == 0) | (sample_indices == len(coupling_values) - 2))
    )
    bootstrap_interval_censored = boundary_fraction > 0.025
    if bootstrap_interval_censored:
        bootstrap_low: float | None = None
        bootstrap_high: float | None = None
        seed_interval_width: float | None = None
    else:
        percentile_bounds = np.percentile(midpoints, [2.5, 97.5])
        bootstrap_low, bootstrap_high = (
            float(percentile_bounds[0]),
            float(percentile_bounds[1]),
        )
        seed_interval_width = float(bootstrap_high - bootstrap_low)
    if boundary:
        censoring = "left" if index == 0 else "right"
        estimate: float | None = None
    else:
        censoring = "none"
        estimate = float((bracket_low + bracket_high) / 2.0)
    return {
        "transition_estimate_rad_per_s": estimate,
        "censoring": censoring,
        "grid_bracket_low_rad_per_s": bracket_low,
        "grid_bracket_high_rad_per_s": bracket_high,
        "bootstrap_low_rad_per_s": bootstrap_low,
        "bootstrap_high_rad_per_s": bootstrap_high,
        "bootstrap_boundary_fraction": boundary_fraction,
        "bootstrap_interval_censored": bootstrap_interval_censored,
        "seed_interval_width_rad_per_s": seed_interval_width,
        "selected_grid_interval_width_rad_per_s": bracket_high - bracket_low,
    }


def _curve_mean_and_responses(
    frequencies: np.ndarray,
    phases: np.ndarray,
    size: int,
    factors: np.ndarray,
    critical_value: float,
    times: np.ndarray,
    analysis_mask: np.ndarray,
    integration: dict[str, Any],
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    curve = np.empty(len(factors), dtype=np.float64)
    rows: list[dict[str, Any]] = []
    for coupling_index, factor in enumerate(factors):
        coupling = float(factor * critical_value)
        result = simulate_kuramoto(
            frequencies[:size],
            phases[:size],
            coupling,
            times,
            method=integration["method"],
            relative_tolerance=integration["relative_tolerance"],
            absolute_tolerance=integration["absolute_tolerance"],
            maximum_step_s=integration["maximum_step_seconds"],
        )
        early_mask = (times >= 80.0) & (times < 100.0)
        late_mask = (times >= 100.0) & (times <= 120.0)
        early_mean = float(np.mean(result.coherence[early_mask]))
        late_mean = float(np.mean(result.coherence[late_mask]))
        mean_coherence = float(np.mean(result.coherence[analysis_mask]))
        curve[coupling_index] = mean_coherence
        rows.append(
            {
                "coupling_factor_of_Kc": float(factor),
                "coupling_rad_per_s": coupling,
                "mean_coherence": mean_coherence,
                "early_window_mean_coherence": early_mean,
                "late_window_mean_coherence": late_mean,
                "late_minus_early_coherence": late_mean - early_mean,
                "final_coherence": float(result.coherence[-1]),
                "function_evaluations": result.function_evaluations,
            }
        )
    return curve, rows


def _save_figures(
    output: Path,
    summary_rows: list[dict[str, Any]],
    transition_rows: list[dict[str, Any]],
    numerical_rows: list[dict[str, Any]],
    critical_value: float,
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for axis, normalized in zip(axes, (False, True), strict=True):
        for size in POPULATION_SIZES:
            selected = [row for row in summary_rows if row["population_size"] == size]
            selected.sort(key=lambda row: row["coupling_rad_per_s"])
            x = np.array(
                [
                    row["coupling_factor_of_Kc"]
                    if normalized
                    else row["coupling_rad_per_s"]
                    for row in selected
                ]
            )
            y = np.array([row["mean_coherence"] for row in selected])
            e = np.array([row["standard_error_coherence"] for row in selected])
            axis.plot(x, y, marker="o", markersize=3, label=f"N={size}")
            axis.fill_between(x, y - e, y + e, alpha=0.14)
        axis.axvline(
            1.0 if normalized else critical_value, color="black", linestyle="--"
        )
        axis.set_xlabel("K/Kc" if normalized else "Coupling K (rad/s)")
        axis.grid(alpha=0.25)
        axis.legend()
    axes[0].set_ylabel("Mean order parameter R (80–120 s)")
    axes[0].set_title("Absolute coupling")
    axes[1].set_title("Coupling normalized by continuum Kc")
    figure.suptitle("Truncated-Lorentzian finite-size synchronization")
    figure.tight_layout()
    figure.savefig(
        output / "coherence_by_coupling_and_size.png",
        dpi=140,
        metadata={"Software": "Newton Lab"},
    )
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(8, 4.8))
    for row in transition_rows:
        if (
            row["transition_estimate_rad_per_s"] is None
            or row["bootstrap_low_rad_per_s"] is None
            or row["bootstrap_high_rad_per_s"] is None
        ):
            continue
        estimate = row["transition_estimate_rad_per_s"]
        low = row["bootstrap_low_rad_per_s"]
        high = row["bootstrap_high_rad_per_s"]
        axis.vlines(row["population_size"], low, high, linewidth=1.5)
        axis.hlines(
            [low, high],
            row["population_size"] * 0.97,
            row["population_size"] * 1.03,
            linewidth=1.5,
        )
        axis.plot(
            row["population_size"],
            estimate,
            "o",
            label=f"N={row['population_size']}",
        )
    axis.axhline(critical_value, color="black", linestyle="--", label="continuum Kc")
    axis.set_xscale("log", base=2)
    axis.set_xticks(POPULATION_SIZES, labels=[str(size) for size in POPULATION_SIZES])
    axis.set_xlabel("Population size N")
    axis.set_ylabel("Operational transition K (rad/s)")
    axis.set_title("Seed-bootstrap interval; censored estimates omitted")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(
        output / "transition_vs_population_size.png",
        dpi=140,
        metadata={"Software": "Newton Lab"},
    )
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(8, 4.5))
    axis.plot(
        [row["coupling_factor_of_Kc"] for row in numerical_rows],
        [row["primary_mean_coherence"] for row in numerical_rows],
        marker="o",
        label="primary settings",
    )
    axis.plot(
        [row["coupling_factor_of_Kc"] for row in numerical_rows],
        [row["refined_mean_coherence"] for row in numerical_rows],
        marker="x",
        label="refined settings",
    )
    axis.set_xlabel("K/Kc")
    axis.set_ylabel("Mean coherence (four paired seeds)")
    axis.set_title("Largest-population numerical-resolution check")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(
        output / "numerical_resolution_check.png",
        dpi=140,
        metadata={"Software": "Newton Lab"},
    )
    plt.close(figure)


def run_transition_refinement_study(output_directory: Path) -> dict[str, Any]:
    """Run the frozen Phase 61 study without modifying existing outputs."""
    output = Path(output_directory)
    protocol_path = output / PROTOCOL_NAME
    if not protocol_path.is_file():
        raise FileNotFoundError(f"frozen Phase 61 protocol not found: {protocol_path}")
    unexpected = [path.name for path in output.iterdir() if path.name != PROTOCOL_NAME]
    if unexpected:
        raise FileExistsError(
            "refusing to overwrite Phase 61 outputs: " + ", ".join(sorted(unexpected))
        )
    protocol, protocol_hash = _load_protocol(protocol_path)
    critical_value = truncated_lorentzian_critical_coupling()
    if not math.isclose(
        critical_value,
        protocol["theoretical_reference"]["continuum_Kc_rad_per_second"],
        rel_tol=1e-14,
        abs_tol=0.0,
    ):
        raise ScientificValidationError("protocol critical coupling is inconsistent")

    integration = protocol["integration"]
    times = np.linspace(
        0.0,
        integration["duration_seconds"],
        int(integration["duration_seconds"] / integration["output_interval_seconds"])
        + 1,
    )
    analysis_mask = (times >= integration["analysis_window_start_seconds"]) & (
        times <= integration["analysis_window_end_seconds"]
    )
    factors = np.asarray(protocol["coupling_factors_of_continuum_Kc"], dtype=np.float64)
    seeds = tuple(int(seed) for seed in protocol["replicate_seeds"])
    sizes = tuple(int(size) for size in protocol["finite_population_sizes"])
    maximum_size = max(sizes)
    replicate_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    transition_rows: list[dict[str, Any]] = []
    curves_by_size: dict[int, np.ndarray] = {}
    integration_count = 0

    for size in sizes:
        curves = np.empty((len(seeds), len(factors)), dtype=np.float64)
        for seed_index, seed in enumerate(seeds):
            frequencies, phases = _nested_population_sample(seed, maximum_size)
            curve, observations = _curve_mean_and_responses(
                frequencies,
                phases,
                size,
                factors,
                critical_value,
                times,
                analysis_mask,
                integration,
            )
            integration_count += len(factors)
            curves[seed_index, :] = curve
            for observation in observations:
                replicate_rows.append(
                    {
                        "population_size": size,
                        "replicate_seed": seed,
                        **observation,
                    }
                )
        curves_by_size[size] = curves

        coupling_values = factors * critical_value
        for index, factor in enumerate(factors):
            values = curves[:, index]
            summary_rows.append(
                {
                    "population_size": size,
                    "coupling_factor_of_Kc": float(factor),
                    "coupling_rad_per_s": float(coupling_values[index]),
                    "continuum_Kc_rad_per_s": critical_value,
                    "mean_coherence": float(np.mean(values)),
                    "sd_across_seeds": float(np.std(values, ddof=1)),
                    "standard_error_coherence": float(
                        np.std(values, ddof=1) / math.sqrt(len(values))
                    ),
                    "replicate_count": len(seeds),
                }
            )

        estimate = _estimate_transition(
            curves,
            coupling_values,
            bootstrap_seed=protocol["transition_estimator"]["bootstrap_seed"] + size,
        )
        estimate_value = estimate["transition_estimate_rad_per_s"]
        offset = (
            None if estimate_value is None else float(estimate_value - critical_value)
        )
        transition_rows.append(
            {
                "population_size": size,
                "continuum_Kc_rad_per_s": critical_value,
                **estimate,
                "offset_from_continuum_rad_per_s": offset,
                "offset_from_continuum_fraction": (
                    None if offset is None else float(offset / critical_value)
                ),
                "replicate_count": len(seeds),
            }
        )

    numerical_protocol = protocol["numerical_resolution_check"]
    numerical_factors = np.asarray(
        numerical_protocol["coupling_factors_of_continuum_Kc"], dtype=np.float64
    )
    numerical_sizes = [int(value) for value in numerical_protocol["replicate_seeds"]]
    numerical_size = int(numerical_protocol["population_size"])
    numerical_rows: list[dict[str, Any]] = []
    primary_numerical_curves = np.empty(
        (len(numerical_sizes), len(numerical_factors)), dtype=np.float64
    )
    refined_numerical_curves = np.empty_like(primary_numerical_curves)
    index_by_factor = {float(factor): index for index, factor in enumerate(factors)}
    refined_integration = {
        "method": numerical_protocol["method"],
        "relative_tolerance": numerical_protocol["relative_tolerance"],
        "absolute_tolerance": numerical_protocol["absolute_tolerance"],
        "maximum_step_seconds": numerical_protocol["maximum_step_seconds"],
    }
    for seed_index, seed in enumerate(numerical_sizes):
        frequencies, phases = _nested_population_sample(seed, maximum_size)
        primary_curve = curves_by_size[numerical_size][seed_index, :]
        for factor_index, factor in enumerate(numerical_factors):
            primary_numerical_curves[seed_index, factor_index] = primary_curve[
                index_by_factor[float(factor)]
            ]
        refined_curve, _observations = _curve_mean_and_responses(
            frequencies,
            phases,
            numerical_size,
            numerical_factors,
            critical_value,
            times,
            analysis_mask,
            refined_integration,
        )
        refined_numerical_curves[seed_index, :] = refined_curve
        integration_count += len(numerical_factors)
    for index, factor in enumerate(numerical_factors):
        numerical_rows.append(
            {
                "coupling_factor_of_Kc": float(factor),
                "coupling_rad_per_s": float(factor * critical_value),
                "primary_mean_coherence": float(
                    np.mean(primary_numerical_curves[:, index])
                ),
                "refined_mean_coherence": float(
                    np.mean(refined_numerical_curves[:, index])
                ),
                "absolute_mean_coherence_difference": float(
                    abs(
                        np.mean(refined_numerical_curves[:, index])
                        - np.mean(primary_numerical_curves[:, index])
                    )
                ),
                "paired_seed_count": len(numerical_sizes),
            }
        )

    primary_numerical_transition = _estimate_transition(
        primary_numerical_curves,
        numerical_factors * critical_value,
        bootstrap_seed=616262,
    )
    refined_numerical_transition = _estimate_transition(
        refined_numerical_curves,
        numerical_factors * critical_value,
        bootstrap_seed=616262,
    )
    primary_estimate = primary_numerical_transition["transition_estimate_rad_per_s"]
    refined_estimate = refined_numerical_transition["transition_estimate_rad_per_s"]
    if primary_estimate is None or refined_estimate is None:
        numerical_shift: float | None = None
    else:
        numerical_shift = abs(refined_estimate - primary_estimate)

    _write_csv(
        output / "replicate_results.csv",
        [
            "population_size",
            "replicate_seed",
            "coupling_factor_of_Kc",
            "coupling_rad_per_s",
            "mean_coherence",
            "early_window_mean_coherence",
            "late_window_mean_coherence",
            "late_minus_early_coherence",
            "final_coherence",
            "function_evaluations",
        ],
        replicate_rows,
    )
    _write_csv(
        output / "synchronization_summary.csv",
        [
            "population_size",
            "coupling_factor_of_Kc",
            "coupling_rad_per_s",
            "continuum_Kc_rad_per_s",
            "mean_coherence",
            "sd_across_seeds",
            "standard_error_coherence",
            "replicate_count",
        ],
        summary_rows,
    )
    _write_csv(
        output / "transition_estimates.csv",
        [
            "population_size",
            "continuum_Kc_rad_per_s",
            "transition_estimate_rad_per_s",
            "censoring",
            "grid_bracket_low_rad_per_s",
            "grid_bracket_high_rad_per_s",
            "bootstrap_low_rad_per_s",
            "bootstrap_high_rad_per_s",
            "bootstrap_boundary_fraction",
            "bootstrap_interval_censored",
            "seed_interval_width_rad_per_s",
            "selected_grid_interval_width_rad_per_s",
            "offset_from_continuum_rad_per_s",
            "offset_from_continuum_fraction",
            "replicate_count",
        ],
        transition_rows,
    )
    _write_csv(
        output / "numerical_resolution_check.csv",
        [
            "coupling_factor_of_Kc",
            "coupling_rad_per_s",
            "primary_mean_coherence",
            "refined_mean_coherence",
            "absolute_mean_coherence_difference",
            "paired_seed_count",
        ],
        numerical_rows,
    )

    fine_step = 0.025 * critical_value
    numerical_censored = primary_estimate is None or refined_estimate is None
    numerical_resolution_shift = (
        None if numerical_shift is None else float(numerical_shift / fine_step)
    )
    transition_512 = next(
        row for row in transition_rows if row["population_size"] == 512
    )
    transition_1024 = next(
        row for row in transition_rows if row["population_size"] == 1024
    )
    high_size_rows = (transition_512, transition_1024)
    intervals_cover_reference = all(
        row["censoring"] == "none"
        and not row["bootstrap_interval_censored"]
        and row["bootstrap_low_rad_per_s"] is not None
        and row["bootstrap_high_rad_per_s"] is not None
        and row["seed_interval_width_rad_per_s"] is not None
        and row["bootstrap_low_rad_per_s"]
        <= critical_value
        <= row["bootstrap_high_rad_per_s"]
        and row["seed_interval_width_rad_per_s"] <= 0.10 * critical_value
        for row in high_size_rows
    )
    high_size_disagree_same_side = all(
        row["censoring"] == "none"
        and not row["bootstrap_interval_censored"]
        and row["bootstrap_low_rad_per_s"] is not None
        and row["bootstrap_high_rad_per_s"] is not None
        and (
            row["bootstrap_high_rad_per_s"] < critical_value
            or row["bootstrap_low_rad_per_s"] > critical_value
        )
        for row in high_size_rows
    ) and (
        (transition_512["offset_from_continuum_rad_per_s"] or 0.0)
        * (transition_1024["offset_from_continuum_rad_per_s"] or 0.0)
        > 0.0
    )
    point_difference_high_sizes = (
        None
        if transition_512["transition_estimate_rad_per_s"] is None
        or transition_1024["transition_estimate_rad_per_s"] is None
        else abs(
            transition_512["transition_estimate_rad_per_s"]
            - transition_1024["transition_estimate_rad_per_s"]
        )
    )
    disagreement_stable_sizes = (
        high_size_disagree_same_side
        and point_difference_high_sizes is not None
        and point_difference_high_sizes <= fine_step
    )
    numerically_stable = (
        not numerical_censored
        and numerical_shift is not None
        and numerical_shift <= fine_step
    )
    if intervals_cover_reference and numerically_stable:
        conclusion = "agreement"
    elif disagreement_stable_sizes and numerically_stable:
        conclusion = "persistent_disagreement"
    else:
        conclusion = "inconclusive"

    metadata = {
        "phase": 61,
        "study": "kuramoto_truncated_lorentzian_finite_size_resolution",
        "protocol_sha256": protocol_hash,
        "continuum_critical_coupling_rad_per_s": critical_value,
        "primary_integrations": len(sizes) * len(seeds) * len(factors),
        "numerical_resolution_integrations": len(numerical_sizes)
        * len(numerical_factors),
        "total_successful_integrations": integration_count,
        "replicate_rows": len(replicate_rows),
        "summary_rows": len(summary_rows),
        "transition_rows": len(transition_rows),
        "numerical_resolution_rows": len(numerical_rows),
        "numerical_transition_shift_rad_per_s": numerical_shift,
        "numerical_transition_shift_in_fine_grid_intervals": numerical_resolution_shift,
        "conclusion_under_frozen_criteria": conclusion,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "matplotlib_version": matplotlib.__version__,
        "phase32_status": "BLOCKED_AUTHORIZATION",
        "phase52_status": "NO-GO",
    }
    (output / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _save_figures(output, summary_rows, transition_rows, numerical_rows, critical_value)

    manifest: dict[str, str] = {}
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "sha256_manifest.json":
            manifest[path.name] = hashlib.sha256(path.read_bytes()).hexdigest().upper()
    (output / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata
