"""Frozen Phase 60 experiment for finite-size Kuramoto synchronization."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import shutil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.kuramoto import KuramotoResult, simulate_kuramoto

FROZEN_PROTOCOL_SHA256 = (
    "49A87B327663E08686E4B618CA01B4D02C8A23EA2B3121BFBB2063819A507A86"
)
GAUSSIAN_ID = "gaussian_sigma_1"
LORENTZIAN_ID = "truncated_lorentzian_gamma_1_cutoff_5"
MODEL_IDS = (GAUSSIAN_ID, LORENTZIAN_ID)
POPULATION_SIZES = (64, 256)
REPLICATE_SEEDS = tuple(range(60001, 60009))
COUPLING_FACTORS = (0.6, 0.68, 0.76, 0.84, 0.92, 1.0, 1.08, 1.16, 1.24, 1.32, 1.4)
BOOTSTRAP_REPLICATES = 2000


def critical_coupling(population: str) -> float:
    """Return continuum Kc for the frozen Gaussian or Lorentzian density.

    The relation Kc = 2/(pi*g(0)) is the infinite-population result for a
    centered, symmetric unimodal density under the standard all-to-all
    sinusoidal mean-field coupling. It is not an exact finite-N threshold.
    """
    if population == GAUSSIAN_ID:
        sigma = 1.0
        density_at_zero = 1.0 / (math.sqrt(2.0 * math.pi) * sigma)
    elif population == LORENTZIAN_ID:
        half_width = 1.0
        cutoff = 5.0
        normalization_angle = math.atan(cutoff / half_width)
        density_at_zero = 1.0 / (2.0 * half_width * normalization_angle)
    else:
        raise ScientificValidationError(f"unknown frequency population: {population}")
    return 2.0 / (math.pi * density_at_zero)


def _validate_protocol(protocol_path: Path) -> tuple[dict[str, Any], str]:
    raw = protocol_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest().upper()
    if digest != FROZEN_PROTOCOL_SHA256:
        raise ScientificValidationError("Phase 60 frozen protocol hash mismatch")
    protocol = json.loads(raw)
    if (
        protocol.get("phase") != 60
        or protocol.get("protocol_version") != "1.1"
        or protocol.get("protocol_status") != "frozen_before_final_evaluation"
        or protocol.get("phase32_status") != "BLOCKED_AUTHORIZATION"
        or protocol.get("phase52_status") != "NO-GO"
    ):
        raise ScientificValidationError("Phase 60 protocol status is invalid")
    return protocol, digest


def _frequency_sample(
    population: str, population_size: int, seed: int, population_index: int
) -> tuple[np.ndarray, np.ndarray]:
    frequency_rng = np.random.default_rng(
        np.random.SeedSequence([seed, population_index, population_size, 1])
    )
    phase_rng = np.random.default_rng(
        np.random.SeedSequence([seed, population_index, population_size, 2])
    )
    if population == GAUSSIAN_ID:
        frequencies = frequency_rng.normal(0.0, 1.0, population_size)
    else:
        cutoff_angle = math.atan(5.0)
        angles = frequency_rng.uniform(-cutoff_angle, cutoff_angle, population_size)
        frequencies = np.tan(angles)
    phases = phase_rng.uniform(-math.pi, math.pi, population_size)
    return np.asarray(frequencies, dtype=np.float64), np.asarray(
        phases, dtype=np.float64
    )


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _transition_estimate(
    replicate_curves: np.ndarray,
    coupling_values: np.ndarray,
    *,
    bootstrap_seed: int,
) -> tuple[float, float, float, bool]:
    """Estimate steepest-rise interval midpoint with paired-seed bootstrap CI."""
    mean_curve = np.mean(replicate_curves, axis=0)
    slopes = np.diff(mean_curve) / np.diff(coupling_values)
    selected_index = int(np.argmax(slopes))
    estimate = float(
        (coupling_values[selected_index] + coupling_values[selected_index + 1]) / 2
    )
    rng = np.random.default_rng(bootstrap_seed)
    sample_indices = rng.integers(
        0,
        replicate_curves.shape[0],
        size=(BOOTSTRAP_REPLICATES, replicate_curves.shape[0]),
    )
    bootstrap_curves = np.mean(replicate_curves[sample_indices, :], axis=1)
    bootstrap_slopes = np.diff(bootstrap_curves, axis=1) / np.diff(coupling_values)
    bootstrap_indices = np.argmax(bootstrap_slopes, axis=1)
    bootstrap_estimates = (
        coupling_values[bootstrap_indices] + coupling_values[bootstrap_indices + 1]
    ) / 2.0
    lower, upper = np.percentile(bootstrap_estimates, [2.5, 97.5])
    at_boundary = selected_index in (0, len(coupling_values) - 2)
    return estimate, float(lower), float(upper), at_boundary


def _save_figures(
    output_directory: Path,
    summary_rows: list[dict[str, Any]],
    timeseries_rows: list[dict[str, Any]],
    phase_rows: list[dict[str, Any]],
    protocol: dict[str, Any],
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for axis, population in zip(axes, MODEL_IDS, strict=True):
        predicted = protocol["frequency_populations"][population][
            "continuum_critical_coupling_rad_per_second"
        ]
        for size in POPULATION_SIZES:
            selected = [
                row
                for row in summary_rows
                if row["frequency_population"] == population
                and row["population_size"] == size
            ]
            selected.sort(key=lambda row: row["coupling_rad_per_s"])
            coupling = np.array([row["coupling_rad_per_s"] for row in selected])
            means = np.array([row["mean_coherence"] for row in selected])
            errors = np.array([row["standard_error_coherence"] for row in selected])
            axis.plot(coupling, means, marker="o", label=f"N={size}")
            axis.fill_between(coupling, means - errors, means + errors, alpha=0.18)
        axis.axvline(predicted, color="black", linestyle="--", label="continuum Kc")
        axis.set_title(population.replace("_", " "))
        axis.set_xlabel("Coupling K (rad/s)")
        axis.grid(alpha=0.25)
        axis.legend()
    axes[0].set_ylabel("Mean order parameter R (analysis window)")
    figure.suptitle("Finite-population synchronization across coupling strength")
    figure.tight_layout()
    figure.savefig(
        output_directory / "order_parameter_vs_coupling.png",
        dpi=140,
        metadata={"Software": "Newton Lab"},
    )
    plt.close(figure)

    representative = [
        row for row in timeseries_rows if row["frequency_population"] == GAUSSIAN_ID
    ]
    figure, axis = plt.subplots(figsize=(8, 4.5))
    for factor in protocol["representative_dynamics"]["coupling_factors_of_Kc"]:
        selected = [
            row for row in representative if row["coupling_factor_of_Kc"] == factor
        ]
        selected.sort(key=lambda row: row["time_s"])
        label = "below predicted Kc" if factor < 1.0 else "above predicted Kc"
        axis.plot(
            [row["time_s"] for row in selected],
            [row["coherence"] for row in selected],
            label=f"{label}, K/Kc={factor:.2f}",
        )
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Order parameter R(t)")
    axis.set_title("Representative Gaussian population trajectories (N=256)")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(
        output_directory / "representative_synchronization_traces.png",
        dpi=140,
        metadata={"Software": "Newton Lab"},
    )
    plt.close(figure)

    selected_phases = [
        row for row in phase_rows if row["frequency_population"] == GAUSSIAN_ID
    ]
    figure, axes = plt.subplots(
        1, 2, figsize=(10, 4), sharey=True, constrained_layout=True
    )
    for axis, factor in zip(
        axes,
        protocol["representative_dynamics"]["coupling_factors_of_Kc"],
        strict=True,
    ):
        phases = [
            row["phase_rad"]
            for row in selected_phases
            if row["coupling_factor_of_Kc"] == factor
        ]
        axis.hist(phases, bins=24, range=(-math.pi, math.pi), color="#376a9f")
        axis.set_title(f"K/Kc={factor:.2f}")
        axis.set_xlabel("Phase (rad)")
    axes[0].set_ylabel("Oscillator count")
    figure.suptitle("Final phase snapshots: one representative realization")
    figure.savefig(
        output_directory / "representative_phase_distribution.png",
        dpi=140,
        metadata={"Software": "Newton Lab"},
    )
    plt.close(figure)


def run_synchronization_study(output_directory: Path) -> dict[str, Any]:
    """Run the frozen Phase 60 study in a directory containing protocol.json.

    The output directory may contain only its frozen protocol. Existing results
    are never overwritten; use a new directory for a reproducibility rerun.
    """
    destination = Path(output_directory)
    protocol_path = destination / "protocol.json"
    if not protocol_path.is_file():
        raise FileNotFoundError(f"frozen protocol not found: {protocol_path}")
    unexpected = [
        child.name for child in destination.iterdir() if child.name != "protocol.json"
    ]
    if unexpected:
        raise FileExistsError(
            "refusing to overwrite Phase 60 outputs: " + ", ".join(sorted(unexpected))
        )
    protocol, protocol_hash = _validate_protocol(protocol_path)

    integration = protocol["integration"]
    output_times = np.linspace(
        0.0,
        integration["duration_seconds"],
        int(integration["duration_seconds"] / integration["output_interval_seconds"])
        + 1,
    )
    analysis_mask = (output_times >= integration["analysis_window_start_seconds"]) & (
        output_times <= integration["analysis_window_end_seconds"]
    )
    transient = protocol["transient_diagnostic"]
    earlier_mask = (output_times >= transient["earlier_window_seconds"][0]) & (
        output_times < transient["earlier_window_seconds"][1]
    )
    later_mask = (output_times >= transient["later_window_seconds"][0]) & (
        output_times <= transient["later_window_seconds"][1]
    )
    factors = np.asarray(protocol["coupling_sweep_factors_of_continuum_Kc"])
    seeds = tuple(int(value) for value in protocol["replicate_seeds"])
    sizes = tuple(int(value) for value in protocol["finite_population_sizes"])

    replicate_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    transition_rows: list[dict[str, Any]] = []
    timeseries_rows: list[dict[str, Any]] = []
    phase_rows: list[dict[str, Any]] = []
    integration_count = 0

    for population_index, population in enumerate(MODEL_IDS):
        critical_value = critical_coupling(population)
        for size in sizes:
            replicate_curves = np.empty((len(seeds), len(factors)), dtype=np.float64)
            for replicate_index, seed in enumerate(seeds):
                frequencies, initial_phases = _frequency_sample(
                    population, size, seed, population_index
                )
                for coupling_index, factor in enumerate(factors):
                    coupling = float(factor * critical_value)
                    result: KuramotoResult = simulate_kuramoto(
                        frequencies,
                        initial_phases,
                        coupling,
                        output_times,
                        method=integration["method"],
                        relative_tolerance=integration["relative_tolerance"],
                        absolute_tolerance=integration["absolute_tolerance"],
                        maximum_step_s=integration["maximum_step_seconds"],
                    )
                    integration_count += 1
                    mean_coherence = float(np.mean(result.coherence[analysis_mask]))
                    earlier_mean = float(np.mean(result.coherence[earlier_mask]))
                    later_mean = float(np.mean(result.coherence[later_mask]))
                    replicate_curves[replicate_index, coupling_index] = mean_coherence
                    replicate_rows.append(
                        {
                            "frequency_population": population,
                            "population_size": size,
                            "replicate_seed": seed,
                            "coupling_factor_of_Kc": float(factor),
                            "coupling_rad_per_s": coupling,
                            "mean_coherence": mean_coherence,
                            "earlier_half_mean_coherence": earlier_mean,
                            "later_half_mean_coherence": later_mean,
                            "late_minus_early_coherence": later_mean - earlier_mean,
                            "temporal_sd_coherence": float(
                                np.std(result.coherence[analysis_mask], ddof=1)
                            ),
                            "final_coherence": float(result.coherence[-1]),
                            "function_evaluations": result.function_evaluations,
                        }
                    )
                    representative = protocol["representative_dynamics"]
                    if (
                        population == representative["frequency_population"]
                        and size == representative["population_size"]
                        and seed == representative["replicate_seed"]
                        and float(factor) in representative["coupling_factors_of_Kc"]
                    ):
                        for time, coherence in zip(
                            result.time_s, result.coherence, strict=True
                        ):
                            timeseries_rows.append(
                                {
                                    "frequency_population": population,
                                    "population_size": size,
                                    "replicate_seed": seed,
                                    "coupling_factor_of_Kc": float(factor),
                                    "coupling_rad_per_s": coupling,
                                    "time_s": float(time),
                                    "coherence": float(coherence),
                                }
                            )
                        for phase in result.phases_rad[-1]:
                            phase_rows.append(
                                {
                                    "frequency_population": population,
                                    "population_size": size,
                                    "replicate_seed": seed,
                                    "coupling_factor_of_Kc": float(factor),
                                    "phase_rad": float(
                                        (phase + math.pi) % (2 * math.pi) - math.pi
                                    ),
                                }
                            )

            coupling_values = factors * critical_value
            for coupling_index, factor in enumerate(factors):
                values = replicate_curves[:, coupling_index]
                summary_rows.append(
                    {
                        "frequency_population": population,
                        "population_size": size,
                        "coupling_factor_of_Kc": float(factor),
                        "coupling_rad_per_s": float(coupling_values[coupling_index]),
                        "continuum_Kc_rad_per_s": critical_value,
                        "mean_coherence": float(np.mean(values)),
                        "sd_across_replicates": float(np.std(values, ddof=1)),
                        "standard_error_coherence": float(
                            np.std(values, ddof=1) / math.sqrt(len(values))
                        ),
                        "replicate_count": len(values),
                    }
                )

            estimate, lower, upper, at_boundary = _transition_estimate(
                replicate_curves,
                coupling_values,
                bootstrap_seed=606060 + population_index * 100 + size,
            )
            transition_rows.append(
                {
                    "frequency_population": population,
                    "population_size": size,
                    "continuum_Kc_rad_per_s": critical_value,
                    "estimated_transition_rad_per_s": estimate,
                    "bootstrap_95_percentile_low_rad_per_s": lower,
                    "bootstrap_95_percentile_high_rad_per_s": upper,
                    "grid_interval_width_rad_per_s": float(np.diff(coupling_values)[0]),
                    "steepest_rise_at_sweep_boundary": at_boundary,
                    "replicate_count": len(seeds),
                    "operational_definition": "midpoint of max-slope coupling interval",
                }
            )

    _write_csv(
        destination / "replicate_results.csv",
        [
            "frequency_population",
            "population_size",
            "replicate_seed",
            "coupling_factor_of_Kc",
            "coupling_rad_per_s",
            "mean_coherence",
            "earlier_half_mean_coherence",
            "later_half_mean_coherence",
            "late_minus_early_coherence",
            "temporal_sd_coherence",
            "final_coherence",
            "function_evaluations",
        ],
        replicate_rows,
    )
    _write_csv(
        destination / "synchronization_summary.csv",
        [
            "frequency_population",
            "population_size",
            "coupling_factor_of_Kc",
            "coupling_rad_per_s",
            "continuum_Kc_rad_per_s",
            "mean_coherence",
            "sd_across_replicates",
            "standard_error_coherence",
            "replicate_count",
        ],
        summary_rows,
    )
    _write_csv(
        destination / "transition_estimates.csv",
        [
            "frequency_population",
            "population_size",
            "continuum_Kc_rad_per_s",
            "estimated_transition_rad_per_s",
            "bootstrap_95_percentile_low_rad_per_s",
            "bootstrap_95_percentile_high_rad_per_s",
            "grid_interval_width_rad_per_s",
            "steepest_rise_at_sweep_boundary",
            "replicate_count",
            "operational_definition",
        ],
        transition_rows,
    )
    _write_csv(
        destination / "representative_timeseries.csv",
        [
            "frequency_population",
            "population_size",
            "replicate_seed",
            "coupling_factor_of_Kc",
            "coupling_rad_per_s",
            "time_s",
            "coherence",
        ],
        timeseries_rows,
    )
    _write_csv(
        destination / "representative_phase_snapshot.csv",
        [
            "frequency_population",
            "population_size",
            "replicate_seed",
            "coupling_factor_of_Kc",
            "phase_rad",
        ],
        phase_rows,
    )
    _save_figures(destination, summary_rows, timeseries_rows, phase_rows, protocol)

    metadata = {
        "phase": 60,
        "study": "finite_size_kuramoto_synchronization",
        "protocol_sha256": protocol_hash,
        "successful_integrations": integration_count,
        "expected_integrations": len(MODEL_IDS)
        * len(sizes)
        * len(seeds)
        * len(factors),
        "replicate_rows": len(replicate_rows),
        "summary_rows": len(summary_rows),
        "transition_rows": len(transition_rows),
        "representative_time_rows": len(timeseries_rows),
        "representative_phase_rows": len(phase_rows),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "matplotlib_version": matplotlib.__version__,
        "phase32_status": "BLOCKED_AUTHORIZATION",
        "phase52_status": "NO-GO",
    }
    (destination / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    manifest: dict[str, str] = {}
    for path in sorted(destination.iterdir()):
        if path.is_file() and path.name != "sha256_manifest.json":
            manifest[path.name] = hashlib.sha256(path.read_bytes()).hexdigest().upper()
    (destination / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata


def copy_protocol_for_rerun(source_directory: Path, rerun_directory: Path) -> None:
    """Copy the frozen protocol into a fresh directory for a safe rerun."""
    destination = Path(rerun_directory)
    destination.mkdir(parents=True, exist_ok=True)
    if any(destination.iterdir()):
        raise FileExistsError(f"rerun directory is not empty: {destination}")
    shutil.copyfile(
        Path(source_directory) / "protocol.json", destination / "protocol.json"
    )
