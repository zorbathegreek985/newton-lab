"""Phase 49: compare slowing recovery in heat diffusion and a pendulum."""

from __future__ import annotations

import csv
import json
import math
import platform
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy  # type: ignore[import-untyped]
from scipy.signal import find_peaks  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.heat_diffusion_case_study import (
    HeatDiffusionModel,
    SpatialMode,
    evaluate_heat_diffusion,
)
from newton_lab.modal_research import heat_mode_decay_rate
from newton_lab.pendulum import simulate_damped_pendulum

DISTANCES = (0.025, 0.05, 0.1, 0.2, 0.4)
PERTURBATION_AMPLITUDES = (0.01, 0.3)
SOLVER_SETTINGS = ((1e-10, 1e-12), (1e-12, 1e-14))
MODE_NUMBER = 1
OBSERVATION_SAMPLES = 401
RECOVERY_E_FOLDS = 12.0
LOCAL_RATE_RELATIVE_ERROR_LIMIT = 0.02
AMPLITUDE_SENSITIVITY_LIMIT = 0.02
EXPONENT_ACCEPTANCE = (0.95, 1.05)
EXPONENT_R_SQUARED_MINIMUM = 0.999


@dataclass(frozen=True)
class RecoveryCondition:
    """One model, parameter, perturbation, and solver setting."""

    model: str
    dimensionless_distance: float
    physical_parameter_value: float
    perturbation_amplitude: float
    solver_rtol: float | None
    solver_atol: float | None
    output_sample_count: int
    equilibrium: str
    stability: str
    predicted_rate_per_s: float
    measured_rate_per_s: float
    relative_rate_error: float
    fit_r_squared: float
    measurement_points: int
    measurement_window_s: float
    numerical_convergence_difference: float | None
    local_approximation_accepted: bool
    diagnostics: str


@dataclass(frozen=True)
class ModelSummary:
    """Frozen-range log-log fit and model-level measurement summary."""

    model: str
    fitted_power_law_exponent: float
    exponent_r_squared: float
    predeclared_distance_min: float
    predeclared_distance_max: float
    maximum_primary_rate_relative_error: float
    maximum_amplitude_sensitivity_relative_difference: float
    all_primary_rates_within_tolerance: bool
    amplitude_sensitivity_within_tolerance: bool
    exponent_within_predeclared_interval: bool
    exponent_r_squared_accepted: bool
    qualitative_slowing_observed: bool
    model_prediction_supported: bool


def _fit_log_rate(times: np.ndarray, amplitudes: np.ndarray) -> tuple[float, float]:
    """Fit log amplitude against time; return decay rate and R-squared."""
    t = np.asarray(times, dtype=np.float64)
    a = np.asarray(amplitudes, dtype=np.float64)
    if (
        t.ndim != 1
        or a.ndim != 1
        or t.size != a.size
        or t.size < 3
        or not np.isfinite(t).all()
        or not np.isfinite(a).all()
        or np.any(np.diff(t) <= 0.0)
        or np.any(a <= 0.0)
    ):
        raise ScientificValidationError("rate fit requires positive finite amplitudes")
    log_amplitude = np.log(a / a[0])
    slope, intercept = np.polyfit(t, log_amplitude, 1)
    predicted = slope * t + intercept
    residual_sum = float(np.sum((log_amplitude - predicted) ** 2))
    total_sum = float(np.sum((log_amplitude - np.mean(log_amplitude)) ** 2))
    if slope >= 0.0 or total_sum == 0.0:
        raise ScientificValidationError("recovery fit did not resolve positive decay")
    return float(-slope), float(1.0 - residual_sum / total_sum)


def pendulum_linearized_spectrum(
    damping_ratio: float,
    *,
    natural_frequency_rad_per_s: float = 1.0,
) -> tuple[complex, complex, float, float]:
    """Return eigenvalues, envelope rate and damped frequency near equilibrium.

    ``damping_ratio`` is ``b/(2*m*L**2*omega0)``. This helper deliberately
    supports only the underdamped regime used in this study.
    """
    if not math.isfinite(damping_ratio) or damping_ratio <= 0.0 or damping_ratio >= 1.0:
        raise ScientificValidationError("damping_ratio must be between zero and one")
    if (
        not math.isfinite(natural_frequency_rad_per_s)
        or natural_frequency_rad_per_s <= 0.0
    ):
        raise ScientificValidationError("natural frequency must be positive")
    gamma = damping_ratio * natural_frequency_rad_per_s
    damped_frequency = natural_frequency_rad_per_s * math.sqrt(1.0 - damping_ratio**2)
    return (
        complex(-gamma, damped_frequency),
        complex(-gamma, -damped_frequency),
        gamma,
        damped_frequency,
    )


def _heat_condition(distance: float, amplitude: float) -> RecoveryCondition:
    length_m = 1.0
    reference_frequency = 1.0
    diffusivity = distance * length_m**2 * reference_frequency
    predicted_rate = heat_mode_decay_rate(diffusivity, MODE_NUMBER, length_m)
    duration = RECOVERY_E_FOLDS / predicted_rate
    times = np.linspace(0.0, duration, OBSERVATION_SAMPLES)
    model = HeatDiffusionModel(
        length_m=length_m,
        left_temperature_k=300.0,
        right_temperature_k=300.0,
        diffusivity_m2_per_s=diffusivity,
        modes=(SpatialMode(mode_number=MODE_NUMBER, amplitude=amplitude),),
    )
    evaluation = evaluate_heat_diffusion(
        model,
        tuple(float(value) for value in np.linspace(0.0, length_m, 101)),
        tuple(float(value) for value in times),
    )
    modal_values = np.asarray(evaluation.modal_amplitudes, dtype=np.float64)[:, 0]
    measured_rate, r_squared = _fit_log_rate(times, modal_values)
    error = abs(measured_rate - predicted_rate) / predicted_rate
    return RecoveryCondition(
        model="fixed_end_heat_mode_n1",
        dimensionless_distance=distance,
        physical_parameter_value=diffusivity,
        perturbation_amplitude=amplitude,
        solver_rtol=None,
        solver_atol=None,
        output_sample_count=OBSERVATION_SAMPLES,
        equilibrium="uniform steady profile T_base(x)=300 K for 300 K endpoints",
        stability="asymptotically stable for alpha>0; non-attracting at alpha=0",
        predicted_rate_per_s=predicted_rate,
        measured_rate_per_s=measured_rate,
        relative_rate_error=error,
        fit_r_squared=r_squared,
        measurement_points=OBSERVATION_SAMPLES,
        measurement_window_s=duration,
        numerical_convergence_difference=None,
        local_approximation_accepted=error <= LOCAL_RATE_RELATIVE_ERROR_LIMIT,
        diagnostics="sampled exact single-mode evaluator; no transient PDE solver used",
    )


def _pendulum_condition(
    distance: float,
    amplitude: float,
    solver_rtol: float,
    solver_atol: float,
) -> RecoveryCondition:
    mass_kg = 1.0
    length_m = 1.0
    gravity_m_per_s2 = 1.0
    natural_frequency = math.sqrt(gravity_m_per_s2 / length_m)
    damping_coefficient = 2.0 * distance * mass_kg * length_m**2 * natural_frequency
    _eigenvalue_1, _eigenvalue_2, predicted_rate, damped_frequency = (
        pendulum_linearized_spectrum(
            distance, natural_frequency_rad_per_s=natural_frequency
        )
    )
    duration = RECOVERY_E_FOLDS / predicted_rate
    period = 2.0 * math.pi / damped_frequency
    sample_interval = min(0.01, period / 100.0)
    sample_count = int(math.ceil(duration / sample_interval)) + 1
    trace = simulate_damped_pendulum(
        mass_kg,
        length_m,
        damping_coefficient,
        gravity_m_per_s2,
        amplitude,
        0.0,
        duration,
        num_points=sample_count,
        method="DOP853",
        relative_tolerance=solver_rtol,
        absolute_tolerance=solver_atol,
    )
    peak_indices, _properties = find_peaks(trace.angle_rad)
    peak_times = trace.time_s[peak_indices]
    peak_amplitudes = trace.angle_rad[peak_indices]
    if trace.angle_rad[0] > trace.angle_rad[1]:
        peak_times = np.insert(peak_times, 0, trace.time_s[0])
        peak_amplitudes = np.insert(peak_amplitudes, 0, trace.angle_rad[0])
    keep = peak_amplitudes >= amplitude * 1e-6
    peak_times = peak_times[keep]
    peak_amplitudes = peak_amplitudes[keep]
    measured_rate, r_squared = _fit_log_rate(peak_times, peak_amplitudes)
    error = abs(measured_rate - predicted_rate) / predicted_rate
    return RecoveryCondition(
        model="nonlinear_damped_pendulum",
        dimensionless_distance=distance,
        physical_parameter_value=damping_coefficient,
        perturbation_amplitude=amplitude,
        solver_rtol=solver_rtol,
        solver_atol=solver_atol,
        output_sample_count=sample_count,
        equilibrium="downward rest state (theta, omega)=(0 rad, 0 rad/s)",
        stability="locally asymptotically stable for b>0; center at b=0",
        predicted_rate_per_s=predicted_rate,
        measured_rate_per_s=measured_rate,
        relative_rate_error=error,
        fit_r_squared=r_squared,
        measurement_points=int(peak_times.size),
        measurement_window_s=duration,
        numerical_convergence_difference=None,
        local_approximation_accepted=error <= LOCAL_RATE_RELATIVE_ERROR_LIMIT,
        diagnostics=(
            "log-linear fit of positive angular peaks; peak floor is 1e-6 "
            "times the initial angle"
        ),
    )


def _log_log_fit(distances: np.ndarray, rates: np.ndarray) -> tuple[float, float]:
    slope, intercept = np.polyfit(np.log(distances), np.log(rates), 1)
    fitted = slope * np.log(distances) + intercept
    residual = float(np.sum((np.log(rates) - fitted) ** 2))
    total = float(np.sum((np.log(rates) - np.mean(np.log(rates))) ** 2))
    return float(slope), float(1.0 - residual / total)


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _summarize_model(
    model_name: str, conditions: list[RecoveryCondition]
) -> ModelSummary:
    primary = sorted(
        (
            item
            for item in conditions
            if item.perturbation_amplitude == PERTURBATION_AMPLITUDES[0]
            and (item.solver_rtol is None or item.solver_rtol == SOLVER_SETTINGS[0][0])
        ),
        key=lambda item: item.dimensionless_distance,
    )
    distances = np.array([item.dimensionless_distance for item in primary])
    rates = np.array([item.measured_rate_per_s for item in primary])
    exponent, r_squared = _log_log_fit(distances, rates)
    sensitivity_differences: list[float] = []
    for distance in DISTANCES:
        amp_rates = [
            item.measured_rate_per_s
            for item in conditions
            if item.dimensionless_distance == distance
            and item.perturbation_amplitude == PERTURBATION_AMPLITUDES[1]
            and (item.solver_rtol is None or item.solver_rtol == SOLVER_SETTINGS[0][0])
        ]
        small_rates = [
            item.measured_rate_per_s
            for item in conditions
            if item.dimensionless_distance == distance
            and item.perturbation_amplitude == PERTURBATION_AMPLITUDES[0]
            and (item.solver_rtol is None or item.solver_rtol == SOLVER_SETTINGS[0][0])
        ]
        sensitivity_differences.append(
            abs(amp_rates[0] - small_rates[0]) / small_rates[0]
        )
    max_error = max(item.relative_rate_error for item in primary)
    slope_positive = bool(np.all(np.diff(rates) > 0.0))
    exponent_in_range = EXPONENT_ACCEPTANCE[0] <= exponent <= EXPONENT_ACCEPTANCE[1]
    r_squared_accepted = r_squared >= EXPONENT_R_SQUARED_MINIMUM
    primary_rates_accepted = all(item.local_approximation_accepted for item in primary)
    amplitude_accepted = max(sensitivity_differences) <= AMPLITUDE_SENSITIVITY_LIMIT
    return ModelSummary(
        model=model_name,
        fitted_power_law_exponent=exponent,
        exponent_r_squared=r_squared,
        predeclared_distance_min=float(distances[0]),
        predeclared_distance_max=float(distances[-1]),
        maximum_primary_rate_relative_error=max_error,
        maximum_amplitude_sensitivity_relative_difference=max(sensitivity_differences),
        all_primary_rates_within_tolerance=primary_rates_accepted,
        amplitude_sensitivity_within_tolerance=amplitude_accepted,
        exponent_within_predeclared_interval=exponent_in_range,
        exponent_r_squared_accepted=r_squared_accepted,
        qualitative_slowing_observed=slope_positive,
        model_prediction_supported=(
            exponent_in_range
            and r_squared_accepted
            and primary_rates_accepted
            and amplitude_accepted
            and slope_positive
        ),
    )


def _make_plots(
    output: Path,
    conditions: list[RecoveryCondition],
    summaries: list[ModelSummary],
) -> None:
    names = ("fixed_end_heat_mode_n1", "nonlinear_damped_pendulum")
    colors = ("#2463a6", "#c8582b")
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for name, color in zip(names, colors, strict=True):
        selected = sorted(
            (
                row
                for row in conditions
                if row.model == name
                and row.perturbation_amplitude == PERTURBATION_AMPLITUDES[0]
                and (
                    row.solver_rtol is None or row.solver_rtol == SOLVER_SETTINGS[0][0]
                )
            ),
            key=lambda row: row.dimensionless_distance,
        )
        ax.loglog(
            [row.dimensionless_distance for row in selected],
            [row.measured_rate_per_s for row in selected],
            marker="o",
            color=color,
            label=f"{name}: measured",
        )
        ax.loglog(
            [row.dimensionless_distance for row in selected],
            [row.predicted_rate_per_s for row in selected],
            linestyle="--",
            color=color,
            alpha=0.75,
            label=f"{name}: linearized prediction",
        )
    ax.set_xlabel("Dimensionless distance from zero-dissipation boundary")
    ax.set_ylabel("Recovery rate (s⁻¹)")
    ax.set_title("Local recovery slows as dissipation approaches zero")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "recovery_rate_vs_boundary_distance.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, name, color in zip(axes, names, colors, strict=True):
        selected = [
            row
            for row in conditions
            if row.model == name
            and row.perturbation_amplitude == PERTURBATION_AMPLITUDES[0]
            and (row.solver_rtol is None or row.solver_rtol == SOLVER_SETTINGS[0][0])
        ]
        ax.scatter(
            [row.predicted_rate_per_s for row in selected],
            [row.measured_rate_per_s for row in selected],
            color=color,
        )
        limits = [
            min(row.predicted_rate_per_s for row in selected),
            max(row.predicted_rate_per_s for row in selected),
        ]
        ax.plot(limits, limits, "k--", label="identity")
        ax.set_title(name)
        ax.set_xlabel("Linearized predicted rate (s⁻¹)")
        ax.set_ylabel("Measured recovery rate (s⁻¹)")
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.tight_layout()
    fig.savefig(output / "predicted_vs_measured_rates.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 5))
    for summary, color in zip(summaries, colors, strict=True):
        selected = sorted(
            (
                row
                for row in conditions
                if row.model == summary.model
                and row.perturbation_amplitude == PERTURBATION_AMPLITUDES[0]
                and (
                    row.solver_rtol is None or row.solver_rtol == SOLVER_SETTINGS[0][0]
                )
            ),
            key=lambda row: row.dimensionless_distance,
        )
        ax.loglog(
            [row.dimensionless_distance for row in selected],
            [row.measured_rate_per_s for row in selected],
            marker="o",
            color=color,
            label=(
                f"{summary.model}, fitted exponent="
                f"{summary.fitted_power_law_exponent:.3f}"
            ),
        )
        anchor = selected[2].measured_rate_per_s
        anchor_distance = selected[2].dimensionless_distance
        reference_x = np.array(
            [selected[0].dimensionless_distance, selected[-1].dimensionless_distance]
        )
        ax.loglog(
            reference_x,
            anchor * (reference_x / anchor_distance),
            linestyle="--",
            color=color,
            alpha=0.65,
            label=f"{summary.model}: exponent 1 reference",
        )
    ax.set_xlabel("Dimensionless distance from boundary")
    ax.set_ylabel("Measured recovery rate (s⁻¹)")
    ax.set_title("Predeclared log-log scaling fit across five distances")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(output / "recovery_rate_scaling.png", dpi=160)
    plt.close(fig)


def run_critical_slowing_down_study(output_directory: Path) -> dict[str, Any]:
    """Run the frozen Phase 49 study and write tables, plots, and metadata."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    conditions: list[RecoveryCondition] = []
    for distance in DISTANCES:
        for amplitude in PERTURBATION_AMPLITUDES:
            conditions.append(_heat_condition(distance, amplitude))
            for rtol, atol in SOLVER_SETTINGS:
                conditions.append(_pendulum_condition(distance, amplitude, rtol, atol))
    pendulum_pairs: dict[tuple[float, float], list[RecoveryCondition]] = {}
    for row in conditions:
        if row.model == "nonlinear_damped_pendulum":
            pendulum_pairs.setdefault(
                (row.dimensionless_distance, row.perturbation_amplitude), []
            ).append(row)
    conditions = [
        replace(
            row,
            numerical_convergence_difference=(
                abs(
                    pendulum_pairs[
                        (row.dimensionless_distance, row.perturbation_amplitude)
                    ][0].measured_rate_per_s
                    - pendulum_pairs[
                        (row.dimensionless_distance, row.perturbation_amplitude)
                    ][1].measured_rate_per_s
                )
                if row.model == "nonlinear_damped_pendulum"
                else None
            ),
        )
        for row in conditions
    ]
    summaries = [
        _summarize_model(
            "fixed_end_heat_mode_n1",
            [row for row in conditions if row.model == "fixed_end_heat_mode_n1"],
        ),
        _summarize_model(
            "nonlinear_damped_pendulum",
            [row for row in conditions if row.model == "nonlinear_damped_pendulum"],
        ),
    ]
    _write_rows(
        output / "per_condition_results.csv", [asdict(row) for row in conditions]
    )
    _write_rows(output / "model_summary.csv", [asdict(row) for row in summaries])
    _make_plots(output, conditions, summaries)
    metadata: dict[str, Any] = {
        "study": "Phase 49 critical slowing down and cross-system stability signatures",
        "models": {
            "heat": {
                "equation": "u_t=alpha*u_xx; u(0,t)=u(L,t)=0",
                "equilibrium": (
                    "uniform steady profile T_base=300 K, L=1 m, homogeneous fixed ends"
                ),
                "parameter": "dimensionless distance alpha/(L^2*omega_ref)",
                "boundary": "alpha=0; all profiles become stationary",
                "rate": "alpha*(n*pi/L)^2 for mode n=1",
                "measurement": (
                    "log-linear fit to sampled modal amplitudes from exact evaluator"
                ),
            },
            "pendulum": {
                "equation": "theta_ddot + b/(m*L^2)*theta_dot + (g/L)*sin(theta)=0",
                "equilibrium": "downward rest; m=L=g=1 SI",
                "parameter": "damping ratio b/(2*m*L^2*omega0)",
                "boundary": "zero damping; center with eigenvalues +/- i*omega0",
                "linearized_eigenvalues": "-zeta*omega0 +/- i*omega0*sqrt(1-zeta^2)",
                "rate": "envelope rate zeta*omega0 for underdamped zeta<1",
                "measurement": (
                    "log-linear fit to positive angular peaks from the existing "
                    "nonlinear solver"
                ),
            },
        },
        "dimensionless_distances": list(DISTANCES),
        "perturbation_amplitudes": list(PERTURBATION_AMPLITUDES),
        "solver": "existing nonlinear pendulum solve_ivp API with DOP853",
        "solver_settings": [
            {"rtol": rtol, "atol": atol} for rtol, atol in SOLVER_SETTINGS
        ],
        "heat_observation_samples": OBSERVATION_SAMPLES,
        "recovery_window_e_folds": RECOVERY_E_FOLDS,
        "primary_exponent_fit": (
            "OLS log(measured rate) vs log(dimensionless distance), all five "
            "frozen values, small perturbation"
        ),
        "criteria": {
            "primary_rate_relative_error_max": LOCAL_RATE_RELATIVE_ERROR_LIMIT,
            "amplitude_sensitivity_relative_difference_max": (
                AMPLITUDE_SENSITIVITY_LIMIT
            ),
            "power_law_exponent_interval": list(EXPONENT_ACCEPTANCE),
            "power_law_fit_r_squared_minimum": EXPONENT_R_SQUARED_MINIMUM,
            "qualitative_slowing_rule": (
                "measured rate increases at each adjacent distance"
            ),
        },
        "random_seed": None,
        "deterministic": True,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "conditions": len(conditions),
        "model_summaries": [asdict(summary) for summary in summaries],
    }
    (output / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return metadata
