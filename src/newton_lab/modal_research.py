"""Focused Phase 22 analysis of oscillator and heat modal evolution.

This research helper uses Newton Lab's existing oscillator and Phase 19 heat
APIs. It is a deterministic controlled comparison, not a generic discovery
framework or a claim that the physical systems are equivalent.
"""

from __future__ import annotations

from math import pi, sqrt

import numpy as np
from pydantic import Field

from newton_lab.dynamics import simulate_damped_oscillator
from newton_lab.exceptions import ScientificValidationError
from newton_lab.heat_diffusion_case_study import (
    HeatDiffusionSmoothingCaseStudy,
    SpatialMode,
    evaluate_heat_diffusion,
    run_heat_diffusion_smoothing_case_study,
    steady_temperature_profile,
)
from newton_lab.simulation import SimulationModel


class OscillatorRateResult(SimulationModel):
    """Analytical oscillator properties and a sampled envelope estimate."""

    mass_kg: float
    damping_coefficient_kg_per_s: float
    stiffness_n_per_m: float
    damping_ratio: float
    natural_frequency_rad_per_s: float
    envelope_decay_rate_per_s: float
    damped_frequency_rad_per_s: float
    sampled_peak_count: int = Field(ge=0)
    sampled_decay_rate_per_s: float | None
    maximum_analytical_displacement_error_m: float


class OscillatorSamplingResult(SimulationModel):
    """Effect of finite output interval or observation duration on rate fitting."""

    output_point_count: int = Field(ge=2)
    duration_s: float
    sample_interval_s: float
    absolute_peak_count: int = Field(ge=0)
    fitted_decay_rate_per_s: float | None


class HeatModeResult(SimulationModel):
    """Analytical modal attenuation recorded at the Phase 19 output times."""

    mode_number: int = Field(ge=1)
    initial_amplitude_k: float
    decay_rate_per_s: float
    times_s: tuple[float, ...]
    attenuation_factors: tuple[float, ...]
    modal_amplitudes_k: tuple[float, ...]


class HeatRateSensitivity(SimulationModel):
    """One-factor-at-a-time analytical heat modal-rate calculation."""

    varied_parameter: str
    values: tuple[float, ...]
    units: str
    decay_rates_per_s: tuple[float, ...]


class ModalDynamicsResearchResult(SimulationModel):
    """Reproducible numerical evidence and analytical checks for Phase 22."""

    heat_case_study_id: str
    steady_bvp_numerically_accepted: bool
    steady_bvp_max_profile_error_k: float
    heat_length_m: float
    heat_diffusivity_m2_per_s: float
    heat_times_s: tuple[float, ...]
    heat_modes: tuple[HeatModeResult, ...]
    heat_sensitivities: tuple[HeatRateSensitivity, ...]
    heat_slow_mode_fit_per_s: float
    heat_fast_mode_fit_above_resolution_per_s: float
    heat_fast_mode_at_final_time_k: float
    declared_observation_floor_k: float
    combined_signal_position_m: float
    combined_signal_times_s: tuple[float, ...]
    combined_signal_perturbation_k: tuple[float, ...]
    oscillator_rates: tuple[OscillatorRateResult, ...]
    oscillator_dense_sampling: OscillatorSamplingResult
    oscillator_sparse_sampling: OscillatorSamplingResult
    oscillator_short_window: OscillatorSamplingResult
    assumptions: tuple[str, ...] = Field(min_length=1)
    limitations: tuple[str, ...] = Field(min_length=1)


def heat_mode_decay_rate(
    diffusivity_m2_per_s: float, mode_number: int, length_m: float
) -> float:
    """Return ``alpha*(n*pi/L)^2`` for homogeneous fixed-end heat modes."""
    if not np.isfinite(diffusivity_m2_per_s) or diffusivity_m2_per_s <= 0.0:
        raise ScientificValidationError("diffusivity must be finite and positive")
    if (
        isinstance(mode_number, bool)
        or not isinstance(mode_number, int)
        or mode_number < 1
    ):
        raise ScientificValidationError("mode_number must be a positive integer")
    if not np.isfinite(length_m) or length_m <= 0.0:
        raise ScientificValidationError("length_m must be finite and positive")
    return diffusivity_m2_per_s * (mode_number * pi / length_m) ** 2


def oscillator_modal_properties(
    mass_kg: float, damping_coefficient_kg_per_s: float, stiffness_n_per_m: float
) -> tuple[float, float, float, float]:
    """Return ``(zeta, omega_n, gamma, omega_d)`` for an underdamped oscillator."""
    values = (mass_kg, damping_coefficient_kg_per_s, stiffness_n_per_m)
    if not np.isfinite(values).all() or mass_kg <= 0.0 or stiffness_n_per_m <= 0.0:
        raise ScientificValidationError(
            "mass and stiffness must be finite and positive"
        )
    if damping_coefficient_kg_per_s < 0.0:
        raise ScientificValidationError("damping must be finite and non-negative")
    omega_n = sqrt(stiffness_n_per_m / mass_kg)
    zeta = damping_coefficient_kg_per_s / (2.0 * sqrt(mass_kg * stiffness_n_per_m))
    if zeta >= 1.0:
        raise ScientificValidationError(
            "oscillator_modal_properties requires an underdamped oscillator"
        )
    gamma = damping_coefficient_kg_per_s / (2.0 * mass_kg)
    omega_d = omega_n * sqrt(1.0 - zeta**2)
    return zeta, omega_n, gamma, omega_d


def _estimate_log_decay_rate(
    times_s: tuple[float, ...], amplitudes: tuple[float, ...]
) -> float:
    """Fit a log slope only to finite, strictly positive modal amplitudes."""
    times = np.asarray(times_s, dtype=np.float64)
    values = np.asarray(amplitudes, dtype=np.float64)
    if (
        times.ndim != 1
        or values.ndim != 1
        or len(times) != len(values)
        or len(times) < 3
        or not np.isfinite(times).all()
        or not np.isfinite(values).all()
        or np.any(values <= 0.0)
        or np.any(np.diff(times) <= 0.0)
    ):
        raise ScientificValidationError(
            "log-rate fitting requires at least three ordered times and "
            "strictly positive finite amplitudes"
        )
    slope = float(np.polyfit(times, np.log(values), 1)[0])
    return -slope


def _absolute_peaks(
    time_s: np.ndarray, values: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    magnitudes = np.abs(values)
    indices = (
        np.flatnonzero(
            (magnitudes[1:-1] >= magnitudes[:-2]) & (magnitudes[1:-1] >= magnitudes[2:])
        )
        + 1
    )
    if len(magnitudes) > 1 and magnitudes[0] >= magnitudes[1]:
        indices = np.insert(indices, 0, 0)
    if len(magnitudes) > 1 and magnitudes[-1] >= magnitudes[-2]:
        indices = np.append(indices, len(magnitudes) - 1)
    indices = np.unique(indices)
    keep = magnitudes[indices] > 0.0
    return time_s[indices[keep]], magnitudes[indices[keep]]


def _sampled_rate(
    mass_kg: float,
    damping_kg_per_s: float,
    stiffness_n_per_m: float,
    duration_s: float,
    point_count: int,
) -> OscillatorSamplingResult:
    trace = simulate_damped_oscillator(
        mass_kg,
        damping_kg_per_s,
        stiffness_n_per_m,
        0.1,
        0.0,
        duration_s,
        num_points=point_count,
        method="RK45",
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
    )
    peak_times, peak_values = _absolute_peaks(trace.time_s, trace.displacement_m)
    rate = (
        _estimate_log_decay_rate(
            tuple(float(value) for value in peak_times),
            tuple(float(value) for value in peak_values),
        )
        if len(peak_values) >= 3
        else None
    )
    return OscillatorSamplingResult(
        output_point_count=point_count,
        duration_s=duration_s,
        sample_interval_s=duration_s / (point_count - 1),
        absolute_peak_count=len(peak_values),
        fitted_decay_rate_per_s=rate,
    )


def _oscillator_rate_result(
    mass_kg: float,
    damping_kg_per_s: float,
    stiffness_n_per_m: float,
    *,
    duration_s: float,
    point_count: int,
) -> OscillatorRateResult:
    zeta, omega_n, gamma, omega_d = oscillator_modal_properties(
        mass_kg, damping_kg_per_s, stiffness_n_per_m
    )
    trace = simulate_damped_oscillator(
        mass_kg,
        damping_kg_per_s,
        stiffness_n_per_m,
        0.1,
        0.0,
        duration_s,
        num_points=point_count,
        method="RK45",
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
    )
    peak_times, peak_values = _absolute_peaks(trace.time_s, trace.displacement_m)
    fitted = (
        _estimate_log_decay_rate(
            tuple(float(value) for value in peak_times),
            tuple(float(value) for value in peak_values),
        )
        if len(peak_values) >= 3
        else None
    )
    damped_coefficient = gamma / omega_d
    analytic = (
        0.1
        * np.exp(-gamma * trace.time_s)
        * (
            np.cos(omega_d * trace.time_s)
            + damped_coefficient * np.sin(omega_d * trace.time_s)
        )
    )
    return OscillatorRateResult(
        mass_kg=mass_kg,
        damping_coefficient_kg_per_s=damping_kg_per_s,
        stiffness_n_per_m=stiffness_n_per_m,
        damping_ratio=zeta,
        natural_frequency_rad_per_s=omega_n,
        envelope_decay_rate_per_s=gamma,
        damped_frequency_rad_per_s=omega_d,
        sampled_peak_count=len(peak_values),
        sampled_decay_rate_per_s=fitted,
        maximum_analytical_displacement_error_m=float(
            np.max(np.abs(trace.displacement_m - analytic))
        ),
    )


def run_modal_dynamics_research(
    heat_case: HeatDiffusionSmoothingCaseStudy | None = None,
) -> ModalDynamicsResearchResult:
    """Run deterministic Phase 22 checks using existing physical model APIs.

    If no Phase 19 result is supplied, its established case-study API is run
    once to provide the declared BVP provenance and analytical modal outputs.
    The oscillator checks use the existing SI oscillator solver. No market or
    external data are used.
    """
    case = heat_case or run_heat_diffusion_smoothing_case_study(sample_count=101)
    model = case.diffusion_model
    times = case.evaluation.times_s
    modes = tuple(
        HeatModeResult(
            mode_number=mode.mode_number,
            initial_amplitude_k=mode.amplitude,
            decay_rate_per_s=heat_mode_decay_rate(
                model.diffusivity_m2_per_s, mode.mode_number, model.length_m
            ),
            times_s=times,
            attenuation_factors=tuple(
                row[index] for row in case.evaluation.attenuation_factors
            ),
            modal_amplitudes_k=tuple(
                row[index] for row in case.evaluation.modal_amplitudes
            ),
        )
        for index, mode in enumerate(model.modes)
    )
    sensitivity_rows = (
        HeatRateSensitivity(
            varied_parameter="diffusivity",
            values=(
                0.5 * model.diffusivity_m2_per_s,
                model.diffusivity_m2_per_s,
                2.0 * model.diffusivity_m2_per_s,
            ),
            units="m^2/s",
            decay_rates_per_s=tuple(
                heat_mode_decay_rate(value, 1, model.length_m)
                for value in (
                    0.5 * model.diffusivity_m2_per_s,
                    model.diffusivity_m2_per_s,
                    2.0 * model.diffusivity_m2_per_s,
                )
            ),
        ),
        HeatRateSensitivity(
            varied_parameter="mode_number",
            values=(1.0, 2.0, 4.0),
            units="dimensionless",
            decay_rates_per_s=tuple(
                heat_mode_decay_rate(model.diffusivity_m2_per_s, n, model.length_m)
                for n in (1, 2, 4)
            ),
        ),
        HeatRateSensitivity(
            varied_parameter="domain_length",
            values=(0.5 * model.length_m, model.length_m, 2.0 * model.length_m),
            units="m",
            decay_rates_per_s=tuple(
                heat_mode_decay_rate(model.diffusivity_m2_per_s, 1, length)
                for length in (
                    0.5 * model.length_m,
                    model.length_m,
                    2.0 * model.length_m,
                )
            ),
        ),
    )

    floor_k = 1e-8
    slow = modes[0]
    fast = max(modes, key=lambda item: item.mode_number)
    slow_fit = _estimate_log_decay_rate(
        slow.times_s, tuple(abs(v) for v in slow.modal_amplitudes_k)
    )
    visible = tuple(
        (time, abs(amplitude))
        for time, amplitude in zip(fast.times_s, fast.modal_amplitudes_k, strict=True)
        if abs(amplitude) >= floor_k
    )
    fast_fit = _estimate_log_decay_rate(
        tuple(item[0] for item in visible), tuple(item[1] for item in visible)
    )

    combined_model = model.model_copy(
        update={
            "modes": (
                SpatialMode(mode_number=1, amplitude=0.5),
                SpatialMode(mode_number=2, amplitude=-1.0),
            )
        }
    )
    combined_times = (0.0, 100.0, 300.0, 400.0, 600.0)
    position = 0.25 * model.length_m
    combined_eval = evaluate_heat_diffusion(
        combined_model, (position, 0.75 * model.length_m), combined_times
    )
    base = steady_temperature_profile(combined_model, (position,))[0]
    combined_perturbation = tuple(row[0] - base for row in combined_eval.temperature_k)

    oscillator_rates = tuple(
        _oscillator_rate_result(1.0, damping, 4.0, duration_s=10.0, point_count=20001)
        for damping in (0.05, 0.1, 0.2)
    ) + tuple(
        _oscillator_rate_result(mass, 0.1, 4.0, duration_s=10.0, point_count=20001)
        for mass in (0.5, 2.0)
    )
    return ModalDynamicsResearchResult(
        heat_case_study_id=case.case_study_id,
        steady_bvp_numerically_accepted=case.source_bvp_numerically_accepted,
        steady_bvp_max_profile_error_k=case.source_bvp_max_profile_error_k,
        heat_length_m=model.length_m,
        heat_diffusivity_m2_per_s=model.diffusivity_m2_per_s,
        heat_times_s=times,
        heat_modes=modes,
        heat_sensitivities=sensitivity_rows,
        heat_slow_mode_fit_per_s=slow_fit,
        heat_fast_mode_fit_above_resolution_per_s=fast_fit,
        heat_fast_mode_at_final_time_k=abs(fast.modal_amplitudes_k[-1]),
        declared_observation_floor_k=floor_k,
        combined_signal_position_m=position,
        combined_signal_times_s=combined_times,
        combined_signal_perturbation_k=combined_perturbation,
        oscillator_rates=oscillator_rates,
        oscillator_dense_sampling=_sampled_rate(1.0, 0.1, 4.0, 10.0, 20001),
        oscillator_sparse_sampling=_sampled_rate(1.0, 0.1, 4.0, 10.0, 13),
        oscillator_short_window=_sampled_rate(1.0, 0.1, 4.0, 1.5, 3001),
        assumptions=(
            "Heat: homogeneous 1D constant-property medium; fixed endpoints; "
            "no source; homogeneous Dirichlet perturbation modes.",
            "Oscillator: linear spring and viscous damping, constant positive "
            "mass/stiffness, no force, and underdamped parameters.",
            "Oscillator rate fits use log absolute local peaks, never signed "
            "displacement samples or zero crossings.",
            "The 1e-8 K heat floor is an analysis threshold, not measured "
            "instrument resolution.",
        ),
        limitations=(
            "Heat modal attenuation and spatial samples use the exact finite-mode "
            "formula; the steady BVP is a separate numerical experiment.",
            "Oscillator displacement is numerically integrated and sampled; its "
            "error is measured against the analytical solution.",
            "Envelope fits require resolved peaks and a sufficient finite window; "
            "sparse output can bias peak locations or prevent fitting.",
            "This parameter grid is illustrative, not universal coverage.",
            "The comparison shows conditional structure, not physical equivalence "
            "or application and financial validity.",
        ),
    )


def render_modal_dynamics_report(result: ModalDynamicsResearchResult) -> str:
    """Render compact deterministic tables of computed Phase 22 results."""
    lines = [
        "# Phase 22 computed modal results",
        "",
        f"Heat case `{result.heat_case_study_id}`: L={result.heat_length_m:g} m; "
        f"alpha={result.heat_diffusivity_m2_per_s:.6g} m^2/s; "
        f"sample times={result.heat_times_s} s. Steady BVP acceptance="
        f"{result.steady_bvp_numerically_accepted}; max profile error="
        f"{result.steady_bvp_max_profile_error_k:.6g} K.",
        "",
        "| Mode | Initial (K) | Rate (s^-1) | Attenuation factors |",
        "|---:|---:|---:|---|",
    ]
    for mode in result.heat_modes:
        factors = ", ".join(f"{value:.6g}" for value in mode.attenuation_factors)
        lines.append(
            f"| n={mode.mode_number} | {mode.initial_amplitude_k:.6g} | "
            f"{mode.decay_rate_per_s:.9g} | {factors} |"
        )
    lines.extend(
        [
            "",
            "| Heat rate factor | Values | Rates (s^-1) |",
            "|---|---|---|",
        ]
    )
    for row in result.heat_sensitivities:
        lines.append(
            f"| {row.varied_parameter} ({row.units}) | {row.values} | "
            f"{tuple(round(value, 9) for value in row.decay_rates_per_s)} |"
        )
    lines.extend(
        [
            "",
            "| m (kg) | c (kg/s) | zeta | gamma (s^-1) | "
            "omega_d (rad/s) | peak fit (s^-1) | max error (m) |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for item in result.oscillator_rates:
        fitted = (
            "unavailable"
            if item.sampled_decay_rate_per_s is None
            else f"{item.sampled_decay_rate_per_s:.9g}"
        )
        lines.append(
            f"| {item.mass_kg:g} | {item.damping_coefficient_kg_per_s:g} | "
            f"{item.damping_ratio:.6g} | {item.envelope_decay_rate_per_s:.6g} | "
            f"{item.damped_frequency_rad_per_s:.9g} | {fitted} | "
            f"{item.maximum_analytical_displacement_error_m:.3g} |"
        )
    for label, sample in (
        ("dense", result.oscillator_dense_sampling),
        ("sparse", result.oscillator_sparse_sampling),
        ("short window", result.oscillator_short_window),
    ):
        lines.append(
            f"- {label}: points={sample.output_point_count}, "
            f"dt={sample.sample_interval_s:.6g} s, "
            f"peaks={sample.absolute_peak_count}, "
            f"fitted rate={sample.fitted_decay_rate_per_s}."
        )
    combined_values = tuple(
        round(value, 6) for value in result.combined_signal_perturbation_k
    )
    lines.extend(
        [
            "",
            f"Heat combined perturbation at x={result.combined_signal_position_m:g} m "
            f"and times {result.combined_signal_times_s} s: "
            f"{combined_values} K.",
            "",
        ]
    )
    return "\n".join(lines)
