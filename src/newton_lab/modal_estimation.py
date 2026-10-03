"""Phase 23 synthetic-observation experiments for estimating modal decay rates.

Observation records intentionally contain no ground-truth rate. Generators
return truth separately, and estimator functions accept observation records
only. This keeps evaluation post hoc and prevents parameter leakage.
"""

from __future__ import annotations

from enum import StrEnum
from math import isfinite, pi, sin, sqrt
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from pydantic import Field, field_validator, model_validator
from scipy.optimize import least_squares  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.modal_research import heat_mode_decay_rate, oscillator_modal_properties
from newton_lab.simulation import NamedValue, SimulationModel


class EstimateStatus(StrEnum):
    """Interpretation statuses kept separate from optimizer convergence."""

    INTERPRETABLE = "available_interpretable_under_tested_assumptions"
    FIT_FAILED = "fit_failed"
    INSUFFICIENT = "insufficient_observations"
    NOT_IDENTIFIABLE = "not_identifiable_at_declared_observation_resolution"
    UNRELIABLE = "available_but_unreliable_under_declared_diagnostics"


class ObservationSeries(SimulationModel):
    """Observed values and observation-process metadata, with no hidden truth."""

    scenario_id: str = Field(min_length=1)
    observable: str = Field(min_length=1)
    times: tuple[float, ...] = Field(min_length=1)
    values: tuple[float, ...] = Field(min_length=1)
    time_unit: str = Field(min_length=1)
    value_unit: str = Field(min_length=1)
    noise_model: Literal["none", "iid_gaussian_additive"]
    noise_standard_deviation: float = Field(ge=0.0)
    noise_seed: int | None = Field(default=None, ge=0)

    @model_validator(mode="after")
    def aligned_finite_samples(self) -> ObservationSeries:
        if len(self.times) != len(self.values):
            raise ValueError("observation times and values must have equal length")
        if not np.isfinite(self.times).all() or not np.isfinite(self.values).all():
            raise ValueError("observations must contain finite times and values")
        if len(self.times) > 1 and not np.all(np.diff(self.times) > 0.0):
            raise ValueError("observation times must be strictly increasing")
        if self.noise_model == "iid_gaussian_additive" and self.noise_seed is None:
            raise ValueError("Gaussian synthetic noise must record its fixed seed")
        if self.noise_model == "none" and self.noise_standard_deviation != 0.0:
            raise ValueError("noise-free observations must declare zero noise")
        return self

    @field_validator("noise_standard_deviation")
    @classmethod
    def finite_noise(cls, value: float) -> float:
        if not isfinite(value):
            raise ValueError("noise standard deviation must be finite")
        return value


class SyntheticGroundTruth(SimulationModel):
    """Generator-side truth stored separately from estimator inputs."""

    scenario_id: str
    model_description: str
    component_rates_per_s: tuple[float, ...] = Field(min_length=1)
    component_amplitudes: tuple[float, ...] = Field(min_length=1)
    amplitude_unit: str
    target_rate_per_s: float | None
    model_parameters: tuple[NamedValue, ...] = ()
    oscillation_frequency_rad_per_s: float | None = None

    @model_validator(mode="after")
    def aligned_components(self) -> SyntheticGroundTruth:
        if len(self.component_rates_per_s) != len(self.component_amplitudes):
            raise ValueError("ground-truth component rates and amplitudes must align")
        if any(
            rate <= 0.0 or not isfinite(rate) for rate in self.component_rates_per_s
        ):
            raise ValueError("ground-truth decay rates must be finite and positive")
        if not np.isfinite(self.component_amplitudes).all():
            raise ValueError("ground-truth amplitudes must be finite")
        return self


class SyntheticObservationCase(SimulationModel):
    """Synthetic observations paired with separately held generator truth."""

    observations: ObservationSeries
    ground_truth: SyntheticGroundTruth

    @model_validator(mode="after")
    def matching_identity(self) -> SyntheticObservationCase:
        if self.observations.scenario_id != self.ground_truth.scenario_id:
            raise ValueError("observation and truth scenario IDs must match")
        return self


class ExcludedObservation(SimulationModel):
    """One sample omitted from an estimator, with its explicit reason."""

    sample_index: int = Field(ge=0)
    reason: str = Field(min_length=1)


class DecayRateEstimate(SimulationModel):
    """Rate fit, convergence evidence, exclusions, and residual diagnostics."""

    status: EstimateStatus
    method: str
    estimated_rate_per_s: float | None
    estimated_amplitude: float | None
    fit_converged: bool
    convergence_message: str
    input_sample_count: int = Field(ge=0)
    used_sample_indices: tuple[int, ...] = ()
    excluded_observations: tuple[ExcludedObservation, ...] = ()
    fitted_values: tuple[float, ...] = ()
    observation_rmse: float | None
    observation_rmse_unit: str
    log_rmse: float | None
    local_rate_noise_sensitivity_per_s: float | None
    jacobian_condition_number: float | None
    diagnostics: tuple[str, ...] = ()


class RateEvaluation(SimulationModel):
    """Post-fit comparison to truth; never passed into an estimator."""

    method: str
    target_rate_per_s: float | None
    absolute_error_per_s: float | None
    relative_error: float | None
    nearest_component_rate_per_s: float | None
    nearest_component_gap_per_s: float | None
    within_declared_10_percent_tolerance: bool | None
    interpretation: str


class EstimateOutcome(SimulationModel):
    """Estimator result and independent post hoc ground-truth evaluation."""

    estimate: DecayRateEstimate
    evaluation: RateEvaluation


class EstimationExperiment(SimulationModel):
    """One controlled condition, retaining observations and separate truth."""

    observations: ObservationSeries
    ground_truth: SyntheticGroundTruth
    condition_summary: str
    outcomes: tuple[EstimateOutcome, ...] = Field(min_length=1)


class Phase23ResearchResults(SimulationModel):
    """Reproducible cadence, duration, noise, oscillator, and mixture results."""

    experiments: tuple[EstimationExperiment, ...] = Field(min_length=1)
    seed_policy: str
    interpretation_policy: tuple[str, ...] = Field(min_length=1)


def _validated_times(times: tuple[float, ...]) -> tuple[float, ...]:
    values = tuple(float(value) for value in times)
    if not values or not np.isfinite(values).all():
        raise ScientificValidationError("times must be a non-empty finite sequence")
    if len(values) > 1 and not np.all(np.diff(values) > 0.0):
        raise ScientificValidationError("times must be strictly increasing")
    if values[0] < 0.0:
        raise ScientificValidationError("times must be non-negative")
    return values


def _synthetic_noise(
    values: np.ndarray, noise_standard_deviation: float, seed: int
) -> np.ndarray:
    if noise_standard_deviation == 0.0:
        return values.copy()
    return values + np.random.default_rng(seed).normal(
        0.0, noise_standard_deviation, size=len(values)
    )


def generate_exponential_observations(
    *,
    scenario_id: str,
    amplitude: float,
    rate_per_s: float,
    times_s: tuple[float, ...],
    value_unit: str,
    noise_standard_deviation: float = 0.0,
    noise_seed: int = 2301,
    observable: str = "scalar_exponential_decay",
) -> SyntheticObservationCase:
    """Generate ``A*exp(-lambda*t)`` observations and hold truth separately."""
    if not isfinite(amplitude) or amplitude <= 0.0:
        raise ScientificValidationError("amplitude must be finite and positive")
    if not isfinite(rate_per_s) or rate_per_s <= 0.0:
        raise ScientificValidationError("rate_per_s must be finite and positive")
    if not isfinite(noise_standard_deviation) or noise_standard_deviation < 0.0:
        raise ScientificValidationError("noise standard deviation must be non-negative")
    if isinstance(noise_seed, bool) or noise_seed < 0:
        raise ScientificValidationError("noise_seed must be a non-negative integer")
    times = _validated_times(times_s)
    exact_values = amplitude * np.exp(-rate_per_s * np.asarray(times))
    observed = _synthetic_noise(exact_values, noise_standard_deviation, noise_seed)
    observation = ObservationSeries(
        scenario_id=scenario_id,
        observable=observable,
        times=times,
        values=tuple(float(value) for value in observed),
        time_unit="s",
        value_unit=value_unit,
        noise_model=(
            "iid_gaussian_additive" if noise_standard_deviation > 0.0 else "none"
        ),
        noise_standard_deviation=noise_standard_deviation,
        noise_seed=(noise_seed if noise_standard_deviation > 0.0 else None),
    )
    truth = SyntheticGroundTruth(
        scenario_id=scenario_id,
        model_description="single exponential A*exp(-lambda*t)",
        component_rates_per_s=(rate_per_s,),
        component_amplitudes=(amplitude,),
        amplitude_unit=value_unit,
        target_rate_per_s=rate_per_s,
    )
    return SyntheticObservationCase(observations=observation, ground_truth=truth)


def generate_heat_mode_observations(
    *,
    scenario_id: str,
    diffusivity_m2_per_s: float,
    mode_number: int,
    length_m: float,
    initial_amplitude_k: float,
    times_s: tuple[float, ...],
    noise_standard_deviation_k: float = 0.0,
    noise_seed: int = 2302,
) -> SyntheticObservationCase:
    """Generate one supported Phase 19 heat mode at exact requested times."""
    rate = heat_mode_decay_rate(diffusivity_m2_per_s, mode_number, length_m)
    case = generate_exponential_observations(
        scenario_id=scenario_id,
        amplitude=initial_amplitude_k,
        rate_per_s=rate,
        times_s=times_s,
        value_unit="K",
        noise_standard_deviation=noise_standard_deviation_k,
        noise_seed=noise_seed,
        observable=f"heat_mode_n{mode_number}_temperature_amplitude",
    )
    truth = case.ground_truth.model_copy(
        update={
            "model_description": (
                "Phase 19 homogeneous fixed-end 1D heat mode; lambda=alpha*(n*pi/L)^2"
            ),
            "model_parameters": (
                NamedValue(
                    name="diffusivity",
                    value=diffusivity_m2_per_s,
                    unit="m^2/s",
                ),
                NamedValue(name="mode_number", value=float(mode_number), unit="1"),
                NamedValue(name="domain_length", value=length_m, unit="m"),
            ),
        }
    )
    return case.model_copy(update={"ground_truth": truth})


def generate_oscillator_observations(
    *,
    scenario_id: str,
    mass_kg: float,
    damping_coefficient_kg_per_s: float,
    stiffness_n_per_m: float,
    initial_displacement_m: float,
    initial_velocity_m_per_s: float,
    times_s: tuple[float, ...],
    noise_standard_deviation_m: float = 0.0,
    noise_seed: int = 2303,
) -> SyntheticObservationCase:
    """Generate exact underdamped Phase 18 displacement observations."""
    zeta, omega_n, gamma, omega_d = oscillator_modal_properties(
        mass_kg, damping_coefficient_kg_per_s, stiffness_n_per_m
    )
    if not np.isfinite((initial_displacement_m, initial_velocity_m_per_s)).all():
        raise ScientificValidationError("oscillator initial conditions must be finite")
    times = _validated_times(times_s)
    coefficient = (initial_velocity_m_per_s + gamma * initial_displacement_m) / omega_d
    time_values = np.asarray(times, dtype=np.float64)
    exact_values = np.exp(-gamma * time_values) * (
        initial_displacement_m * np.cos(omega_d * time_values)
        + coefficient * np.sin(omega_d * time_values)
    )
    if np.all(exact_values == 0.0):
        raise ScientificValidationError(
            "initial conditions produce no observable motion"
        )
    if not isfinite(noise_standard_deviation_m) or noise_standard_deviation_m < 0.0:
        raise ScientificValidationError("noise standard deviation must be non-negative")
    observed = _synthetic_noise(exact_values, noise_standard_deviation_m, noise_seed)
    observation = ObservationSeries(
        scenario_id=scenario_id,
        observable="signed_oscillator_displacement",
        times=times,
        values=tuple(float(value) for value in observed),
        time_unit="s",
        value_unit="m",
        noise_model=(
            "iid_gaussian_additive" if noise_standard_deviation_m > 0.0 else "none"
        ),
        noise_standard_deviation=noise_standard_deviation_m,
        noise_seed=(noise_seed if noise_standard_deviation_m > 0.0 else None),
    )
    truth = SyntheticGroundTruth(
        scenario_id=scenario_id,
        model_description=(
            "Phase 18 linear unforced underdamped oscillator analytical solution"
        ),
        component_rates_per_s=(gamma,),
        component_amplitudes=(sqrt(initial_displacement_m**2 + coefficient**2),),
        amplitude_unit="m",
        target_rate_per_s=gamma,
        model_parameters=(
            NamedValue(name="mass", value=mass_kg, unit="kg"),
            NamedValue(
                name="damping_coefficient",
                value=damping_coefficient_kg_per_s,
                unit="kg/s",
            ),
            NamedValue(name="stiffness", value=stiffness_n_per_m, unit="N/m"),
            NamedValue(
                name="initial_displacement",
                value=initial_displacement_m,
                unit="m",
            ),
            NamedValue(
                name="initial_velocity",
                value=initial_velocity_m_per_s,
                unit="m/s",
            ),
        ),
        oscillation_frequency_rad_per_s=omega_d,
    )
    return SyntheticObservationCase(observations=observation, ground_truth=truth)


def generate_mixed_heat_observations(
    *,
    scenario_id: str,
    mode_numbers: tuple[int, ...],
    effective_amplitudes_k: tuple[float, ...],
    diffusivity_m2_per_s: float,
    length_m: float,
    times_s: tuple[float, ...],
    position_fraction: float = 0.125,
    noise_standard_deviation_k: float = 0.0,
    noise_seed: int = 2304,
) -> SyntheticObservationCase:
    """Generate a fixed-position sum of Phase 19 spatial diffusion modes."""
    if len(mode_numbers) < 2 or len(mode_numbers) != len(effective_amplitudes_k):
        raise ScientificValidationError("mixed modes require aligned components")
    if not 0.0 < position_fraction < 1.0:
        raise ScientificValidationError("position_fraction must lie strictly in (0, 1)")
    if not isfinite(noise_standard_deviation_k) or noise_standard_deviation_k < 0.0:
        raise ScientificValidationError(
            "noise standard deviation must be finite and non-negative"
        )
    if (
        isinstance(noise_seed, bool)
        or not isinstance(noise_seed, int)
        or noise_seed < 0
    ):
        raise ScientificValidationError("noise_seed must be a non-negative integer")
    if any(
        isinstance(n, bool) or not isinstance(n, int) or n < 1 for n in mode_numbers
    ):
        raise ScientificValidationError("mode numbers must be positive integers")
    rates = tuple(
        heat_mode_decay_rate(diffusivity_m2_per_s, n, length_m) for n in mode_numbers
    )
    sine_factors = tuple(sin(n * pi * position_fraction) for n in mode_numbers)
    if any(abs(value) < 1e-12 for value in sine_factors):
        raise ScientificValidationError(
            "observation position suppresses a selected mode"
        )
    times = _validated_times(times_s)
    components = tuple(
        amplitude * np.exp(-rate * np.asarray(times))
        for amplitude, rate in zip(effective_amplitudes_k, rates, strict=True)
    )
    exact_values = np.sum(np.vstack(components), axis=0)
    values = _synthetic_noise(exact_values, noise_standard_deviation_k, noise_seed)
    observation = ObservationSeries(
        scenario_id=scenario_id,
        observable="mixed_heat_modal_temperature_perturbation",
        times=times,
        values=tuple(float(value) for value in values),
        time_unit="s",
        value_unit="K",
        noise_model=(
            "iid_gaussian_additive" if noise_standard_deviation_k > 0.0 else "none"
        ),
        noise_standard_deviation=noise_standard_deviation_k,
        noise_seed=noise_seed if noise_standard_deviation_k > 0.0 else None,
    )
    truth = SyntheticGroundTruth(
        scenario_id=scenario_id,
        model_description=(
            "sum of fixed-position Phase 19 sine modes with homogeneous "
            "fixed-end boundary conditions"
        ),
        component_rates_per_s=rates,
        component_amplitudes=effective_amplitudes_k,
        amplitude_unit="K",
        target_rate_per_s=None,
        model_parameters=(
            NamedValue(name="diffusivity", value=diffusivity_m2_per_s, unit="m^2/s"),
            NamedValue(name="domain_length", value=length_m, unit="m"),
            NamedValue(
                name="position_fraction",
                value=position_fraction,
                unit="dimensionless",
            ),
        ),
    )
    return SyntheticObservationCase(observations=observation, ground_truth=truth)


def _excluded_indices(
    indices: tuple[int, ...], reason: str
) -> tuple[ExcludedObservation, ...]:
    return tuple(
        ExcludedObservation(sample_index=index, reason=reason) for index in indices
    )


def _status_for_fit(
    *,
    rate: float,
    noise_sensitivity: float | None,
    condition_number: float | None,
    log_rmse: float | None,
    excluded_fraction: float,
) -> tuple[EstimateStatus, tuple[str, ...]]:
    diagnostics: list[str] = []
    if condition_number is not None and condition_number > 1e8:
        return EstimateStatus.NOT_IDENTIFIABLE, (
            "scaled fit Jacobian is ill-conditioned",
        )
    if noise_sensitivity is not None and abs(rate) <= 3.0 * noise_sensitivity:
        return EstimateStatus.NOT_IDENTIFIABLE, (
            "estimated rate does not exceed three times its local "
            "noise-sensitivity scale; "
            "this is a declared screening rule, not a confidence interval",
        )
    if rate <= 0.0:
        return EstimateStatus.UNRELIABLE, ("fitted rate is not positive decay",)
    if log_rmse is not None and log_rmse > 0.15:
        diagnostics.append(
            "log residual exceeds the declared 0.15 model-mismatch screen"
        )
    if excluded_fraction > 0.5:
        diagnostics.append("more than half of input samples were excluded")
    if diagnostics:
        return EstimateStatus.UNRELIABLE, tuple(diagnostics)
    return EstimateStatus.INTERPRETABLE, (
        "estimate passes the declared observation and residual screens only",
    )


def estimate_exponential_rate(
    observations: ObservationSeries,
    *,
    method: Literal["log_linear", "nonlinear"] = "log_linear",
) -> DecayRateEstimate:
    """Estimate scalar exponential decay using observations only.

    Log-linear fitting excludes nonpositive values and, for declared additive
    Gaussian noise, positive values no greater than three noise standard
    deviations. Nonlinear least squares fits all signed observations to
    ``A*exp(-rate*t)`` with positive amplitude and rate bounds.
    """
    times = np.asarray(observations.times, dtype=np.float64)
    values = np.asarray(observations.values, dtype=np.float64)
    count = len(times)
    if count < 3:
        return DecayRateEstimate(
            status=EstimateStatus.INSUFFICIENT,
            method=method,
            estimated_rate_per_s=None,
            estimated_amplitude=None,
            fit_converged=False,
            convergence_message="at least three observations are required",
            input_sample_count=count,
            observation_rmse=None,
            observation_rmse_unit=observations.value_unit,
            log_rmse=None,
            local_rate_noise_sensitivity_per_s=None,
            jacobian_condition_number=None,
        )
    if method == "nonlinear":
        return _estimate_nonlinear(observations)

    excluded: list[ExcludedObservation] = []
    included: list[int] = []
    floor = 3.0 * observations.noise_standard_deviation
    for index, value in enumerate(values):
        if value <= 0.0:
            excluded.append(
                ExcludedObservation(
                    sample_index=index,
                    reason="nonpositive value cannot be log-transformed",
                )
            )
        elif observations.noise_standard_deviation > 0.0 and value <= floor:
            excluded.append(
                ExcludedObservation(
                    sample_index=index,
                    reason="value is at or below the declared three-sigma noise floor",
                )
            )
        else:
            included.append(index)
    if len(included) < 3:
        return DecayRateEstimate(
            status=EstimateStatus.INSUFFICIENT,
            method="log_linear",
            estimated_rate_per_s=None,
            estimated_amplitude=None,
            fit_converged=False,
            convergence_message="fewer than three usable positive observations",
            input_sample_count=count,
            used_sample_indices=tuple(included),
            excluded_observations=tuple(excluded),
            observation_rmse=None,
            observation_rmse_unit=observations.value_unit,
            log_rmse=None,
            local_rate_noise_sensitivity_per_s=None,
            jacobian_condition_number=None,
            diagnostics=("insufficient positive above-floor samples",),
        )
    used_times = times[included]
    used_values = values[included]
    design = np.column_stack((np.ones(len(included)), -used_times))
    log_values = np.log(used_values)
    if observations.noise_standard_deviation > 0.0:
        log_sigma = observations.noise_standard_deviation / used_values
        weights = 1.0 / np.square(log_sigma)
    else:
        weights = np.ones(len(included), dtype=np.float64)
    weighted_design = design * np.sqrt(weights)[:, None]
    weighted_values = log_values * np.sqrt(weights)
    beta, _, rank, _ = np.linalg.lstsq(weighted_design, weighted_values, rcond=None)
    if rank < 2:
        return DecayRateEstimate(
            status=EstimateStatus.NOT_IDENTIFIABLE,
            method="log_linear",
            estimated_rate_per_s=None,
            estimated_amplitude=None,
            fit_converged=True,
            convergence_message="weighted design matrix has rank below two",
            input_sample_count=count,
            used_sample_indices=tuple(included),
            excluded_observations=tuple(excluded),
            observation_rmse=None,
            observation_rmse_unit=observations.value_unit,
            log_rmse=None,
            local_rate_noise_sensitivity_per_s=None,
            jacobian_condition_number=None,
            diagnostics=("rate and intercept are not separately determined",),
        )
    amplitude = float(np.exp(beta[0]))
    rate = float(beta[1])
    fitted = amplitude * np.exp(-rate * used_times)
    log_residual = log_values - (beta[0] - rate * used_times)
    log_rmse = float(np.sqrt(np.mean(np.square(log_residual))))
    rmse = float(np.sqrt(np.mean(np.square(used_values - fitted))))
    condition = float(np.linalg.cond(weighted_design))
    sensitivity: float | None = None
    if observations.noise_standard_deviation > 0.0:
        covariance = np.linalg.pinv(weighted_design.T @ weighted_design)
        sensitivity = float(sqrt(max(0.0, covariance[1, 1])))
    status, diagnostics = _status_for_fit(
        rate=rate,
        noise_sensitivity=sensitivity,
        condition_number=condition,
        log_rmse=log_rmse,
        excluded_fraction=len(excluded) / count,
    )
    return DecayRateEstimate(
        status=status,
        method="weighted_log_linear" if sensitivity is not None else "log_linear",
        estimated_rate_per_s=rate,
        estimated_amplitude=amplitude,
        fit_converged=True,
        convergence_message="weighted linear least squares converged",
        input_sample_count=count,
        used_sample_indices=tuple(included),
        excluded_observations=tuple(excluded),
        fitted_values=tuple(float(value) for value in fitted),
        observation_rmse=rmse,
        observation_rmse_unit=observations.value_unit,
        log_rmse=log_rmse,
        local_rate_noise_sensitivity_per_s=sensitivity,
        jacobian_condition_number=condition,
        diagnostics=diagnostics,
    )


def _estimate_nonlinear(observations: ObservationSeries) -> DecayRateEstimate:
    times = np.asarray(observations.times, dtype=np.float64)
    values = np.asarray(observations.values, dtype=np.float64)
    span = float(times[-1] - times[0])
    if span <= 0.0:
        return DecayRateEstimate(
            status=EstimateStatus.INSUFFICIENT,
            method="nonlinear_exponential_least_squares",
            estimated_rate_per_s=None,
            estimated_amplitude=None,
            fit_converged=False,
            convergence_message="observation times have no span",
            input_sample_count=len(times),
            observation_rmse=None,
            observation_rmse_unit=observations.value_unit,
            log_rmse=None,
            local_rate_noise_sensitivity_per_s=None,
            jacobian_condition_number=None,
        )
    positive = np.flatnonzero(values > 3.0 * observations.noise_standard_deviation)
    if len(positive) >= 2:
        estimate = -float(np.polyfit(times[positive], np.log(values[positive]), 1)[0])
        rate_initial = max(1e-10, estimate)
        amplitude_initial = max(1e-12, float(values[positive[0]]))
    else:
        rate_initial = 1.0 / span
        amplitude_initial = max(1e-12, float(np.max(np.abs(values))))

    def model(parameters: NDArray[np.float64]) -> NDArray[np.float64]:
        amplitude, rate = float(parameters[0]), float(parameters[1])
        return amplitude * np.exp(-rate * times)

    def jacobian(parameters: NDArray[np.float64]) -> NDArray[np.float64]:
        amplitude, rate = float(parameters[0]), float(parameters[1])
        decay = np.exp(-rate * times)
        return np.column_stack((decay, -times * amplitude * decay))

    try:
        fit = least_squares(
            lambda parameters: model(parameters) - values,
            np.asarray((amplitude_initial, rate_initial)),
            jac=jacobian,
            bounds=((1e-15, 1e-15), (np.inf, np.inf)),
            x_scale="jac",
            max_nfev=10000,
        )
    except (ValueError, FloatingPointError, RuntimeError) as exc:
        return DecayRateEstimate(
            status=EstimateStatus.FIT_FAILED,
            method="nonlinear_exponential_least_squares",
            estimated_rate_per_s=None,
            estimated_amplitude=None,
            fit_converged=False,
            convergence_message=f"{type(exc).__name__}: {exc}",
            input_sample_count=len(times),
            observation_rmse=None,
            observation_rmse_unit=observations.value_unit,
            log_rmse=None,
            local_rate_noise_sensitivity_per_s=None,
            jacobian_condition_number=None,
        )
    amplitude, rate = (float(value) for value in fit.x)
    fitted = model(fit.x)
    residual = values - fitted
    rmse = float(np.sqrt(np.mean(np.square(residual))))
    condition = _scaled_condition_number(fit.jac)
    sensitivity: float | None = None
    if observations.noise_standard_deviation > 0.0:
        covariance = np.linalg.pinv(fit.jac.T @ fit.jac)
        sensitivity = float(
            observations.noise_standard_deviation * sqrt(max(0.0, covariance[1, 1]))
        )
    status, diagnostics = _status_for_fit(
        rate=rate,
        noise_sensitivity=sensitivity,
        condition_number=condition,
        log_rmse=None,
        excluded_fraction=0.0,
    )
    return DecayRateEstimate(
        status=status if fit.success else EstimateStatus.FIT_FAILED,
        method="nonlinear_exponential_least_squares",
        estimated_rate_per_s=rate,
        estimated_amplitude=amplitude,
        fit_converged=bool(fit.success),
        convergence_message=str(fit.message),
        input_sample_count=len(times),
        used_sample_indices=tuple(range(len(times))),
        fitted_values=tuple(float(value) for value in fitted),
        observation_rmse=rmse,
        observation_rmse_unit=observations.value_unit,
        log_rmse=None,
        local_rate_noise_sensitivity_per_s=sensitivity,
        jacobian_condition_number=condition,
        diagnostics=diagnostics,
    )


def _scaled_condition_number(jacobian: np.ndarray) -> float | None:
    norms = np.linalg.norm(jacobian, axis=0)
    if np.any(norms == 0.0) or not np.isfinite(norms).all():
        return float("inf")
    scaled = jacobian / norms
    value = float(np.linalg.cond(scaled))
    return value if isfinite(value) else float("inf")


def estimate_oscillator_envelope(
    observations: ObservationSeries,
) -> DecayRateEstimate:
    """Fit log absolute local displacement peaks; never log signed samples."""
    values = np.asarray(observations.values, dtype=np.float64)
    times = np.asarray(observations.times, dtype=np.float64)
    magnitudes = np.abs(values)
    peak_indices = (
        np.flatnonzero(
            (magnitudes[1:-1] >= magnitudes[:-2]) & (magnitudes[1:-1] >= magnitudes[2:])
        )
        + 1
    )
    if len(magnitudes) > 1 and magnitudes[0] >= magnitudes[1]:
        peak_indices = np.insert(peak_indices, 0, 0)
    if len(magnitudes) > 1 and magnitudes[-1] >= magnitudes[-2]:
        peak_indices = np.append(peak_indices, len(magnitudes) - 1)
    peak_indices = np.unique(peak_indices)
    exclusions = [
        ExcludedObservation(
            sample_index=index,
            reason="not a local maximum of absolute sampled displacement",
        )
        for index in range(len(values))
        if index not in set(int(value) for value in peak_indices)
    ]
    accepted_peaks: list[int] = []
    for index in peak_indices:
        if magnitudes[index] <= 0.0:
            exclusions.append(
                ExcludedObservation(
                    sample_index=int(index),
                    reason="zero peak cannot be log-transformed",
                )
            )
        elif (
            observations.noise_standard_deviation > 0.0
            and magnitudes[index] <= 3.0 * observations.noise_standard_deviation
        ):
            exclusions.append(
                ExcludedObservation(
                    sample_index=int(index),
                    reason="peak is at or below three-sigma floor",
                )
            )
        else:
            accepted_peaks.append(int(index))
    if len(accepted_peaks) < 4:
        return DecayRateEstimate(
            status=EstimateStatus.INSUFFICIENT,
            method="absolute_sample_peak_log_linear",
            estimated_rate_per_s=None,
            estimated_amplitude=None,
            fit_converged=False,
            convergence_message="fewer than four usable absolute local peaks",
            input_sample_count=len(values),
            used_sample_indices=tuple(accepted_peaks),
            excluded_observations=tuple(
                sorted(exclusions, key=lambda item: item.sample_index)
            ),
            observation_rmse=None,
            observation_rmse_unit=observations.value_unit,
            log_rmse=None,
            local_rate_noise_sensitivity_per_s=None,
            jacobian_condition_number=None,
            diagnostics=(
                "absolute-peak envelope fit requires at least four above-floor peaks",
            ),
        )
    peak_series = ObservationSeries(
        scenario_id=observations.scenario_id,
        observable="absolute_local_displacement_peak",
        times=tuple(float(times[index]) for index in accepted_peaks),
        values=tuple(float(magnitudes[index]) for index in accepted_peaks),
        time_unit=observations.time_unit,
        value_unit=observations.value_unit,
        noise_model="none",
        noise_standard_deviation=0.0,
    )
    fit = estimate_exponential_rate(peak_series)
    intervals = np.diff(times)
    median_interval = float(np.median(intervals))
    estimated_period = 2.0 * float(np.median(np.diff(times[accepted_peaks])))
    samples_per_period = estimated_period / median_interval
    sampling_diagnostics: tuple[str, ...] = ()
    if samples_per_period < 16.0:
        sampling_diagnostics = (
            f"only {samples_per_period:.3g} samples per estimated oscillation period; "
            "sampled-peak rate is cadence-sensitive",
        )
    mapped_exclusions = tuple(
        sorted(
            exclusions
            + [
                ExcludedObservation(
                    sample_index=accepted_peaks[item.sample_index], reason=item.reason
                )
                for item in fit.excluded_observations
            ],
            key=lambda item: item.sample_index,
        )
    )
    return fit.model_copy(
        update={
            "method": "absolute_sample_peak_log_linear",
            "input_sample_count": len(values),
            "used_sample_indices": tuple(
                accepted_peaks[index] for index in fit.used_sample_indices
            ),
            "excluded_observations": mapped_exclusions,
            "status": (
                EstimateStatus.UNRELIABLE
                if sampling_diagnostics and fit.estimated_rate_per_s is not None
                else fit.status
            ),
            "diagnostics": fit.diagnostics
            + (
                "local peak selection assumes sampled extrema approximate "
                "envelope peaks; "
                "sparse cadence can bias peak timing and height",
            )
            + sampling_diagnostics,
        }
    )


def evaluate_estimate(
    estimate: DecayRateEstimate, truth: SyntheticGroundTruth
) -> RateEvaluation:
    """Compare an estimate with separately held truth after fitting is complete."""
    value = estimate.estimated_rate_per_s
    target = truth.target_rate_per_s
    if value is None:
        return RateEvaluation(
            method=estimate.method,
            target_rate_per_s=target,
            absolute_error_per_s=None,
            relative_error=None,
            nearest_component_rate_per_s=None,
            nearest_component_gap_per_s=None,
            within_declared_10_percent_tolerance=None,
            interpretation="no estimate available for post-fit comparison",
        )
    if target is not None:
        absolute = abs(value - target)
        relative = absolute / abs(target)
        return RateEvaluation(
            method=estimate.method,
            target_rate_per_s=target,
            absolute_error_per_s=absolute,
            relative_error=relative,
            nearest_component_rate_per_s=target,
            nearest_component_gap_per_s=absolute,
            within_declared_10_percent_tolerance=relative <= 0.1,
            interpretation=(
                "within the study's 10% descriptive error tolerance"
                if relative <= 0.1
                else "outside the study's 10% descriptive error tolerance"
            ),
        )
    nearest = min(truth.component_rates_per_s, key=lambda rate: abs(value - rate))
    gap = abs(value - nearest)
    relative_gap = gap / abs(nearest)
    return RateEvaluation(
        method=estimate.method,
        target_rate_per_s=None,
        absolute_error_per_s=None,
        relative_error=None,
        nearest_component_rate_per_s=nearest,
        nearest_component_gap_per_s=gap,
        within_declared_10_percent_tolerance=relative_gap <= 0.1,
        interpretation=(
            "single-exponential estimate is an effective mixture rate; "
            "nearest-component gap is descriptive, not component identification"
        ),
    )


def _outcome(
    case: SyntheticObservationCase,
    condition_summary: str,
    methods: tuple[DecayRateEstimate, ...],
) -> EstimationExperiment:
    return EstimationExperiment(
        observations=case.observations,
        ground_truth=case.ground_truth,
        condition_summary=condition_summary,
        outcomes=tuple(
            EstimateOutcome(
                estimate=estimate,
                evaluation=evaluate_estimate(estimate, case.ground_truth),
            )
            for estimate in methods
        ),
    )


def _uniform_times(duration_s: float, interval_s: float) -> tuple[float, ...]:
    count = int(round(duration_s / interval_s)) + 1
    return tuple(float(value) for value in np.linspace(0.0, duration_s, count))


def _window_observations(
    observations: ObservationSeries, start_s: float, end_s: float, scenario_id: str
) -> ObservationSeries:
    selected = tuple(
        index
        for index, time in enumerate(observations.times)
        if start_s <= time <= end_s
    )
    return ObservationSeries(
        scenario_id=scenario_id,
        observable=observations.observable,
        times=tuple(observations.times[index] for index in selected),
        values=tuple(observations.values[index] for index in selected),
        time_unit=observations.time_unit,
        value_unit=observations.value_unit,
        noise_model=observations.noise_model,
        noise_standard_deviation=observations.noise_standard_deviation,
        noise_seed=observations.noise_seed,
    )


def run_phase23_research() -> Phase23ResearchResults:
    """Run fixed-seed observation experiments for exponential, oscillator, heat."""
    experiments: list[EstimationExperiment] = []

    for index, interval in enumerate((0.05, 0.25, 0.8), start=1):
        case = generate_exponential_observations(
            scenario_id=f"scalar_cadence_{index}",
            amplitude=1.0,
            rate_per_s=0.3,
            times_s=_uniform_times(4.0, interval),
            value_unit="arbitrary signal units",
            noise_standard_deviation=0.02,
            noise_seed=2310,
        )
        fits = (
            estimate_exponential_rate(case.observations),
            estimate_exponential_rate(case.observations, method="nonlinear"),
        )
        experiments.append(
            _outcome(case, f"cadence dt={interval:g} s; duration=4 s", fits)
        )

    for index, duration in enumerate((1.0, 4.0, 10.0), start=1):
        case = generate_exponential_observations(
            scenario_id=f"scalar_duration_{index}",
            amplitude=1.0,
            rate_per_s=0.3,
            times_s=_uniform_times(duration, 0.25),
            value_unit="arbitrary signal units",
            noise_standard_deviation=0.02,
            noise_seed=2320,
        )
        fits = (
            estimate_exponential_rate(case.observations),
            estimate_exponential_rate(case.observations, method="nonlinear"),
        )
        experiments.append(
            _outcome(case, f"duration={duration:g} s; cadence dt=0.25 s", fits)
        )

    for index, noise in enumerate((0.001, 0.01, 0.05), start=1):
        case = generate_exponential_observations(
            scenario_id=f"scalar_noise_{index}",
            amplitude=1.0,
            rate_per_s=0.3,
            times_s=_uniform_times(12.0, 0.25),
            value_unit="arbitrary signal units",
            noise_standard_deviation=noise,
            noise_seed=2330,
        )
        fits = (
            estimate_exponential_rate(case.observations),
            estimate_exponential_rate(case.observations, method="nonlinear"),
        )
        experiments.append(
            _outcome(case, f"additive Gaussian sigma={noise:g}; seed=2330", fits)
        )

    period = 2.0 * pi / sqrt(4.0 - 0.05**2)
    for index, interval in enumerate((0.05, 0.4, 10.0 / 12.0, 10.0 / 6.0), start=1):
        points = int(round(10.0 / interval)) + 1
        actual_interval = 10.0 / (points - 1)
        case = generate_oscillator_observations(
            scenario_id=f"oscillator_cadence_{index}",
            mass_kg=1.0,
            damping_coefficient_kg_per_s=0.1,
            stiffness_n_per_m=4.0,
            initial_displacement_m=0.1,
            initial_velocity_m_per_s=0.0,
            times_s=_uniform_times(10.0, actual_interval),
        )
        estimate = estimate_oscillator_envelope(case.observations)
        experiments.append(
            _outcome(
                case,
                f"oscillator dt={actual_interval:.6g} s; duration=10 s; "
                f"period={period:.6g} s; samples/period={period / actual_interval:.4g}",
                (estimate,),
            )
        )

    for index, duration in enumerate((1.5, 3.0, 10.0), start=1):
        case = generate_oscillator_observations(
            scenario_id=f"oscillator_duration_{index}",
            mass_kg=1.0,
            damping_coefficient_kg_per_s=0.1,
            stiffness_n_per_m=4.0,
            initial_displacement_m=0.1,
            initial_velocity_m_per_s=0.0,
            times_s=_uniform_times(duration, 0.05),
        )
        estimate = estimate_oscillator_envelope(case.observations)
        experiments.append(
            _outcome(
                case,
                f"oscillator duration={duration:g} s; cadence dt=0.05 s; "
                f"period={period:.6g} s",
                (estimate,),
            )
        )

    alpha = 1e-4
    length = 1.0
    for mode, amplitude, duration, interval in (
        (1, 1.0, 1200.0, 60.0),
        (4, 0.5, 300.0, 15.0),
    ):
        case = generate_heat_mode_observations(
            scenario_id=f"heat_single_mode_{mode}",
            diffusivity_m2_per_s=alpha,
            mode_number=mode,
            length_m=length,
            initial_amplitude_k=amplitude,
            times_s=_uniform_times(duration, interval),
        )
        experiments.append(
            _outcome(
                case,
                f"individual heat mode n={mode}; L=1 m; alpha=1e-4 m^2/s; "
                f"dt={interval:g} s; duration={duration:g} s",
                (estimate_exponential_rate(case.observations),),
            )
        )

    mixture_specs = (
        ("mix_separated_equal", (1, 4), (1.0, 0.5), (0.0, 1200.0, 20.0)),
        ("mix_weak_fast", (1, 4), (1.0, 0.01), (0.0, 1200.0, 20.0)),
        ("mix_strong_fast", (1, 4), (0.01, 1.0), (0.0, 1200.0, 20.0)),
        ("mix_close_rates", (10, 11), (1.0, 1.0), (0.0, 80.0, 2.0)),
    )
    for scenario_id, mode_numbers, amplitudes, window in mixture_specs:
        start, end, interval = window
        case = generate_mixed_heat_observations(
            scenario_id=scenario_id,
            mode_numbers=mode_numbers,
            effective_amplitudes_k=amplitudes,
            diffusivity_m2_per_s=alpha,
            length_m=length,
            times_s=_uniform_times(end, interval),
        )
        full = estimate_exponential_rate(case.observations)
        early_obs = _window_observations(
            case.observations,
            start,
            min(start + end * 0.1, end),
            f"{scenario_id}_early",
        )
        late_obs = _window_observations(
            case.observations, end * 0.5, end, f"{scenario_id}_late"
        )
        estimates = (
            full,
            estimate_exponential_rate(early_obs),
            estimate_exponential_rate(late_obs),
        )
        experiments.append(
            _outcome(
                case,
                f"heat mode mixture n={mode_numbers}, effective amplitudes "
                f"K={amplitudes}; "
                f"early/full/late windows; dt={interval:g} s",
                estimates,
            )
        )

    noisy_mixture = generate_mixed_heat_observations(
        scenario_id="mix_fast_component_below_noise_floor",
        mode_numbers=(1, 4),
        effective_amplitudes_k=(0.01, 1.0),
        diffusivity_m2_per_s=alpha,
        length_m=length,
        times_s=_uniform_times(1200.0, 20.0),
        noise_standard_deviation_k=0.005,
        noise_seed=2340,
    )
    noisy_early = _window_observations(
        noisy_mixture.observations,
        0.0,
        120.0,
        "mix_fast_component_below_noise_floor_early",
    )
    noisy_late = _window_observations(
        noisy_mixture.observations,
        600.0,
        1200.0,
        "mix_fast_component_below_noise_floor_late",
    )
    experiments.append(
        _outcome(
            noisy_mixture,
            "strong-fast n=(1,4) mixture; sigma=0.005 K; seed=2340; "
            "fast-mode tail crosses the three-sigma floor; dt=20 s",
            (
                estimate_exponential_rate(noisy_mixture.observations),
                estimate_exponential_rate(noisy_early),
                estimate_exponential_rate(noisy_late),
            ),
        )
    )

    return Phase23ResearchResults(
        experiments=tuple(experiments),
        seed_policy=(
            "Synthetic additive Gaussian noise uses NumPy default_rng with fixed "
            "scenario seeds recorded in each observation series. Noise-free cases "
            "are analytical samples."
        ),
        interpretation_policy=(
            "Estimators receive ObservationSeries only; ground truth is stored "
            "separately.",
            "INTERPRETABLE, FIT_FAILED, INSUFFICIENT, NOT_IDENTIFIABLE, and "
            "UNRELIABLE are distinct statuses.",
            "A local noise-sensitivity threshold of three times the linearized "
            "rate scale is a study screen, not a confidence interval or universal "
            "identifiability threshold.",
            "A 10% rate error is used only as a descriptive post-fit score in "
            "this tested design.",
            "Single-exponential estimates of mixtures are effective "
            "window-dependent summaries, not identified physical mode rates.",
        ),
    )


def render_phase23_results(results: Phase23ResearchResults) -> str:
    """Render compact deterministic tables for the Phase 23 research report."""
    lines = [
        "# Phase 23 computed modal-estimation results",
        "",
        "| Scenario | Condition | Method | Status | Estimate (s^-1) | "
        "Truth / nearest mode (s^-1) | Abs. error/gap (s^-1) | Rel. error | "
        "n used / excluded | RMSE (observation units) |",
        "|---|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for experiment in results.experiments:
        for outcome in experiment.outcomes:
            estimate = outcome.estimate
            evaluation = outcome.evaluation
            rate = (
                "unavailable"
                if estimate.estimated_rate_per_s is None
                else f"{estimate.estimated_rate_per_s:.7g}"
            )
            reference = (
                evaluation.target_rate_per_s
                if evaluation.target_rate_per_s is not None
                else evaluation.nearest_component_rate_per_s
            )
            reference_text = "n/a" if reference is None else f"{reference:.7g}"
            gap = (
                evaluation.absolute_error_per_s
                if evaluation.absolute_error_per_s is not None
                else evaluation.nearest_component_gap_per_s
            )
            gap_text = "n/a" if gap is None else f"{gap:.4g}"
            relative = (
                "n/a"
                if evaluation.relative_error is None
                else f"{evaluation.relative_error:.4g}"
            )
            rmse = (
                "n/a"
                if estimate.observation_rmse is None
                else f"{estimate.observation_rmse:.4g}"
            )
            n_used = len(estimate.used_sample_indices)
            n_excluded = len(estimate.excluded_observations)
            lines.append(
                f"| {experiment.observations.scenario_id} | "
                f"{experiment.condition_summary} | "
                f"{estimate.method} | {estimate.status.value} | {rate} | "
                f"{reference_text} | {gap_text} | {relative} | "
                f"{n_used}/{n_excluded} | {rmse} {estimate.observation_rmse_unit} |"
            )
    lines.extend(("", "## Estimation diagnostics", ""))
    for experiment in results.experiments:
        for outcome in experiment.outcomes:
            estimate = outcome.estimate
            reason_indices: dict[str, list[int]] = {}
            for item in estimate.excluded_observations:
                reason_indices.setdefault(item.reason, []).append(item.sample_index)
            exclusions = (
                "; ".join(
                    f"{len(indices)} x {reason} (indices {indices[:8]}"
                    f"{' ...' if len(indices) > 8 else ''})"
                    for reason, indices in reason_indices.items()
                )
                or "none"
            )
            lines.append(
                f"- `{experiment.observations.scenario_id}` / `{estimate.method}`: "
                f"converged={estimate.fit_converged}; "
                f"message={estimate.convergence_message}; "
                f"rate-noise-sensitivity="
                f"{estimate.local_rate_noise_sensitivity_per_s}; "
                f"Jacobian condition={estimate.jacobian_condition_number}; "
                f"log RMSE={estimate.log_rmse}; exclusions={exclusions}. "
                f"{'; '.join(estimate.diagnostics)}"
            )
    lines.extend(("", results.seed_policy, ""))
    return "\n".join(lines)
