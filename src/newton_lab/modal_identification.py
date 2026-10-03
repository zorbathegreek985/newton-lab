"""Phase 24 controlled robustness and two-mode identification experiments.

Estimators receive observation records and explicit initial guesses only.
Generator truth is retained in a separate case record and is used only after
model fitting for post-fit evaluation.
"""

from __future__ import annotations

from enum import StrEnum
from math import isfinite
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray
from pydantic import Field, model_validator
from scipy.optimize import least_squares  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.modal_estimation import (
    DecayRateEstimate,
    ObservationSeries,
    SyntheticGroundTruth,
    estimate_exponential_rate,
    estimate_oscillator_envelope,
    generate_oscillator_observations,
)
from newton_lab.simulation import SimulationModel

BIC_MODEL_THRESHOLD = 6.0
"""Predeclared BIC difference screen; a study rule, not a universal cutoff."""

RATE_RECOVERY_RELATIVE_TOLERANCE = 0.10
"""Descriptive post-fit per-mode tolerance, not a universal standard."""


class IdentificationStatus(StrEnum):
    """Post-fit outcomes after estimates are evaluated against hidden truth."""

    SUPPORTED = "supported_under_tested_conditions"
    INCONCLUSIVE = "mode_identification_inconclusive"
    INSUFFICIENT = "insufficient_observations"
    FIT_FAILED = "numerical_fitting_failure"
    INCONSISTENT = "fitted_model_inconsistent_with_synthetic_truth"
    UNRELIABLE = "unreliable_fit_or_diagnostics"


class ModelSelection(StrEnum):
    """Selection based on the known-noise Gaussian BIC difference."""

    SINGLE = "single_mode"
    TWO = "two_modes"
    INCONCLUSIVE = "inconclusive"
    UNAVAILABLE = "unavailable"


class ModeObservationCase(SimulationModel):
    """Observed data, exact noise-free signal, and separately held truth."""

    observations: ObservationSeries
    noise_free_signal_values: tuple[float, ...]
    ground_truth: SyntheticGroundTruth

    @model_validator(mode="after")
    def validate_case_alignment(self) -> ModeObservationCase:
        if len(self.observations.values) != len(self.noise_free_signal_values):
            raise ValueError("noise-free and observed samples must be aligned")
        if not np.isfinite(self.noise_free_signal_values).all():
            raise ValueError("noise-free signal values must be finite")
        if self.observations.scenario_id != self.ground_truth.scenario_id:
            raise ValueError("observation and truth scenario IDs must match")
        return self


class StartFitRecord(SimulationModel):
    """Outcome from one explicit two-mode nonlinear initialization."""

    initial_amplitudes: tuple[float, float]
    initial_rates_per_s: tuple[float, float]
    converged: bool
    convergence_message: str
    residual_rmse: float | None
    fitted_amplitudes: tuple[float, float] | None
    fitted_rates_per_s: tuple[float, float] | None


class TwoModeFit(SimulationModel):
    """Best converged positive two-exponential fit and diagnostics."""

    converged: bool
    convergence_message: str
    fitted_amplitudes: tuple[float, float] | None
    fitted_rates_per_s: tuple[float, float] | None
    fitted_values: tuple[float, ...] = ()
    residual_rmse: float | None
    residual_rmse_unit: str
    jacobian_condition_number: float | None
    start_records: tuple[StartFitRecord, ...] = ()
    free_parameter_count: int = 4


class ModelComparison(SimulationModel):
    """Single/two-mode fit comparison using known-variance Gaussian BIC."""

    single_mode_fit: DecayRateEstimate
    two_mode_fit: TwoModeFit
    single_free_parameter_count: int = 2
    two_free_parameter_count: int = 4
    single_bic: float | None
    two_bic: float | None
    bic_difference_single_minus_two: float | None
    selected_model: ModelSelection
    comparison_assumption: str


class MatchedModeError(SimulationModel):
    """Post-fit rate errors after ascending-rate mode matching."""

    true_rate_per_s: float
    estimated_rate_per_s: float
    absolute_error_per_s: float
    relative_error: float


class IdentificationEvaluation(SimulationModel):
    """Truth-aware evaluation, kept separate from the fitting result."""

    status: IdentificationStatus
    true_mode_count: int
    selected_model: ModelSelection
    matched_rate_errors: tuple[MatchedModeError, ...] = ()
    matched_amplitude_errors: tuple[float, ...] = ()
    interpretation: str


class IdentificationRun(SimulationModel):
    """One generated case, comparison, and separate truth evaluation."""

    case: ModeObservationCase
    design_name: str
    initial_guess_design: str
    comparison: ModelComparison
    evaluation: IdentificationEvaluation


class SamplingCoverage(SimulationModel):
    """Observed timing coverage for one sampling design."""

    design_name: str
    sample_count: int
    duration_s: float
    minimum_interval_s: float
    median_interval_s: float
    maximum_interval_s: float


class RepeatedNoiseSummary(SimulationModel):
    """Descriptive across-seed rate and identification summaries."""

    condition_name: str
    realization_count: int
    successful_fit_count: int
    failed_fit_count: int
    unreliable_count: int
    inconclusive_count: int
    supported_count: int
    inconsistent_count: int
    per_mode_summaries: tuple[ModeRateDistribution, ...]
    seeds: tuple[int, ...]
    interpretation: str


class ModeRateDistribution(SimulationModel):
    """Descriptive estimate and absolute/relative error summaries by mode."""

    mode_name: str
    estimate_median_per_s: float
    estimate_minimum_per_s: float
    estimate_maximum_per_s: float
    estimate_population_sd_per_s: float
    signed_error_mean_per_s: float
    absolute_error_median_per_s: float
    absolute_error_minimum_per_s: float
    absolute_error_maximum_per_s: float
    absolute_error_population_sd_per_s: float
    relative_error_median: float
    relative_error_minimum: float
    relative_error_maximum: float
    relative_error_population_sd: float


class Phase24ResearchResults(SimulationModel):
    """Repeated-noise, sampling-design, and mode-matrix results."""

    repeated_noise_runs: tuple[IdentificationRun, ...]
    repeated_noise_summaries: tuple[RepeatedNoiseSummary, ...]
    mode_matrix_runs: tuple[IdentificationRun, ...]
    sampling_design_runs: tuple[IdentificationRun, ...]
    sampling_coverages: tuple[SamplingCoverage, ...]
    oscillator_sampling_summary: tuple[tuple[str, str, float | None, int, float], ...]
    initial_guess_comparison: tuple[
        tuple[
            tuple[float, float],
            tuple[float, float],
            bool,
            str,
            tuple[float | None, float | None],
        ],
        ...,
    ]
    selection_assumptions: tuple[str, ...] = Field(min_length=1)


def _validated_times(times_s: tuple[float, ...]) -> tuple[float, ...]:
    times = tuple(float(value) for value in times_s)
    if len(times) < 3 or not np.isfinite(times).all():
        raise ScientificValidationError(
            "at least three finite observation times are required"
        )
    if times[0] < 0.0 or not np.all(np.diff(times) > 0.0):
        raise ScientificValidationError(
            "observation times must be non-negative and strictly increasing"
        )
    return times


def generate_mode_observations(
    *,
    scenario_id: str,
    amplitudes: tuple[float, ...],
    rates_per_s: tuple[float, ...],
    times_s: tuple[float, ...],
    noise_standard_deviation: float,
    noise_seed: int,
) -> ModeObservationCase:
    """Generate a single- or positive two-mode exponential observation case.

    The estimator receives only ``case.observations``. Noise-free values and
    generator parameters stay in separate fields. Gaussian noise is iid,
    additive, and seeded independently for each realization.
    """
    if len(amplitudes) not in (1, 2) or len(amplitudes) != len(rates_per_s):
        raise ScientificValidationError("one or two aligned modal terms are required")
    if any(not isfinite(value) or value <= 0.0 for value in amplitudes):
        raise ScientificValidationError("amplitudes must be finite and positive")
    if any(not isfinite(value) or value <= 0.0 for value in rates_per_s):
        raise ScientificValidationError("rates must be finite and positive")
    if len(rates_per_s) == 2 and rates_per_s[0] == rates_per_s[1]:
        raise ScientificValidationError("two-mode rates must be distinct")
    if not isfinite(noise_standard_deviation) or noise_standard_deviation < 0.0:
        raise ScientificValidationError(
            "noise standard deviation must be finite and non-negative"
        )
    if (
        isinstance(noise_seed, bool)
        or not isinstance(noise_seed, int)
        or noise_seed < 0
    ):
        raise ScientificValidationError("noise_seed must be a non-negative integer")
    times = _validated_times(times_s)
    time_array = np.asarray(times, dtype=np.float64)
    noise_free = np.sum(
        np.vstack(
            [
                amplitude * np.exp(-rate * time_array)
                for amplitude, rate in zip(amplitudes, rates_per_s, strict=True)
            ]
        ),
        axis=0,
    )
    observed = noise_free.copy()
    if noise_standard_deviation > 0.0:
        observed += np.random.default_rng(noise_seed).normal(
            0.0, noise_standard_deviation, size=len(times)
        )
    observation = ObservationSeries(
        scenario_id=scenario_id,
        observable="positive_modal_signal",
        times=times,
        values=tuple(float(value) for value in observed),
        time_unit="s",
        value_unit="signal units",
        noise_model=(
            "iid_gaussian_additive" if noise_standard_deviation > 0.0 else "none"
        ),
        noise_standard_deviation=noise_standard_deviation,
        noise_seed=noise_seed if noise_standard_deviation > 0.0 else None,
    )
    truth = SyntheticGroundTruth(
        scenario_id=scenario_id,
        model_description=f"sum of {len(amplitudes)} positive exponential mode(s)",
        component_rates_per_s=rates_per_s,
        component_amplitudes=amplitudes,
        amplitude_unit="signal units",
        target_rate_per_s=rates_per_s[0] if len(rates_per_s) == 1 else None,
    )
    return ModeObservationCase(
        observations=observation,
        noise_free_signal_values=tuple(float(value) for value in noise_free),
        ground_truth=truth,
    )


def _two_mode_values(
    parameters: NDArray[np.float64], times: NDArray[np.float64]
) -> NDArray[np.float64]:
    amplitude_1, rate_1, amplitude_2, rate_2 = (
        float(parameters[index]) for index in range(4)
    )
    return amplitude_1 * np.exp(-rate_1 * times) + amplitude_2 * np.exp(-rate_2 * times)


def _two_mode_jacobian(
    parameters: NDArray[np.float64], times: NDArray[np.float64]
) -> NDArray[np.float64]:
    amplitude_1, rate_1, amplitude_2, rate_2 = (
        float(parameters[index]) for index in range(4)
    )
    mode_1 = np.exp(-rate_1 * times)
    mode_2 = np.exp(-rate_2 * times)
    return np.column_stack(
        (mode_1, -amplitude_1 * times * mode_1, mode_2, -amplitude_2 * times * mode_2)
    )


def default_two_mode_initial_guesses(
    observations: ObservationSeries,
) -> tuple[tuple[float, float, float, float], ...]:
    """Build deterministic starts from observations and time span only."""
    times = np.asarray(observations.times, dtype=np.float64)
    span = float(times[-1] - times[0])
    if span <= 0.0:
        return ()
    scale = max(float(np.max(np.abs(observations.values))), 1e-12)
    base = 1.0 / span
    return (
        (0.5 * scale, 0.5 * base, 0.5 * scale, 4.0 * base),
        (0.8 * scale, 0.2 * base, 0.2 * scale, 8.0 * base),
        (0.2 * scale, 3.0 * base, 0.8 * scale, 0.5 * base),
        (0.5 * scale, 1.0 * base, 0.5 * scale, 1.1 * base),
        (0.9 * scale, 0.1 * base, 0.1 * scale, 12.0 * base),
    )


def fit_two_mode_exponential(
    observations: ObservationSeries,
    *,
    initial_guesses: tuple[tuple[float, float, float, float], ...] | None = None,
) -> TwoModeFit:
    """Fit ``A1 exp(-r1*t) + A2 exp(-r2*t)`` from observations alone.

    Amplitudes and rates are constrained positive. Rates are sorted only after
    fitting, so the symmetry under swapping modes cannot change the fitted
    curve. Multiple deterministic initializations address, but do not remove,
    local-minimum and weak-identification risks.
    """
    times = np.asarray(observations.times, dtype=np.float64)
    values = np.asarray(observations.values, dtype=np.float64)
    span = float(times[-1] - times[0])
    if len(times) < 8 or span <= 0.0:
        return TwoModeFit(
            converged=False,
            convergence_message=(
                "at least eight samples and positive time span required"
            ),
            fitted_amplitudes=None,
            fitted_rates_per_s=None,
            residual_rmse=None,
            residual_rmse_unit=observations.value_unit,
            jacobian_condition_number=None,
        )
    starts = (
        default_two_mode_initial_guesses(observations)
        if initial_guesses is None
        else initial_guesses
    )
    if not starts:
        return TwoModeFit(
            converged=False,
            convergence_message="no initial guesses supplied",
            fitted_amplitudes=None,
            fitted_rates_per_s=None,
            residual_rmse=None,
            residual_rmse_unit=observations.value_unit,
            jacobian_condition_number=None,
        )
    for start in starts:
        if len(start) != 4 or any(
            not isfinite(value) or value <= 0.0 for value in start
        ):
            raise ScientificValidationError(
                "each initial guess must contain positive finite amplitudes and rates"
            )

    amplitude_upper_scale = max(float(np.max(np.abs(values))), 1e-12)
    rate_lower = max(1e-12, 1e-8 / span)
    rate_upper = max(50.0 / span, rate_lower * 100.0)
    records: list[StartFitRecord] = []
    successful: list[tuple[float, Any]] = []
    for start in starts:
        start_array = np.asarray(start, dtype=np.float64)
        initial = np.asarray(
            (
                min(max(start_array[0], 1e-12), amplitude_upper_scale * 1e6),
                min(max(start_array[1], rate_lower * 1.01), rate_upper * 0.99),
                min(max(start_array[2], 1e-12), amplitude_upper_scale * 1e6),
                min(max(start_array[3], rate_lower * 1.01), rate_upper * 0.99),
            ),
            dtype=np.float64,
        )
        try:
            fit = least_squares(
                lambda parameters: _two_mode_values(parameters, times) - values,
                initial,
                jac=lambda parameters: _two_mode_jacobian(parameters, times),
                bounds=(
                    (1e-12, rate_lower, 1e-12, rate_lower),
                    (
                        amplitude_upper_scale * 1e6,
                        rate_upper,
                        amplitude_upper_scale * 1e6,
                        rate_upper,
                    ),
                ),
                x_scale="jac",
                max_nfev=5000,
                ftol=1e-11,
                xtol=1e-11,
                gtol=1e-11,
            )
            parameters = np.asarray(fit.x, dtype=np.float64)
            residual = values - _two_mode_values(parameters, times)
            rmse = float(np.sqrt(np.mean(np.square(residual))))
            converged = bool(fit.success and np.isfinite(parameters).all())
            amplitudes = (float(parameters[0]), float(parameters[2]))
            rates = (float(parameters[1]), float(parameters[3]))
            if converged:
                successful.append((float(np.sum(np.square(residual))), fit))
            records.append(
                StartFitRecord(
                    initial_amplitudes=(float(start[0]), float(start[2])),
                    initial_rates_per_s=(float(start[1]), float(start[3])),
                    converged=converged,
                    convergence_message=str(fit.message),
                    residual_rmse=rmse,
                    fitted_amplitudes=amplitudes,
                    fitted_rates_per_s=rates,
                )
            )
        except (ValueError, FloatingPointError, RuntimeError) as exc:
            records.append(
                StartFitRecord(
                    initial_amplitudes=(float(start[0]), float(start[2])),
                    initial_rates_per_s=(float(start[1]), float(start[3])),
                    converged=False,
                    convergence_message=f"{type(exc).__name__}: {exc}",
                    residual_rmse=None,
                    fitted_amplitudes=None,
                    fitted_rates_per_s=None,
                )
            )
    if not successful:
        return TwoModeFit(
            converged=False,
            convergence_message="all two-mode initializations failed",
            fitted_amplitudes=None,
            fitted_rates_per_s=None,
            residual_rmse=None,
            residual_rmse_unit=observations.value_unit,
            jacobian_condition_number=None,
            start_records=tuple(records),
        )
    _, best = min(successful, key=lambda item: item[0])
    parameters = np.asarray(best.x, dtype=np.float64)
    pairs = sorted(
        (
            (float(parameters[1]), float(parameters[0])),
            (float(parameters[3]), float(parameters[2])),
        )
    )
    rates = (pairs[0][0], pairs[1][0])
    amplitudes = (pairs[0][1], pairs[1][1])
    fitted = _two_mode_values(parameters, times)
    rmse = float(np.sqrt(np.mean(np.square(values - fitted))))
    jacobian = np.asarray(best.jac, dtype=np.float64)
    norms = np.linalg.norm(jacobian, axis=0)
    condition = (
        float(np.linalg.cond(jacobian / norms))
        if np.all(norms > 0.0) and np.isfinite(norms).all()
        else float("inf")
    )
    return TwoModeFit(
        converged=True,
        convergence_message=str(best.message),
        fitted_amplitudes=amplitudes,
        fitted_rates_per_s=rates,
        fitted_values=tuple(float(value) for value in fitted),
        residual_rmse=rmse,
        residual_rmse_unit=observations.value_unit,
        jacobian_condition_number=condition,
        start_records=tuple(records),
    )


def compare_single_and_two_mode_models(
    observations: ObservationSeries,
    *,
    initial_guesses: tuple[tuple[float, float, float, float], ...] | None = None,
) -> ModelComparison:
    """Compare fits with BIC under the declared known-variance Gaussian model.

    For positive declared noise sigma, the Gaussian log likelihood has known
    variance. Terms shared by models are omitted, yielding
    ``BIC = RSS / sigma**2 + k*log(n)``. With noise-free or undeclared-scale
    observations the criterion is unavailable; a lower training residual alone
    never selects the more complex model.
    """
    single = estimate_exponential_rate(observations, method="nonlinear")
    two = fit_two_mode_exponential(observations, initial_guesses=initial_guesses)
    single_bic: float | None = None
    two_bic: float | None = None
    delta: float | None = None
    selection = ModelSelection.UNAVAILABLE
    assumptions = (
        "BIC comparison assumes independent additive Gaussian errors with the "
        "positive observation noise standard deviation declared in the input; "
        "the variance is treated as known. BIC is not a mode-identifiability "
        "proof, and no model is selected from residual reduction alone."
    )
    if (
        observations.noise_standard_deviation > 0.0
        and single.fit_converged
        and two.converged
    ):
        sigma = observations.noise_standard_deviation
        n = len(observations.values)
        single_rss = sum(
            (actual - fitted) ** 2
            for actual, fitted in zip(
                observations.values, single.fitted_values, strict=True
            )
        )
        two_rss = sum(
            (actual - fitted) ** 2
            for actual, fitted in zip(
                observations.values, two.fitted_values, strict=True
            )
        )
        single_bic = single_rss / (sigma * sigma) + 2 * np.log(n)
        two_bic = two_rss / (sigma * sigma) + 4 * np.log(n)
        delta = float(single_bic - two_bic)
        if delta >= BIC_MODEL_THRESHOLD:
            selection = ModelSelection.TWO
        elif delta <= -BIC_MODEL_THRESHOLD:
            selection = ModelSelection.SINGLE
        else:
            selection = ModelSelection.INCONCLUSIVE
    elif not single.fit_converged or not two.converged:
        selection = ModelSelection.UNAVAILABLE
        assumptions = (
            "BIC unavailable because at least one model fit did not converge; "
            "model selection is not inferred from the remaining fit alone."
        )
    elif observations.noise_standard_deviation <= 0.0:
        assumptions = (
            "BIC unavailable because no positive observation noise variance was "
            "declared; residual ranking alone does not select a model."
        )
    return ModelComparison(
        single_mode_fit=single,
        two_mode_fit=two,
        single_bic=None if single_bic is None else float(single_bic),
        two_bic=None if two_bic is None else float(two_bic),
        bic_difference_single_minus_two=delta,
        selected_model=selection,
        comparison_assumption=assumptions,
    )


def evaluate_mode_identification(
    comparison: ModelComparison, truth: SyntheticGroundTruth
) -> IdentificationEvaluation:
    """Compare fitted results to truth only after estimator/model selection."""
    true_count = len(truth.component_rates_per_s)
    fit = comparison.two_mode_fit
    rate_errors: tuple[MatchedModeError, ...] = ()
    amplitude_errors: tuple[float, ...] = ()
    if fit.converged and fit.fitted_rates_per_s is not None:
        estimated_rates = tuple(sorted(fit.fitted_rates_per_s))
        if true_count == 2:
            true_rates = tuple(sorted(truth.component_rates_per_s))
            rate_errors = tuple(
                MatchedModeError(
                    true_rate_per_s=actual,
                    estimated_rate_per_s=estimate,
                    absolute_error_per_s=abs(estimate - actual),
                    relative_error=abs(estimate - actual) / actual,
                )
                for actual, estimate in zip(true_rates, estimated_rates, strict=True)
            )
            estimated_pairs = sorted(
                zip(
                    fit.fitted_rates_per_s,
                    fit.fitted_amplitudes or (),
                    strict=True,
                )
            )
            truth_pairs = sorted(
                zip(
                    truth.component_rates_per_s,
                    truth.component_amplitudes,
                    strict=True,
                )
            )
            amplitude_errors = tuple(
                abs(estimated[1] - actual[1])
                for actual, estimated in zip(truth_pairs, estimated_pairs, strict=True)
            )
    if not fit.converged and "at least eight samples" in fit.convergence_message:
        status = IdentificationStatus.INSUFFICIENT
        interpretation = "insufficient observations for a two-mode fit"
    elif not comparison.single_mode_fit.fit_converged:
        status = IdentificationStatus.FIT_FAILED
        interpretation = "single-mode reference fit did not converge"
    elif not fit.converged:
        status = IdentificationStatus.FIT_FAILED
        interpretation = "all two-mode fitting initializations failed"
    elif comparison.selected_model is ModelSelection.UNAVAILABLE:
        status = IdentificationStatus.INCONCLUSIVE
        interpretation = (
            "model comparison is unavailable without a declared positive noise "
            "scale; residual ranking alone does not select a mode count"
        )
    elif comparison.selected_model is ModelSelection.INCONCLUSIVE:
        status = IdentificationStatus.INCONCLUSIVE
        interpretation = (
            "BIC difference did not reach the predeclared selection threshold"
        )
    elif comparison.selected_model is ModelSelection.SINGLE:
        if true_count == 1:
            single_error = (
                abs(
                    (comparison.single_mode_fit.estimated_rate_per_s or 0.0)
                    - truth.component_rates_per_s[0]
                )
                / truth.component_rates_per_s[0]
            )
            if single_error <= RATE_RECOVERY_RELATIVE_TOLERANCE:
                status = IdentificationStatus.SUPPORTED
                interpretation = (
                    "single-mode selection agrees with synthetic truth and rate "
                    "recovery tolerance"
                )
            else:
                status = IdentificationStatus.INCONSISTENT
                interpretation = (
                    "single-mode selected but its rate misses the recovery tolerance"
                )
        else:
            status = IdentificationStatus.INCONSISTENT
            interpretation = (
                "single-mode selection misses the known second synthetic mode"
            )
    elif true_count == 1:
        status = IdentificationStatus.INCONSISTENT
        interpretation = (
            "two-mode selection overfits a known single-mode synthetic signal"
        )
    elif all(
        item.relative_error <= RATE_RECOVERY_RELATIVE_TOLERANCE for item in rate_errors
    ):
        status = IdentificationStatus.SUPPORTED
        interpretation = (
            "two-mode selection and both matched rates pass the predeclared "
            "10% descriptive recovery tolerance"
        )
    else:
        status = IdentificationStatus.INCONSISTENT
        interpretation = (
            "two-mode model selected, but at least one matched rate misses the "
            "predeclared recovery tolerance"
        )
    if (
        comparison.selected_model is ModelSelection.SINGLE
        and true_count == 1
        and comparison.single_mode_fit.estimated_rate_per_s is not None
    ):
        true_rate = truth.component_rates_per_s[0]
        estimate = comparison.single_mode_fit.estimated_rate_per_s
        rate_errors = (
            MatchedModeError(
                true_rate_per_s=true_rate,
                estimated_rate_per_s=estimate,
                absolute_error_per_s=abs(estimate - true_rate),
                relative_error=abs(estimate - true_rate) / true_rate,
            ),
        )
        amplitude_errors = (
            abs(
                (comparison.single_mode_fit.estimated_amplitude or 0.0)
                - truth.component_amplitudes[0]
            ),
        )
    condition = fit.jacobian_condition_number
    if (
        status is IdentificationStatus.SUPPORTED
        and comparison.selected_model is ModelSelection.TWO
        and condition is not None
        and condition > 1e8
    ):
        status = IdentificationStatus.UNRELIABLE
        interpretation = (
            "fit is poorly conditioned despite passing rate recovery checks"
        )
    return IdentificationEvaluation(
        status=status,
        true_mode_count=true_count,
        selected_model=comparison.selected_model,
        matched_rate_errors=rate_errors,
        matched_amplitude_errors=amplitude_errors,
        interpretation=interpretation,
    )


def _run_case(
    case: ModeObservationCase,
    *,
    design_name: str,
    initial_guess_design: str = "observation-span-derived five-start grid",
    initial_guesses: tuple[tuple[float, float, float, float], ...] | None = None,
) -> IdentificationRun:
    comparison = compare_single_and_two_mode_models(
        case.observations, initial_guesses=initial_guesses
    )
    evaluation = evaluate_mode_identification(comparison, case.ground_truth)
    return IdentificationRun(
        case=case,
        design_name=design_name,
        initial_guess_design=initial_guess_design,
        comparison=comparison,
        evaluation=evaluation,
    )


def _regular_times(duration_s: float, count: int) -> tuple[float, ...]:
    if duration_s <= 0.0 or count < 3:
        raise ScientificValidationError(
            "duration must be positive and count at least three"
        )
    return tuple(float(value) for value in np.linspace(0.0, duration_s, count))


def generate_irregular_times(
    *,
    duration_s: float,
    count: int,
    design: Literal["mild_jitter", "strong_jitter", "clustered", "missing_middle"],
    seed: int,
) -> tuple[float, ...]:
    """Create deterministic irregular times with fixed endpoints and order."""
    if design not in ("mild_jitter", "strong_jitter", "clustered", "missing_middle"):
        raise ScientificValidationError("unknown irregular sampling design")
    regular = np.asarray(_regular_times(duration_s, count), dtype=np.float64)
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ScientificValidationError("seed must be a non-negative integer")
    if design == "clustered":
        early_count = max(3, int(round((count - 2) * 0.6)))
        late_count = count - early_count
        if late_count < 3:
            raise ScientificValidationError(
                "clustered design needs at least six samples"
            )
        early = np.linspace(0.0, 0.2 * duration_s, early_count, endpoint=False)
        late = np.linspace(0.8 * duration_s, duration_s, late_count)
        times = np.concatenate((early, late))
    elif design == "missing_middle":
        times = regular[(regular < 0.4 * duration_s) | (regular > 0.6 * duration_s)]
        times = np.unique(np.concatenate(([0.0], times, [duration_s])))
    else:
        rng = np.random.default_rng(seed)
        step = duration_s / (count - 1)
        fraction = 0.2 if design == "mild_jitter" else 0.4
        jitter = rng.uniform(-fraction * step, fraction * step, size=count)
        jitter[0] = 0.0
        jitter[-1] = 0.0
        times = np.sort(regular + jitter)
    validated = _validated_times(tuple(float(value) for value in times))
    return validated


def _coverage(name: str, times: tuple[float, ...]) -> SamplingCoverage:
    intervals = np.diff(np.asarray(times, dtype=np.float64))
    return SamplingCoverage(
        design_name=name,
        sample_count=len(times),
        duration_s=times[-1] - times[0],
        minimum_interval_s=float(np.min(intervals)),
        median_interval_s=float(np.median(intervals)),
        maximum_interval_s=float(np.max(intervals)),
    )


def _summarize_repetitions(
    condition: str, runs: tuple[IdentificationRun, ...], seeds: tuple[int, ...]
) -> RepeatedNoiseSummary:
    errors_by_mode: dict[int, list[MatchedModeError]] = {}
    for run in runs:
        for index, error in enumerate(run.evaluation.matched_rate_errors):
            errors_by_mode.setdefault(index + 1, []).append(error)
    summaries = tuple(
        ModeRateDistribution(
            mode_name=f"sorted_mode_{mode}",
            estimate_median_per_s=float(
                np.median([item.estimated_rate_per_s for item in values])
            ),
            estimate_minimum_per_s=float(
                np.min([item.estimated_rate_per_s for item in values])
            ),
            estimate_maximum_per_s=float(
                np.max([item.estimated_rate_per_s for item in values])
            ),
            estimate_population_sd_per_s=float(
                np.std([item.estimated_rate_per_s for item in values], ddof=0)
            ),
            signed_error_mean_per_s=float(
                np.mean(
                    [
                        item.estimated_rate_per_s - item.true_rate_per_s
                        for item in values
                    ]
                )
            ),
            absolute_error_median_per_s=float(
                np.median([item.absolute_error_per_s for item in values])
            ),
            absolute_error_minimum_per_s=float(
                np.min([item.absolute_error_per_s for item in values])
            ),
            absolute_error_maximum_per_s=float(
                np.max([item.absolute_error_per_s for item in values])
            ),
            absolute_error_population_sd_per_s=float(
                np.std([item.absolute_error_per_s for item in values], ddof=0)
            ),
            relative_error_median=float(
                np.median([item.relative_error for item in values])
            ),
            relative_error_minimum=float(
                np.min([item.relative_error for item in values])
            ),
            relative_error_maximum=float(
                np.max([item.relative_error for item in values])
            ),
            relative_error_population_sd=float(
                np.std([item.relative_error for item in values], ddof=0)
            ),
        )
        for mode, values in sorted(errors_by_mode.items())
        if values
    )
    statuses = [run.evaluation.status for run in runs]
    return RepeatedNoiseSummary(
        condition_name=condition,
        realization_count=len(runs),
        successful_fit_count=sum(
            run.comparison.single_mode_fit.fit_converged
            and run.comparison.two_mode_fit.converged
            for run in runs
        ),
        failed_fit_count=sum(
            not run.comparison.single_mode_fit.fit_converged
            or not run.comparison.two_mode_fit.converged
            for run in runs
        ),
        unreliable_count=statuses.count(IdentificationStatus.UNRELIABLE),
        inconclusive_count=statuses.count(IdentificationStatus.INCONCLUSIVE),
        supported_count=statuses.count(IdentificationStatus.SUPPORTED),
        inconsistent_count=statuses.count(IdentificationStatus.INCONSISTENT),
        per_mode_summaries=summaries,
        seeds=seeds,
        interpretation=(
            "Descriptive summaries of these fixed seeds and this design only; "
            "counts are not universal probabilities or confidence intervals."
        ),
    )


def _run_oscillator_timing() -> tuple[
    tuple[tuple[str, str, float | None, int, float], ...],
    tuple[SamplingCoverage, ...],
]:
    duration = 10.0
    count = 81
    regular = _regular_times(duration, count)
    designs: tuple[tuple[str, tuple[float, ...]], ...] = (
        ("regular", regular),
        (
            "mild_jitter",
            generate_irregular_times(
                duration_s=duration, count=count, design="mild_jitter", seed=2401
            ),
        ),
        (
            "clustered",
            generate_irregular_times(
                duration_s=duration, count=count, design="clustered", seed=2402
            ),
        ),
        (
            "missing_middle",
            generate_irregular_times(
                duration_s=duration, count=count, design="missing_middle", seed=2403
            ),
        ),
    )
    outcomes: list[tuple[str, str, float | None, int, float]] = []
    coverages: list[SamplingCoverage] = []
    for name, times in designs:
        if name == "missing_middle":
            matched = _regular_times(duration, len(times))
            matched_case = generate_oscillator_observations(
                scenario_id="phase24_oscillator_regular_matched_missing_middle",
                mass_kg=1.0,
                damping_coefficient_kg_per_s=0.1,
                stiffness_n_per_m=4.0,
                initial_displacement_m=0.1,
                initial_velocity_m_per_s=0.0,
                times_s=matched,
            )
            matched_estimate = estimate_oscillator_envelope(matched_case.observations)
            outcomes.append(
                (
                    "regular_matched_missing_middle",
                    matched_estimate.status.value,
                    matched_estimate.estimated_rate_per_s,
                    len(matched),
                    matched[-1] - matched[0],
                )
            )
            coverages.append(
                _coverage("oscillator_regular_matched_missing_middle", matched)
            )
        case = generate_oscillator_observations(
            scenario_id=f"phase24_oscillator_{name}",
            mass_kg=1.0,
            damping_coefficient_kg_per_s=0.1,
            stiffness_n_per_m=4.0,
            initial_displacement_m=0.1,
            initial_velocity_m_per_s=0.0,
            times_s=times,
        )
        estimate = estimate_oscillator_envelope(case.observations)
        outcomes.append(
            (
                name,
                estimate.status.value,
                estimate.estimated_rate_per_s,
                len(times),
                times[-1] - times[0],
            )
        )
        coverages.append(_coverage(f"oscillator_{name}", times))
    for count in (7, 13):
        sparse_designs: tuple[tuple[str, tuple[float, ...]], ...] = (
            ("regular", _regular_times(duration, count)),
            (
                "mild_jitter",
                generate_irregular_times(
                    duration_s=duration,
                    count=count,
                    design="mild_jitter",
                    seed=2490 + count,
                ),
            ),
            (
                "strong_jitter",
                generate_irregular_times(
                    duration_s=duration,
                    count=count,
                    design="strong_jitter",
                    seed=2500 + count,
                ),
            ),
            (
                "clustered",
                generate_irregular_times(
                    duration_s=duration,
                    count=count,
                    design="clustered",
                    seed=2510 + count,
                ),
            ),
        )
        for name, times in sparse_designs:
            case = generate_oscillator_observations(
                scenario_id=f"phase24_oscillator_sparse_{count}_{name}",
                mass_kg=1.0,
                damping_coefficient_kg_per_s=0.1,
                stiffness_n_per_m=4.0,
                initial_displacement_m=0.1,
                initial_velocity_m_per_s=0.0,
                times_s=times,
            )
            estimate = estimate_oscillator_envelope(case.observations)
            outcomes.append(
                (
                    f"sparse_{count}_{name}",
                    estimate.status.value,
                    estimate.estimated_rate_per_s,
                    len(times),
                    times[-1] - times[0],
                )
            )
            coverages.append(_coverage(f"oscillator_sparse_{count}_{name}", times))
    return tuple(outcomes), tuple(coverages)


def run_phase24_research() -> Phase24ResearchResults:
    """Run fixed-seed robustness, timing, and multi-mode identification cases."""
    duration = 15.0
    regular = _regular_times(duration, 61)
    repeated_seeds = tuple(range(2401, 2413))
    repeated_designs = (
        ("single_mode_control", (1.0,), (0.2,), 0.015),
        ("well_separated_comparable", (1.0, 0.6), (0.2, 0.8), 0.015),
        ("moderately_separated", (1.0, 0.8), (0.2, 0.4), 0.015),
        ("close_rates", (1.0, 1.0), (0.2, 0.24), 0.015),
        ("weak_fast_mode", (1.0, 0.03), (0.2, 0.8), 0.015),
    )
    repeated_runs: list[IdentificationRun] = []
    repeated_summaries: list[RepeatedNoiseSummary] = []
    for name, amplitudes, rates, noise in repeated_designs:
        condition_runs: list[IdentificationRun] = []
        for seed in repeated_seeds:
            case = generate_mode_observations(
                scenario_id=f"repeat_{name}_{seed}",
                amplitudes=amplitudes,
                rates_per_s=rates,
                times_s=regular,
                noise_standard_deviation=noise,
                noise_seed=seed,
            )
            run = _run_case(case, design_name=f"regular repeated noise: {name}")
            condition_runs.append(run)
            repeated_runs.append(run)
        repeated_summaries.append(
            _summarize_repetitions(name, tuple(condition_runs), repeated_seeds)
        )

    matrix_specs = (
        ("single_mode_control", (1.0,), (0.2,), 15.0, 61, 0.015, 2420),
        ("well_separated", (1.0, 0.6), (0.2, 0.8), 15.0, 61, 0.015, 2421),
        ("moderate_separation", (1.0, 0.8), (0.2, 0.4), 15.0, 61, 0.015, 2422),
        ("close_separation", (1.0, 1.0), (0.2, 0.24), 15.0, 61, 0.015, 2423),
        ("weak_fast_amplitude", (1.0, 0.03), (0.2, 0.8), 15.0, 61, 0.015, 2424),
        (
            "fast_mode_crosses_noise_floor",
            (1.0, 0.1),
            (0.2, 0.8),
            15.0,
            61,
            0.015,
            2432,
        ),
        ("strong_fast_amplitude", (0.1, 1.0), (0.2, 0.8), 15.0, 61, 0.015, 2425),
        ("short_window", (1.0, 0.6), (0.2, 0.8), 3.0, 13, 0.015, 2426),
        ("too_few_samples", (1.0, 0.6), (0.2, 0.8), 3.0, 7, 0.015, 2431),
        ("sparse_sampling", (1.0, 0.6), (0.2, 0.8), 15.0, 16, 0.015, 2427),
        ("low_noise", (1.0, 0.6), (0.2, 0.8), 15.0, 61, 0.003, 2428),
        ("moderate_noise", (1.0, 0.6), (0.2, 0.8), 15.0, 61, 0.05, 2429),
        ("noise_free_overfit_check", (1.0, 0.6), (0.2, 0.8), 15.0, 61, 0.0, 2430),
    )
    matrix_runs = tuple(
        _run_case(
            generate_mode_observations(
                scenario_id=f"matrix_{name}",
                amplitudes=amplitudes,
                rates_per_s=rates,
                times_s=_regular_times(window, count),
                noise_standard_deviation=noise,
                noise_seed=seed,
            ),
            design_name=name,
        )
        for name, amplitudes, rates, window, count, noise, seed in matrix_specs
    )

    timing_specs: tuple[tuple[str, tuple[float, ...]], ...] = (
        ("regular_baseline", regular),
        (
            "mild_jitter",
            generate_irregular_times(
                duration_s=duration, count=61, design="mild_jitter", seed=2441
            ),
        ),
        (
            "strong_jitter",
            generate_irregular_times(
                duration_s=duration, count=61, design="strong_jitter", seed=2442
            ),
        ),
        (
            "clustered_long_gap",
            generate_irregular_times(
                duration_s=duration, count=61, design="clustered", seed=2443
            ),
        ),
        (
            "missing_middle",
            generate_irregular_times(
                duration_s=duration, count=61, design="missing_middle", seed=2444
            ),
        ),
    )
    timing_runs: list[IdentificationRun] = []
    coverage: list[SamplingCoverage] = []
    for name, times in timing_specs:
        if name != "regular_baseline":
            matched_times = _regular_times(duration, len(times))
            matched_case = generate_mode_observations(
                scenario_id=f"timing_regular_matched_{name}",
                amplitudes=(1.0, 0.6),
                rates_per_s=(0.2, 0.8),
                times_s=matched_times,
                noise_standard_deviation=0.015,
                noise_seed=2440,
            )
            matched_name = f"regular_matched_{name}"
            timing_runs.append(_run_case(matched_case, design_name=matched_name))
            coverage.append(_coverage(matched_name, matched_times))
        case = generate_mode_observations(
            scenario_id=f"timing_{name}",
            amplitudes=(1.0, 0.6),
            rates_per_s=(0.2, 0.8),
            times_s=times,
            noise_standard_deviation=0.015,
            noise_seed=2440,
        )
        timing_runs.append(_run_case(case, design_name=name))
        coverage.append(_coverage(name, times))

    alternative_starts = (
        (0.9, 0.12, 0.7, 0.55),
        (0.7, 0.55, 0.9, 0.12),
        (1.0, 0.3, 0.5, 0.9),
        (0.4, 0.02, 1.1, 0.4),
    )
    start_case = generate_mode_observations(
        scenario_id="initial_guess_robustness",
        amplitudes=(1.0, 0.6),
        rates_per_s=(0.18, 0.65),
        times_s=_regular_times(25.0, 101),
        noise_standard_deviation=0.01,
        noise_seed=2450,
    )
    start_fit = fit_two_mode_exponential(
        start_case.observations, initial_guesses=alternative_starts
    )
    start_comparison = compare_single_and_two_mode_models(
        start_case.observations, initial_guesses=alternative_starts
    )
    start_run = IdentificationRun(
        case=start_case,
        design_name="explicit alternative initializations",
        initial_guess_design="four listed starts; all receive observations only",
        comparison=start_comparison,
        evaluation=evaluate_mode_identification(
            start_comparison, start_case.ground_truth
        ),
    )
    matrix_runs = matrix_runs + (start_run,)
    guess_comparison = tuple(
        (
            record.initial_amplitudes,
            record.initial_rates_per_s,
            record.converged,
            record.convergence_message,
            (
                tuple(sorted(record.fitted_rates_per_s))
                if record.fitted_rates_per_s is not None
                else (None, None)
            ),
        )
        for record in start_fit.start_records
    )
    oscillator_outcomes, oscillator_coverage = _run_oscillator_timing()
    return Phase24ResearchResults(
        repeated_noise_runs=tuple(repeated_runs),
        repeated_noise_summaries=tuple(repeated_summaries),
        mode_matrix_runs=matrix_runs,
        sampling_design_runs=tuple(timing_runs),
        sampling_coverages=tuple(coverage) + oscillator_coverage,
        oscillator_sampling_summary=oscillator_outcomes,
        initial_guess_comparison=guess_comparison,
        selection_assumptions=(
            "BIC uses the iid additive Gaussian likelihood with known positive "
            "noise variance declared in ObservationSeries; only shared constants "
            "are omitted.",
            f"Two-mode selection requires BIC(single)-BIC(two) >= "
            f"{BIC_MODEL_THRESHOLD:g}; a difference <= -{BIC_MODEL_THRESHOLD:g} "
            "selects one mode; intermediate values are inconclusive.",
            "The per-rate recovery tolerance is 10% and is used only after fitting "
            "to summarize agreement with synthetic truth.",
            "Noise-free cases have no BIC selection; lower residual in a "
            "four-parameter "
            "fit is reported but cannot establish an additional mode.",
            "Rate-separation Delta=abs(r2-r1)/max(r1,r2) is descriptive, not a "
            "universal threshold.",
        ),
    )


def render_phase24_results(results: Phase24ResearchResults) -> str:
    """Render deterministic compact result tables for the Phase 24 report."""
    lines = [
        "# Phase 24 computed modal-identification results",
        "",
        "## Repeated noise realizations",
        "",
        "| Condition | Runs | Both fits successful | Failed | Unreliable | "
        "Inconclusive | Supported | Inconsistent | Mode estimate and abs/rel "
        "error summaries (median, min, max, population SD) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for summary in results.repeated_noise_summaries:
        error_text = (
            "; ".join(
                f"{item.mode_name}: rate "
                f"{item.estimate_median_per_s:.4g}/{item.estimate_minimum_per_s:.4g}/"
                f"{item.estimate_maximum_per_s:.4g}/"
                f"{item.estimate_population_sd_per_s:.3g}; signed bias "
                f"{item.signed_error_mean_per_s:.3g}; abs err "
                f"{item.absolute_error_median_per_s:.3g}/"
                f"{item.absolute_error_minimum_per_s:.3g}/"
                f"{item.absolute_error_maximum_per_s:.3g}/"
                f"{item.absolute_error_population_sd_per_s:.3g}; rel err "
                f"{item.relative_error_median:.3g}/"
                f"{item.relative_error_minimum:.3g}/"
                f"{item.relative_error_maximum:.3g}/"
                f"{item.relative_error_population_sd:.3g}"
                for item in summary.per_mode_summaries
            )
            or "no matched rate estimates"
        )
        lines.append(
            f"| {summary.condition_name} | {summary.realization_count} | "
            f"{summary.successful_fit_count} | {summary.failed_fit_count} | "
            f"{summary.unreliable_count} | {summary.inconclusive_count} | "
            f"{summary.supported_count} | "
            f"{summary.inconsistent_count} | {error_text} |"
        )
        lines.append(f"  Seeds: `{summary.seeds}`. {summary.interpretation}")
    lines.extend(
        (
            "",
            "## Mode and sampling matrix",
            "",
            "| Design | True modes | N | Span (s) | Sigma | A true | "
            "Rates true (s^-1) | A/rate single-fit | A two-fit | "
            "Rates two-fit (s^-1) | "
            "Afast/Aslow | Delta rate | RMSE single | RMSE two | "
            "Scaled Jacobian cond. | "
            "k single/two | BIC single | BIC two | "
            "Delta BIC | Selection | Evaluation | Rate abs/rel errors | "
            "Amplitude abs errors |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"
            "---|---|---|---|---|---|---|---|---|",
        )
    )
    all_runs = results.mode_matrix_runs + results.sampling_design_runs
    for run in all_runs:
        comp = run.comparison
        fit = comp.two_mode_fit
        truth = run.case.ground_truth
        observed = run.case.observations
        fitted_rates = (
            "failed"
            if fit.fitted_rates_per_s is None
            else ", ".join(f"{rate:.5g}" for rate in fit.fitted_rates_per_s)
        )
        fitted_text = f"({fitted_rates})"
        fitted_amplitudes = (
            "failed"
            if fit.fitted_amplitudes is None
            else ", ".join(f"{value:.5g}" for value in fit.fitted_amplitudes)
        )
        true_amplitudes = ", ".join(
            f"{value:.5g}" for value in truth.component_amplitudes
        )
        true_text = ", ".join(f"{rate:.5g}" for rate in truth.component_rates_per_s)
        single_amplitude = comp.single_mode_fit.estimated_amplitude
        single_rate = comp.single_mode_fit.estimated_rate_per_s
        single_fit_text = (
            "failed"
            if single_amplitude is None or single_rate is None
            else f"{single_amplitude:.5g}/{single_rate:.5g}"
        )
        delta_rate = None
        amplitude_ratio = None
        if len(truth.component_rates_per_s) == 2:
            ordered_truth = sorted(
                zip(
                    truth.component_rates_per_s,
                    truth.component_amplitudes,
                    strict=True,
                )
            )
            r1, a1 = ordered_truth[0]
            r2, a2 = ordered_truth[1]
            delta_rate = abs(r2 - r1) / max(r1, r2)
            amplitude_ratio = a2 / a1

        def format_optional(value: float | None) -> str:
            return "n/a" if value is None else f"{value:.5g}"

        single_rmse = comp.single_mode_fit.observation_rmse
        amplitude_errors = (
            ", ".join(
                f"{value:.4g}" for value in run.evaluation.matched_amplitude_errors
            )
            or "n/a"
        )
        rate_errors = (
            "; ".join(
                f"{item.absolute_error_per_s:.4g}/{item.relative_error:.4g}"
                for item in run.evaluation.matched_rate_errors
            )
            or "n/a"
        )
        lines.append(
            f"| {run.design_name} | {len(truth.component_rates_per_s)} | "
            f"{len(observed.values)} | "
            f"{observed.times[-1] - observed.times[0]:.5g} | "
            f"{observed.noise_standard_deviation:.4g} | ({true_amplitudes}) | "
            f"({true_text}) | {single_fit_text} | "
            f"({fitted_amplitudes}) | {fitted_text} | "
            f"{format_optional(amplitude_ratio)} | "
            f"{format_optional(delta_rate)} | "
            f"{format_optional(single_rmse)} | {format_optional(fit.residual_rmse)} | "
            f"{format_optional(fit.jacobian_condition_number)} | "
            f"{comp.single_free_parameter_count}/{comp.two_free_parameter_count} | "
            f"{format_optional(comp.single_bic)} | {format_optional(comp.two_bic)} | "
            f"{format_optional(comp.bic_difference_single_minus_two)} | "
            f"{comp.selected_model.value} | {run.evaluation.status.value} | "
            f"{rate_errors} | {amplitude_errors} |"
        )
    lines.extend(("", "## Sampling coverage", ""))
    lines.append(
        "| Design | Count | Span (s) | Min dt (s) | Median dt (s) | Max dt (s) |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|")
    for item in results.sampling_coverages:
        lines.append(
            f"| {item.design_name} | {item.sample_count} | {item.duration_s:.5g} | "
            f"{item.minimum_interval_s:.5g} | {item.median_interval_s:.5g} | "
            f"{item.maximum_interval_s:.5g} |"
        )
    lines.extend(("", "## Oscillator timing comparison", ""))
    for (
        name,
        status,
        rate,
        sample_count,
        duration,
    ) in results.oscillator_sampling_summary:
        rate_text = "unavailable" if rate is None else f"{rate:.7g} s^-1"
        lines.append(
            f"- `{name}`: n={sample_count}, span={duration:g} s, {status}; "
            f"estimated envelope rate {rate_text}."
        )
    lines.extend(("", "## Initial-guess sensitivity", ""))
    for (
        amplitudes,
        initial_rates,
        converged,
        message,
        rates,
    ) in results.initial_guess_comparison:
        formatted_rates = ", ".join(
            "unavailable" if rate is None else f"{rate:.7g}" for rate in rates
        )
        lines.append(
            f"- start A=({amplitudes[0]:g}, {amplitudes[1]:g}), "
            f"rates=({initial_rates[0]:g}, {initial_rates[1]:g}) s^-1: "
            f"converged={converged}; {message}; fitted rates "
            f"{formatted_rates} s^-1."
        )
    lines.extend(("", "## Selection assumptions", ""))
    lines.extend(f"- {item}" for item in results.selection_assumptions)
    return "\n".join(lines)
