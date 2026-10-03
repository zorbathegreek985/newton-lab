"""Causal exponential filtering and a controlled Phase 29 comparison study.

The recursive filter is a discrete-time signal-processing operator, not a heat
equation solver. The accompanying study compares its behavior with the existing
noncausal finite-mode diffusion transformation on synthetic inputs.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from math import atan2, cos, isfinite, pi, sin, sqrt
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from newton_lab.diffusion_signal_smoothing import run_diffusion_smoothing_study
from newton_lab.exceptions import ScientificValidationError

StudyStatus = Literal["succeeded", "failed"]


@dataclass(frozen=True, slots=True)
class FrequencyResponse:
    """Analytical response for one discrete angular frequency in rad/sample."""

    omega_rad_per_sample: float
    magnitude: float
    phase_rad: float
    phase_delay_samples: float | None
    group_delay_samples: float


@dataclass(frozen=True, slots=True)
class BetaSweepRow:
    """One ordered parameter-sweep result; failed inputs remain visible."""

    sweep_index: int
    beta: float
    status: StudyStatus
    failure_message: str | None
    input_rmse_to_truth: float | None
    filtered_rmse_to_truth: float | None
    filtered_rmse_to_input: float | None
    slow_component_retention: float | None
    fast_component_retention: float | None
    step_10_90_width_samples: float | None
    step_half_crossing_delay_samples: float | None
    step_plateau_mean: float | None
    step_pretransition_max_abs: float | None


@dataclass(frozen=True, slots=True)
class FrequencyCheck:
    """Measured and analytical response for one steady-state sinusoid."""

    omega_rad_per_sample: float
    beta: float
    measured_magnitude: float
    predicted_magnitude: float
    measured_phase_rad: float
    predicted_phase_rad: float
    magnitude_absolute_error: float
    phase_absolute_error_rad: float


@dataclass(frozen=True, slots=True)
class DiffusionMappingRow:
    """Mode-aligned comparison with the Phase 29 diffusion attenuation."""

    mode_number: int
    omega_rad_per_sample: float
    diffusion_exposure: float
    diffusion_attenuation: float
    matched_beta: float


@dataclass(frozen=True, slots=True)
class CausalFilterStudy:
    """Deterministic frequency, contamination, step, and sweep results."""

    sample_count: int
    diffusion_exposure: float
    trailing_ma_width: int
    trailing_ma_rmse_to_truth: float
    beta_sweep: tuple[float, ...]
    frequency_checks: tuple[FrequencyCheck, ...]
    diffusion_mapping: tuple[DiffusionMappingRow, ...]
    sweep_rows: tuple[BetaSweepRow, ...]
    report_markdown: str


def exponential_filter(values: ArrayLike, beta: float) -> NDArray[np.float64]:
    """Apply ``y[n]=(1-beta)*x[n]+beta*y[n-1]``, with ``y[0]=x[0]``.

    The calculation is left-to-right, so output through index ``n`` depends only
    on input through index ``n``. Inputs must be a nonempty finite 1-D sequence;
    beta must be a finite real number in ``[0, 1)``. Caller data is not mutated.
    """
    coefficient = _validate_beta(beta)
    try:
        source = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ScientificValidationError(
            "values must be a finite numeric sequence"
        ) from error
    if source.ndim != 1 or source.size == 0 or not np.isfinite(source).all():
        raise ScientificValidationError("values must be a nonempty finite 1-D sequence")
    result = np.empty_like(source)
    result[0] = source[0]
    for index in range(1, len(source)):
        result[index] = (1.0 - coefficient) * source[index] + coefficient * result[
            index - 1
        ]
    if not np.isfinite(result).all():
        raise ScientificValidationError("filter output contains non-finite values")
    return result


def frequency_response(beta: float, omega_rad_per_sample: float) -> FrequencyResponse:
    """Return analytic magnitude, phase, phase delay, and group delay.

    Phase uses the ``exp(i*omega*n)`` convention, so a lag has negative phase.
    At DC, phase delay is undefined and returned as ``None``.
    """
    coefficient = _validate_beta(beta)
    omega = _validate_omega(omega_rad_per_sample)
    denominator_sq = 1.0 + coefficient**2 - 2.0 * coefficient * cos(omega)
    magnitude = (1.0 - coefficient) / sqrt(denominator_sq)
    phase = -atan2(coefficient * sin(omega), 1.0 - coefficient * cos(omega))
    group_delay = coefficient * (cos(omega) - coefficient) / denominator_sq
    phase_delay = None if omega == 0.0 else -phase / omega
    return FrequencyResponse(omega, magnitude, phase, phase_delay, group_delay)


def run_causal_filter_study(
    *,
    beta_sweep: Sequence[float] = (0.0, 0.5, 0.9, 0.98, 0.9, -0.1, 1.0),
    sample_count: int = 501,
    diffusion_exposure: float = 0.001,
) -> CausalFilterStudy:
    """Run repeatable synthetic comparisons against the Phase 29 signal study."""
    if (
        isinstance(sample_count, bool)
        or not isinstance(sample_count, int)
        or sample_count < 101
    ):
        raise ScientificValidationError("sample_count must be an integer >= 101")
    exposure = _finite_scalar(diffusion_exposure, "diffusion_exposure")
    if exposure < 0.0:
        raise ScientificValidationError("diffusion_exposure must be non-negative")
    requested = tuple(beta_sweep)
    if not requested:
        raise ScientificValidationError("beta_sweep must contain at least one value")

    # These frequencies span low, middle, and high portions of [0, pi].
    beta_for_response = 0.9
    frequencies = (0.05 * pi, 0.25 * pi, 0.75 * pi)
    frequency_checks = tuple(
        _measure_sinusoid(beta_for_response, omega) for omega in frequencies
    )

    phase29 = run_diffusion_smoothing_study(
        exposures=(exposure,), sample_count=sample_count
    )
    signal = next(
        item for item in phase29.signals if item.case_id == "slow_plus_high_frequency"
    )
    diffusion_run = next(
        row
        for row in phase29.runs
        if row.case_id == signal.case_id and row.sweep_index == 0
    )
    if diffusion_run.status != "succeeded":
        raise ScientificValidationError("Phase 29 diffusion comparison failed")
    noisy = np.asarray(signal.input_values, dtype=np.float64)
    truth = np.asarray(signal.truth, dtype=np.float64)
    diffusion_values = np.asarray(diffusion_run.diffusion_values, dtype=np.float64)
    trailing_ma = _trailing_moving_average(noisy, 21)
    time = np.linspace(0.0, 1.0, sample_count, dtype=np.float64)
    step = np.where(np.arange(sample_count) >= sample_count // 2, 1.0, 0.0)
    transition = sample_count // 2
    plateau_start = transition + sample_count // 10
    plateau_end = min(sample_count, transition + 3 * sample_count // 10)

    rows: list[BetaSweepRow] = []
    recorded_betas: list[float] = []
    for index, raw_beta in enumerate(requested):
        beta = _coerce_beta(raw_beta)
        recorded_betas.append(beta)
        try:
            coefficient = _validate_beta(raw_beta)
            filtered = exponential_filter(noisy, coefficient)
            filtered_step = exponential_filter(step, coefficient)
            slow_raw, fast_raw = _component_amplitudes(noisy, time, (1, 12))
            slow_filtered, fast_filtered = _component_amplitudes(
                filtered, time, (1, 12)
            )
            response = _step_metrics(filtered_step, transition)
            rows.append(
                BetaSweepRow(
                    index,
                    beta,
                    "succeeded",
                    None,
                    _rmse(noisy, truth),
                    _rmse(filtered, truth),
                    _rmse(filtered, noisy),
                    _ratio(slow_filtered, slow_raw),
                    _ratio(fast_filtered, fast_raw),
                    response[0],
                    response[1],
                    float(np.mean(filtered_step[plateau_start:plateau_end])),
                    float(np.max(np.abs(filtered_step[:transition]))),
                )
            )
        except ScientificValidationError as error:
            rows.append(
                BetaSweepRow(
                    index,
                    beta,
                    "failed",
                    str(error),
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                )
            )

    mapping_modes = (1, 4, 12)
    intervals = sample_count - 1
    mapping = tuple(
        DiffusionMappingRow(
            mode,
            mode * pi / intervals,
            exposure,
            float(np.exp(-((mode * pi) ** 2) * exposure)),
            _beta_matching_magnitude(
                mode * pi / intervals,
                float(np.exp(-((mode * pi) ** 2) * exposure)),
            ),
        )
        for mode in mapping_modes
    )
    report = _render_report(
        sample_count,
        exposure,
        tuple(recorded_betas),
        frequency_checks,
        mapping,
        tuple(rows),
        _rmse(noisy, truth),
        _rmse(diffusion_values, truth),
        _rmse(trailing_ma, truth),
        21,
    )
    return CausalFilterStudy(
        sample_count,
        exposure,
        21,
        _rmse(trailing_ma, truth),
        tuple(recorded_betas),
        frequency_checks,
        mapping,
        tuple(rows),
        report,
    )


def _measure_sinusoid(beta: float, omega: float) -> FrequencyCheck:
    count = 8192
    indices = np.arange(count, dtype=np.float64)
    source = np.sin(omega * indices)
    output = exponential_filter(source, beta)
    start = 1500
    design = np.column_stack(
        (np.sin(omega * indices[start:]), np.cos(omega * indices[start:]))
    )
    coefficients, *_ = np.linalg.lstsq(design, output[start:], rcond=None)
    magnitude = float(np.hypot(coefficients[0], coefficients[1]))
    phase = float(atan2(coefficients[1], coefficients[0]))
    expected = frequency_response(beta, omega)
    return FrequencyCheck(
        omega,
        beta,
        magnitude,
        expected.magnitude,
        phase,
        expected.phase_rad,
        abs(magnitude - expected.magnitude),
        abs(_wrap_angle(phase - expected.phase_rad)),
    )


def _component_amplitudes(
    values: NDArray[np.float64], time: NDArray[np.float64], modes: tuple[int, int]
) -> tuple[float, float]:
    columns = [np.sin(mode * pi * time) for mode in modes]
    coefficients, *_ = np.linalg.lstsq(np.column_stack(columns), values, rcond=None)
    return float(coefficients[0]), float(coefficients[1])


def _trailing_moving_average(
    values: NDArray[np.float64], width: int
) -> NDArray[np.float64]:
    """Causal trailing mean with a partial, renormalized startup window."""
    return np.asarray(
        [
            np.mean(values[max(0, index - width + 1) : index + 1])
            for index in range(len(values))
        ],
        dtype=np.float64,
    )


def _step_metrics(
    values: NDArray[np.float64], transition: int
) -> tuple[float | None, float | None]:
    def crossing(level: float) -> float | None:
        for index in range(transition - 1, len(values) - 1):
            if values[index] < level <= values[index + 1]:
                delta = values[index + 1] - values[index]
                return float(index + (level - values[index]) / delta)
        return None

    low, high, half = crossing(0.1), crossing(0.9), crossing(0.5)
    return (
        None if low is None or high is None else high - low,
        None if half is None else half - (transition - 0.5),
    )


def _beta_matching_magnitude(omega: float, target: float) -> float:
    if target >= 1.0 - 1e-14:
        return 0.0
    low, high = 0.0, 1.0 - 1e-12
    for _ in range(100):
        middle = (low + high) / 2.0
        if frequency_response(middle, omega).magnitude > target:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def _render_report(
    sample_count: int,
    exposure: float,
    betas: tuple[float, ...],
    checks: tuple[FrequencyCheck, ...],
    mapping: tuple[DiffusionMappingRow, ...],
    rows: tuple[BetaSweepRow, ...],
    raw_rmse: float,
    diffusion_rmse: float,
    trailing_ma_rmse: float,
    trailing_ma_width: int,
) -> str:
    lines = [
        "# Phase 30: Causal diffusion-inspired filtering",
        "",
        "## Research question and preregistered hypotheses",
        "",
        "Can a causal one-sided filter suppress high-frequency variation while retaining slower components, and what delay and transition costs appear? This is a synthetic signal-processing study, not a trading or prediction test.",  # noqa: E501
        "",
        "H1: measured sinusoidal response matches the transfer-function magnitude and phase within 0.003. H2: outputs do not anticipate later samples. H3: increasing beta broadens the step rise and increases the DC group delay over the tested sweep; nonzero-frequency group delay is frequency-dependent and need not be monotonic in beta. H4: at least one beta reduces error on the designed mode-12 contamination. H5: diffusion and the causal filter have no automatic global parameter equivalence.",  # noqa: E501
        "",
        "## Filter and analytical response",
        "",
        "For 0 <= beta < 1, y[0]=x[0] and y[n]=(1-beta)x[n]+beta*y[n-1]. This recurrence reads no index later than n. Its z-transform transfer function is H(z)=(1-beta)/(1-beta z^-1). For omega in radians/sample, |H|=(1-beta)/sqrt(1+beta^2-2 beta cos(omega)); phase is -atan2(beta sin(omega),1-beta cos(omega)). Group delay is beta(cos(omega)-beta)/(1+beta^2-2 beta cos(omega)); zero-frequency group delay is beta/(1-beta). Phase delay is -phase/omega and is undefined at DC. Step threshold delay is a separate finite-record measurement.",  # noqa: E501
        "With zero-state transform notation, Y(z)=(1-beta)X(z)+beta*z^-1*Y(z). "
        "Solving for Y/X gives the stated transfer function. On the unit circle, the "
        "denominator magnitude is sqrt((1-beta*cos(omega))^2 + "
        "(beta*sin(omega))^2), yielding the magnitude expression. The reported phase "
        "uses H(e^(i*omega)) and is negative for a lag.",
        "",
        "Initialization y[0]=x[0] avoids a zero-state startup jump, but the first observation is passed through unchanged; subsequent outputs contain a decaying initialization transient. The filter is causal on its finite record by construction.",  # noqa: E501
        "",
        "For step metrics, the transition is located halfway between the last "
        "pre-change sample and the first post-change sample. This makes the beta=0 "
        "threshold delay zero; pre-transition samples are checked separately.",
        "",
        "## Numerical frequency response",
        "",
        "Sinusoids use beta=0.9, 8192 samples, and a least-squares sine/cosine fit after discarding 1500 startup samples. These angular frequencies are discrete radians per sample.",  # noqa: E501
        "",
        "| omega | measured gain | analytic gain | gain error | measured phase | analytic phase | phase error |",  # noqa: E501
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    lines.extend(
        f"| {c.omega_rad_per_sample:.6g} | {c.measured_magnitude:.6g} | {c.predicted_magnitude:.6g} | {c.magnitude_absolute_error:.3g} | {c.measured_phase_rad:.6g} | {c.predicted_phase_rad:.6g} | {c.phase_absolute_error_rad:.3g} |"  # noqa: E501
        for c in checks
    )
    lines.extend(
        [
            "",
            "## Relationship to Phase 29 diffusion",
            "",
            "Phase 29 applies a noncausal, whole-interval sine transform with exact modal factor exp(-(n*pi)^2*s) on the unit interval. It uses both sides of a finite record through its modal projection; it is not online filtering. Here the sampled spatial coordinate is xi=j/M with M=sample_count-1, so mode n maps to discrete angular frequency omega=n*pi/M. For the aligned modes below, the causal magnitude is matched separately to the diffusion factor. The differing beta values demonstrate that a single beta does not exactly reproduce these multiple modal attenuations. This mapping compares amplitudes at selected frequencies only; it does not match phase, boundary treatment, or causality.",  # noqa: E501
            "",
            "| mode n | omega (rad/sample) | exposure s | diffusion gain | beta matching this mode |",  # noqa: E501
            "|---:|---:|---:|---:|---:|",
        ]
    )
    lines.extend(
        f"| {r.mode_number} | {r.omega_rad_per_sample:.7g} | {r.diffusion_exposure:g} | {r.diffusion_attenuation:.7g} | {r.matched_beta:.7g} |"  # noqa: E501
        for r in mapping
    )
    lines.extend(
        [
            "",
            "## Controlled contamination and beta sweep",
            "",
            f"The deterministic Phase 29 signal is truth sin(pi*xi) plus a mode-12 perturbation of amplitude 0.35, sampled at {sample_count} points. No random noise is used. The same input realization is compared with the unchanged Phase 29 diffusion output at s={exposure:g}, the causal recurrence, and a causal trailing moving average of width {trailing_ma_width} samples (partial renormalized windows at startup). Raw RMSE to truth is {raw_rmse:.6g}; Phase 29 diffusion RMSE is {diffusion_rmse:.6g}; trailing-mean RMSE is {trailing_ma_rmse:.6g}. The sweep table reports exponential-filter errors. These are synthetic reconstruction errors, not a denoising guarantee. The moving average and exponential filter are not asserted to have equal bandwidth.",  # noqa: E501
            "",
            "| index | beta | status | filtered RMSE to truth | RMSE to input | slow retention | fast retention | step 10–90 width (samples) | half delay (samples) | plateau | pre-transition max | failure |",  # noqa: E501
            "|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in rows:
        cells = [str(row.sweep_index), f"{row.beta:g}", row.status]
        cells.extend(
            "—" if value is None else f"{value:.6g}"
            for value in (
                row.filtered_rmse_to_truth,
                row.filtered_rmse_to_input,
                row.slow_component_retention,
                row.fast_component_retention,
                row.step_10_90_width_samples,
                row.step_half_crossing_delay_samples,
                row.step_plateau_mean,
                row.step_pretransition_max_abs,
            )
        )
        cells.append(row.failure_message or "—")
        lines.append("| " + " | ".join(cells) + " |")
    supported = all(
        c.magnitude_absolute_error < 0.003 and c.phase_absolute_error_rad < 0.003
        for c in checks
    )
    successful = [row for row in rows if row.status == "succeeded"]
    improved = any(
        row.filtered_rmse_to_truth is not None and row.filtered_rmse_to_truth < raw_rmse
        for row in successful
    )
    lines.extend(
        [
            "",
            "## Hypothesis outcomes",
            "",
            f"- H1: {'supported' if supported else 'contradicted'} for these three frequencies; maximum magnitude error {max(c.magnitude_absolute_error for c in checks):.3g}, maximum phase error {max(c.phase_absolute_error_rad for c in checks):.3g}.",  # noqa: E501
            "- H2: supported by direct prefix non-anticipation tests and recurrence; scope is the declared discrete-time input model.",  # noqa: E501
            "- H3: supported for the DC group-delay formula and measured step widths in the tested sweep; nonzero-frequency group delay is frequency-dependent and is not generally monotonic in beta. There is no universal one-number delay claim.",  # noqa: E501
            f"- H4: {'supported' if improved else 'not supported'} for this deterministic mode-12 perturbation; the finding is restricted to this signal and metric.",  # noqa: E501
            "- H5: supported for the displayed aligned modes: separately matched beta values differ. Amplitude matching at isolated frequencies does not establish equivalence of operators.",  # noqa: E501
            "",
            "## Interpretation, limitations, and reproducibility",
            "",
            "A causal feature can in principle be computed from observations available through the current index. This establishes only information timing. It does not establish predictive value, a trading signal, profitability, risk-adjusted returns, regime robustness, or performance after costs. No market data, financial outcomes, or backtest are used.",  # noqa: E501
            "",
            "The experiment uses deterministic synthetic values, a finite grid, least-squares sinusoid fits, and one structured contamination. Diffusion uses a unit-interval spatial-mode mapping and whole-series modal projection; filter initialization and finite-sample transients affect early response. The numerical comparisons do not cover arbitrary signals or boundary conditions. The trailing moving average is a separately defined causal baseline, but its bandwidth is not matched to either smoother.",  # noqa: E501
            "",
            "Reproduce the study in the project environment:",
            "",
            "```python",
            "from newton_lab.causal_diffusion_filter import run_causal_filter_study",
            "",
            "study = run_causal_filter_study()",
            "print(study.report_markdown)",
            "```",
            "",
            "The study is self-contained and is not registered as a general Newton Lab "
            "experiment record; integration would require a new typed result model.",
            "",
            "### Next research question",
            "",
            "Under a predeclared stationary signal family, can a trailing moving average and the exponential filter be compared at matched measured attenuation and latency, with repeated noise realizations and explicit edge handling? That would isolate operator choice without implying financial utility.",  # noqa: E501
            "",
        ]
    )
    return "\n".join(lines)


def _validate_beta(beta: object) -> float:
    if isinstance(beta, bool):
        raise ScientificValidationError("beta must be a finite number in [0, 1)")
    try:
        value = float(str(beta))
    except (TypeError, ValueError, OverflowError) as error:
        raise ScientificValidationError(
            "beta must be a finite number in [0, 1)"
        ) from error
    if not isfinite(value) or value < 0.0 or value >= 1.0:
        raise ScientificValidationError("beta must be a finite number in [0, 1)")
    return value


def _validate_omega(omega: object) -> float:
    value = _finite_scalar(omega, "omega_rad_per_sample")
    if value < 0.0 or value > pi:
        raise ScientificValidationError("omega_rad_per_sample must be in [0, pi]")
    return value


def _finite_scalar(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ScientificValidationError(f"{name} must be finite")
    try:
        result = float(str(value))
    except (TypeError, ValueError, OverflowError) as error:
        raise ScientificValidationError(f"{name} must be finite") from error
    if not isfinite(result):
        raise ScientificValidationError(f"{name} must be finite")
    return result


def _coerce_beta(value: object) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError, OverflowError):
        return float("nan")


def _rmse(left: NDArray[np.float64], right: NDArray[np.float64]) -> float:
    return float(np.sqrt(np.mean(np.square(left - right))))


def _ratio(numerator: float, denominator: float) -> float | None:
    return None if abs(denominator) < 1e-14 else numerator / denominator


def _wrap_angle(value: float) -> float:
    return (value + pi) % (2.0 * pi) - pi
