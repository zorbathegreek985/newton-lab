"""Controlled tests of noncausal diffusion-inspired signal smoothing.

Signal time is mapped to position on a fixed unit interval. The existing exact
finite-mode heat solution then attenuates sine components by ``exp(-(n*pi)^2*s)``.
This is a whole-series transform, not a causal online filter or a financial model.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from math import exp, isfinite, pi
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from newton_lab.exceptions import ScientificValidationError
from newton_lab.heat_diffusion_case_study import (
    HeatDiffusionModel,
    SpatialMode,
    evaluate_heat_diffusion,
)

RunStatus = Literal["succeeded", "failed"]

DEFAULT_EXPOSURES: tuple[float, ...] = (
    0.0,
    1e-5,
    5e-5,
    1e-4,
    5e-4,
    1e-3,
    2e-3,
)
MOVING_AVERAGE_WIDTH = 21
STEP_ON = 0.35
STEP_OFF = 0.70
STEP_PLATEAU = (0.45, 0.60)
STEP_PREWINDOW = (0.25, STEP_ON)


@dataclass(frozen=True, slots=True)
class SyntheticSignal:
    """One deterministic sampled input, truth and finite sine representation."""

    case_id: str
    description: str
    truth: tuple[float, ...]
    input_values: tuple[float, ...]
    modes: tuple[SpatialMode, ...]
    slow_mode: int
    fast_mode: int
    transition_on: float | None = None
    transition_off: float | None = None


@dataclass(frozen=True, slots=True)
class SmoothingRun:
    """One exposure result; invalid settings remain in their original position."""

    case_id: str
    sweep_index: int
    exposure: float
    status: RunStatus
    failure_message: str | None
    diffusion_values: tuple[float, ...]
    moving_average_values: tuple[float, ...]
    raw_rmse_to_truth: float | None
    diffusion_rmse_to_truth: float | None
    moving_average_rmse_to_truth: float | None
    diffusion_rmse_to_input: float | None
    moving_average_rmse_to_input: float | None
    projection_rmse_to_input: float | None
    slow_mode_retention: float | None
    fast_mode_retention: float | None
    moving_average_slow_mode_retention: float | None
    moving_average_fast_mode_retention: float | None
    analytical_modal_error_max: float | None
    expected_slow_attenuation: float | None
    expected_fast_attenuation: float | None
    step_rising_delay: float | None
    step_falling_delay: float | None
    step_plateau_retention: float | None
    moving_average_step_rising_delay: float | None
    moving_average_step_falling_delay: float | None
    moving_average_step_plateau_retention: float | None
    pre_transition_max_influence: float | None
    pre_transition_change_from_projection: float | None
    moving_average_pre_transition_influence: float | None
    step_rising_10_90_width: float | None
    step_falling_90_10_width: float | None
    moving_average_rising_10_90_width: float | None
    moving_average_falling_90_10_width: float | None


@dataclass(frozen=True, slots=True)
class DiffusionSmoothingStudy:
    """Reproducible synthetic cases, ordered sweep, controls and report."""

    normalized_signal_time: tuple[float, ...]
    exposures: tuple[float, ...]
    moving_average_width_samples: int
    modal_truncation: int
    signals: tuple[SyntheticSignal, ...]
    runs: tuple[SmoothingRun, ...]
    report_markdown: str


def _rmse(left: NDArray[np.float64], right: NDArray[np.float64]) -> float:
    return float(np.sqrt(np.mean(np.square(left - right))))


def _sine_coefficient(
    values: NDArray[np.float64], normalized_time: NDArray[np.float64], mode: int
) -> float:
    """Discrete sine coefficient on an endpoint-inclusive uniform grid (DST-I)."""
    intervals = len(values) - 1
    basis = np.sin(mode * pi * normalized_time[1:-1])
    return float((2.0 / intervals) * np.dot(values[1:-1], basis))


def _moving_average(values: NDArray[np.float64], width: int) -> NDArray[np.float64]:
    """Centered finite window; clip and renormalize the window at each edge."""
    half_width = width // 2
    return np.asarray(
        [
            np.mean(
                values[max(0, i - half_width) : min(len(values), i + half_width + 1)]
            )
            for i in range(len(values))
        ],
        dtype=np.float64,
    )


def _make_signals(
    normalized_time: NDArray[np.float64], modal_truncation: int
) -> tuple[SyntheticSignal, ...]:
    two_frequency_modes = (
        SpatialMode(mode_number=1, amplitude=1.0),
        SpatialMode(mode_number=4, amplitude=0.5),
    )
    two_frequency = np.sin(pi * normalized_time) + 0.5 * np.sin(
        4.0 * pi * normalized_time
    )
    slow_truth = np.sin(pi * normalized_time)
    high_frequency_noise = 0.35 * np.sin(12.0 * pi * normalized_time)
    noisy_input = slow_truth + high_frequency_noise

    coefficients = tuple(
        SpatialMode(
            mode_number=mode,
            amplitude=(
                2.0
                * (np.cos(mode * pi * STEP_ON) - np.cos(mode * pi * STEP_OFF))
                / (mode * pi)
            ),
        )
        for mode in range(1, modal_truncation + 1)
    )
    pulse_truth = np.where(
        (normalized_time >= STEP_ON) & (normalized_time < STEP_OFF), 1.0, 0.0
    )
    return (
        SyntheticSignal(
            "two_frequency",
            "Ground truth is mode 1 plus mode 4 with amplitudes 1 and 0.5.",
            tuple(float(value) for value in two_frequency),
            tuple(float(value) for value in two_frequency),
            two_frequency_modes,
            1,
            4,
        ),
        SyntheticSignal(
            "slow_plus_high_frequency",
            "Mode 1 truth plus a deterministic mode 12 perturbation of amplitude "
            "0.35; not a random-noise realization.",
            tuple(float(value) for value in slow_truth),
            tuple(float(value) for value in noisy_input),
            (
                SpatialMode(mode_number=1, amplitude=1.0),
                SpatialMode(mode_number=12, amplitude=0.35),
            ),
            1,
            12,
        ),
        SyntheticSignal(
            "abrupt_pulse",
            f"Unit pulse on [{STEP_ON:g}, {STEP_OFF:g}), represented by "
            f"{modal_truncation} sine modes.",
            tuple(float(value) for value in pulse_truth),
            tuple(float(value) for value in pulse_truth),
            coefficients,
            1,
            12,
            STEP_ON,
            STEP_OFF,
        ),
    )


def _nearest_crossing(
    values: NDArray[np.float64],
    normalized_time: NDArray[np.float64],
    threshold: float,
    interval: tuple[float, float],
    *,
    rising: bool,
) -> float | None:
    candidates: list[float] = []
    for index in range(len(values) - 1):
        left = values[index] - threshold
        right = values[index + 1] - threshold
        crossed = left <= 0.0 < right if rising else left >= 0.0 > right
        if crossed:
            fraction = -left / (right - left)
            crossing = float(
                normalized_time[index]
                + fraction * (normalized_time[index + 1] - normalized_time[index])
            )
            if interval[0] <= crossing <= interval[1]:
                candidates.append(crossing)
    if not candidates:
        return None
    reference = (interval[0] + interval[1]) / 2.0
    return min(candidates, key=lambda point: abs(point - reference))


def _step_metrics(
    values: NDArray[np.float64], normalized_time: NDArray[np.float64]
) -> tuple[float | None, float | None, float, float, float | None, float | None]:
    rising = _nearest_crossing(
        values,
        normalized_time,
        0.5,
        (STEP_ON - 0.15, STEP_ON + 0.15),
        rising=True,
    )
    falling = _nearest_crossing(
        values,
        normalized_time,
        0.5,
        (STEP_OFF - 0.15, STEP_OFF + 0.15),
        rising=False,
    )
    rising_10 = _nearest_crossing(
        values,
        normalized_time,
        0.1,
        (STEP_ON - 0.15, STEP_ON + 0.15),
        rising=True,
    )
    rising_90 = _nearest_crossing(
        values,
        normalized_time,
        0.9,
        (STEP_ON - 0.15, STEP_ON + 0.15),
        rising=True,
    )
    falling_90 = _nearest_crossing(
        values,
        normalized_time,
        0.9,
        (STEP_OFF - 0.15, STEP_OFF + 0.15),
        rising=False,
    )
    falling_10 = _nearest_crossing(
        values,
        normalized_time,
        0.1,
        (STEP_OFF - 0.15, STEP_OFF + 0.15),
        rising=False,
    )
    plateau = (normalized_time >= STEP_PLATEAU[0]) & (
        normalized_time <= STEP_PLATEAU[1]
    )
    prewindow = (normalized_time >= STEP_PREWINDOW[0]) & (
        normalized_time < STEP_PREWINDOW[1]
    )
    plateau_retention = float(np.mean(values[plateau]))
    pre_transition = float(np.max(np.abs(values[prewindow])))
    return (
        None if rising is None else rising - STEP_ON,
        None if falling is None else falling - STEP_OFF,
        plateau_retention,
        pre_transition,
        None if rising_10 is None or rising_90 is None else rising_90 - rising_10,
        None if falling_90 is None or falling_10 is None else falling_10 - falling_90,
    )


def run_diffusion_smoothing_study(
    *,
    exposures: Sequence[float] = DEFAULT_EXPOSURES,
    sample_count: int = 501,
    modal_truncation: int = 80,
    moving_average_width: int = MOVING_AVERAGE_WIDTH,
) -> DiffusionSmoothingStudy:
    """Run deterministic synthetic cases through exact diffusion and a MA baseline.

    Normalized signal time xi maps to position ``x=L*xi``. Diffusion length is
    ``L=1 m`` and diffusivity ``alpha=1 m^2/s``, so the dimensionless exposure
    ``s=alpha*t/L^2`` numerically equals ``t`` when expressed in seconds. The
    moving average is centered and noncausal too; no bandwidth equivalence is implied.

    Invalid exposure values are retained as failed rows in input order, including
    repeated values. This function never reads external or market data.
    """
    if (
        isinstance(sample_count, bool)
        or not isinstance(sample_count, int)
        or sample_count < 101
    ):
        raise ScientificValidationError("sample_count must be an integer >= 101")
    if (
        isinstance(modal_truncation, bool)
        or not isinstance(modal_truncation, int)
        or modal_truncation < 1
        or modal_truncation >= sample_count - 1
    ):
        raise ScientificValidationError(
            "modal_truncation must be positive and less than sample_count - 1"
        )
    if (
        isinstance(moving_average_width, bool)
        or not isinstance(moving_average_width, int)
        or moving_average_width < 1
        or moving_average_width % 2 == 0
        or moving_average_width > sample_count
    ):
        raise ScientificValidationError(
            "moving_average_width must be odd and in [1, sample_count]"
        )
    requested_exposures = tuple(exposures)
    if not requested_exposures:
        raise ScientificValidationError("exposures must contain at least one value")

    normalized_time = np.linspace(0.0, 1.0, sample_count, dtype=np.float64)
    positions_m = tuple(float(value) for value in normalized_time)
    signals = _make_signals(normalized_time, modal_truncation)
    runs: list[SmoothingRun] = []
    for signal in signals:
        raw = np.asarray(signal.input_values, dtype=np.float64)
        truth = np.asarray(signal.truth, dtype=np.float64)
        moving_average = _moving_average(raw, moving_average_width)
        model = HeatDiffusionModel(
            length_m=1.0,
            left_temperature_k=0.0,
            right_temperature_k=0.0,
            diffusivity_m2_per_s=1.0,
            modes=signal.modes,
        )
        initial_modal_values = np.zeros_like(normalized_time)
        for mode in signal.modes:
            initial_modal_values += mode.amplitude * np.sin(
                mode.mode_number * pi * normalized_time
            )
        projection_error = _rmse(initial_modal_values, raw)
        raw_coefficients = {
            mode: _sine_coefficient(raw, normalized_time, mode)
            for mode in (signal.slow_mode, signal.fast_mode)
        }
        raw_slow = raw_coefficients[signal.slow_mode]
        raw_fast = raw_coefficients[signal.fast_mode]
        ma_slow = _sine_coefficient(moving_average, normalized_time, signal.slow_mode)
        ma_fast = _sine_coefficient(moving_average, normalized_time, signal.fast_mode)
        ma_step = (
            _step_metrics(moving_average, normalized_time)
            if signal.transition_on is not None
            else (None, None, None, None, None, None)
        )

        for sweep_index, requested in enumerate(requested_exposures):
            exposure = _coerce_exposure(requested)
            if isinstance(requested, bool) or not isfinite(exposure) or exposure < 0.0:
                runs.append(_failed_run(signal, sweep_index, exposure, moving_average))
                continue
            try:
                evaluation = evaluate_heat_diffusion(model, positions_m, (exposure,))
                filtered = np.asarray(evaluation.temperature_k[0], dtype=np.float64)
                modal_error = max(
                    abs(factor - exp(-((mode.mode_number * pi) ** 2) * exposure))
                    for mode, factor in zip(
                        signal.modes,
                        evaluation.attenuation_factors[0],
                        strict=True,
                    )
                )
                if not np.isfinite(filtered).all() or not isfinite(modal_error):
                    raise ScientificValidationError(
                        "diffusion evaluation produced non-finite values"
                    )
                filtered_slow = _sine_coefficient(
                    filtered, normalized_time, signal.slow_mode
                )
                filtered_fast = _sine_coefficient(
                    filtered, normalized_time, signal.fast_mode
                )
                slow_retention = (
                    None if abs(raw_slow) <= 1e-14 else filtered_slow / raw_slow
                )
                fast_retention = (
                    None if abs(raw_fast) <= 1e-14 else filtered_fast / raw_fast
                )
                diffusion_step = (
                    _step_metrics(filtered, normalized_time)
                    if signal.transition_on is not None
                    else (None, None, None, None, None, None)
                )
                expected_slow = exp(-((signal.slow_mode * pi) ** 2) * exposure)
                expected_fast = exp(-((signal.fast_mode * pi) ** 2) * exposure)
                runs.append(
                    SmoothingRun(
                        signal.case_id,
                        sweep_index,
                        exposure,
                        "succeeded",
                        None,
                        tuple(float(value) for value in filtered),
                        tuple(float(value) for value in moving_average),
                        _rmse(raw, truth),
                        _rmse(filtered, truth),
                        _rmse(moving_average, truth),
                        _rmse(filtered, raw),
                        _rmse(moving_average, raw),
                        projection_error,
                        slow_retention,
                        fast_retention,
                        None if abs(raw_slow) <= 1e-14 else ma_slow / raw_slow,
                        None if abs(raw_fast) <= 1e-14 else ma_fast / raw_fast,
                        modal_error,
                        expected_slow,
                        expected_fast,
                        diffusion_step[0],
                        diffusion_step[1],
                        diffusion_step[2],
                        ma_step[0],
                        ma_step[1],
                        ma_step[2],
                        diffusion_step[3],
                        float(
                            np.max(
                                np.abs(filtered - initial_modal_values)[
                                    (normalized_time >= STEP_PREWINDOW[0])
                                    & (normalized_time < STEP_PREWINDOW[1])
                                ]
                            )
                        )
                        if signal.transition_on is not None
                        else None,
                        ma_step[3],
                        diffusion_step[4],
                        diffusion_step[5],
                        ma_step[4],
                        ma_step[5],
                    )
                )
            except (ScientificValidationError, ValueError, OverflowError) as error:
                runs.append(
                    _failed_run(
                        signal,
                        sweep_index,
                        exposure,
                        moving_average,
                        f"{type(error).__name__}: {error}",
                    )
                )

    run_tuple = tuple(runs)
    exposures_recorded = tuple(_coerce_exposure(value) for value in requested_exposures)
    report = render_smoothing_report(
        tuple(float(value) for value in normalized_time),
        exposures_recorded,
        moving_average_width,
        modal_truncation,
        signals,
        run_tuple,
    )
    return DiffusionSmoothingStudy(
        tuple(float(value) for value in normalized_time),
        exposures_recorded,
        moving_average_width,
        modal_truncation,
        signals,
        run_tuple,
        report,
    )


def _failed_run(
    signal: SyntheticSignal,
    sweep_index: int,
    exposure: float,
    moving_average: NDArray[np.float64],
    message: str = (
        "exposure must be numeric, finite, and non-negative; booleans are invalid"
    ),
) -> SmoothingRun:
    """Keep an invalid sweep item visible without fabricating metrics."""
    return SmoothingRun(
        case_id=signal.case_id,
        sweep_index=sweep_index,
        exposure=exposure,
        status="failed",
        failure_message=message,
        diffusion_values=(),
        moving_average_values=tuple(float(value) for value in moving_average),
        raw_rmse_to_truth=None,
        diffusion_rmse_to_truth=None,
        moving_average_rmse_to_truth=None,
        diffusion_rmse_to_input=None,
        moving_average_rmse_to_input=None,
        projection_rmse_to_input=None,
        slow_mode_retention=None,
        fast_mode_retention=None,
        moving_average_slow_mode_retention=None,
        moving_average_fast_mode_retention=None,
        analytical_modal_error_max=None,
        expected_slow_attenuation=None,
        expected_fast_attenuation=None,
        step_rising_delay=None,
        step_falling_delay=None,
        step_plateau_retention=None,
        moving_average_step_rising_delay=None,
        moving_average_step_falling_delay=None,
        moving_average_step_plateau_retention=None,
        pre_transition_max_influence=None,
        pre_transition_change_from_projection=None,
        moving_average_pre_transition_influence=None,
        step_rising_10_90_width=None,
        step_falling_90_10_width=None,
        moving_average_rising_10_90_width=None,
        moving_average_falling_90_10_width=None,
    )


def _coerce_exposure(value: object) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError, OverflowError):
        return float("nan")


def render_smoothing_report(
    normalized_time: tuple[float, ...],
    exposures: tuple[float, ...],
    moving_average_width: int,
    modal_truncation: int,
    signals: tuple[SyntheticSignal, ...],
    runs: tuple[SmoothingRun, ...],
) -> str:
    """Render a deterministic report from the recorded study results."""
    lines = [
        "# Phase 29: Diffusion-inspired signal smoothing",
        "",
        "## Research question and hypothesis",
        "",
        "Can transient diffusion attenuate rapidly varying synthetic components",
        "more than slow components, while making responsiveness loss measurable?",
        "The hypothesis is limited to the declared modal model and synthetic",
        "signals; it makes no claim about price prediction or trading utility.",
        "",
        "Hypotheses specified before measurement: H1, evaluated attenuation",
        "matches the closed form; H2, higher modes attenuate more; H3, removing",
        "controlled high-frequency contamination can reduce error to clean truth",
        "for that case but need not improve every signal; H4, pulse smoothing",
        "broadens the transition and creates pre-edge influence.",
        "",
        "## Mathematical mapping and assumptions",
        "",
        "For fixed Dirichlet endpoints, `du/dt = alpha*d2u/dx2` has modal",
        "solution `a_n sin(n*pi*x/L) exp(-alpha*(n*pi/L)^2*t)`. With",
        "`x=L*xi`, xi maps to normalized signal time in `[0,1]`; this is a",
        "mathematical coordinate mapping, not physical space in a market. With",
        "`s=alpha*t/L^2`, modal retention is `exp(-(n*pi)^2*s)`. The length",
        "is 1 m and alpha is 1 m^2/s. Exposure is dimensionless; its numeric",
        "value equals diffusion time in seconds for these selected units. Signal",
        "values use arbitrary units, not literal kelvin temperatures.",
        "",
        "This is a whole-interval sine-series transform. Coefficients depend on",
        "the complete interval, including observations after a given xi; the",
        "transformation is noncausal. The comparator is also noncausal: a centered",
        f"{moving_average_width}-sample moving average, with clipped windows",
        "renormalized at both endpoints. No equal-bandwidth claim is made.",
        "The heat model assumes constant positive diffusivity, no source, and",
        "fixed zero endpoint perturbations. These are modeling assumptions, not",
        "properties attributed to a financial time series.",
        "",
        "## Cases and metrics",
        "",
        "Synthetic cases:",
        "",
        "RMSE is `sqrt(mean((estimate-reference)^2))` on the shared sampled xi grid.",
        "Input RMSE is measured to the known clean truth. Change-to-input RMSE",
        "measures how much a transformation alters the raw input. Sine-mode",
        "retention is the discrete DST-I coefficient after/before ratio.",
        "Modal formula error is max absolute difference between evaluated and",
        "closed-form attenuation factors. Step delay is the nearest linearly",
        "interpolated 0.5 crossing within +/-0.15 xi of each known edge, minus",
        "that edge. Transition width is the 10-to-90% (rising) or 90-to-10%",
        "(falling) crossing distance. Plateau retention is mean output on",
        "xi=[0.45,0.60].",
        "Pre-transition influence is max absolute output on xi=[0.25,0.35),",
        "divided by unit step height. Additional diffusion influence is the max",
        "absolute change from the zero-exposure modal projection in that window.",
        "These are operational pulse metrics, not universal responsiveness measures.",
        "",
        f"Dimensionless exposures retain order and repeats: {exposures}.",
        f"Sample count={len(normalized_time)}; modal truncation={modal_truncation};",
        f"moving-average width={moving_average_width} samples.",
        "The moving average uses the same raw input realization per case and is",
        "independent of diffusion exposure. At edges, its available window is",
        "clipped and renormalized. The unsmoothed input is the third baseline.",
        "",
        "## Results",
        "",
        "| "
        + " | ".join(
            (
                "Case",
                "s",
                "status",
                "projection RMSE",
                "raw RMSE truth",
                "diffusion RMSE truth",
                "MA RMSE truth",
                "diffusion change to input",
                "MA change to input",
                "formula error",
            )
        )
        + " |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for signal in signals:
        lines.append(f"- `{signal.case_id}`: {signal.description}")
    lines.extend(
        (
            "The pulse uses a truncated sine representation; its zero-exposure",
            "projection error is reported separately from smoothing error.",
            "For pulse edges a and b, coefficients are",
            "`b_n=2*(cos(n*pi*a)-cos(n*pi*b))/(n*pi)`.",
            "",
        )
    )

    def fmt(value: float | None) -> str:
        return "n/a" if value is None else f"{value:.5g}"

    for run in runs:
        lines.append(
            "| "
            + " | ".join(
                (
                    run.case_id,
                    fmt(run.exposure),
                    run.status,
                    fmt(run.projection_rmse_to_input),
                    fmt(run.raw_rmse_to_truth),
                    fmt(run.diffusion_rmse_to_truth),
                    fmt(run.moving_average_rmse_to_truth),
                    fmt(run.diffusion_rmse_to_input),
                    fmt(run.moving_average_rmse_to_input),
                    fmt(run.analytical_modal_error_max),
                )
            )
            + " |"
        )
    lines.extend(
        (
            "",
            "Mode attenuation and retention; analytic factors are shown next to",
            "the sampled coefficient ratios. Moving-average ratios are measured on",
            "the same input and do not share the diffusion transfer function.",
            "",
            "| "
            + " | ".join(
                (
                    "Case",
                    "s",
                    "analytic slow H",
                    "diffusion slow",
                    "MA slow",
                    "analytic fast H",
                    "diffusion fast",
                    "MA fast",
                )
            )
            + " |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        )
    )
    for run in runs:
        lines.append(
            "| "
            + " | ".join(
                (
                    run.case_id,
                    fmt(run.exposure),
                    fmt(run.expected_slow_attenuation),
                    fmt(run.slow_mode_retention),
                    fmt(run.moving_average_slow_mode_retention),
                    fmt(run.expected_fast_attenuation),
                    fmt(run.fast_mode_retention),
                    fmt(run.moving_average_fast_mode_retention),
                )
            )
            + " |"
        )
    lines.extend(
        (
            "",
            "Pulse responsiveness metrics use the declared edge, threshold and",
            "windows. Pre-projection influence at zero exposure is Fourier/Gibbs",
            "ringing; the additional-change column isolates change caused by the",
            "diffusion evolution relative to that truncated initial projection.",
            "",
            "| "
            + " | ".join(
                (
                    "s",
                    "diffusion rise/fall delay",
                    "diffusion rise/fall width",
                    "MA rise/fall delay",
                    "MA rise/fall width",
                    "diffusion/MA plateau",
                    "diffusion pre / added / MA influence",
                )
            )
            + " |",
            "|---:|---|---|---|---|---|---|",
        )
    )
    for run in runs:
        if run.case_id != "abrupt_pulse":
            continue
        delay = f"{fmt(run.step_rising_delay)}/{fmt(run.step_falling_delay)}"
        width = (
            f"{fmt(run.step_rising_10_90_width)}/{fmt(run.step_falling_90_10_width)}"
        )
        ma_delay = (
            f"{fmt(run.moving_average_step_rising_delay)}/"
            f"{fmt(run.moving_average_step_falling_delay)}"
        )
        ma_width = (
            f"{fmt(run.moving_average_rising_10_90_width)}/"
            f"{fmt(run.moving_average_falling_90_10_width)}"
        )
        plateau = (
            f"{fmt(run.step_plateau_retention)}/"
            f"{fmt(run.moving_average_step_plateau_retention)}"
        )
        influence = (
            f"{fmt(run.pre_transition_max_influence)}/"
            f"{fmt(run.pre_transition_change_from_projection)}/"
            f"{fmt(run.moving_average_pre_transition_influence)}"
        )
        lines.append(
            "| "
            + " | ".join(
                (
                    fmt(run.exposure),
                    delay,
                    width,
                    ma_delay,
                    ma_width,
                    plateau,
                    influence,
                )
            )
            + " |"
        )
    valid_runs = [run for run in runs if run.status == "succeeded"]
    formula_supported = bool(valid_runs) and all(
        run.analytical_modal_error_max is not None
        and run.analytical_modal_error_max <= 1e-14
        for run in valid_runs
    )
    positive_runs = [run for run in valid_runs if run.exposure > 0.0]
    order_supported = bool(positive_runs) and all(
        run.expected_slow_attenuation is not None
        and run.expected_fast_attenuation is not None
        and run.expected_fast_attenuation < run.expected_slow_attenuation
        for run in positive_runs
    )
    noise_improvement = any(
        run.raw_rmse_to_truth is not None
        and run.diffusion_rmse_to_truth is not None
        and run.diffusion_rmse_to_truth < run.raw_rmse_to_truth
        for run in runs
        if run.case_id == "slow_plus_high_frequency"
    )
    noise_tested = any(
        run.exposure > 0.0 for run in runs if run.case_id == "slow_plus_high_frequency"
    )
    signal_worsening = any(
        run.exposure > 0.0
        and run.raw_rmse_to_truth is not None
        and run.diffusion_rmse_to_truth is not None
        and run.diffusion_rmse_to_truth > run.raw_rmse_to_truth
        for run in runs
        if run.case_id == "two_frequency"
    )
    widths = [
        run.step_rising_10_90_width
        for run in runs
        if run.case_id == "abrupt_pulse" and run.step_rising_10_90_width is not None
    ]
    pulse_broadening = len(widths) > 1 and max(widths) > min(widths)
    h1 = "supported" if formula_supported else "inconclusive"
    h2 = "supported" if order_supported else "inconclusive"
    h3 = (
        "supported for the controlled high-frequency case"
        if noise_improvement
        else "not supported for its tested positive exposures"
        if noise_tested
        else "inconclusive; no positive exposure was tested"
    )
    h4 = (
        "supported by increased measured rise width"
        if pulse_broadening
        else "inconclusive"
    )
    lines.extend(
        (
            "",
            "## Hypothesis assessment for this sweep",
            "",
            f"- H1 attenuation formula agreement: {h1}.",
            "  Successful runs have maximum formula error <= 1e-14.",
            f"- H2 higher-mode attenuation ordering: {h2} for positive exposures.",
            f"- H3 reconstruction improvement: {h3}.",
            f"  Two-frequency detail loss observed: {signal_worsening}.",
            f"- H4 pulse responsiveness cost: {h4}.",
            "  Midpoint delay and pre-edge effects depend on pulse and truncation.",
        )
    )
    succeeded = sum(run.status == "succeeded" for run in runs)
    failed = len(runs) - succeeded
    lines.extend(
        (
            "",
            f"Recorded rows: {len(runs)}; succeeded: {succeeded}; failed: {failed}.",
            "Failures retain their sweep index and exposure; failed metrics remain",
            "undefined rather than being replaced with zero.",
            "",
            "## Findings by evidence class",
            "",
            "**Established mathematical result.** For this fixed-boundary model,",
            "mode n has attenuation `exp(-(n*pi)^2*s)`, so larger n has a faster",
            "decay rate. Agreement between the evaluator and this expression checks",
            "the modal scaling and mode indexing. This residual is an algebraic",
            "cross-check, not an independent PDE discretization or convergence test.",
            "",
            "**Experimentally supported result.** The table reports measured modal",
            "retention, clean-signal error, input change, moving-average comparison,",
            "and pulse metrics only for the finite cases and exposure values shown.",
            "Lower high-frequency energy does not by itself imply lower error to",
            "truth; compare both metrics for each exposure.",
            "",
            "**Candidate cross-domain application.** Diffusion-inspired smoothing is",
            "a candidate whole-series preprocessing or feature-construction method.",
            "Smoothing modifies the input and cannot restore the original signal",
            "automatically. The centered transform can let later values influence",
            "earlier outputs, as quantified by the pulse pre-influence metric.",
            "",
            "**Unvalidated financial application.** No financial observations were",
            "downloaded or tested. This study establishes no predictive power,",
            "trading signal, profitability, excess returns, risk-adjusted performance,",
            "asset/regime robustness, or executable performance after costs. A",
            "causal redesign and separate out-of-sample validation would be required",
            "before any real-time financial feature could be assessed.",
            "",
            "## Newton Lab discovery integration",
            "",
            "Phase 19 and the existing knowledge integration already record the",
            "diffusion-to-synthetic-signal relationship as a synthetic analogy.",
            "Those APIs do not consume this study's per-exposure result rows; changing",
            "their artifact schema would be broader than this focused experiment.",
            "This report therefore preserves the evidence classes locally rather",
            "than promoting the measurements into a registry relationship.",
            "",
            "## Reproducibility and limitations",
            "",
            "The synthetic inputs, coefficients, sample grid, moving-average window,",
            "modal cutoff and ordered exposures are deterministic. The exact modal",
            "evaluator is reused; no finite-difference PDE solver is introduced.",
            "The modal formula residual checks the dimensionless scaling against the",
            "existing analytical evaluator; this is not a solver accuracy bound.",
            "For the pulse, finite truncation causes Gibbs/projection error even at",
            "zero exposure. The discrete DST coefficient and sampled metrics have",
            "finite-grid effects. The moving average and diffusion operator have",
            "different transfer functions and boundary behavior, so this is not an",
            "equal-bandwidth comparison. Neither method is a causal feature here.",
            "",
            "## Conclusion and next question",
            "",
            "The study can test attenuation ordering and synthetic trade-offs under",
            "a declared mathematical mapping. It cannot establish that smoothing",
            "improves real data or predicts future outcomes. A next scientific step",
            "would be to derive and validate a causal diffusion-inspired filter,",
            "then compare its transfer function and delay with causal baselines",
            "before considering any real-data study.",
            "",
        )
    )
    return "\n".join(lines)
