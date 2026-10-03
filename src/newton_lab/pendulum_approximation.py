"""Phase 25 study of the undamped pendulum's small-angle approximation.

The study reuses the production full-sine pendulum solver and keeps exact
period predictions, numerical trajectories, and post hoc estimates separate.
"""

from __future__ import annotations

from math import pi, sqrt

import numpy as np
from pydantic import Field
from scipy.special import ellipk  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.pendulum import DampedPendulumResult, simulate_damped_pendulum
from newton_lab.simulation import SimulationModel

AMPLITUDES_RAD = (0.01, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
BASE_RELATIVE_TOLERANCE = 1e-10
BASE_ABSOLUTE_TOLERANCE = 1e-12
BASE_SAMPLE_COUNT = 12001
OBSERVATION_SMALL_ANGLE_PERIODS = 12
PERIOD_RELATIVE_TOLERANCE = 0.01
DISPLACEMENT_TOLERANCE_RAD = 0.05
DISPLACEMENT_WINDOW_SMALL_ANGLE_PERIODS = 3
PHASE_TOLERANCE_RAD = 0.1
PHASE_WINDOW_SMALL_ANGLE_PERIODS = 5


class PendulumWindowError(SimulationModel):
    """Trajectory and phase comparison over one requested horizon."""

    window_small_angle_periods: int = Field(ge=1)
    window_duration_s: float
    maximum_absolute_angle_error_rad: float
    accumulated_phase_difference_rad: float


class PendulumAmplitudeComparison(SimulationModel):
    """Analytical predictions and numerical estimates for one release angle."""

    initial_angle_rad: float
    initial_angle_deg: float
    small_angle_period_s: float
    nonlinear_analytical_period_s: float
    nonlinear_numerical_period_s: float
    numerical_period_absolute_error_s: float
    small_angle_period_difference_s: float
    relative_period_error_denominator_nonlinear: float
    relative_period_error_denominator_small_angle: float
    windows: tuple[PendulumWindowError, ...]
    solver_method: str
    relative_tolerance: float
    absolute_tolerance: float
    sample_count: int
    duration_s: float
    function_evaluations: int


class PendulumNumericalSensitivity(SimulationModel):
    """Focused tolerance and output-grid sensitivity for one amplitude."""

    initial_angle_rad: float
    comparison_name: str
    base_sample_count: int
    comparison_sample_count: int
    maximum_shared_sample_angle_difference_rad: float
    base_numerical_period_s: float
    comparison_numerical_period_s: float
    numerical_period_difference_s: float


class PendulumApproximationResearch(SimulationModel):
    """Reproducible amplitude sweep and numerical reliability checks."""

    mass_kg: float
    length_m: float
    gravity_m_per_s2: float
    damping_coefficient_kg_m2_per_s: float
    initial_angular_velocity_rad_per_s: float
    small_angle_angular_frequency_rad_per_s: float
    small_angle_period_s: float
    amplitude_comparisons: tuple[PendulumAmplitudeComparison, ...]
    numerical_sensitivities: tuple[PendulumNumericalSensitivity, ...]
    assumptions: tuple[str, ...] = Field(min_length=1)
    limitations: tuple[str, ...] = Field(min_length=1)


def nonlinear_period_reference(
    initial_angle_rad: float, *, length_m: float, gravity_m_per_s2: float
) -> float:
    """Return the ideal undamped period using SciPy's elliptic parameter m.

    SciPy ``ellipk(m)`` accepts the parameter ``m`` in
    ``integral (1 - m sin(t)^2)^(-1/2) dt``. Here
    ``m = sin(initial_angle/2)^2``; it is not the elliptic modulus.
    """
    if not np.isfinite((initial_angle_rad, length_m, gravity_m_per_s2)).all():
        raise ScientificValidationError("period inputs must be finite")
    if not 0.0 < initial_angle_rad < pi:
        raise ScientificValidationError("initial angle must be between 0 and pi rad")
    if length_m <= 0.0 or gravity_m_per_s2 <= 0.0:
        raise ScientificValidationError("length and gravity must be positive")
    parameter = float(np.sin(initial_angle_rad / 2.0) ** 2)
    return 4.0 * sqrt(length_m / gravity_m_per_s2) * float(ellipk(parameter))


def _estimate_period_from_downward_zero_crossings(
    trajectory: DampedPendulumResult,
) -> float:
    """Estimate period from linearly interpolated positive-to-negative zeros.

    Linear interpolation limits grid bias but does not remove it; this remains
    an estimate from the solver's finite output samples.
    """
    time = trajectory.time_s
    angle = trajectory.angle_rad
    crossing_times: list[float] = []
    for index in range(len(angle) - 1):
        left, right = float(angle[index]), float(angle[index + 1])
        if left > 0.0 and right <= 0.0:
            fraction = left / (left - right)
            crossing_times.append(
                float(time[index] + fraction * (time[index + 1] - time[index]))
            )
    if len(crossing_times) < 2:
        raise ScientificValidationError(
            "period estimation requires at least two downward zero crossings"
        )
    return float(np.mean(np.diff(crossing_times)))


def _phase_difference(
    trajectory: DampedPendulumResult, *, omega0: float, initial_angle: float
) -> np.ndarray:
    """Return unwrapped linear-coordinate phase lag in radians.

    State phase is ``atan2(-theta_dot/(omega0*theta0), theta/theta0)``.
    The linear model has phase ``omega0*t``. Their unwrapped difference is
    reported, avoiding a wrapped phase that would hide accumulated cycle lag.
    """
    raw_phase = np.arctan2(
        -trajectory.angular_velocity_rad_per_s / (omega0 * initial_angle),
        trajectory.angle_rad / initial_angle,
    )
    nonlinear_phase = np.unwrap(raw_phase)
    return omega0 * trajectory.time_s - nonlinear_phase


def _window_comparisons(
    trajectory: DampedPendulumResult,
    *,
    initial_angle: float,
    omega0: float,
    small_angle_period: float,
) -> tuple[PendulumWindowError, ...]:
    linear_angle = initial_angle * np.cos(omega0 * trajectory.time_s)
    phase_difference = _phase_difference(
        trajectory, omega0=omega0, initial_angle=initial_angle
    )
    records: list[PendulumWindowError] = []
    for cycles in (1, 3, 5, 10):
        duration = cycles * small_angle_period
        mask = trajectory.time_s <= duration + 1e-12
        if not np.any(mask):
            continue
        records.append(
            PendulumWindowError(
                window_small_angle_periods=cycles,
                window_duration_s=duration,
                maximum_absolute_angle_error_rad=float(
                    np.max(np.abs(trajectory.angle_rad[mask] - linear_angle[mask]))
                ),
                accumulated_phase_difference_rad=float(
                    np.interp(duration, trajectory.time_s, phase_difference)
                ),
            )
        )
    return tuple(records)


def _compare_amplitude(
    initial_angle: float,
    *,
    mass: float,
    length: float,
    gravity: float,
    omega0: float,
    small_period: float,
    relative_tolerance: float,
    absolute_tolerance: float,
    sample_count: int,
) -> tuple[PendulumAmplitudeComparison, DampedPendulumResult]:
    duration = OBSERVATION_SMALL_ANGLE_PERIODS * small_period
    trajectory = simulate_damped_pendulum(
        mass,
        length,
        0.0,
        gravity,
        initial_angle,
        0.0,
        duration,
        num_points=sample_count,
        method="DOP853",
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )
    numerical_period = _estimate_period_from_downward_zero_crossings(trajectory)
    analytical_period = nonlinear_period_reference(
        initial_angle, length_m=length, gravity_m_per_s2=gravity
    )
    windows = _window_comparisons(
        trajectory,
        initial_angle=initial_angle,
        omega0=omega0,
        small_angle_period=small_period,
    )
    return (
        PendulumAmplitudeComparison(
            initial_angle_rad=initial_angle,
            initial_angle_deg=float(np.degrees(initial_angle)),
            small_angle_period_s=small_period,
            nonlinear_analytical_period_s=analytical_period,
            nonlinear_numerical_period_s=numerical_period,
            numerical_period_absolute_error_s=abs(numerical_period - analytical_period),
            small_angle_period_difference_s=analytical_period - small_period,
            relative_period_error_denominator_nonlinear=(
                abs(analytical_period - small_period) / analytical_period
            ),
            relative_period_error_denominator_small_angle=(
                abs(analytical_period - small_period) / small_period
            ),
            windows=windows,
            solver_method=trajectory.method,
            relative_tolerance=trajectory.relative_tolerance,
            absolute_tolerance=trajectory.absolute_tolerance,
            sample_count=sample_count,
            duration_s=duration,
            function_evaluations=trajectory.function_evaluations,
        ),
        trajectory,
    )


def _numerical_sensitivity(
    initial_angle: float,
    *,
    mass: float,
    length: float,
    gravity: float,
    small_period: float,
    name: str,
    relative_tolerance: float,
    absolute_tolerance: float,
    sample_count: int,
) -> PendulumNumericalSensitivity:
    base, base_trajectory = _compare_amplitude(
        initial_angle,
        mass=mass,
        length=length,
        gravity=gravity,
        omega0=sqrt(gravity / length),
        small_period=small_period,
        relative_tolerance=BASE_RELATIVE_TOLERANCE,
        absolute_tolerance=BASE_ABSOLUTE_TOLERANCE,
        sample_count=BASE_SAMPLE_COUNT,
    )
    other, other_trajectory = _compare_amplitude(
        initial_angle,
        mass=mass,
        length=length,
        gravity=gravity,
        omega0=sqrt(gravity / length),
        small_period=small_period,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        sample_count=sample_count,
    )
    if len(base_trajectory.time_s) == len(other_trajectory.time_s):
        difference = float(
            np.max(np.abs(base_trajectory.angle_rad - other_trajectory.angle_rad))
        )
    elif (len(base_trajectory.time_s) - 1) % (len(other_trajectory.time_s) - 1) == 0:
        stride = (len(base_trajectory.time_s) - 1) // (len(other_trajectory.time_s) - 1)
        difference = float(
            np.max(
                np.abs(base_trajectory.angle_rad[::stride] - other_trajectory.angle_rad)
            )
        )
    else:
        difference = float(
            np.max(
                np.abs(
                    base_trajectory.angle_rad
                    - np.interp(
                        base_trajectory.time_s,
                        other_trajectory.time_s,
                        other_trajectory.angle_rad,
                    )
                )
            )
        )
    return PendulumNumericalSensitivity(
        initial_angle_rad=initial_angle,
        comparison_name=name,
        base_sample_count=BASE_SAMPLE_COUNT,
        comparison_sample_count=sample_count,
        maximum_shared_sample_angle_difference_rad=difference,
        base_numerical_period_s=base.nonlinear_numerical_period_s,
        comparison_numerical_period_s=other.nonlinear_numerical_period_s,
        numerical_period_difference_s=(
            other.nonlinear_numerical_period_s - base.nonlinear_numerical_period_s
        ),
    )


def run_pendulum_approximation_research() -> PendulumApproximationResearch:
    """Run the fixed Phase 25 amplitude, window, and numerical checks."""
    mass, length, gravity = 1.0, 1.0, 9.81
    omega0 = sqrt(gravity / length)
    small_period = 2.0 * pi / omega0
    comparisons = tuple(
        _compare_amplitude(
            angle,
            mass=mass,
            length=length,
            gravity=gravity,
            omega0=omega0,
            small_period=small_period,
            relative_tolerance=BASE_RELATIVE_TOLERANCE,
            absolute_tolerance=BASE_ABSOLUTE_TOLERANCE,
            sample_count=BASE_SAMPLE_COUNT,
        )[0]
        for angle in AMPLITUDES_RAD
    )
    sensitivities = tuple(
        _numerical_sensitivity(
            angle,
            mass=mass,
            length=length,
            gravity=gravity,
            small_period=small_period,
            name=name,
            relative_tolerance=rtol,
            absolute_tolerance=atol,
            sample_count=count,
        )
        for angle, name, rtol, atol, count in (
            (0.1, "looser_tolerance_same_sampling", 1e-8, 1e-10, BASE_SAMPLE_COUNT),
            (
                1.0,
                "denser_sampling_same_tolerance",
                BASE_RELATIVE_TOLERANCE,
                BASE_ABSOLUTE_TOLERANCE,
                24001,
            ),
            (
                2.0,
                "coarser_sampling_same_tolerance",
                BASE_RELATIVE_TOLERANCE,
                BASE_ABSOLUTE_TOLERANCE,
                601,
            ),
        )
    )
    return PendulumApproximationResearch(
        mass_kg=mass,
        length_m=length,
        gravity_m_per_s2=gravity,
        damping_coefficient_kg_m2_per_s=0.0,
        initial_angular_velocity_rad_per_s=0.0,
        small_angle_angular_frequency_rad_per_s=omega0,
        small_angle_period_s=small_period,
        amplitude_comparisons=comparisons,
        numerical_sensitivities=sensitivities,
        assumptions=(
            "Point-mass bob, massless rigid rod, fixed length and pivot, "
            "uniform gravity.",
            "Undamped, unforced full-sine model; release from rest at positive "
            "amplitude.",
            "Small-angle comparison is the exact solution of the linearized equation.",
        ),
        limitations=(
            "All evidence is numerical or analytical within the idealized "
            "equations; no real measurements are used.",
            "Period estimates use finite uniform samples and linearly "
            "interpolated downward zero crossings.",
            "Tolerance-based amplitude limits are grid-limited and depend on "
            "the stated metric and horizon.",
        ),
    )


def render_pendulum_approximation_report(
    results: PendulumApproximationResearch,
) -> str:
    """Render reproducible result tables for the Phase 25 research report."""
    lines = [
        "## Computed amplitude-sweep results",
        "",
        "| Initial angle (rad) | Initial angle (deg) | Solver nfev | "
        "T small-angle (s) | "
        "T nonlinear analytical (s) | T nonlinear sampled estimate (s) | "
        "Estimate abs error (s) | Delta T (s) | Relative error / T nonlinear | "
        "Relative error / T small-angle | Max angle error at 1, 3, 5, 10 T0 "
        "(rad) | Phase lag at 1, 3, 5, 10 T0 (rad) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for item in results.amplitude_comparisons:
        displacement = ", ".join(
            f"{window.maximum_absolute_angle_error_rad:.6g}" for window in item.windows
        )
        phase = ", ".join(
            f"{window.accumulated_phase_difference_rad:.6g}" for window in item.windows
        )
        lines.append(
            f"| {item.initial_angle_rad:g} | {item.initial_angle_deg:.4f} | "
            f"{item.function_evaluations} | "
            f"{item.small_angle_period_s:.9g} | "
            f"{item.nonlinear_analytical_period_s:.9g} | "
            f"{item.nonlinear_numerical_period_s:.9g} | "
            f"{item.numerical_period_absolute_error_s:.3g} | "
            f"{item.small_angle_period_difference_s:.6g} | "
            f"{item.relative_period_error_denominator_nonlinear:.6g} | "
            f"{item.relative_period_error_denominator_small_angle:.6g} | "
            f"{displacement} | {phase} |"
        )
    lines.extend(
        (
            "",
            "### Numerical setting comparisons",
            "",
            "| Amplitude (rad) | Comparison | Samples (base/comparison) | "
            "Maximum shared-sample angle difference (rad) | "
            "Period estimate base (s) | Period estimate comparison (s) | "
            "Period difference (s) |",
            "|---:|---|---:|---:|---:|---:|---:|",
        )
    )
    for sensitivity in results.numerical_sensitivities:
        lines.append(
            f"| {sensitivity.initial_angle_rad:g} | {sensitivity.comparison_name} | "
            f"{sensitivity.base_sample_count}/{sensitivity.comparison_sample_count} | "
            f"{sensitivity.maximum_shared_sample_angle_difference_rad:.6g} | "
            f"{sensitivity.base_numerical_period_s:.9g} | "
            f"{sensitivity.comparison_numerical_period_s:.9g} | "
            f"{sensitivity.numerical_period_difference_s:.3g} |"
        )
    lines.extend(("", "### Tolerance-based grid limits", ""))

    def period_passes(row: PendulumAmplitudeComparison) -> bool:
        return (
            row.relative_period_error_denominator_nonlinear <= PERIOD_RELATIVE_TOLERANCE
        )

    def displacement_passes(row: PendulumAmplitudeComparison) -> bool:
        error = next(
            window.maximum_absolute_angle_error_rad
            for window in row.windows
            if window.window_small_angle_periods == 3
        )
        return error <= DISPLACEMENT_TOLERANCE_RAD

    def phase_passes(row: PendulumAmplitudeComparison) -> bool:
        lag = next(
            window.accumulated_phase_difference_rad
            for window in row.windows
            if window.window_small_angle_periods == 5
        )
        return lag <= PHASE_TOLERANCE_RAD

    thresholds = (
        (
            "Relative period error <= 1% (denominator: nonlinear analytical period)",
            period_passes,
        ),
        (
            "Maximum angle error <= 0.05 rad over 3 T0",
            displacement_passes,
        ),
        (
            "Accumulated phase lag <= 0.1 rad at 5 T0",
            phase_passes,
        ),
    )
    for label, predicate in thresholds:
        passing = [
            row.initial_angle_rad
            for row in results.amplitude_comparisons
            if predicate(row)
        ]
        first_failed = next(
            (
                row.initial_angle_rad
                for row in results.amplitude_comparisons
                if not predicate(row)
            ),
            None,
        )
        lines.append(
            f"- {label}: largest tested passing amplitude `"
            f"{max(passing):g} rad` ({np.degrees(max(passing)):.4f} deg)"
            + (
                f"; next tested failure at `{first_failed:g} rad`."
                if first_failed is not None
                else "; all tested amplitudes pass."
            )
        )
    lines.extend(("", "### Numerical settings", ""))
    for item in results.amplitude_comparisons[:1]:
        lines.append(
            f"- Main sweep uses `{item.solver_method}`, "
            f"`rtol={item.relative_tolerance:g}`, `atol={item.absolute_tolerance:g}`, "
            f"{item.sample_count} uniform requested output points over "
            f"{item.duration_s:.6g} s "
            f"({OBSERVATION_SMALL_ANGLE_PERIODS} small-angle periods)."
        )
        lines.append(
            f"- Solver function evaluations in the first amplitude case: "
            f"{item.function_evaluations}. "
            "Adaptive internal steps differ from the uniform output grid."
        )
    return "\n".join(lines) + "\n"
