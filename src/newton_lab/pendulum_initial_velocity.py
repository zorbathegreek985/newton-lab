"""Phase 26 pendulum study with nonzero initial angular velocity.

This module reuses Newton Lab's full-sine pendulum solver. Energy determines
the libration turning amplitude; only librating states use the elliptic period
reference. Rotations and the separatrix are explicitly classified instead.
"""

from __future__ import annotations

from enum import StrEnum
from math import pi, sqrt

import numpy as np
from pydantic import Field, model_validator
from scipy.special import ellipk  # type: ignore[import-untyped]

from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.pendulum import DampedPendulumResult, simulate_damped_pendulum
from newton_lab.simulation import SimulationModel

ROUND_OFF_ENERGY_TOLERANCE = 64.0 * np.finfo(np.float64).eps
NEAR_SEPARATRIX_ENERGY_GAP = 0.02
BASE_RELATIVE_TOLERANCE = 1e-10
BASE_ABSOLUTE_TOLERANCE = 1e-12
BASE_SAMPLE_COUNT = 24001
PERIOD_RELATIVE_TOLERANCE = 0.01
ANGLE_ERROR_TOLERANCE_RAD = 0.05
ANGLE_ERROR_WINDOW_LINEAR_PERIODS = 3
PHASE_ERROR_TOLERANCE_RAD = 0.1
PHASE_ERROR_WINDOW_LINEAR_PERIODS = 5


class InitialStateClass(StrEnum):
    """Energy-based state classification under the principal-angle convention."""

    LIBRATING = "librating"
    NEAR_SEPARATRIX = "near_separatrix_libration"
    SEPARATRIX = "separatrix"
    ROTATING = "rotating"
    EQUILIBRIUM = "stable_equilibrium"
    INVALID = "invalid_initial_condition"


class CaseStatus(StrEnum):
    """Outcome of classifying, simulating, and measuring one sweep case."""

    ANALYZED = "analyzed"
    SEPARATRIX_EXCLUDED = "separatrix_period_not_finite"
    ROTATION_EXCLUDED = "rotation_excluded_from_libration_period"
    EQUILIBRIUM_EXCLUDED = "equilibrium_has_no_observed_period"
    INVALID_INPUT = "invalid_input"
    SOLVER_FAILURE = "solver_failure"
    INSUFFICIENT_DATA = "insufficient_period_data"


class InitialConditionSpec(SimulationModel):
    """One ordered initial-condition pair, with radians and radians/second."""

    case_id: str
    initial_angle_rad: float
    initial_angular_velocity_rad_per_s: float


class InitialStateClassification(SimulationModel):
    """Dimensionless energy and turning amplitude for an initial state."""

    classification: InitialStateClass
    normalized_energy: float
    distance_below_separatrix: float
    turning_point_amplitude_rad: float | None
    reason: str


class AccuracyCriteria(SimulationModel):
    """Configurable, observable-specific Phase 26 screening tolerances."""

    relative_period_error: float = Field(default=0.01, gt=0.0)
    maximum_angle_error_rad: float = Field(default=0.05, gt=0.0)
    angle_window_linear_periods: int = Field(default=3, ge=1)
    maximum_absolute_phase_difference_rad: float = Field(default=0.1, gt=0.0)
    phase_window_linear_periods: int = Field(default=5, ge=1)

    @model_validator(mode="after")
    def validate_supported_windows(self) -> AccuracyCriteria:
        if self.angle_window_linear_periods not in (1, 3, 5, 10):
            raise ValueError(
                "angle window must be one of 1, 3, 5, or 10 linear periods"
            )
        if self.phase_window_linear_periods not in (1, 3, 5, 10):
            raise ValueError(
                "phase window must be one of 1, 3, 5, or 10 linear periods"
            )
        return self


class PeriodEstimate(SimulationModel):
    """Same-direction zero-crossing period estimate from finite samples."""

    period_s: float | None
    crossing_direction: str | None
    crossing_count: int
    status: str
    explanation: str


class WindowMetrics(SimulationModel):
    """Sampled angle and unwrapped phase discrepancies at one horizon."""

    window_linear_periods: int = Field(ge=1)
    duration_s: float
    maximum_sampled_absolute_angle_error_rad: float
    signed_accumulated_phase_difference_rad: float


class InitialConditionResult(SimulationModel):
    """Traceable outcome for one input pair, including explicit exclusions."""

    case_id: str
    sweep_index: int = Field(ge=0)
    initial_angle_rad: float
    initial_angular_velocity_rad_per_s: float
    initial_linear_phase_rad: float | None
    classification: InitialStateClass
    status: CaseStatus
    normalized_energy: float | None
    distance_below_separatrix: float | None
    turning_point_amplitude_rad: float | None
    amplitude_increase_over_initial_angle_rad: float | None
    linear_period_s: float
    nonlinear_analytical_period_s: float | None
    nonlinear_numerical_period_s: float | None
    numerical_period_absolute_error_s: float | None
    signed_period_difference_s: float | None
    signed_relative_period_difference_nonlinear_denominator: float | None
    period_estimate: PeriodEstimate | None
    windows: tuple[WindowMetrics, ...] = ()
    period_criterion_passed: bool | None
    angle_criterion_passed: bool | None
    phase_criterion_passed: bool | None
    solver_method: str
    relative_tolerance: float
    absolute_tolerance: float
    sample_count: int
    duration_s: float | None
    function_evaluations: int | None
    failure_reason: str | None


class ThresholdSlice(SimulationModel):
    """Grid-limited approximation limit conditional on signed initial speed."""

    observable: str
    initial_angular_velocity_rad_per_s: float
    largest_tested_passing_absolute_initial_angle_rad: float | None
    next_higher_tested_absolute_initial_angle_rad: float | None
    next_higher_point_failed: bool | None
    passing_case_ids: tuple[str, ...]
    failing_case_ids: tuple[str, ...]
    interpretation: str


class NumericalSensitivity(SimulationModel):
    """Solver-tolerance or sampling sensitivity for one initial state."""

    case_id: str
    comparison: str
    base_sample_count: int
    comparison_sample_count: int
    maximum_angle_difference_at_shared_samples_rad: float | None
    base_period_estimate_s: float | None
    comparison_period_estimate_s: float | None
    period_estimate_difference_s: float | None
    base_solver_success: bool
    comparison_solver_success: bool


class PendulumInitialVelocityResearch(SimulationModel):
    """Complete Phase 26 ordered sweep, criteria, and sensitivity output."""

    mass_kg: float
    length_m: float
    gravity_m_per_s2: float
    damping_coefficient_kg_m2_per_s: float
    small_angle_natural_frequency_rad_per_s: float
    linear_period_s: float
    criteria: AccuracyCriteria
    sweep_order: tuple[InitialConditionSpec, ...]
    cases: tuple[InitialConditionResult, ...]
    threshold_slices: tuple[ThresholdSlice, ...]
    numerical_sensitivities: tuple[NumericalSensitivity, ...]
    assumptions: tuple[str, ...] = Field(min_length=1)
    limitations: tuple[str, ...] = Field(min_length=1)


def _validate_physical_inputs(
    initial_angle_rad: float,
    initial_angular_velocity_rad_per_s: float,
    length_m: float,
    gravity_m_per_s2: float,
) -> None:
    values = (
        initial_angle_rad,
        initial_angular_velocity_rad_per_s,
        length_m,
        gravity_m_per_s2,
    )
    if not np.isfinite(values).all():
        raise ScientificValidationError(
            "initial conditions and parameters must be finite"
        )
    if length_m <= 0.0 or gravity_m_per_s2 <= 0.0:
        raise ScientificValidationError("length and gravity must be positive")
    if not -pi <= initial_angle_rad <= pi:
        raise ScientificValidationError(
            "initial angle must use the principal convention [-pi, pi] rad"
        )


def classify_initial_condition(
    initial_angle_rad: float,
    initial_angular_velocity_rad_per_s: float,
    *,
    length_m: float,
    gravity_m_per_s2: float,
) -> InitialStateClassification:
    """Classify a state using dimensionless energy relative to the separatrix.

    The initial angle is in the principal range ``[-pi, pi]``. Define
    ``h = 1-cos(theta) + L*omega**2/(2g)``. Libration requires ``h < 2``;
    ``h=2`` is the separatrix and ``h>2`` is rotation. The near-separatrix
    label applies to valid librations with ``2-h <= 0.02``. The energy is
    evaluated with ``2*sin(theta/2)**2`` to avoid cancellation near zero.
    Only a roundoff-sized excess above 2 is treated as separatrix-level; no
    physically meaningful excess is clamped into a libration.
    """
    _validate_physical_inputs(
        initial_angle_rad,
        initial_angular_velocity_rad_per_s,
        length_m,
        gravity_m_per_s2,
    )
    normalized_energy = float(
        2.0 * np.sin(initial_angle_rad / 2.0) ** 2
        + length_m * initial_angular_velocity_rad_per_s**2 / (2.0 * gravity_m_per_s2)
    )
    gap = 2.0 - normalized_energy
    if normalized_energy == 0.0:
        return InitialStateClassification(
            classification=InitialStateClass.EQUILIBRIUM,
            normalized_energy=normalized_energy,
            distance_below_separatrix=2.0,
            turning_point_amplitude_rad=0.0,
            reason="zero angle and zero angular velocity are stable equilibrium",
        )
    if abs(gap) <= ROUND_OFF_ENERGY_TOLERANCE:
        return InitialStateClassification(
            classification=InitialStateClass.SEPARATRIX,
            normalized_energy=normalized_energy,
            distance_below_separatrix=gap,
            turning_point_amplitude_rad=None,
            reason="energy is at the separatrix within floating-point roundoff",
        )
    if gap < 0.0:
        return InitialStateClassification(
            classification=InitialStateClass.ROTATING,
            normalized_energy=normalized_energy,
            distance_below_separatrix=gap,
            turning_point_amplitude_rad=None,
            reason="normalized energy exceeds the separatrix value 2",
        )
    amplitude = float(2.0 * np.arcsin(np.sqrt(normalized_energy / 2.0)))
    near = gap <= NEAR_SEPARATRIX_ENERGY_GAP
    return InitialStateClassification(
        classification=(
            InitialStateClass.NEAR_SEPARATRIX if near else InitialStateClass.LIBRATING
        ),
        normalized_energy=normalized_energy,
        distance_below_separatrix=gap,
        turning_point_amplitude_rad=amplitude,
        reason=(
            "libration is close to the separatrix under the 0.02 energy-gap flag"
            if near
            else "energy is below the separatrix and the state librates"
        ),
    )


def turning_point_amplitude(
    initial_angle_rad: float,
    initial_angular_velocity_rad_per_s: float,
    *,
    length_m: float,
    gravity_m_per_s2: float,
) -> float:
    """Return positive libration turning amplitude in radians.

    The result is ``2*asin(sqrt(h/2))`` from energy conservation. An
    equilibrium returns zero. Separatrix and rotating states are rejected.
    """
    state = classify_initial_condition(
        initial_angle_rad,
        initial_angular_velocity_rad_per_s,
        length_m=length_m,
        gravity_m_per_s2=gravity_m_per_s2,
    )
    if state.classification in (
        InitialStateClass.SEPARATRIX,
        InitialStateClass.ROTATING,
    ):
        raise ScientificValidationError(
            f"turning-point amplitude is not a libration amplitude: {state.reason}"
        )
    assert state.turning_point_amplitude_rad is not None
    return state.turning_point_amplitude_rad


def nonlinear_period_reference(
    initial_angle_rad: float,
    initial_angular_velocity_rad_per_s: float,
    *,
    length_m: float,
    gravity_m_per_s2: float,
) -> float:
    """Return exact libration period from energy-derived turning amplitude.

    SciPy's ``ellipk(m)`` accepts the parameter ``m`` rather than modulus
    ``k``. Here ``m=sin(theta_max/2)**2``. Equilibrium, separatrix, and
    rotation have no nontrivial finite libration period and are rejected.
    """
    state = classify_initial_condition(
        initial_angle_rad,
        initial_angular_velocity_rad_per_s,
        length_m=length_m,
        gravity_m_per_s2=gravity_m_per_s2,
    )
    if state.classification in (
        InitialStateClass.SEPARATRIX,
        InitialStateClass.ROTATING,
        InitialStateClass.EQUILIBRIUM,
    ):
        raise ScientificValidationError(
            f"finite nonlinear libration period is unavailable: {state.reason}"
        )
    amplitude = state.turning_point_amplitude_rad
    assert amplitude is not None
    parameter = float(np.sin(amplitude / 2.0) ** 2)
    return 4.0 * sqrt(length_m / gravity_m_per_s2) * float(ellipk(parameter))


def linear_period_reference(*, length_m: float, gravity_m_per_s2: float) -> float:
    """Return the small-angle period ``2*pi*sqrt(L/g)`` in seconds."""
    if not np.isfinite((length_m, gravity_m_per_s2)).all():
        raise ScientificValidationError("length and gravity must be finite")
    if length_m <= 0.0 or gravity_m_per_s2 <= 0.0:
        raise ScientificValidationError("length and gravity must be positive")
    return 2.0 * pi * sqrt(length_m / gravity_m_per_s2)


def estimate_nonlinear_period(
    trajectory: DampedPendulumResult,
) -> PeriodEstimate:
    """Estimate period from repeated crossings of a common direction.

    Both positive-going and negative-going crossings are collected, including
    a time-zero crossing when the initial state is at equilibrium with nonzero
    velocity. The direction with more complete crossings is used. Successive
    crossings of that direction are one period apart; fewer than two yield an
    explicit insufficient-data result rather than a guessed period.
    """
    time = trajectory.time_s
    angle = trajectory.angle_rad
    if len(time) != len(angle) or len(time) < 2:
        return PeriodEstimate(
            period_s=None,
            crossing_direction=None,
            crossing_count=0,
            status="insufficient_trajectory_data",
            explanation="at least two aligned angle and time samples are required",
        )
    positive_crossings: list[float] = []
    negative_crossings: list[float] = []
    initial_angle = float(angle[0])
    initial_velocity = trajectory.initial_angular_velocity_rad_per_s
    if initial_angle == 0.0 and initial_velocity != 0.0:
        (positive_crossings if initial_velocity > 0 else negative_crossings).append(0.0)
    for index in range(len(angle) - 1):
        left, right = float(angle[index]), float(angle[index + 1])
        if left < 0.0 <= right:
            fraction = -left / (right - left)
            positive_crossings.append(
                float(time[index] + fraction * (time[index + 1] - time[index]))
            )
        elif left > 0.0 >= right:
            fraction = left / (left - right)
            negative_crossings.append(
                float(time[index] + fraction * (time[index + 1] - time[index]))
            )
    direction, crossings = max(
        (("positive", positive_crossings), ("negative", negative_crossings)),
        key=lambda item: len(item[1]),
    )
    if len(crossings) < 2:
        return PeriodEstimate(
            period_s=None,
            crossing_direction=direction if crossings else None,
            crossing_count=len(crossings),
            status="insufficient_complete_cycles",
            explanation=(
                "at least two crossings in the same direction are required; "
                "the record may end before a complete period is observed"
            ),
        )
    intervals = np.diff(np.asarray(crossings, dtype=np.float64))
    return PeriodEstimate(
        period_s=float(np.mean(intervals)),
        crossing_direction=direction,
        crossing_count=len(crossings),
        status="estimated_from_same_direction_crossings",
        explanation=(
            "mean interval between same-direction zero crossings with linear "
            "interpolation of finite output samples"
        ),
    )


def _trajectory_metrics(
    trajectory: DampedPendulumResult,
    *,
    initial_angle_rad: float,
    initial_angular_velocity_rad_per_s: float,
    omega_linear: float,
    linear_period_s: float,
    maximum_linear_periods: int,
) -> tuple[WindowMetrics, ...]:
    times = trajectory.time_s
    linear_angle = initial_angle_rad * np.cos(omega_linear * times) + (
        initial_angular_velocity_rad_per_s / omega_linear
    ) * np.sin(omega_linear * times)
    linear_velocity = -initial_angle_rad * omega_linear * np.sin(
        omega_linear * times
    ) + initial_angular_velocity_rad_per_s * np.cos(omega_linear * times)
    nonlinear_phase = np.unwrap(
        np.arctan2(
            -trajectory.angular_velocity_rad_per_s / omega_linear,
            trajectory.angle_rad,
        )
    )
    linear_phase = np.unwrap(np.arctan2(-linear_velocity / omega_linear, linear_angle))
    # Match phase branches at t=0, including the initial direction of motion.
    linear_phase += (
        2.0
        * pi
        * round((float(nonlinear_phase[0]) - float(linear_phase[0])) / (2.0 * pi))
    )
    phase_difference = linear_phase - nonlinear_phase
    windows: list[WindowMetrics] = []
    for periods in (1, 3, 5, 10):
        if periods > maximum_linear_periods:
            continue
        duration = periods * linear_period_s
        mask = times <= duration + 1e-12
        if not np.any(mask):
            continue
        windows.append(
            WindowMetrics(
                window_linear_periods=periods,
                duration_s=duration,
                maximum_sampled_absolute_angle_error_rad=float(
                    np.max(np.abs(trajectory.angle_rad[mask] - linear_angle[mask]))
                ),
                signed_accumulated_phase_difference_rad=float(
                    np.interp(duration, times, phase_difference)
                ),
            )
        )
    return tuple(windows)


def _initial_linear_phase(
    initial_angle_rad: float,
    initial_angular_velocity_rad_per_s: float,
    omega_linear: float,
) -> float:
    """Return linear-coordinate phase at release, including direction."""
    return float(
        np.arctan2(
            -initial_angular_velocity_rad_per_s / omega_linear, initial_angle_rad
        )
    )


def _failed_case(
    spec: InitialConditionSpec,
    index: int,
    *,
    linear_period_s: float,
    status: CaseStatus,
    classification: InitialStateClass,
    reason: str,
    method: str,
    relative_tolerance: float,
    absolute_tolerance: float,
    sample_count: int,
    normalized_energy: float | None = None,
    energy_gap: float | None = None,
    turning_amplitude: float | None = None,
    analytical_period: float | None = None,
    period_estimate: PeriodEstimate | None = None,
    duration_s: float | None = None,
    function_evaluations: int | None = None,
) -> InitialConditionResult:
    return InitialConditionResult(
        case_id=spec.case_id,
        sweep_index=index,
        initial_angle_rad=spec.initial_angle_rad,
        initial_angular_velocity_rad_per_s=spec.initial_angular_velocity_rad_per_s,
        initial_linear_phase_rad=(
            None
            if classification is InitialStateClass.INVALID
            else _initial_linear_phase(
                spec.initial_angle_rad,
                spec.initial_angular_velocity_rad_per_s,
                2.0 * pi / linear_period_s,
            )
        ),
        classification=classification,
        status=status,
        normalized_energy=normalized_energy,
        distance_below_separatrix=energy_gap,
        turning_point_amplitude_rad=turning_amplitude,
        amplitude_increase_over_initial_angle_rad=(
            None
            if turning_amplitude is None
            else turning_amplitude - abs(spec.initial_angle_rad)
        ),
        linear_period_s=linear_period_s,
        nonlinear_analytical_period_s=analytical_period,
        nonlinear_numerical_period_s=(
            None if period_estimate is None else period_estimate.period_s
        ),
        numerical_period_absolute_error_s=None,
        signed_period_difference_s=(
            None if analytical_period is None else analytical_period - linear_period_s
        ),
        signed_relative_period_difference_nonlinear_denominator=(
            None
            if analytical_period is None
            else (analytical_period - linear_period_s) / analytical_period
        ),
        period_estimate=period_estimate,
        period_criterion_passed=None,
        angle_criterion_passed=None,
        phase_criterion_passed=None,
        solver_method=method,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        sample_count=sample_count,
        duration_s=duration_s,
        function_evaluations=function_evaluations,
        failure_reason=reason,
    )


def analyze_initial_condition(
    spec: InitialConditionSpec,
    *,
    sweep_index: int,
    mass_kg: float,
    length_m: float,
    gravity_m_per_s2: float,
    criteria: AccuracyCriteria,
    sample_count: int = BASE_SAMPLE_COUNT,
    relative_tolerance: float = BASE_RELATIVE_TOLERANCE,
    absolute_tolerance: float = BASE_ABSOLUTE_TOLERANCE,
) -> InitialConditionResult:
    """Classify and, for periodic states, analyze one initial-condition case."""
    if not np.isfinite(mass_kg) or mass_kg <= 0.0:
        raise ScientificValidationError("mass_kg must be finite and positive")
    linear_period = linear_period_reference(
        length_m=length_m, gravity_m_per_s2=gravity_m_per_s2
    )
    omega_linear = sqrt(gravity_m_per_s2 / length_m)
    try:
        classification = classify_initial_condition(
            spec.initial_angle_rad,
            spec.initial_angular_velocity_rad_per_s,
            length_m=length_m,
            gravity_m_per_s2=gravity_m_per_s2,
        )
    except ScientificValidationError as exc:
        return _failed_case(
            spec,
            sweep_index,
            linear_period_s=linear_period,
            status=CaseStatus.INVALID_INPUT,
            classification=InitialStateClass.INVALID,
            reason=str(exc),
            method="DOP853",
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            sample_count=sample_count,
        )
    state = classification.classification
    if state is InitialStateClass.SEPARATRIX:
        return _failed_case(
            spec,
            sweep_index,
            linear_period_s=linear_period,
            status=CaseStatus.SEPARATRIX_EXCLUDED,
            classification=state,
            reason=classification.reason,
            method="DOP853",
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            sample_count=sample_count,
            normalized_energy=classification.normalized_energy,
            energy_gap=classification.distance_below_separatrix,
        )
    if state is InitialStateClass.ROTATING:
        return _failed_case(
            spec,
            sweep_index,
            linear_period_s=linear_period,
            status=CaseStatus.ROTATION_EXCLUDED,
            classification=state,
            reason=classification.reason,
            method="DOP853",
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            sample_count=sample_count,
            normalized_energy=classification.normalized_energy,
            energy_gap=classification.distance_below_separatrix,
        )
    if state is InitialStateClass.EQUILIBRIUM:
        return _failed_case(
            spec,
            sweep_index,
            linear_period_s=linear_period,
            status=CaseStatus.EQUILIBRIUM_EXCLUDED,
            classification=state,
            reason=classification.reason,
            method="DOP853",
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            sample_count=sample_count,
            normalized_energy=classification.normalized_energy,
            energy_gap=classification.distance_below_separatrix,
            turning_amplitude=0.0,
        )
    amplitude = classification.turning_point_amplitude_rad
    assert amplitude is not None
    analytical_period = nonlinear_period_reference(
        spec.initial_angle_rad,
        spec.initial_angular_velocity_rad_per_s,
        length_m=length_m,
        gravity_m_per_s2=gravity_m_per_s2,
    )
    duration = max(12.0 * linear_period, 4.0 * analytical_period)
    try:
        trajectory = simulate_damped_pendulum(
            mass_kg,
            length_m,
            0.0,
            gravity_m_per_s2,
            spec.initial_angle_rad,
            spec.initial_angular_velocity_rad_per_s,
            duration,
            num_points=sample_count,
            method="DOP853",
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )
    except IntegrationError as exc:
        return _failed_case(
            spec,
            sweep_index,
            linear_period_s=linear_period,
            status=CaseStatus.SOLVER_FAILURE,
            classification=state,
            reason=str(exc),
            method="DOP853",
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
            sample_count=sample_count,
            normalized_energy=classification.normalized_energy,
            energy_gap=classification.distance_below_separatrix,
            turning_amplitude=amplitude,
            analytical_period=analytical_period,
            duration_s=duration,
        )
    period_estimate = estimate_nonlinear_period(trajectory)
    windows = _trajectory_metrics(
        trajectory,
        initial_angle_rad=spec.initial_angle_rad,
        initial_angular_velocity_rad_per_s=spec.initial_angular_velocity_rad_per_s,
        omega_linear=omega_linear,
        linear_period_s=linear_period,
        maximum_linear_periods=10,
    )
    if period_estimate.period_s is None:
        return _failed_case(
            spec,
            sweep_index,
            linear_period_s=linear_period,
            status=CaseStatus.INSUFFICIENT_DATA,
            classification=state,
            reason=period_estimate.explanation,
            method=trajectory.method,
            relative_tolerance=trajectory.relative_tolerance,
            absolute_tolerance=trajectory.absolute_tolerance,
            sample_count=sample_count,
            normalized_energy=classification.normalized_energy,
            energy_gap=classification.distance_below_separatrix,
            turning_amplitude=amplitude,
            analytical_period=analytical_period,
            period_estimate=period_estimate,
            duration_s=duration,
            function_evaluations=trajectory.function_evaluations,
        ).model_copy(update={"windows": windows})
    period_error = period_estimate.period_s - analytical_period
    relative_error = abs(analytical_period - linear_period) / analytical_period
    angle_window = next(
        window.maximum_sampled_absolute_angle_error_rad
        for window in windows
        if window.window_linear_periods == criteria.angle_window_linear_periods
    )
    phase_window = next(
        window.signed_accumulated_phase_difference_rad
        for window in windows
        if window.window_linear_periods == criteria.phase_window_linear_periods
    )
    return InitialConditionResult(
        case_id=spec.case_id,
        sweep_index=sweep_index,
        initial_angle_rad=spec.initial_angle_rad,
        initial_angular_velocity_rad_per_s=spec.initial_angular_velocity_rad_per_s,
        initial_linear_phase_rad=_initial_linear_phase(
            spec.initial_angle_rad,
            spec.initial_angular_velocity_rad_per_s,
            omega_linear,
        ),
        classification=state,
        status=CaseStatus.ANALYZED,
        normalized_energy=classification.normalized_energy,
        distance_below_separatrix=classification.distance_below_separatrix,
        turning_point_amplitude_rad=amplitude,
        amplitude_increase_over_initial_angle_rad=amplitude
        - abs(spec.initial_angle_rad),
        linear_period_s=linear_period,
        nonlinear_analytical_period_s=analytical_period,
        nonlinear_numerical_period_s=period_estimate.period_s,
        numerical_period_absolute_error_s=abs(period_error),
        signed_period_difference_s=analytical_period - linear_period,
        signed_relative_period_difference_nonlinear_denominator=(
            (analytical_period - linear_period) / analytical_period
        ),
        period_estimate=period_estimate,
        windows=windows,
        period_criterion_passed=relative_error <= criteria.relative_period_error,
        angle_criterion_passed=angle_window <= criteria.maximum_angle_error_rad,
        phase_criterion_passed=(
            abs(phase_window) <= criteria.maximum_absolute_phase_difference_rad
        ),
        solver_method=trajectory.method,
        relative_tolerance=trajectory.relative_tolerance,
        absolute_tolerance=trajectory.absolute_tolerance,
        sample_count=sample_count,
        duration_s=duration,
        function_evaluations=trajectory.function_evaluations,
        failure_reason=None,
    )


def _default_sweep(gravity: float, length: float) -> tuple[InitialConditionSpec, ...]:
    """Return the declared ordered grid, including excluded state controls."""
    specs = [
        InitialConditionSpec(
            case_id=f"grid_theta_{angle:g}_omega_{speed:+g}",
            initial_angle_rad=angle,
            initial_angular_velocity_rad_per_s=speed,
        )
        for angle in (0.1, 0.2, 0.35, 0.5, 1.2, 2.0)
        for speed in (-3.0, -1.0, 0.0, 1.0, 3.0)
    ]
    for speed in (4.3, -4.3):
        specs.append(
            InitialConditionSpec(
                case_id=f"high_energy_libration_theta_0_1_omega_{speed:+g}",
                initial_angle_rad=0.1,
                initial_angular_velocity_rad_per_s=speed,
            )
        )
    specs.extend(
        (
            InitialConditionSpec(
                case_id="negative_angle_positive_velocity",
                initial_angle_rad=-0.5,
                initial_angular_velocity_rad_per_s=1.0,
            ),
            InitialConditionSpec(
                case_id="negative_angle_negative_velocity",
                initial_angle_rad=-0.5,
                initial_angular_velocity_rad_per_s=-1.0,
            ),
            InitialConditionSpec(
                case_id="near_boundary_3_0_rest",
                initial_angle_rad=3.0,
                initial_angular_velocity_rad_per_s=0.0,
            ),
            InitialConditionSpec(
                case_id="near_boundary_3_12_rest",
                initial_angle_rad=3.12,
                initial_angular_velocity_rad_per_s=0.0,
            ),
            InitialConditionSpec(
                case_id="separatrix_control",
                initial_angle_rad=0.0,
                initial_angular_velocity_rad_per_s=sqrt(4.0 * gravity / length),
            ),
            InitialConditionSpec(
                case_id="rotational_control",
                initial_angle_rad=0.0,
                initial_angular_velocity_rad_per_s=7.0,
            ),
            InitialConditionSpec(
                case_id="equilibrium_control",
                initial_angle_rad=0.0,
                initial_angular_velocity_rad_per_s=0.0,
            ),
            InitialConditionSpec(
                case_id="invalid_principal_angle_control",
                initial_angle_rad=3.3,
                initial_angular_velocity_rad_per_s=0.0,
            ),
        )
    )
    return tuple(specs)


def _threshold_summaries(
    cases: tuple[InitialConditionResult, ...], criteria: AccuracyCriteria
) -> tuple[ThresholdSlice, ...]:
    definitions = (
        ("relative_period_error", "period_criterion_passed"),
        ("sampled_angle_error", "angle_criterion_passed"),
        ("absolute_accumulated_phase_difference", "phase_criterion_passed"),
    )
    speed_values = tuple(
        dict.fromkeys(
            case.initial_angular_velocity_rad_per_s
            for case in cases
            if case.status is CaseStatus.ANALYZED
        )
    )
    summaries: list[ThresholdSlice] = []
    for observable, field_name in definitions:
        for speed in speed_values:
            group = sorted(
                (
                    case
                    for case in cases
                    if case.initial_angular_velocity_rad_per_s == speed
                    and case.status is CaseStatus.ANALYZED
                ),
                key=lambda case: abs(case.initial_angle_rad),
            )
            passing = [case for case in group if getattr(case, field_name) is True]
            failing = [case for case in group if getattr(case, field_name) is False]
            largest = max(
                (abs(case.initial_angle_rad) for case in passing), default=None
            )
            next_case = next(
                (
                    case
                    for case in group
                    if largest is not None and abs(case.initial_angle_rad) > largest
                ),
                None,
            )
            summaries.append(
                ThresholdSlice(
                    observable=observable,
                    initial_angular_velocity_rad_per_s=speed,
                    largest_tested_passing_absolute_initial_angle_rad=largest,
                    next_higher_tested_absolute_initial_angle_rad=(
                        None if next_case is None else abs(next_case.initial_angle_rad)
                    ),
                    next_higher_point_failed=(
                        None
                        if next_case is None
                        else getattr(next_case, field_name) is False
                    ),
                    passing_case_ids=tuple(case.case_id for case in passing),
                    failing_case_ids=tuple(case.case_id for case in failing),
                    interpretation=(
                        "conditional, grid-limited comparison for this signed initial "
                        "angular velocity; not a universal amplitude threshold"
                    ),
                )
            )
    return tuple(summaries)


def _sensitivity_case(
    spec: InitialConditionSpec,
    *,
    length: float,
    gravity: float,
    sample_count: int,
    rtol: float,
    atol: float,
) -> tuple[DampedPendulumResult | None, PeriodEstimate, float | None]:
    state = classify_initial_condition(
        spec.initial_angle_rad,
        spec.initial_angular_velocity_rad_per_s,
        length_m=length,
        gravity_m_per_s2=gravity,
    )
    if state.turning_point_amplitude_rad is None or state.classification in (
        InitialStateClass.EQUILIBRIUM,
        InitialStateClass.SEPARATRIX,
        InitialStateClass.ROTATING,
    ):
        empty = PeriodEstimate(
            period_s=None,
            crossing_direction=None,
            crossing_count=0,
            status="not_a_libration",
            explanation=state.reason,
        )
        return None, empty, None
    tlinear = linear_period_reference(length_m=length, gravity_m_per_s2=gravity)
    tnonlinear = nonlinear_period_reference(
        spec.initial_angle_rad,
        spec.initial_angular_velocity_rad_per_s,
        length_m=length,
        gravity_m_per_s2=gravity,
    )
    duration = max(12.0 * tlinear, 4.0 * tnonlinear)
    try:
        result = simulate_damped_pendulum(
            1.0,
            length,
            0.0,
            gravity,
            spec.initial_angle_rad,
            spec.initial_angular_velocity_rad_per_s,
            duration,
            num_points=sample_count,
            method="DOP853",
            relative_tolerance=rtol,
            absolute_tolerance=atol,
        )
    except IntegrationError:
        return (
            None,
            PeriodEstimate(
                period_s=None,
                crossing_direction=None,
                crossing_count=0,
                status="solver_failure",
                explanation="integration failed under this sensitivity configuration",
            ),
            None,
        )
    return result, estimate_nonlinear_period(result), duration


def _run_sensitivities(
    specs: tuple[InitialConditionSpec, ...], *, length: float, gravity: float
) -> tuple[NumericalSensitivity, ...]:
    selected = tuple(
        next(spec for spec in specs if spec.case_id == case_id)
        for case_id in (
            "grid_theta_0.1_omega_+1",  # small initial displacement
            "grid_theta_0.5_omega_+3",  # substantially increased turning point
            "near_boundary_3_12_rest",
        )
    )
    records: list[NumericalSensitivity] = []
    for spec in selected:
        base, base_period, _ = _sensitivity_case(
            spec,
            length=length,
            gravity=gravity,
            sample_count=BASE_SAMPLE_COUNT,
            rtol=BASE_RELATIVE_TOLERANCE,
            atol=BASE_ABSOLUTE_TOLERANCE,
        )
        for name, points, rtol, atol in (
            ("looser_tolerances", BASE_SAMPLE_COUNT, 1e-7, 1e-9),
            (
                "coarser_output_grid",
                601,
                BASE_RELATIVE_TOLERANCE,
                BASE_ABSOLUTE_TOLERANCE,
            ),
        ):
            comparison, comparison_period, _ = _sensitivity_case(
                spec,
                length=length,
                gravity=gravity,
                sample_count=points,
                rtol=rtol,
                atol=atol,
            )
            if base is None or comparison is None:
                angle_difference = None
            elif points == BASE_SAMPLE_COUNT:
                angle_difference = float(
                    np.max(np.abs(base.angle_rad - comparison.angle_rad))
                )
            else:
                stride = (BASE_SAMPLE_COUNT - 1) // (points - 1)
                angle_difference = float(
                    np.max(np.abs(base.angle_rad[::stride] - comparison.angle_rad))
                )
            period_difference = (
                None
                if base_period.period_s is None or comparison_period.period_s is None
                else comparison_period.period_s - base_period.period_s
            )
            records.append(
                NumericalSensitivity(
                    case_id=spec.case_id,
                    comparison=name,
                    base_sample_count=BASE_SAMPLE_COUNT,
                    comparison_sample_count=points,
                    maximum_angle_difference_at_shared_samples_rad=angle_difference,
                    base_period_estimate_s=base_period.period_s,
                    comparison_period_estimate_s=comparison_period.period_s,
                    period_estimate_difference_s=period_difference,
                    base_solver_success=base is not None,
                    comparison_solver_success=comparison is not None,
                )
            )
    return tuple(records)


def run_initial_velocity_research(
    *, criteria: AccuracyCriteria | None = None
) -> PendulumInitialVelocityResearch:
    """Run the reproducible Phase 26 grid and selected numerical checks."""
    mass, length, gravity = 1.0, 1.0, 9.81
    selected_criteria = criteria or AccuracyCriteria()
    omega_linear = sqrt(gravity / length)
    linear_period = linear_period_reference(length_m=length, gravity_m_per_s2=gravity)
    specs = _default_sweep(gravity, length)
    cases = tuple(
        analyze_initial_condition(
            spec,
            sweep_index=index,
            mass_kg=mass,
            length_m=length,
            gravity_m_per_s2=gravity,
            criteria=selected_criteria,
        )
        for index, spec in enumerate(specs)
    )
    return PendulumInitialVelocityResearch(
        mass_kg=mass,
        length_m=length,
        gravity_m_per_s2=gravity,
        damping_coefficient_kg_m2_per_s=0.0,
        small_angle_natural_frequency_rad_per_s=omega_linear,
        linear_period_s=linear_period,
        criteria=selected_criteria,
        sweep_order=specs,
        cases=cases,
        threshold_slices=_threshold_summaries(cases, selected_criteria),
        numerical_sensitivities=_run_sensitivities(
            specs, length=length, gravity=gravity
        ),
        assumptions=(
            "Ideal point-mass pendulum, fixed rigid length, uniform gravity.",
            "No damping or forcing; all period references assume libration.",
            "Initial angle is represented in the principal interval [-pi, pi] rad.",
        ),
        limitations=(
            "No real pendulum measurements are used; results test only the "
            "stated equations.",
            "Near-separatrix periods grow and require longer trajectories for "
            "repeated-crossing estimates.",
            "Error maxima are sampled and threshold summaries are conditional "
            "on the finite grid.",
        ),
    )


def render_initial_velocity_report(results: PendulumInitialVelocityResearch) -> str:
    """Render result, conditional-threshold, and sensitivity tables as Markdown."""
    case_table_separator = "|" + "|".join(["---"] * 20) + "|"
    lines = [
        "## Computed initial-condition sweep",
        "",
        "| # | Case | Initial angle (rad) | Initial angular velocity (rad/s) | "
        "Initial linear phase (rad) | State class | Outcome | Energy h | "
        "Separatrix gap | Turning amplitude (rad) | Amplitude increase (rad) | "
        "T nonlinear exact (s) | T nonlinear sampled (s) | T linear (s) | "
        "Signed Delta T (s) | Relative error / T nonlinear | Period crossings | "
        "Angle error at 1, 3, 5, 10 Tlinear (rad) | "
        "Phase difference at 1, 3, 5, 10 Tlinear (rad) | Criteria P/A/phase |",
        case_table_separator,
    ]
    for case in results.cases:
        window_by_period = {item.window_linear_periods: item for item in case.windows}
        angle_values = ", ".join(
            _format_window_angle(window_by_period.get(periods))
            for periods in (1, 3, 5, 10)
        )
        phase_values = ", ".join(
            _format_window_phase(window_by_period.get(periods))
            for periods in (1, 3, 5, 10)
        )
        crossings = (
            "n/a"
            if case.period_estimate is None
            else (
                f"{case.period_estimate.crossing_direction}:"
                f"{case.period_estimate.crossing_count}"
            )
        )
        criteria_text = "/".join(
            "n/a" if value is None else ("pass" if value else "fail")
            for value in (
                case.period_criterion_passed,
                case.angle_criterion_passed,
                case.phase_criterion_passed,
            )
        )
        fields = (
            str(case.sweep_index),
            case.case_id,
            f"{case.initial_angle_rad:.5g}",
            f"{case.initial_angular_velocity_rad_per_s:.5g}",
            _format_optional(case.initial_linear_phase_rad),
            case.classification.value,
            case.status.value,
            "n/a"
            if case.normalized_energy is None
            else f"{case.normalized_energy:.7g}",
            "n/a"
            if case.distance_below_separatrix is None
            else f"{case.distance_below_separatrix:.5g}",
            "n/a"
            if case.turning_point_amplitude_rad is None
            else f"{case.turning_point_amplitude_rad:.7g}",
            "n/a"
            if case.amplitude_increase_over_initial_angle_rad is None
            else f"{case.amplitude_increase_over_initial_angle_rad:.5g}",
            "n/a"
            if case.nonlinear_analytical_period_s is None
            else f"{case.nonlinear_analytical_period_s:.8g}",
            "n/a"
            if case.nonlinear_numerical_period_s is None
            else f"{case.nonlinear_numerical_period_s:.8g}",
            f"{case.linear_period_s:.8g}",
            "n/a"
            if case.signed_period_difference_s is None
            else f"{case.signed_period_difference_s:.6g}",
            "n/a"
            if case.signed_relative_period_difference_nonlinear_denominator is None
            else f"{case.signed_relative_period_difference_nonlinear_denominator:.6g}",
            crossings,
            angle_values,
            phase_values,
            criteria_text,
        )
        lines.append("| " + " | ".join(fields) + " |")
        if case.failure_reason is not None:
            lines.append(f"  Reason for `{case.case_id}`: {case.failure_reason}.")
    lines.extend(
        (
            "",
            "## Conditional grid limits",
            "",
            "| Observable | Initial angular velocity (rad/s) | "
            "Largest passing absolute initial angle (rad) | "
            "Next higher tested angle (rad) | Next point failed? | "
            "Passing cases | Failing cases |",
            "|---|---:|---:|---:|---|---|---|",
        )
    )
    for threshold in results.threshold_slices:
        largest_angle = _format_optional(
            threshold.largest_tested_passing_absolute_initial_angle_rad
        )
        next_angle = _format_optional(
            threshold.next_higher_tested_absolute_initial_angle_rad
        )
        lines.append(
            f"| {threshold.observable} | "
            f"{threshold.initial_angular_velocity_rad_per_s:.5g} | "
            f"{largest_angle} | {next_angle} | "
            f"{_format_optional_bool(threshold.next_higher_point_failed)} | "
            f"{', '.join(threshold.passing_case_ids) or 'none'} | "
            f"{', '.join(threshold.failing_case_ids) or 'none'} |"
        )
    lines.extend(
        (
            "",
            "## Numerical sensitivity",
            "",
            "| Case | Comparison | Samples base/other | "
            "Max angle difference at shared samples (rad) | Period base (s) | "
            "Period other (s) | Period estimate difference (s) | "
            "Solver success base/other |",
            "|---|---|---:|---:|---:|---:|---:|---|",
        )
    )
    for sensitivity in results.numerical_sensitivities:
        angle_difference = _format_optional(
            sensitivity.maximum_angle_difference_at_shared_samples_rad
        )
        lines.append(
            f"| {sensitivity.case_id} | {sensitivity.comparison} | "
            f"{sensitivity.base_sample_count}/{sensitivity.comparison_sample_count} | "
            f"{angle_difference} | "
            f"{_format_optional(sensitivity.base_period_estimate_s)} | "
            f"{_format_optional(sensitivity.comparison_period_estimate_s)} | "
            f"{_format_optional(sensitivity.period_estimate_difference_s)} | "
            f"{sensitivity.base_solver_success}/"
            f"{sensitivity.comparison_solver_success} |"
        )
    lines.extend(
        (
            "",
            "## Numerical settings and criteria",
            "",
            f"- Main solver: DOP853, rtol={BASE_RELATIVE_TOLERANCE:g}, "
            f"atol={BASE_ABSOLUTE_TOLERANCE:g}, "
            f"{BASE_SAMPLE_COUNT} uniform requested samples; "
            "duration is max(12 linear periods, 4 exact nonlinear periods).",
            f"- Criteria: relative period error "
            f"<= {results.criteria.relative_period_error:g}; "
            f"sampled angle error <= {results.criteria.maximum_angle_error_rad:g} "
            f"rad over {results.criteria.angle_window_linear_periods} linear periods; "
            f"absolute phase difference <= "
            f"{results.criteria.maximum_absolute_phase_difference_rad:g} rad over "
            f"{results.criteria.phase_window_linear_periods} linear periods.",
        )
    )
    return "\n".join(lines) + "\n"


def _format_optional(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.7g}"


def _format_optional_bool(value: bool | None) -> str:
    return "n/a" if value is None else str(value).lower()


def _format_window_angle(window: WindowMetrics | None) -> str:
    return (
        "n/a"
        if window is None
        else f"{window.maximum_sampled_absolute_angle_error_rad:.5g}"
    )


def _format_window_phase(window: WindowMetrics | None) -> str:
    return (
        "n/a"
        if window is None
        else f"{window.signed_accumulated_phase_difference_rad:.5g}"
    )
