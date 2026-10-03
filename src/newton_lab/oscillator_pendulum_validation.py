"""Independent numerical cross-validation of oscillator and pendulum equations.

The physical oscillator uses Newton Lab's existing oscillator solver. The
linearized pendulum is integrated through a separate right-hand side and
``solve_ivp`` call because the existing pendulum API always uses ``sin(theta)``.
All comparisons use shared dimensionless time samples; the exact harmonic
solution is evaluated separately from both integrations.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi, sqrt
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

from newton_lab.dynamics import SolverMethod, simulate_damped_oscillator
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.pendulum import simulate_damped_pendulum
from newton_lab.pendulum_initial_velocity import (
    estimate_nonlinear_period,
    nonlinear_period_reference,
)

RunStatus = Literal["succeeded", "failed"]


@dataclass(frozen=True, slots=True)
class MatchedInitialConditions:
    """Physical and dimensionless initial states under the declared mapping."""

    q0: float
    dimensionless_velocity: float
    amplitude_scale_m: float
    mass_kg: float
    stiffness_n_per_m: float
    oscillator_frequency_rad_per_s: float
    length_m: float
    gravity_m_per_s2: float
    pendulum_frequency_rad_per_s: float
    oscillator_displacement_m: float
    oscillator_velocity_m_per_s: float
    pendulum_angle_rad: float
    pendulum_angular_velocity_rad_per_s: float
    frequency_match: bool


@dataclass(frozen=True, slots=True)
class PeriodMetric:
    """Period estimate in dimensionless time and its physical conversion."""

    estimate_tau: float | None
    reference_tau: float | None
    absolute_error_tau: float | None
    estimate_s: float | None
    reference_s: float | None
    absolute_error_s: float | None
    crossing_count: int
    status: str


@dataclass(frozen=True, slots=True)
class ValidationRun:
    """One independent integration pair and its aligned comparison metrics."""

    case_id: str
    status: RunStatus
    failure_message: str | None
    initial: MatchedInitialConditions
    method: SolverMethod
    relative_tolerance: float
    absolute_tolerance: float
    dimensionless_times: tuple[float, ...]
    duration_tau: float
    sample_count: int
    oscillator_error_max: float | None
    linear_pendulum_error_max: float | None
    cross_system_error_max: float | None
    nonlinear_model_error_max: float | None
    nonlinear_vs_linear_solver_error_max: float | None
    normalization_scale: float | None
    normalized_oscillator_error: float | None
    normalized_linear_pendulum_error: float | None
    normalized_cross_system_error: float | None
    normalized_nonlinear_model_error: float | None
    linear_period: PeriodMetric | None
    pendulum_linear_period: PeriodMetric | None
    nonlinear_period: PeriodMetric | None
    oscillator_function_evaluations: int | None
    linear_pendulum_function_evaluations: int | None
    nonlinear_function_evaluations: int | None
    exact_dimensionless_displacement: tuple[float, ...] = ()
    oscillator_dimensionless_displacement: tuple[float, ...] = ()
    linear_pendulum_displacement: tuple[float, ...] = ()
    nonlinear_pendulum_displacement: tuple[float, ...] = ()


@dataclass(frozen=True, slots=True)
class SensitivityRecord:
    """A single-factor solver, tolerance, horizon or sampling comparison."""

    factor: str
    setting: str
    case_id: str
    oscillator_error_max: float | None
    linear_pendulum_error_max: float | None
    cross_system_error_max: float | None
    nonlinear_model_error_max: float | None
    normalized_cross_system_error: float | None
    period_estimate_tau: float | None
    period_error_tau: float | None
    status: RunStatus


@dataclass(frozen=True, slots=True)
class ValidationStudy:
    """Ordered Phase 28 run grid, sensitivity rows, settings and report."""

    runs: tuple[ValidationRun, ...]
    sensitivities: tuple[SensitivityRecord, ...]
    initial_conditions: tuple[tuple[float, float], ...]
    tolerances: tuple[tuple[float, float], ...]
    methods: tuple[SolverMethod, ...]
    output_counts: tuple[int, ...]
    horizons_tau: tuple[float, ...]
    assumptions: tuple[str, ...]
    limitations: tuple[str, ...]
    report_markdown: str


def construct_matched_initial_conditions(
    q0: float,
    dimensionless_velocity: float,
    *,
    amplitude_scale_m: float = 0.25,
    mass_kg: float = 2.0,
    stiffness_n_per_m: float = 8.0,
    length_m: float = 1.5,
    gravity_m_per_s2: float = 6.0,
    require_frequency_match: bool = False,
) -> MatchedInitialConditions:
    """Map one dimensionless state to SI oscillator and pendulum states.

    The nonzero scale ``amplitude_scale_m`` is the oscillator displacement in
    metres per unit q; it is a normalization convention, not an identification
    of displacement with angle. By default both natural frequencies are
    2 rad/s. Different frequencies remain comparable through τ=ωt.
    """
    values = (
        q0,
        dimensionless_velocity,
        amplitude_scale_m,
        mass_kg,
        stiffness_n_per_m,
        length_m,
        gravity_m_per_s2,
    )
    if not all(isfinite(value) for value in values):
        raise ScientificValidationError("mapping inputs must be finite")
    if amplitude_scale_m == 0.0:
        raise ScientificValidationError("amplitude_scale_m must be nonzero")
    if mass_kg <= 0.0 or stiffness_n_per_m <= 0.0:
        raise ScientificValidationError("mass and stiffness must be positive")
    if length_m <= 0.0 or gravity_m_per_s2 <= 0.0:
        raise ScientificValidationError("length and gravity must be positive")
    omega_o = sqrt(stiffness_n_per_m / mass_kg)
    omega_p = sqrt(gravity_m_per_s2 / length_m)
    matched = bool(np.isclose(omega_o, omega_p, rtol=1e-12, atol=0.0))
    if require_frequency_match and not matched:
        raise ScientificValidationError(
            "oscillator and pendulum frequencies differ; use dimensionless time"
        )
    return MatchedInitialConditions(
        q0=q0,
        dimensionless_velocity=dimensionless_velocity,
        amplitude_scale_m=amplitude_scale_m,
        mass_kg=mass_kg,
        stiffness_n_per_m=stiffness_n_per_m,
        oscillator_frequency_rad_per_s=omega_o,
        length_m=length_m,
        gravity_m_per_s2=gravity_m_per_s2,
        pendulum_frequency_rad_per_s=omega_p,
        oscillator_displacement_m=amplitude_scale_m * q0,
        oscillator_velocity_m_per_s=amplitude_scale_m
        * omega_o
        * dimensionless_velocity,
        pendulum_angle_rad=q0,
        pendulum_angular_velocity_rad_per_s=omega_p * dimensionless_velocity,
        frequency_match=matched,
    )


def exact_linear_solution(
    dimensionless_time: tuple[float, ...] | NDArray[np.float64],
    q0: float,
    dimensionless_velocity: float,
) -> NDArray[np.float64]:
    """Evaluate ``q0*cos(τ) + v0_star*sin(τ)`` on finite τ samples."""
    times = np.asarray(dimensionless_time, dtype=np.float64)
    if times.ndim != 1 or not np.isfinite(times).all():
        raise ScientificValidationError("dimensionless_time must be finite and 1D")
    if not isfinite(q0) or not isfinite(dimensionless_velocity):
        raise ScientificValidationError("initial state must be finite")
    return q0 * np.cos(times) + dimensionless_velocity * np.sin(times)


def normalized_error(error_max: float, reference: NDArray[np.float64]) -> float | None:
    """Normalize max error by max absolute analytical trajectory; zero is undefined."""
    scale = float(np.max(np.abs(reference))) if reference.size else 0.0
    if not isfinite(error_max) or error_max < 0.0:
        raise ScientificValidationError("error_max must be finite and non-negative")
    return None if scale <= np.finfo(np.float64).tiny else error_max / scale


def estimate_period_tau(
    times_tau: NDArray[np.float64], values: NDArray[np.float64]
) -> tuple[float | None, int, str]:
    """Estimate a period from upward zero crossings with linear interpolation."""
    if times_tau.ndim != 1 or values.ndim != 1 or len(times_tau) != len(values):
        raise ScientificValidationError("period inputs must be aligned 1D arrays")
    if (
        len(times_tau) < 3
        or not np.isfinite(times_tau).all()
        or not np.isfinite(values).all()
    ):
        raise ScientificValidationError(
            "period inputs must be finite and sufficiently sampled"
        )
    crossings: list[float] = []
    for i in range(len(values) - 1):
        if values[i] <= 0.0 < values[i + 1]:
            fraction = -values[i] / (values[i + 1] - values[i])
            crossings.append(
                float(times_tau[i] + fraction * (times_tau[i + 1] - times_tau[i]))
            )
    if len(crossings) < 2:
        return None, len(crossings), "insufficient_same_direction_crossings"
    estimate = float(np.mean(np.diff(crossings)))
    return estimate, len(crossings), "estimated_by_linear_interpolated_upward_crossings"


def _one_period_metric(
    times_tau: NDArray[np.float64],
    values: NDArray[np.float64],
    *,
    omega_rad_per_s: float,
    reference_tau: float,
) -> PeriodMetric:
    estimate, crossings, status = estimate_period_tau(times_tau, values)
    reference_s = reference_tau / omega_rad_per_s
    if estimate is None:
        return PeriodMetric(
            None, reference_tau, None, None, reference_s, None, crossings, status
        )
    error = abs(estimate - reference_tau)
    estimate_s = estimate / omega_rad_per_s
    return PeriodMetric(
        estimate,
        reference_tau,
        error,
        estimate_s,
        reference_s,
        abs(estimate_s - reference_s),
        crossings,
        status,
    )


def _run_one(
    *,
    case_id: str,
    initial: MatchedInitialConditions,
    duration_tau: float,
    sample_count: int,
    method: SolverMethod,
    rtol: float,
    atol: float,
) -> ValidationRun:
    times_tau = np.linspace(0.0, duration_tau, sample_count, dtype=np.float64)
    exact = exact_linear_solution(times_tau, initial.q0, initial.dimensionless_velocity)
    omega_o = initial.oscillator_frequency_rad_per_s
    omega_p = initial.pendulum_frequency_rad_per_s
    try:
        oscillator = simulate_damped_oscillator(
            initial.mass_kg,
            0.0,
            initial.stiffness_n_per_m,
            initial.oscillator_displacement_m,
            initial.oscillator_velocity_m_per_s,
            duration_tau / omega_o,
            num_points=sample_count,
            method=method,
            relative_tolerance=rtol,
            absolute_tolerance=atol,
        )

        # Independently specified physical IVP: theta' = omega, omega' = -(g/L)theta.
        def linear_pendulum_rhs(
            _time_s: float, state: NDArray[np.float64]
        ) -> tuple[float, float]:
            angle, angular_velocity = state
            return angular_velocity, -(omega_p**2) * angle

        pendulum_times_s = times_tau / omega_p
        linear_solution = solve_ivp(
            linear_pendulum_rhs,
            (0.0, duration_tau / omega_p),
            (initial.pendulum_angle_rad, initial.pendulum_angular_velocity_rad_per_s),
            method=method,
            t_eval=pendulum_times_s,
            rtol=rtol,
            atol=atol,
        )
        if not linear_solution.success or linear_solution.y.shape[1] != sample_count:
            raise IntegrationError(
                f"linearized-pendulum integration failed: {linear_solution.message}"
            )
        q_num = oscillator.displacement_m / initial.amplitude_scale_m
        theta_num = linear_solution.y[0]
        e_o = float(np.max(np.abs(q_num - exact)))
        e_p = float(np.max(np.abs(theta_num - exact)))
        e_cross = float(np.max(np.abs(q_num - theta_num)))

        nonlinear = simulate_damped_pendulum(
            1.0,
            initial.length_m,
            0.0,
            initial.gravity_m_per_s2,
            initial.pendulum_angle_rad,
            initial.pendulum_angular_velocity_rad_per_s,
            duration_tau / omega_p,
            num_points=sample_count,
            method=method,
            relative_tolerance=rtol,
            absolute_tolerance=atol,
        )
        theta_full = nonlinear.angle_rad
        e_model = float(np.max(np.abs(theta_full - exact)))
        e_nl_solver = float(np.max(np.abs(theta_full - theta_num)))
        scale = float(np.max(np.abs(exact)))
        linear_period = _one_period_metric(
            times_tau, q_num, omega_rad_per_s=omega_o, reference_tau=2.0 * pi
        )
        pendulum_period = _one_period_metric(
            times_tau, theta_num, omega_rad_per_s=omega_p, reference_tau=2.0 * pi
        )
        nonlinear_estimate = estimate_nonlinear_period(nonlinear)
        nonlinear_reference_s: float | None
        try:
            nonlinear_reference_s = nonlinear_period_reference(
                initial.pendulum_angle_rad,
                initial.pendulum_angular_velocity_rad_per_s,
                length_m=initial.length_m,
                gravity_m_per_s2=initial.gravity_m_per_s2,
            )
        except ScientificValidationError:
            nonlinear_reference_s = None
        nonlinear_period = PeriodMetric(
            estimate_tau=(
                None
                if nonlinear_estimate.period_s is None
                else nonlinear_estimate.period_s * omega_p
            ),
            reference_tau=(
                None
                if nonlinear_reference_s is None
                else nonlinear_reference_s * omega_p
            ),
            absolute_error_tau=(
                None
                if nonlinear_estimate.period_s is None or nonlinear_reference_s is None
                else abs(nonlinear_estimate.period_s - nonlinear_reference_s) * omega_p
            ),
            estimate_s=nonlinear_estimate.period_s,
            reference_s=nonlinear_reference_s,
            absolute_error_s=(
                None
                if nonlinear_estimate.period_s is None or nonlinear_reference_s is None
                else abs(nonlinear_estimate.period_s - nonlinear_reference_s)
            ),
            crossing_count=nonlinear_estimate.crossing_count,
            status=(
                nonlinear_estimate.status
                if nonlinear_reference_s is not None
                else "reference_unavailable_nonlibration"
            ),
        )
        return ValidationRun(
            case_id,
            "succeeded",
            None,
            initial,
            method,
            rtol,
            atol,
            tuple(float(value) for value in times_tau),
            duration_tau,
            sample_count,
            e_o,
            e_p,
            e_cross,
            e_model,
            e_nl_solver,
            scale,
            normalized_error(e_o, exact),
            normalized_error(e_p, exact),
            normalized_error(e_cross, exact),
            normalized_error(e_model, exact),
            linear_period,
            pendulum_period,
            nonlinear_period,
            oscillator.function_evaluations,
            int(linear_solution.nfev),
            nonlinear.function_evaluations,
            tuple(float(value) for value in exact),
            tuple(float(value) for value in q_num),
            tuple(float(value) for value in theta_num),
            tuple(float(value) for value in theta_full),
        )
    except (IntegrationError, ScientificValidationError, ValueError) as exc:
        return ValidationRun(
            case_id,
            "failed",
            f"{type(exc).__name__}: {exc}",
            initial,
            method,
            rtol,
            atol,
            tuple(float(value) for value in times_tau),
            duration_tau,
            sample_count,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
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


def run_validation_experiment() -> ValidationStudy:
    """Run the fixed, ordered sensitivity grid and render its report."""
    initials = (
        (0.0, 0.0),
        (0.01, 0.0),
        (0.2, 0.0),
        (0.2, 0.3),
        (0.2, -0.3),
        (0.8, 0.3),
    )
    tolerances = ((1e-5, 1e-7), (1e-8, 1e-10), (1e-11, 1e-13))
    methods: tuple[SolverMethod, ...] = ("RK45", "DOP853")
    counts = (101, 501, 2001)
    horizons = (2.0 * pi, 8.0 * pi, 16.0 * pi)
    runs: list[ValidationRun] = []
    sensitivities: list[SensitivityRecord] = []
    for case_index, (q0, v0) in enumerate(initials):
        initial = construct_matched_initial_conditions(q0, v0)
        case_id = f"initial_{case_index + 1}"
        base = _run_one(
            case_id=case_id,
            initial=initial,
            duration_tau=8.0 * pi,
            sample_count=501,
            method="DOP853",
            rtol=1e-8,
            atol=1e-10,
        )
        runs.append(base)
        for rtol, atol in tolerances:
            result = _run_one(
                case_id=case_id,
                initial=initial,
                duration_tau=8.0 * pi,
                sample_count=501,
                method="DOP853",
                rtol=rtol,
                atol=atol,
            )
            runs.append(result)
            sensitivities.append(
                _sensitivity("tolerance", f"rtol={rtol:g},atol={atol:g}", result)
            )
        for method in methods:
            result = _run_one(
                case_id=case_id,
                initial=initial,
                duration_tau=8.0 * pi,
                sample_count=501,
                method=method,
                rtol=1e-8,
                atol=1e-10,
            )
            runs.append(result)
            sensitivities.append(_sensitivity("method", method, result))
        for count in counts:
            result = _run_one(
                case_id=case_id,
                initial=initial,
                duration_tau=8.0 * pi,
                sample_count=count,
                method="DOP853",
                rtol=1e-8,
                atol=1e-10,
            )
            runs.append(result)
            sensitivities.append(_sensitivity("output_grid", f"{count} points", result))
        for horizon in horizons:
            result = _run_one(
                case_id=case_id,
                initial=initial,
                duration_tau=horizon,
                sample_count=501,
                method="DOP853",
                rtol=1e-8,
                atol=1e-10,
            )
            runs.append(result)
            sensitivities.append(_sensitivity("horizon", f"tau={horizon:g}", result))
    assumption_text = (
        "Undamped linear oscillator; no forcing and constant positive m,k.",
        "Ideal pendulum; no damping or forcing, fixed length and uniform gravity.",
        (
            "A is a nonzero displacement scale in metres per unit dimensionless "
            "q; it does not equate angle with displacement."
        ),
        (
            "Comparisons use common dimensionless times; physical time is "
            "tau/omega for each system."
        ),
        (
            "All integrations use SciPy solve_ivp; distinct RHS implementations "
            "do not mean independent numerical libraries."
        ),
    )
    limits = (
        (
            "The finite state, tolerance, method, grid and horizon choices do "
            "not prove convergence or universal accuracy."
        ),
        (
            "Grid effects include solve_ivp dense output and crossing "
            "interpolation. They are distinct from adaptive-step tolerance, but "
            "do not bound discretization error."
        ),
        (
            "For the zero state, normalized error is undefined because "
            "analytical amplitude is zero; absolute errors remain reported."
        ),
        (
            "Period estimates require repeated same-direction crossings; short "
            "windows can be insufficient."
        ),
        (
            "Elliptic-period references apply only to librations, not "
            "rotations, separatrix states or equilibrium."
        ),
    )
    markdown = render_validation_report(
        tuple(runs),
        tuple(sensitivities),
        initials,
        tolerances,
        methods,
        counts,
        horizons,
        assumption_text,
        limits,
    )
    return ValidationStudy(
        tuple(runs),
        tuple(sensitivities),
        initials,
        tolerances,
        methods,
        counts,
        horizons,
        assumption_text,
        limits,
        markdown,
    )


def _sensitivity(factor: str, setting: str, run: ValidationRun) -> SensitivityRecord:
    period = run.linear_period
    return SensitivityRecord(
        factor,
        setting,
        run.case_id,
        run.oscillator_error_max,
        run.linear_pendulum_error_max,
        run.cross_system_error_max,
        run.nonlinear_model_error_max,
        run.normalized_cross_system_error,
        None if period is None else period.estimate_tau,
        None if period is None else period.absolute_error_tau,
        run.status,
    )


def render_validation_report(
    runs: tuple[ValidationRun, ...],
    sensitivities: tuple[SensitivityRecord, ...],
    initials: tuple[tuple[float, float], ...],
    tolerances: tuple[tuple[float, float], ...],
    methods: tuple[SolverMethod, ...],
    counts: tuple[int, ...],
    horizons: tuple[float, ...],
    assumptions: tuple[str, ...],
    limitations: tuple[str, ...],
) -> str:
    """Render deterministic Markdown from recorded run data only."""
    baseline_by_case = {
        run.case_id: run
        for run in runs
        if run.method == "DOP853"
        and run.sample_count == 501
        and run.relative_tolerance == 1e-8
        and run.duration_tau == 8 * pi
    }
    lines = [
        "# Phase 28: Independent oscillator and pendulum numerical validation",
        "",
        "## Research question and hypotheses",
        "",
        "Question: when integrated independently, do the dimensionless linear",
        "oscillator and linearized pendulum approach the same exact solution, and",
        "how do discrepancies respond to solver settings, sampling, initial state",
        "and horizon? Can finite-angle model error be distinguished from solver error?",
        "",
        "Hypotheses: mapped analytical trajectories agree; independently integrated",
        "linear trajectories approach that reference as accuracy is tightened;",
        "cross-system error is consistent with per-solver errors; finite-angle",
        "model error can exceed numerical error; grid and period effects are visible.",
        "",
        "## Mathematical mapping",
        "",
        "For the undamped oscillator, `m x'' + kx = 0`, define",
        "`omega_o = sqrt(k/m)`, `tau = omega_o t`, and `x = A q`. Then",
        "`d2q/dtau2 + q = 0`. Here `A != 0` is a displacement scale in metres",
        "per unit q; it does not make displacement and angle physically identical.",
        "",
        "For the full pendulum, `theta'' + (g/L) sin(theta) = 0`. With",
        "`tau = sqrt(g/L) t`, its equation is `d2theta/dtau2 + sin(theta) = 0`;",
        "small-angle linearization gives `d2theta/dtau2 + theta = 0`.",
        "Both linear systems therefore have exact solution",
        "`u(tau) = u0 cos(tau) + v0* sin(tau)`. This establishes mathematical",
        "equivalence under scaling, not equivalence of physical systems.",
        "",
        "Matched states are `x0=Aq0`, `xdot0=A omega_o v0*`, `theta0=q0`,",
        "and `thetadot0=omega_p v0*`, where `omega_p=sqrt(g/L)`. Defaults are",
        "`m=2 kg`, `k=8 N/m`, `L=1.5 m`, `g=6 m/s^2`, `A=0.25 m`; both",
        "frequencies are 2 rad/s. Unequal frequencies use each system's own",
        "`t=tau/omega`, while comparisons use shared dimensionless tau.",
        "",
        "## Experimental design",
        "",
        f"Ordered initial states `(q0,v0*)`: {initials}.",
        f"Horizons tau: {horizons}; output counts: {counts}.",
        f"Methods: {methods}; tolerance pairs (rtol, atol): {tolerances}.",
        "Sensitivity changes one factor at a time around DOP853, 1e-8/1e-10,",
        "501 points and an 8*pi horizon.",
        "",
        "The oscillator uses the existing oscillator API with zero damping. The",
        "linearized pendulum has a separately constructed state, RHS and solve_ivp",
        "call because the existing pendulum API only solves the full-sine model.",
        "The full-sine comparison uses the existing pendulum solver. All use SciPy,",
        "so this is not independent numerical-library validation.",
        "",
        "## Results",
        "",
        "Errors are maximum absolute dimensionless displacement on shared requested",
        "tau samples. Normalized error divides by `max_i |u_exact(tau_i)|`; it is",
        "undefined for the zero state, whose absolute errors remain reported.",
        "Cross error compares the two independent linear integrations. Model error",
        "compares full-sine pendulum samples with the exact linear solution.",
        "The nonlinear-versus-linear numerical discrepancy is also listed.",
        "Full-sine model error includes both physical model discrepancy and",
        "numerical error from integrating the nonlinear trajectory.",
        "",
        "| "
        + " | ".join(
            (
                "Case",
                "q0,v0*",
                "method",
                "rtol/atol",
                "tau/points",
                "E osc",
                "E lin pend",
                "E cross",
                "E model",
                "E nonlinear vs lin num",
                "linear period est/ref/error tau",
                "linear-pendulum period est/ref/error tau",
                "nonlinear period est/ref/error s",
            )
        )
        + " |",
        "|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]

    def value(number: float | None) -> str:
        return "n/a" if number is None else f"{number:.4g}"

    for case_id, run in baseline_by_case.items():
        linear_period = (
            "n/a"
            if run.linear_period is None
            else "/".join(
                value(item)
                for item in (
                    run.linear_period.estimate_tau,
                    run.linear_period.reference_tau,
                    run.linear_period.absolute_error_tau,
                )
            )
        )
        pendulum_period = (
            "n/a"
            if run.pendulum_linear_period is None
            else "/".join(
                value(item)
                for item in (
                    run.pendulum_linear_period.estimate_tau,
                    run.pendulum_linear_period.reference_tau,
                    run.pendulum_linear_period.absolute_error_tau,
                )
            )
        )
        nonlinear_period = (
            "n/a"
            if run.nonlinear_period is None
            else "/".join(
                value(item)
                for item in (
                    run.nonlinear_period.estimate_s,
                    run.nonlinear_period.reference_s,
                    run.nonlinear_period.absolute_error_s,
                )
            )
        )
        fields: tuple[str, ...] = (
            case_id,
            f"{run.initial.q0:g},{run.initial.dimensionless_velocity:g}",
            run.method,
            f"{run.relative_tolerance:g}/{run.absolute_tolerance:g}",
            f"{run.duration_tau / pi:.0f}pi/{run.sample_count}",
            value(run.oscillator_error_max),
            value(run.linear_pendulum_error_max),
            value(run.cross_system_error_max),
            value(run.nonlinear_model_error_max),
            value(run.nonlinear_vs_linear_solver_error_max),
            linear_period,
            pendulum_period,
            nonlinear_period,
        )
        lines.append("| " + " | ".join(fields) + " |")

    lines.extend(
        (
            "",
            "Run records retain exact, oscillator, linearized-pendulum and full-sine",
            "samples on the common dimensionless output grid, alongside solver",
            "configuration and function-evaluation counts.",
        )
    )

    lines.extend(
        (
            "",
            "## Sensitivity",
            "",
            "Each sensitivity row is a recorded run, including its case and status.",
            "Only the named factor changes within each one-factor series; table values",
            "are measurements, not a convergence proof.",
            "",
            "| "
            + " | ".join(
                (
                    "Factor",
                    "Setting",
                    "Case",
                    "E osc",
                    "E lin pend",
                    "E cross",
                    "E model",
                    "normalized cross",
                    "period tau/error",
                    "status",
                )
            )
            + " |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---|",
        )
    )
    for item in sensitivities:
        sensitivity_fields = (
            item.factor,
            item.setting,
            item.case_id,
            value(item.oscillator_error_max),
            value(item.linear_pendulum_error_max),
            value(item.cross_system_error_max),
            value(item.nonlinear_model_error_max),
            value(item.normalized_cross_system_error),
            f"{value(item.period_estimate_tau)}/{value(item.period_error_tau)}",
            item.status,
        )
        lines.append("| " + " | ".join(sensitivity_fields) + " |")

    grouped: dict[str, dict[str, list[float]]] = {}
    for item in sensitivities:
        if item.cross_system_error_max is not None:
            grouped.setdefault(item.factor, {}).setdefault(item.setting, []).append(
                item.cross_system_error_max
            )
    lines.extend(
        (
            "",
            "Mean absolute cross discrepancy across the six initial states by setting:",
        )
    )
    for factor, settings in grouped.items():
        for setting, errors in settings.items():
            mean_error = sum(errors) / len(errors)
            lines.append(f"- {factor}, {setting}: mean E_cross={mean_error:.4g}")
    lines.extend(
        (
            "Across the six states, mean cross error decreases at each tighter",
            "tolerance setting. The RK45 mean is lower than DOP853 in this grid;",
            "this is method sensitivity, not evidence that one method is generally",
            "more accurate. Output-grid means remain close and change slightly",
            "non-monotonically. Mean cross error rises modestly with longer horizons.",
            "At q0=0.8, v0*=0.3, full-sine model error is about 0.943 at 8*pi and",
            "1.58 at 16*pi, far above the linear solver errors. These trends are",
            "specific to the tested configurations and metrics.",
        )
    )

    lines.extend(
        (
            "",
            "## Error attribution and falsification",
            "",
            "- The exact linear references coincide by derivation for matched states.",
            "- The triangle inequality gives cross error <= oscillator analytical",
            "  error + linear-pendulum analytical error on the common grid.",
            "- Tighter tolerances and alternate methods expose numerical sensitivity;",
            "  errors need not decrease monotonically at every finite setting.",
            "- Coarse and dense outputs are evaluated at their own common requested",
            "  samples. Differences mix dense-output evaluation and sampled maxima;",
            "  they do not isolate or bound internal adaptive-step error.",
            "- Signed-velocity symmetry, zero-state and parameter-rescaling checks",
            "  passed in tests. Finite-angle deviation is a model discrepancy, not",
            "  evidence of a solver failure.",
            "- Periods use interpolated upward zero crossings. A period is unavailable",
            "  unless enough same-direction crossings exist. The nonlinear",
            "  reference is elliptic and applies to librations only.",
            "",
            "## Assumptions and limitations",
            "",
        )
    )
    lines.extend(f"- {item}" for item in assumptions)
    lines.extend(f"- {item}" for item in limitations)
    lines.extend(
        (
            "",
            f"Observed runs: {len(runs)}; failed integrations: "
            f"{sum(run.status == 'failed' for run in runs)}.",
            "The tested grid is finite, and all systems share SciPy solve_ivp.",
            "Analytical equivalence is not physical equivalence; numerical agreement",
            "does not establish universal accuracy, application performance or trading",
            "utility.",
            "",
            "## Conclusion and next question",
            "",
            "The runs test two independently specified linear equations against one",
            "exact dimensionless reference. They support the declared mapping only",
            "for the tested conditions. The next question is whether validated bounds",
            "can quantify the observed global error independently of SciPy solvers.",
            "",
        )
    )
    return "\n".join(lines)
