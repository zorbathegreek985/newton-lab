"""Numerical models for simple dynamical systems."""

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.utils.validation import validate_finite_scalar

SolverMethod = Literal["RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"]


@dataclass(frozen=True, slots=True)
class DampedOscillatorResult:
    """Sampled state and numerical settings for an oscillator simulation.

    Arrays contain SI values and have one entry per requested output time.
    """

    time_s: NDArray[np.float64]
    displacement_m: NDArray[np.float64]
    velocity_m_per_s: NDArray[np.float64]
    mass_kg: float
    damping_coefficient_kg_per_s: float
    stiffness_n_per_m: float
    initial_displacement_m: float
    initial_velocity_m_per_s: float
    method: SolverMethod
    relative_tolerance: float
    absolute_tolerance: float
    function_evaluations: int
    maximum_step_s: float | None = None


def simulate_damped_oscillator(
    mass_kg: float,
    damping_coefficient_kg_per_s: float,
    stiffness_n_per_m: float,
    initial_displacement_m: float,
    initial_velocity_m_per_s: float,
    duration_s: float,
    *,
    num_points: int = 1001,
    method: SolverMethod = "RK45",
    relative_tolerance: float = 1e-8,
    absolute_tolerance: float = 1e-10,
    maximum_step_s: float | None = None,
) -> DampedOscillatorResult:
    """Integrate a free, linearly damped harmonic oscillator.

    The model is ``m*x'' + c*x' + k*x = 0`` with constant SI parameters, a
    linear spring, linear viscous damping, and no external force. ``duration_s``
    is the positive end time; output is sampled uniformly from zero to that
    time, including both endpoints. SciPy's adaptive solver controls internal
    steps, while the tolerances control its local error estimate.

    Raises:
        ScientificValidationError: If a parameter, initial condition, or
            numerical setting is outside its valid domain.
        IntegrationError: If the numerical integrator reports failure.
    """
    mass = validate_finite_scalar(mass_kg, name="mass_kg")
    damping = validate_finite_scalar(
        damping_coefficient_kg_per_s, name="damping_coefficient_kg_per_s"
    )
    stiffness = validate_finite_scalar(stiffness_n_per_m, name="stiffness_n_per_m")
    initial_displacement = validate_finite_scalar(
        initial_displacement_m, name="initial_displacement_m"
    )
    initial_velocity = validate_finite_scalar(
        initial_velocity_m_per_s, name="initial_velocity_m_per_s"
    )
    duration = validate_finite_scalar(duration_s, name="duration_s")
    rtol = validate_finite_scalar(relative_tolerance, name="relative_tolerance")
    atol = validate_finite_scalar(absolute_tolerance, name="absolute_tolerance")
    maximum_step = (
        None
        if maximum_step_s is None
        else validate_finite_scalar(maximum_step_s, name="maximum_step_s")
    )

    if mass <= 0.0:
        raise ScientificValidationError("mass_kg must be positive")
    if stiffness <= 0.0:
        raise ScientificValidationError("stiffness_n_per_m must be positive")
    if damping < 0.0:
        raise ScientificValidationError(
            "damping_coefficient_kg_per_s must be non-negative"
        )
    if duration <= 0.0:
        raise ScientificValidationError("duration_s must be positive")
    if (
        isinstance(num_points, bool)
        or not isinstance(num_points, int)
        or num_points < 2
    ):
        raise ScientificValidationError("num_points must be an integer of at least 2")
    if rtol <= 0.0:
        raise ScientificValidationError("relative_tolerance must be positive")
    if atol <= 0.0:
        raise ScientificValidationError("absolute_tolerance must be positive")
    if maximum_step is not None and maximum_step <= 0.0:
        raise ScientificValidationError("maximum_step_s must be positive")
    if method not in ("RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"):
        raise ScientificValidationError(f"unsupported integration method: {method!r}")

    time_s = np.linspace(0.0, duration, num_points, dtype=np.float64)

    def derivative(_time_s: float, state: NDArray[np.float64]) -> tuple[float, float]:
        displacement, velocity = state
        acceleration = -(damping * velocity + stiffness * displacement) / mass
        return velocity, acceleration

    solution = solve_ivp(
        derivative,
        (0.0, duration),
        (initial_displacement, initial_velocity),
        method=method,
        t_eval=time_s,
        rtol=rtol,
        atol=atol,
        max_step=np.inf if maximum_step is None else maximum_step,
    )
    if not solution.success:
        raise IntegrationError(
            f"Damped oscillator integration failed: {solution.message}"
        )

    return DampedOscillatorResult(
        time_s=solution.t.copy(),
        displacement_m=solution.y[0].copy(),
        velocity_m_per_s=solution.y[1].copy(),
        mass_kg=mass,
        damping_coefficient_kg_per_s=damping,
        stiffness_n_per_m=stiffness,
        initial_displacement_m=initial_displacement,
        initial_velocity_m_per_s=initial_velocity,
        method=method,
        relative_tolerance=rtol,
        absolute_tolerance=atol,
        function_evaluations=solution.nfev,
        maximum_step_s=maximum_step,
    )
