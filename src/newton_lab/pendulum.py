"""Numerical simulation of a damped nonlinear pendulum."""

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.utils.validation import validate_finite_scalar

PendulumSolverMethod = Literal["RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"]


@dataclass(frozen=True, slots=True)
class DampedPendulumResult:
    """Sampled pendulum state and metadata, with angles in radians."""

    time_s: NDArray[np.float64]
    angle_rad: NDArray[np.float64]
    angular_velocity_rad_per_s: NDArray[np.float64]
    mass_kg: float
    length_m: float
    damping_coefficient_kg_m2_per_s: float
    gravity_m_per_s2: float
    initial_angle_rad: float
    initial_angular_velocity_rad_per_s: float
    method: PendulumSolverMethod
    relative_tolerance: float
    absolute_tolerance: float
    function_evaluations: int
    maximum_step_s: float | None = None


def simulate_damped_pendulum(
    mass_kg: float,
    length_m: float,
    damping_coefficient_kg_m2_per_s: float,
    gravity_m_per_s2: float,
    initial_angle_rad: float,
    initial_angular_velocity_rad_per_s: float,
    duration_s: float,
    *,
    num_points: int = 1001,
    method: PendulumSolverMethod = "RK45",
    relative_tolerance: float = 1e-8,
    absolute_tolerance: float = 1e-10,
    maximum_step_s: float | None = None,
) -> DampedPendulumResult:
    """Integrate the unforced nonlinear pendulum equation.

    The model is ``theta'' + b/(m*L**2)*theta' + g/L*sin(theta) = 0``.
    All inputs use SI units; angles are radians. The full sine term is retained
    at every amplitude. Output is uniformly sampled from zero through the
    positive ``duration_s`` endpoint, while SciPy adapts internal step sizes.

    Raises:
        ScientificValidationError: If an input is non-finite or outside its
            valid physical or numerical domain.
        IntegrationError: If SciPy reports that integration did not complete.
    """
    mass = validate_finite_scalar(mass_kg, name="mass_kg")
    length = validate_finite_scalar(length_m, name="length_m")
    damping = validate_finite_scalar(
        damping_coefficient_kg_m2_per_s,
        name="damping_coefficient_kg_m2_per_s",
    )
    gravity = validate_finite_scalar(gravity_m_per_s2, name="gravity_m_per_s2")
    initial_angle = validate_finite_scalar(initial_angle_rad, name="initial_angle_rad")
    initial_angular_velocity = validate_finite_scalar(
        initial_angular_velocity_rad_per_s,
        name="initial_angular_velocity_rad_per_s",
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
    if length <= 0.0:
        raise ScientificValidationError("length_m must be positive")
    if gravity <= 0.0:
        raise ScientificValidationError("gravity_m_per_s2 must be positive")
    if damping < 0.0:
        raise ScientificValidationError(
            "damping_coefficient_kg_m2_per_s must be non-negative"
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
        angle, angular_velocity = state
        angular_acceleration = -damping * angular_velocity / (
            mass * length**2
        ) - gravity / length * np.sin(angle)
        return angular_velocity, float(angular_acceleration)

    solution = solve_ivp(
        derivative,
        (0.0, duration),
        (initial_angle, initial_angular_velocity),
        method=method,
        t_eval=time_s,
        rtol=rtol,
        atol=atol,
        max_step=np.inf if maximum_step is None else maximum_step,
    )
    if not solution.success:
        raise IntegrationError(
            f"Damped pendulum integration failed: {solution.message}"
        )

    return DampedPendulumResult(
        time_s=solution.t.copy(),
        angle_rad=solution.y[0].copy(),
        angular_velocity_rad_per_s=solution.y[1].copy(),
        mass_kg=mass,
        length_m=length,
        damping_coefficient_kg_m2_per_s=damping,
        gravity_m_per_s2=gravity,
        initial_angle_rad=initial_angle,
        initial_angular_velocity_rad_per_s=initial_angular_velocity,
        method=method,
        relative_tolerance=rtol,
        absolute_tolerance=atol,
        function_evaluations=solution.nfev,
        maximum_step_s=maximum_step,
    )
