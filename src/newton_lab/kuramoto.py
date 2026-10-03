"""Phase-reduced all-to-all Kuramoto oscillator model."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

from newton_lab.exceptions import IntegrationError, ScientificValidationError

FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]
SolverMethod = Literal["RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"]


@dataclass(frozen=True, slots=True)
class KuramotoResult:
    """Numerical phase trajectory and its collective order parameter."""

    time_s: FloatArray
    phases_rad: FloatArray
    order_parameter: ComplexArray
    coherence: FloatArray
    coupling_rad_per_s: float
    function_evaluations: int
    method: str
    relative_tolerance: float
    absolute_tolerance: float


def order_parameter(phases_rad: FloatArray) -> complex:
    """Return the complex Kuramoto order parameter for one phase snapshot."""
    phases = np.asarray(phases_rad, dtype=np.float64)
    if phases.ndim != 1 or phases.size == 0 or not np.isfinite(phases).all():
        raise ScientificValidationError(
            "phases_rad must be a nonempty, finite one-dimensional array"
        )
    return complex(np.mean(np.exp(1j * phases)))


def simulate_kuramoto(
    natural_frequencies_rad_per_s: FloatArray,
    initial_phases_rad: FloatArray,
    coupling_rad_per_s: float,
    output_times_s: FloatArray,
    *,
    method: SolverMethod = "DOP853",
    relative_tolerance: float = 1e-8,
    absolute_tolerance: float = 1e-10,
    maximum_step_s: float = math.inf,
) -> KuramotoResult:
    """Integrate the all-to-all sinusoidally coupled phase-oscillator model.

    The equations are ``d theta_i/dt = omega_i + K/N * sum_j
    sin(theta_j-theta_i)``. Phases are in radians, natural frequencies and
    coupling in radians per second, and time in seconds. The mean-field
    implementation evaluates coupling in O(N) per derivative call using the
    complex order parameter. Inputs are copied before integration.

    This is a phase-reduced model appropriate when oscillators remain close to
    stable limit cycles and coupling is weak enough for phase reduction. It
    assumes fixed frequencies, identical sinusoidal all-to-all interactions,
    and no forcing, delays, noise, or adaptive links.

    Raises:
        ScientificValidationError: If arrays, coupling, or solver settings
            are invalid.
        IntegrationError: If SciPy reports an unsuccessful integration.
    """
    frequencies = np.array(natural_frequencies_rad_per_s, dtype=np.float64, copy=True)
    phases = np.array(initial_phases_rad, dtype=np.float64, copy=True)
    times = np.array(output_times_s, dtype=np.float64, copy=True)
    coupling = float(coupling_rad_per_s)
    rtol = float(relative_tolerance)
    atol = float(absolute_tolerance)
    maximum_step = float(maximum_step_s)

    if (
        frequencies.ndim != 1
        or frequencies.size < 2
        or not np.isfinite(frequencies).all()
    ):
        raise ScientificValidationError(
            "natural frequencies must be a finite one-dimensional array of size >= 2"
        )
    if (
        phases.ndim != 1
        or phases.size != frequencies.size
        or not np.isfinite(phases).all()
    ):
        raise ScientificValidationError(
            "initial phases must be finite and match the frequency array"
        )
    if (
        times.ndim != 1
        or times.size < 2
        or not np.isfinite(times).all()
        or times[0] != 0.0
        or np.any(np.diff(times) <= 0.0)
    ):
        raise ScientificValidationError(
            "output times must be finite, increasing, and begin at zero"
        )
    if not math.isfinite(coupling) or coupling < 0.0:
        raise ScientificValidationError("coupling_rad_per_s must be finite and >= 0")
    if not math.isfinite(rtol) or rtol <= 0.0:
        raise ScientificValidationError("relative_tolerance must be finite and > 0")
    if not math.isfinite(atol) or atol <= 0.0:
        raise ScientificValidationError("absolute_tolerance must be finite and > 0")
    if math.isnan(maximum_step) or maximum_step <= 0.0:
        raise ScientificValidationError("maximum_step_s must be positive")
    if method not in ("RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"):
        raise ScientificValidationError(f"unsupported integration method: {method!r}")

    population_size = frequencies.size

    def derivative(_time_s: float, state: FloatArray) -> FloatArray:
        collective = np.mean(np.exp(1j * state))
        coupling_term = coupling * np.imag(collective * np.exp(-1j * state))
        return frequencies + coupling_term

    solution = solve_ivp(
        derivative,
        (float(times[0]), float(times[-1])),
        phases,
        method=method,
        t_eval=times,
        rtol=rtol,
        atol=atol,
        max_step=maximum_step,
    )
    if not solution.success:
        raise IntegrationError(f"Kuramoto integration failed: {solution.message}")
    if solution.y.shape != (population_size, times.size):
        raise IntegrationError("Kuramoto integration returned an incomplete trajectory")

    phase_trajectory = np.asarray(solution.y.T, dtype=np.float64).copy()
    collective_trajectory = np.mean(np.exp(1j * phase_trajectory), axis=1)
    return KuramotoResult(
        time_s=solution.t.copy(),
        phases_rad=phase_trajectory,
        order_parameter=np.asarray(collective_trajectory, dtype=np.complex128),
        coherence=np.asarray(np.abs(collective_trajectory), dtype=np.float64),
        coupling_rad_per_s=coupling,
        function_evaluations=int(solution.nfev),
        method=method,
        relative_tolerance=rtol,
        absolute_tolerance=atol,
    )
