"""Tests for the nonlinear damped pendulum model."""

from types import SimpleNamespace

import numpy as np
import pytest

import newton_lab.pendulum as pendulum
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.pendulum import simulate_damped_pendulum


def _mechanical_energy(
    mass_kg: float,
    length_m: float,
    gravity_m_per_s2: float,
    angle_rad: np.ndarray,
    angular_velocity_rad_per_s: np.ndarray,
) -> np.ndarray:
    return (
        0.5 * mass_kg * length_m** 2 * angular_velocity_rad_per_s** 2
        + mass_kg * gravity_m_per_s2 * length_m * (1.0 - np.cos(angle_rad))
    )


def test_valid_simulation_has_expected_shapes_and_metadata() -> None:
    result = simulate_damped_pendulum(1.0, 1.5, 0.1, 9.81, 0.2, 0.0, 4.0, num_points=81)

    assert result.time_s.shape == (81,)
    assert result.angle_rad.shape == (81,)
    assert result.angular_velocity_rad_per_s.shape == (81,)
    assert result.time_s[0] == 0.0
    assert result.time_s[-1] == 4.0
    assert result.mass_kg == 1.0
    assert result.length_m == 1.5
    assert result.gravity_m_per_s2 == 9.81
    assert result.method == "RK45"
    assert result.function_evaluations > 0


def test_initial_angle_and_angular_velocity_are_preserved() -> None:
    result = simulate_damped_pendulum(0.8, 2.0, 0.2, 9.81, -0.4, 0.3, 1.0)

    assert result.angle_rad[0] == -0.4
    assert result.angular_velocity_rad_per_s[0] == 0.3


def test_zero_damping_energy_is_approximately_conserved() -> None:
    mass, length, gravity = 0.8, 1.4, 9.81
    result = simulate_damped_pendulum(
        mass,
        length,
        0.0,
        gravity,
        1.2,
        -0.4,
        20.0,
        num_points=1001,
        method="DOP853",
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
    )
    energy = _mechanical_energy(
        mass,
        length,
        gravity,
        result.angle_rad,
        result.angular_velocity_rad_per_s,
    )

    # Tolerance allows accumulated local integration error over many cycles.
    assert np.ptp(energy) < energy[0] * 2e-8


def test_positive_damping_does_not_increase_energy_beyond_solver_error() -> None:
    mass, length, gravity = 0.8, 1.4, 9.81
    result = simulate_damped_pendulum(
        mass,
        length,
        0.3,
        gravity,
        1.2,
        -0.4,
        10.0,
        num_points=1001,
        method="DOP853",
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
    )
    energy = _mechanical_energy(
        mass,
        length,
        gravity,
        result.angle_rad,
        result.angular_velocity_rad_per_s,
    )

    # Exact model energy derivative is -b*omega**2; this allowance covers
    # small numerical fluctuations from adaptive integration and sampling.
    assert np.diff(energy).max() < energy[0] * 2e-8
    assert energy[-1] < energy[0]


def test_small_angle_trajectory_agrees_with_linear_analytical_solution() -> None:
    mass, length, gravity = 1.0, 1.2, 9.81
    initial_angle = 1e-3
    result = simulate_damped_pendulum(
        mass,
        length,
        0.0,
        gravity,
        initial_angle,
        0.0,
        10.0,
        num_points=1001,
        method="DOP853",
        relative_tolerance=1e-11,
        absolute_tolerance=1e-13,
    )
    natural_frequency = np.sqrt(gravity / length)
    exact_small_angle = initial_angle * np.cos(natural_frequency * result.time_s)
    exact_small_angle_velocity = (
        -initial_angle * natural_frequency * np.sin(natural_frequency * result.time_s)
    )

    # At 1 mrad, the omitted cubic term produces only a few nanoradians
    # of absolute angle difference over this interval.
    np.testing.assert_allclose(
        result.angle_rad, exact_small_angle, rtol=1e-6, atol=1e-8
    )
    np.testing.assert_allclose(
        result.angular_velocity_rad_per_s,
        exact_small_angle_velocity,
        rtol=1e-6,
        atol=1e-8,
    )


@pytest.mark.parametrize(
    ("parameter", "value"),
    [
        ("mass_kg", 0.0),
        ("mass_kg", -1.0),
        ("length_m", 0.0),
        ("length_m", -1.0),
        ("gravity_m_per_s2", 0.0),
        ("gravity_m_per_s2", -9.81),
        ("damping_coefficient_kg_m2_per_s", -0.1),
    ],
)
def test_rejects_invalid_physical_parameters(parameter: str, value: float) -> None:
    parameters: dict[str, float] = {
        "mass_kg": 1.0,
        "length_m": 1.0,
        "damping_coefficient_kg_m2_per_s": 0.1,
        "gravity_m_per_s2": 9.81,
        "initial_angle_rad": 0.0,
        "initial_angular_velocity_rad_per_s": 0.0,
        "duration_s": 1.0,
    }
    parameters[parameter] = value

    with pytest.raises(ScientificValidationError):
        simulate_damped_pendulum(
            mass_kg=parameters["mass_kg"],
            length_m=parameters["length_m"],
            damping_coefficient_kg_m2_per_s=parameters[
                "damping_coefficient_kg_m2_per_s"
            ],
            gravity_m_per_s2=parameters["gravity_m_per_s2"],
            initial_angle_rad=parameters["initial_angle_rad"],
            initial_angular_velocity_rad_per_s=parameters[
                "initial_angular_velocity_rad_per_s"
            ],
            duration_s=parameters["duration_s"],
        )


@pytest.mark.parametrize(
    ("parameter", "value"),
    [
        ("mass_kg", float("nan")),
        ("length_m", float("inf")),
        ("damping_coefficient_kg_m2_per_s", float("nan")),
        ("gravity_m_per_s2", float("inf")),
        ("initial_angle_rad", float("nan")),
        ("initial_angular_velocity_rad_per_s", float("inf")),
        ("duration_s", 0.0),
        ("duration_s", -1.0),
        ("duration_s", float("inf")),
        ("num_points", 1),
        ("num_points", 1.5),
        ("relative_tolerance", 0.0),
        ("absolute_tolerance", -1e-9),
    ],
)
def test_rejects_nonfinite_conditions_and_invalid_settings(
    parameter: str, value: float
) -> None:
    parameters: dict[str, float] = {
        "mass_kg": 1.0,
        "length_m": 1.0,
        "damping_coefficient_kg_m2_per_s": 0.1,
        "gravity_m_per_s2": 9.81,
        "initial_angle_rad": 0.0,
        "initial_angular_velocity_rad_per_s": 0.0,
        "duration_s": 1.0,
    }
    parameters[parameter] = value

    with pytest.raises(ScientificValidationError):
        simulate_damped_pendulum(**parameters)  # type: ignore[arg-type]


def test_reports_solver_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    def failed_solver(*_args: object, **_kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(success=False, message="forced test failure")

    monkeypatch.setattr(pendulum, "solve_ivp", failed_solver)

    with pytest.raises(IntegrationError, match="forced test failure"):
        simulate_damped_pendulum(1.0, 1.0, 0.1, 9.81, 0.5, 0.0, 1.0)
