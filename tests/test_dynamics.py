"""Tests for the damped harmonic oscillator model."""

from types import SimpleNamespace

import numpy as np
import pytest

import newton_lab.dynamics as dynamics
from newton_lab.dynamics import simulate_damped_oscillator
from newton_lab.exceptions import IntegrationError, ScientificValidationError


def test_valid_simulation_has_requested_shapes_and_metadata() -> None:
    result = simulate_damped_oscillator(1.0, 0.2, 4.0, 0.5, 0.0, 3.0, num_points=101)

    assert result.time_s.shape == (101,)
    assert result.displacement_m.shape == (101,)
    assert result.velocity_m_per_s.shape == (101,)
    assert result.time_s[0] == 0.0
    assert result.time_s[-1] == 3.0
    assert result.method == "RK45"
    assert result.function_evaluations > 0


def test_initial_state_is_preserved() -> None:
    result = simulate_damped_oscillator(2.0, 0.1, 3.0, -0.25, 1.5, 1.0)

    assert result.displacement_m[0] == -0.25
    assert result.velocity_m_per_s[0] == 1.5


def test_zero_damping_agrees_with_undamped_solution() -> None:
    result = simulate_damped_oscillator(
        1.0,
        0.0,
        4.0,
        1.0,
        0.0,
        2.0,
        num_points=201,
        method="DOP853",
        relative_tolerance=1e-11,
        absolute_tolerance=1e-13,
    )
    angular_frequency = 2.0

    np.testing.assert_allclose(
        result.displacement_m,
        np.cos(angular_frequency * result.time_s),
        rtol=2e-10,
        atol=2e-11,
    )
    np.testing.assert_allclose(
        result.velocity_m_per_s,
        -angular_frequency * np.sin(angular_frequency * result.time_s),
        rtol=2e-10,
        atol=2e-11,
    )


def test_positive_damping_reduces_late_displacement() -> None:
    result = simulate_damped_oscillator(1.0, 0.5, 4.0, 1.0, 0.0, 5.0)

    assert abs(result.displacement_m[-1]) < abs(result.displacement_m[0])


@pytest.mark.parametrize(
    ("parameter", "value"),
    [
        ("mass_kg", 0.0),
        ("mass_kg", -1.0),
        ("stiffness_n_per_m", 0.0),
        ("stiffness_n_per_m", -1.0),
        ("damping_coefficient_kg_per_s", -0.1),
    ],
)
def test_rejects_invalid_physical_parameters(parameter: str, value: float) -> None:
    parameters = {
        "mass_kg": 1.0,
        "damping_coefficient_kg_per_s": 0.1,
        "stiffness_n_per_m": 1.0,
        "initial_displacement_m": 0.0,
        "initial_velocity_m_per_s": 0.0,
        "duration_s": 1.0,
    }
    parameters[parameter] = value

    with pytest.raises(ScientificValidationError):
        simulate_damped_oscillator(
            mass_kg=parameters["mass_kg"],
            damping_coefficient_kg_per_s=parameters["damping_coefficient_kg_per_s"],
            stiffness_n_per_m=parameters["stiffness_n_per_m"],
            initial_displacement_m=parameters["initial_displacement_m"],
            initial_velocity_m_per_s=parameters["initial_velocity_m_per_s"],
            duration_s=parameters["duration_s"],
        )


@pytest.mark.parametrize(
    ("parameter", "value"),
    [
        ("initial_displacement_m", float("nan")),
        ("initial_velocity_m_per_s", float("inf")),
        ("duration_s", 0.0),
        ("duration_s", -1.0),
        ("duration_s", float("inf")),
        ("num_points", 1),
        ("num_points", 1.5),
        ("relative_tolerance", 0.0),
        ("absolute_tolerance", -1e-8),
    ],
)
def test_rejects_invalid_initial_conditions_and_settings(
    parameter: str, value: float
) -> None:
    parameters: dict[str, float] = {
        "mass_kg": 1.0,
        "damping_coefficient_kg_per_s": 0.1,
        "stiffness_n_per_m": 1.0,
        "initial_displacement_m": 0.0,
        "initial_velocity_m_per_s": 0.0,
        "duration_s": 1.0,
    }
    parameters[parameter] = value

    with pytest.raises(ScientificValidationError):
        simulate_damped_oscillator(**parameters)  # type: ignore[arg-type]


def test_reports_solver_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    def failed_solver(*_args: object, **_kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(success=False, message="forced test failure")

    monkeypatch.setattr(dynamics, "solve_ivp", failed_solver)

    with pytest.raises(IntegrationError, match="forced test failure"):
        simulate_damped_oscillator(1.0, 0.1, 1.0, 1.0, 0.0, 1.0)
