from types import SimpleNamespace

import numpy as np
import pytest

import newton_lab.kuramoto as kuramoto
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.kuramoto import order_parameter, simulate_kuramoto


def test_zero_coupling_matches_independent_phase_rotation_and_preserves_inputs() -> (
    None
):
    frequencies = np.array([-0.7, 0.0, 0.4, 1.2], dtype=np.float64)
    phases = np.array([-1.0, -0.2, 0.5, 2.0], dtype=np.float64)
    times = np.linspace(0.0, 3.0, 31)
    original_frequencies = frequencies.copy()
    original_phases = phases.copy()
    result = simulate_kuramoto(frequencies, phases, 0.0, times)

    expected = phases[None, :] + times[:, None] * frequencies[None, :]
    np.testing.assert_allclose(result.phases_rad, expected, rtol=1e-9, atol=1e-10)
    np.testing.assert_allclose(result.coherence, np.abs(result.order_parameter))
    assert result.phases_rad.shape == (times.size, frequencies.size)
    assert result.coherence.shape == times.shape
    assert result.coherence[0] == pytest.approx(abs(order_parameter(phases)))
    np.testing.assert_array_equal(frequencies, original_frequencies)
    np.testing.assert_array_equal(phases, original_phases)
    np.testing.assert_array_equal(times, np.linspace(0.0, 3.0, 31))


def test_identical_in_phase_oscillators_remain_fully_coherent() -> None:
    phases = np.full(12, 0.3, dtype=np.float64)
    frequencies = np.full(12, 1.1, dtype=np.float64)
    result = simulate_kuramoto(frequencies, phases, 4.0, np.linspace(0.0, 5.0, 51))

    np.testing.assert_allclose(result.coherence, 1.0, atol=1e-12)
    np.testing.assert_allclose(
        result.phases_rad,
        phases[None, :] + result.time_s[:, None] * 1.1,
        atol=1e-8,
    )


def test_phase_shift_and_oscillator_permutation_leave_coherence_unchanged() -> None:
    frequencies = np.array([-0.9, -0.2, 0.1, 0.8])
    phases = np.array([-2.0, -0.4, 0.7, 2.4])
    times = np.linspace(0.0, 2.0, 21)
    baseline = simulate_kuramoto(frequencies, phases, 1.5, times)
    shifted = simulate_kuramoto(frequencies, phases + 1.7, 1.5, times)
    permutation = np.array([2, 0, 3, 1])
    permuted = simulate_kuramoto(
        frequencies[permutation], phases[permutation], 1.5, times
    )

    np.testing.assert_allclose(shifted.coherence, baseline.coherence, atol=1e-8)
    np.testing.assert_allclose(permuted.coherence, baseline.coherence, atol=1e-8)


@pytest.mark.parametrize(
    ("frequencies", "phases", "coupling", "times"),
    [
        ([0.0], [0.0], 1.0, [0.0, 1.0]),
        ([0.0, np.nan], [0.0, 1.0], 1.0, [0.0, 1.0]),
        ([0.0, 1.0], [0.0], 1.0, [0.0, 1.0]),
        ([0.0, 1.0], [0.0, 1.0], -0.1, [0.0, 1.0]),
        ([0.0, 1.0], [0.0, 1.0], 1.0, [0.1, 1.0]),
        ([0.0, 1.0], [0.0, 1.0], 1.0, [0.0, 0.0]),
    ],
)
def test_invalid_model_inputs_are_rejected(
    frequencies: list[float],
    phases: list[float],
    coupling: float,
    times: list[float],
) -> None:
    with pytest.raises(ScientificValidationError):
        simulate_kuramoto(
            np.asarray(frequencies),
            np.asarray(phases),
            coupling,
            np.asarray(times),
        )


def test_integration_failure_is_reported(monkeypatch: pytest.MonkeyPatch) -> None:
    def failed_solver(*_args: object, **_kwargs: object) -> SimpleNamespace:
        return SimpleNamespace(success=False, message="deliberate solver failure")

    monkeypatch.setattr(kuramoto, "solve_ivp", failed_solver)
    with pytest.raises(IntegrationError, match="deliberate solver failure"):
        simulate_kuramoto(
            np.array([-1.0, 1.0]),
            np.array([0.0, 1.0]),
            1.0,
            np.linspace(0.0, 1.0, 3),
        )


def test_order_parameter_is_zero_for_equally_spaced_phases() -> None:
    phases = np.linspace(-np.pi, np.pi, 32, endpoint=False)

    assert abs(order_parameter(phases)) == pytest.approx(0.0, abs=1e-15)


def test_unsupported_solver_method_is_rejected() -> None:
    with pytest.raises(
        ScientificValidationError, match="unsupported integration method"
    ):
        simulate_kuramoto(
            np.array([-1.0, 1.0]),
            np.array([0.0, 1.0]),
            1.0,
            np.linspace(0.0, 1.0, 3),
            method="not-a-solver",  # type: ignore[arg-type]
        )
