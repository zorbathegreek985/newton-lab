"""Tests for the Phase 48 pinned-path modal validation study."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest

from newton_lab import pinned_path_network
from newton_lab.exceptions import ScientificValidationError
from newton_lab.pinned_path_network import (
    GRID_SIZES,
    MODE_NUMBERS,
    continuum_mode_decay_rate,
    discrete_mode_decay_rate,
    estimate_modal_decay_rate,
    initial_sine_mode,
    modal_projection,
    path_laplacian_operator,
    run_pinned_path_study,
    simulate_pinned_path_mode,
    validate_mode,
)


@pytest.mark.parametrize("bad_count", [0, -1, 1.5, True])
def test_agent_count_must_be_positive_integer(bad_count: object) -> None:
    with pytest.raises(ScientificValidationError):
        path_laplacian_operator(cast(int, bad_count))


@pytest.mark.parametrize("mode", [0, -1, 1.2, True])
def test_mode_must_be_positive_integer(mode: object) -> None:
    with pytest.raises(ScientificValidationError):
        validate_mode(5, cast(int, mode))


def test_mode_cannot_exceed_agent_count() -> None:
    with pytest.raises(ScientificValidationError, match="between 1 and agent_count"):
        validate_mode(3, 4)


def test_path_operator_has_fixed_boundary_dirichlet_stencil() -> None:
    np.testing.assert_array_equal(
        path_laplacian_operator(3),
        np.array([[32.0, -16.0, 0.0], [-16.0, 32.0, -16.0], [0.0, -16.0, 32.0]]),
    )
    # For one interior node the two fixed leaders contribute both boundary terms.
    assert path_laplacian_operator(1)[0, 0] == 8.0


@pytest.mark.parametrize("count", [3, 7, 15])
def test_discrete_sine_modes_are_operator_eigenvectors(count: int) -> None:
    operator = path_laplacian_operator(count)
    for mode in range(1, count + 1):
        vector = initial_sine_mode(count, mode)
        rate = discrete_mode_decay_rate(count, mode)
        np.testing.assert_allclose(
            operator @ vector, rate * vector, rtol=2e-12, atol=2e-12
        )


@pytest.mark.parametrize("count", [3, 7, 15])
def test_discrete_eigenvalues_match_independent_matrix_eigensystem(count: int) -> None:
    expected = np.linalg.eigvalsh(path_laplacian_operator(count))
    actual = np.array(
        [discrete_mode_decay_rate(count, mode) for mode in range(1, count + 1)]
    )
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)


def test_continuum_rate_is_independent_heat_api_reference() -> None:
    assert continuum_mode_decay_rate(2) == pytest.approx((2.0 * math.pi) ** 2)


def test_projection_and_log_rate_recovery() -> None:
    vector = np.array([1.0, 2.0, 1.0])
    times = np.linspace(0.0, 2.0, 101)
    states = np.exp(-3.5 * times[:, None]) * vector
    np.testing.assert_allclose(modal_projection(states, vector), np.exp(-3.5 * times))
    assert estimate_modal_decay_rate(
        times, modal_projection(states, vector)
    ) == pytest.approx(3.5)


def test_invalid_projection_and_rate_inputs_are_rejected() -> None:
    with pytest.raises(ScientificValidationError):
        modal_projection(np.ones((3, 2)), np.zeros(2))
    with pytest.raises(ScientificValidationError):
        estimate_modal_decay_rate(np.array([0.0, 1.0]), np.array([1.0, 0.5]))


def test_simulation_returns_numerical_and_exact_discrete_trajectories() -> None:
    run, times, numerical, exact = simulate_pinned_path_mode(7, 2, output_points=101)
    assert run.solver_success
    assert numerical.shape == exact.shape == (101, 7)
    assert times.shape == (101,)
    np.testing.assert_allclose(
        numerical[0], initial_sine_mode(7, 2), rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(numerical, exact, rtol=2e-8, atol=2e-10)
    assert run.discrete_rate_per_tau != run.continuum_rate_per_tau


def test_simulation_is_deterministic_for_identical_settings() -> None:
    first = simulate_pinned_path_mode(15, 1, output_points=81)
    second = simulate_pinned_path_mode(15, 1, output_points=81)
    for left, right in zip(first[1:], second[1:], strict=True):
        np.testing.assert_array_equal(left, right)
    assert first[0] == second[0]


def test_solver_failure_is_reported_clearly(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        pinned_path_network,
        "solve_ivp",
        lambda *_args, **_kwargs: SimpleNamespace(
            success=False, y=np.empty((0, 0)), message="synthetic solver failure"
        ),
    )
    with pytest.raises(RuntimeError, match="synthetic solver failure"):
        simulate_pinned_path_mode(7, 1, output_points=11)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"output_points": 1},
        {"relative_tolerance": 0.0},
        {"absolute_tolerance": float("nan")},
    ],
)
def test_simulation_rejects_invalid_sampling_and_tolerances(
    kwargs: dict[str, Any],
) -> None:
    with pytest.raises(ScientificValidationError):
        simulate_pinned_path_mode(7, 1, **kwargs)


def test_sampling_and_tolerance_settings_do_not_change_exact_reference() -> None:
    loose, _t1, _x1, exact_loose = simulate_pinned_path_mode(
        31, 2, output_points=101, relative_tolerance=1e-9, absolute_tolerance=1e-11
    )
    tight, _t2, _x2, exact_tight = simulate_pinned_path_mode(
        31, 2, output_points=201, relative_tolerance=1e-12, absolute_tolerance=1e-14
    )
    assert loose.discrete_rate_per_tau == tight.discrete_rate_per_tau
    assert loose.continuum_rate_per_tau == tight.continuum_rate_per_tau
    np.testing.assert_array_equal(exact_loose[0], exact_tight[0])


def test_study_outputs_schema_and_frozen_convergence_criteria(tmp_path: Path) -> None:
    metadata = run_pinned_path_study(tmp_path)
    assert metadata["row_count_per_condition_results"] == 60
    assert metadata["row_count_convergence_summary"] == len(MODE_NUMBERS)
    assert metadata["prediction_supported_all_modes"] is True
    expected_files = {
        "per_condition_results.csv",
        "convergence_summary.csv",
        "metadata.json",
        "modal_rate_convergence.png",
        "modal_trajectory_comparison.png",
    }
    assert {path.name for path in tmp_path.iterdir()} == expected_files
    with (tmp_path / "per_condition_results.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        conditions = list(csv.DictReader(stream))
    assert len(conditions) == 60
    assert {int(row["agent_count"]) for row in conditions} == set(GRID_SIZES)
    with (tmp_path / "convergence_summary.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        summaries = list(csv.DictReader(stream))
    assert len(summaries) == 3
    for row in summaries:
        assert row["all_ratios_below_one"] == "True"
        assert row["error_decreases_each_refinement"] == "True"
        assert row["order_within_1_7_to_2_3"] == "True"
    assert json.loads((tmp_path / "metadata.json").read_text(encoding="utf-8"))[
        "prediction_supported_all_modes"
    ]
