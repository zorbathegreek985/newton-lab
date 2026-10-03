"""Tests for the deterministic Phase 42 synthetic benchmark."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.synthetic_timeseries_benchmark import (
    SignalKind,
    generate_signal,
    run_benchmark,
    write_benchmark,
)


def test_generation_is_reproducible_and_keeps_truth_separate() -> None:
    first = generate_signal("mean_shift", seed=7, parameter_value=1.5)
    second = generate_signal("mean_shift", seed=7, parameter_value=1.5)

    assert np.array_equal(first.values, second.values)
    assert np.array_equal(first.truth, second.truth)
    assert first.values.shape == first.truth.shape == (512,)
    assert np.all(first.truth[:256] == 0.0)
    assert np.all(first.truth[256:] == 1.5)


@pytest.mark.parametrize(
    ("kind", "parameter"),
    [
        ("stationary_noise", 1.0),
        ("trend_noise", 0.75),
        ("autocorrelated", 0.5),
        ("damped_oscillation", 0.2),
        ("variance_shift", 1.4),
        ("mean_shift", 1.0),
    ],
)
def test_all_signal_families_are_finite(kind: SignalKind, parameter: float) -> None:
    signal = generate_signal(kind, seed=11, parameter_value=parameter)

    assert signal.values.shape == signal.truth.shape == (512,)
    assert np.isfinite(signal.values).all()
    assert np.isfinite(signal.truth).all()


def test_generator_rejects_invalid_settings() -> None:
    with pytest.raises(ScientificValidationError):
        generate_signal("stationary_noise", sample_count=32)
    with pytest.raises(ScientificValidationError):
        generate_signal("stationary_noise", seed=-1)
    with pytest.raises(ScientificValidationError):
        generate_signal("mean_shift", parameter_value=0.0)
    with pytest.raises(ScientificValidationError):
        generate_signal("mean_shift", parameter_value=True)


def test_benchmark_records_sensitivity_grid_and_intended_detections() -> None:
    rows = run_benchmark()

    assert len(rows) == 48
    assert len({row.seed for row in rows}) == 48
    assert all(row.kind != "nse" for row in rows)
    for kind in ("autocorrelated", "mean_shift", "variance_shift"):
        assert any(row.kind == kind for row in rows)


def test_writer_creates_machine_readable_results_and_figures(tmp_path: Path) -> None:
    rows = write_benchmark(tmp_path / "phase42")
    output = tmp_path / "phase42"

    assert len(rows) == 48
    assert (output / "results.csv").is_file()
    assert (output / "metadata.json").is_file()
    assert (output / "detection_rates.png").is_file()
    assert (output / "diagnostic_overlap.png").is_file()
    assert (output / "filter_error.png").is_file()
