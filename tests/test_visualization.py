"""Headless scientific visualizations and experiment reporting."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from newton_lab.algebraic import AlgebraicSolverConfiguration
from newton_lab.dynamics import simulate_damped_oscillator
from newton_lab.exceptions import ScientificValidationError
from newton_lab.experiment_analysis import (
    analyze_parameter_sweep,
    compare_metric_observations,
    extract_metric_observations,
    project_experiment_runs,
    summarize_experiments,
)
from newton_lab.experiments import (
    AnalyticalReferenceRequest,
    BVPReferenceRequest,
    make_heat_conduction_experiment,
    make_nonlinear_spring_experiment,
    make_oscillator_experiment,
    run_experiment,
)
from newton_lab.visualization import (
    build_experiment_report,
    plot_algebraic_sweep,
    plot_bvp_profile,
    plot_collection_summary,
    plot_metric_distribution,
    plot_metric_sweep,
    plot_ode_trajectory,
    plot_reference_error_sweep,
    write_experiment_report,
)


def test_ode_trajectory_uses_returned_samples_and_keeps_units() -> None:
    result = simulate_damped_oscillator(1.0, 0.2, 4.0, 0.1, 0.0, 2.0, num_points=31)
    time_before = result.time_s.copy()
    displacement_before = result.displacement_m.copy()

    artifact = plot_ode_trajectory(result, components=("displacement_m",))
    axis = artifact.figure.axes[0]

    np.testing.assert_array_equal(axis.lines[0].get_xdata(), result.time_s)
    np.testing.assert_array_equal(axis.lines[0].get_ydata(), result.displacement_m)
    assert axis.get_xlabel() == "Time (s)"
    assert axis.get_ylabel() == "displacement_m (m)"
    assert artifact.source_ids == ("damped_harmonic_oscillator",)
    np.testing.assert_array_equal(result.time_s, time_before)
    np.testing.assert_array_equal(result.displacement_m, displacement_before)
    artifact.close()


def test_sweep_plot_preserves_repeats_failures_and_disconnected_points() -> None:
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="visual_sweep",
            sweep_values=(0.0, -0.1, 0.0),
            metrics=("final_displacement_m",),
        )
    )
    analysis = analyze_parameter_sweep(
        record,
        parameter_name="damping_coefficient_kg_per_s",
        metric_id="final_displacement_m",
    )

    artifact = plot_metric_sweep(analysis)

    assert artifact.source_ids == tuple(f"visual_sweep#{i}" for i in range(3))
    assert len(artifact.figure.axes[0].collections) == 2
    assert len(artifact.figure.axes[0].lines) == 0
    assert "Failed or unavailable" in " ".join(artifact.notes)
    assert [point.parameter_value for point in analysis.points] == [0.0, -0.1, 0.0]
    artifact.close()


def test_algebraic_candidates_are_not_connected_or_called_unique() -> None:
    record = run_experiment(
        make_nonlinear_spring_experiment(
            experiment_id="algebraic_plot",
            sweep_values=(0.8, 1.001, 2.0),
            metrics=("solution:displacement",),
            solver_configuration=AlgebraicSolverConfiguration(),
        )
    )

    artifact = plot_algebraic_sweep(
        record,
        parameter_name="applied_force_n",
        metric_id="solution:displacement",
    )

    assert len(artifact.figure.axes[0].lines) == 0
    assert "root-branch identity" in " ".join(artifact.notes)
    assert "uniqueness" in " ".join(artifact.notes)
    assert artifact.source_ids == tuple(f"algebraic_plot#{i}" for i in range(3))
    artifact.close()


def test_bvp_profile_uses_spatial_mesh_and_only_matching_passed_reference() -> None:
    record = run_experiment(
        make_heat_conduction_experiment(
            experiment_id="heat_profile",
            left_temperature_k=300.0,
            right_temperature_k=400.0,
            sweep_values=None,
            bvp_reference=BVPReferenceRequest(maximum_absolute_tolerance=1e-6),
        )
    )
    projection = project_experiment_runs((record,))[0]
    diagnostics = projection.bvp_diagnostics
    assert diagnostics is not None
    lower, upper = diagnostics.domain
    midpoint = lower + 0.3 * (upper - lower)
    nonuniform = diagnostics.model_copy(
        update={
            "mesh": (lower, midpoint, upper),
            "solution": (
                (300.0, 330.0, 400.0),
                ((400.0 - 300.0) / (upper - lower),) * 3,
            ),
        }
    )
    projection = projection.model_copy(update={"bvp_diagnostics": nonuniform})

    artifact = plot_bvp_profile(
        projection, components=("temperature", "temperature_gradient")
    )

    axis, gradient_axis = artifact.figure.axes
    np.testing.assert_array_equal(axis.lines[0].get_xdata(), nonuniform.mesh)
    np.testing.assert_array_equal(axis.lines[0].get_ydata(), nonuniform.solution[0])
    np.testing.assert_array_equal(
        gradient_axis.lines[0].get_ydata(), nonuniform.solution[1]
    )
    assert nonuniform.mesh[1] != pytest.approx(
        (nonuniform.mesh[0] + nonuniform.mesh[-1]) / 2
    )
    assert len(axis.lines) == 1
    assert "without smoothing" in " ".join(artifact.notes)

    matching_projection = project_experiment_runs((record,))[0]
    reference_artifact = plot_bvp_profile(
        matching_projection,
        components=("temperature",),
        include_reference=True,
    )
    assert len(reference_artifact.figure.axes[0].lines) == 2

    passed_reference = record.runs[0].bvp_reference_comparison
    assert passed_reference is not None
    bad_reference = passed_reference.model_copy(
        update={"maximum_absolute_error": 100.0}
    )
    with pytest.raises(ScientificValidationError, match="metrics do not match"):
        plot_bvp_profile(
            projection,
            components=("temperature",),
            include_reference=True,
            reference_comparison=bad_reference,
        )
    artifact.close()
    reference_artifact.close()


def test_collection_and_metric_views_keep_exclusions_visible() -> None:
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="collection_metrics",
            sweep_values=(0.0, -0.2),
            metrics=("final_displacement_m",),
        )
    )
    summary_artifact = plot_collection_summary(summarize_experiments((record,)))
    observations = extract_metric_observations((record,), "final_displacement_m")
    comparison = compare_metric_observations(observations)
    distribution = plot_metric_distribution(comparison)

    assert summary_artifact.figure.axes[0].patches
    assert distribution.figure.axes[0].patches
    assert distribution.figure.axes[1].patches
    assert "does not assert a probability distribution" in " ".join(distribution.notes)
    summary_artifact.close()
    distribution.close()


def test_reference_sweep_plots_recorded_error_and_not_acceptance() -> None:
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="reference_sweep",
            sweep_values=(0.0, 0.1),
            metrics=("final_displacement_m",),
            analytical_reference=AnalyticalReferenceRequest(
                metric_id="final_displacement_m", absolute_tolerance=1e-6
            ),
        )
    )
    analysis = analyze_parameter_sweep(
        record,
        parameter_name="damping_coefficient_kg_per_s",
        metric_id="final_displacement_m",
    )
    artifact = plot_reference_error_sweep(analysis)

    assert artifact.figure.axes[0].get_ylabel().startswith("reference_absolute_error")
    assert "distinct from solver convergence" in " ".join(artifact.notes)
    artifact.close()


def test_plot_rejects_invalid_components_nonfinite_data_and_output_suffix(
    tmp_path: Path,
) -> None:
    result = simulate_damped_oscillator(1.0, 0.1, 2.0, 0.1, 0.0, 1.0)
    with pytest.raises(ScientificValidationError, match="unknown ODE components"):
        plot_ode_trajectory(result, components=("acceleration_m_per_s2",))

    nonfinite = replace(
        result,
        displacement_m=np.array(result.displacement_m, copy=True),
    )
    nonfinite.displacement_m[2] = np.nan
    with pytest.raises(ScientificValidationError, match="finite"):
        plot_ode_trajectory(nonfinite, components=("displacement_m",))

    with pytest.raises(ScientificValidationError, match=".png or .svg"):
        plot_ode_trajectory(
            result,
            components=("displacement_m",),
            save_path=tmp_path / "trajectory.pdf",
        )


def test_plot_saving_refuses_overwrite_and_report_is_reproducible(
    tmp_path: Path,
) -> None:
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="report_record",
            sweep_values=None,
            metrics=("final_displacement_m",),
        )
    )
    plot_path = tmp_path / "trajectory.png"
    artifact = plot_ode_trajectory(
        simulate_damped_oscillator(1.0, 0.1, 2.0, 0.1, 0.0, 1.0),
        components=("displacement_m",),
        save_path=plot_path,
    )
    assert plot_path.is_file()
    with pytest.raises(FileExistsError):
        plot_ode_trajectory(
            simulate_damped_oscillator(1.0, 0.1, 2.0, 0.1, 0.0, 1.0),
            components=("displacement_m",),
            save_path=plot_path,
        )

    report = build_experiment_report(
        (record,), selected_metrics=("final_displacement_m",), plots=(artifact,)
    )
    markdown = report.render_markdown(relative_to=tmp_path)
    assert markdown.count("# Newton Lab experiment report") == 1
    assert "physical validity, uniqueness" in markdown
    assert "trajectory.png" in markdown
    report_path = write_experiment_report(report, tmp_path / "report.md")
    assert report_path.read_text(encoding="utf-8") == markdown
    with pytest.raises(FileExistsError):
        write_experiment_report(report, report_path)
    artifact.close()


def test_report_includes_failures_provenance_and_is_deterministic() -> None:
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="failed_report",
            sweep_values=(0.0, -0.1),
            metrics=("final_displacement_m",),
        )
    )
    before = record.model_dump(mode="python")

    first = build_experiment_report(
        (record,), selected_metrics=("final_displacement_m",)
    )
    second = build_experiment_report(
        (record,), selected_metrics=("final_displacement_m",)
    )

    assert first.body == second.body
    assert "| 1 | failed |" in first.body
    assert "Parameters and units" in first.body
    assert "Runs without reproducibility metadata: 0 of 2" in first.body
    assert record.model_dump(mode="python") == before


def test_report_rejects_unsaved_plot_artifacts() -> None:
    record = run_experiment(
        make_oscillator_experiment(sweep_values=None, metrics=("final_displacement_m",))
    )
    artifact = plot_ode_trajectory(
        simulate_damped_oscillator(1.0, 0.1, 2.0, 0.1, 0.0, 1.0),
        components=("displacement_m",),
    )

    with pytest.raises(ScientificValidationError, match="only saved plot"):
        build_experiment_report((record,), plots=(artifact,))
    artifact.close()


def test_reports_distinguish_algebraic_and_ode_families() -> None:
    ode = run_experiment(
        make_oscillator_experiment(sweep_values=None, metrics=("final_displacement_m",))
    )
    algebraic = run_experiment(
        make_nonlinear_spring_experiment(
            sweep_values=None,
            metrics=("solution:displacement",),
            solver_configuration=AlgebraicSolverConfiguration(),
        )
    )

    report = build_experiment_report((ode, algebraic))

    assert "Problem family: `ode`" in report.body
    assert "Problem family: `algebraic`" in report.body
    assert "family-specific diagnostics" in report.body
