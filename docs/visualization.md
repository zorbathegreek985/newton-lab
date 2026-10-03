# Scientific visualization and experiment reports

The optional `newton_lab.visualization` module turns existing Newton Lab
results and Phase 14 analysis records into headless figures and deterministic
Markdown. It never reruns a solver. Install Matplotlib with the project and
development dependencies:

```powershell
.venv\Scripts\python.exe -m pip install -e ".[dev,visualization]"
```

Figures are built with Matplotlib's object API and the noninteractive Agg
backend. They can be saved as PNG or SVG. Existing files are protected by
default; pass `overwrite=True` to replace one, and pass
`create_parent_dirs=True` only when the destination folders should be created.
Call `artifact.close()` when a figure is no longer needed.

## Plot types

`plot_ode_trajectory` plots selected oscillator or nonlinear-pendulum state
components against the actual returned time samples. Displacement, velocity,
angle, and angular velocity carry their declared SI or angular units. The
figure does not interpolate between samples, so a sampled maximum can miss an
extremum between output times.

`plot_metric_sweep` accepts a Phase 14 `ParameterSweepAnalysis` and preserves
run order, duplicates, recorded parameter values, and failed or unavailable
metric status. It uses disconnected points; it does not sort, interpolate, or
connect values. The separate `plot_algebraic_sweep` accepts a nonlinear
algebraic `ExperimentRecord` and only plots solution components, residual
components, or `residual_norm`. Candidate points are left disconnected because
root-branch identity is not established by a local solve. Neither convergence
nor a small residual proves uniqueness, completeness, or physical validity.

`plot_bvp_profile` uses the returned spatial mesh and nodal solution without
resampling. An analytical overlay is available only for a matching passed
linear heat-conduction reference comparison whose error metrics are verified
against the plotted profile. The comparison is evaluated at each coordinate
on that same returned mesh.

`plot_reference_error_sweep` plots the recorded scalar or BVP profile absolute
reference error. Reference error is separate from solver convergence and
numerical acceptance. `plot_collection_summary` displays Phase 14 run counts
by family and execution status. `plot_metric_distribution` displays finite
metric observations from successful runs alongside reason counts for excluded
observations; its histogram is descriptive and is not a fitted probability
distribution.

## Example

```python
from newton_lab.experiment_analysis import analyze_parameter_sweep
from newton_lab.experiments import make_oscillator_experiment, run_experiment
from newton_lab.visualization import plot_metric_sweep

record = run_experiment(
    make_oscillator_experiment(
        sweep_values=(0.0, 0.1, 0.2),
        metrics=("final_displacement_m",),
    )
)
analysis = analyze_parameter_sweep(
    record,
    parameter_name="damping_coefficient_kg_per_s",
    metric_id="final_displacement_m",
)
artifact = plot_metric_sweep(
    analysis, save_path="artifacts/damping.png", create_parent_dirs=True
)
print(artifact.saved_path, artifact.source_ids)
artifact.close()
```

The plot's metadata identifies experiment and run IDs, while `notes` records
important interpretation limits. A plot is a view of a particular set of
stored outputs, not an independent validation of the model or solver.

## Markdown reports

`build_experiment_report` accepts experiment records, optional selected metric
IDs, and saved plot artifacts. It uses Phase 14 summaries and projections and
includes run outcomes, solver and numerical acceptance fields, parameters,
failures, requested metrics, reference comparisons, and provenance limits.
It keeps family-specific record data intact and does not combine ODE samples,
algebraic candidates, and BVP profiles into a common result type.

```python
from newton_lab.visualization import build_experiment_report, write_experiment_report

report = build_experiment_report(
    (record,),
    selected_metrics=("final_displacement_m",),
    plots=(artifact,),
)
write_experiment_report(report, "reports/damping.md", create_parent_dirs=True)
```

Report text is deterministic for the same records and options. Plot references
must point to existing saved artifacts. Markdown output also refuses to
overwrite an existing file unless `overwrite=True` is explicit. Reports
summarize recorded evidence; they do not prove causal relationships,
uniqueness, physical validity, or predictive performance.

## Numerical limitations

Visual appearance depends on chosen units, axis ranges, sampling, solver
tolerances, and rendering. A smooth line between ODE samples is a display
connection, not an additional solver evaluation. A BVP line joins returned
nodes, and a histogram depends on its bins and sample size. These plots expose
the supplied observations without estimating unobserved values. Solver
tolerances constrain numerical algorithms but do not quantify model error or
experimental uncertainty. Keep the original records and inspect their
diagnostics when interpreting any figure.
