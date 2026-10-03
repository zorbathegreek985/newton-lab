# Experiment analysis

`newton_lab.experiment_analysis` provides read-only summaries and comparisons
over existing `ExperimentRecord` objects. It does not replace experiment
records, alter solver outcomes, or treat ODE, algebraic, and boundary-value
results as one numerical result type.

## Collection and run summaries

`summarize_experiments(records)` counts records and runs in the supplied
collection. Execution success means exactly `ExperimentRunStatus.SUCCEEDED`;
it is not inferred from metric values. Solver success and numerical acceptance
are counted separately. `None` is counted as unavailable, not false. Family
counts cover ODE, algebraic, boundary-value, or unknown model IDs. Parameter
configuration counts use the exact ordered parameter names, values, and units.

Duplicate experiment IDs are listed with their input record indices and are
counted as separate records. The summary never merges duplicates. The current
record schema has succeeded and failed run statuses, but no third invalid run
status; malformed metrics are classified during metric extraction instead.
Counts by family and model ID count runs (including failed runs), while
parameter configuration counts count exact run configurations.

`project_experiment_runs(records)` creates a common run view. It carries
execution and acceptance states, parameters, reproducibility data, reference
comparisons, and the original typed algebraic or BVP diagnostics. An absent
family-specific diagnostic remains `None`.

```python
from newton_lab.experiment_analysis import (
    analyze_parameter_sweep,
    extract_metric_observations,
    summarize_experiments,
    summarize_metric,
    compare_metric_observations,
)
from newton_lab.experiments import make_oscillator_experiment, run_experiment

record = run_experiment(
    make_oscillator_experiment(
        sweep_values=(0.0, 0.1, 0.2),
        metrics=("final_displacement_m",),
    )
)
summary = summarize_experiments((record,))
observations = extract_metric_observations((record,), "final_displacement_m")
comparison = compare_metric_observations(observations)
statistics = summarize_metric(comparison)
sweep = analyze_parameter_sweep(
    record,
    parameter_name="damping_coefficient_kg_per_s",
    metric_id="final_displacement_m",
)
```

## Metric observations and compatibility

`extract_metric_observations(records, metric_id)` preserves caller record order
and run order. Each observation retains record index and experiment ID, model
and problem IDs, parameters, metric name, unit, adapter-established definition
when known, execution status, solver/acceptance states, failure details,
reproducibility metadata, and applicable family diagnostics and reference
comparisons.

The observation has separate execution and metric states. A failed run remains
an execution failure even if the adapter retained a finite candidate metric.
Missing metrics remain missing, zero remains a measured zero, non-finite values
are marked `non_finite` and not copied as numeric values, and malformed or
duplicate metric entries are marked invalid. Statistics include only finite
metrics from successfully executed runs.

`compare_metric_observations` marks a comparison `direct` only when all
observations share the same adapter-defined model, problem, metric, definition,
and known unit. Matching names or units across different models do not satisfy
this rule. When the records do not establish compatibility, a caller can
provide a `MetricCompatibilityDeclaration` that maps each exact
`(model_id, problem_id, metric_id)` source to its source definition and a
shared definition. The declaration must use one identical unit throughout;
this layer does not convert units. Such a result is `conditional`, and retains
the caller's assumptions. Incomplete declarations, unknown units, and
conflicting definitions produce `not_comparable` with an explanation.

Neither a direct nor conditional comparison declares a smaller value better.
Cross-family comparisons depend on their definitions and stated physical
assumptions; matching units alone are insufficient.

## Sweeps and descriptive statistics

`analyze_parameter_sweep(record, parameter_name, metric_id)` keeps the recorded
run order, reports repeated parameter values, and includes failed or missing
points. It does not sort, interpolate, or fill values. `ordering_reliable` is
false when the record lacks a matching sweep declaration or its stored run
metadata does not match that declaration; observed values and outcomes are
still returned for inspection.

`summarize_metric(comparison)` reports valid and excluded counts, including
execution failures, missing metrics, non-finite values, and invalid metrics.
Exclusion reason counts can overlap when, for example, a failed run also lacks
the selected metric; `excluded_count` counts excluded observations once.
Compatible finite observations receive count, mean, minimum, maximum, and
population standard deviation (`ddof=0`). One valid value has a standard
deviation of zero. No valid values produce `no_valid_observations` and no
statistics. Incompatible observations produce `not_comparable` and no
statistics. These descriptive values establish neither causality nor
predictive power and do not validate the underlying physical model.

Exclusion reason counts are separate indicators and may overlap. For example,
a failed run with no selected metric increments both execution-failure and
missing-metric counts, but increments `excluded_count` only once.

Analysis does not alter source records or their numerical data. It does not
establish solver quality, physical validity, root uniqueness, or superiority
of one solver family over another.
