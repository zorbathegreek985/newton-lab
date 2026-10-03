# Phases 16–17: evidence-aware and trajectory-aware discovery

`newton_lab.discovery.run_discovery` is the first end-to-end discovery
workflow. It accepts existing `ExperimentRecord` objects, reuses Phase 14
metric extraction and sweep analysis, compares the records' associated
equations with the existing knowledge base, and returns typed observations,
structure matches, hypotheses, and a deterministic Markdown report. It does
not run experiments, change the registry, or need a CLI.

## Run the example

Install the project's environment and run:

```powershell
.venv\Scripts\python.exe -c "from newton_lab.discovery import run_discovery; from newton_lab.experiments import make_oscillator_experiment, run_experiment; record = run_experiment(make_oscillator_experiment(experiment_id='damping_discovery', sweep_values=(0.0, 0.1, 0.2), metrics=('maximum_absolute_velocity_m_per_s',))); result = run_discovery((record,)); print(result.report_markdown)"
```

Or use the API in a script:

```python
from newton_lab.discovery import run_discovery
from newton_lab.experiments import make_oscillator_experiment, run_experiment

record = run_experiment(
    make_oscillator_experiment(
        experiment_id="damping_discovery",
        sweep_values=(0.0, 0.1, 0.2),
        metrics=("maximum_absolute_velocity_m_per_s",),
        retain_ode_trajectory=True,
    )
)
result = run_discovery((record,))
print(result.observations)
print(result.structure_matches)
print(result.hypotheses)
print(result.report_markdown)
```

Pass `knowledge_base=...` to compare against another existing workspace and
`selected_metrics=(...)` to restrict the requested metrics. If omitted, the
workflow uses the built-in example knowledge base and every metric requested
by each record. A metric not requested by a record is not fabricated for it.

## What is characterized

Each requested Phase 14 metric becomes a `BehavioralObservation` carrying the
experiment and run identity, exact model ID, parameters and units, metric
definition, equation association method, solver success, numerical acceptance,
and sampling limitations. Successful, finite records are labeled as computed.
Failed, missing, non-finite, or invalid metrics remain unavailable with an
explicit reason; they are never replaced by zero.

For a reliable parameter sweep, the first implementation also records a
derived endpoint change: last metric value minus first metric value in stored
run order. It uses Phase 14's sweep analysis, does not sort or interpolate, and
leaves the descriptor unavailable when ordering or either endpoint is
unreliable. Repeated parameter values and failures remain in provenance and
limitations. This endpoint difference does not establish monotonic response or
causation.

## Optional sampled ODE trajectories

Set `retain_ode_trajectory=True` in an oscillator or pendulum specification to
copy the solver's returned time grid and state arrays into each successful
`ExperimentRun`. The default is `False` to avoid storing those arrays for
experiments that only need scalar metrics. Each snapshot contains the model ID,
state names and units, sample times, sampled state values, and solver method;
the enclosing run retains its parameters, requested output count, solver
tolerances, run status, and equation/source provenance. Algebraic and BVP runs
do not acquire ODE time-series data.

For retained oscillator and pendulum samples, discovery emits three temporal
observations per run:

- `observed_period`: median spacing between linearly interpolated upward
  crossings of the model's zero-equilibrium displacement/angle. At least three
  crossings (two cycles) and 12 returned samples are required. This is a
  sampled estimate, not proof of exact periodicity.
- `amplitude_change`: last-window maximum absolute displacement/angle minus
  the first-window maximum, using equal-count windows spanning 20% of returned
  samples. At least 10 samples are required. This finite-window difference is
  not an exponential fit or damping coefficient.
- `finite_window_return_indicator`: 1 when the last window's maximum absolute
  state is at most 5% of the first window's maximum, and 0 otherwise. The
  threshold is an explicit descriptive convention, not a stability test.

Unavailable estimates have an explicit status and no numeric value. Failed
runs remain represented as failed; sweep observations preserve original run
order, repeated parameter values, and the swept value so temporal descriptors
can be compared without sorting away failures. All temporal results identify
the source run, state, units, sample count, time interval, and equation
association provenance. They describe returned numerical samples only and do
not establish continuous-time behavior, asymptotic stability, physical
validity, or application usefulness. Algebraic solution candidates and BVP
profiles remain separate result types.

BVP and algebraic metrics are retained as their own records; where no equation
ID is declared or exactly matches a registry ID, no structure is guessed from
a model name.

## Structure matching and evidence status

When an experiment has explicit equation record IDs, those IDs are used. If
they are absent, an equation is associated only when `model_id` exactly equals
an existing registry `equation_id`. The matcher reuses
`compare_equations` and the existing relationship proposal API. A match needs
at least two shared, explicitly declared structural properties among
mathematical classification, ODE order, and derivative operators. Equation
expression text, units, and names are not parsed or treated as proof.

Every `DiscoveryStructureMatch` includes the shared properties, supporting
observation IDs, any registry relationship ID, its review state, applicable
conditions, and limitations. Record-level relationships with
`mathematically_derived_relationship` or sourced established-physical evidence
remain established as recorded. A review-graph edge is called established
only when it is accepted and cites source references. The built-in
oscillator/pendulum approximation remains a
`candidate_structural_connection`; its small-angle and coefficient-matching
conditions are shown in the report. Discovery never adds or promotes an edge.

Evidence levels have separate meanings:

- `established_relationship`: qualified relationship record accepted with
  source references.
- `candidate_structural_connection`: metadata-based connection still
  requiring scientific review.
- `proposed_application_hypothesis`: testable cross-domain idea, not validated.
- `validated_application`: reserved for independently documented validation;
  this phase does not emit it. Hypothesis objects enforce `validated_application=False`.

## Hypotheses and limitations

The deliberately narrow application rule emits a quantitative-finance risk
control hypothesis only when a damped harmonic oscillator record has a
reliable damping sweep with at least two distinct parameter values, successful
maximum-absolute-velocity metrics at every run, a negative last-minus-first
endpoint change, exact oscillator equation association, and at least one
structure match. The hypothesis reports that observed difference and asks
whether a damping-like term in a second-order risk-control update can reduce
peak oscillatory adjustment rates. This is a research question, not evidence of
a monotonic relationship. It lists assumptions, similarities, differences,
evidence IDs, and a test plan. The registry contains no finance equation
mapping, and the hypothesis is not a trading signal, financial recommendation,
or claim of profitability.

An appropriate future test would require an explicitly defined target model,
prespecified baselines and rejection criteria, chronological out-of-sample
evaluation, transaction costs and market impact, risk measures, and
overfitting safeguards. A successful physical simulation, shared form, or
plausible plot cannot validate that transfer. Numerical acceptance, physical
validity, mathematical validity, and usefulness in another domain remain
separate questions.

Inspect and challenge a result through its source run IDs, parameter values,
metric definition, association method, match feature list, relationship review
status, conditions, hypothesis assumptions, and test plan. Missing or
insufficient evidence appears in the report instead of being filled by a
guess. Report generation reuses the Phase 15 experiment report and does not
rerun experiments or embed large numerical arrays.
