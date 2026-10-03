# Phase 31: Financial time-series research protocol

## Purpose and scope

This document defines a protocol for a possible future study of the Phase 30
causal exponential filter on a manually supplied, authorized financial dataset.
It does not perform that study. No market data was accessed, downloaded, or
analyzed for Phase 31, and no forecasting or investment performance is claimed.

The initial question is:

> Does a causal diffusion-inspired filter provide measurable value as a
> transformation of financial time-series data, compared with simple causal
> baselines, when evaluated using only information available at each timestamp?

“Value” must be an outcome defined before evaluation. Smooth-looking output is
not sufficient. Keep these questions separate:

1. Does the filter behave as its recurrence and response equations predict?
2. Does it improve a defined signal-processing measure?
3. Does it improve a prespecified forecasting task?
4. Does any forecast improvement have economic value after execution costs and
   risk constraints?
5. Does the result generalize across periods and instruments?

Phase 30 provides synthetic mathematical and signal-processing evidence for the
first two kinds of questions under its stated synthetic conditions. Phase 31
defines how a later financial study could address the remaining questions. It
does not turn the Phase 30 analogy into financial evidence.

## Established properties and protocol decisions

The candidate is the Phase 30 one-sided exponential filter:

\[
y_0=x_0,\qquad y_t=(1-\beta)x_t+\beta y_{t-1},\qquad 0\leq\beta<1.
\]

On an ordered sequence, output through timestamp \(t\) reads observations only
through \(t\). Phase 30 tests this directly by changing the suffix after a fixed
prefix and checking that the outputs on that prefix are unchanged. The standard
linear, time-invariant transfer function and frequency response describe its
steady-state behavior; finite-record initialization affects early samples.
These are mathematical properties, not evidence that a filtered financial series
is more informative.

The actual Phase 30 API is
`exponential_filter(values, beta) -> NDArray[np.float64]`. It accepts a finite,
nonempty, one-dimensional array-like input and returns a new array. Initialization
is `y[0] = x[0]`; it does not take timestamps or implement missing-data or
resampling policy. Its `run_causal_filter_study()` compares deterministic
synthetic signals with Phase 29's whole-interval diffusion transform and a
trailing mean. Phase 29's transform uses an entire finite record and is
noncausal; it must not be used as a real-time feature.

No timestamped split, financial preprocessing, forecasting, or backtesting
component exists in the inspected project. This protocol does not add one.

## Initial study specification: decisions still required

The initial empirical study must be kept to one named market, instrument or
small universe, one observation frequency, and one target. Do not silently fill
these fields with assumed values:

| Study field | Required declaration before analysis | Current status |
|---|---|---|
| Instrument or universe | Instrument identifiers, inclusion/exclusion rules, and point-in-time membership if a universe is used | Unresolved |
| Provider and dataset | Provider, product/dataset identifier, access method, and dataset version | Unresolved |
| Sampling | Bar/event frequency, timestamp meaning, time zone, daylight-saving rule, calendar, and session policy | Unresolved |
| Observation period | Start/end dates and rationale, declared before test inspection | Unresolved |
| Input field | Exact field(s), units, adjustment/vintage status, and availability timestamp | Unresolved |
| Target and task | Signal-processing reference, or an explicitly defined forecast target and horizon | Unresolved |
| Filter | Candidate beta value(s), initialization, and any training/validation selection rule | Unresolved |
| Baselines | No smoothing; causal trailing mean with declared window; persistence only if forecasting | Unresolved |
| Primary and secondary metrics | Named definitions, aggregation, missing-data denominator, and uncertainty method | Unresolved |
| Data-quality exclusions | Duplicate, missing, stale, corrected, and outlier policies | Unresolved |
| Split boundaries | Chronological training, validation, any required gap, and final test dates | Unresolved |

Every unresolved choice that can affect results must be documented before the
evaluation period is inspected. A change after inspection becomes a new
exploratory analysis and cannot replace the originally specified result.

## Data-access permission gate

**Readiness: not yet established.** No provider, permission evidence, dataset,
or data-use terms were supplied or verified in this phase. In particular, no
IEX data was accessed or downloaded. A future study must not proceed to data
retrieval or analysis until an authorized source and intended use are documented.

Before loading real observations, complete this checklist:

- [ ] Name the provider, product/dataset, version, and responsible account or
  license holder.
- [ ] Record evidence that access and the proposed research use are authorized;
  do not infer permission from technical access.
- [ ] Record whether derived results may be retained, shared, or published.
- [ ] Record whether raw data may be stored or redistributed; assume neither
  unless the terms explicitly allow it.
- [ ] Declare instruments/market coverage and any membership or survivorship
  limitations.
- [ ] Declare observation frequency, timestamp convention, time zone, and
  daylight-saving treatment.
- [ ] Declare trading calendar, holidays, sessions, auctions, and overnight
  boundary treatment where applicable.
- [ ] Describe missing observations, duplicate timestamps, stale values, and
  corrections/revisions.
- [ ] Describe corporate-action adjustment for securities or contract-roll
  treatment for futures, where relevant.
- [ ] Describe survivorship and point-in-time constituent issues, where relevant.
- [ ] Record retrieval date, source identifier, file name, and a file-integrity
  digest or equivalent provenance record.
- [ ] Confirm the manually supplied files can be used without automated
  retrieval and can be handled under the provider's terms.

Permission status, data quality, and scientific suitability are separate gates.
Passing one does not establish the others.

## Chronological evaluation design

Use observations in timestamp order. Do not use random train/test splitting as
the default. Before examining evaluation outcomes, declare fixed date boundaries
for:

1. **Training:** fit any model parameters, normalizers, thresholds, or other
   learned transformations. A filter beta may be fixed in advance or selected
   using only training/validation data under the declared rule.
2. **Validation:** choose among the small, declared candidate configurations
   using the prespecified selection metric. Once selected, freeze the choice
   before the final test.
3. **Final test:** one untouched chronological interval used only for the
   prespecified final evaluation. Do not select parameters, metrics, models,
   exclusions, or a favorable reporting window using its results.

The actual dates, boundaries, and any purge/embargo gap remain unresolved until
the data source, frequency, horizon, and observation period are declared. If the
forecast horizon or target construction overlaps a split boundary, exclude or
purge the observations whose target windows cross the boundary. The gap must be
at least sufficient to prevent target overlap under the declared target and
horizon; feature warm-up history is handled separately. Do not discard valid
past observations merely because they precede a test boundary: a causal filter
may carry its state forward from the past, provided that the state is computed
only from observations that would have been available then and the rule is fixed
in advance.

At prediction timestamp \(t\), use only fields whose publication/availability
time is no later than \(t\). Specify whether a prediction is formed before or
after the timestamp's close/value becomes available. That decision determines
whether the current observation \(x_t\) may enter \(y_t\). Align predictions to
the target's realization interval and retain both decision and target timestamps
in any future result record.

For overlapping forecast horizons, disclose overlap and dependence. Do not treat
each overlapping prediction as an independent observation. Choose an aggregation
and uncertainty method suitable for serial dependence before testing. Report
per-window results as well as the declared aggregate; do not hide poor windows in
a single average. A later walk-forward study may repeat training and validation
on expanding or rolling history, but the window schedule and selection rules must
be specified in advance. Chronological splitting does not by itself prevent
leakage through revisions, feature engineering, global preprocessing, target
overlap, or repeated model selection.

## Causal preprocessing and leakage controls

The information rule is: a feature labeled available at timestamp \(t\) must be
computed solely from information actually available at or before \(t\). Apply
the following prohibitions:

- No centered moving average or other transformation using observations after
  the feature timestamp.
- No Phase 29 whole-series diffusion transform as a real-time feature.
- No normalization or imputation parameter fitted on the complete dataset.
- No backfill, interpolation, or other missing-value fill that uses future
  observations. Carrying a past value forward can be causal, but requires an
  explicit staleness limit and must not cross an undeclared session boundary.
- No selecting beta, window length, target, metric, model, exclusions, or report
  window from final test performance.
- No revised, corrected, or subsequently published value treated as known at an
  earlier time unless its original availability and revision timing are
  represented correctly.

Before a real study, record the following for every transformation:

| Transformation | Required declaration |
|---|---|
| Input | Exact fields, units, timestamp and availability time |
| History | Required prior observations/state; for this recurrence, the prior filter state plus current observation is sufficient after initialization |
| Initialization | Phase 30 uses `y[0] = x[0]`; declare the first valid observation and any state reset rule |
| Missing values | Dataset-specific behavior must be decided before analysis; no future fill. The Phase 30 API itself rejects non-finite inputs, so missing-value policy must not be improvised around it |
| Alignment | Output is labeled with the timestamp at which its input became available; any execution/forecast convention is separately specified |
| Warm-up | Declare whether early recursive outputs are retained or excluded and how any exclusion is chosen without test outcomes |
| Causality | Declare yes/no and the information set available at each output time |
| Test | Prefix-invariance test: create two inputs identical through index/timestamp `t`, alter only their suffixes, and verify outputs through `t` are identical |

The recurrence has no fixed finite lookback: its state summarizes prior inputs.
The initialization `y[0]=x[0]` passes the first value through and can affect
early outputs. The warm-up policy is therefore a study decision, not something
to estimate from the final test. Missing timestamps are especially important:
holding state, skipping an update, restarting, or using elapsed-time-adjusted
weights have different meanings. Select and test one policy based on the declared
sampling semantics; do not claim the current API decides it.

Fit scaling, clipping thresholds, feature selection, and learned imputation
parameters on training data only. Apply frozen parameters to validation and
test. Any operation performed separately at each timestamp must still obey the
availability rule. For manually supplied revised data, use point-in-time
vintages or mark the question unanswerable with that dataset.

## Baselines and fairness

Keep the initial comparison small:

1. **No smoothing:** the same input field passed through without transformation.
2. **Phase 30 exponential filter:** declared beta and `y[0]=x[0]` initialization.
3. **Causal trailing moving average:** declared window, startup partial-window
   convention, and timestamp alignment.
4. **Persistence:** include only in a forecasting track, predicting the
   prespecified target from its last available value under the target's exact
   definition.

Do not use Phase 29's centered moving average as a causal baseline. A trailing
mean and exponential filter are not fairly matched because their window and beta
look numerically similar. If the scientific question requires a matched
comparison, specify whether matching means measured frequency attenuation,
latency, or another property; estimate that match using analytical/synthetic
response or training data only. Otherwise report the declared configurations as
separate candidates, with their response and delay, and do not rank them as
equivalent filters. Choose a small candidate set before analysis. If validation
selects a configuration, also disclose the full candidate set and selection
rule. No method is a universal winner by this protocol.

## Two separate evaluation tracks

### Track A: signal-processing properties

Possible measurements include attenuation by frequency, phase/group delay,
step/abrupt-change response, sampling-frequency sensitivity, missing-data
sensitivity, and reconstruction error **only when an independent valid reference
exists**. For a synthetic signal, the generated clean input supplies known truth
and permits controlled reconstruction metrics. For real financial observations,
there is generally no observed noise-free “true signal.” A smoother's output is
not ground truth, and lower variation or a smoother plot alone is not evidence
of improved signal recovery or forecasting.

### Track B: forecasting (future phase only)

Do not implement a forecast in Phase 31. If later approved, first declare a
target such as a precisely defined future value/change, horizon, information
cutoff, timestamp alignment, and primary metric. Compare the same forecasting
procedure with and without the filter, plus appropriate naive baselines, using
chronological validation and a final untouched test. Report uncertainty and
variation across evaluation windows using a method selected before test
inspection. A historical transformation that is smoother has not necessarily
improved a forecast. The protocol makes no forecast claim.

## Multiplicity, experiment log, and research integrity

Keep the confirmatory question narrow and record all tried beta values,
instruments, periods, metrics, and exclusions. Repeated tuning, trying many
instruments/windows, changing metrics after seeing results, and publishing only
favorable configurations create selection bias. Overlapping horizons and serial
dependence reduce effective independent information. Label analyses exploratory
or confirmatory; exploratory results can motivate a later preregistered test but
cannot be silently relabeled as confirmatory.

Use this compact per-study log (a document or table is sufficient; no tracking
platform is proposed):

| Field | Record |
|---|---|
| Experiment ID | Stable human-readable identifier |
| Question and hypothesis | Wording fixed before evaluation |
| Data identifier and permission status | Provider/dataset/version, provenance, and evidence status; never store restricted raw data in the log |
| Configuration | Fields, dates, sampling, split boundaries, beta/windows, preprocessing, target/horizon, baselines, metrics, and selection rule |
| Code/version identifier | Release/package version or a recorded source-tree file digest; no Git operation is required |
| Evaluation dates | Actual train, validation, test, and execution/availability conventions |
| Metrics | Prespecified definitions, aggregate and per-window values, uncertainty method |
| Exclusions and failures | Counts, reasons, and whether decisions were made before outcomes were viewed |
| Study status | Exploratory or confirmatory; protocol deviations and date recorded |

Record failed configurations rather than silently dropping them. Preserve
protocol deviations and report their impact separately.

## Economic evaluation is a separate future question

Phase 31 calculates no investment returns. Even a statistically measurable
forecast improvement would not by itself establish an economically meaningful
result. A later economic study would additionally need a fully specified
decision rule and execution timing, realistic commissions/fees, bid-ask spread,
slippage and market impact where material, position sizing and exposure limits,
appropriate risk measures, turnover/capacity analysis, out-of-sample testing,
robustness across periods/instruments, and an appropriate benchmark. Any
decision must use only data available before execution. This protocol gives no
personalized investment recommendation.

## Readiness and unresolved decisions

**Real-data study status: not ready.** Before any authorized data is loaded, the
permission checklist above must be complete and the study table must have
prespecified values for its material fields. In particular, provider/use rights,
instrument, sampling/timestamps, observation period, target/task, beta/baseline
selection, quality exclusions, split dates, and primary metric remain unresolved.
The missing-value and warm-up policies must be fixed for the chosen data's
semantics. No empirical choices should be inferred from Phase 30's synthetic
defaults.

## Evidence interpretation and limitations

- **Established mathematical properties:** Phase 30 defines a causal recurrence,
  its initialization, and analytical response and includes a direct
  non-anticipation test for finite inputs.
- **Protocol decisions:** use chronological partitions, a sealed final test,
  train-only fitted preprocessing, timestamp availability checks, a small
  causal baseline set, and a compact experiment log.
- **Synthetic evidence:** Phase 30's deterministic synthetic results are limited
  to those signal definitions, parameter values, finite sampling, and metrics.
  They are not financial-data results.
- **Future empirical hypotheses:** whether this transformation helps a defined
  signal-processing or forecast task on authorized observations; whether any
  improvement generalizes or has net economic value. No such hypothesis has
  been tested here.
- **Possible biases:** revisions and survivorship, timestamp errors, session and
  missing-data policies, overlapping targets, tuning across many alternatives,
  dependence between observations, and choosing favorable windows or metrics.

Evidence supporting a future forecasting hypothesis would require improvement
on the prespecified validation/test metrics against the declared causal
baselines, with uncertainty and temporal variation reported and no test-driven
selection. Contradictory evidence would include no improvement, unstable or
reversed effects across the declared windows, leakage on audit, or gains erased
by realistic economic frictions in a separately approved economic evaluation.
Neither outcome would by itself establish universal applicability.

## Reproducibility and next step

This phase adds a protocol only; it does not add code, a data connector, a
template framework, or financial observations. The Phase 30 filter can be
reproduced with the existing synthetic study:

```python
from newton_lab.causal_diffusion_filter import run_causal_filter_study

study = run_causal_filter_study()
print(study.report_markdown)
```

This command reproduces synthetic Phase 30 results, not a financial study. The
next step is to resolve the data provider and permission/use terms, then fill
the initial study declaration before accessing any dataset. If authorization
cannot be established, stop at the protocol; do not retrieve data.
