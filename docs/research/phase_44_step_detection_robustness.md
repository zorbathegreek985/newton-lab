# Phase 44 — Robust Step Detection Under Confounding Signals

> **All records in this study are synthetic. These results do not establish
> real-market forecasting ability, causal identification, or trading
> profitability.** No NSE or external data was used.

## Research question

Can simple trend and oscillation-envelope checks reduce false abrupt-step
classifications while retaining sensitivity to weak true steps? The primary
target is a nonzero abrupt midpoint mean step, a narrower claim than a
difference between the first and second half of a record.

## Phase 43 motivation and baseline

Phase 43 showed that the two-half mean statistic responds to gradual drift and
that the two-half variance statistic responds to damped oscillation amplitude.
The Phase 44 baseline was checked against the existing code and report. It is
the same Phase 43/42 mean rule on 512 samples: calculate the two-half mean
difference divided by the standard error from the within-half sample
variances, then flag when the absolute statistic is at least 3. No formula or
threshold was changed.

## Signal ground truth and evaluation set

The fixed grid has 16 conditions. Each has 40 calibration records and 100
evaluation records. Every evaluation record is paired across all rules: the
same generated observations and seed produce the baseline, veto-rule, and raw
candidate-check results. Calibration seeds use
`44000 + condition_index * 1000 + replicate_index`; evaluation seeds use
`100000 + condition_index * 1000 + replicate_index`. The sets are disjoint.
Exact seed lists, condition definitions, and frozen thresholds are in
[`metadata.json`](../../reports/phase_44_step_detection/metadata.json).

The primary positive class consists of four pure step magnitudes (0.1, 0.2,
0.5, 1.0; unit Gaussian noise) plus a weak 0.2 step combined with either a
linear trend or a damped oscillation. These mixed records are positive because
the known step component is present; they are also reported separately so
their confounding structure remains visible. Negatives are stationary noise,
five trend slopes (0, 0.0005, 0.001, 0.002, 0.004; noise SD 0.75), and four
damping rates (0, 0.001, 0.002, 0.004; oscillation amplitude 1, frequency
0.10 cycles/sample, noise SD 0.2). Slope zero and damping zero provide
no-effect controls for their respective components. Unvaried settings remain
fixed. Phase 44 adds a compositional step input to the existing generator;
Phase 42/43 defaults and generated outputs remain unchanged.

## Candidate checks and frozen rules

The **trend score** is the absolute global OLS slope of value against sample
index. It is a simple full-record estimator. A step can also create a global
slope, so a nonzero score is not evidence that the record contains only a
trend.

The **envelope score** is `1 - late_amplitude / early_amplitude`. Each
amplitude is the magnitude of the sine/cosine coefficients from a least
squares fit with an intercept to one half of the record. The fitted frequency
is fixed at the generator's known 0.10 cycles/sample. This transparent check
assumes the tested oscillatory shape and known frequency; it is not a general
envelope estimator.

Calibration records alone select each score threshold by maximizing Youden's
J for detecting its component label: trend component present, or damped
envelope present. Exact ties select the largest threshold, favoring fewer
vetoes. The frozen trend threshold was `0.000452773` absolute value units per
sample (calibration J `0.63045`); the frozen envelope threshold was `0.156356`
(fractional amplitude reduction; calibration J `0.64792`). Evaluation seeds
were not used in this selection.

Both checks are evaluated as independent vetoes and as a combined veto:

- `baseline`: baseline mean-step flag.
- `trend_veto`: baseline flag **and** trend check is negative.
- `envelope_veto`: baseline flag **and** envelope check is negative.
- `combined_veto`: baseline flag **and** both checks are negative.

Raw trend and envelope flags are also evaluated against their own component
labels. A veto does not add evidence of a step; it only suppresses a baseline
flag when a confounder check fires.

## Primary evaluation results

The primary confusion counts pool 600 positive records and 1,000 negative
records. FPR is FP divided by all negative records; FNR is FN divided by all
positive records. Deltas below are absolute percentage-point changes against
the paired baseline.

| Rule | TP | FP | FN | TN | FPR | FNR | FPR change | FNR change |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Baseline | 339 | 300 | 261 | 700 | 30.0% | 43.5% | — | — |
| Trend veto | 0 | 0 | 600 | 1,000 | 0.0% | 100.0% | -30.0 pp | +56.5 pp |
| Envelope veto | 189 | 174 | 411 | 826 | 17.4% | 68.5% | -12.6 pp | +25.0 pp |
| Combined veto | 0 | 0 | 600 | 1,000 | 0.0% | 100.0% | -30.0 pp | +56.5 pp |

The trend veto removed all 300 baseline false positives, which were all in the
gradual-trend class, but it also vetoed every true positive. The global slope
score could not separate step-induced slope from actual drift under its
calibrated threshold. The combined rule therefore did no better than the
trend-only veto.

The envelope veto reduced false positives from 300 to 174 (126 fewer, a 42%
relative reduction) but lost 150 of the baseline's 339 true detections. The
false-negative rate rose by 25 percentage points. Baseline false positives
occurred in trends; stationary noise and damped oscillation controls did not
trigger the baseline in this evaluation, so the envelope check did not improve
those classes. Its additional effect was to veto genuine steps and some trend
flags.

### Weak steps and mixed cases

| Step condition | Baseline detections | Trend veto | Envelope veto | Combined veto |
|---|---:|---:|---:|---:|
| Pure step 0.1 | 2/100 | 0/100 | 1/100 | 0/100 |
| Pure step 0.2 | 26/100 | 0/100 | 17/100 | 0/100 |
| Pure step 0.5 | 99/100 | 0/100 | 52/100 | 0/100 |
| Pure step 1.0 | 100/100 | 0/100 | 59/100 | 0/100 |
| Step 0.2 + trend 0.001 | 96/100 | 0/100 | 60/100 | 0/100 |
| Step 0.2 + damped oscillation | 16/100 | 0/100 | 0/100 | 0/100 |

Baseline sensitivity was already low for the weakest steps (2% at 0.1 and 26%
at 0.2). The trend veto removed the remaining weak-step detections. The
envelope veto retained one 0.1 step and 17/100 pure 0.2 steps, while reducing
the 0.2 step-plus-damping detection from 16 to zero. These mixed signals have
an unambiguous step component by construction, but their overall statistics
are also shaped by the added component.

### Results by negative confounder

| Negative class | Records | Baseline false flags | Trend-veto false flags | Envelope-veto false flags |
|---|---:|---:|---:|---:|
| Stationary noise | 100 | 0 | 0 | 0 |
| Gradual trends | 500 | 300 | 0 | 174 |
| Damped oscillations | 400 | 0 | 0 | 0 |

The baseline mean rule did not flag the tested pure damped oscillations, so
this grid does not demonstrate an envelope-related baseline false positive to
remove. Raw component-check evaluation provides a separate view: the trend
check detected 462/500 trend-component records (92.4%) with 373/1,100 false
component flags (33.9%); the envelope check detected 397/400 damped-component
records (99.25%) with 457/1,200 false component flags (38.1%). The component
checks can recognize many target components, but their scores also occur in
records without those components. In particular, a step can influence global
OLS slope, and random half-window sinusoid fits can have different amplitudes
even without damping.

## Interpretation, limitations, and reproducibility

The trend check is too indiscriminate for a veto under this configuration; it
eliminates the tested false positives and all true detections. The envelope
check offers a smaller false-positive reduction at a substantial sensitivity
cost, and did not solve a baseline damped-oscillation false-positive problem
because none occurred in this grid. Neither check is a universal improvement.
The two statistics and thresholds are specific to the selected samples,
noise, generator, and known midpoint.

Calibration uses 40 seeds per condition and evaluation uses 100 disjoint seeds
per condition. This gives visible replicate variation but is not definitive
statistical validation; no confidence intervals are claimed. Thresholds are
selected by Youden's J on the same calibration pool across the declared
conditions, then frozen before the disjoint evaluation. Different calibration
mixtures or costs for missed steps versus false alarms could produce different
thresholds and veto trade-offs.

The envelope test uses the generator's known frequency and a sinusoidal model.
The global slope estimate intentionally illustrates that a step can look
trend-like in a whole-record regression. The fixed midpoint and Gaussian noise
omit unknown change locations, heavy tails, missing data, irregular cadence,
and many other processes. This work does not establish causal identification
or real-world performance.

Run the study from the project root with
`.venv\Scripts\python.exe -m newton_lab.step_detection_robustness`. The output
directory contains [per-record data](../../reports/phase_44_step_detection/per_record_results.csv),
[summary metrics](../../reports/phase_44_step_detection/summary.csv),
[metadata](../../reports/phase_44_step_detection/metadata.json), an
[overall rule error plot](../../reports/phase_44_step_detection/rule_error_rates.png),
and a [step sensitivity plot](../../reports/phase_44_step_detection/weak_step_detection.png).
The per-record file includes calibration and evaluation rows, ground truth,
settings, seeds, every diagnostic statistic, and each rule output. Repeated
runner executions reproduced the CSV and metadata SHA-256 hashes exactly.

Phase 32 remains `BLOCKED_AUTHORIZATION`. Phase 42 and 43 outputs and raw NSE
data were not modified.

## Conclusion and next research question

Under this synthetic design, the trend veto removed baseline false positives
but also rejected all true steps. The envelope veto reduced false positives
less and raised missed-step frequency. Mixed records show that real step
components can coexist with trend or damping evidence; a veto can discard a
true change because another structure is present. No tested rule is a clear
winner across these outcomes.

**Next research question:** On a preregistered grid with separate calibration
seeds, can a local before/after level comparison plus explicit slope removal
distinguish a weak step from a smooth trend without sacrificing step-plus-trend
cases?
