# Phase 43 — Diagnostic Robustness and Signal Ambiguity

> **This phase uses synthetic signals only. Its results do not establish real-
> market forecasting ability, trading profitability, or universal diagnostic
> performance.** No NSE or external data was used.

## Research question

Can Phase 42's two-half mean and variance diagnostics distinguish an abrupt
step from other generating processes, or do different processes produce
similar results? This study varies trend slope, step magnitude, damping, and
noise, while retaining fixed sample size and diagnostic thresholds. It reuses
the Phase 42 generator and diagnostic formulas without adding a change-point
algorithm.

## Phase 42 findings under investigation

Phase 42 flagged 8/9 trend records with its mean diagnostic and 5/9 damped
oscillations with its variance diagnostic. Code and result inspection found
no formula defect. The mean statistic is a first-half versus second-half mean
difference divided by an estimated within-half standard error. The variance
statistic is the ratio of the larger to smaller within-half sample variance.
Both use the known midpoint and do not search for an unknown change time.

The formulas measure differences between halves; they do not identify an
abrupt transition or its cause by themselves.

## Controlled design

Every record has 512 samples and every configuration has 20 deterministic
replicates. Mean-study seeds use `43000 + configuration_index * 100 +
replicate_index`; variance-study seeds use `53000 + configuration_index * 100
+ replicate_index`. Exact per-row settings and seeds are in
[`parameter_results.csv`](../../reports/phase_43_diagnostic_robustness/parameter_results.csv).
Thresholds and seed rules are also recorded in
[`metadata.json`](../../reports/phase_43_diagnostic_robustness/metadata.json).

The **mean-step target** is a nonzero abrupt midpoint mean shift. Positive
configurations use step sizes 0.1, 0.2, and 0.5 with unit observation-noise
standard deviation. Negative controls are iid Gaussian noise, AR(1) with phi
0.8 and observation-noise SD 0.5, and linear trends. The trend study varies
slope over 0, 0.0005, 0.001, 0.002, and 0.004 per sample at noise SD 0.75;
slope zero is the no-trend control. A separate noise sweep tests SD 0.25 and
1.5 at slope 0.002; SD 0.75 is represented in the slope sweep.

The **random variance-step target** is a midpoint change in Gaussian noise
standard deviation, from 1 to a multiplier greater than 1. Positive
multipliers are 1.1, 1.4, and 2.0; multiplier 1.0 is the no-change control.
Negative controls include equal-variance iid noise and damped sinusoids. For
the sinusoid, amplitude 1 and frequency 0.10 cycles/sample are fixed. Damping
rates 0, 0.001, 0.002, and 0.004 are tested at noise SD 0.2; noise SD 0.05
and 0.5 are tested at damping 0.002. These are one-factor-at-a-time
comparisons around the Phase 42 settings.

The thresholds were not calibrated or selected on these evaluation runs. The
positive and negative classes are narrowly defined as above. The reported
rates are conditional on these specified controls.

## Classification definitions and results

The mean diagnostic flags when the absolute two-half statistic is >= 3. A true
positive is a flagged midpoint mean-step record; a false positive is a flagged
trend, AR(1), or iid control. The variance diagnostic flags when the two-half
variance ratio is >= 2. A true positive is a flagged random variance step; a
false positive is a flagged equal-variance or damped-sinusoid control. Rates
use FP / all negative controls and FN / all positive cases.

| Target property | TP | FP | FN | TN | False-positive rate | False-negative rate |
|---|---:|---:|---:|---:|---:|---:|
| Abrupt midpoint mean step (60 positive, 180 negative) | 26 | 105 | 34 | 75 | 58.3% | 56.7% |
| Random midpoint variance step (60 positive, 160 negative) | 30 | 59 | 30 | 101 | 36.9% | 50.0% |

### Mean statistic parameter sensitivity

| Configuration | Flagged / 20 |
|---|---:|
| Mean step 0.1 | 2/20 |
| Mean step 0.2 | 4/20 |
| Mean step 0.5 | 20/20 |
| Trend slope 0 (noise SD 0.75) | 0/20 |
| Trend slope 0.0005 | 6/20 |
| Trend slope 0.001 | 19/20 |
| Trend slope 0.002 | 20/20 |
| Trend slope 0.004 | 20/20 |
| Trend noise SD 0.25 at slope 0.002 | 20/20 |
| Trend noise SD 1.5 at slope 0.002 | 17/20 |
| Stationary iid control | 1/20 |
| AR(1) control | 2/20 |

The mean overlap is explained by the statistic responding to a sustained
first-half/second-half level difference, whether caused by a step or gradual
drift. With the Phase 42 slope 0.002 and noise SD 0.25, 0.75, and 1.5, the new
runs flagged 57/60 records (20/20, 20/20, and 17/20). The difference from
Phase 42's 8/9 is consistent with its three seeds and finite-sample variability
at the noisiest setting. There is no implementation defect: the formula does
what it is defined to do. It is too broad to support “abrupt mean step” without
additional shape information.

### Variance statistic parameter sensitivity

| Configuration | Flagged / 20 |
|---|---:|
| Post-change SD multiplier 1.0 (null control) | 0/20 |
| Post-change SD multiplier 1.1 | 0/20 |
| Post-change SD multiplier 1.4 | 10/20 |
| Post-change SD multiplier 2.0 | 20/20 |
| Damping 0 at noise SD 0.2 | 0/20 |
| Damping 0.001 at noise SD 0.2 | 0/20 |
| Damping 0.002 at noise SD 0.2 | 19/20 |
| Damping 0.004 at noise SD 0.2 | 20/20 |
| Noise SD 0.05 at damping 0.002 | 20/20 |
| Noise SD 0.5 at damping 0.002 | 0/20 |
| Stationary iid control | 0/20 |

The damped-oscillation overlap is explained by its decaying deterministic
envelope: early and late half-record variances differ even though random noise
variance stays constant. With damping 0, no run crossed the threshold; at
damping 0.002, 19/20 did. Added noise masks the envelope: at damping 0.002,
detection fell from 20/20 at noise SD 0.05 to 0/20 at 0.5. The ratio reports
unequal observed half variances, but does not identify a random variance
regime change.

Seed sensitivity is visible near thresholds: 19/20 for damping 0.002, 10/20
for variance multiplier 1.4, 6/20 for trend slope 0.0005, and 4/20 for mean
step 0.2. These rates describe only the finite grid, not calibrated error
probabilities for other data or processes.

## Ambiguity and interpretation

The mean statistic detects a difference in half-record averages. A trend and
an abrupt step can have similar values; the statistic alone does not
distinguish their shape or identify a generating mechanism. The variance ratio
detects unequal dispersion between halves. A decaying oscillation envelope and
a stochastic variance step can both create that pattern. Similar outputs
indicate ambiguity in these summaries, not that the processes are identical.

Keep these claims separate:

1. **Pattern detection:** the statistic exceeds its fixed threshold.
2. **Process identification:** evidence distinguishes among trends, steps,
   damped oscillations, dependence, and other candidates.
3. **Causal explanation:** evidence establishes why the pattern arose.

This experiment supports only the first claim for its generated records. It
does not infer real-world causes or universal limits for diagnostics.

## Outputs, limitations, and implementation review

- [`parameter_results.csv`](../../reports/phase_43_diagnostic_robustness/parameter_results.csv)
  contains 460 runs with exact settings, seeds, target labels, statistics, and
  detections.
- [`confusion_summary.csv`](../../reports/phase_43_diagnostic_robustness/confusion_summary.csv)
  contains confusion counts and class-conditional error rates.
- [`mean_diagnostic_sensitivity.png`](../../reports/phase_43_diagnostic_robustness/mean_diagnostic_sensitivity.png)
  plots flag rates by trend slope and mean-step magnitude.
- [`variance_diagnostic_sensitivity.png`](../../reports/phase_43_diagnostic_robustness/variance_diagnostic_sensitivity.png)
  plots flag rates by damping and variance-step magnitude.

The study uses Gaussian noise, fixed midpoint windows, 20 seeds per setting,
and hand-declared thresholds. It does not adjust standard errors for AR
dependence, search unknown change points, or calibrate thresholds on separate
development data. The two Phase 42 overlaps do not reveal an implementation
bug; they show that interpretation must stay narrow. More seeds could refine
these rates. Other processes, noise distributions, windows, sample sizes, and
thresholds could change them.

The Phase 42 generator defaults are unchanged. Phase 42 result files, NSE
files, and the Phase 32 authorization state were not modified. Phase 32 remains
`BLOCKED_AUTHORIZATION`.

## Conclusion and next research question

For these synthetic settings, both statistics flag some designed steps but
also respond to non-step processes with similar between-half differences.
Trend slope, damping, signal-to-noise balance, and finite samples explain the
two Phase 42 overlaps; the diagnostic formulas behaved as implemented. The
flags support claims about measured summaries, not unique process
identification or causality.

**Next research question:** On separate calibration and evaluation seeds, can
adding a simple within-half trend or oscillation-envelope check reduce false
step classifications while preserving detection of weak abrupt changes?
