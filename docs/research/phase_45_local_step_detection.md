# Phase 45 — Local Level-Shift Detection with Trend Adjustment

**Scope:** synthetic methodological experiment only. These results do not establish performance on real-world data, causal identification, or trading profitability. No NSE or external data was used.

## Research question

Can a local before-and-after level comparison, with explicit linear-trend adjustment, separate an abrupt midpoint step from gradual drift more effectively than the Phase 44 two-half mean diagnostic and its global trend veto?

## Phase 44 baseline

Phase 44 uses 512 observations and flags when the absolute two-half mean statistic is at least 3. Its global trend veto retains a baseline flag only when the absolute full-record OLS slope is below a calibration threshold. The envelope veto and combined veto are also carried forward unchanged. The default Phase 44 thresholds reproduced here are `0.0004527730452371508` for absolute slope and `0.1563558789122506` for fractional envelope decay. The Phase 44 report's overall baseline result was reproduced: TP 339, FP 300, FN 261, TN 700; the trend veto yielded TP 0, FP 0, FN 600, TN 1,000.

## Signals and ground truth

The study reused the 16 Phase 44 conditions and its exact signal generator. Each record has 512 samples. The grid contains stationary Gaussian noise; linear-trend-only records with slopes 0, 0.0005, 0.001, 0.002, and 0.004; damped-oscillation-only records with damping rates 0, 0.001, 0.002, and 0.004; pure midpoint steps of 0.1, 0.2, 0.5, and 1.0; and a 0.2 step combined separately with a 0.001 linear trend or 0.002 damping. Noise SD is 1.0 for stationary and pure-step records, 0.75 for trend records, 0.2 for oscillation-only records, and 1.0 for step-plus-confounder records. All settings are in `reports/phase_45_local_step_detection/metadata.json`.

The positive label is the known abrupt midpoint step, including the two step-plus-confounder conditions. Negative labels have no step. A smooth trend is always negative, even when the two halves have different means. This known midpoint is oracle location information: these retrospective comparisons do not search for a change point and are not causal or real-time detectors.

## Methods and calibration

For a record of length `n`, midpoint `m = floor(n/2)`, and window width `w` on each side, the unadjusted statistic is

`D_w = mean(x[m:m+w]) - mean(x[m-w:m])`.

The predeclared widths are 16, 32, and 64 samples on each side. Each local decision uses `abs(D_w) >= threshold_w`, with a separate threshold for each window and method. The trend-adjusted method estimates one full-record ordinary least-squares slope `b` against sample index and uses

`D_adjusted,w = D_w - b*w`.

The window-center separation is `w`, hence the expected linear-trend difference is `b*w`. The step can influence the full-record slope estimate, so the adjustment can subtract part of a real step. Both methods use a two-sided absolute score.

Thresholds were selected using calibration records only by maximizing Youden's J for the known step label; ties follow the existing Phase 44 selector rule. Window sizes were fixed in advance and all are reported; none was selected using evaluation performance. The 40 calibration replicates per condition use seed blocks beginning at 44000, 45000, …, 59000, with 40 seeds per block. The 100 evaluation replicates per condition use blocks beginning at 100000, 101000, …, 115000, with 100 seeds per block. This is 640 calibration and 1,600 evaluation condition-records, with disjoint seed sets; the same evaluation record is paired across methods. The replicate count makes rates interpretable for this small synthetic design but does not provide broad population-level evidence.

| Local rule | Frozen absolute threshold | Calibration Youden J |
|---|---:|---:|
| Raw, w=16 | 0.399685 | 0.3750 |
| Trend adjusted, w=16 | 0.376292 | 0.3608 |
| Raw, w=32 | 0.201934 | 0.5317 |
| Trend adjusted, w=32 | 0.154182 | 0.4917 |
| Raw, w=64 | 0.143648 | 0.4983 |
| Trend adjusted, w=64 | 0.191616 | 0.5192 |

## Held-out results

There are 600 positive and 1,000 negative evaluation records. FPR and FNR are proportions of negatives and positives, respectively. Changes are absolute percentage-point differences from the paired Phase 44 baseline.

| Method | TP | FP | FN | TN | FPR | FNR | Δ FPR | Δ FNR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Phase 44 baseline | 339 | 300 | 261 | 700 | 30.0% | 43.5% | — | — |
| Phase 44 trend veto | 0 | 0 | 600 | 1,000 | 0.0% | 100.0% | -30.0 pp | +56.5 pp |
| Phase 44 envelope veto | 189 | 174 | 411 | 826 | 17.4% | 68.5% | -12.6 pp | +25.0 pp |
| Phase 44 combined veto | 0 | 0 | 600 | 1,000 | 0.0% | 100.0% | -30.0 pp | +56.5 pp |
| Local raw, w=16 | 253 | 135 | 347 | 865 | 13.5% | 57.8% | -16.5 pp | +14.3 pp |
| Local trend adjusted, w=16 | 271 | 156 | 329 | 844 | 15.6% | 54.8% | -14.4 pp | +11.3 pp |
| Local raw, w=32 | 385 | 212 | 215 | 788 | 21.2% | 35.8% | -8.8 pp | -7.7 pp |
| Local trend adjusted, w=32 | 418 | 288 | 182 | 712 | 28.8% | 30.3% | -1.2 pp | -13.2 pp |
| Local raw, w=64 | 440 | 258 | 160 | 742 | 25.8% | 26.7% | -4.2 pp | -16.8 pp |
| Local trend adjusted, w=64 | 347 | 91 | 253 | 909 | 9.1% | 42.2% | -20.9 pp | -1.3 pp |

The local rules produced several trade-offs rather than one uniformly best setting. Raw w=32 and trend-adjusted w=32 lowered both aggregate error rates modestly, with adjusted w=32 reducing FNR most (by 13.2 pp) at only a 1.2 pp FPR reduction. Trend-adjusted w=64 made the strongest FPR reduction and left aggregate FNR close to baseline. These pooled rates conceal meaningful differences by signal class and do not establish a universally preferable threshold.

## Weak steps, trends, and mixed signals

| Condition | Phase 44 baseline detections | Raw w=16 | Adjusted w=16 | Raw w=32 | Adjusted w=32 | Raw w=64 | Adjusted w=64 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Pure step 0.1 | 2/100 | 21/100 | 28/100 | 41/100 | 50/100 | 41/100 | 32/100 |
| Pure step 0.2 | 26/100 | 25/100 | 29/100 | 48/100 | 58/100 | 69/100 | 45/100 |
| Pure step 0.5 | 99/100 | 58/100 | 58/100 | 84/100 | 84/100 | 97/100 | 81/100 |
| Pure step 1.0 | 100/100 | 99/100 | 98/100 | 100/100 | 100/100 | 100/100 | 100/100 |
| Step 0.2 + trend 0.001 | 96/100 | 31/100 | 36/100 | 65/100 | 68/100 | 72/100 | 43/100 |
| Step 0.2 + damping 0.002 | 16/100 | 19/100 | 22/100 | 47/100 | 58/100 | 61/100 | 46/100 |

The local methods improved weak pure-step detections in some configurations, especially w=32 or raw w=64, but generally lost detections on the stronger 0.5 step relative to the baseline. The trend-adjusted w=64 rule detected fewer 0.2 step-plus-trend cases than raw w=64 (43 versus 72), consistent with its full-record slope subtracting step-related variation as well as drift. The trend veto missed every positive in the Phase 44 evaluation, so all local rules retained much greater positive sensitivity than that veto.

On the 500 pure-trend negatives, Phase 44 baseline produced 300 false alarms and trend veto produced none. Local false alarms were 73/500 (raw w=16), 81/500 (adjusted w=16), 157/500 (raw w=32), 218/500 (adjusted w=32), 207/500 (raw w=64), and 70/500 (adjusted w=64). On the 100 stationary-noise records, local false alarms ranged from 21/100 to 65/100, despite zero baseline false alarms. Thus local differences can also respond to random local imbalance. On damped-oscillation negatives (400 records), local false alarms ranged from 0/400 to 44/400.

The 64-sample window reduced false alarms on trends, but increased exposure to other broad local fluctuations and did not preserve the strongest sensitivity on the step-plus-trend case. The step magnitude results and all conditions are in `summary.csv`; condition and record detail is in `per_record_results.csv`.

Noise levels were not varied factorially within the same signal conditions. The experiment includes different fixed noise SDs across condition families, so its results cannot isolate a general noise-level effect. The local-window comparison is direct for this grid; a controlled noise-dose response remains untested.

## Outputs and reproducibility

The reproducible command is:

`.venv\Scripts\python.exe -m newton_lab.local_step_detection`

Environment used: Python 3.12.10 and the existing project virtual environment. Outputs are in `reports/phase_45_local_step_detection/`:

- `per_record_results.csv` — matched method-level records with condition, parameters, seed, truth, threshold, score, and detections.
- `summary.csv` — confusion counts and error rates by method, condition, group, and step magnitude.
- `metadata.json` — condition definitions, calibration/evaluation seeds, frozen thresholds, rules, and Phase 44 thresholds.
- `method_error_rates.png` — held-out overall FPR and FNR by method.
- `step_magnitude_sensitivity.png` — detection rate by step magnitude for preselected methods.

The full experiment was run twice; SHA-256 hashes matched for both CSVs, metadata, and both figures. Matplotlib emitted a non-fatal warning that it could not write its user font cache due to permissions; both plot files were created and matched on repeat.

## Limitations and conclusion

This benchmark uses a known midpoint, simple additive signal generators, Gaussian noise, and a finite predeclared condition grid. Thresholds are calibrated and evaluated on distinct generated seeds but within the same specified generator families. The estimator does not locate unknown change points, account for nonlinear trends, estimate uncertainty for a changing process, or handle broad classes of nonstationarity. The full-record slope can be biased by a true step. Aggregate gains depend on the chosen condition mixture and calibration objective. No formal confidence intervals were calculated; denominators and per-condition replicate rates are provided instead.

Within this synthetic design, local comparisons can retain step detections that a global trend veto suppresses entirely. Window width changes the false-alarm/sensitivity balance, and full-record trend adjustment is not consistently beneficial: it helps the w=64 false-alarm rate on trends but can reduce mixed step detection. Local comparisons are not sufficient by themselves to distinguish all abrupt steps from noise or smooth drift, and statistical step detection does not identify the underlying physical or financial generating process.

**Next research question:** How do these frozen local rules perform when the change location is unknown and the evaluation signals include a separately controlled range of noise levels, with all search and multiplicity handling specified before held-out evaluation?

Phase 32 remains `BLOCKED_AUTHORIZATION`; this synthetic study provides no authorization evidence and does not use market data.
