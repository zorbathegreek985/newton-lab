# Phase 46 — Independent Noise-Level Sensitivity Study

**Scope:** controlled synthetic detection experiment. It uses no external or NSE data and does not establish performance on real financial data, causal identification, forecasting, or profitability.

## Research question

How do independently varied Gaussian noise levels affect false alarms and missed step detections for the fixed Phase 44/45 methods? Phase 45 used noise levels that were coupled to signal families; this experiment crosses noise with each tested family while holding its other parameters fixed.

## Methods and frozen thresholds

The four compared rules are unchanged:

| Method | Statistic and decision rule | Threshold source |
|---|---|---|
| Phase 44 baseline | Absolute two-half mean Z statistic `>= 3` | Fixed Phase 42/44 threshold |
| Phase 45 local raw, w=32 | Absolute post-window mean minus pre-window mean `>= 0.20193448253267382` | Phase 45 calibration, frozen |
| Phase 45 local trend-adjusted, w=32 | Absolute local difference minus full-record OLS slope times 32 `>= 0.15418153275303856` | Phase 45 calibration, frozen |
| Phase 45 local trend-adjusted, w=64 | Absolute local difference minus full-record OLS slope times 64 `>= 0.19161592877807834` | Phase 45 calibration, frozen |

The midpoint is the known sample 256 in records of length 512. The methods evaluate a retrospective comparison at this oracle location; they do not search for a change point. No Phase 46 thresholds were fitted. The calibration partition is recorded separately, but the methods and thresholds stayed fixed from Phases 42–45 throughout the held-out evaluation. Calibration labels and evaluation labels were not used to choose or change any threshold.

## Factorial design

Noise standard deviations were `0.5`, `1.0`, and `1.5`. The same seeded unit-Gaussian residual was scaled to each level within a matched signal configuration, so the latent signal and seed were held fixed while noise changed.

| Family | Non-noise settings | Ground truth |
|---|---|---|
| Stationary null | Zero latent signal | No step |
| Gradual trend | Linear slope `0.001` per sample | No step |
| Pure step | Midpoint steps of `0.1`, `0.2`, `0.5`, and `1.0` | Step present |
| Step plus trend | The same four midpoint step sizes plus linear slope `0.001` | Step present |

This is 10 base signal configurations crossed with 3 noise levels, or 30 factorial cells. Signal duration, midpoint, trend, and step settings remain fixed across noise levels. The implementation reuses the Phase 42 generator's latent signal and unit-noise realization, scaling only its noise residual; no existing generator defaults or Phase 42–45 outputs were changed.

Each cell has 40 calibration and 100 evaluation replicates. That is 1,200 calibration and 3,000 evaluation signal records, each scored by all four methods. Seed IDs are deliberately reused across the three noise levels for matched comparisons: 400 unique calibration seeds and 1,000 unique evaluation seeds. The seed ranges are disjoint. Each evaluation condition therefore has 100 independent seed replicates; pure-step and step-plus-trend family/noise strata pool four step magnitudes and have 400 positive records. Pointwise 95% Wilson intervals describe binomial rate uncertainty. They do not adjust for the many reported comparisons or quantify uncertainty outside this generator design.

The positive label is the known abrupt step, including step-plus-trend cells. Stationary and gradual-trend cells are negative. For each subgroup, FPR is FP / (FP + TN), and FNR is FN / (TP + FN); a rate is undefined where its class denominator is zero. Full TP/FP/FN/TN counts and denominators appear in the summary CSV.

## Held-out error rates by family and noise

The entries show rate (errors / class denominator); family/noise summaries use 100 records for stationary and trend families, and 400 records for each step family/noise cell. Wilson intervals are in `summary.csv`.

### False-positive rate

| Family | Noise SD | Phase 44 baseline | Raw w=32 | Adjusted w=32 | Adjusted w=64 |
|---|---:|---:|---:|---:|---:|
| Stationary null | 0.5 | 2% (2/100) | 5% (5/100) | 20% (20/100) | 0% (0/100) |
| Stationary null | 1.0 | 2% (2/100) | 43% (43/100) | 52% (52/100) | 22% (22/100) |
| Stationary null | 1.5 | 2% (2/100) | 57% (57/100) | 67% (67/100) | 45% (45/100) |
| Gradual trend | 0.5 | 100% (100/100) | 17% (17/100) | 23% (23/100) | 1% (1/100) |
| Gradual trend | 1.0 | 42% (42/100) | 37% (37/100) | 49% (49/100) | 32% (32/100) |
| Gradual trend | 1.5 | 11% (11/100) | 53% (53/100) | 68% (68/100) | 53% (53/100) |

Under stationary nulls, the fixed local thresholds generated more false alarms as noise increased. On smooth trends, the local w=64 adjustment had low false alarms at the lowest noise level but its FPR rose to 53% at the highest. The Phase 44 standardized statistic's response differed: it was nearly unchanged on the stationary null, while its trend false alarms declined as noise obscured the fixed drift.

### False-negative rate

| Family | Noise SD | Phase 44 baseline | Raw w=32 | Adjusted w=32 | Adjusted w=64 |
|---|---:|---:|---:|---:|---:|
| Pure steps, pooled magnitudes | 0.5 | 22.0% (88/400) | 34.0% (136/400) | 27.5% (110/400) | 38.5% (154/400) |
| Pure steps, pooled magnitudes | 1.0 | 44.5% (178/400) | 31.0% (124/400) | 24.0% (96/400) | 32.8% (131/400) |
| Pure steps, pooled magnitudes | 1.5 | 53.2% (213/400) | 27.3% (109/400) | 22.0% (88/400) | 30.2% (121/400) |
| Step plus trend, pooled magnitudes | 0.5 | 0.0% (0/400) | 30.2% (121/400) | 29.5% (118/400) | 38.5% (154/400) |
| Step plus trend, pooled magnitudes | 1.0 | 2.5% (10/400) | 28.7% (115/400) | 24.8% (99/400) | 33.2% (133/400) |
| Step plus trend, pooled magnitudes | 1.5 | 25.2% (101/400) | 24.2% (97/400) | 20.2% (81/400) | 32.0% (128/400) |

These pooled miss rates can hide the step-size effect. The next table reports the held-out detection count out of 100 for each magnitude, noise level, and method; `summary.csv` also contains rates and Wilson intervals.

### Detection by step magnitude

| Family | Noise SD | Step | Baseline | Raw w=32 | Adjusted w=32 | Adjusted w=64 |
|---|---:|---:|---:|---:|---:|---:|
| Pure step | 0.5 | 0.1 | 16 | 17 | 30 | 9 |
| Pure step | 0.5 | 0.2 | 96 | 47 | 60 | 37 |
| Pure step | 0.5 | 0.5 | 100 | 100 | 100 | 100 |
| Pure step | 0.5 | 1.0 | 100 | 100 | 100 | 100 |
| Pure step | 1.0 | 0.1 | 0 | 40 | 60 | 36 |
| Pure step | 1.0 | 0.2 | 23 | 51 | 58 | 45 |
| Pure step | 1.0 | 0.5 | 99 | 85 | 86 | 88 |
| Pure step | 1.0 | 1.0 | 100 | 100 | 100 | 100 |
| Pure step | 1.5 | 0.1 | 0 | 60 | 74 | 47 |
| Pure step | 1.5 | 0.2 | 7 | 58 | 65 | 56 |
| Pure step | 1.5 | 0.5 | 80 | 74 | 76 | 79 |
| Pure step | 1.5 | 1.0 | 100 | 99 | 97 | 97 |
| Step plus trend | 0.5 | 0.1 | 100 | 31 | 34 | 13 |
| Step plus trend | 0.5 | 0.2 | 100 | 48 | 48 | 33 |
| Step plus trend | 0.5 | 0.5 | 100 | 100 | 100 | 100 |
| Step plus trend | 0.5 | 1.0 | 100 | 100 | 100 | 100 |
| Step plus trend | 1.0 | 0.1 | 91 | 45 | 54 | 33 |
| Step plus trend | 1.0 | 0.2 | 99 | 48 | 58 | 43 |
| Step plus trend | 1.0 | 0.5 | 100 | 92 | 89 | 91 |
| Step plus trend | 1.0 | 1.0 | 100 | 100 | 100 | 100 |
| Step plus trend | 1.5 | 0.1 | 37 | 58 | 69 | 46 |
| Step plus trend | 1.5 | 0.2 | 62 | 62 | 69 | 51 |
| Step plus trend | 1.5 | 0.5 | 100 | 84 | 82 | 75 |
| Step plus trend | 1.5 | 1.0 | 100 | 99 | 99 | 100 |

## Interpretation and limits

Increasing noise raised local-method false alarms on stationary noise and generally raised false alarms on trends at the higher noise levels. The Phase 44 baseline behaved differently across families: its stationary-null FPR stayed at 2/100, while its trend FPR fell from 100/100 to 11/100 as the trend became less distinct against noise.

The local thresholds were calibrated in Phase 45 at a different mix of noise levels and are deliberately not retuned here. For several step families their miss rates fell as noise increased. That does not mean more noise improved the underlying information: with a fixed absolute threshold, larger random excursions can cross the detection boundary in both positive and negative records. The high stationary-null FPRs show this directly. The baseline's pure-step and step-plus-trend miss rates generally increased with noise, especially for weaker steps. Step size and method interact strongly with noise: strong steps were usually detected, while the weakest steps had very different outcomes across thresholds and signal families.

The results are descriptive for this specified additive-Gaussian generator, fixed slope, known midpoint, and frozen thresholds. The paired seeds make noise-level comparisons less sensitive to which underlying random draw was used, but this finite grid cannot establish universal detector behavior. Wilson intervals are marginal; no multiple-comparison adjustment or formal paired interval for rate differences was calculated. There are no unknown change locations, non-Gaussian noise, missing observations, parameter estimation errors, real market series, or causal claims in this study.

## Outputs and reproduction

Run from the project root with the existing virtual environment:

`.venv\Scripts\python.exe -m newton_lab.noise_sensitivity`

Files are under `reports/phase_46_noise_sensitivity/`:

- `per_record_results.csv` — calibration and evaluation method predictions with exact signal condition, noise, step, slope, seed, statistic, frozen threshold, source, and label.
- `summary.csv` — held-out confusion counts, class denominators, error rates, and pointwise 95% Wilson intervals by family, noise, and step magnitude.
- `metadata.json` — factorial settings, seed lists and policy, frozen method definitions, and threshold origins.
- `error_rates_by_noise.png` — FPR for negative families and FNR for positive families against noise SD, with Wilson intervals.
- `step_magnitude_sensitivity.png` — detection rates by step size, noise SD, and method for pure and trend-confounded steps.

The full experiment ran twice. SHA-256 hashes matched for all five artifacts. Matplotlib emitted a nonfatal warning that it could not write its user font cache due to permissions; both plot files were created and reproduced identically.

**Conclusion:** Under this experiment's frozen absolute local thresholds, noise scaling materially changed false alarms and missed detections, and the direction and size of the change depended on signal family and step magnitude. The standardized baseline had a different response, especially for pure trends and step-plus-trend signals. No one rule dominated every family/noise/step-size combination. These are controlled synthetic observations, not evidence about real-world or financial data.

Phase 32 remains `BLOCKED_AUTHORIZATION`.
