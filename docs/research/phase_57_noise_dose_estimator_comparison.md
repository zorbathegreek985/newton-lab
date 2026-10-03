# Phase 57 — Noise-Dose and Estimator-Form Comparison

## Executive summary

This frozen synthetic benchmark tested how additive observation-noise dose and fit form affect recovery-rate estimation and exponent identifiability in the Phase 53 fold, transcritical and positive supercritical-pitchfork normal forms. It generated 6,000 held-out traces and made 18,000 estimator attempts. The primary comparison was between two fits using only the sampled observations. A correct-normal-form fit was retained as a separate model-informed diagnostic.

For both observation-only fits, increasing noise sharply worsened rate error at the smallest perturbation. The free-offset exponential remained numerically valid in all 6,000 evaluation attempts, while its mean relative absolute error at amplitude `a=0.0025` rose from below 0.2% in the zero-noise diagnostic to 2.61 (fold), 21.95 (transcritical), and 7.40 (pitchfork) at noise fraction `0.05`. The tail-offset log-linear alternative was frequently invalid at nonzero noise and, when valid, commonly returned rates near zero. It was not an improvement. The exact-form fit had low errors on these exact generator equations, as expected from its model information; this does not establish equal-information superiority or performance under model mismatch.

Only 3 of 36 nonzero-noise, model/amplitude/method exponent conditions met the unchanged Phase 55 p-identifiability criterion. All three were from the free-offset fit at `a=0.04`: fold at noise fraction `0.0025`, and transcritical at `0.0025` and `0.01`. No tail-offset condition met it. Zero-noise exponent rows are diagnostics and are marked not assessed because the seed-labelled trajectories repeat exactly and cannot provide replicate uncertainty.

The findings answer the question conditionally: noise dose and estimator form interacted in the tested design, but the direction and magnitude depended on model and amplitude. Individual-rate accuracy did not ensure exponent identifiability. The experiment is a synthetic method-capability study; its equations encode the scaling laws.

## Research question and hypotheses

**Question:** How do observation-noise level and estimator form interact to affect recovery-rate estimation and exponent identifiability in the three existing normal forms?

The frozen protocol considered increasing dose, estimator-form differences, their interaction, and translation from individual-rate estimates to exponent identification. It did not require either observation-only estimator to win. The complete protocol is [Phase 57 protocol](phase_57_noise_dose_estimator_comparison_protocol.md), with machine-readable configuration in [`protocol.json`](../../reports/phase_57_noise_dose_estimator_comparison/protocol.json).

## Models, estimands and methods

For `μ > 0`, only the specified stable branches were tested:

| Normal form | Selected equilibrium | Theoretical local rate | `(q,p)` for `x* ~ μ^q`, `ρ ~ μ^p` |
|---|---|---|---|
| Fold: `dx/dt = μ − x²` | `+√μ` | `2√μ` | `(0.5, 0.5)` |
| Transcritical: `dx/dt = μx − x²` | `μ` | `μ` | `(1, 1)` |
| Supercritical pitchfork: `dx/dt = μx − x³` | `+√μ` | `2μ` | `(0.5, 1)` |

The margins were `(0.04, 0.09, 0.16, 0.25, 0.36)` and the control settings were `u=1+μ`, with boundary `1`. Initial displacement was set by `x(0)=x*(1+a)` at `a=0.0025` and `0.04`. Observation noise was iid additive Gaussian with `σ=ηx*`, at dose fractions `η=(0, 0.0025, 0.01, 0.05)`. The lower dose makes noise SD equal the smallest initial displacement; `0.01` is the Phase 55 dose, and `0.05` is the Phase 54 high dose. Dynamics did not change with dose. Standard-normal innovations were paired over doses and amplitudes within each seed/model/margin.

All conditions used the same fixed absolute grid from 0 to 75 seconds at 0.1-second intervals (751 samples). This differs from the Phase 56 sketch's `3τ` window: Phase 55's τ-scaled timestamps could reveal the true rate through the sample interval, which is visible to any estimator. The fixed 75-second grid spans 3τ for the slowest included rate and samples the fastest at over eight points per τ. It consequently spans different numbers of τ at other margins; results are reported by model and the five-margin exponent design retains this observation-window heterogeneity.

The two primary methods received only times and observed values:

1. **Free-offset exponential:** Phase 53/55's `b + A exp(−ρt)` fit, estimating its offset and rate. Its initialization scale was derived only from the observed range.
2. **Tail-offset log-linear:** mean of the final 10% of observations as an estimated offset; OLS of log positive residuals against time. Fewer than five positive residuals or a nonnegative slope invalidated the fit. Positive-residual selection is part of this estimator and is recorded; no input rows were pre-cleaned.

The **correct-normal-form diagnostic** fit a named model trajectory and estimated its margin from observations. It received model-family information unavailable to the model-agnostic methods, so it is reported separately and is excluded from primary estimator contrasts and exponent criteria. No method received the true equilibrium, margin, boundary, recovery rate or generator noise SD. The observed-range scale is a numerical initializer, not a noise estimate.

Rate errors were scored against the analytical local rates after each fit. Exponents and the unknown boundary were fitted using five estimated equilibria/rates per seed. Invalid rates and groups remained in attempted denominators. The unchanged Phase 55 p criterion was: at least 80% valid groups, mean absolute p error at most 0.20, and a 95% seed-bootstrap interval covering theoretical p. q is reported separately. Intervals use 2,000 deterministic seed-cluster percentile-bootstrap replicates, are unadjusted across the prespecified family, and are descriptive; non-exclusion of zero does not establish equivalence. No interval or inferential interpretation is assigned to deterministic zero-noise repetitions.

The SNR diagnostic is RMS clean sampled excursion from true equilibrium divided by generator noise SD; it is infinite at zero dose. It uses generator truth, is not an estimator input, and is descriptive rather than causal.

## Evaluation design and reproducibility

Calibration seeds were `58000–58007` and held-out evaluation seeds `58050–58099`; the sets are disjoint. Calibration was limited to procedure and schema checks, with no tuning. There were 960 calibration traces and 6,000 evaluation traces; each trace was sent to all three estimators, for 2,880 and 18,000 fit attempts, respectively. The two primary estimators generated 2,400 evaluation exponent groups, five margins per group.

Run command:

```powershell
.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.noise_dose_estimator_comparison import run_study; print(run_study(Path('reports/phase_57_noise_dose_estimator_comparison')))"
```

The run used Python 3.12.10, NumPy 2.5.3 and SciPy 1.18.1. Prescribed outputs are covered by [`sha256_manifest.json`](../../reports/phase_57_noise_dose_estimator_comparison/sha256_manifest.json). A separate rerun in an isolated temporary directory reproduced every manifested hash. See the result CSVs for all attempts, margin-stratified and aggregated summaries, exponent groups, paired contrasts and SNR.

## Recovery-rate results

All rates and relative absolute errors below aggregate the five margins and held-out seeds. Errors are conditional on valid fits and are accompanied by validity. Full paired bootstrap intervals and margin-level detail are in `rate_condition_summary.csv`, `rate_margin_summary.csv` and `paired_contrasts.csv`.

At the smallest amplitude, the noise-dose pattern was:

| Model | Estimator | η=0 validity / MAE | η=0.0025 validity / MAE | η=0.01 validity / MAE | η=0.05 validity / MAE |
|---|---|---:|---:|---:|---:|
| Fold | Free-offset exponential | 100% / 0.00058 | 100% / 0.920 | 100% / 1.923 | 100% / 2.610 |
| Fold | Tail-offset log-linear | 100% / 0.104 | 66.0% / 0.996 | 54.4% / 0.997 | 50.4% / 0.997 |
| Transcritical | Free-offset exponential | 100% / 0.00133 | 100% / 0.429 | 100% / 4.294 | 100% / 21.951 |
| Transcritical | Tail-offset log-linear | 100% / 0.192 | 90.8% / 0.941 | 58.4% / 0.975 | 45.6% / 0.978 |
| Pitchfork | Free-offset exponential | 100% / 0.00180 | 100% / 0.719 | 100% / 2.882 | 100% / 7.395 |
| Pitchfork | Tail-offset log-linear | 100% / 0.062 | 84.0% / 0.976 | 60.0% / 0.986 | 54.4% / 0.988 |

MAE is mean absolute relative rate error, so `1.0` is 100% error. At amplitude `0.04`, η=`0.01`, free-offset MAE was 0.152 fold, 0.092 transcritical and 0.120 pitchfork; at η=`0.05` it was 1.090, 0.594 and 0.877. Across the full evaluation, the free-offset fit was marked valid in 6,000/6,000 attempts despite very large errors. The tail-offset method was valid in 4,969/6,000; its 1,031 failures were nonnegative or undefined log slopes. A valid optimizer outcome therefore did not guarantee an accurate rate.

The model-informed correct-form fit was valid in 6,000/6,000 attempts. At the smallest amplitude its η=`0.05` mean relative absolute error was 0.00157 fold, 0.00446 transcritical and 0.00459 pitchfork. This is evidence that the named trajectory fit can recover parameters in data generated by that same exact model under these settings. The structural match is built into the diagnostic, so its low error does not validate model selection, mismatched systems, or general superiority.

### Paired estimator and dose contrasts

At `a=0.0025`, the paired dose-by-estimator interaction was supported in the observed direction at both η=`0.01` and `0.05` for each model: the 95% interval for `(tail MAE − exponential MAE at dose) − (same estimator difference at zero)` excluded zero. At η=`0.05`, the estimated differences were −1.605 [−2.747, −0.735] fold (48 paired seeds), −10.563 [−17.979, −4.873] transcritical (45), and −5.181 [−9.550, −2.204] pitchfork (50). Negative values mean the free-offset estimator's error increased more relative to the tail method; they do not mean the tail estimator performed well. The tail method's high invalidity and near-unit relative errors remain visible.

At `a=0.04`, the η=`0.05` interaction was inconclusive for fold (−0.022 [−0.269, 0.198]) and pitchfork (0.072 [−0.194, 0.291]); transcritical was positive, 0.192 [0.069, 0.290]. Noise-vs-zero contrasts and direct tail-minus-exponential comparisons are available for every predeclared model/amplitude/dose cell. There is no universal estimator winner across these conditions.

### SNR diagnostic

At `a=0.0025`, mean generator SNR across the five margins fell from 0.100, 0.238 and 0.169 at η=`0.0025` to 0.0250, 0.0595 and 0.0424 at η=`0.01`, and to 0.0050, 0.0119 and 0.0085 at η=`0.05` for fold, transcritical and pitchfork. At `a=0.04`, corresponding SNR means at η=`0.01` were 0.397, 0.940 and 0.666. These are deterministic generator diagnostics. Their association with errors does not establish a causal mechanism.

## Exponent identifiability — separate outcome

Of 2,400 evaluation exponent groups (50 seeds × 3 models × 2 amplitudes × 4 doses × 2 observation-only estimators), all 1,200 free-offset groups were valid; 777/1,200 tail-offset groups were valid. Group validity required all five rates and equilibria to be usable and the unknown-boundary scaling fit to converge. Invalid groups stayed in each condition denominator.

The unchanged p criterion passed only three of 36 nonzero-noise conditions:

| Model | Estimator | a | η | Valid groups | Mean absolute p error | 95% interval for mean p |
|---|---|---:|---:|---:|---:|---:|
| Fold | Free-offset exponential | 0.04 | 0.0025 | 50/50 | 0.0518 | [0.4814, 0.5139] |
| Transcritical | Free-offset exponential | 0.04 | 0.0025 | 50/50 | 0.0239 | [0.9918, 1.0085] |
| Transcritical | Free-offset exponential | 0.04 | 0.01 | 50/50 | 0.0946 | [0.9749, 1.0426] |

Every other nonzero-dose condition failed at least one of validity, error or interval-coverage requirements. No tail-offset condition and no pitchfork condition passed. The result shows that adequate rate-fit validity is not enough: only a small subset of model/amplitude/dose cells yielded p estimates meeting all frozen criteria. q and estimated-boundary diagnostics are in `exponent_condition_summary.csv` and per-seed `exponent_group_results.csv`.

## What the benchmark establishes and does not establish

**Direct observations:** rate errors and fit-validity changed with dose and estimator in the fixed synthetic design; paired interaction intervals varied by model and amplitude; only three p-identifiability conditions passed; the model-informed exact-form diagnostic had much lower rate error.

**Interpretation:** low perturbation relative to the specified measurement noise is associated with severe rate error in these generators. The dose contrasts quantify results conditional on paired innovations and the specified Gaussian observation process; they do not assign a general causal share to SNR. The tail method illustrates how finite-tail offset bias and positive-residual restrictions can hurt both validity and rate accuracy.

This benchmark does not discover the scaling exponents; the normal forms encode them. It does not establish real-system predictive value, engineering qualification, financial applicability, forecasting ability or profitability. It does not test process noise, drift, missing data, sensor effects, model mismatch, unseen margins, or another noise family. Bootstrap intervals are unadjusted across conditions and zero-noise replicates are deterministic duplicates. The common absolute grid avoids direct rate leakage in timestamp spacing, but creates unequal effective windows in units of τ across margins.

## Next-step recommendation

Do not treat the model-informed fit as a deployable winner. A useful next synthetic question is: **Can an observation-only early-window/local-linear rate estimator reduce the tail-offset bias while retaining better low-SNR behavior than the free-offset fit, when all methods use the same fixed timestamps?** This directly addresses the tail-offset failure observed here and keeps model family and true equilibrium unavailable. It would add evidence about a specific, testable model-agnostic alternative rather than rerun the same noise-dose sweep. No later phase is started by this recommendation.

Phase 32 remains `BLOCKED_AUTHORIZATION`; Phase 52 remains `NO-GO`. No market data or external observations were used.

## Verification record

Using `.venv\Scripts\python.exe` (Python 3.12.10):

- PowerShell `$env:OPENBLAS_NUM_THREADS='1'; .venv\Scripts\python.exe -m pytest tests/test_noise_dose_estimator_comparison.py -q --basetemp .pytest-temp` ? 17 passed.
- PowerShell `$env:OPENBLAS_NUM_THREADS='1'; .venv\Scripts\python.exe -m pytest --basetemp .pytest-temp` ? 570 passed.
- `.venv\Scripts\ruff.exe check .` — passed.
- `.venv\Scripts\ruff.exe format --check .` — passed, 156 files formatted.
- `.venv\Scripts\mypy.exe src` — passed, 50 source files.
- `.venv\Scripts\mypy.exe .` — passed, 95 source files.
- The full frozen experiment was run twice; all nine prescribed file hashes matched. `verify_manifest()` returned true before and after the repeated run.
- Output hashes and final metadata are in `reports/phase_57_noise_dose_estimator_comparison/`.
- Phase 32 and Phase 52 statuses were checked and unchanged. Protected Phase 48–56 artifacts are checked against a pre-edit SHA-256 baseline; the Phase 51 Atlas and roadmap are the only intentionally updated prior-phase files.
- The Phase 51 Atlas JSON and local links are validated after the update. No Git commands or external services were used.
