# Phase 53 — Recovery-Rate Scaling Identifiability Across Canonical Bifurcations

## Executive summary

This phase tested whether finite noisy recovery trajectories from three canonical one-dimensional bifurcations can be identified when equilibrium and distance to the bifurcation are not supplied to the estimators. The frozen benchmark used 3,600 synthetic trajectories, independent calibration/evaluation noise seeds, two noise levels, and two observation windows.

Under the predeclared low-noise, long-window condition (`sigma=0.0005`, window 8), both direct model-informed curve selection and a separate model-agnostic calibrated feature classifier passed their classification criteria. Their results apply to seed variation on the same five control settings. The model-agnostic classifier uses calibration labels and local equilibrium/rate features; it does not establish that fitted scaling exponents alone classify unseen conditions.

The central scaling-identifiability criterion **failed** for all three systems. Although equilibrium and margin recovery from fitting the correct named equation were accurate in the primary condition, free exponential estimates and their fitted exponents showed bias and nonconvergence, especially for low-margin transcritical trajectories. For example, only 26/40 transcritical seed groups produced a valid scaling fit in the primary condition; the equilibrium exponent estimate was `1.228` (95% seed-bootstrap interval `[1.078, 1.364]`) instead of 1. The pitchfork rate exponent was `0.935` (interval `[0.923, 0.950]`) instead of 1. Thus a classifier can separate these generated cases while the quantitative scaling estimand remains poorly identified.

This is a synthetic method-capability study. The equations encode the laws, and no real system, financial data, or application was tested. Phase 32 remains `BLOCKED_AUTHORIZATION`.

## Research question and frozen hypotheses

**Question:** When equilibrium and bifurcation distance must be estimated from finite, noisy trajectories, can recovery-rate scaling distinguish the fold from canonical bifurcations with linear recovery-rate scaling?

The protocol was written to [`reports/phase_53_recovery_scaling/protocol.json`](../../reports/phase_53_recovery_scaling/protocol.json) before implementation and evaluation. It fixes the tested branches, margins, settings, noise, windows, independent seed partitions, fitting methods, rejection thresholds, and success criteria. The protocol SHA-256 is `9af29103a3ad1c7efe69d6daad9c5c775d2903131488cd64434b222c6b1fcbce`.

## Equations, branches, and theoretical signatures

All tested trajectories use positive margin `mu>0` and the stable positive branch. The bifurcation point `mu=0` is not treated as an ordinary stable case.

| Normal form | Equilibria and stability | Tested stable branch; positive recovery rate | Scaling signature |
|---|---|---|---|
| Fold: `dx/dt = mu - x^2` | For `mu>0`, equilibria are `x=+/-sqrt(mu)`. The positive branch has Jacobian `-2sqrt(mu)` and is stable; the negative branch is unstable. For `mu<0` there are no real equilibria. | `x*=+sqrt(mu)`, `rho=2sqrt(mu)`. | Equilibrium exponent `q=1/2`; rate exponent `p=1/2`. |
| Transcritical: `dx/dt = mu*x - x^2` | Equilibria are `x=0` and `x=mu`. Their Jacobians are `mu` and `-mu`, respectively. For `mu>0`, zero is unstable and `x=mu` is stable; for `mu<0`, zero is stable and `x=mu` is unstable. At zero, the exchange point is non-hyperbolic. | `x*=mu`, `rho=mu` for `mu>0`. | `q=1`; `p=1`. |
| Supercritical pitchfork: `dx/dt = mu*x - x^3` | For `mu<0`, zero is stable with eigenvalue `mu`; for `mu>0`, zero is unstable and branches `x=+/-sqrt(mu)` are stable, each with eigenvalue `-2mu`. The origin at zero is non-hyperbolic. | Positive branch `x*=+sqrt(mu)`, `rho=2mu`; the negative branch has the same rate by symmetry. | `q=1/2`; `p=1`. |

Here `rho` is the positive decay rate `-lambda` on the selected stable branch. Both `q` and `p` use the independently observed control setting's distance from an unknown boundary, `u-u_c`; the fitted prefactors absorb the normalized control scale. The transcritical and pitchfork systems both have linear rate exponents, but their equilibrium exponents differ. Fold and pitchfork share a square-root branch exponent but have different rate exponents. These differences are branch- and side-specific; no exponent is asserted for other branches or parameter sides.

## Data-generating design and estimators

The exact equations above generate all records. The latent positive margins are `{0.04, 0.09, 0.16, 0.25, 0.36}` and the observed control settings are `u=1+mu`. The bifurcation setting `u_c=1` and each true margin are used by the generator and evaluator only. Neither estimator receives them. It receives the observed control setting, time samples, measured state values, and declared observation-noise standard deviation. Neither receives true equilibrium or true initial state.

Each trajectory starts at `x(0)=1.2 x*` on the chosen stable branch and is sampled at interval 0.1 over windows 3 and 8 in the benchmark's normalized time coordinate. Additive independent Gaussian observation noise has standard deviation `0.0005` or `0.002` in normalized state units. There is no process noise, drift, missingness, parameter noise, sensor lag, or additional mode. Closed-form trajectories provide truth without numerical integration error.

Calibration uses seeds 41000–41019; evaluation uses disjoint seeds 51000–51039. The five settings occur in both partitions, so this evaluates independent noise-seed generalization at a fixed setting grid, **not extrapolation to unseen margins**. Calibration labels are known because this is a supervised synthetic benchmark.

Two routes are reported separately:

1. **Model-informed trajectory fitting:** for each trace, fit each named model's exact solution with its positive margin unknown; use the first observed sample as the fitted trajectory's initial state. Select by known-noise Gaussian BIC only when the best-to-runner-up gap is at least 6; otherwise reject as indeterminate. The fitted candidate margin determines that candidate's equilibrium and Jacobian rate. The true-model fit is also scored for parameter recovery. The estimator is not passed the true branch point, margin, or equilibrium.
2. **Model-agnostic estimation and classification:** fit each trace to `b + a exp(-rho t)` with free positive amplitude and rate and free offset; `b` estimates equilibrium and `rho` estimates rate. For each seed group spanning the five observed settings, fit equilibrium and rate exponents with a shared unknown boundary. Separately, a diagonal-Gaussian nearest-centroid classifier uses calibration-set labels and log estimated `(b,rho)`, stratified by observed setting, noise, and window. It receives no candidate equation or true parameter values. It rejects when the top posterior is below 0.60 or the log-score gap is below `ln(2)`.

The second classifier is a calibrated feature baseline, not an exponent-only classifier. It is allowed to use calibration labels at the same five settings as evaluation. That choice makes the classification results conditional on this known grid and synthetic calibration design.

The reusable functions and study runner are in [`src/newton_lab/recovery_scaling_identifiability.py`](../../src/newton_lab/recovery_scaling_identifiability.py); focused tests are in [`tests/test_recovery_scaling_identifiability.py`](../../tests/test_recovery_scaling_identifiability.py). The `_s` suffix in the internal time-setting names is bookkeeping only: these canonical states/time are normalized and the values are not measurements in physical seconds.

## Frozen success criteria

The primary condition is `sigma=0.0005`, window 8.0. Before evaluation, the protocol required: (a) model-informed three-class accuracy among non-rejected trajectories at least 0.70 with reject rate at most 0.25; (b) model-agnostic three-class balanced accuracy at least 0.70, reject rate at most 0.25, and fold sensitivity and linear-family specificity each at least 0.75; (c) for each model, median absolute error of mean `q` and `p` estimates at most 0.20 and the 95% seed-bootstrap interval containing both theoretical exponents; and (d) median relative error in model-informed margin and equilibrium estimates at most 0.25 at each tested margin. These thresholds were not changed after seeing evaluation results.

## Results

### Classification and abstention across noise and windows

Each noise/window cell contains 600 evaluation trajectories (3 models × 40 seeds × 5 margins). “Accuracy, all” counts rejects as not correct; false classifications count incorrect non-rejected labels.

| Method | Noise SD | Window | Correct / 600 | False labels | Indeterminate | Accuracy, all | Accuracy among classified |
|---|---:|---:|---:|---:|---:|---:|---:|
| Model-informed | 0.0005 | 3 | 459 | 19 | 122 | 76.5% | 96.0% |
| Model-informed | 0.0005 | 8 | 589 | 3 | 8 | 98.2% | 99.5% |
| Model-informed | 0.002 | 3 | 266 | 24 | 310 | 44.3% | 91.7% |
| Model-informed | 0.002 | 8 | 503 | 34 | 63 | 83.8% | 93.7% |
| Model-agnostic feature baseline | 0.0005 | 3 | 571 | 2 | 27 | 95.2% | 99.7% |
| Model-agnostic feature baseline | 0.0005 | 8 | 586 | 0 | 14 | 97.7% | 100% |
| Model-agnostic feature baseline | 0.002 | 3 | 545 | 5 | 50 | 90.8% | 99.1% |
| Model-agnostic feature baseline | 0.002 | 8 | 578 | 3 | 19 | 96.3% | 99.5% |

At the primary condition, model-informed confusion totals by true class were: fold `200/200` correct; transcritical `189` correct, `3` classified as fold, and `8` rejected; pitchfork `200/200` correct. The feature baseline had fold `200/200` correct, transcritical `186` correct and `14` rejected, and pitchfork `200/200` correct. No false classifications occurred for the feature baseline in this cell.

The distance-stratified confusion and abstention counts are retained in [`classification_confusion.csv`](../../reports/phase_53_recovery_scaling/classification_confusion.csv); rates and denominators are in [`classification_metrics.csv`](../../reports/phase_53_recovery_scaling/classification_metrics.csv). At primary noise/window, model-informed accuracy across the 120 trajectories per margin ranged from 95.0% to 100% (reject rate 0–4.2%); the feature baseline ranged from 88.3% to 100% (reject rate 0–11.7%). At the highest noise and shortest window, model-informed abstention rose to 51.7%, while its accuracy among the remaining classified trajectories was 91.7%. Abstention is therefore material and prevents the high conditional accuracy from hiding low coverage.

Both frozen classification criteria passed in the primary condition. This does not show exponent-only classification: the feature baseline learned class-conditional local equilibrium/rate distributions from labeled calibration seeds at the same settings.

### Scaling exponents and uncertainty

The model-agnostic power fit estimates both branch exponent `q` and recovery exponent `p`, plus the shared boundary. Intervals below are 95% percentile bootstrap intervals over evaluation seed groups; for transcritical rows, invalid/nonconverged groups are excluded and their count is shown explicitly.

| Model | Valid groups / 40 | Mean estimated boundary `u_c` (truth 1) | `q` mean [95% interval], abs mean error (truth) | `p` mean [95% interval], abs mean error (truth) | Exponent criterion |
|---|---:|---|---|---|
| Fold | 40 | 0.9992 | 0.503 [0.500, 0.506], 0.0031 (0.5) | 0.496 [0.493, 0.498], 0.0043 (0.5) | **Fail:** rate interval misses 0.5 |
| Transcritical | 26 | 0.9689 | 1.228 [1.078, 1.364], 0.2282 (1) | 0.994 [0.946, 1.049], 0.0060 (1) | **Fail:** equilibrium exponent error exceeds 0.20; only 26 fits valid |
| Supercritical pitchfork | 40 | 0.9947 | 0.514 [0.503, 0.526], 0.0143 (0.5) | 0.935 [0.923, 0.950], 0.0645 (1) | **Fail:** rate interval misses 1 |

The fold branch exponent is recovered within the frozen error limit, but its rate interval misses the theoretical value. For the pitchfork, the branch exponent is close while the rate interval is below the theoretical value. For the transcritical model, its linear rate exponent is estimated reasonably at this primary condition, while the estimated equilibrium exponent and boundary are biased; just 65% of groups yielded a usable full power fit. The frozen criterion requires both absolute error of each class's mean exponent to be at most 0.20 and both intervals to cover truth. The transcritical equilibrium error exceeds the threshold, and fold/pitchfork intervals miss their rate exponents; therefore the all-class scaling criterion fails.

Across harder settings, accuracy worsened as expected for short and noisy observations. For example, transcritical valid scaling groups ranged from 10/40 at `sigma=0.002`, window 3 to 26/40 at primary settings. The full per-class/per-condition estimates, bootstrap intervals, log-fit residuals and estimated boundaries are in [`scaling_exponent_summary.csv`](../../reports/phase_53_recovery_scaling/scaling_exponent_summary.csv) and [`scaling_exponent_estimates.csv`](../../reports/phase_53_recovery_scaling/scaling_exponent_estimates.csv). The interval describes seed variability conditional on this generator; it is not uncertainty across physical systems.

### Parameter recovery diagnostics

In the primary condition, the true-model-informed fit's median relative errors for margin and branch equilibrium were at most 2.4% across model-by-margin cells; the largest was transcritical at the smallest margin `mu=0.04`. This is conditional on fitting the correct candidate equation and scoring against synthetic truth. It does not mean the model was known in a real application.

The free exponential estimator's largest per-cell median equilibrium error was 10.6% and its largest median rate error was 98.8%, both for low-margin transcritical trajectories. The latter reflects a rate too slow to resolve reliably over the finite window with a small offset/amplitude and noise. This diagnostic explains why a correct class label from the feature classifier does not imply unbiased recovery-rate or exponent estimation. All per-margin diagnostics are in [`parameter_recovery_summary.csv`](../../reports/phase_53_recovery_scaling/parameter_recovery_summary.csv), and individual estimates/standard errors are in [`per_trajectory_estimates.csv`](../../reports/phase_53_recovery_scaling/per_trajectory_estimates.csv).

### Frozen criteria outcome

- Model-informed classification: **passed** (99.5% accuracy among classified; 1.3% rejected).
- Model-agnostic feature classification: **passed** (100% balanced accuracy among classified; 2.3% rejected; fold sensitivity and linear-family specificity 100%).
- Scaling identifiability across all three models: **failed** (none of the three class-specific exponent checks passed all error and interval-coverage conditions).
- Model-informed parameter recovery for the correct candidate: **passed** at all 15 primary model-by-margin cells.
- Overall method-success criterion: **failed**, because exponent identifiability is a required component.

## Interpretation and limitations

For this fixed synthetic setting grid, class-conditioned local features and exact candidate-trajectory fitting often distinguish the generated systems at low noise and longer observation duration. The strict quantitative scaling task is harder: transient curvature, free-asymptote estimation, finite perturbation size, low-margin signal amplitude, and observation noise affect both `b` and `rho`. The particularly poor transcritical low-margin rate estimate demonstrates that finite trajectories can make equilibrium and recovery jointly weakly identified.

The models have known forms in this benchmark, and the full candidate equations were used both to generate data and to define model-informed fits. This is controlled method-capability testing, not independent validation of the normal forms. The model-agnostic classifier is trained with synthetic class labels at the same settings as evaluation; it does not establish performance on unseen settings. There is no process noise, operating-point drift, measurement lag, higher-dimensional mode competition, parameter uncertainty, or model misspecification. The study does not establish a universal bifurcation diagnostic or real-world predictive value.

## Reproducibility and verification

Reproduce from the workspace root with:

```powershell
.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.recovery_scaling_identifiability import run_identifiability_study; run_identifiability_study(Path('reports/phase_53_recovery_scaling'))"
```

The report artifacts are `protocol.json`, six CSV outputs, `metadata.json`, and `sha256_manifest.json`. During verification, the implementation's first criterion calculation was found to use median per-seed error where the frozen protocol specifies error of the mean exponent. The calculation was corrected to match the unchanged protocol and the benchmark was regenerated; no criterion or threshold was changed. The corrected benchmark was run twice with the same protocol. The manifest SHA-256 before and after the repeat was identical: `355F28AEB47A14B77B77244BD74620DA5B0E5B0AFFF0092012266617B5E476CE`.

The focused command `.venv\Scripts\python.exe -m pytest tests\test_recovery_scaling_identifiability.py -q` passed (18 tests). The full test command `.venv\Scripts\python.exe -m pytest --basetemp .pytest-temp` passed (522 tests, 533.31 seconds). `.venv\Scripts\ruff.exe check .`, `.venv\Scripts\ruff.exe format --check .`, `.venv\Scripts\mypy.exe src`, and `.venv\Scripts\mypy.exe .` all passed. Hash comparisons confirmed Phases 48–52 reports, source, tests, and report artifacts were unchanged. Phase 32 remains `BLOCKED_AUTHORIZATION`; no market data or external service was used.
