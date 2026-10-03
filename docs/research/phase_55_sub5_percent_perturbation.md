# Phase 55 — Sub-5% Perturbation Study

## Question and scope

This controlled synthetic study asks whether poor recovery-rate estimation at very small perturbations is more consistent with finite-amplitude transient/model-fit effects, declining observation signal-to-noise, both, or remains indeterminate. It reuses the Phase 53 positive stable branches for the fold (`dx/dt = mu - x^2`, `x*=sqrt(mu)`, `rho=2sqrt(mu)`), transcritical (`dx/dt = mu*x - x^2`, `x*=mu`, `rho=mu`), and selected positive supercritical-pitchfork branch (`dx/dt = mu*x - x^3`, `x*=sqrt(mu)`, `rho=2mu`), all at `mu>0`. Their predicted recovery-rate exponent pairs are respectively `(q,p)=(1/2,1/2)`, `(1,1)`, and `(1/2,1)`.

The equations encode the scaling laws. This is an identifiability and estimator-capability benchmark, not discovery of those laws or validation of a physical, engineering, financial, or market process.

## Frozen design and estimation

Protocol v1.2 was frozen before the held-out evaluation; its executable specification is [`protocol.json`](../../reports/phase_55_sub5_percent_perturbation/protocol.json), and its rationale/amendments are in the [protocol document](phase_55_sub5_percent_perturbation_protocol.md). The estimator sees sampled times, observations, and the declared measurement-scale initialization input; it is not supplied the true equilibrium, margin, rate, or bifurcation boundary. The same estimator and initialization scale (`0.01*x*`) are used in both arms. The diagnostic arm changes only generated observation noise to zero.

The primary amplitudes are `a=delta/x*` = 0.0025, 0.005, 0.01, 0.02, and 0.04, with Phase 54 reference levels 0.05 and 0.20. The noisy arm has independent Gaussian observation noise with standard deviation `0.01*x*`; the paired diagnostic arm has no observation noise. There are five Phase 53 margins, a `3 tau` window, 31 samples at `0.1 tau`, and paired standardized innovations across amplitudes and arms. Eight calibration seeds (56000–56007) were reserved for procedural/schema checks only. Fifty disjoint evaluation seeds (57000–57049) were used, with 2,000 deterministic seed-bootstrap replicates for inferential summaries. No thresholds or criteria were tuned on evaluation results.

Rate estimates come from the existing free-offset exponential fit. The shared-boundary scaling fit estimates equilibrium and rate scaling parameters `q,p` from the five margins. Invalid attempts remain in denominators. The primary p-identifiability rule requires at least 80% valid groups, absolute mean p error at most 0.20, and an interval covering the theoretical p. The zero-noise arm is descriptive and is not scored as an inferential criterion. SNR is a deterministic generator diagnostic; its amplitude trend is reported exactly, without a bootstrap interval or inferential status.

## Run and artifacts

Run from the repository root with the project environment:

```powershell
.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.sub5_percent_perturbation import run_study; print(run_study(Path('reports/phase_55_sub5_percent_perturbation')))"
```

The run used Python 3.12.10, NumPy 2.5.3, and SciPy 1.18.1. It generated 12,180 trajectories total: 1,680 calibration and 10,500 held-out evaluation trajectories. Evaluation contained 2,100 exponent groups. Of evaluation trajectories, 10,348/10,500 (98.55%) had valid rate fits; 1,717/2,100 (81.76%) exponent groups were valid. CSV outputs retain individual attempts, groups, condition summaries, paired contrasts, model contrasts, and the exact amplitude specification. [`sha256_manifest.json`](../../reports/phase_55_sub5_percent_perturbation/sha256_manifest.json) covers the generated artifacts other than itself.

## Findings

### Rate estimation

Under the noisy arm, sub-5% rate error decreased with increasing amplitude in all three models. The seed-bootstrap slope of mean absolute relative rate error against `log10(a)` was fold −1.314 (95% interval −1.753 to −0.933), transcritical −20.697 (−27.879 to −14.297), and pitchfork −6.027 (−8.532 to −3.983). Normalized fit RMSE also fell with amplitude in all models. At the representative margin `mu=0.16`, noisy-arm mean absolute relative rate error at `a=0.0025` versus `a=0.04` was 2.137 versus 0.458 for fold, 17.698 versus 0.382 for transcritical, and 2.005 versus 0.388 for pitchfork. These are large estimator errors at very low signal, including estimates far above the true rate.

In the zero-noise diagnostic, the exact deterministic amplitude trend in signed/absolute rate error was small and positive: 0.0113, 0.0225, and 0.0338 per log10-amplitude unit for fold, transcritical, and pitchfork. Normalized residual slopes were 0.000221, 0.000440, and 0.000667. These small deterministic changes are consistent with a finite-amplitude/model-fit contribution, but do not establish that nonlinear transient bias explains the much larger noisy errors. Zero-noise repeated seeds are identical generator observations, not independent uncertainty replicates.

The deterministic noisy-arm SNR slope per `log10(a)` was 1.253 (fold), 1.245 (transcritical), and 1.238 (pitchfork). At `mu=0.16`, SNR rose from about 0.105 at `a=0.0025` to 1.66–1.68 at `a=0.04`. SNR contains generator truth by design and is descriptive; it is not an estimated predictor nor a causal test.

### Exponent identifiability is distinct from rate accuracy

No sub-5% noisy amplitude met the p-specific identifiability criterion for any model. At `a=0.04`, valid groups were 45/50 fold, 43/50 transcritical, and 47/50 pitchfork, but the mean absolute p errors were 0.231, 0.218, and 0.376, and the p intervals missed theory in all three. At `a=0.20`, fold and transcritical met the p-specific criterion; pitchfork did not because its interval failed to cover p=1. This reference result does not change the failed sub-5% conclusion.

The exponent estimates can remain biased or uncertain even when most individual rates fit successfully. Conversely, rate validity alone does not establish a correctly identified scaling exponent. The per-amplitude table reports p validity, error, interval, coverage, and criterion separately.

### Paired noisy versus zero-noise diagnostics

At each of the three smallest reference points shown (`a=0.0025`, 0.01, 0.04), the noisy-minus-zero-noise mean absolute relative rate-error contrast was positive and its paired seed-bootstrap interval excluded zero for every model. At `a=0.0025`, the contrasts were 2.029 [1.581, 2.562] fold, 22.957 [16.289, 30.511] transcritical, and 7.545 [5.146, 10.336] pitchfork. The contrast remained positive but smaller at `a=0.04`: 0.447 [0.383, 0.530], 0.420 [0.380, 0.463], and 0.381 [0.340, 0.424]. This supports a strong observation-noise/SNR contribution under the specified generator and estimator. It is not a causal decomposition: zero-noise fits are deterministic and the two arms cannot represent unmodeled measurement or process effects.

## Assessment of the competing explanations

The evidence is **more consistent with a substantial low-SNR/measurement-noise contribution, alongside a small finite-amplitude/model-fit contribution**. The noisy rate errors and normalized residuals worsen sharply at small amplitudes, deterministic SNR falls with amplitude, and paired noisy-minus-zero contrasts are large at the smallest level. Meanwhile the zero-noise rate-error and residual trends are nonzero but small. The evidence does not establish unique mechanisms or quantify a general causal share, so conclusions remain conditional on these synthetic models, noise, estimator, margins, and observation window.

The primary research question is answered only for these generator settings. There is no evidence here about real-world systems, forecasting, financial applications, or profitability.

## Limitations and unresolved tests

Only the three canonical normal forms, their selected positive stable branches, five margins, one nonzero noise fraction, one sampling/window condition, and one observation noise family were studied. The time scale is set from true recovery rates; measurement noise is Gaussian and dynamics contain no process noise. The equilibrium-relative amplitude is not the same absolute displacement across models. The experiment does not vary noise magnitude independently, compare the exponential estimator to a directly specified nonlinear/local-linear fit, or test model mismatch, irregular/missing samples, drift, sensor response, branch transfer, or unseen margins. Those are plausible follow-ups, not explanations established by this run.

Phase 32 remains `BLOCKED_AUTHORIZATION`; Phase 52 remains `NO-GO`. This study used no market data and makes no applied claims.

## Verification record

Using `.venv` Python 3.12.10:

- `.venv\Scripts\python.exe -m pytest tests/test_sub5_percent_perturbation.py -q --basetemp .pytest-temp` — 15 passed.
- `$env:OPENBLAS_NUM_THREADS='1'; .venv\Scripts\python.exe -m pytest --basetemp .pytest-temp` — 553 passed.
- `.venv\Scripts\ruff.exe check .` — passed.
- `.venv\Scripts\ruff.exe format --check .` — passed after formatting the edited source helper (151 files already formatted).
- `.venv\Scripts\mypy.exe src` — passed (49 source files).
- `.venv\Scripts\mypy.exe .` — passed (93 source files).

The frozen v1.2 run was repeated; both runs produced the same SHA-256 manifest bytes. The final manifest contains 10 artifact hashes and validates against every covered file; its SHA-256 is `ee0ff92333866200e0b1a1f914a9000747e775b120ab81ebcdf8e49e1127903e`. Metadata's protocol hash matches `protocol.json`. The protected Phase 48–54 baseline had 52 files; all remained byte-identical. Phase 32 remained `BLOCKED_AUTHORIZATION` and Phase 52 remained `NO-GO`. Phase 51 atlas JSON parsed, and internal links in the Phase 51 and Phase 55 reports/protocol resolved. No Git commands were run.
