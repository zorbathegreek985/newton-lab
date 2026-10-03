# Phase 54 — Recovery-Rate Scaling Identifiability Ablation

## Executive summary

This controlled synthetic study varied perturbation amplitude, supplied equilibrium/boundary information, observation noise, and observation-window length across the Phase 53 fold, transcritical, and supercritical-pitchfork models. It generated 30,240 traces in a frozen 3x4x3x2 factorial design with paired innovations and disjoint calibration/evaluation seeds.

The factorial contrasts support a clear within-design pattern: higher observation noise reduced valid exponent-group rates and increased `q` and `p` errors; longer windows improved valid-group rates and reduced exponent errors; larger perturbations improved fitted exponents in this tested range; and oracle equilibrium information improved fit validity and exponent errors. Oracle boundary information improved exponent estimates but had no effect on trajectory-rate fitting, as expected because boundary information enters only at the cross-setting scaling stage. These are effects under the declared generators and estimators, not universal causal claims.

The Phase 53 all-class identifiability criterion did not hold throughout the ablation. Only 3/18 reference-information cells passed for fold, 4/18 for transcritical, and 8/18 for pitchfork. A mid-level reference-design cell (`a=0.20`, noise SD `0.01 x*`, window `3 tau`) passed for all three models; it is an illustrative factorial cell, not a separately predeclared primary condition. Recovery-rate estimation remains a distinct and less stable outcome: valid trajectory fits averaged about 97.9%, but relative-error means were strongly right-skewed, especially for transcritical traces.

This phase identifies factor contrasts and localizes estimator sensitivity; it does not establish that any one factor alone caused Phase 53's failures, discover the encoded scaling laws, or validate a real engineering or financial system. Phase 52 remains `NO-GO`; Phase 32 remains `BLOCKED_AUTHORIZATION`.

## Question, hypotheses, and frozen protocol

**Question:** Which observation and estimation components contribute to Phase 53 exponent-identifiability failures across fold, transcritical, and supercritical-pitchfork branches?

The protocol was frozen before implementation and evaluation in [`phase_54_exponent_identifiability_ablation_protocol.md`](phase_54_exponent_identifiability_ablation_protocol.md) and [`protocol.json`](../../reports/phase_54_exponent_identifiability_ablation/protocol.json). Its SHA-256 is `5b548abfe3ad7cb46dc4e6576daf47d843df4352e97aed507b786c135171cc89`. It specifies no directional hypothesis for factor effects. A factor contrast is marked supported only when its paired-seed bootstrap 95% interval excludes zero and at least half of evaluation seeds match. Intervals containing zero are reported as inconclusive; this design has no equivalence margin with which to claim evidence of no effect.

Phase 53 source was treated as authoritative. It defines, for positive `mu` on the selected stable branches:

| Model | Equation | Stable branch | Exact rate | `(q,p)` |
|---|---|---|---|---|
| Fold | `dx/dt = mu - x^2` | `x*=+sqrt(mu)` | `2 sqrt(mu)` | `(0.5,0.5)` |
| Transcritical | `dx/dt = mu*x - x^2` | `x*=mu` | `mu` | `(1,1)` |
| Supercritical pitchfork | `dx/dt = mu*x - x^3` | `x*=+sqrt(mu)` (the negative branch is also stable) | `2 mu` | `(0.5,1)` |

The bifurcation point is non-hyperbolic and excluded; no other branch or parameter side is included. Phase 54 reuses Phase 53's exact trajectory, branch/rate functions, free exponential estimator, and unknown-boundary scaling estimator. It does not change Phase 53 code or outputs.

## Design and estimator details

The complete factor grid has 72 cells per model:

| Factor | Levels |
|---|---|
| Initial perturbation | `x(0)=x*(1+a)`, `a in {0.05, 0.20, 0.50}` |
| Information supplied | Reference; oracle equilibrium; oracle boundary; both oracles |
| Additive observation-noise SD | `0`, `0.01 x*`, `0.05 x*` |
| Observation duration | `1 tau`, `3 tau`, where `tau=1/rho` is generator-only |

The sampling interval is `0.1 tau` (11 or 31 samples). The five Phase 53 margins and controls remain unchanged. All initial states remain on the positive side and above the selected stable branch; no basin-crossing case was introduced. The raw trajectory table records the absolute displacement, displacement/equilibrium, and displacement/margin. A common standardized Gaussian stream is paired across factor levels within seed/model/margin, and the short trace uses a prefix of the long trace. Noise changes the observations only, not the dynamics.

Eight calibration seeds and 20 disjoint evaluation seeds were used. The calibration partition ran the frozen estimators as a procedural/schema check and did not tune estimators or thresholds. Inference uses evaluation seeds only: 4,320 five-margin exponent groups and 21,600 individual traces, plus the calibration partition, for 30,240 traces total.

Reference rates use Phase 53's free-offset fit `b+a exp(-rho t)`. Oracle-equilibrium conditions fix only `b=x*` and fit amplitude/rate; they do not supply the rate. Oracle-boundary conditions use true margins only in the scaling stage; the trajectory fit remains unchanged. Reference scaling reuses Phase 53's estimated-boundary fit; known-boundary cases regress log values on log true margins. At zero observation noise, observations are exactly the clean trajectory; a `1e-12` numerical scale floor is supplied to the positive-scale Phase 53 estimator input only. No artificial noise is added.

A trajectory rate fit is valid only when the optimizer converges and returns a finite positive rate. A five-margin group is valid only if all its trajectory rate fits and its scaling fit are valid with finite exponents. Invalid groups remain in denominators and receive no exponent score. Each cell's bootstrap interval resamples evaluation seeds among valid groups and reports the valid denominator. The unchanged Phase 53 exponent threshold is absolute error of each mean `q`/`p` at most `0.20` and both 95% intervals covering theory; Phase 54 additionally requires at least 80% valid groups. Oracle-both `q` is a construction upper-bound check, not observation-based identifiability.

Paired-seed bootstrap 95% intervals estimate all main contrasts and three preselected interactions: noise-by-window, amplitude-by-noise, and oracle-both-versus-reference by window. Contrasts are model-specific and equally marginalize balanced factor contexts. Conditional exponent-error comparisons retain paired valid cells and report matched-seed counts. No evaluation labels or outcomes were used for calibration or model selection. Classification was not repeated; it is not the primary or a proxy outcome here.

## Primary exponent-identifiability results

There is no single scalar “overall pass”: the criterion is evaluated for every model/factor cell. Under the reference-information estimator, only a minority of the 18 cells per model passed:

| Model | Reference cells passing / 18 | Mean valid-group rate across reference cells | Mean absolute `q` error across reference cells | Mean absolute `p` error across reference cells |
|---|---:|---:|---:|---:|
| Fold | 3 | 76.4% | 0.309 | 0.157 |
| Transcritical | 4 | 78.3% | 0.121 | 0.138 |
| Supercritical pitchfork | 8 | 78.6% | 0.089 | 0.140 |

For illustration, the middle-amplitude, middle-noise, long-window reference cell had 20/20 valid groups for each model and passed the exponent criterion:

| Model | Mean `q` (95% seed interval) | Mean `p` (95% seed interval) | Mean absolute errors `(q,p)` |
|---|---|---|---|
| Fold | `0.523` (`0.477, 0.570`) | `0.523` (`0.461, 0.595`) | `(0.023, 0.023)` |
| Transcritical | `1.001` (`0.953, 1.050`) | `1.019` (`0.968, 1.069`) | `(0.001, 0.019)` |
| Supercritical pitchfork | `0.495` (`0.467, 0.527`) | `1.004` (`0.948, 1.069`) | `(0.005, 0.004)` |

Across the full design, 3,840 of 4,320 evaluation groups were valid (88.9%). The reference-information cell pass counts show how strongly identifiability depends on observation settings. The information-oracle pass counts are retained in `exponent_condition_summary.csv`; oracle conditions are diagnostics only and must not be read as deployable performance.

## Factor contrasts and interactions

The following are paired marginal contrasts, averaged over the other balanced factors. Values are changes in valid-group probability, absolute exponent error, and mean absolute relative rate error; confidence intervals and exact matched denominators for every model/outcome/contrast are in [`factor_contrasts.csv`](../../reports/phase_54_exponent_identifiability_ablation/factor_contrasts.csv). Negative error contrasts indicate lower error at the named higher/right-hand level.

| Model | Contrast | Valid-group probability | Absolute `q` error | Absolute `p` error |
|---|---|---:|---:|---:|
| Fold | `a=0.50 - 0.05` | +0.146 | -0.137 | -0.490 |
| Fold | noise `0.05 - 0` | -0.271 | +0.220 | +0.662 |
| Fold | `3 tau - 1 tau` | +0.158 | -0.049 | -0.166 |
| Transcritical | `a=0.50 - 0.05` | +0.183 | -0.110 | -0.381 |
| Transcritical | noise `0.05 - 0` | -0.254 | +0.169 | +0.492 |
| Transcritical | `3 tau - 1 tau` | +0.122 | -0.037 | -0.112 |
| Pitchfork | `a=0.50 - 0.05` | +0.200 | -0.100 | -0.438 |
| Pitchfork | noise `0.05 - 0` | -0.246 | +0.152 | +0.580 |
| Pitchfork | `3 tau - 1 tau` | +0.103 | -0.039 | -0.129 |

All listed amplitude, noise, and window contrasts met the frozen supported-effect rule for all three models. The paired bootstrap intervals exclude zero; for example, the noise contrast's valid-group intervals were fold `[-0.300,-0.238]`, transcritical `[-0.288,-0.221]`, and pitchfork `[-0.275,-0.217]`. The long-window effects are consistent with better observability in this design; they do not show that window length alone explains every Phase 53 discrepancy.

Information contrasts were also model-specific. Supplying the equilibrium improved valid-group probability by `0.214-0.236` and reduced both exponent errors with supported intervals for each model. Supplying the boundary left rate-fit validity unchanged (a structural consequence of the estimator split), but reduced `q` error in all models; its `p`-error contrast was supported for fold/transcritical and inconclusive for pitchfork. Supplying both oracles improved exponent errors, but these are upper-bound diagnostics using synthetic truth.

The three preselected interactions were estimated separately for each model. Noise-by-window and amplitude-by-noise contrasts were supported for valid-group probability and both exponent errors in all three models. Oracle-information-by-window contrasts were likewise supported for those outcomes. This indicates that factor effects are not additive over the tested grid. Relative-rate-error interaction results were less stable: several intervals included zero, especially for transcritical and fold, and should be treated as inconclusive. Full intervals are retained in the contrast artifact; no unlisted higher-order interactions are interpreted.

## Recovery-rate estimation (separate secondary outcome)

Across evaluation trajectories, rate-fit validity averaged 97.9% for fold, 97.7% for transcritical, and 98.1% for pitchfork. Rate accuracy did not track exponent accuracy one-for-one. Across condition-level summaries, mean relative rate error was much larger than the median summary for fold (`1.10` vs `0.17`), transcritical (`7.13` vs `0.26`), and pitchfork (`0.83` vs `0.38`), indicating heavy right tails and occasional very poor rate estimates. Noise `0.05 x*` increased mean relative rate error with supported paired intervals in each model; window effects on rate error were supported for pitchfork but inconclusive for fold and transcritical. Oracle equilibrium reduced rate error for transcritical and pitchfork; the fold contrast was inconclusive. Oracle boundary correctly had no direct effect on trajectory-rate fits.

`rate_condition_summary.csv` reports attempted and valid rate fits, mean/median relative rate error, and equilibrium error by model, margin, and all four factors. It preserves rate estimation as its own outcome rather than treating it as exponent estimation. The large means are sensitive to the right-tailed errors; medians and per-condition denominators should be considered alongside them.

## Interpretation, uncertainty, and limitations

Within this synthetic grid, higher observation noise worsened exponent errors and fit validity; longer observation windows and larger tested perturbations improved them. Oracle equilibrium information materially improved estimation, and oracle boundary information particularly reduced equilibrium-exponent error. The paired factorial design and seed-bootstrap intervals support these contrasts for these generators and levels. They do not identify a unique internal mechanism: changing perturbation amplitude also changes the nonlinear transient, and oracle information is intentionally unavailable in practical estimation.

The finite window, finite amplitude, free asymptote, low-margin signal, and noise mechanisms discussed in Phase 53 are consistent with the results, but Phase 54 does not isolate all of them beyond its specified factor interventions. Rate-error tails limit mean-based interpretation. Only 20 evaluation seeds support each cell, and invalid-group bootstrap intervals are conditional on valid groups. Effects with intervals spanning zero are inconclusive, not evidence of no effect. The design uses Gaussian observation noise only, a known generator timescale to make windows comparable, and the same five controls across calibration/evaluation. It excludes process noise, drift, sensor effects, irregular/missing observations, model mismatch, multimode competition, and unseen-margin transfer.

The equations themselves encode the theoretical exponents. Phase 54 is method-capability and identifiability research; it is not law discovery, independent real-system validation, engineering qualification, finance research, or profitability evidence.

## Reproducibility and verification

Run from the workspace root:

```powershell
.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.exponent_identifiability_ablation import run_ablation; run_ablation(Path('reports/phase_54_exponent_identifiability_ablation'))"
```

The study was run twice using the same final code, protocol, and environment. All seven manifested files (protocol, five CSVs, and metadata) were byte-identical between runs. The final manifest SHA-256 is `9a8e584e542aa06b4f023813d82c29ee7f965dc259d7617f5e783b9e51e2694b`.

Phase 32 remains `BLOCKED_AUTHORIZATION`. Phase 52 remains `NO-GO`. The Phase 48-50 reports/artifacts and Phase 53 source, tests, report, and outputs were hash-checked unchanged. No external data or services were used.

## Supported next question

The ablation explains which tested information and observation factors move exponent error and fit validity, but nonlinear finite-amplitude bias and rate-error tails remain intertwined. A justified next methods question is: **At fixed noise, window, and information condition, how does a progressively smaller perturbation trade nonlinear transient bias against signal-to-noise and exponent-fit validity across these same branches?** This would narrow the amplitude mechanism rather than introduce another model family. It remains synthetic and should be pursued only as a focused follow-up; independent application work still requires a qualified target under Phase 52's gates.
