# Phase 54 — Recovery-Rate Scaling Identifiability Ablation Protocol

**Status: frozen before implementation and evaluation.** The executable copy is [`reports/phase_54_exponent_identifiability_ablation/protocol.json`](../../reports/phase_54_exponent_identifiability_ablation/protocol.json). Any implementation defect requiring a design change must increment the protocol version and trigger complete regeneration.

## Question and estimands

Which observation and estimation components—perturbation amplitude, equilibrium/boundary information, observation noise, and window length—contribute to errors in equilibrium exponent `q` and recovery-rate exponent `p` for the Phase 53 fold, transcritical, and supercritical-pitchfork branches?

For each model and factorial condition, estimate mean `q` and `p`, absolute error from exact Phase 53 values, 95% seed-bootstrap intervals, interval widths and truth coverage, and valid-group probability. A replicate group is one seed/model/factor condition spanning all five controls. Invalid groups stay in denominators. Secondary outcomes are per-trajectory relative recovery-rate error and rate-fit validity by model and margin. Classification is not rerun because Phase 54 studies scaling/rate estimation, not model selection.

## Models and fixed settings

Reuse Phase 53's `stable_equilibrium`, `local_recovery_rate`, `exact_trajectory`, `estimate_exponential`, and `estimate_scaling`. Do not duplicate or modify them. The tested positive branches and exact `(q,p)` are fold `(0.5,0.5)`, transcritical `(1,1)`, and supercritical pitchfork `(0.5,1)`, for positive margins `(0.04, 0.09, 0.16, 0.25, 0.36)` and controls `u=1+mu`, with true boundary `u_c=1` hidden from reference estimates.

The design uses sampling interval `0.1 tau`, with `tau=1/rho` computed from the exact model rate for design only. Windows are `1 tau` (11 samples) and `3 tau` (31 samples). Estimators receive sample times and states, not `tau`, `rho`, or margin except where explicitly supplied by an oracle-information condition. The short window is a deliberate stress case and still provides enough points for the exponential estimator.

## Full factorial and pairing

| Factor | Levels and definition |
|---|---|
| Perturbation fraction | `0.05`, `0.20`, `0.50`; initialize at `x(0)=x*(1+a)`. Record displacement `a*x*` and displacement/margin. These starts remain positive and above the chosen branch; no basin crossing is included. |
| Information | `reference`, `oracle_equilibrium`, `oracle_boundary`, `oracle_both`. Reference estimates exponential offset and shared boundary. Oracle equilibrium fixes only the exponential offset to true `x*`; oracle boundary fixes only `u_c=1` in the scaling fit. Both fixes both values. Oracles are diagnostic upper bounds, not deployable estimators. |
| Relative noise SD | `0`, `0.01`, `0.05` times equilibrium. Independent additive Gaussian observation noise only; zero means exactly noiseless observations. Underlying dynamics do not change. |
| Window | `1 tau`, `3 tau`; sample interval fixed at `0.1 tau`, giving 11/31 points. |

The full `3 x 4 x 3 x 2` factorial has 72 cells per model. It identifies all main effects and preselected interactions: noise-by-window, amplitude-by-noise, and information-by-window. Higher-order interactions are not interpreted. Eight calibration seeds (`54000-54007`) and 20 disjoint evaluation seeds (`55000-55019`) span all models, conditions, and five margins: 1,728 calibration groups, 4,320 evaluation groups, 30,240 trajectories total.

Seeds are paired across factor cells by model and margin. A common standardized Gaussian stream is reused across amplitude, information, noise, and window; the short window uses the prefix of the long-window stream. This matches observation noise while allowing the clean trajectory to change with amplitude. Calibration is a procedural/schema check only: no estimator, threshold, factor, or criterion is tuned from it. Inferential results use evaluation seeds only.

## Estimation and validity

Reference and oracle-boundary rate fits reuse Phase 53's free exponential `a exp(-rho t)+b`. Oracle-equilibrium fixes `b=x*` and fits only positive amplitude and rate by nonlinear least squares. At zero noise, the Phase 53 estimator receives a `1e-12` scale floor only for its positive-scale validation/initialization; observed values remain exact. The fixed-offset fitter does not receive the true rate.

Reference scaling reuses Phase 53 `estimate_scaling`. Oracle-equilibrium substitutes true equilibria while boundary remains estimated. Oracle-boundary performs OLS on log values versus log true margins while equilibrium remains estimated. Oracle-both supplies true equilibria and margins; `q` is then a deterministic construction check, while `p` still depends on estimated rates. A group is valid only if all five rate estimates are finite, positive, and converged and the required scaling fit converges with finite exponents. Invalid groups are counted, never assigned zero error. Intervals bootstrap replicate seeds among valid groups and always accompany validity denominators.

Identifiability passes only if at least 80% of groups are valid, absolute error of both mean exponents is at most the Phase 53 threshold `0.20`, and both 95% percentile intervals cover theory. For oracle-both `q`, report the zero-uncertainty upper-bound result separately; it is not evidence of observation-based identifiability.

## Frozen statistical comparisons

For each model, report all cells and model-specific balanced marginal contrasts: amplitude `0.50-0.05`; each oracle-information level versus reference; noise `0.05-0`; window `3 tau-1 tau`. Preselected difference-in-differences are noise-by-window, amplitude-by-noise, and oracle-both-versus-reference by window. Apply contrasts to valid-group probability, conditional absolute `q` error, conditional absolute `p` error, and mean relative rate error. Use paired replicate-seed bootstrap intervals; exponent-error contrasts include only seeds valid in both paired conditions and report their matched denominator. Never pool models as a substitute for model-specific results.

An effect is called supported only when the paired 95% interval excludes zero and at least half the evaluation seeds are matched. Otherwise report inconclusive or not detected; non-detection is not evidence of no effect. If fewer than two valid groups support a cell interval, report the interval as undefined. Overall identifiability is not declared if any model's reference condition fails. No evaluation-tuned threshold is allowed.

## Limitations frozen in advance

All trajectories are synthetic and generated by the equations under study. Noise is observation-only and Gaussian, scaled to equilibrium; process noise, drift, sensor lag, missingness, irregular sampling, model mismatch, and multimode dynamics are excluded. Oracle conditions are unattainable diagnostic upper bounds. The generator uses true `rho` to define comparable windows, though estimators only receive timestamps and observations. Calibration and evaluation share the five control settings; no unseen-margin extrapolation is tested.
