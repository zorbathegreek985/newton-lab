# Phase 55 â€” Sub-5% Perturbation Study Protocol

**Protocol version:** 1.2
**Status:** Amended and frozen before the final held-out evaluation.
The executable protocol is [`protocol.json`](../../reports/phase_55_sub5_percent_perturbation/protocol.json). Version 1.0 and 1.1 outputs are superseded and excluded from conclusions. Amendments are documented below; all final results are regenerated under version 1.2.

### Version 1.1 amendment

The initial v1.0 run revealed that the zero-noise arm passed a `1e-12` scale floor to the Phase 53 exponential fitter, while the noisy arm passed `0.01*x*`. Because that argument controls the fit's starting-amplitude initialization, the diagnostic had changed estimator setup as well as observation noise. Version 1.1 fixes the estimator's supplied scale at the primary `0.01*x*` value in both arms. In the zero-noise arm, no random noise is added to observations; only the generated observations change. This preserves the same estimator initialization and thresholds across arms. The diagnostic remains noncausal because the deterministic zero-noise arm cannot represent sampling uncertainty.

Version 1.2 records that SNR is deterministic for each model/margin/amplitude setting, not a stochastic estimate over replicate seeds. Its amplitude trend is reported exactly and descriptively, with no bootstrap interval or inferential-effect label. This changes no design factor or success criterion.

## Question and hypotheses

Across the Phase 53 fold, transcritical, and supercritical-pitchfork normal forms, how does reducing perturbation amplitude below 5% affect recovery-rate and recovery-rate-exponent estimation, and what evidence distinguishes nonlinear transient bias from declining signal-to-noise (SNR)?

- **H1, nonlinear transient bias:** finite perturbations can make the exponential recovery fit systematically amplitude-dependent. Support would include amplitude-dependent signed rate bias or normalized fit residual in the zero-noise arm as well as in the noisy arm. This is consistent with, but does not prove, nonlinear-dynamics bias.
- **H2, declining SNR:** at fixed nonzero noise, smaller perturbations reduce the clean recovery signal relative to noise and can increase uncertainty/error or reduce fit validity. Support would include worsening noisy-arm outcomes as the predeclared SNR falls and a larger noisy-versus-zero-noise penalty at smaller amplitudes.
- **Indeterminate:** if intervals are too wide, fits are invalid, or the noisy/zero-noise contrast does not isolate the explanations, report that the mechanisms were not distinguished.

The exact Phase 53 equations and positive stable branches are reused for `mu > 0`: fold `xdot=mu-x^2`, `x*=sqrt(mu)`, `rho=2sqrt(mu)`, `(q,p)=(1/2,1/2)`; transcritical `xdot=mu*x-x^2`, `x*=mu`, `rho=mu`, `(q,p)=(1,1)`; supercritical pitchfork `xdot=mu*x-x^3`, selected branch `x*=sqrt(mu)`, `rho=2mu`, `(q,p)=(1/2,1)`. The bifurcation point and other branches/sides are excluded.

## Frozen settings

Margins are Phase 53's `(0.04, 0.09, 0.16, 0.25, 0.36)`, with controls `u=1+mu` and true boundary 1. The amplitude fraction is `a=delta/x*`, exactly Phase 54's normalization; the initial state is `x(0)=x*(1+a)`. Evaluate `a=(0.0025, 0.005, 0.01, 0.02, 0.04)` as the sub-5% range. Retain `0.05` and `0.20` as previously tested Phase 54 references. There are no post-result amplitude additions. Every model/margin/amplitude displacement is `delta=a*x*`; also record `delta/mu`. Because `x*` differs by normal form, equal `a` is a relative-to-equilibrium comparison, not equal displacement or equal displacement-to-margin across models. All starts lie above the selected positive equilibrium and remain on its positive basin side.

The primary arm fixes additive, independent Gaussian observation noise at `sigma=0.01*x*`, the Phase 54 middle noise level. The fixed window is `3 tau` with `dt=0.1 tau` (31 samples), where `tau=1/rho` is generator-only as in Phase 54. The information condition is Phase 54 `reference`: the free-offset exponential estimates equilibrium and rate, and the existing shared-boundary scaling fit estimates q and p. Estimators receive only sampled times and observations plus the declared measurement-scale input; they do not receive true equilibrium, margin, rate, or boundary as model parameters. A paired zero-noise diagnostic changes only the generated observation noise to zero. To hold estimator setup exactly fixed, the same `0.01*x*` scale input is passed to the fitter in both arms; it is used only as the existing initialization scale. No artificial noise or `1e-12` floor is added to zero-noise observations. It is a diagnostic contrast, not a second primary condition.

Calibration seeds are `56000â€“56007`; held-out evaluation seeds are `57000â€“57049`. These sets are disjoint. Calibration is limited to deterministic execution/schema/edge-case checks and cannot tune estimators, levels, thresholds, or criteria. Evaluation results are not inspected before the implementation, focused tests, and protocol checks are complete. Within each model/margin/seed, standardized Gaussian innovations are reused across all amplitudes and both noise arms; the 31-point sample times are identical. Thus amplitude/noise comparisons are paired within a model, while cross-model trend contrasts use independent seed bootstrap samples.

## Outcomes, validity, and uncertainty

Rate fitting reuses Phase 53 `estimate_exponential`. A rate fit is valid only if the optimizer converges and the rate is finite and positive. The exact theoretical rate is the reference. Preserve every attempt; invalid fits stay in denominators and are never assigned zero error. Report rate estimates and signed/absolute errors by model, margin, amplitude, and arm, plus valid counts, distribution summaries, fitted RMSE, and RMSE divided by the true initial displacement.

The predeclared SNR is `RMS_t[x_clean(t)-x*] / sigma`, evaluated over the sampled window using generator truth solely as a diagnostic. At zero noise report SNR as undefined/infinite (no measurement noise), not as a finite value. This SNR summarizes signal magnitude relative to the specified observation noise; it is not an estimator input.

For each seed/model/amplitude/arm, fit Phase 54's unknown-boundary scaling model across the five margins. A scaling group is valid only if all five rate fits are valid and the scaling optimizer converges with finite q and p. Primary exponent is p; q is a prespecified equilibrium-scaling diagnostic. Report attempted/valid groups, means, bias and absolute error, seed-bootstrap 95% intervals, interval width, theory coverage, and the Phase 54 joint criterion. Also report a p-specific criterion: at least 80% valid groups, mean absolute p error no greater than 0.20, and its 95% interval covers theory. Neither criterion turns invalid groups into successes.

Use 2,000 deterministic seed-bootstrap replicates for cell intervals and paired contrasts. Rate and exponent error summaries are conditional on valid fits and always accompanied by validity denominators. Paired noisy-minus-zero contrasts use seeds valid in both arms and state matched denominators; validity contrasts use all paired seeds. For trend summaries, regress seed-level mean outcomes on `log10(a)` over the five sub-5% levels only. Rate validity and group validity slopes use every seed. Conditional error/width slopes use seeds with at least three valid amplitude levels and report how many contribute. Bootstrap seeds, not individual margins or amplitudes, to preserve repeated-measure dependence. Cross-model differences in slopes use independent model-specific seed resampling. Report SNR exactly from the generator and noise configuration; do not assign it a sampling interval or inferential effect status. No significance is interpreted as equivalence or proof of no effect.

An amplitude effect is called supported only if its two-sided 95% seed-bootstrap interval excludes zero; otherwise it is inconclusive/not detected. An identifiability success is the cell criterion above, not a hypothesis test. Do not pool models in place of model-specific results.

## Planned comparisons and stopping rules

1. Within each normal form, estimate sub-5% log-amplitude trends in absolute relative rate error, signed relative rate bias, valid rate probability, normalized fit RMSE, SNR, absolute p error, p-validity, p-interval width, and coverage.
2. Report the same outcomes at the `0.05` and `0.20` Phase 54 reference levels without including them in the sub-5% trend fit.
3. At each amplitude, compare noisy and zero-noise arms with paired seed intervals for rate error, valid rate probability, exponent error, exponent validity, and interval width. Report coverage separately; do not subtract arms and label the difference a causal decomposition.
4. Compare model-specific sub-5% trend slopes using independent bootstrap samples.
5. Treat zero-noise amplitude trends in rate bias and normalized residual as evidence consistent with finite-amplitude/model-fit effects. Treat noisy-only degradation aligned with falling SNR and a larger noisy-vs-zero penalty at low amplitude as evidence consistent with measurement-noise limitation. These signatures do not uniquely identify causes.

Stop or label a cell/trend indeterminate when it has fewer than 2 valid exponent groups for an interval, fewer than 25 contributing seeds for a conditional trend, or optimizer/serialization failure. Do not change thresholds or levels after evaluation. If the two signatures overlap or intervals fail to resolve them, recommend a follow-up with an independently varied noise scale or a directly specified nonlinear-vs-local-linear fit comparison; do not claim either mechanism was isolated.

## Limitations

The equations encode the theoretical exponents. This is a synthetic method-capability study, not discovery of scaling laws or engineering, financial, or real-world validation. Gaussian observation noise is the only stochastic observation process; the dynamics themselves are deterministic. The true-rate timescale defines comparable windows. Only five fixed margins and positive branches are tested. Equilibrium-relative perturbation fractions do not equate absolute amplitudes across models. Deterministic zero-noise seed replicates do not represent sampling uncertainty; the estimator receives the fixed primary noise-scale input solely for its initialization in both arms. No process noise, drift, sensor effects, missing/irregular sampling, model mismatch, branch crossing, or unseen-margin transfer is tested.
