# Phase 58 — Research Thread Closure and Next-Direction Decision

## Decision

**Close the synthetic recovery-rate scaling thread (Phases 48–57) and pause new experiments while consolidating its strongest results.** This is a priority decision, not a claim that recovery-rate estimation, bifurcation scaling, or cross-domain transfer is fully solved. Another adjacent synthetic estimator/noise ablation is not the default next step. The next scientific thread should begin only when it tests a genuinely distinct mechanism or has independent measurements that can challenge the current synthetic conclusions.

The immediate consolidation milestone is bounded: prepare one concise research brief from the extant reports and outputs, with each quantitative statement linked to its derivation or result table, reproducibility command, assumptions, and claim boundary. No new simulation, model, dataset, or dependency is needed for that milestone.

## Scope and evidence reviewed

The workspace contains the Phase 51 Atlas and roadmap; Phase 52 engineering qualification report; Phase 53 report; Phase 54 and 55 reports and protocols; Phase 56 synthesis and evidence matrix; Phase 57 report and protocol; and the Phase 53–57 report directories. The listed Phase 53, 54, 55 and 57 manifests each validated against their listed files. Phase 56's matrix covers Phases 48–55; Phase 57 is the subsequent experiment and is recorded in the current Atlas. The Phase 51 Atlas, roadmap, and this closure report are the decision records updated for this phase. Historical reports, protocols, datasets and generated results were not edited.

The reviewed evidence is local to Newton Lab. Its categories are mathematical derivation, numerical verification, synthetic empirical results, and method diagnostics. No independent real-world measurements validating these recovery or scaling claims were found in the audited project. Phase 56 has an evidence matrix and synthesis report but no separate SHA-256 output manifest; its JSON matrix was parsed during inspection.

## What the thread established

### Phase 48: a specific second-order convergence result

For a uniform nearest-neighbor path with fixed zero leaders, the discrete modal decay rate is

`λₙ(h) = 4 h⁻² sin²(nπh/2)`.

Relative to the fixed-end continuum heat-mode rate `(nπ)²`, this gives the ratio `sinc²(nπh/2)`. For modes 1–3 over the tested path sizes, the finest-four-grid convergence orders were **1.99942, 1.99767, and 1.99476**. All tested discrete rates were below their continuum references. Across 60 ODE integrations, the largest fitted-rate relative error against the exact discrete eigenvalue was `4.14e-11`. This supports second-order spatial convergence and accurate temporal fitting for this deliberately matched uniform path/stencil and fixed-boundary setup. It does not establish the same rates for arbitrary graphs, different update rules, or a measured thermal medium.

### Phases 49–50: different model-specific recovery relationships

Phase 49 found exponent-one recovery-rate scaling for the selected first heat mode and the underdamped nonlinear pendulum envelope over their respective five-point normalized dissipation sweeps; fitted exponents were **1.000000** and **1.000029**. The heat rate came from an exact mode evaluator, not an independent PDE solve. Pendulum rates came from full-sine numerical trajectories near rest. These are separate model-specific derivations and numerical checks, not evidence of one shared physical mechanism.

Phase 50 compared the canonical fold law `xdot=μ−x²`, `μ>0`, on its stable branch `x*=+√μ`, where the local recovery rate is `ρ=2√μ`. On five held-out lower margins, its exponent-1/2 predictor had mean absolute log error **0.004010**, versus **0.918405** for transferring exponent one. The simulated generator was the fold equation, and equilibrium/margin were supplied. This demonstrates discrimination under that encoded synthetic model and design; it is not validation on an independent controller or device.

### Phases 53–55: classification, rate estimation and exponent identification differ

Phase 53 tested fold, transcritical and positive supercritical-pitchfork branches using finite noisy trajectories. In its primary synthetic classification condition, the model-informed classifier labelled 589/600 correctly with 3 false labels and 8 indeterminate cases; its feature baseline labelled 586/600 correctly with no false labels and 14 indeterminate cases. Both classification criteria passed while all-class exponent-identifiability criteria failed. The transcritical exponent fit had only 26/40 valid groups and a biased estimated boundary. Thus successful class prediction did not imply identifiable exponents or reliable recovery-rate parameters.

Phase 54 varied amplitude, information condition, noise and observation window in a frozen factorial design. Reference-information exponent criteria passed only 3/18 fold, 4/18 transcritical and 8/18 pitchfork cells. Higher tested noise generally worsened validity and exponent errors; larger tested perturbations and longer windows generally helped. Oracle-equilibrium and oracle-boundary modes were diagnostic access conditions, not available-data performance. The factorial evidence localizes sensitivity to its tested interventions; it does not uniquely attribute error to any one mechanism.

Phase 55 narrowed the amplitude range. Under its v1.2 design, no sub-5% condition met the p-identifiability criterion. At very small amplitudes noisy rate errors were large, while zero-noise finite-amplitude signatures were much smaller. This was consistent with an important observation-noise contribution plus a possible finite-amplitude/model-fit contribution, but its single nonzero noise level and single estimator could not assign unique causes. Its SNR used generator truth and was descriptive, not a causal explanation.

### Phase 57: additional dose and estimator evidence

Phase 57 compared noise fractions `(0, 0.0025, 0.01, 0.05)` of equilibrium at amplitudes `0.0025` and `0.04`, on the three normal forms. It used paired innovations, 8 calibration seeds and 50 evaluation seeds. To avoid a hidden rate-information channel in timestamps, its final protocol used a fixed absolute grid, which gives different effective windows in units of the true recovery time across margins.

At amplitude `0.0025`, the observation-only free-offset exponential fit was numerically valid in every evaluation attempt, but its mean relative absolute rate error rose from below 0.2% at zero noise to **2.61 fold, 21.95 transcritical, and 7.40 pitchfork** at dose `0.05`. The tail-offset log-linear alternative was valid in **4,969/6,000** evaluation attempts and often had large error when it was valid. Its dose-by-estimator interaction with the exponential method varied by model and amplitude; neither method is a universal winner. Only **3/36** nonzero-noise model/amplitude/method exponent conditions met the unchanged Phase 55 p criterion. Rate validity did not ensure exponent identification.

The correct-normal-form diagnostic had low rate error on its matching generated equations, but it received the declared model family and estimates under that structural assumption. It is an explicitly model-informed diagnostic, not an equal-information result and not evidence that an unknown real system can be assigned the right normal form. The named equations encode their scaling laws. No unique causal decomposition of estimator error, noise, finite-window effects, and finite-amplitude effects follows from these experiments.

## What remains unresolved

The following evidence levels must stay distinct:

- **Demonstrated within a synthetic model:** the Phase 48 discrete spectrum/convergence, Phase 49–50 model-specific rates and predictor comparison, and Phases 53–57 estimator/identifiability behaviors hold under their stated equations, branches, grids, noise processes and criteria.
- **Potentially generalizable, not independently validated:** local Jacobian eigenvalues, modal decay, and sensitivity of estimation to signal size and observation design are useful mathematical structures. Their practical persistence under different mechanisms, measurement processes, parameter mappings and model misspecification is untested here.
- **Unresolved hypotheses:** whether finite-amplitude nonlinear transient bias materially explains rate error beyond the specified noise and fit settings; how results change under process noise, correlated/non-Gaussian noise, drift, sensor dynamics, missingness or irregular sampling; whether another observation-only estimator improves the tested trade-offs; and whether estimated margins/boundaries remain useful on unseen settings.
- **External-measurement-supported application:** none was established by Phases 48–57. No synthetic result demonstrates real-world predictive performance, physical equivalence across domains, engineering qualification, financial usefulness or profitability.

The existing work already explored the canonical branches, observation noise, window and amplitude effects, information access, classifier behavior, estimator forms and p-identifiability to a substantial degree. More adjacent variations within the same exact generators have diminishing value unless they expose a clearly distinct, falsifiable method question. Revisit the synthetic thread only if independent evidence challenges an assumption or an estimator claim, a predeclared failure mode remains consequential, or a new result would change the scientific decision. Avoid tuning more variants to the same finite grid.

## Wider direction review

| Direction | Scientific value | Present feasibility and decision |
|---|---|---|
| **New physical phenomenon** | Could expand Newton Lab beyond recovery scaling if tied to a distinct mechanism and a meaningful falsification test. | No candidate was qualified in the reviewed backlog with an independent target and measurements. Selecting another canonical ODE now risks extending the synthetic sequence without resolving external validity. Defer until a candidate and bounded test are identified. |
| **Cross-domain mathematical connection** | Can clarify when two systems share a mathematical structure and where their physical meanings diverge. | The recovery sequence already contains synthetic cross-domain analogies and model-family comparisons. Without independently specified target evidence, another analogy would not test physical equivalence or application value. Defer unless it yields a new testable prediction with independent target variables. |
| **Independent empirical validation** | Would most directly test whether a mathematical structure survives real measurement, uncertainty and model error. | Not currently executable from workspace evidence. Phase 52 is `NO-GO`: no named qualified engineering system, documented mechanism/margin, measured recovery data, provenance/terms and feasible validation protocol are jointly established. Phase 32 separately remains `BLOCKED_AUTHORIZATION`; this is not a reason to reopen financial-data acquisition. Reconsider only with a concrete system, accessible measurements, usage terms and protocol. |
| **Pause and consolidate** | Improves the reliability and communication of the strongest evidence and makes unsupported extrapolation less likely. | **Selected now.** It is feasible with existing artifacts, needs no new data or model and prevents the project from treating a longer sequence of synthetic ablations as increasing real-world evidence. |

## Bounded next direction and stopping rule

**Immediate direction:** pause new experiments and consolidate Phases 48–57 into a concise, traceable research brief. The first milestone is one brief that (1) links each quantitative claim to its derivation or output table, (2) states the model, branch, assumptions and evidence category, (3) supplies reproduction commands and current manifest status for the key computational results, and (4) presents the unresolved application/authorization gates without implying they are resolved. This milestone is documentation and audit work; it creates no new experiment.

**Consolidation question:** Can a reader reconstruct the strongest Phase 48–57 claims from the cited equations, reports and manifested outputs without relying on undocumented assumptions or treating synthetic results as external validation?

**Comparison and completion criteria:** audit each quantitative claim in the brief against its source equation/report and, where generated, its existing output; require a source link and scope statement for every included claim, no unsupported applied claim, resolving links, and valid listed manifests. If a material number cannot be traced, mark it unresolved and omit it from the summary rather than infer it. The milestone is complete when the brief and audit checklist satisfy those conditions and independent review can follow the evidence chain. There is no calibration/evaluation split because no model or data are fitted.

**Stop/change rule:** stop the consolidation milestone after the evidence chain is complete; do not extend it into new simulation work. If unsupported or irreproducible claims are found, correct the brief or explicitly flag the artifact for a separately authorized repair decision. Reopen scientific experiments only when one of these gates is met: (a) a concrete independently measured candidate system, documented source and permitted-use terms, sufficient observations and a feasible held-out protocol are available; or (b) a distinct physical mechanism is identified with a bounded model and a result that can falsify the proposed correspondence. Negative or inconclusive results should close the proposed question or redirect to the specific failed prerequisite, not automatically spawn another adjacent benchmark.

**Later alternative:** independent engineering validation near a documented stability boundary could be pursued later because it would add evidence unavailable from the synthetic sequence. The blocking prerequisite is a named, accessible, permitted system with documented mechanism, measured variables, sampling, uncertainty and perturbation/recovery observations. A public report or a plausible source page alone does not satisfy this gate.

The other options are deferred for defined reasons: another physical model is not selected until its target and test add information beyond a new normal form; another cross-domain analogy is not selected without an independently specified target prediction; empirical validation is not currently qualified. This is a bounded pause, not abandonment of Newton Lab's scientific goal.

## Preserved status and constraints

Phase 32 remains `BLOCKED_AUTHORIZATION`. Phase 52 remains `NO-GO`. This phase uses no market data, acquires no data, contacts no source, and performs no experiment. It does not claim independent validation, physical equivalence, forecasting ability or profitability.

## Links and inspection record

- [Phase 51 Research Atlas and roadmap](phase_51_research_atlas_and_roadmap.md)
- [Phase 52 engineering qualification](phase_52_engineering_system_qualification.md)
- [Phase 53 identifiability report](phase_53_recovery_scaling_identifiability.md)
- [Phase 54 report and protocol](phase_54_exponent_identifiability_ablation.md), [protocol details](phase_54_exponent_identifiability_ablation_protocol.md)
- [Phase 55 report and protocol](phase_55_sub5_percent_perturbation.md), [protocol details](phase_55_sub5_percent_perturbation_protocol.md)
- [Phase 56 synthesis](phase_56_evidence_synthesis_and_research_decision.md), [evidence matrix](../../reports/phase_56_evidence_synthesis/evidence_matrix.json)
- [Phase 57 report and protocol](phase_57_noise_dose_estimator_comparison.md), [protocol details](phase_57_noise_dose_estimator_comparison_protocol.md)
- [Phase 51 machine-readable Atlas](../../reports/phase_51_research_atlas/research_atlas.json)

The Phase 53, 54, 55 and 57 SHA-256 manifests were present and all listed entries validated during inspection. The Phase 56 matrix was present and contained eight evidence rows (Phases 48–55); Phase 57 is documented in its own report and the current Atlas. No Phase 48 report was missing. All claims above are scoped to the cited local records; the Atlas also contains research areas and decisions beyond this recovery thread.
