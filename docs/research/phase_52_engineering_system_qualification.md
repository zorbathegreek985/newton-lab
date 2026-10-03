# Phase 52 — Real Engineering System Qualification

**Decision: NO-GO.** No candidate in the inspected workspace currently passes the qualification gates. The strongest candidate concept is a measured electrical or process-control system near a documented stability limit, but no named plant, measurement dataset, access terms, or disturbance/recovery records are present. The available synthetic fold study cannot substitute for an independent target.

## Executive summary and connection to Phase 51

[Phase 51](phase_51_research_atlas_and_roadmap.md) identified the lack of independently measured application evidence as Newton Lab's main research gap. This phase assessed whether the current evidence can support an engineering test without first implementing another model. The answer is no: the local workspace contains the Phase 48–50 synthetic/numerical study artifacts and two user-supplied NSE report files, but no engineering telemetry or laboratory observations.

The assessed hypothesis is the Phase 50 fold-specific relationship: near an independently established saddle-node stability boundary, the local recovery rate is predicted to scale with the square root of the independently defined positive operating margin. Phase 50 derives this result for its canonical model `dx/dt = mu - x^2`, where the stable equilibrium is `sqrt(mu)` and the positive local rate is `2 sqrt(mu)`. Its target equation was chosen to embody that fold relationship. It tests the estimator and predictor under a controlled synthetic process; it is not evidence from a real plant.

Phase 49 offers a secondary exponent-one baseline from the selected heat and pendulum models, but their exponents arise from their own parameterizations. Phase 50's local margin and equilibrium are supplied exactly. Neither study tests whether those quantities can be inferred independently from partially observed engineering measurements.

## Scope and inspected evidence

Verified local paths were inspected before use:

- [Phase 51 atlas and roadmap](phase_51_research_atlas_and_roadmap.md) and [machine-readable atlas](../../reports/phase_51_research_atlas/research_atlas.json).
- [Phase 49 report](phase_49_critical_slowing_down.md), [source](../../src/newton_lab/critical_slowing_down.py), and its focused tests.
- [Phase 50 report](phase_50_cross_domain_validation.md), [source](../../src/newton_lab/fold_control_validation.py), focused tests, frozen protocol, metadata, and result files.
- [Phase 32 authorization gate](phase_32_study_readiness.md), [Phase 33 source audit](phase_33_data_source_audit.md), [Phase 37 acquisition plan](phase_37_historical_data_acquisition_plan.md), [Phase 39 file audit](phase_39_historical_dataset_audit.md), [Phase 40 coverage/schema plan](phase_40_historical_coverage_and_schema_plan.md), and [Phase 41 authorization checkpoint](phase_41_authorization_evidence_checkpoint.md).
- The workspace `data/` file inventory. It contains the previously documented NSE ZIP/CSVs; no engineering-system data file was found.

No external research tools, websites, data downloads, or services were accessed in Phase 52. The official source pages mentioned in Phases 33 and 37 are referenced here only through those existing local reports; they were not rechecked for this qualification. No external engineering dataset, source title, URL, version, access term, or measurement is asserted as independently verified in this phase.

## Candidate systems assessed

The candidates below are target-system concepts, not claims that a particular plant or dataset has been located. “Source evidence” means what the inspected project records establish; a candidate with no identified dataset is not treated as publicly available or authorized.

| Candidate and possible existing hypothesis | Measurements and source evidence actually available | Equilibrium and margin testability | Disturbance, recovery, and independence | Access, baseline, and main limitations | Gate result |
|---|---|---|---|---|---|
| **Electrical power-system voltage-control equipment near a stability limit**; assess the Phase 50 fold-rate hypothesis only if an independent source establishes a fold and defines an operating margin. | Phase 33/37 discuss CEA/CERC monthly market-monitoring publications, but describe them as market reports/summaries, not plant-level voltage, controller-state, perturbation, or recovery telemetry. No specific system or time series is identified locally. | No local observation of the operating point, state variables, or stability margin. A load/voltage margin cannot be assumed from a report label or chosen parameter. An independent engineering definition and uncertainty estimate are required. | No disturbance-and-return trace is available. It is unknown whether the target data arose from measurements, simulation, or analytic generation because no candidate dataset has been identified. | Applicable data rights and output terms for a suitable engineering dataset are not documented. Potential baseline: constant local rate and exponent-one law, fit only on calibration operating points, against the independently derived fold prediction using identical observations. Risks include multidimensional modes, protection/control actions, changing inputs, and non-fold limits. | **FAIL**: no system, data, terms, or response records. |
| **Instrumented industrial process-control plant near a documented loss of operating point**; test the same fold-rate hypothesis if the plant's mechanism is established independently. | No named plant, source, report, historian export, variables, sensor descriptions, or sampling specification was found in the workspace. Phase 50's canonical scalar normal form is a simulator, not a plant source. | Not assessable. Need operating inputs and state measurements across stable operating margins, an equilibrium estimation procedure independent of recovery outcomes, and evidence that the parameter measures distance to the relevant boundary. | No empirical disturbance/recovery sequence is available. A simulator based on `mu-x^2` would reproduce the tested law and count only as simulation validation. | Access, research-use permission, retention and redistribution terms are absent. Same-information baselines could include constant-rate and exponent-one predictors trained on the same high-margin records. Risks include unobserved controls, saturation, delays, drift, and unsafe or unrepresentative disturbances. | **FAIL**: target identity, measurements, permission, and mechanism absent. |
| **Measured mechanical/electromechanical oscillator with adjustable damping**; test the Phase 49 envelope-rate-versus-damping relationship on an actual device. | Newton Lab's pendulum/oscillator APIs and reports are equation-based; no instrumented apparatus, device records, sensor data, calibration records, or data source is present. | Equilibrium could potentially be estimated from quiet baseline observations, and damping could potentially be independently varied/measured, but no such records or validated mapping from a control setting to stability distance exist here. | No measured ring-down traces or sampling/duration details are available. A controlled physical ring-down could be independent of the simulated Phase 49 trajectory if its observations were separately collected and documented. | No source/access terms are available because no dataset is identified. Baseline could be a direct local exponential envelope fit and a model-specific eigenvalue prediction on equivalent data. Oscillation frequency must be kept separate from envelope decay; sensor noise, forcing, nonlinear amplitude and changing damping are risks. | **FAIL**: a plausible experiment design, but no actual system or evidence to qualify. |

The most promising *concept* for the selected fold hypothesis is the first or second candidate, provided its stability mechanism is independently documented and suitable response measurements exist. The mechanical candidate is a more direct route to a real measurement of Phase 49's damping law, but it would not test the Phase 50 fold law unless its actual stability boundary were independently shown to be a fold. No preference between an actual plant and dataset can be made from the current evidence.

## Testability and qualification gates

The Phase 50 hypothesis is stated precisely as follows: for a target whose independently documented local normal form is a saddle-node fold, the positive local recovery rate should scale as `rho = C mu^(1/2)` on the stable side near the boundary, with the exponent predicted from that target's Jacobian. The positive margin `mu` must be defined in the target's physical variables and estimated without using held-out recovery rates. If the target does not have a supported fold mechanism, the hypothesis does not apply.

| Qualification gate | Workspace evidence | Status |
|---|---|---|
| Physical identity and operating context documented | No named plant, component, apparatus, or operating regime identified. | **NOT MET** |
| Measurement provenance and relevant variable definitions documented | No engineering measurements or engineering data source in the workspace. Existing CEA/CERC descriptions concern market-monitoring reports, not the required plant telemetry. | **NOT MET** |
| Access and intended-use terms sufficiently clear | No engineering dataset is identified, so no applicable access, research-use, retention, or redistribution terms can be assessed. | **NOT MET** |
| Sampling, temporal resolution, duration, units, and measurement procedure known | No engineering response record exists locally. | **NOT MET** |
| Equilibrium/operating point estimable with uncertainty | No state observations or estimation protocol available. A theoretical equilibrium supplied by a simulator is not an independent measurement. | **NOT MET** |
| Margin to a relevant boundary independently definable and estimable | No physical parameter map, boundary evidence, or independent margin estimate exists. A convenient setpoint is not automatically distance to a stability boundary. | **NOT MET** |
| Disturbance and recovery measured independently | No empirical perturbation/recovery trace is available. | **NOT MET** |
| Recovery rate estimable without imposing the exponent or law | Feasible as a future method if raw response traces and a prespecified observation window exist; not testable on current data. For oscillatory systems, fit envelope decay separately from frequency. | **NOT MET** |
| Target data independent of the law being tested | Phase 50 is a canonical simulation whose equation encodes the fold law. No independent measured target data identified. | **NOT MET** |
| Fair baseline and falsifiable held-out evaluation possible | Phase 50 provides a useful protocol pattern, but no engineering records exist for equivalent-information calibration and evaluation. | **NOT MET** |

**Qualification outcome: no candidate passes. Final status: `NO-GO`.** This is a decision about evidence currently available in the inspected workspace, not a conclusion that suitable engineering datasets cannot exist elsewhere.

## Exact missing requirements and next action

Before implementation or experiment design, identify one actual candidate and provide or locate evidence for all of the following:

1. Named physical system, owner/operator, operating context, and independent description of its stability boundary/mechanism.
2. Dataset/source identity, original measurement provenance, acquisition method, version/date, and whether records are empirical measurements rather than simulation output.
3. Written terms or permission covering local research, storage/retention, derived outputs, and any intended publication or sharing.
4. Field dictionary with units, sensor locations/accuracy, timestamps/time zone, sampling cadence, missingness, calibration, and operating-condition/control inputs.
5. Multiple stable operating conditions approaching the boundary, with equilibrium estimable from observations and an independently defined margin with uncertainty.
6. Disturbance and subsequent recovery observations of sufficient duration and temporal resolution. Disturbances must not be chosen using the held-out response.
7. A preregistered calibration/test split by operating condition, rate estimator that does not impose the exponent, target-specific prediction, and same-information baselines (at minimum a constant-rate baseline and the exponent-one transfer baseline where meaningful).
8. A stop rule for unsupported fold mechanism, unresolvable margin, inadequate signal-to-noise or cadence, absent permission, or failure to distinguish the proposed relation from plausible alternatives.

The immediate next action is a **target-readiness dossier** for one named system, containing the documents and data dictionary above before any model implementation. If no independently documented and permitted system can be identified, stop the engineering application effort and record the access/measurement gap. Do not substitute Phase 50's simulated trajectories or the NSE cross-sectional files.

## Authorization status, files, and verification

Phase 32 remains **`BLOCKED_AUTHORIZATION`**. Nothing in this phase supplies relevant permission evidence for financial data or changes that gate. No dataset was downloaded, copied, opened for analysis, or uploaded; no external source was accessed; no model or experiment was implemented or run; and no prior phase finding was changed.

Created for Phase 52:

- This report.
- [`candidate_register.json`](../../reports/phase_52_engineering_qualification/candidate_register.json), which records the three candidates and their unmet gates.

Verification checked local path existence, report links, candidate-register JSON parsing/path references, the Phase 32 status, and hashes of the Phase 48–50 reports, source files, tests, and artifact directories. Tests were not run because this phase changes documentation and JSON only. No Git commands were used.
