# Persistent discovery artifacts

`DiscoveryArtifactStore` saves completed Phase 18 and Phase 19 case studies as
versioned JSON. The default directory is `.newton_lab/artifacts` under the
current working directory; saving creates it when needed. Tests and explicit
workflows can pass another directory.

## Format and identity

Schema version 1 uses an envelope with `schema_version`, `artifact_id`,
`artifact_type`, a complete structured `payload`, and a SHA-256 content
integrity digest. The artifact ID is the original case-study ID and determines
the JSON filename. The digest detects payload changes; it does not validate
scientific claims. JSON output is stable, rejects non-finite numbers, and
preserves numeric precision and tuple order (including repeated sweep values
and failed sweep points).

The two case-study classes have distinct payload schemas. All validated model
fields are retained, including experiments and run outcomes, parameters and
units, analytical values, sampled observations, assumptions, evidence
classifications, limitations, and embedded report text. Loading validates the
envelope and the corresponding concrete case-study model. Unsupported schema
versions, unknown types, duplicate JSON properties, invalid payloads, and
identity conflicts raise `PersistenceError`.
No artifact timestamp is created because the source case-study models do not
carry a trustworthy creation time.

Phase 19's source BVP model contains Python callbacks. Artifacts do not contain
functions or executable expressions. For the built-in
`steady_linear_heat_conduction` problem, the store records the two explicit
callback IDs `steady_heat_ode_v1` and `fixed_endpoint_temperatures_v1`; loading
rebinds them from Newton Lab's built-in factory. Other callback-bearing
problems are rejected. No experiment or solver runs during loading.

## Save, load, and integrate

```python
from newton_lab.cross_domain import run_damping_transfer_case_study
from newton_lab.heat_diffusion_case_study import (
    run_heat_diffusion_smoothing_case_study,
)
from newton_lab.knowledge.artifacts import (
    DiscoveryArtifactStore,
    build_knowledge_integration_from_artifacts,
)

store = DiscoveryArtifactStore()  # .newton_lab/artifacts

# Initial case creation runs the original APIs once.
oscillator = run_damping_transfer_case_study()
heat = run_heat_diffusion_smoothing_case_study()
store.save(oscillator)
store.save(heat)

# These steps only read, validate, integrate, and render stored evidence.
loaded_oscillator = store.load(oscillator.case_study_id)
loaded_heat = store.load(heat.case_study_id)
index = build_knowledge_integration_from_artifacts(store)
markdown = index.report_markdown
```

`store.list_artifacts()` returns IDs in stable order. Existing artifacts are
never silently replaced; pass `overwrite=True` to `save` to explicitly replace
one. Writes use a temporary sibling file and atomic creation or replacement to
avoid exposing partial JSON. Artifact IDs are restricted to lowercase letters,
digits, and underscores, and resolved paths must remain inside the selected
artifact directory.

The Phase 20 integration API receives the loaded Phase 18 and Phase 19 model
objects and builds its existing registry links, retrospective associations,
relationship graph, and Markdown report. If either required artifact is
missing or invalid, integration raises an error instead of silently producing
a partial report. Retrospective equation associations remain separate from
original experiment metadata. Evidence classifications and limitations are
preserved: the oscillator's financial hypothesis is unvalidated, the shared
oscillator/diffusion edge remains a candidate, and synthetic smoothing is not
evidence of practical prediction benefit.

Artifact persistence preserves recorded data; it does not establish the
validity of the underlying models or scientific claims. This format is local
JSON rather than a migration or revision-management system. Schema changes
require an explicit supported version or migration in a later implementation.
