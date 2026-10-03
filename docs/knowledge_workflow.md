# Knowledge ingestion, review, and persistence

Newton Lab's knowledge workflow keeps submitted records separate from the
accepted equation registry. It accepts manually supplied structured mappings,
records already present in the project example registry, and one local JSON
object at a time. It does not fetch web pages, crawl sources, or execute an
equation expression.

## Equation and source records

An `EquationRecord` describes one equation or model using a stable ID, its
mathematical expression, assumptions, taxonomy placement, related records,
evidence claims, and source references. A source can identify authors or a
responsible organization, publication date and edition, source type, publisher,
DOI, URL, ISBN or other identifier, a page or section locator, and access date.
Reliability notes and limitations capture provenance concerns. A source ID
must refer to the same metadata wherever it is reused in a workspace.

These fields document provenance; they do not certify source quality. Claims
retain their evidence category, and software tests do not establish scientific
validity.

## Ingest candidates

```python
from newton_lab.knowledge import candidate_from_json_file, candidate_from_mapping

candidate = candidate_from_mapping(
    {
        "equation_id": "linear_example",
        "canonical_name": "Linear example",
        "branch_id": "classical_mechanics",
        "mathematical_expression": "y = a x",
    }
)
candidate_from_disk = candidate_from_json_file("candidate.json")
```

Manual mappings are validated against the equation schema. The local JSON
importer requires a single object, rejects duplicate keys and malformed JSON,
and retains the submitted text with the candidate. Invalid submissions raise
`CandidateIngestionError` with actionable validation detail and, when
available, the raw submission. Ingested project examples also begin as drafts;
existing records are not silently accepted.

## Human review

The states and permitted transitions are:

```text
draft -> needs_review -> accepted -> superseded
             |
             +-> rejected -> draft
                    +------> superseded
             +-> draft
```

`needs_review` can return to draft, be rejected, or be accepted. Rejected
records can return to draft for revision or be superseded. Accepted records can
only be superseded; superseded records are terminal. Every transition records
an actor, timezone-aware timestamp, notes, a snapshot of the record, and any
review checklist, source IDs, and edit description.

Acceptance requires at least one source, inspection of every attached source,
review of the mathematical form and assumptions, a reviewer identity, and
notes. Records with empirical or independent-validation claims also require an
explicit empirical evidence assessment. This is a human editorial gate, not a
change to `verification_status` and not a claim that the result has been
scientifically proven.

```python
from newton_lab.knowledge import KnowledgeBase, ReviewChecklist, ReviewState

workspace = KnowledgeBase().add_candidate(candidate)
workspace = workspace.transition_entry(
    candidate.record.equation_id,
    ReviewState.NEEDS_REVIEW,
    actor_id="submitter:researcher",
)
# Acceptance is possible after a reviewer has inspected every referenced source
# and completed all required review dimensions.
```

Only accepted entries appear in `accepted_registry()`. That method validates
the usual taxonomy, source, and relationship constraints, including that a
related record is also accepted.

## JSON persistence

`JsonKnowledgeStore` stores the taxonomy, candidates, source metadata, original
submissions, current states, and complete review history in a JSON envelope
with `schema_version: 1`. Writes are deterministic and use a temporary file
plus atomic creation/replacement. Existing files are preserved unless
`save(..., overwrite=True)` is explicitly requested. Load rejects malformed,
duplicate-key, invalid, or unsupported-version documents. No automatic
migration is attempted: unsupported versions need an explicit migration before
they can be loaded.

```python
from newton_lab.knowledge import JsonKnowledgeStore

store = JsonKnowledgeStore("research/knowledge.json")
store.save(workspace)  # refuses to overwrite an existing file
restored_workspace = store.load()
```

## Numerical and scientific limits

Ingestion checks structure and known references, not the truth of an equation.
Review relies on people inspecting the cited material. A source may be
incomplete, inaccessible, or interpreted incorrectly. The registry stores
equation expressions as text and does not parse or solve them. Persistence
validates against its declared schema but is not a substitute for backups or
independent review. Future schema changes require explicit versioned
migrations.
