"""Tests for the physics taxonomy and equation registry."""

import pytest
from pydantic import ValidationError

from newton_lab.exceptions import RegistryError
from newton_lab.knowledge.examples import build_example_registry
from newton_lab.knowledge.registry import (
    ApplicationMapping,
    ApplicationMappingLevel,
    EquationRecord,
    EquationRegistry,
    EvidenceCategory,
    EvidenceClaim,
    ImplementationStatus,
    MathematicalClassification,
    SourceReference,
    VerificationStatus,
)
from newton_lab.knowledge.taxonomy import (
    PhysicsBranch,
    PhysicsTaxonomy,
    Subfield,
    default_physics_taxonomy,
)


def _minimal_record(**overrides: object) -> EquationRecord:
    values: dict[str, object] = {
        "equation_id": "test_equation",
        "canonical_name": "Test equation",
        "branch_id": "classical_mechanics",
        "subfield_id": "newtonian_mechanics",
        "mathematical_expression": "x = y",
    }
    values.update(overrides)
    return EquationRecord(**values)


def test_default_taxonomy_covers_major_branches_without_claiming_exhaustiveness() -> (
    None
):
    taxonomy = default_physics_taxonomy()
    branch_names = {branch.name for branch in taxonomy.branches}

    assert len(taxonomy.branches) >= 14
    assert "Classical mechanics" in branch_names
    assert "Quantum physics" in branch_names
    assert "Computational and interdisciplinary physics" in branch_names


def test_taxonomy_can_be_extended_and_resolves_subfields_by_branch() -> None:
    taxonomy = default_physics_taxonomy()
    added_branch = PhysicsBranch(
        id="acoustofluidics",
        name="Acoustofluidics",
        subfields=(
            Subfield(id="microfluidic_acoustics", name="Microfluidic acoustics"),
        ),
    )
    extended = PhysicsTaxonomy(branches=(*taxonomy.branches, added_branch))

    assert extended.branch("acoustofluidics") == added_branch
    assert extended.contains_subfield("acoustofluidics", "microfluidic_acoustics")
    assert not extended.contains_subfield(
        "classical_mechanics", "microfluidic_acoustics"
    )
    custom_record = _minimal_record(
        equation_id="custom_branch_equation",
        branch_id="acoustofluidics",
        subfield_id="microfluidic_acoustics",
    )
    registry = EquationRegistry(taxonomy=extended, equations=(custom_record,))
    assert registry.get("custom_branch_equation") == custom_record


def test_taxonomy_rejects_duplicate_identifiers() -> None:
    branch = PhysicsBranch(id="same_branch", name="First")

    with pytest.raises(ValidationError, match="branch identifiers must be unique"):
        PhysicsTaxonomy(branches=(branch, branch.model_copy(update={"name": "Second"})))


def test_record_supports_multiple_classifications_and_optional_metadata() -> None:
    record = _minimal_record(
        mathematical_classifications=(
            MathematicalClassification.ODE,
            MathematicalClassification.NONLINEAR,
        ),
        assumptions=("A stated idealization.",),
    )

    assert len(record.mathematical_classifications) == 2
    assert record.source_references == ()
    assert record.boundary_conditions == ()
    assert record.assumptions == ("A stated idealization.",)


def test_record_rejects_invalid_stable_identifier_and_missing_required_fields() -> None:
    with pytest.raises(ValidationError):
        _minimal_record(equation_id="Equation 01")

    with pytest.raises(ValidationError):
        EquationRecord.model_validate(
            {"equation_id": "missing_expression", "canonical_name": "Incomplete"}
        )


def test_established_claim_requires_provenance() -> None:
    with pytest.raises(ValidationError, match="require a source reference"):
        EvidenceClaim(
            category=EvidenceCategory.ESTABLISHED_PHYSICAL_RESULT,
            statement="Established claim without traceable source.",
        )


def test_established_application_field_rejects_untested_hypothesis() -> None:
    with pytest.raises(ValidationError, match="must use established physical-result"):
        _minimal_record(
            established_applications=(
                {
                    "category": "plausible_untested_application_hypothesis",
                    "statement": "A proposal that has not been tested.",
                },
            )
        )


def test_hypothesis_and_analogy_are_distinct_and_mapping_is_grade_checked() -> None:
    hypothesis = EvidenceClaim(
        category=EvidenceCategory.APPLICATION_HYPOTHESIS,
        statement="A possible future application to an untested system.",
    )
    mapping = ApplicationMapping(
        domain="engineering",
        application="control model comparison",
        level=ApplicationMappingLevel.HYPOTHESIS,
        evidence=hypothesis,
    )
    assert mapping.evidence.category == EvidenceCategory.APPLICATION_HYPOTHESIS

    with pytest.raises(ValidationError, match="requires evidence category"):
        ApplicationMapping(
            domain="engineering",
            application="structural comparison",
            level=ApplicationMappingLevel.STRUCTURAL_SIMILARITY,
            evidence=hypothesis,
        )


def test_registry_checks_branch_and_subfield_relationships() -> None:
    record = _minimal_record(subfield_id="special_relativity")

    with pytest.raises(ValidationError, match="outside branch"):
        EquationRegistry(equations=(record,))


def test_registry_rejects_duplicate_equation_ids() -> None:
    record = _minimal_record()

    with pytest.raises(ValidationError, match="equation IDs must be unique"):
        EquationRegistry(equations=(record, record))


def test_registry_rejects_invalid_equation_references() -> None:
    record = _minimal_record(
        related_equations=(
            {
                "target_equation_id": "not_in_registry",
                "relationship": "derived_from",
                "evidence": {
                    "category": "mathematically_derived_relationship",
                    "statement": "A derivation relationship.",
                },
            },
        )
    )

    with pytest.raises(ValidationError, match="unknown equation"):
        EquationRegistry(equations=(record,))


def test_registry_rejects_claims_with_unresolved_source_references() -> None:
    record = _minimal_record(
        physical_interpretations=(
            {
                "category": "mathematically_derived_relationship",
                "statement": "A derivation.",
                "source_reference_ids": ("missing_source",),
            },
        )
    )

    with pytest.raises(ValidationError, match="unknown source IDs"):
        EquationRegistry(equations=(record,))


def test_registry_round_trip_serialization_preserves_records() -> None:
    registry = build_example_registry()
    restored = EquationRegistry.model_validate_json(registry.model_dump_json())

    assert restored == registry
    assert restored.get("maxwell_equations_vacuum_si").source_references[0].doi == (
        "10.1098/rstl.1865.0008"
    )


def test_registry_lookup_search_and_immutable_registration() -> None:
    registry = build_example_registry()
    original_count = len(registry.equations)
    new_record = _minimal_record(
        equation_id="additional_test_record",
        canonical_name="Additional test record",
    )
    extended = registry.register(new_record)

    assert len(registry.equations) == original_count
    assert extended.get("additional_test_record") == new_record
    assert registry.search("Schrödinger")
    assert registry.search("not present") == ()
    with pytest.raises(RegistryError, match="duplicate equation ID"):
        extended.register(new_record)
    with pytest.raises(RegistryError, match="unknown equation ID"):
        registry.get("not_present")


def test_source_reference_requires_stable_id_and_title() -> None:
    with pytest.raises(ValidationError):
        SourceReference(source_id="", title="A source")


def test_implemented_simulators_are_linked_but_not_scientifically_verified() -> None:
    registry = build_example_registry()

    for equation_id in ("damped_harmonic_oscillator", "damped_nonlinear_pendulum"):
        record = registry.get(equation_id)
        assert record.implementation_status == ImplementationStatus.TESTED
        assert record.verification_status == VerificationStatus.UNASSESSED
        assert record.related_equations[0].target_equation_id == "newton_second_law"
