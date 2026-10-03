"""Typed equation records and a validated, serializable knowledge registry."""

from __future__ import annotations

from datetime import date
from enum import StrEnum
from typing import Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

from newton_lab.exceptions import RegistryError
from newton_lab.knowledge.taxonomy import PhysicsTaxonomy, default_physics_taxonomy

_IDENTIFIER_PATTERN = r"^[a-z0-9]+(?:_[a-z0-9]+)*$"


class MathematicalClassification(StrEnum):
    """Mathematical forms or properties an equation may have simultaneously."""

    ALGEBRAIC = "algebraic"
    TRANSCENDENTAL = "transcendental"
    ODE = "ordinary_differential_equation"
    PDE = "partial_differential_equation"
    INTEGRAL_EQUATION = "integral_equation"
    DIFFERENCE_EQUATION = "difference_or_recurrence_equation"
    STOCHASTIC_EQUATION = "stochastic_equation"
    VARIATIONAL = "variational_formulation"
    OPTIMIZATION = "optimization_formulation"
    CONSERVATION_LAW = "conservation_law"
    EIGENVALUE_PROBLEM = "spectral_or_eigenvalue_problem"
    TENSOR_FORMULATION = "tensor_formulation"
    VECTOR_FORMULATION = "vector_formulation"
    OPERATOR_FORMULATION = "operator_formulation"
    FIELD_FORMULATION = "field_formulation"
    NETWORK_FORMULATION = "network_or_graph_model"
    LINEAR = "linear_system"
    NONLINEAR = "nonlinear_system"
    DETERMINISTIC = "deterministic_system"
    STOCHASTIC = "stochastic_system"


class EvidenceCategory(StrEnum):
    """Evidence status for a scientific or application claim."""

    ESTABLISHED_PHYSICAL_RESULT = "established_physical_result_or_application"
    MATHEMATICALLY_DERIVED = "mathematically_derived_relationship"
    KNOWN_ANALOGY = "known_analogy_or_structural_similarity"
    APPLICATION_HYPOTHESIS = "plausible_untested_application_hypothesis"
    EMPIRICALLY_INVESTIGATED = "empirically_investigated_application"
    INDEPENDENTLY_VALIDATED = "independently_validated_result"


class VerificationStatus(StrEnum):
    """Review status for a record's scientific source and evidence trail."""

    UNASSESSED = "unassessed"
    SOURCES_REVIEWED = "sources_reviewed"
    INDEPENDENTLY_REVIEWED = "independently_reviewed"
    EXPERIMENTALLY_VALIDATED = "experimentally_validated"


class SourceType(StrEnum):
    """Kinds of source material that can support a knowledge claim."""

    PRIMARY_RESEARCH = "primary_research"
    TEXTBOOK = "textbook"
    REVIEW_ARTICLE = "review_article"
    DATASET = "dataset"
    TECHNICAL_STANDARD = "technical_standard"
    ARCHIVAL_SOURCE = "archival_source"
    SOFTWARE_DOCUMENTATION = "software_documentation"
    OTHER = "other"


class ImplementationStatus(StrEnum):
    """Software implementation status, distinct from scientific verification."""

    NOT_IMPLEMENTED = "not_implemented"
    SPECIFIED = "specified"
    IMPLEMENTED = "implemented"
    TESTED = "tested"


class RelationshipType(StrEnum):
    """Supported kinds of equation-to-equation links."""

    DERIVED_FROM = "derived_from"
    LIMIT_OF = "limit_of"
    SPECIAL_CASE_OF = "special_case_of"
    COUPLED_WITH = "coupled_with"
    SHARES_MATHEMATICAL_STRUCTURE = "shares_mathematical_structure"
    SHARES_ASSUMPTIONS = "shares_assumptions"
    SHARES_SYMMETRY = "shares_symmetry"
    SHARES_CONSERVATION_LAW = "shares_conservation_law"
    SHARES_VARIABLES = "shares_variables"
    SHARES_SOLUTION_METHOD = "shares_solution_method"
    APPROXIMATION_OF = "approximation_of"
    EQUIVALENT_UNDER_ASSUMPTIONS = "equivalent_under_stated_assumptions"
    TRANSFORMABLE_INTO = "transformable_into"
    SHARES_MATHEMATICAL_OPERATOR = "shares_mathematical_operator"
    ANALOGY_WITH = "analogy_with"
    CANDIDATE_RELATIONSHIP = "candidate_relationship_requiring_investigation"


class RegistryModel(BaseModel):
    """Common strict, immutable settings for registry records."""

    model_config = ConfigDict(frozen=True, extra="forbid", str_strip_whitespace=True)


class Linearity(StrEnum):
    """Known linearity of an equation in its dependent variables."""

    LINEAR = "linear"
    NONLINEAR = "nonlinear"


class MathematicalStructure(RegistryModel):
    """Explicitly assessed features; ``None`` means unknown or unassessed.

    Empty tuples represent an assessed absence. Features are metadata only and
    do not form or parse a symbolic expression.
    """

    ode_order: int | None = Field(default=None, ge=1)
    pde_order: int | None = Field(default=None, ge=1)
    linearity: Linearity | None = None
    autonomous: bool | None = None
    explicit_time_dependence: bool | None = None
    derivative_operators: tuple[str, ...] | None = None
    integral_operators: tuple[str, ...] | None = None
    state_variables: tuple[str, ...] | None = None
    parameters: tuple[str, ...] | None = None
    forcing_terms: tuple[str, ...] | None = None
    conservation_laws: tuple[str, ...] | None = None
    symmetries: tuple[str, ...] | None = None
    analytical_solution_families: tuple[str, ...] | None = None
    applicable_numerical_methods: tuple[str, ...] | None = None

    @model_validator(mode="after")
    def unique_features(self) -> Self:
        for field_name in (
            "derivative_operators",
            "integral_operators",
            "state_variables",
            "parameters",
            "forcing_terms",
            "conservation_laws",
            "symmetries",
            "analytical_solution_families",
            "applicable_numerical_methods",
        ):
            values = getattr(self, field_name)
            if values is not None and len(values) != len(set(values)):
                raise ValueError(f"{field_name} must not contain duplicates")
        if self.ode_order is not None and self.pde_order is not None:
            raise ValueError("a structure cannot declare both ODE and PDE order")
        return self


class ApplicationMappingLevel(StrEnum):
    """Strength of a proposed or evidenced cross-domain mapping."""

    HYPOTHESIS = "hypothesis"
    STRUCTURAL_SIMILARITY = "structural_similarity"
    MATHEMATICAL_MAPPING = "mathematical_mapping"
    PHYSICALLY_JUSTIFIED_MODEL = "physically_justified_model"
    EMPIRICALLY_INVESTIGATED = "empirically_investigated"
    INDEPENDENTLY_VALIDATED = "independently_validated"


class QuantityDefinition(RegistryModel):
    """A symbol or quantity with optional unit and dimensional metadata."""

    symbol: str = Field(min_length=1)
    name: str = Field(min_length=1)
    meaning: str = Field(min_length=1)
    si_unit: str | None = None
    dimensions: str | None = None


class SourceReference(RegistryModel):
    """Traceable bibliographic, archival, or dataset provenance."""

    source_id: str = Field(pattern=_IDENTIFIER_PATTERN)
    title: str = Field(min_length=1)
    authors: tuple[str, ...] = ()
    responsible_organization: str | None = None
    publication_date: date | None = None
    year: int | None = Field(default=None, ge=0)
    edition: str | None = None
    source_type: SourceType = SourceType.OTHER
    publisher_or_journal: str | None = None
    doi: str | None = None
    url: str | None = None
    isbn: str | None = None
    other_identifier: str | None = None
    locator: str | None = None
    accessed_date: date | None = None
    reliability_notes: str | None = None
    limitations: str | None = None


class EvidenceClaim(RegistryModel):
    """A claim paired with its explicit evidence category and provenance."""

    category: EvidenceCategory
    statement: str = Field(min_length=1)
    source_reference_ids: tuple[str, ...] = ()

    @model_validator(mode="after")
    def strong_evidence_has_provenance(self) -> Self:
        needs_source = {
            EvidenceCategory.ESTABLISHED_PHYSICAL_RESULT,
            EvidenceCategory.EMPIRICALLY_INVESTIGATED,
            EvidenceCategory.INDEPENDENTLY_VALIDATED,
        }
        if self.category in needs_source and not self.source_reference_ids:
            raise ValueError(
                f"claims categorized as {self.category.value} require "
                "a source reference"
            )
        if len(self.source_reference_ids) != len(set(self.source_reference_ids)):
            raise ValueError("source reference IDs in a claim must be unique")
        return self


class EquationRelation(RegistryModel):
    """An evidence-tagged link to another equation record."""

    target_equation_id: str = Field(pattern=_IDENTIFIER_PATTERN)
    relationship: RelationshipType
    evidence: EvidenceClaim


class ApplicationMapping(RegistryModel):
    """An explicitly graded relation between an equation and another domain."""

    domain: str = Field(min_length=1)
    application: str = Field(min_length=1)
    level: ApplicationMappingLevel
    evidence: EvidenceClaim

    @model_validator(mode="after")
    def evidence_matches_mapping_level(self) -> Self:
        category_by_level = {
            ApplicationMappingLevel.HYPOTHESIS: (
                EvidenceCategory.APPLICATION_HYPOTHESIS
            ),
            ApplicationMappingLevel.STRUCTURAL_SIMILARITY: (
                EvidenceCategory.KNOWN_ANALOGY
            ),
            ApplicationMappingLevel.MATHEMATICAL_MAPPING: (
                EvidenceCategory.MATHEMATICALLY_DERIVED
            ),
            ApplicationMappingLevel.PHYSICALLY_JUSTIFIED_MODEL: (
                EvidenceCategory.ESTABLISHED_PHYSICAL_RESULT
            ),
            ApplicationMappingLevel.EMPIRICALLY_INVESTIGATED: (
                EvidenceCategory.EMPIRICALLY_INVESTIGATED
            ),
            ApplicationMappingLevel.INDEPENDENTLY_VALIDATED: (
                EvidenceCategory.INDEPENDENTLY_VALIDATED
            ),
        }
        required_category = category_by_level[self.level]
        if self.evidence.category != required_category:
            raise ValueError(
                f"mapping level {self.level.value} requires evidence category "
                f"{required_category.value}"
            )
        return self


class EquationRecord(RegistryModel):
    """Metadata for one equation, model, or explicitly grouped equation system."""

    equation_id: str = Field(pattern=_IDENTIFIER_PATTERN)
    canonical_name: str = Field(min_length=1)
    alternative_names: tuple[str, ...] = ()
    branch_id: str = Field(pattern=_IDENTIFIER_PATTERN)
    subfield_id: str | None = Field(default=None, pattern=_IDENTIFIER_PATTERN)
    mathematical_expression: str = Field(min_length=1)
    mathematical_structure: MathematicalStructure | None = None
    description: str | None = None
    symbol_definitions: tuple[QuantityDefinition, ...] = ()
    parameters: tuple[QuantityDefinition, ...] = ()
    state_variables: tuple[QuantityDefinition, ...] = ()
    assumptions: tuple[str, ...] = ()
    validity_domain: str | None = None
    applicability_limits: tuple[str, ...] = ()
    initial_conditions: tuple[str, ...] = ()
    boundary_conditions: tuple[str, ...] = ()
    mathematical_classifications: tuple[MathematicalClassification, ...] = ()
    related_equations: tuple[EquationRelation, ...] = ()
    physical_interpretations: tuple[EvidenceClaim, ...] = ()
    established_applications: tuple[EvidenceClaim, ...] = ()
    application_mappings: tuple[ApplicationMapping, ...] = ()
    source_references: tuple[SourceReference, ...] = ()
    verification_status: VerificationStatus = VerificationStatus.UNASSESSED
    implementation_status: ImplementationStatus = ImplementationStatus.NOT_IMPLEMENTED
    numerical_methods: tuple[str, ...] = ()
    tags: tuple[str, ...] = ()

    @model_validator(mode="after")
    def local_identifiers_are_unique(self) -> Self:
        if len(self.mathematical_classifications) != len(
            set(self.mathematical_classifications)
        ):
            raise ValueError("mathematical classifications must be unique")
        source_ids = [source.source_id for source in self.source_references]
        if len(source_ids) != len(set(source_ids)):
            raise ValueError("source reference IDs must be unique within an equation")
        relation_ids = [
            relation.target_equation_id for relation in self.related_equations
        ]
        if len(relation_ids) != len(set(relation_ids)):
            raise ValueError("related equation IDs must be unique within an equation")
        if self.equation_id in relation_ids:
            raise ValueError("an equation cannot declare a relation to itself")
        if any(
            claim.category != EvidenceCategory.ESTABLISHED_PHYSICAL_RESULT
            for claim in self.established_applications
        ):
            raise ValueError(
                "established_applications must use established physical-result evidence"
            )
        return self


class EquationRegistry(RegistryModel):
    """Immutable registry validating taxonomy, equation, and source references."""

    taxonomy: PhysicsTaxonomy = Field(default_factory=default_physics_taxonomy)
    equations: tuple[EquationRecord, ...] = ()

    @model_validator(mode="after")
    def references_are_consistent(self) -> Self:
        by_id = {record.equation_id: record for record in self.equations}
        if len(by_id) != len(self.equations):
            raise ValueError("equation IDs must be unique within a registry")

        for record in self.equations:
            branch = self.taxonomy.branch(record.branch_id)
            if branch is None:
                raise ValueError(
                    f"equation {record.equation_id!r} references unknown branch "
                    f"{record.branch_id!r}"
                )
            if record.subfield_id is not None and not self.taxonomy.contains_subfield(
                record.branch_id, record.subfield_id
            ):
                raise ValueError(
                    f"equation {record.equation_id!r} references subfield "
                    f"{record.subfield_id!r} outside branch {record.branch_id!r}"
                )

            source_ids = {source.source_id for source in record.source_references}
            claims = (
                *record.physical_interpretations,
                *record.established_applications,
                *(relation.evidence for relation in record.related_equations),
                *(mapping.evidence for mapping in record.application_mappings),
            )
            for claim in claims:
                missing_sources = set(claim.source_reference_ids) - source_ids
                if missing_sources:
                    raise ValueError(
                        f"equation {record.equation_id!r} has claims referencing "
                        f"unknown source IDs: {sorted(missing_sources)}"
                    )
            for relation in record.related_equations:
                if relation.target_equation_id not in by_id:
                    raise ValueError(
                        f"equation {record.equation_id!r} references unknown equation "
                        f"{relation.target_equation_id!r}"
                    )
        sources_by_id: dict[str, SourceReference] = {}
        for record in self.equations:
            for source in record.source_references:
                existing_source = sources_by_id.get(source.source_id)
                if existing_source is not None and existing_source != source:
                    raise ValueError(
                        f"source ID {source.source_id!r} has conflicting metadata"
                    )
                sources_by_id[source.source_id] = source
        return self

    def get(self, equation_id: str) -> EquationRecord:
        """Return an equation by stable ID, or raise ``RegistryError``."""
        for equation in self.equations:
            if equation.equation_id == equation_id:
                return equation
        raise RegistryError(f"unknown equation ID: {equation_id}")

    def search(self, query: str) -> tuple[EquationRecord, ...]:
        """Find records whose names, expression, description, or tags match text."""
        normalized_query = query.strip().casefold()
        if not normalized_query:
            return ()
        return tuple(
            equation
            for equation in self.equations
            if normalized_query
            in " ".join(
                (
                    equation.canonical_name,
                    *equation.alternative_names,
                    equation.mathematical_expression,
                    equation.description or "",
                    *equation.tags,
                )
            ).casefold()
        )

    def register(self, equation: EquationRecord) -> EquationRegistry:
        """Return a validated registry containing one additional equation."""
        if any(item.equation_id == equation.equation_id for item in self.equations):
            raise RegistryError(f"duplicate equation ID: {equation.equation_id}")
        return EquationRegistry(
            taxonomy=self.taxonomy,
            equations=(*self.equations, equation),
        )
