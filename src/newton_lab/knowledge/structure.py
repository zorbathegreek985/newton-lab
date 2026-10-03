"""Interpretable structure comparison and conservative relationship discovery.

This module compares declared metadata only. It does not parse expressions,
derive identities, or infer physical equivalence.
"""

from __future__ import annotations

from enum import StrEnum
from itertools import combinations
from typing import Literal

from newton_lab.knowledge.registry import EquationRecord, RelationshipType
from newton_lab.knowledge.workflow import (
    EquationRelationship,
    KnowledgeBase,
    RelationshipReviewStatus,
    WorkflowModel,
)


class FeatureComparisonStatus(StrEnum):
    """Whether an explicitly recorded feature is shared, different, or unknown."""

    SHARED = "shared"
    DIFFERENT = "different"
    PARTIALLY_SHARED = "partially_shared"
    UNKNOWN = "unknown_or_unassessed"


class FeatureComparison(WorkflowModel):
    """One named feature comparison; values are rendered as stable primitives."""

    feature_name: str
    status: FeatureComparisonStatus
    left_value: str | bool | int | tuple[str, ...] | None
    right_value: str | bool | int | tuple[str, ...] | None
    shared_values: tuple[str, ...] = ()
    left_only_values: tuple[str, ...] = ()
    right_only_values: tuple[str, ...] = ()


class StructureComparison(WorkflowModel):
    """Transparent comparison with no scalar similarity or equivalence verdict."""

    left_equation_id: str
    right_equation_id: str
    feature_comparisons: tuple[FeatureComparison, ...]
    shared_assumptions: tuple[str, ...]
    differing_assumptions: tuple[str, ...]
    shared_conditions: tuple[str, ...]
    differing_conditions: tuple[str, ...]
    documented_relationships: tuple[EquationRelationship, ...]
    comparison_kind: Literal["structural_metadata_comparison"] = (
        "structural_metadata_comparison"
    )
    physical_equivalence_concluded: Literal[False] = False

    @property
    def shared_features(self) -> tuple[str, ...]:
        """Names of features explicitly shared by both records."""
        return tuple(
            item.feature_name
            for item in self.feature_comparisons
            if item.status
            in {
                FeatureComparisonStatus.SHARED,
                FeatureComparisonStatus.PARTIALLY_SHARED,
            }
        )

    @property
    def different_features(self) -> tuple[str, ...]:
        """Names of features explicitly different between records."""
        return tuple(
            item.feature_name
            for item in self.feature_comparisons
            if item.status
            in {
                FeatureComparisonStatus.DIFFERENT,
                FeatureComparisonStatus.PARTIALLY_SHARED,
            }
        )

    @property
    def unknown_features(self) -> tuple[str, ...]:
        """Names of features unknown or unassessed on either record."""
        return tuple(
            item.feature_name
            for item in self.feature_comparisons
            if item.status == FeatureComparisonStatus.UNKNOWN
        )


_STRUCTURE_FIELDS = (
    "classifications",
    "ode_order",
    "pde_order",
    "linearity",
    "autonomous",
    "explicit_time_dependence",
    "derivative_operators",
    "integral_operators",
    "state_variables",
    "parameters",
    "forcing_terms",
    "conservation_laws",
    "symmetries",
    "analytical_solution_families",
    "applicable_numerical_methods",
)


def _feature_value(record: EquationRecord, field_name: str) -> object:
    if field_name == "classifications":
        if not record.mathematical_classifications:
            return None
        return tuple(
            classification.value
            for classification in record.mathematical_classifications
        )
    structure = record.mathematical_structure
    if structure is None:
        return None
    value = getattr(structure, field_name)
    if value is None:
        return None
    if isinstance(value, tuple):
        return tuple(item.value if hasattr(item, "value") else item for item in value)
    return value.value if hasattr(value, "value") else value


def _comparable_value(value: object) -> str | bool | int | tuple[str, ...] | None:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, tuple) and all(isinstance(item, str) for item in value):
        return value
    return str(value)


def _ordered_set_comparison(
    left_values: tuple[str, ...], right_values: tuple[str, ...]
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    left = set(left_values)
    right = set(right_values)
    return tuple(sorted(left & right)), tuple(sorted(left ^ right))


def compare_equations(
    left: EquationRecord,
    right: EquationRecord,
    *,
    knowledge_base: KnowledgeBase | None = None,
) -> StructureComparison:
    """Compare declared structural features, assumptions, conditions, and links.

    A missing structure or ``None`` feature is unassessed. An empty feature
    tuple is an explicit assessed absence. No expression-text matching occurs.
    """
    if left.equation_id == right.equation_id:
        raise ValueError("compare_equations requires two distinct equation IDs")

    comparisons: list[FeatureComparison] = []
    for field_name in _STRUCTURE_FIELDS:
        left_raw = _feature_value(left, field_name)
        right_raw = _feature_value(right, field_name)
        if left_raw is None or right_raw is None:
            status = FeatureComparisonStatus.UNKNOWN
            shared_values: tuple[str, ...] = ()
            left_only_values: tuple[str, ...] = ()
            right_only_values: tuple[str, ...] = ()
        elif isinstance(left_raw, tuple) and isinstance(right_raw, tuple):
            left_set = {str(value) for value in left_raw}
            right_set = {str(value) for value in right_raw}
            shared_values = tuple(sorted(left_set & right_set))
            left_only_values = tuple(sorted(left_set - right_set))
            right_only_values = tuple(sorted(right_set - left_set))
            if shared_values and not left_only_values and not right_only_values:
                status = FeatureComparisonStatus.SHARED
            elif shared_values:
                status = FeatureComparisonStatus.PARTIALLY_SHARED
            else:
                status = FeatureComparisonStatus.DIFFERENT
        elif left_raw == right_raw:
            status = FeatureComparisonStatus.SHARED
            shared_values = ()
            left_only_values = ()
            right_only_values = ()
        else:
            status = FeatureComparisonStatus.DIFFERENT
            shared_values = ()
            left_only_values = ()
            right_only_values = ()
        comparisons.append(
            FeatureComparison(
                feature_name=field_name,
                status=status,
                left_value=_comparable_value(left_raw),
                right_value=_comparable_value(right_raw),
                shared_values=shared_values,
                left_only_values=left_only_values,
                right_only_values=right_only_values,
            )
        )

    shared_assumptions, differing_assumptions = _ordered_set_comparison(
        left.assumptions, right.assumptions
    )
    left_conditions = (*left.initial_conditions, *left.boundary_conditions)
    right_conditions = (*right.initial_conditions, *right.boundary_conditions)
    shared_conditions, differing_conditions = _ordered_set_comparison(
        left_conditions, right_conditions
    )
    documented: tuple[EquationRelationship, ...] = ()
    if knowledge_base is not None:
        pair = {left.equation_id, right.equation_id}
        documented = tuple(
            relationship
            for relationship in knowledge_base.relationships
            if {relationship.source_equation_id, relationship.target_equation_id}
            == pair
        )
    return StructureComparison(
        left_equation_id=left.equation_id,
        right_equation_id=right.equation_id,
        feature_comparisons=tuple(comparisons),
        shared_assumptions=shared_assumptions,
        differing_assumptions=differing_assumptions,
        shared_conditions=shared_conditions,
        differing_conditions=differing_conditions,
        documented_relationships=documented,
    )


def _slug(value: str) -> str:
    return "_".join(part for part in value.lower().replace("-", " ").split() if part)


def propose_relationships(
    knowledge_base: KnowledgeBase,
) -> tuple[EquationRelationship, ...]:
    """Suggest review-required edges from explicit shared classifications/operators.

    Suggestions carry their matched metadata in the justification. This
    deliberately conservative first pass never creates sources or acceptance.
    """
    existing = {
        (
            relation.source_equation_id,
            relation.target_equation_id,
            relation.relationship_type,
        )
        if relation.directional
        else (
            *sorted((relation.source_equation_id, relation.target_equation_id)),
            relation.relationship_type,
        )
        for relation in knowledge_base.relationships
    }
    proposed: list[EquationRelationship] = []
    ordered_records = sorted(
        (entry.record for entry in knowledge_base.entries),
        key=lambda record: record.equation_id,
    )
    for left, right in combinations(ordered_records, 2):
        comparison = compare_equations(left, right)
        shared_classes = next(
            item.shared_values
            for item in comparison.feature_comparisons
            if item.feature_name == "classifications"
        )
        shared_operators = next(
            item.shared_values
            for item in comparison.feature_comparisons
            if item.feature_name == "derivative_operators"
        )
        evidence: list[tuple[RelationshipType, tuple[str, ...]]] = []
        if isinstance(shared_classes, tuple) and shared_classes:
            evidence.append(
                (RelationshipType.SHARES_MATHEMATICAL_STRUCTURE, shared_classes)
            )
        if isinstance(shared_operators, tuple) and shared_operators:
            evidence.append(
                (RelationshipType.SHARES_MATHEMATICAL_OPERATOR, shared_operators)
            )
        for relation_type, matched in evidence:
            identity = (left.equation_id, right.equation_id, relation_type)
            if identity in existing:
                continue
            token = _slug("_".join(matched))
            relationship_id = (
                f"candidate_{left.equation_id}_{right.equation_id}_"
                f"{relation_type.value}_{token}"
            )
            proposed.append(
                EquationRelationship(
                    relationship_id=relationship_id,
                    source_equation_id=left.equation_id,
                    target_equation_id=right.equation_id,
                    relationship_type=relation_type,
                    description="Candidate shared mathematical structure for review.",
                    justification=(
                        "Explicitly recorded features match: " + ", ".join(matched)
                    ),
                    directional=False,
                    review_status=RelationshipReviewStatus.CANDIDATE,
                )
            )
    return tuple(proposed)
