"""Synthetic tests for schema-aware market CSV exploration."""

from pathlib import Path

import pytest

from newton_lab.exceptions import ScientificValidationError
from newton_lab.market_data_exploration import (
    analyze_market_csv,
    create_exploration_plots,
)


@pytest.fixture
def synthetic_market_csv(tmp_path: Path) -> Path:
    path = tmp_path / "synthetic_market.csv"
    path.write_text(
        "TradDt,FinInstrmId,ISIN,TckrSymb,Category,Measure,Other,Blank\n"
        "2026-01-01,001,AA,A,red,0,1,\n"
        "2026-01-01,002,BB,B,blue,-2,-3,NA\n"
        "2026-01-01,002,BB,B,blue,4,9,\n"
        "2026-01-01,002,BB,B,blue,4,9,\n",
        encoding="utf-8",
    )
    return path


def test_profile_reports_raw_shape_categories_numeric_and_identifiers(
    synthetic_market_csv: Path,
) -> None:
    profile = analyze_market_csv(synthetic_market_csv)

    assert (profile.row_count, profile.column_count) == (4, 8)
    assert profile.duplicate_row_count == 1
    assert (
        next(item for item in profile.columns if item.name == "Blank").missing_count
        == 4
    )
    assert (
        next(item for item in profile.columns if item.name == "TradDt").distinct_count
        == 1
    )

    category = next(
        item for item in profile.categorical_summaries if item.name == "Category"
    )
    assert [(item.value, item.count) for item in category.value_counts] == [
        ("blue", 3),
        ("red", 1),
    ]
    numeric = next(item for item in profile.numeric_summaries if item.name == "Measure")
    assert numeric.count == 4
    assert numeric.zero_count == 1
    assert numeric.negative_count == 1
    assert numeric.positive_count == 2
    assert numeric.mean == 1.5

    candidate = next(
        item for item in profile.identifier_candidates if item.name == "FinInstrmId"
    )
    assert candidate.non_missing_count == 4
    assert candidate.distinct_count == 2
    assert candidate.repeated_value_count == 2
    assert not candidate.unique_within_file

    pair = next(
        item
        for item in profile.numeric_column_correlations
        if (item.left_column, item.right_column) == ("Measure", "Other")
    )
    assert pair.pair_count == 4
    assert pair.pearson_r == pytest.approx(1.0)


def test_pairwise_association_uses_row_aligned_complete_pairs(tmp_path: Path) -> None:
    path = tmp_path / "missing_pairs.csv"
    path.write_text("A,B\n1,2\nNA,4\n3,6\n", encoding="utf-8")

    pair = analyze_market_csv(path).numeric_column_correlations[0]

    assert pair.pair_count == 2
    assert pair.pearson_r == pytest.approx(1.0)


def test_rejects_rows_with_inconsistent_width(tmp_path: Path) -> None:
    path = tmp_path / "bad_width.csv"
    path.write_text("A,B\n1,2,3\n", encoding="utf-8")

    with pytest.raises(ScientificValidationError, match="header width"):
        analyze_market_csv(path)


def test_plots_are_created_and_existing_outputs_are_not_overwritten(
    synthetic_market_csv: Path, tmp_path: Path
) -> None:
    pytest.importorskip("matplotlib")
    output = tmp_path / "figures"
    profile = analyze_market_csv(synthetic_market_csv)

    plots = create_exploration_plots(profile, output, source_path=synthetic_market_csv)

    assert plots
    assert all(path.is_file() and path.stat().st_size > 0 for path in plots)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        create_exploration_plots(profile, output, source_path=synthetic_market_csv)
