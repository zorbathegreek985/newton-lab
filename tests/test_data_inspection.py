"""Tests for local CSV inspection using synthetic fixtures only."""

from pathlib import Path

import pytest

from newton_lab.data_inspection import inspect_local_csv
from newton_lab.exceptions import ScientificValidationError


def test_inspect_local_csv_reports_shape_quality_and_frequency(tmp_path: Path) -> None:
    csv_path = tmp_path / "synthetic.csv"
    csv_path.write_text(
        "timestamp,value\n"
        "2024-01-01T00:00:00,1\n"
        "2024-01-01T01:00:00,NA\n"
        "2024-01-01T02:00:00,3\n"
        "2024-01-01T02:00:00,3\n",
        encoding="utf-8",
    )

    report = inspect_local_csv(csv_path, timestamp_column="timestamp")

    assert report.columns == ("timestamp", "value")
    assert report.row_count == 4
    assert report.missing_values_by_column == {"timestamp": 0, "value": 1}
    assert report.duplicate_row_count == 1
    assert report.duplicate_timestamp_count == 1
    assert report.invalid_timestamp_count == 0
    assert report.date_start == "2024-01-01T00:00:00"
    assert report.date_end == "2024-01-01T02:00:00"
    assert report.inferred_frequency == "1h"


def test_inspect_local_csv_marks_invalid_dates_and_irregular_frequency(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "irregular.csv"
    csv_path.write_text(
        "when,value\n2024-01-01,1\nnot-a-date,2\n2024-01-03,3\n",
        encoding="utf-8",
    )

    report = inspect_local_csv(csv_path, timestamp_column="when")

    assert report.invalid_timestamp_count == 1
    assert report.date_start == "2024-01-01T00:00:00"
    assert report.date_end == "2024-01-03T00:00:00"
    assert report.inferred_frequency is None


@pytest.mark.parametrize(
    ("path", "message"),
    [
        ("https://example.invalid/data.csv", "existing local CSV"),
        ("missing.csv", "existing local CSV"),
    ],
)
def test_inspect_local_csv_rejects_nonlocal_or_missing_path(
    path: str, message: str
) -> None:
    with pytest.raises(ScientificValidationError, match=message):
        inspect_local_csv(path, timestamp_column="timestamp")


def test_inspect_local_csv_requires_timestamp_column(tmp_path: Path) -> None:
    csv_path = tmp_path / "no_timestamp.csv"
    csv_path.write_text("value\n1\n", encoding="utf-8")

    with pytest.raises(ScientificValidationError, match="not present"):
        inspect_local_csv(csv_path, timestamp_column="timestamp")
