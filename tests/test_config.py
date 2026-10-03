"""Tests for application configuration."""

import pytest
from pydantic import ValidationError

from newton_lab.config import AppConfig


def test_default_config() -> None:
    config = AppConfig()

    assert config.app_name == "Newton Lab"
    assert config.environment == "development"
    assert config.log_level == "INFO"
    assert config.random_seed == 0
    assert config.deterministic is True


def test_valid_custom_config() -> None:
    config = AppConfig(
        app_name="Research run",
        environment="test",
        log_level="DEBUG",
        random_seed=42,
        deterministic=False,
    )

    assert config.app_name == "Research run"
    assert config.environment == "test"
    assert config.log_level == "DEBUG"
    assert config.random_seed == 42
    assert config.deterministic is False


@pytest.mark.parametrize(
    "values",
    [
        {"environment": "staging"},
        {"log_level": "TRACE"},
        {"random_seed": -1},
        {"app_name": ""},
        {"unknown_setting": True},
    ],
)
def test_invalid_config(values: dict[str, object]) -> None:
    with pytest.raises(ValidationError):
        AppConfig(**values)
