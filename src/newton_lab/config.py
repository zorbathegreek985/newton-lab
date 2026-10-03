"""Validated application configuration."""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class AppConfig(BaseModel):
    """Settings shared by Newton Lab command-line and research workflows."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    app_name: str = Field(default="Newton Lab", min_length=1)
    environment: Literal["development", "test", "production"] = "development"
    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO"
    random_seed: int | None = Field(default=0, ge=0)
    deterministic: bool = True
