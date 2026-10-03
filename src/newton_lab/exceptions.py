"""Project-specific exception types."""


class NewtonLabError(Exception):
    """Base exception for errors raised by Newton Lab."""


class ScientificValidationError(NewtonLabError, ValueError):
    """Raised when data fails a scientific input validation check."""


class IntegrationError(NewtonLabError, RuntimeError):
    """Raised when a numerical integrator cannot complete a simulation."""


class RegistryError(NewtonLabError, ValueError):
    """Raised when a registry operation cannot be completed consistently."""


class CandidateIngestionError(RegistryError):
    """Raised for unreadable candidate inputs, retaining submitted text if known."""

    def __init__(self, message: str, *, raw_submission: str | None = None) -> None:
        super().__init__(message)
        self.raw_submission = raw_submission


class PersistenceError(RegistryError):
    """Raised when a knowledge store cannot be safely read or written."""
