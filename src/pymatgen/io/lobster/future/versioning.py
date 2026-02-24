from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

F = TypeVar("F", bound=Callable)


def version_processor(min_version: str = "0.0", max_version: str | None = None) -> Callable[[F], F]:
    """Decorator to mark a method as a version processor.

    Args:
        min_version (str): Minimum version for which the processor is valid.
        max_version (str | None): Maximum version for which the processor is valid.

    Returns:
        Callable[[F], F]: Decorator for versioned processor methods.
    """

    def decorator(func: F) -> F:
        setattr(func, "version_info", (min_version, max_version))  # NOQA: B010

        return func

    return decorator
