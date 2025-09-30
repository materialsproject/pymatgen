from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable


def version_processor(
    min_version: str = "0.0", max_version: str | None = None
) -> Callable[[Callable], Callable]:
    """Decorator to mark a method as a version processor.

    Args:
        min_version (str): Minimum version for which the processor is valid.
        max_version (str | None): Maximum version for which the processor is valid.

    Returns:
        Callable[[Callable], Callable]: Decorator for versioned processor methods.
    """

    def decorator(func: Callable) -> Callable:
        setattr(func, "version_info", (min_version, max_version))  # NOQA: B010

        return func

    return decorator
