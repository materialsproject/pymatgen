"""Store generic minimization settings read from a JDFTx out file.

This module contains the JMinSettings class for storing generic minimization
and mutants for storing specific minimization settings read from a JDFTx out
file.
"""

from __future__ import annotations

import pprint
from dataclasses import dataclass
from typing import Any


@dataclass
class JMinSettings:
    """Store generic minimization settings read from a JDFTx out file.

    Store generic minimization settings read from a JDFTx out file.
    """

    params: dict[str, Any] | None = None

    def __init__(
        self,
        params: dict[str, Any] | None = None,
    ) -> None:
        """Initialize a generic JMinSettings class.

        Args:
            params (dict[str, Any] | None): A dictionary of minimization settings.
        """
        self.params = None if params is None else dict(params)

    def __str__(self) -> str:
        """Return a string representation of the minimization settings."""
        return pprint.pformat(self)


@dataclass
class JMinSettingsElectronic(JMinSettings):
    """JMInSettings mutant for electronic minimization settings.

    A class for storing electronic minimization settings read from a JDFTx out file.
    """

    start_flag: str = "electronic-minimize"

    def __init__(
        self,
        params: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(params=params)


@dataclass
class JMinSettingsFluid(JMinSettings):
    """JMInSettings mutant for fluid minimization settings.

    A class for storing fluid minimization settings read from a JDFTx out file.
    """

    start_flag: str = "fluid-minimize"

    def __init__(
        self,
        params: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(params=params)


@dataclass
class JMinSettingsLattice(JMinSettings):
    """JMInSettings mutant for lattice minimization settings.

    A class for storing lattice minimization settings read from a JDFTx out file.
    """

    start_flag: str = "lattice-minimize"

    def __init__(
        self,
        params: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(params=params)


@dataclass
class JMinSettingsIonic(JMinSettings):
    """JMInSettings mutant for ionic minimization settings.

    A class for storing ionic minimization settings read from a JDFTx out file.
    """

    start_flag: str = "ionic-minimize"

    def __init__(
        self,
        params: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(params=params)
