from __future__ import annotations

import pytest

from pymatgen.analysis.chemenv.utils.chemenv_errors import (
    AbstractChemenvError,
    ChemenvError,
    EquivalentSiteSearchError,
    NeighborsNotComputedChemenvError,
    SolidAngleError,
)


@pytest.mark.parametrize(
    ("error", "expected"),
    [
        (AbstractChemenvError("MyClass", "my_method", "oops"), "MyClass: my_method\n'oops'"),
        (ChemenvError("MyClass", "my_method", "oops"), "MyClass: my_method\noops"),
        (
            NeighborsNotComputedChemenvError("site0"),
            "The neighbors were not computed for the following site : \nsite0",
        ),
        (EquivalentSiteSearchError("site1"), "Equivalent site could not be found for the following site : site1"),
        (
            SolidAngleError(1.5),
            "Value of cosinus (1.5) from which an angle should be retrieved is not between -1.0 and 1.0",
        ),
    ],
)
def test_error_str(error, expected):
    assert str(error) == expected
    with pytest.raises(type(error)):
        raise error


def test_error_attributes():
    err = ChemenvError("cls", "meth", "msg")
    assert (err.cls, err.method, err.msg) == ("cls", "meth", "msg")
    assert NeighborsNotComputedChemenvError("s").site == "s"
    assert EquivalentSiteSearchError("s").site == "s"
    assert SolidAngleError(-2).cosinus == -2
    assert issubclass(SolidAngleError, AbstractChemenvError)
