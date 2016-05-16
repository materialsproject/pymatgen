#!/usr/bin/env python
from __future__ import print_function
from pymatgen.io.abinit.xc import XcFunctional

ixc11 = XcFunctional.from_abinit_ixc("11")
assert ixc11.name == "PBE"
assert ixc11.type == "GGA"

d = {ixc11: ixc11.name}
print(d)
assert "PBE" in d
assert ixc11 in d

ixc_101130 = XcFunctional.from_abinit_ixc("-101130")
assert ixc_101130.name == "PBE"
assert ixc_101130 == ixc11

gga_pbe = XcFunctional.from_name("PBE")
assert gga_pbe.name == "PBE"
#gga_pbe_libxc = XcFunctional.from_type_name("GGA", "GGA_X_PBE+GGA_C_PBE")

assert ixc11 == gga_pbe
#assert ixc11 == gga_pbe_libxc
