from __future__ import annotations

import pytest
from pytest import approx

from pymatgen.core import Structure
from pymatgen.io.res import AirssProvider, ResParseError, ResWriter
from pymatgen.util.testing import TEST_FILES_DIR

res_coc = f"{TEST_FILES_DIR}/res/coc-115925-9326-14.res"


@pytest.mark.parametrize("provider", [AirssProvider.from_file(res_coc, "strict")])
class TestAirssProvider:
    def test_cggf(self, provider: AirssProvider):
        entry = provider.get_cut_grid_gmax_fsbc()
        assert entry is not None
        cut, gs, gm, fsbc = entry
        assert cut == approx(326.5366)
        assert gs == approx(1.75)
        assert gm == approx(16.201)
        assert fsbc == "automatic"

    def test_moks(self, provider: AirssProvider):
        entry = provider.get_mpgrid_offset_nkpts_spacing()
        assert entry is not None
        grid, offset, nkpts, spacing = entry
        assert grid == (6, 6, 8)
        assert offset == (0, 0, 0)
        assert nkpts == 144
        assert spacing == 0.05

    def test_pspots(self, provider: AirssProvider):
        ps_pots = provider.get_pspots()
        assert ps_pots["Co"] == "3|2.2|2.0|1.0|8|9|10|40:41:32(qc=5)"
        assert ps_pots["C"] == "2|1.4|8|9|10|20:21(qc=5)"

    def test_titl(self, provider: AirssProvider):
        assert provider.seed == "coc-115925-9326-14"
        assert provider.energy == approx(-3.90427411e003)
        assert provider.spacegroup_label == "R3"
        assert provider.pressure == approx(15.0252)
        assert provider.volume == approx(57.051984)

    def test_lattice(self, provider: AirssProvider):
        assert provider.lattice.lengths == approx((5.07144, 5.07144, 3.89024))
        assert provider.lattice.angles == approx((49.32125, 49.32125, 60))

    def test_misc(self, provider: AirssProvider):
        rs_info = provider.get_run_start_info()
        assert rs_info is not None
        date, path = rs_info
        assert path == "/path/to/airss/run"
        assert date.day == 16

        castep_v = provider.get_castep_version()
        assert castep_v is not None
        assert castep_v == "19.11"

        frd = provider.get_func_rel_disp()
        assert frd is not None
        f, r, d = frd
        assert f == "Perdew Burke Ernzerhof"

        airss_v = provider.get_airss_version()
        assert airss_v is not None
        assert airss_v[0] == "0.9.1"
        assert airss_v[1].year == 2018

    def test_entry(self, provider: AirssProvider):
        entry1 = provider.entry
        assert entry1.energy == provider.energy

        string1 = ResWriter(entry1).string
        entry2 = AirssProvider.from_str(string1).entry
        string2 = ResWriter(entry2).string
        assert entry1 == entry2
        assert string1 == string2

    def test_raise(self, provider: AirssProvider):
        entry = provider.entry
        string = ResWriter(entry).string
        string_strip = "\n".join(line for line in string.splitlines() if "REM" not in line)
        prov = AirssProvider.from_str(string_strip, "strict")
        assert entry.structure == prov.entry.structure
        with pytest.raises(ResParseError, match="No CASTEP version found in REM"):
            prov.get_castep_version()

    def test_as_dict(self, provider: AirssProvider):
        verbose_dict = provider.as_dict(verbose=True)

        assert sorted(verbose_dict) == ["@class", "@module", "@version", "parse_rems", "res"]

        # test round-trip serialization/deserialization gives same dict
        assert AirssProvider.from_dict(verbose_dict).as_dict() == verbose_dict

        # non-verbose case
        dct = provider.as_dict(verbose=False)
        assert sorted(dct) == [
            "appearances",
            "energy",
            "integrated_absolute_spin_density",
            "integrated_spin_density",
            "pressure",
            "rems",
            "seed",
            "spacegroup_label",
            "structure",
            "volume",
        ]
        assert dct["seed"] == "coc-115925-9326-14"
        assert dct["energy"] == approx(-3904.2741)
        assert dct["spacegroup_label"] == "R3"
        assert dct["pressure"] == approx(15.0252)
        assert dct["volume"] == approx(57.051984)


class TestSpin:
    def test_read_spin(self):
        with open(res_coc) as f:
            lines = f.readlines()
        # add spin to a line
        lines[25] = lines[25][:-1] + " -1.4\n"
        contents = "".join(lines)
        provider = AirssProvider.from_str(contents)

        for site in provider.structure:
            if site.properties["magmom"] is not None:
                assert site.properties.get("magmom") == approx(-1.4)
                return
        pytest.fail("valid 'magmom' not found in any site properties")

    def test_gh_2938_example(self):
        res_spin_file = f"{TEST_FILES_DIR}/res/spins-in-last-col.res"
        with open(res_spin_file) as res_file:
            contents = res_file.read()

        provider = AirssProvider.from_str(contents)

        for site in provider.structure:
            if site.properties["magmom"] is not None:
                assert site.properties.get("magmom") in (3.31, 4.12, -0.01, -0.04, -0.17)


class TestStructureModule:
    def test_structure_from_file(self):
        structure: Structure = Structure.from_file(res_coc)
        # just check that we read it
        assert structure.lattice.angles == approx((49.32125, 49.32125, 60))
