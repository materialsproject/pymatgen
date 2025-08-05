from __future__ import annotations

import io
from unittest.mock import patch

import pytest
import requests

from pymatgen.io.vasp.help import VaspDoc

BeautifulSoup = pytest.importorskip("bs4").BeautifulSoup


try:
    response = requests.get("https://www.vasp.at/wiki/index.php/The_VASP_Manual", timeout=5)
    website_down = response.status_code != 200
    reason = f"VASP manual returned status {response.status_code}"
except (requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout) as e:
    website_down = True
    reason = f"VASP manual not reachable: {type(e).__name__}"

if website_down:
    pytest.skip(reason, allow_module_level=True)


class TestVaspDoc:
    @pytest.mark.parametrize("tag", ["ISYM"])
    def test_print_help(self, tag):
        vasp_doc = VaspDoc()
        with patch("sys.stdout", new=io.StringIO()) as fake_stdout:
            vasp_doc.print_help(tag)
            output = fake_stdout.getvalue()

        assert tag in output

    @pytest.mark.parametrize("tag", ["ISYM"])
    def test_get_help(self, tag):
        docstr = VaspDoc.get_help(tag)

        assert tag in docstr

    def test_get_incar_tags(self):
        """Get all INCAR tags and check incar_parameters.json file."""
        incar_tags_wiki = VaspDoc.get_incar_tags()
        assert isinstance(incar_tags_wiki, list)

        known_incar_tags = ("ENCUT", "ISMEAR")
        for tag in known_incar_tags:
            assert tag in incar_tags_wiki
