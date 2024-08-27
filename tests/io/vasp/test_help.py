from __future__ import annotations

import io
from unittest.mock import patch

import pytest
import requests

from pymatgen.io.vasp.help import VaspDoc

try:
    from bs4 import BeautifulSoup
except ImportError:
    BeautifulSoup = None


try:
    website_down = requests.get("https://www.vasp.at", timeout=60).status_code != 200
except requests.exceptions.ConnectionError:
    website_down = True


@pytest.mark.skipif(BeautifulSoup is None, reason="BeautifulSoup4 is not installed")
@pytest.mark.skipif(website_down, reason="www.vasp.at is down")
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
        vasp_doc = VaspDoc()
        docstr = vasp_doc.get_help(tag)

        assert tag in docstr

    def test_get_incar_tags(self):
        incar_tags = VaspDoc.get_incar_tags()
        assert incar_tags
