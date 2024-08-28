from __future__ import annotations

import io
import json
from pathlib import Path
from unittest.mock import patch

import pytest
import requests

from pymatgen.io.vasp.help import VaspDoc

BeautifulSoup = pytest.importorskip("bs4").BeautifulSoup


try:
    website_down = requests.get("https://www.vasp.at", timeout=5).status_code != 200
except (requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout):
    website_down = True


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
        docstr = VaspDoc.get_help(tag)

        assert tag in docstr

    def test_get_incar_tags(self):
        """Get all INCAR tags and check incar_parameters.json file."""
        incar_tags_wiki: set[str] = {tag.replace(" ", "_") for tag in VaspDoc.get_incar_tags()}

        incar_tags_json_file = (
            Path(__file__).resolve().parent.parent.parent.parent / "src/pymatgen/io/vasp/incar_parameters.json"
        )
        with open(incar_tags_json_file, encoding="utf-8") as file:
            incar_tags_json: set[str] = set(json.load(file).keys())

        if tags_wiki_only := incar_tags_wiki.difference(incar_tags_json):
            raise RuntimeError(f"{len(tags_wiki_only)} INCAR tags are missing in json file: {tags_wiki_only}")

        # The following is commented out because we want to make
        # the json more permissive and still include deprecated tags
        # if tags_json_only := incar_tags_json.difference(incar_tags_wiki):
        #     raise RuntimeError(
        #     f"{len(tags_json_only)} INCAR tags might have been removed from VASP wiki: {tags_json_only}"
        # )
