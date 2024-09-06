"""Get help with VASP parameters from VASP wiki."""

from __future__ import annotations

import re

import requests
from monty.dev import requires

try:
    from bs4 import BeautifulSoup
except ImportError:
    BeautifulSoup = None


class VaspDoc:
    """A VASP documentation helper."""

    @requires(BeautifulSoup, "BeautifulSoup4 must be installed to fetch from the VASP wiki.")
    def __init__(self) -> None:
        """Init for VaspDoc."""
        self.url_template = "https://www.vasp.at/wiki/index.php/%s"

    def print_help(self, tag: str) -> None:
        """
        Print the help for a TAG.

        Args:
            tag (str): Tag used in VASP.
        """
        print(self.get_help(tag))

    def print_jupyter_help(self, tag: str) -> None:
        """
        Display HTML help in ipython notebook.

        Args:
            tag (str): Tag used in VASP.
        """
        html_str = self.get_help(tag, "html")
        from IPython.core.display import HTML, display

        display(HTML(html_str))

    @classmethod
    def get_help(cls, tag: str, fmt: str = "text") -> str:
        """Get help on a VASP tag.

        Args:
            tag (str): VASP tag, e.g. ISYM.

        Returns:
            Help text.
        """
        tag = tag.upper()
        response = requests.get(
            f"https://www.vasp.at/wiki/index.php/{tag}",
            timeout=60,
        )
        soup = BeautifulSoup(response.text, features="html.parser")
        main_doc = soup.find(id="mw-content-text")
        if fmt == "text":
            output = main_doc.text
            return re.sub("\n{2,}", "\n\n", output)

        return str(main_doc)

    @classmethod
    def get_incar_tags(cls) -> list[str]:
        """Get a list of all INCAR tags from the VASP wiki."""
        tags = []
        for url in (
            "https://www.vasp.at/wiki/index.php/Category:INCAR_tag",
            "https://www.vasp.at/wiki/index.php?title=Category:INCAR_tag&pagefrom=LREAL#mw-pages",
            "https://www.vasp.at/wiki/index.php?title=Category:INCAR_tag&pagefrom=Profiling#mw-pages",
        ):
            response = requests.get(url, timeout=60)
            soup = BeautifulSoup(response.text, features="html.parser")
            for div in soup.findAll("div", {"class": "mw-category-group"}):
                children = div.findChildren("li")
                for child in children:
                    tags.append(child.text.strip())
        return tags
