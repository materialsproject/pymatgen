"""Get help with VASP parameters from VASP wiki."""

from __future__ import annotations

import re

import orjson
import requests
from monty.dev import requires

try:
    from bs4 import BeautifulSoup
except ImportError:
    BeautifulSoup = None  # type:ignore[assignment]


class VaspDoc:
    """A VASP documentation helper."""

    @requires(BeautifulSoup is not None, "BeautifulSoup4 must be installed to fetch from the VASP wiki.")
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
        # Use Mediawiki API as documented in
        # https://www.vasp.at/wiki/api.php?action=help&modules=query
        url = (
            "https://www.vasp.at/wiki/api.php?"
            "action=query&list=categorymembers"
            "&cmtitle=Category:INCAR_tag"
            "&cmlimit=500&format=json"
        )
        response = requests.get(url, timeout=60)
        response_dict = orjson.loads(response.text)

        def extract_titles(data):
            """Extract keywords from from Wikimedia response data.
            See https://www.vasp.at/wiki/api.php?action=help&modules=query%2Bcategorymembers
            Returns: List of keywords as strings.
            """
            return [category_data["title"] for category_data in data["query"]["categorymembers"]]

        tags = extract_titles(response_dict)

        # If there are more than 500 items in the response, we will
        # get 'continue' field in the response
        # See https://www.mediawiki.org/wiki/API:Continue
        while "continue" in response_dict:
            response = requests.get(url + f"&cmcontinue={response_dict['continue']['cmcontinue']}", timeout=60)
            response_dict = orjson.loads(response.text)
            tags = tags + extract_titles(response_dict)

        return tags
