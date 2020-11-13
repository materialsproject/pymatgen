"""
Get help with VASP parameters from VASP wiki.
"""

import re

import requests
import urllib3
from bs4 import BeautifulSoup


class VaspDoc:
    """
    A VASP documentation helper.
    """

    def __init__(self):
        """
        Init for VaspDoc.
        """
        self.url_template = "http://www.vasp.at/wiki/index.php/%s"
        urllib3.disable_warnings()

    def print_help(self, tag):
        """
        Print the help for a TAG.

        Args:
            tag (str): Tag used in VASP.
        """
        print(self.get_help(tag))

    def print_jupyter_help(self, tag):
        """
        Display HTML help in ipython notebook.

        Args:
            tag (str): Tag used in VASP.
        """
        help = self.get_help(tag, "html")
        from IPython.core.display import display, HTML
        display(HTML(help))

    @classmethod
    def get_help(cls, tag, fmt="text"):
        """
        Get help on a VASP tag.

        Args:
            tag (str): VASP tag, e.g., ISYM.

        Returns:
            Help text.
        """
        tag = tag.upper()
        r = requests.get("http://www.vasp.at/wiki/index.php/%s" % tag, verify=False)
        soup = BeautifulSoup(r.text)
        main_doc = soup.find(id="mw-content-text")
        if fmt == "text":
            output = main_doc.text
            output = re.sub("\n{2,}", "\n\n", output)
        else:
            output = str(main_doc)

        return output

    @classmethod
    def get_incar_tags(cls):
        """
        Returns: All incar tags
        """
        tags = []
        for page in ["http://www.vasp.at/wiki/index.php/Category:INCAR",
                     "http://www.vasp.at/wiki/index.php?title=Category:INCAR&pagefrom=ML+FF+LCONF+DISCARD#mw-pages"]:
            r = requests.get(page, verify=False)
            soup = BeautifulSoup(r.text)
            for div in soup.findAll('div', {'class': 'mw-category-group'}):
                children = div.findChildren('li')
                for child in children:
                    tags.append(child.text.strip())
        return tags


if __name__ == "__main__":
    doc = VaspDoc()
    doc.print_help("ISYM")
    print(doc.get_incar_tags())
