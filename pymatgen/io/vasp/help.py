"""
Get help with VASP parameters from VASP wiki.
"""

import requests
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

    def print_help(self, tag):
        """
        Print the help for a TAG.

        Args:
            tag (str): Tag used in VASP.
        """
        print(self.get_help(tag))

    def get_help(self, tag):
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
        main_doc = soup.find(id="bodyContent")
        contents = main_doc.text
        return contents

    def get_incar_tags(self):
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
                    tags.append(child.text)
        return tags


if __name__ == "__main__":
    doc = VaspDoc()
    doc.print_help("ISYM")
    print(doc.get_incar_tags())
