"""
Optimade support.
"""

from collections import namedtuple
from typing import Dict

import requests

from pymatgen.core.periodic_table import DummySpecies
from pymatgen.core.structure import Structure
from pymatgen.util.sequence import PBar


class OptimadeRester:
    """
    Class to call OPTIMADE-compliant APIs, see optimade.org

    This class is ready to use but considered in-development and subject to change
    until the OPTIMADE paper is published.
    """

    # regenerate on-demand from official providers.json using OptimadeRester.refresh_aliases()
    # these aliases are provided as a convenient shortcut for users of the OptimadeRester class
    aliases = {
        "aflow": "http://aflow.org/API/optimade/",
        "cod": "https://www.crystallography.net/cod/optimade",
        "mcloud.2dstructures": "https://aiida.materialscloud.org/2dstructures/optimade",
        "mcloud.2dtopo": "https://aiida.materialscloud.org/2dtopo/optimade",
        "mcloud.curated-cofs": "https://aiida.materialscloud.org/curated-cofs/optimade",
        "mcloud.li-ion-conductors": "https://aiida.materialscloud.org/li-ion-conductors/optimade",
        "mcloud.optimade-sample": "https://aiida.materialscloud.org/optimade-sample/optimade",
        "mcloud.pyrene-mofs": "https://aiida.materialscloud.org/pyrene-mofs/optimade",
        "mcloud.scdm": "https://aiida.materialscloud.org/autowannier/optimade",
        "mcloud.sssp": "https://aiida.materialscloud.org/sssplibrary/optimade",
        "mcloud.stoceriaitf": "https://aiida.materialscloud.org/stoceriaitf/optimade",
        "mcloud.tc-applicability": "https://aiida.materialscloud.org/tc-applicability/optimade",
        "mcloud.threedd": "https://aiida.materialscloud.org/3dd/optimade",
        "mp": "https://optimade.materialsproject.org",
        "odbx": "https://optimade.odbx.science",
        "omdb.omdb_production": "http://optimade.openmaterialsdb.se",
        "oqmd": "http://oqmd.org/optimade/",
        "tcod": "https://www.crystallography.net/tcod/optimade",
    }

    def __init__(self, alias_or_structure_resource_url="mp"):
        """
        OPTIMADE is an effort to provide a standardized interface to retrieve information
        from many different materials science databases.

        This is a client to retrieve structures from OPTIMADE v1 compliant endpoints. It
        does not yet support all features of the OPTIMADE v1 specification but is intended
        as a way to quickly search an endpoint in a way familiar to users of pymatgen without
        needing to know the full OPTIMADE specification.

        For advanced usage, please see the OPTIMADE documentation at optimate.org and
        consider calling the APIs directly.

        For convenience, known OPTIMADE endpoints have been given aliases in pymatgen to save
        typing the full URL. The current list of aliases is:

        aflow, cod, mcloud.sssp, mcloud.2dstructures, mcloud.2dtopo, mcloud.tc-applicability,
        mcloud.threedd, mcloud.scdm, mcloud.curated-cofs, mcloud.optimade-sample, mcloud.stoceriaitf,
        mcloud.pyrene-mofs, mcloud.li-ion-conductors, mp, odbx, omdb.omdb_production, oqmd, tcod

        To refresh this list of aliases, generated from the current list of OPTIMADE providers
        at optimade.org, call the refresh_aliases() method.

        Args:
            alias_or_structure_resource_url: the alias or structure resource URL
        """

        # TODO: maybe we should use the nice pydantic models from optimade-python-tools
        #  for response validation, and use the Lark parser for filter validation
        self.session = requests.Session()
        self._timeout = 10  # seconds

        if alias_or_structure_resource_url in self.aliases:
            self.resource = self.aliases[alias_or_structure_resource_url]
        else:
            self.resource = alias_or_structure_resource_url

    @staticmethod
    def _build_filter(
        elements=None, nelements=None, nsites=None, chemical_formula_anonymous=None, chemical_formula_hill=None
    ):
        """
        Convenience method to build an OPTIMADE filter.
        """

        filters = []

        if elements:
            if isinstance(elements, str):
                elements = [elements]
            elements_str = ", ".join([f'"{el}"' for el in elements])
            filters.append(f"(elements HAS ALL {elements_str})")

        if nsites:
            if isinstance(nsites, (list, tuple)):
                filters.append(f"(nsites>={min(nsites)} AND nsites<={max(nsites)})")
            else:
                filters.append(f"(nsites={int(nsites)})")

        if nelements:
            if isinstance(nelements, (list, tuple)):
                filters.append(f"(nelements>={min(nelements)} AND nelements<={max(nelements)})")
            else:
                filters.append(f"(nelements={int(nelements)})")

        if chemical_formula_anonymous:
            filters.append(f'(chemical_formula_anonymous="{chemical_formula_anonymous}")')

        if chemical_formula_hill:
            filters.append(f'(chemical_formula_hill="{chemical_formula_anonymous}")')

        return " AND ".join(filters)

    def get_structures(
        self,
        elements=None,
        nelements=None,
        nsites=None,
        chemical_formula_anonymous=None,
        chemical_formula_hill=None,
    ) -> Dict[str, Structure]:
        """
        Retrieve structures from the OPTIMADE database.

        Not all functionality of OPTIMADE is currently exposed in this convenience method. To
        use a custom filter, call get_structures_with_filter().

        Args:
            elements: List of elements
            nelements: Number of elements, e.g. 4 or [2, 5] for the range >=2 and <=5
            nsites: Number of sites, e.g. 4 or [2, 5] for the range >=2 and <=5
            chemical_formula_anonymous: Anonymous chemical formula
            chemical_formula_hill: Chemical formula following Hill convention

        Returns: Dict of Structures keyed by that database's id system
        """

        optimade_filter = self._build_filter(
            elements=elements,
            nelements=nelements,
            nsites=nsites,
            chemical_formula_anonymous=chemical_formula_anonymous,
            chemical_formula_hill=chemical_formula_hill,
        )

        return self.get_structures_with_filter(optimade_filter)

    def get_structures_with_filter(self, optimade_filter: str) -> Dict[str, Structure]:
        """
        Get structures satisfying a given OPTIMADE filter.

        Args:
            filter: An OPTIMADE-compliant filter

        Returns: Dict of Structures keyed by that database's id system
        """

        fields = "response_fields=lattice_vectors,cartesian_site_positions,species,species_at_sites"

        url = f"{self.resource}/v1/structures?filter={optimade_filter}&fields={fields}"

        json = self.session.get(url, timeout=self._timeout).json()

        structures = self._get_structures_from_resource(json)

        if "next" in json["links"] and json["links"]["next"]:
            pbar = PBar(total=json["meta"].get("data_returned"))
            while "next" in json["links"] and json["links"]["next"]:
                json = self.session.get(json["links"]["next"], timeout=self._timeout).json()
                structures.update(self._get_structures_from_resource(json))
                pbar.update(len(structures))

        return structures

    @staticmethod
    def _get_structures_from_resource(json):

        structures = {}

        def _sanitize_symbol(symbol):
            if symbol == "vacancy":
                symbol = DummySpecies("X_vacancy", oxidation_state=None)
            elif symbol == "X":
                symbol = DummySpecies("X", oxidation_state=None)
            return symbol

        def _get_comp(sp_dict):
            return {
                _sanitize_symbol(symbol): conc
                for symbol, conc in zip(sp_dict["chemical_symbols"], sp_dict["concentration"])
            }

        for data in json["data"]:

            # TODO: check the spec! and remove this try/except (are all providers following spec?)

            try:
                # e.g. COD
                structure = Structure(
                    lattice=data["attributes"]["lattice_vectors"],
                    species=[_get_comp(d) for d in data["attributes"]["species"]],
                    coords=data["attributes"]["cartesian_site_positions"],
                    coords_are_cartesian=True,
                )
                structures[data["id"]] = structure
            except Exception:

                try:
                    # e.g. MP (all ordered, no vacancies)
                    structure = Structure(
                        lattice=data["attributes"]["lattice_vectors"],
                        species=data["attributes"]["species_at_sites"],
                        coords=data["attributes"]["cartesian_site_positions"],
                        coords_are_cartesian=True,
                    )
                    structures[data["id"]] = structure
                except Exception:
                    pass

        return structures

    def refresh_aliases(self, providers_url="https://providers.optimade.org/providers.json"):
        """
        Updates available OPTIMADE structure resources based on the current list of OPTIMADE
        providers.
        """
        json = self.session.get(url=providers_url, timeout=self._timeout).json()
        providers_from_url = {
            entry["id"]: entry["attributes"]["base_url"] for entry in json["data"] if entry["attributes"]["base_url"]
        }

        providers = {}
        for provider, link in providers_from_url.items():
            try:
                providers[provider] = self.session.get(f"{link}/v1/links", timeout=self._timeout).json()
            except Exception as exc:
                print(f"Failed to parse {provider} at {link}: {exc}")

        # TODO: importing optimade-python-tool's data structures will make more sense
        Provider = namedtuple("Provider", ["name", "base_url", "description", "homepage"])

        def _parse_provider_link(provider, provider_link_json):
            """No validation attempted."""
            ps = {}
            try:
                d = [d for d in provider_link_json["data"] if d["attributes"]["link_type"] == "child"]
                for link in d:
                    key = f"{provider}.{link['id']}" if provider != link["id"] else provider
                    if link["attributes"]["base_url"]:
                        ps[key] = Provider(
                            name=link["attributes"]["name"],
                            base_url=link["attributes"]["base_url"],
                            description=link["attributes"]["description"],
                            homepage=link["attributes"].get("homepage"),
                        )
            except Exception:
                # print(f"Failed to parse {provider}: {exc}")
                # Not all providers parse yet.
                pass
            return ps

        structure_providers = {}
        for provider, provider_link_json in providers.items():
            structure_providers.update(_parse_provider_link(provider, provider_link_json))

        self.aliases = {alias: provider.base_url for alias, provider in structure_providers.items()}

    def __enter__(self):
        """
        Support for "with" context.
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Support for "with" context.
        """
        self.session.close()
