"""Optimade support."""

from __future__ import annotations

import logging
import sys
from typing import TYPE_CHECKING, NamedTuple
from urllib.parse import urljoin, urlparse

import requests
from tqdm import tqdm

from pymatgen.core import DummySpecies, Structure
from pymatgen.util.due import Doi, due
from pymatgen.util.provenance import StructureNL

if TYPE_CHECKING:
    from typing import ClassVar

    from typing_extensions import Self

    from pymatgen.core.structure import Molecule

_logger = logging.getLogger(__name__)
_handler = logging.StreamHandler(sys.stdout)
_logger.addHandler(_handler)
_logger.setLevel(logging.WARNING)


class Provider(NamedTuple):
    """TODO: Import optimade-python-tool's data structures will make more sense."""

    name: str
    base_url: str
    description: str
    homepage: str
    prefix: str


@due.dcite(
    Doi("10.1038/s41597-021-00974-z"),
    description="OPTIMADE, an API for exchanging materials data",
)
class OptimadeRester:
    """Call OPTIMADE-compliant APIs, see https://optimade.org and [1].

    This class is ready to use but considered in-development and subject to change.

    Please also consider using the client in "OPTIMADE Python tools":

    https://www.optimade.org/optimade-python-tools/latest/getting_started/client/

    The "OPTIMADE Python tools" client is less integrated with pymatgen, but
    more actively developed for the latest OPTIMADE features.

    [1] Andersen, C.W., *et al*.
        OPTIMADE, an API for exchanging materials data.
        Sci Data 8, 217 (2021). https://doi.org/10.1038/s41597-021-00974-z
    """

    # regenerate on-demand from official providers.json using OptimadeRester.refresh_aliases()
    # these aliases are provided as a convenient shortcut for users of the OptimadeRester class
    aliases: ClassVar[dict[str, str]] = {
        "aflow": "https://aflow.org/API/optimade/",
        "alexandria": "https://alexandria.icams.rub.de/pbe",
        "alexandria.pbe": "https://alexandria.icams.rub.de/pbe",
        "alexandria.pbesol": "https://alexandria.icams.rub.de/pbesol",
        "cod": "https://www.crystallography.net/cod/optimade",
        "cmr": "https://cmr-optimade.fysik.dtu.dk",
        "mcloud.mc3d": "https://aiida.materialscloud.org/mc3d/optimade",
        "mcloud.mc2d": "https://aiida.materialscloud.org/mc2d/optimade",
        "mcloud.2dtopo": "https://aiida.materialscloud.org/2dtopo/optimade",
        "mcloud.tc-applicability": "https://aiida.materialscloud.org/tc-applicability/optimade",
        "mcloud.pyrene-mofs": "https://aiida.materialscloud.org/pyrene-mofs/optimade",
        "mcloud.curated-cofs": "https://aiida.materialscloud.org/curated-cofs/optimade",
        "mcloud.stoceriaitf": "https://aiida.materialscloud.org/stoceriaitf/optimade",
        "mcloud.scdm": "https://aiida.materialscloud.org/autowannier/optimade",
        "mcloud.tin-antimony-sulfoiodide": "https://aiida.materialscloud.org/tin-antimony-sulfoiodide/optimade",
        "mcloud.optimade-sample": "https://aiida.materialscloud.org/optimade-sample/optimade",
        "mp": "https://optimade.materialsproject.org",
        "mpdd": "http://mpddoptimade.phaseslab.org",
        "mpds": "https://api.mpds.io",
        "nmd": "https://nomad-lab.eu/prod/rae/optimade/",
        "odbx": "https://optimade.odbx.science",
        "odbx.odbx_misc": "https://optimade-misc.odbx.science",
        "odbx.gnome": "https://optimade-gnome.odbx.science",
        "omdb.omdb_production": "https://optimade.openmaterialsdb.se",
        "oqmd": "https://oqmd.org/optimade/",
        "jarvis": "https://jarvis.nist.gov/optimade/jarvisdft",
        "tcod": "https://www.crystallography.net/tcod/optimade",
        "twodmatpedia": "http://optimade.2dmatpedia.org",
    }

    # The set of OPTIMADE fields that are required to define a `pymatgen.core.Structure`
    mandatory_response_fields = (
        "lattice_vectors",
        "cartesian_site_positions",
        "species",
        "species_at_sites",
    )

    def __init__(
        self,
        aliases_or_resource_urls: str | list[str] | None = None,
        refresh_aliases: bool = False,
        timeout: int = 5,
    ):
        """OPTIMADE is an effort to provide a standardized interface to retrieve information
        from many different materials science databases.

        This is a client to retrieve structures from OPTIMADE v1 compliant endpoints. It
        does not yet support all features of the OPTIMADE v1 specification but is intended
        as a way to quickly search an endpoint in a way familiar to users of pymatgen without
        needing to know the full OPTIMADE specification.

        For advanced usage, please see the OPTIMADE documentation at optimade.org and
        consider calling the APIs directly.

        For convenience, known OPTIMADE endpoints have been given aliases in pymatgen to save
        typing the full URL.

        To get an up-to-date list aliases, generated from the current list of OPTIMADE providers
        at optimade.org, call the refresh_aliases() method or pass refresh_aliases=True when
        creating instances of this class.

        Args:
            aliases_or_resource_urls: the alias or structure resource URL or a list of
            aliases or resource URLs, if providing the resource URL directly it should not
            be an index, this interface can only currently access the "v1/structures"
            information from the specified resource URL
            refresh_aliases: if True, use an up-to-date list of providers/aliases from the live
            list of OPTIMADE providers hosted at https://providers.optimade.org.
            timeout: number of seconds before an attempted request is abandoned, a good
            timeout is useful when querying many providers, some of which may be offline
        """
        # TODO: maybe we should use the nice pydantic models from optimade-python-tools
        #  for response validation, and use the Lark parser for filter validation
        self.session = requests.Session()
        self._timeout = timeout  # seconds

        # Optionally refresh the aliases before interpreting those provided by the user
        # or using potentially outdated set provided in the code
        if refresh_aliases:
            _logger.warning("Refreshing OPTIMADE provider aliases from https://providers.optimade.org")
            self.refresh_aliases()

        if isinstance(aliases_or_resource_urls, str):
            aliases_or_resource_urls = [aliases_or_resource_urls]

        # this stores a dictionary with keys provider id (in the same format as the aliases)
        # and values as the corresponding URL
        self.resources = {}

        # preprocess aliases to ensure they have a trailing slash where appropriate
        for alias, url in self.aliases.items():
            if urlparse(url).path is not None and not url.endswith("/"):
                self.aliases[alias] += "/"

        if not aliases_or_resource_urls:
            aliases_or_resource_urls = list(self.aliases)
            _logger.warning(
                "Connecting to all known OPTIMADE providers, this will be slow. Please connect to only the "
                f"OPTIMADE providers you want to query. Choose from: {', '.join(self.aliases)}"
            )

        for alias_or_resource_url in aliases_or_resource_urls:
            if alias_or_resource_url in self.aliases:
                self.resources[alias_or_resource_url] = self.aliases[alias_or_resource_url]

            elif self._validate_provider(alias_or_resource_url):
                # TODO: unclear what the key should be here, the "prefix" is for the root provider,
                # may need to walk back to the index for the given provider to find the correct identifier

                self.resources[alias_or_resource_url] = alias_or_resource_url

            else:
                _logger.error(f"The following is not a known alias or a valid url: {alias_or_resource_url}")

        self._providers = {url: self._validate_provider(provider_url=url) for url in self.resources.values()}

    def __repr__(self):
        return f"OptimadeRester connected to: {', '.join(self.resources.values())}"

    def __str__(self):
        return self.describe()

    def describe(self):
        """Human-readable information about the resources being searched by the OptimadeRester."""
        provider_text = "\n".join(map(str, (provider for provider in self._providers.values() if provider)))
        return f"OptimadeRester connected to:\n{provider_text}"

    def _get_json(self, url):
        """Retrieve and returns JSON resource from given url."""
        return self.session.get(url, timeout=self._timeout).json()

    @staticmethod
    def _build_filter(
        elements: str | list[str] | None = None,
        nelements: int | None = None,
        nsites: int | None = None,
        chemical_formula_anonymous: str | None = None,
        chemical_formula_hill: str | None = None,
    ):
        """Convenience method to build an OPTIMADE filter."""
        filters = []

        if elements:
            if isinstance(elements, str):
                elements = [elements]
            elements_str = ", ".join(f'"{el}"' for el in elements)
            filters.append(f"(elements HAS ALL {elements_str})")

        if nsites:
            if isinstance(nsites, list | tuple):
                filters.append(f"(nsites>={min(nsites)} AND nsites<={max(nsites)})")
            else:
                filters.append(f"({nsites=})")

        if nelements:
            if isinstance(nelements, list | tuple):
                filters.append(f"(nelements>={min(nelements)} AND nelements<={max(nelements)})")
            else:
                filters.append(f"({nelements=})")

        if chemical_formula_anonymous:
            filters.append(f"({chemical_formula_anonymous=})")

        if chemical_formula_hill:
            filters.append(f"({chemical_formula_hill=})")

        return " AND ".join(filters)

    def get_structures(
        self,
        elements: list[str] | str | None = None,
        nelements: int | None = None,
        nsites: int | None = None,
        chemical_formula_anonymous: str | None = None,
        chemical_formula_hill: str | None = None,
    ) -> dict[str, dict[str, Structure | Molecule]]:
        """Retrieve Structures from OPTIMADE providers.

        Not all functionality of OPTIMADE is currently exposed in this convenience method. To
        use a custom filter, call get_structures_with_filter().

        Args:
            elements: List of elements
            nelements: Number of elements, e.g. 4 or [2, 5] for the range >=2 and <=5
            nsites: Number of sites, e.g. 4 or [2, 5] for the range >=2 and <=5
            chemical_formula_anonymous: The desired chemical formula in OPTIMADE anonymous formula format
            (NB. The ordering is reversed from the pymatgen format, e.g. pymatgen "ABC2" should become "A2BC").
            chemical_formula_hill: The desired chemical formula in the OPTIMADE take on the Hill formula format.
            (NB. Again, this is different from the pymatgen format, as the OPTIMADE version is a reduced chemical
            formula simply using the IUPAC/Hill ordering.)

        Returns:
            dict[str, Structure]: keyed by that database provider's id system
        """
        optimade_filter = self._build_filter(
            elements=elements,
            nelements=nelements,
            nsites=nsites,
            chemical_formula_anonymous=chemical_formula_anonymous,
            chemical_formula_hill=chemical_formula_hill,
        )

        return self.get_structures_with_filter(optimade_filter)

    def get_snls(
        self,
        elements: list[str] | str | None = None,
        nelements: int | None = None,
        nsites: int | None = None,
        chemical_formula_anonymous: str | None = None,
        chemical_formula_hill: str | None = None,
        additional_response_fields: str | list[str] | set[str] | None = None,
    ) -> dict[str, dict[str, StructureNL]]:
        """Retrieve StructureNL from OPTIMADE providers.

        A StructureNL is an object provided by pymatgen which combines Structure with
        associated metadata, such as the URL is was downloaded from and any additional namespaced
        data.

        Not all functionality of OPTIMADE is currently exposed in this convenience method. To
        use a custom filter, call get_structures_with_filter().

        Args:
            elements: List of elements
            nelements: Number of elements, e.g. 4 or [2, 5] for the range >=2 and <=5
            nsites: Number of sites, e.g. 4 or [2, 5] for the range >=2 and <=5
            chemical_formula_anonymous: The desired chemical formula in OPTIMADE anonymous formula format
            (NB. The ordering is reversed from the pymatgen format, e.g. pymatgen "ABC2" should become "A2BC").
            chemical_formula_hill: The desired chemical formula in the OPTIMADE take on the Hill formula format.
            (NB. Again, this is different from the pymatgen format, as the OPTIMADE version is a reduced chemical
            formula simply using the IUPAC/Hill ordering.)
            additional_response_fields: Any additional fields desired from the OPTIMADE API,
            these will be stored under the `'_optimade'` key in each `StructureNL.data` dictionary.

        Returns:
            dict[str, StructureNL]: keyed by that database provider's id system
        """
        optimade_filter = self._build_filter(
            elements=elements,
            nelements=nelements,
            nsites=nsites,
            chemical_formula_anonymous=chemical_formula_anonymous,
            chemical_formula_hill=chemical_formula_hill,
        )

        return self.get_snls_with_filter(optimade_filter, additional_response_fields=additional_response_fields)

    def get_structures_with_filter(self, optimade_filter: str) -> dict[str, dict[str, Structure | Molecule]]:
        """Get structures satisfying a given OPTIMADE filter.

        Args:
            optimade_filter: An OPTIMADE-compliant filter

        Returns:
            dict[str, Structure]: keyed by that database provider's id system
        """
        all_snls = self.get_snls_with_filter(optimade_filter)
        all_structures = {}

        for identifier, snls_dict in all_snls.items():
            all_structures[identifier] = {k: snl.structure for k, snl in snls_dict.items()}

        return all_structures

    def get_snls_with_filter(
        self,
        optimade_filter: str,
        additional_response_fields: str | list[str] | set[str] | None = None,
    ) -> dict[str, dict[str, StructureNL]]:
        """Get structures satisfying a given OPTIMADE filter.

        Args:
            optimade_filter: An OPTIMADE-compliant filter
            additional_response_fields: Any additional fields desired from the OPTIMADE API,

        Returns:
            dict[str, Structure]: keyed by that database provider's id system
        """
        all_snls = {}

        response_fields = self._handle_response_fields(additional_response_fields)

        for identifier, resource in self.resources.items():
            url = urljoin(resource, f"v1/structures?filter={optimade_filter}&{response_fields=!s}")

            try:
                json = self._get_json(url)

                structures = self._get_snls_from_resource(json, url, identifier)

                pbar = tqdm(
                    total=json["meta"].get("data_returned", 0),
                    desc=identifier,
                    initial=len(structures),
                )

                # TODO: check spec for `more_data_available` boolean, may simplify this conditional
                while next_link := json.get("links", {}).get("next"):
                    if isinstance(next_link, dict) and "href" in next_link:
                        next_link = next_link["href"]
                    json = self._get_json(next_link)
                    additional_structures = self._get_snls_from_resource(json, url, identifier)
                    structures.update(additional_structures)
                    pbar.update(len(additional_structures))

                if structures:
                    all_snls[identifier] = structures

            except Exception:
                # TODO: manually inspect failures to either (a) correct a bug or (b) raise more appropriate error

                _logger.exception(f"Could not retrieve required information from provider {identifier} and {url=}")

        return all_snls

    @staticmethod
    def _get_snls_from_resource(json, url, identifier) -> dict[str, StructureNL]:
        snls = {}

        exceptions = set()

        def _sanitize_symbol(symbol):
            if symbol == "vacancy":
                symbol = DummySpecies("X_vacancy", oxidation_state=None)
            elif symbol == "X":
                symbol = DummySpecies("X", oxidation_state=None)
            return symbol

        def _get_comp(sp_dict):
            return {
                _sanitize_symbol(symbol): conc
                for symbol, conc in zip(sp_dict["chemical_symbols"], sp_dict["concentration"], strict=True)
            }

        for data in json["data"]:
            # TODO: check the spec! and remove this try/except (are all providers following spec?)
            # e.g. can check data["type"] == "structures"

            try:
                # e.g. COD
                structure = Structure(
                    lattice=data["attributes"]["lattice_vectors"],
                    species=[_get_comp(d) for d in data["attributes"]["species"]],
                    coords=data["attributes"]["cartesian_site_positions"],
                    coords_are_cartesian=True,
                )
                # Grab any custom fields or non-mandatory fields if they were requested
                namespaced_data = {
                    k: v
                    for k, v in data["attributes"].items()
                    if k.startswith("_") or k not in {"lattice_vectors", "species", "cartesian_site_positions"}
                }

                # TODO: follow `references` to add reference information here
                snl = StructureNL(
                    structure,
                    authors=[],
                    history=[
                        {
                            "name": identifier,
                            "url": url,
                            "description": {"id": data["id"]},
                        }
                    ],
                    data={"_optimade": namespaced_data},
                )

                snls[data["id"]] = snl

            # TODO: bare exception, remove...
            except Exception:
                try:
                    # e.g. MP (all ordered, no vacancies)
                    structure = Structure(
                        lattice=data["attributes"]["lattice_vectors"],
                        species=data["attributes"]["species_at_sites"],
                        coords=data["attributes"]["cartesian_site_positions"],
                        coords_are_cartesian=True,
                    )
                    # Grab any custom fields or non-mandatory fields if they were requested
                    namespaced_data = {
                        k: v
                        for k, v in data["attributes"].items()
                        if k.startswith("_")
                        or k
                        not in {
                            "lattice_vectors",
                            "species",
                            "cartesian_site_positions",
                        }
                    }

                    # TODO: follow `references` to add reference information here
                    snl = StructureNL(
                        structure,
                        authors=[],
                        history=[
                            {
                                "name": identifier,
                                "url": url,
                                "description": {"id": data["id"]},
                            }
                        ],
                        data={"_optimade": namespaced_data},
                    )

                    snls[data["id"]] = snl

                except Exception as exc:
                    if str(exc) not in exceptions:
                        exceptions.add(str(exc))

        if exceptions:
            _logger.error(f"Failed to parse returned data for {url}: {', '.join(exceptions)}")

        return snls

    def _validate_provider(self, provider_url) -> Provider | None:
        """Check that a given URL is indeed an OPTIMADE provider,
        returning None if it is not a provider, or the provider
        prefix if it is.

        TODO: careful reading of OPTIMADE specification required
        TODO: add better exception handling, intentionally permissive currently
        """
        # Add trailing slash to all URLs if missing; prevents urljoin from scrubbing
        # sections of the path
        if urlparse(provider_url).path is not None and not provider_url.endswith("/"):
            provider_url += "/"

        def is_url(url) -> bool:
            """Basic URL validation thanks to https://stackoverflow.com/a/52455972."""
            try:
                result = urlparse(url)
                return all([result.scheme, result.netloc])
            except ValueError:
                return False

        if not is_url(provider_url):
            _logger.warning(f"An invalid url was supplied: {provider_url}")
            return None

        url = None

        try:
            url = urljoin(provider_url, "v1/info")
            provider_info_json = self._get_json(url)
        except Exception as exc:
            _logger.warning(f"Failed to parse {url} when validating: {exc}")
            return None

        try:
            return Provider(
                name=provider_info_json["meta"].get("provider", {}).get("name", "Unknown"),
                base_url=provider_url,
                description=provider_info_json["meta"].get("provider", {}).get("description", "Unknown"),
                homepage=provider_info_json["meta"].get("provider", {}).get("homepage"),
                prefix=provider_info_json["meta"].get("provider", {}).get("prefix", "Unknown"),
            )
        except Exception as exc:
            _logger.warning(f"Failed to extract required information from {url}: {exc}")
            return None

    def _parse_provider(self, provider: str, provider_url: str) -> dict[str, Provider]:
        """Used internally to update the list of providers or to
        check a given URL is valid.

        It does not raise exceptions but will instead _logger.warning and provide
        an empty dictionary in the case of invalid data.

        In future, when the specification is sufficiently well adopted,
        we might be more strict here.

        Args:
            provider: the provider prefix
            provider_url: An OPTIMADE provider URL

        Returns:
            dict: keys (in format of "provider.database") mapped to Provider objects.
        """
        # Add trailing slash to all URLs if missing; prevents urljoin from scrubbing
        if urlparse(provider_url).path is not None and not provider_url.endswith("/"):
            provider_url += "/"

        url = None

        try:
            url = urljoin(provider_url, "v1/links")
            provider_link_json = self._get_json(url)
        except Exception:
            _logger.exception(f"Failed to parse {url} when following links")
            return {}

        def _parse_provider_link(provider, provider_link_json):
            """No validation attempted."""
            ps = {}
            try:
                data = [dct for dct in provider_link_json["data"] if dct["attributes"]["link_type"] == "child"]
                for link in data:
                    key = f"{provider}.{link['id']}" if provider != link["id"] else provider
                    if link["attributes"]["base_url"]:
                        ps[key] = Provider(
                            name=link["attributes"]["name"],
                            base_url=link["attributes"]["base_url"],
                            description=link["attributes"]["description"],
                            homepage=link["attributes"].get("homepage"),
                            prefix=link["attributes"].get("prefix"),
                        )
            except Exception:
                # print(f"Failed to parse {provider}: {exc}")
                # Not all providers parse yet.
                pass
            return ps

        return _parse_provider_link(provider, provider_link_json)

    def _handle_response_fields(self, additional_response_fields: str | list[str] | set[str] | None = None) -> str:
        """Used internally to handle the mandatory and additional response fields.

        Args:
            additional_response_fields: A set of additional fields to request.

        Returns:
            A string of comma-separated OPTIMADE response fields.
        """
        if isinstance(additional_response_fields, str):
            additional_response_fields = {additional_response_fields}
        if not additional_response_fields:
            additional_response_fields = set()
        return ",".join({*additional_response_fields, *self.mandatory_response_fields})

    def refresh_aliases(self, providers_url="https://providers.optimade.org/providers.json"):
        """Update available OPTIMADE structure resources based on the current list of OPTIMADE
        providers.
        """
        json = self._get_json(providers_url)
        providers_from_url = {
            entry["id"]: entry["attributes"]["base_url"] for entry in json["data"] if entry["attributes"]["base_url"]
        }

        structure_providers = {}
        for provider, provider_link in providers_from_url.items():
            structure_providers.update(self._parse_provider(provider, provider_link))

        self.aliases = {alias: provider.base_url for alias, provider in structure_providers.items()}

        # Add missing trailing slashes to any aliases with a path that need them
        for alias, url in self.aliases.items():
            if urlparse(url).path is not None and not url.endswith("/"):
                self.aliases[alias] += "/"

    # TODO: revisit context manager logic here and in MPRester
    def __enter__(self) -> Self:
        """Support for "with" context."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Support for "with" context."""
        self.session.close()
