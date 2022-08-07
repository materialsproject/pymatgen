# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Determine functional groups present in a Molecule.
"""


import copy

from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.core.structure import Molecule
from pymatgen.io.babel import BabelMolAdaptor

try:
    import networkx as nx
    import networkx.algorithms.isomorphism as iso
except ImportError:
    raise ImportError("pymatgen.analysis.functional_groups requires the NetworkX graph library to be installed.")

__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "July 2018"
__credit__ = "Peiyuan Yu"


class FunctionalGroupExtractor:
    """
    This class is used to algorithmically parse a molecule (represented by an
    instance of pymatgen.analysis.graphs.MoleculeGraph) and determine arbitrary
    functional groups.
    """

    def __init__(self, molecule, optimize=False):
        """
        Instantiation method for FunctionalGroupExtractor.

        :param molecule: Either a filename, a pymatgen.core.structure.Molecule
            object, or a pymatgen.analysis.graphs.MoleculeGraph object.
        :param optimize: Default False. If True, then the input molecule will be
            modified, adding Hydrogens, performing a simple conformer search,
            etc.
        """

        self.molgraph = None

        if isinstance(molecule, str):
            try:
                if optimize:
                    obmol = BabelMolAdaptor.from_file(molecule, file_format="mol")
                    # OBMolecule does not contain pymatgen Molecule information
                    # So, we need to wrap the obmol in a BabelMolAdapter
                    obmol.add_hydrogen()
                    obmol.make3d()
                    obmol.localopt()
                    self.molecule = obmol.pymatgen_mol
                else:
                    self.molecule = Molecule.from_file(molecule)
            except OSError:
                raise ValueError("Input must be a valid molecule file, a Molecule object, or a MoleculeGraph object.")

        elif isinstance(molecule, Molecule):
            if optimize:
                obmol = BabelMolAdaptor(molecule)
                obmol.add_hydrogen()
                obmol.make3d()
                obmol.localopt()

                self.molecule = obmol.pymatgen_mol
            else:
                self.molecule = molecule

        elif isinstance(molecule, MoleculeGraph):
            if optimize:
                obmol = BabelMolAdaptor(molecule.molecule)
                obmol.add_hydrogen()
                obmol.make3d()
                obmol.localopt()

                self.molecule = obmol.pymatgen_mol

            else:
                self.molecule = molecule.molecule
                self.molgraph = molecule

        else:
            raise ValueError("Input to FunctionalGroupExtractor must be str, Molecule, or MoleculeGraph.")

        if self.molgraph is None:
            self.molgraph = MoleculeGraph.with_local_env_strategy(self.molecule, OpenBabelNN())

        # Assign a specie and coordinates to each node in the graph,
        # corresponding to the Site in the Molecule object
        self.molgraph.set_node_attributes()

        self.species = nx.get_node_attributes(self.molgraph.graph, "specie")

    def get_heteroatoms(self, elements=None):
        """
        Identify non-H, non-C atoms in the MoleculeGraph, returning a list of
        their node indices.

        :param elements: List of elements to identify (if only certain
            functional groups are of interest).
        :return: set of ints representing node indices
        """

        heteroatoms = set()

        for node in self.molgraph.graph.nodes():
            if elements is not None:
                if str(self.species[node]) in elements:
                    heteroatoms.add(node)
            else:
                if str(self.species[node]) not in ["C", "H"]:
                    heteroatoms.add(node)

        return heteroatoms

    def get_special_carbon(self, elements=None):
        """
        Identify Carbon atoms in the MoleculeGraph that fit the characteristics
        defined Ertl (2017), returning a list of their node indices.

        The conditions for marking carbon atoms are (quoted from Ertl):
            "- atoms connected by non-aromatic double or triple bond to any
            heteroatom
            - atoms in nonaromatic carbon-carbon double or triple bonds
            - acetal carbons, i.e. sp3 carbons connected to two or more oxygens,
            nitrogens or sulfurs; these O, N or S atoms must have only single bonds
            - all atoms in oxirane, aziridine and thiirane rings"

        :param elements: List of elements that will qualify a carbon as special
            (if only certain functional groups are of interest).
            Default None.
        :return: set of ints representing node indices
        """

        specials = set()

        # For this function, only carbons are considered
        carbons = [n for n in self.molgraph.graph.nodes if str(self.species[n]) == "C"]

        # Condition one: double/triple bonds to heteroatoms
        for node in carbons:
            neighbors = self.molgraph.graph[node]

            for neighbor, attributes in neighbors.items():
                if elements is not None:
                    if str(self.species[neighbor]) in elements and int(attributes[0]["weight"]) in [2, 3]:
                        specials.add(node)
                else:
                    if str(self.species[neighbor]) not in ["C", "H"] and int(attributes[0]["weight"]) in [2, 3]:
                        specials.add(node)

        # Condition two: carbon-carbon double & triple bonds
        for node in carbons:
            neighbors = self.molgraph.graph[node]

            for neighbor, attributes in neighbors.items():
                if str(self.species[neighbor]) == "C" and int(attributes[0]["weight"]) in [2, 3]:
                    specials.add(node)
                    specials.add(neighbor)

        # Condition three: Acetal carbons
        for node in carbons:
            neighbors = self.molgraph.graph[node]

            neighbor_spec = [str(self.species[n]) for n in neighbors]

            ons = len([n for n in neighbor_spec if n in ["O", "N", "S"]])

            if len(neighbors) == 4 and ons >= 2:
                specials.add(node)

        # Condition four: oxirane/aziridine/thiirane rings
        rings = self.molgraph.find_rings()
        rings_indices = [set(sum(ring, ())) for ring in rings]

        for ring in rings_indices:
            ring_spec = sorted(str(self.species[node]) for node in ring)
            # All rings of interest are three-member rings
            if len(ring) == 3 and ring_spec in [
                ["C", "C", "O"],
                ["C", "C", "N"],
                ["C", "C", "S"],
            ]:
                for node in ring:
                    if node in carbons:
                        specials.add(node)

        return specials

    def link_marked_atoms(self, atoms):
        """
        Take a list of marked "interesting" atoms (heteroatoms, special carbons)
        and attempt to connect them, returning a list of disjoint groups of
        special atoms (and their connected hydrogens).

        :param atoms: set of marked "interesting" atoms, presumably identified
            using other functions in this class.
        :return: list of sets of ints, representing groups of connected atoms
        """

        # We will add hydrogens to functional groups
        hydrogens = {n for n in self.molgraph.graph.nodes if str(self.species[n]) == "H"}

        # Graph representation of only marked atoms
        subgraph = self.molgraph.graph.subgraph(list(atoms)).to_undirected()

        func_grps = []
        for func_grp in nx.connected_components(subgraph):
            grp_hs = set()
            for node in func_grp:
                neighbors = self.molgraph.graph[node]
                for neighbor in neighbors:
                    # Add all associated hydrogens into the functional group
                    if neighbor in hydrogens:
                        grp_hs.add(neighbor)
            func_grp = func_grp.union(grp_hs)

            func_grps.append(func_grp)

        return func_grps

    def get_basic_functional_groups(self, func_groups=None):
        """
        Identify functional groups that cannot be identified by the Ertl method
        of get_special_carbon and get_heteroatoms, such as benzene rings, methyl
        groups, and ethyl groups.

        TODO: Think of other functional groups that are important enough to be
        added (ex: do we need ethyl, butyl, propyl?)

        :param func_groups: List of strs representing the functional groups of
            interest. Default to None, meaning that all of the functional groups
            defined in this function will be sought.
        :return: list of sets of ints, representing groups of connected atoms
        """

        strat = OpenBabelNN()

        hydrogens = {n for n in self.molgraph.graph.nodes if str(self.species[n]) == "H"}

        carbons = [n for n in self.molgraph.graph.nodes if str(self.species[n]) == "C"]

        if func_groups is None:
            func_groups = ["methyl", "phenyl"]

        results = []

        if "methyl" in func_groups:
            for node in carbons:
                neighbors = strat.get_nn_info(self.molecule, node)
                hs = {n["site_index"] for n in neighbors if n["site_index"] in hydrogens}
                # Methyl group is CH3, but this will also catch methane
                if len(hs) >= 3:
                    hs.add(node)
                    results.append(hs)

        if "phenyl" in func_groups:
            rings_indices = [set(sum(ring, ())) for ring in self.molgraph.find_rings()]

            possible_phenyl = [r for r in rings_indices if len(r) == 6]

            for ring in possible_phenyl:
                # Phenyl group should have only one (0 for benzene) member whose
                # neighbors are not two carbons and one hydrogen
                num_deviants = 0
                for node in ring:
                    neighbors = strat.get_nn_info(self.molecule, node)
                    neighbor_spec = sorted(str(self.species[n["site_index"]]) for n in neighbors)
                    if neighbor_spec != ["C", "C", "H"]:
                        num_deviants += 1

                if num_deviants <= 1:
                    for node in ring:
                        ring_group = copy.deepcopy(ring)
                        neighbors = self.molgraph.graph[node]

                        # Add hydrogens to the functional group
                        for neighbor in neighbors:
                            if neighbor in hydrogens:
                                ring_group.add(neighbor)

                    results.append(ring_group)

        return results

    def get_all_functional_groups(self, elements=None, func_groups=None, catch_basic=True):
        """
        Identify all functional groups (or all within a certain subset) in the
        molecule, combining the methods described above.

        :param elements: List of elements that will qualify a carbon as special
            (if only certain functional groups are of interest).
            Default None.
        :param func_groups: List of strs representing the functional groups of
            interest. Default to None, meaning that all of the functional groups
            defined in this function will be sought.
        :param catch_basic: bool. If True, use get_basic_functional_groups and
            other methods
        :return: list of sets of ints, representing groups of connected atoms
        """

        heteroatoms = self.get_heteroatoms(elements=elements)
        special_cs = self.get_special_carbon(elements=elements)
        groups = self.link_marked_atoms(heteroatoms.union(special_cs))

        if catch_basic:
            groups += self.get_basic_functional_groups(func_groups=func_groups)

        return groups

    def categorize_functional_groups(self, groups):
        """
        Determine classes of functional groups present in a set.

        :param groups: Set of functional groups.
        :return: dict containing representations of the groups, the indices of
            where the group occurs in the MoleculeGraph, and how many of each
            type of group there is.
        """

        categories = {}

        em = iso.numerical_edge_match("weight", 1)  # pylint: disable=E1102
        nm = iso.categorical_node_match("specie", "C")

        for group in groups:
            atoms = [self.molecule[a] for a in group]
            species = [a.specie for a in atoms]
            coords = [a.coords for a in atoms]

            adaptor = BabelMolAdaptor(Molecule(species, coords))
            # Use Canonical SMILES to ensure uniqueness
            smiles = adaptor.pybel_mol.write("can").strip()

            if smiles in categories:
                this_subgraph = self.molgraph.graph.subgraph(list(group)).to_undirected()
                for other in categories[smiles]["groups"]:
                    other_subgraph = self.molgraph.graph.subgraph(list(other)).to_undirected()

                    if not nx.is_isomorphic(this_subgraph, other_subgraph, edge_match=em, node_match=nm):
                        break

                    if group not in categories[smiles]["groups"]:
                        categories[smiles]["groups"].append(group)
                        categories[smiles]["count"] += 1

            else:
                categories[smiles] = {"groups": [group], "count": 1}

        return categories
