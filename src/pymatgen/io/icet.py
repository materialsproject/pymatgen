from __future__ import annotations

import multiprocessing as multiproc
import warnings
from string import ascii_uppercase
from time import time
from typing import TYPE_CHECKING

from pymatgen.command_line.mcsqs_caller import Sqs
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

try:
    from icet import ClusterSpace
    from icet.tools import enumerate_structures
    from icet.tools.structure_generation import _get_sqs_cluster_vector, _validate_concentrations, generate_sqs
    from mchammer.calculators import compare_cluster_vectors
except ImportError:
    ClusterSpace = _validate_concentrations = _get_sqs_cluster_vector = compare_cluster_vectors = generate_sqs = (
        enumerate_structures
    ) = None


if TYPE_CHECKING:
    from typing import Any, ClassVar

    from _icet import _ClusterSpace
    from ase import Atoms


class IcetSQS:
    """Interface to the Icet library of SQS structure generation tools.

    https://icet.materialsmodeling.org
    """

    sqs_kwarg_names: ClassVar[dict[str, tuple[str, ...]]] = {
        "monte_carlo": (
            "include_smaller_cells",
            "pbc",
            "T_start",
            "T_stop",
            "n_steps",
            "optimality_weight",
            "random_seed",
            "tol",
        ),
        "enumeration": ("include_smaller_cells", "pbc", "optimality_weight", "tol"),
    }
    _sqs_kwarg_defaults: ClassVar[dict[str, Any]] = {
        "optimality_weight": None,
        "tol": 1.0e-5,
        "include_smaller_cells": False,  # for consistency with ATAT
        "pbc": (True, True, True),
    }
    sqs_methods: tuple[str, ...] = ("enumeration", "monte_carlo")

    def __init__(
        self,
        structure: Structure,
        scaling: int,
        instances: int | None,
        cluster_cutoffs: dict[int, float],
        sqs_method: str | None = None,
        sqs_kwargs: dict | None = None,
    ) -> None:
        """
        Instantiate an IcetSQS interface.

        Args:
            structure (Structure): disordered structure to compute SQS
            scaling (int): SQS supercell contains scaling * len(structure) sites
            instances (int): number of parallel SQS jobs to run
            cluster_cutoffs (dict): dict of cluster size (pairs, triplets, ...) and
                the size of the cluster
        Kwargs:
            sqs_method (str or None): if a str, one of ("enumeration", "monte_carlo")
                If None, default to "enumeration" for a supercell of < 24 sites, and
                "monte carlo" otherwise.
            sqs_kwargs (dict): kwargs to pass to the icet SQS generators.
                See self.sqs_kwarg_names for possible options.
        """
        if ClusterSpace is None:
            raise ImportError("IcetSQS requires the icet package. Use `pip install icet`")

        self._structure = structure
        self.scaling = scaling
        self.instances = instances or multiproc.cpu_count()

        self._get_site_composition()

        # The peculiar way that icet works requires a copy of the
        # disordered structure, but without any fractionally-occupied sites
        # Essentially the host structure
        _ordered_structure = structure.copy()

        original_composition = _ordered_structure.composition.as_dict()
        dummy_comp = next(iter(_ordered_structure.composition))
        _ordered_structure.replace_species(
            {species: dummy_comp for species in original_composition if species != dummy_comp}
        )
        self._ordered_atoms = AseAtomsAdaptor.get_atoms(_ordered_structure)

        self.cutoffs_list = []
        for i in range(2, max(cluster_cutoffs.keys()) + 1):
            if i not in cluster_cutoffs:
                # pad missing non-sequential values
                cluster_cutoffs[i] = 0.0
            self.cutoffs_list.append(cluster_cutoffs[i])

        # For safety, enumeration works well on 1 core for ~< 24 sites/cell
        # The bottleneck is **generation** of the structures via enumeration,
        # less checking their SQS objective.
        # Beyond ~24 sites/cell, monte carlo is more efficient
        sqs_method = sqs_method or ("enumeration" if self.scaling * len(self._structure) < 24 else "monte_carlo")

        # Default sqs_kwargs
        self.sqs_kwargs = self._sqs_kwarg_defaults.copy()
        self.sqs_kwargs.update(sqs_kwargs or {})

        unrecognized_kwargs = {key for key in self.sqs_kwargs if key not in self.sqs_kwarg_names[sqs_method]}
        if len(unrecognized_kwargs) > 0:
            warnings.warn(f"Ignoring unrecognized icet {sqs_method} kwargs: {', '.join(unrecognized_kwargs)}")

        self.sqs_kwargs = {
            key: value for key, value in self.sqs_kwargs.items() if key in self.sqs_kwarg_names[sqs_method]
        }

        if sqs_method == "monte_carlo":
            self.sqs_getter = self.monte_carlo_sqs_structures
            if self.sqs_kwargs.get("random_seed") is None:
                self.sqs_kwargs["random_seed"] = int(1e6 * time())

        elif sqs_method == "enumeration":
            self.sqs_getter = self.enumerate_sqs_structures

        else:
            raise ValueError(f"Unknown {sqs_method=}! Must be one of {self.sqs_methods}")

        self._sqs_obj_kwargs = {}
        for key in ("optimality_weight", "tol"):
            if value := self.sqs_kwargs.get(key, self._sqs_kwarg_defaults[key]):
                self._sqs_obj_kwargs[key] = value

        cluster_space = self._get_cluster_space()
        self.target_concentrations = _validate_concentrations(
            concentrations=self.composition, cluster_space=cluster_space
        )
        self.sqs_vector = _get_sqs_cluster_vector(
            cluster_space=cluster_space,
            target_concentrations=self.target_concentrations,
        )

    def run(self) -> Sqs:
        """
        Run the SQS search with icet.

        Returns:
            pymatgen Sqs object
        """
        sqs_structures = self.sqs_getter()
        for idx in range(len(sqs_structures)):
            sqs_structures[idx]["structure"] = AseAtomsAdaptor.get_structure(sqs_structures[idx]["structure"])
        sqs_structures = sorted(sqs_structures, key=lambda entry: entry["objective_function"])

        return Sqs(
            bestsqs=sqs_structures[0]["structure"],
            objective_function=sqs_structures[0]["objective_function"],
            allsqs=sqs_structures,
            directory="./",
            clusters=str(self._get_cluster_space()),
        )

    def _get_site_composition(self) -> dict[str, dict]:
        """Get Icet-format composition from structure.

        Returns:
            Dict with sublattice compositions specified by uppercase letters,
                e.g. In_x Ga_1-x As becomes: {
                    "A": {"In": x, "Ga": 1 - x},
                    "B": {"As": 1}
                }
        """
        uppercase_letters = list(ascii_uppercase)
        self.composition: dict[str, dict] = {}
        for idx, site in enumerate(self._structure):
            site_comp = site.species.as_dict()
            if site_comp not in self.composition.values():
                self.composition[uppercase_letters[idx]] = site_comp

        return self.composition

    def _get_cluster_space(self) -> ClusterSpace:
        """Generate the ClusterSpace object for icet."""
        chemical_symbols = [list(site.species.as_dict()) for site in self._structure]
        return ClusterSpace(
            structure=self._ordered_atoms,
            cutoffs=self.cutoffs_list,
            chemical_symbols=chemical_symbols,
        )

    def get_icet_sqs_obj(self, material: Atoms | Structure, cluster_space: _ClusterSpace | None = None) -> float:
        """Get the SQS objective function.

        Args:
            material (pymatgen.Structure | ase.Atoms): structure to compute SQS objective function for.
            cluster_space (ClusterSpace): ClusterSpace of the SQS search.

        Returns:
            float: the SQS objective function
        """
        if isinstance(material, Structure):
            material = AseAtomsAdaptor.get_atoms(material)

        cluster_space = cluster_space or self._get_cluster_space()
        return compare_cluster_vectors(
            cv_1=cluster_space.get_cluster_vector(material),
            cv_2=self.sqs_vector,
            orbit_data=cluster_space.orbit_data,
            **self._sqs_obj_kwargs,
        )

    def enumerate_sqs_structures(self, cluster_space: _ClusterSpace | None = None) -> list:
        """Generate an SQS by enumeration of all possible arrangements.

        Adapted from icet.tools.structure_generation.generate_sqs_by_enumeration
        to accommodate multiprocessing.

        Args:
            cluster_space (ClusterSpace) : ClusterSpace of the SQS search.

        Returns:
            list: dicts of the form: {
                    "structure": SQS structure,
                    "objective_function": SQS objective function,
                }
        """
        # Translate concentrations to the format required for concentration
        # restricted enumeration
        cr: dict[str, tuple] = {}
        cluster_space = cluster_space or self._get_cluster_space()
        sub_lattices = cluster_space.get_sublattices(cluster_space.primitive_structure)
        for sl in sub_lattices:
            mult_factor = len(sl.indices) / len(cluster_space.primitive_structure)
            if sl.symbol in self.target_concentrations:
                sl_conc = self.target_concentrations[sl.symbol]
            else:
                sl_conc = {sl.chemical_symbols[0]: 1.0}
            for species, value in sl_conc.items():
                c = value * mult_factor
                if species in cr:
                    cr[species] = (cr[species][0] + c, cr[species][1] + c)
                else:
                    cr[species] = (c, c)

        # Check to be sure...
        c_sum = sum(c[0] for c in cr.values())
        if abs(c_sum - 1) >= self.sqs_kwargs["tol"]:
            raise ValueError(f"Site occupancies sum to {abs(c_sum - 1)} instead of 1!")

        sizes = list(range(1, self.scaling + 1)) if self.sqs_kwargs["include_smaller_cells"] else [self.scaling]

        # Prepare primitive structure with the right boundary conditions
        prim = cluster_space.primitive_structure
        prim.set_pbc(self.sqs_kwargs["pbc"])

        structures = enumerate_structures(prim, sizes, cluster_space.chemical_symbols, concentration_restrictions=cr)
        chunks: list[list[Atoms]] = [[] for _ in range(self.instances)]
        proc_idx = 0
        for structure in structures:
            chunks[proc_idx].append(structure)
            proc_idx = (proc_idx + 1) % self.instances

        manager = multiproc.Manager()
        working_list = manager.list()
        processes = []
        for proc_idx in range(self.instances):
            process = multiproc.Process(
                target=self._get_best_sqs_from_list,
                args=(chunks[proc_idx], working_list),
            )
            processes.append(process)
            process.start()

        for process in processes:
            process.join()

        return list(working_list)

    def _get_best_sqs_from_list(self, structures: list[Atoms], output_list: list[dict]) -> dict[str, Any]:
        """Find best SQS structure from list of SQS structures.

        Args:
            structures (list of ase Atoms) : list of SQS structures
            output_list (list of dicts) : shared list between
                multiprocessing processes to store best SQS objects.
        """
        best_sqs: dict[str, Any] = {"structure": None, "objective_function": 1.0e20}
        cluster_space = self._get_cluster_space()
        for structure in structures:
            objective = self.get_icet_sqs_obj(structure, cluster_space=cluster_space)
            if objective < best_sqs["objective_function"]:
                best_sqs = {"structure": structure, "objective_function": objective}
        output_list.append(best_sqs)

        return best_sqs

    def _single_monte_carlo_sqs_run(self):
        """Run a single Monte Carlo SQS search with Icet."""
        cluster_space = self._get_cluster_space()
        sqs_structure = generate_sqs(
            cluster_space=cluster_space,
            max_size=self.scaling,
            target_concentrations=self.target_concentrations,
            **self.sqs_kwargs,
        )
        return {
            "structure": sqs_structure,
            "objective_function": self.get_icet_sqs_obj(sqs_structure, cluster_space=cluster_space),
        }

    def monte_carlo_sqs_structures(self) -> list:
        """Run `self.instances` Monte Carlo SQS search with Icet."""
        with multiproc.Pool(self.instances) as pool:
            return pool.starmap(self._single_monte_carlo_sqs_run, [() for _ in range(self.instances)])
