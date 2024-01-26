from __future__ import annotations

import multiprocessing as multiproc
from importlib.metadata import PackageNotFoundError
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

    loaded_icet = True

except ImportError:
    loaded_icet = False


if TYPE_CHECKING:
    from typing import Any

    from ase import Atoms


class IcetSQS:
    sqs_kwarg_names: dict[str, tuple[str, ...]] = {
        "monte_carlo": ("T_start", "T_stop", "n_steps", "optimality_weight", "random_seed", "tol"),
        "enumeration": ("include_smaller_cells", "pbc", "optimality_weight", "tol"),
    }
    _sqs_kwarg_defaults: dict[str, Any] = {
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
        instances: int,
        sqs_method: str,
        cluster_cutoffs: dict[int, float],
        sqs_kwargs: dict | None = None,
    ) -> None:
        if not loaded_icet:
            raise PackageNotFoundError("IcetSQS requires the icet package. Use `pip install icet`.")

        self._structure = structure
        self.scaling = scaling
        self.instances = instances

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
        # Beyond that, monte carlo is more efficient
        sqs_method = sqs_method or ("enumeration" if self.scaling * len(self._structure) < 24 else "monte_carlo")

        # Default sqs_kwargs
        self.sqs_kwargs = self._sqs_kwarg_defaults.copy()
        self.sqs_kwargs.update(sqs_kwargs or {})
        self.sqs_kwargs = {k: v for k, v in self.sqs_kwargs.items() if k in self.sqs_kwarg_names[sqs_method]}

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
            cluster_space=cluster_space, target_concentrations=self.target_concentrations
        )

    def run(self) -> Sqs:
        sqs_structures = self.sqs_getter()
        for ientry in range(len(sqs_structures)):
            sqs_structures[ientry]["structure"] = AseAtomsAdaptor.get_structure(sqs_structures[ientry]["structure"])
        sqs_structures = sorted(sqs_structures, key=lambda entry: entry["objective_function"])

        return Sqs(
            bestsqs=sqs_structures[0]["structure"],
            objective_function=sqs_structures[0]["objective_function"],
            allsqs=sqs_structures,
            directory="./",
            clusters=str(self._get_cluster_space()),
        )

    def _get_site_composition(self) -> None:
        """
        Get Icet-format composition from structure.

        Sublattice compositions are specified with new uppercase letters,
        e.g., In_x Ga_1-x As becomes:
        {
            "A": {"In": x, "Ga": 1 - x},
            "B": {"As": 1}
        }
        """
        uppercase_letters = list(ascii_uppercase)
        iletter = 0
        self.composition = {}
        for isite in range(len(self._structure)):
            site_comp = self._structure.sites[isite].species.as_dict()
            if site_comp not in self.composition.values():
                self.composition[uppercase_letters[iletter]] = site_comp
                iletter += 1

    def _get_cluster_space(self) -> ClusterSpace:
        chemical_symbols = [
            list(self._structure.sites[isite].species.as_dict())
            for isite in range(self._structure.num_sites)
        ]
        return ClusterSpace(structure=self._ordered_atoms, cutoffs=self.cutoffs_list, chemical_symbols=chemical_symbols)

    def get_icet_sqs_obj(self, material: Atoms | Structure, cluster_space: ClusterSpace | None = None) -> float:
        if isinstance(material, Structure):
            material = AseAtomsAdaptor.get_atoms(material)

        cluster_space = cluster_space or self._get_cluster_space()
        return compare_cluster_vectors(
            cv_1=cluster_space.get_cluster_vector(material),
            cv_2=self.sqs_vector,
            orbit_data=cluster_space.orbit_data,
            **self._sqs_obj_kwargs,
        )

    def enumerate_sqs_structures(self, cluster_space: ClusterSpace | None = None) -> list:
        """
        Adapted from icet.tools.structure_generation.generate_sqs_by_enumeration.

        Given a ``cluster_space``, generate a special quasirandom structure
        (SQS), i.e., a structure that for a given supercell size provides
        the best possible approximation to a random alloy [ZunWeiFer90]_.

        In the present case, this means that the generated structure will
        have a cluster vector that as closely as possible matches the
        cluster vector of an infintely large randomly occupied supercell.
        Internally the function uses a simulated annealing algorithm and the
        difference between two cluster vectors is calculated with the
        measure suggested by A. van de Walle et al. in Calphad **42**, 13-18
        (2013) [WalTiwJon13]_ (for more information, see
        :class:`mchammer.calculators.TargetVectorCalculator`).

        This functions generates SQS cells by exhaustive enumeration, which
        means that the generated SQS cell is guaranteed to be optimal with
        regard to the specified measure and cell size.

        Parameters
        ----------
        cluster_space : ClusterSpace | None = None
            a cluster space defining the lattice to be occupied
        scaling
            maximum supercell size
        target_concentrations
            concentration of each species in the target structure, per
            sublattice (for example ``{'Au': 0.5, 'Pd': 0.5}`` for a
            single sublattice Au-Pd structure, or
            ``{'A': {'Au': 0.5, 'Pd': 0.5}, 'B': {'H': 0.25, 'X': 0.75}}``
            for a system with two sublattices.
            The symbols defining sublattices ('A', 'B' etc) can be
            found by printing the `cluster_space`
        include_smaller_cells
            if True, search among all supercell sizes including
            ``max_size``, else search only among those exactly matching
            ``max_size``
        pbc
            Periodic boundary conditions for each direction, e.g.,
            ``(True, True, False)``. The axes are defined by
            the cell of ``cluster_space.primitive_structure``.
            Default is periodic boundary in all directions.
        optimality_weight
            controls weighting :math:`L` of perfect correlations, see
            :class:`mchammer.calculators.TargetVectorCalculator`
        tol
            Numerical tolerance
        """

        # Translate concentrations to the format required for concentration
        # restricted enumeration
        cr: dict[str, tuple] = {}
        cluster_space = cluster_space or self._get_cluster_space()
        sublattices = cluster_space.get_sublattices(cluster_space.primitive_structure)
        for sl in sublattices:
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
            raise ValueError(f"This should never happen, but {abs(c_sum - 1)}")

        sizes = list(range(1, self.scaling + 1)) if self.sqs_kwargs["include_smaller_cells"] else [self.scaling]

        # Prepare primitive structure with the right boundary conditions
        prim = cluster_space.primitive_structure
        prim.set_pbc(self.sqs_kwargs["pbc"])

        structures = enumerate_structures(prim, sizes, cluster_space.chemical_symbols, concentration_restrictions=cr)
        chunks = [[] for _ in range(self.instances)]
        iproc = 0
        for structure in structures:
            chunks[iproc].append(structure)
            iproc = (iproc + 1) % self.instances

        manager = multiproc.Manager()
        working_list = manager.list()
        processes = []
        for iproc in range(self.instances):
            process = multiproc.Process(
                target=self._get_best_sqs_from_list,
                args=(
                    chunks[iproc],
                    working_list,
                ),
            )
            processes.append(process)
            process.start()

        for process in processes:
            process.join()

        return list(working_list)

    def _get_best_sqs_from_list(self, structures: list[Atoms], output_list: list[dict]) -> None:
        best_sqs = {"structure": None, "objective_function": 1.0e20}
        cluster_space = self._get_cluster_space()
        for structure in structures:
            objective = self.get_icet_sqs_obj(structure, cluster_space=cluster_space)
            if objective < best_sqs["objective_function"]:
                best_sqs = {"structure": structure, "objective_function": objective}
        output_list.append(best_sqs)

    def monte_carlo_sqs_structures(self) -> list:
        with multiproc.Pool(self.instances) as pool:
            return pool.starmap(self._single_monte_carlo_sqs_run, [() for _ in range(self.instances)])

    def _single_monte_carlo_sqs_run(self):
        cluster_space = self._get_cluster_space()
        sqs_structure = generate_sqs(
            cluster_space=cluster_space,
            scaling=self.scaling,
            target_concentrations=self.target_concentrations,
            **self.sqs_kwargs,
        )
        return {
            "structure": sqs_structure,
            "objective_function": self.get_icet_sqs_obj(sqs_structure, cluster_space=cluster_space),
        }
