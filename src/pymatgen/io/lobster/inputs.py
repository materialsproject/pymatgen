"""Module for reading LOBSTER input files.
For more information on LOBSTER see www.cohp.de.

If you use this module, please cite:
J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier,
"Automated Bonding Analysis with Crystal Orbital Hamilton Populations",
ChemPlusChem 2022, e202200123,
DOI: 10.1002/cplu.202200123.
"""

from __future__ import annotations

import itertools
import os
import re
import warnings
from collections import UserDict
from typing import TYPE_CHECKING

import numpy as np
import spglib
from monty.io import zopen
from monty.json import MSONable
from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Vasprun
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing import Any, ClassVar, Literal

    from typing_extensions import Self

    from pymatgen.core.composition import Composition
    from pymatgen.util.typing import PathLike, Tuple3Ints

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

__author__ = "Janine George, Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__maintainer__ = "Janine George"
__email__ = "janinegeorge.ulfen@gmail.com"
__date__ = "Dec 13, 2017"


due.cite(
    Doi("10.1002/cplu.202200123"),
    description="Automated Bonding Analysis with Crystal Orbital Hamilton Populations",
)


class Lobsterin(UserDict, MSONable):
    """Handle and generate lobsterin files.
    Furthermore, it can also modify INCAR files for LOBSTER, generate KPOINTS files for fatband calculations in LOBSTER,
    and generate the standard primitive cells in a POSCAR file that are needed for the fatband calculations.
    There are also several standard lobsterin files that can be easily generated.

    Reminder: lobsterin keywords are not case sensitive.
    """

    # These keywords need an additional float suffix
    _FLOAT_KEYWORDS: tuple[str, ...] = (
        "COHPstartEnergy",
        "COHPendEnergy",
        "gaussianSmearingWidth",
        "useDecimalPlaces",
        "COHPSteps",
        "basisRotation",
    )

    # These keywords need an additional string suffix
    _STRING_KEYWORDS: tuple[str, ...] = (
        "basisSet",
        "cohpGenerator",
        "realspaceHamiltonian",
        "realspaceOverlap",
        "printPAWRealSpaceWavefunction",
        "printLCAORealSpaceWavefunction",
        "kSpaceCOHP",
        "EwaldSum",
    )

    # The keywords themselves (without suffix) can trigger additional functionalities
    _BOOLEAN_KEYWORDS: tuple[str, ...] = (
        "saveProjectionToFile",
        "skipdos",
        "skipcohp",
        "skipcoop",
        "skipcobi",
        "skipMadelungEnergy",
        "loadProjectionFromFile",
        "printTotalSpilling",
        "forceEnergyRange",
        "DensityOfEnergy",
        "BWDF",
        "BWDFCOHP",
        "skipPopulationAnalysis",
        "skipGrossPopulation",
        "userecommendedbasisfunctions",
        "skipProjection",
        "printLmosOnAtoms",
        "writeBasisFunctions",
        "writeMatricesToFile",
        "noFFTforVisualization",
        "RMSp",
        "onlyReadVasprun.xml",
        "noMemoryMappedFiles",
        "skipPAWOrthonormalityTest",
        "doNotIgnoreExcessiveBands",
        "doNotUseAbsoluteSpilling",
        "skipReOrthonormalization",
        "forceV1HMatrix",
        "useOriginalTetrahedronMethod",
        "forceEnergyRange",
        "bandwiseSpilling",
        "kpointwiseSpilling",
        "LSODOS",
        "autoRotate",
        "doNotOrthogonalizeBasis",
    )

    # These keywords need additional string suffixes.
    # They could be used multiple times within one lobsterin.
    _LIST_KEYWORDS: tuple[str, ...] = (
        "basisfunctions",
        "cohpbetween",
        "createFatband",
        "customSTOforAtom",
        "cobiBetween",
    )

    # Generate {lowered: original} mappings
    FLOAT_KEYWORDS: ClassVar[dict[str, str]] = {key.lower(): key for key in _FLOAT_KEYWORDS}
    STRING_KEYWORDS: ClassVar[dict[str, str]] = {key.lower(): key for key in _STRING_KEYWORDS}
    BOOLEAN_KEYWORDS: ClassVar[dict[str, str]] = {key.lower(): key for key in _BOOLEAN_KEYWORDS}
    LIST_KEYWORDS: ClassVar[dict[str, str]] = {key.lower(): key for key in _LIST_KEYWORDS}

    # All known keywords
    AVAILABLE_KEYWORDS: ClassVar[dict[str, str]] = {
        **FLOAT_KEYWORDS,
        **STRING_KEYWORDS,
        **BOOLEAN_KEYWORDS,
        **LIST_KEYWORDS,
    }

    def __init__(self, settingsdict: dict) -> None:
        """
        Args:
            settingsdict: dict to initialize Lobsterin.
        """
        super().__init__()

        # Check for duplicates from case sensitivity
        keys = tuple(map(str.lower, settingsdict.keys()))
        if len(keys) != len(set(keys)):
            raise KeyError("There are duplicates for the keywords!")

        self.update(settingsdict)

    def __setitem__(self, key: str, val: Any) -> None:
        """
        Necessary due to the missing case sensitivity of lobsterin
        keywords. Also clean the keys and values by stripping white spaces.

        Raises:
            KeyError: if keyword is not available.
        """
        key = key.strip().lower()

        if key not in type(self).AVAILABLE_KEYWORDS:
            raise KeyError(f"Key {key} is currently not available")

        super().__setitem__(key, val.strip() if isinstance(val, str) else val)

    def __getitem__(self, key: str) -> Any:
        """To avoid cases sensitivity problems."""
        try:
            return super().__getitem__(key.strip().lower())

        except KeyError as exc:
            raise KeyError(f"{key=} is not available") from exc

    def __contains__(self, key: str) -> bool:
        """To avoid cases sensitivity problems."""
        return super().__contains__(key.lower().strip())

    def __delitem__(self, key: str) -> None:
        """To avoid cases sensitivity problems."""
        super().__delitem__(key.lower().strip())

    def diff(self, other: Self) -> dict[str, dict[str, Any]]:
        """Compare two Lobsterin and find which parameters are the same.
        Similar to the diff method of Incar.

        Args:
            other (Lobsterin): Lobsterin object to compare to.

        Returns:
            dict: {"Same": same_params, "Different": diff_params}
        """
        same_params = {}
        diff_params = {}

        # Check self
        for k1, v1 in self.items():
            if k1 not in other:
                diff_params[k1] = {"lobsterin1": v1, "lobsterin2": None}

            # String keywords
            elif isinstance(v1, str):
                if v1 != other[k1]:
                    diff_params[k1] = {"lobsterin1": v1, "lobsterin2": other[k1]}
                else:
                    same_params[k1] = v1

            # List keywords
            elif isinstance(v1, list):
                new_set1 = {value.strip().lower() for value in v1}
                new_set2 = {value.strip().lower() for value in other[k1]}
                if new_set1 != new_set2:
                    diff_params[k1] = {"lobsterin1": v1, "lobsterin2": other[k1]}

            # Float/boolean keywords
            elif v1 != other[k1]:
                diff_params[k1] = {"lobsterin1": v1, "lobsterin2": other[k1]}
            else:
                same_params[k1] = v1

        # Check other
        for k2, v2 in other.items():
            if k2 not in self and k2 not in diff_params:
                diff_params[k2] = {"lobsterin1": None, "lobsterin2": v2}

        return {"Same": same_params, "Different": diff_params}

    def write_lobsterin(
        self,
        path: PathLike = "lobsterin",
        overwritedict: dict | None = None,
    ) -> None:
        """Write a lobsterin file, and recover keys to Camel case.

        Args:
            path (str): filename of the output lobsterin file
            overwritedict (dict): dict that can be used to update lobsterin, e.g. {"skipdos": True}
        """
        # Update previous entries
        if overwritedict is not None:
            self.update(overwritedict)

        with open(path, mode="w", encoding="utf-8") as file:
            for key in self:
                if key in type(self).FLOAT_KEYWORDS or key in type(self).STRING_KEYWORDS:
                    file.write(f"{type(self).AVAILABLE_KEYWORDS[key]} {self.get(key)}\n")

                elif key in type(self).BOOLEAN_KEYWORDS:
                    file.write(f"{type(self).BOOLEAN_KEYWORDS[key]}\n")

                elif key in type(self).LIST_KEYWORDS:
                    for value in self.get(key):  # type: ignore[union-attr]
                        file.write(f"{type(self).LIST_KEYWORDS[key]} {value}\n")

    def as_dict(self) -> dict:
        """MSONable dict."""
        dct = dict(self)
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            Lobsterin
        """
        return cls({key: val for key, val in dct.items() if key not in {"@module", "@class"}})

    def _get_nbands(self, structure: Structure) -> int:
        """Get number of bands."""
        if self.get("basisfunctions") is None:
            raise ValueError("No basis functions are provided. The program cannot calculate nbands.")

        basis_functions: list[str] = []
        for string_basis in self["basisfunctions"]:
            string_basis_raw = string_basis.strip().split(" ")
            while "" in string_basis_raw:
                string_basis_raw.remove("")
            for _idx in range(int(structure.composition.element_composition[string_basis_raw[0]])):
                basis_functions.extend(string_basis_raw[1:])

        num_basis_functions = 0
        for basis in basis_functions:
            if "s" in basis:
                num_basis_functions += 1
            elif "p" in basis:
                num_basis_functions += 3
            elif "d" in basis:
                num_basis_functions += 5
            elif "f" in basis:
                num_basis_functions += 7

        return int(num_basis_functions)

    def write_INCAR(
        self,
        incar_input: PathLike = "INCAR",
        incar_output: PathLike = "INCAR.lobster",
        poscar_input: PathLike = "POSCAR",
        isym: Literal[-1, 0] = 0,
        further_settings: dict | None = None,
    ) -> None:
        """Write INCAR file. Will only make the run static, insert NBANDS,
        set ISYM=0, LWAVE=True and you have to check for the rest.

        Args:
            incar_input (PathLike): path to input INCAR
            incar_output (PathLike): path to output INCAR
            poscar_input (PathLike): path to input POSCAR
            isym (-1 | 0): ISYM value.
            further_settings (dict): A dict can be used to include further settings, e.g. {"ISMEAR":-5}
        """
        # Read INCAR from file, which will be modified
        incar = Incar.from_file(incar_input)
        warnings.warn("Please check your incar_input before using it. This method only changes three settings!")
        if isym in {-1, 0}:
            incar["ISYM"] = isym
        else:
            raise ValueError(f"Got {isym=}, must be -1 or 0")

        incar["NSW"] = 0
        incar["LWAVE"] = True

        # Get NBANDS from _get_nbands (use basis set that is inserted)
        incar["NBANDS"] = self._get_nbands(Structure.from_file(poscar_input))
        if further_settings is not None:
            for key, item in further_settings.items():
                incar[key] = item
        incar.write_file(incar_output)

    @staticmethod
    def get_basis(
        structure: Structure,
        potcar_symbols: list[str],
        address_basis_file: PathLike | None = None,
    ) -> list[str]:
        """Get the basis functions from given potcar_symbols, e.g., ["Fe_pv", "Si"].

        Args:
            structure (Structure): Structure object
            potcar_symbols: list of potcar symbols
            address_basis_file (PathLike): path to the basis file

        Returns:
            basis
        """
        if address_basis_file is None:
            address_basis_file = f"{MODULE_DIR}/lobster_basis/BASIS_PBE_54_standard.yaml"

        atom_types_potcar = [name.split("_")[0] for name in potcar_symbols]

        if set(structure.symbol_set) != set(atom_types_potcar):
            raise ValueError("Your POSCAR does not correspond to your POTCAR!")

        basis = loadfn(address_basis_file)["BASIS"]

        basis_functions = []
        list_forin = []
        for idx, name in enumerate(potcar_symbols):
            if name not in basis:
                raise ValueError(
                    f"Missing basis information for POTCAR symbol: {name}. Please provide the basis manually."
                )
            basis_functions.append(basis[name].split())
            list_forin.append(f"{atom_types_potcar[idx]} {basis[name]}")
        return list_forin

    @staticmethod
    def get_all_possible_basis_functions(
        structure: Structure,
        potcar_symbols: list[str],
        address_basis_file_min: PathLike | None = None,
        address_basis_file_max: PathLike | None = None,
    ) -> list[dict]:
        """
        Args:
            structure: Structure object
            potcar_symbols: list of the potcar symbols
            address_basis_file_min: path to file with the minimum required basis by the POTCAR
            address_basis_file_max: path to file with the largest possible basis of the POTCAR.

        Returns:
            list[dict]: Can be used to create new Lobsterin objects in
                standard_calculations_from_vasp_files as dict_for_basis
        """
        max_basis = Lobsterin.get_basis(
            structure=structure,
            potcar_symbols=potcar_symbols,
            address_basis_file=address_basis_file_max or f"{MODULE_DIR}/lobster_basis/BASIS_PBE_54_max.yaml",
        )
        min_basis = Lobsterin.get_basis(
            structure=structure,
            potcar_symbols=potcar_symbols,
            address_basis_file=address_basis_file_min or f"{MODULE_DIR}/lobster_basis/BASIS_PBE_54_min.yaml",
        )
        all_basis = get_all_possible_basis_combinations(min_basis=min_basis, max_basis=max_basis)
        list_basis_dict = []
        for basis in all_basis:
            basis_dict = {}

            for elba in basis:
                basplit = elba.split()
                basis_dict[basplit[0]] = " ".join(basplit[1:])
            list_basis_dict.append(basis_dict)
        return list_basis_dict

    @staticmethod
    def write_POSCAR_with_standard_primitive(
        POSCAR_input: PathLike = "POSCAR",
        POSCAR_output: PathLike = "POSCAR.lobster",
        symprec: float = 0.01,
    ) -> None:
        """Write a POSCAR with the standard primitive cell.
        This is needed to arrive at the correct kpath.

        Args:
            POSCAR_input (PathLike): Input POSCAR file
            POSCAR_output (PathLike): Output POSCAR file
            symprec (float): precision to find symmetry
        """
        structure = Structure.from_file(POSCAR_input)
        kpath = HighSymmKpath(structure, symprec=symprec)
        new_structure = kpath.prim
        new_structure.to(fmt="POSCAR", filename=POSCAR_output)

    @staticmethod
    def write_KPOINTS(
        POSCAR_input: PathLike = "POSCAR",
        KPOINTS_output: PathLike = "KPOINTS.lobster",
        reciprocal_density: int = 100,
        isym: Literal[-1, 0] = 0,
        from_grid: bool = False,
        input_grid: Tuple3Ints = (5, 5, 5),
        line_mode: bool = True,
        kpoints_line_density: int = 20,
        symprec: float = 0.01,
    ) -> None:
        """Write a gamma-centered KPOINTS file for LOBSTER.

        Args:
            POSCAR_input (PathLike): path to POSCAR
            KPOINTS_output (PathLike): path to output KPOINTS
            reciprocal_density (int): Grid density
            isym (-1 | 0): ISYM value.
            from_grid (bool): If True KPOINTS will be generated with the help of a grid given in input_grid.
                Otherwise, they will be generated from the reciprocal_density
            input_grid (tuple): grid to generate the KPOINTS file
            line_mode (bool): If True, band structure will be generated
            kpoints_line_density (int): density of the lines in the band structure
            symprec (float): precision to determine symmetry
        """
        structure = Structure.from_file(POSCAR_input)
        if not from_grid:
            kpoint_grid = Kpoints.automatic_density_by_vol(structure, reciprocal_density).kpts
            mesh = kpoint_grid[0]
        else:
            mesh = input_grid

        # The following code is taken from SpacegroupAnalyzer
        # We need to switch off symmetry here
        matrix = structure.lattice.matrix
        positions = structure.frac_coords
        unique_species: list[Composition] = []
        zs = []
        magmoms = []

        for species, group in itertools.groupby(structure, key=lambda s: s.species):
            if species in unique_species:
                ind = unique_species.index(species)
                zs.extend([ind + 1] * len(tuple(group)))
            else:
                unique_species.append(species)
                zs.extend([len(unique_species)] * len(tuple(group)))

        for site in structure:
            if hasattr(site, "magmom"):
                magmoms.append(site.magmom)
            elif site.is_ordered and hasattr(site.specie, "spin"):
                magmoms.append(site.specie.spin)
            else:
                magmoms.append(0)

        # For now, we are setting MAGMOM to zero. (Taken from INCAR class)
        cell = matrix, positions, zs, magmoms
        # TODO: what about this shift?
        mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])

        # Get the KPOINTS for the grid
        if isym == -1:
            kpts = []
            weights = []
            all_labels = []
            for gp in grid:
                kpts.append(gp.astype(float) / mesh)
                weights.append(float(1))
                all_labels.append("")

        # Time reversal symmetry: k and -k are equivalent
        elif isym == 0:
            kpts = []
            weights = []
            all_labels = []
            newlist = [list(gp) for gp in list(grid)]
            mapping = []
            for gp in newlist:
                minus_gp = [-k for k in gp]
                if minus_gp in newlist and minus_gp != [0, 0, 0]:
                    mapping.append(newlist.index(minus_gp))
                else:
                    mapping.append(newlist.index(gp))

            for igp, gp in enumerate(newlist):
                if mapping[igp] > igp:
                    kpts.append(np.array(gp).astype(float) / mesh)
                    weights.append(float(2))
                    all_labels.append("")
                elif mapping[igp] == igp:
                    kpts.append(np.array(gp).astype(float) / mesh)
                    weights.append(float(1))
                    all_labels.append("")

        else:
            raise ValueError(f"Got {isym=}, must be -1 or 0")

        # Line mode
        if line_mode:
            kpath = HighSymmKpath(structure, symprec=symprec)
            if not np.allclose(kpath.prim.lattice.matrix, structure.lattice.matrix):
                raise ValueError(
                    "You are not using the standard primitive cell. The k-path is not correct. Please generate a "
                    "standard primitive cell first."
                )

            frac_k_points, labels = kpath.get_kpoints(line_density=kpoints_line_density, coords_are_cartesian=False)

            for k, f in enumerate(frac_k_points):
                kpts.append(f)
                weights.append(0.0)
                all_labels.append(labels[k])
        comment = f"{isym=}, grid: {mesh} plus kpoint path" if line_mode else f"{isym=}, grid: {mesh}"

        kpoints_instance = Kpoints(
            comment=comment,
            style=Kpoints.supported_modes.Reciprocal,
            num_kpts=len(kpts),
            kpts=tuple(kpts),
            kpts_weights=weights,
            labels=all_labels,
        )

        kpoints_instance.write_file(filename=KPOINTS_output)

    @classmethod
    def from_file(cls, lobsterin: PathLike) -> Self:
        """Create Lobsterin from lobsterin file.

        Args:
            lobsterin (PathLike): path to lobsterin.

        Returns:
            Lobsterin object
        """
        with zopen(lobsterin, mode="rt") as file:
            lines = file.read().split("\n")
        if not lines:
            raise RuntimeError("lobsterin file contains no data.")

        lobsterin_dict: dict[str, Any] = {}
        for line in lines:
            # Remove comment lines and in-line comments
            if line := re.split(r"[!#//]", line)[0].strip():
                # Extract keywords
                line_parts = line.replace("\t", " ").strip().split()
                if line_parts:
                    key = line_parts[0].lower()
                else:
                    continue

                # Avoid duplicates for float/string keywords
                if (key in cls.FLOAT_KEYWORDS or key in cls.STRING_KEYWORDS) and key in lobsterin_dict:
                    raise ValueError(f"Same keyword {key} twice!")

                # Parse by keyword type
                if key in cls.BOOLEAN_KEYWORDS:
                    lobsterin_dict[key] = True

                elif key in cls.FLOAT_KEYWORDS:
                    lobsterin_dict[key] = float(line_parts[1])

                elif key in cls.STRING_KEYWORDS:
                    lobsterin_dict[key] = " ".join(line_parts[1:])

                elif key in cls.LIST_KEYWORDS:
                    if key in lobsterin_dict:
                        lobsterin_dict[key].append(" ".join(line_parts[1:]))
                    else:
                        lobsterin_dict[key] = [" ".join(line_parts[1:])]

                else:
                    raise ValueError(f"Invalid {key=}.")

        return cls(lobsterin_dict)

    @staticmethod
    def _get_potcar_symbols(POTCAR_input: PathLike) -> list[str]:
        """Get the name of the species in the POTCAR.

        Args:
            POTCAR_input (PathLike): path to potcar file

        Returns:
            list[str]: names of the species
        """
        potcar = Potcar.from_file(POTCAR_input)
        for pot in potcar:
            if pot.potential_type != "PAW":
                raise ValueError("Lobster only works with PAW! Use different POTCARs")

        # Warning about a bug in LOBSTER-4.1.0
        with zopen(POTCAR_input, mode="r") as file:
            data = file.read()

        if isinstance(data, bytes):
            data = data.decode("utf-8")

        if "SHA256" in data or "COPYR" in data:
            warnings.warn(
                "These POTCARs are not compatible with "
                "Lobster up to version 4.1.0."
                "\n The keywords SHA256 and COPYR "
                "cannot be handled by Lobster"
                " \n and will lead to wrong results."
            )

        if potcar.functional != "PBE":
            raise RuntimeError("We only have BASIS options for PBE so far")

        return [name["symbol"] for name in potcar.spec]

    @classmethod
    def standard_calculations_from_vasp_files(
        cls,
        POSCAR_input: PathLike = "POSCAR",
        INCAR_input: PathLike = "INCAR",
        POTCAR_input: PathLike | None = None,
        Vasprun_output: PathLike = "vasprun.xml",
        dict_for_basis: dict | None = None,
        option: str = "standard",
    ) -> Self:
        """Generate lobsterin with standard settings.

        Args:
            POSCAR_input (PathLike): path to POSCAR
            INCAR_input (PathLike): path to INCAR
            POTCAR_input (PathLike): path to POTCAR
            Vasprun_output (PathLike): path to vasprun.xml
            dict_for_basis (dict): can be provided: it should look the following:
                dict_for_basis={"Fe":'3p 3d 4s 4f', "C": '2s 2p'} and will overwrite all settings from POTCAR_input
            option (str): 'standard' will start a normal LOBSTER run where COHPs, COOPs, DOS, CHARGE etc. will be
                calculated
                'standard_with_energy_range_from_vasprun' will start a normal LOBSTER run for entire energy range
                of VASP static run. vasprun.xml file needs to be in current directory.
                'standard_from_projection' will start a normal LOBSTER run from a projection
                'standard_with_fatband' will do a fatband calculation, run over all orbitals
                'onlyprojection' will only do a projection
                'onlydos' will only calculate a projected dos
                'onlycohp' will only calculate cohp
                'onlycoop' will only calculate coop
                'onlycohpcoop' will only calculate cohp and coop

        Returns:
            Lobsterin with standard settings
        """
        warnings.warn(
            "Always check and test the provided basis functions. The spilling of your Lobster calculation might help"
        )

        if option not in {
            "standard",
            "standard_from_projection",
            "standard_with_fatband",
            "standard_with_energy_range_from_vasprun",
            "onlyprojection",
            "onlydos",
            "onlycohp",
            "onlycoop",
            "onlycobi",
            "onlycohpcoop",
            "onlycohpcoopcobi",
            "onlymadelung",
        }:
            raise ValueError("The option is not valid!")

        lobsterin_dict: dict[str, Any] = {
            # This basis set covers elements up to Lr (Z = 103)
            "basisSet": "pbeVaspFit2015",
            # Energies around e-fermi
            "COHPstartEnergy": -35.0,
            "COHPendEnergy": 5.0,
        }

        if option in {
            "standard",
            "standard_with_energy_range_from_vasprun",
            "onlycohp",
            "onlycoop",
            "onlycobi",
            "onlycohpcoop",
            "onlycohpcoopcobi",
            "standard_with_fatband",
        }:
            # Every interaction with a distance of 6.0 Å is checked
            lobsterin_dict["cohpGenerator"] = "from 0.1 to 6.0 orbitalwise"
            # Save the projection
            lobsterin_dict["saveProjectionToFile"] = True

        if option == "standard_from_projection":
            lobsterin_dict["cohpGenerator"] = "from 0.1 to 6.0 orbitalwise"
            lobsterin_dict["loadProjectionFromFile"] = True

        elif option == "standard_with_energy_range_from_vasprun":
            vasp_run = Vasprun(Vasprun_output)
            lobsterin_dict["COHPstartEnergy"] = round(
                min(vasp_run.complete_dos.energies - vasp_run.complete_dos.efermi), 4
            )
            lobsterin_dict["COHPendEnergy"] = round(
                max(vasp_run.complete_dos.energies - vasp_run.complete_dos.efermi), 4
            )
            lobsterin_dict["COHPSteps"] = len(vasp_run.complete_dos.energies)

        # TODO: add COBI here! might be relevant LOBSTER version
        elif option == "onlycohp":
            lobsterin_dict["skipdos"] = True
            lobsterin_dict["skipcoop"] = True
            lobsterin_dict["skipPopulationAnalysis"] = True
            lobsterin_dict["skipGrossPopulation"] = True
            # LOBSTER-4.1.0
            lobsterin_dict["skipcobi"] = True
            lobsterin_dict["skipMadelungEnergy"] = True

        elif option == "onlycoop":
            lobsterin_dict["skipdos"] = True
            lobsterin_dict["skipcohp"] = True
            lobsterin_dict["skipPopulationAnalysis"] = True
            lobsterin_dict["skipGrossPopulation"] = True
            # LOBSTER-4.1.0
            lobsterin_dict["skipcobi"] = True
            lobsterin_dict["skipMadelungEnergy"] = True

        elif option == "onlycohpcoop":
            lobsterin_dict["skipdos"] = True
            lobsterin_dict["skipPopulationAnalysis"] = True
            lobsterin_dict["skipGrossPopulation"] = True
            # LOBSTER-4.1.0
            lobsterin_dict["skipcobi"] = True
            lobsterin_dict["skipMadelungEnergy"] = True

        elif option == "onlycohpcoopcobi":
            lobsterin_dict["skipdos"] = True
            lobsterin_dict["skipPopulationAnalysis"] = True
            lobsterin_dict["skipGrossPopulation"] = True
            lobsterin_dict["skipMadelungEnergy"] = True

        elif option == "onlydos":
            lobsterin_dict["skipcohp"] = True
            lobsterin_dict["skipcoop"] = True
            lobsterin_dict["skipPopulationAnalysis"] = True
            lobsterin_dict["skipGrossPopulation"] = True
            # LOBSTER-4.1.0
            lobsterin_dict["skipcobi"] = True
            lobsterin_dict["skipMadelungEnergy"] = True

        elif option == "onlyprojection":
            lobsterin_dict["skipdos"] = True
            lobsterin_dict["skipcohp"] = True
            lobsterin_dict["skipcoop"] = True
            lobsterin_dict["skipPopulationAnalysis"] = True
            lobsterin_dict["skipGrossPopulation"] = True
            lobsterin_dict["saveProjectionToFile"] = True
            # LOBSTER-4.1.0
            lobsterin_dict["skipcobi"] = True
            lobsterin_dict["skipMadelungEnergy"] = True

        elif option == "onlycobi":
            lobsterin_dict["skipdos"] = True
            lobsterin_dict["skipcohp"] = True
            lobsterin_dict["skipPopulationAnalysis"] = True
            lobsterin_dict["skipGrossPopulation"] = True
            # LOBSTER-4.1.0
            lobsterin_dict["skipcobi"] = True
            lobsterin_dict["skipMadelungEnergy"] = True

        elif option == "onlymadelung":
            lobsterin_dict["skipdos"] = True
            lobsterin_dict["skipcohp"] = True
            lobsterin_dict["skipcoop"] = True
            lobsterin_dict["skipPopulationAnalysis"] = True
            lobsterin_dict["skipGrossPopulation"] = True
            lobsterin_dict["saveProjectionToFile"] = True
            # LOBSTER-4.1.0
            lobsterin_dict["skipcobi"] = True

        incar = Incar.from_file(INCAR_input)

        if incar["ISMEAR"] == 0:
            lobsterin_dict["gaussianSmearingWidth"] = incar["SIGMA"]

        if incar["ISMEAR"] != 0 and option == "standard_with_fatband":
            raise ValueError("ISMEAR has to be 0 for a fatband calculation with Lobster")

        if dict_for_basis is not None:
            # dict_for_basis = {"Fe":"3p 3d 4s 4f", "C": "2s 2p"}
            # Will just insert this basis and not check with poscar
            basis = [f"{key} {value}" for key, value in dict_for_basis.items()]
        elif POTCAR_input is not None:
            # Get basis functions from POTCAR
            potcar_names = cls._get_potcar_symbols(POTCAR_input=POTCAR_input)

            basis = cls.get_basis(structure=Structure.from_file(POSCAR_input), potcar_symbols=potcar_names)
        else:
            raise ValueError("basis cannot be generated")

        lobsterin_dict["basisfunctions"] = basis

        if option == "standard_with_fatband":
            lobsterin_dict["createFatband"] = basis

        return cls(lobsterin_dict)


def get_all_possible_basis_combinations(min_basis: list, max_basis: list) -> list[list[str]]:
    """Get all possible basis combinations.

    Args:
        min_basis: list of basis entries: e.g., ["Si 3p 3s"]
        max_basis: list of basis entries: e.g., ["Si 3p 3s"].

    Returns:
        list[list[str]]: all possible combinations of basis functions, e.g. [["Si 3p 3s"]]
    """
    max_basis_lists = [x.split() for x in max_basis]
    min_basis_lists = [x.split() for x in min_basis]

    # Get all possible basis functions
    basis_dict: dict[str, dict] = {}
    for iel, el in enumerate(max_basis_lists):
        basis_dict[el[0]] = {"fixed": [], "variable": [], "combinations": []}
        for basis in el[1:]:
            if basis in min_basis_lists[iel]:
                basis_dict[el[0]]["fixed"].append(basis)
            if basis not in min_basis_lists[iel]:
                basis_dict[el[0]]["variable"].append(basis)
        for L in range(len(basis_dict[el[0]]["variable"]) + 1):
            for subset in itertools.combinations(basis_dict[el[0]]["variable"], L):
                basis_dict[el[0]]["combinations"].append(" ".join([el[0]] + basis_dict[el[0]]["fixed"] + list(subset)))

    list_basis = [item["combinations"] for item in basis_dict.values()]

    # Get all combinations
    start_basis = list_basis[0]
    if len(list_basis) > 1:
        for el in list_basis[1:]:
            new_start_basis = []
            for elbasis in start_basis:
                for elbasis2 in el:
                    if not isinstance(elbasis, list):
                        new_start_basis.append([elbasis, elbasis2])
                    else:
                        new_start_basis.append([*elbasis.copy(), elbasis2])
            start_basis = new_start_basis
        return start_basis
    return [[basis] for basis in start_basis]
