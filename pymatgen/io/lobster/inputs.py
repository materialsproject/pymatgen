# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License

"""
Module for reading Lobster input files. For more information
on LOBSTER see www.cohp.de.
"""

from __future__ import annotations

import itertools
import os
import warnings
from typing import Any, Sequence

import numpy as np
import spglib
from monty.io import zopen
from monty.json import MSONable
from monty.serialization import loadfn

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Vasprun
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.symmetry.bandstructure import HighSymmKpath

__author__ = "Janine George, Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__maintainer__ = "Janine George, Marco Esters"
__email__ = "janine.george@uclouvain.be, esters@uoregon.edu"
__date__ = "Dec 13, 2017"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class Lobsterin(dict, MSONable):
    """
    This class can handle and generate lobsterin files
    Furthermore, it can also modify INCAR files for lobster, generate KPOINT files for fatband calculations in Lobster,
    and generate the standard primitive cells in a POSCAR file that are needed for the fatband calculations.
    There are also several standard lobsterin files that can be easily generated.
    """

    # reminder: lobster is not case sensitive

    # keyword + one float can be used in file
    FLOATKEYWORDS = [
        "COHPstartEnergy",
        "COHPendEnergy",
        "gaussianSmearingWidth",
        "useDecimalPlaces",
        "COHPSteps",
    ]
    # one of these keywords +endstring can be used in file
    STRINGKEYWORDS = [
        "basisSet",
        "cohpGenerator",
        "realspaceHamiltonian",
        "realspaceOverlap",
        "printPAWRealSpaceWavefunction",
        "printLCAORealSpaceWavefunction",
        "kSpaceCOHP",
        "EwaldSum",
    ]
    # the keyword alone will turn on or off a function
    BOOLEANKEYWORDS = [
        "saveProjectionToFile",
        "skipdos",
        "skipcohp",
        "skipcoop",
        "skipcobi",
        "skipMadelungEnergy",
        "loadProjectionFromFile",
        "forceEnergyRange",
        "DensityOfEnergy",
        "BWDF",
        "BWDFCOHP",
        "skipPopulationAnalysis",
        "skipGrossPopulation",
        "userecommendedbasisfunctions",
        "skipProjection",
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
    ]
    # several of these keywords + ending can be used in a lobsterin file:
    LISTKEYWORDS = ["basisfunctions", "cohpbetween", "createFatband"]

    # all keywords known to this class so far
    AVAILABLEKEYWORDS = FLOATKEYWORDS + STRINGKEYWORDS + BOOLEANKEYWORDS + LISTKEYWORDS

    def __init__(self, settingsdict: dict):
        """
        Args:
            settingsdict: dict to initialize Lobsterin
        """
        super().__init__()
        # check for duplicates
        listkey = [key.lower() for key in settingsdict]
        if len(listkey) != len(list(set(listkey))):
            raise OSError("There are duplicates for the keywords! The program will stop here.")
        self.update(settingsdict)

    def __setitem__(self, key, val):
        """
        Add parameter-val pair to Lobsterin. Warns if parameter is not in list of
        valid lobsterintags. Also cleans the parameter and val by stripping
        leading and trailing white spaces. Similar to INCAR class.
        """
        # due to the missing case sensitivity of lobster, the following code is necessary
        found = False
        for key_here in self:
            if key.strip().lower() == key_here.lower():
                new_key = key_here
                found = True
        if not found:
            new_key = key
        if new_key.lower() not in [element.lower() for element in Lobsterin.AVAILABLEKEYWORDS]:
            raise ValueError("Key is currently not available")

        super().__setitem__(new_key, val.strip() if isinstance(val, str) else val)

    def __getitem__(self, item):
        """
        Implements getitem from dict to avoid problems with cases
        """
        found = False
        for key_here in self:
            if item.strip().lower() == key_here.lower():
                new_key = key_here
                found = True
        if not found:
            new_key = item

        val = dict.__getitem__(self, new_key)
        return val

    def diff(self, other):
        """
        Diff function for lobsterin. Compares two lobsterin and indicates which parameters are the same.
        Similar to the diff in INCAR.

        Args:
            other (Lobsterin): Lobsterin object to compare to

        Returns:
            dict with differences and similarities
        """
        similar_param = {}
        different_param = {}
        key_list_others = [element.lower() for element in other]

        for k1, v1 in self.items():
            k1lower = k1.lower()
            if k1lower not in key_list_others:
                different_param[k1.upper()] = {"lobsterin1": v1, "lobsterin2": None}
            else:
                for key_here in other:
                    if k1.lower() == key_here.lower():
                        new_key = key_here

                if isinstance(v1, str):
                    if v1.strip().lower() != other[new_key].strip().lower():

                        different_param[k1.upper()] = {
                            "lobsterin1": v1,
                            "lobsterin2": other[new_key],
                        }
                    else:
                        similar_param[k1.upper()] = v1
                elif isinstance(v1, list):
                    new_set1 = {element.strip().lower() for element in v1}
                    new_set2 = {element.strip().lower() for element in other[new_key]}
                    if new_set1 != new_set2:
                        different_param[k1.upper()] = {
                            "lobsterin1": v1,
                            "lobsterin2": other[new_key],
                        }
                else:
                    if v1 != other[new_key]:
                        different_param[k1.upper()] = {
                            "lobsterin1": v1,
                            "lobsterin2": other[new_key],
                        }
                    else:
                        similar_param[k1.upper()] = v1

        for k2, v2 in other.items():
            if k2.upper() not in similar_param and k2.upper() not in different_param:
                for key_here in self:
                    if k2.lower() == key_here.lower():
                        new_key = key_here
                    else:
                        new_key = k2
                if new_key not in self:
                    different_param[k2.upper()] = {"lobsterin1": None, "lobsterin2": v2}
        return {"Same": similar_param, "Different": different_param}

    def _get_nbands(self, structure: Structure):
        """
        Get number of bands.
        """
        if self.get("basisfunctions") is None:
            raise OSError("No basis functions are provided. The program cannot calculate nbands.")

        basis_functions: list[str] = []
        for string_basis in self["basisfunctions"]:
            # string_basis.lstrip()
            string_basis_raw = string_basis.strip().split(" ")
            while "" in string_basis_raw:
                string_basis_raw.remove("")
            for _idx in range(0, int(structure.composition.element_composition[string_basis_raw[0]])):
                basis_functions.extend(string_basis_raw[1:])

        no_basis_functions = 0
        for basis in basis_functions:
            if "s" in basis:
                no_basis_functions = no_basis_functions + 1
            elif "p" in basis:
                no_basis_functions = no_basis_functions + 3
            elif "d" in basis:
                no_basis_functions = no_basis_functions + 5
            elif "f" in basis:
                no_basis_functions = no_basis_functions + 7

        return int(no_basis_functions)

    def write_lobsterin(self, path="lobsterin", overwritedict=None):
        """
        Writes a lobsterin file

        Args:
            path (str): filename of the lobsterin file that will be written
            overwritedict (dict): dict that can be used to overwrite lobsterin, e.g. {"skipdos": True}
        """
        # will overwrite previous entries
        # has to search first if entry is already in Lobsterindict (due to case insensitivity)
        if overwritedict is not None:
            for key, entry in overwritedict.items():
                found = False
                for key2 in self:
                    if key.lower() == key2.lower():
                        self[key2] = entry
                        found = True
                if not found:
                    self[key] = entry

        filename = path
        with open(filename, "w") as f:
            for key in Lobsterin.AVAILABLEKEYWORDS:
                if key.lower() in [element.lower() for element in self]:
                    if key.lower() in [element.lower() for element in Lobsterin.FLOATKEYWORDS]:
                        f.write(key + " " + str(self.get(key)) + "\n")
                    elif key.lower() in [element.lower() for element in Lobsterin.BOOLEANKEYWORDS]:
                        # checks if entry is True or False
                        for key_here in self:
                            if key.lower() == key_here.lower():
                                new_key = key_here
                        if self.get(new_key):
                            f.write(key + "\n")
                    elif key.lower() in [element.lower() for element in Lobsterin.STRINGKEYWORDS]:
                        f.write(key + " " + str(self.get(key) + "\n"))
                    elif key.lower() in [element.lower() for element in Lobsterin.LISTKEYWORDS]:
                        for entry in self.get(key):
                            f.write(key + " " + str(entry) + "\n")

    def as_dict(self):
        """
        :return: MSONable dict
        """
        d = dict(self)
        d["@module"] = type(self).__module__
        d["@class"] = type(self).__name__
        return d

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: Lobsterin
        """
        return Lobsterin({k: v for k, v in d.items() if k not in ["@module", "@class"]})

    def write_INCAR(
        self,
        incar_input: str = "INCAR",
        incar_output: str = "INCAR.lobster",
        poscar_input: str = "POSCAR",
        isym: int = -1,
        further_settings: dict | None = None,
    ):
        """
        Will only make the run static, insert nbands, make ISYM=-1, set LWAVE=True and write a new INCAR.
        You have to check for the rest.

        Args:
            incar_input (str): path to input INCAR
            incar_output (str): path to output INCAR
            poscar_input (str): path to input POSCAR
            isym (int): isym equal to -1 or 0 are possible. Current Lobster version only allow -1.
            further_settings (dict): A dict can be used to include further settings, e.g. {"ISMEAR":-5}
        """
        # reads old incar from file, this one will be modified
        incar = Incar.from_file(incar_input)
        warnings.warn("Please check your incar_input before using it. This method only changes three settings!")
        if isym == -1:
            incar["ISYM"] = -1
        elif isym == 0:
            incar["ISYM"] = 0
        else:
            ValueError("isym has to be -1 or 0.")
        incar["NSW"] = 0
        incar["LWAVE"] = True
        # get nbands from _get_nbands (use basis set that is inserted)
        incar["NBANDS"] = self._get_nbands(Structure.from_file(poscar_input))
        if further_settings is not None:
            for key, item in further_settings.items():
                incar[key] = item
        # print it to file
        incar.write_file(incar_output)

    @staticmethod
    def get_basis(
        structure: Structure,
        potcar_symbols: list,
        address_basis_file: str | None = None,
    ):
        """
        Will get the basis from given potcar_symbols (e.g., ["Fe_pv","Si"]
        #include this in lobsterin class

        Args:
            structure (Structure): Structure object
            potcar_symbols: list of potcar symbols

        Returns:
            returns basis
        """
        if address_basis_file is None:
            address_basis_file = os.path.join(MODULE_DIR, "lobster_basis/BASIS_PBE_54_standard.yaml")
        Potcar_names = list(potcar_symbols)

        AtomTypes_Potcar = [name.split("_")[0] for name in Potcar_names]

        AtomTypes = structure.symbol_set

        if set(AtomTypes) != set(AtomTypes_Potcar):
            raise OSError("Your POSCAR does not correspond to your POTCAR!")
        BASIS = loadfn(address_basis_file)["BASIS"]

        basis_functions = []
        list_forin = []
        for itype, type in enumerate(Potcar_names):
            if type not in BASIS:
                raise ValueError(
                    "You have to provide the basis for"
                    + str(type)
                    + "manually. We don't have any information on this POTCAR."
                )
            basis_functions.append(BASIS[type].split())
            tojoin = str(AtomTypes_Potcar[itype]) + " "
            tojoin2 = "".join(str(str(e) + " ") for e in BASIS[type].split())
            list_forin.append(str(tojoin + tojoin2))
        return list_forin

    @staticmethod
    def get_all_possible_basis_functions(
        structure: Structure,
        potcar_symbols: list,
        address_basis_file_min: str | None = None,
        address_basis_file_max: str | None = None,
    ):
        """
        Args:
            structure: Structure object
            potcar_symbols: list of the potcar symbols
            address_basis_file_min: path to file with the minimum required basis by the POTCAR
            address_basis_file_max: path to file with the largest possible basis of the POTCAR

        Returns: List of dictionaries that can be used to create new Lobsterin objects in
        standard_calculations_from_vasp_files as dict_for_basis
        """
        max_basis = Lobsterin.get_basis(
            structure=structure,
            potcar_symbols=potcar_symbols,
            address_basis_file=address_basis_file_max
            or os.path.join(MODULE_DIR, "lobster_basis/BASIS_PBE_54_max.yaml"),
        )
        min_basis = Lobsterin.get_basis(
            structure=structure,
            potcar_symbols=potcar_symbols,
            address_basis_file=address_basis_file_min
            or os.path.join(MODULE_DIR, "lobster_basis/BASIS_PBE_54_min.yaml"),
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
        POSCAR_input="POSCAR", POSCAR_output="POSCAR.lobster", symprec: float = 0.01
    ):
        """
        Writes a POSCAR with the standard primitive cell. This is needed to arrive at the correct kpath

        Args:
            POSCAR_input (str): filename of input POSCAR
            POSCAR_output (str): filename of output POSCAR
            symprec (float): precision to find symmetry
        """
        structure = Structure.from_file(POSCAR_input)
        kpath = HighSymmKpath(structure, symprec=symprec)
        new_structure = kpath.prim
        new_structure.to(fmt="POSCAR", filename=POSCAR_output)

    @staticmethod
    def write_KPOINTS(
        POSCAR_input: str = "POSCAR",
        KPOINTS_output="KPOINTS.lobster",
        reciprocal_density: int = 100,
        isym: int = -1,
        from_grid: bool = False,
        input_grid: Sequence[int] = (5, 5, 5),
        line_mode: bool = True,
        kpoints_line_density: int = 20,
        symprec: float = 0.01,
    ):
        """
        Writes a KPOINT file for lobster (only ISYM=-1 and ISYM=0 are possible), grids are gamma centered

        Args:
            POSCAR_input (str): path to POSCAR
            KPOINTS_output (str): path to output KPOINTS
            reciprocal_density (int): Grid density
            isym (int): either -1 or 0. Current Lobster versions only allow -1.
            from_grid (bool): If True KPOINTS will be generated with the help of a grid given in input_grid. Otherwise,
                they will be generated from the reciprocal_density
            input_grid (list): grid to generate the KPOINTS file
            line_mode (bool): If True, band structure will be generated
            kpoints_line_density (int): density of the lines in the band structure
            symprec (float): precision to determine symmetry
        """
        structure = Structure.from_file(POSCAR_input)
        if not from_grid:
            kpointgrid = Kpoints.automatic_density_by_vol(structure, reciprocal_density).kpts
            mesh = kpointgrid[0]
        else:
            mesh = input_grid

        # The following code is taken from: SpacegroupAnalyzer
        # we need to switch off symmetry here
        latt = structure.lattice.matrix
        positions = structure.frac_coords
        unique_species: list[Composition] = []
        zs = []
        magmoms = []

        for species, g in itertools.groupby(structure, key=lambda s: s.species):
            if species in unique_species:
                ind = unique_species.index(species)
                zs.extend([ind + 1] * len(tuple(g)))
            else:
                unique_species.append(species)
                zs.extend([len(unique_species)] * len(tuple(g)))

        for site in structure:
            if hasattr(site, "magmom"):
                magmoms.append(site.magmom)
            elif site.is_ordered and hasattr(site.specie, "spin"):
                magmoms.append(site.specie.spin)
            else:
                magmoms.append(0)

        # For now, we are setting magmom to zero. (Taken from INCAR class)
        cell = latt, positions, zs, magmoms
        # TODO: what about this shift?
        mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])

        # exit()
        # get the kpoints for the grid
        if isym == -1:
            kpts = []
            weights = []
            all_labels = []
            for gp in grid:
                kpts.append(gp.astype(float) / mesh)
                weights.append(float(1))
                all_labels.append("")
        elif isym == 0:
            # time reversal symmetry: k and -k are equivalent
            kpts = []
            weights = []
            all_labels = []
            newlist = [list(gp) for gp in list(grid)]
            mapping = []
            for gp in newlist:
                minus_gp = [-k for k in gp]
                if minus_gp in newlist and minus_gp not in [[0, 0, 0]]:
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
            ValueError("Only isym=-1 and isym=0 are allowed.")
        # line mode
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
        if isym == -1:
            comment = (
                "ISYM=-1, grid: " + str(mesh) if not line_mode else "ISYM=-1, grid: " + str(mesh) + " plus kpoint path"
            )
        elif isym == 0:
            comment = (
                "ISYM=0, grid: " + str(mesh) if not line_mode else "ISYM=0, grid: " + str(mesh) + " plus kpoint path"
            )

        KpointObject = Kpoints(
            comment=comment,
            style=Kpoints.supported_modes.Reciprocal,
            num_kpts=len(kpts),
            kpts=kpts,
            kpts_weights=weights,
            labels=all_labels,
        )

        KpointObject.write_file(filename=KPOINTS_output)

    @classmethod
    def from_file(cls, lobsterin: str):
        """
        Args:
            lobsterin (str): path to lobsterin

        Returns:
            Lobsterin object
        """
        with zopen(lobsterin, "rt") as f:
            data = f.read().split("\n")
        if len(data) == 0:
            raise OSError("lobsterin file contains no data.")
        Lobsterindict: dict[str, Any] = {}

        for datum in data:
            # will remove all comments to avoid complications
            raw_datum = datum.split("!")[0]
            raw_datum = raw_datum.split("//")[0]
            raw_datum = raw_datum.split("#")[0]
            raw_datum = raw_datum.split(" ")
            while "" in raw_datum:
                raw_datum.remove("")
            if len(raw_datum) > 1:
                # check which type of keyword this is, handle accordingly
                if raw_datum[0].lower() not in [datum2.lower() for datum2 in Lobsterin.LISTKEYWORDS]:
                    if raw_datum[0].lower() not in [datum2.lower() for datum2 in Lobsterin.FLOATKEYWORDS]:
                        if raw_datum[0].lower() not in Lobsterindict:
                            Lobsterindict[raw_datum[0].lower()] = " ".join(raw_datum[1:])
                        else:
                            raise ValueError("Same keyword " + str(raw_datum[0].lower()) + "twice!")
                    else:
                        if raw_datum[0].lower() not in Lobsterindict:
                            Lobsterindict[raw_datum[0].lower()] = float(raw_datum[1])
                        else:
                            raise ValueError("Same keyword " + str(raw_datum[0].lower()) + "twice!")
                else:
                    if raw_datum[0].lower() not in Lobsterindict:
                        Lobsterindict[raw_datum[0].lower()] = [" ".join(raw_datum[1:])]
                    else:
                        Lobsterindict[raw_datum[0].lower()].append(" ".join(raw_datum[1:]))
            elif len(raw_datum) > 0:
                Lobsterindict[raw_datum[0].lower()] = True

        return cls(Lobsterindict)

    @staticmethod
    def _get_potcar_symbols(POTCAR_input: str) -> list:
        """
        Will return the name of the species in the POTCAR

        Args:
            POTCAR_input(str): string to potcar file

        Returns:
            list of the names of the species in string format
        """
        potcar = Potcar.from_file(POTCAR_input)
        for pot in potcar:
            if pot.potential_type != "PAW":
                raise OSError("Lobster only works with PAW! Use different POTCARs")

        # Warning about a bug in lobster-4.1.0
        with zopen(POTCAR_input, "r") as f:
            data = f.read()
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
            raise OSError("We only have BASIS options for PBE so far")

        Potcar_names = [name["symbol"] for name in potcar.spec]
        return Potcar_names

    @classmethod
    def standard_calculations_from_vasp_files(
        cls,
        POSCAR_input: str = "POSCAR",
        INCAR_input: str = "INCAR",
        POTCAR_input: str | None = None,
        Vasprun_output: str = "vasprun.xml",
        dict_for_basis: dict | None = None,
        option: str = "standard",
    ):
        """
        Will generate Lobsterin with standard settings

        Args:
            POSCAR_input(str): path to POSCAR
            INCAR_input(str): path to INCAR
            POTCAR_input (str): path to POTCAR
            dict_for_basis (dict): can be provided: it should look the following:
                dict_for_basis={"Fe":'3p 3d 4s 4f', "C": '2s 2p'} and will overwrite all settings from POTCAR_input

            option (str): 'standard' will start a normal lobster run where COHPs, COOPs, DOS, CHARGE etc. will be
                calculated
                'standard_with_energy_range_from_vasprun' will start a normal lobster run for entire energy range
                of VASP static run. vasprun.xml file needs to be in current directory.
                'standard_from_projection' will start a normal lobster run from a projection
                'standard_with_fatband' will do a fatband calculation, run over all orbitals
                'onlyprojection' will only do a projection
                'onlydos' will only calculate a projected dos
                'onlycohp' will only calculate cohp
                'onlycoop' will only calculate coop
                'onlycohpcoop' will only calculate cohp and coop

        Returns:
            Lobsterin Object with standard settings
        """
        warnings.warn(
            "Always check and test the provided basis functions. The spilling of your Lobster calculation might help"
        )
        # warn that fatband calc cannot be done with tetrahedron method at the moment
        if option not in [
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
        ]:
            raise ValueError("The option is not valid!")

        Lobsterindict: dict[str, Any] = {}
        # this basis set covers most elements
        Lobsterindict["basisSet"] = "pbeVaspFit2015"
        # energies around e-fermi
        Lobsterindict["COHPstartEnergy"] = -35.0
        Lobsterindict["COHPendEnergy"] = 5.0

        if option in [
            "standard",
            "standard_with_energy_range_from_vasprun",
            "onlycohp",
            "onlycoop",
            "onlycobi",
            "onlycohpcoop",
            "onlycohpcoopcobi",
            "standard_with_fatband",
        ]:
            # every interaction with a distance of 6.0 is checked
            Lobsterindict["cohpGenerator"] = "from 0.1 to 6.0 orbitalwise"
            # the projection is saved
            Lobsterindict["saveProjectionToFile"] = True

        if option == "standard_from_projection":
            Lobsterindict["cohpGenerator"] = "from 0.1 to 6.0 orbitalwise"
            Lobsterindict["loadProjectionFromFile"] = True

        if option == "standard_with_energy_range_from_vasprun":
            Vr = Vasprun(Vasprun_output)
            Lobsterindict["COHPstartEnergy"] = round(min(Vr.complete_dos.energies - Vr.complete_dos.efermi), 4)
            Lobsterindict["COHPendEnergy"] = round(max(Vr.complete_dos.energies - Vr.complete_dos.efermi), 4)
            Lobsterindict["COHPSteps"] = len(Vr.complete_dos.energies)

        # TODO: add cobi here! might be relevant lobster version
        if option == "onlycohp":
            Lobsterindict["skipdos"] = True
            Lobsterindict["skipcoop"] = True
            Lobsterindict["skipPopulationAnalysis"] = True
            Lobsterindict["skipGrossPopulation"] = True
            # lobster-4.1.0
            Lobsterindict["skipcobi"] = True
            Lobsterindict["skipMadelungEnergy"] = True

        if option == "onlycoop":
            Lobsterindict["skipdos"] = True
            Lobsterindict["skipcohp"] = True
            Lobsterindict["skipPopulationAnalysis"] = True
            Lobsterindict["skipGrossPopulation"] = True
            # lobster-4.1.0
            Lobsterindict["skipcobi"] = True
            Lobsterindict["skipMadelungEnergy"] = True

        if option == "onlycohpcoop":
            Lobsterindict["skipdos"] = True
            Lobsterindict["skipPopulationAnalysis"] = True
            Lobsterindict["skipGrossPopulation"] = True
            # lobster-4.1.0
            Lobsterindict["skipcobi"] = True
            Lobsterindict["skipMadelungEnergy"] = True

        if option == "onlycohpcoopcobi":
            Lobsterindict["skipdos"] = True
            Lobsterindict["skipPopulationAnalysis"] = True
            Lobsterindict["skipGrossPopulation"] = True
            Lobsterindict["skipMadelungEnergy"] = True

        if option == "onlydos":
            Lobsterindict["skipcohp"] = True
            Lobsterindict["skipcoop"] = True
            Lobsterindict["skipPopulationAnalysis"] = True
            Lobsterindict["skipGrossPopulation"] = True
            # lobster-4.1.0
            Lobsterindict["skipcobi"] = True
            Lobsterindict["skipMadelungEnergy"] = True

        if option == "onlyprojection":
            Lobsterindict["skipdos"] = True
            Lobsterindict["skipcohp"] = True
            Lobsterindict["skipcoop"] = True
            Lobsterindict["skipPopulationAnalysis"] = True
            Lobsterindict["skipGrossPopulation"] = True
            Lobsterindict["saveProjectionToFile"] = True
            # lobster-4.1.0
            Lobsterindict["skipcobi"] = True
            Lobsterindict["skipMadelungEnergy"] = True

        if option == "onlycobi":
            Lobsterindict["skipdos"] = True
            Lobsterindict["skipcohp"] = True
            Lobsterindict["skipPopulationAnalysis"] = True
            Lobsterindict["skipGrossPopulation"] = True
            # lobster-4.1.0
            Lobsterindict["skipcobi"] = True
            Lobsterindict["skipMadelungEnergy"] = True

        if option == "onlymadelung":
            Lobsterindict["skipdos"] = True
            Lobsterindict["skipcohp"] = True
            Lobsterindict["skipcoop"] = True
            Lobsterindict["skipPopulationAnalysis"] = True
            Lobsterindict["skipGrossPopulation"] = True
            Lobsterindict["saveProjectionToFile"] = True
            # lobster-4.1.0
            Lobsterindict["skipcobi"] = True

        incar = Incar.from_file(INCAR_input)
        if incar["ISMEAR"] == 0:
            Lobsterindict["gaussianSmearingWidth"] = incar["SIGMA"]
        if incar["ISMEAR"] != 0 and option == "standard_with_fatband":
            raise ValueError("ISMEAR has to be 0 for a fatband calculation with Lobster")
        if dict_for_basis is not None:
            # dict_for_basis={"Fe":'3p 3d 4s 4f', "C": '2s 2p'}
            # will just insert this basis and not check with poscar
            basis = [key + " " + value for key, value in dict_for_basis.items()]
        elif POTCAR_input is not None:
            # get basis from POTCAR
            potcar_names = Lobsterin._get_potcar_symbols(POTCAR_input=POTCAR_input)

            basis = Lobsterin.get_basis(structure=Structure.from_file(POSCAR_input), potcar_symbols=potcar_names)
        else:
            raise ValueError("basis cannot be generated")
        Lobsterindict["basisfunctions"] = basis
        if option == "standard_with_fatband":
            Lobsterindict["createFatband"] = basis

        return cls(Lobsterindict)


def get_all_possible_basis_combinations(min_basis: list, max_basis: list) -> list:
    """

    Args:
        min_basis: list of basis entries: e.g., ['Si 3p 3s ']
        max_basis: list of basis entries: e.g., ['Si 3p 3s ']

    Returns: all possible combinations of basis functions, e.g. [['Si 3p 3s']]
    """
    max_basis_lists = [x.split() for x in max_basis]
    min_basis_lists = [x.split() for x in min_basis]

    # get all possible basis functions
    basis_dict: dict[str, dict] = {}
    for iel, el in enumerate(max_basis_lists):
        basis_dict[el[0]] = {"fixed": [], "variable": [], "combinations": []}
        for basis in el[1:]:
            if basis in min_basis_lists[iel]:
                basis_dict[el[0]]["fixed"].append(basis)
            if basis not in min_basis_lists[iel]:
                basis_dict[el[0]]["variable"].append(basis)
        for L in range(0, len(basis_dict[el[0]]["variable"]) + 1):
            for subset in itertools.combinations(basis_dict[el[0]]["variable"], L):
                basis_dict[el[0]]["combinations"].append(" ".join([el[0]] + basis_dict[el[0]]["fixed"] + list(subset)))

    list_basis = []
    for item in basis_dict.values():
        list_basis.append(item["combinations"])

    # get all combinations
    start_basis = list_basis[0]
    if len(list_basis) > 1:
        for el in list_basis[1:]:
            new_start_basis = []
            for elbasis in start_basis:
                for elbasis2 in el:
                    if not isinstance(elbasis, list):
                        new_start_basis.append([elbasis, elbasis2])
                    else:
                        new_start_basis.append(elbasis.copy() + [elbasis2])
            start_basis = new_start_basis
        return start_basis
    return [[basis] for basis in start_basis]
