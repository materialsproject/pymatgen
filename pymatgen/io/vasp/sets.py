"""
This module defines the VaspInputSet abstract base class and a concrete implementation for the parameters developed
and tested by the core team of pymatgen, including the Materials Virtual Lab, Materials Project and the MIT high
throughput project. The basic concept behind an input set is to specify a scheme to generate a consistent set of VASP
inputs from a structure without further user intervention. This ensures comparability across runs.

Read the following carefully before implementing new input sets:

1. 99% of what needs to be done can be done by specifying user_incar_settings to override some of the defaults of
   various input sets. Unless there is an extremely good reason to add a new set, **do not** add one. E.g., if you want
   to turn the Hubbard U off, just set "LDAU": False as a user_incar_setting.
2. All derivative input sets should inherit appropriate configurations (e.g., from MPRelaxSet), and more often than
   not, DictSet should be the superclass. Proper superclass delegation should be used where possible. In particular,
   you are not supposed to implement your own as_dict or from_dict for derivative sets unless you know what you are
   doing. Improper overriding the as_dict and from_dict protocols is the major cause of implementation headaches. If
   you need an example, look at how the MPStaticSet is initialized.

The above are recommendations. The following are **UNBREAKABLE** rules:

1. All input sets must take in a structure, list of structures or None as the first argument. If None, the input set
   should perform a stateless initialization and before any output can be written, a structure must be set.
2. user_incar_settings, user_kpoints_settings and user_<whatever>_settings are ABSOLUTE. Any new sets you implement
   must obey this. If a user wants to override your settings, you assume he knows what he is doing. Do not
   magically override user supplied settings. You can issue a warning if you think the user is wrong.
3. All input sets must save all supplied args and kwargs as instance variables. E.g., self.arg = arg and
   self.kwargs = kwargs in the __init__. This ensures the as_dict and from_dict work correctly.
"""

from __future__ import annotations

import abc
import itertools
import os
import re
import shutil
import warnings
from copy import deepcopy
from dataclasses import dataclass, field
from glob import glob
from itertools import chain
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Union
from zipfile import ZipFile

import numpy as np
from monty.dev import deprecated
from monty.io import zopen
from monty.json import MSONable
from monty.serialization import loadfn

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Element, PeriodicSite, SiteCollection, Species, Structure
from pymatgen.io.core import InputGenerator
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pymatgen.core.trajectory import Vector3D

MODULE_DIR = Path(__file__).resolve().parent
# TODO (janosh): replace with following line once PMG is py3.9+ only
# which guarantees absolute path to __file__
# reason: faster and shorter than pathlib (slow import and less OOP overhead)
# MODULE_DIR = os.path.dirname(__file__)


class VaspInputSet(InputGenerator, metaclass=abc.ABCMeta):
    """
    Base class representing a set of VASP input parameters with a structure
    supplied as init parameters. Typically, you should not inherit from this
    class. Start from DictSet or MPRelaxSet or MITRelaxSet.
    """

    _valid_potcars: Sequence[str] | None = None

    @property
    @abc.abstractmethod
    def incar(self):
        """Incar object."""

    @property
    @abc.abstractmethod
    def kpoints(self):
        """Kpoints object."""

    @property
    @abc.abstractmethod
    def poscar(self):
        """Poscar object."""

    @property
    def potcar_symbols(self):
        """List of POTCAR symbols."""

        elements = self.poscar.site_symbols
        potcar_symbols = []
        settings = self._config_dict["POTCAR"]

        if isinstance(settings[elements[-1]], dict):
            for el in elements:
                potcar_symbols.append(settings[el]["symbol"] if el in settings else el)
        else:
            for el in elements:
                potcar_symbols.append(settings.get(el, el))

        return potcar_symbols

    @property
    def potcar(self) -> Potcar:
        """Potcar object."""
        user_potcar_functional = self.user_potcar_functional
        potcar = Potcar(self.potcar_symbols, functional=user_potcar_functional)

        # warn if the selected POTCARs do not correspond to the chosen user_potcar_functional
        for p_single in potcar:
            if user_potcar_functional not in p_single.identify_potcar()[0]:
                warnings.warn(
                    f"POTCAR data with symbol {p_single.symbol} is not known by pymatgen to "
                    f"correspond with the selected {user_potcar_functional=}. This POTCAR "
                    f"is known to correspond with functionals {p_single.identify_potcar(mode='data')[0]}. "
                    "Please verify that you are using the right POTCARs!",
                    BadInputSetWarning,
                )

        return potcar

    @deprecated(message="get_vasp_input will be removed in a future version of pymatgen. Use get_input_set instead.")
    def get_vasp_input(self, structure=None) -> VaspInput:
        """Get a VaspInput object.

        Returns:
            VaspInput.
        """
        return self.get_input_set(structure=structure)

    def get_input_set(self, structure=None) -> VaspInput:
        """Get a VaspInput object.

        Returns:
            VaspInput.
        """
        if structure is not None:
            self.structure = structure
        return VaspInput(incar=self.incar, kpoints=self.kpoints, poscar=self.poscar, potcar=self.potcar)

    def write_input(
        self,
        output_dir: str,
        make_dir_if_not_present: bool = True,
        include_cif: bool = False,
        potcar_spec: bool = False,
        zip_output: bool = False,
    ) -> None:
        """
        Writes a set of VASP input to a directory.

        Args:
            output_dir (str): Directory to output the VASP input files
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            include_cif (bool): Whether to write a CIF file in the output
                directory for easier opening by VESTA.
            potcar_spec (bool): Instead of writing the POTCAR, write a "POTCAR.spec".
                This is intended to help sharing an input set with people who might
                not have a license to specific Potcar files. Given a "POTCAR.spec",
                the specific POTCAR file can be re-generated using pymatgen with the
                "generate_potcar" function in the pymatgen CLI.
            zip_output (bool): If True, output will be zipped into a file with the
                same name as the InputSet (e.g., MPStaticSet.zip)
        """
        if potcar_spec:
            if make_dir_if_not_present:
                os.makedirs(output_dir, exist_ok=True)

            with zopen(f"{output_dir}/POTCAR.spec", mode="wt") as file:
                file.write("\n".join(self.potcar_symbols))

            for key in ["INCAR", "POSCAR", "KPOINTS"]:
                if (val := getattr(self, key.lower())) is not None:
                    with zopen(os.path.join(output_dir, key), mode="wt") as file:
                        file.write(str(val))
        else:
            vasp_input = self.get_vasp_input()
            vasp_input.write_input(output_dir, make_dir_if_not_present=make_dir_if_not_present)

        cif_name = ""
        if include_cif:
            struct = vasp_input["POSCAR"].structure
            cif_name = f"{output_dir}/{struct.formula.replace(' ', '')}.cif"
            struct.to(filename=cif_name)

        if zip_output:
            filename = type(self).__name__ + ".zip"
            with ZipFile(os.path.join(output_dir, filename), mode="w") as zip_file:
                for file in ["INCAR", "POSCAR", "KPOINTS", "POTCAR", "POTCAR.spec", cif_name]:
                    try:
                        zip_file.write(os.path.join(output_dir, file), arcname=file)
                    except FileNotFoundError:
                        pass
                    try:
                        os.remove(os.path.join(output_dir, file))
                    except (FileNotFoundError, PermissionError, IsADirectoryError):
                        pass

    def as_dict(self, verbosity=2):
        """
        Args:
            verbosity: Verbosity for generated dict. If 1, structure is
            excluded.

        Returns:
            MSONable dict
        """
        dct = MSONable.as_dict(self)
        if verbosity == 1:
            dct.pop("structure", None)
        return dct


def _load_yaml_config(fname):
    config = loadfn(MODULE_DIR / (f"{fname}.yaml"))
    if "PARENT" in config:
        parent_config = _load_yaml_config(config["PARENT"])
        for k, v in parent_config.items():
            if k not in config:
                config[k] = v
            elif isinstance(v, dict):
                v_new = config.get(k, {})
                v_new.update(v)
                config[k] = v_new
    return config


UserPotcarFunctional = Union[
    Literal["PBE", "PBE_52", "PBE_54", "LDA", "LDA_52", "LDA_54", "PW91", "LDA_US", "PW91_US"], None
]


@dataclass
class DictSet(VaspInputSet):
    """
    Concrete implementation of VaspInputSet that is initialized from a dict
    settings. This allows arbitrary settings to be input. In general,
    this is rarely used directly unless there is a source of settings in yaml
    format (e.g., from a REST interface). It is typically used by other
    VaspInputSets for initialization.

    Special consideration should be paid to the way the MAGMOM initialization
    for the INCAR is done. The initialization differs depending on the type of
    structure and the configuration settings. The order in which the magmom is
    determined is as follows:

    1. If the site itself has a magmom setting (i.e. site.properties["magmom"] = float),
        that is used. This can be set with structure.add_site_property().
    2. If the species of the site has a spin setting, that is used. This can be set
        with structure.add_spin_by_element().
    3. If the species itself has a particular setting in the config file, that
       is used, e.g., Mn3+ may have a different magmom than Mn4+.
    4. Lastly, the element symbol itself is checked in the config file. If
       there are no settings, a default value of 0.6 is used.

    Args:
        structure (Structure): The Structure to create inputs for. If None, the input
            set is initialized without a Structure but one must be set separately before
            the inputs are generated.
        config_dict (dict): The config dictionary to use.
        files_to_transfer (dict): A dictionary of {filename: filepath}. This allows the
            transfer of files from a previous calculation.
        user_incar_settings (dict): User INCAR settings. This allows a user to override
            INCAR settings, e.g., setting a different MAGMOM for various elements or
            species. Note that in the new scheme, ediff_per_atom and hubbard_u are no
            longer args. Instead, the CONFIG supports EDIFF_PER_ATOM and EDIFF keys.
            The former scales with # of atoms, the latter does not. If both are present,
            EDIFF is preferred. To force such settings, just supply
            user_incar_settings={"EDIFF": 1e-5, "LDAU": False} for example. The keys
            'LDAUU', 'LDAUJ', 'LDAUL' are special cases since pymatgen defines different
            values depending on what anions are present in the structure, so these keys
            can be defined in one of two ways, e.g. either {"LDAUU":{"O":{"Fe":5}}} to
            set LDAUU for Fe to 5 in an oxide, or {"LDAUU":{"Fe":5}} to set LDAUU to 5
            regardless of the input structure. If a None value is given, that key is
            unset. For example, {"ENCUT": None} will remove ENCUT from the
            incar settings. Finally, KSPACING is a special setting and can be set to
            "auto" in which the KSPACING is set automatically based on the band gap.
        user_kpoints_settings (dict or Kpoints): Allow user to override kpoints setting
            by supplying a dict. E.g., {"reciprocal_density": 1000}. User can also
            supply Kpoints object.
        user_potcar_settings (dict): Allow user to override POTCARs. E.g., {"Gd":
            "Gd_3"}. This is generally not recommended.
        constrain_total_magmom (bool): Whether to constrain the total magmom (NUPDOWN in
            INCAR) to be the sum of the expected MAGMOM for all species.
        sort_structure (bool): Whether to sort the structure (using the default sort
            order of electronegativity) before generating input files. Defaults to True,
            the behavior you would want most of the time. This ensures that similar
            atomic species are grouped together.
        user_potcar_functional (str): Functional to use. Default (None) is to use the
            functional in the config dictionary. Valid values: "PBE", "PBE_52",
            "PBE_54", "LDA", "LDA_52", "LDA_54", "PW91", "LDA_US", "PW91_US".
        force_gamma (bool): Force gamma centered kpoint generation. Default (False) is
            to use the Automatic Density kpoint scheme, which will use the Gamma
            centered generation scheme for hexagonal cells, and Monkhorst-Pack otherwise.
        reduce_structure (None/str): Before generating the input files, generate the
            reduced structure. Default (None), does not alter the structure. Valid
            values: None, "niggli", "LLL".
        vdw: Adds default parameters for van-der-Waals functionals supported by VASP to
            INCAR. Supported functionals are: DFT-D2, undamped DFT-D3, DFT-D3 with
            Becke-Jonson damping, Tkatchenko-Scheffler, Tkatchenko-Scheffler with
            iterative Hirshfeld partitioning, MBD@rSC, dDsC, Dion's vdW-DF, DF2, optPBE,
            optB88, optB86b and rVV10.
        use_structure_charge (bool): If set to True, then the overall charge of the
            structure (structure.charge) is used to set the NELECT variable in the
            INCAR. Default is False.
        standardize (float): Whether to standardize to a primitive standard cell.
            Defaults to False.
        sym_prec (float): Tolerance for symmetry finding.
        international_monoclinic (bool): Whether to use international convention (vs
            Curtarolo) for monoclinic. Defaults True.
        validate_magmom (bool): Ensure that the missing magmom values are filled in with
            the VASP default value of 1.0.
        inherit_incar (bool): Whether to inherit INCAR settings from previous
            calculation. This might be useful to port Custodian fixes to child jobs but
            can also be dangerous e.g. when switching from GGA to meta-GGA or relax to
            static jobs. Defaults to True.
        auto_ismear (bool): If true, the values for ISMEAR and SIGMA will be set
            automatically depending on the bandgap of the system. If the bandgap is not
            known (e.g., there is no previous VASP directory) then ISMEAR=0 and
            SIGMA=0.2; if the bandgap is zero (a metallic system) then ISMEAR=2 and
            SIGMA=0.2; if the system is an insulator, then ISMEAR=-5 (tetrahedron
            smearing). Note, this only works when generating the input set from a
            previous VASP directory.
        bandgap_tol (float): Tolerance for determining if a system is metallic when
            KSPACING is set to "auto". If the bandgap is less than this value, the
            system is considered metallic. Defaults to 1e-4 (eV).
        bandgap (float): Used for determining KSPACING if KSPACING == "auto" or
            ISMEAR if auto_ismear == True. Set automatically when using from_prev_calc.
        prev_incar (str or dict): Previous INCAR used for setting parent INCAR when
            inherit_incar == True. Set automatically when using from_prev_calc.
        prev_kpoints (str or Kpoints): Previous Kpoints. Set automatically when using
            from_prev_calc.
    """

    structure: Structure | None = None
    config_dict: dict = field(default_factory=dict)
    files_to_transfer: dict = field(default_factory=dict)
    user_incar_settings: dict = field(default_factory=dict)
    user_kpoints_settings: dict = field(default_factory=dict)
    user_potcar_settings: dict = field(default_factory=dict)
    constrain_total_magmom: bool = False
    sort_structure: bool = True
    user_potcar_functional: UserPotcarFunctional = None
    force_gamma: bool = False
    reduce_structure: Literal["niggli", "LLL"] | None = None
    vdw: str | None = None
    use_structure_charge: bool = False
    standardize: bool = False
    sym_prec: float = 0.1
    international_monoclinic: bool = True
    validate_magmom: bool = True
    inherit_incar: bool | list[str] = False
    auto_ismear: bool = False
    bandgap_tol: float = 1e-4
    bandgap: float | None = None
    prev_incar: str | dict | None = None
    prev_kpoints: str | Kpoints | None = None

    def __post_init__(self):
        """Perform validation"""
        if (valid_potcars := self._valid_potcars) and self.user_potcar_functional not in valid_potcars:
            raise ValueError(f"Invalid {self.user_potcar_functional=}, must be one of {valid_potcars}")

        if hasattr(self, "CONFIG"):
            self.config_dict = self.CONFIG

        self._config_dict = deepcopy(self.config_dict)

        # these have been left to stay consistent with previous API
        self.user_incar_settings = self.user_incar_settings or {}
        self.user_kpoints_settings = self.user_kpoints_settings or {}

        self.vdw = self.vdw.lower() if isinstance(self.vdw, str) else self.vdw
        if self.user_incar_settings.get("KSPACING") and self.user_kpoints_settings is not None:
            warnings.warn(
                "You have specified KSPACING and also supplied kpoints "
                "settings. KSPACING only has effect when there is no "
                "KPOINTS file. Since both settings were given, pymatgen"
                "will generate a KPOINTS file and ignore KSPACING."
                "Remove the `user_kpoints_settings` argument to enable KSPACING.",
                BadInputSetWarning,
            )

        if self.vdw:
            vdw_par = loadfn(str(MODULE_DIR / "vdW_parameters.yaml"))
            try:
                self._config_dict["INCAR"].update(vdw_par[self.vdw])
            except KeyError:
                raise KeyError(
                    f"Invalid or unsupported van-der-Waals functional. Supported functionals are {', '.join(vdw_par)}."
                )
        # 'or' case reads the POTCAR_FUNCTIONAL from the .yaml
        self.user_potcar_functional: UserPotcarFunctional = self.user_potcar_functional or self._config_dict.get(
            "POTCAR_FUNCTIONAL", "PBE"
        )

        # warn if a user is overriding POTCAR_FUNCTIONAL
        if self.user_potcar_functional != self._config_dict.get("POTCAR_FUNCTIONAL", "PBE"):
            warnings.warn(
                "Overriding the POTCAR functional is generally not recommended "
                " as it significantly affect the results of calculations and "
                "compatibility with other calculations done with the same "
                "input set. Note that some POTCAR symbols specified in "
                "the configuration file may not be available in the selected "
                "functional.",
                BadInputSetWarning,
            )

        if self.user_potcar_settings:
            warnings.warn(
                "Overriding POTCARs is generally not recommended as it "
                "significantly affect the results of calculations and "
                "compatibility with other calculations done with the same "
                "input set. In many instances, it is better to write a "
                "subclass of a desired input set and override the POTCAR in "
                "the subclass to be explicit on the differences.",
                BadInputSetWarning,
            )
            for key, val in self.user_potcar_settings.items():
                self._config_dict["POTCAR"][key] = val

        if not isinstance(self.structure, Structure):
            self._structure = None
        else:
            self.structure = self.structure

        if isinstance(self.prev_incar, (Path, str)):
            self.prev_incar = Incar.from_file(self.prev_incar)

        if isinstance(self.prev_kpoints, (Path, str)):
            self.prev_kpoints = Kpoints.from_file(self.prev_kpoints)

        self.prev_vasprun = None
        self.prev_outcar = None

    @property  # type: ignore
    def structure(self) -> Structure:
        """Structure"""
        return self._structure

    @structure.setter
    def structure(self, structure: Structure | None) -> None:
        if not hasattr(self, "_config_dict"):
            self._structure = structure
            return

        if isinstance(structure, SiteCollection):  # could be Structure or Molecule
            if self.user_potcar_functional == "PBE_54" and "W" in structure.symbol_set:
                # when using 5.4 POTCARs, default Tungsten POTCAR to W_Sv but still allow user to override
                self.user_potcar_settings = {"W": "W_sv", **(self.user_potcar_settings or {})}
            if self.reduce_structure:
                structure = structure.get_reduced_structure(self.reduce_structure)
            if self.sort_structure:
                structure = structure.get_sorted_structure()
            if self.validate_magmom:
                get_valid_magmom_struct(structure, spin_mode="auto")

            struct_has_Yb = any(specie.symbol == "Yb" for site in structure for specie in site.species)
            potcar_settings = self._config_dict.get("POTCAR", {})
            if self.user_potcar_settings:
                potcar_settings.update(self.user_potcar_settings)
            uses_Yb_2_psp = potcar_settings.get("Yb", None) == "Yb_2"
            if struct_has_Yb and uses_Yb_2_psp:
                warnings.warn(
                    "The structure contains Ytterbium (Yb) and this InputSet uses the Yb_2 PSP.\n"
                    "Yb_2 is known to often give bad results since Yb has oxidation state 3+ in most compounds.\n"
                    "See https://github.com/materialsproject/pymatgen/issues/2968 for details.",
                    BadInputSetWarning,
                )
            if self.standardize and self.sym_prec:
                structure = standardize_structure(
                    structure,
                    sym_prec=self.sym_prec,
                    international_monoclinic=self.international_monoclinic,
                )
        self._structure = structure

    def get_input_set(
        self,
        structure: Structure | None = None,
        prev_dir: str | Path | None = None,
        potcar_spec: bool = False,
    ) -> VaspInput:
        """
        Get a VASP input set.

        Note, if both ``structure`` and ``prev_dir`` are set, then the structure
        specified will be preferred over the final structure from the last VASP run.

        Args:
            structure (Structure): A structure.
            prev_dir (str or Path): A previous directory to generate the input set from.
            potcar_spec (bool): Instead of generating a Potcar object, use a list of
                potcar symbols. This will be written as a "POTCAR.spec" file. This is
                intended to help sharing an input set with people who might not have a
                license to specific Potcar files. Given a "POTCAR.spec", the specific
                POTCAR file can be re-generated using pymatgen with the
                "generate_potcar" function in the pymatgen CLI.

        Returns:
            VaspInput: A VASP input object.
        """
        if structure is None and prev_dir is None and self.structure is None:
            raise ValueError("Either structure or prev_dir must be set.")

        self._set_previous(prev_dir)

        if structure is not None:
            self.structure = structure

        return VaspInput(
            incar=self.incar,
            kpoints=self.kpoints,
            poscar=self.poscar,
            potcar=self.potcar_symbols if potcar_spec else self.potcar,
        )

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        return {}

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Get updates to the kpoints configuration for this calculation type.

        Note, these updates will be ignored if the user has set user_kpoint_settings.

        Returns:
            dict or Kpoints: A dictionary of updates to apply to the KPOINTS config
                or a Kpoints object.
        """
        return {}

    def _set_previous(self, prev_dir: str | Path | None = None):
        """Load previous calculation outputs."""
        if prev_dir is not None:
            vasprun, outcar = get_vasprun_outcar(prev_dir)
            self.prev_vasprun = vasprun
            self.prev_outcar = outcar
            self.prev_incar = vasprun.incar
            self.prev_kpoints = Kpoints.from_dict(vasprun.kpoints.as_dict())

            if vasprun.efermi is None:
                # VASP doesn't output efermi in vasprun if IBRION = 1
                vasprun.efermi = outcar.efermi

            bs = vasprun.get_band_structure(efermi="smart")
            self.bandgap = 0 if bs.is_metal() else bs.get_band_gap()["energy"]
            self.structure = get_structure_from_prev_run(vasprun, outcar)

    @property
    def incar(self) -> Incar:
        """Get the INCAR."""
        if self.structure is None:
            raise RuntimeError("No structure is associated with the input set!")

        prev_incar: dict[str, Any] = {}
        if self.inherit_incar is True and self.prev_incar:
            prev_incar = self.prev_incar  # type: ignore
        elif isinstance(self.inherit_incar, (list, tuple)) and self.prev_incar:
            prev_incar = {k: self.prev_incar[k] for k in self.inherit_incar if k in self.prev_incar}  # type: ignore

        incar_updates = self.incar_updates
        settings = dict(self._config_dict["INCAR"])
        auto_updates = {}

        # apply updates from input set generator to SETTINGS
        _apply_incar_updates(settings, incar_updates)

        # apply user incar settings to SETTINGS not to INCAR
        _apply_incar_updates(settings, self.user_incar_settings)

        # generate incar
        structure = self.structure
        comp = structure.composition
        elements = sorted((el for el in comp.elements if comp[el] > 0), key=lambda e: e.X)
        most_electro_neg = elements[-1].symbol
        poscar = Poscar(structure)
        hubbard_u = settings.get("LDAU", False)
        incar = Incar()

        for k, v in settings.items():
            if k == "MAGMOM":
                mag = []
                for site in structure:
                    if hasattr(site, "magmom"):
                        mag.append(site.magmom)
                    elif getattr(site.specie, "spin", None) is not None:
                        mag.append(site.specie.spin)
                    elif str(site.specie) in v:
                        if site.specie.symbol == "Co" and v[str(site.specie)] <= 1.0:
                            warnings.warn(
                                "Co without an oxidation state is initialized as low spin by default in Pymatgen. "
                                "If this default behavior is not desired, please set the spin on the magmom on the "
                                "site directly to ensure correct initialization."
                            )
                        mag.append(v.get(str(site.specie)))
                    else:
                        if site.specie.symbol == "Co":
                            warnings.warn(
                                "Co without an oxidation state is initialized as low spin by default in Pymatgen. "
                                "If this default behavior is not desired, please set the spin on the magmom on the "
                                "site directly to ensure correct initialization."
                            )
                        mag.append(v.get(site.specie.symbol, 0.6))
                incar[k] = mag
            elif k in ("LDAUU", "LDAUJ", "LDAUL"):
                if hubbard_u:
                    if hasattr(structure[0], k.lower()):
                        m = {site.specie.symbol: getattr(site, k.lower()) for site in structure}
                        incar[k] = [m[sym] for sym in poscar.site_symbols]
                        # lookup specific LDAU if specified for most_electroneg atom
                    elif most_electro_neg in v and isinstance(v[most_electro_neg], dict):
                        incar[k] = [v[most_electro_neg].get(sym, 0) for sym in poscar.site_symbols]
                        # else, use fallback LDAU value if it exists
                    else:
                        incar[k] = [
                            v.get(sym, 0) if isinstance(v.get(sym, 0), (float, int)) else 0
                            for sym in poscar.site_symbols
                        ]
            elif k.startswith("EDIFF") and k != "EDIFFG":
                if "EDIFF" not in settings and k == "EDIFF_PER_ATOM":
                    incar["EDIFF"] = float(v) * len(structure)
                else:
                    incar["EDIFF"] = float(settings["EDIFF"])
            elif k == "KSPACING" and v == "auto":
                # default to metal if no prev calc available
                bandgap = 0 if self.bandgap is None else self.bandgap
                incar[k] = auto_kspacing(bandgap, self.bandgap_tol)
            else:
                incar[k] = v
        has_u = hubbard_u and sum(incar["LDAUU"]) > 0
        if not has_u:
            for key in list(incar):
                if key.startswith("LDAU"):
                    del incar[key]

        # Modify LMAXMIX if you have d or f electrons present. Note that if the user
        # explicitly sets LMAXMIX in settings it will override this logic.
        # Previously, this was only set if Hubbard U was enabled as per the VASP manual
        # but following an investigation it was determined that this would lead to a
        # significant difference between SCF -> NonSCF even without Hubbard U enabled.
        # Thanks to Andrew Rosen for investigating and reporting.
        if "LMAXMIX" not in settings:
            # contains f-electrons
            if any(el.Z > 56 for el in structure.composition):
                incar["LMAXMIX"] = 6
            # contains d-electrons
            elif any(el.Z > 20 for el in structure.composition):
                incar["LMAXMIX"] = 4

        # Warn user about LASPH for +U, meta-GGAs, hybrids, and vdW-DF
        if not incar.get("LASPH", False) and (
            incar.get("METAGGA")
            or incar.get("LHFCALC", False)
            or incar.get("LDAU", False)
            or incar.get("LUSE_VDW", False)
        ):
            warnings.warn("LASPH = True should be set for +U, meta-GGAs, hybrids, and vdW-DFT", BadInputSetWarning)

        # apply previous incar settings, be careful not to override user_incar_settings
        # or the settings updates from the specific input set implementations
        # also skip LDAU/MAGMOM as structure may have changed.
        skip = list(self.user_incar_settings) + list(incar_updates)
        skip += ["MAGMOM", "NUPDOWN", "LDAUU", "LDAUL", "LDAUJ"]
        _apply_incar_updates(incar, prev_incar, skip=skip)

        if self.constrain_total_magmom:
            nupdown = sum(mag if abs(mag) > 0.6 else 0 for mag in incar["MAGMOM"])
            if abs(nupdown - round(nupdown)) > 1e-5:
                warnings.warn(
                    "constrain_total_magmom was set to True, but the sum of MAGMOM "
                    "values is not an integer. NUPDOWN is meant to set the spin "
                    "multiplet and should typically be an integer. You are likely "
                    "better off changing the values of MAGMOM or simply setting "
                    "NUPDOWN directly in your INCAR settings.",
                    UserWarning,
                    stacklevel=1,
                )
            auto_updates["NUPDOWN"] = nupdown

        if self.use_structure_charge:
            auto_updates["NELECT"] = self.nelect

        # Check that ALGO is appropriate
        if incar.get("LHFCALC", False) is True and incar.get("ALGO", "Normal") not in ["Normal", "All", "Damped"]:
            warnings.warn(
                "Hybrid functionals only support Algo = All, Damped, or Normal.",
                BadInputSetWarning,
            )

        if self.auto_ismear:
            if self.bandgap is None:
                # don't know if we are a metal or insulator so set ISMEAR and SIGMA to
                # be safe with the most general settings
                auto_updates.update(ISMEAR=2, SIGMA=0.2)
            elif self.bandgap <= self.bandgap_tol:
                auto_updates.update(ISMEAR=2, SIGMA=0.2)  # metal
            else:
                auto_updates.update(ISMEAR=-5, SIGMA=0.05)  # insulator

        kpoints = self.kpoints
        if kpoints is not None:
            # unset KSPACING as we are using a KPOINTS file
            incar.pop("KSPACING", None)
        elif "KSPACING" in incar and "KSPACING" not in self.user_incar_settings and "KSPACING" in prev_incar:
            # prefer to inherit KSPACING from previous INCAR if it exists
            # TODO: Is that we actually want to do? Copied from current pymatgen inputsets
            incar["KSPACING"] = prev_incar["KSPACING"]

        # apply updates from auto options, careful not to override user_incar_settings
        _apply_incar_updates(incar, auto_updates, skip=list(self.user_incar_settings))

        # Remove unused INCAR parameters
        _remove_unused_incar_params(incar, skip=list(self.user_incar_settings))

        # Ensure adequate number of KPOINTS are present for the tetrahedron method
        # (ISMEAR=-5). If KSPACING is in the INCAR file the number of kpoints is not
        # known before calling VASP, but a warning is raised when the KSPACING value is
        # > 0.5 (2 reciprocal Angstrom). An error handler in Custodian is available to
        # correct overly large KSPACING values (small number of kpoints) if necessary.
        if kpoints is not None and np.prod(kpoints.kpts) < 4 and incar.get("ISMEAR", 0) == -5:
            incar["ISMEAR"] = 0

        if self.user_incar_settings.get("KSPACING", 0) > 0.5 and incar.get("ISMEAR", 0) == -5:
            warnings.warn(
                "Large KSPACING value detected with ISMEAR = -5. Ensure that VASP "
                "generates an adequate number of KPOINTS, lower KSPACING, or "
                "set ISMEAR = 0",
                BadInputSetWarning,
            )

        if all(k.is_metal for k in structure.composition) and incar.get("NSW", 0) > 0 and incar.get("ISMEAR", 1) < 1:
            warnings.warn(
                "Relaxation of likely metal with ISMEAR < 1 detected. See VASP "
                "recommendations on ISMEAR for metals.",
                BadInputSetWarning,
            )

        return incar

    @property
    def poscar(self) -> Poscar:
        """Poscar"""
        if self.structure is None:
            raise RuntimeError("No structure is associated with the input set!")
        return Poscar(self.structure)

    @property
    def potcar_functional(self) -> UserPotcarFunctional:
        """Returns the functional used for POTCAR generation."""
        return self.user_potcar_functional

    @property
    def nelect(self) -> float:
        """Gets the default number of electrons for a given structure."""
        if self.structure is None:
            raise RuntimeError("No structure is associated with the input set!")

        n_electrons_by_element = {p.element: p.nelectrons for p in self.potcar}
        n_elect = sum(
            num_atoms * n_electrons_by_element[el.symbol] for el, num_atoms in self.structure.composition.items()
        )

        if self.use_structure_charge:
            return n_elect - self.structure.charge
        return n_elect

    @property
    def kpoints(self) -> Kpoints | None:
        """Get the kpoints file."""
        if self.structure is None:
            raise RuntimeError("No structure is associated with the input set!")

        if (
            self.user_incar_settings.get("KSPACING", None) is not None
            or self.incar_updates.get("KSPACING", None) is not None
            or self._config_dict["INCAR"].get("KSPACING", None) is not None
        ) and self.user_kpoints_settings == {}:
            # If KSPACING specified then always use this over k-points
            return None

        # use user setting if set otherwise default to base config settings
        kpoints_updates = self.kpoints_updates
        if self.user_kpoints_settings != {}:
            kconfig = deepcopy(self.user_kpoints_settings)
        elif isinstance(kpoints_updates, Kpoints):
            return kpoints_updates
        elif kpoints_updates != {}:
            kconfig = kpoints_updates
        else:
            kconfig = deepcopy(self._config_dict.get("KPOINTS", {}))

        if isinstance(kconfig, Kpoints):
            return kconfig

        explicit = (
            kconfig.get("explicit")
            or len(kconfig.get("added_kpoints", [])) > 0
            or "zero_weighted_reciprocal_density" in kconfig
            or "zero_weighted_line_density" in kconfig
        )
        # handle length generation first as this doesn't support any additional options
        if kconfig.get("length"):
            if explicit:
                raise ValueError(
                    "length option cannot be used with explicit k-point generation, "
                    "added_kpoints, or zero weighted k-points."
                )
            # If length is in kpoints settings use Kpoints.automatic
            return Kpoints.automatic(kconfig["length"])

        base_kpoints = None
        if kconfig.get("line_density"):
            # handle line density generation
            kpath = HighSymmKpath(self.structure, **kconfig.get("kpath_kwargs", {}))
            frac_k_points, k_points_labels = kpath.get_kpoints(
                line_density=kconfig["line_density"], coords_are_cartesian=False
            )
            base_kpoints = Kpoints(
                comment="Non SCF run along symmetry lines",
                style=Kpoints.supported_modes.Reciprocal,
                num_kpts=len(frac_k_points),
                kpts=frac_k_points,
                labels=k_points_labels,
                kpts_weights=[1] * len(frac_k_points),
            )
        elif kconfig.get("grid_density") or kconfig.get("reciprocal_density"):
            # handle regular weighted k-point grid generation
            if kconfig.get("grid_density"):
                base_kpoints = Kpoints.automatic_density(self.structure, int(kconfig["grid_density"]), self.force_gamma)
            elif kconfig.get("reciprocal_density"):
                density = kconfig["reciprocal_density"]
                base_kpoints = Kpoints.automatic_density_by_vol(self.structure, density, self.force_gamma)
            if explicit:
                sga = SpacegroupAnalyzer(self.structure, symprec=self.sym_prec)
                mesh = sga.get_ir_reciprocal_mesh(base_kpoints.kpts[0])  # type: ignore
                base_kpoints = Kpoints(
                    comment="Uniform grid",
                    style=Kpoints.supported_modes.Reciprocal,
                    num_kpts=len(mesh),
                    kpts=[i[0] for i in mesh],
                    kpts_weights=[i[1] for i in mesh],
                )
            else:
                # if not explicit that means no other options have been specified
                # so we can return the k-points as is
                return base_kpoints

        zero_weighted_kpoints = None
        if kconfig.get("zero_weighted_line_density"):
            # zero_weighted k-points along line mode path
            kpath = HighSymmKpath(self.structure)
            frac_k_points, k_points_labels = kpath.get_kpoints(
                line_density=kconfig["zero_weighted_line_density"],
                coords_are_cartesian=False,
            )
            zero_weighted_kpoints = Kpoints(
                comment="Hybrid run along symmetry lines",
                style=Kpoints.supported_modes.Reciprocal,
                num_kpts=len(frac_k_points),
                kpts=frac_k_points,
                labels=k_points_labels,
                kpts_weights=[0] * len(frac_k_points),
            )
        elif kconfig.get("zero_weighted_reciprocal_density"):
            zero_weighted_kpoints = Kpoints.automatic_density_by_vol(
                self.structure, kconfig["zero_weighted_reciprocal_density"], self.force_gamma
            )
            sga = SpacegroupAnalyzer(self.structure, symprec=self.sym_prec)
            mesh = sga.get_ir_reciprocal_mesh(zero_weighted_kpoints.kpts[0])
            zero_weighted_kpoints = Kpoints(
                comment="Uniform grid",
                style=Kpoints.supported_modes.Reciprocal,
                num_kpts=len(mesh),
                kpts=[i[0] for i in mesh],
                kpts_weights=[0 for i in mesh],
            )

        added_kpoints = None
        if kconfig.get("added_kpoints"):
            points: list = kconfig.get("added_kpoints")  # type: ignore
            added_kpoints = Kpoints(
                comment="Specified k-points only",
                style=Kpoints.supported_modes.Reciprocal,
                num_kpts=len(points),
                kpts=points,
                labels=["user-defined"] * len(points),
                kpts_weights=[0] * len(points),
            )

        if base_kpoints and not (added_kpoints or zero_weighted_kpoints):
            return base_kpoints
        if added_kpoints and not (base_kpoints or zero_weighted_kpoints):
            return added_kpoints

        # do some sanity checking
        if "line_density" in kconfig and zero_weighted_kpoints:
            raise ValueError("Cannot combine line_density and zero weighted k-points options")
        if zero_weighted_kpoints and not base_kpoints:
            raise ValueError("Zero weighted k-points must be used with reciprocal_density or grid_density options")
        if not (base_kpoints or zero_weighted_kpoints or added_kpoints):
            raise ValueError(
                "Invalid k-point generation algo. Supported Keys are 'grid_density' "
                "for Kpoints.automatic_density generation, 'reciprocal_density' for "
                "KPoints.automatic_density_by_vol generation, 'length' for "
                "Kpoints.automatic generation, 'line_density' for line mode generation,"
                " 'added_kpoints' for specific k-points to include, "
                " 'zero_weighted_reciprocal_density' for a zero weighted uniform mesh,"
                " or 'zero_weighted_line_density' for a zero weighted line mode mesh."
            )

        return _combine_kpoints(base_kpoints, zero_weighted_kpoints, added_kpoints)  # type: ignore

    @property
    def potcar(self) -> Potcar:
        """Potcar object"""
        if self.structure is None:
            raise RuntimeError("No structure is associated with the input set!")
        return super().potcar

    def estimate_nbands(self) -> int:
        """
        Estimate the number of bands that VASP will initialize a
        calculation with by default. Note that in practice this
        can depend on # of cores (if not set explicitly).
        Note that this formula is slightly different than the formula on the VASP wiki
        (as of July 2023). This is because the formula in the source code (`main.F`) is
        slightly different than what is on the wiki.
        """
        if self.structure is None:
            raise RuntimeError("No structure is associated with the input set!")

        n_ions = len(self.structure)

        if self.incar["ISPIN"] == 1:  # per the VASP source, if non-spin polarized ignore n_mag
            n_mag = 0
        else:  # otherwise set equal to sum of total magmoms
            n_mag = sum(self.incar["MAGMOM"])
            n_mag = np.floor((n_mag + 1) / 2)

        possible_val_1 = np.floor((self.nelect + 2) / 2) + max(np.floor(n_ions / 2), 3)
        possible_val_2 = np.floor(self.nelect * 0.6)

        n_bands = max(possible_val_1, possible_val_2) + n_mag

        if self.incar.get("LNONCOLLINEAR") is True:
            n_bands = n_bands * 2

        if n_par := self.incar.get("NPAR"):
            n_bands = (np.floor((n_bands + n_par - 1) / n_par)) * n_par

        return int(n_bands)

    def override_from_prev_calc(self, prev_calc_dir="."):
        """
        Update the input set to include settings from a previous calculation.

        Args:
            prev_calc_dir (str): The path to the previous calculation directory.

        Returns:
            The input set with the settings (structure, k-points, incar, etc)
            updated using the previous VASP run.
        """
        self._set_previous(prev_calc_dir)

        if self.standardize:
            warnings.warn(
                "Use of standardize=True with from_prev_run is not "
                "recommended as there is no guarantee the copied "
                "files will be appropriate for the standardized "
                "structure."
            )

        files_to_transfer = {}
        if getattr(self, "copy_chgcar", False):
            chgcars = sorted(glob(str(Path(prev_calc_dir) / "CHGCAR*")))
            if chgcars:
                files_to_transfer["CHGCAR"] = str(chgcars[-1])

        if getattr(self, "copy_wavecar", False):
            for fname in ("WAVECAR", "WAVEDER", "WFULL"):
                wavecar_files = sorted(glob(str(Path(prev_calc_dir) / (fname + "*"))))
                if wavecar_files:
                    if fname == "WFULL":
                        for f in wavecar_files:
                            fname = Path(f).name
                            fname = fname.split(".")[0]
                            files_to_transfer[fname] = f
                    else:
                        files_to_transfer[fname] = str(wavecar_files[-1])

        self.files_to_transfer.update(files_to_transfer)
        return self

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, **kwargs):
        """
        Generate a set of VASP input files for static calculations from a
        directory of previous VASP run.

        Args:
            prev_calc_dir (str): Directory containing the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            **kwargs: All kwargs supported by MPStaticSet, other than prev_incar
                and prev_structure and prev_kpoints which are determined from
                the prev_calc_dir.
        """
        input_set = cls(_dummy_structure, **kwargs)
        return input_set.override_from_prev_calc(prev_calc_dir=prev_calc_dir)

    def __str__(self):
        return type(self).__name__

    def __repr__(self):
        return type(self).__name__

    def write_input(
        self,
        output_dir: str,
        make_dir_if_not_present: bool = True,
        include_cif: bool = False,
        potcar_spec: bool = False,
        zip_output: bool = False,
    ):
        """
        Writes out all input to a directory.

        Args:
            output_dir (str): Directory to output the VASP input files
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            include_cif (bool): Whether to write a CIF file in the output
                directory for easier opening by VESTA.
            potcar_spec (bool): Instead of writing the POTCAR, write a "POTCAR.spec".
                This is intended to help sharing an input set with people who might
                not have a license to specific Potcar files. Given a "POTCAR.spec",
                the specific POTCAR file can be re-generated using pymatgen with the
                "generate_potcar" function in the pymatgen CLI.
            zip_output (bool): Whether to zip each VASP input file written to the output directory.
        """
        super().write_input(
            output_dir=output_dir,
            make_dir_if_not_present=make_dir_if_not_present,
            include_cif=include_cif,
            potcar_spec=potcar_spec,
            zip_output=zip_output,
        )
        for k, v in self.files_to_transfer.items():
            with zopen(v, "rb") as fin, zopen(str(Path(output_dir) / k), "wb") as fout:
                shutil.copyfileobj(fin, fout)

    def calculate_ng(
        self,
        max_prime_factor: int = 7,
        must_inc_2: bool = True,
        custom_encut: float | None = None,
        custom_prec: str | None = None,
    ) -> tuple:
        """
        Calculates the NGX, NGY, and NGZ values using the information available in the INCAR and POTCAR
        This is meant to help with making initial guess for the FFT grid so we can interact with the Charge density API.

        Args:
            max_prime_factor (int): the valid prime factors of the grid size in each direction
                VASP has many different setting for this to handle many compiling options.
                For typical MPI options all prime factors up to 7 are allowed
            must_inc_2 (bool): Whether 2 must be a prime factor of the result. Defaults to True.
            custom_encut (float | None): Calculates the FFT grid parameters using a custom
                ENCUT that may be different from what is generated by the input set. Defaults to None.
                Do *not* use this unless you know what you are doing.
            custom_prec (str | None): Calculates the FFT grid parameters using a custom prec
                that may be different from what is generated by the input set. Defaults to None.
                Do *not* use this unless you know what you are doing.
        """
        # TODO throw error for Ultrasoft potentials

        _RYTOEV = 13.605826
        _AUTOA = 0.529177249

        # TODO Only do this for VASP 6 for now. Older version require more advanced logic

        if custom_encut is not None:
            encut = custom_encut
        else:
            if self.incar.get("ENCUT", 0) > 0:
                encut = self.incar["ENCUT"]  # get the ENCUT val
            else:
                encut = max(i_species.enmax for i_species in self.get_vasp_input()["POTCAR"])

        # PREC=Normal is VASP default
        PREC = self.incar.get("PREC", "Normal") if custom_prec is None else custom_prec

        # Check for unsupported / invalid PREC tags
        if PREC[0].lower() in {"l", "m", "h"}:
            raise NotImplementedError(
                "PREC = LOW/MEDIUM/HIGH from VASP 4.x and not supported, Please use NORMA/SINGLE/ACCURATE"
            )
        if PREC[0].lower() not in {"a", "s", "n", "l", "m", "h"}:
            raise ValueError(f"{PREC=} does not exist. If this is no longer correct, please update this code.")

        CUTOFF = [
            np.sqrt(encut / _RYTOEV) / (2 * np.pi / (anorm / _AUTOA)) for anorm in self.poscar.structure.lattice.abc
        ]

        # TODO This only works in VASP 6.x
        _WFACT = 4 if PREC[0].lower() in {"a", "s"} else 3

        def next_g_size(cur_g_size):
            g_size = int(_WFACT * cur_g_size + 0.5)
            return next_num_with_prime_factors(g_size, max_prime_factor, must_inc_2)

        ng_vec = [*map(next_g_size, CUTOFF)]

        # TODO This works for VASP 5.x and 6.x
        finer_g_scale = 2 if PREC[0].lower() in {"a", "n"} else 1

        return ng_vec, [ng_ * finer_g_scale for ng_ in ng_vec]


# Helper functions to determine valid FFT grids for VASP
def next_num_with_prime_factors(n: int, max_prime_factor: int, must_inc_2: bool = True) -> int:
    """
    Return the next number greater than or equal to n that only has the desired prime factors.

    Args:
        n (int): Initial guess at the grid density
        max_prime_factor (int): the maximum prime factor
        must_inc_2 (bool): 2 must be a prime factor of the result

    Returns:
        int: first product of the prime_factors that is >= n
    """
    if max_prime_factor < 2:
        raise ValueError("Must choose a maximum prime factor greater than 2")
    prime_factors = primes_less_than(max_prime_factor)
    for new_val in itertools.count(start=n):
        if must_inc_2 and new_val % 2 != 0:
            continue
        cur_val_ = new_val
        for j in prime_factors:
            while cur_val_ % j == 0:
                cur_val_ //= j
        if cur_val_ == 1:
            return new_val
    raise ValueError("No factorable number found, not possible.")


def primes_less_than(max_val: int) -> list[int]:
    """Get the primes less than or equal to the max value."""
    res = []
    for i in range(2, max_val + 1):
        for j in range(2, i):
            if i % j == 0:
                break
        else:
            res.append(i)
    return res


@due.dcite(
    Doi("10.1016/j.commatsci.2011.02.023"),
    description="A high-throughput infrastructure for density functional theory calculations",
)
@dataclass
class MITRelaxSet(DictSet):
    """
    Standard implementation of VaspInputSet utilizing parameters in the MIT
    High-throughput project.
    The parameters are chosen specifically for a high-throughput project,
    which means in general pseudopotentials with fewer electrons were chosen.

    Args:
        structure (Structure): The Structure to create inputs for. If None, the input
            set is initialized without a Structure but one must be set separately before
            the inputs are generated.
        **kwargs: Same as those supported by DictSet.

    Please refer::

        A Jain, G. Hautier, C. Moore, S. P. Ong, C. Fischer, T. Mueller,
        K. A. Persson, G. Ceder. A high-throughput infrastructure for density
        functional theory calculations. Computational Materials Science,
        2011, 50(8), 2295-2310. doi:10.1016/j.commatsci.2011.02.023
    """

    CONFIG = _load_yaml_config("MITRelaxSet")


@dataclass
class MPRelaxSet(DictSet):
    """
    Implementation of VaspInputSet utilizing parameters in the public
    Materials Project. Typically, the pseudopotentials chosen contain more
    electrons than the MIT parameters, and the k-point grid is ~50% more dense.
    The LDAUU parameters are also different due to the different PSPs used,
    which result in different fitted values.

    Args:
        structure (Structure): The Structure to create inputs for. If None, the input
            set is initialized without a Structure but one must be set separately before
            the inputs are generated.
        **kwargs: Same as those supported by DictSet.
    """

    CONFIG = _load_yaml_config("MPRelaxSet")


@due.dcite(
    Doi("10.1021/acs.jpclett.0c02405"),
    description="AccurAccurate and Numerically Efficient r2SCAN Meta-Generalized Gradient Approximation",
)
@due.dcite(
    Doi("10.1103/PhysRevLett.115.036402"),
    description="Strongly Constrained and Appropriately Normed Semilocal Density Functional",
)
@due.dcite(
    Doi("10.1103/PhysRevB.93.155109"),
    description="Efficient generation of generalized Monkhorst-Pack grids through the use of informatics",
)
@dataclass
class MPScanRelaxSet(DictSet):
    """
    Class for writing a relaxation input set using the accurate and numerically
    efficient r2SCAN variant of the Strongly Constrained and Appropriately Normed
    (SCAN) metaGGA density functional.

    Notes:
        1. This functional is officially supported in VASP 6.0.0 and above. On older version,
        source code may be obtained by contacting the authors of the referenced manuscript.
        The original SCAN functional, available from VASP 5.4.3 onwards, maybe used instead
        by passing `user_incar_settings={"METAGGA": "SCAN"}` when instantiating this InputSet.
        r2SCAN and SCAN are expected to yield very similar results.

        2. Meta-GGA calculations require POTCAR files that include
        information on the kinetic energy density of the core-electrons,
        i.e. "PBE_52" or "PBE_54". Make sure the POTCARs include the
        following lines (see VASP wiki for more details):

            $ grep kinetic POTCAR
            kinetic energy-density
            mkinetic energy-density pseudized
            kinetic energy density (partial)

    Args:
        bandgap (float): Bandgap of the structure in eV. The bandgap is used to
            compute the appropriate k-point density and determine the
            smearing settings.
            Metallic systems (default, bandgap = 0) use a KSPACING value of 0.22
            and Methfessel-Paxton order 2 smearing (ISMEAR=2, SIGMA=0.2).
            Non-metallic systems (bandgap > 0) use the tetrahedron smearing
            method (ISMEAR=-5, SIGMA=0.05). The KSPACING value is
            calculated from the bandgap via Eqs. 25 and 29 of Wisesa, McGill,
            and Mueller [1] (see References). Note that if 'user_incar_settings'
            or 'user_kpoints_settings' override KSPACING, the calculation from
            bandgap is not performed.
        vdw (str): set "rVV10" to enable SCAN+rVV10, which is a versatile
            van der Waals density functional by combing the SCAN functional
            with the rVV10 non-local correlation functional. rvv10 is the only
            dispersion correction available for SCAN at this time.
        **kwargs: Same as those supported by DictSet.

        References:
            [1] P. Wisesa, K.A. McGill, T. Mueller, Efficient generation of
            generalized Monkhorst-Pack grids through the use of informatics,
            Phys. Rev. B. 93 (2016) 1-10. doi:10.1103/PhysRevB.93.155109.

    References:
        James W. Furness, Aaron D. Kaplan, Jinliang Ning, John P. Perdew, and Jianwei Sun.
        Accurate and Numerically Efficient r2SCAN Meta-Generalized Gradient Approximation.
        The Journal of Physical Chemistry Letters 0, 11 DOI: 10.1021/acs.jpclett.0c02405
    """

    bandgap: float | None = None
    user_potcar_functional: UserPotcarFunctional = "PBE_54"
    auto_ismear: bool = True
    CONFIG = _load_yaml_config("MPSCANRelaxSet")
    _valid_potcars = ("PBE_52", "PBE_54")

    def __post_init__(self):
        super().__post_init__()
        if self.vdw and self.vdw != "rvv10":
            warnings.warn("Use of van der waals functionals other than rVV10 with SCAN is not supported at this time. ")
            # delete any vdw parameters that may have been added to the INCAR
            vdw_par = loadfn(str(MODULE_DIR / "vdW_parameters.yaml"))
            for k in vdw_par[self.vdw]:
                self._config_dict["INCAR"].pop(k, None)

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        return {"KSPACING": "auto"}  # enable automatic KSPACING based on band gap


@dataclass
class MPMetalRelaxSet(DictSet):
    """
    Implementation of VaspInputSet utilizing parameters in the public
    Materials Project, but with tuning for metals. Key things are a denser
    k point density, and a.
    """

    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        return {"ISMEAR": 1, "SIGMA": 0.2}

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Get updates to the kpoints configuration for this calculation type."""
        return {"reciprocal_density": 200}


@dataclass
class MPHSERelaxSet(DictSet):
    """Same as the MPRelaxSet, but with HSE parameters."""

    CONFIG = _load_yaml_config("MPHSERelaxSet")


@dataclass
class MPStaticSet(DictSet):
    """Creates input files for a static calculation.

    Args:
        structure (Structure): Structure from previous run.
        lepsilon (bool): Whether to add static dielectric calculation
        lcalcpol (bool): Whether to turn on evaluation of the Berry phase approximations
            for electronic polarization
        reciprocal_density (int): For static calculations, we usually set the
            reciprocal density by volume. This is a convenience arg to change
            that, rather than using user_kpoints_settings. Defaults to 100,
            which is ~50% more than that of standard relaxation calculations.
        small_gap_multiply ([float, float]): If the gap is less than
            1st index, multiply the default reciprocal_density by the 2nd
            index.
        **kwargs: kwargs supported by MPRelaxSet.
    """

    lepsilon: bool = False
    lcalcpol: bool = False
    reciprocal_density: int = 100
    small_gap_multiply: tuple[float, float] | None = None
    inherit_incar: bool = True
    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: dict[str, Any] = {"NSW": 0, "ISMEAR": -5, "LCHARG": True, "LORBIT": 11, "LREAL": False}
        if self.lepsilon:
            # LPEAD=T: numerical evaluation of overlap integral prevents LRF_COMMUTATOR
            # errors and can lead to better expt. agreement but produces slightly
            # different results
            updates.update({"IBRION": 8, "LEPSILON": True, "LPEAD": True, "NSW": 1, "EDIFF": 1e-5})

        if self.lcalcpol:
            updates["LCALCPOL"] = True
        return updates

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Get updates to the kpoints configuration for this calculation type."""
        factor = 1.0
        if self.bandgap is not None and self.small_gap_multiply and self.bandgap <= self.small_gap_multiply[0]:
            factor = self.small_gap_multiply[1]

        # prefer to use k-point scheme from previous run unless lepsilon = True is specified
        if self.prev_kpoints and self.prev_kpoints.style == Kpoints.supported_modes.Monkhorst and not self.lepsilon:  # type: ignore
            kpoints = Kpoints.automatic_density_by_vol(
                self.structure,  # type: ignore
                int(self.reciprocal_density * factor),
                self.force_gamma,
            )
            k_div = [kp + 1 if kp % 2 == 1 else kp for kp in kpoints.kpts[0]]  # type: ignore
            return Kpoints.monkhorst_automatic(k_div)  # type: ignore

        return {"reciprocal_density": self.reciprocal_density * factor}


@dataclass
class MatPESStaticSet(DictSet):
    """Creates input files for a MatPES static calculation.

    The goal of MatPES is to generate potential energy surface data. This is a distinctly different
    from the objectives of the MP static calculations, which aims to obtain primarily accurate
    energies and also electronic structure (DOS). For PES data, force accuracy (and to some extent,
    stress accuracy) is of paramount importance.

    The default POTCAR versions have been updated to PBE_54 from the old PBE set used in the
    MPStaticSet. However, **U values** are still based on PBE. The implicit assumption here is that
    the PBE_54 and PBE POTCARs are sufficiently similar that the U values fitted to the old PBE
    functional still applies.

    Args:
        structure (Structure): The Structure to create inputs for. If None, the input
            set is initialized without a Structure but one must be set separately before
            the inputs are generated.
        xc_functional ('R2SCAN'|'PBE'): Exchange-correlation functional to use. Defaults to 'PBE'.
        **kwargs: Same as those supported by DictSet.
    """

    xc_functional: Literal["R2SCAN", "PBE", "PBE+U"] = "PBE"
    prev_incar: dict | str | None = None
    # These are parameters that we will inherit from any previous INCAR supplied. They are mostly parameters related
    # to symmetry and convergence set by Custodian when errors are encountered in a previous run. Given that our goal
    # is to have a strictly homogeneous PES data, all other parameters (e.g., ISMEAR, ALGO, etc.) are not inherited.
    inherit_incar: list[str] | bool = (  # type: ignore
        "LPEAD",
        "NGX",
        "NGY",
        "NGZ",
        "SYMPREC",
        "IMIX",
        "LMAXMIX",
        "KGAMMA",
        "ISYM",
        "NCORE",
        "NPAR",
        "NELMIN",
        "IOPT",
        "NBANDS",
        "KPAR",
        "AMIN",
        "NELMDL",
        "BMIX",
        "AMIX_MAG",
        "BMIX_MAG",
    )
    CONFIG = _load_yaml_config("MatPESStaticSet")

    def __post_init__(self):
        """Validate inputs"""
        super().__post_init__()
        valid_xc_functionals = ("R2SCAN", "PBE", "PBE+U")
        if self.xc_functional.upper() not in valid_xc_functionals:
            raise ValueError(
                f"Unrecognized xc_functional='{self.xc_functional}'. "
                f"Supported exchange-correlation functionals are {valid_xc_functionals}"
            )

        default_potcars = self.CONFIG["PARENT"].replace("PBE", "PBE_").replace("BASE", "")  # PBE64BASE -> PBE_64
        self.user_potcar_functional = self.user_potcar_functional or default_potcars
        if self.user_potcar_functional.upper() != default_potcars:
            warnings.warn(
                f"{self.user_potcar_functional=} is inconsistent with the recommended {default_potcars}.", UserWarning
            )

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: dict[str, Any] = {}
        if self.xc_functional.upper() == "R2SCAN":
            updates.update({"METAGGA": "R2SCAN", "ALGO": "ALL", "GGA": None})
        if self.xc_functional.upper().endswith("+U"):
            updates["LDAU"] = True
        return updates


@dataclass
class MPScanStaticSet(MPScanRelaxSet):
    """
    Creates input files for a static calculation using the accurate and numerically
    efficient r2SCAN variant of the Strongly Constrained and Appropriately Normed
    (SCAN) metaGGA functional.

    Args:
        structure (Structure): Structure from previous run.
        bandgap (float): Bandgap of the structure in eV. The bandgap is used to
            compute the appropriate k-point density and determine the smearing settings.
        lepsilon (bool): Whether to add static dielectric calculation
        lcalcpol (bool): Whether to turn on evaluation of the Berry phase approximations
            for electronic polarization.
        **kwargs: kwargs supported by MPScanRelaxSet.
    """

    lepsilon: bool = False
    lcalcpol: bool = False
    inherit_incar: bool = True

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: dict[str, Any] = {
            "LREAL": False,
            "NSW": 0,
            "LORBIT": 11,
            "LVHAR": True,
            "ISMEAR": -5,
            "KSPACING": "auto",
        }

        if self.lepsilon:
            # LPEAD=T: numerical evaluation of overlap integral prevents
            # LRF_COMMUTATOR errors and can lead to better expt. agreement
            # but produces slightly different results
            updates.update({"IBRION": 8, "LEPSILON": True, "LPEAD": True, "NSW": 1, "NPAR": None})

        if self.lcalcpol:
            updates["LCALCPOL"] = True

        return updates


@dataclass
class MPHSEBSSet(DictSet):
    """
    Implementation of a VaspInputSet for HSE band structure computations.

    Remember that HSE band structures must be self-consistent in VASP. A band structure
    along symmetry lines for instance needs BOTH a uniform grid with appropriate weights
    AND a path along the lines with weight 0.

    Thus, the "uniform" mode is just like regular static SCF but allows adding custom
    kpoints (e.g., corresponding to known VBM/CBM) to the uniform grid that have zero
    weight (e.g., for better gap estimate).

    The "gap" mode behaves just like the "uniform" mode, however, if starting from a
    previous calculation, the VBM and CBM k-points will automatically be added to
    ``added_kpoints``.

    The "line" mode is just like "uniform" mode, but additionally adds k-points along
    symmetry lines with zero weight.

    The "uniform_dense" mode is like "uniform" mode but additionally adds a denser
    uniform mesh with zero weight. This can be useful when calculating Fermi surfaces
    or BoltzTraP/AMSET electronic transport using hybrid DFT.

    Args:
        structure (Structure): Structure to compute
        added_kpoints (list): a list of kpoints (list of 3 number list) added to the
            run. The k-points are in fractional coordinates
        mode (str): "Line" - generate k-points along symmetry lines for bandstructure.
            "Uniform" - generate uniform k-points grid.
        reciprocal_density (int): k-point density to use for uniform mesh.
        copy_chgcar (bool): Whether to copy the CHGCAR of a previous run.
        kpoints_line_density (int): k-point density for high symmetry lines
        dedos (float): Energy difference used to set NEDOS, based on the total energy
            range.
        optics (bool): Whether to add LOPTICS (used for calculating optical response).
        nbands_factor (float): Multiplicative factor for NBANDS when starting from a
            previous calculation. Choose a higher number if you are doing an LOPTICS
            calculation.
        **kwargs (dict): Any other parameters to pass into DictSet.
    """

    added_kpoints: list[Vector3D] = field(default_factory=list)
    mode: str = "gap"
    reciprocal_density: float = 50
    copy_chgcar: bool = True
    kpoints_line_density: float = 20
    nbands_factor: float = 1.2
    zero_weighted_reciprocal_density: float = 100
    dedos: float = 0.02
    optics: bool = False
    CONFIG = MPHSERelaxSet.CONFIG

    def __post_init__(self) -> None:
        """Ensure mode is set correctly."""
        super().__post_init__()

        if "reciprocal_density" in self.user_kpoints_settings:
            self.reciprocal_density = self.user_kpoints_settings["reciprocal_density"]

        self.mode = self.mode.lower()
        supported_modes = ("line", "uniform", "gap", "uniform_dense")
        if self.mode not in supported_modes:
            raise ValueError(f"Supported modes are: {', '.join(supported_modes)}")

    @property
    def kpoints_updates(self) -> dict:
        """Get updates to the kpoints configuration for this calculation type."""
        kpoints: dict[str, Any] = {"reciprocal_density": self.reciprocal_density, "explicit": True}

        if self.mode == "line":
            # add line_density on top of reciprocal density
            kpoints["zero_weighted_line_density"] = self.kpoints_line_density

        elif self.mode == "uniform_dense":
            kpoints["zero_weighted_reciprocal_density"] = self.zero_weighted_reciprocal_density

        added_kpoints = deepcopy(self.added_kpoints)
        if self.prev_vasprun is not None and self.mode == "gap":
            bs = self.prev_vasprun.get_band_structure()
            if not bs.is_metal():
                added_kpoints.append(bs.get_vbm()["kpoint"].frac_coords)
                added_kpoints.append(bs.get_cbm()["kpoint"].frac_coords)

        kpoints["added_kpoints"] = added_kpoints

        return kpoints

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates = dict(NSW=0, ISMEAR=0, SIGMA=0.05, ISYM=3, LCHARG=False, NELMIN=5)

        if self.mode == "uniform" and len(self.added_kpoints) == 0:
            # automatic setting of nedos using the energy range and the energy step
            nedos = _get_nedos(self.prev_vasprun, self.dedos)

            # use tetrahedron method for DOS and optics calculations
            updates.update({"ISMEAR": -5, "NEDOS": nedos})

        else:
            # if line mode or explicit k-points (gap) can't use ISMEAR=-5
            # use small sigma to avoid partial occupancies for small band gap materials
            updates.update({"ISMEAR": 0, "SIGMA": 0.01})

        if self.prev_vasprun is not None:
            # set nbands
            nbands = int(np.ceil(self.prev_vasprun.parameters["NBANDS"] * self.nbands_factor))
            updates["NBANDS"] = nbands

        if self.optics:
            # LREAL not supported with LOPTICS
            updates.update({"LOPTICS": True, "LREAL": False, "CSHIFT": 1e-5})

        if self.prev_vasprun is not None and self.prev_outcar is not None:
            # turn off spin when magmom for every site is smaller than 0.02.
            updates["ISPIN"] = _get_ispin(self.prev_vasprun, self.prev_outcar)

        return updates


@dataclass
class MPNonSCFSet(DictSet):
    """
    Init a MPNonSCFSet. Typically, you would use the classmethod
    from_prev_calc to initialize from a previous SCF run.

    Args:
        structure (Structure): Structure to compute
        mode (str): Line, Uniform or Boltztrap mode supported.
        nedos (int): nedos parameter. Default to 2001.
        dedos (float): setting nedos=0 and uniform mode in from_prev_calc,
            an automatic nedos will be calculated using the total energy range
            divided by the energy step dedos
        reciprocal_density (int): density of k-mesh by reciprocal
            volume (defaults to 100)
        kpoints_line_density (int): Line density for Line mode.
        optics (bool): whether to add dielectric function
        copy_chgcar: Whether to copy the old CHGCAR when starting from a
            previous calculation.
        nbands_factor (float): Multiplicative factor for NBANDS when starting
            from a previous calculation. Choose a higher number if you are
            doing an LOPTICS calculation.
        small_gap_multiply ([float, float]): When starting from a previous
            calculation, if the gap is less than 1st index, multiply the default
            reciprocal_density by the 2nd index.
        **kwargs: kwargs supported by MPRelaxSet.
    """

    mode: str = "line"
    nedos: int = 2001
    dedos: float = 0.005
    reciprocal_density: float = 100
    kpoints_line_density: float = 20
    optics: bool = False
    copy_chgcar: bool = True
    nbands_factor: float = 1.2
    small_gap_multiply: tuple[float, float] | None = None
    inherit_incar: bool = True
    CONFIG = MPRelaxSet.CONFIG

    def __post_init__(self):
        """Perform inputset validation."""
        super().__post_init__()

        self.mode = self.mode.lower()

        if self.mode not in ["line", "uniform", "boltztrap"]:
            raise ValueError("Supported modes for NonSCF runs are 'line', 'uniform' and 'boltztrap!")

        if (self.mode.lower() != "uniform" or self.nedos < 2000) and self.optics:
            warnings.warn("It is recommended to use Uniform mode with a high NEDOS for optics calculations.")

        if self.standardize:
            warnings.warn(
                "Use of standardize=True with from_prev_run is not "
                "recommended as there is no guarantee the copied "
                "files will be appropriate for the standardized"
                " structure. copy_chgcar is enforced to be false."
            )
            self.copy_chgcar = False

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: dict[str, Any] = {"LCHARG": False, "LORBIT": 11, "LWAVE": False, "NSW": 0, "ISYM": 0, "ICHARG": 11}

        if self.prev_vasprun is not None:
            # set NBANDS
            n_bands = int(np.ceil(self.prev_vasprun.parameters["NBANDS"] * self.nbands_factor))
            updates["NBANDS"] = n_bands

        # automatic setting of NEDOS using the energy range and the energy step
        nedos = _get_nedos(self.prev_vasprun, self.dedos) if self.nedos == 0 else self.nedos

        if self.mode == "uniform":
            # use tetrahedron method for DOS and optics calculations
            updates.update({"ISMEAR": -5, "ISYM": 2, "NEDOS": nedos})

        elif self.mode in ("line", "boltztrap"):
            # if line mode or explicit k-points (boltztrap) can't use ISMEAR=-5
            # use small sigma to avoid partial occupancies for small band gap materials
            # use a larger sigma if the material is a metal
            sigma = 0.2 if self.bandgap == 0 or self.bandgap is None else 0.01
            updates.update({"ISMEAR": 0, "SIGMA": sigma})

        if self.optics:
            # LREAL not supported with LOPTICS = True; automatic NEDOS usually
            # underestimates, so set it explicitly
            updates.update({"LOPTICS": True, "LREAL": False, "CSHIFT": 1e-5, "NEDOS": nedos})

        if self.prev_vasprun is not None and self.prev_outcar is not None:
            # turn off spin when magmom for every site is smaller than 0.02.
            updates["ISPIN"] = _get_ispin(self.prev_vasprun, self.prev_outcar)

        updates["MAGMOM"] = None
        return updates

    @property
    def kpoints_updates(self) -> dict:
        """Get updates to the kpoints configuration for this calculation type."""
        factor = 1.0
        if self.bandgap is not None and self.small_gap_multiply and self.bandgap <= self.small_gap_multiply[0]:
            factor = self.small_gap_multiply[1]

        if self.mode == "line":
            return {"line_density": self.kpoints_line_density * factor}

        if self.mode == "boltztrap":
            return {"explicit": True, "reciprocal_density": self.reciprocal_density * factor}

        return {"reciprocal_density": self.reciprocal_density * factor}


@dataclass
class MPSOCSet(DictSet):
    """An input set for running spin-orbit coupling (SOC) calculations.

    Args:
        structure (Structure): the structure must have the 'magmom' site
            property and each magnetic moment value must have 3
            components. eg: ``magmom = [[0,0,2], ...]``
        saxis (tuple): magnetic moment orientation
        copy_chgcar: Whether to copy the old CHGCAR. Defaults to True.
        nbands_factor (float): Multiplicative factor for NBANDS. Choose a
            higher number if you are doing an LOPTICS calculation.
        reciprocal_density (int): density of k-mesh by reciprocal volume.
        small_gap_multiply ([float, float]): If the gap is less than
            1st index, multiply the default reciprocal_density by the 2nd
            index.
        lepsilon (bool): Whether to add static dielectric calculation
        lcalcpol (bool): Whether to turn on evaluation of the Berry phase approximations
            for electronic polarization
        magmom (list[list[float]]): Override for the structure magmoms.
        **kwargs: kwargs supported by DictSet.
    """

    saxis: tuple[int, int, int] = (0, 0, 1)
    nbands_factor: float = 1.2
    lepsilon: bool = False
    lcalcpol: bool = False
    reciprocal_density: float = 100
    small_gap_multiply: tuple[float, float] | None = None
    magmom: list[Vector3D] | None = None
    inherit_incar: bool = True
    copy_chgcar: bool = True
    CONFIG = MPRelaxSet.CONFIG

    def __post_init__(self):
        super().__post_init__()
        if (
            self.structure
            and not hasattr(self.structure[0], "magmom")
            and not isinstance(self.structure[0].magmom, list)
        ):
            raise ValueError(
                "The structure must have the 'magmom' site "
                "property and each magnetic moment value must have 3 "
                "components. eg:- magmom = [0,0,2]"
            )

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates = {
            "ISYM": -1,
            "LSORBIT": "T",
            "ICHARG": 11,
            "SAXIS": list(self.saxis),
            "NSW": 0,
            "ISMEAR": -5,
            "LCHARG": True,
            "LORBIT": 11,
            "LREAL": False,
        }

        if self.lepsilon:
            # LPEAD=T: numerical evaluation of overlap integral prevents LRF_COMMUTATOR
            # errors and can lead to better expt. agreement but produces slightly
            # different results
            updates.update({"IBRION": 8, "LEPSILON": True, "LPEAD": True, "NSW": 1})

        if self.lcalcpol:
            updates["LCALCPOL"] = True

        if self.prev_vasprun is not None:
            # set NBANDS
            n_bands = int(np.ceil(self.prev_vasprun.parameters["NBANDS"] * self.nbands_factor))
            updates["NBANDS"] = n_bands
        return updates

    @property
    def kpoints_updates(self) -> dict:
        """Get updates to the kpoints configuration for this calculation type."""
        factor = 1.0
        if self.bandgap is not None and self.small_gap_multiply and self.bandgap <= self.small_gap_multiply[0]:
            factor = self.small_gap_multiply[1]
        return {"reciprocal_density": self.reciprocal_density * factor}

    @DictSet.structure.setter  # type: ignore
    def structure(self, structure: Structure | None) -> None:
        if structure is not None:
            if self.magmom:
                structure = structure.copy(site_properties={"magmom": self.magmom})

            # magmom has to be 3D for SOC calculation.
            if hasattr(structure[0], "magmom"):
                if not isinstance(structure[0].magmom, list):
                    structure = structure.copy(site_properties={"magmom": [[0, 0, site.magmom] for site in structure]})
            else:
                raise ValueError("Neither the previous structure has magmom property nor magmom provided")

        DictSet.structure.fset(self, structure)  # type: ignore


@dataclass
class MPNMRSet(DictSet):
    """Init a MPNMRSet.

    Args:
        structure (Structure): Structure from previous run.
        mode (str): The NMR calculation to run
            "cs": for Chemical Shift
            "efg" for Electric Field Gradient
        isotopes (list): list of Isotopes for quadrupole moments
        reciprocal_density (int): density of k-mesh by reciprocal volume.
        lepsilon (bool): Whether to add static dielectric calculation
        lcalcpol (bool): Whether to turn on evaluation of the Berry phase approximations
            for electronic polarization
        reciprocal_density (int): For static calculations, we usually set the
            reciprocal density by volume. This is a convenience arg to change
            that, rather than using user_kpoints_settings. Defaults to 100,
            which is ~50% more than that of standard relaxation calculations.
        small_gap_multiply ([float, float]): If the gap is less than
            1st index, multiply the default reciprocal_density by the 2nd
            index.
        **kwargs: kwargs supported by MPRelaxSet.
    """

    mode: Literal["cs", "efg"] = "cs"
    isotopes: list = field(default_factory=list)
    reciprocal_density: int = 100
    small_gap_multiply: tuple[float, float] | None = None
    inherit_incar: bool = True
    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates: dict[str, Any] = {"NSW": 0, "ISMEAR": -5, "LCHARG": True, "LORBIT": 11, "LREAL": False}
        if self.mode.lower() == "cs":
            updates.update(
                LCHIMAG=True,
                EDIFF=-1.0e-10,
                ISYM=0,
                LCHARG=False,
                LNMR_SYM_RED=True,
                NELMIN=10,
                NLSPLINE=True,
                PREC="ACCURATE",
                SIGMA=0.01,
            )
        elif self.mode.lower() == "efg":
            isotopes = {ist.split("-")[0]: ist for ist in self.isotopes}
            quad_efg = [
                float(Species(s.name).get_nmr_quadrupole_moment(isotopes.get(s.name)))
                for s in self.structure.species  # type: ignore
            ]
            updates.update(
                ALGO="FAST",
                EDIFF=-1.0e-10,
                ISYM=0,
                LCHARG=False,
                LEFG=True,
                QUAD_EFG=quad_efg,
                NELMIN=10,
                PREC="ACCURATE",
                SIGMA=0.01,
            )
        return updates

    @property
    def kpoints_updates(self) -> dict:
        """Get updates to the kpoints configuration for this calculation type."""
        factor = 1.0
        if self.bandgap is not None and self.small_gap_multiply and self.bandgap <= self.small_gap_multiply[0]:
            factor = self.small_gap_multiply[1]
        return {"reciprocal_density": self.reciprocal_density * factor}


@due.dcite(
    Doi("10.1149/2.0061602jes"),
    description="Elastic Properties of Alkali Superionic Conductor Electrolytes from First Principles Calculations",
)
class MVLElasticSet(DictSet):
    """
    MVL denotes VASP input sets that are implemented by the Materials Virtual
    Lab (http://materialsvirtuallab.org) for various research.

    This input set is used to calculate elastic constants in VASP. It is used
    in the following work::

        Z. Deng, Z. Wang, I.-H. Chu, J. Luo, S. P. Ong.
        “Elastic Properties of Alkali Superionic Conductor Electrolytes
        from First Principles Calculations”, J. Electrochem. Soc.
        2016, 163(2), A67-A74. doi: 10.1149/2.0061602jes

    To read the elastic constants, you may use the Outcar class which parses the
    elastic constants.

    Args:
        structure (pymatgen.Structure): Input structure.
        potim (float): POTIM parameter. The default of 0.015 is usually fine,
            but some structures may require a smaller step.
        **kwargs: Parameters supported by MPRelaxSet.
    """

    potim: float = 0.015
    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        return {"IBRION": 6, "NFREE": 2, "POTIM": self.potim, "NPAR": None}


@dataclass
class MVLGWSet(DictSet):
    """
    MVL denotes VASP input sets that are implemented by the Materials Virtual
    Lab (http://materialsvirtuallab.org) for various research. This is a
    flexible input set for GW calculations.

    Note that unlike all other input sets in this module, the PBE_54 series of
    functional is set as the default. These have much improved performance for
    GW calculations.

    A typical sequence is mode="STATIC" -> mode="DIAG" -> mode="GW" ->
    mode="BSE". For all steps other than the first one (static), the
    recommendation is to use from_prev_calculation on the preceding run in
    the series.

    Args:
        structure (Structure): Input structure.
        mode (str): Supported modes are "STATIC" (default), "DIAG", "GW",
            and "BSE".
        nbands (int): For subsequent calculations, it is generally
            recommended to perform NBANDS convergence starting from the
            NBANDS of the previous run for DIAG, and to use the exact same
            NBANDS for GW and BSE. This parameter is used by
            from_previous_calculation to set nband.
        copy_wavecar: Whether to copy the old WAVECAR, WAVEDER and associated
            files when starting from a previous calculation.
        nbands_factor (int): Multiplicative factor for NBANDS when starting
            from a previous calculation. Only applies if mode=="DIAG".
            Need to be tested for convergence.
        reciprocal_density (int): Density of k-mesh by reciprocal atom. Only
            applies if mode=="STATIC". Defaults to 100.
        ncores (int): Numbers of cores used for the calculation. VASP will alter
            NBANDS if it was not dividable by ncores. Only applies if
            mode=="DIAG".
        **kwargs: All kwargs supported by DictSet. Typically,
            user_incar_settings is a commonly used option.
    """

    reciprocal_density: float = 100
    mode: str = "STATIC"
    copy_wavecar: bool = True
    nbands_factor: int = 5
    ncores: int = 16
    nbands: int | None = None
    force_gamma: bool = True
    inherit_incar: bool = True  # inherit incar from previous run if available
    SUPPORTED_MODES = ("DIAG", "GW", "STATIC", "BSE")
    CONFIG = _load_yaml_config("MVLGWSet")

    def __post_init__(self):
        """Validate input settings."""
        super().__post_init__()
        self.mode = self.mode.upper()

        if self.mode not in MVLGWSet.SUPPORTED_MODES:
            raise ValueError(f"{self.mode} not one of the support modes : {MVLGWSet.SUPPORTED_MODES}")

    @property
    def kpoints_updates(self) -> dict:
        """Get updates to the kpoints configuration for this calculation type."""
        # Generate gamma center k-points mesh grid for GW calc, which is requested
        # by GW calculation.
        return {"reciprocal_density": self.reciprocal_density}

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates = {}
        nbands = int(self.prev_vasprun.parameters["NBANDS"]) if self.prev_vasprun is not None else None

        if self.mode == "DIAG":
            # Default parameters for diagonalization calculation.
            updates.update({"ALGO": "Exact", "NELM": 1, "LOPTICS": True, "LPEAD": True})
            if nbands:
                nbands = int(np.ceil(nbands * self.nbands_factor / self.ncores) * self.ncores)

        elif self.mode == "GW":
            # Default parameters for GW calculation.
            updates.update(
                {"ALGO": "GW0", "NELM": 1, "NOMEGA": 80, "ENCUTGW": 250, "EDIFF": None, "LOPTICS": None, "LPEAD": None}
            )
        elif self.mode == "BSE":
            # Default parameters for BSE calculation.
            updates.update({"ALGO": "BSE", "ANTIRES": 0, "NBANDSO": 20, "NBANDSV": 20})

        if nbands:
            updates["NBANDS"] = nbands

        return updates

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, mode="DIAG", **kwargs):
        """
        Generate a set of VASP input files for GW or BSE calculations from a
        directory of previous Exact Diag VASP run.

        Args:
            prev_calc_dir (str): The directory contains the outputs(
                vasprun.xml of previous vasp run.
            mode (str): Supported modes are "STATIC", "DIAG" (default), "GW",
                and "BSE".
            **kwargs: All kwargs supported by MVLGWSet, other than structure,
                prev_incar and mode, which are determined from the
                prev_calc_dir.
        """
        input_set = cls(_dummy_structure, mode=mode, **kwargs)
        return input_set.override_from_prev_calc(prev_calc_dir=prev_calc_dir)


@dataclass
class MVLSlabSet(DictSet):
    """
    Class for writing a set of slab vasp runs,
    including both slabs (along the c direction) and orient unit cells (bulk),
    to ensure the same KPOINTS, POTCAR and INCAR criterion.

    Args:
        structure: Structure
        k_product: default to 50, kpoint number * length for a & b
            directions, also for c direction in bulk calculations
        bulk:
        auto_dipole:
        set_mix:
        sort_structure:
        **kwargs: Other kwargs supported by DictSet.
    """

    k_product: int = 50
    bulk: bool = False
    auto_dipole: bool = False
    set_mix: bool = True
    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates = {"EDIFF": 1e-4, "EDIFFG": -0.02, "ENCUT": 400, "ISMEAR": 0, "SIGMA": 0.05, "ISIF": 3}
        if not self.bulk:
            updates.update({"ISIF": 2, "LVTOT": True, "NELMIN": 8})
            if self.set_mix:
                updates.update({"AMIN": 0.01, "AMIX": 0.2, "BMIX": 0.001})
            if self.auto_dipole:
                weights = [s.species.weight for s in self.structure]  # type: ignore
                center_of_mass = np.average(self.structure.frac_coords, weights=weights, axis=0)  # type: ignore
                updates.update({"IDIPOL": 3, "LDIPOL": True, "DIPOL": center_of_mass})
        return updates

    @property
    def kpoints_updates(self):
        """Get updates to the kpoints configuration for this calculation type.

        k_product, default to 50, is kpoint number * length for a & b
            directions, also for c direction in bulk calculations
        Automatic mesh & Gamma is the default setting.
        """
        # To get input sets, the input structure has to has the same number
        # of required parameters as a Structure object (ie. 4). Slab
        # attributes aren't going to affect the VASP inputs anyways so
        # converting the slab into a structure should not matter
        # use k_product to calculate kpoints, k_product = kpts[0][0] * a
        lattice_abc = self.structure.lattice.abc
        kpt_calc = [
            int(self.k_product / lattice_abc[0] + 0.5),
            int(self.k_product / lattice_abc[1] + 0.5),
            1,
        ]

        # calculate kpts (c direction) for bulk. (for slab, set to 1)
        if self.bulk:
            kpt_calc[2] = int(self.k_product / lattice_abc[2] + 0.5)

        return Kpoints(comment="Generated by pymatgen's MVLGBSet", style=Kpoints.supported_modes.Gamma, kpts=[kpt_calc])

    def as_dict(self, verbosity=2):
        """
        Args:
            verbosity (int): Verbosity of dict. E.g., whether to include Structure.

        Returns:
            MSONable dict
        """
        dct = MSONable.as_dict(self)
        if verbosity == 1:
            dct.pop("structure", None)
        return dct


@dataclass
class MVLGBSet(DictSet):
    """
    Class for writing a vasp input files for grain boundary calculations, slab or bulk.

    Args:
        structure(Structure): provide the structure
        k_product: Kpoint number * length for a & b directions, also for c direction in
            bulk calculations. Default to 40.
        slab_mode (bool): Defaults to False. Use default (False) for a bulk supercell.
            Use True if you are performing calculations on a slab-like (i.e., surface)
            of the GB, for example, when you are calculating the work of separation.
        is_metal (bool): Defaults to True. This determines whether an ISMEAR of 1 is
            used (for metals) or not (for insulators and semiconductors) by default.
            Note that it does *not* override user_incar_settings, which can be set by
            the user to be anything desired.
        **kwargs:
            Other kwargs supported by MPRelaxSet.
    """

    k_product: int = 40
    slab_mode: bool = False
    is_metal: bool = True
    CONFIG = MPRelaxSet.CONFIG

    @property
    def kpoints_updates(self):
        """
        k_product, default to 40, is kpoint number * length for a & b
        directions, also for c direction in bulk calculations
        Automatic mesh & Gamma is the default setting.
        """
        # use k_product to calculate kpoints, k_product = kpts[0][0] * a
        lengths = self.structure.lattice.abc
        kpt_calc = [
            int(self.k_product / lengths[0] + 0.5),
            int(self.k_product / lengths[1] + 0.5),
            int(self.k_product / lengths[2] + 0.5),
        ]

        if self.slab_mode:
            kpt_calc[2] = 1

        return Kpoints(comment="Generated by pymatgen's MVLGBSet", style=Kpoints.supported_modes.Gamma, kpts=[kpt_calc])

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        # The default incar setting is used for metallic system, for
        # insulator or semiconductor, ISMEAR need to be changed.
        updates = dict(LCHARG=False, NELM=60, PREC="Normal", EDIFFG=-0.02, ICHARG=0, NSW=200, EDIFF=0.0001)

        if self.is_metal:
            updates["ISMEAR"] = 1
            updates["LDAU"] = False

        if self.slab_mode:
            # for clean grain boundary and bulk relaxation, full optimization
            # relaxation (ISIF=3) is used. For slab relaxation (ISIF=2) is used.
            updates["ISIF"] = 2
            updates["NELMIN"] = 8

        return updates


@dataclass
class MVLRelax52Set(DictSet):
    """
    Implementation of VaspInputSet utilizing the public Materials Project
    parameters for INCAR & KPOINTS and VASP's recommended PAW potentials for
    POTCAR.

    Keynotes from VASP manual:
        1. Recommended potentials for calculations using vasp.5.2+
        2. If dimers with short bonds are present in the compound (O2, CO,
            N2, F2, P2, S2, Cl2), it is recommended to use the h potentials.
            Specifically, C_h, O_h, N_h, F_h, P_h, S_h, Cl_h
        3. Released on Oct 28, 2018 by VASP. Please refer to VASP
            Manual 1.2, 1.3 & 10.2.1 for more details.

    Args:
        structure (Structure): input structure.
        user_potcar_functional (str): choose from "PBE_52" and "PBE_54".
        **kwargs: Other kwargs supported by DictSet.
    """

    user_potcar_functional: UserPotcarFunctional = "PBE_52"
    CONFIG = _load_yaml_config("MVLRelax52Set")
    _valid_potcars = ("PBE_52", "PBE_54")


class MITNEBSet(DictSet):
    """
    Class for writing NEB inputs.

    Note that EDIFF is not on a per atom basis for this input set.
    """

    def __init__(self, structures, unset_encut=False, **kwargs):
        """
        Args:
            structures: List of Structure objects.
            unset_encut (bool): Whether to unset ENCUT.
            **kwargs: Other kwargs supported by DictSet.
        """
        if len(structures) < 3:
            raise ValueError("You need at least 3 structures for an NEB.")
        kwargs["sort_structure"] = False
        super().__init__(structures[0], MITRelaxSet.CONFIG, **kwargs)
        self.structures = self._process_structures(structures)

        self.unset_encut = False
        if unset_encut:
            self._config_dict["INCAR"].pop("ENCUT", None)

        if "EDIFF" not in self._config_dict["INCAR"]:
            self._config_dict["INCAR"]["EDIFF"] = self._config_dict["INCAR"].pop("EDIFF_PER_ATOM")

        # NEB specific defaults
        defaults = {"IMAGES": len(structures) - 2, "IBRION": 1, "ISYM": 0, "LCHARG": False, "LDAU": False}
        self._config_dict["INCAR"].update(defaults)

    @property
    def poscar(self):
        """Poscar for structure of first end point."""
        return Poscar(self.structures[0])

    @property
    def poscars(self):
        """List of Poscars."""
        return [Poscar(s) for s in self.structures]

    @staticmethod
    def _process_structures(structures):
        """Remove any atom jumps across the cell."""
        input_structures = structures
        structures = [input_structures[0]]
        for s in input_structures[1:]:
            prev = structures[-1]
            for i, site in enumerate(s):
                t = np.round(prev[i].frac_coords - site.frac_coords)
                if np.any(np.abs(t) > 0.5):
                    s.translate_sites([i], t, to_unit_cell=False)
            structures.append(s)
        return structures

    def write_input(
        self,
        output_dir,
        make_dir_if_not_present=True,
        write_cif=False,
        write_path_cif=False,
        write_endpoint_inputs=False,
    ):
        """
        NEB inputs has a special directory structure where inputs are in 00,
        01, 02, ....

        Args:
            output_dir (str): Directory to output the VASP input files
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            write_cif (bool): If true, writes a cif along with each POSCAR.
            write_path_cif (bool): If true, writes a cif for each image.
            write_endpoint_inputs (bool): If true, writes input files for
                running endpoint calculations.
        """
        output_dir = Path(output_dir)
        if make_dir_if_not_present and not output_dir.exists():
            output_dir.mkdir(parents=True)
        self.incar.write_file(str(output_dir / "INCAR"))
        self.kpoints.write_file(str(output_dir / "KPOINTS"))
        self.potcar.write_file(str(output_dir / "POTCAR"))

        for i, p in enumerate(self.poscars):
            d = output_dir / str(i).zfill(2)
            if not d.exists():
                d.mkdir(parents=True)
            p.write_file(str(d / "POSCAR"))
            if write_cif:
                p.structure.to(filename=str(d / f"{i}.cif"))
        if write_endpoint_inputs:
            end_point_param = MITRelaxSet(self.structures[0], user_incar_settings=self.user_incar_settings)

            for image in ["00", str(len(self.structures) - 1).zfill(2)]:
                end_point_param.incar.write_file(str(output_dir / image / "INCAR"))
                end_point_param.kpoints.write_file(str(output_dir / image / "KPOINTS"))
                end_point_param.potcar.write_file(str(output_dir / image / "POTCAR"))
        if write_path_cif:
            sites = set()
            lat = self.structures[0].lattice
            for site in chain(*(struct for struct in self.structures)):
                sites.add(PeriodicSite(site.species, site.frac_coords, lat))
            nebpath = Structure.from_sites(sorted(sites))
            nebpath.to(filename=str(output_dir / "path.cif"))


@dataclass
class MITMDSet(DictSet):
    """
    Class for writing a vasp md run. This DOES NOT do multiple stage runs.

    Args:
        structure (Structure): Input structure.
        start_temp (float): Starting temperature.
        end_temp (float): Final temperature.
        nsteps (int): Number of time steps for simulations. NSW parameter.
        time_step (float): The time step for the simulation. The POTIM
            parameter. Defaults to 2fs.
        spin_polarized (bool): Whether to do spin polarized calculations.
            The ISPIN parameter. Defaults to False.
        **kwargs: Other kwargs supported by DictSet.
    """

    structure: Structure | None = None
    start_temp: float = 0.0
    end_temp: float = 300.0
    nsteps: int = 1000
    time_step: float = 2
    spin_polarized: bool = False
    CONFIG = MITRelaxSet.CONFIG

    @property
    def incar_updates(self):
        """Get updates to the INCAR config for this calculation type."""
        # MD default settings
        return {
            "TEBEG": self.start_temp,
            "TEEND": self.end_temp,
            "NSW": self.nsteps,
            "EDIFF_PER_ATOM": 0.000001,
            "LSCALU": False,
            "LCHARG": False,
            "LPLANE": False,
            "LWAVE": True,
            "ISMEAR": 0,
            "NELMIN": 4,
            "LREAL": True,
            "BMIX": 1,
            "MAXMIX": 20,
            "NELM": 500,
            "NSIM": 4,
            "ISYM": 0,
            "ISIF": 0,
            "IBRION": 0,
            "NBLOCK": 1,
            "KBLOCK": 100,
            "SMASS": 0,
            "POTIM": self.time_step,
            "PREC": "Low",
            "ISPIN": 2 if self.spin_polarized else 1,
            "LDAU": False,
            "ENCUT": None,
        }

    @property
    def kpoints_updates(self) -> Kpoints | dict:
        """Get updates to the kpoints configuration for this calculation type."""
        return Kpoints.gamma_automatic()


@dataclass
class MPMDSet(DictSet):
    """
    This a modified version of the old MITMDSet pre 2018/03/12.

    This set serves as the basis for the amorphous skyline paper.

    (1) Aykol, M.; Dwaraknath, S. S.; Sun, W.; Persson, K. A. Thermodynamic
        Limit for Synthesis of Metastable Inorganic Materials. Sci. Adv. 2018,
        4 (4).

    Class for writing a vasp md run. This DOES NOT do multiple stage runs.
    Precision remains normal, to increase accuracy of stress tensor.

    Args:
        structure (Structure): Input structure.
        start_temp (int): Starting temperature.
        end_temp (int): Final temperature.
        nsteps (int): Number of time steps for simulations. NSW parameter.
        time_step (float): The time step for the simulation. The POTIM
            parameter. Defaults to None, which will set it automatically
            to 2.0 fs for non-hydrogen containing structures and 0.5 fs
            for hydrogen containing structures.
        spin_polarized (bool): Whether to do spin polarized calculations.
            The ISPIN parameter. Defaults to False.
        **kwargs: Other kwargs supported by DictSet.
    """

    start_temp: float = 0.0
    end_temp: float = 300.0
    nsteps: int = 1000
    time_step: float | None = None
    spin_polarized: bool = False
    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates = {
            "TEBEG": self.start_temp,
            "TEEND": self.end_temp,
            "NSW": self.nsteps,
            "EDIFF_PER_ATOM": 0.00001,
            "LSCALU": False,
            "LCHARG": False,
            "LPLANE": False,
            "LWAVE": True,
            "ISMEAR": 0,
            "NELMIN": 4,
            "LREAL": True,
            "BMIX": 1,
            "MAXMIX": 20,
            "NELM": 500,
            "NSIM": 4,
            "ISYM": 0,
            "ISIF": 0,
            "IBRION": 0,
            "NBLOCK": 1,
            "KBLOCK": 100,
            "SMASS": 0,
            "PREC": "Normal",
            "ISPIN": 2 if self.spin_polarized else 1,
            "LDAU": False,
            "ADDGRID": True,
            "ENCUT": None,
        }
        if not self.spin_polarized:
            updates["MAGMOM"] = None

        if self.time_step is None:
            if Element("H") in self.structure.species:  # type: ignore
                updates.update({"POTIM": 0.5, "NSW": self.nsteps * 4})
            else:
                updates["POTIM"] = 2.0
        else:
            updates["POTIM"] = self.time_step

        return updates

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Get updates to the kpoints configuration for this calculation type."""
        return Kpoints.gamma_automatic()


@dataclass
class MVLNPTMDSet(DictSet):
    """
    Class for writing a vasp md run in NPT ensemble.

    Notes:
        To eliminate Pulay stress, the default ENCUT is set to a rather large
        value of ENCUT, which is 1.5 * ENMAX.
    """

    start_temp: float = 0.0
    end_temp: float = 300.0
    nsteps: int = 1000
    time_step: float = 2
    spin_polarized: bool = False
    CONFIG = MITRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        # NPT-AIMD default settings
        updates = {
            "ALGO": "Fast",
            "ISIF": 3,
            "LANGEVIN_GAMMA": [10] * self.structure.ntypesp,  # type: ignore
            "LANGEVIN_GAMMA_L": 1,
            "MDALGO": 3,
            "PMASS": 10,
            "PSTRESS": 0,
            "SMASS": 0,
            "TEBEG": self.start_temp,
            "TEEND": self.end_temp,
            "NSW": self.nsteps,
            "EDIFF_PER_ATOM": 0.000001,
            "LSCALU": False,
            "LCHARG": False,
            "LPLANE": False,
            "LWAVE": True,
            "ISMEAR": 0,
            "NELMIN": 4,
            "LREAL": True,
            "BMIX": 1,
            "MAXMIX": 20,
            "NELM": 500,
            "NSIM": 4,
            "ISYM": 0,
            "IBRION": 0,
            "NBLOCK": 1,
            "KBLOCK": 100,
            "POTIM": self.time_step,
            "PREC": "Low",
            "ISPIN": 2 if self.spin_polarized else 1,
            "LDAU": False,
        }
        # Set NPT-AIMD ENCUT = 1.5 * VASP_default
        enmax = [self.potcar[i].keywords["ENMAX"] for i in range(self.structure.ntypesp)]  # type: ignore
        updates["ENCUT"] = max(enmax) * 1.5
        return updates

    @property
    def kpoints_updates(self) -> Kpoints | dict:
        """Get updates to the kpoints configuration for this calculation type."""
        return Kpoints.gamma_automatic()


@dataclass
class MVLScanRelaxSet(DictSet):
    """
    Class for writing a relax input set using Strongly Constrained and
    Appropriately Normed (SCAN) semilocal density functional.

    Notes:
        1. This functional is only available from VASP.5.4.3 upwards.

        2. Meta-GGA calculations require POTCAR files that include
        information on the kinetic energy density of the core-electrons,
        i.e. "PBE_52" or "PBE_54". Make sure the POTCAR including the
        following lines (see VASP wiki for more details):

            $ grep kinetic POTCAR
            kinetic energy-density
            mkinetic energy-density pseudized
            kinetic energy density (partial)

    Args:
        structure (Structure): input structure.
        vdw (str): set "rVV10" to enable SCAN+rVV10, which is a versatile
            van der Waals density functional by combing the SCAN functional
            with the rVV10 non-local correlation functional.
        **kwargs: Other kwargs supported by DictSet.
    """

    user_potcar_functional: UserPotcarFunctional = "PBE_52"
    _valid_potcars = ("PBE_52", "PBE_54")
    CONFIG = MPRelaxSet.CONFIG

    def __post_init__(self):
        super().__post_init__()
        if self.user_potcar_functional not in ("PBE_52", "PBE_54"):
            raise ValueError("SCAN calculations required PBE_52 or PBE_54!")

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates = {
            "ADDGRID": True,
            "EDIFF": 1e-5,
            "EDIFFG": -0.05,
            "LASPH": True,
            "LDAU": False,
            "METAGGA": "SCAN",
            "NELM": 200,
        }
        if self.vdw and self.vdw.lower() == "rvv10":
            updates["BPARAM"] = 15.7  # This is the correct BPARAM for SCAN+rVV10
        return updates


@dataclass
class LobsterSet(DictSet):
    """Input set to prepare VASP runs that can be digested by Lobster (See cohp.de).

    Args:
        structure (Structure): input structure.
        isym (int): ISYM entry for INCAR, only isym=-1 and isym=0 are allowed
        ismear (int): ISMEAR entry for INCAR, only ismear=-5 and ismear=0 are allowed
        reciprocal_density (int): density of k-mesh by reciprocal volume
        user_supplied_basis (dict): dict including basis functions for all elements in
            structure, e.g. {"Fe": "3d 3p 4s", "O": "2s 2p"}; if not supplied, a
            standard basis is used
        address_basis_file (str): address to a file similar to
            "BASIS_PBE_54_standard.yaml" in pymatgen.io.lobster.lobster_basis
        user_potcar_settings (dict): dict including potcar settings for all elements in
            structure, e.g. {"Fe": "Fe_pv", "O": "O"}; if not supplied, a standard basis
            is used.
        **kwargs: Other kwargs supported by DictSet.
    """

    isym: int = 0
    ismear: int = -5
    reciprocal_density: int | None = None
    address_basis_file: str | None = None
    user_supplied_basis: dict | None = None

    # newest potcars are preferred
    # Choose PBE_54 unless the user specifies a different potcar_functional
    user_potcar_functional: UserPotcarFunctional = "PBE_54"

    CONFIG = MPRelaxSet.CONFIG
    _valid_potcars = ("PBE_52", "PBE_54")

    def __post_init__(self):
        super().__post_init__()
        warnings.warn("Make sure that all parameters are okay! This is a brand new implementation.")

        if self.isym not in (-1, 0):
            raise ValueError("Lobster cannot digest WAVEFUNCTIONS with symmetry. isym must be -1 or 0")
        if self.ismear not in (-5, 0):
            raise ValueError("Lobster usually works with ismear=-5 or ismear=0")

        self._config_dict["POTCAR"]["W"] = "W_sv"

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Get updates to the kpoints configuration for this calculation type."""
        # test, if this is okay
        return {"reciprocal_density": self.reciprocal_density if self.reciprocal_density else 310}

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        from pymatgen.io.lobster import Lobsterin

        potcar_symbols = self.poscar.site_symbols

        # predefined basis! Check if the basis is okay! (charge spilling and bandoverlaps!)
        if self.user_supplied_basis is None and self.address_basis_file is None:
            basis = Lobsterin.get_basis(structure=self.structure, potcar_symbols=potcar_symbols)  # type: ignore
        elif self.address_basis_file is not None:
            basis = Lobsterin.get_basis(
                structure=self.structure,  # type: ignore
                potcar_symbols=potcar_symbols,
                address_basis_file=self.address_basis_file,
            )
        elif self.user_supplied_basis is not None:
            # test if all elements from structure are in user_supplied_basis
            for atom_type in self.structure.symbol_set:  # type: ignore
                if atom_type not in self.user_supplied_basis:
                    raise ValueError(f"There are no basis functions for the atom type {atom_type}")
            basis = [f"{key} {value}" for key, value in self.user_supplied_basis.items()]

        lobsterin = Lobsterin(settingsdict={"basisfunctions": basis})
        nbands = lobsterin._get_nbands(structure=self.structure)  # type: ignore

        return {
            "EDIFF": 1e-6,
            "NSW": 0,
            "LWAVE": True,
            "ISYM": self.isym,
            "NBANDS": nbands,
            "IBRION": -1,
            "ISMEAR": self.ismear,
            "LORBIT": 11,
            "ICHARG": 0,
            "ALGO": "Normal",
        }


def get_vasprun_outcar(path: str | Path, parse_dos: bool = True, parse_eigen: bool = True) -> tuple[Vasprun, Outcar]:
    """
    Get a Vasprun and Outcar from a directory.

    Args:
        path: Path to get the vasprun.xml and OUTCAR.
        parse_dos: Whether to parse dos. Defaults to True.
        parse_eigen: Whether to parse eigenvalue. Defaults to True.

    Returns:
        Vasprun and Outcar files.
    """
    path = Path(path)
    vruns = list(glob(str(path / "vasprun.xml*")))
    outcars = list(glob(str(path / "OUTCAR*")))

    if len(vruns) == 0 or len(outcars) == 0:
        raise ValueError(f"Unable to get vasprun.xml/OUTCAR from prev calculation in {path}")
    vsfile_fullpath = str(path / "vasprun.xml")
    outcarfile_fullpath = str(path / "OUTCAR")
    vsfile = vsfile_fullpath if vsfile_fullpath in vruns else sorted(vruns)[-1]
    outcarfile = outcarfile_fullpath if outcarfile_fullpath in outcars else sorted(outcars)[-1]
    return (
        Vasprun(vsfile, parse_dos=parse_dos, parse_eigen=parse_eigen),
        Outcar(outcarfile),
    )


def get_structure_from_prev_run(vasprun, outcar=None):
    """
    Process structure from previous run.

    Args:
        vasprun (Vasprun): Vasprun that contains the final structure from previous run.
        outcar (Outcar): Outcar that contains the magnetization info from previous run.

    Returns:
        The magmom-decorated structure that can be passed to get VASP input files, e.g.
        get_kpoints.
    """
    structure = vasprun.final_structure

    site_properties = {}
    # magmom
    if vasprun.is_spin:
        if outcar and outcar.magnetization:
            site_properties["magmom"] = [i["tot"] for i in outcar.magnetization]
        else:
            site_properties["magmom"] = vasprun.parameters["MAGMOM"]
    # ldau
    if vasprun.parameters.get("LDAU", False):
        for k in ("LDAUU", "LDAUJ", "LDAUL"):
            vals = vasprun.incar[k]
            m = {}
            l_val = []
            s = 0
            for site in structure:
                if site.specie.symbol not in m:
                    m[site.specie.symbol] = vals[s]
                    s += 1
                l_val.append(m[site.specie.symbol])
            if len(l_val) == len(structure):
                site_properties.update({k.lower(): l_val})
            else:
                raise ValueError(f"length of list {l_val} not the same as structure")

    return structure.copy(site_properties=site_properties)


def standardize_structure(structure, sym_prec=0.1, international_monoclinic=True):
    """
    Get the symmetrically standardized structure.

    Args:
        structure (Structure): The structure.
        sym_prec (float): Tolerance for symmetry finding for standardization.
        international_monoclinic (bool): Whether to use international
            convention (vs Curtarolo) for monoclinic. Defaults True.

    Returns:
        The symmetrized structure.
    """
    sym_finder = SpacegroupAnalyzer(structure, symprec=sym_prec)
    new_structure = sym_finder.get_primitive_standard_structure(international_monoclinic=international_monoclinic)

    # the primitive structure finding has had several bugs in the past
    # defend through validation
    vpa_old = structure.volume / len(structure)
    vpa_new = new_structure.volume / len(new_structure)

    if abs(vpa_old - vpa_new) / vpa_old > 0.02:
        raise ValueError(f"Standardizing cell failed! VPA old: {vpa_old}, VPA new: {vpa_new}")

    matcher = StructureMatcher()
    if not matcher.fit(structure, new_structure):
        raise ValueError("Standardizing cell failed! Old structure doesn't match new.")

    return new_structure


class BadInputSetWarning(UserWarning):
    """Warning class for bad but legal VASP inputs."""


def batch_write_input(
    structures,
    vasp_input_set=MPRelaxSet,
    output_dir=".",
    make_dir_if_not_present=True,
    subfolder=None,
    sanitize=False,
    include_cif=False,
    potcar_spec=False,
    zip_output=False,
    **kwargs,
):
    """
    Batch write vasp input for a sequence of structures to
    output_dir, following the format output_dir/{group}/{formula}_{number}.

    Args:
        structures ([Structure]): Sequence of Structures.
        vasp_input_set (VaspInputSet): VaspInputSet class that creates
            vasp input files from structures. Note that a class should be
            supplied. Defaults to MPRelaxSet.
        output_dir (str): Directory to output files. Defaults to current
            directory ".".
        make_dir_if_not_present (bool): Create the directory if not present.
            Defaults to True.
        subfolder (callable): Function to create subdirectory name from
            structure. Defaults to simply "formula_count".
        sanitize (bool): Boolean indicating whether to sanitize the
            structure before writing the VASP input files. Sanitized output
            are generally easier for viewing and certain forms of analysis.
            Defaults to False.
        include_cif (bool): Whether to output a CIF as well. CIF files are
            generally better supported in visualization programs.
        potcar_spec (bool): Instead of writing the POTCAR, write a "POTCAR.spec".
                This is intended to help sharing an input set with people who might
                not have a license to specific Potcar files. Given a "POTCAR.spec",
                the specific POTCAR file can be re-generated using pymatgen with the
                "generate_potcar" function in the pymatgen CLI.
        zip_output (bool): If True, output will be zipped into a file with the
            same name as the InputSet (e.g., MPStaticSet.zip)
        **kwargs: Additional kwargs are passed to the vasp_input_set class
            in addition to structure.
    """
    output_dir = Path(output_dir)
    for idx, site in enumerate(structures):
        formula = re.sub(r"\s+", "", site.formula)
        if subfolder is not None:
            subdir = subfolder(site)
            d = output_dir / subdir
        else:
            d = output_dir / f"{formula}_{idx}"
        if sanitize:
            site = site.copy(sanitize=True)
        v = vasp_input_set(site, **kwargs)
        v.write_input(
            str(d),
            make_dir_if_not_present=make_dir_if_not_present,
            include_cif=include_cif,
            potcar_spec=potcar_spec,
            zip_output=zip_output,
        )


_dummy_structure = Structure(
    [1, 0, 0, 0, 1, 0, 0, 0, 1],
    ["I"],
    [[0, 0, 0]],
    site_properties={"magmom": [[0, 0, 1]]},
)


def get_valid_magmom_struct(structure, inplace=True, spin_mode="auto"):
    """
    Make sure that the structure has valid magmoms based on the kind of calculation.

    Fill in missing Magmom values.

    Args:
        structure: The input structure
        inplace: True: edit magmoms of the input structure; False: return new structure
        spin_mode: "scalar"/"vector"/"none"/"auto" only first letter (s/v/n) is needed.
            dictates how the spin configuration will be determined.

            - auto: read the existing magmom values and decide
            - scalar: use a single scalar value (for spin up/down)
            - vector: use a vector value for spin-orbit systems
            - none: Remove all the magmom information

    Returns:
        New structure if inplace is False
    """
    default_values = {"s": 1.0, "v": [1.0, 1.0, 1.0], "n": None}
    if spin_mode[0].lower() == "a":
        mode = "n"
        for site in structure:
            if "magmom" not in site.properties or site.properties["magmom"] is None:
                pass
            elif isinstance(site.properties["magmom"], (float, int)):
                if mode == "v":
                    raise TypeError("Magmom type conflict")
                mode = "s"
                if isinstance(site.properties["magmom"], int):
                    site.properties["magmom"] = float(site.properties["magmom"])
            elif len(site.properties["magmom"]) == 3:
                if mode == "s":
                    raise TypeError("Magmom type conflict")
                mode = "v"
            else:
                raise TypeError("Unrecognized Magmom Value")
    else:
        mode = spin_mode[0].lower()

    new_struct = structure.copy() if not inplace else structure
    for site in new_struct:
        if mode == "n":
            if "magmom" in site.properties:
                site.properties.pop("magmom")
        elif "magmom" not in site.properties or site.properties["magmom"] is None:
            site.properties["magmom"] = default_values[mode]

    if not inplace:
        return new_struct
    return None


@dataclass
class MPAbsorptionSet(DictSet):
    """
    MP input set for generating frequency dependent dielectrics.

    Two modes are supported: "IPA" or "RPA".
    A typical sequence is mode="STATIC" -> mode="IPA" -> mode="RPA"(optional)
    For all steps other than the first one (static), the
    recommendation is to use from_prev_calculation on the preceding run in
    the series. It is important to ensure Gamma centred kpoints for the RPA step.

    Args:
        structure (Structure): Input structure.
        mode (str): Supported modes are "IPA", "RPA"
        copy_wavecar (bool): Whether to copy the WAVECAR from a previous run.
            Defaults to True.
        nbands (int): For subsequent calculations, it is generally
            recommended to perform NBANDS convergence starting from the
            NBANDS of the previous run for DIAG, and to use the exact same
            NBANDS for RPA. This parameter is used by
            from_previous_calculation to set nband.
        nbands_factor (int): Multiplicative factor for NBANDS when starting
            from a previous calculation. Only applies if mode=="IPA".
            Need to be tested for convergence.
        reciprocal_density: the k-points density
        nkred: the reduced number of kpoints to calculate, equal to the k-mesh.
            Only applies in "RPA" mode because of the q->0 limit.
        nedos: the density of DOS, default: 2001.
        **kwargs: All kwargs supported by DictSet. Typically, user_incar_settings is a
            commonly used option.
    """

    # CONFIG = _load_yaml_config("MPAbsorptionSet")

    mode: str = "IPA"
    copy_wavecar: bool = True
    nbands_factor: float = 2
    reciprocal_density: float = 400
    nkred: tuple[int, int, int] | None = None
    nedos: int = 2001
    inherit_incar: bool = True
    force_gamma: bool = True
    CONFIG = MPRelaxSet.CONFIG
    nbands: int | None = None
    SUPPORTED_MODES = ("IPA", "RPA")

    def __post_init__(self):
        """Validate settings"""
        super().__post_init__()
        self.mode = self.mode.upper()
        if self.mode not in MPAbsorptionSet.SUPPORTED_MODES:
            raise ValueError(f"{self.mode} not one of the support modes : {MPAbsorptionSet.SUPPORTED_MODES}")

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Get updates to the kpoints configuration for this calculation type.

        Generate gamma center k-points mesh grid for optical calculation. It is not
        mandatory for 'ALGO = Exact', but is requested by 'ALGO = CHI' calculation.
        """
        return {"reciprocal_density": self.reciprocal_density}

    @property
    def incar_updates(self) -> dict:
        """Get updates to the INCAR config for this calculation type."""
        updates = {
            "ALGO": "Exact",
            "EDIFF": 1.0e-8,
            "IBRION": -1,
            "ICHARG": 1,
            "ISMEAR": 0,
            "SIGMA": 0.01,
            "LWAVE": True,
            "LREAL": False,  # for small cell it's more efficient to use reciprocal
            "NELM": 100,
            "NSW": 0,
            "LOPTICS": True,
            "CSHIFT": 0.1,
            "NEDOS": self.nedos,
        }

        if self.mode == "RPA":
            # Default parameters for the response function calculation. NELM has to be
            # set to 1. NOMEGA is set to 1000 in order to get smooth spectrum
            updates.update({"ALGO": "CHI", "NELM": 1, "NOMEGA": 1000, "EDIFF": None, "LOPTICS": None, "LWAVE": None})

        if self.prev_vasprun is not None and self.mode == "IPA":
            prev_nbands = int(self.prev_vasprun.parameters["NBANDS"]) if self.nbands is None else self.nbands
            updates["NBANDS"] = int(np.ceil(prev_nbands * self.nbands_factor))

        if self.prev_vasprun is not None and self.mode == "RPA":
            # Since in the optical calculation, only the q->0 transition is of interest,
            # we can reduce the number of q by the factor of the number of kpoints in
            # each corresponding x, y, z directions. This will reduce the computational
            # work by factor of 1/nkredx*nkredy*nkredz. An isotropic NKRED can be used
            # for cubic lattices, but using NKREDX, NKREDY, NKREDZ are more sensible for
            # other lattices.
            self.nkred = self.prev_vasprun.kpoints.kpts[0] if self.nkred is None else self.nkred
            updates.update({"NKREDX": self.nkred[0], "NKREDY": self.nkred[1], "NKREDZ": self.nkred[2]})

        return updates


def _get_ispin(vasprun: Vasprun | None, outcar: Outcar | None) -> int:
    """Get value of ISPIN depending on the magnetisation in the OUTCAR and vasprun."""
    if outcar is not None and outcar.magnetization is not None:
        # Turn off spin when magmom for every site is smaller than 0.02.
        site_magmom = np.array([i["tot"] for i in outcar.magnetization])
        return 2 if np.any(np.abs(site_magmom) > 0.02) else 1
    if vasprun is not None:
        return 2 if vasprun.is_spin else 1
    return 2


def _combine_kpoints(*kpoints_objects: Kpoints) -> Kpoints:
    """Combine k-points files together."""
    labels = []
    kpoints = []
    weights = []

    for kpoints_object in filter(None, kpoints_objects):
        if kpoints_object.style != Kpoints.supported_modes.Reciprocal:
            raise ValueError("Can only combine kpoints with style=Kpoints.supported_modes.Reciprocal")
        if kpoints_object.labels is None:
            labels.append([""] * len(kpoints_object.kpts))
        else:
            labels.append(kpoints_object.labels)

        weights.append(kpoints_object.kpts_weights)
        kpoints.append(kpoints_object.kpts)

    labels = np.concatenate(labels).tolist()
    weights = np.concatenate(weights).tolist()
    kpoints = np.concatenate(kpoints)
    return Kpoints(
        comment="Combined k-points",
        style=Kpoints.supported_modes.Reciprocal,
        num_kpts=len(kpoints),
        kpts=kpoints,
        labels=labels,
        kpts_weights=weights,
    )


def _apply_incar_updates(incar, updates, skip: Sequence[str] = ()) -> None:
    """
    Apply updates to an INCAR file.

    Args:
        incar (dict): An incar.
        updates (dict): Updates to apply.
        skip (list of str): Keys to skip.
    """
    for k, v in updates.items():
        if k in skip:
            continue

        if v is None:
            incar.pop(k, None)
        else:
            incar[k] = v


def _remove_unused_incar_params(incar, skip: Sequence[str] = ()) -> None:
    """
    Remove INCAR parameters that are not actively used by VASP.

    Args:
        incar (dict): An incar.
        skip (list of str): Keys to skip.
    """
    # Turn off IBRION/ISIF/POTIM if NSW = 0
    opt_flags = ["EDIFFG", "IBRION", "ISIF", "POTIM"]
    if incar.get("NSW", 0) == 0:
        for opt_flag in opt_flags:
            if opt_flag not in skip:
                incar.pop(opt_flag, None)

    # Remove MAGMOMs if they aren't used
    if incar.get("ISPIN", 1) == 1 and "MAGMOM" not in skip:
        incar.pop("MAGMOM", None)

    # Turn off +U flags if +U is not even used
    ldau_flags = ["LDAUU", "LDAUJ", "LDAUL", "LDAUTYPE"]
    if incar.get("LDAU", False) is False:
        for ldau_flag in ldau_flags:
            if ldau_flag not in skip:
                incar.pop(ldau_flag, None)


def _get_nedos(vasprun: Vasprun | None, dedos: float) -> int:
    """Automatic setting of nedos using the energy range and the energy step."""
    if vasprun is None:
        return 2000

    emax = max(eigs.max() for eigs in vasprun.eigenvalues.values())
    emin = min(eigs.min() for eigs in vasprun.eigenvalues.values())
    return int((emax - emin) / dedos)


def auto_kspacing(bandgap, bandgap_tol):
    """Automatically set kspacing based on the bandgap"""
    if bandgap is None or bandgap <= bandgap_tol:  # metallic
        return 0.22

    rmin = max(1.5, 25.22 - 2.87 * bandgap)  # Eq. 25
    kspacing = 2 * np.pi * 1.0265 / (rmin - 1.0183)  # Eq. 29

    # cap kspacing at a max of 0.44, per internal benchmarking
    return min(kspacing, 0.44)
