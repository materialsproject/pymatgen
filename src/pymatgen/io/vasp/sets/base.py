# ruff: noqa: PGH003
"""
This module defines the VaspInputSet abstract base class and a concrete implementation for the parameters developed
and tested by the core team of pymatgen, including the Materials Virtual Lab, Materials Project and the MIT high
throughput project. The basic concept behind an input set is to specify a scheme to generate a consistent set of VASP
inputs from a structure without further user intervention. This ensures comparability across runs.

Read the following carefully before implementing new input sets:

1. 99% of what needs to be done can be done by specifying user_incar_settings to override some of the defaults of
   various input sets. Unless there is an extremely good reason to add a new set, **do not** add one. e.g. if you want
   to turn the Hubbard U off, just set "LDAU": False as a user_incar_setting.
2. All derivative input sets should inherit appropriate configurations (e.g., from MPRelaxSet), and more often than
   not, VaspInputSet should be the superclass. Superclass delegation should be used where possible. In particular,
   you are not supposed to implement your own as_dict or from_dict for derivative sets unless you know what you are
   doing. Improper overriding the as_dict and from_dict protocols is the major cause of implementation headaches. If
   you need an example, look at how the MPStaticSet is initialized.

The above are recommendations. The following are **UNBREAKABLE** rules:

1. All input sets must take in a structure, list of structures or None as the first argument. If None, the input set
   should perform a stateless initialization and before any output can be written, a structure must be set.
2. user_incar_settings, user_kpoints_settings and user_<whatever>_settings are ABSOLUTE. Any new sets you implement
   must obey this. If a user wants to override your settings, you assume he knows what he is doing. Do not
   magically override user supplied settings. You can issue a warning if you think the user is wrong.
3. All input sets must save all supplied args and kwargs as instance variables. e.g. self.arg = arg and
   self.kwargs = kwargs in the __init__. This ensures the as_dict and from_dict work correctly.
"""

from __future__ import annotations

import abc
import itertools
import re
import warnings
from collections.abc import Sequence
from copy import deepcopy
from dataclasses import dataclass, field
from glob import glob
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Union, cast

import numpy as np
from monty.dev import deprecated
from monty.json import MSONable
from monty.serialization import loadfn

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import SiteCollection, Structure
from pymatgen.io.core import InputGenerator
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.util.typing import Kpoint

if TYPE_CHECKING:
    from typing import Any

    from typing_extensions import Self


UserPotcarFunctional = Union[
    Literal["PBE", "PBE_52", "PBE_54", "LDA", "LDA_52", "LDA_54", "PW91", "LDA_US", "PW91_US"], None
]
MODULE_DIR = Path(__file__).resolve().parent.parent


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


@dataclass
class VaspInputSet(InputGenerator, abc.ABC):
    """
    Base class representing a set of VASP input parameters with a structure
    supplied as init parameters and initialized from a dict of settings.
    This allows arbitrary settings to be input. In general,
    this is rarely used directly unless there is a source of settings in yaml
    format (e.g., from a REST interface). It is typically used by other
    VaspInputSets for initialization.

    Special consideration should be paid to the way the MAGMOM initialization
    for the INCAR is done. The initialization differs depending on the type of
    structure and the configuration settings. The order in which the magmom is
    determined is as follows:

    1. If the site is specified in user_incar_settings, use that setting.
    2. If the site itself has a magmom setting (i.e. site.properties["magmom"] = float),
        that is used. This can be set with structure.add_site_property().
    3. If the species of the site has a spin setting, that is used. This can be set
        with structure.add_spin_by_element().
    4. If the species itself has a particular setting in the config file, that
       is used, e.g. Mn3+ may have a different magmom than Mn4+.
    5. Lastly, the element symbol itself is checked in the config file. If
       there are no settings, a default value of 0.6 is used.

    Args:
        structure (Structure): The Structure to create inputs for. If None, the input
            set is initialized without a Structure but one must be set separately before
            the inputs are generated.
        config_dict (dict): The config dictionary to use.
        files_to_transfer (dict): A dictionary of {filename: filepath}. This allows the
            transfer of files from a previous calculation.
        user_incar_settings (dict): User INCAR settings. This allows a user to override
            INCAR settings, e.g. setting a different MAGMOM for various elements or
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
            by supplying a dict. e.g. {"reciprocal_density": 1000}. User can also
            supply Kpoints object.
        user_potcar_settings (dict): Allow user to override POTCARs. e.g. {"Gd":
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
        auto_kspacing (bool): If true, determines the value of KSPACING from the bandgap
            of a previous calculation.
        auto_ismear (bool): If true, the values for ISMEAR and SIGMA will be set
            automatically depending on the bandgap of the system. If the bandgap is not
            known (e.g., there is no previous VASP directory) then ISMEAR=0 and
            SIGMA=0.2; if the bandgap is zero (a metallic system) then ISMEAR=2 and
            SIGMA=0.2; if the system is an insulator, then ISMEAR=-5 (tetrahedron
            smearing). Note, this only works when generating the input set from a
            previous VASP directory.
        auto_ispin (bool) = False:
            If generating input set from a previous calculation, this controls whether
            to disable magnetisation (ISPIN = 1) if the absolute value of all magnetic
            moments are less than 0.02.
        auto_lreal (bool) = False:
            If True, automatically use the VASP recommended LREAL based on cell size.
        auto_metal_kpoints
            If true and the system is metallic, try and use ``reciprocal_density_metal``
            instead of ``reciprocal_density`` for metallic systems. Note, this only works
            if the bandgap is not None.
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
    auto_kspacing: bool = False
    auto_ismear: bool = False
    auto_ispin: bool = False
    auto_lreal: bool = False
    auto_metal_kpoints: bool = False
    bandgap_tol: float = 1e-4
    bandgap: float | None = None
    prev_incar: str | dict | None = None
    prev_kpoints: str | Kpoints | None = None
    _valid_potcars: Sequence[str] | None = None

    def __post_init__(self):
        """Perform validation."""
        user_potcar_functional = self.user_potcar_functional
        if (valid_potcars := self._valid_potcars) and user_potcar_functional not in valid_potcars:
            raise ValueError(f"Invalid {user_potcar_functional=}, must be one of {valid_potcars}")

        if hasattr(self, "CONFIG"):
            self.config_dict = self.CONFIG

        self._config_dict = deepcopy(self.config_dict)

        # these have been left to stay consistent with previous API
        self.user_incar_settings = self.user_incar_settings or {}
        self.user_kpoints_settings = self.user_kpoints_settings or {}

        self.vdw = self.vdw.lower() if isinstance(self.vdw, str) else self.vdw
        if self.user_incar_settings.get("KSPACING") and self.user_kpoints_settings:
            # self.user_kpoints_settings will never be `None` because it is set to
            # an empty dict if it is `None`.
            warnings.warn(
                "You have specified KSPACING and also supplied KPOINTS "
                "settings. KSPACING only has effect when there is no "
                "KPOINTS file. Since both settings were given, pymatgen"
                "will generate a KPOINTS file and ignore KSPACING."
                "Remove the `user_kpoints_settings` argument to enable KSPACING.",
                BadInputSetWarning,
            )

        if self.vdw:
            vdw_par = loadfn(MODULE_DIR / "vdW_parameters.yaml")
            if vdw_param := vdw_par.get(self.vdw):
                self._config_dict["INCAR"].update(vdw_param)
            else:
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
        self._ispin = None

    @deprecated(message="get_vasp_input will be removed in a future version of pymatgen. Use get_input_set instead.")
    def get_vasp_input(self, structure=None) -> VaspInput:
        """Get a VaspInput object.

        Returns:
            VaspInput.
        """
        return self.get_input_set(structure=structure)

    def write_input(
        self,
        output_dir: str,
        make_dir_if_not_present: bool = True,
        include_cif: bool | str = False,
        potcar_spec: bool = False,
        zip_output: bool | str = False,
    ) -> None:
        """Write a set of VASP input to a directory.

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
                same name as the InputSet (e.g., MPStaticSet.zip).
        """
        vasp_input = self.get_input_set(potcar_spec=potcar_spec)

        cif_name = None
        if include_cif:
            struct = vasp_input["POSCAR"].structure
            cif_name = f"{output_dir}/{struct.formula.replace(' ', '')}.cif"

        vasp_input.write_input(
            output_dir=output_dir,
            make_dir_if_not_present=make_dir_if_not_present,
            cif_name=cif_name,
            zip_name=f"{type(self).__name__}.zip" if zip_output else None,
            files_to_transfer=self.files_to_transfer,
        )

    def as_dict(self, verbosity=2):
        """
        Args:
            verbosity: Verbosity for generated dict. If 1, structure is
            excluded.

        Returns:
            dict: MSONable VaspInputSet representation.
        """
        dct = MSONable.as_dict(self)
        if verbosity == 1:
            dct.pop("structure", None)
        return dct

    @property  # type: ignore
    def structure(self) -> Structure:  # noqa: F811
        """Structure."""
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
        """Get a VASP input set.

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
            raise ValueError("Either structure or prev_dir must be set")

        self._set_previous(prev_dir)

        if structure is not None:
            self.structure = structure

        return VaspInput(
            incar=self.incar,
            kpoints=self.kpoints,
            poscar=self.poscar,
            potcar="\n".join(self.potcar_symbols) if potcar_spec else self.potcar,
            potcar_spec=potcar_spec,
        )

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
        return {}

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Updates to the kpoints configuration for this calculation type.

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
            if self.auto_ispin:
                # turn off spin when magmom for every site is smaller than 0.02.
                self._ispin = _get_ispin(vasprun, outcar)

            self.structure = get_structure_from_prev_run(vasprun, outcar)

    @property
    def incar(self) -> Incar:
        """The INCAR."""
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
        if self.auto_ispin and (self._ispin is not None):
            auto_updates["ISPIN"] = self._ispin

        # breaking change - order in which settings applied inconsistent with atomate2
        # apply updates from input set generator to SETTINGS
        # _apply_incar_updates(settings, incar_updates)

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

        for key, setting in settings.items():
            if key == "MAGMOM":
                mag = []
                for site in structure:
                    if uic_magmom := self.user_incar_settings.get("MAGMOM", {}).get(site.species_string):
                        mag.append(uic_magmom)
                    elif hasattr(site, "magmom"):
                        mag.append(site.magmom)
                    elif getattr(site.specie, "spin", None) is not None:
                        mag.append(site.specie.spin)
                    elif str(site.specie) in setting:
                        if site.specie.symbol == "Co" and setting[str(site.specie)] <= 1.0:
                            warnings.warn(
                                "Co without an oxidation state is initialized as low spin by default in Pymatgen. "
                                "If this default behavior is not desired, please set the spin on the magmom on the "
                                "site directly to ensure correct initialization."
                            )
                        mag.append(setting.get(str(site.specie)))
                    else:
                        if site.specie.symbol == "Co":
                            warnings.warn(
                                "Co without an oxidation state is initialized as low spin by default in Pymatgen. "
                                "If this default behavior is not desired, please set the spin on the magmom on the "
                                "site directly to ensure correct initialization."
                            )
                        mag.append(setting.get(site.specie.symbol, 0.6))
                incar[key] = mag
            elif key in ("LDAUU", "LDAUJ", "LDAUL"):
                if hubbard_u:
                    if hasattr(structure[0], key.lower()):
                        m = {site.specie.symbol: getattr(site, key.lower()) for site in structure}
                        incar[key] = [m[sym] for sym in poscar.site_symbols]
                        # lookup specific LDAU if specified for most_electroneg atom
                    elif most_electro_neg in setting and isinstance(setting[most_electro_neg], dict):
                        incar[key] = [setting[most_electro_neg].get(sym, 0) for sym in poscar.site_symbols]
                        # else, use fallback LDAU value if it exists
                    else:
                        incar[key] = [
                            setting.get(sym, 0) if isinstance(setting.get(sym, 0), (float, int)) else 0
                            for sym in poscar.site_symbols
                        ]
            elif key.startswith("EDIFF") and key != "EDIFFG":
                if "EDIFF" not in settings and key == "EDIFF_PER_ATOM":
                    incar["EDIFF"] = float(setting) * len(structure)
                else:
                    incar["EDIFF"] = float(settings["EDIFF"])
            elif key == "KSPACING" and self.auto_kspacing:
                # default to metal if no prev calc available
                bandgap = 0 if self.bandgap is None else self.bandgap
                incar[key] = auto_kspacing(bandgap, self.bandgap_tol)
            else:
                incar[key] = setting
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
                auto_updates.update(ISMEAR=0, SIGMA=0.2)
            elif self.bandgap <= self.bandgap_tol:
                auto_updates.update(ISMEAR=2, SIGMA=0.2)  # metal
            else:
                auto_updates.update(ISMEAR=-5, SIGMA=0.05)  # insulator

        if self.auto_lreal:
            auto_updates.update(LREAL=_get_recommended_lreal(structure))

        # apply updates from auto options, careful not to override user_incar_settings
        _apply_incar_updates(incar, auto_updates, skip=list(self.user_incar_settings))

        # apply updates from input set generator to INCAR
        _apply_incar_updates(incar, incar_updates, skip=list(self.user_incar_settings))

        # Finally, re-apply `self.user_incar_settings` to make sure any accidentally
        # overwritten settings are changed back to the intended values.
        # skip dictionary parameters to avoid dictionaries appearing in the INCAR
        _apply_incar_updates(incar, self.user_incar_settings, skip=["LDAUU", "LDAUJ", "LDAUL", "MAGMOM"])

        # Remove unused INCAR parameters
        _remove_unused_incar_params(incar, skip=list(self.user_incar_settings))

        kpoints = self.kpoints
        if kpoints is not None:
            # unset KSPACING as we are using a KPOINTS file
            incar.pop("KSPACING", None)
        elif "KSPACING" in incar and "KSPACING" not in self.user_incar_settings and "KSPACING" in prev_incar:
            # prefer to inherit KSPACING from previous INCAR if it exists
            # TODO: Is that we actually want to do? Copied from current pymatgen inputsets
            incar["KSPACING"] = prev_incar["KSPACING"]

        # Ensure adequate number of KPOINTS are present for the tetrahedron method
        # (ISMEAR=-5). If KSPACING is in the INCAR file the number of kpoints is not
        # known before calling VASP, but a warning is raised when the KSPACING value is
        # > 0.5 (2 reciprocal Angstrom). An error handler in Custodian is available to
        # correct overly large KSPACING values (small number of kpoints) if necessary.
        if kpoints is not None and np.prod(kpoints.kpts) < 4 and incar.get("ISMEAR", 0) == -5:
            incar["ISMEAR"] = 0

        if incar.get("KSPACING", 0) > 0.5 and incar.get("ISMEAR", 0) == -5:
            warnings.warn(
                "Large KSPACING value detected with ISMEAR = -5. Ensure that VASP "
                "generates an adequate number of KPOINTS, lower KSPACING, or "
                "set ISMEAR = 0",
                BadInputSetWarning,
            )

        ismear = incar.get("ISMEAR", 1)
        sigma = incar.get("SIGMA", 0.2)
        if (
            all(elem.is_metal for elem in structure.composition)
            and incar.get("NSW", 0) > 0
            and (ismear < 0 or (ismear == 0 and sigma > 0.05))
        ):
            msg = ""
            if ismear < 0:
                msg = f"Relaxation of likely metal with ISMEAR < 0 ({ismear})."
            elif ismear == 0 and sigma > 0.05:
                msg = f"ISMEAR = 0 with a small SIGMA ({sigma}) detected."
            warnings.warn(
                f"{msg} See VASP recommendations on ISMEAR for metals (https://www.vasp.at/wiki/index.php/ISMEAR).",
                BadInputSetWarning,
                stacklevel=1,
            )

        return incar

    @property
    def poscar(self) -> Poscar:
        """Poscar."""
        if self.structure is None:
            raise RuntimeError("No structure is associated with the input set!")
        site_properties = self.structure.site_properties
        return Poscar(
            self.structure,
            velocities=site_properties.get("velocities"),
            predictor_corrector=site_properties.get("predictor_corrector"),
            predictor_corrector_preamble=self.structure.properties.get("predictor_corrector_preamble"),
            lattice_velocities=self.structure.properties.get("lattice_velocities"),
        )

    @property
    def potcar_functional(self) -> UserPotcarFunctional:
        """The functional used for POTCAR generation."""
        return self.user_potcar_functional

    @property
    def nelect(self) -> float:
        """The default number of electrons for a given structure."""
        if self.structure is None:
            raise RuntimeError("No structure is associated with the input set!")

        n_electrons_by_element = {p.element: p.nelectrons for p in self.potcar}
        n_elect = sum(
            num_atoms * n_electrons_by_element[el.symbol] for el, num_atoms in self.structure.composition.items()
        )

        return n_elect - (self.structure.charge if self.use_structure_charge else 0)

    @property
    def kpoints(self) -> Kpoints | None:
        """The kpoints file."""
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
                    kpts=tuple(i[0] for i in mesh),
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
            mesh = sga.get_ir_reciprocal_mesh(zero_weighted_kpoints.kpts[0])  # type: ignore[arg-type]
            zero_weighted_kpoints = Kpoints(
                comment="Uniform grid",
                style=Kpoints.supported_modes.Reciprocal,
                num_kpts=len(mesh),
                kpts=tuple(i[0] for i in mesh),
                kpts_weights=[0 for _ in mesh],
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
        """The input set's POTCAR."""
        if self.structure is None:
            raise RuntimeError("No structure is associated with the input set!")

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

    def estimate_nbands(self) -> int:
        """Estimate the number of bands that VASP will initialize a
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
            VaspInputSet: A new input set with settings (Structure, k-points, incar, etc)
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
                wavecar_files = sorted(glob(str(Path(prev_calc_dir) / (f"{fname}*"))))
                if wavecar_files:
                    if fname == "WFULL":
                        for wavecar_file in wavecar_files:
                            fname = Path(wavecar_file).name
                            fname = fname.split(".")[0]
                            files_to_transfer[fname] = wavecar_file
                    else:
                        files_to_transfer[fname] = str(wavecar_files[-1])

        self.files_to_transfer.update(files_to_transfer)
        return self

    @classmethod
    def from_prev_calc(cls, prev_calc_dir: str, **kwargs) -> Self:
        """Generate a set of VASP input files for static calculations from a
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

    def __str__(self) -> str:
        return type(self).__name__

    def __repr__(self) -> str:
        return type(self).__name__

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
        elif self.incar.get("ENCUT", 0) > 0:
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

    @staticmethod
    def from_directory(directory: str | Path, optional_files: dict | None = None) -> VaspInput:
        """Load a set of VASP inputs from a directory.

        Note that only the standard INCAR, POSCAR, POTCAR and KPOINTS files are read
        unless optional_filenames is specified.

        Parameters
        ----------
        directory
            Directory to read VASP inputs from.
        optional_files
            Optional files to read in as well as a dict of {filename: Object class}.
            Object class must have a static/class method from_file.
        """
        directory = Path(directory)
        objs = {"INCAR": Incar, "KPOINTS": Kpoints, "POSCAR": Poscar, "POTCAR": Potcar}

        inputs = {}
        for name, obj in objs.items():
            if (directory / name).exists():
                inputs[name.upper()] = obj.from_file(directory / name)  # type: ignore[attr-defined]
            else:
                # handle the case where there is no KPOINTS file
                inputs[name.upper()] = None

        optional_inputs = {}
        if optional_files is not None:
            for name, obj in optional_files.items():
                optional_inputs[str(name)] = obj.from_file(directory / name)  # type: ignore[attr-defined]

        return VaspInput(
            incar=inputs["INCAR"],
            kpoints=inputs["KPOINTS"],
            poscar=inputs["POSCAR"],
            potcar=inputs["POTCAR"],
            optional_files=optional_inputs,  # type: ignore[arg-type]
        )

    def _get_nedos(self, dedos: float) -> int:
        """Automatic setting of nedos using the energy range and the energy step."""
        if self.prev_vasprun is None:
            return 2000

        emax = max(eigs.max() for eigs in self.prev_vasprun.eigenvalues.values())
        emin = min(eigs.min() for eigs in self.prev_vasprun.eigenvalues.values())
        return int((emax - emin) / dedos)


# create VaspInputGenerator alias to follow atomate2 terminology
VaspInputGenerator = VaspInputSet


class DictSet(VaspInputSet):
    """Alias for VaspInputSet."""

    def __post_init__(self):
        super().__post_init__()
        warnings.warn(
            "DictSet is deprecated, and will be removed on 2025-12-31. Use VaspInputSet",
            category=FutureWarning,
            stacklevel=2,
        )


# Helper functions to determine valid FFT grids for VASP
def next_num_with_prime_factors(n: int, max_prime_factor: int, must_inc_2: bool = True) -> int:
    """Get the next number greater than or equal to n that only has the desired prime factors.

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


def get_vasprun_outcar(path: str | Path, parse_dos: bool = True, parse_eigen: bool = True) -> tuple[Vasprun, Outcar]:
    """Get a Vasprun and Outcar from a directory.

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
    outcarfile_fullpath = str(path / "OUTCAR.gz")
    vsfile = vsfile_fullpath if vsfile_fullpath in vruns else max(vruns)
    outcarfile = outcarfile_fullpath if outcarfile_fullpath in outcars else max(outcars)
    return Vasprun(vsfile, parse_dos=parse_dos, parse_eigen=parse_eigen), Outcar(outcarfile)


def get_structure_from_prev_run(vasprun, outcar=None) -> Structure:
    """
    Process structure from previous run.

    Args:
        vasprun (Vasprun): Vasprun that contains the final structure from previous run.
        outcar (Outcar): Outcar that contains the magnetization info from previous run.

    Returns:
        Structure: The magmom-decorated structure that can be passed to get VASP input files, e.g.
            get_kpoints().
    """
    structure = vasprun.final_structure

    site_properties = {}
    # magmom
    if vasprun.is_spin:
        if outcar and outcar.magnetization:
            site_properties["magmom"] = [i["tot"] for i in outcar.magnetization]
        else:
            site_properties["magmom"] = vasprun.parameters["MAGMOM"]
    # LDAU
    if vasprun.parameters.get("LDAU", False):
        for key in ("LDAUU", "LDAUJ", "LDAUL"):
            vals = vasprun.incar[key]
            m = {}
            l_val = []
            s = 0
            for site in structure:
                if site.specie.symbol not in m:
                    m[site.specie.symbol] = vals[s]
                    s += 1
                l_val.append(m[site.specie.symbol])
            if len(l_val) == len(structure):
                site_properties.update({key.lower(): l_val})
            else:
                raise ValueError(f"length of list {l_val} not the same as structure")

    return structure.copy(site_properties=site_properties)


def standardize_structure(structure, sym_prec=0.1, international_monoclinic=True):
    """Get the symmetrically standardized structure.

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
    vasp_input_set=None,
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
    if vasp_input_set is None:
        from pymatgen.io.vasp.sets.mp import MPRelaxSet

        vasp_input_set = MPRelaxSet

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


def get_valid_magmom_struct(structure: Structure, inplace: bool = True, spin_mode: str = "auto") -> Structure:
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

    ret_struct = structure if inplace else structure.copy()
    for site in ret_struct:
        if mode == "n":
            if "magmom" in site.properties:
                site.properties.pop("magmom")
        elif "magmom" not in site.properties or site.properties["magmom"] is None:
            site.properties["magmom"] = default_values[mode]

    return ret_struct


def _get_ispin(vasprun: Vasprun | None, outcar: Outcar | None) -> int:
    """Get value of ISPIN depending on the magnetisation in the OUTCAR and vasprun."""
    if outcar is not None and outcar.magnetization is not None:
        # Turn off spin when magmom for every site is smaller than 0.02.
        site_magmom = np.array([i["tot"] for i in outcar.magnetization])
        return 2 if np.any(np.abs(site_magmom) > 0.02) else 1
    if vasprun is not None:
        return 2 if vasprun.is_spin else 1
    return 2


def _get_recommended_lreal(structure: Structure) -> str | bool:
    """Get recommended LREAL flag based on the structure."""
    return "Auto" if structure.num_sites > 16 else False


def _combine_kpoints(*kpoints_objects: Sequence[Kpoints]) -> Kpoints:
    """Combine multiple Kpoints objects."""
    _labels: list[list[str]] = []
    _kpoints: list[Sequence[Kpoint]] = []
    _weights = []

    kpoints_obj: Kpoints
    for kpoints_obj in kpoints_objects:  # type: ignore[assignment]
        if kpoints_obj is None:
            continue
        if kpoints_obj.style != Kpoints.supported_modes.Reciprocal:
            raise ValueError("Can only combine kpoints with style=Kpoints.supported_modes.Reciprocal")
        if kpoints_obj.labels is None:
            _labels.append([""] * len(kpoints_obj.kpts))
        else:
            _labels.append(kpoints_obj.labels)

        _kpoints.append(kpoints_obj.kpts)
        _weights.append(kpoints_obj.kpts_weights)

    labels = np.concatenate(_labels).tolist()
    kpoints = np.concatenate(_kpoints).tolist()
    weights = np.concatenate(_weights).tolist()

    return Kpoints(
        comment="Combined k-points",
        style=Kpoints.supported_modes.Reciprocal,
        num_kpts=len(kpoints),
        kpts=cast(Sequence[Kpoint], kpoints),
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
    """Get NEDOS using the energy range and the energy step,
    defaults to 2000.
    """
    if vasprun is None:
        return 2000

    if vasprun.eigenvalues is None:
        raise RuntimeError("eigenvalues cannot be None.")

    emax = max(eigs.max() for eigs in vasprun.eigenvalues.values())
    emin = min(eigs.min() for eigs in vasprun.eigenvalues.values())
    return int((emax - emin) / dedos)


def auto_kspacing(bandgap, bandgap_tol):
    """Set kspacing based on the bandgap."""
    if bandgap is None or bandgap <= bandgap_tol:  # metallic
        return 0.22

    rmin = max(1.5, 25.22 - 2.87 * bandgap)  # Eq. 25
    kspacing = 2 * np.pi * 1.0265 / (rmin - 1.0183)  # Eq. 29

    # cap kspacing at a max of 0.44, per internal benchmarking
    return min(kspacing, 0.44)
