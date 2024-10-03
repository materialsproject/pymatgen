"""
This module defines input sets for CP2K and is a work in progress. The structure/philosophy
of this module is based on the Vasp input sets in Pymatgen. These sets are meant to contain
tested parameters that will result in successful, reproducible, consistent calculations without
need for intervention 99% of the time. 99% of the time, you only need to provide a pymatgen
structure object and let the defaults take over from there.

The sets are intended to be very general, e.g. a set for geometry relaxation, and so most of the
time, if you have specific needs, you can simply specify them via the keyword argument
override_default_params (see Section.update() method). If you have the need to create a new input
set (say for a standardized high throughput calculation) then you can create a new child of the
Cp2kInputSet class.

In order to implement a new Set within the current code structure, follow this 3 step flow:
    (1) Inherit from Cp2kInputSet or one of its children and call the super() constructor
    (2) Create the new sections and insert them into self and its subsections as needed
    (3) Call self.update(override_default_params) in order to allow user settings.
"""

from __future__ import annotations

import itertools
import os
import warnings
from typing import TYPE_CHECKING

import numpy as np
from ruamel.yaml import YAML

from pymatgen.core import SETTINGS
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Element, Molecule, Structure
from pymatgen.io.cp2k.inputs import (
    DOS,
    LDOS,
    PBE,
    PDOS,
    QS,
    BandStructure,
    BasisFile,
    BasisInfo,
    BrokenSymmetry,
    Cell,
    Coord,
    Cp2kInput,
    Dft,
    EDensityCube,
    ForceEval,
    GaussianTypeOrbitalBasisSet,
    Global,
    GthPotential,
    Keyword,
    Kind,
    Kpoints,
    Mgrid,
    MOCubes,
    OrbitalTransformation,
    PotentialFile,
    PotentialInfo,
    Scf,
    Section,
    SectionList,
    Smear,
    Subsys,
    VHartreeCube,
    XCFunctional,
)
from pymatgen.io.cp2k.utils import get_truncated_coulomb_cutoff, get_unique_site_indices
from pymatgen.io.vasp.inputs import Kpoints as VaspKpoints
from pymatgen.io.vasp.inputs import KpointsSupportedModes

if TYPE_CHECKING:
    from typing import Literal

__author__ = "Nicholas Winner"
__version__ = "2.0"
__email__ = "nwinner@berkeley.edu"
__date__ = "September 2022"


class DftSet(Cp2kInput):
    """
    Base for an input set using the Quickstep module (i.e. a DFT calculation). The DFT section is
    pretty vast in CP2K, so this set hopes to make the DFT setup fairly simple. The provided
    parameters are pretty conservative, and so they should not need to be changed very often.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        project_name: str = "CP2K",
        basis_and_potential: dict | None = None,
        xc_functionals: list | str | None = None,
        multiplicity: int = 0,
        ot: bool = True,
        energy_gap: float = -1,
        qs_method: str = "GPW",
        eps_default: float = 1e-12,
        eps_scf: float = 1e-6,
        max_scf: int | None = None,
        minimizer: str = "DIIS",
        preconditioner: str = "FULL_SINGLE_INVERSE",
        algorithm: str = "IRAC",
        linesearch: str = "2PNT",
        rotation: bool = True,
        occupation_preconditioner: bool = False,
        cutoff: float | None = None,
        rel_cutoff: int = 50,
        ngrids: int = 5,
        progression_factor: int = 3,
        override_default_params: dict | None = None,
        wfn_restart_file_name: str | None = None,
        kpoints: VaspKpoints | None = None,
        smearing: bool = False,
        **kwargs,
    ) -> None:
        """
        Args:
            structure: Pymatgen structure or molecule object
            ot (bool): Whether or not to use orbital transformation method for matrix
                diagonalization. OT is the flagship scf solver of CP2K, and will provide
                speed-ups for this part of the calculation, but the system must have a band gap
                for OT to be used (higher band-gap --> faster convergence).
            energy_gap (float): Estimate of energy gap for pre-conditioner. Default is -1, leaving
                it up to CP2K.
            eps_default (float): Replaces all EPS_XX Keywords in the DFT section value, ensuring
                an overall accuracy of at least this much.
            eps_scf (float): The convergence criteria for leaving the SCF loop. Default is 1e-6.
                Should ensure reasonable results, but is not applicable to all situations.
                    Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
                    it is the largest gradient of the energy with respect to changing any of the
                    molecular orbital coefficients. For diagonalization, it is the largest change
                    in the density matrix from the last step.
            max_scf (int): The max number of SCF cycles before terminating the solver. NOTE: With
                the OT solver, this corresponds to the max number of INNER scf loops, and then
                the outer loops are set with outer_max_scf, while with diagonalization it
                corresponds to the overall (INNER*OUTER) number of SCF steps, with the
                inner loop limit set by
            minimizer (str): The minimization scheme. DIIS can be as much as 50% faster than the
                more robust conjugate gradient method, and so it is chosen as default. Switch to CG
                if dealing with a difficult system.
            preconditioner (str): Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
                robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
                considered "best" but needs algorithm to be set to STRICT. Only change from these
                two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
                preferred.
            algorithm (str): Algorithm for the OT method. STRICT assumes that the orbitals are
                strictly orthogonal to each other, which works well for wide gap ionic systems,
                but can diverge for systems with small gaps, fractional occupations, and some
                other cases. IRAC (iterative refinement of the approximate congruency)
                transformation is not analytically correct and uses a truncated polynomial
                expansion, but is robust to the problems with STRICT, and so is the default.
            linesearch (str): Linesearch method for CG. 2PNT is the default, and is the fastest,
                but is not as robust as 3PNT. 2PNT is required as of CP2K v9.1 for compatibility
                with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
                gapped transition metal systems as an alternative.
            rotation (bool): Whether or not to allow for rotation of the orbitals in the OT method.
                This equates to allowing for fractional occupations in the calculation.
            occupation_preconditioner (bool): Whether or not to account for fractional occupations
                in the preconditioner. This method is not fully integrated as of CP2K v9.1 and is
                set to false by default.
            cutoff (int): Cutoff energy (in Ry) for the finest level of the multigrid. A high
                cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
                is appropriate. By default cutoff is set to 0, leaving it up to the set.
            rel_cutoff (int): This cutoff decides how the Gaussians are mapped onto the different
                levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
                CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
                and thus the effective integration grid for the calculation may still be too
                coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
                large enough.
            ngrids (int): number of multi-grids to use. CP2K default is 4, but the molopt basis
                files recommend 5.
            progression_factor (int): Divisor of CUTOFF to get the cutoff for the next level of
                the multigrid.
            wfn_restart_file_name (str): RESTART file for the initial wavefunction guess.
            kpoints (Kpoints): kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
                CP2K runs with gamma point only.
            smearing (bool): whether or not to activate smearing (should be done for systems
                containing no (or a very small) band gap.
        """
        super().__init__(name="CP2K_INPUT", subsections={})

        self.structure = structure
        self.basis_and_potential = basis_and_potential or {}
        self.project_name = project_name
        self.charge = int(structure.charge)
        if not multiplicity and isinstance(self.structure, Molecule):
            self.multiplicity = self.structure.spin_multiplicity
        else:
            self.multiplicity = multiplicity
        self.ot = ot
        self.qs_method = qs_method
        self.energy_gap = energy_gap
        self.eps_default = eps_default
        self.eps_scf = eps_scf
        self.max_scf = max_scf
        self.minimizer = minimizer
        self.preconditioner = preconditioner
        self.algorithm = algorithm
        self.rotation = rotation
        self.occupation_preconditioner = occupation_preconditioner
        self.linesearch = linesearch
        self.cutoff = cutoff
        self.rel_cutoff = rel_cutoff
        self.ngrids = ngrids
        self.progression_factor = progression_factor
        self.override_default_params = override_default_params or {}
        self.wfn_restart_file_name = wfn_restart_file_name
        self.kpoints = kpoints
        self.smearing = smearing
        self.kwargs = kwargs

        # Enable force and energy evaluations (most calculations)
        self.insert(ForceEval())

        if self.kpoints:
            # As of CP2K v2022.1 kpoint module is not fully integrated, so even specifying
            # 0,0,0 will disable certain features. So, you have to drop it all together to
            # get full support
            if (
                self.kpoints.style in [KpointsSupportedModes.Gamma, KpointsSupportedModes.Monkhorst]
                and np.array_equal(self.kpoints.kpts[0], (1, 1, 1))
            ) or (
                self.kpoints.style in [KpointsSupportedModes.Reciprocal, KpointsSupportedModes.Cartesian]
                and np.array_equal(self.kpoints.kpts[0], (0, 0, 0))
            ):
                self.kpoints = None
            if ot and self.kpoints:
                warnings.warn("As of 2022.1, kpoints not supported with OT. Defaulting to diagonalization")
                ot = False

        # Build the global section
        global_sec = Global(
            project_name=self.kwargs.get("project_name", "CP2K"),
            run_type=self.kwargs.get("run_type", "ENERGY_FORCE"),
        )
        self.insert(global_sec)

        # Build the QS Section
        qs = QS(method=self.qs_method, eps_default=eps_default, eps_pgf_orb=kwargs.get("eps_pgf_orb"))
        max_scf = max_scf or 20 if ot else 400  # If ot, max_scf is for inner loop
        scf = Scf(eps_scf=eps_scf, max_scf=max_scf, subsections={})

        if ot:
            scf.insert(
                OrbitalTransformation(
                    minimizer=minimizer,
                    preconditioner=preconditioner,
                    energy_gap=energy_gap,
                    algorithm=algorithm,
                    linesearch=linesearch,
                    rotation=rotation,
                    occupation_preconditioner=occupation_preconditioner,
                )
            )
            scf.insert(
                Section(
                    "OUTER_SCF",
                    subsections={},
                    keywords={
                        "MAX_SCF": Keyword("MAX_SCF", kwargs.get("outer_max_scf", 20)),
                        "EPS_SCF": Keyword("EPS_SCF", kwargs.get("outer_eps_scf", eps_scf)),
                    },
                )
            )
        else:
            scf.insert(Section("DIAGONALIZATION", subsections={}))
            mixing_kwds = {
                "ALPHA": Keyword("ALPHA", kwargs.get("alpha", 0.05)),
                "BETA": Keyword("BETA", kwargs.get("beta", 0.01)),
                "NBUFFER": Keyword("NBUFFER", kwargs.get("nbuffer", 10)),
                "N_SIMPLE_MIX": Keyword("N_SIMPLE_MIX", kwargs.get("n_simple_mix", 3)),
                "METHOD": Keyword("METHOD", kwargs.get("mixing_method", "BROYDEN_MIXING")),
            }
            mixing = Section("MIXING", keywords=mixing_kwds, subsections=None)
            scf.insert(mixing)
            scf["MAX_DIIS"] = Keyword("MAX_DIIS", 15)

        # Get basis, potential, and xc info
        self.basis_and_potential = DftSet.get_basis_and_potential(self.structure, self.basis_and_potential)
        self.basis_set_file_names = self.basis_and_potential.get("basis_filenames")
        self.potential_file_name = self.basis_and_potential.get("potential_filename")
        self.xc_functionals = DftSet.get_xc_functionals(
            xc_functionals=xc_functionals
        )  # kwargs.get("xc_functional", "PBE"))

        # create the subsys (structure)
        self.create_subsys(self.structure)

        # Create the multigrid for FFTs
        if not self.cutoff:
            basis_sets = [self.basis_and_potential[el].get("basis") for el in self.structure.symbol_set]
            self.cutoff = DftSet.get_cutoff_from_basis(basis_sets=basis_sets, rel_cutoff=rel_cutoff)

        mgrid = Mgrid(
            cutoff=self.cutoff,
            rel_cutoff=rel_cutoff,
            ngrids=ngrids,
            progression_factor=progression_factor,
        )

        if smearing and not ot:
            scf["ADDED_MOS"] = Keyword("ADDED_MOS", -1, -1) if self.kwargs.get("spin_polarized", True) else -1
            scf.insert(Smear(elec_temp=kwargs.get("elec_temp", 300)))

        # Set the DFT calculation with global parameters
        dft = Dft(
            MULTIPLICITY=self.multiplicity,
            CHARGE=self.charge,
            uks=self.kwargs.get("spin_polarized", True),
            basis_set_filenames=self.basis_set_file_names or [],
            potential_filename=self.potential_file_name,
            subsections={"QS": qs, "SCF": scf, "MGRID": mgrid},
            wfn_restart_file_name=wfn_restart_file_name,
        )

        # Set kpoints only if user supplies them
        if self.kpoints:
            dft.insert(Kpoints.from_kpoints(self.kpoints, structure=self.structure))

        # Create subsections and insert into them
        self["FORCE_EVAL"].insert(dft)

        xc_functional = XCFunctional(functionals=self.xc_functionals)
        xc = Section("XC", subsections={"XC_FUNCTIONAL": xc_functional})
        self["FORCE_EVAL"]["DFT"].insert(xc)
        self["FORCE_EVAL"]["DFT"].insert(Section("PRINT", subsections={}))

        if isinstance(structure, Molecule):
            self.activate_nonperiodic()

        if kwargs.get("print_forces", True):
            self.print_forces()
        if kwargs.get("print_dos", True):
            self.print_dos()
        if kwargs.get("print_pdos", True):
            self.print_pdos()
        if kwargs.get("print_ldos", False):
            self.print_ldos()
        if kwargs.get("print_mo_cubes", True):
            self.print_mo_cubes()
        if kwargs.get("print_v_hartree", True):
            self.print_v_hartree()
        if kwargs.get("print_e_density", True):
            self.print_e_density()
        if kwargs.get("print_bandstructure", False):
            self.print_bandstructure(kwargs.get("kpoints_line_density", 20))

        self.update(self.override_default_params)
        if kwargs.get("validate", True):
            self.validate()

    @staticmethod
    def get_basis_and_potential(structure, basis_and_potential):
        """Get a dictionary of basis and potential info for constructing the input file.

        data in basis_and_potential argument can be specified in several ways:

            Strategy 1: Element-specific info (takes precedence)

                1. Provide a basis and potential object:

                    el: {'basis': obj, 'potential': obj}

                2. Provide a hash of the object that matches the keys in the pmg configured CP2K data files.

                    el: {'basis': hash, 'potential': hash}

                3. Provide the name of the basis and potential AND the basis_filenames and potential_filename
                keywords specifying where to find these objects

                    el: {
                        'basis': name, 'potential': name, 'basis_filenames': [filenames],
                        'potential_filename': filename
                    }

            Strategy 2: global descriptors

                In this case, any elements not present in the argument will be dealt with by searching the pmg
                configured CP2K data files to find a objects matching your requirements.

                    - functional: Find potential and basis that have been optimized for a specific functional like PBE.
                        Can be None if you do not require them to match.
                    - basis_type: type of basis to search for (e.g. DZVP-MOLOPT).
                    - aux_basis_type: type of basis to search for (e.g. pFIT). Some elements do not have all aux types
                        available. Use aux_basis_type that is universal to avoid issues, or avoid using this strategy.
                    - potential_type: "Pseudopotential" or "All Electron"

                ***BE WARNED*** CP2K data objects can have the same name, this will sort those and choose the first one
                that matches.

        Will raise an error if no basis/potential info can be found according to the input.
        """
        data = {"basis_filenames": []}
        functional = basis_and_potential.get("functional", SETTINGS.get("PMG_DEFAULT_CP2K_FUNCTIONAL"))
        basis_type = basis_and_potential.get("basis_type", SETTINGS.get("PMG_DEFAULT_CP2K_BASIS_TYPE"))
        potential_type = basis_and_potential.get(
            "potential_type", SETTINGS.get("PMG_DEFAULT_POTENTIAL_TYPE", "Pseudopotential")
        )
        aux_basis_type = basis_and_potential.get("aux_basis_type", SETTINGS.get("PMG_DEFAULT_CP2K_AUX_BASIS_TYPE"))

        for el in structure.symbol_set:
            possible_basis_sets = []
            possible_potentials = []
            basis, aux_basis, potential, DATA = None, None, None, None
            desired_basis, desired_aux_basis, desired_potential = None, None, None
            have_element_file = os.path.isfile(os.path.join(SETTINGS.get("PMG_CP2K_DATA_DIR", "."), el))

            # Necessary if matching data to CP2K data files
            if have_element_file:
                with open(os.path.join(SETTINGS.get("PMG_CP2K_DATA_DIR", "."), el), encoding="utf-8") as file:
                    yaml = YAML(typ="unsafe", pure=True)
                    DATA = yaml.load(file)
                    if not DATA.get("basis_sets"):
                        raise ValueError(f"No standard basis sets available in data directory for {el}")
                    if not DATA.get("potentials"):
                        raise ValueError(f"No standard potentials available in data directory for {el}")

            # Give precedence to explicitly listing info for the element
            if el in basis_and_potential:
                _basis = basis_and_potential[el].get("basis")
                if isinstance(_basis, GaussianTypeOrbitalBasisSet):
                    possible_basis_sets.append(_basis)
                elif have_element_file:
                    if _basis in DATA["basis_sets"]:
                        possible_basis_sets.append(GaussianTypeOrbitalBasisSet.from_dict(DATA["basis_sets"][_basis]))
                    elif _basis:
                        desired_basis = GaussianTypeOrbitalBasisSet(name=_basis)

                _aux_basis = basis_and_potential[el].get("aux_basis")
                if isinstance(_aux_basis, GaussianTypeOrbitalBasisSet):
                    aux_basis = _aux_basis
                elif have_element_file:
                    if _aux_basis in DATA["basis_sets"]:
                        aux_basis = GaussianTypeOrbitalBasisSet.from_dict(DATA["basis_sets"][_aux_basis])
                    elif _aux_basis:
                        desired_aux_basis = GaussianTypeOrbitalBasisSet(name=_aux_basis)

                _potential = basis_and_potential[el].get("potential")
                if isinstance(_potential, GthPotential):
                    possible_potentials.append(_potential)
                elif have_element_file:
                    if _potential in DATA["potentials"]:
                        possible_potentials.append(GthPotential.from_dict(DATA["potentials"][_potential]))
                    elif _potential:
                        desired_potential = GthPotential(name=_potential)

                if basis_and_potential[el].get("basis_filename"):
                    data["basis_filenames"].append(basis_and_potential[el].get("basis_filename"))
                pfn1 = basis_and_potential[el].get("potential_filename")
                pfn2 = data.get("potential_filename")
                if pfn1 and pfn2 and pfn1 != pfn2:
                    raise ValueError(
                        "Provided potentials have more than one corresponding file."
                        "CP2K does not support multiple potential filenames"
                    )
                data["potential_filename"] = basis_and_potential[el].get("potential_filename")

            # Else, if functional/basis settings are provided, create objects to search the data files
            else:
                if basis_type and have_element_file:
                    desired_basis = GaussianTypeOrbitalBasisSet(
                        element=Element(el),
                        potential=potential_type,
                        info=BasisInfo.from_str(f"{basis_type}-{functional}"),
                    )
                    desired_potential = GthPotential(
                        element=Element(el), potential=potential_type, info=PotentialInfo(xc=functional)
                    )
                if aux_basis_type and have_element_file:
                    desired_aux_basis = GaussianTypeOrbitalBasisSet(info=BasisInfo.from_str(aux_basis_type))

            # If basis/potential are not explicit, match the desired ones to available ones in the element file
            if desired_basis:
                for _possible_basis in DATA.get("basis_sets").values():
                    possible_basis = GaussianTypeOrbitalBasisSet.from_dict(_possible_basis)
                    if desired_basis.softmatch(possible_basis):
                        possible_basis_sets.append(possible_basis)
            if desired_aux_basis:
                for _possible_basis in DATA.get("basis_sets").values():
                    possible_basis = GaussianTypeOrbitalBasisSet.from_dict(_possible_basis)
                    if desired_aux_basis.softmatch(possible_basis):
                        aux_basis = possible_basis
                        data["basis_filenames"].append(aux_basis.filename)
                        break
            if desired_potential:
                for _possible_potential in DATA.get("potentials").values():
                    possible_potential = GthPotential.from_dict(_possible_potential)
                    if desired_potential.softmatch(possible_potential):
                        possible_potentials.append(possible_potential)

            possible_basis_sets = sorted(
                filter(lambda x: x.info.electrons, possible_basis_sets), key=lambda x: x.info.electrons, reverse=True
            )
            possible_potentials = sorted(
                filter(lambda x: x.info.electrons, possible_potentials), key=lambda x: x.info.electrons, reverse=True
            )

            def match_elecs(x):
                for p in possible_potentials:
                    if x.info.electrons == p.info.electrons:
                        return p
                return None

            for b in possible_basis_sets:
                fb = match_elecs(b)
                if fb is not None:
                    basis = b
                    potential = fb
                    break

            if basis is None:
                if not basis_and_potential.get(el, {}).get("basis"):
                    raise ValueError(f"No explicit basis found for {el} and matching has failed.")
                warnings.warn(f"Unable to validate basis for {el}. Exact name provided will be put in input file.")
                basis = basis_and_potential[el].get("basis")

            if aux_basis is None and basis_and_potential.get(el, {}).get("aux_basis"):
                warnings.warn(
                    f"Unable to validate auxiliary basis for {el}. Exact name provided will be put in input file."
                )
                aux_basis = basis_and_potential[el].get("aux_basis")

            if potential is None:
                if basis_and_potential.get(el, {}).get("potential"):
                    warnings.warn(
                        f"Unable to validate potential for {el}. Exact name provided will be put in input file."
                    )
                    potential = basis_and_potential.get(el, {}).get("potential")
                else:
                    raise ValueError("No explicit potential found and matching has failed.")

            if hasattr(basis, "filename"):
                data["basis_filenames"].append(basis.filename)
            pfn1 = data.get("potential_filename")
            pfn2 = potential.filename
            if pfn1 and pfn2 and pfn1 != pfn2:
                raise ValueError(
                    "Provided potentials have more than one corresponding file."
                    "CP2K does not support multiple potential filenames"
                )
            data["potential_filename"] = pfn2

            data[el] = {"basis": basis, "aux_basis": aux_basis, "potential": potential}
        return data

    @staticmethod
    def get_cutoff_from_basis(basis_sets, rel_cutoff) -> float:
        """Given a basis and a relative cutoff. Determine the ideal cutoff variable."""
        for basis in basis_sets:
            if not basis.exponents:
                raise ValueError(f"Basis set {basis} contains missing exponent info. Please specify cutoff manually")

        def get_soft_exponents(b):
            if b.potential == "All Electron":
                radius = 1.2 if b.element == Element("H") else 1.512  # Hard radius defaults for gapw
                threshold = 1e-4  # Default for gapw
                max_lshell = max(shell for shell in b.lmax)
                exponent = np.log(radius**max_lshell / threshold) / radius**2
                return [[exponent]]
            return b.exponents

        exponents = [get_soft_exponents(b) for b in basis_sets if b.exponents]
        exponents = list(itertools.chain.from_iterable(exponents))
        return np.ceil(max(itertools.chain.from_iterable(exponents))) * rel_cutoff

    @staticmethod
    def get_xc_functionals(xc_functionals: list | str | None = None) -> list:
        """Get XC functionals. If simplified names are provided in kwargs, they
        will be expanded into their corresponding X and C names.
        """
        names = xc_functionals or SETTINGS.get("PMG_DEFAULT_CP2K_FUNCTIONAL")
        if not names:
            raise ValueError(
                "No XC functional provided. Specify kwarg xc_functional or configure PMG_DEFAULT_FUNCTIONAL "
                "in your .pmgrc.yaml file"
            )
        if isinstance(names, str):
            names = [names]
        names = [n.upper() for n in names]
        cp2k_names = []
        for name in names:
            if name in ["LDA", "LSDA"]:
                cp2k_names.append("PADE")
            elif name == "SCAN":
                cp2k_names.extend(["MGGA_X_SCAN", "MGGA_C_SCAN"])
            elif name == "SCANL":
                cp2k_names.extend(["MGGA_X_SCANL", "MGGA_C_SCANL"])
            elif name == "R2SCAN":
                cp2k_names.extend(["MGGA_X_R2SCAN", "MGGA_C_R2SCAN"])
            elif name == "R2SCANL":
                cp2k_names.extend(["MGGA_X_R2SCANL", "MGGA_C_R2SCANL"])
            else:
                cp2k_names.append(name)

        return cp2k_names

    def write_basis_set_file(self, basis_sets, fn="BASIS") -> None:
        """Write the basis sets to a file."""
        BasisFile(objects=basis_sets).write_file(fn)

    def write_potential_file(self, potentials, fn="POTENTIAL") -> None:
        """Write the potentials to a file."""
        PotentialFile(objects=potentials).write_file(fn)

    def print_forces(self) -> None:
        """Print out the forces and stress during calculation."""
        self["FORCE_EVAL"].insert(Section("PRINT", subsections={}))
        self["FORCE_EVAL"]["PRINT"].insert(Section("FORCES", subsections={}))
        self["FORCE_EVAL"]["PRINT"].insert(Section("STRESS_TENSOR", subsections={}))

    def print_dos(self, ndigits=6) -> None:
        """
        Activate printing of the overall DOS file.

        Note: As of 2022.1, ndigits needs to be set to a sufficient value to ensure data is not lost.
        Note: As of 2022.1, can only be used with a k-point calculation.
        """
        if self.kpoints:
            self["force_eval"]["dft"]["print"].insert(DOS(ndigits))

    def print_pdos(self, nlumo: int = -1) -> None:
        """
        Activate creation of the PDOS file.

        Args:
            nlumo (int): Number of virtual orbitals to be added to the MO set (-1=all).
                CAUTION: Setting this value to be higher than the number of states present may
                cause a Cholesky error.
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/PDOS") and not self.kpoints:
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(PDOS(nlumo=nlumo))

    def print_ldos(self, nlumo: int = -1) -> None:
        """
        Activate the printing of LDOS files, printing one for each atom kind by default.

        Args:
            nlumo (int): Number of virtual orbitals to be added to the MO set (-1=all).
                CAUTION: Setting this value to be higher than the number of states present may
                cause a Cholesky error.
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/PDOS"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(PDOS(nlumo=nlumo))
        for idx in range(len(self.structure)):
            self["FORCE_EVAL"]["DFT"]["PRINT"]["PDOS"].insert(LDOS(idx + 1, alias=f"LDOS {idx + 1}", verbose=False))

    def print_mo_cubes(self, write_cube: bool = False, nlumo: int = -1, nhomo: int = -1) -> None:
        """
        Activate printing of molecular orbitals.

        Args:
            write_cube (bool): whether to write cube file for the MOs instead of out file
            nlumo (int): Controls the number of lumos printed and dumped as a cube (-1=all)
            nhomo (int): Controls the number of homos printed and dumped as a cube (-1=all)
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/MO_CUBES"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(MOCubes(write_cube=write_cube, nlumo=nlumo, nhomo=nhomo))

    def print_mo(self) -> None:
        """Print molecular orbitals when running non-OT diagonalization."""
        raise NotImplementedError

    def print_v_hartree(self, stride=(2, 2, 2)) -> None:
        """
        Controls the printing of a cube file with eletrostatic potential generated by the
            total density (electrons+ions). It is valid only for QS with GPW formalism.
        Note that by convention the potential has opposite sign than the expected physical one.
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(VHartreeCube(keywords={"STRIDE": Keyword("STRIDE", *stride)}))

    def print_e_density(self, stride=(2, 2, 2)) -> None:
        """Controls the printing of cube files with electronic density and, for UKS, the spin density."""
        if not self.check("FORCE_EVAL/DFT/PRINT/E_DENSITY_CUBE"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(EDensityCube(keywords={"STRIDE": Keyword("STRIDE", *stride)}))

    def print_bandstructure(self, kpoints_line_density: int = 20) -> None:
        """
        Attaches a non-scf band structure calc the end of an SCF loop.

        This requires a kpoint calculation, which is not always default in CP2K.

        Args:
            kpoints_line_density: number of kpoints along each branch in line-mode calc.
        """
        if not self.kpoints:
            raise ValueError("Kpoints must be provided to enable band structure printing")

        bs = BandStructure.from_kpoints(
            self.kpoints,
            kpoints_line_density=kpoints_line_density,
        )
        self["force_eval"]["dft"]["print"].insert(bs)

    def print_hirshfeld(self, on=True) -> None:
        """Activate or deactivate printing of Hirshfeld charges."""
        section = Section("HIRSHFELD", section_parameters=["ON" if on else "OFF"])
        self["force_eval"]["dft"]["print"].insert(section)

    def print_mulliken(self, on=False) -> None:
        """Activate or deactivate printing of Mulliken charges."""
        section = Section("MULLIKEN", section_parameters=["ON" if on else "OFF"])
        self["force_eval"]["dft"]["print"].insert(section)

    def set_charge(self, charge: int) -> None:
        """Set the overall charge of the simulation cell."""
        self["FORCE_EVAL"]["DFT"]["CHARGE"] = Keyword("CHARGE", int(charge))

    def activate_hybrid(
        self,
        hybrid_functional: str = "PBE0",
        hf_fraction: float = 0.25,
        gga_x_fraction: float = 0.75,
        gga_c_fraction: float = 1,
        max_memory: int = 2000,
        cutoff_radius: float = 8.0,
        potential_type: str | None = None,
        omega: float = 0.11,
        scale_coulomb: float = 1,
        scale_gaussian: float = 1,
        scale_longrange: float = 1,
        admm: bool = True,
        admm_method: str = "BASIS_PROJECTION",
        admm_purification_method: str = "NONE",
        admm_exch_correction_func: str = "DEFAULT",
        eps_schwarz: float = 1e-7,
        eps_schwarz_forces: float = 1e-6,
        screen_on_initial_p: bool = True,
        screen_p_forces: bool = True,
    ) -> None:
        """
        Basic set for activating hybrid DFT calculation using Auxiliary Density Matrix Method.

        Note 1: When running ADMM with CP2K, memory is very important. If the memory requirements
        exceed what is available (see max_memory), then CP2K will have to calculate the 4-electron
        integrals for HFX during each step of the SCF cycle. ADMM provides a huge speed up by
        making the memory requirements *feasible* to fit into RAM, which means you only need to
        calculate the integrals once each SCF cycle. But, this only works if it fits into memory.
        When setting up ADMM calculations, we recommend doing whatever is possible to fit all the
        4EI into memory.

        Note 2: This set is designed for reliable high-throughput calculations, NOT for extreme
        accuracy. Please review the in-line comments in this method if you want more control.

        Args:
            hybrid_functional (str): Type of hybrid functional. This set supports HSE (screened)
                and PBE0 (truncated). Default is PBE0, which converges easier in the GPW basis
                used by CP2K.
            hf_fraction (float): fraction of exact HF exchange energy to mix. Default: 0.25
            gga_x_fraction (float): fraction of gga exchange energy to retain. Default: 0.75
            gga_c_fraction (float): fraction of gga correlation energy to retain. Default: 1.0
            max_memory (int): Maximum memory available to each MPI process (in Mb) in the
                calculation. Most modern computing nodes will have ~2Gb per core, or 2048 Mb,
                but check for your specific system. This value should be as large as possible
                while still leaving some memory for the other parts of CP2K. Important: If
                this value is set larger than the memory limits, CP2K will likely seg-fault.
                Default: 2000
            cutoff_radius (float): for truncated hybrid functional (i.e. PBE0), this is the cutoff
                radius. The default is selected as that which generally gives convergence, but
                maybe too low (if you want very high accuracy) or too high (if you want a quick
                screening). Default: 8 angstroms
            potential_type (str): what interaction potential to use for HFX. Available in CP2K are
                COULOMB, GAUSSIAN, IDENTITY, LOGRANGE, MIX_CL, MIX_CL_TRUNC, MIX_LG, SHORTRANGE,
                and TRUNCATED. Default is None, and it will be set automatically depending on the
                named hybrid_functional that you use, but setting it to one of the acceptable
                values will constitute a user-override.
            omega (float): For HSE, this specifies the screening parameter. HSE06 sets this as
                0.2, which is the default.
            scale_coulomb: Scale for the coulomb operator if using a range separated functional
            scale_gaussian: Scale for the gaussian operator (if applicable)
            scale_longrange: Scale for the coulomb operator if using a range separated functional
            admm: Whether or not to use the auxiliary density matrix method for the exact
                HF exchange contribution. Highly recommended. Speed ups between 10x and 1000x are
                possible when compared to non ADMM hybrid calculations.
            admm_method: Method for constructing the auxiliary basis
            admm_purification_method: Method for purifying the auxiliary density matrix so as to
                preserve properties, such as idempotency. May lead to shifts in the
                eigenvalues.
            admm_exch_correction_func: Which functional to use to calculate the exchange correction
                E_x(primary) - E_x(aux)
            eps_schwarz: Screening threshold for HFX, in Ha. Contributions smaller than
                this will be screened. The smaller the value, the more accurate, but also the more
                costly. Default value is 1e-7. 1e-6 works in a large number of cases, but is
                quite aggressive, which can lead to convergence issues.
            eps_schwarz_forces: Same as for eps_schwarz, but for screening contributions to
                forces. Convergence is not as sensitive with respect to eps_schwarz forces as
                compared to eps_schwarz, and so 1e-6 should be good default.
            screen_on_initial_p: If an initial density matrix is provided, in the form of a
                CP2K wfn restart file, then this initial density will be used for screening. This
                is generally very computationally efficient, but, as with eps_schwarz, can lead to
                instabilities if the initial density matrix is poor.
            screen_p_forces: Same as screen_on_initial_p, but for screening of forces.
        """
        if not admm:
            for k, v in self.basis_and_potential.items():
                if "aux_basis" in v:
                    del self.basis_and_potential[k]["aux_basis"]
            del self["force_eval"]["subsys"]
            self.create_subsys(self.structure)

        else:
            aux_matrix_params = {
                "ADMM_PURIFICATION_METHOD": Keyword("ADMM_PURIFICATION_METHOD", admm_purification_method),
                "METHOD": Keyword("METHOD", admm_method),
                "EXCH_CORRECTION_FUNC": Keyword("EXCH_CORRECTION_FUNC", admm_exch_correction_func),
            }
            aux_matrix = Section(
                "AUXILIARY_DENSITY_MATRIX_METHOD",
                keywords=aux_matrix_params,
                subsections={},
            )
            self.subsections["FORCE_EVAL"]["DFT"].insert(aux_matrix)

        # Define the GGA functional as PBE
        screening = Section(
            "SCREENING",
            subsections={},
            keywords={
                "EPS_SCHWARZ": Keyword("EPS_SCHWARZ", eps_schwarz),
                "EPS_SCHWARZ_FORCES": Keyword("EPS_SCHWARZ_FORCES", eps_schwarz_forces),
                "SCREEN_ON_INITIAL_P": Keyword("SCREEN_ON_INITIAL_P", screen_on_initial_p),
                "SCREEN_P_FORCES": Keyword("SCREEN_P_FORCES", screen_p_forces),
            },
        )

        # Check for unphysical cutoff radius
        if isinstance(self.structure, Structure):
            max_cutoff_radius = get_truncated_coulomb_cutoff(self.structure)
            if max_cutoff_radius < cutoff_radius:
                warnings.warn(
                    "Provided cutoff radius exceeds half the minimum"
                    " distance between atoms. I hope you know what you're doing."
                )

        ip_keywords: dict[str, Keyword] = {}
        if hybrid_functional == "HSE06":
            pbe = PBE("ORIG", scale_c=1, scale_x=0)
            xc_functional = XCFunctional(functionals=[], subsections={"PBE": pbe})

            potential_type = potential_type or "SHORTRANGE"
            xc_functional.insert(
                Section(
                    "XWPBE",
                    subsections={},
                    keywords={
                        "SCALE_X0": Keyword("SCALE_X0", 1),
                        "SCALE_X": Keyword("SCALE_X", -0.25),
                        "OMEGA": Keyword("OMEGA", 0.11),
                    },
                )
            )
            ip_keywords |= {
                "POTENTIAL_TYPE": Keyword("POTENTIAL_TYPE", potential_type),
                "OMEGA": Keyword("OMEGA", 0.11),
                "CUTOFF_RADIUS": Keyword("CUTOFF_RADIUS", cutoff_radius),
            }
        elif hybrid_functional == "PBE0":
            pbe = PBE("ORIG", scale_c=1, scale_x=1 - hf_fraction)
            xc_functional = XCFunctional(functionals=[], subsections={"PBE": pbe})

            if isinstance(self.structure, Molecule):
                potential_type = "COULOMB"
            else:
                potential_type = "TRUNCATED"
                xc_functional.insert(
                    Section(
                        "PBE_HOLE_T_C_LR",
                        subsections={},
                        keywords={
                            "CUTOFF_RADIUS": Keyword("CUTOFF_RADIUS", cutoff_radius),
                            "SCALE_X": Keyword("SCALE_X", hf_fraction),
                        },
                    )
                )
                ip_keywords["CUTOFF_RADIUS"] = Keyword("CUTOFF_RADIUS", cutoff_radius)
                ip_keywords["T_C_G_DATA"] = Keyword("T_C_G_DATA", "t_c_g.dat")

            ip_keywords["POTENTIAL_TYPE"] = Keyword("POTENTIAL_TYPE", potential_type)

        elif hybrid_functional == "RSH":
            # Activates range separated functional using mixing of the truncated
            # coulomb operator and the long range operator using scale_longrange,
            # scale_coulomb, cutoff_radius, and omega.
            pbe = PBE("ORIG", scale_c=1, scale_x=0)
            xc_functional = XCFunctional(functionals=[], subsections={"PBE": pbe})

            potential_type = potential_type or "MIX_CL_TRUNC"
            hf_fraction = 1
            ip_keywords |= {
                "POTENTIAL_TYPE": Keyword("POTENTIAL_TYPE", potential_type),
                "CUTOFF_RADIUS": Keyword("CUTOFF_RADIUS", cutoff_radius),
                "T_C_G_DATA": Keyword("T_C_G_DATA", "t_c_g.dat"),
                "OMEGA": Keyword("OMEGA", omega),
                "SCALE_COULOMB": Keyword("SCALE_COULOMB", scale_coulomb),
                "SCALE_LONGRANGE": Keyword("SCALE_LONGRANGE", scale_longrange - scale_coulomb),
            }
            xc_functional.insert(
                Section(
                    "XWPBE",
                    subsections={},
                    keywords={
                        "SCALE_X0": Keyword("SCALE_X0", 1 - scale_longrange),
                        "SCALE_X": Keyword("SCALE_X", scale_longrange - scale_coulomb),
                        "OMEGA": Keyword("OMEGA", omega),
                    },
                )
            )
            xc_functional.insert(
                Section(
                    "PBE_HOLE_T_C_LR",
                    subsections={},
                    keywords={
                        "CUTOFF_RADIUS": Keyword("CUTOFF_RADIUS", cutoff_radius),
                        "SCALE_X": Keyword("SCALE_X", scale_longrange),
                    },
                )
            )
        else:
            warnings.warn(
                "Unknown hybrid functional. Using PBE base functional and overriding all "
                "settings manually. Proceed with caution."
            )
            pbe = PBE("ORIG", scale_c=gga_c_fraction, scale_x=gga_x_fraction)
            xc_functional = XCFunctional(functionals=[], subsections={"PBE": pbe})

            ip_keywords |= {
                "POTENTIAL_TYPE": Keyword("POTENTIAL_TYPE", potential_type),
                "CUTOFF_RADIUS": Keyword("CUTOFF_RADIUS", cutoff_radius),
                "T_C_G_DATA": Keyword("T_C_G_DATA", "t_c_g.dat"),
                "SCALE_COULOMB": Keyword("SCALE_COULOMB", scale_coulomb),
                "SCALE_GAUSSIAN": Keyword("SCALE_GAUSSIAN", scale_gaussian),
                "SCALE_LONGRANGE": Keyword("SCALE_LONGRANGE", scale_longrange),
                "OMEGA": Keyword("OMEGA", omega),
            }

        interaction_potential = Section("INTERACTION_POTENTIAL", subsections={}, keywords=ip_keywords)

        # Unlikely for users to override
        load_balance = Section(
            "LOAD_BALANCE",
            keywords={"RANDOMIZE": Keyword("RANDOMIZE", True)},  # noqa: FBT003
            subsections={},
        )

        # EPS_STORAGE_SCALING squashes the integrals for efficient storage
        # Unlikely for users to override.
        memory = Section(
            "MEMORY",
            subsections={},
            keywords={
                "EPS_STORAGE_SCALING": Keyword("EPS_STORAGE_SCALING", 0.1),
                "MAX_MEMORY": Keyword("MAX_MEMORY", max_memory),
            },
        )
        hf = Section(
            "HF",
            keywords={"FRACTION": Keyword("FRACTION", hf_fraction)},
            subsections={
                "SCREENING": screening,
                "INTERACTION_POTENTIAL": interaction_potential,
                "LOAD_BALANCE": load_balance,
                "MEMORY": memory,
            },
        )
        xc = Section("XC", subsections={"XC_FUNCTIONAL": xc_functional, "HF": hf})

        self.subsections["FORCE_EVAL"]["DFT"].insert(xc)

    def activate_motion(
        self,
        max_drift: float = 3e-3,
        rms_drift: float = 1.5e-3,
        max_force: float = 4.5e-4,
        rms_force: float = 3e-4,
        max_iter: int = 200,
        optimizer: str = "BFGS",
        trust_radius: float = 0.25,
        line_search: str = "2PNT",
        ensemble: str = "NVE",
        temperature: float = 300,
        timestep: float = 0.5,
        nsteps: int = 3,
        thermostat: str = "NOSE",
        nproc_rep: int = 1,
    ) -> None:
        """
        Turns on the motion section for GEO_OPT, CELL_OPT, etc. calculations.
        Will turn on the printing subsections and also bind any constraints
        to their respective atoms.
        """
        if not self.check("MOTION"):
            self.insert(Section("MOTION", subsections={}))

        run_type = self["global"].get("run_type", Keyword("run_type", "energy")).values[0].upper()
        run_type = {"GEOMETRY_OPTIMIZATION": "GEO_OPT", "MOLECULAR_DYNAMICS": "MD"}.get(run_type, run_type)

        self["MOTION"].insert(Section("PRINT", subsections={}))
        self["MOTION"]["PRINT"].insert(Section("TRAJECTORY", section_parameters=["ON"], subsections={}))
        self["MOTION"]["PRINT"].insert(Section("CELL", subsections={}))
        self["MOTION"]["PRINT"].insert(Section("FORCES", subsections={}))
        self["MOTION"]["PRINT"].insert(Section("STRESS", subsections={}))

        # ACTIVATE RELAX IF REQUESTED
        if run_type in ["GEO_OPT", "CELL_OPT"]:
            opt_params = {
                "MAX_DR": Keyword("MAX_DR", max_drift),
                "MAX_FORCE": Keyword("MAX_FORCE", max_force),
                "RMS_DR": Keyword("RMS_DR", rms_drift),
                "RMS_FORCE": Keyword("RMS_FORCE", rms_force),
                "MAX_ITER": Keyword("MAX_ITER", max_iter),
                "OPTIMIZER": Keyword("OPTIMIZER", optimizer),
            }
            opt = Section(run_type, subsections={}, keywords=opt_params)
            if optimizer.upper() == "CG":
                ls = Section("LINE_SEARCH", subsections={}, keywords={"TYPE": Keyword("TYPE", line_search)})
                cg = Section("CG", subsections={"LINE_SEARCH": ls}, keywords={})
                opt.insert(cg)
            elif optimizer.upper() == "BFGS":
                bfgs = Section("BFGS", subsections={}, keywords={"TRUST_RADIUS": Keyword("TRUST_RADIUS", trust_radius)})
                opt.insert(bfgs)

            self["MOTION"].insert(opt)

        # ACTIVATE MD IF REQUESTED
        elif run_type == "MD":
            md_keywords = {
                "ENSEMBLE": Keyword("ENSEMBLE", ensemble),
                "TEMPERATURE": Keyword("TEMPERATURE", temperature),
                "TIMESTEP": Keyword("TIMESTEP", timestep),
                "STEPS": Keyword("STEPS", nsteps),
            }
            thermostat = Section("THERMOSTAT", keywords={"TYPE": thermostat})
            md = Section("MD", subsections={"THERMOSTAT": thermostat}, keywords=md_keywords)
            self["MOTION"].insert(md)

        elif run_type == "BAND":
            convergence_control_params = {
                "MAX_DR": Keyword("MAX_DR", max_drift),
                "MAX_FORCE": Keyword("MAX_FORCE", max_force),
                "RMS_DR": Keyword("RMS_DR", rms_drift),
                "RMS_FORCE": Keyword("RMS_FORCE", rms_force),
            }
            band_kwargs = {
                "BAND_TYPE": Keyword("BAND_TYPE", "IT-NEB", description="Improved tangent NEB"),
                "NUMBER_OF_REPLICA": Keyword("NUMBER_OF_REPLICA"),
                "NPROC_REP": Keyword("NPROC_REP", nproc_rep),
            }
            band = Section("BAND", keywords=band_kwargs)
            band.insert(Section("CONVERGENCE_CONTROL", keywords=convergence_control_params))
            self["MOTION"].insert(band)

        self.modify_dft_print_iters(0, add_last="numeric")

        if "fix" in self.structure.site_properties:
            self["motion"].insert(Section("CONSTRAINT"))

            i = 0
            components = []
            tuples = []
            while i < len(self.structure):
                end = i + sum(
                    1
                    for j in itertools.takewhile(
                        lambda x: x == self.structure.site_properties["fix"][i],
                        self.structure.site_properties["fix"][i:],
                    )
                )
                components.append(self.structure.site_properties["fix"][i])
                tuples.append((i + 1, end))
                i = end
            self["motion"]["constraint"].insert(
                SectionList(
                    sections=[
                        Section(
                            "FIXED_ATOMS",
                            keywords={
                                "COMPONENTS_TO_FIX": Keyword("COMPONENTS_TO_FIX", c),
                                "LIST": Keyword("LIST", f"{t[0]}..{t[1]}"),
                            },
                        )
                        for t, c in zip(tuples, components, strict=True)
                        if c
                    ]
                )
            )

    def activate_tddfpt(self, **kwargs) -> None:
        """Activate TDDFPT for calculating excited states. Only works with GPW. Supports hfx."""
        if not self.check("force_eval/properties"):
            self["FORCE_EVAL"].insert(Section("PROPERTIES"))
        self["FORCE_EVAL"]["PROPERTIES"].insert(Section("TDDFPT", **kwargs))

    def activate_epr(self, **kwargs) -> None:
        """Calculate g-tensor. Requires localize. Suggested with GAPW."""
        if not self.check("force_eval/properties/linres/localize"):
            self.activate_localize()
        self["FORCE_EVAL"]["PROPERTIES"]["LINRES"].insert(Section("EPR", **kwargs))
        self["FORCE_EVAL"]["PROPERTIES"]["LINRES"]["EPR"] |= {"PRINT": {"G_TENSOR": {}}}

    def activate_nmr(self, **kwargs) -> None:
        """Calculate nmr shifts. Requires localize. Suggested with GAPW."""
        if not self.check("force_eval/properties/linres/localize"):
            self.activate_localize()
        self["FORCE_EVAL"]["PROPERTIES"]["LINRES"].insert(Section("NMR", **kwargs))
        self["FORCE_EVAL"]["PROPERTIES"]["LINRES"]["NMR"] |= {"PRINT": {"CHI_TENSOR": {}, "SHIELDING_TENSOR": {}}}

    def activate_spinspin(self, **kwargs) -> None:
        """Calculate spin-spin coupling tensor. Requires localize."""
        if not self.check("force_eval/properties/linres/localize"):
            self.activate_localize()
        self["FORCE_EVAL"]["PROPERTIES"]["LINRES"].insert(Section("SPINSPIN", **kwargs))

    def activate_polar(self, **kwargs) -> None:
        """Calculate polarizations (including raman)."""
        if not self.check("force_eval/properties"):
            self["FORCE_EVAL"].insert(Section("PROPERTIES"))
        if not self.check("force_eval/properties/linres"):
            self["FORCE_EVAL"]["PROPERTIES"].insert("LINRES")
        self["FORCE_EVAL"]["PROPERTIES"]["LINRES"].insert(Section("POLAR", **kwargs))

    def activate_hyperfine(self) -> None:
        """Print the hyperfine coupling constants."""
        self["FORCE_EVAL"]["DFT"]["PRINT"].insert(Section("HYPERFINE_COUPLING_TENSOR", FILENAME="HYPERFINE"))

    def activate_localize(self, states="OCCUPIED", preconditioner="FULL_ALL", restart=False) -> None:
        """
        Activate calculation of the maximally localized wannier functions.

        Args:
            states: Which states to calculate. occupied, unoccupied, mixed states, or all states. At
                present, unoccupied orbitals are only implemented for GPW.
            preconditioner: Preconditioner to use for optimize
            restart: Initialize from the localization restart file
        """
        if not self.check("force_eval/properties"):
            self["FORCE_EVAL"].insert(Section("PROPERTIES"))
        if not self.check("force_eval/properties/linres"):
            self["FORCE_EVAL"]["PROPERTIES"].insert(Section("LINRES"))

        self["FORCE_EVAL"]["PROPERTIES"]["LINRES"].insert(
            Section("LOCALIZE", PRECONDITIONER=preconditioner, STATES=states, RESTART=restart)
        )
        self["FORCE_EVAL"]["PROPERTIES"]["LINRES"]["LOCALIZE"].insert(Section("PRINT"))
        self["FORCE_EVAL"]["PROPERTIES"]["LINRES"]["LOCALIZE"]["PRINT"].insert(Section("LOC_RESTART"))

    def activate_vdw_potential(
        self,
        dispersion_functional: str,
        potential_type: str,
    ) -> None:
        """
        Activate van der Waals dispersion corrections.

        Args:
            dispersion_functional: Type of dispersion functional.
                Options: pair_potential or non_local
            potential_type: What type of potential to use, given a dispersion functional type
                Options: DFTD2, DFTD3, DFTD3(BJ), DRSLL, LMKLL, RVV10
        """
        vdw = Section(
            "VDW_POTENTIAL", keywords={"DISPERSION_FUNCTIONAL": Keyword("DISPERSION_FUNCTIONAL", dispersion_functional)}
        )
        keywords = {"TYPE": Keyword("TYPE", potential_type)}
        if dispersion_functional.upper() == "PAIR_POTENTIAL":
            reference_functional = self.xc_functionals[0]
            warnings.warn(
                "Reference functional will not be checked for validity. "
                "Calculation will fail if the reference functional "
                "does not exist in the dftd3 reference data"
            )
            keywords["PARAMETER_FILE_NAME"] = Keyword("PARAMETER_FILE_NAME", "dftd3.dat")
            keywords["REFERENCE_FUNCTIONAL"] = Keyword("REFERENCE_FUNCTIONAL", reference_functional)

        vdw.insert(Section(dispersion_functional, keywords=keywords))
        self["FORCE_EVAL"]["DFT"]["XC"].insert(vdw)

    def activate_fast_minimization(self, on) -> None:
        """Modify the set to use fast SCF minimization."""
        if on:
            ot = OrbitalTransformation(
                minimizer="DIIS",
                preconditioner="FULL_ALL",
                algorithm="IRAC",
                linesearch="2PNT",
            )
            self.update({"FORCE_EVAL": {"DFT": {"SCF": {"OT": ot}}}})  # type: ignore[assignment]

    def activate_robust_minimization(self) -> None:
        """Modify the set to use more robust SCF minimization technique."""
        ot = OrbitalTransformation(
            minimizer="CG",
            preconditioner="FULL_ALL",
            algorithm="STRICT",
            linesearch="3PNT",
        )
        self.update({"FORCE_EVAL": {"DFT": {"SCF": {"OT": ot}}}})  # type: ignore[assignment]

    def activate_very_strict_minimization(self) -> None:
        """Method to modify the set to use very strict SCF minimization scheme."""
        ot = OrbitalTransformation(
            minimizer="CG",
            preconditioner="FULL_ALL",
            algorithm="STRICT",
            linesearch="GOLD",
        )
        self.update({"FORCE_EVAL": {"DFT": {"SCF": {"OT": ot}}}})  # type: ignore[assignment]

    def activate_nonperiodic(self, solver="ANALYTIC") -> None:
        """
        Activates a calculation with non-periodic calculations by turning of PBC and
        changing the poisson solver. Still requires a CELL to put the atoms.
        """
        kwds = {
            "POISSON_SOLVER": Keyword("POISSON_SOLVER", solver),
            "PERIODIC": Keyword("PERIODIC", "NONE"),
        }
        self["FORCE_EVAL"]["DFT"].insert(Section("POISSON", subsections={}, keywords=kwds))

    def create_subsys(self, structure: Structure | Molecule) -> None:
        """Create the structure for the input."""
        subsys = Subsys()
        if isinstance(structure, Structure):
            subsys.insert(Cell(structure.lattice))
        else:
            x = max(*structure.cart_coords[:, 0], 1)
            y = max(*structure.cart_coords[:, 1], 1)
            z = max(*structure.cart_coords[:, 2], 1)
            cell = Cell(lattice=Lattice([[10 * x, 0, 0], [0, 10 * y, 0], [0, 0, 10 * z]]))
            cell.add(Keyword("PERIODIC", "NONE"))
            subsys.insert(cell)

        # Insert atom kinds by identifying the unique sites (unique element and site properties)
        unique_kinds = get_unique_site_indices(structure)
        for k, v in unique_kinds.items():
            kind = k.split("_")[0]
            kwargs = {}

            _ox = (
                self.structure.site_properties["oxi_state"][v[0]]
                if "oxi_state" in self.structure.site_properties
                else 0
            )
            _sp = self.structure.site_properties["spin"][v[0]] if "spin" in self.structure.site_properties else 0

            bs = BrokenSymmetry.from_el(kind, _ox, _sp) if _ox else None

            if "magmom" in self.structure.site_properties and not bs:
                kwargs["magnetization"] = self.structure.site_properties["magmom"][v[0]]

            if "ghost" in self.structure.site_properties:
                kwargs["ghost"] = self.structure.site_properties["ghost"][v[0]]

            if "basis_set" in self.structure.site_properties:
                basis_set = self.structure.site_properties["basis_set"][v[0]]
            else:
                basis_set = self.basis_and_potential[kind].get("basis")

            if "potential" in self.structure.site_properties:
                potential = self.structure.site_properties["potential"][v[0]]
            else:
                potential = self.basis_and_potential[kind].get("potential")

            if "aux_basis" in self.structure.site_properties:
                kwargs["aux_basis"] = self.structure.site_properties["aux_basis"][v[0]]
            elif self.basis_and_potential[kind].get("aux_basis"):
                kwargs["aux_basis"] = self.basis_and_potential[kind].get("aux_basis")

            _kind = Kind(
                kind,
                alias=k,
                basis_set=basis_set,
                potential=potential,
                subsections={"BS": bs} if bs else {},
                **kwargs,
            )
            if self.qs_method.upper() == "GAPW":
                _kind.add(Keyword("RADIAL_GRID", 200))
                _kind.add(Keyword("LEBEDEV_GRID", 80))

            subsys.insert(_kind)

        coord = Coord(structure, aliases=unique_kinds)
        subsys.insert(coord)
        self["FORCE_EVAL"].insert(subsys)

    def modify_dft_print_iters(self, iters, add_last: Literal["no", "numeric", "symbolic"] = "no"):
        """
        Modify all DFT print iterations at once. Common use is to set iters to the max
        number of iterations + 1 and then set add_last to numeric. This would have the
        effect of printing only the first and last iteration, which might be useful for
        speeding up/saving space on GEO_OPT or MD runs where you don't need the intermediate
        values.

        Args:
            iters (int): print each "iters" iterations.
            add_last (str): Whether to explicitly include the last iteration, and how to mark it.
                numeric: mark last iteration with the iteration number
                symbolic: mark last iteration with the letter "l"
                no: do not explicitly include the last iteration
        """
        if add_last.lower() not in {"no", "numeric", "symbolic"}:
            raise ValueError(f"add_list should be no/numeric/symbolic, got {add_last.lower()}")

        run_type = self["global"].get("run_type", Keyword("run_type", "energy")).values[0].upper()
        if run_type not in ["ENERGY_FORCE", "ENERGY", "WAVEFUNCTION_OPTIMIZATION", "WFN_OPT"] and self.check(
            "FORCE_EVAL/DFT/PRINT"
        ):
            for v in self["force_eval"]["dft"]["print"].subsections.values():
                if v.name.upper() in [
                    "ACTIVE_SPACE",
                    "BAND_STRUCTURE",
                    "GAPW",
                    "IMPLICIT_PSOLVER",
                    "SCCS",
                    "WFN_MIX",
                ]:
                    continue

                v.insert(Section("EACH", subsections=None, keywords={run_type: Keyword(run_type, iters)}))
                v.keywords["ADD_LAST"] = Keyword("ADD_LAST", add_last)

    def validate(self):
        """Implements a few checks for a valid input set."""
        if self.check("force_eval/dft/kpoints") and self.check("force_eval/dft/xc/hf"):
            raise Cp2kValidationError("Does not support hartree fock with kpoints")

        for val in self["force_eval"]["subsys"].subsections.values():
            if (
                val.name.upper() == "KIND"
                and val["POTENTIAL"].values[0].upper() == "ALL"
                and self["force_eval"]["dft"]["qs"]["method"].values[0].upper() != "GAPW"
            ):
                raise Cp2kValidationError("All electron basis sets require GAPW method")


class StaticSet(DftSet):
    """Quick Constructor for static calculations."""

    def __init__(self, **kwargs) -> None:
        super().__init__(run_type="ENERGY_FORCE", **kwargs)


class RelaxSet(DftSet):
    """Quick Constructor for geometry relaxation."""

    def __init__(self, **kwargs) -> None:
        super().__init__(run_type="GEO_OPT", **kwargs)
        self.activate_motion()


class CellOptSet(DftSet):
    """Quick Constructor for cell optimization relaxation."""

    def __init__(self, **kwargs) -> None:
        super().__init__(run_type="CELL_OPT", **kwargs)
        self.activate_motion()


class HybridStaticSet(DftSet):
    """Quick Constructor for static calculations."""

    def __init__(self, **kwargs) -> None:
        super().__init__(run_type="ENERGY_FORCE", **kwargs)
        self.activate_hybrid(**kwargs.get("activate_hybrid", {}))


class HybridRelaxSet(DftSet):
    """Quick Constructor for hybrid geometry relaxation."""

    def __init__(self, **kwargs) -> None:
        super().__init__(run_type="GEO_OPT", **kwargs)
        self.activate_hybrid(**kwargs.get("activate_hybrid", {}))


class HybridCellOptSet(DftSet):
    """Quick Constructor for hybrid cell optimization relaxation."""

    def __init__(self, **kwargs) -> None:
        super().__init__(run_type="CELL_OPT", **kwargs)
        self.activate_hybrid(**kwargs.get("activate_hybrid", {}))


class Cp2kValidationError(Exception):
    """
    CP2K validation exception. Not exhausted. May raise validation
    errors for features which actually do work if using a newer version
    of CP2K.
    """

    CP2K_VERSION = "v2022.1"

    def __init__(self, message) -> None:
        message = f"CP2K {self.CP2K_VERSION}: {message}"
        super().__init__(message)
