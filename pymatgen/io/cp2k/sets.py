"""
This module defines input sets for CP2K and is a work in progress. The structure/philosophy
of this module is based on the VASP input sets in Pymatgen. These sets are meant to contain
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
from pathlib import Path

from ruamel.yaml import YAML

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cp2k.inputs import (
    LDOS,
    PBE,
    PDOS,
    QS,
    BrokenSymmetry,
    Cell,
    Coord,
    Cp2kInput,
    Dft,
    E_Density_Cube,
    ForceEval,
    Global,
    Keyword,
    KeywordList,
    Kind,
    Kpoints,
    Mgrid,
    MO_Cubes,
    OrbitalTransformation,
    Scf,
    Section,
    SectionList,
    Smear,
    Subsys,
    V_Hartree_Cube,
    Xc_Functional,
)
from pymatgen.io.cp2k.utils import (
    get_aux_basis,
    get_basis_and_potential,
    get_cutoff_from_basis,
    get_truncated_coulomb_cutoff,
    get_unique_site_indices,
    get_xc_functionals,
)

__author__ = "Nicholas Winner"
__version__ = "1.0"
__email__ = "nwinner@berkeley.edu"
__date__ = "March 2022"

MODULE_DIR = Path(__file__).resolve().parent

with open(os.path.join(MODULE_DIR, "settings.yaml")) as f:
    yaml = YAML(typ="unsafe", pure=True)
    SETTINGS = yaml.load(f)


class Cp2kInputSet(Cp2kInput):

    """
    The basic representation of a CP2K input set as a collection of "sections" defining the simulation
    connected to a structure object. At the most basis level, CP2K requires a &GLOBAL section and
    &FORCE_EVAL section. Global sets parameters like "RUN_TYPE" or the overall verbosity. FORCE_EVAL is
    the largest section usually, containing the cell and coordinates of atoms, the DFT settings, and more.
    This top level input set is meant to initialize GLOBAL and FORCE_EVAL based on a structure object and
    and sections that the user provides.

    Like everything that goes into a cp2k input file, this base input set is essentially a section object.
    These sets are distinguished by saving default settings for easy implementation of calculations such
    as relaxation and static calculations. This base set is here to transfer a pymatgen structure object
    into the input format for cp2k and associate the basis set and pseudopotential to use with each
    element in the structure.

    Generally, this class will not be used directly, and instead one of
    its child-classes will be used, which contain more predefined initializations of various sections, and,
    if modifications are required, the user can specify override_default_settings.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        basis_and_potential: dict | str = "preferred",
        multiplicity: int = 0,
        project_name: str = "CP2K",
        override_default_params: dict = None,
        **kwargs,
    ):
        """
        Args:
            structure: (Structure or Molecule) pymatgen structure or molecule object used to define
                the lattice, coordinates, and elements. This structure object cannot contain "special"
                species like the Dummy species, e.g. X, or fractional occupations, e.g. Fe0.2, etc.
            potential_and_basis: (dict) Specifies what basis set and potential to use. Specify these
                as a dict of the form:
                    { element: {'cardinality': __, 'sr': __, 'q': __},
                     'cardinality': __, 'functional': __}
                Where cardinality and functional are overall specifications (for all elements), while
                <key='element'> specifies the overrides for a specific element. Currently the following
                conventions must be followed:
                    (a) All species of a particular element must have the same potential/basis
            multiplicity: (int) Specify the system's multiplicity if appropriate
            project_name: (str) Specify the project name. This will be used to name the output files
                from a CP2K calculation
            override_default_params: (dict) Specifies user-defined settings to override the settings of any
                input set (See Section.update())
        """
        super().__init__(name="CP2K_INPUT", subsections={})

        # Important CP2K set parameters
        self.structure = structure
        self.charge = int(structure.charge)
        self.basis_and_potential = basis_and_potential
        self.multiplicity = multiplicity  # spin multiplicity = 2s+1
        self.override_default_params = override_default_params or {}
        self.project_name = project_name
        self.kwargs = kwargs

        self.insert(ForceEval())  # always present in cp2k
        self.basis_set_file_names = None  # need for dft
        self.potential_file_name = None  # need for dft
        self.create_subsys(self.structure)  # assemble structure with atom types and pseudopotentials assigned

        if self.kwargs.get("print_forces", True):
            self.print_forces()

        self.update(self.override_default_params)

    def create_subsys(self, structure: Structure | Molecule):
        """
        Create the structure for the input
        """
        subsys = Subsys()
        if isinstance(structure, Structure):
            subsys.insert(Cell(structure.lattice))

        # Decide what basis sets/pseudopotentials to use
        basis_and_potential = get_basis_and_potential(structure.symbol_set, self.basis_and_potential)

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
                basis_set = basis_and_potential[kind]["basis"]

            if "potential" in self.structure.site_properties:
                potential = self.structure.site_properties["potential"][v[0]]
            else:
                potential = basis_and_potential[kind]["potential"]

            if "aux_basis" in self.structure.site_properties:
                kwargs["aux_basis"] = self.structure.site_properties["aux_basis"][v[0]]
                if kwargs["aux_basis"] == "match":
                    kwargs["aux_basis"] = basis_set

            _kind = Kind(
                kind, alias=k, basis_set=basis_set, potential=potential, subsections={"BS": bs} if bs else {}, **kwargs
            )

            subsys.insert(_kind)

        coord = Coord(structure, aliases=unique_kinds)
        subsys.insert(coord)
        self["FORCE_EVAL"].insert(subsys)
        self.basis_set_file_names = basis_and_potential["basis_filenames"]
        self.potential_file_name = basis_and_potential["potential_filename"]

    def print_forces(self):
        """
        Print out the forces and stress during calculation
        """
        self["FORCE_EVAL"].insert(Section("PRINT", subsections={}))
        self["FORCE_EVAL"]["PRINT"].insert(Section("FORCES", subsections={}))
        self["FORCE_EVAL"]["PRINT"].insert(Section("STRESS_TENSOR", subsections={}))


class DftSet(Cp2kInputSet):
    """
    Base for an input set using the Quickstep module (i.e. a DFT calculation). The DFT section is pretty vast
    in CP2K, so this set hopes to make the DFT setup fairly simple. The provided parameters are pretty conservative,
    and so they should not need to be changed very often.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        ot: bool = True,
        energy_gap: float = -1,
        eps_default: float = 1e-12,
        eps_scf: float = 1e-6,
        max_scf: int | None = None,
        minimizer: str = "DIIS",
        preconditioner: str = "FULL_SINGLE_INVERSE",
        algorithm: str = "IRAC",
        linesearch: str = "2PNT",
        rotation: bool = True,
        occupation_preconditioner: bool = False,
        cutoff: int = 0,
        rel_cutoff: int = 50,
        ngrids: int = 5,
        progression_factor: int = 3,
        override_default_params: dict = None,
        wfn_restart_file_name: str = None,
        kpoints: Kpoints | None = None,
        smearing: bool = False,
        **kwargs,
    ):
        """
        Args:
            structure: Pymatgen structure or molecule object
            ot (bool): Whether or not to use orbital transformation method for matrix diagonalization. OT is the
                flagship scf solver of CP2K, and will provide huge speed-ups for this part of the calculation,
                but the system must have a band gap for OT to be used (higher band-gap --> faster convergence).
                Band gap is also used by the preconditioner for OT, and should be set as a value SMALLER than the true
                band gap to get good efficiency. Generally, this parameter does not need to be changed from
                default of 0.01
            energy_gap (float): Estimate of energy gap for preconditioner. Default is -1, leaving it up to cp2k.
            eps_default (float): Replaces all EPS_XX Keywords in the DFT section (NOT its subsections!) to have this
                value, ensuring an overall accuracy of at least this much.
            eps_scf (float): The convergence criteria for leaving the SCF loop. Default is 1e-6. Should
                ensure reasonable results, but is not applicable to all situations.
                    Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method, it is the largest
                    gradient of the energy with respect to changing any of the molecular orbital coefficients. For
                    diagonalization, it is the largest change in the density matrix from the last step.
            max_scf (int): The max number of SCF cycles before terminating the solver. NOTE: With the OT solver, this
                corresponds to the max number of INNER scf loops, and then the outer loops are set with outer_max_scf,
                while with diagonalization it corresponds to the overall (INNER*OUTER) number of SCF steps, with the
                inner loop limit set by
            minimizer (str): The minimization scheme. DIIS can be as much as 50% faster than the more robust conjugate
                gradient method, and so it is chosen as default. Switch to CG if dealing with a difficult system.
            preconditioner (str): Preconditioner for the OT method. FULL_SINGLE_INVERSE is very robust and compatible
                with non-integer occupations from IRAC+rotation. FULL_ALL is considered "best" but needs algorithm to
                be set to STRICT. Only change from these two when simulation cell gets to be VERY large, in which case
                FULL_KINETIC might be preferred.
            algorithm (str): Algorithm for the OT method. STRICT assumes that the orbitals are strictly orthogonal to
                each other, which works well for wide gap ionic systems, but can diverge for systems with small gaps,
                fractional occupations, and some other cases. IRAC (iterative refinement of the approximate congruency)
                transformation is not analytically correct and uses a truncated polynomial expansion, but is robust to
                the problems with STRICT, and so is the default.
            linesearch (str): Linesearch method for CG. 2PNT is the default, and is the fastest, but is not as robust
                as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility with irac+rotation. This may be upgraded
                in the future. 3PNT can be good for wide gapped transition metal systems as an alternative.
            rotation (bool): Whether or not to allow for rotation of the orbitals in the OT method. This equates to
                allowing for fractional occupations in the calculation.
            occupation_preconditioner (bool): Whether or not to account for fractional occupations in the
                preconditioner. This method is not fully integrated as of cp2k v9.1 and is set to false by default.
            cutoff (int): Cutoff energy (in Ry) for the finest level of the multigrid. A high cutoff will allow you to
                have very accurate calculations PROVIDED that REL_CUTOFF is appropriate. By default cutoff is set to 0,
                which will assign it to be the largest exponent of your basis times the rel_cutoff.
            rel_cutoff (int): This cutoff decides how the Gaussians are mapped onto the different levels of the
                multigrid. If REL_CUTOFF is too low, then even if you have a high CUTOFF, all Gaussians will be
                mapped onto the coarsest level of the multi-grid, and thus the effective integration grid for
                the calculation may still be too coarse. By default 50Ry is chosen, which should be sufficient
                given the cutoff is large enough.
                    From CP2K manual: A Gaussian is mapped onto the coarsest level of the multi-grid, on which the
                    function will cover number of grid points greater than or equal to the number of grid points
                    will cover on a reference grid defined by REL_CUTOFF.
            ngrids (int): number of multi-grids to use. CP2K default is 4, but the molopt basis files recommend 5.
            progression_factor (int): Divisor of CUTOFF to get the cutoff for the next level of the multigrid.
            wfn_restart_file_name (str): RESTART file for the initial wavefunction guess.
            kpoints (Kpoints): kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default, CP2K runs with gamma
                point only.
            smearing (bool): whether or not to activate smearing (should be done for systems containing no (or a very
                small) band gap.
        """

        super().__init__(structure, **kwargs)

        self.structure = structure
        self.ot = ot
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

        # Build the QS Section
        qs = QS(eps_default=eps_default)
        max_scf = max_scf or 20 if ot else 400  # If ot, max_scf is for inner loop
        scf = Scf(eps_scf=eps_scf, max_scf=max_scf, subsections={})

        # If there's a band gap, use OT, else use diagonalization
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

        # Create the multigrid for FFTs
        if not cutoff:
            basis_and_potential = get_basis_and_potential(structure.symbol_set, self.basis_and_potential)
            cutoff = get_cutoff_from_basis(
                els=self.structure.symbol_set,
                bases=[basis_and_potential[s]["basis"] for s in self.structure.symbol_set],
                rel_cutoff=rel_cutoff,
            )
        mgrid = Mgrid(
            cutoff=cutoff,
            rel_cutoff=rel_cutoff,
            ngrids=ngrids,
            progression_factor=progression_factor,
        )

        # Set the DFT calculation with global parameters
        dft = Dft(
            MULTIPLICITY=self.multiplicity,
            CHARGE=self.charge,
            basis_set_filenames=self.basis_set_file_names,
            potential_filename=self.potential_file_name,
            subsections={"QS": qs, "SCF": scf, "MGRID": mgrid},
            wfn_restart_file_name=wfn_restart_file_name,
        )

        # Set kpoints only if user supplies them
        if kpoints:
            dft.insert(Kpoints.from_kpoints(kpoints, structure=self.structure, reduce=True))
        if smearing:
            scf.kwargs["ADDED_MOS"] = 500
            scf["ADDED_MOS"] = 500  # TODO: how to grab the appropriate number?
            scf.insert(Smear(elec_temp=kwargs.get("elec_temp", 300)))

        # Create subsections and insert into them
        self["FORCE_EVAL"].insert(dft)
        xc_functionals = get_xc_functionals(kwargs.get("xc_functional", "PBE"))
        xc_functional = Xc_Functional(functionals=xc_functionals)
        xc = Section("XC", subsections={"XC_FUNCTIONAL": xc_functional})
        self["FORCE_EVAL"]["DFT"].insert(xc)
        self["FORCE_EVAL"]["DFT"].insert(Section("PRINT", subsections={}))

        if isinstance(structure, Molecule):
            self.activate_nonperiodic()

        if kwargs.get("print_pdos", True):
            self.print_pdos()
        if kwargs.get("print_ldos", False):
            self.print_ldos()
        if kwargs.get("print_mo_cubes", True):
            self.print_mo_cubes()
        if kwargs.get("print_hartree_potential", True):
            self.print_hartree_potential()
        if kwargs.get("print_e_density", True):
            self.print_e_density()

        self.update(self.override_default_params)

    def print_pdos(self, nlumo=-1):
        """
        Activate creation of the PDOS file.

        Args:
            nlumo (int): Number of virtual orbitals to be added to the MO set (-1=all).
                CAUTION: Setting this value to be higher than the number of states present may cause a Cholesky error.
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/PDOS"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(PDOS(nlumo=nlumo))

    def print_ldos(self, nlumo=-1):
        """
        Activate the printing of LDOS files, printing one for each atom kind by default

        Args:
            nlumo (int): Number of virtual orbitals to be added to the MO set (-1=all).
                CAUTION: Setting this value to be higher than the number of states present may cause a Cholesky error.
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/PDOS"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(PDOS(nlumo=nlumo))
        for i in range(self.structure.num_sites):
            self["FORCE_EVAL"]["DFT"]["PRINT"]["PDOS"].insert(LDOS(i + 1, alias=f"LDOS {i + 1}", verbose=False))

    def print_mo_cubes(self, write_cube=False, nlumo=-1, nhomo=-1):
        """
        Activate printing of molecular orbitals.

        Args:
            write_cube (bool): whether to write cube file for the MOs (setting false will just print levels in out file)
            nlumo (int): Controls the number of lumos that are printed and dumped as a cube (-1=all)
            nhomo (int): Controls the number of homos that are printed and dumped as a cube (-1=all)
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/MO_CUBES"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(MO_Cubes(write_cube=write_cube, nlumo=nlumo, nhomo=nhomo))

    def print_mo(self):
        """
        Print molecular orbitals when running non-OT diagonalization
        """
        raise NotImplementedError

    def print_hartree_potential(self, stride=(2, 2, 2)):
        """
        Controls the printing of a cube file with eletrostatic potential generated by the total density
        (electrons+ions). It is valid only for QS with GPW formalism.
        Note that by convention the potential has opposite sign than the expected physical one.
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(V_Hartree_Cube(keywords={"STRIDE": Keyword("STRIDE", *stride)}))

    def print_e_density(self, stride=(2, 2, 2)):
        """
        Controls the printing of cube files with the electronic density and, for LSD calculations, the spin density
        """
        if not self.check("FORCE_EVAL/DFT/PRINT/E_DENSITY_CUBE"):
            self["FORCE_EVAL"]["DFT"]["PRINT"].insert(E_Density_Cube(keywords={"STRIDE": Keyword("STRIDE", *stride)}))

    def set_charge(self, charge: int):
        """
        Set the overall charge of the simulation cell
        """
        self["FORCE_EVAL"]["DFT"]["CHARGE"] = Keyword("CHARGE", int(charge))

    def activate_hybrid(
        self,
        hybrid_functional: str = "PBE0",
        hf_fraction: float = 0.25,
        gga_x_fraction: float = 0.75,
        gga_c_fraction: float = 1,
        max_memory: int = 2000,
        cutoff_radius: float = 8.0,
        potential_type: str = None,
        scale_coulomb: float = 1,
        scale_gaussian: float = 1,
        scale_longrange: float = 1,
        omega: float = 0.11,
        aux_basis: dict | None = None,
        admm: bool = True,
        eps_schwarz: float = 1e-7,
        eps_schwarz_forces: float = 1e-6,
        screen_on_initial_p: bool = True,
        screen_p_forces: bool = True,
    ):
        """
        Basic set for activating hybrid DFT calculation using Auxiliary Density Matrix Method.

        Note 1: When running ADMM with cp2k, memory is very important. If the memory requirements exceed
        what is available (see max_memory), then CP2K will have to calculate the 4-electron integrals
        for HFX during each step of the SCF cycle. ADMM provides a huge speed up by making the memory
        requirements *feasible* to fit into RAM, which means you only need to calculate the integrals
        once each SCF cycle. But, this only works if it fits into memory. When setting up ADMM
        calculations, we recommend doing whatever is possible to fit all the 4EI into memory.

        Note 2: This set is designed for reliable high-throughput calculations, NOT for extreme
        accuracy. Please review the in-line comments in this method if you want more control.

        Args:
            hybrid_functional (str): Type of hybrid functional. This set supports HSE (screened) and PBE0
                (truncated). Default is PBE0, which converges easier in the GPW basis used by
                cp2k.
            hf_fraction (float): fraction of exact HF exchange energy to mix. Default: 0.25
            gga_x_fraction (float): fraction of gga exchange energy to retain. Default: 0.75
            gga_c_fraction (float): fraction of gga correlation energy to retain. Default: 1.0
            max_memory (int): Maximum memory available to each MPI process (in Mb) in the calculation.
                Most modern computing nodes will have ~2Gb per core, or 2048 Mb, but check for
                your specific system. This value should be as large as possible while still leaving
                some memory for the other parts of cp2k. Important: If this value is set larger
                than the memory limits, CP2K will likely seg-fault.
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
            aux_basis (dict): If you want to specify the aux basis to use, specify it as a dict of
                the form {'specie_1': 'AUX_BASIS_1', 'specie_2': 'AUX_BASIS_2'}
            admm (bool): Whether or not to use the auxiliary density matrix method for the exact
                HF exchange contribution. Highly recommended. Speed ups between 10x and aaa1000x are
                possible when compared to non ADMM hybrid calculations. Default: True
            eps_schwarz (float): Screening threshold for HFX, in Ha. Contributions smaller than this
                will be screened. The smaller the value, the more accurate, but also the more
                costly. Default value is 1e-7. 1e-6 works in a large number of cases, but is
                quite aggressive, which can lead to convergence issues.
            eps_schwarz_forces (float): Same as for eps_schwarz, but for screening contributions to
                forces. Convergence is not as sensitive with respect to eps_schwarz forces as
                compared to eps_schwarz, and so 1e-6 should be good default.
            screen_on_initial_p (bool): If an initial density matrix is provided, in the form of a
                CP2K wfn restart file, then this initial density will be used for screening. This
                is generally very computationally efficient, but, as with eps_schwarz, can lead to
                instabilities if the initial density matrix is poor.
            screen_p_forces (bool): Same as screen_on_initial_p, but for screening of forces.
        """
        if admm:
            aux_basis = aux_basis or {}
            aux_basis = {s: aux_basis[s] if s in aux_basis else None for s in self.structure.symbol_set}
            basis = get_aux_basis(basis_type=aux_basis)
            if isinstance(self["FORCE_EVAL"]["DFT"]["BASIS_SET_FILE_NAME"], KeywordList):
                self["FORCE_EVAL"]["DFT"]["BASIS_SET_FILE_NAME"].extend(
                    [Keyword("BASIS_SET_FILE_NAME", k) for k in ["BASIS_ADMM", "BASIS_ADMM_MOLOPT"]],
                )

            for k, v in self["FORCE_EVAL"]["SUBSYS"].subsections.items():
                if v.name.upper() == "KIND":
                    if isinstance(v.keywords["BASIS_SET"], KeywordList) and any(
                        "AUX_FIT" in k.values for k in v.keywords["BASIS_SET"]
                    ):
                        continue
                    if any(k.upper() == "AUX_FIT" for k in v.keywords["BASIS_SET"].values):
                        continue
                    kind = v["ELEMENT"].values[0]
                    v.keywords["BASIS_SET"] += Keyword("BASIS_SET", "AUX_FIT", basis[kind])

            # Don't change unless you know what you're doing
            # Use NONE for accurate eigenvalues (static calcs)
            aux_matrix_params = {
                "ADMM_PURIFICATION_METHOD": Keyword("ADMM_PURIFICATION_METHOD", "NONE"),
                "METHOD": Keyword("METHOD", "BASIS_PROJECTION"),
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

        ip_keywords = {}
        if hybrid_functional == "HSE06":
            pbe = PBE("ORIG", scale_c=1, scale_x=0)
            xc_functional = Xc_Functional(functionals=[], subsections={"PBE": pbe})

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
            ip_keywords.update(
                {
                    "POTENTIAL_TYPE": Keyword("POTENTIAL_TYPE", potential_type),
                    "OMEGA": Keyword("OMEGA", 0.11),
                    "CUTOFF_RADIUS": Keyword("CUTOFF_RADIUS", cutoff_radius),
                }
            )
        elif hybrid_functional == "PBE0":
            pbe = PBE("ORIG", scale_c=1, scale_x=1 - hf_fraction)
            xc_functional = Xc_Functional(functionals=[], subsections={"PBE": pbe})
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

            if isinstance(self.structure, Molecule):
                potential_type = "COULOMB"
            else:
                potential_type = "TRUNCATED"
            ip_keywords.update(
                {
                    "POTENTIAL_TYPE": Keyword("POTENTIAL_TYPE", potential_type),
                    "CUTOFF_RADIUS": Keyword("CUTOFF_RADIUS", cutoff_radius),
                    "T_C_G_DATA": Keyword("T_C_G_DATA", "t_c_g.dat"),
                }
            )
        elif hybrid_functional == "RSH":
            """
            Activates range separated functional using mixing of the truncated
            coulomb operator and the long range operator using scale_longrange,
            scale_coulomb, cutoff_radius, and omega.
            """
            pbe = PBE("ORIG", scale_c=1, scale_x=0)
            xc_functional = Xc_Functional(functionals=[], subsections={"PBE": pbe})

            potential_type = potential_type or "MIX_CL_TRUNC"
            hf_fraction = 1
            ip_keywords.update(
                {
                    "POTENTIAL_TYPE": Keyword("POTENTIAL_TYPE", potential_type),
                    "CUTOFF_RADIUS": Keyword("CUTOFF_RADIUS", cutoff_radius),
                    "T_C_G_DATA": Keyword("T_C_G_DATA", "t_c_g.dat"),
                    "OMEGA": Keyword("OMEGA", omega),
                    "SCALE_COULOMB": Keyword("SCALE_COULOMB", scale_coulomb),
                    "SCALE_LONGRANGE": Keyword("SCALE_LONGRANGE", scale_longrange - scale_coulomb),
                }
            )
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
                "Uknown hybrid functional. Using PBE base functional"
                " and overriding all settings manually. Proceed with"
                " caution."
            )
            pbe = PBE("ORIG", scale_c=gga_c_fraction, scale_x=gga_x_fraction)
            xc_functional = Xc_Functional(functionals=[], subsections={"PBE": pbe})

            ip_keywords.update(
                {
                    "POTENTIAL_TYPE": Keyword("POTENTIAL_TYPE", potential_type),
                    "CUTOFF_RADIUS": Keyword("CUTOFF_RADIUS", cutoff_radius),
                    "T_C_G_DATA": Keyword("T_C_G_DATA", "t_c_g.dat"),
                    "SCALE_COULOMB": Keyword("SCALE_COULOMB", scale_coulomb),
                    "SCALE_GAUSSIAN": Keyword("SCALE_GAUSSIAN", scale_gaussian),
                    "SCALE_LONGRANGE": Keyword("SCALE_LONGRANGE", scale_longrange),
                    "OMEGA": Keyword("OMEGA", omega),
                }
            )

        interaction_potential = Section("INTERACTION_POTENTIAL", subsections={}, keywords=ip_keywords)

        # Unlikely for users to override
        load_balance = Section(
            "LOAD_BALANCE",
            keywords={"RANDOMIZE": Keyword("RANDOMIZE", True)},
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

    def activate_motion(self):
        """
        Turns on the motion section for GEO_OPT, CELL_OPT, etc. calculations.
        Will turn on the printing subsections and also bind any constraints
        to their respective atoms.
        """
        if not self.check("MOTION"):
            self.insert(Section("MOTION", subsections={}))

        self["MOTION"].insert(Section("PRINT", subsections={}))
        self["MOTION"]["PRINT"].insert(Section("TRAJECTORY", section_parameters=["ON"], subsections={}))
        self["MOTION"]["PRINT"].insert(Section("CELL", subsections={}))
        self["MOTION"]["PRINT"].insert(Section("FORCES", subsections={}))
        self["MOTION"]["PRINT"].insert(Section("STRESS", subsections={}))

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
                        for t, c in zip(tuples, components)
                        if c
                    ]
                )
            )

    def activate_fast_minimization(self, on):
        """
        Method to modify the set to use fast SCF minimization.
        """
        if on:
            ot = OrbitalTransformation(
                minimizer="DIIS",
                preconditioner="FULL_ALL",
                algorithm="IRAC",
                linesearch="2PNT",
            )
            self.update({"FORCE_EVAL": {"DFT": {"SCF": {"OT": ot}}}})

    def activate_robust_minimization(self):
        """
        Method to modify the set to use more robust SCF minimization technique
        """
        ot = OrbitalTransformation(
            minimizer="CG",
            preconditioner="FULL_ALL",
            algorithm="STRICT",
            linesearch="3PNT",
        )
        self.update({"FORCE_EVAL": {"DFT": {"SCF": {"OT": ot}}}})

    def activate_very_strict_minimization(self):
        """
        Method to modify the set to use very strict SCF minimization scheme
        :return:
        """
        ot = OrbitalTransformation(
            minimizer="CG",
            preconditioner="FULL_ALL",
            algorithm="STRICT",
            linesearch="GOLD",
        )
        self.update({"FORCE_EVAL": {"DFT": {"SCF": {"OT": ot}}}})

    def activate_nonperiodic(self, solver="MT"):
        """
        Activates a calculation with non-periodic calculations by turning of PBC and
        changing the poisson solver. Still requires a CELL to put the atoms
        """
        if not self.check("FORCE_EVAL/SUBSYS/CELL"):
            x = max(s.coords[0] for s in self.structure.sites)
            y = max(s.coords[1] for s in self.structure.sites)
            z = max(s.coords[2] for s in self.structure.sites)
            self["FORCE_EVAL"]["SUBSYS"].insert(Cell(lattice=Lattice([[x, 0, 0], [0, y, 0], [0, 0, z]])))
        self["FORCE_EVAL"]["SUBSYS"]["CELL"].add(Keyword("PERIODIC", "NONE"))
        kwds = {"POISSON_SOLVER": Keyword("POISSON_SOLVER", solver), "PERIODIC": Keyword("PERIODIC", "NONE")}
        self["FORCE_EVAL"]["DFT"].insert(Section("POISSON", subsections={}, keywords=kwds))

    def modify_dft_print_iters(self, iters, add_last="no"):
        """
        Modify all DFT print iterations at once. Common use is to set iters to the max
        number of iterations + 1 and then set add_last to numeric. This would have the
        effect of printing only the first and last iteration, which might be useful for
        speeding up/saving space on GEO_OPT or MD runs where you don't need the intermediate
        values.

        Args
            iters (int): print each "iters" iterations.
            add_last (str): Whether to explicitly include the last iteration, and how to mark it.
                numeric: mark last iteration with the iteration number
                symbolic: mark last iteration with the letter "l"
                no: do not explicitly include the last iteration
        """
        assert add_last.lower() in ["no", "numeric", "symbolic"]
        if self.check("FORCE_EVAL/DFT/PRINT"):
            run_type = self["global"].get("run_type", Keyword("run_type", "energy")).values[0]
            for v in self["force_eval"]["dft"]["print"].subsections.values():

                if v.name.upper() in ["ACTIVE_SPACE", "BAND_STRUCTURE", "GAPW", "IMPLICIT_PSOLVER", "SCCS", "WFN_MIX"]:
                    continue

                v.insert(Section("EACH", subsections=None, keywords={run_type: Keyword(run_type, iters)}))
                v.keywords["ADD_LAST"] = Keyword("ADD_LAST", add_last)


class StaticSet(DftSet):

    """
    Basic static energy calculation. Turns on Quickstep module, sets the run_type in global,
    and uses structure object to build the subsystem.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        project_name: str = "Static",
        run_type: str = "ENERGY_FORCE",
        override_default_params: dict = None,
        **kwargs,
    ):
        """
        Args:
            structure: Pymatgen structure object
            project_name (str): What to name this cp2k project (controls naming of files printed out)
            run_type (str): Run type. As a static set it should be one of the static aliases, like 'ENERGY_FORCE'
        """
        super().__init__(structure, **kwargs)
        global_section = Global(project_name=project_name, run_type=run_type)
        self.structure = structure
        self.project_name = project_name
        self.run_type = run_type
        self.override_default_params = override_default_params or {}
        self.insert(global_section)
        self.update(self.override_default_params)
        self.kwargs = kwargs


class RelaxSet(DftSet):

    """
    CP2K input set containing the basic settings for performing geometry optimization. Values are all cp2k
    defaults, and should be good for most systems of interest.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        max_drift: float = 3e-3,
        rms_drift: float = 1.5e-3,
        max_force: float = 4.5e-4,
        rms_force: float = 3e-4,
        max_iter: int = 200,
        project_name: str = "Relax",
        optimizer: str = "BFGS",
        override_default_params: dict = None,
        **kwargs,
    ):
        """
        Args:
            structure: Pymatgen structure object
            max_drift: Convergence criterion for the maximum geometry change between the current and the
                last optimizer iteration. This keyword cannot be repeated and it expects precisely one real.
                Default value: 1.5.00000000E-003
                Default unit: [bohr]
            rms_drift: Convergence criterion for the RMS geometry change between the current and the
                last optimizer iteration. This keyword cannot be repeated and it expects precisely one real.
                Default value: 1.00000000E-003
                Default unit: [bohr]
            max_force (float): Convergence criterion for the maximum force component of the current configuration.
                This keyword cannot be repeated and it expects precisely one real.
                Default value: 1e-3
                Default unit: [bohr^-1*hartree]
            rms_force (float): Convergence criterion for the RMS force component of the current configuration.
                This keyword cannot be repeated and it expects precisely one real.
                Default value: 1e-3
                Default unit: [bohr^-1*hartree]
            max_iter (int): Specifies the maximum number of geometry optimization steps.
                One step might imply several force evaluations for the CG and LBFGS optimizers.
                This keyword cannot be repeated and it expects precisely one integer.
                Default value: 200
            optimizer (str): Specify which method to use to perform a geometry optimization.
                This keyword cannot be repeated and it expects precisely one keyword. BFGS is a
                quasi-newtonian method, and will best for "small" systems near the minimum. LBFGS
                is a limited memory version that can be used for "large" (>1000 atom) systems when
                efficiency outweighs robustness. CG is more robust, especially when you are far from
                the minimum, but it slower.
                Default value: BFGS
        """
        super().__init__(structure, **kwargs)

        self.structure = structure
        self.max_drift = max_drift
        self.rms_drift = rms_drift
        self.max_force = max_force
        self.rms_force = rms_force
        self.max_iter = max_iter
        self.project_name = project_name
        self.optimizer = optimizer
        self.override_default_params = override_default_params or {}
        self.kwargs = kwargs

        global_section = Global(project_name=project_name, run_type="GEO_OPT")

        geo_opt_params = {
            "TYPE": Keyword("TYPE", "MINIMIZATION"),
            "MAX_DR": Keyword("MAX_DR", max_drift),
            "MAX_FORCE": Keyword("MAX_FORCE", max_force),
            "RMS_DR": Keyword("RMS_DR", rms_drift),
            "RMS_FORCE": Keyword("RMS_FORCE", rms_force),
            "MAX_ITER": Keyword("MAX_ITER", max_iter),
            "OPTIMIZER": Keyword("OPTIMIZER", optimizer),
        }
        geo_opt = Section("GEO_OPT", subsections={}, keywords=geo_opt_params)

        if optimizer.upper() == "CG":
            ls = Section("LINE_SEARCH", subsections={}, keywords={"TYPE": Keyword("TYPE", "2PNT")})
            cg = Section("CG", subsections={"LINE_SEARCH": ls}, keywords={})
            geo_opt.insert(cg)

        self.activate_motion()
        self["MOTION"].insert(geo_opt)
        self.insert(global_section)
        self.modify_dft_print_iters(0, add_last="numeric")
        self.update(self.override_default_params)


# TODO Add cell opt convergence criteria
class CellOptSet(DftSet):

    """
    CP2K input set containing the basic settings for performing geometry optimization. Values are all cp2k
    defaults, and should be good for most systems of interest.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        project_name: str = "CellOpt",
        override_default_params: dict = None,
        **kwargs,
    ):
        """
        Args:
            structure: Pymatgen structure object
            max_drift: Convergence criterion for the maximum geometry change between the current and the
                last optimizer iteration. This keyword cannot be repeated and it expects precisely one real.
                Default value: 3.00000000E-003
                Default unit: [bohr]
            max_force (float): Convergence criterion for the maximum force component of the current configuration.
                This keyword cannot be repeated and it expects precisely one real.
                Default value: 4.50000000E-004
                Default unit: [bohr^-1*hartree]
            max_iter (int): Specifies the maximum number of geometry optimization steps.
                One step might imply several force evaluations for the CG and LBFGS optimizers.
                This keyword cannot be repeated and it expects precisely one integer.
                Default value: 200
            optimizer (str): Specify which method to use to perform a geometry optimization.
                This keyword cannot be repeated and it expects precisely one keyword. BFGS is a
                quasi-newtonian method, and will best for "small" systems near the minimum. LBFGS
                is a limited memory version that can be used for "large" (>1000 atom) systems when
                efficiency outweighs robustness. CG is more robust, especially when you are far from
                the minimum, but it slower.
                Default value: BFGS
        """
        super().__init__(structure, **kwargs)

        self.structure = structure
        self.project_name = project_name
        self.override_default_params = override_default_params or {}
        self.kwargs = kwargs
        global_section = Global(project_name=project_name, run_type="CELL_OPT")
        self.insert(global_section)
        self.activate_motion()
        self.modify_dft_print_iters(0, add_last="numeric")
        self.update(self.override_default_params)


class HybridStaticSet(StaticSet):

    """
    Static calculation using hybrid DFT with the ADMM formalism in Cp2k.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        hybrid_functional: str = "PBE0",
        hf_fraction: float = 0.25,
        project_name: str = "Hybrid-Static",
        gga_x_fraction: float = 0.75,
        gga_c_fraction: float = 1,
        potential_type: str = None,
        scale_coulomb: float = 1,
        scale_gaussian: float = 1,
        scale_longrange: float = 1,
        override_default_params: dict = None,
        max_memory: int = 2000,
        cutoff_radius: float = 8.0,
        omega: float = 0.2,
        aux_basis: dict | None = None,
        admm: bool = True,
        eps_schwarz: float = 1e-6,
        eps_schwarz_forces: float = 1e-6,
        screen_on_initial_p: bool = True,
        screen_p_forces: bool = True,
        **kwargs,
    ):
        """
        Args:
            structure: pymatgen structure object
            method: hybrid dft method to use (currently select between HSE06 and PBE0)
            hf_fraction: percentage of exact HF to mix-in
            project_name: what to call this project
            gga_x_fraction: percentage of gga exchange to use
            gga_c_fraction: percentage of gga correlation to use
            override_default_params: override settings (see above).
        """
        super().__init__(structure, project_name=project_name, **kwargs)

        self.structure = structure
        self.hybrid_functional = hybrid_functional
        self.hf_fraction = hf_fraction
        self.project_name = project_name
        self.potential_type = potential_type
        self.gga_x_fraction = gga_x_fraction
        self.gga_c_fraction = gga_c_fraction
        self.scale_coulomb = scale_coulomb
        self.scale_gaussian = scale_gaussian
        self.scale_longrange = scale_longrange
        self.override_default_params = override_default_params or {}
        self.max_memory = max_memory
        self.cutoff_radius = cutoff_radius
        self.omega = omega
        self.aux_basis = aux_basis
        self.admm = admm
        self.eps_schwarz = eps_schwarz
        self.eps_schwarz_forces = eps_schwarz_forces
        self.screen_on_initial_p = screen_on_initial_p
        self.screen_p_forces = screen_p_forces
        self.kwargs = kwargs

        self.activate_hybrid(
            hybrid_functional=self.hybrid_functional,
            hf_fraction=self.hf_fraction,
            potential_type=self.potential_type,
            gga_x_fraction=self.gga_x_fraction,
            gga_c_fraction=self.gga_c_fraction,
            scale_coulomb=self.scale_coulomb,
            scale_gaussian=self.scale_gaussian,
            scale_longrange=self.scale_longrange,
            max_memory=self.max_memory,
            cutoff_radius=cutoff_radius,
            omega=self.omega,
            aux_basis=self.aux_basis,
            admm=self.admm,
            eps_schwarz=self.eps_schwarz,
            eps_schwarz_forces=self.eps_schwarz_forces,
            screen_on_initial_p=self.screen_on_initial_p,
            screen_p_forces=self.screen_p_forces,
        )
        self.update(self.override_default_params)


class HybridRelaxSet(RelaxSet):

    """
    Static calculation using hybrid DFT with the ADMM formalism in Cp2k.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        hybrid_functional: str = "PBE0",
        hf_fraction: float = 0.25,
        project_name: str = "Hybrid-Static",
        potential_type: str = None,
        gga_x_fraction: float = 0.75,
        gga_c_fraction: float = 1,
        scale_coulomb: float = 1,
        scale_gaussian: float = 1,
        scale_longrange: float = 1,
        override_default_params: dict = None,
        max_memory: int = 2000,
        cutoff_radius: float = 8.0,
        omega: float = 0.2,
        aux_basis: dict | None = None,
        admm: bool = True,
        eps_schwarz: float = 1e-6,
        eps_schwarz_forces: float = 1e-6,
        screen_on_initial_p: bool = True,
        screen_p_forces: bool = True,
        **kwargs,
    ):
        """
        Args:
            structure: pymatgen structure object
            method: hybrid dft method to use (currently select between HSE06 and PBE0)
            hf_fraction: percentage of exact HF to mix-in
            project_name: what to call this project
            gga_x_fraction: percentage of gga exchange to use
            gga_c_fraction: percentage of gga correlation to use
            override_default_params: override settings (see above).
        """
        super().__init__(structure, project_name=project_name, **kwargs)

        self.structure = structure
        self.hybrid_functional = hybrid_functional
        self.hf_fraction = hf_fraction
        self.project_name = project_name
        self.potential_type = potential_type
        self.gga_x_fraction = gga_x_fraction
        self.gga_c_fraction = gga_c_fraction
        self.scale_coulomb = scale_coulomb
        self.scale_gaussian = scale_gaussian
        self.scale_longrange = scale_longrange
        self.override_default_params = override_default_params or {}
        self.max_memory = max_memory
        self.cutoff_radius = cutoff_radius
        self.omega = omega
        self.aux_basis = aux_basis
        self.admm = admm
        self.eps_schwarz = eps_schwarz
        self.eps_schwarz_forces = eps_schwarz_forces
        self.screen_on_initial_p = screen_on_initial_p
        self.screen_p_forces = screen_p_forces
        self.kwargs = kwargs

        self.activate_hybrid(
            hybrid_functional=self.hybrid_functional,
            hf_fraction=self.hf_fraction,
            potential_type=self.potential_type,
            gga_x_fraction=self.gga_x_fraction,
            gga_c_fraction=self.gga_c_fraction,
            scale_coulomb=self.scale_coulomb,
            scale_gaussian=self.scale_gaussian,
            scale_longrange=self.scale_longrange,
            max_memory=self.max_memory,
            cutoff_radius=cutoff_radius,
            omega=self.omega,
            aux_basis=self.aux_basis,
            admm=self.admm,
            eps_schwarz=self.eps_schwarz,
            eps_schwarz_forces=self.eps_schwarz_forces,
            screen_on_initial_p=self.screen_on_initial_p,
            screen_p_forces=self.screen_p_forces,
        )
        self.update(self.override_default_params)


class HybridCellOptSet(CellOptSet):

    """
    Static calculation using hybrid DFT with the ADMM formalism in Cp2k.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        hybrid_functional: str = "PBE0",
        hf_fraction: float = 0.25,
        project_name: str = "Hybrid-Static",
        potential_type: str = None,
        gga_x_fraction: float = 0.75,
        gga_c_fraction: float = 1,
        scale_coulomb: float = 1,
        scale_gaussian: float = 1,
        scale_longrange: float = 1,
        override_default_params: dict = None,
        max_memory: int = 2000,
        cutoff_radius: float = 8.0,
        omega: float = 0.2,
        aux_basis: dict | None = None,
        admm: bool = True,
        eps_schwarz: float = 1e-6,
        eps_schwarz_forces: float = 1e-6,
        screen_on_initial_p: bool = True,
        screen_p_forces: bool = True,
        **kwargs,
    ):
        """
        Args:
            structure: pymatgen structure object
            method: hybrid dft method to use (currently select between HSE06 and PBE0)
            hf_fraction: percentage of exact HF to mix-in
            project_name: what to call this project
            gga_x_fraction: percentage of gga exchange to use
            gga_c_fraction: percentage of gga correlation to use
            override_default_params: override settings (see above).
        """
        super().__init__(structure, project_name=project_name, **kwargs)

        self.structure = structure
        self.hybrid_functional = hybrid_functional
        self.hf_fraction = hf_fraction
        self.project_name = project_name
        self.potential_type = potential_type
        self.gga_x_fraction = gga_x_fraction
        self.gga_c_fraction = gga_c_fraction
        self.scale_coulomb = scale_coulomb
        self.scale_gaussian = scale_gaussian
        self.scale_longrange = scale_longrange
        self.override_default_params = override_default_params or {}
        self.max_memory = max_memory
        self.cutoff_radius = cutoff_radius
        self.omega = omega
        self.aux_basis = aux_basis
        self.admm = admm
        self.eps_schwarz = eps_schwarz
        self.eps_schwarz_forces = eps_schwarz_forces
        self.screen_on_initial_p = screen_on_initial_p
        self.screen_p_forces = screen_p_forces
        self.kwargs = kwargs

        self.activate_hybrid(
            hybrid_functional=self.hybrid_functional,
            hf_fraction=self.hf_fraction,
            potential_type=self.potential_type,
            gga_x_fraction=self.gga_x_fraction,
            gga_c_fraction=self.gga_c_fraction,
            scale_coulomb=self.scale_coulomb,
            scale_gaussian=self.scale_gaussian,
            scale_longrange=self.scale_longrange,
            max_memory=self.max_memory,
            cutoff_radius=cutoff_radius,
            omega=self.omega,
            aux_basis=self.aux_basis,
            admm=self.admm,
            eps_schwarz=self.eps_schwarz,
            eps_schwarz_forces=self.eps_schwarz_forces,
            screen_on_initial_p=self.screen_on_initial_p,
            screen_p_forces=self.screen_p_forces,
        )
        self.update(self.override_default_params)
