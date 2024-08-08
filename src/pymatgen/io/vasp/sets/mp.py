# ruff: noqa: PGH003

from __future__ import annotations

import warnings
from copy import deepcopy
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from monty.serialization import loadfn

from pymatgen.core import Element, Species, Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets.base import (
    MODULE_DIR,
    UserPotcarFunctional,
    VaspInputSet,
    _get_ispin,
    _get_nedos,
    _load_yaml_config,
)
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, Literal

    from pymatgen.util.typing import Vector3D


@dataclass
class MPRelaxSet(VaspInputSet):
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
        **kwargs: Keywords supported by VaspInputSet.
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
class MPScanRelaxSet(VaspInputSet):
    """Write a relaxation input set using the accurate and numerically
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
        **kwargs: Keywords supported by VaspInputSet.

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
    auto_kspacing: bool = True
    user_potcar_functional: UserPotcarFunctional = "PBE_54"
    auto_ismear: bool = True
    CONFIG = _load_yaml_config("MPSCANRelaxSet")
    _valid_potcars: Sequence[str] | None = ("PBE_52", "PBE_54")

    def __post_init__(self):
        super().__post_init__()
        if self.vdw and self.vdw != "rvv10":
            warnings.warn("Use of van der waals functionals other than rVV10 with SCAN is not supported at this time")
            # delete any vdw parameters that may have been added to the INCAR
            vdw_par = loadfn(str(MODULE_DIR / "vdW_parameters.yaml"))
            for key in vdw_par[self.vdw]:
                self._config_dict["INCAR"].pop(key, None)


@dataclass
class MPMetalRelaxSet(VaspInputSet):
    """
    Implementation of VaspInputSet utilizing parameters in the public
    Materials Project, but with tuning for metals. Key things are a denser
    k point density, and a.
    """

    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
        return {"ISMEAR": 1, "SIGMA": 0.2}

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Updates to the kpoints configuration for this calculation type."""
        return {"reciprocal_density": 200}


@dataclass
class MPHSERelaxSet(VaspInputSet):
    """Same as the MPRelaxSet, but with HSE parameters."""

    CONFIG = _load_yaml_config("MPHSERelaxSet")


@dataclass
class MPStaticSet(VaspInputSet):
    """Create input files for a static calculation.

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
        **kwargs: Keywords supported by MPRelaxSet.
    """

    lepsilon: bool = False
    lcalcpol: bool = False
    reciprocal_density: int = 100
    small_gap_multiply: tuple[float, float] | None = None
    inherit_incar: bool = True
    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
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
        """Updates to the kpoints configuration for this calculation type."""
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
class MPScanStaticSet(MPScanRelaxSet):
    """Create input files for a static calculation using the accurate and numerically
    efficient r2SCAN variant of the Strongly Constrained and Appropriately Normed
    (SCAN) metaGGA functional.

    Args:
        structure (Structure): Structure from previous run.
        bandgap (float): Bandgap of the structure in eV. The bandgap is used to
            compute the appropriate k-point density and determine the smearing settings.
        lepsilon (bool): Whether to add static dielectric calculation
        lcalcpol (bool): Whether to turn on evaluation of the Berry phase approximations
            for electronic polarization.
        **kwargs: Keywords supported by MPScanRelaxSet.
    """

    lepsilon: bool = False
    lcalcpol: bool = False
    inherit_incar: bool = True
    auto_kspacing: bool = True

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
        updates: dict[str, Any] = {
            "LREAL": False,
            "NSW": 0,
            "LORBIT": 11,
            "LVHAR": True,
            "ISMEAR": -5,
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
class MPHSEBSSet(VaspInputSet):
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
        **kwargs: Keywords supported by VaspInputSet.
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
        """Updates to the kpoints configuration for this calculation type."""
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
        """Updates to the INCAR config for this calculation type."""
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
class MPNonSCFSet(VaspInputSet):
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
        **kwargs: Keywords supported by MPRelaxSet.
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

        mode = self.mode = self.mode.lower()

        valid_modes = ("line", "uniform", "boltztrap")
        if mode not in valid_modes:
            raise ValueError(
                f"Invalid {mode=}. Supported modes for NonSCF runs are {', '.join(map(repr, valid_modes))}"
            )

        if (mode != "uniform" or self.nedos < 2000) and self.optics:
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
        """Updates to the INCAR config for this calculation type."""
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
        """Updates to the kpoints configuration for this calculation type."""
        factor = 1.0
        if self.bandgap is not None and self.small_gap_multiply and self.bandgap <= self.small_gap_multiply[0]:
            factor = self.small_gap_multiply[1]

        if self.mode == "line":
            return {"line_density": self.kpoints_line_density * factor}

        if self.mode == "boltztrap":
            return {"explicit": True, "reciprocal_density": self.reciprocal_density * factor}

        return {"reciprocal_density": self.reciprocal_density * factor}


@dataclass
class MPSOCSet(VaspInputSet):
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
        **kwargs: Keywords supported by VaspInputSet.
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
                "The structure must have the 'magmom' site property and each magnetic "
                "moment value must have 3 components. e.g. magmom = [0,0,2]"
            )

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
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
        """Updates to the kpoints configuration for this calculation type."""
        factor = 1.0
        if self.bandgap is not None and self.small_gap_multiply and self.bandgap <= self.small_gap_multiply[0]:
            factor = self.small_gap_multiply[1]
        return {"reciprocal_density": self.reciprocal_density * factor}

    @VaspInputSet.structure.setter  # type: ignore
    def structure(self, structure: Structure | None) -> None:
        if structure is not None:
            if self.magmom:
                structure = structure.copy(site_properties={"magmom": self.magmom})

            # magmom has to be 3D for SOC calculation.
            if hasattr(structure[0], "magmom"):
                if not isinstance(structure[0].magmom, list):
                    # project magmom to z-axis
                    structure = structure.copy(site_properties={"magmom": [[0, 0, site.magmom] for site in structure]})
            else:
                raise ValueError("Neither the previous structure has magmom property nor magmom provided")

        VaspInputSet.structure.fset(self, structure)  # type: ignore


@dataclass
class MPNMRSet(VaspInputSet):
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
        **kwargs: Keywords supported by MPRelaxSet.
    """

    mode: Literal["cs", "efg"] = "cs"
    isotopes: list = field(default_factory=list)
    reciprocal_density: int = 100
    small_gap_multiply: tuple[float, float] | None = None
    inherit_incar: bool = True
    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
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
        """Updates to the kpoints configuration for this calculation type."""
        factor = 1.0
        if self.bandgap is not None and self.small_gap_multiply and self.bandgap <= self.small_gap_multiply[0]:
            factor = self.small_gap_multiply[1]
        return {"reciprocal_density": self.reciprocal_density * factor}


@dataclass
class MPMDSet(VaspInputSet):
    """
    This a modified version of the old MITMDSet pre 2018/03/12.

    This set serves as the basis for the amorphous skyline paper.

    (1) Aykol, M.; Dwaraknath, S. S.; Sun, W.; Persson, K. A. Thermodynamic
        Limit for Synthesis of Metastable Inorganic Materials. Sci. Adv. 2018,
        4 (4).

    Class for writing a VASP MD run. This DOES NOT do multiple stage runs.
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
        **kwargs: Other kwargs supported by VaspInputSet.
    """

    start_temp: float = 0.0
    end_temp: float = 300.0
    nsteps: int = 1000
    time_step: float | None = None
    spin_polarized: bool = False
    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
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
        """Updates to the kpoints configuration for this calculation type."""
        return Kpoints.gamma_automatic()


@dataclass
class MPAbsorptionSet(VaspInputSet):
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
        **kwargs: All kwargs supported by VaspInputSet. Typically, user_incar_settings is a
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
        """Validate settings."""
        super().__post_init__()
        self.mode = self.mode.upper()
        if self.mode not in MPAbsorptionSet.SUPPORTED_MODES:
            raise ValueError(f"{self.mode} not one of the support modes : {MPAbsorptionSet.SUPPORTED_MODES}")

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Updates to the kpoints configuration for this calculation type.

        Generate gamma center k-points mesh grid for optical calculation. It is not
        mandatory for 'ALGO = Exact', but is requested by 'ALGO = CHI' calculation.
        """
        return {"reciprocal_density": self.reciprocal_density}

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
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
