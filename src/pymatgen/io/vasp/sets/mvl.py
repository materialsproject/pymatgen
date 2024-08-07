# ruff: noqa: PGH003

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from monty.json import MSONable

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets.base import UserPotcarFunctional, VaspInputSet, _load_yaml_config
from pymatgen.io.vasp.sets.mit import MITRelaxSet
from pymatgen.io.vasp.sets.mp import MPRelaxSet
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self


@dataclass
class MVLGWSet(VaspInputSet):
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
        **kwargs: All kwargs supported by VaspInputSet. Typically,
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
        self.mode = mode = self.mode.upper()

        if mode not in MVLGWSet.SUPPORTED_MODES:
            raise ValueError(f"Invalid {mode=}, supported modes are {', '.join(map(repr, MVLGWSet.SUPPORTED_MODES))}")

    @property
    def kpoints_updates(self) -> dict:
        """Updates to the kpoints configuration for this calculation type."""
        # Generate gamma center k-points mesh grid for GW calc, which is requested
        # by GW calculation.
        return {"reciprocal_density": self.reciprocal_density}

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
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
    def from_prev_calc(cls, prev_calc_dir: str, mode: str = "DIAG", **kwargs) -> Self:
        """Generate a set of VASP input files for GW or BSE calculations from a
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
        input_set = cls(None, mode=mode, **kwargs)
        return input_set.override_from_prev_calc(prev_calc_dir=prev_calc_dir)


@dataclass
class MVLSlabSet(VaspInputSet):
    """Write a set of slab vasp runs, including both slabs (along the c direction)
    and orient unit cells (bulk), to ensure the same KPOINTS, POTCAR and INCAR criterion.

    Args:
        structure: Structure
        k_product: default to 50, kpoint number * length for a & b
            directions, also for c direction in bulk calculations
        bulk:
        auto_dipole:
        set_mix:
        sort_structure:
        **kwargs: Other kwargs supported by VaspInputSet.
    """

    k_product: int = 50
    bulk: bool = False
    auto_dipole: bool = False
    set_mix: bool = True
    CONFIG = MPRelaxSet.CONFIG

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
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
        """Updates to the kpoints configuration for this calculation type.

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
            verbosity (int): Verbosity of dict. e.g. whether to include Structure.

        Returns:
            dict: MSONable MVLGBSet representation.
        """
        dct = MSONable.as_dict(self)
        if verbosity == 1:
            dct.pop("structure", None)
        return dct


@dataclass
class MVLGBSet(VaspInputSet):
    """Write a vasp input files for grain boundary calculations, slab or bulk.

    Args:
        structure (Structure): provide the structure
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
        """k_product is kpoint number * length for a & b directions, also for c direction
        in bulk calculations Automatic mesh & Gamma is the default setting.
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
        """Updates to the INCAR config for this calculation type."""
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
class MVLRelax52Set(VaspInputSet):
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
        **kwargs: Other kwargs supported by VaspInputSet.
    """

    user_potcar_functional: UserPotcarFunctional = "PBE_52"
    CONFIG = _load_yaml_config("MVLRelax52Set")
    _valid_potcars: Sequence[str] | None = ("PBE_52", "PBE_54")


@due.dcite(
    Doi("10.1149/2.0061602jes"),
    description="Elastic Properties of Alkali Superionic Conductor Electrolytes from First Principles Calculations",
)
class MVLElasticSet(VaspInputSet):
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
        """Updates to the INCAR config for this calculation type."""
        return {"IBRION": 6, "NFREE": 2, "POTIM": self.potim, "NPAR": None}


@dataclass
class MVLNPTMDSet(VaspInputSet):
    """Write a VASP MD run in NPT ensemble.

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
        """Updates to the INCAR config for this calculation type."""
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
        enmax = [self.potcar[i].keywords["ENMAX"] for i in range(self.structure.ntypesp)]  # type: ignore[union-attr]
        updates["ENCUT"] = max(enmax) * 1.5
        return updates

    @property
    def kpoints_updates(self) -> Kpoints | dict:
        """Updates to the kpoints configuration for this calculation type."""
        return Kpoints.gamma_automatic()


@dataclass
class MVLScanRelaxSet(VaspInputSet):
    """Write a relax input set using Strongly Constrained and
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
        **kwargs: Other kwargs supported by VaspInputSet.
    """

    user_potcar_functional: UserPotcarFunctional = "PBE_52"
    _valid_potcars: Sequence[str] | None = ("PBE_52", "PBE_54")
    CONFIG = MPRelaxSet.CONFIG

    def __post_init__(self):
        super().__post_init__()
        if self.user_potcar_functional not in ("PBE_52", "PBE_54"):
            raise ValueError("SCAN calculations require PBE_52 or PBE_54!")

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
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
