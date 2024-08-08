from __future__ import annotations

from dataclasses import dataclass
from itertools import chain
from pathlib import Path

import numpy as np

from pymatgen.core import Structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.inputs import Kpoints, Poscar
from pymatgen.io.vasp.sets.base import VaspInputSet, _load_yaml_config
from pymatgen.util.due import Doi, due


@due.dcite(
    Doi("10.1016/j.commatsci.2011.02.023"),
    description="A high-throughput infrastructure for density functional theory calculations",
)
@dataclass
class MITRelaxSet(VaspInputSet):
    """
    Standard implementation of VaspInputSet utilizing parameters in the MIT
    High-throughput project.
    The parameters are chosen specifically for a high-throughput project,
    which means in general pseudopotentials with fewer electrons were chosen.

    Args:
        structure (Structure): The Structure to create inputs for. If None, the input
            set is initialized without a Structure but one must be set separately before
            the inputs are generated.
        **kwargs: Keywords supported by VaspInputSet.

    Please refer::

        A Jain, G. Hautier, C. Moore, S. P. Ong, C. Fischer, T. Mueller,
        K. A. Persson, G. Ceder. A high-throughput infrastructure for density
        functional theory calculations. Computational Materials Science,
        2011, 50(8), 2295-2310. doi:10.1016/j.commatsci.2011.02.023
    """

    CONFIG = _load_yaml_config("MITRelaxSet")


class MITNEBSet(VaspInputSet):
    """Write NEB inputs.

    Note that EDIFF is not on a per atom basis for this input set.
    """

    def __init__(self, structures, unset_encut=False, **kwargs) -> None:
        """
        Args:
            structures: List of Structure objects.
            unset_encut (bool): Whether to unset ENCUT.
            **kwargs: Other kwargs supported by VaspInputSet.
        """
        if len(structures) < 3:
            raise ValueError(f"You need at least 3 structures for an NEB, got {len(structures)}")
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
            for idx, site in enumerate(s):
                translate = np.round(prev[idx].frac_coords - site.frac_coords)
                if np.any(np.abs(translate) > 0.5):
                    s.translate_sites([idx], translate, to_unit_cell=False)
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

        for idx, poscar in enumerate(self.poscars):
            d = output_dir / str(idx).zfill(2)
            if not d.exists():
                d.mkdir(parents=True)
            poscar.write_file(str(d / "POSCAR"))
            if write_cif:
                poscar.structure.to(filename=str(d / f"{idx}.cif"))
        if write_endpoint_inputs:
            end_point_param = MITRelaxSet(self.structures[0], user_incar_settings=self.user_incar_settings)

            for image in ["00", str(len(self.structures) - 1).zfill(2)]:
                end_point_param.incar.write_file(str(output_dir / image / "INCAR"))
                end_point_param.kpoints.write_file(str(output_dir / image / "KPOINTS"))
                end_point_param.potcar.write_file(str(output_dir / image / "POTCAR"))
        if write_path_cif:
            sites = {
                PeriodicSite(site.species, site.frac_coords, self.structures[0].lattice)
                for site in chain(*(struct for struct in self.structures))
            }
            neb_path = Structure.from_sites(sorted(sites))
            neb_path.to(filename=f"{output_dir}/path.cif")


@dataclass
class MITMDSet(VaspInputSet):
    """Write a VASP MD run. This DOES NOT do multiple stage runs.

    Args:
        structure (Structure): Input structure.
        start_temp (float): Starting temperature.
        end_temp (float): Final temperature.
        nsteps (int): Number of time steps for simulations. NSW parameter.
        time_step (float): The time step for the simulation. The POTIM
            parameter. Defaults to 2fs.
        spin_polarized (bool): Whether to do spin polarized calculations.
            The ISPIN parameter. Defaults to False.
        **kwargs: Other kwargs supported by VaspInputSet.
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
        """Updates to the INCAR config for this calculation type."""
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
        """Updates to the kpoints configuration for this calculation type."""
        return Kpoints.gamma_automatic()
