from __future__ import annotations

import linecache
from abc import ABC, abstractmethod
from collections import Counter
from typing import TYPE_CHECKING

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core import Lattice, Structure

if TYPE_CHECKING:
    from pymatgen.util.typing import PathLike


class LineLocator(MSONable):
    @staticmethod
    def locate_all_lines(file_path: PathLike, content: str):
        """
        Description:
            Locate the line where a certain paragraph of text is located (return all lines)

        Args:
            file_path: str
                Absolute path of file.
            2. content: str
                Contents that needs to be located.
        """
        row_idxs_lst: list[int] = []  # starts from 1, be compatible to linecache package
        row_no: int = 0
        with zopen(file_path, "rt") as f:
            for row_content in f:
                row_no += 1
                if content.upper() in row_content.upper():
                    row_idxs_lst.append(row_no)
        return row_idxs_lst


class ListLocator(MSONable):
    @staticmethod
    def locate_all_lines(strs_lst: list[str], content: str):
        """
        Description:
            Locate the indices in list where a certain paragraph of text is located (return all indices)

        Args:
            strs_lst: List[str]
            content: str
                Contents that needs to be located.
        """
        str_idxs_lst: list[int] = []  # starts from 0, be compatible to list in python
        str_no = -1
        for tmp_str in strs_lst:
            str_no += 1
            if content.upper() in tmp_str.upper():
                str_idxs_lst.append(str_no)
        return str_idxs_lst


class ACExtractorBase(ABC):
    @abstractmethod
    def get_num_atoms(self):
        pass

    @abstractmethod
    def get_lattice(self):
        pass

    @abstractmethod
    def get_types(self):
        pass

    @abstractmethod
    def get_coords(self):
        """
        Returns:
            Fractional coordinates
        """

    @abstractmethod
    def get_magmoms(self):
        pass


class ACExtractor(ACExtractorBase):
    """
    Description:
        Extract information contained in atom.config : number of atoms, lattice, types, frac_coords, magmoms

    @Author: Hanyu Liu
    @email:  domainofbuaa@gmail.com
    """

    def __init__(self, file_path: PathLike):
        """
        Args
            file_path: str
                The absolute path of `atom.config` file
        """
        self.atom_config_path = file_path
        self.num_atoms = self.get_num_atoms()
        self.lattice = self.get_lattice()
        self.types = self.get_types()
        self.coords = self.get_coords()
        self.magmoms = self.get_magmoms()

    def get_num_atoms(self):
        """
        Returns:
            num_atoms: int
        """
        first_row = linecache.getline(self.atom_config_path, 1)
        return int(first_row.split()[0])

    def get_lattice(self):
        """
        Returns:
            lattice: np.ndarray, shape = (9,)
        """
        basis_vectors: list[float] = []
        content: str = "LATTICE"
        idx_row: int = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[0]
        for row_idx in [idx_row + 1, idx_row + 2, idx_row + 3]:
            row_content: list[str] = linecache.getline(self.atom_config_path, row_idx).split()[:3]
            for value in row_content:
                basis_vectors.append(float(value))

        return np.array(basis_vectors)

    def get_types(self):
        """
        Returns:
            atomic_numbers : np.ndarray
                Atomic numbers in order corresponding to sites
        """
        content = "POSITION"
        idx_row = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[0]
        with open(self.atom_config_path) as f:
            atom_config_content = f.readlines()
        atomic_numbers_content = atom_config_content[idx_row : idx_row + self.num_atoms]
        atomic_numbers_lst = [int(row.split()[0]) for row in atomic_numbers_content]  # convert str to int
        return np.array(atomic_numbers_lst)

    def get_coords(self):
        """
        Returns:
            coords: np.ndarray
                Fractional coordinates.
        """
        coords_lst: list[float] = []
        content: str = "POSITION"
        idx_row: int = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[0]
        with open(self.atom_config_path) as f:
            atom_config_content = f.readlines()
        """
        row_content
        -----------
            '29         0.377262291145329         0.128590184800933         0.257759805813488     1  1  1'
        """
        for row_content in atom_config_content[idx_row : idx_row + self.num_atoms]:
            row_content_lst = row_content.split()
            coord_tmp = [float(value) for value in row_content_lst[1:4]]  # convert str to float.
            coords_lst.append(np.array(coord_tmp))
        return np.array(coords_lst).reshape(-1)

    def get_magmoms(self):
        content = "MAGNETIC"
        magnetic_moments_lst = []
        try:  # Error: not containing magmoms info.
            idx_row = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[-1]

            with open(self.atom_config_path) as f:
                atom_config_content = f.readlines()

            magnetic_moments_content = atom_config_content[idx_row : idx_row + self.num_atoms]
            # MAGNETIC
            # 3 0.0 # atomic_number magmom
            # ...
            magnetic_moments_lst = [
                float(tmp_magnetic_moment.split()[-1]) for tmp_magnetic_moment in magnetic_moments_content
            ]
        except Exception:
            magnetic_moments_lst = [0 for _ in range(self.num_atoms)]
        return magnetic_moments_lst


class ACstrExtractor(ACExtractorBase):
    """
    Description:
        Extract information from atom.config file. You can get str by slicing the MOVEMENT.

    @Author: Hanyu Liu
    @email:  domainofbuaa@gmail.com
    """

    def __init__(self, atom_config_str: str):
        self.atom_config_str = atom_config_str
        self.strs_lst = self.atom_config_str.split("\n")
        self.num_atoms = self.get_num_atoms()

    def get_num_atoms(self):
        return int(self.strs_lst[0].split()[0].strip())

    def get_lattice(self):
        """
        Returns:
            lattice: np.ndarray, shape=(9,)
        """
        basis_vectors_lst = []
        aim_content = "LATTICE"
        aim_idx = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)[0]
        for idx_str in [aim_idx + 1, aim_idx + 2, aim_idx + 3]:
            # ['0.8759519000E+01', '0.0000000000E+00', '0.0000000000E+00']
            str_lst = self.strs_lst[idx_str].split()[:3]
            for tmp_str in str_lst:
                basis_vectors_lst.append(float(tmp_str))  # convert str to float
        return np.array(basis_vectors_lst)

    def get_types(self):
        """
        Returns:
            types: np.ndarray
        """
        aim_content = "POSITION"
        aim_idx = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)[0]
        strs_lst = self.strs_lst[aim_idx + 1 : aim_idx + self.num_atoms + 1]
        atomic_numbers_lst = [int(entry.split()[0]) for entry in strs_lst]
        return np.array(atomic_numbers_lst)

    def get_coords(self):
        """
        Returns:
            coords: np.ndarray, shape=(num_atoms*3, )
        """
        coords_lst = []
        aim_content = "POSITION"
        aim_idx = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)[0]
        for tmp_str in self.strs_lst[aim_idx + 1 : aim_idx + self.num_atoms + 1]:
            # ['14', '0.751401861790384', '0.501653718883189', '0.938307102003243', '1', '1', '1']
            tmp_strs_lst = tmp_str.split()
            tmp_coord = [float(value) for value in tmp_strs_lst[1:4]]  # convert str to float
            coords_lst.append(tmp_coord)
        return np.array(coords_lst).reshape(-1)

    def get_magmoms(self):
        """
        Returns:
            magmoms: np.ndarray
        """
        magnetic_moments_lst = []
        aim_content = "MAGNETIC"
        aim_idxs = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)
        if len(aim_idxs) == 0:
            magnetic_moments_lst = [0 for _ in range(self.num_atoms)]
        else:
            aim_idx = aim_idxs[0]
            magnetic_moments_content = self.strs_lst[aim_idx + 1 : aim_idx + self.num_atoms + 1]
            magnetic_moments_lst = [
                float(tmp_magnetic_moment.split()[-1]) for tmp_magnetic_moment in magnetic_moments_content
            ]
        return np.array(magnetic_moments_lst)

    def get_etot(self):
        """
        [' 216 atoms', 'Iteration (fs) =    0.3000000000E+01',
        ' Etot', 'Ep', 'Ek (eV) =   -0.2831881714E+05  -0.2836665392E+05   0.4783678177E+02',
        ' SCF =     7']
        """
        strs_lst = self.strs_lst[0].split(",")
        aim_index = ListLocator.locate_all_lines(strs_lst=strs_lst, content="EK (EV) =")[0]
        # strs_lst[aim_index].split() :
        #   ['Ek', '(eV)', '=', '-0.2831881714E+05', '-0.2836665392E+05', '0.4783678177E+02']
        return np.array([float(strs_lst[aim_index].split()[3].strip())])

    def get_eatoms(self):
        """
        Description:
            When turn on `ENERGY DEPOSITION`, PWmat will output energy per atom.
        Returns:
            eatoms: np.ndarray | None
        """
        eatoms_lst = []
        aim_content = "Atomic-Energy, ".upper()
        aim_idxs = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)
        if len(aim_idxs) == 0:
            return None
        aim_idx = aim_idxs[0]
        for tmp_str in self.strs_lst[aim_idx + 1 : aim_idx + self.num_atoms + 1]:
            """
            Atomic-Energy, Etot(eV),E_nonloc(eV),Q_atom:dE(eV)=  -0.1281163115E+06
            14   0.6022241483E+03    0.2413350871E+02    0.3710442365E+01
            """
            eatoms_lst.append(float(tmp_str.split()[1]))
        return np.array(eatoms_lst)

    def get_fatoms(self):
        """
        Returns:
            fatoms: np.ndarray
        """
        forces_lst = []
        aim_content = "Force".upper()
        aim_idx = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)[0]

        for tmp_str in self.strs_lst[aim_idx + 1 : aim_idx + self.num_atoms + 1]:
            # ['14', '0.089910342901203', '0.077164252174742', '0.254144099204679']
            tmp_str_lst = tmp_str.split()
            tmp_forces = [float(value) for value in tmp_str_lst[1:4]]
            forces_lst.append(tmp_forces)
        return -np.array(forces_lst).reshape(-1)

    def get_virial(self):
        """
        Returns:
            1. virial_tensor: np.ndarray | None
        """
        virial_tensor = []

        ### Step 1. 得到所有原子的原子序数、坐标
        aim_content = "LATTICE"
        aim_idx = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)[0]

        for tmp_idx in [aim_idx + 1, aim_idx + 2, aim_idx + 3]:
            # tmp_strs_lst =
            # ['0.8759519000E+01', '0.0000000000E+00', '0.0000000000E+00',
            # 'stress', '(eV):', '0.115558E+02', '0.488108E+01', '0.238778E+01']
            tmp_strs_lst = self.strs_lst[tmp_idx].split()
            tmp_aim_row_lst = ListLocator.locate_all_lines(strs_lst=tmp_strs_lst, content="STRESS")
            if len(tmp_aim_row_lst) == 0:
                return None

        ### Step 2. 获取维里张量
        for tmp_idx in [aim_idx + 1, aim_idx + 2, aim_idx + 3]:
            # tmp_str_lst = ['0.120972E+02', '0.483925E+01', '0.242063E+01']
            tmp_str_lst = self.strs_lst[tmp_idx].split()[-3:]

            virial_tensor.append(float(tmp_str_lst[0]))
            virial_tensor.append(float(tmp_str_lst[1]))
            virial_tensor.append(float(tmp_str_lst[2]))

        return np.array(virial_tensor)


class AtomConfig(MSONable):
    def __init__(self, structure: Structure, sort_structure: bool = False):
        """
        Args:
            structure (Structure): Structure object
            sort_structure (Optional[bool]): Whether to sort the structure. Useful if species
                are not grouped properly together. Defaults to False.

        Returns:
            void

        @Author: Hanyu Liu
        @email:  domainofbuaa@gmail.com
        """
        self.structure: Structure = structure
        if sort_structure:
            pass
        elements_counter = dict(sorted(Counter(self.structure.species).items()))
        true_names = [f"{tmp_key}{tmp_value}" for (tmp_key, tmp_value) in elements_counter.items()]
        self.true_names = "".join(true_names)

    def __repr__(self):
        return self.get_str()

    def __str__(self):
        return self.get_str()

    @classmethod
    def from_str(cls, data: str, mag: bool = True):
        """
        Reads a atom.config from a string

        Args:
            data: string containing atom.config data

        Returns:
            AtomConfig object
        """
        acextractor = ACstrExtractor(atom_config_str=data)
        properties: dict[str, float] = {}
        structure = Structure(
            lattice=acextractor.get_lattice(),
            species=acextractor.get_types(),
            coords=acextractor.get_coords().reshape(-1, 3),
            coords_are_cartesian=False,
            properties=properties,
        )
        if mag:
            magmoms = acextractor.get_magmoms()
            for ii, tmp_site in enumerate(structure.sites):
                tmp_site.properties.update({"magmom": magmoms[ii]})

        return cls(
            structure,
        )

    @classmethod
    def from_file(cls, filename: str, mag: bool = True):
        with zopen(filename, "rt") as f:
            return cls.from_str(data=f.read(), mag=mag)

    @classmethod
    def from_dict(cls, d: dict):
        """
        Args:
            d: dict
        """
        return cls(
            Structure.from_dict(d["structure"]),
        )

    def get_str(self):
        """
        Returns:
            String representation of atom.config
        """
        # This corrects for VASP really annoying bug of crashing on lattices
        # which have triple product < 0. We will just invert the lattice
        # vectors.
        latt = self.structure.lattice
        if np.linalg.det(latt.matrix) < 0:
            latt = Lattice(-latt.matrix)

        lines: list[str] = []
        lines.append(f"\t{self.structure.num_sites} atoms\n")
        lines.append("Lattice vector\n")
        for ii in range(3):
            lines.append(f"{latt.matrix[ii][0]:>15f}{latt.matrix[ii][1]:>15f}{latt.matrix[ii][2]:>15f}\n")
        lines.append("Position, move_x, move_y, move_z\n")
        for ii in range(self.structure.num_sites):
            lines.append(f"{int(self.structure.species[ii].Z):>4d}")
            lines.append(f"{self.structure.frac_coords[ii][0]:>15f}")
            lines.append(f"{self.structure.frac_coords[ii][1]:>15f}")
            lines.append(f"{self.structure.frac_coords[ii][2]:>15f}")
            lines.append("   1   1   1\n")
        if "magmom" in self.structure.sites[0].properties:
            lines.append("MAGNETIC\n")
            for _, tmp_site in enumerate(self.structure.sites):
                lines.append(f"{int(tmp_site.specie.Z):>4d}{tmp_site.properties['magmom']:>15f}\n")
        return "".join(lines)

    def write_file(self, filename: PathLike, **kwargs):
        with zopen(filename, "wt") as f:
            f.write(self.get_str(**kwargs))

    def as_dict(self):
        """
        Returns:
            dict
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "structure": self.structure.as_dict(),
            "true_names": self.true_names,
        }
