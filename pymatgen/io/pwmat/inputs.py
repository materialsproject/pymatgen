from __future__ import annotations

import linecache
from abc import ABC, abstractmethod
from collections import Counter
from typing import TYPE_CHECKING

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.kpath import KPathSeek

if TYPE_CHECKING:
    from pymatgen.util.typing import PathLike

__author__ = "Hanyu Liu"
__email__ = "domainofbuaa@gmail.com"
__date__ = "2023-12-28"


class LineLocator(MSONable):
    @staticmethod
    def locate_all_lines(file_path: PathLike, content: str) -> list[int]:
        """Locate the line where a certain paragraph of text is located (return all lines)

        Args:
            file_path (str): Absolute path to file.
            content (str): Contents that needs to be located.
        """
        row_idxs: list[int] = []  # starts from 1 to be compatible with linecache package
        row_no: int = 0
        with zopen(file_path, "rt") as f:
            for row_content in f:
                row_no += 1
                if content.upper() in row_content.upper():
                    row_idxs.append(row_no)
        return row_idxs


class ListLocator(MSONable):
    @staticmethod
    def locate_all_lines(strs_lst: list[str], content: str) -> list[int]:
        """Locate the indices in list where a certain paragraph of text is located (return all indices)

        Args:
            strs_lst (list[str]): List of strings.
            content (str): Contents that needs to be located.
        """
        str_idxs: list[int] = []  # starts from 0 to be compatible with list
        str_no: int = -1
        for tmp_str in strs_lst:
            str_no += 1
            if content.upper() in tmp_str.upper():
                str_idxs.append(str_no)
        return str_idxs


class ACExtractorBase(ABC):
    @abstractmethod
    def get_n_atoms(self):
        pass

    @abstractmethod
    def get_lattice(self):
        pass

    @abstractmethod
    def get_types(self):
        pass

    @abstractmethod
    def get_coords(self) -> np.ndarray:  # fractional coordinates
        pass

    @abstractmethod
    def get_magmoms(self):
        pass


class ACExtractor(ACExtractorBase):
    """Extract information contained in atom.config : number of atoms, lattice, types, frac_coords, magmoms"""

    def __init__(self, file_path: PathLike):
        """
        Args
            file_path: str
                The absolute path of atom.config file
        """
        self.atom_config_path = file_path
        self.n_atoms = self.get_n_atoms()
        self.lattice = self.get_lattice()
        self.types = self.get_types()
        self.coords = self.get_coords()
        self.magmoms = self.get_magmoms()

    def get_n_atoms(self) -> int:
        """
        Returns:
            int: number of atoms
        """
        first_row = linecache.getline(str(self.atom_config_path), 1)
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

    def get_types(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: Atomic numbers in order corresponding to sites
        """
        content = "POSITION"
        idx_row = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[0]
        with open(self.atom_config_path) as f:
            atom_config_content = f.readlines()
        atomic_numbers_content = atom_config_content[idx_row : idx_row + self.n_atoms]
        atomic_numbers_lst = [int(row.split()[0]) for row in atomic_numbers_content]  # convert str to int
        return np.array(atomic_numbers_lst)

    def get_coords(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: Fractional coordinates.
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
        for row_content in atom_config_content[idx_row : idx_row + self.n_atoms]:
            row_content_lst = row_content.split()
            coord_tmp = [float(value) for value in row_content_lst[1:4]]  # convert str to float.
            coords_lst.append(np.array(coord_tmp))
        return np.array(coords_lst).reshape(-1)

    def get_magmoms(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: The magnetic moments of individual atoms.
        """
        content = "MAGNETIC"
        magnetic_moments_lst = []
        try:  # Error: not containing magmoms info.
            idx_row = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[-1]

            with open(self.atom_config_path) as f:
                atom_config_content = f.readlines()

            magnetic_moments_content = atom_config_content[idx_row : idx_row + self.n_atoms]
            # MAGNETIC
            # 3 0.0 # atomic_number magmom
            # ...
            magnetic_moments_lst = [
                float(tmp_magnetic_moment.split()[-1]) for tmp_magnetic_moment in magnetic_moments_content
            ]
        except Exception:
            magnetic_moments_lst = [0 for _ in range(self.n_atoms)]
        return magnetic_moments_lst


class ACstrExtractor(ACExtractorBase):
    """Extract information from atom.config file. You can get str by slicing the MOVEMENT."""

    def __init__(self, atom_config_str: str):
        self.atom_config_str = atom_config_str
        self.strs_lst = self.atom_config_str.split("\n")
        self.num_atoms = self.get_n_atoms()

    def get_n_atoms(self) -> int:
        """
        Returns:
            int: The number of atoms
        """
        return int(self.strs_lst[0].split()[0].strip())

    def get_lattice(self):
        """
        Returns:
            np.ndarray: Lattice basis vectors of shape=(9,)
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

    def get_types(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: Types of elements.
        """
        aim_content = "POSITION"
        aim_idx = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)[0]
        strs_lst = self.strs_lst[aim_idx + 1 : aim_idx + self.num_atoms + 1]
        atomic_numbers_lst = [int(entry.split()[0]) for entry in strs_lst]
        return np.array(atomic_numbers_lst)

    def get_coords(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: Fractional coordinates of atoms of shape=(num_atoms*3,)
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

    def get_magmoms(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: Atomic magnetic moments.
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

    def get_etot(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: The total energy of the material system.
        """
        # strs_lst:
        #   [' 216 atoms', 'Iteration (fs) =    0.3000000000E+01',
        #    ' Etot', 'Ep', 'Ek (eV) =   -0.2831881714E+05  -0.2836665392E+05   0.4783678177E+02',
        #    ' SCF =     7']
        strs_lst = self.strs_lst[0].split(",")
        aim_index = ListLocator.locate_all_lines(strs_lst=strs_lst, content="EK (EV) =")[0]
        # strs_lst[aim_index].split() :
        #   ['Ek', '(eV)', '=', '-0.2831881714E+05', '-0.2836665392E+05', '0.4783678177E+02']
        return np.array([float(strs_lst[aim_index].split()[3].strip())])

    def get_eatoms(self) -> np.ndarray | None:
        """
        Returns:
            np.ndarray | None : The energies of individual atoms within the system.
            
        Description:
            When turn on `ENERGY DEPOSITION`, PWmat will output energy per atom.
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

    def get_fatoms(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: Forces acting on individual atoms of shape=(num_atoms*3,)
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
            np.ndarray | None: Virial tensor of shape=(9,)
        """
        virial_tensor = []

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

        for tmp_idx in [aim_idx + 1, aim_idx + 2, aim_idx + 3]:
            # tmp_str_lst = ['0.120972E+02', '0.483925E+01', '0.242063E+01']
            tmp_str_lst = self.strs_lst[tmp_idx].split()[-3:]

            virial_tensor.append(float(tmp_str_lst[0]))
            virial_tensor.append(float(tmp_str_lst[1]))
            virial_tensor.append(float(tmp_str_lst[2]))

        return np.array(virial_tensor)


class AtomConfig(MSONable):
    """Object for representing the data in a atom.config or final.config file."""
    def __init__(self, structure: Structure, sort_structure: bool = False):
        """
        Args:
            structure (Structure): Structure object
            sort_structure (bool, optional): Whether to sort the structure. Useful if species
                are not grouped properly together. Defaults to False.
        """
        self.structure: Structure = structure
        if sort_structure:
            self.structure = self.structure.get_sorted_structure()
        elements_counter = dict(sorted(Counter(self.structure.species).items()))
        true_names = [f"{tmp_key}{tmp_value}" for (tmp_key, tmp_value) in elements_counter.items()]
        self.true_names = "".join(true_names)

    def __repr__(self):
        return self.get_str()

    def __str__(self):
        return self.get_str()

    @classmethod
    def from_str(cls, data: str, mag: bool = True):
        """Reads a atom.config from a string

        Args:
            data: string containing atom.config data

        Returns:
            AtomConfig object
        """
        ac_extractor = ACstrExtractor(atom_config_str=data)
        properties: dict[str, float] = {}
        structure = Structure(
            lattice=ac_extractor.get_lattice(),
            species=ac_extractor.get_types(),
            coords=ac_extractor.get_coords().reshape(-1, 3),
            coords_are_cartesian=False,
            properties=properties,
        )
        if mag:
            magmoms = ac_extractor.get_magmoms()
            for idx, tmp_site in enumerate(structure):
                tmp_site.properties.update({"magmom": magmoms[idx]})

        return cls(structure)

    @classmethod
    def from_file(cls, filename: str, mag: bool = True) -> AtomConfig:
        """Reads a AtomConfig from a file

        Args:
            filename (str): File name containing AtomConfig data
            mag (bool, optional): Whether to read magnetic moments. Defaults to True.

        Returns:
            AtomConfig object.
        """
        with zopen(filename, "rt") as file:
            return cls.from_str(data=file.read(), mag=mag)

    @classmethod
    def from_dict(cls, dct: dict) -> AtomConfig:
        """
        Args:
            dct: dict containing atom.config data
        """
        return cls(Structure.from_dict(dct["structure"]))

    def get_str(self) -> str:
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
        """Writes AtomConfig to a file."""
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


class GenKpt(MSONable):
    """GenKpt object for reading and writing gen.kpt. This file just generate line-mode kpoints."""

    def __init__(
        self, reciprocal_lattice: np.array, kpoints: dict[str, np.array], path: list[list[str]], density: float = 0.01
    ):
        """
        Args:
            reciprocal_lattice (np.array): Reciprocal lattice with factor of 2*pi
            kpoints (dict[str, np.array]): Kpoints and their corresponding fractional coordinates.
            kpath (list[list[str]]): All kpaths, with each list representing one kpath.
            density (float): The density of kpoints mesh with factor of 2*pi.
        """
        self._reciprocal_lattice: np.array = reciprocal_lattice
        self._kpath: dict = {}
        self._kpath.update({"kpoints": kpoints})
        self._kpath.update({"path": path})
        self._density = density

    @staticmethod
    def from_structure(structure: Structure, dim: int, density: float = 0.01):
        """
        Args:
            strutcure (Structure)
            dim (int): The dimension of material (2 or 3).
            density (float): Kpoints mesh without factor with 2*pi. Program will 
                automatically convert it with 2*pi.
        """
        kpath_set = KPathSeek(structure)
        if dim == 2:
            kpts_2d: dict[str, np.array] = {}
            for tmp_name, tmp_kpt in kpath_set.kpath["kpoints"].items():
                if (tmp_kpt[2]) == 0:
                    kpts_2d.update({tmp_name: tmp_kpt})

            path_2d: list[list[str]] = []
            for tmp_path in kpath_set.kpath["path"]:
                tmp_path_2d: list[str] = []
                for tmp_hsp in tmp_path:
                    if tmp_hsp in kpts_2d:
                        tmp_path_2d.append(tmp_hsp)
                if len(tmp_path_2d) > 1:
                    path_2d.append(tmp_path_2d)
            kpts: dict[str, np.array] = kpts_2d
            path: list[list[str]] = path_2d
        else:
            kpts = kpath_set.kpath["kpoints"]
            path = kpath_set.kpath["path"]
        rec_lattice: np.array = structure.lattice.reciprocal_lattice.matrix  # with 2*pi
        return GenKpt(rec_lattice, kpts, path, density * 2 * np.pi)

    def get_str(self):
        """Returns a string to be written as a gen.kpt file."""
        def calc_distance(hsp1: str, hsp2: str):
            """_summary_
            Returns:
                distance (float). With factor of 2*pi
            """
            hsp1_coord: np.array = np.dot(
                np.array(self._kpath["kpoints"][hsp1]).reshape(1, 3), self._reciprocal_lattice
            )
            hsp2_coord: np.array = np.dot(
                np.array(self._kpath["kpoints"][hsp2]).reshape(1, 3), self._reciprocal_lattice
            )
            return float(np.linalg.norm(hsp2_coord - hsp1_coord))

        discontinue_pairs: list[list[int]] = []
        for ii in range(len(self._kpath["path"]) - 1):
            discontinue_pairs.append([self._kpath["path"][ii][-1], self._kpath["path"][ii + 1][0]])

        flatten_paths: list[str] = [tmp_hsp for tmp_path in self._kpath["path"] for tmp_hsp in tmp_path]
        genkpt_str: str = f"Generated by pymatgen. density={self._density/(2*np.pi)}, "
        genkpt_str += "fractional coordinates in reciprocal lattice.\n"
        for ii in range(len(flatten_paths) - 1):
            if [flatten_paths[ii], flatten_paths[ii + 1]] not in discontinue_pairs:
                genkpt_str += f"{np.ceil(calc_distance(flatten_paths[ii], flatten_paths[ii+1])/self._density)}\n"
                genkpt_str += f"    {self._kpath['kpoints'][flatten_paths[ii]][0]:>12.6f}\t"
                genkpt_str += f"{self._kpath['kpoints'][flatten_paths[ii]][1]:>12.6f}\t"
                genkpt_str += f"{self._kpath['kpoints'][flatten_paths[ii]][2]:>12.6f}\t"
                genkpt_str += f"{flatten_paths[ii]}\n"
                genkpt_str += f"    {self._kpath['kpoints'][flatten_paths[ii+1]][0]:>12.6f}\t"
                genkpt_str += f"{self._kpath['kpoints'][flatten_paths[ii+1]][1]:>12.6f}\t"
                genkpt_str += f"{self._kpath['kpoints'][flatten_paths[ii+1]][2]:>12.6f}\t"
                genkpt_str += f"{flatten_paths[ii+1]}\n"
        return genkpt_str

    def write_file(self, filename: str):
        """Writes gen.kpt to a file."""
        with zopen(filename, "wt") as f:
            f.write(self.get_str())


class HighSymmetryPoint(MSONable):
    """HighSymmetryPoint object for reading and writing HIGH_SYMMETRY_POINTS file which generate line-mode kpoints."""

    def __init__(self, reciprocal_lattice: np.array, kpts: dict[str, list], path: list[list[str]], density: float):
        """
        Args:
            reciprocal_lattice (np.array): Reciprocal lattice.
            kpts (dict[str, list[float]]): Kpoints and their corresponding fractional coordinates.
            path (list[list[str]]): All k-paths, with each list representing one k-path.
            density (float): Density of kpoints mesh with factor of 2*pi.
        """
        self._reciprocal_lattice: np.array = reciprocal_lattice
        self._kpath: dict = {}
        self._kpath.update({"kpoints": kpts})
        self._kpath.update({"path": path})
        self._density = density

    @staticmethod
    def from_structure(structure: Structure, dim: int, density: float = 0.01):
        """_summary_

        Args:
            structure (Structure)
            dim (int): Dimension of material, 2 or 3.
            density (float, optional): Density of kpoints mesh without factor of 2*pi.. Defaults to 0.01.
        """
        reciprocal_lattice: np.array = structure.lattice.reciprocal_lattice.matrix
        gen_kpt = GenKpt.from_structure(structure=structure, dim=dim, density=density)
        return HighSymmetryPoint(
            reciprocal_lattice, gen_kpt._kpath["kpoints"], gen_kpt._kpath["path"], density * 2 * np.pi
        )

    def get_str(self):
        """Returns a string representation of the HIGH_SYMMETRY_POINTS."""
        def calc_distance(hsp1: str, hsp2: str):
            """
            Returns:
                distance (float). With factor of 2*pi
            """
            hsp1_coord: np.array = np.dot(
                np.array(self._kpath["kpoints"][hsp1]).reshape(1, 3), self._reciprocal_lattice
            )
            hsp2_coord: np.array = np.dot(
                np.array(self._kpath["kpoints"][hsp2]).reshape(1, 3), self._reciprocal_lattice
            )
            return float(np.linalg.norm(hsp2_coord - hsp1_coord))

        def get_hsp_row_str(label: str, index: int, coordinate: float):
            if label == "GAMMA":
                return f"G            {index:>4d}         {coordinate:>.6f}\n"
            return f"{label}            {index:>4d}         {coordinate:>.6f}\n"

        discontinue_pairs: list[list[str]] = []
        for ii in range(len(self._kpath["path"]) - 1):
            discontinue_pairs.append([self._kpath["path"][ii][-1], self._kpath["path"][ii + 1][0]])
        flatten_paths: list[str] = [tmp_hsp for tmp_path in self._kpath["path"] for tmp_hsp in tmp_path]
        # flatten_paths = [hsp.replace("GAMMA", "G") for hsp in flatten_paths]

        index: int = 1
        coordinate: float = 0.0
        hsp_str: str = "Label       Index       Coordinate\n"
        hsp_str += get_hsp_row_str(flatten_paths[0], index, coordinate)
        for ii in range(1, len(flatten_paths)):
            if [flatten_paths[ii - 1], flatten_paths[ii]] not in discontinue_pairs:
                coordinate += calc_distance(flatten_paths[ii - 1], flatten_paths[ii]) / (2 * np.pi)
                index += int(np.ceil(calc_distance(flatten_paths[ii - 1], flatten_paths[ii]) / self._density + 1))
            else:
                coordinate += 0.0
                index += 1
            hsp_str += get_hsp_row_str(flatten_paths[ii], index, coordinate)
        return hsp_str

    def write_file(self, filename: str):
        """Write HighSymmetryPoint to a file."""
        with zopen(filename, "wt") as f:
            f.write(self.get_str())
