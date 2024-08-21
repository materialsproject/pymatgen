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
    from typing_extensions import Self

    from pymatgen.util.typing import PathLike

__author__ = "Hanyu Liu"
__email__ = "domainofbuaa@gmail.com"
__date__ = "2024-1-16"


class LineLocator(MSONable):
    """Find the line indices (starts from 1) of a certain paragraph of text from the file."""

    @staticmethod
    def locate_all_lines(file_path: PathLike, content: str, exclusion: str = "") -> list[int]:
        """Locate the line in file where a certain paragraph of text is located (return all indices).

        Args:
            file_path (PathLike): Absolute path to file.
            content (str): Certain paragraph of text that needs to be located.
            exclusion (str): Certain paragraph of text that is excluded.
        """
        row_idxs: list[int] = []  # starts from 1 to be compatible with linecache package
        row_no: int = 0
        with zopen(file_path, mode="rt") as file:
            for row_content in file:
                row_no += 1
                if content.upper() in row_content.upper() and (
                    not exclusion or exclusion.upper() not in row_content.upper()
                ):
                    row_idxs.append(row_no)
        return row_idxs


class ListLocator(MSONable):
    """Find the element indices (starts from 0) of a certain paragraph of text from the list."""

    @staticmethod
    def locate_all_lines(strs_lst: list[str], content: str, exclusion: str = "") -> list[int]:
        """Locate the elements in list where a certain paragraph of text is located (return all indices).

        Args:
            strs_lst (list[str]): List of strings.
            content (str): Certain paragraph of text that needs to be located.
            exclusion (str): Certain paragraph of text that is excluded.
        """
        str_idxs: list[int] = []  # starts from 0 to be compatible with list
        str_no: int = -1
        for tmp_str in strs_lst:
            str_no += 1
            if (content.upper() in tmp_str.upper()) and (not exclusion or exclusion.upper() not in tmp_str.upper()):
                str_idxs.append(str_no)
        return str_idxs


class ACExtractorBase(ABC):
    """A parent class of ACExtractor and ACstrExtractor, ensuring that they are as consistent as possible."""

    @abstractmethod
    def get_n_atoms(self) -> int:
        """Get the number of atoms in structure defined by atom.config file."""

    @abstractmethod
    def get_lattice(self) -> np.ndarray:
        """Get the lattice of structure defined by atom.config file."""

    @abstractmethod
    def get_types(self) -> np.ndarray:
        """Get atomic number of atoms in structure defined by atom.config file."""

    @abstractmethod
    def get_coords(self) -> np.ndarray:
        """Get fractional coordinates of atoms in structure defined by atom.config file."""

    @abstractmethod
    def get_magmoms(self) -> np.ndarray:
        """Get atomic magmoms of atoms in structure defined by atom.config file."""


class ACExtractor(ACExtractorBase):
    """Extract information contained in atom.config : number of atoms, lattice, types, frac_coords, magmoms."""

    def __init__(self, file_path: PathLike) -> None:
        """Initialization function.

        Args:
            file_path (str): The absolute path of atom.config file.
        """
        self.atom_config_path = file_path
        self.n_atoms = self.get_n_atoms()
        self.lattice = self.get_lattice()
        self.types = self.get_types()
        self.coords = self.get_coords()
        self.magmoms = self.get_magmoms()

    def get_n_atoms(self) -> int:
        """Return the number of atoms in the structure."""
        first_row = linecache.getline(str(self.atom_config_path), 1)
        return int(first_row.split()[0])

    def get_lattice(self) -> np.ndarray:
        """Return the lattice of structure.

        Returns:
            lattice: np.ndarray, shape = (9,)
        """
        basis_vectors: list[float] = []
        content: str = "LATTICE"
        idx_row: int = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[0]
        for row_idx in [idx_row + 1, idx_row + 2, idx_row + 3]:
            row_content: list[str] = linecache.getline(str(self.atom_config_path), row_idx).split()[:3]
            for value in row_content:
                basis_vectors.append(float(value))

        return np.array(basis_vectors)

    def get_types(self) -> np.ndarray:
        """Return the atomic number of atoms in structure.

        Returns:
            np.ndarray: Atomic numbers in order corresponding to sites
        """
        content = "POSITION"
        idx_row = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[0]
        with open(self.atom_config_path) as file:
            atom_config_content = file.readlines()
        atomic_numbers_content = atom_config_content[idx_row : idx_row + self.n_atoms]
        atomic_numbers_lst = [int(row.split()[0]) for row in atomic_numbers_content]  # convert str to int
        return np.array(atomic_numbers_lst)

    def get_coords(self) -> np.ndarray:
        """Return the fractional coordinates in structure.

        Returns:
            np.ndarray: Fractional coordinates.
        """
        coords_lst: list[np.ndarray] = []
        content: str = "POSITION"
        idx_row: int = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[0]
        with open(self.atom_config_path) as file:
            atom_config_content = file.readlines()
        """
        row_content:
            '29         0.377262291145329         0.128590184800933         0.257759805813488     1  1  1'
        """
        for row_content in atom_config_content[idx_row : idx_row + self.n_atoms]:
            row_content_lst = row_content.split()
            coord_tmp = [float(value) for value in row_content_lst[1:4]]  # convert str to float.
            coords_lst.append(np.array(coord_tmp))
        return np.array(coords_lst).reshape(-1)

    def get_magmoms(self) -> np.ndarray:
        """Return the magenetic moments of atoms in structure.

        Returns:
            np.ndarray: The magnetic moments of individual atoms.
        """
        content: str = "MAGNETIC"
        magnetic_moments_lst: list[float] = []
        try:  # Error: not containing magmoms info.
            idx_row = LineLocator.locate_all_lines(file_path=self.atom_config_path, content=content)[-1]

            with open(self.atom_config_path) as file:
                atom_config_content = file.readlines()

            magnetic_moments_content = atom_config_content[idx_row : idx_row + self.n_atoms]
            # MAGNETIC
            # 3 0.0 # atomic_number magmom
            # ...
            magnetic_moments_lst = [
                float(tmp_magnetic_moment.split()[-1]) for tmp_magnetic_moment in magnetic_moments_content
            ]
        except Exception:
            magnetic_moments_lst = [0 for _ in range(self.n_atoms)]
        return np.array(magnetic_moments_lst)


class ACstrExtractor(ACExtractorBase):
    """Extract information from atom.config file. You can get str by slicing the MOVEMENT."""

    def __init__(self, atom_config_str: str):
        """Initialization function.

        Args:
            atom_config_str (str): A string describing the structure in atom.config file.
        """
        self.atom_config_str = atom_config_str
        self.strs_lst = self.atom_config_str.split("\n")
        self.num_atoms = self.get_n_atoms()

    def get_n_atoms(self) -> int:
        """Return the number of atoms in structure.

        Returns:
            int: The number of atoms
        """
        return int(self.strs_lst[0].split()[0].strip())

    def get_lattice(self) -> np.ndarray:
        """Return the lattice of structure.

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
        """Return the atomic number of atoms in structure.

        Returns:
            np.ndarray: Types of elements.
        """
        aim_content = "POSITION"
        aim_idx = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)[0]
        strs_lst = self.strs_lst[aim_idx + 1 : aim_idx + self.num_atoms + 1]
        atomic_numbers_lst = [int(entry.split()[0]) for entry in strs_lst]
        return np.array(atomic_numbers_lst)

    def get_coords(self) -> np.ndarray:
        """Return the fractional coordinate of atoms in structure.

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
        """Return the magnetic moments of atoms in structure.

        Returns:
            np.ndarray: Atomic magnetic moments.
        """
        magnetic_moments_lst: list[float] = []
        aim_content: str = "MAGNETIC"
        aim_idxs: list[int] = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content)
        if len(aim_idxs) == 0:
            magnetic_moments_lst = [0.0 for _ in range(self.num_atoms)]
        else:
            aim_idx = aim_idxs[0]
            magnetic_moments_content = self.strs_lst[aim_idx + 1 : aim_idx + self.num_atoms + 1]
            magnetic_moments_lst = [
                float(tmp_magnetic_moment.split()[-1]) for tmp_magnetic_moment in magnetic_moments_content
            ]
        return np.array(magnetic_moments_lst)

    def get_e_tot(self) -> np.ndarray:
        """Return the total energy of structure.

        Returns:
            np.ndarray: The total energy of the material system.
        """
        # strs_lst:
        #   [' 216 atoms', 'Iteration (fs) =    0.3000000000E+01',
        #    ' Etot', 'Ep', 'Ek =   -0.2831881714E+05  -0.2836665392E+05   0.4783678177E+02',
        #    ' SCF =     7']
        strs_lst = self.strs_lst[0].split(",")
        aim_index = ListLocator.locate_all_lines(strs_lst=strs_lst, content="EK")[0]
        # strs_lst[aim_index].split() :
        #   ['Ek', '(eV)', '=', '-0.2831881714E+05', '-0.2836665392E+05', '0.4783678177E+02']
        return np.array([float(strs_lst[aim_index].split("=")[1].split()[0].strip())])

    def get_atom_energies(self) -> np.ndarray | None:
        """Return the energies of individual atoms in material system.

        When turning on `ENERGY DEPOSITION`, PWmat will output energy per atom.

        Returns:
            np.ndarray | None: The energies of individual atoms within the material system.
        """
        energies = []
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
            energies.append(float(tmp_str.split()[1]))
        return np.array(energies)

    def get_atom_forces(self) -> np.ndarray:
        """Return the force on atoms in material system.

        Returns:
            np.ndarray: Forces acting on individual atoms of shape=(num_atoms*3,)
        """
        forces = []
        aim_content = "Force".upper()
        aim_idx = ListLocator.locate_all_lines(strs_lst=self.strs_lst, content=aim_content, exclusion="average")[0]
        for line in self.strs_lst[aim_idx + 1 : aim_idx + self.num_atoms + 1]:
            # ['14', '0.089910342901203', '0.077164252174742', '0.254144099204679']
            forces.append([float(val) for val in line.split()[1:4]])
        return -np.array(forces).reshape(-1)

    def get_virial(self) -> np.ndarray | None:
        """Return the virial tensor of material system.

        Returns:
            np.ndarray | None: Virial tensor of shape=(9,)
        """
        virial_tensor: list[float] = []
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
            virial_tensor += (float(tmp_str_lst[0]), float(tmp_str_lst[1]), float(tmp_str_lst[2]))

        return np.array(virial_tensor)


class AtomConfig(MSONable):
    """Object for representing the data in a atom.config or final.config file."""

    def __init__(self, structure: Structure, sort_structure: bool = False):
        """Initialization function.

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
    def from_str(cls, data: str, mag: bool = False) -> Self:
        """Reads a atom.config from a string.

        Args:
            data (str): string containing atom.config data
            mag (bool, optional): Whether to read magnetic moment information.

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
                tmp_site.properties |= {"magmom": magmoms[idx]}

        return cls(structure)

    @classmethod
    def from_file(cls, filename: PathLike, mag: bool = False) -> Self:
        """Get a AtomConfig from a file.

        Args:
            filename (PathLike): File name containing AtomConfig data
            mag (bool, optional): Whether to read magnetic moments. Defaults to True.

        Returns:
            AtomConfig object.
        """
        with zopen(filename, "rt") as file:
            return cls.from_str(data=file.read(), mag=mag)

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Get a AtomConfig object from a dictionary.

        Args:
            dct: dict containing atom.config data

        Returns:
            AtomConfig object.
        """
        return cls(Structure.from_dict(dct["structure"]))

    def get_str(self) -> str:
        """Return a string describing the structure in atom.config format.

        Returns:
            str: String representation of atom.config
        """
        # This corrects for VASP really annoying bug of crashing on lattices
        # which have triple product < 0. We will just invert the lattice vectors.
        lattice = self.structure.lattice
        if np.linalg.det(lattice.matrix) < 0:
            lattice = Lattice(-lattice.matrix)

        lines: list[str] = []
        lines += (f"\t{self.structure.num_sites} atoms\n", "Lattice vector\n")
        for idx in range(3):
            lines.append(f"{lattice.matrix[idx][0]:>15f}{lattice.matrix[idx][1]:>15f}{lattice.matrix[idx][2]:>15f}\n")
        lines.append("Position, move_x, move_y, move_z\n")
        for site_idx in range(len(self.structure)):
            lines += (
                f"{int(self.structure.species[site_idx].Z):>4d}",
                f"{self.structure.frac_coords[site_idx][0]:>15f}",
                f"{self.structure.frac_coords[site_idx][1]:>15f}",
                f"{self.structure.frac_coords[site_idx][2]:>15f}",
                "   1   1   1\n",
            )
        if "magmom" in self.structure.sites[0].properties:
            lines.append("MAGNETIC\n")
            for _, tmp_site in enumerate(self.structure.sites):
                lines.append(f"{int(tmp_site.specie.Z):>4d}{tmp_site.properties['magmom']:>15f}\n")
        return "".join(lines)

    def write_file(self, filename: PathLike, **kwargs):
        """Write AtomConfig to a file."""
        with zopen(filename, "wt") as file:
            file.write(self.get_str(**kwargs))

    def as_dict(self):
        """
        Returns:
            dict.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "structure": self.structure.as_dict(),
            "true_names": self.true_names,
        }


class GenKpt(MSONable):
    """Read and write gen.kpt. This file just generates line-mode kpoints."""

    def __init__(
        self,
        reciprocal_lattice: np.ndarray,
        kpoints: dict[str, np.ndarray],
        path: list[list[str]],
        density: float = 0.01,
    ) -> None:
        """Initialization function.

        Args:
            reciprocal_lattice (np.array): Reciprocal lattice with factor of 2*pi.
            kpoints (dict[str, np.array]): Kpoints and their corresponding fractional coordinates.
            kpath (list[list[str]]): All kpaths, with each list representing one kpath.
            density (float): The density of kpoints mesh with factor of 2*pi.
        """
        self.reciprocal_lattice: np.ndarray = reciprocal_lattice
        self.kpath: dict = {}
        self.kpath |= {"kpoints": kpoints}
        self.kpath |= {"path": path}
        self.density = density

    @classmethod
    def from_structure(cls, structure: Structure, dim: int, density: float = 0.01) -> Self:
        """Obtain a AtomConfig object from Structure object.

        Args:
            structure (Structure): A structure object.
            dim (int): The dimension of the material system (2 or 3).
            density (float): Kpoints mesh without factor with 2*pi. Program will
                automatically convert it with 2*pi.
        """
        kpath_set = KPathSeek(structure)
        if dim == 2:
            kpts_2d: dict[str, np.ndarray] = {}
            for tmp_name, tmp_kpt in kpath_set.kpath["kpoints"].items():
                if (tmp_kpt[2]) == 0:
                    kpts_2d |= {tmp_name: tmp_kpt}

            path_2d: list[list[str]] = []
            for tmp_path in kpath_set.kpath["path"]:
                tmp_path_2d: list[str] = []
                for tmp_hsp in tmp_path:
                    if tmp_hsp in kpts_2d:
                        tmp_path_2d.append(tmp_hsp)
                if len(tmp_path_2d) > 1:
                    path_2d.append(tmp_path_2d)
            kpts: dict[str, np.ndarray] = kpts_2d
            path: list[list[str]] = path_2d
        else:
            kpts = kpath_set.kpath["kpoints"]
            path = kpath_set.kpath["path"]
        rec_lattice: np.ndarray = structure.lattice.reciprocal_lattice.matrix  # with 2*pi
        return cls(rec_lattice, kpts, path, density * 2 * np.pi)

    def get_str(self):
        """Get a string to be written as a gen.kpt file."""

        def calc_distance(hsp1: str, hsp2: str) -> float:
            """Calculate the distance between two high symmetry points.

            Args:
                hsp1 (str): The name of the first high symmetry point.
                hsp2 (str): The name of the second high symmetry point.

            Returns:
                float: The distance between two high symmetry points. With factor of 2*pi.
            """
            hsp1_coord: np.ndarray = np.dot(
                np.array(self.kpath["kpoints"][hsp1]).reshape(1, 3), self.reciprocal_lattice
            )
            hsp2_coord: np.ndarray = np.dot(
                np.array(self.kpath["kpoints"][hsp2]).reshape(1, 3), self.reciprocal_lattice
            )
            return float(np.linalg.norm(hsp2_coord - hsp1_coord))

        discontinue_pairs: list[list[int]] = []
        for idx in range(len(self.kpath["path"]) - 1):
            discontinue_pairs.append([self.kpath["path"][idx][-1], self.kpath["path"][idx + 1][0]])

        flatten_paths: list[str] = [tmp_hsp for tmp_path in self.kpath["path"] for tmp_hsp in tmp_path]
        gen_kpt_str: str = f"Generated by pymatgen. density={self.density / (2 * np.pi)}, "
        gen_kpt_str += "fractional coordinates in reciprocal lattice.\n"
        for idx in range(len(flatten_paths) - 1):
            if [flatten_paths[idx], flatten_paths[idx + 1]] not in discontinue_pairs:
                gen_kpt_str += f"{np.ceil(calc_distance(flatten_paths[idx], flatten_paths[idx + 1]) / self.density)}\n"
                gen_kpt_str += f"    {self.kpath['kpoints'][flatten_paths[idx]][0]:>12.6f}\t"
                gen_kpt_str += f"{self.kpath['kpoints'][flatten_paths[idx]][1]:>12.6f}\t"
                gen_kpt_str += f"{self.kpath['kpoints'][flatten_paths[idx]][2]:>12.6f}\t"
                gen_kpt_str += f"{flatten_paths[idx]}\n"
                gen_kpt_str += f"    {self.kpath['kpoints'][flatten_paths[idx + 1]][0]:>12.6f}\t"
                gen_kpt_str += f"{self.kpath['kpoints'][flatten_paths[idx + 1]][1]:>12.6f}\t"
                gen_kpt_str += f"{self.kpath['kpoints'][flatten_paths[idx + 1]][2]:>12.6f}\t"
                gen_kpt_str += f"{flatten_paths[idx + 1]}\n"
        return gen_kpt_str

    def write_file(self, filename: PathLike):
        """Write gen.kpt to a file.

        Args:
            filename (PathLike): The absolute path of file to be written.
        """
        with zopen(filename, "wt") as file:
            file.write(self.get_str())


class HighSymmetryPoint(MSONable):
    """Read and write HIGH_SYMMETRY_POINTS file which generate line-mode kpoints."""

    def __init__(self, reciprocal_lattice: np.ndarray, kpts: dict[str, list], path: list[list[str]], density: float):
        """Initialization function.

        Args:
            reciprocal_lattice (np.array): Reciprocal lattice.
            kpts (dict[str, list[float]]): Kpoints and their corresponding fractional coordinates.
            path (list[list[str]]): All k-paths, with each list representing one k-path.
            density (float): Density of kpoints mesh with factor of 2*pi.
        """
        self.reciprocal_lattice: np.ndarray = reciprocal_lattice
        self.kpath: dict = {}
        self.kpath |= {"kpoints": kpts}
        self.kpath |= {"path": path}
        self.density = density

    @classmethod
    def from_structure(cls, structure: Structure, dim: int, density: float = 0.01) -> Self:
        """Obtain HighSymmetry object from Structure object.

        Args:
            structure (Structure): A structure object.
            dim (int): Dimension of the material system (2 or 3).
            density (float, optional): Density of kpoints mesh without factor of 2*pi. Defaults to 0.01.
                The program will automatically convert it to with factor of 2*pi.
        """
        reciprocal_lattice: np.ndarray = structure.lattice.reciprocal_lattice.matrix
        gen_kpt = GenKpt.from_structure(structure=structure, dim=dim, density=density)
        return cls(reciprocal_lattice, gen_kpt.kpath["kpoints"], gen_kpt.kpath["path"], density * 2 * np.pi)

    def get_str(self) -> str:
        """Get a string describing high symmetry points in HIGH_SYMMETRY_POINTS format."""

        def calc_distance(hsp1: str, hsp2: str) -> float:
            """Calculate the distance of two high symmetry points.

            Returns:
                float: The distance between two high symmetry points with factor of 2*pi.
            """
            hsp1_coord: np.ndarray = np.dot(
                np.array(self.kpath["kpoints"][hsp1]).reshape(1, 3), self.reciprocal_lattice
            )
            hsp2_coord: np.ndarray = np.dot(
                np.array(self.kpath["kpoints"][hsp2]).reshape(1, 3), self.reciprocal_lattice
            )
            return float(np.linalg.norm(hsp2_coord - hsp1_coord))

        def get_hsp_row_str(label: str, index: int, coordinate: float) -> str:
            """
            Return string containing name, index, coordinate of the certain high symmetry point
            in HIGH_SYMMETRY_POINTS format.

            Args:
                label (str): Name of the high symmetry point.
                index (int): Index of the high symmetry point.
                coordinate (float): Coordinate in bandstructure of the high symmetry point.

            Returns:
                str: String containing name, index, coordinate of the certain high symmetry point
                    in HIGH_SYMMETRY_POINTS format.
            """
            if label == "GAMMA":
                return f"G            {index:>4d}         {coordinate:>.6f}\n"
            return f"{label}            {index:>4d}         {coordinate:>.6f}\n"

        discontinue_pairs: list[list[str]] = []
        for ii in range(len(self.kpath["path"]) - 1):
            discontinue_pairs.append([self.kpath["path"][ii][-1], self.kpath["path"][ii + 1][0]])
        flatten_paths: list[str] = [tmp_hsp for tmp_path in self.kpath["path"] for tmp_hsp in tmp_path]
        # flatten_paths = [hsp.replace("GAMMA", "G") for hsp in flatten_paths]

        index: int = 1
        coordinate: float = 0.0
        hsp_str: str = "Label       Index       Coordinate\n"
        hsp_str += get_hsp_row_str(flatten_paths[0], index, coordinate)
        for ii in range(1, len(flatten_paths)):
            if [flatten_paths[ii - 1], flatten_paths[ii]] not in discontinue_pairs:
                coordinate += calc_distance(flatten_paths[ii - 1], flatten_paths[ii]) / (2 * np.pi)
                index += int(np.ceil(calc_distance(flatten_paths[ii - 1], flatten_paths[ii]) / self.density + 1))
            else:
                coordinate += 0.0
                index += 1
            hsp_str += get_hsp_row_str(flatten_paths[ii], index, coordinate)
        return hsp_str

    def write_file(self, filename: PathLike):
        """Write HighSymmetryPoint to a file."""
        with zopen(filename, "wt") as file:
            file.write(self.get_str())
