# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines utility classes and functions.
"""

from __future__ import annotations

import os
import tempfile
from shutil import which
from subprocess import PIPE, Popen

import numpy as np
from monty.dev import deprecated
from monty.tempfile import ScratchDir

from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Molecule
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.io.packmol import PackmolBoxGen
from pymatgen.util.coord import get_angle

try:
    from openbabel import pybel as pb
except ImportError:
    pb = None

__author__ = "Kiran Mathew, Brandon Wood, Michael Humbert"
__email__ = "kmathew@lbl.gov"


class Polymer:
    """
    Generate polymer chain via Random walk. At each position there are
    a total of 5 possible moves(excluding the previous direction).
    """

    def __init__(
        self,
        start_monomer,
        s_head,
        s_tail,
        monomer,
        head,
        tail,
        end_monomer,
        e_head,
        e_tail,
        n_units,
        link_distance=1.0,
        linear_chain=False,
    ):
        """
        Args:
            start_monomer (Molecule): Starting molecule
            s_head (int): starting atom index of the start_monomer molecule
            s_tail (int): tail atom index of the start_monomer
            monomer (Molecule): The monomer
            head (int): index of the atom in the monomer that forms the head
            tail (int): tail atom index. monomers will be connected from
                tail to head
            end_monomer (Molecule): Terminal molecule
            e_head (int): starting atom index of the end_monomer molecule
            e_tail (int): tail atom index of the end_monomer
            n_units (int): number of monomer units excluding the start and
                terminal molecules
            link_distance (float): distance between consecutive monomers
            linear_chain (bool): linear or random walk polymer chain
        """
        self.start = s_head
        self.end = s_tail
        self.monomer = monomer
        self.n_units = n_units
        self.link_distance = link_distance
        self.linear_chain = linear_chain
        # translate monomers so that head atom is at the origin
        start_monomer.translate_sites(range(len(start_monomer)), -monomer.cart_coords[s_head])
        monomer.translate_sites(range(len(monomer)), -monomer.cart_coords[head])
        end_monomer.translate_sites(range(len(end_monomer)), -monomer.cart_coords[e_head])
        self.mon_vector = monomer.cart_coords[tail] - monomer.cart_coords[head]
        self.moves = {
            1: [1, 0, 0],
            2: [0, 1, 0],
            3: [0, 0, 1],
            4: [-1, 0, 0],
            5: [0, -1, 0],
            6: [0, 0, -1],
        }
        self.prev_move = 1
        # places the start monomer at the beginning of the chain
        self.molecule = start_monomer.copy()
        self.length = 1
        # create the chain
        self._create(self.monomer, self.mon_vector)
        # terminate the chain with the end_monomer
        self.n_units += 1
        end_mon_vector = end_monomer.cart_coords[e_tail] - end_monomer.cart_coords[e_head]
        self._create(end_monomer, end_mon_vector)
        self.molecule = Molecule.from_sites(self.molecule.sites)

    def _create(self, monomer, mon_vector):
        """
        create the polymer from the monomer

        Args:
            monomer (Molecule)
            mon_vector (numpy.array): molecule vector that starts from the
                start atom index to the end atom index
        """
        while self.length != (self.n_units - 1):
            if self.linear_chain:
                move_direction = np.array(mon_vector) / np.linalg.norm(mon_vector)
            else:
                move_direction = self._next_move_direction()
            self._add_monomer(monomer.copy(), mon_vector, move_direction)

    def _next_move_direction(self):
        """
        pick a move at random from the list of moves
        """
        nmoves = len(self.moves)
        move = np.random.randint(1, nmoves + 1)
        while self.prev_move == (move + 3) % nmoves:
            move = np.random.randint(1, nmoves + 1)
        self.prev_move = move
        return np.array(self.moves[move])

    def _align_monomer(self, monomer, mon_vector, move_direction):
        """
        rotate the monomer so that it is aligned along the move direction

        Args:
            monomer (Molecule)
            mon_vector (numpy.array): molecule vector that starts from the
                start atom index to the end atom index
            move_direction (numpy.array): the direction of the polymer chain
                extension
        """
        axis = np.cross(mon_vector, move_direction)
        origin = monomer[self.start].coords
        angle = get_angle(mon_vector, move_direction)
        op = SymmOp.from_origin_axis_angle(origin, axis, angle)
        monomer.apply_operation(op)

    def _add_monomer(self, monomer, mon_vector, move_direction):
        """
        extend the polymer molecule by adding a monomer along mon_vector direction

        Args:
            monomer (Molecule): monomer molecule
            mon_vector (numpy.array): monomer vector that points from head to tail.
            move_direction (numpy.array): direction along which the monomer
                will be positioned
        """
        translate_by = self.molecule.cart_coords[self.end] + self.link_distance * move_direction
        monomer.translate_sites(range(len(monomer)), translate_by)
        if not self.linear_chain:
            self._align_monomer(monomer, mon_vector, move_direction)
        # add monomer if there are no crossings
        does_cross = False
        for i, site in enumerate(monomer):
            try:
                self.molecule.append(site.specie, site.coords, properties=site.properties)
            except Exception:
                does_cross = True
                polymer_length = len(self.molecule)
                self.molecule.remove_sites(range(polymer_length - i, polymer_length))
                break
        if not does_cross:
            self.length += 1
            self.end += len(self.monomer)


@deprecated(PackmolBoxGen, "PackmolRunner is being phased out in favor of the packmol I/O class.")
class PackmolRunner:
    """
    Wrapper for the Packmol software that can be used to pack various types of
    molecules into a one single unit.
    """

    def __init__(
        self,
        mols,
        param_list,
        input_file="pack.inp",
        tolerance=2.0,
        filetype="xyz",
        control_params=None,
        auto_box=True,
        output_file="packed.xyz",
        bin="packmol",
    ):
        """
        Args:
            mols:
                list of Molecules to pack
            input_file:
                name of the packmol input file
            tolerance:
                min distance between the atoms
            filetype:
                input/output structure file type
            control_params:
                packmol control parameters dictionary. Basically
                all parameters other than structure/atoms
            param_list:
                list of parameters containing dicts for each molecule
            auto_box:
                put the molecule assembly in a box
            output_file:
                output file name. The extension will be adjusted
                according to the filetype
        """
        self.packmol_bin = bin.split()
        if not which(self.packmol_bin[-1]):
            raise RuntimeError(
                "PackmolRunner requires the executable 'packmol' to be in "
                "the path. Please download packmol from "
                "https://github.com/leandromartinez98/packmol "
                "and follow the instructions in the README to compile. "
                "Don't forget to add the packmol binary to your path"
            )
        self.mols = mols
        self.param_list = param_list
        self.input_file = input_file
        self.boxit = auto_box
        self.control_params = control_params or {"maxit": 20, "nloop": 600}
        if not self.control_params.get("tolerance"):
            self.control_params["tolerance"] = tolerance
        if not self.control_params.get("filetype"):
            self.control_params["filetype"] = filetype
        if not self.control_params.get("output"):
            self.control_params["output"] = f"{output_file.split('.')[0]}.{self.control_params['filetype']}"
        if self.boxit:
            self._set_box()

    @staticmethod
    def _format_param_val(param_val):
        """
        Internal method to format values in the packmol parameter dictionaries

        Args:
            param_val:
                Some object to turn into String

        Returns:
            String representation of the object
        """
        if isinstance(param_val, list):
            return " ".join(str(x) for x in param_val)
        return str(param_val)

    def _set_box(self):
        """
        Set the box size for the molecular assembly
        """
        net_volume = 0.0
        for idx, mol in enumerate(self.mols):
            length = max(np.max(mol.cart_coords[:, i]) - np.min(mol.cart_coords[:, i]) for i in range(3)) + 2.0
            net_volume += (length**3.0) * float(self.param_list[idx]["number"])
        length = net_volume ** (1.0 / 3.0)
        for idx, _mol in enumerate(self.mols):
            self.param_list[idx]["inside box"] = f"0.0 0.0 0.0 {length} {length} {length}"

    def _write_input(self, input_dir="."):
        """
        Write the packmol input file to the input directory.

        Args:
            input_dir (str): path to the input directory
        """
        with open(os.path.join(input_dir, self.input_file), "w", encoding="utf-8") as inp:
            for k, v in self.control_params.items():
                inp.write(f"{k} {self._format_param_val(v)}\n")
            # write the structures of the constituent molecules to file and set
            # the molecule id and the corresponding filename in the packmol
            # input file.
            for idx, mol in enumerate(self.mols):
                filename = os.path.join(input_dir, f"{idx}.{self.control_params['filetype']}")
                # pdb
                if self.control_params["filetype"] == "pdb":
                    self.write_pdb(mol, filename, num=idx + 1)
                # all other filetypes
                else:
                    a = BabelMolAdaptor(mol)
                    pm = pb.Molecule(a.openbabel_mol)
                    pm.write(
                        self.control_params["filetype"],
                        filename=filename,
                        overwrite=True,
                    )

                inp.write("\n")
                inp.write(f"structure {os.path.join(input_dir, str(idx))}.{self.control_params['filetype']}\n")
                for k, v in self.param_list[idx].items():
                    inp.write(f"  {k} {self._format_param_val(v)}\n")
                inp.write("end structure\n")

    def run(self, site_property=None):
        """
        Write the input file to the scratch directory, run packmol and return
        the packed molecule to the current working directory.

        Args:
            site_property (str): if set then the specified site property
                for the final packed molecule will be restored.

        Returns:
                Molecule object
        """
        with tempfile.TemporaryDirectory() as scratch_dir:
            self._write_input(input_dir=scratch_dir)
            with open(os.path.join(scratch_dir, self.input_file)) as packmol_input:
                with Popen(self.packmol_bin, stdin=packmol_input, stdout=PIPE, stderr=PIPE) as p:
                    (stdout, stderr) = p.communicate()
            output_file = os.path.join(self.control_params["output"])
            if os.path.isfile(output_file):
                packed_mol = BabelMolAdaptor.from_file(output_file, self.control_params["filetype"])
                packed_mol = packed_mol.pymatgen_mol
                print(f"packed molecule written to {self.control_params['output']}")
                if site_property:
                    packed_mol = self.restore_site_properties(site_property=site_property, filename=output_file)
                return packed_mol
            raise RuntimeError(f"Packmol execution failed. {stdout}\n{stderr}")

    @staticmethod
    def write_pdb(mol, filename, name=None, num=None):
        """
        dump the molecule into pdb file with custom residue name and number.
        """
        # ugly hack to get around the openbabel issues with inconsistent
        # residue labelling.
        scratch = tempfile.gettempdir()
        with ScratchDir(scratch, copy_to_current_on_exit=False) as _:
            mol.to(fmt="pdb", filename="tmp.pdb")
            bma = BabelMolAdaptor.from_file("tmp.pdb", "pdb")

        num = num or 1
        name = name or f"ml{num}"

        # bma = BabelMolAdaptor(mol)
        pbm = pb.Molecule(bma._obmol)
        for x in pbm.residues:
            x.OBResidue.SetName(name)
            x.OBResidue.SetNum(num)

        pbm.write(format="pdb", filename=filename, overwrite=True)

    def _set_residue_map(self):
        """
        map each residue to the corresponding molecule.
        """
        self.map_residue_to_mol = {}
        lookup = {}
        for idx, mol in enumerate(self.mols):
            if mol.formula not in lookup:
                mol.translate_sites(indices=range(len(mol)), vector=-mol.center_of_mass)
                lookup[mol.formula] = mol.copy()
            self.map_residue_to_mol[f"ml{idx + 1}"] = lookup[mol.formula]

    def convert_obatoms_to_molecule(self, atoms, residue_name=None, site_property="ff_map"):
        """
        Convert list of openbabel atoms to MOlecule.

        Args:
            atoms ([OBAtom]): list of OBAtom objects
            residue_name (str): the key in self.map_residue_to_mol. Usec to
                restore the site properties in the final packed molecule.
            site_property (str): the site property to be restored.

        Returns:
            Molecule object
        """
        restore_site_props = residue_name is not None

        if restore_site_props and not hasattr(self, "map_residue_to_mol"):
            self._set_residue_map()

        coords = []
        zs = []
        for atm in atoms:
            coords.append(list(atm.coords))
            zs.append(atm.atomicnum)

        mol = Molecule(zs, coords)

        if restore_site_props:

            props = []

            ref = self.map_residue_to_mol[residue_name].copy()

            # sanity check
            assert len(mol) == len(ref)
            assert ref.formula == mol.formula

            # the packed molecules have the atoms in the same order..sigh!
            for i, site in enumerate(mol):
                assert site.specie.symbol == ref[i].specie.symbol
                props.append(getattr(ref[i], site_property))

            mol.add_site_property(site_property, props)

        return mol

    def restore_site_properties(self, site_property="ff_map", filename=None):
        """
        Restore the site properties for the final packed molecule.

        Args:
            site_property (str):
            filename (str): path to the final packed molecule.

        Returns:
            Molecule
        """
        # only for pdb
        if not self.control_params["filetype"] == "pdb":
            raise ValueError()

        filename = filename or self.control_params["output"]
        bma = BabelMolAdaptor.from_file(filename, "pdb")
        pbm = pb.Molecule(bma._obmol)

        assert len(pbm.residues) == sum(x["number"] for x in self.param_list)

        packed_mol = self.convert_obatoms_to_molecule(
            pbm.residues[0].atoms,
            residue_name=pbm.residues[0].name,
            site_property=site_property,
        )

        for resid in pbm.residues[1:]:
            mol = self.convert_obatoms_to_molecule(resid.atoms, residue_name=resid.name, site_property=site_property)
            for site in mol:
                packed_mol.append(site.species, site.coords, properties=site.properties)

        return packed_mol


class LammpsRunner:
    """
    LAMMPS wrapper
    """

    def __init__(self, input_filename="lammps.in", bin="lammps"):
        """
        Args:
            input_filename (str): input file name
            bin (str): command to run, excluding the input file name
        """
        self.lammps_bin = bin.split()
        if not which(self.lammps_bin[-1]):
            raise RuntimeError(
                f"LammpsRunner requires the executable {self.lammps_bin[-1]} to be in the path. "
                "Please download and install LAMMPS from "
                "https://www.lammps.org/. "
                "Don't forget to add the binary to your path"
            )
        self.input_filename = input_filename

    def run(self):
        """
        Write the input/data files and run LAMMPS.
        """
        lammps_cmd = self.lammps_bin + ["-in", self.input_filename]
        print(f"Running: {' '.join(lammps_cmd)}")
        with Popen(lammps_cmd, stdout=PIPE, stderr=PIPE) as p:
            (stdout, stderr) = p.communicate()
        return stdout, stderr


if __name__ == "__main__":
    ethanol_coords = [
        [0.00720, -0.56870, 0.00000],
        [-1.28540, 0.24990, 0.00000],
        [1.13040, 0.31470, 0.00000],
        [0.03920, -1.19720, 0.89000],
        [0.03920, -1.19720, -0.89000],
        [-1.31750, 0.87840, 0.89000],
        [-1.31750, 0.87840, -0.89000],
        [-2.14220, -0.42390, -0.00000],
        [1.98570, -0.13650, -0.00000],
    ]
    ethanol = Molecule(["C", "C", "O", "H", "H", "H", "H", "H", "H"], ethanol_coords)
    water_coords = [
        [9.626, 6.787, 12.673],
        [9.626, 8.420, 12.673],
        [10.203, 7.604, 12.673],
    ]
    water = Molecule(["H", "H", "O"], water_coords)
    pmr = PackmolRunner(
        [ethanol, water],
        [
            {"number": 1, "fixed": [0, 0, 0, 0, 0, 0], "centerofmass": ""},
            {"number": 15, "inside sphere": [0, 0, 0, 5]},
        ],
        input_file="packmol_input.inp",
        tolerance=2.0,
        filetype="xyz",
        control_params={"nloop": 1000},
        auto_box=False,
        output_file="cocktail.xyz",
    )
    s = pmr.run()
