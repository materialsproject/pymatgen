"""
OpenBabel interface module, which opens up access to the hundreds of file
formats supported by OpenBabel. Requires openbabel with python bindings to be
installed. Please consult the openbabel docs https://openbabel.org.
"""

from __future__ import annotations

import copy
import warnings
from typing import TYPE_CHECKING

from monty.dev import requires

from pymatgen.core.structure import IMolecule, Molecule

try:
    from openbabel import openbabel, pybel
except Exception:
    openbabel = pybel = None

if TYPE_CHECKING:
    from typing_extensions import Self

    from pymatgen.analysis.graphs import MoleculeGraph


__author__ = "Shyue Ping Ong, Qi Wang"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 28, 2012"


needs_openbabel = requires(
    openbabel,
    "BabelMolAdaptor requires openbabel to be installed with Python bindings. "
    "Please get it at https://openbabel.org (version >=3.0.0).",
)


class BabelMolAdaptor:
    """
    Adaptor serves as a bridge between OpenBabel's Molecule and pymatgen's
    Molecule.
    """

    @needs_openbabel
    def __init__(self, mol: Molecule | openbabel.OBMol | pybel.Molecule) -> None:
        """Initialize with pymatgen Molecule or OpenBabel's OBMol.

        Args:
            mol: pymatgen's Molecule/IMolecule or OpenBabel OBMol
        """
        if isinstance(mol, IMolecule):
            if not mol.is_ordered:
                raise ValueError("OpenBabel Molecule only supports ordered molecules.")

            # For some reason, manually adding atoms does not seem to create
            # the correct OBMol representation to do things like force field
            # optimization. So we go through the indirect route of creating
            # an XYZ file and reading in that file.
            ob_mol = openbabel.OBMol()
            ob_mol.BeginModify()
            for site in mol:
                coords = list(site.coords)
                atom_no = site.specie.Z
                ob_atom = openbabel.OBAtom()
                ob_atom.thisown = 0
                ob_atom.SetAtomicNum(atom_no)
                ob_atom.SetVector(*map(float, coords))
                ob_mol.AddAtom(ob_atom)
                del ob_atom
            ob_mol.ConnectTheDots()
            ob_mol.PerceiveBondOrders()
            ob_mol.SetTotalSpinMultiplicity(mol.spin_multiplicity)
            ob_mol.SetTotalCharge(int(mol.charge))
            ob_mol.Center()
            ob_mol.EndModify()
            self._ob_mol = ob_mol
        elif isinstance(mol, openbabel.OBMol):
            self._ob_mol = mol
        elif isinstance(mol, pybel.Molecule):
            self._ob_mol = mol.OBMol
        else:
            raise TypeError(f"Unsupported input type {type(mol)}, must be Molecule, openbabel.OBMol or pybel.Molecule")

    @property
    def pymatgen_mol(self) -> Molecule:
        """Pymatgen Molecule object."""
        sp = []
        coords = []
        for atom in openbabel.OBMolAtomIter(self._ob_mol):
            sp.append(atom.GetAtomicNum())
            coords.append([atom.GetX(), atom.GetY(), atom.GetZ()])
        return Molecule(sp, coords)

    @property
    def openbabel_mol(self):
        """OpenBabel's OBMol."""
        return self._ob_mol

    def localopt(self, forcefield: str = "mmff94", steps: int = 500) -> None:
        """
        A wrapper to pybel's localopt method to optimize a Molecule.

        Args:
            forcefield: Default is mmff94. Options are 'gaff', 'ghemical',
                'mmff94', 'mmff94s', and 'uff'.
            steps: Default is 500.
        """
        pybelmol = pybel.Molecule(self._ob_mol)
        pybelmol.localopt(forcefield=forcefield, steps=steps)
        self._ob_mol = pybelmol.OBMol

    def make3d(self, forcefield: str = "mmff94", steps: int = 50) -> None:
        """
        A wrapper to pybel's make3D method generate a 3D structure from a
        2D or 0D structure.
        The 3D structure is made very quickly using a combination of rules
        (e.g. sp3 atoms should have four bonds arranged in a tetrahedron) and
        ring templates (e.g. cyclohexane is shaped like a chair). Once 3D
        coordinates are generated, hydrogens are added and a quick local
        optimization is carried out as default.

        The generated 3D structure can have clashes or have high energy
        structures due to some strain. Please consider to use the conformer
        search or geometry optimization to further optimize the structure.

        Args:
            forcefield: Default is mmff94. Options are 'gaff', 'ghemical',
                'mmff94', 'mmff94s', and 'uff'.
            steps: Default is 50.
        """
        pybelmol = pybel.Molecule(self._ob_mol)
        pybelmol.make3D(forcefield=forcefield, steps=steps)
        self._ob_mol = pybelmol.OBMol

    def add_hydrogen(self) -> None:
        """Add hydrogens (make all hydrogen explicit)."""
        self._ob_mol.AddHydrogens()

    def remove_bond(self, idx1: int, idx2: int) -> None:
        """
        Remove a bond from an openbabel molecule.

        Args:
            idx1: The atom index of one of the atoms participating the in bond
            idx2: The atom index of the other atom participating in the bond
        """
        for obbond in openbabel.OBMolBondIter(self._ob_mol):
            if (obbond.GetBeginAtomIdx() == idx1 and obbond.GetEndAtomIdx() == idx2) or (
                obbond.GetBeginAtomIdx() == idx2 and obbond.GetEndAtomIdx() == idx1
            ):
                self._ob_mol.DeleteBond(obbond)

    def rotor_conformer(self, *rotor_args, algo: str = "WeightedRotorSearch", forcefield: str = "mmff94") -> None:
        """
        Conformer search based on several Rotor Search algorithms of openbabel.
        If the input molecule is not 3D, make3d will be called (generate 3D
        structure, add hydrogen, a quick localopt). All hydrogen atoms need
        to be made explicit.

        Args:
            rotor_args: pass args to Rotor Search in openbabel.
                for "WeightedRotorSearch": (conformers, geomSteps,
                sampleRingBonds-default False)
                for "SystematicRotorSearch": (geomSteps-default 2500,
                sampleRingBonds-default False)
                for "RandomRotorSearch": (conformers, geomSteps-default 2500,
                sampleRingBonds-default False)
            algo (str): Default is "WeightedRotorSearch". Options are
                "SystematicRotorSearch", "RandomRotorSearch", and
                "WeightedRotorSearch".
            forcefield (str): Default is mmff94. Options are 'gaff', 'ghemical',
                'mmff94', 'mmff94s', and 'uff'.
        """
        if self._ob_mol.GetDimension() != 3:
            self.make3d()
        else:
            self.add_hydrogen()

        ff = openbabel.OBForceField.FindType(forcefield)
        if ff == 0:
            warnings.warn(
                f"This input {forcefield=} is not supported "
                "in openbabel. The forcefield will be reset as "
                "default 'mmff94' for now.",
                stacklevel=2,
            )
            ff = openbabel.OBForceField.FindType("mmff94")

        try:
            rotor_search = getattr(ff, algo)
        except AttributeError:
            warnings.warn(
                f"This input conformer search algorithm {algo} is not "
                "supported in openbabel. Options are "
                "'SystematicRotorSearch', 'RandomRotorSearch' "
                "and 'WeightedRotorSearch'. "
                "The algorithm will be reset as default "
                "'WeightedRotorSearch' for now.",
                stacklevel=2,
            )
            rotor_search = ff.WeightedRotorSearch
        rotor_search(*rotor_args)
        ff.GetConformers(self._ob_mol)

    def gen3d_conformer(self) -> None:
        """
        A combined method to first generate 3D structures from 0D or 2D
        structures and then find the minimum energy conformer:

        1. Use OBBuilder to create a 3D structure using rules and ring templates
        2. Do 250 steps of a steepest descent geometry optimization with the
           MMFF94 forcefield
        3. Do 200 iterations of a Weighted Rotor conformational search
           (optimizing each conformer with 25 steps of a steepest descent)
        4. Do 250 steps of a conjugate gradient geometry optimization.

        Warning from openbabel docs:
        For many applications where 100s if not 1000s of molecules need to be
        processed, gen3d is rather SLOW. Sometimes this function can cause a
        segmentation fault.
        A future version of Open Babel will provide options for slow/medium/fast
        3D structure generation which will involve different compromises
        between speed and finding the global energy minimum.
        """
        gen3d = openbabel.OBOp.FindType("Gen3D")
        gen3d.Do(self._ob_mol)

    def confab_conformers(
        self,
        forcefield: str = "mmff94",
        freeze_atoms: list[int] | None = None,
        rmsd_cutoff: float = 0.5,
        energy_cutoff: float = 50.0,
        conf_cutoff: int = 100000,
        verbose: bool = False,
    ) -> list[Molecule]:
        """
        Conformer generation based on Confab to generate all diverse low-energy
        conformers for molecules. This is different from rotor_conformer or
        gen3d_conformer as it aims to not simply to find a low energy
        conformation but to generate several different conformations.

        Args:
            forcefield (str): Default is mmff94. Options are 'gaff', 'ghemical',
                'mmff94', 'mmff94s', and 'uff'.
            freeze_atoms ([int]): index of atoms to be freezed when performing
                conformer search, default is None.
            rmsd_cutoff (float): rmsd_cufoff, default is 0.5 Angstrom.
            energy_cutoff (float): energy_cutoff, default is 50.0 kcal/mol.
            conf_cutoff (int): max number of conformers to test,
                default is 1 million.
            verbose (bool): whether to display information on torsions found,
                default is False.

        Returns:
            list[Molecule]: Molecule objects for generated conformers.
        """
        if self._ob_mol.GetDimension() != 3:
            self.make3d()
        else:
            self.add_hydrogen()

        ff = openbabel.OBForceField.FindType(forcefield)
        if ff == 0:
            print(f"Could not find {forcefield=} in openbabel, the forcefield will be reset as default 'mmff94'")
            ff = openbabel.OBForceField.FindType("mmff94")

        if freeze_atoms:
            print(f"{len(freeze_atoms)} atoms will be freezed")
            constraints = openbabel.OBFFConstraints()

            for atom in openbabel.OBMolAtomIter(self._ob_mol):
                atom_id = atom.GetIndex() + 1
                if id in freeze_atoms:
                    constraints.AddAtomConstraint(atom_id)
            ff.SetConstraints(constraints)

        # Confab conformer generation
        ff.DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, verbose)
        ff.GetConformers(self._ob_mol)

        # Number of conformers generated by Confab conformer generation
        conformer_num = self._ob_mol.NumConformers()

        conformers = []
        for i in range(conformer_num):
            self._ob_mol.SetConformer(i)
            conformer = copy.deepcopy(BabelMolAdaptor(self._ob_mol).pymatgen_mol)
            conformers.append(conformer)
        self._ob_mol.SetConformer(0)
        return conformers

    @property
    def pybel_mol(self) -> Molecule:
        """Pybel's Molecule object."""
        return pybel.Molecule(self._ob_mol)

    def write_file(self, filename: str, file_format: str = "xyz") -> None:
        """
        Uses OpenBabel to output all supported formats.

        Args:
            filename: Filename of file to output
            file_format: String specifying any OpenBabel supported formats.
        """
        mol = pybel.Molecule(self._ob_mol)
        mol.write(file_format, filename, overwrite=True)

    @classmethod
    def from_file(
        cls, filename: str, file_format: str = "xyz", return_all_molecules: bool = False
    ) -> Self | list[Self]:
        """
        Uses OpenBabel to read a molecule from a file in all supported formats.

        Args:
            filename: Filename of input file
            file_format: String specifying any OpenBabel supported formats.
            return_all_molecules: If ``True``, will return a list of
                ``BabelMolAdaptor`` instances, one for each molecule found in
                the file. If ``False``, will return only the first molecule.

        Returns:
            BabelMolAdaptor object or list thereof
        """
        mols = pybel.readfile(str(file_format), str(filename))
        if return_all_molecules:
            return [cls(mol.OBMol) for mol in mols]

        return cls(next(mols).OBMol)

    @classmethod
    def from_molecule_graph(cls, mol: MoleculeGraph) -> Self:
        """
        Read a molecule from a pymatgen MoleculeGraph object.

        Args:
            mol: pymatgen MoleculeGraph object.

        Returns:
            BabelMolAdaptor object
        """
        return cls(mol.molecule)

    @classmethod
    @needs_openbabel
    def from_str(cls, string_data: str, file_format: str = "xyz") -> Self:
        """
        Uses OpenBabel to read a molecule from a string in all supported
        formats.

        Args:
            string_data: String containing molecule data.
            file_format: String specifying any OpenBabel supported formats.

        Returns:
            BabelMolAdaptor object
        """
        mols = pybel.readstring(file_format, string_data)
        return cls(mols.OBMol)
