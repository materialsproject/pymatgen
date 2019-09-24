# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import warnings
import copy
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from monty.dev import requires

try:
    import openbabel as ob
    import pybel as pb
except Exception:
    pb = None
    ob = None

"""
OpenBabel interface module, which opens up access to the hundreds of file
formats supported by OpenBabel. Requires openbabel with python bindings to be
installed. Please consult the
`openbabel documentation <http://openbabel.org/wiki/Main_Page>`_.
"""

__author__ = "Shyue Ping Ong, Qi Wang"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 28, 2012"


class BabelMolAdaptor:
    """
    Adaptor serves as a bridge between OpenBabel's Molecule and pymatgen's
    Molecule.
    """

    @requires(pb and ob,
              "BabelMolAdaptor requires openbabel to be installed with "
              "Python bindings. Please get it at http://openbabel.org.")
    def __init__(self, mol):
        """
        Initializes with pymatgen Molecule or OpenBabel"s OBMol.

        Args:
            mol: pymatgen's Molecule or OpenBabel OBMol
        """
        if isinstance(mol, Molecule):
            if not mol.is_ordered:
                raise ValueError("OpenBabel Molecule only supports ordered "
                                 "molecules.")

            # For some reason, manually adding atoms does not seem to create
            # the correct OBMol representation to do things like force field
            # optimization. So we go through the indirect route of creating
            # an XYZ file and reading in that file.
            obmol = ob.OBMol()
            obmol.BeginModify()
            for site in mol:
                coords = [c for c in site.coords]
                atomno = site.specie.Z
                obatom = ob.OBAtom()
                obatom.thisown = 0
                obatom.SetAtomicNum(atomno)
                obatom.SetVector(*coords)
                obmol.AddAtom(obatom)
                del obatom
            obmol.ConnectTheDots()
            obmol.PerceiveBondOrders()
            obmol.SetTotalSpinMultiplicity(mol.spin_multiplicity)
            obmol.SetTotalCharge(int(mol.charge))
            obmol.Center()
            obmol.Kekulize()
            obmol.EndModify()
            self._obmol = obmol
        elif isinstance(mol, ob.OBMol):
            self._obmol = mol

    @property
    def pymatgen_mol(self):
        """
        Returns pymatgen Molecule object.
        """
        sp = []
        coords = []
        for atom in ob.OBMolAtomIter(self._obmol):
            sp.append(atom.GetAtomicNum())
            coords.append([atom.GetX(), atom.GetY(), atom.GetZ()])
        return Molecule(sp, coords)

    @property
    def openbabel_mol(self):
        """
        Returns OpenBabel's OBMol.
        """
        return self._obmol

    def localopt(self, forcefield='mmff94', steps=500):
        """
        A wrapper to pybel's localopt method to optimize a Molecule.

        Args:
            forcefield: Default is mmff94. Options are 'gaff', 'ghemical',
                'mmff94', 'mmff94s', and 'uff'.
            steps: Default is 500.
        """
        pbmol = pb.Molecule(self._obmol)
        pbmol.localopt(forcefield=forcefield, steps=steps)
        self._obmol = pbmol.OBMol

    def make3d(self, forcefield="mmff94", steps=50):
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
        pbmol = pb.Molecule(self._obmol)
        pbmol.make3D(forcefield=forcefield, steps=steps)
        self._obmol = pbmol.OBMol

    def add_hydrogen(self):
        """
        Add hydrogens (make all hydrogen explicit).
        """
        self._obmol.AddHydrogens()

    def remove_bond(self, idx1, idx2):
        """
        Remove a bond from an openbabel molecule

        Args:
            idx1: The atom index of one of the atoms participating the in bond
            idx2: The atom index of the other atom participating in the bond
        """
        for obbond in ob.OBMolBondIter(self._obmol):
            if (obbond.GetBeginAtomIdx() == idx1 and obbond.GetEndAtomIdx() == idx2) or (
                    obbond.GetBeginAtomIdx() == idx2 and obbond.GetEndAtomIdx() == idx1):
                self._obmol.DeleteBond(obbond)

    def rotor_conformer(self, *rotor_args, algo="WeightedRotorSearch",
                        forcefield="mmff94"):
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
        if self._obmol.GetDimension() != 3:
            self.make3d()
        else:
            self.add_hydrogen()

        ff = ob.OBForceField_FindType(forcefield)
        if ff == 0:
            warnings.warn("This input forcefield {} is not supported "
                          "in openbabel. The forcefield will be reset as "
                          "default 'mmff94' for now.".format(forcefield))
            ff = ob.OBForceField_FindType("mmff94")

        try:
            rotor_search = getattr(ff, algo)
        except AttributeError:
            warnings.warn("This input conformer search algorithm {} is not "
                          "supported in openbabel. Options are "
                          "'SystematicRotorSearch', 'RandomRotorSearch' "
                          "and 'WeightedRotorSearch'. "
                          "The algorithm will be reset as default "
                          "'WeightedRotorSearch' for now.".format(algo))
            rotor_search = ff.WeightedRotorSearch
        rotor_search(*rotor_args)
        ff.GetConformers(self._obmol)

    def gen3d_conformer(self):
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
        gen3d = ob.OBOp.FindType("Gen3D")
        gen3d.Do(self._obmol)

    def confab_conformers(self, forcefield="mmff94", freeze_atoms=None,
                          rmsd_cutoff=0.5, energy_cutoff=50.0,
                          conf_cutoff=100000, verbose=False):
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
            conf_cutoff (float): max number of conformers to test,
                default is 1 million.
            verbose (bool): whether to display information on torsions found,
                default is False.

        Returns:
             (list): list of pymatgen Molecule objects for generated conformers.
        """
        if self._obmol.GetDimension() != 3:
            self.make3d()
        else:
            self.add_hydrogen()

        ff = ob.OBForceField_FindType(forcefield)
        if ff == 0:
            print("Could not find forcefield {} in openbabel, the forcefield "
                  "will be reset as default 'mmff94'".format(forcefield))
            ff = ob.OBForceField_FindType("mmff94")

        if freeze_atoms:
            print('{} atoms will be freezed'.format(len(freeze_atoms)))
            constraints = ob.OBFFConstraints()

            for atom in ob.OBMolAtomIter(self._obmol):
                atom_id = atom.GetIndex() + 1
                if id in freeze_atoms:
                    constraints.AddAtomConstraint(atom_id)
            ff.SetConstraints(constraints)

        # Confab conformer generation
        ff.DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff,
                          verbose)
        ff.GetConformers(self._obmol)

        # Number of conformers generated by Confab conformer generation
        conformer_num = self._obmol.NumConformers()

        conformers = []
        for i in range(conformer_num):
            self._obmol.SetConformer(i)
            conformer = copy.deepcopy(BabelMolAdaptor(self._obmol).pymatgen_mol)
            conformers.append(conformer)
        self._obmol.SetConformer(0)
        return conformers

    @property
    def pybel_mol(self):
        """
        Returns Pybel's Molecule object.
        """
        return pb.Molecule(self._obmol)

    def write_file(self, filename, file_format="xyz"):
        """
        Uses OpenBabel to output all supported formats.

        Args:
            filename: Filename of file to output
            file_format: String specifying any OpenBabel supported formats.
        """
        mol = pb.Molecule(self._obmol)
        return mol.write(file_format, filename, overwrite=True)

    @staticmethod
    def from_file(filename, file_format="xyz", return_all_molecules=False):
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
        mols = pb.readfile(str(file_format), str(filename))
        if return_all_molecules:
            return [BabelMolAdaptor(mol.OBMol) for mol in mols]
        else:
            return BabelMolAdaptor(next(mols).OBMol)

    @staticmethod
    def from_molecule_graph(mol):
        """
        Read a molecule from a pymatgen MoleculeGraph object.

        Args:
            mol: pymatgen MoleculeGraph object.

        Returns:
            BabelMolAdaptor object
        """
        if isinstance(mol, MoleculeGraph):
            return BabelMolAdaptor(mol.molecule)

    @staticmethod
    def from_string(string_data, file_format="xyz"):
        """
        Uses OpenBabel to read a molecule from a string in all supported
        formats.

        Args:
            string_data: String containing molecule data.
            file_format: String specifying any OpenBabel supported formats.

        Returns:
            BabelMolAdaptor object
        """
        mols = pb.readstring(str(file_format), str(string_data))
        return BabelMolAdaptor(mols.OBMol)
