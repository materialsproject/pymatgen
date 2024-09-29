# from __future__ import annotations

# from pathlib import Path
# from typing import TYPE_CHECKING, Any

# from pymatgen.io.jdftx.jdftxoutfile import JDFTXOutfile
# import numpy as np

# from os.path import exists as ope

# if TYPE_CHECKING:
#     from pymatgen.electronic_structure.bandstructure import BandStructure


# class JDFTxOutput:
#     """
#     Class for parsing output of JDFTx calculations.
#     """

#     calc_dir: Path | None = None
#     outfile_path: Path | None = None
#     outfile: JDFTXOutfile | None = None
#     set_prefix: str | None = None
#     fprefix: str | None = None
#     #####
#     # pymatgen objects
#     #####
#     bandstruc: BandStructure | None = None
#     #####
#     #

#     @classmethod
#     def from_out_file(cls, filename: str) -> JDFTxOutput:
#         """
#         Initializes a JDFTxOutput object from an output file.

#         Args:
#             filename: Filename of the JDFTx output file.

#         Returns:
#             JDFTxOutput object.
#         """
#         instance = cls()
#         instance._set_outfile_path(Path(filename))
#         instance.calc_dir = Path(filename).parent
#         return cls._from_calc_dir(instance)

#     @classmethod
#     def from_calc_dir(cls, calc_dir: str, set_prefix: str | None = None) -> JDFTxOutput:
#         """
#         Initializes a JDFTxOutput object from a calculation directory.

#         Args:
#             calc_dir: Directory containing the JDFTx output file.

#         Returns:
#             JDFTxOutput object.
#         """
#         instance = cls()
#         instance.set_prefix = set_prefix
#         instance.calc_dir = Path(calc_dir)
#         if not instance.calc_dir.is_dir():
#             raise ValueError(f"{calc_dir} is not a directory. To initialize from an out file, use from_out_file.")
#         return cls._from_calc_dir(instance)

#     @classmethod
#     def _from_calc_dir(cls, instance: JDFTxOutput) -> JDFTxOutput:
#         """
#         Initializes a JDFTxOutput object from a calculation directory.

#         Args:
#             calc_dir: Directory containing the JDFTx output file.

#         Returns:
#             JDFTxOutput object.
#         """
#         if instance.outfile_path is None:
#             instance._set_outfile_path(instance._find_outfile_path())
#         instance.outfile = JDFTXOutfile.from_file(instance.outfile_path)
#         return instance

#     def _set_elec_data(self) -> None:
#         """ Set electronic data.

#         Set electronic data from the JDFTx output file.
#         """
#         self.bandstruc = self._get_bandstruc()

#     def _set_bandstruc(self) -> BandStructure | None:
#         """
#         """
#         eigvals_filename = self._find_generic_path("eigenvals")
#         if eigvals_filename is None:
#             return None
#         eigvals = np.fromfile(eigvals_filename)
#         kfolding = outfile[-1].kgrid
#         kpts_filename = self._find_generic_path("kPts")


#     def _set_outfile_path(self, outfile_path: str) -> None:
#         """ Set the output file path.

#         Set the output file path. Uses the name of the file to set the prefix and fprefix attributes, if not
#         explicitly set by the user.

#         Parameters:
#         ----------
#         outfile_path: str
#             Path to the JDFTx output file.
#         """
#         self.outfile_path = Path(outfile_path)
#         if self.set_prefix is None:
#             outfname = self.outfile_path.name
#             if outfname == "out":
#                 self.set_prefix = ""
#                 self.fprefix = ""
#             else:
#                 self.set_prefix = outfname.split(".")[0]
#                 self.fprefix = f"{self.set_prefix}."

#     def _find_outfile_path(self) -> Path:
#         """
#         Finds the output file in a calculation directory.

#         Args:
#             calc_dir: Directory containing the JDFTx output file.

#         Returns:
#             Path to the output file.
#         """
#         # outfile_path recognized for files name 'out' or (prefix).out
#         calc_dir = self.calc_dir
#         if calc_dir is None:
#             raise ValueError("calc_dir not set.")
#         _fs: list[Path] = [f for f in calc_dir.iterdir() if f.is_file() and f.name.endswith("out")]
#         fs = [f for f in _fs if "." in f.name and f.name.split(".")[-1] == "out"]
#         fs += [f for f in _fs if f.name == "out"]
#         if not len(fs):
#             raise FileNotFoundError(f"No out file found in {calc_dir}")
#         if len(fs) > 1:
#             raise ValueError(
#                 f"Multiple out files found in {calc_dir}. Please specify the out file by "
#                 "initializing with the from_out_file method, or by cleaning up the directory."
#             )
#         return fs[0]

#     def _find_generic_path(self, filesuffix: str) -> Path | None:
#         """ Find a file of a given suffix in a calculation directory.

#         Find a file of a given suffix in a calculation directory.
#         If multiple files are found, the dump prefix will be extracted the the outfile object to narrow it down.

#         Parameters:
#         ----------
#         filesuffix: str
#             Suffix of the file to find.

#         Returns:
#             Path to the file.
#         """
#         calc_dir = self.calc_dir
#         if calc_dir is None:
#             raise ValueError("calc_dir not set.")
#         fs = [f for f in calc_dir.iterdir() if f.is_file() and f.name.endswith(filesuffix)]
#         if not len(fs):
#             return None
#             #raise FileNotFoundError(f"No file with suffix {filesuffix} found in {calc_dir}")
#         elif len(fs) == 1:
#             return fs[0]
#         fname_excpected = f"{self.fprefix}{filesuffix}"
#         fs_exact = [f for f in fs if f.name == fname_excpected]
#         if len(fs_exact):
#             return fs_exact[0]
#         else:
#             raise ValueError(
#                 f"Multiple files with suffix {filesuffix} found in {calc_dir}, "
#                 f"None of which match with the expected file name {fname_excpected}."
#                 "Either specify the calculation prefix upon initialization, or remove the extra files."
#             )

#     def __getattr__(self, name: str) -> Any:
#         """Return attribute.

#         Return attribute. Explicitly defined to allow for passed access from certain contained objects.

#         Parameters:
#         ----------
#         name: str
#             Attribute name.

#         Returns:
#         -------
#         Any
#             Attribute value.
#         """
#         if hasattr(self.outfile, name):
#             return getattr(self.outfile, name)
#         raise AttributeError(f"{name} not found in JDFTxOutput.")
