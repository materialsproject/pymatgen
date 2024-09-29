from __future__ import annotations

from pathlib import Path

from pymatgen.io.jdftx.jdftxoutfile import JDFTXOutfile


class JDFTxOutput:
    """
    Class for parsing output of JDFTx calculations.
    """

    calc_dir: Path | None = None
    outfile_path: Path | None = None
    outfile: JDFTXOutfile | None = None

    @classmethod
    def from_out_file(cls, filename: str) -> JDFTxOutput:
        """
        Initializes a JDFTxOutput object from an output file.

        Args:
            filename: Filename of the JDFTx output file.

        Returns:
            JDFTxOutput object.
        """
        instance = cls()
        instance.outfile_path = Path(filename)
        instance.calc_dir = Path(filename).parent
        return cls._from_calc_dir(instance)

    @classmethod
    def from_calc_dir(cls, calc_dir: str) -> JDFTxOutput:
        """
        Initializes a JDFTxOutput object from a calculation directory.

        Args:
            calc_dir: Directory containing the JDFTx output file.

        Returns:
            JDFTxOutput object.
        """
        instance = cls()
        instance.calc_dir = Path(calc_dir)
        if not instance.calc_dir.is_dir():
            raise ValueError(f"{calc_dir} is not a directory. To initialize from an out file, use from_out_file.")
        return cls._from_calc_dir(instance)

    @classmethod
    def _from_calc_dir(cls, instance: JDFTxOutput) -> JDFTxOutput:
        """
        Initializes a JDFTxOutput object from a calculation directory.

        Args:
            calc_dir: Directory containing the JDFTx output file.

        Returns:
            JDFTxOutput object.
        """
        if instance.outfile_path is None:
            instance.outfile_path = instance._find_outfile_path()
        instance.outfile = JDFTXOutfile.from_file(instance.outfile_path)
        return instance

    def _find_outfile_path(self) -> Path:
        """
        Finds the output file in a calculation directory.

        Args:
            calc_dir: Directory containing the JDFTx output file.

        Returns:
            Path to the output file.
        """
        # outfile_path recognized for files name 'out' or (prefix).out
        calc_dir = self.calc_dir
        if calc_dir is None:
            raise ValueError("calc_dir not set.")
        _fs: list[Path] = [f for f in calc_dir.iterdir() if f.is_file() and f.name.endswith("out")]
        fs = [f for f in _fs if "." in f.name and f.name.split(".")[-1] == "out"]
        fs += [f for f in _fs if f.name == "out"]
        if not len(fs):
            raise FileNotFoundError(f"No out file found in {calc_dir}")
        if len(fs) > 1:
            raise ValueError(
                f"Multiple out files found in {calc_dir}. Please specify the out file by "
                "initializing with the from_out_file method, or by cleaning up the directory."
            )
        return fs[0]
