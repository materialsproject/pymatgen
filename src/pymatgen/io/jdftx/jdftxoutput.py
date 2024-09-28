from __future__ import annotations

from pathlib import Path

from jdftxoutfile import JDFTXOutfile


class JDFTxOutput:
    """
    Class for parsing output of JDFTx calculations.
    """

    calc_dir: Path | None = None
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
        instance.outfile = JDFTXOutfile.from_file(filename)
        instance.calc_dir = Path(filename).parent
        return instance
