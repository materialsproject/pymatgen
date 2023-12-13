from __future__ import annotations

import copy
from typing import Union

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.io.pwmat.inputs import ACstrExtractor, AtomConfig, LineLocator
from pymatgen.util.typing import PathLike


class Movement(MSONable):
    def __init__(
        self, filename: PathLike, ionic_step_skip: Union[int, None] = None, ionic_step_offset: Union[int, None] = None
    ):
        """
        Description:
            Extract information from MOVEMENT file which records trajectory during MD.

        @Author: Hanyu Liu
        @email:  domainofbuaa@gmail.com
        """
        self.filename = filename
        self.ionic_step_skip = ionic_step_skip
        self.ionic_step_offset = ionic_step_offset
        self.split_mark = "--------------------------------------"

        self.chunksizes, self.chunkstarts = self._get_chunkinfo()

        self.nionic_steps: int = len(self.chunksizes)
        self.ionic_steps: list[dict] = self._parse_sefv()

    def _get_chunkinfo(self):
        """
        Description:
            Split MOVEMENT into many chunks, so that process it chunk by chunk.
            Chunk contains the row of "--------------------------------------" .

        Returns:
            Tuple[List[int], List[int]]
                chunksizes: List[int]
                chunkstarts: List[int]
        """
        chunksizes: list[int] = []
        row_idxs: list[int] = LineLocator.locate_all_lines(self.filename, self.split_mark)
        chunksizes.append(row_idxs[0])
        for ii in range(1, len(row_idxs)):
            chunksizes.append(row_idxs[ii] - row_idxs[ii - 1])
        chunksizes_bak: list[int] = copy.deepcopy(chunksizes)
        chunksizes_bak.insert(0, 0)
        chunkstarts: list[int] = np.cumsum(chunksizes_bak).tolist()
        chunkstarts.pop(-1)
        return chunksizes, chunkstarts

    @property
    def atom_configs(self):
        """
        Returns:
            atom_configs: List[Structure]
        """
        return [step["atom_config"] for _, step in enumerate(self.ionic_steps)]

    @property
    def etots(self):
        """
        Returns:
            etots: np.ndarray
        """
        return np.array([step["etot"] for _, step in enumerate(self.ionic_steps)])

    @property
    def fatoms(self):
        """
        Returns:
            fatoms: np.ndarray
        """
        return [step["fatoms"] for _, step in enumerate(self.ionic_steps)]

    @property
    def eatoms(self):
        """
        Returns:
            eatoms: np.ndarray
        """
        try:
            return [step["eatoms"] for _, step in enumerate(self.ionic_steps)]
        except KeyError as e:
            print(e)

    @property
    def virials(self):
        try:
            return [step["virial"] for _, step in enumerate(self.ionic_steps)]
        except KeyError as e:
            print(e)

    def _parse_sefv(self):
        ionic_steps: list[dict] = []
        with zopen(self.filename, "rt") as mvt:
            tmp_step: dict = {}
            for ii in range(self.nionic_steps):
                tmp_chunk: str = ""
                for jj in range(self.chunksizes[ii]):
                    tmp_chunk += mvt.readline()
                tmp_step.update({"atom_config": AtomConfig.from_str(tmp_chunk)})
                tmp_step.update({"etot": ACstrExtractor(tmp_chunk).get_etot()[0]})
                tmp_step.update({"fatoms": ACstrExtractor(tmp_chunk).get_fatoms().reshape(-1, 3)})
                try:
                    tmp_step.update({"eatoms": ACstrExtractor(tmp_chunk).get_eatoms()})
                except AttributeError():
                    print(f"Ionic step #{ii} : Energy deposition is turn down.")
                try:
                    tmp_step.update({"virial": ACstrExtractor(tmp_chunk).get_virial().reshape(3, 3)})
                except:
                    print(f"Ionic step #{ii} : No virial infomation.")
                ionic_steps.append(tmp_step)
        return ionic_steps
