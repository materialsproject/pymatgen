from __future__ import annotations

import copy
import linecache
from typing import TYPE_CHECKING
from io import StringIO

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.io.pwmat.inputs import (
    ACstrExtractor, 
    AtomConfig, 
    LineLocator
)

if TYPE_CHECKING:
    from pymatgen.util.typing import PathLike


class Movement(MSONable):
    def __init__(self, filename: PathLike, ionic_step_skip: int | None = None, ionic_step_offset: int | None = None):
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
        return np.array([step["fatoms"] for _, step in enumerate(self.ionic_steps)])

    @property
    def eatoms(self):
        """
        Returns:
            eatoms: np.ndarray
        """
        return np.array([step["eatoms"] for _, step in enumerate(self.ionic_steps) if ("eatoms" in step)])

    @property
    def virials(self):
        return np.array([step["virial"] for _, step in enumerate(self.ionic_steps) if ("virial" in step)])

    def _parse_sefv(self):
        ionic_steps: list[dict] = []
        with zopen(self.filename, "rt") as mvt:
            tmp_step: dict = {}
            for ii in range(self.nionic_steps):
                tmp_chunk: str = ""
                for _ in range(self.chunksizes[ii]):
                    tmp_chunk += mvt.readline()
                tmp_step.update({"atom_config": AtomConfig.from_str(tmp_chunk)})
                tmp_step.update({"etot": ACstrExtractor(tmp_chunk).get_etot()[0]})
                tmp_step.update({"fatoms": ACstrExtractor(tmp_chunk).get_fatoms().reshape(-1, 3)})
                eatoms: np.ndarray | None = ACstrExtractor(tmp_chunk).get_fatoms()
                if eatoms is not None:
                    tmp_step.update({"eatoms": ACstrExtractor(tmp_chunk).get_eatoms()})
                else:
                    print(f"Ionic step #{ii} : Energy deposition is turn down.")
                virial: np.ndarray | None = ACstrExtractor(tmp_chunk).get_virial()
                if virial is not None:
                    tmp_step.update({"virial": virial.reshape(3, 3)})
                else:
                    print(f"Ionic step #{ii} : No virial information.")
                ionic_steps.append(tmp_step)
        return ionic_steps



class OutFermi(MSONable):
    """
    Description:
        Extract fermi energy (eV) from OUT.FERMI

    @Author: Hanyu Liu
    @Email: domainofbuaa@gmail.com
    """
    def __init__(self, filename:str):
        self.filename:str = filename
        with zopen(self.filename, "rt") as f:
            self._efermi:float = np.round(
                float(f.readline().split()[-2].strip()), 3)
    
    @property
    def efermi(self):
        return self._efermi



class Report(MSONable):
    """
    Description:
        Extract information of spin, kpoints, bands, eigenvalues from REPORT file.
        
    @Author: Hanyu Liu
    @Email: domainofbuaa@gmail.com
    """
    def __init__(self, filename:str):
        self.filename = filename
        self._spin, self._num_kpts, self._num_bands = self._parse_band()
        self._eigenvalues = self._parse_eigen()
        self._kpts, self._kpts_weight, self._hsps = self._parse_kpt()
    
    
    def _parse_band(self):
        """_summary_
        

        Returns:
            spin: int, Whether turn on spin or not.
                1 : turn down the spin
                2 : turn on the spin
            num_kpts: int, The number of kpoints.
            num_bands: int, The number of bands.
        """
        content = "SPIN"
        row_idx:int = LineLocator.locate_all_lines(file_path=self.filename, content=content)[0]
        spin = int( linecache.getline(self.filename, row_idx).split()[-1].strip() )
        
        content:str = "NUM_KPT"
        row_idx:int = LineLocator.locate_all_lines(file_path=self.filename, content=content)[0]
        num_kpts = int( linecache.getline(self.filename, row_idx).split()[-1].strip() )
        
        content:str = "NUM_BAND"
        row_idx:int = LineLocator.locate_all_lines(file_path=self.filename, content=content)[0]
        num_bands = int( linecache.getline(self.filename, row_idx).split()[-1].strip() )
        return spin, num_kpts, num_bands
    
    
    def _parse_eigen(self):
        """
        Return:
            eigenvales: np.array
                shape = (1 or 2, self.num_kpts, self.num_bands)
        """
        num_rows:int = int( np.ceil(self._num_bands / 5) )
        content:str = "eigen energies, in eV"
        rows_lst:list[int] = LineLocator.locate_all_lines(file_path=self.filename, content=content)
        rows_array:np.array = np.array(rows_lst).reshape(self._spin, -1)
        eigenvalues:np.array = np.zeros((self._spin, self._num_kpts, self._num_bands))
        
        for ii in range(self._spin):
            for jj in range(self._num_kpts):
                tmp_eigenvalues_str = ""
                for kk in range(num_rows):
                    tmp_eigenvalues_str += linecache.getline(self.filename, rows_array[ii][jj]+kk+1)
                tmp_eigenvalues_array = np.array( 
                    [float(eigen_value) for eigen_value in tmp_eigenvalues_str.split()] )        
                for kk in range(self._num_bands):
                    eigenvalues[ii][jj][kk] = tmp_eigenvalues_array[kk]
        return eigenvalues

    
    def _parse_kpt(self):
        """
        Returns:
            kpts: np.array, The fractional coordinates of KPoints
            kpts_weight: np.array, The weight of KPoints
            hsps: dict[str, np.array], The name and coordinates of high symmetric points
        """
        num_rows:int = int( self._num_kpts )
        content:str = "total number of K-point:"
        row_idx:int = LineLocator.locate_all_lines(self.filename, content)[0]
        kpts:np.array = np.zeros((self._num_kpts, 3))
        kpts_weight:np.array = np.zeros((self._num_kpts))
        hsps:dict[str, np.array] = {}
        for ii in range(num_rows):
            #  0.00000     0.00000    0.00000     0.03704           G
            tmp_row_lst:list[str] = linecache.getline(self.filename, row_idx+ii+1).split();
            for jj in range(3):
                kpts[ii][jj] = float( tmp_row_lst[jj].strip() )
            kpts_weight[ii] = float( tmp_row_lst[3].strip() )
            
            if len(tmp_row_lst) == 5:
                hsps.update(
                    {   tmp_row_lst[4]: np.array([
                            float(tmp_row_lst[0]), 
                            float(tmp_row_lst[1]), 
                            float(tmp_row_lst[2])]
                        )
                    }
                )
        return kpts, kpts_weight, hsps
    
    
    @property
    def spin(self):
        return self._spin
    
    @property
    def num_kpts(self):
        return self._num_kpts
    
    @property
    def num_bands(self):
        return self._num_bands
    
    @property
    def eigenvalues(self):
        return self._eigenvalues
    
    @property
    def kpts(self):
        return self._kpts
    
    @property
    def kpts_weight(self):
        return self._kpts_weight
    
    @property
    def hsps(self):
        return self._hsps


class Dosspin(MSONable):
    """_summary_
    Description:
        Extract information of DOS from :
            - DOS.totalspin, DOS.totalspin_projected
            - DOS.spinup, DOS.spinup_projected
            - DOS.spindown, DOS.spindown_projected

    @Author: Hanyu Liu
    @Email: domainofbuaa@gmail.com
    """
    def __init__(self, filename:str):
        self.filename = filename
        self._labels, self._dos = self._parse()
    
    
    def _parse(self):
        """
        Returns:
            labels: list[str], The label of DOS, e.g. Total, Cr-3S, ...
            dos: np.array
        """
        labels:list[str] = []
        labels = linecache.getline(self.filename, 1).split()[1:]
        dos_str:str = ""
        with zopen(self.filename) as f:
            f.readline()
            dos_str = f.read()
        dos:np.array = np.loadtxt(StringIO(dos_str))
        return labels, dos
    
    
    @property
    def labels(self):
        return self._labels
    
    
    @property
    def dos(self):
        return self._dos
    
    
    def get_partial_dos(self, part:str):
        """
        Description:
            Get partial dos for give element or orbital.
        
        Args:
            part: str, 
                e.g. 'Energy', 'Total', 'Cr-3S', 'Cr-3P', 
                     'Cr-4S', 'Cr-3D', 'I-4D', 'I-5S', 'I-5P'
                e.g. 'Energy', 'Total', 'Cr-3S', 'Cr-3Pz', 
                     'Cr-3Px', 'Cr-3Py', 'Cr-4S', 'Cr-3Dz2', 
                     'Cr-3Dxz', 'Cr-3Dyz', 'Cr-3D(x^2-y^2)', 
                     'Cr-3Dxy', 'I-4Dz2', 'I-4Dxz', 'I-4Dyz', 
                     'I-4D(x^2-y^2)', 'I-4Dxy', 'I-5S', 'I-5Pz', 
                     'I-5Px', 'I-5Py'
        
        Returns:
            partial_dos: np.array
        """
        part:str = part.upper()
        labels_upper:list[str] = [tmp_label.upper() for tmp_label in self._labels]
        idx_dos = labels_upper.index(part) 
        return self._dos[:, idx_dos]