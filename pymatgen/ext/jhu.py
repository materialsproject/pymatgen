# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module creates an interface to the JHU kpoints servlet,
see http://muellergroup.jhu.edu/K-Points.html.
"""

import os
import shutil
import tempfile

import requests

from pymatgen.io.vasp.inputs import Kpoints

__author__ = "Joseph Montoya"
__copyright__ = "Copyright 2017, The Materials Project"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__date__ = "June 22, 2017"


def get_kpoints(
    structure,
    min_distance=0,
    min_total_kpoints=1,
    kppra=None,
    gap_distance=7,
    remove_symmetry=None,
    include_gamma="auto",
    header="simple",
    incar=None,
):
    """
    Get kpoints object from JHU servlet, per Wisesa-McGill-Mueller
    methodology.  Refer to http://muellergroup.jhu.edu/K-Points.html
    and P. Wisesa, K. A. McGill, T. Mueller, Phys. Rev. B 93,
    155109 (2016)

    Args:
        structure (Structure): structure object
        min_distance (float): The minimum allowed distance
            between lattice points on the real-space superlattice
        min_total_kpoints (int): The minimum allowed number of
            total k-points in the Brillouin zone.
        kppra (float): minimum k-points per reciprocal atom.
        gap_distance (float): auto-detection threshold for
            non-periodicity (in slabs, nanowires, etc.)
        remove_symmetry (string): optional flag to control
            symmetry options, can be none, structural,
            time_reversal, or all
        include_gamma (string or bool): whether to include
            gamma point
        header (string): "verbose" or "simple", denotes
            the verbosity of the header
        incar (Incar): incar object to upload
    """
    config = locals()
    config.pop("structure", "incar")

    # Generate PRECALC string
    precalc = "".join(["{}={}\n".format(k, v) for k, v in config.items()])
    precalc = precalc.replace("_", "").upper()
    precalc = precalc.replace("REMOVESYMMETRY", "REMOVE_SYMMETRY")
    precalc = precalc.replace("TIMEREVERSAL", "TIME_REVERSAL")
    url = "http://muellergroup.jhu.edu:8080/PreCalcServer/PreCalcServlet"
    temp_dir_name = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(temp_dir_name)

    precalc_file = open("PRECALC", "w+")
    poscar_file = open("POSCAR", "w+")
    incar_file = open("INCAR", "w+")

    precalc_file.write(precalc)
    poscar_file.write(structure.to("POSCAR"))
    files = [("fileupload", precalc_file), ("fileupload", poscar_file)]
    if incar:
        incar_file.write(incar.get_string())
        files.append(("fileupload", incar_file))

    precalc_file.seek(0)
    poscar_file.seek(0)
    incar_file.seek(0)

    r = requests.post(url, files=files)

    precalc_file.close()
    poscar_file.close()
    incar_file.close()
    kpoints = Kpoints.from_string(r.text)
    os.chdir(cwd)

    shutil.rmtree(temp_dir_name)

    return kpoints
