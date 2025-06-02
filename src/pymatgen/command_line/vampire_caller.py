"""This module implements an interface to the VAMPIRE code for atomistic
simulations of magnetic materials.

This module depends on a compiled vampire executable available in the path.
Please download at https://vampire.york.ac.uk/download/ and
follow the instructions to compile the executable.

If you use this module, please cite:

"Atomistic spin model simulations of magnetic nanomaterials."
R. F. L. Evans, W. J. Fan, P. Chureemart, T. A. Ostler, M. O. A. Ellis
and R. W. Chantrell. J. Phys.: Condens. Matter 26, 103202 (2014)
"""

from __future__ import annotations

import logging
import subprocess
from shutil import which

import pandas as pd
from monty.dev import requires
from monty.json import MSONable

from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper

__author__ = "ncfrey"
__version__ = "0.1"
__maintainer__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"
__status__ = "Development"
__date__ = "June 2019"

logger = logging.getLogger(__name__)

VAMP_EXE = which("vampire-serial")


class VampireCaller:
    """Run Vampire on a material with magnetic ordering and exchange parameter
    information to compute the critical temperature with classical Monte Carlo.

    Attributes:
            sgraph (StructureGraph): Ground state graph.
            unique_site_ids (dict): Maps each site to its unique identifier
            nn_interactions (dict): {i: j} pairs of NN interactions
                between unique sites.
            ex_params (dict): Exchange parameter values (meV/atom)
            mft_t (float): Mean field theory estimate of critical T
            mat_name (str): Formula unit label for input files
            mat_id_dict (dict): Maps sites to material id # for vampire
                indexing.
    """

    @requires(
        VAMP_EXE,  # type: ignore[arg-type]
        "VampireCaller requires vampire-serial to be in the path."
        "Please follow the instructions at https://vampire.york.ac.uk/download/.",
    )
    def __init__(
        self,
        ordered_structures=None,
        energies=None,
        mc_box_size=4.0,
        equil_timesteps=2000,
        mc_timesteps=4000,
        save_inputs=False,
        hm=None,
        avg=True,
        user_input_settings=None,
    ):
        """user_input_settings is a dictionary that can contain:
        * start_t (int): Start MC sim at this temp, defaults to 0 K.
        * end_t (int): End MC sim at this temp, defaults to 1500 K.
        * temp_increment (int): Temp step size, defaults to 25 K.

        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Energies of each relaxed magnetic structure.
            mc_box_size (float): x=y=z dimensions (nm) of MC simulation box
            equil_timesteps (int): number of MC steps for equilibrating
            mc_timesteps (int): number of MC steps for averaging
            save_inputs (bool): if True, save scratch dir of vampire input files
            hm (HeisenbergModel): object already fit to low energy
                magnetic orderings.
            avg (bool): If True, simply use <J> exchange parameter estimate.
                If False, attempt to use NN, NNN, etc. interactions.
            user_input_settings (dict): optional commands for VAMPIRE Monte Carlo

        Todo:
            * Create input files in a temp folder that gets cleaned up after run terminates
        """
        self.mc_box_size = mc_box_size
        self.equil_timesteps = equil_timesteps
        self.mc_timesteps = mc_timesteps
        self.save_inputs = save_inputs
        self.avg = avg

        if not user_input_settings:  # set to empty dict
            self.user_input_settings = {}
        else:
            self.user_input_settings = user_input_settings

        # Get exchange parameters and set instance variables
        if not hm:
            hmapper = HeisenbergMapper(ordered_structures, energies, cutoff=3.0, tol=0.02)

            hm = hmapper.get_heisenberg_model()

        # Attributes from HeisenbergModel
        self.hm = hm
        self.structure = hm.structures[0]  # ground state
        self.sgraph = hm.sgraphs[0]  # ground state graph
        self.unique_site_ids = hm.unique_site_ids
        self.nn_interactions = hm.nn_interactions
        self.dists = hm.dists
        self.tol = hm.tol
        self.ex_params = hm.ex_params
        self.javg = hm.javg

        # Full structure name before reducing to only magnetic ions
        self.mat_name = hm.formula

        # Switch to scratch dir which automatically cleans up vampire inputs files unless user specifies to save them
        # with ScratchDir(
        #     "/scratch", copy_from_current_on_enter=self.save_inputs, copy_to_current_on_exit=self.save_inputs
        # ):

        # Create input files
        self._create_mat()
        self._create_input()
        self._create_ucf()

        # Call Vampire
        with subprocess.Popen([VAMP_EXE], stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
            _stdout, stderr = process.communicate()
            stdout: str = _stdout.decode()

        if stderr:
            van_helsing = stderr.decode()
            if len(van_helsing) > 27:  # Suppress blank warning msg
                logger.warning(van_helsing)

        if process.returncode != 0:
            raise RuntimeError(f"Vampire exited with return code {process.returncode}.")

        self._stdout = stdout
        self._stderr = stderr

        # Process output
        n_mats = max(self.mat_id_dict.values())
        parsed_out, critical_temp = VampireCaller.parse_stdout("output", n_mats)
        self.output = VampireOutput(parsed_out, n_mats, critical_temp)

    def _create_mat(self):
        structure = self.structure
        mat_name = self.mat_name
        magmoms = structure.site_properties["magmom"]

        # Maps sites to material id for vampire inputs
        mat_id_dict = {}

        n_mats = 0
        for key in self.unique_site_ids:
            spin_up, spin_down = False, False
            n_mats += 1  # at least 1 mat for each unique site

            # Check which spin sublattices exist for this site id
            for site in key:
                if magmoms[site] > 0:
                    spin_up = True
                if magmoms[site] < 0:
                    spin_down = True

            # Assign material id for each site
            for site in key:
                if spin_up and not spin_down:
                    mat_id_dict[site] = n_mats
                if spin_down and not spin_up:
                    mat_id_dict[site] = n_mats
                if spin_up and spin_down:
                    # Check if spin up or down shows up first
                    m0 = magmoms[key[0]]
                    if magmoms[site] > 0 and m0 > 0:
                        mat_id_dict[site] = n_mats
                    if magmoms[site] < 0 and m0 < 0:
                        mat_id_dict[site] = n_mats
                    if magmoms[site] > 0 > m0:
                        mat_id_dict[site] = n_mats + 1
                    if magmoms[site] < 0 < m0:
                        mat_id_dict[site] = n_mats + 1

            # Increment index if two sublattices
            if spin_up and spin_down:
                n_mats += 1

        mat_file = [f"material:num-materials={n_mats}"]

        for key in self.unique_site_ids:
            i = self.unique_site_ids[key]  # unique site id

            for site in key:
                mat_id = mat_id_dict[site]

                # Only positive magmoms allowed
                m_magnitude = abs(magmoms[site])

                if magmoms[site] > 0:
                    spin = 1
                elif magmoms[site] < 0:
                    spin = -1
                else:
                    spin = 0

                atom = structure[i].species.reduced_formula

                mat_file += [f"material[{mat_id}]:material-element={atom}"]
                mat_file += [
                    f"material[{mat_id}]:damping-constant=1.0",
                    f"material[{mat_id}]:uniaxial-anisotropy-constant=1.0e-24",  # xx - do we need this?
                    f"material[{mat_id}]:atomic-spin-moment={m_magnitude:.2f} !muB",
                    f"material[{mat_id}]:initial-spin-direction=0,0,{spin}",
                ]

        mat_file = "\n".join(mat_file)
        mat_file_name = f"{mat_name}.mat"

        self.mat_id_dict = mat_id_dict

        with open(mat_file_name, mode="w", encoding="utf-8") as file:
            file.write(mat_file)

    def _create_input(self):
        structure = self.structure
        mc_box_size = self.mc_box_size
        equil_timesteps = self.equil_timesteps
        mc_timesteps = self.mc_timesteps
        mat_name = self.mat_name

        input_script = [f"material:unit-cell-file={mat_name}.ucf"]
        input_script += [f"material:file={mat_name}.mat"]

        # Specify periodic boundary conditions
        input_script += [
            "create:periodic-boundaries-x",
            "create:periodic-boundaries-y",
            "create:periodic-boundaries-z",
        ]

        # Unit cell size in Angstrom
        abc = structure.lattice.abc
        ucx, ucy, ucz = abc[0], abc[1], abc[2]

        input_script += [f"dimensions:unit-cell-size-x = {ucx:.10f} !A"]
        input_script += [f"dimensions:unit-cell-size-y = {ucy:.10f} !A"]
        input_script += [f"dimensions:unit-cell-size-z = {ucz:.10f} !A"]

        # System size in nm
        input_script += [
            f"dimensions:system-size-x = {mc_box_size:.1f} !nm",
            f"dimensions:system-size-y = {mc_box_size:.1f} !nm",
            f"dimensions:system-size-z = {mc_box_size:.1f} !nm",
        ]

        # Critical temperature Monte Carlo calculation
        input_script += [
            "sim:integrator = monte-carlo",
            "sim:program = curie-temperature",
        ]

        # Default Monte Carlo params
        input_script += [
            f"sim:equilibration-time-steps = {equil_timesteps}",
            f"sim:loop-time-steps = {mc_timesteps}",
            "sim:time-steps-increment = 1",
        ]

        # Set temperature range and step size of simulation
        start_t = self.user_input_settings.get("start_t", 0)

        end_t = self.user_input_settings.get("end_t", 1500)

        temp_increment = self.user_input_settings.get("temp_increment", 25)

        input_script += [
            f"sim:minimum-temperature = {start_t}",
            f"sim:maximum-temperature = {end_t}",
            f"sim:temperature-increment = {temp_increment}",
        ]

        # Output to save
        input_script += [
            "output:temperature",
            "output:mean-magnetisation-length",
            "output:material-mean-magnetisation-length",
            "output:mean-susceptibility",
        ]

        input_script = "\n".join(input_script)

        with open("input", mode="w", encoding="utf-8") as file:
            file.write(input_script)

    def _create_ucf(self):
        structure = self.structure
        mat_name = self.mat_name

        abc = structure.lattice.abc
        ucx, ucy, ucz = abc[0], abc[1], abc[2]

        ucf = ["# Unit cell size:"]
        ucf += [f"{ucx:.10f} {ucy:.10f} {ucz:.10f}"]

        ucf += ["# Unit cell lattice vectors:"]
        a1 = list(structure.lattice.matrix[0])
        ucf += [f"{a1[0]:.10f} {a1[1]:.10f} {a1[2]:.10f}"]
        a2 = list(structure.lattice.matrix[1])
        ucf += [f"{a2[0]:.10f} {a2[1]:.10f} {a2[2]:.10f}"]
        a3 = list(structure.lattice.matrix[2])
        ucf += [f"{a3[0]:.10f} {a3[1]:.10f} {a3[2]:.10f}"]

        nmats = max(self.mat_id_dict.values())

        ucf += ["# Atoms num_materials; id cx cy cz mat cat hcat"]
        ucf += [f"{len(structure)} {nmats}"]

        # Fractional coordinates of atoms
        for site, r in enumerate(structure.frac_coords):
            # Back to 0 indexing for some reason...
            mat_id = self.mat_id_dict[site] - 1
            ucf += [f"{site} {r[0]:.10f} {r[1]:.10f} {r[2]:.10f} {mat_id} 0 0"]

        # J_ij exchange interaction matrix
        sgraph = self.sgraph
        n_inter = 0
        for idx in range(len(sgraph.graph.nodes)):
            n_inter += sgraph.get_coordination_of_site(idx)

        ucf += ["# Interactions"]
        ucf += [f"{n_inter} isotropic"]

        iid = 0  # counts number of interaction
        for idx in range(len(sgraph.graph.nodes)):
            connections = sgraph.get_connected_sites(idx)
            for c in connections:
                jimage = c[1]  # relative integer coordinates of atom j
                dx = jimage[0]
                dy = jimage[1]
                dz = jimage[2]
                j = c[2]  # index of neighbor
                dist = round(c[-1], 2)

                # Look up J_ij between the sites
                # if case: Just use <J> estimate
                j_exc = self.hm.javg if self.avg is True else self.hm._get_j_exc(idx, j, dist)

                # Convert J_ij from meV to Joules
                j_exc *= 1.6021766e-22

                j_exc = str(j_exc)  # otherwise this rounds to 0

                ucf += [f"{iid} {idx} {j} {dx} {dy} {dz} {j_exc}"]
                iid += 1

        ucf = "\n".join(ucf)
        ucf_file_name = f"{mat_name}.ucf"

        with open(ucf_file_name, mode="w", encoding="utf-8") as file:
            file.write(ucf)

    @staticmethod
    def parse_stdout(vamp_stdout, n_mats: int) -> tuple:
        """Parse stdout from Vampire.

        Args:
            vamp_stdout (txt file): Vampire 'output' file.
            n_mats (int): Number of materials in Vampire simulation.

        Returns:
            parsed_out (DataFrame): MSONable vampire output.
            critical_temp (float): Calculated critical temp.
        """
        names = [
            "T",
            "m_total",
            *[f"m_{idx + 1}" for idx in range(n_mats)],
            "X_x",
            "X_y",
            "X_z",
            "X_m",
            "nan",
        ]

        # Parsing vampire MC output
        df_stdout = pd.read_csv(vamp_stdout, sep="\t", skiprows=9, header=None, names=names).drop("nan", axis=1)

        parsed_out = df_stdout.to_json()

        # Max of susceptibility <-> critical temp
        critical_temp = df_stdout.iloc[df_stdout.X_m.idxmax()]["T"]

        return parsed_out, critical_temp


class VampireOutput(MSONable):
    """This class processes results from a Vampire Monte Carlo simulation
    and parses the critical temperature.
    """

    def __init__(self, parsed_out=None, nmats=None, critical_temp=None):
        """
        Args:
            parsed_out (str): JSON rep of parsed stdout DataFrame.
            nmats (int): Number of distinct materials (1 for each specie and up/down spin).
            critical_temp (float): Monte Carlo Tc result.
        """
        self.parsed_out = parsed_out
        self.nmats = nmats
        self.critical_temp = critical_temp
