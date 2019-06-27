# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import subprocess
import logging
import numpy as np
import pandas as pd

from monty.dev import requires
from monty.os.path import which

from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper
from pymatgen.analysis.magnetism.analyzer import CollinearMagneticStructureAnalyzer

"""
This module implements an interface to the VAMPIRE code for atomistic
simulations of magnetic materials.

This module depends on a compiled vampire executable available in the path.
Please download at https://vampire.york.ac.uk/download/ and
follow the instructions to compile the executable.

If you use this module, please cite the following:
    
"Atomistic spin model simulations of magnetic nanomaterials."
R. F. L. Evans, W. J. Fan, P. Chureemart, T. A. Ostler, M. O. A. Ellis 
and R. W. Chantrell. J. Phys.: Condens. Matter 26, 103202 (2014)
"""

__author__ = "ncfrey"
__version__ = "0.1"
__maintainer__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"
__status__ = "Development"
__date__ = "June 2019"

VAMPEXE = which("vampire-serial")

class VampireCaller:
        
    @requires(VAMPEXE, "VampireCaller requires vampire-serial to be in the path."
    "Please follow the instructions at https://vampire.york.ac.uk/download/.")
    
    def __init__(self, ordered_structures, energies, 
        mc_box_size=4.0, equil_timesteps=10000, mc_timesteps=10000, hm=None):
        
        """
        Run Vampire on a material with magnetic ordering and exchange parameter information to compute the critical temperature with classical Monte Carlo.
        
        Args:
            ordered_structures (list): Structure objects with magmoms.
            energies (list): Energies of each relaxed magnetic structure.
            mc_box_size (float): x=y=z dimensions (nm) of MC simulation box
            equil_timesteps (int): number of MC steps for equilibrating
            mc_timesteps (int): number of MC steps for averaging
            hm (HeisenbergMapper): object already fit to low energy
                magnetic orderings.

        Parameters:
            sgraph (StructureGraph): Ground state graph.
            unique_site_ids (dict): Maps each site to its unique identifier
            nn_interacations (dict): {i: j} pairs of NN interactions
                between unique sites.
            ex_params (dict): Exchange parameter values (meV/atom)
            mft_t (float): Mean field theory estimate of critical T
            mat_name (str): Formula unit label for input files

        """
    
        self.mc_box_size = mc_box_size
        self.equil_timesteps = equil_timesteps
        self.mc_timesteps = mc_timesteps

        # Sort by energy if not already sorted
        ordered_structures = [s for _, s in
                              sorted(zip(energies, ordered_structures), reverse=False)]

        energies = sorted(energies, reverse=False)

        # Get exchange parameters and set instance variables
        if not hm:
            hm = HeisenbergMapper(ordered_structures, energies, cutoff=7.5,
                tol=0.02)

        # Instance attributes from HeisenbergMapper
        self.hm = hm
        self.sgraph = hm.sgraphs[0]  # ground state graph
        self.unique_site_ids = hm.unique_site_ids
        self.nn_interactions = hm.nn_interactions
        self.dists = hm.dists
        self.tol = hm.tol
        self.ex_params = hm.get_exchange()

        # Get mean field estimate of critical temp
        j_avg = hm.estimate_exchange()
        self.mft_t = hm.get_mft_temperature(j_avg)

        # Create input files
        structure = ordered_structures[0]  # ground state
        self.mat_name = str(structure.composition.reduced_formula)
        self._create_mat(structure)
        self._create_input(structure)
        self._create_ucf(structure)
        
        # Call Vampire
        process = subprocess.Popen(['vampire-serial'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        stdout = stdout.decode()

        if stderr:
            stderr = stderr.decode()
            logging.warning(stderr)
        
        if process.returncode != 0:
            raise RuntimeError("Vampire exited with return code {}.".format(process.returncode))
        
        self._stdout = stdout
        self._stderr = stderr
    
        # Process output
        self.output = VampireOutput('output')
       
    def _create_mat(self, structure):
        """Todo: 
            * make different materials (sublattices) for up and down spin
            * if up/down spins are on inequivalent sites this is already
              taken care of
        """

        mat_name = self.mat_name
        nmats = len(self.unique_site_ids)
        
        mat_file = ['material:num-materials=%d' % (nmats)]
        
        for key in self.unique_site_ids:
            for site in key:
                i = self.unique_site_ids[key]
                mat_id = i + 1  # Starts at 1
                atom = structure[i].species.reduced_formula

                # Only positive magmoms are allowed
                mag_mom = abs(structure.site_properties['magmom'][site])
                mat_file += ['material[%d]:material-element=%s' % (mat_id, atom)]
                mat_file += ['material[%d]:damping-constant=1.0' % (mat_id),
                                 'material[%d]:uniaxial-anisotropy-constant=0.0' % (mat_id),
                                 'material[%d]:atomic-spin-moment=%.2f !muB' % (mat_id, mag_mom)]
            
        mat_file = "\n".join(mat_file)
        mat_file_name = mat_name + '.mat'
        
        with open(mat_file_name, 'w') as f:
            f.write(mat_file)
        
    def _create_input(self, structure):
        """Todo:
            * Does temperature increment need to evenly divide the
              interval?
        """
        
        mcbs = self.mc_box_size
        equil_timesteps = self.equil_timesteps
        mc_timesteps = self.mc_timesteps
        mat_name = self.mat_name

        input_script = ['material:unit-cell-file=%s.ucf' % (mat_name)]
        input_script += ['material:file=%s.mat' % (mat_name)]
        
        # Specify periodic boundary conditions
        input_script += ['create:periodic-boundaries-x',
                         'create:periodic-boundaries-y',
                         'create:periodic-boundaries-z']
        
        # Unit cell size in Angstrom
        abc = structure.lattice.abc
        ucx, ucy, ucz = abc[0], abc[1], abc[2]
        
        input_script += ['dimensions:unit-cell-size-x = %.10f !A' % (ucx)]
        input_script += ['dimensions:unit-cell-size-y = %.10f !A' % (ucy)]
        input_script += ['dimensions:unit-cell-size-z = %.10f !A' % (ucz)]

        # System size in nm
        input_script += ['dimensions:system-size-x = %.1f !nm' % (mcbs),
                         'dimensions:system-size-y = %.1f !nm' % (mcbs),
                         'dimensions:system-size-z = %.1f !nm' % (mcbs)]
        
        # Critical temperature Monte Carlo calculation
        input_script += ['sim:integrator = monte-carlo',
                         'sim:program = curie-temperature']
        
        # Default Monte Carlo params
        input_script += ['sim:equilibration-time-steps = %d' % (equil_timesteps),
                         'sim:loop-time-steps = %d' % (mc_timesteps),
                         'sim:time-steps-increment = 1']
        
        # Do Monte Carlo between +- 400 K from MFT estimate
        mft_t = self.mft_t
        delta_t = 400
        max_t = round(mft_t + delta_t)
        if mft_t - delta_t > 0:
            min_t = round(mft_t - delta_t)
        else:
            min_t = 0

        # Temperature range of simulation
        input_script += ['sim:minimum-temperature = %d' % (min_t),
                         'sim:maximum-temperature = %d' % (max_t),
                         'sim:temperature-increment = 25']
        
        # Output to save
        input_script += ['output:temperature',
                         'output:mean-magnetisation-length']
        
        input_script = "\n".join(input_script)
        
        with open('input', 'w') as f:
            f.write(input_script)
        
    def _create_ucf(self, structure):
        
        mat_name = self.mat_name
        tol = self.tol
        dists = self.dists

        abc = structure.lattice.abc
        ucx, ucy, ucz = abc[0], abc[1], abc[2]
        
        ucf = ['# Unit cell size:']
        ucf += ['%.10f %.10f %.10f' % (ucx, ucy, ucz)]
        
        ucf += ['# Unit cell lattice vectors:']
        a1 = list(structure.lattice.matrix[0])
        ucf += ['%.10f %.10f %.10f' % (a1[0], a1[1], a1[2])]
        a2 = list(structure.lattice.matrix[1])
        ucf += ['%.10f %.10f %.10f' % (a2[0], a2[1], a2[2])]
        a3 = list(structure.lattice.matrix[2])
        ucf += ['%.10f %.10f %.10f' % (a3[0], a3[1], a3[2])]
        
        nmats = len(self.unique_site_ids)
        
        ucf += ['# Atoms num_materials; id cx cy cz mat cat hcat']
        ucf += ['%d %d' % (len(structure), nmats)]
        
        # Fractional coordinates of atoms
        for i, r in enumerate(structure.frac_coords):
            for key in self.unique_site_ids:
                if i in key:
                    mat_id = self.unique_site_ids[key]
            ucf += ['%d %.10f %.10f %.10f %d 0 0' % (i, r[0], r[1], r[2], mat_id)]
            
        # J_ij exchange interaction matrix
        sgraph = self.sgraph
        ninter = 0
        for i, node in enumerate(sgraph.graph.nodes):
            ninter += sgraph.get_coordination_of_site(i)

        ucf += ['# Interactions']
        ucf += ['%d isotropic' % (ninter)]
        
        iid = 0  # counts number of interaction
        for i, node in enumerate(sgraph.graph.nodes):
            connections = sgraph.get_connected_sites(i)
            for c in connections:
                jimage = c[1]  # relative integer coordinates of atom j
                dx = jimage[0]
                dy = jimage[1]
                dz = jimage[2]
                j = c[2]  # index of neighbor
                dist = round(c[-1], 2)

                # Look up J_ij between the sites
                j_exc = self.hm._get_j_exc(i, j, dist)

                # Convert J_ij from meV to Joules
                j_exc *= 1.6021766e-22
                j_exc = str(j_exc)  # otherwise this rounds to 0

                ucf += ['%d %d %d %d %d %d %s' % (iid, i, j, dx, dy, dz, j_exc)]
                iid += 1
        
        ucf = '\n'.join(ucf)
        ucf_file_name = mat_name + '.ucf'
        
        with open(ucf_file_name, 'w') as f:
            f.write(ucf)

  
class VampireOutput:
    
    def __init__(self, vamp_stdout):
        """
        This class processes results from a Vampire Monte Carlo simulation
        and returns the critical temperature.
        
        Args:
            vamp_stdout (txt file): stdout from running vampire-serial.

        Parameters:
            critical_temp (float): Monte Carlo Tc result.
        """
        
        self.vamp_stdout = vamp_stdout
        self.critical_temp = np.nan
        self._parse_stdout(vamp_stdout)
        
    def _parse_stdout(self, vamp_stdout):
        
        # Parsing vampire MC output
        df = pd.read_csv(vamp_stdout,sep='\t',skiprows=8,header=None,
                         names=['T','m','nan'])
        df.drop('nan', axis=1, inplace=True)

        # Compute derivative to get critical temp
        grad = np.gradient(df['m'], df['T'])
        T_crit = df.iloc[grad.argmin()]['T']
        
        self.critical_temp = T_crit