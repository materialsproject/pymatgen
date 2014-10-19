"""
DFT Surface Energy Error Correction for LDA and GGA exchange-correlation
Functionals. Ref: Phys. Rev. B 73, 195123, 2006.
"""

from __future__ import division, unicode_literals

__author__ = "Bharat Medasani"
__version__ = "0.1"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Alpha"
__date__ = "10/5/2014"

import math

ergpercmsq_to_evperangsq = 6.24150934e-5

unit_surf_energy_error = {
        "LDA":{"A":448.454, "B":-55.845}, 
        "PW91":{"A":1577.2, "B":-231.29}, 
        "PBE":{"A":1193.7, "B":-174.37}}
bohr_rad = 5.2917721092e-1

def surf_energy_error(xc_func, atom_no, valence, volume):
    """
    Computes the Unit Surface Energy Error for given functional and atom.
    
    Args:
        xc_func: Exchange Correlation Functional.
            Options are "LDA", "PBE", "PW91"
        atom_no: Number of atoms in cell
        valence: Valence of each atom in cell
        Volume: Volume of the cell
    
    Returns:
        Unit Surface Energy Error in eV/Ang^2
    """

    if xc_func not in unit_surf_energy_error.keys():
        raise ValueError("Surface energy error for functional not defined.")

    blk_electron_den = valence*atom_no/volume
    rs = (3/(4*math.pi*blk_electron_den))**(1/3.0)
    rs = rs/bohr_rad
    rspa = rs**(-5.0/2)
    rspb = rs**(-3.0/2)
    a = unit_surf_energy_error[functional]['A']
    b = unit_surf_energy_error[functional]['B']
    corr = (a*rspa + b*rspb) 
    corr = corr * ergpercmsq_to_evperangsq
    return corr

