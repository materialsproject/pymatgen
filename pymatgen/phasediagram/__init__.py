"""
The phasediagram package implements the analysis tools to perform phase
stability analyses, including the constructing of phase diagrams, determination
of decomposition products, etc. The package is designed to be fairly modular
and standalone.
"""

__author__ = "Shyue"
__date__ = "Mar 28 2013"

from pymatgen.phasediagram.pdmaker import PhaseDiagram, \
    GrandPotentialPhaseDiagram
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.phasediagram.plotter import PDPlotter
