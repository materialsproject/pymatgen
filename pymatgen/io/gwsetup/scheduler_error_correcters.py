"""
Error handlers for errors originating from the Submission systems.
"""

from __future__ import division

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "Mar 24, 2014"


import re
from abc import ABCMeta, abstractmethod, abstractproperty
from pymatgen.io.gwsetup.scheduler_error_parsers import AbstractError
from pymatgen.io.gwsetup.scheduler_error_parsers import TimeCancelError, MemoryCancelError, NodeFailureError



