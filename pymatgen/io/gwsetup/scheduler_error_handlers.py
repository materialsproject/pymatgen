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

from pymatgen.io.gwsetup.scheduler_error_parsers import get_parser
from pymatgen.io.gwsetup.scheduler_error_correcters import get_correcter
from custodian.custodian import ErrorHandler


class SchedulerErrorHandler(ErrorHandler):
    """
    Custodian error handler for scheduler related errors
      scheduler takes the schduler, currently: slurm is supported
    """
    def __init__(self, scheduler, code, err_file='queue.err', out_file='queue.out', run_err_file='run.err',
                 batch_err_file='batch.err'):
        self.scheduler = scheduler
        self.code = code
        self.err_file = err_file
        self.out_file = out_file
        self.run_err_file = run_err_file
        self.batch_err_file = batch_err_file
        self.errors = []
        self.corrections = {}

    def check(self):
        parser = get_parser(self.scheduler, err_file=self.err_file, out_file=self.out_file,
                            run_err_file=self.run_err_file, batch_err_file=self.batch_err_file)
        parser.parse()
        self.errors = parser.errors
        if len(self.errors) == 0:
            return False
        else:
            return True

    def correct(self):
        correcter = get_correcter(self.scheduler, self.code)
        correcter.correct()
        self.corrections = correcter.corrections
        return self.corrections