# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

"""
Error handlers for errors originating from the Submission systems.
"""

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "May 2014"

from pymatgen.io.abinit.scheduler_error_parsers import get_parser
try:
    from custodian.custodian import ErrorHandler
except ImportError:
    ErrorHandler = object


# TODO (from SP): Pls move this somewhere else. Custodian and Workflow stuff
# really shouldn't be in pymatgen.
class SchedulerErrorHandler(ErrorHandler):
    """
    Custodian error handler for scheduler related errors
      scheduler_adapter takes the scheduler, it should at least provide a .name attribute indentifying the scheduler,
      currently 'slurm' is supported.
      If the scheduler adapter also provides the methods defined in CorrectorProtocolScheduler, problems can also be
      fixed by .apply_corrections.
      If a application_adapter is also provided and it provides the methods defined in CorrectorProtocolApplication
      problems can also be fixed a the level of the application, e.g. making the application require less memory.
    """
    def __init__(self, scheduler_adapter, application_adapter=None, err_file='queue.err', out_file='queue.out',
                 run_err_file='run.err', batch_err_file='batch.err'):
        self.scheduler_adapter = scheduler_adapter
        self.application_adapter = application_adapter
        self.err_file = err_file
        self.out_file = out_file
        self.run_err_file = run_err_file
        self.batch_err_file = batch_err_file
        self.errors = []
        self.corrections = {}

    def check(self):
        """
        Check for the defined errors, put all found errors in self.errors, return True if any were found False if no
        errors were found
        """
        parser = get_parser(self.scheduler_adapter.name, err_file=self.err_file, out_file=self.out_file,
                            run_err_file=self.run_err_file, batch_err_file=self.batch_err_file)
        parser.parse()
        self.errors = parser.errors
        if len(self.errors) == 0:
            return False
        else:
            return True

    def correct(self):
        """
        For custodian compatibility
        """
        self.return_corrections()

    def return_corrections(self):

        for error in self.errors:
            self.corrections.update({error: {'scheduler_adapter_solutions': [], 'aplication_adapter_solutions': []}})
            self.corrections[error]['scheduler_adapter_solutions'].append(error.scheduler_adapter_solutions)
            self.corrections[error]['application_adapter_solutions'].append(error.application_adapter_solutions)
        return self.corrections

    def apply_corrections(self):
        """
        Method to directly apply the corrections.
        """
        for error in self.errors:
            for solution in error.scheduler_adapter_solutions:
                if self.scheduler_adapter is not None:
                    if self.scheduler_adapter.__getattribut__(solution[0].__name__)(solution[1]):
                        return True
            for solution in error.application_adapter_solutions:
                if self.application_adapter is not None:
                    if self.application_adapter.__getattribut__(solution[0].__name__)(solution[1]):
                        return True
        return False
