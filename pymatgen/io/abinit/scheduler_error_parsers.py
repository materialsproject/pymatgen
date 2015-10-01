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

__all_errors__ = ['SubmitError', 'FullQueueError', 'DiskError', 'TimeCancelError', 'MemoryCancelError',
                  'NodeFailureError']

import re
import abc
import six

from abc import ABCMeta, abstractproperty, abstractmethod

@six.add_metaclass(ABCMeta)
class CorrectorProtocolScheduler(object):
    """
    Abstract class to define the protocol / interface for correction operators. The client code quadapters / submission
    script generator method / ... should implement these methods.
    """

    @abstractproperty
    def name(self):
        return str()

    @abstractmethod
    def exclude_nodes(self, nodes):
        """
        Method to exclude certain nodes from being used in the calculation. It is called when a calculation seemed to
        have been crashed due to a hardware failure at the nodes specified.

            nodes: list of node numbers that were found to cause problems

        returns True is the memory could be increased False otherwise
        """
        return bool

    @abstractmethod
    def increase_mem(self):
        """
        Method to increase then memory in the calculation. It is called when a calculation seemed to have been crashed
        due to a insufficient memory.

        returns True is the memory could be increased False otherwise
        """
        return bool

    @abstractmethod
    def increase_time(self):
        """
        Method to increase te time for the calculation. It is called when a calculation seemed to
        have been crashed due to a time limit.

        returns True is the memory could be increased False otherwise
        """
        return bool

    @abstractmethod
    def increase_cpus(self):
        """
        Method to increse the number of cpus being used in the calculation. It is called when a calculation seemed to
        have been crashed due to time or memory limits being broken.

        returns True is the memory could be increased False otherwise
        """
        return bool


@six.add_metaclass(ABCMeta)
class CorrectorProtocolApplication(object):
    """
    Abstract class to define the protocol / interface for correction operators. The client code quadapters / submission
    script generator method / ... should implement these methods.
    """

    @abstractproperty
    def name(self):
        return str()

    @abstractmethod
    def decrease_mem(self):
        """
        Method to increase then memory in the calculation. It is called when a calculation seemed to have been crashed
        due to a insufficient memory.

        returns True is the memory could be increased False otherwise
        """
        return bool

    @abstractmethod
    def speed_up(self):
        """
        Method to speed_up the calculation. It is called when a calculation seemed to time limits being broken.

        returns True is the memory could be increased False otherwise
        """
        return bool


@six.add_metaclass(ABCMeta)
class AbstractError(object):
    """
    Error base class
    """

    def __init__(self, errmsg, meta_data):
        self.errmsg = errmsg
        self.meta_data = meta_data if meta_data is not None else {}

    def __str__(self):
        _message = '%s  %s\n' \
                   '  error message : %s \n' \
                   '  meta data     : %s' % (self.name, self.__doc__, self.errmsg, str(self.meta_data))
        return _message

    @property
    def name(self):
        return self.__class__.__name__

    @property
    def scheduler_adapter_solutions(self):
        """
        to be implemented by concrete errors returning a list of tuples defining corrections. The First element of the
          tuple should be a string of one of the methods in CorrectorProtocolScheduler, the second element should
          contain the arguments.
        """
        return []

    @property
    def application_adapter_solutions(self):
        """
        to be implemented by concrete errors returning a list of tuples defining corrections. The First element of the
          tuple should be a string of one of the methods in CorrectorProtocolApplication, the second element should
          contain the arguments.
        """
        return []

    def last_resort_solution(self):
        """
        what to do if every thing else fails...
        """
        print('non of the defined solutions for %s returned success...' % self.name)
        return


class SubmitError(AbstractError):
    """
    Errors occurring at submission. The limits on the cluster may have changed.
    """


class FullQueueError(AbstractError):
    """
    Errors occurring at submission. To many jobs in the queue / total cpus / .. .
    """


class DiskError(AbstractError):
    """
    Errors involving problems writing to disk.
    """


class TimeCancelError(AbstractError):
    """
    Error due to exceeding the time limit for the job.
      .limit will return a list of limits that were broken, None if it could not be determined.
    """

    @property
    def limit(self):
        return self.meta_data.get('broken_limit')

    @property
    def scheduler_adapter_solutions(self):
        return [(CorrectorProtocolScheduler.increase_time,)]

    @property
    def application_adapter_solutions(self):
        return [(CorrectorProtocolApplication.speed_up,)]


class MemoryCancelError(AbstractError):
    """
    Error due to exceeding the memory limit for the job.
      .limit will return a list of limits that were broken, None if it could not be determined.
    """

    @property
    def limit(self):
        return self.meta_data.get('broken_limit')

    @property
    def scheduler_adapter_solutions(self):
        return [(CorrectorProtocolScheduler.increase_mem,)]

    @property
    def application_adapter_solutions(self):
        return [(CorrectorProtocolApplication.decrease_mem,)]


class NodeFailureError(AbstractError):
    """
    Error due the hardware failure of a specific node.
     .node will return a list of problematic nodes, None if it could not be determined.
    """

    @property
    def nodes(self):
        return self.meta_data.get('nodes')

    @property
    def scheduler_adapter_solutions(self):
        return [(CorrectorProtocolScheduler.exclude_nodes, [self.nodes])]


@six.add_metaclass(ABCMeta)
class AbstractErrorParser(object):
    """
    Abstract class for parsing errors originating from the scheduler system and error that are not reported by the
    program itself, i.e. segmentation faults.

    A concrete implementation of this class for a specific scheduler needs a class attribute ERRORS for containing a
    dictionary specifying error:

    ERRORS = {ErrorClass: {
                'file_specifier' : {
                    'string': "the string to be looked for",
                    'meta_filter': "string specifing the regular expression to obtain the meta data"
                    }
                }

    """
    def __init__(self, err_file, out_file=None, run_err_file=None, batch_err_file=None):
        self.files = {'err': err_file, 'out': out_file, 'run_err': run_err_file, 'batch_err': batch_err_file}
        self.errors = []
        return

    @abc.abstractproperty
    def error_definitions(self):
        return {}

    @staticmethod
    def extract_metadata(lines, meta_filter):
        meta_dict = {}
        for key in meta_filter.keys():
            values = []
            for line in lines:
                match = re.match(meta_filter[key][0], line)
                if match is not None:
                    values.append(re.match(meta_filter[key][0], line).group(meta_filter[key][1]))
            values = sorted(set(values))
            meta_dict.update({key: values})
        return meta_dict

    def parse_single(self, errmsg):
        """
        Parse the provided files for the corresponding strings.
        """
        found = False
        message = None
        metadata = None
        for k in errmsg.keys():
            if self.files[k] is not None:
                # print 'parsing ', self.files[k], ' for ', errmsg[k]['string']
                try:
                    with open(self.files[k], mode='r') as f:
                        lines = f.read().split('\n')
                    for line in lines:
                        if errmsg[k]['string'] in line:
                            message = line
                            found = True
                    if found:
                        metadata = self.extract_metadata(lines, errmsg[k]['meta_filter'])
                except (IOError, OSError):
                    print(self.files[k], 'not found')
                    pass
                except TypeError:
                    print('type error', self.files[k], ' has type ', self.files[k].cls(), ' should be string.')
                    pass

        return found, message, metadata

    def parse(self):
        """
        Parse for the occurens of all errors defined in ERRORS
        """
        for error in self.error_definitions:
            result = self.parse_single(self.error_definitions[error])
            if result[0]:
                self.errors.append(error(result[1], result[2]))
        if len(self.errors) > 0:
            print('QUEUE_ERROR FOUND')
            for error in self.errors:
                print(error)


class SlurmErrorParser(AbstractErrorParser):
    """
    Implementation of the error definitions for the Slurm scheduler
    """

    @property
    def error_definitions(self):
        return {
            SubmitError: {
                'batch_err': {
                    'string': "Batch job submission failed",
                    'meta_filter': {}
                }
            },
            FullQueueError: {
                'batch_err': {
                    'string': "sbatch: error: Batch job submission failed: Job violates accounting/QOS policy",
                    'meta_filter': {}
                }
            },
            MemoryCancelError: {
                'err': {
                    'string': "Exceeded job memory limit",
                    'meta_filter': {}
                }
            },
            TimeCancelError: {
                'err': {
                    'string': "DUE TO TIME LIMIT",
                    'meta_filter': {
                        'time_of_cancel': [r"JOB (\d+) CANCELLED AT (\S*) DUE TO TIME LIMIT", 1]
                    }
                }
            },
            NodeFailureError: {
                'run_err': {
                    'string': "can't open /dev/ipath, network down",
                    'meta_filter': {
                        'nodes': [r"node(\d+)\.(\d+)can't open (\S*), network down \(err=26\)", 1]
                    }
                }
            },
            AbstractError: {
                'out': {
                    'string': "a string to be found",
                    'meta_filter': {}
                }
            }
        }


class PBSErrorParser(AbstractErrorParser):
    """
    Implementation for the PBS scheduler
    """

#=>> PBS: job killed: walltime 932 exceeded limit 900
#=>> PBS: job killed: walltime 46 exceeded limit 30
    @property
    def error_definitions(self):
        return {
            TimeCancelError: {
                'out': {
                    'string': "job killed: walltime",
                    'meta_filter': {
                        'broken_limit': [r"job killed: walltime (\d+) exceeded limit (\d+)", 1]
                    }
                }
            },
            AbstractError: {
                'out': {
                    'string': "a string to be found",
                    'meta_filter': {}
                }
            }
        }


ALL_PARSERS = {'slurm': SlurmErrorParser, 'pbspro': PBSErrorParser, 'torque': PBSErrorParser}


def get_parser(scheduler, err_file, out_file=None, run_err_file=None, batch_err_file=None):
    """
    Factory function to provide the parser for the specified scheduler. If the scheduler is not implemented None is
    returned. The files, string, correspond to file names of the out and err files:
    err_file        stderr of the scheduler
    out_file        stdout of the scheduler
    run_err_file    stderr of the application
    batch_err_file  stderr of the submission

    Returns:
        None if scheduler is not supported.
    """
    cls = ALL_PARSERS.get(scheduler)
    return cls if cls is None else cls(err_file, out_file, run_err_file, batch_err_file)


if __name__ == "__main__":
    my_parser = get_parser('pbs', err_file='queue.err', out_file='queue.out', run_err_file='run.err',
                           batch_err_file='sbatch.err')
    my_parser.parse()
    print('parser.errors', my_parser.errors)
    for my_error in my_parser.errors:
        print(my_error)
