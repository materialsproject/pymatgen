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

__all_errors__ = ['SubmitError', 'FullQueueError', 'DiskError', 'TimeCancelError', 'MemoryCancelError',
                  'NodeFailureError']

import re
import abc
from abc import ABCMeta


class AbstractError():
    """
    Error base class
    """
    __metaclass__ = ABCMeta

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


class MemoryCancelError(AbstractError):
    """
    Error due to exceeding the memory limit for the job.
      .limit will return a list of limits that were broken, None if it could not be determined.
    """

    @property
    def limit(self):
        return self.meta_data.get('broken_limit')


class NodeFailureError(AbstractError):
    """
    Error due the hardware failure of a specific node.
     .node will return a list of problematic nodes, None if it could not be determined.
    """

    @property
    def node(self):
        return self.meta_data.get('node')


class AbstractErrorParser():
    """
    Abstract class for parsing errors originating from the scheduler system and error that are not reported by the
    program itself, i.e. segmentation faults.

    A concrete implementation of this class for a specific sceduler needs a class attribute ERRORS for containing a
    dictionary specifying error:

    ERRORS = {ErrorClass: {'filespecifier' : {'string': "the string to be lookedfor",
                                              'metafilter': "string specifing the regular exprecion to obtian the meta data"}}

    """
    __metaclass__ = ABCMeta

    def __init__(self, err_file, out_file=None, run_err_file=None, batch_err_file=None):
        self.files = {'err': err_file, 'out': out_file, 'run_err': run_err_file, 'batch_err': batch_err_file}
        self.errors = []
        return

    @abc.abstractproperty
    def error_definitions(self):
        return {}

    @staticmethod
    def extract_metadata(lines, metafilter):
        meta_dict = {}
        for key in metafilter.keys():
            values = []
            for line in lines:
                match = re.match(metafilter[key][0], line)
                if match is not None:
                    values.append(re.match(metafilter[key][0], line).group(metafilter[key][1]))
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
            print 'parsing ', self.files[k], ' for ', errmsg[k]['string']
            try:
                with open(self.files[k], mode='r') as f:
                    lines = f.read().split('\n')
                for line in lines:
                    if errmsg[k]['string'] in line:
                        message = line
                        found = True
                if found:
                    metadata = self.extract_metadata(lines, errmsg[k]['metafilter'])
            except (IOError, OSError):
                print self.files[k], 'not found'

        return found, message, metadata

    def parse(self):
        """
        Parse for the occurens of all errors defined in ERRORS
        """
        for my_error in self.error_definitions:
            result = self.parse_single(self.error_definitions[my_error])
            if result[0]:
                self.errors.append(my_error(result[1], result[2]))


class SlurmErrorParser(AbstractErrorParser):
    """
    Implementation for the Slurm scheduler
    """

    @property
    def error_definitions(self):
        return {
            SubmitError: {
                'batch_err': {
                    'string': "sbatch: error: Batch job submission failed:",
                    'metafilter': {}
                }
            },
            FullQueueError: {
                'batch_err': {
                    'string': "sbatch: error: Batch job submission failed: Job violates accounting/QOS policy",
                    'metafilter': {}
                }
            },
            MemoryCancelError: {
                'err': {
                    'string': "Exceeded job memory limit",
                    'metafilter': {}
                }
            },
            TimeCancelError: {
                'err': {
                    'string': "Exceeded job time limit",
                    'metafilter': {}
                }
            },
            NodeFailureError: {
                'run_err': {
                    'string': "can't open /dev/ipath, network down",
                    'metafilter': {
                        'node': [r"node(\d+)\.(\d+)can't open (\S*), network down \(err=26\)", 1]
                    }
                }
            },
            AbstractError: {
                'out': {
                    'string': "a string to be found",
                    'metafilter': {}
                }
            }
        }


class PBSErrorParse(AbstractErrorParser):
    """
    Implementation for the PBS scheduler
    """


ALL_PARSERS = {'slurm': SlurmErrorParser, 'pbs': PBSErrorParse}


def get_parser(scheduler, err_file, out_file=None, run_err_file=None, batch_err_file=None):
    """
    Factory function to provide the parser for the specified scheduler. If the scheduler is not implemented None is
    returned. The files correspond to file names of the out and err files:
    err_file        stderr of the scheduler
    out_file        stdout of the scheduler
    run_err_file    stderr of the application
    batch_err_file  stderr of the submission
    """
    cls = ALL_PARSERS.get(scheduler)

    return cls(err_file, out_file, run_err_file, batch_err_file)


if __name__ == "__main__":
    my_parser = get_parser('slurm', err_file='queue.err', out_file='queue.out', run_err_file='run.err', batch_err_file='sbatch.err')
    my_parser.parse()
    print 'parser.errors', my_parser.errors
    for error in my_parser.errors:
        print error
        print error.node