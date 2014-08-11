"""
This module contains useful decorators for a variety of functions.
"""
from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Dec 31, 2011"

import logging

from functools import wraps


def enable_logging(main):
    """
    This decorator is used to decorate main functions.
    It adds the initialization of the logger and an argument parser that allows one to select the loglevel.
    Useful if we are writing simple main functions that call libraries where the logging module is used

    Args:
        main:
            main function.
    """
    @wraps(main)
    def wrapper(*args, **kwargs):
        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('--loglevel', default="ERROR", type=str,
                            help="set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG")

        options = parser.parse_args()

        # loglevel is bound to the string value obtained from the command line argument.
        # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
        numeric_level = getattr(logging, options.loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % options.loglevel)
        logging.basicConfig(level=numeric_level)

        retcode = main(*args, **kwargs)
        return retcode

    return wrapper



