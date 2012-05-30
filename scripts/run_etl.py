#!/usr/bin/env python
"""
Perform Extract-Transform-Load operations.
"""
__author__ = "Dan Gunter"
__copyright__ = "Copyright 2012, The Materials Project"
__maintainer__ = "Dan Gunter"
__email__ = "dkgunter@lbl.gov"
__date__ = "25 May 2012"

# System imports
import argparse
import os
import sys

# Package imports
from pymatgen.etl import base as etl

## Functions

def _build_conf(args):
    """Build a configuration dictionary from the command-line
    arguments.
    """
    conf = { 
        "sources" : [ {
            "collection" : args.src_coll,
            "module" : args.src_mod,
            "class" : args.src_class,
        } ],
        "target" : {
            "port" : args.tgt_port,
            "db" : args.tgt_db,
            "collection" : args.tgt_coll,
        }
    }
    return conf

def main(cmdline=None):
    """Program entry point.
    """
    if cmdline is None:
        cmdline = sys.argv[1:]
    desc = __doc__
    parser = argparse.ArgumentParser(
                    usage="\n  %(prog)s -c config-file\n"
                    "  %(prog)s -s/--source COLL -t/--target COLL [options]",
                    description=desc)
    parser.add_argument("-c", "--config", dest="cfg_file", default=None,
                       metavar="FILE", type=argparse.FileType('r'), 
                       help="Configuration file (*)")
    parser.add_argument("-d", "--db", dest="tgt_db",
                        metavar="NAME", default="test",
                        help="Database name (%(default)s)")
    parser.add_argument("-k", "--class", dest="src_class", default="ETL",
                        metavar="CLASS", help="Class to instantiate "
                        "in source module (%(default)s)")
    parser.add_argument("-m", "--module", dest="src_mod",
                        metavar="MOD", help="Source module (*)")
    parser.add_argument("-p", "--eport", dest="tgt_port", type=int,
                        metavar="NUM", default=27017,
                        help="Database port (%(default)d)")
    parser.add_argument("-s", "--source", dest="src_coll", default=None,
                        metavar="COLL", help="Source collection (*)")
    parser.add_argument("-t", "--target", dest="tgt_coll", default=None,
                        metavar="COLL", help="Target collection (*)")
    args = parser.parse_args(cmdline)
    # Config. file
    if args.cfg_file is not None:
        try:
            conf = open(args.cfg_file, "r")
        except IOError, err:
            parser.error("Failed to open '{f}': {e}".format(
                            f=args.cfg_file, e=err))
    elif args.src_coll is None or args.tgt_coll is None or \
         args.src_mod is None:
        parser.error("missing required arguments")
    else:
        conf = _build_conf(args)
    # Run
    runner = etl.ETLRunner(conf)
    result = 0
    try:
        runner.run()
    except etl.ETLError, err:
        print("Error: {e}".format(e=err)) # XXX: should log
        result = 1
    # Done
    return result

if __name__ == "__main__":
    sys.exit(main())