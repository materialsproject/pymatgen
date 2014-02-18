import time
import socket
import os
import sys


def now():
    """
    helper to return a time string
    """
    return time.strftime("%H:%M:%S %d/%m/%Y")


def get_cluster_params():
    """
    task independent function to obtain cluster specific values for the execution of VASP.
    """
    hostname = socket.gethostname()
    try:
        vasp_exe = os.environ['VASPEXE']
    except KeyError:
        print 'Error, could not find VASPEXE environment variable on cluster "{hn}" ... exit'.format(hn=hostname)
        sys.exit()
    try:
        vasp_modules = os.environ['VASP_MODULES'].split()
    except KeyError:
        print 'Error, could not find VASP_MODULES environment variable on cluster "{hn}" ... exit'.format(hn=hostname)
        sys.exit()
    try:
        vasp_mpi_executer = os.environ['VASP_MPI_EXECUTER']
    except KeyError:
        print 'Error, could not find VASP_MPI_EXECUTER environment variable on cluster "{hn}" ... exit'.format(hn=hostname)
        sys.exit()
    return {'vasp_exe': vasp_exe, 'vasp_modules': vasp_modules, 'vasp_mpi_executer': vasp_mpi_executer}