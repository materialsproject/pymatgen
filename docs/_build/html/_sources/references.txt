References
==========

Some of pymatgen's functionality is based on scientific advances / principles
developed by various scientists. If you use some of these functionality in
your research, you may wish to consider citing the following works:

pymatgen.io.vaspio_set
----------------------

The MIT parameter sets, which are optimized for high-throughput computing, are
outlined the following work:

    A. Jain, G. Hautier, C. Moore, S. P. Ong, C. C. Fischer, T. Mueller,
    K. A. Persson, and G. Ceder. *A high-throughput infrastructure for density
    functional theory calculations.* Computational Materials Science, 2011,
    50(8), 2295-2310. `doi:10.1016/j.commatsci.2011.02.023
    <http://dx.doi.org/10.1016/j.commatsci.2011.02.023>`_

pymatgen.phasediagram
---------------------

The phase diagram code, in particular the grand canonical phase diagram
analysis, is based on the work of Ong et al. and are used in following works:

    S. P. Ong, L. Wang, B. Kang, and G. Ceder. *Li-Fe-P-O2 Phase Diagram from
    First Principles Calculations.* Chemistry of Materials, 2008, 20(5),
    1798-1807. `doi:10.1021/cm702327g <http://dx.doi.org/10.1021/cm702327g>`_

    S. P. Ong, A. Jain, G. Hautier, B. Kang, and G. Ceder. *Thermal stabilities
    of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
    principles calculations.* Electrochemistry Communications, 2010, 12(3),
    427-430. `doi:10.1016/j.elecom.2010.01.010
    <http://dx.doi.org/10.1016/j.elecom.2010.01.010>`_

pymatgen.entries.compatibility
------------------------------

The compatibility processing, which allows mixing of GGA and GGA+U runs that
have been calculated using the MaterialsProjectVaspInputSet or MITVaspInputSet,
is based on the following work:

    A. Jain, G. Hautier, S. P. Ong, C. Moore, C. C. Fischer, K. A. Persson, and
    G. Ceder. *Formation enthalpies by mixing GGA and GGA + U calculations.*
    Physical Review B, 2011, 84(4), 045115. `doi:10.1103/PhysRevB.84.045115
    <http://dx.doi.org/10.1103/PhysRevB.84.045115>`_

pymatgen.matproj
----------------

The matproj package contains an interface to the `Materials Project REST API
<http://www.materialsproject.org/open>`_ (Materials API). If you use data
from the Materials Project, please cite the following works:

    A. Jain, G. Hautier, C. Moore, S. P. Ong, C. Fischer, T. Mueller,
    K. Persson, G. Ceder. *A high-throughput infrastructure for density
    functional theory calculations.* Computational Materials Science, 2011,
    50(8), 2295–2310. `doi:10 .1016/j.commatsci.2011.02.023
    <http://dx.doi.org/10.1016/j.commatsci.2011.02.023>`_

    S. P. Ong, A. Jain, G. Hautier, M. Kocher, S. Cholia, D. Gunter, D. Bailey,
    D. Skinner, K. Persson, G. Ceder. *The Materials Project.*
    http://materialsproject.org/

pymatgen.symmetry
-----------------

The symmetry package is based on the excellent spglib developed by Atz Togo. For
more information, please refer to Atz Togo's site at
http://spglib.sourceforge.net/.

pymatgen.command_line.bader_caller
----------------------------------

This module implements an interface to the Henkelmann et al.'s excellent
Fortran code for calculating a Bader charge analysis. Please cite the
following:

    G. Henkelman, A. Arnaldsson, and H. Jonsson, "A fast and robust algorithm
    for Bader decomposition of charge density", Comput. Mater. Sci. 36,
    254-360 (2006).

pymatgen.io.feffio
------------------

This module implements an io interface for FEFF calculations. Please
acknowledge the contribution of Alan Dozier, UKY.

pymatgen.io.zeoio
-----------------

This implements an interface to the excellent Zeo++ code base. Please
consider citing the following publications:

    T.F. Willems, C.H. Rycroft, M. Kazi, J.C. Meza, and M. Haranczyk,
    Algorithms and tools for high-throughput geometry- based analysis of
    crystalline porous materials, Microporous and Mesoporous Materials,
    149 (2012) 134-141, `doi:10.1016/j.micromeso.2011.08.020
    <http://dx.doi.org/10.1016/j.micromeso.2011.08.020>`_.

    R.L. Martin, B. Smit, and M. Haranczyk, Addressing challenges of
    identifying geometrically diverse sets of crystalline porous materials,
    J. Chem. Information and Modelling, `doi:10.1021/ci200386x
    <http://dx.doi.org/10.1021/ci200386x>`_.
