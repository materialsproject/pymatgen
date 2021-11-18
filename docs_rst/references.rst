References
==========

Some of pymatgen's functionality is based on scientific advances / principles
developed by various scientists. If you use some of these functionality in
your research, you may wish to consider citing the following works:

pymatgen.analysis.path_finder
-----------------------------

The path finder code, which finds diffusion paths through a structure based on
a given potential field, is written by the Ceder group at UC Berkley::

    Rong, Z., Kitchaev, D., Canepa, P., Huang, W., & Ceder, G. (2016).
    An efficient algorithm for finding the minimum energy path for cation
    migration in ionic materials. The Journal of Chemical Physics, 145(7),
    74112. doi:10.1063/1.4960790

pymatgen.core.surface and pymatgen.analysis.wulff
-------------------------------------------------

The surface generation code, which can automatically generate surfaces based
on any crystal, and the Wulff code, which plots the Wulff shape given a
crystal and surface energies, are written by the Materials Virtual Lab::

    Tran, R.; Xu, Z.; Radhakrishnan, B.; Winston, D.; Sun, W.; Persson, K. A.;
    Ong, S. P. Surface energies of elemental crystals, Sci. Data, 2016, 3,
    160080, doi:10.1038/sdata.2016.80.

and contains elements from the following publication::

    Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
    Surface Science, 2013, 617, 53–59, doi:10.1016/j.susc.2013.05.016.

pymatgen.io.vasp.sets
---------------------

The MIT parameter sets, which are optimized for high-throughput computing, are
outlined the following work::

    Jain, A.; Hautier, G.; Moore, C. J.; Ong, S. P.; Fischer, C. C.;
    Mueller, T.; Persson, K. A.; Ceder, G. A high-throughput infrastructure for
    density functional theory calculations, Comput. Mater. Sci., 2011, 50,
    2295–2310, doi:10.1016/j.commatsci.2011.02.023.

pymatgen.phasediagram
---------------------

The phase diagram code, in particular the grand canonical phase diagram
analysis, is based on the work of Ong et al. and are used in following works::

    Ong, S. P.; Wang, L.; Kang, B.; Ceder, G. Li−Fe−P−O2 Phase Diagram from
    First Principles Calculations, Chem. Mater., 2008, 20, 1798–1807,
    doi:10.1021/cm702327g.

    Ong, S. P.; Jain, A.; Hautier, G.; Kang, B.; Ceder, G. Thermal stabilities
    of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
    principles calculations, Electrochem. commun., 2010, 12, 427–430,
    doi:10.1016/j.elecom.2010.01.010.

pymatgen.entries.compatibility
------------------------------

The compatibility processing, which allows mixing of GGA and GGA+U runs that
have been calculated using the MaterialsProjectVaspInputSet or MITVaspInputSet,
is based on the following work::

    Jain, A.; Hautier, G.; Ong, S. P.; Moore, C. J.; Fischer, C. C.;
    Persson, K. A.; Ceder, G. Formation enthalpies by mixing GGA and GGA+U
    calculations, Phys. Rev. B, 2011, 84, 45115, doi:10.1103/PhysRevB.84.045115.

pymatgen.matproj
----------------

The matproj package contains an interface to the `Materials Project REST API
<http://www.materialsproject.org/open>`_ (Materials API). If you use data
from the Materials Project, please cite the following works::

    Jain, A.; Ong, S. P.; Hautier, G.; Chen, W.; Richards, W. D.; Dacek,
    S.; Cholia, S.; Gunter, D.; Skinner, D.; Ceder, G.; Persson, K. A.
    Commentary: The Materials Project: A materials genome approach to
    accelerating materials innovation, APL Mater., 2013, 1, 11002,
    doi:10.1063/1.4812323.

    Ong, S. P.; Cholia, S.; Jain, A.; Brafman, M.; Gunter, D.; Ceder, G.;
    Persson, K. a. The Materials Application Programming Interface (API): A
    simple, flexible and efficient API for materials data based on
    REpresentational State Transfer (REST) principles, Comput. Mater. Sci.,
    2015, 97, 209–215, doi:10.1016/j.commatsci.2014.10.037.

pymatgen.symmetry
-----------------

The symmetry package is based on the excellent spglib developed by Atz Togo. For
more information, please refer to Atz Togo's site at
http://spglib.sourceforge.net/.

pymatgen.command_line.bader_caller
----------------------------------

This module implements an interface to the Henkelmann et al.'s excellent
Fortran code for calculating a Bader charge analysis. Please cite the
following::

    Henkelman, G., Arnaldsson, A., & Jónsson, H. (2006). A fast and robust
    algorithm for Bader decomposition of charge density. Computational
    Materials Science, 36(3), 354–360. doi:10.1016/j.commatsci.2005.04.010

pymatgen.io.feff
----------------

This module implements an io interface for FEFF calculations. Please
acknowledge the contribution of Alan Dozier, UKY.

pymatgen.io.zeo
---------------

This implements an interface to the excellent Zeo++ code base. Please
consider citing the following publications::

    T.F. Willems, C.H. Rycroft, M. Kazi, J.C. Meza, and M. Haranczyk,
    Algorithms and tools for high-throughput geometry- based analysis of
    crystalline porous materials, Microporous and Mesoporous Materials,
    149 (2012) 134-141, `doi:10.1016/j.micromeso.2011.08.020
    <https://doi.org/10.1016/j.micromeso.2011.08.020>`_.

    R.L. Martin, B. Smit, and M. Haranczyk, Addressing challenges of
    identifying geometrically diverse sets of crystalline porous materials,
    J. Chem. Information and Modelling, `doi:10.1021/ci200386x
    <https://doi.org/10.1021/ci200386x>`_.
