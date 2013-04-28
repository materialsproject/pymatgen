# eosfit.py fits E(V) data to a Birch-Murnaghan equation of state. 
# Current version: 1.1
#
# Copyright (C) 2012 Kurt Lejaeghere <Kurt.Lejaeghere@UGent.be>, Center for
# Molecular Modeling (CMM), Ghent University, Ghent, Belgium
#
# eosfit.py is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# eosfit.py is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for 
# more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with eosfit.py; if not, see <http://www.gnu.org/licenses/>.

# The following code is based on the source code of eos.py from the Atomic 
# Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/>.

# Python and numpy are required to use this script.

import numpy as np
import sys 

def BM(energies):

    fitdata = np.polyfit(energies[:,0]**(-2./3.), energies[:,1], 3, full=True)
    ssr = fitdata[1]
    sst = np.sum((energies[:,1] - np.average(energies[:,1]))**2.)
    residuals0 = ssr/sst
    deriv0 = np.poly1d(fitdata[0])
    deriv1 = np.polyder(deriv0, 1)
    deriv2 = np.polyder(deriv1, 1)
    deriv3 = np.polyder(deriv2, 1)

    volume0 = 0
    x = 0
    for x in np.roots(deriv1):
        if x > 0 and deriv2(x) > 0:
            volume0 = x**(-3./2.)
            break

    if volume0 == 0:
        print('Error: No minimum could be found')
        exit()
    
    derivV2 = 4./9. * x**5. * deriv2(x)
    derivV3 = (-20./9. * x**(13./2.) * deriv2(x) -
        8./27. * x**(15./2.) * deriv3(x))
    bulk_modulus0 = derivV2 / x**(3./2.)
    bulk_deriv0 = -1 - x**(-3./2.) * derivV3 / derivV2

    return volume0, bulk_modulus0, bulk_deriv0, residuals0

if __name__ == "__main__":
    usage = '''\
    Use: python eosfit.py filename
        calculates the Birch-Murnaghan equation of state from a given file, 
        containing in its columns the volumes in A^3/atom and energies in eV/atom,
        respectively
        output is printed in filename.eosout
    --help gives an overview of all options
    '''

    if len(sys.argv) != 2:
        print 'Error: Wrong number of arguments'
        exit()
    else:
        if sys.argv[1] == '--help':
            print usage            
            exit()

    infile = sys.argv[1]

    data = np.loadtxt(infile)
    volume, bulk_modulus, bulk_deriv, residuals = BM(data)

    echarge = 1.60217733e-19

    outstr = '''\
    Equation Of State parameters - least squares fit of a Birch Murnaghan curve

    %.5f \t %.5f \t %.3f

       V0 \t \t  B0 \t \t  BP
    [A^3/at] \t [GPa] \t \t [--] 

    1-R^2: %f

    ''' % (volume, (bulk_modulus * echarge * 1.0e21), bulk_deriv, residuals[0])

    outfile = open(infile+'.eosout', 'w')
    outfile.write(outstr)

    outfile.close()
