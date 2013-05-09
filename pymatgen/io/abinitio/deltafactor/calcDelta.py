#!/usr/bin/env python
#
# calcDelta.py determines the Delta factor of a code. Current version: 1.1
#
# Copyright (C) 2012 Kurt Lejaeghere <Kurt.Lejaeghere@UGent.be>, Center for
# Molecular Modeling (CMM), Ghent University, Ghent, Belgium
#
# calcDelta.py is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your
# option) any later version.
#
# In addition to the regulations of the GNU Lesser General Public License,
# publications and communications based in parts on this program or on
# parts of this program are required to cite the following articles:
#
# "Error estimates for solid-state density-functional theory predictions: an
# overview by means of the ground-state elemental crystals", K. Lejaeghere,
# V. Van Speybroeck, G. Van Oost, and S. Cottenier (2012), to be published.
#
# calcDelta.py is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with calcDelta.py; if not, see <http://www.gnu.org/licenses/>.

# Python and numpy are required to use this script.

import numpy as np
from sys import argv, stdout, exit


if __name__ == "__main__":
    usage = '''\
    calcDelta.py -- Support script to calculate the Delta factor of a code

    Use: python calcDelta.py infile [reffile] [--stdout]
         where reffile and infile refer to files containing the element name, V0,
         B0, and B1 information (V0 in A^3/atom, B0 in GPa, and B1 dimensionless)
         in columns
         This command calculates the Delta factor of the code in infile compared
         to the one in reffile. When reffile is not explicitly given, WIEN2k.txt
         is used by default.
         Additional output is printed in Delta-out.txt. The option --stdout can be
         used to explicitly print all elements to standard output (on screen)
         instead.
         Attention: the presence of WIEN2k.txt in the same folder is required for
         full functionality!
    python calcDelta.py --help displays the current instructions
    '''

    if len(argv) not in [2, 3, 4]:
        print 'Error: Wrong number of arguments'
        exit()
    if argv[1] == '--help':
        print usage
        exit()

    reffile = 'WIEN2k.txt'
    if len(argv) > 2 and argv[2] != '--stdout':
        reffile = argv[2]

    data_f = np.loadtxt(argv[1], dtype={'names':('element', 'V0', 'B0', 'BP'),
        'formats':('S2', np.float, np.float, np.float)})
    data_w = np.loadtxt(reffile, dtype={'names':('element', 'V0', 'B0', 'BP'),
        'formats':('S2', np.float, np.float, np.float)})

    try:
        len(data_f['element'])
    except TypeError:
        print 'Error: ' + argv[1] + ': at least two elements required'
        exit()
    eloverlap = list(set(data_f['element']) & set(data_w['element']))

    v0w = np.zeros(len(eloverlap))
    b0w = np.zeros(len(eloverlap))
    b1w = np.zeros(len(eloverlap))

    v0f = np.zeros(len(eloverlap))
    b0f = np.zeros(len(eloverlap))
    b1f = np.zeros(len(eloverlap))

    elw = list(data_w['element'])
    elf = list(data_f['element'])

    for i in range(len(eloverlap)):
        searchnr = elw.index(eloverlap[i])
        v0w[i] = data_w['V0'][searchnr]
        b0w[i] = data_w['B0'][searchnr] * 10.**9. / 1.602176565e-19 / 10.**30.
        b1w[i] = data_w['BP'][searchnr]

        searchnr = elf.index(eloverlap[i])
        v0f[i] = data_f['V0'][searchnr]
        b0f[i] = data_f['B0'][searchnr] * 10.**9. / 1.602176565e-19 / 10.**30.
        b1f[i] = data_f['BP'][searchnr]

    Vi = 0.94 * v0w
    Vf = 1.06 * v0w

    a3f = 9. * v0f**3. * b0f / 16. * (b1f - 4.)
    a2f = 9. * v0f**(7./3.) * b0f / 16. * (14. - 3. * b1f)
    a1f = 9. * v0f**(5./3.) * b0f / 16. * (3. * b1f - 16.)
    a0f = 9. * v0f * b0f / 16. * (6. - b1f)

    a3w = 9. * v0w**3. * b0w / 16. * (b1w - 4.)
    a2w = 9. * v0w**(7./3.) * b0w / 16. * (14. - 3. * b1w)
    a1w = 9. * v0w**(5./3.) * b0w / 16. * (3. * b1w - 16.)
    a0w = 9. * v0w * b0w / 16. * (6. - b1w)

    x = [0, 0, 0, 0, 0, 0, 0]

    x[0] = (a0f - a0w)**2
    x[1] = 6. * (a1f - a1w) * (a0f - a0w)
    x[2] = -3. * (2. * (a2f - a2w) * (a0f - a0w) + (a1f - a1w)**2.)
    x[3] = -2. * (a3f - a3w) * (a0f - a0w) - 2. * (a2f - a2w) * (a1f - a1w)
    x[4] = -3./5. * (2. * (a3f - a3w) * (a1f - a1w) + (a2f - a2w)**2.)
    x[5] = -6./7. * (a3f - a3w) * (a2f - a2w)
    x[6] = -1./3. * (a3f - a3w)**2.

    Fi = np.zeros_like(Vi)
    Ff = np.zeros_like(Vf)

    for n in range(7):
        Fi = Fi + x[n] * Vi**(-(2.*n-3.)/3.)
        Ff = Ff + x[n] * Vf**(-(2.*n-3.)/3.)

    Delta1 = 1000. * np.sqrt((Ff - Fi) / (0.12 * v0w))
    Dmax = Delta1.argmax()
    Dmin = Delta1.argmin()
    Delta = Delta1.mean()
    Dstdev = Delta1.std()
    total = len(v0w)

    outfile = stdout
    if argv[-1] != '--stdout':
        outfile = open('Delta-out.txt', 'w')
    elementlist = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na',
                   'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',
                   'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
                   'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
                   'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
                   'I', 'Xe', 'Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
                   'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Rn']

    outfile.write('--------------------\n')
    outfile.write('# Delta-factor of ' + argv[1] + ' with respect to ' +
        reffile + ' (in meV/atom)\n')
    outfile.write('# (%i elements of %i included)\n' % (total, len(elementlist)))
    outfile.write('# calculated with calcDelta.py version 1.1 \n')
    outfile.write('--------------------\n')

    for el in elementlist:
        while True:
            try:
                i = eloverlap.index(el)
                outfile.write(eloverlap[i] + '\t %.3f \n' % Delta1[i])
                break
            except ValueError:
                outfile.write(el + '\t N/A \n')
            break

    outfile.write('--------------------\n')
    outfile.write('np.mean  %.3f\n' % Delta)
    outfile.write('np.std   %.3f\n' % Dstdev)
    outfile.write('np.max   %.3f  (%.2s)\n' % (Delta1[Dmax], eloverlap[Dmax]))
    outfile.write('np.min   %.3f  (%.2s)\n' % (Delta1[Dmin], eloverlap[Dmin]))
    outfile.write('--------------------\n')

    outfile.close()
