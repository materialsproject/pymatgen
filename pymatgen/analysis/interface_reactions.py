from __future__ import division

"""
This module provides InterfacialReactivity class for analyzing interfacial reactions 
and plotting reaction energy vs. mixing ratio between two reactants.

Class Methods:

    -- get_kinks: Gets information on critical reactions: [index of kink, critical
                  mixing ratio x, critical energy, balanced reaction].
    
    -- get_products: Gets all potential reaction products in interfacial reactions.
    
    -- labels: Dictionary version of class method get_kinks.

    -- plot: Plots the reaction energy as a function of mixing ratios of reaction 
             between reactants. Kinks are labeled with indices.

    -- minimum: Gets the minimum reaction energy E_min and corresponding mixing ratio
                X_min: (x_min, energy_min).

    -- get_no_mixing_energy: Gets only the reaction energies of reactants from 
                             equilibration with the external reservoir.

"""

__author__ = "Yihan Xiao"
__email__ = "eric.xyh2011@gmail.com"
__date__ = "Aug 15 2017"

import numpy as np
import matplotlib.pylab as plt

from pymatgen.phasediagram.maker import GrandPotentialPhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen import Composition

PDAnalyzer.numerical_tol = 1e-12

class InterfacialReactivity:
    """
    A class for analyzing interfacial reactions and plotting reaction energy vs.
    mixing ratio between two reactants.
    """

    def __init__(self, c1, c2, pd, norm=1, include_no_mixing_energy=0,
                 pd_non_grand=None):
        '''
        Initiates the class.

        :param c1, c2: Composition objects for reactants.
        :param pd: PhaseDiagram object or GrandPotentialPhaseDiagram object built
                from elements in c1 and c2.
        :param norm: Whether or not the total number of atoms in each reactant
                will be normalized to 1 for reaction energy calculations.
        :param include_no_mixing_energy: No_mixing_energy for a reactant is the
                opposite number of its energy above grand potential convex hull. In
                cases where reactions involve elements reservoir, this param
                determines whether no_mixing_energy of reactants will be included
                in the final reaction energy calculation. By definition, if pd is
                not a GrandPotentialPhaseDiagram object, this param is 0.
        :param pd_non_grand: PhaseDiagram object but not GrandPotentialPhaseDiagram
                object, built from elements in c1 and c2.
        '''
        self.grand = isinstance(pd,GrandPotentialPhaseDiagram)

        # if pd is not a GrandPotentialPhaseDiagram object, include_no_mixing_energy should be 0.
        assert not (include_no_mixing_energy and not self.grand), \
            'Please provide grand phase diagram to compute no_mixing_energy for reactants!'
        # if include_no_mixing_energy = 1, pd should be a GrandPotentialPhaseDiagram object,
        # pd_non_grand should be a normal PhaseDiagram object.
        assert not (include_no_mixing_energy and not pd_non_grand),\
            'Please also provide non-grand phase diagram to compute no_mixing_energy for reactants!'

        # Keeps copy of original compositions.
        self.c1_original = c1
        self.c2_original = c2

        # Two sets of composition attributes for different processing conditions:
        # normalization with and without exluding element(s) from reservoir.
        self.c1 = c1
        self.c2 = c2
        self.comp1 = c1
        self.comp2 = c2

        self.norm = norm
        self.pd = pd
        self.pda = PDAnalyzer(pd)
        if pd_non_grand:
            self.pd_non_grand = pd_non_grand

        # Factor is the compositional ratio between composition self.c1 and processed
        #  composition self.comp1. E.g., factor for Composition('SiO2') and Composition('O')
        # is 2. This factor will be used to convert mixing ratio in self.comp1 - self.comp2
        # tie line to that in self.c1 - self.c2 tie line.
        self.factor1, self.factor2 = 1, 1

        if self.grand:
            # Excludes element(s) from reservoir.
            self.comp1 = Composition({k:v for k,v in c1.items() if k not in pd.chempots})
            self.comp2 = Composition({k:v for k,v in c2.items() if k not in pd.chempots})
            # Calculate the factor for in case where self.grand =1 and self.norm =1.
            factor1 = self.comp1.num_atoms / c1.num_atoms
            factor2 = self.comp2.num_atoms / c2.num_atoms

        if self.norm:
                self.c1 = c1.fractional_composition
                self.c2 = c2.fractional_composition
                self.comp1 = self.comp1.fractional_composition
                self.comp2 = self.comp2.fractional_composition
                if self.grand:
                    # Only when self.grand =1 and self.norm =1 will self.factor be updated.
                    self.factor1 = factor1
                    self.factor2 = factor2

        # Computes energies for reactants in different scenarios.
        if not self.grand:
            # Use entry energy as reactant energy if no reservoir is present.
            self.e1 = self._get_entry_energy(self.pd, self.comp1)
            self.e2 = self._get_entry_energy(self.pd, self.comp2)
        else:
            if include_no_mixing_energy:
                # Computing grand potentials needs compositions containing element(s)
                # from reservoir, so self.c1 and self.c2 are used.
                self.e1 = self._get_grand_potential(self.c1)
                self.e2 = self._get_grand_potential(self.c2)
            else:
                self.e1 = self.pda.get_hull_energy(self.comp1)
                self.e2 = self.pda.get_hull_energy(self.comp2)


    def _get_entry_energy(self,pd,composition):
        """
        Returns the lowest entry energy for entries matching the COMPOSITION. Entries
        with non-negative formation energies are excluded.
        :param pd: PhaseDiagram object.
        :param composition: Composition object that the target entry should match.
        :return: The lowest entry energy among entries matching the COMPOSITION.
        """
        candidate = [i.energy_per_atom for i in pd.qhull_entries if
                     i.composition.fractional_composition == composition.fractional_composition]
        assert candidate != [], 'The reactant {} has no matching entry with negative' \
                                ' formation energy!'.format(composition.reduced_formula)
        min_entry_energy = min(candidate)
        return min_entry_energy * composition.num_atoms

    def _get_grand_potential(self,composition):
        '''
        Returns the grand potential Phi at a given composition and chemical potential(s).
         E.g., Phi[c, mu_{Li}]= E_{hull}[c] - n_{Li}[c]mu_{Li}.
        :param composition: Composition object.
        :return: Grand potential at a given composition and chemical potential(s).
        '''
        grand_potential = self._get_entry_energy(self.pd_non_grand, composition) - \
                          sum([composition[e] * mu for e, mu in self.pd.chempots.items()])
        if self.norm:
            # Normalizes energy to the composition excluding element(s) from reservoir.
            grand_potential /= (1 - sum([composition.get_atomic_fraction(e.symbol)
                                         for e, mu in self.pd.chempots.items()]))
        return grand_potential


    def _get_energy(self, x):
        '''
        Returns reaction energy at mixing ratio x:(1-x) for self.comp1:self.comp2.
        :param x: Mixing ratio x.
        :return: Reaction energy.
        '''
        return self.pda.get_hull_energy(self.comp1 * x + self.comp2 * (1-x)) - \
               self.e1 * x - self.e2 * (1-x)

    def _get_reaction(self, x, normalize=0):
        '''
        Returns balanced reaction at mixing ratio x:(1-x) for self.comp1:self.comp2.
        :param x: Mixing ratio x.
        :param normalize: Whether or not to normalize the sum of coefficients of reactants to 1.
                For unnormalized case, use original reactant compositions in reaction for clarity.
        :return: Reaction object.
        '''
        mix_comp = self.comp1 * x + self.comp2 * (1-x)
        decomp = self.pda.get_decomposition(mix_comp)

        if normalize:
            reactant = list(set([self.c1, self.c2]))
        else:
            # Uses original composition for reactants.
            reactant = list(set([self.c1_original, self.c2_original]))

        if self.grand:
            reactant += [Composition(e.symbol) for e, v in self.pd.chempots.items()]

        product = [Composition(k.name) for k, v in decomp.items()]
        reaction = Reaction(reactant, product)

        if normalize:
            x = self._convert(x,self.factor1,self.factor2)
            if x == 1:
                reaction.normalize_to(self.c1, x)
            else:
                reaction.normalize_to(self.c2, 1-x)
        return reaction

    def get_products(self):
        """
        Returns all potential products in interfacial reactions.
        :return: List of product formulas. E.g., ['Li', 'O2','Mn'].
        """
        products = set()
        for _,_,_,react in self.get_kinks():
            products = products.union(set([k.reduced_formula for k in react.products]))
        return list(products)

    def _convert(self,x, factor1,factor2):
        """
        Converts mixing ratio x in comp1 - comp2 tie line to that in c1 - c2 tie line.
        :param x: Mixing ratio x in comp1 - comp2 tie line.
        :param factor1: Compositional ratio between composition c1 and processed composition comp1.
                E.g., factor for Composition('SiO2') and Composition('O') is 2.
        :param factor2: Compositional ratio between composition c2 and processed composition comp2.
        :return: Mixing ratio in c1 - c2 tieline.
        """
        return x * factor2 / ((1-x) * factor1 + x * factor2)

    def _reverse_convert(self,x, factor1,factor2):
        """
        Converts mixing ratio x in c1 - c2 tie line to that in comp1 - comp2 tie line.
        :param x: Mixing ratio x in c1 - c2 tie line.
        :param factor1: Compositional ratio between composition c1 and processed composition comp1.
                E.g., factor for Composition('SiO2') and Composition('O') is 2.
        :param factor2: Compositional ratio between composition c2 and processed composition comp2.
        :return: Mixing ratio in comp1 - comp2 tie line.
        """
        return x * factor1 / ((1-x) * factor2 + x * factor1)

    def get_kinks(self):
        '''
        Finds all the kinks in mixing ratio where reaction products changes along the tie line of
        composition self.c1 and composition self.c2.
        :return: Zip object of tuples (index, mixing ratio, reaction energy, reaction).
        '''
        c1_coord = self.pd.pd_coords(self.comp1)
        c2_coord = self.pd.pd_coords(self.comp2)
        n1 = self.comp1.num_atoms
        n2 = self.comp2.num_atoms
        critical_comp = self.pda.get_critical_compositions(self.comp1,self.comp2)
        x_kink, energy_kink, react_kink = [], [], []
        if all(c1_coord == c2_coord):
            x_kink = [0,1]
            energy_kink = [self._get_energy(x) for x in x_kink]
            react_kink = [self._get_reaction(x) for x in x_kink]
        else:
            for i in reversed(critical_comp):
                # Gets mixing ratio x at kinks.
                c = self.pd.pd_coords(i)
                x = np.linalg.norm(c-c2_coord) / np.linalg.norm(c1_coord-c2_coord)
                # Modifies mixing ratio in case compositions self.comp1 and self.comp2 are not normalized.
                x = x * n2 / (n1 + x * (n2 - n1))
                # Converts mixing ratio in comp1 - comp2 tie line to that in c1 - c2 tie line.
                x_converted = self._convert(x,self.factor1,self.factor2)
                x_kink.append(x_converted)
                # Gets reaction energy at kinks
                energy_kink.append(self._get_energy(x))
                # Gets balanced reaction at kinks
                react_kink.append(self._get_reaction(x))
        index_kink = range(1, len(critical_comp)+1)
        return zip(index_kink, x_kink, energy_kink, react_kink)

    def labels(self):
        '''
        Returns a dictionary containing kink information:
        {index: 'x= mixing_ratio energy= reaction_energy reaction_equation'}.
        :return: Dictionary object. E.g., {1: 'x= 0.0 energy = 0.0 Mn -> Mn',
        2: 'x= 0.5 energy = -15.0 O2 + Mn -> MnO2', 3: 'x= 1.0 energy = 0.0 O2 -> O2'}.
        '''
        return {j:'x= '+str(round(x,4))+' energy = '+str(round(energy,4))+' '+str(reaction)
                for j, x, energy, reaction in self.get_kinks()}

    def plot(self):
        '''
        Plots reaction energy as a function of mixing ratio x in self.c1 - self.c2 tie line using pylab.
        :return: Plot of reaction energy as a function of mixing ratio x.
        '''
        plt.rcParams['xtick.major.pad'] = '6'
        plt.rcParams['ytick.major.pad'] = '6'
        plt.rcParams['axes.linewidth'] = 2
        npoint = 1000
        xs = np.linspace(0, 1, npoint)

        # Converts sampling points in self.c1 - self.c2 tie line to those in self.comp1 - self.comp2 tie line.
        xs_reverse_converted = self._reverse_convert(xs, self.factor1,self.factor2)
        energies = [self._get_energy(x) for x in xs_reverse_converted]
        plt.plot(xs, energies, 'k-')

        # Marks kinks and minimum energy point.
        kinks = self.get_kinks()
        _, x_kink, energy_kink, _, = zip(*kinks)
        plt.scatter(x_kink, energy_kink, marker='o', c='blue', s=20)
        plt.scatter(self.minimum()[0], self.minimum()[1], marker='*', c='red', s=300)

        # Labels kinks with indices. Labels are made draggable in case of overlapping.
        for index, x, energy, reaction in kinks:
            plt.annotate(
                index,
                xy=(x, energy), xytext=(5, 30),
                textcoords='offset points', ha='right', va='bottom',
                arrowprops= dict(arrowstyle = '->', connectionstyle='arc3,rad=0')).draggable()
        plt.xlim([-0.05,1.05])
        if self.norm:
            plt.ylabel('Energy (eV/atom)')
        else:
            plt.ylabel('Energy (eV/f.u.)')
        plt.xlabel('$x$ in $x$ {} + $(1-x)$ {}'.format(self.c1.reduced_formula,self.c2.reduced_formula))
        return plt

    def minimum(self):
        '''
        Returns the minimum reaction energy E_min and corresponding mixing ratio x_min.
        :return: Tuple (x_min, E_min).
        '''
        return min([(x,energy) for _,x,energy,_ in self.get_kinks()],key=lambda i:i[1])

    def get_no_mixing_energy(self):
        '''
        Returns the opposite number of energy above grand potential convex hull for both reactants.
        :return: [(reactant1, no_mixing_energy1),(reactant2,no_mixing_energy2)].
        '''
        assert self.grand == 1, 'Please provide grand potential phase diagram for computing no_mixing_energy!'

        energy1 = self.pda.get_hull_energy(self.comp1) - self._get_grand_potential(self.c1)
        energy2 = self.pda.get_hull_energy(self.comp2) - self._get_grand_potential(self.c2)
        unit = 'eV/f.u.'
        if self.norm:
            unit = 'eV/atom'
        return [(self.c1_original.reduced_formula + ' ({0})'.format(unit), energy1),
                (self.c2_original.reduced_formula + ' ({0})'.format(unit), energy2)]
