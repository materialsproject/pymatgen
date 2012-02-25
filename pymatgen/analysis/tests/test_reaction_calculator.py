import unittest

from pymatgen.core.structure import Composition
from pymatgen.analysis.reaction_calculator import Reaction, BalancedReaction, ReactionError

class ReactionTest(unittest.TestCase):

    #def setUp(self):
    #    pass
    
    def test_init(self):

        reactants = [Composition.from_formula("Fe"), Composition.from_formula("O2")]
        products = [Composition.from_formula("Fe2O3")]
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "2.000 Fe + 1.500 O2 -> 1.000 Fe2O3", "Wrong reaction obtained!")
        self.assertEquals(rxn.normalized_repr, "4 Fe + 3 O2 -> 2 Fe2O3", "Wrong normalized reaction obtained!")
        
        reactants = [Composition.from_formula("Fe"), Composition.from_formula("O"), Composition.from_formula("Mn"),Composition.from_formula("P")]
        products = [Composition.from_formula("FeP"),Composition.from_formula("MnO")]
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "1.000 Fe + 0.500 O2 + 1.000 Mn + 1.000 P -> 1.000 FeP + 1.000 MnO", "Wrong reaction obtained!")
        self.assertEquals(rxn.normalized_repr, "2 Fe + O2 + 2 Mn + 2 P -> 2 FeP + 2 MnO", "Wrong normalized reaction obtained!")
        reactants = [Composition.from_formula("FePO4")]
        products = [Composition.from_formula("FePO4")]
        
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "1.000 FePO4 -> 1.000 FePO4", "Wrong reaction obtained!")
        
        products = [Composition.from_formula("Ti2 O4"), Composition.from_formula("O1")]
        reactants = [Composition.from_formula("Ti1 O2")]
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "2.000 TiO2 -> 2.000 TiO2", "Wrong reaction obtained!")
        
        reactants = [Composition.from_formula("FePO4"), Composition.from_formula("Li")]
        products = [Composition.from_formula("LiFePO4")]
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "1.000 FePO4 + 1.000 Li -> 1.000 LiFePO4", "Wrong reaction obtained!")
        
        reactants = [Composition.from_formula("MgO")]
        products = [Composition.from_formula("MgO")]
        
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "1.000 MgO -> 1.000 MgO", "Wrong reaction obtained!")
        
        reactants = [Composition.from_formula("Mg")]
        products = [Composition.from_formula("Mg")]
        
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "1.000 Mg -> 1.000 Mg", "Wrong reaction obtained!")
        
        reactants = [Composition.from_formula("FePO4"),Composition.from_formula("LiPO3")]
        products = [Composition.from_formula("LiFeP2O7")]
        
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "1.000 FePO4 + 1.000 LiPO3 -> 1.000 LiFeP2O7", "Wrong reaction obtained!")
        
        reactants = [Composition.from_formula("Na"), Composition.from_formula("K2O")]
        products = [Composition.from_formula("Na2O"), Composition.from_formula("K")]
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "1.000 Na + 0.500 K2O -> 0.500 Na2O + 1.000 K", "Wrong reaction obtained!")
        
        # Test for an old bug which has a problem when excess product is defined.
        products = [Composition.from_formula("FePO4"), Composition.from_formula("O")]
        reactants = [Composition.from_formula("FePO4")]
        rxn = Reaction(reactants, products)
        
        self.assertEquals(str(rxn), "1.000 FePO4 -> 1.000 FePO4", "Wrong reaction obtained!")
        
        products = map(Composition.from_formula, ['La8Ti8O12', 'O2', 'LiCrO2'])
        reactants = [Composition.from_formula('LiLa3Ti3CrO12')]
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "1.000 LiLa3Ti3CrO12 -> 1.500 La2Ti2O3 + 2.750 O2 + 1.000 LiCrO2", "Wrong reaction obtained!")
         
    def test_normalize_to(self):
        products = [Composition.from_formula("Fe"), Composition.from_formula("O2")]
        reactants = [Composition.from_formula("Fe2O3")]
        rxn = Reaction(reactants, products)
        rxn.normalize_to(Composition.from_formula("Fe"), 3)
        self.assertEquals(str(rxn), "1.500 Fe2O3 -> 3.000 Fe + 2.250 O2", "Wrong reaction obtained!")
    
    def test_calculate_energy(self):
               
        reactants = [Composition.from_formula("MgO"), Composition.from_formula("Al2O3")]
        products = [Composition.from_formula("MgAl2O4")]
        energies = {Composition.from_formula("MgO"):-0.1, Composition.from_formula("Al2O3"):-0.2, Composition.from_formula("MgAl2O4"):-0.5}      
        rxn = Reaction(reactants, products)
        self.assertEquals(str(rxn), "1.000 MgO + 1.000 Al2O3 -> 1.000 MgAl2O4", "Wrong reaction obtained!")
        self.assertEquals(rxn.normalized_repr, "MgO + Al2O3 -> MgAl2O4", "Wrong normalized reaction obtained!")
        self.assertAlmostEquals(rxn.calculate_energy(energies), -0.2,5)
    
    
    def test_products_reactants(self):    
        reactants = [Composition.from_formula("Li3Fe2(PO4)3"), Composition.from_formula("Fe2O3"), Composition.from_formula("O2")]
        products = [Composition.from_formula("LiFePO4")]
        energies = {Composition.from_formula("Li3Fe2(PO4)3"):-0.1, Composition.from_formula("Fe2O3"):-0.2, Composition.from_formula("O2"):-0.2, Composition.from_formula("LiFePO4"):-0.5}
        rxn = Reaction(reactants, products)

        self.assertIn(Composition.from_formula("O2"), rxn.products, "O not in products!")
        self.assertIn(Composition.from_formula("Li3Fe2(PO4)3"), rxn.reactants, "Li3Fe2(PO4)4 not in reactants!")
        self.assertEquals(str(rxn), "0.333 Li3Fe2(PO4)3 + 0.167 Fe2O3 -> 0.250 O2 + 1.000 LiFePO4", "Wrong reaction obtained!")
        self.assertEquals(rxn.normalized_repr, "4 Li3Fe2(PO4)3 + 2 Fe2O3 -> 3 O2 + 12 LiFePO4", "Wrong normalized reaction obtained!")
        self.assertAlmostEquals(rxn.calculate_energy(energies), -0.48333333,5)


class BalancedReactionTest(unittest.TestCase):

    def test_init(self):
        rct = {Composition.from_formula('K2SO4'):3, Composition.from_formula('Na2S'):1, Composition.from_formula('Li'):24}
        prod = {Composition.from_formula('KNaS'): 2, Composition.from_formula('K2S'):2, Composition.from_formula('Li2O'):12}
        rxn = BalancedReaction(rct, prod)
        self.assertEquals(str(rxn), "24.000 Li + 1.000 Na2S + 3.000 K2SO4 -> 12.000 Li2O + 2.000 K2S + 2.000 KNaS")

        #Test unbalanced exception
        rct = {Composition.from_formula('K2SO4'):1, Composition.from_formula('Na2S'):1, Composition.from_formula('Li'):24}
        prod = {Composition.from_formula('KNaS'): 2, Composition.from_formula('K2S'):2, Composition.from_formula('Li2O'):12}
        self.assertRaises(ReactionError, BalancedReaction, rct, prod)
        
if __name__ == '__main__':
    unittest.main()

