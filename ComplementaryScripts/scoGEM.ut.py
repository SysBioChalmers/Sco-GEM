import unittest
import cobra

class TestScoGem(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.model = cobra.io.read_sbml_model("../ModelFiles/xml/scoGEM.xml")

    def test_read(self):
        self.assertTrue(bool(self.model))

    def test_metabolite_and_reactions_names(self):
        for r in self.model.reactions:
            self.assertIsInstance(r.name, str, msg = "Reaction {0} has no name".format(r.id))
        for m in self.model.metabolites:
            self.assertIsInstance(m.name, str, msg = "Metabolite {0} has no name".format(m.id))

    def test_growth(self):
        solution = self.model.optimize()
        self.assertAlmostEqual(solution.objective_value, 0.075, delta = 0.002)

if __name__ == '__main__':
    unittest.main()