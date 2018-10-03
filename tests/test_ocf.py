import unittest
import numpy as np
import pykbi

class TestOCF(unittest.TestCase):

    def setUp(self):
        n = 10
        self.rdf = pykbi.RDF(np.linspace(0.1,1.1, n), np.ones(n), closed=False)
        self.rdf.Integrate()


    def test_extrapolation(self):
        self.rdf.FindValues()
        self.assertAlmostEqual(self.rdf.ReturnKBI(),0.0)


    def test_not_extrapolated(self):
        # if the FindValues function has not been called, we should get None
        self.assertIsNone(self.rdf.ReturnKBI())


    def TearDown(self):
        self.rdf = None


class TestRDF_Closed(unittest.TestCase):

    def setUp(self):
        n = 10
        self.rdf = pykbi.RDF(np.linspace(0.1, 1.1, n), np.ones(n), closed=True)
        self.rdf.Integrate()

    def test_extrapolation(self):
        self.rdf.FindValues([2.0, 3.0])
        self.assertAlmostEqual(self.rdf.ReturnKBI(), 0.0)

    def test_not_extrapolated(self):
        # if the FindValues function has not been called, we should get None
        self.assertIsNone(self.rdf.ReturnKBI())

    def TearDown(self):
        self.rdf = None



if __name__ == "__main__":
    unittest.main()
