import unittest
import numpy as np
import pykbi

class TestRDF_Open(unittest.TestCase):

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

    def test_integral_type(self):
        self.assertEqual(self.rdf.integral_type, "open")

    def test_eqint(self):
        self.assertIsNone(self.rdf.eqint)

    def test_lt(self):
        self.assertIsNone(self.rdf.lt)

    def test_volume(self):
        self.assertIsNone(self.rdf.volume)

    def test_AddingBox(self):
        self.rdf.AddBoxSize(10)
        self.assertAlmostEqual(self.rdf.volume,1000)


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

    def test_integral_type(self):
        self.assertEqual(self.rdf.integral_type, "closed")


    def test_eqint(self):
        self.assertIsNone(self.rdf.eqint)

    def test_lt(self):
        self.assertIsNone(self.rdf.lt)

    def test_volume(self):
        self.assertIsNone(self.rdf.volume)

    def test_AddingBox(self):
        self.rdf.AddBoxSize(10)
        self.assertAlmostEqual(self.rdf.volume,1000)

    def TearDown(self):
        self.rdf = None


class TestRDF_Initiating(unittest.TestCase):

    def setUp(self):
        pass

    def test_wrong_input(self):
        self.assertRaises(TypeError, lambda: pykbi.RDF(np.array([]), 1.0))
        self.assertRaises(TypeError, lambda: pykbi.RDF(1.0, np.array([])))

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
