import unittest
import numpy as np
import pykbi

class TestODF(unittest.TestCase):

    def setUp(self):
        n = 250
        radial = np.linspace(0.1, 1000.0, n)
        self.ocf = pykbi.odf(radial, 1.0)


    def test_long_distance(self):
        self.assertAlmostEqual(self.ocf[-1],1.0)

    def test_short_distance(self):
        self.assertAlmostEqual(self.ocf[0],0.0)

    def tearDown(self):
        self.ocf = None


if __name__ == "__main__":
    unittest.main()
