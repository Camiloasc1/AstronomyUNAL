import unittest

from CelestialMechanics.canonical_units.canonical_units import *


class MyTestCase(unittest.TestCase):
    def test_something(self):
        m2 = 0.25
        L = lagrange_points(m2)
        self.assertAlmostEqual(-1.10317, L[0], places=5)
        self.assertAlmostEqual(0.36074, L[1], places=5)
        self.assertAlmostEqual(1.26585, L[2], places=4)

        self.assertAlmostEqual(687.16, T(4317.6), places=1)


if __name__ == '__main__':
    unittest.main()
