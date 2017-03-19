import unittest

from CelestialMechanics.orbits import parable


class MyTestCase(unittest.TestCase):
    def test_parable(self):
        self.assertAlmostEqual(0.603086, parable.r(0.301543, 90.), places=6)


if __name__ == '__main__':
    unittest.main()