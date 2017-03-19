import unittest

from CelestialMechanics.orbits import ellipse


class MyTestCase(unittest.TestCase):
    def test_ellipse(self):
        r = ellipse.r(1.5236164, 0.0932802, 32.)
        self.assertAlmostEqual(1.3996391, r, places=7)
        a, e = ellipse.ae(0.275, 1.168)
        self.assertAlmostEqual(0.722, a, places=3)
        self.assertAlmostEqual(0.618, e, places=2)


if __name__ == '__main__':
    unittest.main()
