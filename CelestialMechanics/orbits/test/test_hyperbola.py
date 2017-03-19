import unittest

from CelestialMechanics.orbits.hyperbola import angle_asymptotic


class MyTestCase(unittest.TestCase):
    def test_something(self):
        angle1, angle2 = angle_asymptotic(1.5)
        self.assertAlmostEqual(131.81031489577862, angle1)
        self.assertAlmostEqual(228.18968510422138, angle2)


if __name__ == '__main__':
    unittest.main()
