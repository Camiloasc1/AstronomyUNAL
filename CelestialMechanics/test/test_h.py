import unittest
import numpy as np
from numpy.linalg import norm

from CelestialMechanics.h import h_and_angle


class MyTestCase(unittest.TestCase):
    def test_h_and_angle(self):
        r = [-0.05956988, -0.46122922, -0.03221964]
        self.assertAlmostEqual(0.46617493, norm(r), places=8)

        r1 = [0.02225882, -0.00217097, -0.00221965]
        self.assertAlmostEqual(0.02247432, norm(r1), places=8)

        h, angle1, angle2, angle_quadrant = h_and_angle(r, r1)
        self.assertAlmostEqual(0.00095382, h[0], places=8)
        self.assertAlmostEqual(-0.00084940, h[1], places=8)
        self.assertAlmostEqual(0.01039574, h[2], places=8)

        self.assertAlmostEqual(0.01047390, norm(h), places=7)
        self.assertAlmostEqual(88.61416619, angle1, places=2)
        self.assertAlmostEqual(91.38583381, angle2, places=2)
        self.assertAlmostEqual(91.38583381, angle_quadrant, places=2)

        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
