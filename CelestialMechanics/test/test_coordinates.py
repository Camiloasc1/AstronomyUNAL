import unittest

from astropy import units as u

from CelestialMechanics.coordinates import eclip2ecu, ecu2eclip


class MyTestCase(unittest.TestCase):
    def test_something(self):
        x = 1.0559276 * u.au
        y = 1.0054148 * u.au
        z = -0.0048456 * u.au
        x, y, z = eclip2ecu(x, y, z)
        self.assertAlmostEqual(1.0559276, x.value)
        self.assertAlmostEqual(0.924377593, y.value)
        self.assertAlmostEqual(0.3954851, z.value)

        x, y, z = ecu2eclip(x, y, z)
        self.assertAlmostEqual(1.0559276, x.value)
        self.assertAlmostEqual(1.0054148, y.value)
        self.assertAlmostEqual(-0.0048456, z.value)


if __name__ == '__main__':
    unittest.main()
