import unittest

from astropy import units as u
from astropy import constants as astroconst
from CelestialMechanics.orbits import parable


class MyTestCase(unittest.TestCase):
    def test_parable(self):
        self.assertAlmostEqual(0.603086, parable.r(0.301543, 90.), places=6)
        self.assertAlmostEqual(120, parable.angle(0.125 * u.au, 0.5 * u.au)[0].value)

        # 4.15
        a = astroconst.R_earth + 560 * u.km
        self.assertAlmostEqual(10704.7, parable.v(a, astroconst.M_earth, 0).value, delta=20)


if __name__ == '__main__':
    unittest.main()
