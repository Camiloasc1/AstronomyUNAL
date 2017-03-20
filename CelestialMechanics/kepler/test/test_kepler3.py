import unittest

from astropy import constants as astroconst
from astropy import units as u
from CelestialMechanics.kepler import constants
from CelestialMechanics.kepler.kepler3 import *


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertAlmostEqual(10759.22, T_sun(9.5549, 1. / constants.Saturn), delta=30)  # from wikipedia
        self.assertAlmostEqual(10713.8599, T_sun(9.512027, 1. / constants.Saturn), places=4)
        self.assertAlmostEqual(0.03360135, n_sun(9.512027, 1. / constants.Saturn), places=7)
        self.assertAlmostEqual(0.03360135, n(T_sun(9.512027, 1. / constants.Saturn)), places=7)
        self.assertAlmostEqual(n_sun(9.512027, 1. / constants.Saturn), n(T_sun(9.512027, 1. / constants.Saturn)))

        T_4_11 = 70.85 * u.a.to(u.d)
        self.assertAlmostEqual(25878.413, T_4_11, delta=1)
        self.assertAlmostEqual(17.122, a_sun(T_4_11, 0), places=3)

        a_4_12 = a((21.0874 * u.d).to(u.s), 1.1 * astroconst.M_sun.to(u.kg), 18 * astroconst.M_jup.to(u.kg)).value
        self.assertAlmostEqual(2.322E10, a_4_12, delta=1E10)
        self.assertAlmostEqual(0.155, a_4_12 * u.m.to(u.au), places=3)


if __name__ == '__main__':
    unittest.main()
