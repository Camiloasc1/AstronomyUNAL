import unittest

from astropy import constants as astroconst
from astropy import units as u
from CelestialMechanics.orbits import ellipse
from CelestialMechanics.kepler import constants


class MyTestCase(unittest.TestCase):
    def test_ellipse(self):
        r = ellipse.r(1.5236164, 0.0932802, 32.)
        self.assertAlmostEqual(1.3996391, r, places=7)
        a, e = ellipse.ae(0.275, 1.168)
        self.assertAlmostEqual(0.722, a, places=3)
        self.assertAlmostEqual(0.618, e, places=2)

        sun = astroconst.M_sun
        mercury = astroconst.M_sun / constants.Mercury
        energy = ellipse.E((0.38709 * u.au).to(u.m), sun, mercury)
        self.assertAlmostEqual(-3.817E32, energy.value, delta=1E32)

        # 4.14
        a = 17.8 * u.au
        e = 0.967
        q, Q = ellipse.qQ(a, e)
        self.assertAlmostEqual(0.031478, ellipse.v_sun(q, a, 0).value, places=5)
        self.assertAlmostEqual(0.000528, ellipse.v_sun(Q, a, 0).value, places=5)
        self.assertAlmostEqual(54.50, ellipse.v_sun(q, a, 0).to(u.km / u.s).value, places=2)
        self.assertAlmostEqual(0.91, ellipse.v_sun(Q, a, 0).to(u.km / u.s).value, places=2)

        # 4.15
        a = astroconst.R_earth + 560 * u.km
        self.assertAlmostEqual(7569.5, ellipse.v(a, a, astroconst.M_earth, 0).value, delta=20)


if __name__ == '__main__':
    unittest.main()
