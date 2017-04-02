import unittest

from astropy import constants as astroconst
from astropy import units as u
from astropy.time import Time

from CelestialMechanics.kepler import constants
from CelestialMechanics.orbits import ellipse
from CelestialMechanics.orbits.ellipse import delta_t_t0_aeangle


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
        vq, vQ = ellipse.vqQ_sun(a, e, 0)
        self.assertAlmostEqual(0.031478, vq.value, places=2)
        self.assertAlmostEqual(0.000528, vQ.value, places=2)
        self.assertAlmostEqual(54.50, vq.to(u.km / u.s).value, places=2)
        self.assertAlmostEqual(0.91, vQ.to(u.km / u.s).value, places=2)

        # 4.15
        a = astroconst.R_earth + 560 * u.km
        self.assertAlmostEqual(7569.5, ellipse.v(a, a, astroconst.M_earth, 0).value, delta=20)

    def test_chapter_5(self):
        from CelestialMechanics.mu import mu_sun

        # 5.5
        t0 = Time('2014-01-03T00:00:00Z', format='isot', scale='utc').jd * u.d + 0.633 * u.d
        t1 = Time('2014-04-03T00:00:00Z', format='isot', scale='utc').jd * u.d + 0.9 * u.d
        t2 = Time('2014-10-05T00:00:00Z', format='isot', scale='utc').jd * u.d + 0.5 * u.d

        a = 1 * u.au
        e = 0.01669
        r = 1 * u.au
        mu = mu_sun(1 / constants.Earth_Moon)

        angles = ellipse.angle(a, e, r)
        self.assertAlmostEqual(90.9563109612867, angles[0].value)
        self.assertAlmostEqual(269.0436890387133, angles[1].value)

        delta_t_t0 = delta_t_t0_aeangle(a, e, angles[0], mu)
        self.assertAlmostEqual((t1 - t0).value, delta_t_t0.value, delta=0.1)
        delta_t_t0 = delta_t_t0_aeangle(a, e, angles[1], mu)
        self.assertAlmostEqual((t2 - t0).value, delta_t_t0.value, delta=0.1)


if __name__ == '__main__':
    unittest.main()
