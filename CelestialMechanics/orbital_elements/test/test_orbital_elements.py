import unittest

import numpy as np

from astropy import units as u
from astropy.time import Time

from CelestialMechanics.mu import mu_gm1m2, mu_sun
from CelestialMechanics.orbital_elements.orbital_elements import SolveE, solve_ellipse, solve, solve_hyperbola, \
    solve_parable
import CelestialMechanics.orbits.hyperbola as hyperbola


class MyTestCase(unittest.TestCase):
    def test_E(self):
        # 5.1
        E = SolveE(30.191502 * u.au, 0.010746, 302.56659 * u.deg, None, 0, 0)
        self.assertAlmostEqual(302.044701, E.value, places=6)
        # 5.2
        E = SolveE(1.6734521 * u.au, 0.2907116, 22.199696 * u.deg, None, 0, 0)
        self.assertAlmostEqual(30.704771, E.value, places=5)
        # 5.3
        E = SolveE(0.7233289 * u.au, 0.0067692, 176.349528 * u.deg, None, 0, 0)
        # self.assertAlmostEqual(175.380763, E.value, places=5)

        # 5.4
        # t_r = Time('2013-03-10T03:56:31.2Z', format='isot', scale='utc').jd * u.d
        t_r = Time('2013-03-10T00:00:00Z', format='isot', scale='utc').jd * u.d + 0.16425 * u.d
        self.assertAlmostEqual(Time('2013-03-10T03:56:31.2Z', format='isot', scale='utc').jd, t_r.value)
        t = Time('2013-12-28T00:00:00Z', format='isot', scale='utc').jd * u.d
        q = 0.3016486 * u.au
        e = 1.2998762
        F = SolveE(hyperbola.a(q, e), e, None, mu_sun(0), t_r, t)
        self.assertAlmostEqual(2.445043, F.value, places=6)

    def test_ellipse(self):
        # 5.1
        a = 30.191502 * u.au
        e = 0.010746
        i = 0 * u.deg
        W = 0 * u.deg
        w = 0 * u.deg
        M_r = 302.370534 * u.deg
        mu = mu_sun(1. / 19412.26)
        t_r = Time('2010-05-30T00:00:00Z', format='isot', scale='utc').jd * u.d
        t = Time('2010-07-02T00:00:00Z', format='isot', scale='utc').jd * u.d
        r, angle, r1, r_angle1 = solve_ellipse(a, e, M_r, mu, t_r, t)
        self.assertAlmostEqual(30.019361, r.value, places=3)
        self.assertAlmostEqual(-58.478691, angle.to(u.deg).value, places=4)

        (x, y, z), (x1, y1, z1) = solve(a, e, W, w, i, M_r, mu, t_r, t)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)

    def test_hyperbola(self):
        # 5.4
        q = 0.3016486 * u.au
        e = 1.2998762
        a = hyperbola.a(q, e)
        i = 0 * u.deg
        W = 0 * u.deg
        w = 0 * u.deg
        mu = mu_sun(0)
        t_r = Time('2013-03-10T00:00:00Z', format='isot', scale='utc').jd * u.d + 0.16425 * u.d
        t = Time('2013-12-28T00:00:00Z', format='isot', scale='utc').jd * u.d
        r, angle, r1, r_angle1 = solve_hyperbola(a, e, mu, t_r, t)
        self.assertAlmostEqual(6.589548, r.value, places=5)
        self.assertAlmostEqual(133.496551, angle.to(u.deg).value, places=5)

        (x, y, z), (x1, y1, z1) = solve(a, e, W, w, i, None, mu, t_r, t)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)

    def test_parable(self):
        # 5.5
        q = 5.1723060 * u.au
        i = 0 * u.deg
        W = 0 * u.deg
        w = 0 * u.deg
        mu = mu_sun(0)
        t_r = Time('2015-01-24T00:00:00Z', format='isot', scale='utc').jd * u.d + 0.08039 * u.d
        t = Time('2014-11-28T00:00:00Z', format='isot', scale='utc').jd * u.d
        r, angle, r1, r_angle1 = solve_parable(q, mu, t_r, t)
        self.assertAlmostEqual(5.190286, r.value, places=5)
        self.assertAlmostEqual(-6.748472, angle.to(u.deg).value, places=5)

        (x, y, z), (x1, y1, z1) = solve(a, e, W, w, i, None, mu, t_r, t)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)


if __name__ == '__main__':
    unittest.main()
