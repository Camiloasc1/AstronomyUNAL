import unittest

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from CelestialMechanics.mu import mu_gm1m2, mu_sun
from CelestialMechanics.orbital_elements.orbital_elements import SolveE, solve_ellipse, solve, solve_hyperbola, \
    solve_parable, eclip2ecu
import CelestialMechanics.orbits.hyperbola as hyperbola
from CelestialMechanics.orbits import ellipse


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
        from CelestialMechanics.kepler import constants
        mu = mu_sun(1. / constants.Neptune)
        t_r = Time('2010-05-30T00:00:00Z', format='isot', scale='utc').jd * u.d
        t = Time('2010-07-02T00:00:00Z', format='isot', scale='utc').jd * u.d
        r, angle, r1, r_angle1 = solve_ellipse(a, e, M_r, mu, t_r, t)
        self.assertAlmostEqual(30.019361, r.value, places=3)
        self.assertAlmostEqual(-58.478691, angle.to(u.deg).value, places=4)
        (x, y, z), (x1, y1, z1) = solve(a, e, W, w, i, M_r, mu, t_r, t)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)

        # 5.2
        a = 1.6734521 * u.au
        e = 0.2907116
        i = 0 * u.deg
        W = 0 * u.deg
        w = 0 * u.deg
        M_r = 336.21581 * u.deg
        mu = mu_sun(0)
        t_r = Time('2013-04-18T00:00:00Z', format='isot', scale='utc').jd * u.d
        t = Time('2013-07-28T00:00:00Z', format='isot', scale='utc').jd * u.d
        r, angle, r1, r_angle1 = solve_ellipse(a, e, M_r, mu, t_r, t)
        self.assertAlmostEqual(1.2551615, r.value, places=4)
        self.assertAlmostEqual(40.645550, angle.to(u.deg).value, places=2)
        (x, y, z), (x1, y1, z1) = solve(a, e, W, w, i, M_r, mu, t_r, t)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)

        # 5.3
        a = 0.7233289 * u.au
        e = 0.0067692
        i = 0 * u.deg
        W = 0 * u.deg
        w = 0 * u.deg
        M_r = ellipse.angle_Mr(273.435890 * u.deg, 131.833001 * u.deg)
        self.assertAlmostEqual(141.602889, M_r.value)
        from CelestialMechanics.kepler import constants
        mu = mu_sun(1. / constants.Venus)
        t1 = Time('2013-09-30T20:30:15Z', format='isot', scale='utc').tt.jd * u.d + 5 * u.h
        t = Time('2013-10-01T01:30:15Z', format='isot', scale='utc').tt.jd * u.d
        self.assertAlmostEqual(t.value, t1.value)
        self.assertAlmostEqual(2456566.563451, t.value, places=6)
        self.assertAlmostEqual(2456566.563451, t1.value, places=6)
        t_r = Time('2013-09-10T00:00:00', format='isot', scale='tt').jd * u.d
        self.assertAlmostEqual(21.0634512, (t - t_r).value, places=7)
        r, angle, r1, r_angle1 = solve_ellipse(a, e, M_r, mu, t_r, t)
        self.assertAlmostEqual(0.728209, r.value, places=6)
        # self.assertAlmostEqual(176.411893, angle.to(u.deg).value, places=6)  # error in the book
        self.assertAlmostEqual(175.411893, angle.to(u.deg).value, places=4)  # error in the book
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
        self.assertAlmostEqual(-6.748472, angle.to(u.deg).value, places=3)

        (x, y, z), (x1, y1, z1) = solve(q, 1, W, w, i, None, mu, t_r, t)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)

    def test_chapter_5(self):
        # 5.10
        a = 0.2714161 * u.au
        e = 1.1418219
        i = 0 * u.deg
        W = 0 * u.deg
        w = 0 * u.deg
        mu = mu_sun(0)
        t_r = Time('2009-03-02T00:00:00Z', format='isot', scale='utc').jd * u.d + 0.5791 * u.d
        t = Time('2009-06-03T00:00:00Z', format='isot', scale='utc').jd * u.d
        r, angle, r1, r_angle1 = solve_hyperbola(a, e, mu, t_r, t)
        self.assertAlmostEqual(3.670381, r.value, places=5)
        self.assertAlmostEqual(148.883814, angle.to(u.deg).value, places=5)
        (x, y, z), (x1, y1, z1) = solve(a, e, W, w, i, None, mu, t_r, t)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)

    def test_chapter_6(self):
        from CelestialMechanics.kepler import constants
        # 5.6
        t_r = Time('2013-08-16T00:00:00Z', format='isot', scale='utc').jd * u.d
        t = Time('2013-12-01T00:00:00Z', format='isot', scale='utc').jd * u.d

        # Venus
        a = 0.72332709 * u.au
        e = 0.00676854
        i = 3.39464972 * u.deg
        W = 76.64385978 * u.deg
        w = 55.16203171 * u.deg
        M_r = 101.57663166 * u.deg
        mu = mu_sun(1. / constants.Venus)
        r, angle, r1, r_angle1 = solve_ellipse(a, e, M_r, mu, t_r, t)
        self.assertAlmostEqual(0.7231033934034238, r.value)
        self.assertAlmostEqual(-87.7684496261713, angle.to(u.deg).value)
        (x, y, z), (x1, y1, z1) = solve(a, e, W, w, i, M_r, mu, t_r, t)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)
        venus = np.array([x.to(u.au).value, y.to(u.au).value, z.to(u.au).value]) * u.au

        # Earth-Moon
        a = 1.00000225 * u.au
        e = 0.01667663
        i = 0.00176991 * u.deg
        W = 178.07928611 * u.deg
        w = 284.87149742 * u.deg
        M_r = 221.40784487 * u.deg
        mu = mu_sun(1. / constants.Earth_Moon)
        r, angle, r1, r_angle1 = solve_ellipse(a, e, M_r, mu, t_r, t)
        self.assertAlmostEqual(0.9861218868042078, r.value)
        self.assertAlmostEqual(-34.195316641554065, angle.to(u.deg).value)
        (x, y, z), (x1, y1, z1) = solve(a, e, W, w, i, M_r, mu, t_r, t)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)
        earth_moon = np.array([x.to(u.au).value, y.to(u.au).value, z.to(u.au).value]) * u.au

        (x, y, z) = venus - earth_moon
        x, y, z = eclip2ecu(x, y, z)
        self.assertAlmostEqual(0.44725, np.sqrt(x * x + y * y + z * z).value, places=5)  # rho
        coord = SkyCoord(x=x, y=y, z=z, unit='au', representation='cartesian')
        coord.representation = 'spherical'

        spected = SkyCoord('19h33m52.34s', '-24d40m42.43s', distance=0.44725 * u.au)
        self.assertAlmostEqual(spected.ra.value, coord.ra.value, places=5)
        self.assertAlmostEqual(spected.ra.value, coord.spherical.lon.value, places=5)
        self.assertAlmostEqual(spected.dec.value, coord.dec.value, places=5)
        self.assertAlmostEqual(spected.dec.value, coord.spherical.lat.value, places=5)
        self.assertAlmostEqual(spected.distance.value, coord.spherical.distance.value, places=5)  # rho

        # print(coord)
        # print(coord.to_string('hmsdms'))

    def test_eclip2ecu(self):
        x = 1.0559276 * u.au
        y = 1.0054148 * u.au
        z = -0.0048456 * u.au
        x, y, z = eclip2ecu(x, y, z)
        self.assertAlmostEqual(1.0559276, x.value)
        self.assertAlmostEqual(0.924377593, y.value)
        self.assertAlmostEqual(0.3954851, z.value)


if __name__ == '__main__':
    unittest.main()
