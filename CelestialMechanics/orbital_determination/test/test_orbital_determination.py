import unittest

import numpy as np
from astropy import units as u
from astropy.time import Time

from CelestialMechanics.coordinates import ecu2eclip_
from CelestialMechanics.mu import mu_sun
from CelestialMechanics.kepler import constants
from CelestialMechanics.orbital_determination.orbital_determination import solve, solve_gauss


class MyTestCase(unittest.TestCase):
    def test_examples(self):
        # 8.1
        r = [2.77904683, -4.28963554, -0.04438092] * u.au
        r1 = [0.00624498, 0.00446529, -0.00015828] * u.au / u.d
        t_r = Time('2009-01-09T00:00:00', format='isot', scale='tt').jd * u.d
        mu = mu_sun(1. / constants.Jupiter)

        a, e, i, W, w, M_r, t0 = solve(r, r1, mu, t_r)
        self.assertAlmostEqual(1.30376234, i.to(u.deg).value, places=8)
        self.assertAlmostEqual(100.50895502, W.to(u.deg).value, places=8)
        self.assertAlmostEqual(274.07925551, w.to(u.deg).value % 360, places=8)
        self.assertAlmostEqual(293.61066092, M_r.to(u.deg).value % 360, places=8)
        # print(Time(t0, format='jd', scale='utc').isot)

        # 8.2
        r = [0.0429740, 3.5483648, -5.0009781] * u.au
        r1 = [0.0069528, -0.000767, 0.0068981] * u.au / u.d
        t_r = Time('2014-02-15T00:00:00', format='isot', scale='tt').jd * u.d
        mu = mu_sun(0)

        a, e, i, W, w, M_r, t0 = solve(r, r1, mu, t_r)
        self.assertAlmostEqual(3.8289407, a.to(u.au).value, places=7)
        self.assertAlmostEqual(121.2623712, i.to(u.deg).value, places=6)
        self.assertAlmostEqual(30.4818530, W.to(u.deg).value, places=7)
        self.assertAlmostEqual(3.02425566, w.to(u.deg).value % 360, places=8)
        self.assertAlmostEqual(2457277.004, t0.to(u.d).value, places=3)

        # 8.3
        r = [-0.5316809, 0.8283019, 0] * u.au
        r1 = [-0.0147583, -0.0093581, 0] * u.au / u.d
        t_r = Time('2014-01-23T00:00:00', format='isot', scale='tt').jd * u.d
        mu = mu_sun(1 / constants.Earth_Moon)

        a, e, i, W, w, M_r, t0 = solve(r, r1, mu, t_r)
        self.assertAlmostEqual(1.0000185, a.to(u.au).value, places=5)
        self.assertAlmostEqual(0.0169947, e.value, places=3)
        self.assertAlmostEqual(0, i.to(u.deg).value)
        self.assertAlmostEqual(0, W.to(u.deg).value)
        self.assertAlmostEqual(100.3653616, w.to(u.deg).value, delta=3)  # Book error
        self.assertAlmostEqual(21.6299438, M_r.to(u.deg).value, delta=3)  # Book error

        # 8.4
        r = [2.5, 0, 0.1] * u.au
        r1 = [0.006, 0, 0] * u.au / u.d
        t_r = 0 * u.d
        mu = mu_sun(0)

        a, e, i, W, w, M_r, t0 = solve(r, r1, mu, t_r)
        self.assertAlmostEqual(1.4755725, a.to(u.au).value, places=7)
        self.assertAlmostEqual(0.9995876, e.value, places=6)
        self.assertAlmostEqual(90, i.to(u.deg).value)
        self.assertAlmostEqual(180, W.to(u.deg).value)
        self.assertAlmostEqual(358.4061828, w.to(u.deg).value % 360, places=5)
        self.assertAlmostEqual(92.9695608, M_r.to(u.deg).value, places=4)

        from CelestialMechanics.orbital_elements import orbital_elements
        r, angle, r1, r_angle1 = orbital_elements.solve_ellipse(a, e, M_r, mu, t_r, 100 * u.d)
        self.assertAlmostEqual(2.8924662, r.value, places=2)
        self.assertAlmostEqual(179.7667993, angle.to(u.deg).value, places=2)
        self.assertAlmostEqual(0.0020063, r1.value, places=4)
        self.assertAlmostEqual(0.0002075, r_angle1.to(u.au / u.d).value, places=6)
        (x, y, z), (x1, y1, z1) = orbital_elements.solve(a, e, W, w, i, M_r, mu, t_r, 100 * u.d)
        self.assertAlmostEqual(r.value, np.sqrt(x * x + y * y + z * z).value)
        self.assertAlmostEqual(2.8909957, x.value, places=1)
        self.assertAlmostEqual(0, y.value)
        self.assertAlmostEqual(0.0922178, z.value, places=3)
        self.assertAlmostEqual(0.00201190, x1.value, places=4)
        self.assertAlmostEqual(0, y1.value)
        self.assertAlmostEqual(-0.0001434, z1.value, places=5)

    def test_problems(self):
        # 8.1
        r_ = [-2.32791156, -0.80227612, -0.35673637] * u.au
        r_ = ecu2eclip_(r_)
        r1_ = [0.00554700, -0.00883579, -0.00261369] * u.au / u.d
        r1_ = ecu2eclip_(r1_)
        t_r = Time('2015-06-26T00:00:00', format='isot', scale='tt').jd * u.d
        mu = mu_sun(0)

        a, e, i, W, w, M_r, t0 = solve(r_, r1_, mu, t_r)
        self.assertAlmostEqual(2.42152141, a.to(u.au).value, places=5)
        self.assertAlmostEqual(0.18479305, e.value, places=6)
        self.assertAlmostEqual(202.44598740, W.to(u.deg).value, places=2)
        self.assertAlmostEqual(107.13869188, w.to(u.deg).value, places=3)
        self.assertAlmostEqual(6.02979307, i.to(u.deg).value, places=3)
        self.assertAlmostEqual(271.92847594, M_r.to(u.deg).value, places=3)
        # print(Time(t0, format='jd', scale='utc').isot)

        # 8.2
        r_ = [-2.57961310, -1.46709088, -1.23199012] * u.au
        r1_ = [-0.00850280, 0.01015010, 0.00297724] * u.au / u.d
        t_r = Time('2005-08-20T00:00:00', format='isot', scale='tt').jd * u.d
        mu = mu_sun(0)

        a, e, i, W, w, M_r, t0 = solve(r_, r1_, mu, t_r)
        self.assertAlmostEqual(3.19393775, a.to(u.au).value)
        self.assertAlmostEqual(1, e.value)
        self.assertAlmostEqual(155.85899889, W.to(u.deg).value, places=7)
        self.assertAlmostEqual(294.20696215, w.to(u.deg).value, places=7)
        self.assertAlmostEqual(152.76699862, i.to(u.deg).value, places=7)
        self.assertAlmostEqual(2453565.9999, t0.to(u.d).value, places=2)
        # print(Time(t0, format='jd', scale='utc').isot)

    def test_gauss(self):
        t1 = Time('2013-04-10T00:00:00', format='isot', scale='tt').jd * u.d

        r_1 = [2.1233484, -0.7935677, 0.4372058] * u.au
        r_2 = [2.1450400, -0.6739860, 0.4174729] * u.au
        r_3 = [2.1555124, -0.6014181, 0.4051345] * u.au

        a, e, i, W, w, M_r, t0 = solve_gauss(r_1, r_2, r_3, mu_sun(0), t1)
        self.assertAlmostEqual(215.4785322, W.to(u.deg).value, places=3)
        self.assertAlmostEqual(13.1011075, i.to(u.deg).value, places=4)
        self.assertAlmostEqual(180.4021798, w.to(u.deg).value, delta=1)
        self.assertAlmostEqual(0.2476931, e.value, places=2)
        self.assertAlmostEqual(2.7898982, a.to(u.au).value, places=1)
        self.assertAlmostEqual(324.3914010, M_r.to(u.deg).value, delta=1)
        self.assertAlmostEqual(2454858.7869853, t0.to(u.d).value, delta=11)


if __name__ == '__main__':
    unittest.main()
