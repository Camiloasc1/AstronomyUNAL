import unittest

from astropy import units as u
from astropy import constants as astroconst


class MyTestCase(unittest.TestCase):
    def test_2(self):
        from numpy.linalg import norm
        from CelestialMechanics.physics import GravitationalForce
        pos = [-4.11415091, 3.40901224, 0.07790249] * u.au
        force = -GravitationalForce(pos, astroconst.M_sun, astroconst.M_jup).to(u.kg * u.m / u.s ** 2)
        acceleration = force / astroconst.M_jup
        acceleration = acceleration.to(u.au / u.d ** 2)
        self.assertAlmostEqual(norm([0.00000791, -0.00000662, -0.00000015]), norm(acceleration), places=3)
        self.assertAlmostEqual(0.00000791, acceleration[0].value, places=3)
        self.assertAlmostEqual(-0.00000662, acceleration[1].value, places=3)
        self.assertAlmostEqual(-0.00000015, acceleration[2].value, places=3)

    def test_3(self):
        from CelestialMechanics.orbits import parable
        q = 0.06727709 * u.au
        r = 2.33987471 * u.au
        angle1, angle2 = parable.angle(q, r)
        self.assertAlmostEqual(199.5251056, angle2.value, places=7)

    def test_4(self):
        from CelestialMechanics.orbits import ellipse
        v_q = 56.76 * u.km / u.s
        v_Q = 20.67 * u.km / u.s
        a, e = ellipse.ae(v_q, v_Q)
        self.assertAlmostEqual(0.4661, e.value, places=4)
        self.assertAlmostEqual(0.75616, a.value, places=5)

    def test_5(self):
        from astropy.constants import R_sun
        from CelestialMechanics.kepler import kepler3
        T = 25.05 * u.d
        a = kepler3.a_sun(T, 0) - R_sun
        self.assertAlmostEqual(0.163, a.value, places=2)


if __name__ == '__main__':
    unittest.main()
