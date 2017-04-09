import numpy as np
from astropy import units as u

from CelestialMechanics.orbits import ellipse, parable, hyperbola
from CelestialMechanics import util


def solve(r_, r1_, mu, t):
    h_ = util.cross(r_, r1_)

    h = util.norm(h_)
    r = util.norm(r_)
    v = util.norm(r1_)

    D = D_rvmu(r, v, mu)

    if np.isclose(D.value, 2):  # Parable
        a = h * h / 2 / mu
        e_ = -util.cross(h_, r1_) / mu - r_ / r
        e = 1
    elif D.value < 2:  # Ellipse
        a = r / (2 - D)
        e_ = -util.cross(h_, r1_) / mu - r_ / r
        assert np.isclose(0, np.dot(h_, e_))
        e = np.sqrt(1 - h * h / mu / a)
        assert np.isclose(e, util.norm(e_))
    elif D.value > 2:  # Hyperbola
        a = -r / (2 - D)
        e_ = -util.cross(h_, r1_) / mu - r_ / r
        assert np.isclose(0, np.dot(h_, e_))
        e = np.sqrt(1 + h * h / mu / a)
        assert np.isclose(e, util.norm(e_))

    # angle = util.angle_orbit(a, e, r, r_, r1_)

    if e < 1:  # Ellipse
        angles = ellipse.angles(a, e, r)
        angle = util.angle_cuadrant(angles, r_, r1_)
        E = ellipse.E_angle(e, angle)
        M_r = ellipse.angle_M_eE(e, E)
        t0 = ellipse.t0(M_r, t, ellipse.n(a, mu))
    elif e > 1:  # Hyperbola
        angles = hyperbola.angles(a, e, r)
        angle = util.angle_cuadrant(angles, r_, r1_)
        F = hyperbola.F_angle(e, angle)
        Mh = hyperbola.angle_Mh_eF(e, F)
        t0 = hyperbola.t0(Mh, t, ellipse.n(a, mu) / u.rad)
        M_r = None
    elif e == 1:  # Parable
        angles = parable.angles(a, r)
        angle = util.angle_cuadrant(angles, r_, r1_)
        FI = parable.FI_angle(angle)
        S = parable.S_FI(FI)
        t0 = parable.t0(S, t, ellipse.n(a, mu) / u.rad)
        M_r = None

    if np.isclose(h_[2].value, h.value):  # Earth
        i = 0 * u.deg
        W = 0 * u.deg
        w = np.arctan2(r_[1], r_[0]) - angle
    else:
        i = np.arccos(h_[2] / h)
        W = np.arctan2(h_[0], -h_[1])
        w = np.arctan2(e_[2] * h, -e_[0] * h_[1] + e_[1] * h_[0])
    return a, e, i, W, w, M_r, t0


def D_rvmu(r: float, v: float, mu: float) -> float:
    """
    D = r * v * v / mu

    :param r: radius vector
    :type r: float
    :param v: velocity
    :type v: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :return: D
    :rtype: float
    """
    return r * v * v / mu
