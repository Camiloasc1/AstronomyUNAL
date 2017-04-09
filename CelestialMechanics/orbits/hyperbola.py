from typing import Tuple

import numpy as np
from astropy import units as u


def r(a: float, e: float, angle: float) -> float:
    """
    Calculate r = (a * (e * e - 1.)) / (1. + e * np.cos(angle))

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param angle: theta angle
    :type angle: float
    :return: radius vector
    :rtype: float
    """
    angle = np.deg2rad(angle)

    r = (a * (e * e - 1.)) / (1. + e * np.cos(angle))

    r = float(r)
    return r


def angles(a: float, e: float, r: float) -> Tuple[float, float]:
    """
    theta = arccos((a * (e * e - 1) - r) / (e * r))

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param r: radius vector
    :type r: float
    :return: the two theta angles
    :rtype: (float, float)
    """
    angle = a * (e * e - 1) - r
    angle /= e * r
    angle = np.arccos(angle)
    angle = angle.to(u.deg)
    return angle, (360. * u.deg) - angle


def a(q: float, e: float) -> float:
    """
    a = q / (e - 1)

    :param q: pericentric distance
    :type q: float
    :param e: eccentricity
    :type e: float
    :return: semi-major axis
    :rtype: float
    """
    return q / (e - 1)


def q(a: float, e: float) -> float:
    """
    q = a * (e - 1)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :return: pericentric distance
    :rtype: float
    """
    return a * (e - 1)


def angle_asymptotic(e: float) -> Tuple[float, float]:
    angle = np.arccos(-1. / e)
    angle = np.rad2deg(angle)

    angle = float(angle)
    return angle, 360. - angle


def e(a: float, b: float) -> float:
    """
    e = sqrt(1 + (b * b) / (a * a))

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: eccentricity
    :rtype: float
    """
    return np.sqrt(1 + (b * b) / (a * a))


def c(a: float, b: float) -> float:
    """
    c = sqrt(a * a + b * b)

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: linear eccentricity
    :rtype: float
    """
    return np.sqrt(a * a + b * b)


def l(a: float, b: float) -> float:
    """
    l = b * b / a

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: semi-latus rectum
    :rtype: float
    """
    return b * b / a


def p(a: float, b: float) -> float:
    """
    p = b * b / c(a, b)

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: focal parameter
    :rtype: float
    """
    return b * b / c(a, b)


def E(a: float, m1: float, m2: float) -> float:
    """
    E = G * m1 * m2 / (2 * a)

    :param a: semi-major axis
    :type a: float
    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: energy of the 2-body system
    :rtype: float
    """
    from astropy.constants import G

    return G * m1 * m2 / (2 * a)


def v(r: float, a: float, m1: float, m2: float) -> float:
    """
    v = sqrt(G * (m1 + m2) * (2 / r + 1 / a))

    :param r: radius vector
    :type r: float
    :param a: semi-major axis
    :type a: float
    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: velocity
    :rtype: float
    """
    from CelestialMechanics.mu import mu_gm1m2

    v = mu_gm1m2(m1, m2) * (2 / r + 1 / a)
    v = np.sqrt(v)

    return v


def Mh(a: float, mu: float, t_r: float, t: float) -> float:
    """
    Mh = sqrt(mu / a ** 3) * (t - t_r)

    :param a: semi-major axis
    :type a: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :param t_r: reference time
    :type t_r: float
    :param t: time
    :type t: float
    :return: Mh angle
    :rtype: float
    """
    Mh = mu / a ** 3
    Mh = np.sqrt(Mh)
    Mh *= t - t_r
    return Mh


def angle_Mh_eF(e: float, F: float) -> float:
    """
    Mh = e * sinh(F) - F

    :param e: eccentricity
    :type e: float
    :param F: eccentric anomaly
    :type F: float
    :return: Mh angle
    :rtype: float
    """
    e = e * u.rad
    return e * np.sinh(F) - F


def solve_F(Mh: float, e: float, ROUNDS: float) -> float:
    Mh = Mh.value
    F = np.arcsinh(1 / e * (2 * Mh))
    for i in range(ROUNDS):
        F = np.arcsinh(1 / e * (Mh + F))
    return F * u.rad


def r_F(a: float, e: float, F: float) -> float:
    """
    r = a * (e * np.cosh(F) - 1)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param F: eccentric anomaly
    :type F: float
    :return: radius vector
    :rtype: float
    """
    return a * (e * np.cosh(F) - 1)


def angle_F(e: float, F: float) -> float:
    """
    theta = 2 * np.arctan(np.sqrt((e + 1) / (e - 1)) * np.tanh(F / 2))

    :param e: eccentricity
    :type e: float
    :param F: eccentric anomaly
    :type F: float
    :return: theta angle
    :rtype: float
    """
    return 2 * np.arctan(np.sqrt((e + 1) / (e - 1)) * np.tanh(F / 2))


def F_angle(e: float, angle: float) -> float:
    """
    F = 2 * np.arctan(np.sqrt((e - 1) / (e + 1)) * np.tanh(angle / 2))

    :param e: eccentricity
    :type e: float
    :param angle: theta angle
    :type angle: float
    :return: eccentric anomaly
    :rtype: float
    """
    return 2 * np.arctan(np.sqrt((e - 1) / (e + 1)) * np.tanh(angle / 2))


def r1(a: float, e: float, angle: float, mu: float) -> float:
    """
    r. = sqrt(mu / a / (e * e - 1)) * e * sin(angle)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param angle: theta angle
    :type angle: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :return: r.
    :rtype: float
    """
    r1 = mu / a / (e * e - 1)
    r1 = np.sqrt(r1)
    r1 *= e * np.sin(angle)
    return r1


def r_angle1(a: float, e: float, r: float, mu: float) -> float:
    """
    rthetha.= sqrt(mu * a * (e * e - 1)) / r

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :param r: radius vector
    :type r: float
    :param mu: G * (m1 + m2)
    :type m1: float
    :return: r theta.
    :rtype: float
    """
    r_angle1 = mu * a * (e * e - 1)
    r_angle1 = np.sqrt(r_angle1)
    r_angle1 /= r
    return r_angle1


def t0(Mh_r: float, t_r: float, n: float) -> float:
    """
    t0 = t_r - M_r / n

    :param Mh_r: mean anomaly at t_r
    :type Mh_r: float
    :param t_r: reference time
    :type t_r: float
    :param n: mean movement
    :type n: float
    :return: t0
    :rtype: float
    """
    return t_r - Mh_r / n
