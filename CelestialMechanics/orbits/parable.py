from typing import Tuple

import numpy as np
from astropy import units as u


def r(q: float, angle: float) -> float:
    """
    r = 2 * q / (1 + np.cos(angle))

    :param q: pericentric distance
    :type q: float
    :param angle: theta angle
    :type angle: float
    :return: radius vector
    :rtype: float
    """
    angle = np.deg2rad(angle)

    r = 2 * q / (1 + np.cos(angle))

    r = float(r)
    return r


def angle(q: float, r: float) -> Tuple[float, float]:
    """
    cos(angle) = 2 * q / r - 1
    angle = arccos(2 * q / r - 1)

    :param q: pericentric distance
    :type q: float
    :param r: radius vector
    :type r: float
    :return: angle
    :rtype: float
    """
    angle = np.rad2deg(np.arccos(2 * q / r - 1))
    return angle, (360. * u.deg) - angle


def e() -> float:
    """
    e = 1

    :return: eccentricity
    :rtype: float
    """
    return 1.


def c(a: float) -> float:
    """
    c = -

    :param a: semi-major axis
    :type a: float
    :return: linear eccentricity
    :rtype: float
    """
    return None


def l(a: float) -> float:
    """
    l = 2 * a

    :param a: semi-major axis
    :type a: float
    :return: semi-latus rectum
    :rtype: float
    """
    return 2. * a


def p(a: float) -> float:
    """
    p = 2 * a

    :param a: semi-major axis
    :type a: float
    :return: focal parameter
    :rtype: float
    """
    return 2. * a


def E() -> float:
    """
    E = 0

    :return: energy of the 2-body system
    :rtype: float
    """

    return 0


def v(r: float, m1: float, m2: float) -> float:
    """
    v = sqrt(2 * G * (m1 + m2) / r)

    :param r: radius vector
    :type r: float
    :param m1: mass 1
    :type m1: float
    :param m2: mass 2
    :type m2: float
    :return: velocity
    :rtype: float
    """
    from CelestialMechanics.mu import mu_gm1m2

    v = 2 * mu_gm1m2(m1, m2) / r
    v = np.sqrt(v)

    return v


def r1(q: float, angle: float, mu: float) -> float:
    """
    r. = sqrt(mu / 2 / q) * e * sin(angle)

    :param q: pericentric distance
    :type q: float
    :param angle: theta angle
    :type angle: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :return: r.
    :rtype: float
    """
    r1 = mu / 2 / q
    r1 = np.sqrt(r1)
    r1 *= np.sin(angle)
    return r1


def r_angle1(q: float, r: float, mu: float) -> float:
    """
    rthetha.= sqrt(2 * q * mu) / r

    :param q: pericentric distance
    :type q: float
    :param r: radius vector
    :type r: float
    :param mu: G * (m1 + m2)
    :type mu: float
    :return: r.
    :rtype: float
    """
    r_angle1 = 2 * q * mu
    r_angle1 = np.sqrt(r_angle1)
    r_angle1 /= r
    return r_angle1
