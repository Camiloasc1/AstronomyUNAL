import numpy as np
from typing import Tuple


def r(a: float, e: float, angle: float) -> float:
    """
    r = (a * (1 - e * e)) / (1 + e * cos(angle))

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

    r = (a * (1. - e * e)) / (1. + e * np.cos(angle))

    r = float(r)
    return r


def ae(q: float, Q: float) -> Tuple[float, float]:
    """
    a = (Q + q) / 2
    e = (Q - q) / (2 * a)

    :param q: pericentric distance
    :type q: float
    :param Q: apocentric distance
    :type Q: float
    :return: semi-major axis, eccentricity
    :rtype: float
    """
    a = (Q + q) / 2
    return a, (Q - q) / (2 * a)


def b(a: float, e: float) -> float:
    """
    b = a * a * (1 - e * e)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :return: semi-minor axis
    :rtype: float
    """
    return a * a * (1 - e * e)


def qQ(a: float, e: float) -> Tuple[float, float]:
    """
    q = a * (1 - e)
    Q = a * (1 + e)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :return: pericentric distance, apocentric distance
    :rtype: (float, float)
    """
    return q(a, e), Q(a, e)


def q(a: float, e: float) -> float:
    """
    q = a * (1 - e)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :return: pericentric distance
    :rtype: float
    """
    return a * (1 - e)


def Q(a: float, e: float) -> float:
    """
    Q = a * (1 + e)

    :param a: semi-major axis
    :type a: float
    :param e: eccentricity
    :type e: float
    :return: apocentric distance
    :rtype: float
    """
    return a * (1 + e)


def e(a: float, b: float) -> float:
    """
    e = sqrt(1 - (b * b) / (a * a))

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: eccentricity
    :rtype: float
    """
    return np.sqrt(1 - (b * b) / (a * a))


def c(a: float, b: float) -> float:
    """
    c = sqrt(a * a - b * b)

    :param a: semi-major axis
    :type a: float
    :param b: semi-minor axis
    :type b: float
    :return: linear eccentricity
    :rtype: float
    """
    return np.sqrt(a * a - b * b)


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
    E = -G * m1 * m2 / (2 * a)

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

    return -G * m1 * m2 / (2 * a)
