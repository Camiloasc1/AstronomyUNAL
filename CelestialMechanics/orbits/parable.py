import numpy as np


def r(q: float, angle: float) -> float:
    """
    Calculate r = (a * (1 - e * e)) / (1 + e * cos(angle))

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


def e() -> float:
    """
    e = 1

    :return: eccentricity
    :rtype: float
    """
    return 1


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
    l = 2*a

    :param a: semi-major axis
    :type a: float
    :return: semi-latus rectum
    :rtype: float
    """
    return 2 * a


def p(a: float) -> float:
    """
    p = b * b / c(a, b)

    :param a: semi-major axis
    :type a: float
    :return: focal parameter
    :rtype: float
    """
    return 2 * a
