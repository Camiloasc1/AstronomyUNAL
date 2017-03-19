import numpy as np


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


def angle(q: float, r: float) -> float:
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
    return np.rad2deg(np.arccos(2 * q / r - 1))


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
